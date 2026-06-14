#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature18_MGI.py

Feature 18: MGI mouse knockout / phenotype severity features.

Purpose
-------
Build gene-level mouse knockout / mouse phenotype features for human genes.

Main idea:
    human HGNC gene
        -> mouse ortholog from MGI homology reports
        -> mouse phenotype annotations from MGI
        -> lethality / viability / phenotype burden features

Useful features:
    - has mouse ortholog
    - one-to-one mouse/human orthology flag
    - mouse phenotype annotation count
    - unique Mammalian Phenotype term count
    - lethality flag
    - embryonic lethality flag
    - prenatal/perinatal/postnatal lethality flags
    - abnormal survival/viability flag
    - growth/development/reproduction/nervous/immune/metabolism phenotype counts
    - PubMed evidence count

Inputs
------
HGNC:
    databases/HGNC/hgnc_complete_set.txt

MGI downloads:
    feature_databases/MGI/HOM_MouseHumanSequence.rpt
    feature_databases/MGI/HMD_HumanPhenotype.rpt
    feature_databases/MGI/MGI_GenePheno.rpt

Outputs
-------
feature18_mgi/
    downloads/HOM_MouseHumanSequence.rpt
    downloads/HMD_HumanPhenotype.rpt
    downloads/MGI_GenePheno.rpt
    processed/feature18_mgi_human_mouse_orthologs.csv
    processed/feature18_mgi_mouse_gene_phenotype_long.csv
    processed/feature18_mgi_gene_features.csv
    processed/feature18_mgi_hgnc_merged.csv
    feature18_mgi_summary.txt
    feature18_mgi_run_metadata.json

Run
---
    python Feature18_MGI.py --download

Fast test:
    python Feature18_MGI.py --download --max-pheno-rows 200000

Leakage policy
--------------
Safe. Uses mouse orthology and mouse phenotype/knockout phenotype annotations.

Excluded:
    human drug labels
    clinical target labels
    ChEMBL
    DrugBank
    DGIdb
    Open Targets
    Pharos
    known drug-target labels
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import shutil
import sys
import time
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pandas as pd


# =============================================================================
# DEFAULTS
# =============================================================================

DEFAULT_HGNC_PATH = Path("databases") / "HGNC" / "hgnc_complete_set.txt"
DEFAULT_DBDIR = Path("feature_databases") / "MGI"
DEFAULT_OUTDIR = Path("feature18_mgi")

MGI_BASE_URL = "https://www.informatics.jax.org/downloads/reports"

MGI_URLS = {
    "HOM_MouseHumanSequence.rpt": f"{MGI_BASE_URL}/HOM_MouseHumanSequence.rpt",
    "HMD_HumanPhenotype.rpt": f"{MGI_BASE_URL}/HMD_HumanPhenotype.rpt",
    "MGI_GenePheno.rpt": f"{MGI_BASE_URL}/MGI_GenePheno.rpt",
}

HUMAN_TAX_ID = "9606"
MOUSE_TAX_ID = "10090"

# High-level MP term keywords/classes.
# These are broad biological summaries, not drug labels.
PHENOTYPE_KEYWORDS = {
    "lethality_any": [
        "lethality", "lethal", "death", "mortality", "dies", "dead",
        "embryonic lethality", "prenatal lethality", "perinatal lethality",
        "postnatal lethality", "preweaning lethality",
    ],
    "embryonic_lethality": [
        "embryonic lethality", "embryonic lethal", "embryo lethal",
        "lethality during fetal growth", "prenatal lethality",
        "prenatal lethal",
    ],
    "perinatal_lethality": [
        "perinatal lethality", "perinatal lethal", "neonatal lethality",
        "neonatal lethal",
    ],
    "postnatal_lethality": [
        "postnatal lethality", "postnatal lethal", "preweaning lethality",
        "preweaning lethal", "premature death",
    ],
    "abnormal_survival": [
        "abnormal survival", "decreased survival", "reduced survival",
        "increased mortality", "viability", "abnormal viability",
    ],
    "growth_size": [
        "growth", "body size", "body weight", "small size", "increased body size",
        "decreased body weight", "dwarf", "growth retardation",
    ],
    "development": [
        "development", "morphogenesis", "embryogenesis", "organogenesis",
        "developmental", "differentiation",
    ],
    "reproduction_fertility": [
        "fertility", "reproductive", "reproduction", "infertility",
        "sterility", "gametogenesis", "spermatogenesis", "oogenesis",
    ],
    "nervous_behavior": [
        "nervous system", "brain", "neuron", "behavior", "behaviour",
        "seizure", "learning", "memory", "motor", "locomotor",
    ],
    "immune_hematopoietic": [
        "immune", "hematopoietic", "haematopoietic", "blood", "lymphocyte",
        "t cell", "b cell", "macrophage", "inflammation", "inflammatory",
    ],
    "cardiovascular": [
        "cardiovascular", "heart", "cardiac", "blood vessel", "vascular",
        "circulatory",
    ],
    "respiratory": [
        "respiratory", "lung", "breathing", "pulmonary",
    ],
    "metabolism_homeostasis": [
        "metabolism", "homeostasis", "glucose", "insulin", "lipid",
        "cholesterol", "energy balance", "metabolic",
    ],
    "renal_urinary": [
        "renal", "kidney", "urinary", "bladder",
    ],
    "liver_digestive": [
        "liver", "hepatic", "digestive", "intestine", "stomach", "pancreas",
    ],
    "skeletal_muscle": [
        "skeleton", "skeletal", "bone", "cartilage", "muscle",
    ],
    "craniofacial": [
        "craniofacial", "face", "skull", "jaw", "tooth", "teeth",
    ],
    "vision_eye": [
        "eye", "vision", "retina", "lens", "cornea",
    ],
    "hearing_ear": [
        "ear", "hearing", "auditory", "cochlea",
    ],
    "skin_hair": [
        "skin", "hair", "coat", "pigmentation",
    ],
    "tumorigenesis": [
        "tumor", "tumour", "neoplasm", "cancer", "carcinoma",
    ],
}

# MP root/high-level IDs that are common in MGI high-level phenotype report.
# Exact file may only contain IDs, so we can still count broad classes from IDs if present.
MP_LETHALITY_IDS = {
    "MP:0010768",  # mortality/aging
    "MP:0011100",  # preweaning lethality, complete penetrance
    "MP:0011101",  # preweaning lethality, incomplete penetrance
    "MP:0011102",  # prenatal lethality, complete penetrance
    "MP:0011103",  # prenatal lethality, incomplete penetrance
    "MP:0011104",  # perinatal lethality, complete penetrance
    "MP:0011105",  # perinatal lethality, incomplete penetrance
    "MP:0011106",  # postnatal lethality, complete penetrance
    "MP:0011107",  # postnatal lethality, incomplete penetrance
}


# =============================================================================
# HELPERS
# =============================================================================

def log(msg: str) -> None:
    print(msg, flush=True)


def now_iso() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def clean_text(x: Any) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip()
    if s.lower() in {"nan", "none", "na", "null"}:
        return ""
    return s


def normalize_symbol(x: Any) -> str:
    s = clean_text(x).upper()
    s = re.sub(r"\s+", "", s)
    return s


def clean_entrez(x: Any) -> str:
    s = clean_text(x)
    if not s:
        return ""
    return s.split(".")[0]


def clean_mgi_id(x: Any) -> str:
    s = clean_text(x)
    if not s:
        return ""
    if s.startswith("MGI:"):
        return s
    if re.match(r"^\d+$", s):
        return f"MGI:{s}"
    return s


def split_pipe(x: Any) -> List[str]:
    s = clean_text(x)
    if not s:
        return []
    parts = re.split(r"[|,;]+", s)
    out = []
    for p in parts:
        p = clean_text(p)
        if p and p not in out:
            out.append(p)
    return out


def split_mp_ids(x: Any) -> List[str]:
    s = clean_text(x)
    if not s:
        return []
    ids = re.findall(r"MP:\d+", s)
    if ids:
        return sorted(set(ids))
    return [p for p in split_pipe(s) if p.startswith("MP:")]


def split_pubmed_ids(x: Any) -> List[str]:
    s = clean_text(x)
    if not s:
        return []
    out = []
    for p in re.split(r"[|,; ]+", s):
        p = clean_text(p)
        if p and re.match(r"^\d+$", p):
            out.append(p)
    return sorted(set(out))


def shannon_entropy(items: Iterable[str]) -> float:
    items = [clean_text(x) for x in items if clean_text(x)]
    if not items:
        return 0.0
    c = Counter(items)
    n = sum(c.values())
    ent = 0.0
    for _, count in c.items():
        p = count / n
        if p > 0:
            ent -= p * math.log2(p)
    return float(ent)


def download_file(url: str, outpath: Path, force: bool = False) -> None:
    mkdir(outpath.parent)

    if outpath.exists() and outpath.stat().st_size > 0 and not force:
        log(f"[DOWNLOAD SKIP] Already exists: {outpath}")
        return

    if outpath.exists() and force:
        log(f"[REMOVE EXISTING] {outpath}")
        outpath.unlink()

    tmp = outpath.with_suffix(outpath.suffix + ".tmp")

    log(f"[DOWNLOAD] {url}")
    log(f"[TO]       {outpath}")

    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})

    with urllib.request.urlopen(req, timeout=180) as response:
        total = response.length
        downloaded = 0
        last_print = time.time()

        with open(tmp, "wb") as f:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)

                if time.time() - last_print > 5:
                    if total:
                        pct = 100 * downloaded / total
                        log(f"[DOWNLOAD] {downloaded / 1024 / 1024:.1f} MB / {total / 1024 / 1024:.1f} MB ({pct:.1f}%)")
                    else:
                        log(f"[DOWNLOAD] {downloaded / 1024 / 1024:.1f} MB")
                    last_print = time.time()

    with open(tmp, "rb") as f:
        head = f.read(200).lower()

    if b"<html" in head or b"<!doctype html" in head:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"Downloaded HTML/error page instead of report: {url}")

    tmp.rename(outpath)
    log(f"[DOWNLOAD DONE] {outpath} ({outpath.stat().st_size / 1024 / 1024:.2f} MB)")


def find_local_file(dbdir: Path, filename: str) -> Optional[Path]:
    if not dbdir.exists():
        return None

    exact = dbdir / filename
    if exact.exists() and exact.stat().st_size > 0:
        return exact

    hits = list(dbdir.rglob(filename))
    hits = [p for p in hits if p.exists() and p.stat().st_size > 0]
    if hits:
        return sorted(hits, key=lambda p: p.stat().st_size, reverse=True)[0]

    stem = filename.replace(".rpt", "").lower()
    hits = [p for p in dbdir.rglob("*") if p.is_file() and stem in p.name.lower()]
    hits = [p for p in hits if p.stat().st_size > 0]
    if hits:
        return sorted(hits, key=lambda p: p.stat().st_size, reverse=True)[0]

    return None


def read_rpt_no_header(path: Path, names: List[str], max_rows: int = 0) -> pd.DataFrame:
    nrows = max_rows if max_rows and max_rows > 0 else None
    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        header=None,
        names=names,
        low_memory=False,
        nrows=nrows,
        comment="#",
    )
    return df


def read_rpt_auto(path: Path, expected_names: List[str], max_rows: int = 0) -> pd.DataFrame:
    """
    MGI .rpt files generally have no header. This reader checks whether the first
    row looks like a header and falls back to fixed names.
    """
    nrows = max_rows if max_rows and max_rows > 0 else None

    # Try no-header first because MGI reports usually do not include a true header.
    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        header=None,
        low_memory=False,
        nrows=nrows,
        comment="#",
    )

    if df.shape[1] <= len(expected_names):
        df.columns = expected_names[:df.shape[1]]
    else:
        extra = [f"extra_{i}" for i in range(df.shape[1] - len(expected_names))]
        df.columns = expected_names + extra

    return df


# =============================================================================
# HGNC
# =============================================================================

def load_hgnc(hgnc_path: Path, protein_coding_only: bool = True, limit_genes: int = 0) -> pd.DataFrame:
    if not hgnc_path.exists():
        raise FileNotFoundError(f"HGNC file not found: {hgnc_path}")

    log("=" * 100)
    log(f"[READ HGNC] {hgnc_path}")

    hgnc = pd.read_csv(hgnc_path, sep="\t", dtype=str, low_memory=False)
    hgnc.columns = [c.strip() for c in hgnc.columns]

    if "symbol" not in hgnc.columns:
        raise RuntimeError("HGNC file must contain symbol column.")

    before = len(hgnc)

    if protein_coding_only:
        if "locus_group" in hgnc.columns:
            hgnc = hgnc[hgnc["locus_group"].astype(str).str.lower().eq("protein-coding gene")].copy()
        elif "locus_type" in hgnc.columns:
            hgnc = hgnc[hgnc["locus_type"].astype(str).str.lower().str.contains("protein", na=False)].copy()

    hgnc["gene_symbol"] = hgnc["symbol"].map(normalize_symbol)

    for col in ["hgnc_id", "alias_symbol", "prev_symbol", "entrez_id", "ensembl_gene_id", "uniprot_ids", "name"]:
        if col not in hgnc.columns:
            hgnc[col] = ""

    hgnc["entrez_id_clean"] = hgnc["entrez_id"].map(clean_entrez)
    hgnc = hgnc.drop_duplicates(subset=["gene_symbol"], keep="first").copy()

    if limit_genes and limit_genes > 0:
        hgnc = hgnc.head(limit_genes).copy()

    log(f"[HGNC ROWS] {before} -> {len(hgnc)}")
    log(f"[PROTEIN CODING ONLY] {protein_coding_only}")

    return hgnc


def build_hgnc_maps(hgnc: pd.DataFrame) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str]]:
    symbol_map: Dict[str, str] = {}
    entrez_map: Dict[str, str] = {}
    hgnc_id_map: Dict[str, str] = {}

    for _, row in hgnc.iterrows():
        gene = normalize_symbol(row.get("gene_symbol", ""))

        if not gene:
            continue

        symbol_map[gene] = gene

        entrez = clean_entrez(row.get("entrez_id_clean", ""))
        if entrez:
            entrez_map[entrez] = gene

        hgnc_id = clean_text(row.get("hgnc_id", ""))
        if hgnc_id:
            hgnc_id_map[hgnc_id] = gene

        for field in ["alias_symbol", "prev_symbol"]:
            raw = clean_text(row.get(field, ""))
            if not raw:
                continue
            for part in re.split(r"[|,;]+", raw):
                s = normalize_symbol(part)
                if s and s not in symbol_map:
                    symbol_map[s] = gene

    log("=" * 100)
    log("[HGNC MAPS]")
    log(f"[SYMBOLS + ALIASES] {len(symbol_map)}")
    log(f"[ENTREZ IDS]         {len(entrez_map)}")
    log(f"[HGNC IDS]           {len(hgnc_id_map)}")

    return symbol_map, entrez_map, hgnc_id_map


# =============================================================================
# MGI ORTHOLOGY
# =============================================================================

def load_hom_mouse_human_sequence(path: Path, max_rows: int = 0) -> pd.DataFrame:
    names = [
        "homology_class_key",
        "common_organism_name",
        "ncbi_taxon_id",
        "symbol",
        "entrez_gene_id",
        "mouse_mgi_id",
        "hgnc_id",
        "omim_gene_id",
        "genetic_location",
        "genome_coordinates",
        "nucleotide_refseq_ids",
        "protein_refseq_ids",
        "swissprot_ids",
    ]

    df = read_rpt_auto(path, names, max_rows=max_rows)
    for c in df.columns:
        df[c] = df[c].map(clean_text)

    log("=" * 100)
    log(f"[READ HOM_MouseHumanSequence] {path}")
    log(f"[SHAPE] {df.shape}")
    return df


def build_human_mouse_orthologs(
    hom: pd.DataFrame,
    hgnc: pd.DataFrame,
    symbol_map: Dict[str, str],
    entrez_map: Dict[str, str],
    hgnc_id_map: Dict[str, str],
    processed_dir: Path,
) -> pd.DataFrame:
    log("=" * 100)
    log("[BUILD HUMAN ? MOUSE ORTHOLOG MAP]")

    rows = []

    by_class = dict(tuple(hom.groupby("homology_class_key")))

    for class_key, g in by_class.items():
        human_rows = g[g["ncbi_taxon_id"].eq(HUMAN_TAX_ID)].copy()
        mouse_rows = g[g["ncbi_taxon_id"].eq(MOUSE_TAX_ID)].copy()

        if human_rows.empty or mouse_rows.empty:
            continue

        human_genes = []

        for _, hrow in human_rows.iterrows():
            human_gene = ""

            h_symbol = normalize_symbol(hrow.get("symbol", ""))
            h_entrez = clean_entrez(hrow.get("entrez_gene_id", ""))
            h_hgnc_id = clean_text(hrow.get("hgnc_id", ""))

            if h_symbol and h_symbol in symbol_map:
                human_gene = symbol_map[h_symbol]
            elif h_entrez and h_entrez in entrez_map:
                human_gene = entrez_map[h_entrez]
            elif h_hgnc_id and h_hgnc_id in hgnc_id_map:
                human_gene = hgnc_id_map[h_hgnc_id]

            if human_gene:
                human_genes.append(
                    {
                        "human_gene_symbol": human_gene,
                        "human_symbol_raw": h_symbol,
                        "human_entrez_gene_id": h_entrez,
                        "human_hgnc_id": h_hgnc_id,
                    }
                )

        if not human_genes:
            continue

        mouse_genes = []

        for _, mrow in mouse_rows.iterrows():
            mouse_symbol = normalize_symbol(mrow.get("symbol", ""))
            mouse_entrez = clean_entrez(mrow.get("entrez_gene_id", ""))
            mouse_mgi_id = clean_mgi_id(mrow.get("mouse_mgi_id", ""))

            if mouse_symbol or mouse_mgi_id:
                mouse_genes.append(
                    {
                        "mouse_gene_symbol": mouse_symbol,
                        "mouse_entrez_gene_id": mouse_entrez,
                        "mouse_mgi_id": mouse_mgi_id,
                    }
                )

        if not mouse_genes:
            continue

        is_one_to_one = int(len(human_genes) == 1 and len(mouse_genes) == 1)

        for h in human_genes:
            for m in mouse_genes:
                rows.append(
                    {
                        "gene_symbol": h["human_gene_symbol"],
                        "homology_class_key": class_key,
                        "human_symbol_raw": h["human_symbol_raw"],
                        "human_entrez_gene_id": h["human_entrez_gene_id"],
                        "human_hgnc_id": h["human_hgnc_id"],
                        "mouse_gene_symbol": m["mouse_gene_symbol"],
                        "mouse_entrez_gene_id": m["mouse_entrez_gene_id"],
                        "mouse_mgi_id": m["mouse_mgi_id"],
                        "is_one_to_one_ortholog": is_one_to_one,
                        "human_genes_in_homology_class": len(human_genes),
                        "mouse_genes_in_homology_class": len(mouse_genes),
                    }
                )

    orth = pd.DataFrame(rows)

    if not orth.empty:
        orth = orth.drop_duplicates()

    outpath = processed_dir / "feature18_mgi_human_mouse_orthologs.csv"
    orth.to_csv(outpath, index=False)

    log(f"[SAVED ORTHOLOGS] {outpath}")
    log(f"[ORTHOLOG SHAPE] {orth.shape}")
    log(f"[HUMAN GENES WITH MOUSE ORTHOLOG] {orth['gene_symbol'].nunique() if not orth.empty else 0}")

    return orth


# =============================================================================
# HMD HIGH-LEVEL PHENOTYPE
# =============================================================================

def load_hmd_human_phenotype(path: Optional[Path], max_rows: int = 0) -> pd.DataFrame:
    if path is None or not path.exists():
        log("[HMD] Not found. Skipping high-level phenotype report.")
        return pd.DataFrame()

    names = [
        "human_marker_symbol",
        "human_entrez_gene_id",
        "mouse_marker_symbol",
        "mouse_mgi_id",
        "high_level_mp_ids",
    ]

    df = read_rpt_auto(path, names, max_rows=max_rows)
    for c in df.columns:
        df[c] = df[c].map(clean_text)

    log("=" * 100)
    log(f"[READ HMD_HumanPhenotype] {path}")
    log(f"[SHAPE] {df.shape}")

    return df


def map_hmd_to_human_genes(
    hmd: pd.DataFrame,
    symbol_map: Dict[str, str],
    entrez_map: Dict[str, str],
) -> pd.DataFrame:
    if hmd.empty:
        return hmd

    rows = []

    for _, row in hmd.iterrows():
        gene = ""

        hsym = normalize_symbol(row.get("human_marker_symbol", ""))
        hentrez = clean_entrez(row.get("human_entrez_gene_id", ""))

        if hsym and hsym in symbol_map:
            gene = symbol_map[hsym]
        elif hentrez and hentrez in entrez_map:
            gene = entrez_map[hentrez]

        if not gene:
            continue

        mp_ids = split_mp_ids(row.get("high_level_mp_ids", ""))

        rows.append(
            {
                "gene_symbol": gene,
                "mouse_gene_symbol": normalize_symbol(row.get("mouse_marker_symbol", "")),
                "mouse_mgi_id": clean_mgi_id(row.get("mouse_mgi_id", "")),
                "hmd_high_level_mp_ids": "|".join(mp_ids),
                "hmd_high_level_mp_count": len(mp_ids),
                "hmd_has_lethality_mp_id": int(any(mp in MP_LETHALITY_IDS for mp in mp_ids)),
            }
        )

    out = pd.DataFrame(rows)

    if not out.empty:
        out = out.drop_duplicates()

    return out


# =============================================================================
# GENE PHENOTYPE REPORT
# =============================================================================

def load_mgi_gene_pheno(path: Path, max_rows: int = 0) -> pd.DataFrame:
    names = [
        "allelic_composition",
        "allele_symbols",
        "allele_ids",
        "genetic_background",
        "mp_id",
        "pubmed_ids",
        "mouse_mgi_ids",
        "mgi_genotype_id",
    ]

    df = read_rpt_auto(path, names, max_rows=max_rows)

    for c in df.columns:
        df[c] = df[c].map(clean_text)

    log("=" * 100)
    log(f"[READ MGI_GenePheno] {path}")
    log(f"[SHAPE] {df.shape}")

    return df


def classify_phenotype_from_text(mp_id: str, text: str) -> Set[str]:
    t = f"{mp_id} {text}".lower()
    cats = set()

    for cat, kws in PHENOTYPE_KEYWORDS.items():
        if any(kw.lower() in t for kw in kws):
            cats.add(cat)

    if mp_id in MP_LETHALITY_IDS:
        cats.add("lethality_any")

    if not cats and (mp_id or clean_text(text)):
        cats.add("other")

    return cats


def infer_knockout_like(allelic_composition: str, allele_symbols: str) -> int:
    """
    Broad knockout/loss-of-function proxy from allele/genotype text.

    This is intentionally conservative but flexible:
      - targeted mutation
      - null allele
      - knockout/KO
      - tm1/tm2 style targeted alleles
      - deletion
    """
    text = f"{allelic_composition} {allele_symbols}".lower()

    kws = [
        "knockout", "ko", "null", "targeted", "targeted mutation",
        "deletion", "deleted", "tm1", "tm2", "tm3", "em1", "em2",
        "loss-of-function", "loss of function",
    ]

    return int(any(k in text for k in kws))


def build_mouse_gene_phenotype_long(
    gene_pheno: pd.DataFrame,
    orth: pd.DataFrame,
    hmd_mapped: pd.DataFrame,
    processed_dir: Path,
) -> pd.DataFrame:
    log("=" * 100)
    log("[BUILD HUMAN-GENE MGI PHENOTYPE LONG]")

    # mouse MGI ID -> human genes
    mgi_to_human: Dict[str, Set[str]] = defaultdict(set)
    mouse_symbol_to_human: Dict[str, Set[str]] = defaultdict(set)

    if not orth.empty:
        for _, row in orth.iterrows():
            h = normalize_symbol(row.get("gene_symbol", ""))
            mid = clean_mgi_id(row.get("mouse_mgi_id", ""))
            msym = normalize_symbol(row.get("mouse_gene_symbol", ""))

            if h and mid:
                mgi_to_human[mid].add(h)
            if h and msym:
                mouse_symbol_to_human[msym].add(h)

    rows = []
    skipped_unmapped = 0

    for idx, row in gene_pheno.iterrows():
        mouse_mgi_ids = [clean_mgi_id(x) for x in split_pipe(row.get("mouse_mgi_ids", ""))]
        mouse_mgi_ids = [x for x in mouse_mgi_ids if x]

        human_genes = set()

        for mid in mouse_mgi_ids:
            human_genes.update(mgi_to_human.get(mid, set()))

        if not human_genes:
            skipped_unmapped += 1
            continue

        mp_id = clean_text(row.get("mp_id", ""))
        allele_symbols = clean_text(row.get("allele_symbols", ""))
        allelic_composition = clean_text(row.get("allelic_composition", ""))
        genotype_id = clean_text(row.get("mgi_genotype_id", ""))

        pubmed_ids = split_pubmed_ids(row.get("pubmed_ids", ""))
        cats = classify_phenotype_from_text(mp_id, f"{allele_symbols} {allelic_composition}")
        ko_like = infer_knockout_like(allelic_composition, allele_symbols)

        for gene in human_genes:
            rows.append(
                {
                    "gene_symbol": gene,
                    "mouse_mgi_ids": "|".join(mouse_mgi_ids),
                    "mp_id": mp_id,
                    "phenotype_categories": "|".join(sorted(cats)),
                    "allelic_composition": allelic_composition,
                    "allele_symbols": allele_symbols,
                    "mgi_genotype_id": genotype_id,
                    "pubmed_ids": "|".join(pubmed_ids),
                    "pubmed_count": len(pubmed_ids),
                    "is_knockout_like": ko_like,
                    "is_lethality_related": int("lethality_any" in cats),
                    "is_embryonic_lethality_related": int("embryonic_lethality" in cats),
                    "is_perinatal_lethality_related": int("perinatal_lethality" in cats),
                    "is_postnatal_lethality_related": int("postnatal_lethality" in cats),
                    "is_abnormal_survival_related": int("abnormal_survival" in cats),
                }
            )

        if (idx + 1) % 500000 == 0:
            log(f"[MGI PHENO] processed={idx + 1:,} mapped_rows={len(rows):,}")

    long_df = pd.DataFrame(rows)

    if not long_df.empty:
        long_df = long_df.drop_duplicates()

    # Add HMD high-level rows as phenotype-level evidence if available.
    # This gives broad phenotype coverage even if MGI_GenePheno lacks term names.
    if not hmd_mapped.empty:
        hmd_rows = []
        for _, row in hmd_mapped.iterrows():
            mp_ids = split_mp_ids(row.get("hmd_high_level_mp_ids", ""))
            for mp_id in mp_ids:
                cats = classify_phenotype_from_text(mp_id, "")
                hmd_rows.append(
                    {
                        "gene_symbol": normalize_symbol(row.get("gene_symbol", "")),
                        "mouse_mgi_ids": clean_mgi_id(row.get("mouse_mgi_id", "")),
                        "mp_id": mp_id,
                        "phenotype_categories": "|".join(sorted(cats)),
                        "allelic_composition": "",
                        "allele_symbols": "",
                        "mgi_genotype_id": "",
                        "pubmed_ids": "",
                        "pubmed_count": 0,
                        "is_knockout_like": 0,
                        "is_lethality_related": int("lethality_any" in cats or row.get("hmd_has_lethality_mp_id", 0) == 1),
                        "is_embryonic_lethality_related": int("embryonic_lethality" in cats),
                        "is_perinatal_lethality_related": int("perinatal_lethality" in cats),
                        "is_postnatal_lethality_related": int("postnatal_lethality" in cats),
                        "is_abnormal_survival_related": int("abnormal_survival" in cats),
                        "source": "HMD_HumanPhenotype",
                    }
                )

        hmd_df = pd.DataFrame(hmd_rows)
        if not hmd_df.empty:
            if "source" not in long_df.columns:
                long_df["source"] = "MGI_GenePheno"
            long_df = pd.concat([long_df, hmd_df], ignore_index=True, sort=False)
            long_df = long_df.drop_duplicates()

    if "source" not in long_df.columns:
        long_df["source"] = "MGI_GenePheno"

    outpath = processed_dir / "feature18_mgi_mouse_gene_phenotype_long.csv"
    long_df.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED MGI PHENOTYPE LONG]")
    log(f"[PATH] {outpath}")
    log(f"[SHAPE] {long_df.shape}")
    log(f"[GENES] {long_df['gene_symbol'].nunique() if not long_df.empty else 0}")
    log(f"[SKIPPED UNMAPPED GENE PHENO ROWS] {skipped_unmapped}")

    return long_df


# =============================================================================
# AGGREGATION
# =============================================================================

def aggregate_mgi_features(
    hgnc: pd.DataFrame,
    orth: pd.DataFrame,
    pheno_long: pd.DataFrame,
    processed_dir: Path,
) -> pd.DataFrame:
    log("=" * 100)
    log("[AGGREGATE MGI FEATURES PER HUMAN GENE]")

    genes = sorted(hgnc["gene_symbol"].unique().tolist())

    orth_by_gene = dict(tuple(orth.groupby("gene_symbol"))) if not orth.empty else {}
    pheno_by_gene = dict(tuple(pheno_long.groupby("gene_symbol"))) if not pheno_long.empty else {}

    phenotype_categories = sorted(set(PHENOTYPE_KEYWORDS.keys()) | {"other"})

    records = []

    for gene in genes:
        og = orth_by_gene.get(gene, pd.DataFrame())
        pg = pheno_by_gene.get(gene, pd.DataFrame())

        rec = {
            "gene_symbol": gene,
            "feature18_mgi_has_mouse_ortholog": int(not og.empty),
            "feature18_mgi_mouse_ortholog_count": int(og["mouse_gene_symbol"].nunique()) if not og.empty and "mouse_gene_symbol" in og.columns else 0,
            "feature18_mgi_has_one_to_one_mouse_ortholog": int((pd.to_numeric(og["is_one_to_one_ortholog"], errors="coerce").fillna(0) == 1).any()) if not og.empty and "is_one_to_one_ortholog" in og.columns else 0,
            "feature18_mgi_mouse_homology_class_count": int(og["homology_class_key"].nunique()) if not og.empty and "homology_class_key" in og.columns else 0,
        }

        if pg.empty:
            rec.update(
                {
                    "feature18_mgi_has_mouse_phenotype": 0,
                    "feature18_mgi_mouse_phenotype_annotation_count": 0,
                    "feature18_mgi_unique_mp_term_count": 0,
                    "feature18_mgi_unique_genotype_count": 0,
                    "feature18_mgi_unique_pubmed_count": 0,
                    "feature18_mgi_has_knockout_like_phenotype": 0,
                    "feature18_mgi_knockout_like_annotation_count": 0,
                    "feature18_mgi_lethality_flag": 0,
                    "feature18_mgi_embryonic_lethality_flag": 0,
                    "feature18_mgi_perinatal_lethality_flag": 0,
                    "feature18_mgi_postnatal_lethality_flag": 0,
                    "feature18_mgi_abnormal_survival_flag": 0,
                }
            )

            for cat in phenotype_categories:
                rec[f"feature18_mgi_pheno_{cat}_count"] = 0

            records.append(rec)
            continue

        pubmeds = []
        if "pubmed_ids" in pg.columns:
            for val in pg["pubmed_ids"].fillna("").astype(str).tolist():
                pubmeds.extend(split_pubmed_ids(val))

        categories = []
        if "phenotype_categories" in pg.columns:
            for val in pg["phenotype_categories"].fillna("").astype(str).tolist():
                categories.extend([x for x in val.split("|") if x])

        cat_counter = Counter(categories)

        ko_vals = pd.to_numeric(pg.get("is_knockout_like", pd.Series(dtype=float)), errors="coerce").fillna(0)
        leth_vals = pd.to_numeric(pg.get("is_lethality_related", pd.Series(dtype=float)), errors="coerce").fillna(0)
        emb_vals = pd.to_numeric(pg.get("is_embryonic_lethality_related", pd.Series(dtype=float)), errors="coerce").fillna(0)
        peri_vals = pd.to_numeric(pg.get("is_perinatal_lethality_related", pd.Series(dtype=float)), errors="coerce").fillna(0)
        post_vals = pd.to_numeric(pg.get("is_postnatal_lethality_related", pd.Series(dtype=float)), errors="coerce").fillna(0)
        surv_vals = pd.to_numeric(pg.get("is_abnormal_survival_related", pd.Series(dtype=float)), errors="coerce").fillna(0)

        rec.update(
            {
                "feature18_mgi_has_mouse_phenotype": 1,
                "feature18_mgi_mouse_phenotype_annotation_count": int(len(pg)),
                "feature18_mgi_unique_mp_term_count": int(pg["mp_id"].replace("", np.nan).dropna().nunique()) if "mp_id" in pg.columns else 0,
                "feature18_mgi_unique_genotype_count": int(pg["mgi_genotype_id"].replace("", np.nan).dropna().nunique()) if "mgi_genotype_id" in pg.columns else 0,
                "feature18_mgi_unique_pubmed_count": int(pd.Series(pubmeds).nunique()) if pubmeds else 0,
                "feature18_mgi_pubmed_mention_count": int(len(pubmeds)),
                "feature18_mgi_phenotype_category_count": int(len(set(categories))),
                "feature18_mgi_phenotype_category_entropy": shannon_entropy(categories),
                "feature18_mgi_mp_term_entropy": shannon_entropy(pg["mp_id"].tolist()) if "mp_id" in pg.columns else 0.0,

                "feature18_mgi_has_knockout_like_phenotype": int(ko_vals.sum() > 0),
                "feature18_mgi_knockout_like_annotation_count": int(ko_vals.sum()),

                "feature18_mgi_lethality_flag": int(leth_vals.sum() > 0),
                "feature18_mgi_lethality_annotation_count": int(leth_vals.sum()),
                "feature18_mgi_embryonic_lethality_flag": int(emb_vals.sum() > 0),
                "feature18_mgi_embryonic_lethality_annotation_count": int(emb_vals.sum()),
                "feature18_mgi_perinatal_lethality_flag": int(peri_vals.sum() > 0),
                "feature18_mgi_perinatal_lethality_annotation_count": int(peri_vals.sum()),
                "feature18_mgi_postnatal_lethality_flag": int(post_vals.sum() > 0),
                "feature18_mgi_postnatal_lethality_annotation_count": int(post_vals.sum()),
                "feature18_mgi_abnormal_survival_flag": int(surv_vals.sum() > 0),
                "feature18_mgi_abnormal_survival_annotation_count": int(surv_vals.sum()),
            }
        )

        # Knockout-specific lethality signals.
        if "is_knockout_like" in pg.columns and "is_lethality_related" in pg.columns:
            ko_lethal = ((ko_vals == 1) & (leth_vals == 1))
            rec["feature18_mgi_knockout_lethality_flag"] = int(ko_lethal.sum() > 0)
            rec["feature18_mgi_knockout_lethality_annotation_count"] = int(ko_lethal.sum())
        else:
            rec["feature18_mgi_knockout_lethality_flag"] = 0
            rec["feature18_mgi_knockout_lethality_annotation_count"] = 0

        for cat in phenotype_categories:
            rec[f"feature18_mgi_pheno_{cat}_count"] = int(cat_counter.get(cat, 0))

        records.append(rec)

    features = pd.DataFrame(records)

    # Derived transforms.
    for col in [
        "feature18_mgi_mouse_phenotype_annotation_count",
        "feature18_mgi_unique_mp_term_count",
        "feature18_mgi_unique_genotype_count",
        "feature18_mgi_unique_pubmed_count",
        "feature18_mgi_lethality_annotation_count",
        "feature18_mgi_knockout_like_annotation_count",
    ]:
        if col in features.columns:
            features[f"{col}_log1p"] = np.log1p(pd.to_numeric(features[col], errors="coerce").fillna(0))

    if (
        "feature18_mgi_unique_mp_term_count_log1p" in features.columns
        and "feature18_mgi_unique_pubmed_count_log1p" in features.columns
    ):
        features["feature18_mgi_evidence_weighted_phenotype_burden"] = (
            features["feature18_mgi_unique_mp_term_count_log1p"].fillna(0)
            * features["feature18_mgi_unique_pubmed_count_log1p"].fillna(0)
        )

    # Fill counts/flags.
    for c in features.columns:
        if c.startswith("feature18_") and (
            "_count" in c
            or "_flag" in c
            or "_has_" in c
            or c.endswith("_log1p")
            or c.endswith("_burden")
        ):
            if features[c].dtype != object:
                features[c] = features[c].fillna(0)

    outpath = processed_dir / "feature18_mgi_gene_features.csv"
    features.to_csv(outpath, index=False)

    log("[SAVED MGI GENE FEATURES]")
    log(f"[PATH] {outpath}")
    log(f"[SHAPE] {features.shape}")
    log(f"[ORTHOLOG COVERAGE] {features['feature18_mgi_has_mouse_ortholog'].mean():.4f}")
    log(f"[PHENOTYPE COVERAGE] {features['feature18_mgi_has_mouse_phenotype'].mean():.4f}")
    log(f"[LETHALITY FLAG COUNT] {int(features['feature18_mgi_lethality_flag'].sum())}")

    return features


# =============================================================================
# MERGE / SUMMARY
# =============================================================================

def merge_with_hgnc(hgnc: pd.DataFrame, features: pd.DataFrame, processed_dir: Path) -> pd.DataFrame:
    merged = hgnc.merge(features, on="gene_symbol", how="left")

    feature_cols = [c for c in features.columns if c != "gene_symbol"]
    merged["feature18_mgi_has_any_feature"] = merged[feature_cols].notna().any(axis=1).astype(int)
    merged["feature18_mgi_n_nonmissing_features"] = merged[feature_cols].notna().sum(axis=1)

    outpath = processed_dir / "feature18_mgi_hgnc_merged.csv"
    merged.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED HGNC MERGED]")
    log(f"[PATH] {outpath}")
    log(f"[SHAPE] {merged.shape}")

    return merged


def write_summary(
    outdir: Path,
    hgnc: pd.DataFrame,
    orth: pd.DataFrame,
    hmd: pd.DataFrame,
    gene_pheno: pd.DataFrame,
    pheno_long: pd.DataFrame,
    features: pd.DataFrame,
    merged: pd.DataFrame,
    file_paths: Dict[str, Path],
    args: argparse.Namespace,
) -> None:
    path = outdir / "feature18_mgi_summary.txt"

    lines = []
    lines.append("Feature 18: MGI mouse knockout / phenotype severity features")
    lines.append("=" * 80)
    lines.append(f"Created at: {now_iso()}")
    lines.append("")
    lines.append("Inputs:")
    lines.append(f"HGNC: {args.hgnc}")
    for k, p in file_paths.items():
        lines.append(f"{k}: {p}")
    lines.append("")
    lines.append("MGI URLs:")
    for k, u in MGI_URLS.items():
        lines.append(f"{k}: {u}")
    lines.append("")
    lines.append("Rows:")
    lines.append(f"HGNC rows: {hgnc.shape[0]}")
    lines.append(f"Ortholog rows: {orth.shape[0]}")
    lines.append(f"HMD rows: {hmd.shape[0]}")
    lines.append(f"MGI_GenePheno rows: {gene_pheno.shape[0]}")
    lines.append(f"Phenotype long rows: {pheno_long.shape[0]}")
    lines.append(f"Gene feature rows: {features.shape[0]}")
    lines.append(f"Gene feature columns: {features.shape[1]}")
    lines.append(f"Merged rows: {merged.shape[0]}")
    lines.append(f"Merged columns: {merged.shape[1]}")
    lines.append("")
    lines.append("Coverage:")
    for col in [
        "feature18_mgi_has_mouse_ortholog",
        "feature18_mgi_has_one_to_one_mouse_ortholog",
        "feature18_mgi_has_mouse_phenotype",
        "feature18_mgi_has_knockout_like_phenotype",
        "feature18_mgi_lethality_flag",
        "feature18_mgi_embryonic_lethality_flag",
        "feature18_mgi_knockout_lethality_flag",
    ]:
        if col in features.columns:
            lines.append(f"{col}: mean={features[col].fillna(0).mean():.4f}, count={int(features[col].fillna(0).sum())}")
    lines.append("")
    lines.append("Feature summaries:")
    for col in [
        "feature18_mgi_mouse_phenotype_annotation_count",
        "feature18_mgi_unique_mp_term_count",
        "feature18_mgi_unique_genotype_count",
        "feature18_mgi_unique_pubmed_count",
        "feature18_mgi_lethality_annotation_count",
        "feature18_mgi_knockout_like_annotation_count",
        "feature18_mgi_evidence_weighted_phenotype_burden",
    ]:
        if col in features.columns:
            x = pd.to_numeric(features[col], errors="coerce")
            lines.append(
                f"{col}: median={x.median(skipna=True):.4f}, "
                f"max={x.max(skipna=True):.4f}, "
                f"nonmissing={int(x.notna().sum())}"
            )
    lines.append("")
    lines.append("Interpretation:")
    lines.append("These features capture mouse orthology and mouse phenotype/knockout severity evidence.")
    lines.append("Lethality flags indicate mouse genotype phenotype annotations related to lethality/survival.")
    lines.append("They are not direct human clinical druggability labels.")
    lines.append("")
    lines.append("Leakage policy:")
    lines.append("Included: MGI orthology, Mammalian Phenotype annotations, PubMed count, knockout-like allele text flags.")
    lines.append("Excluded: clinical target labels, known drug labels, ChEMBL, DrugBank, DGIdb, Open Targets, Pharos.")

    path.write_text("\n".join(lines) + "\n")
    log(f"[SUMMARY] {path}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description="Build Feature 18 MGI mouse knockout / phenotype features.")
    parser.add_argument("--hgnc", default=str(DEFAULT_HGNC_PATH), help="Path to HGNC complete set.")
    parser.add_argument("--dbdir", default=str(DEFAULT_DBDIR), help="MGI database directory.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory.")
    parser.add_argument("--download", action="store_true", help="Download MGI reports.")
    parser.add_argument("--force-download", action="store_true", help="Force re-download.")
    parser.add_argument("--hom-file", default="", help="Manual HOM_MouseHumanSequence.rpt file.")
    parser.add_argument("--hmd-file", default="", help="Manual HMD_HumanPhenotype.rpt file.")
    parser.add_argument("--gene-pheno-file", default="", help="Manual MGI_GenePheno.rpt file.")
    parser.add_argument("--max-hom-rows", type=int, default=0, help="Debug: first N homology rows.")
    parser.add_argument("--max-hmd-rows", type=int, default=0, help="Debug: first N HMD rows.")
    parser.add_argument("--max-pheno-rows", type=int, default=0, help="Debug: first N gene phenotype rows.")
    parser.add_argument("--limit-genes", type=int, default=0, help="Debug: first N HGNC genes.")
    parser.add_argument("--all-hgnc-genes", action="store_true", help="Use all HGNC genes, not only protein-coding.")
    args = parser.parse_args()

    dbdir = mkdir(Path(args.dbdir))
    outdir = mkdir(Path(args.outdir))
    downloads_dir = mkdir(outdir / "downloads")
    processed_dir = mkdir(outdir / "processed")

    log("=" * 100)
    log("FEATURE 18: MGI MOUSE KNOCKOUT / PHENOTYPE SEVERITY")
    log("=" * 100)
    log(f"[HGNC]           {args.hgnc}")
    log(f"[DBDIR]          {dbdir.resolve()}")
    log(f"[OUTDIR]         {outdir.resolve()}")
    log(f"[DOWNLOAD]       {args.download}")
    log(f"[FORCE DOWNLOAD] {args.force_download}")
    log(f"[LIMIT GENES]    {args.limit_genes if args.limit_genes else 'none'}")
    log("=" * 100)

    if args.download or args.force_download:
        for fname, url in MGI_URLS.items():
            download_file(url, dbdir / fname, force=args.force_download)

    hom_file = Path(args.hom_file) if args.hom_file else find_local_file(dbdir, "HOM_MouseHumanSequence.rpt")
    hmd_file = Path(args.hmd_file) if args.hmd_file else find_local_file(dbdir, "HMD_HumanPhenotype.rpt")
    gene_pheno_file = Path(args.gene_pheno_file) if args.gene_pheno_file else find_local_file(dbdir, "MGI_GenePheno.rpt")

    if hom_file is None or not hom_file.exists():
        raise FileNotFoundError("HOM_MouseHumanSequence.rpt not found. Run with --download or provide --hom-file.")

    if gene_pheno_file is None or not gene_pheno_file.exists():
        raise FileNotFoundError("MGI_GenePheno.rpt not found. Run with --download or provide --gene-pheno-file.")

    file_paths = {
        "HOM_MouseHumanSequence": hom_file,
        "HMD_HumanPhenotype": hmd_file if hmd_file else Path(""),
        "MGI_GenePheno": gene_pheno_file,
    }

    # Copy downloaded inputs into output folder for reproducibility.
    for p in [hom_file, hmd_file, gene_pheno_file]:
        if p and p.exists():
            try:
                shutil.copy2(p, downloads_dir / p.name)
            except Exception:
                pass

    hgnc = load_hgnc(
        Path(args.hgnc),
        protein_coding_only=not args.all_hgnc_genes,
        limit_genes=args.limit_genes,
    )

    symbol_map, entrez_map, hgnc_id_map = build_hgnc_maps(hgnc)

    hom = load_hom_mouse_human_sequence(
        hom_file,
        max_rows=args.max_hom_rows,
    )

    orth = build_human_mouse_orthologs(
        hom=hom,
        hgnc=hgnc,
        symbol_map=symbol_map,
        entrez_map=entrez_map,
        hgnc_id_map=hgnc_id_map,
        processed_dir=processed_dir,
    )

    hmd = load_hmd_human_phenotype(
        hmd_file,
        max_rows=args.max_hmd_rows,
    )

    hmd_mapped = map_hmd_to_human_genes(
        hmd=hmd,
        symbol_map=symbol_map,
        entrez_map=entrez_map,
    )

    if not hmd_mapped.empty:
        hmd_mapped_path = processed_dir / "feature18_mgi_hmd_human_phenotype_mapped.csv"
        hmd_mapped.to_csv(hmd_mapped_path, index=False)
        log(f"[SAVED HMD MAPPED] {hmd_mapped_path}")
        log(f"[HMD MAPPED SHAPE] {hmd_mapped.shape}")

    gene_pheno = load_mgi_gene_pheno(
        gene_pheno_file,
        max_rows=args.max_pheno_rows,
    )

    pheno_long = build_mouse_gene_phenotype_long(
        gene_pheno=gene_pheno,
        orth=orth,
        hmd_mapped=hmd_mapped,
        processed_dir=processed_dir,
    )

    features = aggregate_mgi_features(
        hgnc=hgnc,
        orth=orth,
        pheno_long=pheno_long,
        processed_dir=processed_dir,
    )

    merged = merge_with_hgnc(
        hgnc=hgnc,
        features=features,
        processed_dir=processed_dir,
    )

    metadata = {
        "created_at": now_iso(),
        "hgnc": str(Path(args.hgnc)),
        "dbdir": str(dbdir.resolve()),
        "outdir": str(outdir.resolve()),
        "hom_file": str(hom_file),
        "hmd_file": str(hmd_file) if hmd_file else "",
        "gene_pheno_file": str(gene_pheno_file),
        "mgi_urls": MGI_URLS,
        "protein_coding_only": not args.all_hgnc_genes,
        "limit_genes": args.limit_genes,
        "max_hom_rows": args.max_hom_rows,
        "max_hmd_rows": args.max_hmd_rows,
        "max_pheno_rows": args.max_pheno_rows,
        "outputs": {
            "orthologs": str(processed_dir / "feature18_mgi_human_mouse_orthologs.csv"),
            "phenotype_long": str(processed_dir / "feature18_mgi_mouse_gene_phenotype_long.csv"),
            "gene_features": str(processed_dir / "feature18_mgi_gene_features.csv"),
            "hgnc_merged": str(processed_dir / "feature18_mgi_hgnc_merged.csv"),
            "summary": str(outdir / "feature18_mgi_summary.txt"),
        },
        "leakage_policy": {
            "included": [
                "MGI human-mouse orthology",
                "MGI Mammalian Phenotype annotations",
                "mouse genotype phenotype burden",
                "mouse lethality/survival flags",
                "knockout-like allele text flags",
                "PubMed evidence count",
            ],
            "excluded": [
                "clinical target labels",
                "known drug-target labels",
                "ChEMBL",
                "DrugBank",
                "DGIdb",
                "Open Targets",
                "Pharos",
            ],
        },
    }

    with open(outdir / "feature18_mgi_run_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    write_summary(
        outdir=outdir,
        hgnc=hgnc,
        orth=orth,
        hmd=hmd,
        gene_pheno=gene_pheno,
        pheno_long=pheno_long,
        features=features,
        merged=merged,
        file_paths=file_paths,
        args=args,
    )

    log("=" * 100)
    log("[DONE]")
    log(f"[ORTHOLOGS]     {processed_dir / 'feature18_mgi_human_mouse_orthologs.csv'}")
    log(f"[PHENO LONG]    {processed_dir / 'feature18_mgi_mouse_gene_phenotype_long.csv'}")
    log(f"[GENE FEATURES] {processed_dir / 'feature18_mgi_gene_features.csv'}")
    log(f"[HGNC MERGED]   {processed_dir / 'feature18_mgi_hgnc_merged.csv'}")
    log("=" * 100)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        log("\n[STOPPED] Interrupted by user.")
        sys.exit(130)
    except Exception as exc:
        log("=" * 100)
        log("[ERROR]")
        log(str(exc))
        log("=" * 100)
        raise