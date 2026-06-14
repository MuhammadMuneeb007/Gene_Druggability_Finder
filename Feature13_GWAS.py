#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature13_GWAS.py

Feature 13: GWAS Catalog gene-association burden / pleiotropy features.

Purpose
-------
Build leakage-aware, gene-level GWAS Catalog features.

This feature captures whether a gene has been repeatedly implicated near
published GWAS associations across human traits.

Important interpretation
------------------------
This is NOT mutation burden.
This is NOT causal proof.
This is NOT druggability evidence.

It is a human genetic-association burden / pleiotropy feature.

Inputs
------
HGNC:
    databases/HGNC/hgnc_complete_set.txt

GWAS Catalog associations:
    feature_databases/GWASCatalog/gwas_catalog_associations.tsv

or auto-downloaded from:
    https://www.ebi.ac.uk/gwas/api/search/downloads/associations/v1.0?split=false

Outputs
-------
feature13_gwas_catalog/
    downloads/gwas_catalog_associations.tsv
    processed/feature13_gwas_gene_association_long.csv
    processed/feature13_gwas_catalog_gene_features.csv
    processed/feature13_gwas_catalog_hgnc_merged.csv
    feature13_gwas_catalog_summary.txt
    feature13_gwas_catalog_run_metadata.json

Run
---
    python Feature13_GWAS.py

Force download:
    python Feature13_GWAS.py --download --force-download

Fast test:
    python Feature13_GWAS.py --max-rows 100000

Main fix
--------
GWAS Catalog download can be:
    1. plain TSV
    2. gzip-compressed TSV
    3. ZIP archive containing TSV, even if saved as .tsv

This script detects file type using magic bytes instead of file extension.
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
import zipfile
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set

import numpy as np
import pandas as pd


# =============================================================================
# DEFAULTS
# =============================================================================

DEFAULT_HGNC_PATH = Path("databases") / "HGNC" / "hgnc_complete_set.txt"
DEFAULT_DBDIR = Path("feature_databases") / "GWASCatalog"
DEFAULT_OUTDIR = Path("feature13_gwas_catalog")

GWAS_ASSOCIATIONS_URL = "https://www.ebi.ac.uk/gwas/api/search/downloads/associations/v1.0?split=false"


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


def clean_ensembl_gene_id(x: Any) -> str:
    s = clean_text(x)

    if not s:
        return ""

    return s.replace("gene:", "").split(".")[0]


def safe_float(x: Any) -> float:
    try:
        if x is None or pd.isna(x):
            return np.nan

        s = str(x).strip()

        if s in {"", ".", "NA", "NaN", "nan", "None"}:
            return np.nan

        return float(s)

    except Exception:
        return np.nan


def neglog10_pvalue(p: Any) -> float:
    p = safe_float(p)

    if pd.isna(p) or p <= 0:
        return np.nan

    return float(-math.log10(p))


def shannon_entropy(items: Iterable[str]) -> float:
    items = [clean_text(x) for x in items if clean_text(x)]

    if not items:
        return 0.0

    c = Counter(items)
    n = sum(c.values())

    ent = 0.0

    for count in c.values():
        p = count / n
        ent -= p * math.log2(p)

    return float(ent)


def is_gzip_file(path: Path) -> bool:
    """
    Detect gzip by magic bytes.

    gzip starts with:
        1f 8b
    """
    try:
        with open(path, "rb") as f:
            return f.read(2) == b"\x1f\x8b"

    except Exception:
        return False


def is_zip_file(path: Path) -> bool:
    """
    Detect ZIP by magic bytes.

    ZIP usually starts with:
        PK
    """
    try:
        with open(path, "rb") as f:
            return f.read(2) == b"PK"

    except Exception:
        return False


def inspect_archive_or_text(path: Path) -> None:
    """
    Inspect the downloaded GWAS Catalog file.

    Handles:
        - plain TSV
        - gzip TSV
        - ZIP archive containing TSV
    """
    if not path.exists() or path.stat().st_size == 0:
        raise RuntimeError(f"File is missing or empty: {path}")

    log(f"[FILE SIZE] {path.stat().st_size / 1024 / 1024:.2f} MB")

    if is_zip_file(path):
        log("[FILE TYPE] ZIP archive detected by magic bytes")

        with zipfile.ZipFile(path, "r") as z:
            names = z.namelist()
            log(f"[ZIP CONTENTS] {names[:10]}")

            tsv_names = [
                n for n in names
                if n.lower().endswith((".tsv", ".txt", ".csv"))
            ]

            if not tsv_names:
                raise RuntimeError(f"ZIP file does not contain a TSV/TXT/CSV file: {path}")

            inner = tsv_names[0]

            with z.open(inner) as f:
                header = f.readline().decode("utf-8", errors="replace").strip()

            log(f"[ZIP INNER FILE] {inner}")
            log(f"[HEADER PREVIEW] {header[:200]}")

        return

    if is_gzip_file(path):
        log("[FILE TYPE] gzip-compressed content detected by magic bytes")

        with gzip.open(path, "rt", errors="replace") as f:
            header = f.readline().strip()

        log(f"[HEADER PREVIEW] {header[:200]}")

        return

    log("[FILE TYPE] plain text content detected")

    with open(path, "rt", errors="replace") as f:
        header = f.readline().strip()

    log(f"[HEADER PREVIEW] {header[:200]}")


def download_file(url: str, outpath: Path, force: bool = False) -> None:
    mkdir(outpath.parent)

    if outpath.exists() and outpath.stat().st_size > 0 and not force:
        log(f"[DOWNLOAD SKIP] Already exists: {outpath}")
        inspect_archive_or_text(outpath)
        return

    if outpath.exists() and force:
        log(f"[REMOVE EXISTING] {outpath}")
        outpath.unlink()

    tmp = outpath.with_suffix(outpath.suffix + ".tmp")

    log(f"[DOWNLOAD] {url}")
    log(f"[TO]       {outpath}")

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Mozilla/5.0 Feature13_GWASCatalog",
            "Accept": "*/*",
        },
    )

    with urllib.request.urlopen(request) as response:
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
                        log(
                            f"[DOWNLOAD] {downloaded / 1024 / 1024:.1f} MB / "
                            f"{total / 1024 / 1024:.1f} MB ({pct:.1f}%)"
                        )
                    else:
                        log(f"[DOWNLOAD] {downloaded / 1024 / 1024:.1f} MB")

                    last_print = time.time()

    tmp.rename(outpath)

    log(f"[DOWNLOAD DONE] {outpath} ({outpath.stat().st_size / 1024 / 1024:.2f} MB)")
    inspect_archive_or_text(outpath)


def find_local_gwas_file(dbdir: Path) -> Optional[Path]:
    if not dbdir.exists():
        return None

    candidates = []

    for pattern in [
        "*associations*.tsv",
        "*associations*.tsv.gz",
        "*associations*.zip",
        "*gwas*catalog*.tsv",
        "*gwas*catalog*.tsv.gz",
        "*gwas*catalog*.zip",
        "*.tsv",
        "*.tsv.gz",
        "*.zip",
    ]:
        candidates.extend(list(dbdir.rglob(pattern)))

    candidates = [p for p in candidates if p.exists() and p.stat().st_size > 0]

    if not candidates:
        return None

    ranked = []

    for p in candidates:
        name = p.name.lower()
        score = 0

        if "association" in name:
            score += 20
        if "gwas" in name:
            score += 10
        if "catalog" in name:
            score += 10
        if name.endswith(".tsv") or name.endswith(".tsv.gz") or name.endswith(".zip"):
            score += 5
        if is_zip_file(p):
            score += 4
        if is_gzip_file(p):
            score += 3

        ranked.append((score, p.stat().st_size, p))

    ranked.sort(reverse=True)

    return ranked[0][2]


def find_col(
    df: pd.DataFrame,
    candidates: List[str],
    contains: Optional[List[str]] = None,
) -> Optional[str]:
    lower_map = {str(c).strip().lower(): c for c in df.columns}

    for c in candidates:
        key = c.lower().strip()

        if key in lower_map:
            return lower_map[key]

    if contains:
        for c in df.columns:
            low = str(c).lower()

            if all(x.lower() in low for x in contains):
                return c

    return None


# =============================================================================
# HGNC
# =============================================================================

def load_hgnc(
    hgnc_path: Path,
    protein_coding_only: bool = True,
    limit_genes: int = 0,
) -> pd.DataFrame:
    if not hgnc_path.exists():
        raise FileNotFoundError(f"HGNC file not found: {hgnc_path}")

    log("=" * 100)
    log(f"[READ HGNC] {hgnc_path}")

    hgnc = pd.read_csv(hgnc_path, sep="\t", dtype=str, low_memory=False)
    hgnc.columns = [c.strip() for c in hgnc.columns]

    if "symbol" not in hgnc.columns:
        raise RuntimeError("HGNC file must contain a symbol column.")

    before = len(hgnc)

    if protein_coding_only:
        if "locus_group" in hgnc.columns:
            hgnc = hgnc[
                hgnc["locus_group"].astype(str).str.lower().eq("protein-coding gene")
            ].copy()

        elif "locus_type" in hgnc.columns:
            hgnc = hgnc[
                hgnc["locus_type"].astype(str).str.lower().str.contains("protein", na=False)
            ].copy()

    hgnc["gene_symbol"] = hgnc["symbol"].map(normalize_symbol)

    for col in ["alias_symbol", "prev_symbol", "ensembl_gene_id", "entrez_id", "name"]:
        if col not in hgnc.columns:
            hgnc[col] = ""

    hgnc["ensembl_gene_id_clean"] = hgnc["ensembl_gene_id"].map(clean_ensembl_gene_id)
    hgnc = hgnc.drop_duplicates(subset=["gene_symbol"], keep="first").copy()

    if limit_genes and limit_genes > 0:
        hgnc = hgnc.head(limit_genes).copy()

    log(f"[HGNC ROWS] {before} -> {len(hgnc)}")
    log(f"[PROTEIN CODING ONLY] {protein_coding_only}")

    return hgnc


def build_symbol_alias_map(hgnc: pd.DataFrame) -> Dict[str, str]:
    """
    Map current symbols, aliases, and previous symbols back to current HGNC symbol.
    """
    symbol_map: Dict[str, str] = {}

    for _, row in hgnc.iterrows():
        current = normalize_symbol(row.get("gene_symbol", ""))

        if not current:
            continue

        symbol_map[current] = current

        for field in ["alias_symbol", "prev_symbol"]:
            raw = clean_text(row.get(field, ""))

            if not raw:
                continue

            for part in re.split(r"[|,;]+", raw):
                s = normalize_symbol(part)

                if s and s not in symbol_map:
                    symbol_map[s] = current

    log("=" * 100)
    log("[HGNC SYMBOL/ALIAS MAP]")
    log(f"[CURRENT GENES] {hgnc['gene_symbol'].nunique()}")
    log(f"[TOTAL SYMBOLS/ALIASES] {len(symbol_map)}")

    return symbol_map


# =============================================================================
# GWAS GENE PARSING
# =============================================================================

def split_gene_field(x: Any) -> List[str]:
    """
    Parse GWAS Catalog REPORTED GENE(S) / MAPPED_GENE fields.
    """
    s = clean_text(x)

    if not s:
        return []

    bad = {
        "NR",
        "N/R",
        "NA",
        "NONE",
        "NULL",
        "INTERGENIC",
        "NOT REPORTED",
        "NOT MAPPED",
        "UNKNOWN",
    }

    if s.upper() in bad:
        return []

    s = re.sub(r"\([^)]*\)", " ", s)

    s = s.replace(";", ",")
    s = s.replace("|", ",")
    s = s.replace("/", ",")
    s = s.replace(" and ", ",")
    s = s.replace(" AND ", ",")
    s = re.sub(r"\s+-\s+", ",", s)

    parts = re.split(r"[,\s]+", s)

    out = []

    for p in parts:
        p = normalize_symbol(p)

        if not p:
            continue

        if p.upper() in bad:
            continue

        if p in {"UPSTREAM", "DOWNSTREAM", "NEAR", "BETWEEN", "GENE", "GENES"}:
            continue

        if not re.match(r"^[A-Z0-9][A-Z0-9\.-]*$", p):
            continue

        out.append(p)

    seen = set()
    final = []

    for x in out:
        if x not in seen:
            final.append(x)
            seen.add(x)

    return final


def map_gene_symbols(symbols: List[str], symbol_map: Dict[str, str]) -> List[str]:
    mapped = []

    for s in symbols:
        ss = normalize_symbol(s)

        if ss in symbol_map:
            mapped.append(symbol_map[ss])

    seen = set()
    out = []

    for x in mapped:
        if x not in seen:
            out.append(x)
            seen.add(x)

    return out


def classify_trait(text: str) -> Set[str]:
    """
    Simple keyword-based trait category assignment.
    This is not used as causal evidence; it is a broad pleiotropy summary.
    """
    t = clean_text(text).lower()

    cats = set()

    patterns = {
        "anthropometric": [
            "height",
            "body mass",
            "bmi",
            "obesity",
            "weight",
            "waist",
            "hip",
            "anthropometric",
            "lean mass",
            "fat mass",
        ],
        "metabolic": [
            "diabetes",
            "insulin",
            "glucose",
            "metabolic",
            "lipid",
            "cholesterol",
            "triglyceride",
            "hdl",
            "ldl",
            "adiponectin",
        ],
        "cardiovascular": [
            "blood pressure",
            "hypertension",
            "coronary",
            "cardiovascular",
            "heart",
            "atrial",
            "stroke",
            "myocardial",
            "artery",
        ],
        "immune_inflammatory": [
            "immune",
            "autoimmune",
            "inflammatory",
            "asthma",
            "allergy",
            "arthritis",
            "lupus",
            "crohn",
            "colitis",
            "psoriasis",
            "eczema",
        ],
        "neurological_psychiatric": [
            "brain",
            "neuro",
            "alzheimer",
            "parkinson",
            "migraine",
            "epilepsy",
            "schizophrenia",
            "depression",
            "bipolar",
            "autism",
            "cognitive",
            "intelligence",
        ],
        "cancer": [
            "cancer",
            "carcinoma",
            "tumor",
            "tumour",
            "melanoma",
            "leukemia",
            "lymphoma",
            "glioma",
            "neoplasm",
        ],
        "hematological": [
            "blood cell",
            "hemoglobin",
            "haemoglobin",
            "platelet",
            "erythrocyte",
            "leukocyte",
            "white blood",
            "red blood",
            "hematological",
            "haematological",
        ],
        "renal_liver": [
            "kidney",
            "renal",
            "creatinine",
            "liver",
            "hepatic",
            "bilirubin",
            "alanine aminotransferase",
            "aspartate aminotransferase",
        ],
        "reproductive": [
            "menarche",
            "menopause",
            "fertility",
            "reproductive",
            "pregnancy",
            "birth weight",
            "testosterone",
        ],
        "biomarker": [
            "protein level",
            "serum",
            "plasma",
            "biomarker",
            "metabolite",
            "metabolites",
            "enzyme level",
            "hormone",
            "cytokine",
        ],
        "drug_response": [
            "drug response",
            "treatment response",
            "adverse drug",
            "pharmacogen",
            "medication",
            "warfarin",
            "statin",
            "chemotherapy",
        ],
        "infectious": [
            "infection",
            "viral",
            "bacterial",
            "covid",
            "hiv",
            "tuberculosis",
            "malaria",
            "hepatitis",
        ],
    }

    for cat, kws in patterns.items():
        if any(k in t for k in kws):
            cats.add(cat)

    if not cats and t:
        cats.add("other")

    return cats


# =============================================================================
# GWAS READING
# =============================================================================

def read_gwas_associations(path: Path, max_rows: int = 0) -> pd.DataFrame:
    """
    Robust GWAS Catalog reader.

    Handles:
        1. plain TSV
        2. gzip-compressed TSV saved as .tsv or .tsv.gz
        3. ZIP archive saved as .tsv containing the real TSV
    """
    if not path.exists():
        raise FileNotFoundError(f"GWAS association file not found: {path}")

    log("=" * 100)
    log(f"[READ GWAS CATALOG ASSOCIATIONS] {path}")

    nrows = max_rows if max_rows and max_rows > 0 else None

    inspect_archive_or_text(path)

    if is_zip_file(path):
        log("[PANDAS INPUT] reading TSV from inside ZIP archive")

        with zipfile.ZipFile(path, "r") as z:
            names = z.namelist()

            tsv_names = [
                n for n in names
                if n.lower().endswith((".tsv", ".txt", ".csv"))
            ]

            if not tsv_names:
                raise RuntimeError(f"No TSV/TXT/CSV file found inside ZIP archive: {path}")

            inner = tsv_names[0]
            log(f"[ZIP READ INNER FILE] {inner}")

            with z.open(inner) as f:
                try:
                    df = pd.read_csv(
                        f,
                        sep="\t",
                        dtype=str,
                        low_memory=False,
                        nrows=nrows,
                        encoding="utf-8",
                        encoding_errors="replace",
                    )
                except TypeError:
                    df = pd.read_csv(
                        f,
                        sep="\t",
                        dtype=str,
                        low_memory=False,
                        nrows=nrows,
                        encoding="utf-8",
                    )

    elif is_gzip_file(path):
        log("[PANDAS COMPRESSION] gzip")

        try:
            df = pd.read_csv(
                path,
                sep="\t",
                dtype=str,
                low_memory=False,
                nrows=nrows,
                compression="gzip",
                encoding="utf-8",
                encoding_errors="replace",
            )
        except TypeError:
            df = pd.read_csv(
                path,
                sep="\t",
                dtype=str,
                low_memory=False,
                nrows=nrows,
                compression="gzip",
                encoding="utf-8",
            )

    else:
        log("[PANDAS COMPRESSION] none")

        try:
            df = pd.read_csv(
                path,
                sep="\t",
                dtype=str,
                low_memory=False,
                nrows=nrows,
                compression=None,
                encoding="utf-8",
                encoding_errors="replace",
            )
        except TypeError:
            df = pd.read_csv(
                path,
                sep="\t",
                dtype=str,
                low_memory=False,
                nrows=nrows,
                compression=None,
                encoding="utf-8",
            )

    df.columns = [str(c).strip() for c in df.columns]

    log(f"[GWAS SHAPE] {df.shape}")
    log(f"[GWAS COLUMNS] {list(df.columns)[:20]} ...")

    return df


# =============================================================================
# BUILD LONG TABLE
# =============================================================================

def build_gwas_gene_long(
    gwas: pd.DataFrame,
    symbol_map: Dict[str, str],
    processed_dir: Path,
) -> pd.DataFrame:
    log("=" * 100)
    log("[BUILD GWAS GENE LONG TABLE]")

    reported_gene_col = find_col(
        gwas,
        ["REPORTED GENE(S)", "Reported Gene(s)", "reported genes", "reported_gene"],
        contains=["reported", "gene"],
    )

    mapped_gene_col = find_col(
        gwas,
        ["MAPPED_GENE", "MAPPED GENE", "MAPPED GENE(S)", "mapped_gene", "mapped genes"],
        contains=["mapped", "gene"],
    )

    trait_col = find_col(
        gwas,
        ["DISEASE/TRAIT", "Disease/Trait", "disease_trait", "trait"],
        contains=["trait"],
    )

    efo_col = find_col(
        gwas,
        ["MAPPED_TRAIT", "MAPPED TRAIT", "mapped_trait", "EFO trait", "MAPPED_TRAIT_URI"],
        contains=["mapped", "trait"],
    )

    pvalue_col = find_col(
        gwas,
        ["P-VALUE", "P VALUE", "pvalue", "p_value", "P"],
        contains=["p"],
    )

    snp_col = find_col(
        gwas,
        ["SNPS", "SNP_ID_CURRENT", "snp", "variant_id"],
        contains=["snp"],
    )

    study_col = find_col(
        gwas,
        ["STUDY ACCESSION", "STUDY_ACCESSION", "study_accession"],
        contains=["study", "accession"],
    )

    pubmed_col = find_col(
        gwas,
        ["PUBMEDID", "PUBMED ID", "pubmedid", "pmid"],
        contains=["pubmed"],
    )

    date_col = find_col(
        gwas,
        ["DATE", "date"],
        contains=["date"],
    )

    initial_sample_col = find_col(
        gwas,
        ["INITIAL SAMPLE SIZE", "initial_sample_size"],
        contains=["initial", "sample"],
    )

    replication_sample_col = find_col(
        gwas,
        ["REPLICATION SAMPLE SIZE", "replication_sample_size"],
        contains=["replication", "sample"],
    )

    log("[DETECTED COLUMNS]")
    log(f"  reported_gene_col      = {reported_gene_col}")
    log(f"  mapped_gene_col        = {mapped_gene_col}")
    log(f"  trait_col              = {trait_col}")
    log(f"  efo_col                = {efo_col}")
    log(f"  pvalue_col             = {pvalue_col}")
    log(f"  snp_col                = {snp_col}")
    log(f"  study_col              = {study_col}")
    log(f"  pubmed_col             = {pubmed_col}")
    log(f"  date_col               = {date_col}")
    log(f"  initial_sample_col     = {initial_sample_col}")
    log(f"  replication_sample_col = {replication_sample_col}")

    if reported_gene_col is None and mapped_gene_col is None:
        raise RuntimeError("Could not find REPORTED GENE(S) or MAPPED_GENE columns.")

    rows = []

    for idx, row in gwas.iterrows():
        trait = clean_text(row.get(trait_col, "")) if trait_col else ""
        mapped_trait = clean_text(row.get(efo_col, "")) if efo_col else ""
        pval = safe_float(row.get(pvalue_col, np.nan)) if pvalue_col else np.nan
        nlogp = neglog10_pvalue(pval)
        snp = clean_text(row.get(snp_col, "")) if snp_col else ""
        study = clean_text(row.get(study_col, "")) if study_col else ""
        pubmed = clean_text(row.get(pubmed_col, "")) if pubmed_col else ""
        date = clean_text(row.get(date_col, "")) if date_col else ""
        initial_sample = clean_text(row.get(initial_sample_col, "")) if initial_sample_col else ""
        replication_sample = clean_text(row.get(replication_sample_col, "")) if replication_sample_col else ""

        trait_for_cat = " ".join([trait, mapped_trait])
        categories = classify_trait(trait_for_cat)

        reported_raw = split_gene_field(row.get(reported_gene_col, "")) if reported_gene_col else []
        mapped_raw = split_gene_field(row.get(mapped_gene_col, "")) if mapped_gene_col else []

        reported_genes = map_gene_symbols(reported_raw, symbol_map)
        mapped_genes = map_gene_symbols(mapped_raw, symbol_map)

        for gene in reported_genes:
            rows.append(
                {
                    "gene_symbol": gene,
                    "gwas_gene_source": "reported",
                    "raw_gene_symbols": ",".join(reported_raw),
                    "trait": trait,
                    "mapped_trait": mapped_trait,
                    "trait_categories": ",".join(sorted(categories)),
                    "pvalue": pval,
                    "neglog10_pvalue": nlogp,
                    "snp": snp,
                    "study_accession": study,
                    "pubmed_id": pubmed,
                    "date": date,
                    "initial_sample_size_text": initial_sample,
                    "replication_sample_size_text": replication_sample,
                }
            )

        for gene in mapped_genes:
            rows.append(
                {
                    "gene_symbol": gene,
                    "gwas_gene_source": "mapped",
                    "raw_gene_symbols": ",".join(mapped_raw),
                    "trait": trait,
                    "mapped_trait": mapped_trait,
                    "trait_categories": ",".join(sorted(categories)),
                    "pvalue": pval,
                    "neglog10_pvalue": nlogp,
                    "snp": snp,
                    "study_accession": study,
                    "pubmed_id": pubmed,
                    "date": date,
                    "initial_sample_size_text": initial_sample,
                    "replication_sample_size_text": replication_sample,
                }
            )

        combined = sorted(set(reported_genes) | set(mapped_genes))

        for gene in combined:
            rows.append(
                {
                    "gene_symbol": gene,
                    "gwas_gene_source": "any",
                    "raw_gene_symbols": ",".join(sorted(set(reported_raw + mapped_raw))),
                    "trait": trait,
                    "mapped_trait": mapped_trait,
                    "trait_categories": ",".join(sorted(categories)),
                    "pvalue": pval,
                    "neglog10_pvalue": nlogp,
                    "snp": snp,
                    "study_accession": study,
                    "pubmed_id": pubmed,
                    "date": date,
                    "initial_sample_size_text": initial_sample,
                    "replication_sample_size_text": replication_sample,
                }
            )

        if (idx + 1) % 100000 == 0:
            log(f"[GWAS LONG] processed={idx + 1:,} rows={len(rows):,}")

    long_df = pd.DataFrame(rows)

    if not long_df.empty:
        long_df = long_df.drop_duplicates()

    outpath = processed_dir / "feature13_gwas_gene_association_long.csv"
    long_df.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED GWAS GENE LONG]")
    log(f"[PATH]  {outpath}")
    log(f"[SHAPE] {long_df.shape}")
    log(f"[GENES] {long_df['gene_symbol'].nunique() if not long_df.empty else 0}")

    return long_df


# =============================================================================
# AGGREGATION
# =============================================================================

def count_by_source(g: pd.DataFrame, source: str) -> pd.DataFrame:
    if g.empty or "gwas_gene_source" not in g.columns:
        return pd.DataFrame()

    return g[g["gwas_gene_source"].eq(source)].copy()


def aggregate_gwas_features(
    hgnc: pd.DataFrame,
    long_df: pd.DataFrame,
    processed_dir: Path,
    include_drug_response_category: bool,
) -> pd.DataFrame:
    log("=" * 100)
    log("[AGGREGATE GWAS FEATURES PER GENE]")

    genes = sorted(hgnc["gene_symbol"].unique().tolist())
    by_gene = dict(tuple(long_df.groupby("gene_symbol"))) if not long_df.empty else {}

    records = []

    trait_categories_all = [
        "anthropometric",
        "metabolic",
        "cardiovascular",
        "immune_inflammatory",
        "neurological_psychiatric",
        "cancer",
        "hematological",
        "renal_liver",
        "reproductive",
        "biomarker",
        "infectious",
        "other",
    ]

    if include_drug_response_category:
        trait_categories_all.append("drug_response")

    for gene in genes:
        g = by_gene.get(gene, pd.DataFrame())

        if g.empty:
            rec = {
                "gene_symbol": gene,
                "feature13_gwas_has_any_association": 0,
                "feature13_gwas_association_count_any": 0,
                "feature13_gwas_association_count_reported": 0,
                "feature13_gwas_association_count_mapped": 0,
                "feature13_gwas_unique_snp_count": 0,
                "feature13_gwas_unique_trait_count": 0,
                "feature13_gwas_unique_mapped_trait_count": 0,
                "feature13_gwas_unique_all_trait_count": 0,
                "feature13_gwas_unique_study_count": 0,
                "feature13_gwas_unique_pubmed_count": 0,
                "feature13_gwas_trait_entropy": 0.0,
                "feature13_gwas_mapped_trait_entropy": 0.0,
                "feature13_gwas_all_trait_entropy": 0.0,
                "feature13_gwas_genomewide_sig_count_p_le_5e_8": 0,
                "feature13_gwas_strong_sig_count_p_le_5e_9": 0,
                "feature13_gwas_suggestive_count_p_le_1e_5": 0,
                "feature13_gwas_min_pvalue": np.nan,
                "feature13_gwas_median_pvalue": np.nan,
                "feature13_gwas_mean_pvalue": np.nan,
                "feature13_gwas_max_neglog10_pvalue": np.nan,
                "feature13_gwas_mean_neglog10_pvalue": np.nan,
                "feature13_gwas_median_neglog10_pvalue": np.nan,
                "feature13_gwas_first_association_year": np.nan,
                "feature13_gwas_last_association_year": np.nan,
                "feature13_gwas_association_year_span": np.nan,
                "feature13_gwas_trait_category_count": 0,
                "feature13_gwas_trait_category_entropy": 0.0,
            }

            for cat in trait_categories_all:
                rec[f"feature13_gwas_traitcat_{cat}_count"] = 0

            records.append(rec)
            continue

        any_g = count_by_source(g, "any")
        rep_g = count_by_source(g, "reported")
        map_g = count_by_source(g, "mapped")

        if any_g.empty:
            any_g = g.copy()

        pvals = (
            pd.to_numeric(any_g["pvalue"], errors="coerce").dropna()
            if "pvalue" in any_g.columns
            else pd.Series(dtype=float)
        )

        nlogp = (
            pd.to_numeric(any_g["neglog10_pvalue"], errors="coerce").dropna()
            if "neglog10_pvalue" in any_g.columns
            else pd.Series(dtype=float)
        )

        traits = (
            any_g["trait"].dropna().astype(str).replace("", np.nan).dropna().tolist()
            if "trait" in any_g.columns
            else []
        )

        mapped_traits = (
            any_g["mapped_trait"].dropna().astype(str).replace("", np.nan).dropna().tolist()
            if "mapped_trait" in any_g.columns
            else []
        )

        all_trait_names = [x for x in traits + mapped_traits if clean_text(x)]

        rec = {
            "gene_symbol": gene,
            "feature13_gwas_has_any_association": 1,

            "feature13_gwas_association_count_any": int(len(any_g)),
            "feature13_gwas_association_count_reported": int(len(rep_g)),
            "feature13_gwas_association_count_mapped": int(len(map_g)),

            "feature13_gwas_unique_snp_count": int(any_g["snp"].replace("", np.nan).dropna().nunique()) if "snp" in any_g.columns else 0,
            "feature13_gwas_unique_trait_count": int(pd.Series(traits).nunique()) if traits else 0,
            "feature13_gwas_unique_mapped_trait_count": int(pd.Series(mapped_traits).nunique()) if mapped_traits else 0,
            "feature13_gwas_unique_all_trait_count": int(pd.Series(all_trait_names).nunique()) if all_trait_names else 0,
            "feature13_gwas_unique_study_count": int(any_g["study_accession"].replace("", np.nan).dropna().nunique()) if "study_accession" in any_g.columns else 0,
            "feature13_gwas_unique_pubmed_count": int(any_g["pubmed_id"].replace("", np.nan).dropna().nunique()) if "pubmed_id" in any_g.columns else 0,

            "feature13_gwas_trait_entropy": shannon_entropy(traits),
            "feature13_gwas_mapped_trait_entropy": shannon_entropy(mapped_traits),
            "feature13_gwas_all_trait_entropy": shannon_entropy(all_trait_names),

            "feature13_gwas_genomewide_sig_count_p_le_5e_8": int((pvals <= 5e-8).sum()) if len(pvals) else 0,
            "feature13_gwas_strong_sig_count_p_le_5e_9": int((pvals <= 5e-9).sum()) if len(pvals) else 0,
            "feature13_gwas_suggestive_count_p_le_1e_5": int((pvals <= 1e-5).sum()) if len(pvals) else 0,

            "feature13_gwas_min_pvalue": float(pvals.min()) if len(pvals) else np.nan,
            "feature13_gwas_median_pvalue": float(pvals.median()) if len(pvals) else np.nan,
            "feature13_gwas_mean_pvalue": float(pvals.mean()) if len(pvals) else np.nan,

            "feature13_gwas_max_neglog10_pvalue": float(nlogp.max()) if len(nlogp) else np.nan,
            "feature13_gwas_mean_neglog10_pvalue": float(nlogp.mean()) if len(nlogp) else np.nan,
            "feature13_gwas_median_neglog10_pvalue": float(nlogp.median()) if len(nlogp) else np.nan,
        }

        years = []

        if "date" in any_g.columns:
            for d in any_g["date"].dropna().astype(str).tolist():
                m = re.search(r"(19|20)\d{2}", d)

                if m:
                    years.append(int(m.group(0)))

        if years:
            rec["feature13_gwas_first_association_year"] = int(min(years))
            rec["feature13_gwas_last_association_year"] = int(max(years))
            rec["feature13_gwas_association_year_span"] = int(max(years) - min(years))

        else:
            rec["feature13_gwas_first_association_year"] = np.nan
            rec["feature13_gwas_last_association_year"] = np.nan
            rec["feature13_gwas_association_year_span"] = np.nan

        cat_counter = Counter()

        if "trait_categories" in any_g.columns:
            for val in any_g["trait_categories"].dropna().astype(str):
                for cat in val.split(","):
                    cat = cat.strip()

                    if not cat:
                        continue

                    if cat == "drug_response" and not include_drug_response_category:
                        continue

                    cat_counter[cat] += 1

        for cat in trait_categories_all:
            rec[f"feature13_gwas_traitcat_{cat}_count"] = int(cat_counter.get(cat, 0))

        rec["feature13_gwas_trait_category_count"] = int(
            sum(1 for c in trait_categories_all if cat_counter.get(c, 0) > 0)
        )

        rec["feature13_gwas_trait_category_entropy"] = shannon_entropy(
            [cat for cat, count in cat_counter.items() for _ in range(count)]
        )

        records.append(rec)

    features = pd.DataFrame(records)

    features["feature13_gwas_log1p_association_count_any"] = np.log1p(
        pd.to_numeric(features["feature13_gwas_association_count_any"], errors="coerce").fillna(0)
    )

    features["feature13_gwas_log1p_unique_trait_count"] = np.log1p(
        pd.to_numeric(features["feature13_gwas_unique_all_trait_count"], errors="coerce").fillna(0)
    )

    features["feature13_gwas_log1p_unique_study_count"] = np.log1p(
        pd.to_numeric(features["feature13_gwas_unique_study_count"], errors="coerce").fillna(0)
    )

    features["feature13_gwas_pleiotropy_index"] = (
        features["feature13_gwas_log1p_association_count_any"].fillna(0)
        * features["feature13_gwas_log1p_unique_trait_count"].fillna(0)
    )

    for c in features.columns:
        if c.startswith("feature13_") and (
            "_count" in c
            or "_has_" in c
            or c.endswith("_index")
            or c.startswith("feature13_gwas_log1p")
        ):
            if features[c].dtype != object:
                features[c] = features[c].fillna(0)

    outpath = processed_dir / "feature13_gwas_catalog_gene_features.csv"
    features.to_csv(outpath, index=False)

    log("[SAVED GWAS GENE FEATURES]")
    log(f"[PATH]  {outpath}")
    log(f"[SHAPE] {features.shape}")
    log(f"[COVERAGE] {features['feature13_gwas_has_any_association'].mean():.4f}")

    return features


# =============================================================================
# MERGE
# =============================================================================

def merge_with_hgnc(
    hgnc: pd.DataFrame,
    features: pd.DataFrame,
    processed_dir: Path,
) -> pd.DataFrame:
    merged = hgnc.merge(features, on="gene_symbol", how="left")

    feature_cols = [c for c in features.columns if c != "gene_symbol"]

    merged["feature13_gwas_catalog_has_any_feature"] = (
        merged[feature_cols].notna().any(axis=1).astype(int)
    )

    merged["feature13_gwas_catalog_n_nonmissing_features"] = (
        merged[feature_cols].notna().sum(axis=1)
    )

    outpath = processed_dir / "feature13_gwas_catalog_hgnc_merged.csv"
    merged.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED HGNC MERGED]")
    log(f"[PATH]  {outpath}")
    log(f"[SHAPE] {merged.shape}")

    return merged


# =============================================================================
# SUMMARY
# =============================================================================

def write_summary(
    outdir: Path,
    hgnc: pd.DataFrame,
    gwas_file: Path,
    gwas: pd.DataFrame,
    long_df: pd.DataFrame,
    features: pd.DataFrame,
    merged: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    path = outdir / "feature13_gwas_catalog_summary.txt"

    lines = []
    lines.append("Feature 13: GWAS Catalog gene-association burden / pleiotropy")
    lines.append("=" * 80)
    lines.append(f"Created at: {now_iso()}")
    lines.append("")
    lines.append("Inputs:")
    lines.append(f"HGNC: {args.hgnc}")
    lines.append(f"GWAS Catalog association file: {gwas_file}")
    lines.append(f"GWAS Catalog download URL: {GWAS_ASSOCIATIONS_URL}")
    lines.append(f"GWAS file is ZIP: {is_zip_file(gwas_file)}")
    lines.append(f"GWAS file is gzip: {is_gzip_file(gwas_file)}")
    lines.append("")
    lines.append("Rows:")
    lines.append(f"HGNC rows: {hgnc.shape[0]}")
    lines.append(f"GWAS association rows read: {gwas.shape[0]}")
    lines.append(f"GWAS gene-long rows: {long_df.shape[0]}")
    lines.append(f"Genes with GWAS association: {int(features['feature13_gwas_has_any_association'].sum())}")
    lines.append(f"GWAS coverage: {features['feature13_gwas_has_any_association'].mean():.4f}")
    lines.append(f"Gene feature rows: {features.shape[0]}")
    lines.append(f"Gene feature columns: {features.shape[1]}")
    lines.append(f"Merged rows: {merged.shape[0]}")
    lines.append(f"Merged columns: {merged.shape[1]}")
    lines.append("")
    lines.append("Important feature summaries:")

    for col in [
        "feature13_gwas_association_count_any",
        "feature13_gwas_unique_all_trait_count",
        "feature13_gwas_unique_study_count",
        "feature13_gwas_max_neglog10_pvalue",
        "feature13_gwas_pleiotropy_index",
    ]:
        if col in features.columns:
            x = pd.to_numeric(features[col], errors="coerce")
            lines.append(
                f"{col}: median={x.median(skipna=True):.4f}, max={x.max(skipna=True):.4f}"
            )

    lines.append("")
    lines.append("Interpretation:")
    lines.append("These features represent GWAS Catalog association burden and pleiotropy for each gene.")
    lines.append("They should not be interpreted as mutation burden or causal target evidence.")
    lines.append("Mapped genes may reflect nearest/overlapping genes rather than experimentally proven causal genes.")
    lines.append("")
    lines.append("Leakage policy:")
    lines.append("Included: GWAS Catalog association counts, traits, variants, studies, publications, p-value burden.")
    lines.append("Excluded: Open Targets association scores, known drug labels, clinical target labels, ChEMBL, DrugBank, DGIdb, Pharos.")

    path.write_text("\n".join(lines) + "\n")

    log(f"[SUMMARY] {path}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build Feature 13 GWAS Catalog gene-association burden features."
    )

    parser.add_argument("--hgnc", default=str(DEFAULT_HGNC_PATH), help="Path to HGNC complete set.")
    parser.add_argument("--dbdir", default=str(DEFAULT_DBDIR), help="GWAS Catalog database directory.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory.")
    parser.add_argument("--gwas-file", default="", help="Manual GWAS Catalog associations TSV/TSV.GZ/ZIP file.")
    parser.add_argument("--download", action="store_true", help="Download GWAS Catalog association file.")
    parser.add_argument("--force-download", action="store_true", help="Re-download GWAS Catalog association file.")
    parser.add_argument("--max-rows", type=int, default=0, help="Debug only: read first N GWAS rows.")
    parser.add_argument("--limit-genes", type=int, default=0, help="Debug only: first N HGNC genes.")
    parser.add_argument("--all-hgnc-genes", action="store_true", help="Use all HGNC genes, not only protein-coding.")
    parser.add_argument(
        "--include-drug-response-category",
        action="store_true",
        help="Include drug-response trait-category feature. Default excludes it to be conservative.",
    )

    args = parser.parse_args()

    dbdir = mkdir(Path(args.dbdir))
    outdir = mkdir(Path(args.outdir))
    downloads_dir = mkdir(outdir / "downloads")
    processed_dir = mkdir(outdir / "processed")

    default_db_gwas_file = dbdir / "gwas_catalog_associations.tsv"
    out_download_file = downloads_dir / "gwas_catalog_associations.tsv"

    log("=" * 100)
    log("FEATURE 13: GWAS CATALOG GENE-ASSOCIATION BURDEN")
    log("=" * 100)
    log(f"[HGNC]                {args.hgnc}")
    log(f"[DBDIR]               {dbdir.resolve()}")
    log(f"[OUTDIR]              {outdir.resolve()}")
    log(f"[GWAS FILE]           {args.gwas_file if args.gwas_file else 'auto'}")
    log(f"[DOWNLOAD]            {args.download}")
    log(f"[FORCE DOWNLOAD]      {args.force_download}")
    log(f"[MAX ROWS]            {args.max_rows if args.max_rows else 'none'}")
    log(f"[LIMIT GENES]         {args.limit_genes if args.limit_genes else 'none'}")
    log(f"[DRUG RESPONSE CAT]   {args.include_drug_response_category}")
    log("=" * 100)

    if args.gwas_file:
        gwas_file = Path(args.gwas_file)

    else:
        gwas_file = find_local_gwas_file(dbdir)

        if args.download or args.force_download or gwas_file is None:
            download_file(
                GWAS_ASSOCIATIONS_URL,
                default_db_gwas_file,
                force=args.force_download,
            )
            gwas_file = default_db_gwas_file

    if not gwas_file.exists():
        raise FileNotFoundError(
            "GWAS Catalog association file not found. Run with --download or provide --gwas-file."
        )

    inspect_archive_or_text(gwas_file)

    try:
        if gwas_file.resolve() != out_download_file.resolve():
            if not out_download_file.exists() or args.force_download:
                shutil.copy2(gwas_file, out_download_file)

    except Exception as exc:
        log(f"[WARN] Could not copy GWAS file to output downloads directory: {exc}")

    hgnc = load_hgnc(
        Path(args.hgnc),
        protein_coding_only=not args.all_hgnc_genes,
        limit_genes=args.limit_genes,
    )

    symbol_map = build_symbol_alias_map(hgnc)

    gwas = read_gwas_associations(
        path=gwas_file,
        max_rows=args.max_rows,
    )

    long_df = build_gwas_gene_long(
        gwas=gwas,
        symbol_map=symbol_map,
        processed_dir=processed_dir,
    )

    features = aggregate_gwas_features(
        hgnc=hgnc,
        long_df=long_df,
        processed_dir=processed_dir,
        include_drug_response_category=args.include_drug_response_category,
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
        "gwas_file": str(gwas_file.resolve()),
        "gwas_file_is_zip_by_magic_bytes": is_zip_file(gwas_file),
        "gwas_file_is_gzip_by_magic_bytes": is_gzip_file(gwas_file),
        "gwas_download_url": GWAS_ASSOCIATIONS_URL,
        "protein_coding_only": not args.all_hgnc_genes,
        "max_rows": args.max_rows,
        "limit_genes": args.limit_genes,
        "include_drug_response_category": args.include_drug_response_category,
        "outputs": {
            "gwas_gene_long": str(processed_dir / "feature13_gwas_gene_association_long.csv"),
            "gene_features": str(processed_dir / "feature13_gwas_catalog_gene_features.csv"),
            "hgnc_merged": str(processed_dir / "feature13_gwas_catalog_hgnc_merged.csv"),
            "summary": str(outdir / "feature13_gwas_catalog_summary.txt"),
        },
        "leakage_policy": {
            "included": [
                "GWAS Catalog association burden",
                "reported and mapped gene counts",
                "trait diversity",
                "study/publication/variant counts",
                "p-value burden",
                "broad trait categories",
            ],
            "excluded": [
                "Open Targets genetic association scores",
                "known drug labels",
                "clinical target labels",
                "ChEMBL",
                "DrugBank",
                "DGIdb",
                "Pharos",
                "target tractability labels",
            ],
        },
    }

    metadata_path = outdir / "feature13_gwas_catalog_run_metadata.json"

    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=2)

    log(f"[METADATA] {metadata_path}")

    write_summary(
        outdir=outdir,
        hgnc=hgnc,
        gwas_file=gwas_file,
        gwas=gwas,
        long_df=long_df,
        features=features,
        merged=merged,
        args=args,
    )

    log("=" * 100)
    log("[DONE]")
    log(f"[GWAS LONG]     {processed_dir / 'feature13_gwas_gene_association_long.csv'}")
    log(f"[GENE FEATURES] {processed_dir / 'feature13_gwas_catalog_gene_features.csv'}")
    log(f"[HGNC MERGED]   {processed_dir / 'feature13_gwas_catalog_hgnc_merged.csv'}")
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