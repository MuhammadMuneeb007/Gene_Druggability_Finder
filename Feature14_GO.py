#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature14_GeneOntology.py

Feature 14: Gene Ontology annotation features.

Purpose
-------
Build leakage-aware gene-level Gene Ontology features from GOA human annotations.

This feature captures:
    - number of GO annotations per gene
    - number of unique GO terms per gene
    - Biological Process / Molecular Function / Cellular Component counts
    - evidence-code burden
    - annotation diversity
    - broad GO keyword category features
    - optional ontology-depth features from go-basic.obo

Inputs
------
HGNC:
    databases/HGNC/hgnc_complete_set.txt

GO annotation:
    feature_databases/GeneOntology/goa_human.gaf.gz

GO ontology:
    feature_databases/GeneOntology/go-basic.obo

Outputs
-------
feature14_gene_ontology/
    downloads/goa_human.gaf.gz
    downloads/go-basic.obo
    processed/feature14_go_annotation_long.csv
    processed/feature14_go_gene_features.csv
    processed/feature14_go_hgnc_merged.csv
    feature14_go_summary.txt
    feature14_go_run_metadata.json

Run
---
    python Feature14_GeneOntology.py --download

Fast test:
    python Feature14_GeneOntology.py --download --max-gaf-lines 200000

Leakage policy
--------------
Included:
    GO biological process / molecular function / cellular component annotations
    GO evidence code summaries
    GO term diversity
    broad functional keyword categories

Excluded:
    drug labels
    clinical target labels
    ChEMBL
    DrugBank
    DGIdb
    Open Targets
    Pharos
    target tractability labels
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
from collections import Counter, defaultdict, deque
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pandas as pd


# =============================================================================
# DEFAULTS
# =============================================================================

DEFAULT_HGNC_PATH = Path("databases") / "HGNC" / "hgnc_complete_set.txt"
DEFAULT_DBDIR = Path("feature_databases") / "GeneOntology"
DEFAULT_OUTDIR = Path("feature14_gene_ontology")

GOA_HUMAN_GAF_URL = "https://current.geneontology.org/annotations/goa_human.gaf.gz"
GO_BASIC_OBO_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"

# GAF 2.2 columns.
GAF_COLUMNS = [
    "DB",
    "DB_Object_ID",
    "DB_Object_Symbol",
    "Qualifier",
    "GO_ID",
    "DB_Reference",
    "Evidence_Code",
    "With_From",
    "Aspect",
    "DB_Object_Name",
    "DB_Object_Synonym",
    "DB_Object_Type",
    "Taxon",
    "Date",
    "Assigned_By",
    "Annotation_Extension",
    "Gene_Product_Form_ID",
]

ASPECT_MAP = {
    "P": "biological_process",
    "F": "molecular_function",
    "C": "cellular_component",
}

EVIDENCE_GROUPS = {
    "experimental": {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"},
    "phylogenetic": {"IBA", "IBD", "IKR", "IRD"},
    "computational": {"ISS", "ISO", "ISA", "ISM", "IGC", "RCA"},
    "author_statement": {"TAS", "NAS"},
    "curator_statement": {"IC", "ND"},
    "electronic": {"IEA"},
}

# Broad GO keyword categories. These are intentionally generic functional summaries,
# not drug-target labels.
GO_KEYWORD_CATEGORIES = {
    "kinase_phosphorylation": [
        "kinase", "phosphorylation", "phosphotransferase", "protein phosphorylation",
    ],
    "phosphatase_dephosphorylation": [
        "phosphatase", "dephosphorylation",
    ],
    "receptor_signaling": [
        "receptor", "signal transduction", "signaling pathway", "cell surface receptor",
    ],
    "transporter_channel": [
        "transporter", "transport", "channel", "ion transmembrane", "solute",
    ],
    "enzyme_metabolism": [
        "metabolic process", "catalytic activity", "oxidoreductase", "transferase",
        "hydrolase", "ligase", "lyase", "isomerase",
    ],
    "transcription_gene_regulation": [
        "transcription", "gene expression", "dna-binding transcription factor",
        "regulation of transcription", "rna polymerase",
    ],
    "rna_processing": [
        "rna processing", "rna splicing", "mrna", "ribosome", "translation",
    ],
    "dna_repair_replication": [
        "dna repair", "dna replication", "double-strand break", "chromosome segregation",
    ],
    "cell_cycle": [
        "cell cycle", "mitosis", "meiosis", "checkpoint", "cytokinesis",
    ],
    "apoptosis_cell_death": [
        "apoptosis", "programmed cell death", "cell death",
    ],
    "immune_inflammatory": [
        "immune", "inflammatory", "cytokine", "interleukin", "antigen", "complement",
    ],
    "development_differentiation": [
        "development", "differentiation", "morphogenesis", "organogenesis",
    ],
    "neuronal_synaptic": [
        "neuron", "synapse", "synaptic", "neurotransmitter", "axon", "dendrite",
    ],
    "extracellular_matrix_adhesion": [
        "extracellular matrix", "cell adhesion", "collagen", "integrin", "matrix",
    ],
    "membrane_cell_surface": [
        "membrane", "plasma membrane", "cell surface", "integral component of membrane",
    ],
    "mitochondrial": [
        "mitochondrion", "mitochondrial", "respiratory chain", "oxidative phosphorylation",
    ],
    "nuclear_chromatin": [
        "nucleus", "chromatin", "histone", "nucleosome",
    ],
    "secreted_extracellular": [
        "extracellular region", "secreted", "extracellular space",
    ],
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


def clean_uniprot(x: Any) -> str:
    s = clean_text(x)
    if not s:
        return ""
    s = s.replace("UniProtKB:", "")
    s = s.replace("UniProt:", "")
    s = s.replace("sp|", "")
    s = s.replace("tr|", "")
    s = s.split("|")[0]
    s = s.split("-")[0]
    s = s.split(".")[0]
    return s.strip()


def split_uniprot_ids(x: Any) -> List[str]:
    s = clean_text(x)
    if not s:
        return []

    out = []
    for part in re.split(r"[|;, ]+", s):
        acc = clean_uniprot(part)
        if acc and acc not in out:
            out.append(acc)
    return out


def clean_ensembl_gene_id(x: Any) -> str:
    s = clean_text(x)
    if not s:
        return ""
    return s.replace("gene:", "").split(".")[0]


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


def open_text_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", errors="ignore")
    return open(path, "rt", errors="ignore")


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

    with urllib.request.urlopen(url) as response:
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

    tmp.rename(outpath)
    log(f"[DOWNLOAD DONE] {outpath} ({outpath.stat().st_size / 1024 / 1024:.2f} MB)")


def find_local_file(dbdir: Path, patterns: List[str]) -> Optional[Path]:
    if not dbdir.exists():
        return None

    files = []
    for pat in patterns:
        files.extend(list(dbdir.rglob(pat)))

    files = [p for p in files if p.exists() and p.stat().st_size > 0]
    if not files:
        return None

    files = sorted(files, key=lambda p: p.stat().st_size, reverse=True)
    return files[0]


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
            hgnc = hgnc[hgnc["locus_group"].astype(str).str.lower().eq("protein-coding gene")].copy()
        elif "locus_type" in hgnc.columns:
            hgnc = hgnc[hgnc["locus_type"].astype(str).str.lower().str.contains("protein", na=False)].copy()

    hgnc["gene_symbol"] = hgnc["symbol"].map(normalize_symbol)

    for col in ["uniprot_ids", "alias_symbol", "prev_symbol", "ensembl_gene_id", "entrez_id", "name"]:
        if col not in hgnc.columns:
            hgnc[col] = ""

    hgnc["ensembl_gene_id_clean"] = hgnc["ensembl_gene_id"].map(clean_ensembl_gene_id)
    hgnc = hgnc.drop_duplicates(subset=["gene_symbol"], keep="first").copy()

    if limit_genes and limit_genes > 0:
        hgnc = hgnc.head(limit_genes).copy()

    log(f"[HGNC ROWS] {before} -> {len(hgnc)}")
    log(f"[PROTEIN CODING ONLY] {protein_coding_only}")
    log(f"[GENES WITH UNIPROT] {(hgnc['uniprot_ids'].astype(str).str.len() > 0).sum()}")

    return hgnc


def build_maps(hgnc: pd.DataFrame) -> Tuple[Dict[str, Set[str]], Dict[str, str]]:
    uniprot_to_genes: Dict[str, Set[str]] = defaultdict(set)
    symbol_alias_to_gene: Dict[str, str] = {}

    for _, row in hgnc.iterrows():
        gene = normalize_symbol(row["gene_symbol"])
        if not gene:
            continue

        symbol_alias_to_gene[gene] = gene

        for field in ["alias_symbol", "prev_symbol"]:
            raw = clean_text(row.get(field, ""))
            if raw:
                for part in re.split(r"[|,;]+", raw):
                    s = normalize_symbol(part)
                    if s and s not in symbol_alias_to_gene:
                        symbol_alias_to_gene[s] = gene

        for acc in split_uniprot_ids(row.get("uniprot_ids", "")):
            uniprot_to_genes[acc].add(gene)

    log("=" * 100)
    log("[ID MAPS]")
    log(f"[UNIPROT IDS] {len(uniprot_to_genes)}")
    log(f"[SYMBOLS + ALIASES] {len(symbol_alias_to_gene)}")

    return uniprot_to_genes, symbol_alias_to_gene


# =============================================================================
# OBO PARSER
# =============================================================================

def parse_go_obo(obo_path: Optional[Path]) -> Dict[str, Dict[str, Any]]:
    """
    Parse go-basic.obo.

    Returns:
        go_id -> {
            name,
            namespace,
            parents,
            depth
        }
    """
    if obo_path is None or not obo_path.exists():
        log("[OBO] Not found. GO term names/depths will be missing.")
        return {}

    log("=" * 100)
    log(f"[READ GO OBO] {obo_path}")

    terms: Dict[str, Dict[str, Any]] = {}
    current: Dict[str, Any] = {}
    in_term = False

    with open_text_maybe_gzip(obo_path) as f:
        for line in f:
            line = line.rstrip("\n")

            if line == "[Term]":
                if in_term and current.get("id") and not current.get("is_obsolete", False):
                    terms[current["id"]] = current
                current = {"parents": []}
                in_term = True
                continue

            if line.startswith("[") and line != "[Term]":
                if in_term and current.get("id") and not current.get("is_obsolete", False):
                    terms[current["id"]] = current
                current = {}
                in_term = False
                continue

            if not in_term:
                continue

            if line.startswith("id: "):
                current["id"] = line.replace("id: ", "").strip()
            elif line.startswith("name: "):
                current["name"] = line.replace("name: ", "").strip()
            elif line.startswith("namespace: "):
                current["namespace"] = line.replace("namespace: ", "").strip()
            elif line.startswith("is_a: "):
                parent = line.replace("is_a: ", "").split()[0].strip()
                if parent:
                    current.setdefault("parents", []).append(parent)
            elif line.startswith("is_obsolete: true"):
                current["is_obsolete"] = True

    if in_term and current.get("id") and not current.get("is_obsolete", False):
        terms[current["id"]] = current

    # Compute simple depth from is_a roots.
    children = defaultdict(list)
    indeg = defaultdict(int)

    for go_id, rec in terms.items():
        indeg.setdefault(go_id, 0)
        for parent in rec.get("parents", []):
            if parent in terms:
                children[parent].append(go_id)
                indeg[go_id] += 1

    roots = [go_id for go_id in terms if indeg[go_id] == 0]
    depth = {go_id: 0 for go_id in roots}
    q = deque(roots)

    while q:
        parent = q.popleft()
        for child in children.get(parent, []):
            d = depth[parent] + 1
            if child not in depth or d > depth[child]:
                depth[child] = d
            q.append(child)

    for go_id in terms:
        terms[go_id]["depth"] = int(depth.get(go_id, 0))

    log(f"[OBO TERMS] {len(terms)}")

    return terms


# =============================================================================
# GAF PARSER
# =============================================================================

def qualifier_is_negative(qualifier: str) -> bool:
    q = clean_text(qualifier).upper()
    return "NOT" in q.split("|") or q == "NOT"


def map_gaf_row_to_genes(
    row: List[str],
    uniprot_to_genes: Dict[str, Set[str]],
    symbol_alias_to_gene: Dict[str, str],
) -> Set[str]:
    genes = set()

    db_object_id = clean_uniprot(row[1]) if len(row) > 1 else ""
    db_object_symbol = normalize_symbol(row[2]) if len(row) > 2 else ""
    synonyms = clean_text(row[10]) if len(row) > 10 else ""

    if db_object_id in uniprot_to_genes:
        genes.update(uniprot_to_genes[db_object_id])

    if db_object_symbol in symbol_alias_to_gene:
        genes.add(symbol_alias_to_gene[db_object_symbol])

    if synonyms:
        for s in re.split(r"[|,;]+", synonyms):
            ss = normalize_symbol(s)
            if ss in symbol_alias_to_gene:
                genes.add(symbol_alias_to_gene[ss])

    return genes


def build_go_long_table(
    gaf_path: Path,
    go_terms: Dict[str, Dict[str, Any]],
    uniprot_to_genes: Dict[str, Set[str]],
    symbol_alias_to_gene: Dict[str, str],
    processed_dir: Path,
    max_gaf_lines: int = 0,
    include_iea: bool = True,
    include_negative_not: bool = False,
) -> pd.DataFrame:
    if not gaf_path.exists():
        raise FileNotFoundError(f"GAF file not found: {gaf_path}")

    log("=" * 100)
    log(f"[READ GOA HUMAN GAF] {gaf_path}")

    rows = []
    n_seen = 0
    n_data = 0
    n_mapped = 0
    n_skip_iea = 0
    n_skip_not = 0

    t0 = time.time()

    with open_text_maybe_gzip(gaf_path) as f:
        for line in f:
            if line.startswith("!"):
                continue

            n_seen += 1

            if max_gaf_lines and n_seen > max_gaf_lines:
                break

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 15:
                continue

            n_data += 1

            qualifier = clean_text(parts[3])
            go_id = clean_text(parts[4])
            evidence = clean_text(parts[6])
            aspect = clean_text(parts[8])
            object_type = clean_text(parts[11])
            date = clean_text(parts[13])
            assigned_by = clean_text(parts[14])

            if not include_negative_not and qualifier_is_negative(qualifier):
                n_skip_not += 1
                continue

            if not include_iea and evidence == "IEA":
                n_skip_iea += 1
                continue

            genes = map_gaf_row_to_genes(parts, uniprot_to_genes, symbol_alias_to_gene)

            if not genes:
                continue

            term = go_terms.get(go_id, {})
            go_name = clean_text(term.get("name", ""))
            namespace = clean_text(term.get("namespace", ""))
            depth = term.get("depth", np.nan)

            if not namespace:
                namespace = ASPECT_MAP.get(aspect, "")

            ev_group = "other"
            for group, codes in EVIDENCE_GROUPS.items():
                if evidence in codes:
                    ev_group = group
                    break

            text_for_keywords = f"{go_name} {namespace}".lower()
            keyword_cats = []
            for cat, kws in GO_KEYWORD_CATEGORIES.items():
                if any(kw.lower() in text_for_keywords for kw in kws):
                    keyword_cats.append(cat)

            for gene in genes:
                rows.append(
                    {
                        "gene_symbol": gene,
                        "go_id": go_id,
                        "go_name": go_name,
                        "go_namespace": namespace,
                        "go_aspect": aspect,
                        "go_depth": depth,
                        "evidence_code": evidence,
                        "evidence_group": ev_group,
                        "qualifier": qualifier,
                        "db": clean_text(parts[0]),
                        "db_object_id": clean_uniprot(parts[1]),
                        "db_object_symbol": normalize_symbol(parts[2]),
                        "db_object_type": object_type,
                        "assigned_by": assigned_by,
                        "annotation_date": date,
                        "go_keyword_categories": ",".join(sorted(keyword_cats)),
                    }
                )
                n_mapped += 1

            if n_data % 500000 == 0:
                elapsed = (time.time() - t0) / 60
                log(f"[GAF] data_lines={n_data:,} mapped_rows={n_mapped:,} elapsed={elapsed:.1f} min")

    long_df = pd.DataFrame(rows)

    if not long_df.empty:
        long_df = long_df.drop_duplicates()

    outpath = processed_dir / "feature14_go_annotation_long.csv"
    long_df.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED GO LONG]")
    log(f"[PATH] {outpath}")
    log(f"[SHAPE] {long_df.shape}")
    log(f"[GENES] {long_df['gene_symbol'].nunique() if not long_df.empty else 0}")
    log(f"[SKIP IEA] {n_skip_iea}")
    log(f"[SKIP NOT] {n_skip_not}")

    return long_df


# =============================================================================
# AGGREGATION
# =============================================================================

def aggregate_go_features(
    hgnc: pd.DataFrame,
    long_df: pd.DataFrame,
    processed_dir: Path,
) -> pd.DataFrame:
    log("=" * 100)
    log("[AGGREGATE GO FEATURES PER GENE]")

    genes = sorted(hgnc["gene_symbol"].unique().tolist())
    by_gene = dict(tuple(long_df.groupby("gene_symbol"))) if not long_df.empty else {}

    namespaces = [
        "biological_process",
        "molecular_function",
        "cellular_component",
    ]

    evidence_groups = sorted(set(EVIDENCE_GROUPS.keys()) | {"other"})
    evidence_codes = sorted(set(code for codes in EVIDENCE_GROUPS.values() for code in codes) | {"IEA", "TAS", "NAS", "IC", "ND"})

    records = []

    for gene in genes:
        g = by_gene.get(gene, pd.DataFrame())

        if g.empty:
            rec = {
                "gene_symbol": gene,
                "feature14_go_has_annotation": 0,
                "feature14_go_annotation_count": 0,
                "feature14_go_unique_term_count": 0,
            }

            for ns in namespaces:
                rec[f"feature14_go_namespace_{ns}_annotation_count"] = 0
                rec[f"feature14_go_namespace_{ns}_unique_term_count"] = 0

            for group in evidence_groups:
                rec[f"feature14_go_evidence_group_{group}_count"] = 0

            for cat in GO_KEYWORD_CATEGORIES:
                rec[f"feature14_go_keyword_{cat}_count"] = 0

            records.append(rec)
            continue

        rec = {
            "gene_symbol": gene,
            "feature14_go_has_annotation": 1,
            "feature14_go_annotation_count": int(len(g)),
            "feature14_go_unique_term_count": int(g["go_id"].nunique()) if "go_id" in g.columns else 0,
            "feature14_go_unique_name_count": int(g["go_name"].replace("", np.nan).dropna().nunique()) if "go_name" in g.columns else 0,
            "feature14_go_unique_evidence_code_count": int(g["evidence_code"].replace("", np.nan).dropna().nunique()) if "evidence_code" in g.columns else 0,
            "feature14_go_unique_assigned_by_count": int(g["assigned_by"].replace("", np.nan).dropna().nunique()) if "assigned_by" in g.columns else 0,
            "feature14_go_term_entropy": shannon_entropy(g["go_id"].tolist()) if "go_id" in g.columns else 0.0,
            "feature14_go_namespace_entropy": shannon_entropy(g["go_namespace"].tolist()) if "go_namespace" in g.columns else 0.0,
            "feature14_go_evidence_entropy": shannon_entropy(g["evidence_code"].tolist()) if "evidence_code" in g.columns else 0.0,
        }

        # Namespace features.
        for ns in namespaces:
            sub = g[g["go_namespace"].eq(ns)] if "go_namespace" in g.columns else pd.DataFrame()
            rec[f"feature14_go_namespace_{ns}_annotation_count"] = int(len(sub))
            rec[f"feature14_go_namespace_{ns}_unique_term_count"] = int(sub["go_id"].nunique()) if not sub.empty and "go_id" in sub.columns else 0

        # Aspect fallback counts.
        for aspect, ns in ASPECT_MAP.items():
            sub = g[g["go_aspect"].eq(aspect)] if "go_aspect" in g.columns else pd.DataFrame()
            rec[f"feature14_go_aspect_{aspect}_annotation_count"] = int(len(sub))
            rec[f"feature14_go_aspect_{aspect}_unique_term_count"] = int(sub["go_id"].nunique()) if not sub.empty and "go_id" in sub.columns else 0

        # Evidence groups.
        if "evidence_group" in g.columns:
            group_counts = Counter(g["evidence_group"].dropna().astype(str).tolist())
        else:
            group_counts = Counter()

        for group in evidence_groups:
            rec[f"feature14_go_evidence_group_{group}_count"] = int(group_counts.get(group, 0))

        # Common evidence codes.
        if "evidence_code" in g.columns:
            ev_counts = Counter(g["evidence_code"].dropna().astype(str).tolist())
        else:
            ev_counts = Counter()

        for code in evidence_codes:
            rec[f"feature14_go_evidence_code_{code}_count"] = int(ev_counts.get(code, 0))

        # Depth features.
        if "go_depth" in g.columns:
            depths = pd.to_numeric(g["go_depth"], errors="coerce").dropna()
            if len(depths):
                rec["feature14_go_depth_mean"] = float(depths.mean())
                rec["feature14_go_depth_median"] = float(depths.median())
                rec["feature14_go_depth_max"] = float(depths.max())
                rec["feature14_go_depth_min"] = float(depths.min())
                rec["feature14_go_deep_term_count_depth_ge_5"] = int((depths >= 5).sum())
                rec["feature14_go_deep_term_count_depth_ge_7"] = int((depths >= 7).sum())
            else:
                rec["feature14_go_depth_mean"] = np.nan
                rec["feature14_go_depth_median"] = np.nan
                rec["feature14_go_depth_max"] = np.nan
                rec["feature14_go_depth_min"] = np.nan
                rec["feature14_go_deep_term_count_depth_ge_5"] = 0
                rec["feature14_go_deep_term_count_depth_ge_7"] = 0

        # Keyword categories.
        keyword_counter = Counter()
        if "go_keyword_categories" in g.columns:
            for val in g["go_keyword_categories"].dropna().astype(str).tolist():
                for cat in val.split(","):
                    cat = cat.strip()
                    if cat:
                        keyword_counter[cat] += 1

        for cat in GO_KEYWORD_CATEGORIES:
            rec[f"feature14_go_keyword_{cat}_count"] = int(keyword_counter.get(cat, 0))

        rec["feature14_go_keyword_category_count"] = int(sum(1 for cat in GO_KEYWORD_CATEGORIES if keyword_counter.get(cat, 0) > 0))
        rec["feature14_go_keyword_category_entropy"] = shannon_entropy(
            [cat for cat, count in keyword_counter.items() for _ in range(count)]
        )

        # Date features.
        years = []
        if "annotation_date" in g.columns:
            for d in g["annotation_date"].dropna().astype(str).tolist():
                m = re.search(r"(19|20)\d{2}", d)
                if m:
                    years.append(int(m.group(0)))

        if years:
            rec["feature14_go_first_annotation_year"] = int(min(years))
            rec["feature14_go_last_annotation_year"] = int(max(years))
            rec["feature14_go_annotation_year_span"] = int(max(years) - min(years))
        else:
            rec["feature14_go_first_annotation_year"] = np.nan
            rec["feature14_go_last_annotation_year"] = np.nan
            rec["feature14_go_annotation_year_span"] = np.nan

        records.append(rec)

    features = pd.DataFrame(records)

    # Derived normalized features.
    features["feature14_go_log1p_annotation_count"] = np.log1p(
        pd.to_numeric(features["feature14_go_annotation_count"], errors="coerce").fillna(0)
    )
    features["feature14_go_log1p_unique_term_count"] = np.log1p(
        pd.to_numeric(features["feature14_go_unique_term_count"], errors="coerce").fillna(0)
    )
    features["feature14_go_annotation_density_index"] = (
        features["feature14_go_log1p_annotation_count"].fillna(0)
        * features["feature14_go_log1p_unique_term_count"].fillna(0)
    )

    # Fill counts and flags only.
    for c in features.columns:
        if c.startswith("feature14_") and (
            "_count" in c
            or "_has_" in c
            or c.startswith("feature14_go_log1p")
            or c.endswith("_index")
        ):
            if features[c].dtype != object:
                features[c] = features[c].fillna(0)

    outpath = processed_dir / "feature14_go_gene_features.csv"
    features.to_csv(outpath, index=False)

    log("[SAVED GO GENE FEATURES]")
    log(f"[PATH] {outpath}")
    log(f"[SHAPE] {features.shape}")
    log(f"[COVERAGE] {features['feature14_go_has_annotation'].mean():.4f}")

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

    merged["feature14_go_has_any_feature"] = merged[feature_cols].notna().any(axis=1).astype(int)
    merged["feature14_go_n_nonmissing_features"] = merged[feature_cols].notna().sum(axis=1)

    outpath = processed_dir / "feature14_go_hgnc_merged.csv"
    merged.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED HGNC MERGED]")
    log(f"[PATH] {outpath}")
    log(f"[SHAPE] {merged.shape}")

    return merged


# =============================================================================
# SUMMARY
# =============================================================================

def write_summary(
    outdir: Path,
    hgnc: pd.DataFrame,
    gaf_path: Path,
    obo_path: Optional[Path],
    long_df: pd.DataFrame,
    features: pd.DataFrame,
    merged: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    path = outdir / "feature14_go_summary.txt"

    lines = []
    lines.append("Feature 14: Gene Ontology annotation features")
    lines.append("=" * 80)
    lines.append(f"Created at: {now_iso()}")
    lines.append("")
    lines.append("Inputs:")
    lines.append(f"HGNC: {args.hgnc}")
    lines.append(f"GOA human GAF: {gaf_path}")
    lines.append(f"GO ontology OBO: {obo_path if obo_path else 'not used'}")
    lines.append(f"GOA URL: {GOA_HUMAN_GAF_URL}")
    lines.append(f"OBO URL: {GO_BASIC_OBO_URL}")
    lines.append("")
    lines.append("Options:")
    lines.append(f"include_iea: {args.include_iea}")
    lines.append(f"include_negative_not: {args.include_negative_not}")
    lines.append(f"max_gaf_lines: {args.max_gaf_lines}")
    lines.append("")
    lines.append("Rows:")
    lines.append(f"HGNC rows: {hgnc.shape[0]}")
    lines.append(f"GO long rows: {long_df.shape[0]}")
    lines.append(f"Genes with GO annotations: {int(features['feature14_go_has_annotation'].sum())}")
    lines.append(f"GO coverage: {features['feature14_go_has_annotation'].mean():.4f}")
    lines.append(f"Gene feature rows: {features.shape[0]}")
    lines.append(f"Gene feature columns: {features.shape[1]}")
    lines.append(f"Merged rows: {merged.shape[0]}")
    lines.append(f"Merged columns: {merged.shape[1]}")
    lines.append("")
    lines.append("Feature summaries:")

    for col in [
        "feature14_go_annotation_count",
        "feature14_go_unique_term_count",
        "feature14_go_namespace_biological_process_unique_term_count",
        "feature14_go_namespace_molecular_function_unique_term_count",
        "feature14_go_namespace_cellular_component_unique_term_count",
        "feature14_go_depth_mean",
        "feature14_go_annotation_density_index",
    ]:
        if col in features.columns:
            x = pd.to_numeric(features[col], errors="coerce")
            lines.append(f"{col}: median={x.median(skipna=True):.4f}, max={x.max(skipna=True):.4f}")

    lines.append("")
    lines.append("Interpretation:")
    lines.append("These features describe functional annotation burden and diversity across GO namespaces.")
    lines.append("They are not direct druggability labels.")
    lines.append("")
    lines.append("Leakage policy:")
    lines.append("Included: GO term counts, namespace counts, evidence-code counts, term-depth summaries, broad GO keyword categories.")
    lines.append("Excluded: clinical target labels, known drug labels, ChEMBL, DrugBank, DGIdb, Open Targets, Pharos.")

    path.write_text("\n".join(lines) + "\n")
    log(f"[SUMMARY] {path}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description="Build Feature 14 Gene Ontology features.")
    parser.add_argument("--hgnc", default=str(DEFAULT_HGNC_PATH), help="Path to HGNC complete set.")
    parser.add_argument("--dbdir", default=str(DEFAULT_DBDIR), help="GeneOntology database directory.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory.")
    parser.add_argument("--gaf-file", default="", help="Manual GOA human GAF file.")
    parser.add_argument("--obo-file", default="", help="Manual go-basic.obo file.")
    parser.add_argument("--download", action="store_true", help="Download GOA human GAF and go-basic.obo.")
    parser.add_argument("--force-download", action="store_true", help="Force re-download.")
    parser.add_argument("--include-iea", action="store_true", help="Include IEA electronic annotations. Default excludes IEA.")
    parser.add_argument("--include-negative-not", action="store_true", help="Include NOT-qualified annotations. Default excludes NOT.")
    parser.add_argument("--max-gaf-lines", type=int, default=0, help="Debug only: first N GAF data lines.")
    parser.add_argument("--limit-genes", type=int, default=0, help="Debug only: first N HGNC genes.")
    parser.add_argument("--all-hgnc-genes", action="store_true", help="Use all HGNC genes, not only protein-coding.")
    args = parser.parse_args()

    dbdir = mkdir(Path(args.dbdir))
    outdir = mkdir(Path(args.outdir))
    downloads_dir = mkdir(outdir / "downloads")
    processed_dir = mkdir(outdir / "processed")

    default_gaf = dbdir / "goa_human.gaf.gz"
    default_obo = dbdir / "go-basic.obo"

    log("=" * 100)
    log("FEATURE 14: GENE ONTOLOGY ANNOTATION FEATURES")
    log("=" * 100)
    log(f"[HGNC]              {args.hgnc}")
    log(f"[DBDIR]             {dbdir.resolve()}")
    log(f"[OUTDIR]            {outdir.resolve()}")
    log(f"[GAF FILE]          {args.gaf_file if args.gaf_file else 'auto'}")
    log(f"[OBO FILE]          {args.obo_file if args.obo_file else 'auto'}")
    log(f"[DOWNLOAD]          {args.download}")
    log(f"[FORCE DOWNLOAD]    {args.force_download}")
    log(f"[INCLUDE IEA]       {args.include_iea}")
    log(f"[INCLUDE NOT]       {args.include_negative_not}")
    log(f"[MAX GAF LINES]     {args.max_gaf_lines if args.max_gaf_lines else 'none'}")
    log(f"[LIMIT GENES]       {args.limit_genes if args.limit_genes else 'none'}")
    log("=" * 100)

    if args.download or args.force_download:
        download_file(GOA_HUMAN_GAF_URL, default_gaf, force=args.force_download)
        download_file(GO_BASIC_OBO_URL, default_obo, force=args.force_download)

    gaf_path = Path(args.gaf_file) if args.gaf_file else find_local_file(
        dbdir,
        ["goa_human.gaf.gz", "*human*.gaf.gz", "*.gaf.gz", "*.gaf"],
    )

    obo_path = Path(args.obo_file) if args.obo_file else find_local_file(
        dbdir,
        ["go-basic.obo", "*.obo"],
    )

    if gaf_path is None or not gaf_path.exists():
        raise FileNotFoundError("GOA human GAF not found. Run with --download or provide --gaf-file.")

    if obo_path is None or not obo_path.exists():
        log("[WARNING] go-basic.obo not found. Run with --download for GO names/depths.")
        obo_path = None

    # Copy downloads into output folder for reproducibility.
    try:
        shutil.copy2(gaf_path, downloads_dir / gaf_path.name)
    except Exception:
        pass

    if obo_path:
        try:
            shutil.copy2(obo_path, downloads_dir / obo_path.name)
        except Exception:
            pass

    hgnc = load_hgnc(
        Path(args.hgnc),
        protein_coding_only=not args.all_hgnc_genes,
        limit_genes=args.limit_genes,
    )

    uniprot_to_genes, symbol_alias_to_gene = build_maps(hgnc)

    go_terms = parse_go_obo(obo_path)

    long_df = build_go_long_table(
        gaf_path=gaf_path,
        go_terms=go_terms,
        uniprot_to_genes=uniprot_to_genes,
        symbol_alias_to_gene=symbol_alias_to_gene,
        processed_dir=processed_dir,
        max_gaf_lines=args.max_gaf_lines,
        include_iea=args.include_iea,
        include_negative_not=args.include_negative_not,
    )

    features = aggregate_go_features(
        hgnc=hgnc,
        long_df=long_df,
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
        "gaf_path": str(gaf_path),
        "obo_path": str(obo_path) if obo_path else "",
        "goa_human_gaf_url": GOA_HUMAN_GAF_URL,
        "go_basic_obo_url": GO_BASIC_OBO_URL,
        "include_iea": args.include_iea,
        "include_negative_not": args.include_negative_not,
        "max_gaf_lines": args.max_gaf_lines,
        "limit_genes": args.limit_genes,
        "protein_coding_only": not args.all_hgnc_genes,
        "outputs": {
            "go_long": str(processed_dir / "feature14_go_annotation_long.csv"),
            "gene_features": str(processed_dir / "feature14_go_gene_features.csv"),
            "hgnc_merged": str(processed_dir / "feature14_go_hgnc_merged.csv"),
            "summary": str(outdir / "feature14_go_summary.txt"),
        },
        "leakage_policy": {
            "included": [
                "GO annotation counts",
                "GO namespace counts",
                "GO evidence-code counts",
                "GO term diversity",
                "GO ontology depth summaries",
                "broad GO keyword categories",
            ],
            "excluded": [
                "clinical target labels",
                "known drug-target labels",
                "ChEMBL",
                "DrugBank",
                "DGIdb",
                "Open Targets",
                "Pharos",
                "target tractability labels",
            ],
        },
    }

    with open(outdir / "feature14_go_run_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    write_summary(
        outdir=outdir,
        hgnc=hgnc,
        gaf_path=gaf_path,
        obo_path=obo_path,
        long_df=long_df,
        features=features,
        merged=merged,
        args=args,
    )

    log("=" * 100)
    log("[DONE]")
    log(f"[GO LONG]       {processed_dir / 'feature14_go_annotation_long.csv'}")
    log(f"[GENE FEATURES] {processed_dir / 'feature14_go_gene_features.csv'}")
    log(f"[HGNC MERGED]   {processed_dir / 'feature14_go_hgnc_merged.csv'}")
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