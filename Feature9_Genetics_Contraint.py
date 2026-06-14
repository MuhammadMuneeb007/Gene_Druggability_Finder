#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature9_DisorderConstraint.py

Feature 9: Protein disorder + genetic constraint / intolerance.

Purpose
-------
Build leakage-safe gene-level features that capture:

    1. Intrinsic protein disorder
       - DisProt curated disorder regions if available
       - MobiDB disorder summary if available
       - AlphaFold pLDDT / low-confidence proxy from Feature 4 if available

    2. Gene constraint / intolerance
       - gnomAD pLI
       - gnomAD LOEUF
       - observed / expected LoF
       - observed / expected missense
       - observed / expected synonymous

Why Feature 9?
--------------
This captures your idea of genes/proteins that are "constant", constrained, or hard
to perturb. In standard genomics language, this is evolutionary constraint /
loss-of-function intolerance. In structural biology, the complementary concept is
protein disorder / flexibility.

Default inputs
--------------
HGNC:
    databases/HGNC/hgnc_complete_set.txt

Optional local disorder/constraint databases:
    feature_databases/DisProt/
    feature_databases/MobiDB/
    feature_databases/gnomAD/
    feature_databases/GnomAD/
    feature_databases/Constraint/

Optional existing feature table:
    feature4_structure/processed/feature4_structure_hgnc_merged.csv

Outputs
-------
feature9_disorder_constraint/
    processed/feature9_disprot_long.csv
    processed/feature9_mobidb_long.csv
    processed/feature9_constraint_gene_table.csv
    processed/feature9_disorder_constraint_gene_features.csv
    processed/feature9_disorder_constraint_hgnc_merged.csv
    feature9_disorder_constraint_summary.txt
    feature9_disorder_constraint_run_metadata.json

Run
---
    python Feature9_DisorderConstraint.py

Fast test:
    python Feature9_DisorderConstraint.py --limit-genes 500

If you want to try automatic downloads:
    python Feature9_DisorderConstraint.py --download

Important
---------
This script is leakage-safe. It does NOT use:
    ChEMBL
    DrugBank
    DGIdb
    Open Targets knownDrugs
    Open Targets tractability
    Pharos Tclin/Tchem
    clinical target labels
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pandas as pd


try:
    import requests
except Exception:
    requests = None


# =============================================================================
# DEFAULTS
# =============================================================================

DEFAULT_HGNC_PATH = Path("databases") / "HGNC" / "hgnc_complete_set.txt"
DEFAULT_FEATURE4_PATH = Path("feature4_structure") / "processed" / "feature4_structure_hgnc_merged.csv"
DEFAULT_DBDIR = Path("feature_databases")
DEFAULT_OUTDIR = Path("feature9_disorder_constraint")

# gnomAD v2.1.1 is stable and commonly used for LOEUF / pLI style features.
# v4 constraint exists, but gnomAD says v4 gene constraint is still experimental/beta.
GNOMAD_V211_CONSTRAINT_URL = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/"
    "constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
)

# DisProt download endpoints can change. The script tries local files first.
# These are best-effort fallbacks only.
DISPROT_DOWNLOAD_CANDIDATES = [
    "https://disprot.org/api/search?release=current&format=tsv",
    "https://disprot.org/api/search?format=tsv",
]


# =============================================================================
# BASIC HELPERS
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
    return clean_text(x).upper()


def clean_uniprot(x: Any) -> str:
    s = clean_text(x)
    if not s:
        return ""
    return s.split("-")[0].split(".")[0]


def clean_ensembl_gene_id(x: Any) -> str:
    s = clean_text(x)
    if not s:
        return ""
    s = s.replace("gene:", "")
    return s.split(".")[0]


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


def safe_int(x: Any) -> int:
    try:
        if x is None or pd.isna(x):
            return 0
        s = str(x).strip()
        if s in {"", ".", "NA", "NaN", "nan", "None"}:
            return 0
        return int(float(s))
    except Exception:
        return 0


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


def open_text_maybe_gzip(path: Path):
    name = str(path).lower()
    if name.endswith(".gz") or name.endswith(".bgz"):
        return gzip.open(path, "rt", errors="ignore")
    return open(path, "rt", errors="ignore")


def read_table_flexible(path: Path) -> pd.DataFrame:
    """
    Read TSV/CSV/TXT/GZ/BGZ robustly.

    Important:
    gnomAD .bgz files are gzip-compatible, but pandas compression='infer'
    does not always infer .bgz correctly. Therefore we explicitly open
    .gz/.bgz files with gzip.open.
    """
    if not path.exists():
        raise FileNotFoundError(path)

    name = path.name.lower()

    # Explicit gzip/BGZF handling.
    if name.endswith(".gz") or name.endswith(".bgz"):
        with gzip.open(path, "rt", errors="ignore") as f:
            if ".csv" in name:
                return pd.read_csv(f, dtype=str, low_memory=False)
            return pd.read_csv(f, sep="\t", dtype=str, low_memory=False)

    # Plain CSV.
    if name.endswith(".csv"):
        return pd.read_csv(path, dtype=str, low_memory=False)

    # Default plain TSV/TXT.
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)


def standardize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    return df


def find_first_existing_file(
    roots: List[Path],
    patterns: List[str],
    prefer_keywords: Optional[List[str]] = None,
) -> Optional[Path]:
    prefer_keywords = prefer_keywords or []
    files = []

    for root in roots:
        if not root.exists():
            continue
        for pattern in patterns:
            files.extend(list(root.rglob(pattern)))

    files = [p for p in files if p.exists() and p.stat().st_size > 0]

    if not files:
        return None

    ranked = []
    for p in files:
        name = p.name.lower()
        score = 0
        for kw in prefer_keywords:
            if kw.lower() in name:
                score += 10
        if name.endswith(".tsv") or name.endswith(".tsv.gz"):
            score += 5
        if name.endswith(".txt") or name.endswith(".txt.gz") or name.endswith(".bgz"):
            score += 3
        ranked.append((score, p.stat().st_size, p))

    ranked.sort(reverse=True)
    return ranked[0][2]


def interval_union_length(intervals: List[Tuple[int, int]]) -> int:
    clean = []

    for a, b in intervals:
        a = safe_int(a)
        b = safe_int(b)

        if a <= 0 or b <= 0:
            continue

        if b < a:
            a, b = b, a

        clean.append((a, b))

    if not clean:
        return 0

    clean.sort()

    total = 0
    cur_a, cur_b = clean[0]

    for a, b in clean[1:]:
        if a <= cur_b + 1:
            cur_b = max(cur_b, b)
        else:
            total += cur_b - cur_a + 1
            cur_a, cur_b = a, b

    total += cur_b - cur_a + 1

    return int(total)


# =============================================================================
# DOWNLOAD HELPERS
# =============================================================================

def download_file(url: str, outpath: Path, timeout: int = 60) -> bool:
    if requests is None:
        log("[DOWNLOAD SKIP] requests is not installed.")
        return False

    mkdir(outpath.parent)

    if outpath.exists() and outpath.stat().st_size > 0:
        log(f"[DOWNLOAD SKIP] Already exists: {outpath}")
        return True

    log(f"[DOWNLOAD] {url}")
    log(f"[TO]       {outpath}")

    try:
        with requests.get(url, stream=True, timeout=timeout) as r:
            if r.status_code != 200:
                log(f"[DOWNLOAD FAILED] HTTP {r.status_code}")
                return False

            with open(outpath, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        f.write(chunk)

        ok = outpath.exists() and outpath.stat().st_size > 0
        log(f"[DOWNLOAD {'OK' if ok else 'FAILED'}] {outpath}")
        return ok

    except Exception as exc:
        log(f"[DOWNLOAD ERROR] {exc}")
        return False


def maybe_download_databases(dbdir: Path) -> None:
    """
    Best-effort downloader for Feature 9 databases.

    Downloads:
        1. gnomAD constraint table
        2. DisProt disorder annotations
        3. MobiDB human TSV, if endpoint works

    Local files are preferred. Existing files are skipped.
    """
    log("=" * 100)
    log("[OPTIONAL DOWNLOAD MODE]")

    gnomad_dir = mkdir(dbdir / "gnomAD")
    disprot_dir = mkdir(dbdir / "DisProt")
    mobidb_dir = mkdir(dbdir / "MobiDB")

    # -------------------------------------------------------------------------
    # gnomAD v2.1.1 constraint
    # -------------------------------------------------------------------------
    download_file(
        GNOMAD_V211_CONSTRAINT_URL,
        gnomad_dir / "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
        timeout=120,
    )

    # -------------------------------------------------------------------------
    # DisProt current TSV
    # -------------------------------------------------------------------------
    disprot_candidates = [
        "https://disprot.org/api/search?release=current&format=tsv",
        "https://disprot.org/api/search?format=tsv",
    ]

    for i, url in enumerate(disprot_candidates, start=1):
        ok = download_file(
            url,
            disprot_dir / f"disprot_download_candidate_{i}.tsv",
            timeout=120,
        )
        if ok:
            break

    # -------------------------------------------------------------------------
    # MobiDB human TSV - optional / best effort
    # -------------------------------------------------------------------------
    mobidb_candidates = [
        "https://mobidb.org/api/download?format=tsv&organism=9606",
        "https://mobidb.bio.unipd.it/api/download?format=tsv&organism=9606",
    ]

    for i, url in enumerate(mobidb_candidates, start=1):
        ok = download_file(
            url,
            mobidb_dir / f"mobidb_human_candidate_{i}.tsv",
            timeout=120,
        )

        if ok:
            outpath = mobidb_dir / f"mobidb_human_candidate_{i}.tsv"

            # Avoid accepting an HTML error page as a database.
            try:
                head = outpath.read_text(errors="ignore")[:300].lower()
                if "<html" in head or "<!doctype html" in head:
                    log(f"[MOBIDB WARNING] Download looks like HTML. Deleting: {outpath}")
                    outpath.unlink()
                    continue
            except Exception:
                pass

            break

    log("[OPTIONAL DOWNLOAD MODE DONE]")


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
        raise RuntimeError("HGNC file must contain a symbol column.")

    before = len(hgnc)

    if protein_coding_only:
        if "locus_group" in hgnc.columns:
            hgnc = hgnc[hgnc["locus_group"].astype(str).str.lower().eq("protein-coding gene")].copy()
        elif "locus_type" in hgnc.columns:
            hgnc = hgnc[hgnc["locus_type"].astype(str).str.lower().str.contains("protein", na=False)].copy()

    hgnc["gene_symbol"] = hgnc["symbol"].map(normalize_symbol)

    if "ensembl_gene_id" not in hgnc.columns:
        hgnc["ensembl_gene_id"] = ""
    if "entrez_id" not in hgnc.columns:
        hgnc["entrez_id"] = ""
    if "uniprot_ids" not in hgnc.columns:
        hgnc["uniprot_ids"] = ""
    if "name" not in hgnc.columns:
        hgnc["name"] = ""

    hgnc["ensembl_gene_id_clean"] = hgnc["ensembl_gene_id"].map(clean_ensembl_gene_id)
    hgnc = hgnc.drop_duplicates(subset=["gene_symbol"], keep="first").copy()

    if limit_genes and limit_genes > 0:
        hgnc = hgnc.head(limit_genes).copy()

    log(f"[HGNC ROWS] {before} -> {len(hgnc)}")
    log(f"[PROTEIN CODING ONLY] {protein_coding_only}")
    log(f"[GENES WITH UNIPROT] {(hgnc['uniprot_ids'].astype(str).str.len() > 0).sum()}")
    log(f"[GENES WITH ENSEMBL] {(hgnc['ensembl_gene_id_clean'].astype(str).str.len() > 0).sum()}")

    return hgnc


def build_maps(hgnc: pd.DataFrame) -> Tuple[Dict[str, Set[str]], Dict[str, str], Dict[str, str]]:
    """
    Returns:
        uniprot_to_genes
        ensembl_to_gene
        symbol_to_gene
    """
    uniprot_to_genes: Dict[str, Set[str]] = defaultdict(set)
    ensembl_to_gene: Dict[str, str] = {}
    symbol_to_gene: Dict[str, str] = {}

    for _, row in hgnc.iterrows():
        gene = row["gene_symbol"]
        symbol_to_gene[gene] = gene

        ens = clean_ensembl_gene_id(row.get("ensembl_gene_id_clean", ""))
        if ens:
            ensembl_to_gene[ens] = gene

        for acc in split_uniprot_ids(row.get("uniprot_ids", "")):
            uniprot_to_genes[acc].add(gene)

    log("=" * 100)
    log("[ID MAPS]")
    log(f"[UNIPROT IDS] {len(uniprot_to_genes)}")
    log(f"[ENSEMBL IDS] {len(ensembl_to_gene)}")
    log(f"[SYMBOLS]     {len(symbol_to_gene)}")

    return uniprot_to_genes, ensembl_to_gene, symbol_to_gene


# =============================================================================
# FILE DETECTION
# =============================================================================

def find_disprot_file(dbdir: Path) -> Optional[Path]:
    roots = [
        dbdir / "DisProt",
        dbdir / "disprot",
        dbdir / "Disorder",
        dbdir,
    ]
    return find_first_existing_file(
        roots,
        patterns=[
            "*disprot*.tsv",
            "*disprot*.tsv.gz",
            "*disprot*.txt",
            "*disprot*.txt.gz",
            "*DisProt*.tsv",
            "*DisProt*.tsv.gz",
        ],
        prefer_keywords=["disprot", "download", "current"],
    )


def find_mobidb_file(dbdir: Path) -> Optional[Path]:
    roots = [
        dbdir / "MobiDB",
        dbdir / "mobidb",
        dbdir / "Disorder",
        dbdir,
    ]
    return find_first_existing_file(
        roots,
        patterns=[
            "*mobidb*.tsv",
            "*mobidb*.tsv.gz",
            "*mobidb*.txt",
            "*mobidb*.txt.gz",
            "*MobiDB*.tsv",
            "*MobiDB*.tsv.gz",
        ],
        prefer_keywords=["mobidb", "human", "9606"],
    )


def find_gnomad_constraint_file(dbdir: Path) -> Optional[Path]:
    roots = [
        dbdir / "gnomAD",
        dbdir / "GnomAD",
        dbdir / "gnomad",
        dbdir / "Constraint",
        dbdir / "constraint",
        dbdir,
    ]
    return find_first_existing_file(
        roots,
        patterns=[
            "*gnomad*constraint*.tsv",
            "*gnomad*constraint*.tsv.gz",
            "*gnomad*constraint*.txt",
            "*gnomad*constraint*.txt.gz",
            "*gnomad*lof_metrics*.txt",
            "*gnomad*lof_metrics*.txt.gz",
            "*gnomad*lof_metrics*.bgz",
            "*lof_metrics*.txt",
            "*lof_metrics*.txt.gz",
            "*lof_metrics*.bgz",
            "*constraint*.txt",
            "*constraint*.txt.gz",
            "*constraint*.tsv",
            "*constraint*.tsv.gz",
        ],
        prefer_keywords=["gnomad", "lof", "constraint", "gene"],
    )


# =============================================================================
# DISPROT PARSING
# =============================================================================

def find_col(df: pd.DataFrame, candidates: List[str], contains: Optional[List[str]] = None) -> Optional[str]:
    lower_map = {str(c).lower().strip(): c for c in df.columns}

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


def parse_disprot(
    disprot_file: Optional[Path],
    uniprot_to_genes: Dict[str, Set[str]],
    symbol_to_gene: Dict[str, str],
    processed_dir: Path,
) -> pd.DataFrame:
    outpath = processed_dir / "feature9_disprot_long.csv"

    if disprot_file is None or not disprot_file.exists():
        log("[DISPROT] No local DisProt file found. Writing empty table.")
        out = pd.DataFrame(columns=["gene_symbol"])
        out.to_csv(outpath, index=False)
        return out

    log("=" * 100)
    log(f"[READ DISPROT] {disprot_file}")

    try:
        df = read_table_flexible(disprot_file)
    except Exception as exc:
        log(f"[DISPROT ERROR] Could not read file: {exc}")
        out = pd.DataFrame(columns=["gene_symbol"])
        out.to_csv(outpath, index=False)
        return out

    df = standardize_columns(df)

    acc_col = find_col(
        df,
        ["UniProt ACC", "UniProt", "uniprot", "accession", "acc", "uniprot_acc", "UniProt accession"],
        contains=["uniprot"],
    )
    gene_col = find_col(
        df,
        ["Gene name", "gene_name", "Gene", "gene", "symbol"],
        contains=["gene"],
    )
    seq_len_col = find_col(
        df,
        ["Sequence length", "sequence_length", "length", "Protein length"],
        contains=["length"],
    )
    disorder_content_col = find_col(
        df,
        ["Protein Disorder Content", "disorder_content", "Disorder content", "protein_disorder_content"],
        contains=["disorder", "content"],
    )
    start_col = find_col(
        df,
        ["Start", "start", "region_start"],
        contains=["start"],
    )
    end_col = find_col(
        df,
        ["End", "end", "region_end"],
        contains=["end"],
    )
    region_id_col = find_col(
        df,
        ["Region ID", "region_id", "Region"],
        contains=["region"],
    )
    term_col = find_col(
        df,
        ["Term name", "term_name", "Term", "term"],
        contains=["term"],
    )
    af_vlow_col = find_col(
        df,
        ["Alphafold Very Low confidence content", "AlphaFold Very Low confidence content"],
        contains=["alphafold", "low"],
    )

    log(f"[DISPROT COL] accession={acc_col}")
    log(f"[DISPROT COL] gene={gene_col}")
    log(f"[DISPROT COL] seq_len={seq_len_col}")
    log(f"[DISPROT COL] disorder_content={disorder_content_col}")
    log(f"[DISPROT COL] start={start_col}")
    log(f"[DISPROT COL] end={end_col}")

    rows = []

    for _, row in df.iterrows():
        genes = set()

        acc = clean_uniprot(row.get(acc_col, "")) if acc_col else ""
        if acc and acc in uniprot_to_genes:
            genes.update(uniprot_to_genes[acc])

        sym = normalize_symbol(row.get(gene_col, "")) if gene_col else ""
        if sym and sym in symbol_to_gene:
            genes.add(sym)

        if not genes:
            continue

        start = safe_int(row.get(start_col, 0)) if start_col else 0
        end = safe_int(row.get(end_col, 0)) if end_col else 0
        region_len = abs(end - start) + 1 if start > 0 and end > 0 else np.nan

        for gene in genes:
            rows.append(
                {
                    "gene_symbol": gene,
                    "uniprot_accession": acc,
                    "disprot_gene_name": sym,
                    "disprot_sequence_length": safe_float(row.get(seq_len_col, np.nan)) if seq_len_col else np.nan,
                    "disprot_disorder_content": safe_float(row.get(disorder_content_col, np.nan)) if disorder_content_col else np.nan,
                    "disprot_region_start": start,
                    "disprot_region_end": end,
                    "disprot_region_length": region_len,
                    "disprot_region_id": clean_text(row.get(region_id_col, "")) if region_id_col else "",
                    "disprot_term_name": clean_text(row.get(term_col, "")) if term_col else "",
                    "disprot_alphafold_very_low_conf_content": safe_float(row.get(af_vlow_col, np.nan)) if af_vlow_col else np.nan,
                }
            )

    out = pd.DataFrame(rows)

    if not out.empty:
        out = out.drop_duplicates()

    out.to_csv(outpath, index=False)

    log(f"[DISPROT MATCHED ROWS] {out.shape[0]}")
    log(f"[SAVED] {outpath}")

    return out


def aggregate_disprot(disprot_long: pd.DataFrame, genes: List[str]) -> pd.DataFrame:
    records = []
    by_gene = dict(tuple(disprot_long.groupby("gene_symbol"))) if not disprot_long.empty else {}

    for gene in genes:
        g = by_gene.get(gene, pd.DataFrame())

        if g.empty:
            records.append(
                {
                    "gene_symbol": gene,
                    "feature9_disprot_has_annotation": 0,
                    "feature9_disprot_region_count": 0,
                }
            )
            continue

        intervals = []
        for _, row in g.iterrows():
            s = safe_int(row.get("disprot_region_start", 0))
            e = safe_int(row.get("disprot_region_end", 0))
            if s > 0 and e > 0:
                intervals.append((s, e))

        region_lengths = pd.to_numeric(g.get("disprot_region_length", pd.Series(dtype=float)), errors="coerce").dropna()
        disorder_content = pd.to_numeric(g.get("disprot_disorder_content", pd.Series(dtype=float)), errors="coerce").dropna()
        seq_len = pd.to_numeric(g.get("disprot_sequence_length", pd.Series(dtype=float)), errors="coerce").dropna()
        af_vlow = pd.to_numeric(g.get("disprot_alphafold_very_low_conf_content", pd.Series(dtype=float)), errors="coerce").dropna()

        union_len = interval_union_length(intervals)
        max_seq_len = float(seq_len.max()) if len(seq_len) else np.nan
        region_fraction = float(union_len / max_seq_len) if max_seq_len and max_seq_len > 0 else np.nan

        records.append(
            {
                "gene_symbol": gene,
                "feature9_disprot_has_annotation": 1,
                "feature9_disprot_row_count": int(len(g)),
                "feature9_disprot_region_count": int(len(intervals)),
                "feature9_disprot_unique_uniprot_count": int(g["uniprot_accession"].replace("", np.nan).dropna().nunique()) if "uniprot_accession" in g.columns else 0,
                "feature9_disprot_region_union_length": int(union_len),
                "feature9_disprot_longest_region_length": float(region_lengths.max()) if len(region_lengths) else np.nan,
                "feature9_disprot_mean_region_length": float(region_lengths.mean()) if len(region_lengths) else np.nan,
                "feature9_disprot_median_region_length": float(region_lengths.median()) if len(region_lengths) else np.nan,
                "feature9_disprot_max_sequence_length": max_seq_len,
                "feature9_disprot_region_fraction_of_protein": region_fraction,
                "feature9_disprot_disorder_content_max": float(disorder_content.max()) if len(disorder_content) else np.nan,
                "feature9_disprot_disorder_content_mean": float(disorder_content.mean()) if len(disorder_content) else np.nan,
                "feature9_disprot_alphafold_very_low_conf_content_max": float(af_vlow.max()) if len(af_vlow) else np.nan,
                "feature9_disprot_alphafold_very_low_conf_content_mean": float(af_vlow.mean()) if len(af_vlow) else np.nan,
            }
        )

    return pd.DataFrame(records)


# =============================================================================
# MOBIDB PARSING
# =============================================================================

def parse_mobidb(
    mobidb_file: Optional[Path],
    uniprot_to_genes: Dict[str, Set[str]],
    symbol_to_gene: Dict[str, str],
    processed_dir: Path,
) -> pd.DataFrame:
    outpath = processed_dir / "feature9_mobidb_long.csv"

    if mobidb_file is None or not mobidb_file.exists():
        log("[MOBIDB] No local MobiDB file found. Writing empty table.")
        out = pd.DataFrame(columns=["gene_symbol"])
        out.to_csv(outpath, index=False)
        return out

    log("=" * 100)
    log(f"[READ MOBIDB] {mobidb_file}")

    try:
        df = read_table_flexible(mobidb_file)
    except Exception as exc:
        log(f"[MOBIDB ERROR] Could not read file: {exc}")
        out = pd.DataFrame(columns=["gene_symbol"])
        out.to_csv(outpath, index=False)
        return out

    df = standardize_columns(df)

    acc_col = find_col(
        df,
        ["acc", "accession", "uniprot", "uniprot_acc", "UniProt", "UniProt ACC"],
        contains=["acc"],
    )
    gene_col = find_col(
        df,
        ["gene", "gene_name", "Gene", "Gene name", "symbol"],
        contains=["gene"],
    )

    disorder_cols = []
    for c in df.columns:
        low = str(c).lower()
        if "disorder" in low or "idr" in low:
            disorder_cols.append(c)

    start_col = find_col(df, ["start", "region_start", "Start"], contains=["start"])
    end_col = find_col(df, ["end", "region_end", "End"], contains=["end"])
    length_col = find_col(df, ["length", "sequence_length", "Sequence length"], contains=["length"])

    log(f"[MOBIDB COL] accession={acc_col}")
    log(f"[MOBIDB COL] gene={gene_col}")
    log(f"[MOBIDB COL] disorder_cols={disorder_cols[:10]}")

    rows = []

    for _, row in df.iterrows():
        genes = set()

        acc = clean_uniprot(row.get(acc_col, "")) if acc_col else ""
        if acc and acc in uniprot_to_genes:
            genes.update(uniprot_to_genes[acc])

        sym = normalize_symbol(row.get(gene_col, "")) if gene_col else ""
        if sym and sym in symbol_to_gene:
            genes.add(sym)

        if not genes:
            continue

        start = safe_int(row.get(start_col, 0)) if start_col else 0
        end = safe_int(row.get(end_col, 0)) if end_col else 0
        region_len = abs(end - start) + 1 if start > 0 and end > 0 else np.nan

        disorder_numeric_values = []
        for c in disorder_cols:
            v = safe_float(row.get(c, np.nan))
            if not pd.isna(v):
                disorder_numeric_values.append(v)

        disorder_max = max(disorder_numeric_values) if disorder_numeric_values else np.nan
        disorder_mean = float(np.mean(disorder_numeric_values)) if disorder_numeric_values else np.nan

        for gene in genes:
            rows.append(
                {
                    "gene_symbol": gene,
                    "uniprot_accession": acc,
                    "mobidb_gene_name": sym,
                    "mobidb_sequence_length": safe_float(row.get(length_col, np.nan)) if length_col else np.nan,
                    "mobidb_region_start": start,
                    "mobidb_region_end": end,
                    "mobidb_region_length": region_len,
                    "mobidb_disorder_numeric_max": disorder_max,
                    "mobidb_disorder_numeric_mean": disorder_mean,
                }
            )

    out = pd.DataFrame(rows)

    if not out.empty:
        out = out.drop_duplicates()

    out.to_csv(outpath, index=False)

    log(f"[MOBIDB MATCHED ROWS] {out.shape[0]}")
    log(f"[SAVED] {outpath}")

    return out


def aggregate_mobidb(mobidb_long: pd.DataFrame, genes: List[str]) -> pd.DataFrame:
    records = []
    by_gene = dict(tuple(mobidb_long.groupby("gene_symbol"))) if not mobidb_long.empty else {}

    for gene in genes:
        g = by_gene.get(gene, pd.DataFrame())

        if g.empty:
            records.append(
                {
                    "gene_symbol": gene,
                    "feature9_mobidb_has_annotation": 0,
                    "feature9_mobidb_region_count": 0,
                }
            )
            continue

        intervals = []
        for _, row in g.iterrows():
            s = safe_int(row.get("mobidb_region_start", 0))
            e = safe_int(row.get("mobidb_region_end", 0))
            if s > 0 and e > 0:
                intervals.append((s, e))

        region_lengths = pd.to_numeric(g.get("mobidb_region_length", pd.Series(dtype=float)), errors="coerce").dropna()
        disorder_max = pd.to_numeric(g.get("mobidb_disorder_numeric_max", pd.Series(dtype=float)), errors="coerce").dropna()
        disorder_mean = pd.to_numeric(g.get("mobidb_disorder_numeric_mean", pd.Series(dtype=float)), errors="coerce").dropna()
        seq_len = pd.to_numeric(g.get("mobidb_sequence_length", pd.Series(dtype=float)), errors="coerce").dropna()

        union_len = interval_union_length(intervals)
        max_seq_len = float(seq_len.max()) if len(seq_len) else np.nan
        region_fraction = float(union_len / max_seq_len) if max_seq_len and max_seq_len > 0 else np.nan

        records.append(
            {
                "gene_symbol": gene,
                "feature9_mobidb_has_annotation": 1,
                "feature9_mobidb_row_count": int(len(g)),
                "feature9_mobidb_region_count": int(len(intervals)),
                "feature9_mobidb_unique_uniprot_count": int(g["uniprot_accession"].replace("", np.nan).dropna().nunique()) if "uniprot_accession" in g.columns else 0,
                "feature9_mobidb_region_union_length": int(union_len),
                "feature9_mobidb_longest_region_length": float(region_lengths.max()) if len(region_lengths) else np.nan,
                "feature9_mobidb_mean_region_length": float(region_lengths.mean()) if len(region_lengths) else np.nan,
                "feature9_mobidb_region_fraction_of_protein": region_fraction,
                "feature9_mobidb_disorder_numeric_max": float(disorder_max.max()) if len(disorder_max) else np.nan,
                "feature9_mobidb_disorder_numeric_mean": float(disorder_mean.mean()) if len(disorder_mean) else np.nan,
            }
        )

    return pd.DataFrame(records)


# =============================================================================
# GNOMAD CONSTRAINT PARSING
# =============================================================================

def parse_gnomad_constraint(
    constraint_file: Optional[Path],
    hgnc: pd.DataFrame,
    ensembl_to_gene: Dict[str, str],
    symbol_to_gene: Dict[str, str],
    processed_dir: Path,
) -> pd.DataFrame:
    outpath = processed_dir / "feature9_constraint_gene_table.csv"

    if constraint_file is None or not constraint_file.exists():
        log("[GNOMAD] No local gnomAD constraint file found. Writing empty table.")
        out = pd.DataFrame(columns=["gene_symbol"])
        out.to_csv(outpath, index=False)
        return out

    log("=" * 100)
    log(f"[READ GNOMAD CONSTRAINT] {constraint_file}")

    try:
        df = read_table_flexible(constraint_file)
    except Exception as exc:
        log(f"[GNOMAD ERROR] Could not read file: {exc}")
        out = pd.DataFrame(columns=["gene_symbol"])
        out.to_csv(outpath, index=False)
        return out

    df = standardize_columns(df)

    gene_col = find_col(
        df,
        ["gene", "gene_symbol", "symbol", "Gene", "Gene Symbol"],
        contains=["gene"],
    )
    ens_col = find_col(
        df,
        ["gene_id", "ensembl_gene_id", "ensg", "transcript"],
        contains=["gene", "id"],
    )

    log(f"[GNOMAD COL] gene={gene_col}")
    log(f"[GNOMAD COL] ensembl={ens_col}")

    # Candidate constraint columns. The script only keeps columns that exist.
    candidate_numeric_cols = [
        "pLI",
        "pli",
        "pRec",
        "prec",
        "pNull",
        "pnull",
        "oe_lof",
        "oe_lof_lower",
        "oe_lof_upper",
        "oe_lof_upper_rank",
        "oe_mis",
        "oe_mis_lower",
        "oe_mis_upper",
        "oe_syn",
        "oe_syn_lower",
        "oe_syn_upper",
        "lof.oe",
        "lof.oe_ci.lower",
        "lof.oe_ci.upper",
        "mis.oe",
        "mis.oe_ci.lower",
        "mis.oe_ci.upper",
        "syn.oe",
        "syn.oe_ci.lower",
        "syn.oe_ci.upper",
        "obs_lof",
        "exp_lof",
        "obs_mis",
        "exp_mis",
        "obs_syn",
        "exp_syn",
        "lof_z",
        "mis_z",
        "syn_z",
        "constraint_flag",
    ]

    existing_numeric_cols = []
    lower_to_col = {c.lower(): c for c in df.columns}

    for c in candidate_numeric_cols:
        if c.lower() in lower_to_col:
            existing_numeric_cols.append(lower_to_col[c.lower()])

    # Also collect LOEUF if present under any name.
    for c in df.columns:
        low = c.lower()
        if "loeuf" in low and c not in existing_numeric_cols:
            existing_numeric_cols.append(c)
        if "constraint" in low and c not in existing_numeric_cols:
            # only numeric later
            existing_numeric_cols.append(c)

    rows = []

    for _, row in df.iterrows():
        gene = ""

        sym = normalize_symbol(row.get(gene_col, "")) if gene_col else ""
        ens = clean_ensembl_gene_id(row.get(ens_col, "")) if ens_col else ""

        if sym in symbol_to_gene:
            gene = sym
        elif ens in ensembl_to_gene:
            gene = ensembl_to_gene[ens]

        if not gene:
            continue

        rec = {"gene_symbol": gene}

        for c in existing_numeric_cols:
            safe_name = re.sub(r"[^A-Za-z0-9]+", "_", c).strip("_").lower()
            rec[f"feature9_gnomad_{safe_name}"] = safe_float(row.get(c, np.nan))

        rows.append(rec)

    out = pd.DataFrame(rows)

    if not out.empty:
        # Merge duplicate rows by gene using median for numeric columns.
        numeric_cols = [c for c in out.columns if c != "gene_symbol"]
        for c in numeric_cols:
            out[c] = pd.to_numeric(out[c], errors="coerce")

        out = out.groupby("gene_symbol", as_index=False)[numeric_cols].median()

        # Derived standard constraint features.
        out["feature9_gnomad_has_constraint"] = 1

        # pLI flag.
        pli_col = None
        for c in out.columns:
            if c.lower() in {"feature9_gnomad_pli", "feature9_gnomad_pli"} or c.lower().endswith("_pli"):
                pli_col = c
                break

        if pli_col:
            out["feature9_gnomad_pli_ge_0_9"] = (pd.to_numeric(out[pli_col], errors="coerce") >= 0.9).astype(int)
            out["feature9_gnomad_pli_ge_0_99"] = (pd.to_numeric(out[pli_col], errors="coerce") >= 0.99).astype(int)

        # LOEUF flag: lower = more constrained.
        loeuf_col = None
        for c in out.columns:
            if "loeuf" in c.lower() or "oe_lof_upper" in c.lower() or "lof_oe_ci_upper" in c.lower():
                loeuf_col = c
                break

        if loeuf_col:
            lo = pd.to_numeric(out[loeuf_col], errors="coerce")
            out["feature9_gnomad_low_loeuf_le_0_35"] = (lo <= 0.35).astype(int)
            out["feature9_gnomad_low_loeuf_le_0_60"] = (lo <= 0.60).astype(int)
            out["feature9_gnomad_high_loeuf_ge_1"] = (lo >= 1.0).astype(int)
    else:
        out = pd.DataFrame(columns=["gene_symbol"])

    out.to_csv(outpath, index=False)

    log(f"[GNOMAD MATCHED GENES] {out.shape[0]}")
    log(f"[SAVED] {outpath}")

    return out


def add_empty_constraint_features(genes: List[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": genes,
            "feature9_gnomad_has_constraint": 0,
        }
    )


# =============================================================================
# ALPHAFOLD / FEATURE4 PROXY
# =============================================================================

def extract_feature4_alphafold_proxy(feature4_path: Path, genes: List[str]) -> pd.DataFrame:
    """
    Extract existing AlphaFold confidence features from Feature 4 if available.

    This is useful as a backup disorder proxy:
        lower pLDDT = more flexible/disordered/uncertain regions.

    The function is intentionally flexible and looks for feature4 columns that contain:
        alphafold, plddt, confidence, bfactor, radius, compactness
    """
    base = pd.DataFrame({"gene_symbol": genes})

    if not feature4_path.exists():
        log(f"[FEATURE4] Not found: {feature4_path}")
        base["feature9_afproxy_has_feature4"] = 0
        return base

    log("=" * 100)
    log(f"[READ FEATURE4 PROXY] {feature4_path}")

    df = pd.read_csv(feature4_path, low_memory=False)

    if "gene_symbol" not in df.columns:
        if "symbol" in df.columns:
            df["gene_symbol"] = df["symbol"].map(normalize_symbol)
        else:
            log("[FEATURE4] No gene_symbol/symbol column. Skipping.")
            base["feature9_afproxy_has_feature4"] = 0
            return base

    df["gene_symbol"] = df["gene_symbol"].map(normalize_symbol)

    keep = ["gene_symbol"]
    for c in df.columns:
        low = c.lower()
        if not c.startswith("feature4_"):
            continue

        # Keep numeric AlphaFold/structure-confidence proxies only.
        if any(k in low for k in ["alphafold", "plddt", "confidence", "bfactor", "radius", "compactness", "gyration"]):
            keep.append(c)

    keep = list(dict.fromkeys(keep))

    if len(keep) == 1:
        base["feature9_afproxy_has_feature4"] = 0
        return base

    sub = df[keep].drop_duplicates(subset=["gene_symbol"], keep="first").copy()

    rename = {}
    for c in sub.columns:
        if c == "gene_symbol":
            continue
        rename[c] = "feature9_afproxy_" + c.replace("feature4_", "")

    sub = sub.rename(columns=rename)

    for c in sub.columns:
        if c != "gene_symbol":
            sub[c] = pd.to_numeric(sub[c], errors="coerce")

    sub["feature9_afproxy_has_feature4"] = 1

    out = base.merge(sub, on="gene_symbol", how="left")
    out["feature9_afproxy_has_feature4"] = out["feature9_afproxy_has_feature4"].fillna(0).astype(int)

    # Try to derive low-confidence flags if mean pLDDT is present.
    plddt_cols = [c for c in out.columns if "plddt" in c.lower() and pd.api.types.is_numeric_dtype(out[c])]
    mean_cols = [c for c in plddt_cols if "mean" in c.lower() or "avg" in c.lower()]

    if mean_cols:
        c = mean_cols[0]
        out["feature9_afproxy_mean_plddt_lt_70"] = (pd.to_numeric(out[c], errors="coerce") < 70).astype(int)
        out["feature9_afproxy_mean_plddt_lt_50"] = (pd.to_numeric(out[c], errors="coerce") < 50).astype(int)

    return out


# =============================================================================
# FEATURE MERGING
# =============================================================================

def merge_feature_blocks(
    hgnc: pd.DataFrame,
    disprot_features: pd.DataFrame,
    mobidb_features: pd.DataFrame,
    constraint_features: pd.DataFrame,
    afproxy_features: pd.DataFrame,
    processed_dir: Path,
) -> pd.DataFrame:
    genes = pd.DataFrame({"gene_symbol": sorted(hgnc["gene_symbol"].unique().tolist())})

    merged = genes.copy()

    for block in [disprot_features, mobidb_features, constraint_features, afproxy_features]:
        if block is None or block.empty or "gene_symbol" not in block.columns:
            continue
        merged = merged.merge(block, on="gene_symbol", how="left")

    # Add combined disorder features.
    disorder_has_cols = [
        c for c in merged.columns
        if c in {
            "feature9_disprot_has_annotation",
            "feature9_mobidb_has_annotation",
            "feature9_afproxy_has_feature4",
        }
    ]

    if disorder_has_cols:
        tmp = merged[disorder_has_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        merged["feature9_any_disorder_source_available"] = (tmp.sum(axis=1) > 0).astype(int)
    else:
        merged["feature9_any_disorder_source_available"] = 0

    # Combined high-disorder proxies.
    candidate_disorder_fraction_cols = [
        c for c in merged.columns
        if any(x in c.lower() for x in ["disorder_content", "region_fraction", "very_low_conf"])
    ]

    numeric_fraction_cols = []
    for c in candidate_disorder_fraction_cols:
        merged[c] = pd.to_numeric(merged[c], errors="coerce")
        if merged[c].notna().sum() > 0:
            numeric_fraction_cols.append(c)

    if numeric_fraction_cols:
        merged["feature9_combined_max_disorder_fraction_proxy"] = merged[numeric_fraction_cols].max(axis=1, skipna=True)
        merged["feature9_combined_mean_disorder_fraction_proxy"] = merged[numeric_fraction_cols].mean(axis=1, skipna=True)
        merged["feature9_combined_high_disorder_proxy_ge_0_30"] = (
            merged["feature9_combined_max_disorder_fraction_proxy"] >= 0.30
        ).astype(int)
        merged["feature9_combined_high_disorder_proxy_ge_0_50"] = (
            merged["feature9_combined_max_disorder_fraction_proxy"] >= 0.50
        ).astype(int)
    else:
        merged["feature9_combined_max_disorder_fraction_proxy"] = np.nan
        merged["feature9_combined_mean_disorder_fraction_proxy"] = np.nan
        merged["feature9_combined_high_disorder_proxy_ge_0_30"] = 0
        merged["feature9_combined_high_disorder_proxy_ge_0_50"] = 0

    # Combined constraint flags.
    constraint_cols = [c for c in merged.columns if c.startswith("feature9_gnomad_")]
    if constraint_cols:
        merged["feature9_any_constraint_feature_available"] = (
            merged[constraint_cols].notna().sum(axis=1) > 0
        ).astype(int)
    else:
        merged["feature9_any_constraint_feature_available"] = 0

    # Fill count/flag columns only.
    for c in merged.columns:
        if c.startswith("feature9_") and (
            "_count" in c
            or "_has_" in c
            or c.startswith("feature9_any_")
            or "_ge_" in c
            or "_le_" in c
            or c.endswith("_available")
        ):
            if merged[c].dtype != object:
                merged[c] = merged[c].fillna(0)

    outpath = processed_dir / "feature9_disorder_constraint_gene_features.csv"
    merged.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED FEATURE9 GENE FEATURES]")
    log(f"[PATH]  {outpath}")
    log(f"[SHAPE] {merged.shape}")

    return merged


def merge_with_hgnc(hgnc: pd.DataFrame, features: pd.DataFrame, processed_dir: Path) -> pd.DataFrame:
    merged = hgnc.merge(features, on="gene_symbol", how="left")

    feature_cols = [c for c in features.columns if c != "gene_symbol"]
    merged["feature9_disorder_constraint_has_any_feature"] = merged[feature_cols].notna().any(axis=1).astype(int)
    merged["feature9_disorder_constraint_n_nonmissing_features"] = merged[feature_cols].notna().sum(axis=1)

    outpath = processed_dir / "feature9_disorder_constraint_hgnc_merged.csv"
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
    disprot_file: Optional[Path],
    mobidb_file: Optional[Path],
    constraint_file: Optional[Path],
    feature4_path: Path,
    disprot_long: pd.DataFrame,
    mobidb_long: pd.DataFrame,
    constraint_features: pd.DataFrame,
    features: pd.DataFrame,
    merged: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    path = outdir / "feature9_disorder_constraint_summary.txt"

    lines = []
    lines.append("Feature 9: Protein disorder + gene constraint / intolerance")
    lines.append("=" * 80)
    lines.append(f"Created at: {now_iso()}")
    lines.append("")
    lines.append("Inputs:")
    lines.append(f"HGNC: {args.hgnc}")
    lines.append(f"DisProt file: {disprot_file if disprot_file else 'not found / skipped'}")
    lines.append(f"MobiDB file: {mobidb_file if mobidb_file else 'not found / skipped'}")
    lines.append(f"gnomAD constraint file: {constraint_file if constraint_file else 'not found / skipped'}")
    lines.append(f"Feature4 AlphaFold proxy file: {feature4_path if feature4_path.exists() else 'not found / skipped'}")
    lines.append("")
    lines.append("Rows:")
    lines.append(f"HGNC rows: {hgnc.shape[0]}")
    lines.append(f"DisProt long rows: {disprot_long.shape[0]}")
    lines.append(f"MobiDB long rows: {mobidb_long.shape[0]}")
    lines.append(f"gnomAD constraint matched genes: {constraint_features.shape[0] if not constraint_features.empty else 0}")
    lines.append(f"Feature9 gene feature rows: {features.shape[0]}")
    lines.append(f"Feature9 gene feature columns: {features.shape[1]}")
    lines.append(f"Merged rows: {merged.shape[0]}")
    lines.append(f"Merged columns: {merged.shape[1]}")
    lines.append("")
    lines.append("Coverage:")

    for col, label in [
        ("feature9_disprot_has_annotation", "DisProt annotation coverage"),
        ("feature9_mobidb_has_annotation", "MobiDB annotation coverage"),
        ("feature9_gnomad_has_constraint", "gnomAD constraint coverage"),
        ("feature9_afproxy_has_feature4", "AlphaFold/Feature4 proxy coverage"),
        ("feature9_any_disorder_source_available", "Any disorder proxy coverage"),
        ("feature9_any_constraint_feature_available", "Any constraint feature coverage"),
    ]:
        if col in features.columns:
            lines.append(f"{label}: {features[col].fillna(0).mean():.4f}")

    if "feature9_combined_max_disorder_fraction_proxy" in features.columns:
        lines.append(
            "Median combined disorder fraction proxy: "
            f"{features['feature9_combined_max_disorder_fraction_proxy'].median(skipna=True):.4f}"
        )

    lines.append("")
    lines.append("Leakage policy:")
    lines.append("Included: DisProt/MobiDB disorder, AlphaFold confidence proxy, gnomAD gene constraint.")
    lines.append("Excluded: clinical target labels, ChEMBL, DrugBank, DGIdb, Open Targets knownDrugs/tractability, Pharos TDL.")
    lines.append("")
    lines.append("Interpretation:")
    lines.append("High disorder means the protein may contain flexible/intrinsically disordered regions.")
    lines.append("Low LOEUF / high pLI means the gene is more loss-of-function constrained.")
    lines.append("These are not drug labels; they describe protein flexibility and genetic intolerance.")

    path.write_text("\n".join(lines) + "\n")
    log(f"[SUMMARY] {path}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description="Build Feature 9 disorder + constraint features.")
    parser.add_argument("--hgnc", default=str(DEFAULT_HGNC_PATH), help="Path to HGNC complete set.")
    parser.add_argument("--dbdir", default=str(DEFAULT_DBDIR), help="Root feature_databases directory.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory.")
    parser.add_argument("--feature4", default=str(DEFAULT_FEATURE4_PATH), help="Feature4 structure HGNC merged table.")
    parser.add_argument("--disprot-file", default="", help="Manual DisProt TSV/TXT file.")
    parser.add_argument("--mobidb-file", default="", help="Manual MobiDB TSV/TXT file.")
    parser.add_argument("--gnomad-file", default="", help="Manual gnomAD constraint file.")
    parser.add_argument("--all-hgnc-genes", action="store_true", help="Use all HGNC genes, not only protein-coding.")
    parser.add_argument("--limit-genes", type=int, default=0, help="Debug only: first N HGNC genes.")
    parser.add_argument("--download", action="store_true", help="Try to download gnomAD constraint and DisProt into feature_databases.")
    args = parser.parse_args()

    outdir = mkdir(Path(args.outdir))
    processed_dir = mkdir(outdir / "processed")
    dbdir = Path(args.dbdir)
    feature4_path = Path(args.feature4)

    log("=" * 100)
    log("FEATURE 9: PROTEIN DISORDER + GNOMAD CONSTRAINT")
    log("=" * 100)
    log(f"[HGNC]       {args.hgnc}")
    log(f"[DBDIR]      {dbdir.resolve()}")
    log(f"[OUTDIR]     {outdir.resolve()}")
    log(f"[FEATURE4]   {feature4_path}")
    log(f"[DOWNLOAD]   {args.download}")
    log(f"[LIMIT]      {args.limit_genes if args.limit_genes else 'none'}")
    log("=" * 100)

    if args.download:
        maybe_download_databases(dbdir)

    hgnc = load_hgnc(
        Path(args.hgnc),
        protein_coding_only=not args.all_hgnc_genes,
        limit_genes=args.limit_genes,
    )

    genes = sorted(hgnc["gene_symbol"].unique().tolist())
    uniprot_to_genes, ensembl_to_gene, symbol_to_gene = build_maps(hgnc)

    disprot_file = Path(args.disprot_file) if args.disprot_file else find_disprot_file(dbdir)
    mobidb_file = Path(args.mobidb_file) if args.mobidb_file else find_mobidb_file(dbdir)
    constraint_file = Path(args.gnomad_file) if args.gnomad_file else find_gnomad_constraint_file(dbdir)

    log("=" * 100)
    log("[DETECTED FILES]")
    log(f"[DISPROT] {disprot_file if disprot_file else 'not found'}")
    log(f"[MOBIDB]  {mobidb_file if mobidb_file else 'not found'}")
    log(f"[GNOMAD]  {constraint_file if constraint_file else 'not found'}")
    log("=" * 100)

    disprot_long = parse_disprot(
        disprot_file=disprot_file,
        uniprot_to_genes=uniprot_to_genes,
        symbol_to_gene=symbol_to_gene,
        processed_dir=processed_dir,
    )

    mobidb_long = parse_mobidb(
        mobidb_file=mobidb_file,
        uniprot_to_genes=uniprot_to_genes,
        symbol_to_gene=symbol_to_gene,
        processed_dir=processed_dir,
    )

    disprot_features = aggregate_disprot(disprot_long, genes)
    mobidb_features = aggregate_mobidb(mobidb_long, genes)

    constraint_features = parse_gnomad_constraint(
        constraint_file=constraint_file,
        hgnc=hgnc,
        ensembl_to_gene=ensembl_to_gene,
        symbol_to_gene=symbol_to_gene,
        processed_dir=processed_dir,
    )

    if constraint_features.empty or "gene_symbol" not in constraint_features.columns:
        constraint_features = add_empty_constraint_features(genes)

    afproxy_features = extract_feature4_alphafold_proxy(
        feature4_path=feature4_path,
        genes=genes,
    )

    features = merge_feature_blocks(
        hgnc=hgnc,
        disprot_features=disprot_features,
        mobidb_features=mobidb_features,
        constraint_features=constraint_features,
        afproxy_features=afproxy_features,
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
        "feature4": str(feature4_path),
        "disprot_file": str(disprot_file) if disprot_file else "",
        "mobidb_file": str(mobidb_file) if mobidb_file else "",
        "gnomad_constraint_file": str(constraint_file) if constraint_file else "",
        "protein_coding_only": not args.all_hgnc_genes,
        "limit_genes": args.limit_genes,
        "download_attempted": args.download,
        "outputs": {
            "disprot_long": str(processed_dir / "feature9_disprot_long.csv"),
            "mobidb_long": str(processed_dir / "feature9_mobidb_long.csv"),
            "constraint_gene_table": str(processed_dir / "feature9_constraint_gene_table.csv"),
            "gene_features": str(processed_dir / "feature9_disorder_constraint_gene_features.csv"),
            "hgnc_merged": str(processed_dir / "feature9_disorder_constraint_hgnc_merged.csv"),
        },
        "leakage_policy": {
            "included": [
                "DisProt disorder annotations",
                "MobiDB disorder annotations",
                "AlphaFold pLDDT / confidence proxy from Feature4",
                "gnomAD constraint metrics",
            ],
            "excluded": [
                "clinical target labels",
                "known drug target labels",
                "ChEMBL",
                "DrugBank",
                "DGIdb",
                "Open Targets knownDrugs",
                "Open Targets tractability",
                "Pharos TDL",
            ],
        },
    }

    with open(outdir / "feature9_disorder_constraint_run_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    write_summary(
        outdir=outdir,
        hgnc=hgnc,
        disprot_file=disprot_file,
        mobidb_file=mobidb_file,
        constraint_file=constraint_file,
        feature4_path=feature4_path,
        disprot_long=disprot_long,
        mobidb_long=mobidb_long,
        constraint_features=constraint_features,
        features=features,
        merged=merged,
        args=args,
    )

    log("=" * 100)
    log("[DONE]")
    log(f"[GENE FEATURES] {processed_dir / 'feature9_disorder_constraint_gene_features.csv'}")
    log(f"[HGNC MERGED]   {processed_dir / 'feature9_disorder_constraint_hgnc_merged.csv'}")
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