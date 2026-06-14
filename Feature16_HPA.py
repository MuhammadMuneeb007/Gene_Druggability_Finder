#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Feature16_HPA.py

Feature 16: Human Protein Atlas (HPA) subcellular location, RNA expression,
protein tissue expression, immune/single-cell expression, and protein atlas
summary features for gene druggability modelling.

Purpose
-------
Build leakage-aware gene-level HPA features that are useful for estimating
whether a gene/protein may be druggable, without using known drug-target labels.

Included feature groups
-----------------------
1. Subcellular localization
   - plasma membrane / membrane / secreted / extracellular flags
   - nuclear, cytosolic, mitochondrial, ER, Golgi, lysosome, vesicle flags
   - location count and reliability score

2. RNA tissue expression
   - consensus tissue RNA expression
   - HPA tissue RNA expression
   - GTEx tissue RNA expression
   - expression max/mean/median/sum
   - detected tissue count at multiple thresholds
   - tissue specificity metrics: tau and entropy
   - tissue-group expression summaries for brain, liver, kidney, heart,
     lung, intestine, immune/blood, reproductive, endocrine, muscle, skin

3. Protein expression in normal tissues by IHC
   - protein expression max/mean/median
   - high/medium/low/not-detected counts
   - tissue and cell-type diversity

4. Immune-cell and single-cell RNA expression when files are available
   - expression max/mean/median
   - detected cell-type count
   - tau/entropy across cell types

5. Protein atlas summary file, if available
   - selected HPA annotation columns converted into safe categorical/count flags

Leakage policy
--------------
Included:
    HPA subcellular localization
    HPA/GTEx/FANTOM-style expression summaries
    normal tissue protein expression
    immune/single-cell expression
    HPA protein class keywords only if present in proteinatlas.tsv

Excluded by default:
    cancer prognostic/survival data
    TCGA survival outcomes
    disease-outcome labels
    drug/chemical interaction files
    known drug-target labels
    ChEMBL, DrugBank, DGIdb, Open Targets, Pharos labels

Default inputs
--------------
HGNC:
    databases/HGNC/hgnc_complete_set.txt

HPA local directory:
    feature_databases/HPA/

Automatically downloaded core files:
    subcellular_location.tsv.zip
    rna_tissue_consensus.tsv.zip
    rna_tissue_hpa.tsv.zip
    rna_tissue_gtex.tsv.zip
    normal_ihc_data.tsv.zip
    rna_immune_cell.tsv.zip
    rna_single_cell_type.tsv.zip
    proteinatlas.tsv.zip

Outputs
-------
feature16_hpa/
    processed/feature16_hpa_gene_features.csv
    processed/feature16_hpa_hgnc_merged.csv
    processed/feature16_hpa_subcellular_features.csv
    processed/feature16_hpa_rna_consensus_features.csv
    processed/feature16_hpa_rna_hpa_features.csv
    processed/feature16_hpa_rna_gtex_features.csv
    processed/feature16_hpa_normal_ihc_features.csv
    processed/feature16_hpa_immune_cell_features.csv
    processed/feature16_hpa_single_cell_features.csv
    processed/feature16_hpa_proteinatlas_features.csv
    feature16_hpa_summary.txt
    feature16_hpa_run_metadata.json

Run
---
    python Feature16_HPA.py

Force download:
    python Feature16_HPA.py --download --force-download

Fast test:
    python Feature16_HPA.py --limit-genes 1000 --max-rows 200000

If HPC blocks download, manually place ZIP files in:
    feature_databases/HPA/
then rerun:
    python Feature16_HPA.py
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
import zipfile
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


# =============================================================================
# DEFAULTS
# =============================================================================

DEFAULT_HGNC_PATH = Path("databases") / "HGNC" / "hgnc_complete_set.txt"
DEFAULT_DBDIR = Path("feature_databases") / "HPA"
DEFAULT_OUTDIR = Path("feature16_hpa")

# Current HPA TSV download base uses /download/tsv/.
HPA_TSV_BASE_URL = "https://www.proteinatlas.org/download/tsv"

# Core, druggability-relevant, leakage-safe files.
# Some HPA filenames can change between versions; downloader also tries aliases.
HPA_FILE_SPECS: Dict[str, Dict[str, Any]] = {
    "subcellular": {
        "filename": "subcellular_location.tsv.zip",
        "aliases": ["subcellular_location.tsv.zip"],
        "required": True,
        "description": "Subcellular protein location per gene",
    },
    "rna_consensus": {
        "filename": "rna_tissue_consensus.tsv.zip",
        "aliases": ["rna_tissue_consensus.tsv.zip"],
        "required": True,
        "description": "Consensus tissue RNA expression per gene",
    },
    "rna_hpa": {
        "filename": "rna_tissue_hpa.tsv.zip",
        "aliases": ["rna_tissue_hpa.tsv.zip"],
        "required": False,
        "description": "HPA tissue RNA expression per gene",
    },
    "rna_gtex": {
        "filename": "rna_tissue_gtex.tsv.zip",
        "aliases": ["rna_tissue_gtex.tsv.zip"],
        "required": False,
        "description": "GTEx tissue RNA expression per gene",
    },
    "normal_ihc": {
        "filename": "normal_ihc_data.tsv.zip",
        "aliases": ["normal_ihc_data.tsv.zip", "normal_tissue.tsv.zip"],
        "required": False,
        "description": "Normal tissue protein expression by IHC",
    },
    "rna_immune_cell": {
        "filename": "rna_immune_cell.tsv.zip",
        "aliases": ["rna_immune_cell.tsv.zip", "rna_blood_cell.tsv.zip"],
        "required": False,
        "description": "Immune-cell RNA expression per gene",
    },
    "rna_single_cell": {
        "filename": "rna_single_cell_type.tsv.zip",
        "aliases": ["rna_single_cell_type.tsv.zip", "rna_single_cell.tsv.zip"],
        "required": False,
        "description": "Single-cell type RNA expression per gene",
    },
    "proteinatlas": {
        "filename": "proteinatlas.tsv.zip",
        "aliases": ["proteinatlas.tsv.zip"],
        "required": False,
        "description": "HPA protein atlas summary table",
    },
}

# Large / potentially leaky disease-outcome files deliberately not downloaded by default.
# They can be enabled with --include-large-risky-hpa, but the generated features are not parsed
# unless you extend the script. Keeping this explicit helps avoid accidental leakage.
HPA_RISKY_OR_LARGE_FILES = {
    "rna_cancer_sample": "rna_cancer_sample.tsv.zip",
    "pathology": "pathology.tsv.zip",
    "cancer_ihc": "pathology_ihc_data.tsv.zip",
}


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
    if x is None:
        return ""
    try:
        if pd.isna(x):
            return ""
    except Exception:
        pass
    s = str(x).strip()
    if s.lower() in {"nan", "none", "na", "null", "-"}:
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
    s = s.replace("gene:", "")
    s = s.split(".")[0]
    return s


def safe_float(x: Any) -> float:
    try:
        if x is None or pd.isna(x):
            return np.nan
        s = str(x).strip()
        if s in {"", ".", "NA", "NaN", "nan", "None", "null", "-"}:
            return np.nan
        return float(s)
    except Exception:
        return np.nan


def safe_int(x: Any) -> int:
    try:
        if x is None or pd.isna(x):
            return 0
        s = str(x).strip()
        if s in {"", ".", "NA", "NaN", "nan", "None", "null", "-"}:
            return 0
        return int(float(s))
    except Exception:
        return 0


def sanitize_feature_name(x: str) -> str:
    x = str(x).strip().lower()
    x = x.replace("/", "_")
    x = x.replace("+", "_plus_")
    x = x.replace("-", "_")
    x = re.sub(r"[^a-z0-9]+", "_", x)
    x = re.sub(r"_+", "_", x).strip("_")
    return x


def shannon_entropy(values: Sequence[float]) -> float:
    arr = np.asarray([float(v) for v in values if not pd.isna(v) and float(v) > 0], dtype=float)
    if arr.size == 0:
        return 0.0
    p = arr / arr.sum()
    p = p[p > 0]
    if p.size == 0:
        return 0.0
    return float(-(p * np.log2(p)).sum())


def tau_specificity(values: Sequence[float]) -> float:
    """
    Tissue specificity tau. 0 = broad/even expression, 1 = restricted expression.
    """
    arr = np.asarray([float(v) for v in values if not pd.isna(v) and float(v) >= 0], dtype=float)
    if arr.size <= 1:
        return np.nan
    maxv = arr.max()
    if maxv <= 0:
        return 0.0
    return float(np.sum(1.0 - (arr / maxv)) / (arr.size - 1))


def find_col(df: pd.DataFrame, candidates: List[str], contains: Optional[List[str]] = None) -> Optional[str]:
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


def detect_expression_col(df: pd.DataFrame) -> Optional[str]:
    candidates = [
        "nTPM", "NX", "TPM", "pTPM", "Value", "value", "Expression", "expression",
        "Normalized expression", "normalized expression", "FPKM",
    ]
    col = find_col(df, candidates)
    if col:
        return col
    for c in df.columns:
        low = str(c).lower()
        if any(tok in low for tok in ["ntpm", "nx", "tpm", "expression", "value"]):
            return c
    return None


def open_text_from_path(path: Path):
    name = path.name.lower()
    if name.endswith(".gz"):
        return gzip.open(path, "rt", errors="ignore")
    return open(path, "rt", errors="ignore")


def read_table_flexible(path: Path, max_rows: int = 0) -> pd.DataFrame:
    """
    Read HPA TSV/CSV/TXT/GZ/ZIP robustly.
    ZIP files are expected to contain one TSV/TXT/CSV member.
    """
    if not path.exists():
        raise FileNotFoundError(path)
    nrows = max_rows if max_rows and max_rows > 0 else None
    name = path.name.lower()

    if name.endswith(".zip"):
        with zipfile.ZipFile(path, "r") as z:
            names = [n for n in z.namelist() if not n.endswith("/")]
            preferred = [n for n in names if n.lower().endswith((".tsv", ".txt", ".csv"))]
            if not preferred:
                preferred = names
            if not preferred:
                raise RuntimeError(f"ZIP is empty: {path}")
            member = sorted(preferred, key=lambda x: (0 if x.lower().endswith('.tsv') else 1, len(x)))[0]
            log(f"[ZIP MEMBER] {path.name} -> {member}")
            with z.open(member) as f:
                sep = "," if member.lower().endswith(".csv") else "\t"
                return pd.read_csv(f, sep=sep, dtype=str, low_memory=False, nrows=nrows)

    if name.endswith(".gz"):
        with gzip.open(path, "rt", errors="ignore") as f:
            sep = "," if name.endswith(".csv.gz") else "\t"
            return pd.read_csv(f, sep=sep, dtype=str, low_memory=False, nrows=nrows)

    sep = "," if name.endswith(".csv") else "\t"
    return pd.read_csv(path, sep=sep, dtype=str, low_memory=False, nrows=nrows)


# =============================================================================
# DOWNLOAD HELPERS
# =============================================================================


def looks_like_html_error(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return True
    try:
        with open(path, "rb") as f:
            head = f.read(1000).lower()
        markers = [b"<html", b"<!doctype html", b"404 not found", b"403 forbidden", b"access denied"]
        return any(m in head for m in markers)
    except Exception:
        return False


def is_valid_zip(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        with zipfile.ZipFile(path, "r") as z:
            if not z.namelist():
                return False
            bad = z.testzip()
            if bad is not None:
                log(f"[ZIP INVALID] {path}: first bad member={bad}")
                return False
        return True
    except Exception as exc:
        log(f"[ZIP INVALID] {path}: {exc}")
        return False


def download_with_python(url: str, outpath: Path, timeout: int = 180) -> bool:
    tmp = outpath.with_suffix(outpath.suffix + ".tmp")
    headers = {
        "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 Chrome/120 Safari/537.36",
        "Accept": "application/zip,application/octet-stream,text/tab-separated-values,*/*",
        "Accept-Language": "en-US,en;q=0.9",
    }
    req = urllib.request.Request(url, headers=headers)
    log(f"[DOWNLOAD PYTHON] {url}")
    log(f"[TO]                {outpath}")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as response:
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
                            log(f"[DOWNLOAD PYTHON] {downloaded/1024/1024:.1f} / {total/1024/1024:.1f} MB ({pct:.1f}%)")
                        else:
                            log(f"[DOWNLOAD PYTHON] {downloaded/1024/1024:.1f} MB")
                        last_print = time.time()
        tmp.rename(outpath)
        if looks_like_html_error(outpath):
            log("[DOWNLOAD PYTHON FAILED] HTML/error page downloaded.")
            outpath.unlink(missing_ok=True)
            return False
        if outpath.suffix.lower() == ".zip" and not is_valid_zip(outpath):
            log("[DOWNLOAD PYTHON FAILED] Invalid ZIP.")
            outpath.unlink(missing_ok=True)
            return False
        log(f"[DOWNLOAD PYTHON OK] {outpath} ({outpath.stat().st_size/1024/1024:.2f} MB)")
        return True
    except Exception as exc:
        log(f"[DOWNLOAD PYTHON ERROR] {exc}")
        tmp.unlink(missing_ok=True)
        return False


def download_with_command(url: str, outpath: Path, tool: str) -> bool:
    tmp = outpath.with_suffix(outpath.suffix + f".{tool}.tmp")
    if tool == "wget":
        cmd = ["wget", "--user-agent=Mozilla/5.0", "--tries=3", "--timeout=90", "-O", str(tmp), url]
    elif tool == "curl":
        cmd = ["curl", "-L", "--retry", "3", "--connect-timeout", "90", "-A", "Mozilla/5.0", "-o", str(tmp), url]
    else:
        return False
    log(f"[DOWNLOAD {tool.upper()}] {url}")
    log(f"[CMD] {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        if result.returncode != 0:
            log(f"[DOWNLOAD {tool.upper()} ERROR] returncode={result.returncode}")
            if result.stderr:
                log(result.stderr[-2000:])
            tmp.unlink(missing_ok=True)
            return False
        tmp.rename(outpath)
        if looks_like_html_error(outpath):
            log(f"[DOWNLOAD {tool.upper()} FAILED] HTML/error page downloaded.")
            outpath.unlink(missing_ok=True)
            return False
        if outpath.suffix.lower() == ".zip" and not is_valid_zip(outpath):
            log(f"[DOWNLOAD {tool.upper()} FAILED] Invalid ZIP.")
            outpath.unlink(missing_ok=True)
            return False
        log(f"[DOWNLOAD {tool.upper()} OK] {outpath} ({outpath.stat().st_size/1024/1024:.2f} MB)")
        return True
    except FileNotFoundError:
        log(f"[DOWNLOAD {tool.upper()} SKIP] {tool} not found.")
        return False
    except Exception as exc:
        log(f"[DOWNLOAD {tool.upper()} ERROR] {exc}")
        tmp.unlink(missing_ok=True)
        return False


def download_file(urls: Sequence[str], outpath: Path, force: bool = False, required: bool = False) -> bool:
    mkdir(outpath.parent)
    if outpath.exists() and outpath.stat().st_size > 0 and not force:
        if outpath.suffix.lower() == ".zip" and not is_valid_zip(outpath):
            log(f"[LOCAL INVALID] Removing invalid file: {outpath}")
            outpath.unlink(missing_ok=True)
        elif looks_like_html_error(outpath):
            log(f"[LOCAL INVALID] Removing HTML/error file: {outpath}")
            outpath.unlink(missing_ok=True)
        else:
            log(f"[DOWNLOAD SKIP] Existing file: {outpath}")
            return True
    if outpath.exists() and force:
        log(f"[REMOVE EXISTING] {outpath}")
        outpath.unlink(missing_ok=True)

    for url in urls:
        log("=" * 100)
        log(f"[HPA DOWNLOAD ATTEMPT] {url}")
        if download_with_python(url, outpath):
            return True
        if download_with_command(url, outpath, "wget"):
            return True
        if download_with_command(url, outpath, "curl"):
            return True

    msg = f"Could not download {outpath.name}."
    if required:
        raise RuntimeError(msg)
    log(f"[DOWNLOAD OPTIONAL FAILED] {msg} Continuing without this file.")
    return False


def candidate_urls_for_file(filename: str, aliases: Sequence[str]) -> List[str]:
    urls = []
    for fn in [filename] + list(aliases):
        if not fn:
            continue
        urls.append(f"{HPA_TSV_BASE_URL}/{fn}")
        # Older HPA pages sometimes used /download/<file> without /tsv/.
        urls.append(f"https://www.proteinatlas.org/download/{fn}")
    # De-duplicate while preserving order.
    out = []
    for u in urls:
        if u not in out:
            out.append(u)
    return out


def maybe_download_hpa_files(dbdir: Path, force: bool = False, include_large_risky_hpa: bool = False) -> Dict[str, Optional[Path]]:
    mkdir(dbdir)
    paths: Dict[str, Optional[Path]] = {}
    for key, spec in HPA_FILE_SPECS.items():
        filename = spec["filename"]
        outpath = dbdir / filename
        urls = candidate_urls_for_file(filename, spec.get("aliases", []))
        required = bool(spec.get("required", False))
        ok = download_file(urls, outpath, force=force, required=required)
        paths[key] = outpath if ok and outpath.exists() else None

    if include_large_risky_hpa:
        log("=" * 100)
        log("[WARNING] --include-large-risky-hpa was set. Downloading large/risky files, but this script does not parse them by default.")
        for key, filename in HPA_RISKY_OR_LARGE_FILES.items():
            outpath = dbdir / filename
            urls = candidate_urls_for_file(filename, [filename])
            ok = download_file(urls, outpath, force=force, required=False)
            paths[key] = outpath if ok and outpath.exists() else None

    return paths


def find_local_file(dbdir: Path, aliases: Sequence[str]) -> Optional[Path]:
    if not dbdir.exists():
        return None
    for name in aliases:
        p = dbdir / name
        if p.exists() and p.stat().st_size > 0 and not looks_like_html_error(p):
            return p
    # Fallback recursive search.
    files = []
    for name in aliases:
        files.extend(list(dbdir.rglob(name)))
    files = [p for p in files if p.exists() and p.stat().st_size > 0 and not looks_like_html_error(p)]
    if not files:
        return None
    files.sort(key=lambda p: p.stat().st_size, reverse=True)
    return files[0]


def detect_or_download_inputs(dbdir: Path, download: bool, force_download: bool, include_large_risky_hpa: bool) -> Dict[str, Optional[Path]]:
    paths: Dict[str, Optional[Path]] = {}
    for key, spec in HPA_FILE_SPECS.items():
        aliases = [spec["filename"]] + spec.get("aliases", [])
        paths[key] = find_local_file(dbdir, aliases)

    missing_required = [k for k, spec in HPA_FILE_SPECS.items() if spec.get("required", False) and paths.get(k) is None]
    missing_any = [k for k, p in paths.items() if p is None]

    if download or force_download or missing_required or missing_any:
        log("=" * 100)
        if missing_required:
            log(f"[HPA] Missing required files: {missing_required}. Attempting download.")
        elif missing_any:
            log(f"[HPA] Some optional files are missing: {missing_any}. Attempting download/reuse.")
        else:
            log("[HPA] Download requested. Attempting download/reuse.")
        downloaded = maybe_download_hpa_files(dbdir, force=force_download, include_large_risky_hpa=include_large_risky_hpa)
        for key, p in downloaded.items():
            if p is not None:
                paths[key] = p

    return paths


# =============================================================================
# HGNC
# =============================================================================


def load_hgnc(hgnc_path: Path, protein_coding_only: bool = True, limit_genes: int = 0) -> pd.DataFrame:
    if not hgnc_path.exists():
        raise FileNotFoundError(f"HGNC file not found: {hgnc_path}")
    log("=" * 100)
    log(f"[READ HGNC] {hgnc_path}")
    hgnc = pd.read_csv(hgnc_path, sep="\t", dtype=str, low_memory=False)
    hgnc.columns = [str(c).strip() for c in hgnc.columns]
    if "symbol" not in hgnc.columns:
        raise RuntimeError("HGNC file must contain symbol column.")
    before = len(hgnc)
    if protein_coding_only:
        if "locus_group" in hgnc.columns:
            hgnc = hgnc[hgnc["locus_group"].astype(str).str.lower().eq("protein-coding gene")].copy()
        elif "locus_type" in hgnc.columns:
            hgnc = hgnc[hgnc["locus_type"].astype(str).str.lower().str.contains("protein", na=False)].copy()
    hgnc["gene_symbol"] = hgnc["symbol"].map(normalize_symbol)
    for col in ["alias_symbol", "prev_symbol", "entrez_id", "ensembl_gene_id", "uniprot_ids", "name"]:
        if col not in hgnc.columns:
            hgnc[col] = ""
    hgnc["ensembl_gene_id_clean"] = hgnc["ensembl_gene_id"].map(clean_ensembl_gene_id)
    hgnc = hgnc.drop_duplicates(subset=["gene_symbol"], keep="first").copy()
    if limit_genes and limit_genes > 0:
        hgnc = hgnc.head(limit_genes).copy()
    log(f"[HGNC ROWS] {before} -> {len(hgnc)}")
    log(f"[PROTEIN CODING ONLY] {protein_coding_only}")
    log(f"[GENES WITH ENSEMBL] {(hgnc['ensembl_gene_id_clean'].astype(str).str.len() > 0).sum()}")
    return hgnc


def build_hgnc_maps(hgnc: pd.DataFrame) -> Tuple[Dict[str, str], Dict[str, str]]:
    symbol_map: Dict[str, str] = {}
    ensembl_map: Dict[str, str] = {}
    for _, row in hgnc.iterrows():
        gene = normalize_symbol(row.get("gene_symbol", ""))
        if not gene:
            continue
        symbol_map[gene] = gene
        ens = clean_ensembl_gene_id(row.get("ensembl_gene_id_clean", ""))
        if ens:
            ensembl_map[ens] = gene
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
    log(f"[ENSEMBL IDS]       {len(ensembl_map)}")
    return symbol_map, ensembl_map


def map_df_to_gene(df: pd.DataFrame, symbol_map: Dict[str, str], ensembl_map: Dict[str, str]) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    gene_col = find_col(df, ["Gene", "gene", "Ensembl", "Ensembl gene", "Gene ID", "Gene id"], contains=["gene"])
    gene_name_col = find_col(df, ["Gene name", "Gene name ", "gene_name", "Gene symbol", "Symbol", "gene_symbol"], contains=["gene", "name"])
    symbols = []
    ensgs = []
    for _, row in df.iterrows():
        gene = ""
        raw_ens = clean_ensembl_gene_id(row.get(gene_col, "")) if gene_col else ""
        raw_sym = normalize_symbol(row.get(gene_name_col, "")) if gene_name_col else ""
        if raw_ens and raw_ens in ensembl_map:
            gene = ensembl_map[raw_ens]
        elif raw_sym and raw_sym in symbol_map:
            gene = symbol_map[raw_sym]
        elif raw_ens and raw_ens.startswith("ENSG"):
            gene = ensembl_map.get(raw_ens, "")
        symbols.append(gene)
        ensgs.append(raw_ens)
    df["gene_symbol"] = symbols
    df["feature16_hpa_raw_ensembl_gene_id"] = ensgs
    df = df[df["gene_symbol"].astype(str).str.len() > 0].copy()
    return df


# =============================================================================
# FEATURE BUILDERS
# =============================================================================


def empty_features(hgnc: pd.DataFrame, prefix: str, has_col: str) -> pd.DataFrame:
    out = pd.DataFrame({"gene_symbol": sorted(hgnc["gene_symbol"].unique())})
    out[has_col] = 0
    return out


def parse_location_text(row: pd.Series, location_cols: List[str]) -> str:
    parts = []
    for c in location_cols:
        v = clean_text(row.get(c, ""))
        if v:
            parts.append(v)
    return " | ".join(parts)


def reliability_score(x: Any) -> float:
    s = clean_text(x).lower()
    if not s:
        return np.nan
    if "enhanced" in s:
        return 4.0
    if "supported" in s:
        return 3.0
    if "approved" in s:
        return 2.0
    if "uncertain" in s:
        return 1.0
    if "not" in s:
        return 0.0
    return np.nan


def build_subcellular_features(path: Optional[Path], hgnc: pd.DataFrame, symbol_map: Dict[str, str], ensembl_map: Dict[str, str], processed_dir: Path, max_rows: int = 0) -> pd.DataFrame:
    outpath = processed_dir / "feature16_hpa_subcellular_features.csv"
    if path is None or not path.exists():
        log("[SUBCELLULAR] Missing. Writing empty features.")
        out = empty_features(hgnc, "feature16_hpa_subcellular", "feature16_hpa_subcellular_has_data")
        out.to_csv(outpath, index=False)
        return out
    log("=" * 100)
    log(f"[READ HPA SUBCELLULAR] {path}")
    df = read_table_flexible(path, max_rows=max_rows)
    df = map_df_to_gene(df, symbol_map, ensembl_map)
    log(f"[SUBCELLULAR MAPPED] {df.shape}")

    location_cols = []
    for c in df.columns:
        low = str(c).lower()
        if any(tok in low for tok in ["location", "main", "additional", "extracellular"]):
            if "rna" not in low and "url" not in low:
                location_cols.append(c)
    location_cols = list(dict.fromkeys(location_cols))

    rel_col = find_col(df, ["Reliability", "reliability"])

    keyword_groups = {
        "plasma_membrane": ["plasma membrane", "cell membrane", "cell surface"],
        "membrane_any": ["membrane"],
        "secreted_or_extracellular": ["secreted", "extracellular", "extracellular matrix", "extracellular region", "secretory"],
        "cell_surface_accessible": ["plasma membrane", "cell surface", "extracellular", "secreted"],
        "nucleus": ["nucleus", "nucleoplasm", "nucleoli", "nuclear"],
        "cytosol": ["cytosol", "cytoplasm"],
        "mitochondria": ["mitochondria", "mitochondrion"],
        "endoplasmic_reticulum": ["endoplasmic reticulum"],
        "golgi": ["golgi"],
        "lysosome": ["lysosome", "lysosomes"],
        "peroxisome": ["peroxisome", "peroxisomes"],
        "vesicles": ["vesicle", "vesicles", "endosome", "endosomes"],
        "cytoskeleton": ["cytoskeleton", "actin", "microtubule", "intermediate filament"],
        "centrosome": ["centrosome", "centriolar", "centrosomal"],
        "cell_junctions": ["cell junction", "focal adhesion", "adherens", "tight junction"],
    }

    records = []
    for gene, g in df.groupby("gene_symbol"):
        texts = [parse_location_text(row, location_cols) for _, row in g.iterrows()]
        all_text = " | ".join(texts).lower()
        split_locs = []
        for t in texts:
            for part in re.split(r"[;,|]+", t):
                p = clean_text(part)
                if p:
                    split_locs.append(p.lower())
        unique_locs = sorted(set(split_locs))
        rec = {
            "gene_symbol": gene,
            "feature16_hpa_subcellular_has_data": 1,
            "feature16_hpa_subcellular_row_count": int(len(g)),
            "feature16_hpa_subcellular_location_count": int(len(unique_locs)),
        }
        for group, kws in keyword_groups.items():
            rec[f"feature16_hpa_subcellular_is_{group}"] = int(any(kw in all_text for kw in kws))
        if rel_col:
            rel_scores = [reliability_score(x) for x in g[rel_col].tolist()]
            rel_scores = [x for x in rel_scores if not pd.isna(x)]
            rec["feature16_hpa_subcellular_reliability_score_max"] = float(max(rel_scores)) if rel_scores else np.nan
            rec["feature16_hpa_subcellular_reliability_score_mean"] = float(np.mean(rel_scores)) if rel_scores else np.nan
        # Accessibility score for classical modalities: antibodies/biologics prefer surface/secreted.
        access = 0
        access += 3 * rec.get("feature16_hpa_subcellular_is_plasma_membrane", 0)
        access += 3 * rec.get("feature16_hpa_subcellular_is_secreted_or_extracellular", 0)
        access += 1 * rec.get("feature16_hpa_subcellular_is_membrane_any", 0)
        rec["feature16_hpa_subcellular_accessibility_score"] = int(access)
        records.append(rec)

    out = pd.DataFrame(records)
    genes = pd.DataFrame({"gene_symbol": sorted(hgnc["gene_symbol"].unique())})
    out = genes.merge(out, on="gene_symbol", how="left")
    flag_cols = [c for c in out.columns if c.startswith("feature16_hpa_subcellular_is_") or c.endswith("_count") or c.endswith("_score") or c.endswith("has_data")]
    for c in flag_cols:
        if c in out.columns and out[c].dtype != object:
            out[c] = out[c].fillna(0)
    out.to_csv(outpath, index=False)
    log(f"[SAVED SUBCELLULAR FEATURES] {outpath} {out.shape}")
    return out


def expression_group_features(g: pd.DataFrame, value_col: str, label_col: Optional[str], prefix: str) -> Dict[str, Any]:
    vals = pd.to_numeric(g[value_col], errors="coerce").dropna()
    rec: Dict[str, Any] = {}
    rec[f"{prefix}_row_count"] = int(len(g))
    if vals.empty:
        rec[f"{prefix}_max"] = np.nan
        rec[f"{prefix}_mean"] = np.nan
        rec[f"{prefix}_median"] = np.nan
        rec[f"{prefix}_sum"] = np.nan
        rec[f"{prefix}_detected_count_ge_1"] = 0
        rec[f"{prefix}_detected_count_ge_5"] = 0
        rec[f"{prefix}_detected_count_ge_10"] = 0
        rec[f"{prefix}_tau"] = np.nan
        rec[f"{prefix}_entropy"] = 0.0
        return rec
    arr = vals.to_numpy(dtype=float)
    rec[f"{prefix}_max"] = float(np.nanmax(arr))
    rec[f"{prefix}_mean"] = float(np.nanmean(arr))
    rec[f"{prefix}_median"] = float(np.nanmedian(arr))
    rec[f"{prefix}_sum"] = float(np.nansum(arr))
    rec[f"{prefix}_std"] = float(np.nanstd(arr))
    rec[f"{prefix}_detected_count_gt_0"] = int(np.sum(arr > 0))
    rec[f"{prefix}_detected_count_ge_1"] = int(np.sum(arr >= 1))
    rec[f"{prefix}_detected_count_ge_5"] = int(np.sum(arr >= 5))
    rec[f"{prefix}_detected_count_ge_10"] = int(np.sum(arr >= 10))
    rec[f"{prefix}_fraction_detected_ge_1"] = float(np.mean(arr >= 1))
    rec[f"{prefix}_tau"] = tau_specificity(arr)
    rec[f"{prefix}_entropy"] = shannon_entropy(arr)
    rec[f"{prefix}_broadly_expressed_ge_1_frac_ge_0_8"] = int(np.mean(arr >= 1) >= 0.8)
    rec[f"{prefix}_tissue_restricted_tau_ge_0_85"] = int((rec[f"{prefix}_tau"] if not pd.isna(rec[f"{prefix}_tau"]) else 0) >= 0.85)

    if label_col and label_col in g.columns:
        tmp = g[[label_col, value_col]].copy()
        tmp[value_col] = pd.to_numeric(tmp[value_col], errors="coerce")
        tmp = tmp.dropna(subset=[value_col])
        if not tmp.empty:
            top = tmp.sort_values(value_col, ascending=False).iloc[0]
            rec[f"{prefix}_top_label"] = clean_text(top[label_col])
            rec[f"{prefix}_unique_label_count"] = int(tmp[label_col].replace("", np.nan).dropna().nunique())
    return rec


def add_tissue_group_features(g: pd.DataFrame, value_col: str, label_col: Optional[str], prefix: str, rec: Dict[str, Any]) -> None:
    if not label_col or label_col not in g.columns:
        return
    tmp = g[[label_col, value_col]].copy()
    tmp[value_col] = pd.to_numeric(tmp[value_col], errors="coerce")
    tmp[label_col] = tmp[label_col].astype(str).str.lower()
    groups = {
        "brain": ["brain", "cerebell", "cortex", "hippocampus", "amygdala", "basal ganglia"],
        "liver": ["liver"],
        "kidney": ["kidney"],
        "heart": ["heart", "myocard"],
        "lung": ["lung", "bronch"],
        "intestine": ["intestine", "colon", "rectum", "duodenum", "ileum", "jejunum"],
        "stomach": ["stomach", "gastric"],
        "immune_blood": ["blood", "spleen", "lymph", "tonsil", "bone marrow", "thymus", "pbmc", "immune"],
        "reproductive": ["testis", "ovary", "prostate", "endometrium", "fallopian", "placenta", "seminal"],
        "endocrine": ["thyroid", "adrenal", "pituitary", "pancreas", "islet"],
        "muscle": ["muscle", "skeletal muscle", "smooth muscle"],
        "skin": ["skin", "epiderm"],
        "adipose": ["adipose", "fat"],
    }
    for group, kws in groups.items():
        mask = tmp[label_col].apply(lambda x: any(kw in x for kw in kws))
        vals = tmp.loc[mask, value_col].dropna().to_numpy(dtype=float)
        rec[f"{prefix}_{group}_max"] = float(np.max(vals)) if vals.size else np.nan
        rec[f"{prefix}_{group}_mean"] = float(np.mean(vals)) if vals.size else np.nan
        rec[f"{prefix}_{group}_detected_ge_1_count"] = int(np.sum(vals >= 1)) if vals.size else 0


def build_rna_features(path: Optional[Path], hgnc: pd.DataFrame, symbol_map: Dict[str, str], ensembl_map: Dict[str, str], processed_dir: Path, source_name: str, max_rows: int = 0) -> pd.DataFrame:
    prefix = f"feature16_hpa_{source_name}"
    outpath = processed_dir / f"feature16_hpa_{source_name}_features.csv"
    if path is None or not path.exists():
        log(f"[{source_name.upper()}] Missing. Writing empty features.")
        out = empty_features(hgnc, prefix, f"{prefix}_has_data")
        out.to_csv(outpath, index=False)
        return out
    log("=" * 100)
    log(f"[READ HPA {source_name.upper()}] {path}")
    df = read_table_flexible(path, max_rows=max_rows)
    df = map_df_to_gene(df, symbol_map, ensembl_map)
    value_col = detect_expression_col(df)
    label_col = find_col(df, ["Tissue", "tissue", "Sample", "sample", "Brain region", "Immune cell", "Cell type", "Cell line"], contains=["tissue"])
    if label_col is None:
        for candidate in ["Sample", "sample", "Immune cell", "Cell type", "Cell line", "Brain region"]:
            if candidate in df.columns:
                label_col = candidate
                break
    log(f"[{source_name.upper()} COL] value={value_col} label={label_col}")
    if value_col is None:
        log(f"[{source_name.upper()}] Could not detect expression value column. Empty features.")
        out = empty_features(hgnc, prefix, f"{prefix}_has_data")
        out.to_csv(outpath, index=False)
        return out

    records = []
    for gene, g in df.groupby("gene_symbol"):
        rec = {"gene_symbol": gene, f"{prefix}_has_data": 1}
        rec.update(expression_group_features(g, value_col=value_col, label_col=label_col, prefix=prefix))
        add_tissue_group_features(g, value_col=value_col, label_col=label_col, prefix=prefix, rec=rec)
        records.append(rec)
    out = pd.DataFrame(records)
    genes = pd.DataFrame({"gene_symbol": sorted(hgnc["gene_symbol"].unique())})
    out = genes.merge(out, on="gene_symbol", how="left")
    for c in out.columns:
        if c.startswith(prefix) and ("count" in c or "has_data" in c or "flag" in c or "broadly" in c or "restricted" in c):
            if out[c].dtype != object:
                out[c] = out[c].fillna(0)
    out.to_csv(outpath, index=False)
    log(f"[SAVED {source_name.upper()} FEATURES] {outpath} {out.shape}")
    return out


def ihc_level_score(x: Any) -> float:
    s = clean_text(x).lower()
    if not s:
        return np.nan
    if "not detected" in s or s in {"negative", "none"}:
        return 0.0
    if "low" in s or "weak" in s:
        return 1.0
    if "medium" in s or "moderate" in s:
        return 2.0
    if "high" in s or "strong" in s:
        return 3.0
    return np.nan


def build_ihc_features(path: Optional[Path], hgnc: pd.DataFrame, symbol_map: Dict[str, str], ensembl_map: Dict[str, str], processed_dir: Path, max_rows: int = 0) -> pd.DataFrame:
    prefix = "feature16_hpa_normal_ihc"
    outpath = processed_dir / "feature16_hpa_normal_ihc_features.csv"
    if path is None or not path.exists():
        log("[NORMAL IHC] Missing. Writing empty features.")
        out = empty_features(hgnc, prefix, f"{prefix}_has_data")
        out.to_csv(outpath, index=False)
        return out
    log("=" * 100)
    log(f"[READ HPA NORMAL IHC] {path}")
    df = read_table_flexible(path, max_rows=max_rows)
    df = map_df_to_gene(df, symbol_map, ensembl_map)
    tissue_col = find_col(df, ["Tissue", "tissue"])
    cell_col = find_col(df, ["Cell type", "cell type", "Cell type ", "cell_type"], contains=["cell"])
    level_col = find_col(df, ["Level", "level", "Expression", "expression", "Intensity", "Staining"], contains=["level"])
    if level_col is None:
        for c in df.columns:
            low = str(c).lower()
            if any(tok in low for tok in ["level", "expression", "intensity", "staining"]):
                level_col = c
                break
    rel_col = find_col(df, ["Reliability", "reliability"])
    log(f"[NORMAL IHC COL] tissue={tissue_col} cell={cell_col} level={level_col} reliability={rel_col}")
    if level_col is None:
        out = empty_features(hgnc, prefix, f"{prefix}_has_data")
        out.to_csv(outpath, index=False)
        return out
    records = []
    df["_ihc_score"] = df[level_col].map(ihc_level_score)
    for gene, g in df.groupby("gene_symbol"):
        vals = pd.to_numeric(g["_ihc_score"], errors="coerce").dropna().to_numpy(dtype=float)
        rec = {"gene_symbol": gene, f"{prefix}_has_data": 1, f"{prefix}_row_count": int(len(g))}
        rec[f"{prefix}_max_score"] = float(np.max(vals)) if vals.size else np.nan
        rec[f"{prefix}_mean_score"] = float(np.mean(vals)) if vals.size else np.nan
        rec[f"{prefix}_median_score"] = float(np.median(vals)) if vals.size else np.nan
        rec[f"{prefix}_high_count"] = int(np.sum(vals >= 3)) if vals.size else 0
        rec[f"{prefix}_medium_or_high_count"] = int(np.sum(vals >= 2)) if vals.size else 0
        rec[f"{prefix}_detected_count"] = int(np.sum(vals > 0)) if vals.size else 0
        rec[f"{prefix}_not_detected_count"] = int(np.sum(vals == 0)) if vals.size else 0
        rec[f"{prefix}_detected_fraction"] = float(np.mean(vals > 0)) if vals.size else np.nan
        if tissue_col:
            rec[f"{prefix}_unique_tissue_count"] = int(g[tissue_col].replace("", np.nan).dropna().nunique())
        if cell_col:
            rec[f"{prefix}_unique_cell_type_count"] = int(g[cell_col].replace("", np.nan).dropna().nunique())
        if rel_col:
            rvals = [reliability_score(x) for x in g[rel_col].tolist()]
            rvals = [x for x in rvals if not pd.isna(x)]
            rec[f"{prefix}_reliability_score_max"] = float(max(rvals)) if rvals else np.nan
        records.append(rec)
    out = pd.DataFrame(records)
    genes = pd.DataFrame({"gene_symbol": sorted(hgnc["gene_symbol"].unique())})
    out = genes.merge(out, on="gene_symbol", how="left")
    for c in out.columns:
        if c.startswith(prefix) and ("count" in c or "has_data" in c):
            if out[c].dtype != object:
                out[c] = out[c].fillna(0)
    out.to_csv(outpath, index=False)
    log(f"[SAVED NORMAL IHC FEATURES] {outpath} {out.shape}")
    return out


def build_proteinatlas_summary_features(path: Optional[Path], hgnc: pd.DataFrame, symbol_map: Dict[str, str], ensembl_map: Dict[str, str], processed_dir: Path, max_rows: int = 0) -> pd.DataFrame:
    prefix = "feature16_hpa_proteinatlas"
    outpath = processed_dir / "feature16_hpa_proteinatlas_features.csv"
    if path is None or not path.exists():
        log("[PROTEINATLAS] Missing. Writing empty features.")
        out = empty_features(hgnc, prefix, f"{prefix}_has_data")
        out.to_csv(outpath, index=False)
        return out
    log("=" * 100)
    log(f"[READ HPA PROTEINATLAS SUMMARY] {path}")
    df = read_table_flexible(path, max_rows=max_rows)
    df = map_df_to_gene(df, symbol_map, ensembl_map)
    records = []
    safe_keyword_cols = []
    for c in df.columns:
        low = str(c).lower()
        if any(tok in low for tok in ["protein class", "biological process", "molecular function", "subcellular", "rna tissue specificity", "rna tissue distribution", "secretome", "membrane"]):
            safe_keyword_cols.append(c)
    for gene, g in df.groupby("gene_symbol"):
        rec = {"gene_symbol": gene, f"{prefix}_has_data": 1, f"{prefix}_row_count": int(len(g))}
        text = " | ".join(clean_text(v) for c in safe_keyword_cols for v in g[c].tolist()).lower() if safe_keyword_cols else ""
        kws = {
            "protein_class_transporter": ["transporter", "channel", "solute carrier", "abc transporter"],
            "protein_class_receptor": ["receptor", "gpcr", "g-protein coupled receptor"],
            "protein_class_enzyme": ["enzyme", "kinase", "protease", "phosphatase", "transferase", "hydrolase"],
            "protein_class_kinase": ["kinase"],
            "protein_class_ion_channel": ["ion channel", "channel"],
            "protein_class_secreted": ["secreted", "secretome"],
            "protein_class_membrane": ["membrane"],
            "protein_class_tf": ["transcription factor"],
            "rna_tissue_enriched": ["tissue enriched"],
            "rna_group_enriched": ["group enriched"],
            "rna_tissue_enhanced": ["tissue enhanced"],
            "rna_low_specificity": ["low tissue specificity", "detected in all", "detected in many"],
        }
        for k, terms in kws.items():
            rec[f"{prefix}_{k}"] = int(any(t in text for t in terms))
        records.append(rec)
    out = pd.DataFrame(records)
    genes = pd.DataFrame({"gene_symbol": sorted(hgnc["gene_symbol"].unique())})
    out = genes.merge(out, on="gene_symbol", how="left")
    for c in out.columns:
        if c.startswith(prefix) and out[c].dtype != object:
            out[c] = out[c].fillna(0)
    out.to_csv(outpath, index=False)
    log(f"[SAVED PROTEINATLAS FEATURES] {outpath} {out.shape}")
    return out


# =============================================================================
# MERGE / SUMMARY
# =============================================================================


def merge_feature_tables(hgnc: pd.DataFrame, tables: List[pd.DataFrame], processed_dir: Path) -> pd.DataFrame:
    genes = pd.DataFrame({"gene_symbol": sorted(hgnc["gene_symbol"].unique())})
    features = genes.copy()
    for t in tables:
        if t is None or t.empty:
            continue
        # Avoid duplicate columns except gene_symbol.
        dup = [c for c in t.columns if c in features.columns and c != "gene_symbol"]
        if dup:
            t = t.drop(columns=dup)
        features = features.merge(t, on="gene_symbol", how="left")

    # Combined druggability-oriented safe scores.
    def col(name: str) -> pd.Series:
        return pd.to_numeric(features[name], errors="coerce").fillna(0) if name in features.columns else pd.Series(0, index=features.index)

    accessibility = col("feature16_hpa_subcellular_accessibility_score")
    membrane = col("feature16_hpa_subcellular_is_plasma_membrane")
    secreted = col("feature16_hpa_subcellular_is_secreted_or_extracellular")
    receptor = col("feature16_hpa_proteinatlas_protein_class_receptor")
    transporter = col("feature16_hpa_proteinatlas_protein_class_transporter")
    enzyme = col("feature16_hpa_proteinatlas_protein_class_enzyme")
    consensus_max = col("feature16_hpa_rna_consensus_max")
    consensus_tau = col("feature16_hpa_rna_consensus_tau")
    ihc_detected = col("feature16_hpa_normal_ihc_detected_fraction")

    features["feature16_hpa_surface_secreted_accessibility_score"] = accessibility
    features["feature16_hpa_classical_druggability_proxy_score"] = (
        2.0 * membrane + 2.0 * secreted + 1.5 * receptor + 1.2 * transporter + 1.0 * enzyme
        + np.clip(np.log1p(consensus_max) / 5.0, 0, 1)
        + np.clip(ihc_detected, 0, 1)
    )
    features["feature16_hpa_expression_specificity_proxy_score"] = consensus_tau

    outpath = processed_dir / "feature16_hpa_gene_features.csv"
    features.to_csv(outpath, index=False)
    log("=" * 100)
    log(f"[SAVED FEATURE TABLE] {outpath} {features.shape}")
    return features


def merge_with_hgnc(hgnc: pd.DataFrame, features: pd.DataFrame, processed_dir: Path) -> pd.DataFrame:
    merged = hgnc.merge(features, on="gene_symbol", how="left")
    feature_cols = [c for c in features.columns if c != "gene_symbol"]
    merged["feature16_hpa_has_any_feature"] = merged[feature_cols].notna().any(axis=1).astype(int)
    merged["feature16_hpa_n_nonmissing_features"] = merged[feature_cols].notna().sum(axis=1)
    outpath = processed_dir / "feature16_hpa_hgnc_merged.csv"
    merged.to_csv(outpath, index=False)
    log(f"[SAVED HGNC MERGED] {outpath} {merged.shape}")
    return merged


def write_summary(outdir: Path, hgnc: pd.DataFrame, paths: Dict[str, Optional[Path]], features: pd.DataFrame, merged: pd.DataFrame, args: argparse.Namespace) -> None:
    lines = []
    lines.append("Feature 16: HPA subcellular location + expression features")
    lines.append("=" * 80)
    lines.append(f"Created at: {now_iso()}")
    lines.append("")
    lines.append("Inputs:")
    lines.append(f"HGNC: {args.hgnc}")
    lines.append(f"HPA directory: {args.dbdir}")
    for k, p in paths.items():
        lines.append(f"{k}: {p if p else 'missing'}")
    lines.append("")
    lines.append("Rows:")
    lines.append(f"HGNC rows: {hgnc.shape[0]}")
    lines.append(f"Gene feature rows: {features.shape[0]}")
    lines.append(f"Gene feature columns: {features.shape[1]}")
    lines.append(f"Merged rows: {merged.shape[0]}")
    lines.append(f"Merged columns: {merged.shape[1]}")
    lines.append("")
    lines.append("Coverage:")
    for c in features.columns:
        if c.endswith("_has_data"):
            x = pd.to_numeric(features[c], errors="coerce").fillna(0)
            lines.append(f"{c}: genes={int(x.sum())}, fraction={x.mean():.4f}")
    lines.append("")
    lines.append("Key feature summaries:")
    for c in [
        "feature16_hpa_subcellular_accessibility_score",
        "feature16_hpa_classical_druggability_proxy_score",
        "feature16_hpa_rna_consensus_max",
        "feature16_hpa_rna_consensus_tau",
        "feature16_hpa_normal_ihc_detected_fraction",
    ]:
        if c in features.columns:
            x = pd.to_numeric(features[c], errors="coerce")
            lines.append(f"{c}: median={x.median(skipna=True):.4f}, max={x.max(skipna=True):.4f}, nonmissing={int(x.notna().sum())}")
    lines.append("")
    lines.append("Leakage policy:")
    lines.append("Included: HPA localization, normal tissue expression, normal protein expression, immune/single-cell expression, safe protein class keywords.")
    lines.append("Excluded by default: cancer prognosis/survival, TCGA survival outcomes, drug/chemical interactions, known drug-target labels.")
    path = outdir / "feature16_hpa_summary.txt"
    path.write_text("\n".join(lines) + "\n")
    log(f"[SUMMARY] {path}")


# =============================================================================
# MAIN
# =============================================================================


def main() -> None:
    parser = argparse.ArgumentParser(description="Build Feature 16 HPA localization/expression features.")
    parser.add_argument("--hgnc", default=str(DEFAULT_HGNC_PATH), help="Path to HGNC complete set.")
    parser.add_argument("--dbdir", default=str(DEFAULT_DBDIR), help="HPA database directory.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory.")
    parser.add_argument("--download", action="store_true", help="Download/reuse HPA files. Missing files are downloaded automatically even without this flag.")
    parser.add_argument("--force-download", action="store_true", help="Force re-download existing HPA files.")
    parser.add_argument("--max-rows", type=int, default=0, help="Debug only: read first N rows per HPA file.")
    parser.add_argument("--max-rna-rows", type=int, default=0, help="Alias for --max-rows for compatibility with older script.")
    parser.add_argument("--limit-genes", type=int, default=0, help="Debug only: first N HGNC genes.")
    parser.add_argument("--all-hgnc-genes", action="store_true", help="Use all HGNC genes, not only protein-coding.")
    parser.add_argument("--include-large-risky-hpa", action="store_true", help="Download large/risky HPA cancer/prognosis files. Not parsed by default; not recommended for leakage-safe druggability model.")
    parser.add_argument("--subcellular-file", default="", help="Manual subcellular_location.tsv.zip path.")
    parser.add_argument("--rna-consensus-file", default="", help="Manual rna_tissue_consensus.tsv.zip path.")
    parser.add_argument("--rna-file", default="", help="Alias for --rna-consensus-file.")
    parser.add_argument("--rna-hpa-file", default="", help="Manual rna_tissue_hpa.tsv.zip path.")
    parser.add_argument("--rna-gtex-file", default="", help="Manual rna_tissue_gtex.tsv.zip path.")
    parser.add_argument("--normal-ihc-file", default="", help="Manual normal_ihc_data.tsv.zip path.")
    parser.add_argument("--immune-cell-file", default="", help="Manual rna_immune_cell.tsv.zip path.")
    parser.add_argument("--single-cell-file", default="", help="Manual rna_single_cell_type.tsv.zip path.")
    parser.add_argument("--proteinatlas-file", default="", help="Manual proteinatlas.tsv.zip path.")
    args = parser.parse_args()

    dbdir = mkdir(Path(args.dbdir))
    outdir = mkdir(Path(args.outdir))
    processed_dir = mkdir(outdir / "processed")
    max_rows = args.max_rows or args.max_rna_rows

    log("=" * 100)
    log("FEATURE 16: HPA SUBCELLULAR LOCATION + EXPRESSION FEATURES")
    log("=" * 100)
    log(f"[HGNC]           {args.hgnc}")
    log(f"[DBDIR]          {dbdir.resolve()}")
    log(f"[OUTDIR]         {outdir.resolve()}")
    log(f"[DOWNLOAD]       {args.download}")
    log(f"[FORCE DOWNLOAD] {args.force_download}")
    log(f"[MAX ROWS]       {max_rows if max_rows else 'none'}")
    log(f"[LIMIT GENES]    {args.limit_genes if args.limit_genes else 'none'}")
    log("=" * 100)

    paths = detect_or_download_inputs(
        dbdir=dbdir,
        download=args.download,
        force_download=args.force_download,
        include_large_risky_hpa=args.include_large_risky_hpa,
    )

    # Manual overrides.
    if args.subcellular_file:
        paths["subcellular"] = Path(args.subcellular_file)
    if args.rna_consensus_file or args.rna_file:
        paths["rna_consensus"] = Path(args.rna_consensus_file or args.rna_file)
    if args.rna_hpa_file:
        paths["rna_hpa"] = Path(args.rna_hpa_file)
    if args.rna_gtex_file:
        paths["rna_gtex"] = Path(args.rna_gtex_file)
    if args.normal_ihc_file:
        paths["normal_ihc"] = Path(args.normal_ihc_file)
    if args.immune_cell_file:
        paths["rna_immune_cell"] = Path(args.immune_cell_file)
    if args.single_cell_file:
        paths["rna_single_cell"] = Path(args.single_cell_file)
    if args.proteinatlas_file:
        paths["proteinatlas"] = Path(args.proteinatlas_file)

    hgnc = load_hgnc(Path(args.hgnc), protein_coding_only=not args.all_hgnc_genes, limit_genes=args.limit_genes)
    symbol_map, ensembl_map = build_hgnc_maps(hgnc)

    feature_tables = []
    feature_tables.append(build_subcellular_features(paths.get("subcellular"), hgnc, symbol_map, ensembl_map, processed_dir, max_rows=max_rows))
    feature_tables.append(build_rna_features(paths.get("rna_consensus"), hgnc, symbol_map, ensembl_map, processed_dir, "rna_consensus", max_rows=max_rows))
    feature_tables.append(build_rna_features(paths.get("rna_hpa"), hgnc, symbol_map, ensembl_map, processed_dir, "rna_hpa", max_rows=max_rows))
    feature_tables.append(build_rna_features(paths.get("rna_gtex"), hgnc, symbol_map, ensembl_map, processed_dir, "rna_gtex", max_rows=max_rows))
    feature_tables.append(build_ihc_features(paths.get("normal_ihc"), hgnc, symbol_map, ensembl_map, processed_dir, max_rows=max_rows))
    feature_tables.append(build_rna_features(paths.get("rna_immune_cell"), hgnc, symbol_map, ensembl_map, processed_dir, "rna_immune_cell", max_rows=max_rows))
    feature_tables.append(build_rna_features(paths.get("rna_single_cell"), hgnc, symbol_map, ensembl_map, processed_dir, "rna_single_cell", max_rows=max_rows))
    feature_tables.append(build_proteinatlas_summary_features(paths.get("proteinatlas"), hgnc, symbol_map, ensembl_map, processed_dir, max_rows=max_rows))

    features = merge_feature_tables(hgnc, feature_tables, processed_dir)
    merged = merge_with_hgnc(hgnc, features, processed_dir)

    metadata = {
        "created_at": now_iso(),
        "hgnc": str(Path(args.hgnc)),
        "dbdir": str(dbdir.resolve()),
        "outdir": str(outdir.resolve()),
        "paths": {k: str(v) if v else None for k, v in paths.items()},
        "max_rows": max_rows,
        "limit_genes": args.limit_genes,
        "protein_coding_only": not args.all_hgnc_genes,
        "outputs": {
            "gene_features": str(processed_dir / "feature16_hpa_gene_features.csv"),
            "hgnc_merged": str(processed_dir / "feature16_hpa_hgnc_merged.csv"),
            "summary": str(outdir / "feature16_hpa_summary.txt"),
        },
        "leakage_policy": {
            "included": [
                "HPA subcellular localization",
                "HPA/GTEx tissue RNA expression summaries",
                "normal tissue IHC protein expression",
                "immune-cell and single-cell RNA expression when available",
                "safe HPA protein class keywords when available",
            ],
            "excluded_by_default": [
                "cancer prognosis/survival",
                "TCGA survival outcomes",
                "drug/chemical interactions",
                "known drug-target labels",
                "ChEMBL",
                "DrugBank",
                "DGIdb",
                "Open Targets",
                "Pharos",
            ],
        },
    }
    with open(outdir / "feature16_hpa_run_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    write_summary(outdir, hgnc, paths, features, merged, args)

    log("=" * 100)
    log("[DONE]")
    log(f"[GENE FEATURES] {processed_dir / 'feature16_hpa_gene_features.csv'}")
    log(f"[HGNC MERGED]   {processed_dir / 'feature16_hpa_hgnc_merged.csv'}")
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
