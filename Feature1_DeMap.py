# -*- coding: utf-8 -*-

"""
feature1_dmapp.py

Feature 1: DepMap functional-genomics feature builder for gene druggability modelling.

What this script does
---------------------
1. Reads your existing HGNC gene universe:
      databases/HGNC/hgnc_complete_set.txt

2. Downloads selected safe DepMap Public 26Q1 files using the working DepMap
   file-index API:
      https://depmap.org/portal/api/download/files

3. Stores all raw DepMap files here:
      feature1.dmapp.database/downloads/

4. Converts DepMap matrices into one gene-level feature table:
      feature1.dmapp.database/processed/feature1_depmap_gene_features.csv

5. Merges the DepMap features onto the HGNC protein-coding gene table:
      feature1.dmapp.database/processed/feature1_depmap_hgnc_merged.csv

Important leakage rule
----------------------
This script deliberately avoids compound/drug-target files such as PortalCompounds.csv.
It only uses functional-genomics data: CRISPR gene effect, CRISPR dependency,
expression, mutation matrices, copy number, and optional fusion/subtype metadata.

Run
---
Basic safe run:
    python feature1_dmapp.py

Download only:
    python feature1_dmapp.py --download-only

Include optional omics files:
    python feature1_dmapp.py --include-optional

Use all HGNC genes, not only protein-coding:
    python feature1_dmapp.py --all-hgnc-genes

Custom HGNC path:
    python feature1_dmapp.py --hgnc databases/HGNC/hgnc_complete_set.txt

Custom release:
    python feature1_dmapp.py --release "DepMap Public 26Q1"
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import requests


# =============================================================================
# DEFAULTS
# =============================================================================

DEFAULT_RELEASE = "DepMap Public 26Q1"
DEPMAP_INDEX_URL = "https://depmap.org/portal/api/download/files"

DEFAULT_HGNC_PATH = Path("databases") / "HGNC" / "hgnc_complete_set.txt"
DEFAULT_OUTDIR = Path("feature1.dmapp.database")

# Safe/core files. These are functional genomics and expression features.
REQUIRED_DEPMAP_FILES = [
    "Gene.csv",
    "Model.csv",
    "CRISPRGeneEffect.csv",
    "CRISPRGeneDependency.csv",
    "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
]

# Optional useful cancer-omics features.
OPTIONAL_DEPMAP_FILES = [
    "OmicsSomaticMutationsMatrixHotspot.csv",
    "OmicsSomaticMutationsMatrixDamaging.csv",
    "OmicsCNGeneWGS.csv",
    "OmicsFusionFiltered.csv",
    "SubtypeMatrix.csv",
    "SubtypeTree.csv",
]

# Do not download or use as modelling features because they may encode known drug-target relationships.
LEAKAGE_RISK_FILES = [
    "PortalCompounds.csv",
]

MATRIX_FILES = {
    "CRISPRGeneEffect.csv": "depmap_crispr_gene_effect",
    "CRISPRGeneDependency.csv": "depmap_crispr_gene_dependency",
    "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv": "depmap_expression_tpm_logp1",
    "OmicsSomaticMutationsMatrixHotspot.csv": "depmap_hotspot_mutation",
    "OmicsSomaticMutationsMatrixDamaging.csv": "depmap_damaging_mutation",
    "OmicsCNGeneWGS.csv": "depmap_copy_number_wgs",
}

FUSION_FILE = "OmicsFusionFiltered.csv"


# =============================================================================
# GENERAL HELPERS
# =============================================================================

def now_iso() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str) -> None:
    print(msg, flush=True)


def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def md5_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk_size)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def safe_float_series(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def read_csv_auto(path: Path, nrows: Optional[int] = None) -> pd.DataFrame:
    """
    Read CSV/TSV robustly. HGNC txt is tab-separated; DepMap is CSV.
    """
    if path.suffix.lower() in [".txt", ".tsv"]:
        return pd.read_csv(path, sep="\t", dtype=str, low_memory=False, nrows=nrows)
    return pd.read_csv(path, dtype=str, low_memory=False, nrows=nrows)


def read_numeric_matrix(path: Path) -> pd.DataFrame:
    """
    Read a DepMap matrix where rows are ModelID/DepMap_ID and columns are genes.
    Keeps the first ID column as object and all gene columns as numeric where possible.
    """
    log(f"[READ MATRIX] {path}")
    df = pd.read_csv(path, low_memory=False)
    if df.shape[1] < 2:
        raise RuntimeError(f"Matrix file has too few columns: {path}")
    return df


def normalize_gene_symbol(x: object) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


def parse_depmap_gene_column(col: str) -> Tuple[str, Optional[str], str]:
    """
    DepMap gene matrix columns commonly look like:
        TP53 (7157)
        A1BG (1)
        ENSG... style can also appear in some files.

    Returns:
        symbol_upper, entrez_id_or_none, raw_column
    """
    raw = str(col).strip()

    # Common pattern: SYMBOL (ENTREZ)
    m = re.match(r"^(.+?)\s*\(([^()]+)\)\s*$", raw)
    if m:
        symbol = m.group(1).strip()
        entrez = m.group(2).strip()
        return normalize_gene_symbol(symbol), entrez, raw

    # Fallback: first token as symbol
    # This keeps "TP53" as TP53 and avoids breaking if no Entrez is present.
    symbol = raw.strip()
    return normalize_gene_symbol(symbol), None, raw


def get_first_existing_col(df: pd.DataFrame, candidates: Iterable[str]) -> Optional[str]:
    lower_map = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower_map:
            return lower_map[c.lower()]
    return None


# =============================================================================
# DEPMAP DOWNLOAD
# =============================================================================

def make_session() -> requests.Session:
    s = requests.Session()
    s.headers.update({
        "User-Agent": "Mozilla/5.0 DepMapFeatureDownloader/1.0",
        "Accept": "text/csv,application/json,text/plain,*/*",
    })
    return s


def download_text_or_file(
    session: requests.Session,
    url: str,
    outpath: Path,
    retries: int = 5,
    timeout: int = 300,
) -> bool:
    """
    Download a URL to a file with retries.
    """
    if outpath.exists() and outpath.stat().st_size > 0:
        log(f"[SKIP] Already exists: {outpath} ({outpath.stat().st_size / 1024 / 1024:.2f} MB)")
        return True

    tmp = outpath.with_suffix(outpath.suffix + ".tmp")

    for attempt in range(1, retries + 1):
        try:
            log("=" * 100)
            log(f"[DOWNLOAD] {outpath.name}")
            log(f"[URL]      {url}")
            log(f"[SAVE TO]  {outpath}")

            with session.get(url, stream=True, timeout=timeout, allow_redirects=True) as r:
                if r.status_code != 200:
                    text = ""
                    try:
                        text = r.text[:500]
                    except Exception:
                        pass
                    raise RuntimeError(f"HTTP {r.status_code}: {text}")

                total = int(r.headers.get("content-length", 0))
                downloaded = 0

                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if not chunk:
                            continue
                        f.write(chunk)
                        downloaded += len(chunk)

                        if total > 0:
                            pct = downloaded / total * 100
                            print(
                                f"\r       {downloaded / 1024 / 1024:.2f} MB / "
                                f"{total / 1024 / 1024:.2f} MB ({pct:.1f}%)",
                                end="",
                                flush=True,
                            )
                        else:
                            print(
                                f"\r       {downloaded / 1024 / 1024:.2f} MB",
                                end="",
                                flush=True,
                            )
                print()

            tmp.replace(outpath)

            if outpath.stat().st_size == 0:
                raise RuntimeError("Downloaded file is empty")

            log(f"[DONE] {outpath.name} ({outpath.stat().st_size / 1024 / 1024:.2f} MB)")
            return True

        except Exception as e:
            log(f"[FAILED] Attempt {attempt}/{retries}: {outpath.name}")
            log(f"[ERROR]  {e}")
            if tmp.exists():
                tmp.unlink()
            if attempt < retries:
                wait = 5 * attempt
                log(f"[RETRY] Waiting {wait} seconds...")
                time.sleep(wait)

    return False


def load_or_download_depmap_index(session: requests.Session, downloads_dir: Path, force: bool = False) -> pd.DataFrame:
    index_path = downloads_dir / "depmap_files.csv"

    if force or not index_path.exists() or index_path.stat().st_size == 0:
        ok = download_text_or_file(session, DEPMAP_INDEX_URL, index_path, retries=5, timeout=120)
        if not ok:
            raise RuntimeError("Could not download DepMap file index.")

    df = pd.read_csv(index_path, dtype=str)
    required = {"release", "release_date", "filename", "url", "md5_hash"}
    missing = required.difference(df.columns)
    if missing:
        raise RuntimeError(f"DepMap index missing columns: {sorted(missing)}")

    return df


def select_depmap_files(index_df: pd.DataFrame, release: str, filenames: List[str]) -> pd.DataFrame:
    sub = index_df[index_df["release"].astype(str).eq(release)].copy()

    if sub.empty:
        available = sorted(index_df["release"].dropna().unique().tolist())
        raise RuntimeError(
            f"No files found for release '{release}'. Available releases include:\n"
            + "\n".join(f"  - {x}" for x in available[:50])
        )

    selected = sub[sub["filename"].isin(filenames)].copy()

    found = set(selected["filename"].tolist())
    missing = [f for f in filenames if f not in found]

    log("=" * 100)
    log("[DEPMAP FILE SELECTION]")
    log(f"[RELEASE] {release}")
    for f in filenames:
        if f in found:
            log(f"  [FOUND]   {f}")
        else:
            log(f"  [MISSING] {f}")

    if missing:
        log("[WARNING] Some requested files were not found in the DepMap index.")
        log("[WARNING] The script will continue with available files.")

    return selected


def download_depmap_files(
    release: str,
    downloads_dir: Path,
    include_optional: bool,
    force_index: bool,
    force_files: bool,
) -> Path:
    session = make_session()
    mkdir(downloads_dir)

    index_df = load_or_download_depmap_index(session, downloads_dir, force=force_index)

    wanted = list(REQUIRED_DEPMAP_FILES)
    if include_optional:
        wanted.extend(OPTIONAL_DEPMAP_FILES)

    wanted = [f for f in wanted if f not in LEAKAGE_RISK_FILES]

    selected = select_depmap_files(index_df, release, wanted)

    metadata = {
        "created_at": now_iso(),
        "release": release,
        "index_url": DEPMAP_INDEX_URL,
        "downloads_dir": str(downloads_dir.resolve()),
        "requested_files": wanted,
        "leakage_risk_files_excluded": LEAKAGE_RISK_FILES,
        "files": [],
    }

    failed = []

    for _, row in selected.iterrows():
        fname = row["filename"]
        url = row["url"]
        outpath = downloads_dir / fname

        if force_files and outpath.exists():
            outpath.unlink()

        ok = download_text_or_file(session, url, outpath, retries=5, timeout=600)
        if not ok:
            failed.append(fname)
            continue

        observed_md5 = md5_file(outpath)
        expected_md5 = str(row.get("md5_hash", "") or "").strip()

        md5_ok = None
        if expected_md5 and expected_md5.lower() != "nan":
            md5_ok = observed_md5.lower() == expected_md5.lower()
            if md5_ok:
                log(f"[MD5 OK] {fname}")
            else:
                log(f"[MD5 WARNING] {fname}")
                log(f"  expected: {expected_md5}")
                log(f"  observed: {observed_md5}")

        metadata["files"].append({
            "filename": fname,
            "release": row.get("release"),
            "release_date": row.get("release_date"),
            "source_url": url,
            "path": str(outpath),
            "size_bytes": outpath.stat().st_size,
            "md5_expected": expected_md5,
            "md5_observed": observed_md5,
            "md5_ok": md5_ok,
        })

    metadata["failed"] = failed

    with open(downloads_dir / "feature1_depmap_download_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    if failed:
        log("=" * 100)
        log("[FAILED DOWNLOADS]")
        for f in failed:
            log(f"  - {f}")

    return downloads_dir


# =============================================================================
# HGNC INPUT
# =============================================================================

def load_hgnc(hgnc_path: Path, protein_coding_only: bool = True) -> pd.DataFrame:
    if not hgnc_path.exists():
        raise FileNotFoundError(
            f"HGNC file not found: {hgnc_path}\n"
            "Expected default path: databases/HGNC/hgnc_complete_set.txt"
        )

    log("=" * 100)
    log(f"[READ HGNC] {hgnc_path}")

    hgnc = pd.read_csv(hgnc_path, sep="\t", dtype=str, low_memory=False)
    hgnc.columns = [c.strip() for c in hgnc.columns]

    if "symbol" not in hgnc.columns:
        raise RuntimeError("HGNC file must contain a 'symbol' column.")

    before = len(hgnc)

    if protein_coding_only:
        # Your HGNC file shows locus_group = "protein-coding gene".
        if "locus_group" in hgnc.columns:
            hgnc = hgnc[hgnc["locus_group"].astype(str).str.lower().eq("protein-coding gene")].copy()
        elif "locus_type" in hgnc.columns:
            hgnc = hgnc[hgnc["locus_type"].astype(str).str.lower().str.contains("protein", na=False)].copy()
        else:
            log("[WARNING] Could not find locus_group/locus_type. Keeping all HGNC genes.")

    hgnc["gene_symbol"] = hgnc["symbol"].map(normalize_gene_symbol)

    if "entrez_id" in hgnc.columns:
        hgnc["entrez_id_str"] = hgnc["entrez_id"].astype(str).str.strip().replace({"nan": "", "None": ""})
    else:
        hgnc["entrez_id_str"] = ""

    if "ensembl_gene_id" in hgnc.columns:
        hgnc["ensembl_gene_id_str"] = hgnc["ensembl_gene_id"].astype(str).str.strip().replace({"nan": "", "None": ""})
    else:
        hgnc["ensembl_gene_id_str"] = ""

    # Keep one row per approved symbol.
    hgnc = hgnc.drop_duplicates(subset=["gene_symbol"], keep="first").copy()

    log(f"[HGNC ROWS] {before} -> {len(hgnc)}")
    log(f"[PROTEIN CODING ONLY] {protein_coding_only}")

    return hgnc


# =============================================================================
# DEPMAP FEATURE GENERATION
# =============================================================================

def summarise_matrix_by_gene(
    path: Path,
    prefix: str,
    dependency_thresholds: Optional[List[float]] = None,
    expression_threshold: Optional[float] = None,
    mutation_binary: bool = False,
    copy_number: bool = False,
) -> pd.DataFrame:
    """
    Convert a model x gene matrix into a gene-level table.

    Output has one row per parsed gene column.
    """
    if not path.exists():
        log(f"[SKIP MATRIX] Missing file: {path}")
        return pd.DataFrame()

    df = read_numeric_matrix(path)

    id_col = df.columns[0]
    gene_cols = [c for c in df.columns if c != id_col]

    records = []
    log(f"[MATRIX SHAPE] {path.name}: rows={df.shape[0]}, gene_columns={len(gene_cols)}")

    for i, col in enumerate(gene_cols, start=1):
        symbol, entrez, raw_col = parse_depmap_gene_column(col)
        if not symbol:
            continue

        s = pd.to_numeric(df[col], errors="coerce")
        n_total = int(len(s))
        n_non_missing = int(s.notna().sum())

        rec = {
            "gene_symbol": symbol,
            "depmap_entrez_id": entrez or "",
            "depmap_raw_gene_column": raw_col,
            f"{prefix}_n_models_total": n_total,
            f"{prefix}_n_models_nonmissing": n_non_missing,
            f"{prefix}_missing_fraction": float(1.0 - (n_non_missing / n_total)) if n_total else np.nan,
            f"{prefix}_mean": float(s.mean(skipna=True)) if n_non_missing else np.nan,
            f"{prefix}_median": float(s.median(skipna=True)) if n_non_missing else np.nan,
            f"{prefix}_std": float(s.std(skipna=True)) if n_non_missing > 1 else np.nan,
            f"{prefix}_min": float(s.min(skipna=True)) if n_non_missing else np.nan,
            f"{prefix}_max": float(s.max(skipna=True)) if n_non_missing else np.nan,
            f"{prefix}_q05": float(s.quantile(0.05)) if n_non_missing else np.nan,
            f"{prefix}_q25": float(s.quantile(0.25)) if n_non_missing else np.nan,
            f"{prefix}_q75": float(s.quantile(0.75)) if n_non_missing else np.nan,
            f"{prefix}_q95": float(s.quantile(0.95)) if n_non_missing else np.nan,
        }

        if dependency_thresholds:
            for thr in dependency_thresholds:
                safe_thr = str(thr).replace(".", "p")
                rec[f"{prefix}_fraction_gt_{safe_thr}"] = float((s > thr).mean(skipna=True)) if n_non_missing else np.nan

        if expression_threshold is not None:
            safe_thr = str(expression_threshold).replace(".", "p")
            rec[f"{prefix}_fraction_gt_{safe_thr}"] = float((s > expression_threshold).mean(skipna=True)) if n_non_missing else np.nan
            rec[f"{prefix}_fraction_expressed"] = rec[f"{prefix}_fraction_gt_{safe_thr}"]

        if mutation_binary:
            # Mutation matrices are usually 0/1. This also works if values are boolean-like.
            rec[f"{prefix}_fraction_altered"] = float((s > 0).mean(skipna=True)) if n_non_missing else np.nan
            rec[f"{prefix}_count_altered_models"] = int((s > 0).sum(skipna=True)) if n_non_missing else 0

        if copy_number:
            # Conservative generic CN summaries. Thresholds can be adjusted later.
            rec[f"{prefix}_fraction_gain_gt_0p3"] = float((s > 0.3).mean(skipna=True)) if n_non_missing else np.nan
            rec[f"{prefix}_fraction_loss_lt_minus_0p3"] = float((s < -0.3).mean(skipna=True)) if n_non_missing else np.nan
            rec[f"{prefix}_fraction_amp_gt_1"] = float((s > 1.0).mean(skipna=True)) if n_non_missing else np.nan
            rec[f"{prefix}_fraction_del_lt_minus_1"] = float((s < -1.0).mean(skipna=True)) if n_non_missing else np.nan

        records.append(rec)

        if i % 5000 == 0:
            log(f"  processed {i}/{len(gene_cols)} gene columns for {path.name}")

    out = pd.DataFrame(records)

    if out.empty:
        return out

    # If duplicated symbols exist, keep the first raw gene column and average numeric features.
    # Duplicates are uncommon but can occur with aliases/ambiguous identifiers.
    non_numeric_cols = ["gene_symbol", "depmap_entrez_id", "depmap_raw_gene_column"]
    numeric_cols = [c for c in out.columns if c not in non_numeric_cols]
    agg = {c: "mean" for c in numeric_cols}
    agg["depmap_entrez_id"] = "first"
    agg["depmap_raw_gene_column"] = "first"
    out = out.groupby("gene_symbol", as_index=False).agg(agg)

    log(f"[FEATURES] {path.name}: {out.shape[0]} genes x {out.shape[1]} columns")
    return out


def summarise_fusion_file(path: Path) -> pd.DataFrame:
    """
    Make gene-level fusion features if OmicsFusionFiltered.csv exists.

    This file format can vary. The function looks for likely gene columns and
    counts the number of fusion records per gene.
    """
    if not path.exists():
        log(f"[SKIP FUSION] Missing file: {path}")
        return pd.DataFrame()

    log(f"[READ FUSION] {path}")
    df = pd.read_csv(path, dtype=str, low_memory=False)

    possible_gene_cols = [
        "LeftGene", "RightGene", "left_gene", "right_gene",
        "Gene1", "Gene2", "gene1", "gene2",
        "HugoSymbol", "Hugo_Symbol", "gene", "Gene",
        "gene_symbol", "GeneSymbol",
    ]

    gene_cols = [c for c in possible_gene_cols if c in df.columns]

    if not gene_cols:
        # Fallback: any column with gene in name but avoid IDs if possible.
        gene_cols = [c for c in df.columns if "gene" in c.lower() and "id" not in c.lower()]

    if not gene_cols:
        log("[FUSION WARNING] Could not find gene columns. Skipping fusion features.")
        log(f"[FUSION COLUMNS] {df.columns.tolist()}")
        return pd.DataFrame()

    all_genes = []
    for c in gene_cols:
        vals = df[c].dropna().astype(str).str.split(r"[|,;]", regex=True).explode()
        vals = vals.map(normalize_gene_symbol)
        vals = vals[vals != ""]
        all_genes.append(vals)

    if not all_genes:
        return pd.DataFrame()

    genes = pd.concat(all_genes, ignore_index=True)
    counts = genes.value_counts().rename_axis("gene_symbol").reset_index(name="depmap_fusion_record_count")
    counts["depmap_fusion_has_evidence"] = (counts["depmap_fusion_record_count"] > 0).astype(int)

    log(f"[FUSION FEATURES] {counts.shape[0]} genes")
    return counts


def build_depmap_gene_features(downloads_dir: Path, processed_dir: Path, include_optional: bool) -> pd.DataFrame:
    mkdir(processed_dir)

    feature_tables = []

    # Core matrices.
    for fname, prefix in MATRIX_FILES.items():
        if fname in OPTIONAL_DEPMAP_FILES and not include_optional:
            continue

        path = downloads_dir / fname

        if fname == "CRISPRGeneDependency.csv":
            ft = summarise_matrix_by_gene(
                path,
                prefix=prefix,
                dependency_thresholds=[0.1, 0.2, 0.5, 0.9],
            )
        elif fname == "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv":
            ft = summarise_matrix_by_gene(
                path,
                prefix=prefix,
                expression_threshold=1.0,
            )
        elif fname in ["OmicsSomaticMutationsMatrixHotspot.csv", "OmicsSomaticMutationsMatrixDamaging.csv"]:
            ft = summarise_matrix_by_gene(
                path,
                prefix=prefix,
                mutation_binary=True,
            )
        elif fname == "OmicsCNGeneWGS.csv":
            ft = summarise_matrix_by_gene(
                path,
                prefix=prefix,
                copy_number=True,
            )
        else:
            ft = summarise_matrix_by_gene(
                path,
                prefix=prefix,
            )

        if not ft.empty:
            feature_tables.append(ft)

    # Optional fusion evidence.
    if include_optional:
        fusion_ft = summarise_fusion_file(downloads_dir / FUSION_FILE)
        if not fusion_ft.empty:
            feature_tables.append(fusion_ft)

    if not feature_tables:
        raise RuntimeError(
            "No DepMap feature tables were generated. Check that the files downloaded correctly."
        )

    # Merge all feature tables on gene_symbol.
    merged = feature_tables[0]
    for ft in feature_tables[1:]:
        # Avoid repeated mapping columns except gene_symbol.
        dup_cols = [c for c in ft.columns if c in merged.columns and c != "gene_symbol"]
        if dup_cols:
            ft = ft.drop(columns=dup_cols)
        merged = merged.merge(ft, on="gene_symbol", how="outer")

    merged = merged.sort_values("gene_symbol").reset_index(drop=True)

    outpath = processed_dir / "feature1_depmap_gene_features.csv"
    merged.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED DEPMAP GENE FEATURES]")
    log(f"[PATH]  {outpath}")
    log(f"[SHAPE] {merged.shape[0]} genes x {merged.shape[1]} columns")

    return merged


def add_depmap_gene_mapping_features(downloads_dir: Path, depmap_features: pd.DataFrame) -> pd.DataFrame:
    """
    Add optional mapping columns from DepMap Gene.csv where possible.
    """
    gene_path = downloads_dir / "Gene.csv"
    if not gene_path.exists():
        log("[GENE MAPPING] Gene.csv missing; skipping DepMap gene metadata merge.")
        return depmap_features

    log(f"[READ DEPMAP GENE] {gene_path}")
    gene_df = pd.read_csv(gene_path, dtype=str, low_memory=False)
    gene_df.columns = [c.strip() for c in gene_df.columns]

    # Identify likely symbol column.
    symbol_col = get_first_existing_col(gene_df, [
        "gene_symbol", "GeneSymbol", "symbol", "Symbol", "gene_name", "GeneName"
    ])

    if symbol_col is None:
        # Common DepMap Gene.csv may have "label" or "name" depending on release.
        symbol_col = get_first_existing_col(gene_df, ["label", "Label", "name", "Name"])

    if symbol_col is None:
        log("[GENE MAPPING WARNING] Could not identify symbol column in Gene.csv.")
        log(f"[GENE.CSV COLUMNS] {gene_df.columns.tolist()}")
        return depmap_features

    gene_df["gene_symbol"] = gene_df[symbol_col].map(normalize_gene_symbol)

    keep_cols = ["gene_symbol"]
    for c in gene_df.columns:
        cl = c.lower()
        if c == "gene_symbol":
            continue
        if any(k in cl for k in ["entrez", "ensembl", "uniprot", "gene_id", "symbol"]):
            keep_cols.append(c)

    keep_cols = list(dict.fromkeys(keep_cols))
    gene_small = gene_df[keep_cols].drop_duplicates(subset=["gene_symbol"], keep="first").copy()

    # Prefix non-key columns.
    rename = {c: f"depmap_genecsv_{c}" for c in gene_small.columns if c != "gene_symbol"}
    gene_small = gene_small.rename(columns=rename)

    out = depmap_features.merge(gene_small, on="gene_symbol", how="left")
    return out


def merge_with_hgnc(hgnc: pd.DataFrame, depmap_features: pd.DataFrame, processed_dir: Path) -> pd.DataFrame:
    log("=" * 100)
    log("[MERGE HGNC + DEPMAP FEATURES]")

    depmap_features = depmap_features.copy()
    depmap_features["gene_symbol"] = depmap_features["gene_symbol"].map(normalize_gene_symbol)

    merged = hgnc.merge(depmap_features, on="gene_symbol", how="left")

    feature_cols = [c for c in depmap_features.columns if c != "gene_symbol"]
    merged["feature1_depmap_has_any_feature"] = merged[feature_cols].notna().any(axis=1).astype(int)
    merged["feature1_depmap_n_nonmissing_features"] = merged[feature_cols].notna().sum(axis=1)

    outpath = processed_dir / "feature1_depmap_hgnc_merged.csv"
    merged.to_csv(outpath, index=False)

    coverage = merged["feature1_depmap_has_any_feature"].mean() if len(merged) else np.nan

    log(f"[SAVED] {outpath}")
    log(f"[SHAPE] {merged.shape[0]} genes x {merged.shape[1]} columns")
    log(f"[COVERAGE] HGNC genes with at least one DepMap feature: {coverage:.3f}")

    return merged


# =============================================================================
# REPORTING
# =============================================================================

def write_summary_report(
    outdir: Path,
    hgnc: Optional[pd.DataFrame],
    depmap_features: Optional[pd.DataFrame],
    merged: Optional[pd.DataFrame],
    include_optional: bool,
    release: str,
) -> None:
    report_path = outdir / "feature1_depmap_summary.txt"

    lines = []
    lines.append("Feature 1 DepMap summary")
    lines.append("=" * 80)
    lines.append(f"Created at: {now_iso()}")
    lines.append(f"Release: {release}")
    lines.append(f"Output directory: {outdir.resolve()}")
    lines.append("")
    lines.append("Leakage rule:")
    lines.append("  Included: CRISPR essentiality/dependency, expression, optional mutation/CN/fusion features.")
    lines.append("  Excluded: compound-target/drug-response files such as PortalCompounds.csv.")
    lines.append("")
    lines.append("Required files:")
    for f in REQUIRED_DEPMAP_FILES:
        lines.append(f"  - {f}")
    lines.append("")
    lines.append(f"Optional files included: {include_optional}")
    if include_optional:
        for f in OPTIONAL_DEPMAP_FILES:
            lines.append(f"  - {f}")
    lines.append("")
    lines.append("Excluded leakage-risk files:")
    for f in LEAKAGE_RISK_FILES:
        lines.append(f"  - {f}")
    lines.append("")

    if hgnc is not None:
        lines.append(f"HGNC rows: {hgnc.shape[0]}")
    if depmap_features is not None:
        lines.append(f"DepMap feature rows: {depmap_features.shape[0]}")
        lines.append(f"DepMap feature columns: {depmap_features.shape[1]}")
    if merged is not None:
        lines.append(f"Merged rows: {merged.shape[0]}")
        lines.append(f"Merged columns: {merged.shape[1]}")
        if "feature1_depmap_has_any_feature" in merged.columns:
            lines.append(
                "Merged coverage: "
                f"{merged['feature1_depmap_has_any_feature'].mean():.4f}"
            )

    report_path.write_text("\n".join(lines) + "\n")
    log(f"[REPORT] {report_path}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download DepMap data and build HGNC-level DepMap gene features."
    )
    parser.add_argument("--release", default=DEFAULT_RELEASE, help="DepMap release name.")
    parser.add_argument("--hgnc", default=str(DEFAULT_HGNC_PATH), help="Path to hgnc_complete_set.txt.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory.")
    parser.add_argument("--include-optional", action="store_true", help="Download and process optional mutation/CN/fusion files.")
    parser.add_argument("--all-hgnc-genes", action="store_true", help="Use all HGNC rows instead of protein-coding only.")
    parser.add_argument("--download-only", action="store_true", help="Only download files; do not build features.")
    parser.add_argument("--no-download", action="store_true", help="Do not download; use existing files in downloads folder.")
    parser.add_argument("--force-index", action="store_true", help="Re-download depmap_files.csv index.")
    parser.add_argument("--force-files", action="store_true", help="Re-download DepMap files even if they already exist.")
    args = parser.parse_args()

    outdir = mkdir(Path(args.outdir))
    downloads_dir = mkdir(outdir / "downloads")
    processed_dir = mkdir(outdir / "processed")

    log("=" * 100)
    log("FEATURE 1: DEPMAP FUNCTIONAL-GENOMICS FEATURES")
    log("=" * 100)
    log(f"[RELEASE]          {args.release}")
    log(f"[HGNC]             {args.hgnc}")
    log(f"[OUTDIR]           {outdir.resolve()}")
    log(f"[DOWNLOADS]        {downloads_dir.resolve()}")
    log(f"[PROCESSED]        {processed_dir.resolve()}")
    log(f"[INCLUDE OPTIONAL] {args.include_optional}")
    log(f"[DOWNLOAD ONLY]    {args.download_only}")
    log(f"[NO DOWNLOAD]      {args.no_download}")
    log("=" * 100)

    if not args.no_download:
        download_depmap_files(
            release=args.release,
            downloads_dir=downloads_dir,
            include_optional=args.include_optional,
            force_index=args.force_index,
            force_files=args.force_files,
        )

    if args.download_only:
        write_summary_report(
            outdir=outdir,
            hgnc=None,
            depmap_features=None,
            merged=None,
            include_optional=args.include_optional,
            release=args.release,
        )
        log("[DONE] Download-only mode finished.")
        return

    hgnc = load_hgnc(
        hgnc_path=Path(args.hgnc),
        protein_coding_only=not args.all_hgnc_genes,
    )

    depmap_features = build_depmap_gene_features(
        downloads_dir=downloads_dir,
        processed_dir=processed_dir,
        include_optional=args.include_optional,
    )

    depmap_features = add_depmap_gene_mapping_features(downloads_dir, depmap_features)

    # Save after Gene.csv metadata merge too.
    depmap_feature_path = processed_dir / "feature1_depmap_gene_features.csv"
    depmap_features.to_csv(depmap_feature_path, index=False)

    merged = merge_with_hgnc(
        hgnc=hgnc,
        depmap_features=depmap_features,
        processed_dir=processed_dir,
    )

    write_summary_report(
        outdir=outdir,
        hgnc=hgnc,
        depmap_features=depmap_features,
        merged=merged,
        include_optional=args.include_optional,
        release=args.release,
    )

    log("=" * 100)
    log("[DONE]")
    log(f"[RAW DOWNLOADS]      {downloads_dir.resolve()}")
    log(f"[GENE FEATURES]      {processed_dir / 'feature1_depmap_gene_features.csv'}")
    log(f"[HGNC MERGED TABLE]  {processed_dir / 'feature1_depmap_hgnc_merged.csv'}")
    log("=" * 100)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        log("\n[STOPPED] Interrupted by user.")
        sys.exit(130)
    except Exception as e:
        log("=" * 100)
        log("[ERROR]")
        log(str(e))
        log("=" * 100)
        raise