#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature19_gnomAD_full.py

Feature 19: Full gnomAD gene-constraint features.

Purpose
-------
Build gene-level gnomAD constraint features using the gnomAD constraint file
already downloaded for Feature 9.

This expands Feature 9 by extracting the full constraint set:
    - pLI
    - LOEUF / oe_lof_upper
    - observed/expected LoF
    - observed/expected missense
    - observed/expected synonymous
    - missense Z
    - synonymous Z
    - LoF Z
    - observed/expected counts
    - derived constraint flags

Default expected input locations
--------------------------------
feature_databases/gnomAD/
feature_databases/GnomAD/
feature_databases/Constraint/
feature9_disorder_constraint/

Common file names:
    gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
    gnomad.v2.1.1.lof_metrics.by_gene.txt.gz
    gnomad.v2.1.1.lof_metrics.by_gene.txt
    any file containing "gnomad", "lof_metrics", or "constraint"

HGNC:
    databases/HGNC/hgnc_complete_set.txt

Outputs
-------
feature19_gnomad_full/
    processed/feature19_gnomad_constraint_raw_mapped.csv
    processed/feature19_gnomad_full_gene_features.csv
    processed/feature19_gnomad_full_hgnc_merged.csv
    feature19_gnomad_full_summary.txt
    feature19_gnomad_full_run_metadata.json

Run
---
    python Feature19_gnomAD_full.py

Manual input:
    python Feature19_gnomAD_full.py \
      --gnomad-file feature_databases/gnomAD/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz

Fast test:
    python Feature19_gnomAD_full.py --limit-genes 1000

Leakage policy
--------------
Safe. gnomAD constraint metrics are population genetic constraint features.
They are not drug labels and do not use ChEMBL, DrugBank, DGIdb, Open Targets,
Pharos, or clinical target labels.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# =============================================================================
# DEFAULTS
# =============================================================================

DEFAULT_HGNC_PATH = Path("databases") / "HGNC" / "hgnc_complete_set.txt"
DEFAULT_DBDIR = Path("feature_databases")
DEFAULT_OUTDIR = Path("feature19_gnomad_full")


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
    s = s.replace("gene:", "")
    s = s.split(".")[0]
    return s


def safe_float(x: Any) -> float:
    try:
        if x is None or pd.isna(x):
            return np.nan
        s = str(x).strip()
        if s in {"", ".", "NA", "NaN", "nan", "None", "null"}:
            return np.nan
        return float(s)
    except Exception:
        return np.nan


def open_text_maybe_gzip(path: Path):
    name = str(path).lower()
    if name.endswith(".gz") or name.endswith(".bgz"):
        return gzip.open(path, "rt", errors="ignore")
    return open(path, "rt", errors="ignore")


def read_table_flexible(path: Path, max_rows: int = 0) -> pd.DataFrame:
    """
    Read TSV/CSV/TXT/GZ/BGZ robustly.

    Important:
    gnomAD .bgz files are gzip-compatible, but pandas compression='infer'
    does not always infer .bgz correctly. Therefore we explicitly open
    .gz/.bgz files with gzip.open.
    """
    if not path.exists():
        raise FileNotFoundError(path)

    nrows = max_rows if max_rows and max_rows > 0 else None
    name = path.name.lower()

    # Explicit gzip/BGZF handling.
    if name.endswith(".gz") or name.endswith(".bgz"):
        with gzip.open(path, "rt", errors="ignore") as f:
            if ".csv" in name:
                return pd.read_csv(
                    f,
                    dtype=str,
                    low_memory=False,
                    nrows=nrows,
                )
            return pd.read_csv(
                f,
                sep="\t",
                dtype=str,
                low_memory=False,
                nrows=nrows,
            )

    # Plain CSV.
    if name.endswith(".csv"):
        return pd.read_csv(
            path,
            dtype=str,
            low_memory=False,
            nrows=nrows,
        )

    # Default plain TSV/TXT.
    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        low_memory=False,
        nrows=nrows,
    )


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


def sanitize_feature_name(x: str) -> str:
    x = str(x).strip()
    x = x.replace(".", "_")
    x = re.sub(r"[^A-Za-z0-9]+", "_", x)
    x = x.strip("_").lower()
    return x


# =============================================================================
# FILE DETECTION
# =============================================================================

def find_gnomad_constraint_file(dbdir: Path) -> Optional[Path]:
    roots = [
        dbdir / "gnomAD",
        dbdir / "GnomAD",
        dbdir / "gnomad",
        dbdir / "Constraint",
        dbdir / "constraint",
        Path("feature9_disorder_constraint"),
        Path("feature9_disorder_constraint") / "downloads",
        Path("feature9_disorder_constraint") / "processed",
        dbdir,
    ]

    patterns = [
        "*gnomad*lof_metrics*.txt",
        "*gnomad*lof_metrics*.txt.gz",
        "*gnomad*lof_metrics*.bgz",
        "*lof_metrics*.txt",
        "*lof_metrics*.txt.gz",
        "*lof_metrics*.bgz",
        "*gnomad*constraint*.tsv",
        "*gnomad*constraint*.tsv.gz",
        "*gnomad*constraint*.txt",
        "*gnomad*constraint*.txt.gz",
        "*constraint*.tsv",
        "*constraint*.tsv.gz",
        "*constraint*.txt",
        "*constraint*.txt.gz",
        "*constraint*.bgz",
    ]

    files = []

    for root in roots:
        if not root.exists():
            continue
        for pat in patterns:
            files.extend(list(root.rglob(pat)))

    files = [p for p in files if p.exists() and p.stat().st_size > 0]

    if not files:
        return None

    ranked = []

    for p in files:
        name = p.name.lower()
        score = 0

        if "gnomad" in name:
            score += 30
        if "v2.1.1" in name or "2.1.1" in name:
            score += 30
        if "lof_metrics" in name:
            score += 30
        if "by_gene" in name:
            score += 20
        if name.endswith(".bgz") or name.endswith(".gz"):
            score += 5
        if "feature9" in str(p).lower():
            score += 2

        ranked.append((score, p.stat().st_size, p))

    ranked.sort(reverse=True)

    return ranked[0][2]


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
        raise RuntimeError("HGNC file must contain symbol column.")

    before = len(hgnc)

    if protein_coding_only:
        if "locus_group" in hgnc.columns:
            hgnc = hgnc[hgnc["locus_group"].astype(str).str.lower().eq("protein-coding gene")].copy()
        elif "locus_type" in hgnc.columns:
            hgnc = hgnc[hgnc["locus_type"].astype(str).str.lower().str.contains("protein", na=False)].copy()

    hgnc["gene_symbol"] = hgnc["symbol"].map(normalize_symbol)

    for col in ["alias_symbol", "prev_symbol", "ensembl_gene_id", "entrez_id", "uniprot_ids", "name"]:
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
    """
    Returns:
        symbol_alias_to_current_symbol
        ensembl_to_current_symbol
    """
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


# =============================================================================
# GNOMAD COLUMN NORMALISATION
# =============================================================================

def get_metric_columns(df: pd.DataFrame) -> Dict[str, Optional[str]]:
    """
    Detect common gnomAD v2/v4 constraint columns.

    gnomAD v2.1.1 usually has:
        gene
        gene_id
        transcript
        canonical
        obs_mis
        exp_mis
        oe_mis
        mu_mis
        possible_mis
        obs_mis_pphen
        exp_mis_pphen
        oe_mis_pphen
        possible_mis_pphen
        obs_syn
        exp_syn
        oe_syn
        mu_syn
        possible_syn
        obs_lof
        mu_lof
        possible_lof
        exp_lof
        pLI
        pRec
        pNull
        oe_lof
        oe_syn_lower
        oe_syn_upper
        oe_mis_lower
        oe_mis_upper
        oe_lof_lower
        oe_lof_upper
        constraint_flag
        syn_z
        mis_z
        lof_z
    """
    cols = {str(c).lower(): c for c in df.columns}

    def exact(*names: str) -> Optional[str]:
        for name in names:
            if name.lower() in cols:
                return cols[name.lower()]
        return None

    def contains_all(*tokens: str) -> Optional[str]:
        for c in df.columns:
            low = str(c).lower()
            if all(t.lower() in low for t in tokens):
                return c
        return None

    m = {
        "gene_symbol_raw": exact("gene", "gene_symbol", "symbol"),
        "ensembl_gene_id_raw": exact("gene_id", "ensembl_gene_id", "ensg"),
        "transcript": exact("transcript", "transcript_id"),
        "canonical": exact("canonical"),

        "pli": exact("pLI", "pli", "lof_pLI"),
        "prec": exact("pRec", "prec"),
        "pnull": exact("pNull", "pnull"),

        "oe_lof": exact("oe_lof", "lof.oe", "lof_oe"),
        "oe_lof_lower": exact("oe_lof_lower", "lof.oe_ci.lower", "lof_oe_ci_lower"),
        "oe_lof_upper": exact("oe_lof_upper", "lof.oe_ci.upper", "lof_oe_ci_upper"),
        "loeuf": exact("loeuf", "LOEUF") or contains_all("loeuf"),

        "oe_mis": exact("oe_mis", "mis.oe", "mis_oe"),
        "oe_mis_lower": exact("oe_mis_lower", "mis.oe_ci.lower", "mis_oe_ci_lower"),
        "oe_mis_upper": exact("oe_mis_upper", "mis.oe_ci.upper", "mis_oe_ci_upper"),

        "oe_syn": exact("oe_syn", "syn.oe", "syn_oe"),
        "oe_syn_lower": exact("oe_syn_lower", "syn.oe_ci.lower", "syn_oe_ci_lower"),
        "oe_syn_upper": exact("oe_syn_upper", "syn.oe_ci.upper", "syn_oe_ci_upper"),

        "obs_lof": exact("obs_lof", "lof.obs", "lof_obs"),
        "exp_lof": exact("exp_lof", "lof.exp", "lof_exp"),
        "possible_lof": exact("possible_lof"),
        "mu_lof": exact("mu_lof"),

        "obs_mis": exact("obs_mis", "mis.obs", "mis_obs"),
        "exp_mis": exact("exp_mis", "mis.exp", "mis_exp"),
        "possible_mis": exact("possible_mis"),
        "mu_mis": exact("mu_mis"),

        "obs_mis_pphen": exact("obs_mis_pphen"),
        "exp_mis_pphen": exact("exp_mis_pphen"),
        "oe_mis_pphen": exact("oe_mis_pphen"),
        "possible_mis_pphen": exact("possible_mis_pphen"),

        "obs_syn": exact("obs_syn", "syn.obs", "syn_obs"),
        "exp_syn": exact("exp_syn", "syn.exp", "syn_exp"),
        "possible_syn": exact("possible_syn"),
        "mu_syn": exact("mu_syn"),

        "mis_z": exact("mis_z", "missense_z", "mis_z_score"),
        "syn_z": exact("syn_z", "synonymous_z", "syn_z_score"),
        "lof_z": exact("lof_z", "pLoF_z", "lof_z_score"),

        "constraint_flag": exact("constraint_flag"),
    }

    return m


def standardize_gnomad_constraint(
    df: pd.DataFrame,
    symbol_map: Dict[str, str],
    ensembl_map: Dict[str, str],
) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]

    metric_cols = get_metric_columns(df)

    log("=" * 100)
    log("[DETECTED GNOMAD COLUMNS]")
    for k, v in metric_cols.items():
        log(f"{k:25s} = {v}")

    gene_symbol_col = metric_cols["gene_symbol_raw"]
    ensembl_col = metric_cols["ensembl_gene_id_raw"]

    rows = []

    for _, row in df.iterrows():
        gene = ""

        raw_symbol = normalize_symbol(row.get(gene_symbol_col, "")) if gene_symbol_col else ""
        raw_ens = clean_ensembl_gene_id(row.get(ensembl_col, "")) if ensembl_col else ""

        if raw_symbol and raw_symbol in symbol_map:
            gene = symbol_map[raw_symbol]
        elif raw_ens and raw_ens in ensembl_map:
            gene = ensembl_map[raw_ens]

        if not gene:
            continue

        rec = {
            "gene_symbol": gene,
            "feature19_gnomad_raw_gene_symbol": raw_symbol,
            "feature19_gnomad_raw_ensembl_gene_id": raw_ens,
        }

        # Keep raw text metadata separately.
        for meta_key in ["transcript", "canonical", "constraint_flag"]:
            c = metric_cols.get(meta_key)
            if c:
                rec[f"feature19_gnomad_{meta_key}"] = clean_text(row.get(c, ""))

        # Numeric metrics.
        numeric_keys = [
            "pli", "prec", "pnull",
            "oe_lof", "oe_lof_lower", "oe_lof_upper", "loeuf",
            "oe_mis", "oe_mis_lower", "oe_mis_upper",
            "oe_syn", "oe_syn_lower", "oe_syn_upper",
            "obs_lof", "exp_lof", "possible_lof", "mu_lof",
            "obs_mis", "exp_mis", "possible_mis", "mu_mis",
            "obs_mis_pphen", "exp_mis_pphen", "oe_mis_pphen", "possible_mis_pphen",
            "obs_syn", "exp_syn", "possible_syn", "mu_syn",
            "mis_z", "syn_z", "lof_z",
        ]

        for key in numeric_keys:
            c = metric_cols.get(key)
            if c:
                rec[f"feature19_gnomad_{key}"] = safe_float(row.get(c, np.nan))

        rows.append(rec)

    mapped = pd.DataFrame(rows)

    if mapped.empty:
        return mapped

    # Drop exact duplicate rows.
    mapped = mapped.drop_duplicates()

    return mapped


# =============================================================================
# GENE-LEVEL AGGREGATION
# =============================================================================

def aggregate_to_gene_level(hgnc: pd.DataFrame, mapped: pd.DataFrame, processed_dir: Path) -> pd.DataFrame:
    genes = pd.DataFrame({"gene_symbol": sorted(hgnc["gene_symbol"].unique().tolist())})

    if mapped.empty:
        out = genes.copy()
        out["feature19_gnomad_has_constraint"] = 0
        out["feature19_gnomad_constraint_row_count"] = 0
        out.to_csv(processed_dir / "feature19_gnomad_full_gene_features.csv", index=False)
        return out

    numeric_cols = [
        c for c in mapped.columns
        if c.startswith("feature19_gnomad_")
        and c not in {
            "feature19_gnomad_raw_gene_symbol",
            "feature19_gnomad_raw_ensembl_gene_id",
            "feature19_gnomad_transcript",
            "feature19_gnomad_canonical",
            "feature19_gnomad_constraint_flag",
        }
    ]

    for c in numeric_cols:
        mapped[c] = pd.to_numeric(mapped[c], errors="coerce")

    records = []
    by_gene = dict(tuple(mapped.groupby("gene_symbol")))

    for gene in genes["gene_symbol"]:
        g = by_gene.get(gene, pd.DataFrame())

        if g.empty:
            records.append(
                {
                    "gene_symbol": gene,
                    "feature19_gnomad_has_constraint": 0,
                    "feature19_gnomad_constraint_row_count": 0,
                }
            )
            continue

        # Prefer canonical row if canonical exists and has value true/1.
        canonical_g = g.copy()
        if "feature19_gnomad_canonical" in g.columns:
            canon = g["feature19_gnomad_canonical"].astype(str).str.lower()
            tmp = g[canon.isin(["true", "1", "yes", "y"])]
            if not tmp.empty:
                canonical_g = tmp.copy()

        rec = {
            "gene_symbol": gene,
            "feature19_gnomad_has_constraint": 1,
            "feature19_gnomad_constraint_row_count": int(len(g)),
            "feature19_gnomad_canonical_row_count": int(len(canonical_g)),
            "feature19_gnomad_unique_transcript_count": int(g["feature19_gnomad_transcript"].replace("", np.nan).dropna().nunique()) if "feature19_gnomad_transcript" in g.columns else 0,
        }

        # For gnomAD by-gene table, usually one row per gene/transcript.
        # Use median for duplicate rows; also save min/max where useful.
        for c in numeric_cols:
            vals = pd.to_numeric(canonical_g[c], errors="coerce").dropna()

            if len(vals) == 0:
                vals = pd.to_numeric(g[c], errors="coerce").dropna()

            if len(vals) == 0:
                continue

            rec[c] = float(vals.median())

            # Keep min/max for duplicated transcript cases.
            rec[f"{c}_min"] = float(vals.min())
            rec[f"{c}_max"] = float(vals.max())

        # Constraint flag as text count.
        if "feature19_gnomad_constraint_flag" in g.columns:
            flags = g["feature19_gnomad_constraint_flag"].replace("", np.nan).dropna().astype(str).tolist()
            rec["feature19_gnomad_has_constraint_flag_text"] = int(len(flags) > 0)
            rec["feature19_gnomad_constraint_flag_unique_count"] = int(pd.Series(flags).nunique()) if flags else 0
        else:
            rec["feature19_gnomad_has_constraint_flag_text"] = 0
            rec["feature19_gnomad_constraint_flag_unique_count"] = 0

        records.append(rec)

    features = pd.DataFrame(records)

    # Derived clinically/interpretable constraint flags.
    # Lower LOEUF / oe_lof_upper = stronger LoF constraint.
    loeuf_candidate_cols = [
        "feature19_gnomad_loeuf",
        "feature19_gnomad_oe_lof_upper",
    ]

    loeuf_col = None
    for c in loeuf_candidate_cols:
        if c in features.columns:
            loeuf_col = c
            break

    if loeuf_col:
        x = pd.to_numeric(features[loeuf_col], errors="coerce")
        features["feature19_gnomad_lof_constrained_loeuf_le_0_35"] = (x <= 0.35).astype(int)
        features["feature19_gnomad_lof_constrained_loeuf_le_0_60"] = (x <= 0.60).astype(int)
        features["feature19_gnomad_lof_tolerant_loeuf_ge_1_00"] = (x >= 1.00).astype(int)
        features["feature19_gnomad_lof_constraint_score_inv_loeuf"] = 1.0 / (x + 1e-6)

    if "feature19_gnomad_pli" in features.columns:
        x = pd.to_numeric(features["feature19_gnomad_pli"], errors="coerce")
        features["feature19_gnomad_pli_ge_0_90"] = (x >= 0.90).astype(int)
        features["feature19_gnomad_pli_ge_0_99"] = (x >= 0.99).astype(int)

    if "feature19_gnomad_mis_z" in features.columns:
        x = pd.to_numeric(features["feature19_gnomad_mis_z"], errors="coerce")
        features["feature19_gnomad_missense_constrained_mis_z_ge_3_09"] = (x >= 3.09).astype(int)
        features["feature19_gnomad_missense_constrained_mis_z_ge_2"] = (x >= 2.0).astype(int)

    if "feature19_gnomad_syn_z" in features.columns:
        x = pd.to_numeric(features["feature19_gnomad_syn_z"], errors="coerce")
        features["feature19_gnomad_syn_z_abs"] = x.abs()

    if "feature19_gnomad_lof_z" in features.columns:
        x = pd.to_numeric(features["feature19_gnomad_lof_z"], errors="coerce")
        features["feature19_gnomad_lof_z_ge_2"] = (x >= 2.0).astype(int)

    # Observed/expected derived ratios if missing.
    for prefix in ["lof", "mis", "syn"]:
        obs = f"feature19_gnomad_obs_{prefix}"
        exp = f"feature19_gnomad_exp_{prefix}"
        oe = f"feature19_gnomad_oe_{prefix}"

        if oe not in features.columns and obs in features.columns and exp in features.columns:
            obs_x = pd.to_numeric(features[obs], errors="coerce")
            exp_x = pd.to_numeric(features[exp], errors="coerce")
            features[oe] = obs_x / exp_x.replace(0, np.nan)

    # Combined scores.
    score_parts = []

    if "feature19_gnomad_pli" in features.columns:
        score_parts.append(pd.to_numeric(features["feature19_gnomad_pli"], errors="coerce").fillna(0))

    if loeuf_col:
        loeuf = pd.to_numeric(features[loeuf_col], errors="coerce")
        inv = 1.0 / (loeuf + 1e-6)
        inv = inv.replace([np.inf, -np.inf], np.nan).fillna(0)
        # Scale roughly to 0-1 by dividing by 10 and clipping.
        inv_scaled = np.clip(inv / 10.0, 0, 1)
        score_parts.append(inv_scaled)

    if "feature19_gnomad_mis_z" in features.columns:
        misz = pd.to_numeric(features["feature19_gnomad_mis_z"], errors="coerce").fillna(0)
        misz_scaled = np.clip(misz / 10.0, 0, 1)
        score_parts.append(misz_scaled)

    if score_parts:
        stacked = np.vstack([np.asarray(x, dtype=float) for x in score_parts])
        features["feature19_gnomad_combined_constraint_score"] = np.mean(stacked, axis=0)
    else:
        features["feature19_gnomad_combined_constraint_score"] = np.nan

    # Fill counts/flags.
    for c in features.columns:
        if c.startswith("feature19_") and (
            "_count" in c
            or "_has_" in c
            or "_ge_" in c
            or "_le_" in c
            or "_tolerant_" in c
            or "_constrained_" in c
        ):
            if features[c].dtype != object:
                features[c] = features[c].fillna(0)

    outpath = processed_dir / "feature19_gnomad_full_gene_features.csv"
    features.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED GENE FEATURES]")
    log(f"[PATH]  {outpath}")
    log(f"[SHAPE] {features.shape}")
    log(f"[COVERAGE] {features['feature19_gnomad_has_constraint'].mean():.4f}")

    return features


# =============================================================================
# MERGE / SUMMARY
# =============================================================================

def merge_with_hgnc(hgnc: pd.DataFrame, features: pd.DataFrame, processed_dir: Path) -> pd.DataFrame:
    merged = hgnc.merge(features, on="gene_symbol", how="left")

    feature_cols = [c for c in features.columns if c != "gene_symbol"]

    merged["feature19_gnomad_full_has_any_feature"] = merged[feature_cols].notna().any(axis=1).astype(int)
    merged["feature19_gnomad_full_n_nonmissing_features"] = merged[feature_cols].notna().sum(axis=1)

    outpath = processed_dir / "feature19_gnomad_full_hgnc_merged.csv"
    merged.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED HGNC MERGED]")
    log(f"[PATH]  {outpath}")
    log(f"[SHAPE] {merged.shape}")

    return merged


def write_summary(
    outdir: Path,
    hgnc: pd.DataFrame,
    gnomad_file: Path,
    raw: pd.DataFrame,
    mapped: pd.DataFrame,
    features: pd.DataFrame,
    merged: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    path = outdir / "feature19_gnomad_full_summary.txt"

    lines = []
    lines.append("Feature 19: Full gnomAD gene constraint")
    lines.append("=" * 80)
    lines.append(f"Created at: {now_iso()}")
    lines.append("")
    lines.append("Inputs:")
    lines.append(f"HGNC: {args.hgnc}")
    lines.append(f"gnomAD constraint file: {gnomad_file}")
    lines.append("")
    lines.append("Rows:")
    lines.append(f"HGNC rows: {hgnc.shape[0]}")
    lines.append(f"Raw gnomAD rows read: {raw.shape[0]}")
    lines.append(f"Mapped gnomAD rows: {mapped.shape[0]}")
    lines.append(f"Gene feature rows: {features.shape[0]}")
    lines.append(f"Merged rows: {merged.shape[0]}")
    lines.append(f"Merged columns: {merged.shape[1]}")
    lines.append("")
    lines.append("Coverage:")

    if "feature19_gnomad_has_constraint" in features.columns:
        lines.append(f"gnomAD constraint coverage: {features['feature19_gnomad_has_constraint'].mean():.4f}")

    lines.append("")
    lines.append("Key feature summaries:")

    key_cols = [
        "feature19_gnomad_pli",
        "feature19_gnomad_oe_lof",
        "feature19_gnomad_oe_lof_upper",
        "feature19_gnomad_loeuf",
        "feature19_gnomad_mis_z",
        "feature19_gnomad_syn_z",
        "feature19_gnomad_lof_z",
        "feature19_gnomad_oe_mis",
        "feature19_gnomad_oe_syn",
        "feature19_gnomad_combined_constraint_score",
    ]

    for c in key_cols:
        if c in features.columns:
            x = pd.to_numeric(features[c], errors="coerce")
            lines.append(
                f"{c}: nonmissing={int(x.notna().sum())}, "
                f"median={x.median(skipna=True):.4f}, "
                f"min={x.min(skipna=True):.4f}, "
                f"max={x.max(skipna=True):.4f}"
            )

    lines.append("")
    lines.append("Derived flags:")
    for c in [
        "feature19_gnomad_pli_ge_0_90",
        "feature19_gnomad_pli_ge_0_99",
        "feature19_gnomad_lof_constrained_loeuf_le_0_35",
        "feature19_gnomad_lof_constrained_loeuf_le_0_60",
        "feature19_gnomad_missense_constrained_mis_z_ge_3_09",
    ]:
        if c in features.columns:
            lines.append(f"{c}: {int(features[c].sum())}")

    lines.append("")
    lines.append("Interpretation:")
    lines.append("Lower LOEUF / oe_lof_upper indicates stronger loss-of-function constraint.")
    lines.append("Higher pLI indicates stronger intolerance to protein-truncating variation.")
    lines.append("Higher missense Z indicates missense constraint.")
    lines.append("Synonymous Z is a control-style constraint metric and should be interpreted cautiously.")
    lines.append("")
    lines.append("Leakage policy:")
    lines.append("Included: population genetic constraint metrics from gnomAD.")
    lines.append("Excluded: drug labels, clinical target labels, ChEMBL, DrugBank, DGIdb, Open Targets, Pharos.")

    path.write_text("\n".join(lines) + "\n")
    log(f"[SUMMARY] {path}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description="Build Feature 19 full gnomAD constraint features.")
    parser.add_argument("--hgnc", default=str(DEFAULT_HGNC_PATH), help="Path to HGNC complete set.")
    parser.add_argument("--dbdir", default=str(DEFAULT_DBDIR), help="Root feature_databases directory.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory.")
    parser.add_argument("--gnomad-file", default="", help="Manual gnomAD constraint file.")
    parser.add_argument("--all-hgnc-genes", action="store_true", help="Use all HGNC genes, not only protein-coding.")
    parser.add_argument("--limit-genes", type=int, default=0, help="Debug: first N HGNC genes.")
    parser.add_argument("--max-rows", type=int, default=0, help="Debug: first N gnomAD rows.")
    args = parser.parse_args()

    dbdir = Path(args.dbdir)
    outdir = mkdir(Path(args.outdir))
    processed_dir = mkdir(outdir / "processed")

    log("=" * 100)
    log("FEATURE 19: FULL GNOMAD GENE CONSTRAINT")
    log("=" * 100)
    log(f"[HGNC]        {args.hgnc}")
    log(f"[DBDIR]       {dbdir.resolve()}")
    log(f"[OUTDIR]      {outdir.resolve()}")
    log(f"[GNOMAD FILE] {args.gnomad_file if args.gnomad_file else 'auto-detect'}")
    log(f"[LIMIT GENES] {args.limit_genes if args.limit_genes else 'none'}")
    log(f"[MAX ROWS]    {args.max_rows if args.max_rows else 'none'}")
    log("=" * 100)

    gnomad_file = Path(args.gnomad_file) if args.gnomad_file else find_gnomad_constraint_file(dbdir)

    if gnomad_file is None or not gnomad_file.exists():
        raise FileNotFoundError(
            "Could not find gnomAD constraint file. Provide with --gnomad-file.\n"
            "Expected something like:\n"
            "  feature_databases/gnomAD/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
        )

    log(f"[USING GNOMAD FILE] {gnomad_file}")

    hgnc = load_hgnc(
        Path(args.hgnc),
        protein_coding_only=not args.all_hgnc_genes,
        limit_genes=args.limit_genes,
    )

    symbol_map, ensembl_map = build_hgnc_maps(hgnc)

    raw = read_table_flexible(gnomad_file, max_rows=args.max_rows)
    raw.columns = [str(c).strip() for c in raw.columns]

    log("=" * 100)
    log("[RAW GNOMAD TABLE]")
    log(f"[SHAPE] {raw.shape}")
    log(f"[COLUMNS] {list(raw.columns)}")

    mapped = standardize_gnomad_constraint(
        df=raw,
        symbol_map=symbol_map,
        ensembl_map=ensembl_map,
    )

    raw_mapped_path = processed_dir / "feature19_gnomad_constraint_raw_mapped.csv"
    mapped.to_csv(raw_mapped_path, index=False)

    log("=" * 100)
    log("[SAVED RAW MAPPED]")
    log(f"[PATH]  {raw_mapped_path}")
    log(f"[SHAPE] {mapped.shape}")

    features = aggregate_to_gene_level(
        hgnc=hgnc,
        mapped=mapped,
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
        "gnomad_file": str(gnomad_file),
        "protein_coding_only": not args.all_hgnc_genes,
        "limit_genes": args.limit_genes,
        "max_rows": args.max_rows,
        "outputs": {
            "raw_mapped": str(raw_mapped_path),
            "gene_features": str(processed_dir / "feature19_gnomad_full_gene_features.csv"),
            "hgnc_merged": str(processed_dir / "feature19_gnomad_full_hgnc_merged.csv"),
            "summary": str(outdir / "feature19_gnomad_full_summary.txt"),
        },
        "leakage_policy": {
            "included": [
                "gnomAD pLI",
                "gnomAD LOEUF / oe_lof_upper",
                "gnomAD observed/expected LoF",
                "gnomAD observed/expected missense",
                "gnomAD observed/expected synonymous",
                "gnomAD missense Z",
                "gnomAD synonymous Z",
                "gnomAD LoF Z",
                "derived constraint flags",
            ],
            "excluded": [
                "known drug-target labels",
                "clinical target labels",
                "ChEMBL",
                "DrugBank",
                "DGIdb",
                "Open Targets",
                "Pharos",
            ],
        },
    }

    with open(outdir / "feature19_gnomad_full_run_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    write_summary(
        outdir=outdir,
        hgnc=hgnc,
        gnomad_file=gnomad_file,
        raw=raw,
        mapped=mapped,
        features=features,
        merged=merged,
        args=args,
    )

    log("=" * 100)
    log("[DONE]")
    log(f"[RAW MAPPED]    {raw_mapped_path}")
    log(f"[GENE FEATURES] {processed_dir / 'feature19_gnomad_full_gene_features.csv'}")
    log(f"[HGNC MERGED]   {processed_dir / 'feature19_gnomad_full_hgnc_merged.csv'}")
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