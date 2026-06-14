#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature7_GTEx.py

Feature 7: GTEx tissue-expression features for HGNC protein-coding genes.

Default input:
  databases/HGNC/hgnc_complete_set.txt
  feature_databases/GTEx/

Outputs:
  feature7_gtex/processed/feature7_gtex_expression_long.csv
  feature7_gtex/processed/feature7_gtex_gene_features.csv
  feature7_gtex/processed/feature7_gtex_hgnc_merged.csv
  feature7_gtex/feature7_gtex_summary.txt
  feature7_gtex/feature7_gtex_run_metadata.json

Run:
  python Feature7_GTEx.py

Fast test:
  python Feature7_GTEx.py --limit-genes 500

Manual expression file:
  python Feature7_GTEx.py --expression-file feature_databases/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd


DEFAULT_HGNC_PATH = Path("databases") / "HGNC" / "hgnc_complete_set.txt"
DEFAULT_DBDIR = Path("feature_databases") / "GTEx"
DEFAULT_OUTDIR = Path("feature7_gtex")

EXPRESSION_THRESHOLDS = [0.1, 1.0, 5.0, 10.0]

TISSUE_GROUP_PATTERNS = {
    "brain": ["brain", "cortex", "cerebell", "amygdala", "hippocampus", "hypothalamus", "spinal"],
    "blood_immune": ["blood", "spleen", "lymph", "leukocyte", "whole_blood", "ebv"],
    "heart": ["heart", "atrial", "ventric"],
    "muscle": ["muscle", "skeletal"],
    "adipose": ["adipose", "subcutaneous", "visceral"],
    "skin": ["skin", "sun_exposed", "not_sun_exposed"],
    "digestive": ["colon", "esophagus", "stomach", "small_intestine", "liver", "pancreas"],
    "lung": ["lung"],
    "kidney_urinary": ["kidney", "bladder"],
    "reproductive": ["testis", "ovary", "uterus", "vagina", "prostate", "cervix"],
    "endocrine": ["thyroid", "adrenal", "pituitary"],
    "breast": ["breast", "mammary"],
    "nerve": ["nerve", "tibial"],
    "vascular": ["artery", "aorta", "coronary"],
    "cell_line": ["fibroblast", "lymphocyte", "cultured", "ebv"],
}


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
    if s.lower() in {"nan", "none"}:
        return ""
    return s


def normalize_gene_symbol(x: Any) -> str:
    return clean_text(x).upper()


def clean_ensembl_gene_id(x: Any) -> str:
    s = clean_text(x)
    if not s:
        return ""
    return s.split(".")[0]


def safe_float(x: Any) -> float:
    try:
        if x is None or pd.isna(x):
            return np.nan
        return float(x)
    except Exception:
        return np.nan


def open_text_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", errors="ignore")
    return open(path, "rt", errors="ignore")


def sanitize_tissue_name(x: str) -> str:
    s = clean_text(x)
    s = re.sub(r"[^A-Za-z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s)
    return s.strip("_")


def tissue_to_group(tissue: str) -> List[str]:
    low = tissue.lower()
    groups = []

    for group, pats in TISSUE_GROUP_PATTERNS.items():
        if any(p.lower() in low for p in pats):
            groups.append(group)

    return groups if groups else ["other"]


def tau_specificity(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[~np.isnan(arr)]

    if arr.size <= 1:
        return np.nan

    mx = np.max(arr)

    if mx <= 0:
        return 0.0

    return float(np.sum(1.0 - arr / mx) / (arr.size - 1))


def expression_entropy(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[~np.isnan(arr)]
    arr = np.maximum(arr, 0)

    total = np.sum(arr)

    if arr.size == 0 or total <= 0:
        return 0.0

    p = arr / total
    p = p[p > 0]

    return float(-np.sum(p * np.log2(p)))


def expression_gini(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[~np.isnan(arr)]
    arr = np.maximum(arr, 0)

    if arr.size == 0:
        return np.nan

    if np.sum(arr) == 0:
        return 0.0

    arr = np.sort(arr)
    n = arr.size
    idx = np.arange(1, n + 1)

    return float((2 * np.sum(idx * arr) / (n * np.sum(arr))) - ((n + 1) / n))


def find_gtex_expression_file(dbdir: Path) -> Optional[Path]:
    if not dbdir.exists():
        return None

    files = []

    for pattern in ["*.gct", "*.gct.gz", "*.tsv", "*.tsv.gz", "*.txt", "*.txt.gz"]:
        files.extend(list(dbdir.rglob(pattern)))

    if not files:
        return None

    ranked = []

    for p in files:
        name = p.name.lower()
        score = 0

        if "gtex" in name:
            score += 3
        if "gene" in name:
            score += 5
        if "median" in name:
            score += 20
        if "tpm" in name:
            score += 10
        if ".gct" in name:
            score += 8
        if "read_counts" in name or "count" in name:
            score -= 10
        if "junction" in name or "sample" in name or "annotation" in name:
            score -= 20

        ranked.append((score, -p.stat().st_size, p))

    ranked.sort(reverse=True)

    return ranked[0][2]


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
        raise RuntimeError("HGNC file must contain a 'symbol' column.")

    before = len(hgnc)

    if protein_coding_only:
        if "locus_group" in hgnc.columns:
            hgnc = hgnc[hgnc["locus_group"].astype(str).str.lower().eq("protein-coding gene")].copy()
        elif "locus_type" in hgnc.columns:
            hgnc = hgnc[hgnc["locus_type"].astype(str).str.lower().str.contains("protein", na=False)].copy()

    hgnc["gene_symbol"] = hgnc["symbol"].map(normalize_gene_symbol)

    if "ensembl_gene_id" not in hgnc.columns:
        hgnc["ensembl_gene_id"] = ""
    if "entrez_id" not in hgnc.columns:
        hgnc["entrez_id"] = ""
    if "uniprot_ids" not in hgnc.columns:
        hgnc["uniprot_ids"] = ""

    hgnc["ensembl_gene_id_clean"] = hgnc["ensembl_gene_id"].map(clean_ensembl_gene_id)
    hgnc = hgnc.drop_duplicates(subset=["gene_symbol"], keep="first").copy()

    if limit_genes and limit_genes > 0:
        hgnc = hgnc.head(limit_genes).copy()

    log(f"[HGNC ROWS] {before} -> {len(hgnc)}")
    log(f"[PROTEIN CODING ONLY] {protein_coding_only}")
    log(f"[GENES WITH ENSEMBL] {(hgnc['ensembl_gene_id_clean'].astype(str).str.len() > 0).sum()}")

    return hgnc


def build_ensembl_to_gene(hgnc: pd.DataFrame) -> Dict[str, str]:
    mapping = {}

    for _, row in hgnc.iterrows():
        ens = clean_ensembl_gene_id(row.get("ensembl_gene_id_clean", ""))
        gene = row.get("gene_symbol", "")

        if ens and gene and ens not in mapping:
            mapping[ens] = gene

    return mapping


def read_gct_or_tsv_expression(
    expression_file: Path,
    ensembl_to_gene: Dict[str, str],
    processed_dir: Path,
    save_long: bool = True,
) -> pd.DataFrame:
    if not expression_file.exists():
        raise FileNotFoundError(expression_file)

    log("=" * 100)
    log(f"[READ EXPRESSION FILE] {expression_file}")

    with open_text_maybe_gzip(expression_file) as f:
        first = f.readline().strip()
        second = f.readline().strip()

    skiprows = 0

    if first.startswith("#1.") and re.match(r"^\d+\s+\d+", second.replace("\t", " ")):
        skiprows = 2

    log(f"[GCT SKIPROWS] {skiprows}")

    expr = pd.read_csv(
        expression_file,
        sep="\t",
        skiprows=skiprows,
        dtype=str,
        low_memory=False,
        compression="infer",
    )

    expr.columns = [str(c).strip() for c in expr.columns]

    if expr.empty:
        raise RuntimeError("Expression matrix is empty.")

    gene_id_col = None

    for c in ["Name", "name", "gene_id", "Gene ID", "gene", "Ensembl_ID", "ensembl_gene_id"]:
        if c in expr.columns:
            gene_id_col = c
            break

    if gene_id_col is None:
        gene_id_col = expr.columns[0]

    desc_col = None

    for c in ["Description", "description", "gene_name", "Gene Symbol", "gene_symbol"]:
        if c in expr.columns:
            desc_col = c
            break

    metadata_cols = {gene_id_col}

    if desc_col:
        metadata_cols.add(desc_col)

    tissue_cols = []

    for c in expr.columns:
        if c in metadata_cols:
            continue

        test = pd.to_numeric(expr[c].head(200), errors="coerce")

        if test.notna().sum() > 0:
            tissue_cols.append(c)

    if not tissue_cols:
        raise RuntimeError(
            f"No numeric tissue expression columns detected. Columns: {expr.columns[:20].tolist()}"
        )

    log(f"[GENE ID COLUMN] {gene_id_col}")
    log(f"[DESCRIPTION COLUMN] {desc_col if desc_col else 'none'}")
    log(f"[TISSUE COLUMNS] {len(tissue_cols)}")

    hgnc_symbols = set(ensembl_to_gene.values())
    rows = []
    matched_gene_rows = 0

    for _, row in expr.iterrows():
        ens = clean_ensembl_gene_id(row.get(gene_id_col, ""))
        gene = ensembl_to_gene.get(ens, "")

        if not gene and desc_col:
            candidate = normalize_gene_symbol(row.get(desc_col, ""))

            if candidate in hgnc_symbols:
                gene = candidate

        if not gene:
            continue

        matched_gene_rows += 1

        for tissue in tissue_cols:
            val = safe_float(row.get(tissue))

            if np.isnan(val):
                continue

            rows.append(
                {
                    "gene_symbol": gene,
                    "ensembl_gene_id": ens,
                    "tissue": tissue,
                    "tissue_sanitized": sanitize_tissue_name(tissue),
                    "expression_tpm": val,
                }
            )

    long_df = pd.DataFrame(rows)

    if not long_df.empty:
        long_df = long_df.drop_duplicates()

    if save_long:
        long_path = processed_dir / "feature7_gtex_expression_long.csv"
        long_df.to_csv(long_path, index=False)
        log(f"[SAVED LONG] {long_path}")
        log(f"[LONG SHAPE] {long_df.shape}")

    log(f"[MATCHED GTEX GENE ROWS] {matched_gene_rows}")

    return long_df


def build_gene_features(
    hgnc: pd.DataFrame,
    expr_long: pd.DataFrame,
    processed_dir: Path,
) -> pd.DataFrame:
    genes = pd.DataFrame({"gene_symbol": sorted(hgnc["gene_symbol"].unique().tolist())})
    records = []

    by_gene = dict(tuple(expr_long.groupby("gene_symbol"))) if not expr_long.empty else {}

    for gene in genes["gene_symbol"]:
        g = by_gene.get(gene, pd.DataFrame())

        if g.empty:
            records.append(
                {
                    "gene_symbol": gene,
                    "feature7_gtex_has_expression": 0,
                    "feature7_gtex_tissue_count_measured": 0,
                }
            )
            continue

        g = g.copy()
        g["expression_tpm"] = pd.to_numeric(g["expression_tpm"], errors="coerce")
        g = g[g["expression_tpm"].notna()].copy()

        vals = np.maximum(g["expression_tpm"].to_numpy(dtype=float), 0)

        if vals.size == 0:
            records.append(
                {
                    "gene_symbol": gene,
                    "feature7_gtex_has_expression": 0,
                    "feature7_gtex_tissue_count_measured": 0,
                }
            )
            continue

        top_idx = int(np.argmax(vals))
        top_row = g.iloc[top_idx]

        max_expr = float(np.max(vals))
        mean_expr = float(np.mean(vals))
        median_expr = float(np.median(vals))

        rec = {
            "gene_symbol": gene,
            "feature7_gtex_has_expression": 1,
            "feature7_gtex_tissue_count_measured": int(vals.size),

            "feature7_gtex_mean_tpm": mean_expr,
            "feature7_gtex_median_tpm": median_expr,
            "feature7_gtex_max_tpm": max_expr,
            "feature7_gtex_min_tpm": float(np.min(vals)),
            "feature7_gtex_std_tpm": float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0,
            "feature7_gtex_q05_tpm": float(np.quantile(vals, 0.05)),
            "feature7_gtex_q25_tpm": float(np.quantile(vals, 0.25)),
            "feature7_gtex_q75_tpm": float(np.quantile(vals, 0.75)),
            "feature7_gtex_q95_tpm": float(np.quantile(vals, 0.95)),

            "feature7_gtex_mean_log1p_tpm": float(np.mean(np.log1p(vals))),
            "feature7_gtex_median_log1p_tpm": float(np.median(np.log1p(vals))),
            "feature7_gtex_max_log1p_tpm": float(np.max(np.log1p(vals))),

            "feature7_gtex_tau_specificity": tau_specificity(vals),
            "feature7_gtex_expression_entropy": expression_entropy(vals),
            "feature7_gtex_expression_gini": expression_gini(vals),

            "feature7_gtex_top_tissue": str(top_row["tissue"]),
            "feature7_gtex_top_tissue_sanitized": str(top_row["tissue_sanitized"]),
            "feature7_gtex_top_tissue_tpm": max_expr,
            "feature7_gtex_top_to_mean_ratio": float(max_expr / mean_expr) if mean_expr > 0 else np.nan,
            "feature7_gtex_top_to_median_ratio": float(max_expr / median_expr) if median_expr > 0 else np.nan,
        }

        for thr in EXPRESSION_THRESHOLDS:
            label = str(thr).replace(".", "_")
            rec[f"feature7_gtex_tissue_count_tpm_ge_{label}"] = int(np.sum(vals >= thr))
            rec[f"feature7_gtex_tissue_fraction_tpm_ge_{label}"] = float(np.mean(vals >= thr))

        group_to_values: Dict[str, List[float]] = defaultdict(list)

        for _, row in g.iterrows():
            tissue = row["tissue_sanitized"]
            expr = safe_float(row["expression_tpm"])

            for group in tissue_to_group(tissue):
                group_to_values[group].append(expr)

        all_groups = sorted(list(TISSUE_GROUP_PATTERNS.keys()) + ["other"])

        for group in all_groups:
            group_vals = np.asarray(group_to_values.get(group, []), dtype=float)
            group_vals = group_vals[~np.isnan(group_vals)]

            rec[f"feature7_gtex_group_{group}_tissue_count"] = int(group_vals.size)

            if group_vals.size:
                rec[f"feature7_gtex_group_{group}_mean_tpm"] = float(np.mean(group_vals))
                rec[f"feature7_gtex_group_{group}_max_tpm"] = float(np.max(group_vals))
                rec[f"feature7_gtex_group_{group}_has_tpm_ge_1"] = int(np.max(group_vals) >= 1.0)
                rec[f"feature7_gtex_group_{group}_has_tpm_ge_10"] = int(np.max(group_vals) >= 10.0)
            else:
                rec[f"feature7_gtex_group_{group}_mean_tpm"] = np.nan
                rec[f"feature7_gtex_group_{group}_max_tpm"] = np.nan
                rec[f"feature7_gtex_group_{group}_has_tpm_ge_1"] = 0
                rec[f"feature7_gtex_group_{group}_has_tpm_ge_10"] = 0

        records.append(rec)

    features = pd.DataFrame(records)

    for c in features.columns:
        if c.startswith("feature7_gtex_") and (
            "_count" in c
            or "_has_" in c
            or c.startswith("feature7_gtex_has_")
        ):
            if features[c].dtype != object:
                features[c] = features[c].fillna(0)

    for c in ["feature7_gtex_top_tissue", "feature7_gtex_top_tissue_sanitized"]:
        if c in features.columns:
            features[c] = features[c].fillna("")

    outpath = processed_dir / "feature7_gtex_gene_features.csv"
    features.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED GENE FEATURES]")
    log(f"[PATH]  {outpath}")
    log(f"[SHAPE] {features.shape}")

    return features


def merge_with_hgnc(
    hgnc: pd.DataFrame,
    features: pd.DataFrame,
    processed_dir: Path,
) -> pd.DataFrame:
    merged = hgnc.merge(features, on="gene_symbol", how="left")

    feature_cols = [c for c in features.columns if c != "gene_symbol"]

    merged["feature7_gtex_has_any_feature"] = merged[feature_cols].notna().any(axis=1).astype(int)
    merged["feature7_gtex_n_nonmissing_features"] = merged[feature_cols].notna().sum(axis=1)

    outpath = processed_dir / "feature7_gtex_hgnc_merged.csv"
    merged.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED HGNC MERGED]")
    log(f"[PATH]  {outpath}")
    log(f"[SHAPE] {merged.shape}")

    return merged


def write_summary(
    outdir: Path,
    hgnc: pd.DataFrame,
    expr_long: pd.DataFrame,
    features: pd.DataFrame,
    merged: pd.DataFrame,
    expression_file: Path,
    args: argparse.Namespace,
) -> None:
    path = outdir / "feature7_gtex_summary.txt"

    lines = []
    lines.append("Feature 7: GTEx tissue-expression features")
    lines.append("=" * 80)
    lines.append(f"Created at: {now_iso()}")
    lines.append(f"Expression file: {expression_file}")
    lines.append(f"HGNC rows: {hgnc.shape[0]}")
    lines.append(f"Expression long rows: {expr_long.shape[0]}")
    lines.append(f"Genes with expression rows: {expr_long['gene_symbol'].nunique() if not expr_long.empty else 0}")
    lines.append(f"Unique tissues: {expr_long['tissue'].nunique() if not expr_long.empty else 0}")
    lines.append(f"Gene feature rows: {features.shape[0]}")
    lines.append(f"Gene feature columns: {features.shape[1]}")
    lines.append(f"Merged rows: {merged.shape[0]}")
    lines.append(f"Merged columns: {merged.shape[1]}")
    lines.append("")
    lines.append("Coverage:")

    if "feature7_gtex_has_expression" in features.columns:
        lines.append(f"GTEx expression coverage: {features['feature7_gtex_has_expression'].mean():.4f}")

    if "feature7_gtex_tissue_count_measured" in features.columns:
        lines.append(f"Median tissues measured: {features['feature7_gtex_tissue_count_measured'].median():.2f}")

    if "feature7_gtex_mean_tpm" in features.columns:
        lines.append(f"Median mean TPM: {features['feature7_gtex_mean_tpm'].median():.4f}")

    if "feature7_gtex_tau_specificity" in features.columns:
        lines.append(f"Median tau specificity: {features['feature7_gtex_tau_specificity'].median():.4f}")

    lines.append("")
    lines.append("Leakage policy:")
    lines.append("Included: normal-tissue expression breadth, tissue specificity, expression level and tissue-category expression.")
    lines.append("Excluded: disease associations, drug-target labels, tractability labels, ChEMBL/DrugBank/DGIdb labels.")

    path.write_text("\n".join(lines) + "\n")
    log(f"[SUMMARY] {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build GTEx tissue-expression features for HGNC genes.")
    parser.add_argument("--hgnc", default=str(DEFAULT_HGNC_PATH), help="Path to HGNC complete set.")
    parser.add_argument("--dbdir", default=str(DEFAULT_DBDIR), help="GTEx local database directory.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory.")
    parser.add_argument("--expression-file", default="", help="Direct path to GTEx expression matrix.")
    parser.add_argument("--all-hgnc-genes", action="store_true", help="Use all HGNC rows, not only protein-coding.")
    parser.add_argument("--limit-genes", type=int, default=0, help="Debug only: first N HGNC genes.")
    parser.add_argument("--no-long", action="store_true", help="Do not save long expression table.")
    args = parser.parse_args()

    outdir = mkdir(Path(args.outdir))
    processed_dir = mkdir(outdir / "processed")
    dbdir = Path(args.dbdir)

    log("=" * 100)
    log("FEATURE 7: GTEX TISSUE EXPRESSION")
    log("=" * 100)
    log(f"[HGNC]            {args.hgnc}")
    log(f"[DBDIR]           {dbdir.resolve()}")
    log(f"[OUTDIR]          {outdir.resolve()}")
    log(f"[EXPRESSION FILE] {args.expression_file if args.expression_file else 'auto-detect'}")
    log(f"[LIMIT GENES]     {args.limit_genes if args.limit_genes else 'none'}")
    log("=" * 100)

    expression_file = Path(args.expression_file) if args.expression_file else find_gtex_expression_file(dbdir)

    if expression_file is None or not expression_file.exists():
        raise FileNotFoundError(
            "Could not auto-detect a GTEx expression file. Provide one with --expression-file."
        )

    log(f"[USING EXPRESSION FILE] {expression_file}")

    hgnc = load_hgnc(
        Path(args.hgnc),
        protein_coding_only=not args.all_hgnc_genes,
        limit_genes=args.limit_genes,
    )

    ensembl_to_gene = build_ensembl_to_gene(hgnc)

    expr_long = read_gct_or_tsv_expression(
        expression_file=expression_file,
        ensembl_to_gene=ensembl_to_gene,
        processed_dir=processed_dir,
        save_long=not args.no_long,
    )

    features = build_gene_features(
        hgnc=hgnc,
        expr_long=expr_long,
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
        "expression_file": str(expression_file),
        "protein_coding_only": not args.all_hgnc_genes,
        "limit_genes": args.limit_genes,
        "outputs": {
            "expression_long": str(processed_dir / "feature7_gtex_expression_long.csv"),
            "gene_features": str(processed_dir / "feature7_gtex_gene_features.csv"),
            "hgnc_merged": str(processed_dir / "feature7_gtex_hgnc_merged.csv"),
        },
        "leakage_policy": {
            "included": [
                "normal tissue expression",
                "expression breadth",
                "tissue specificity",
                "tissue-category expression features",
            ],
            "excluded": [
                "disease association evidence",
                "known drug-target labels",
                "approved drug annotations",
                "ChEMBL labels",
                "DrugBank labels",
                "DGIdb labels",
                "Open Targets tractability",
                "Open Targets knownDrugs",
            ],
        },
    }

    with open(outdir / "feature7_gtex_run_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    write_summary(outdir, hgnc, expr_long, features, merged, expression_file, args)

    log("=" * 100)
    log("[DONE]")
    log(f"[GENE FEATURES] {processed_dir / 'feature7_gtex_gene_features.csv'}")
    log(f"[HGNC MERGED]   {processed_dir / 'feature7_gtex_hgnc_merged.csv'}")
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