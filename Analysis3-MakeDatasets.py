#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Analysis3-MakeDatasets.py

Generate train-ready gene-level (X, y) dataset pairs from:

    Dataset/Features_All.csv
    Step0_Output/03_HumanGene_DruggabilityLabels.csv

Design
------
TIER 1 — Individual feature blocks:
    F01, F02, ..., F22

TIER 2 — Cumulative feature blocks:
    CUM01 = F01
    CUM02 = F01 + F02
    ...
    CUM22 = F01 + ... + F22

TIER 3 — Thematic feature groups:
    GRP_Structure, GRP_Network, GRP_Expression, GRP_Constraint,
    GRP_Pathway, GRP_Functional, GRP_Omics, GRP_Annotation,
    GRP_Literature, GRP_NoStructure, GRP_Full

Eight targets
-------------
T1 = clinical_target_label
T2 = clinical_investigation_label
T3 = small_molecule_druggable_label
T4 = chemical_tractable_label
T5 = biologic_druggable_label
T6 = dgidb_interaction_supported_label
T7 = potentially_druggable_category_label
T8 = final_any_druggable_label

Outputs
-------
Datasets/
    dataset_catalogue.csv
    dataset_audit.csv
    feature_label_overlap_summary.csv
    run_metadata.json

    Dataset001_F01_T1/
        X_train.csv
        y_train.csv
        meta.json
        feature_columns.txt
        missingness_summary.csv

Notes
-----
This script does not perform model training, scaling, imputation, or train/test
splitting. Missing values are retained in X_train.csv and should be imputed within
the downstream model-training pipeline using training-fold-only parameters.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# =============================================================================
# DEFAULT PATHS
# =============================================================================

DEFAULT_FEATURES_FILE = Path("Dataset") / "Features_All.csv"
DEFAULT_LABELS_FILE = Path("Step0_Output") / "03_HumanGene_DruggabilityLabels.csv"
DEFAULT_OUTDIR = Path("Datasets")


# =============================================================================
# TARGETS
# =============================================================================

TARGETS: Dict[str, Dict[str, str]] = {
    "T1": {
        "column": "clinical_target_label",
        "name": "Clinical Target",
        "description": "Approved or clinically established drug target; strictest target definition.",
        "leakage_note": "Derived from ChEMBL and Pharos/TCRD clinical target evidence. Direct source columns must not be used as predictors.",
    },
    "T2": {
        "column": "clinical_investigation_label",
        "name": "Clinical Investigation Target",
        "description": "Gene with clinical-investigation evidence, including phase 1--3 clinical trial support.",
        "leakage_note": "Derived from clinical-indication evidence. Direct source columns must not be used as predictors.",
    },
    "T3": {
        "column": "small_molecule_druggable_label",
        "name": "Small-Molecule Target",
        "description": "Gene with small-molecule tractability evidence.",
        "leakage_note": "Derived from ChEMBL and Open Targets small-molecule tractability evidence. Direct source columns must not be used as predictors.",
    },
    "T4": {
        "column": "chemical_tractable_label",
        "name": "Chemical Tractability Target",
        "description": "Gene with broader chemical tractability evidence.",
        "leakage_note": "Derived from chemical tractability evidence. Direct source columns must not be used as predictors.",
    },
    "T5": {
        "column": "biologic_druggable_label",
        "name": "Biologic/Modality Target",
        "description": "Gene with biologic, antibody, targeted-degradation or other modality-based tractability evidence.",
        "leakage_note": "Derived from biologic/modality tractability evidence. Direct source columns must not be used as predictors.",
    },
    "T6": {
        "column": "dgidb_interaction_supported_label",
        "name": "Drug--Gene Interaction Target",
        "description": "Gene supported by DGIdb drug--gene interaction evidence.",
        "leakage_note": "Derived from DGIdb. DGIdb-derived columns must not be used as predictors.",
    },
    "T7": {
        "column": "potentially_druggable_category_label",
        "name": "Potentially Druggable Family Target",
        "description": "Gene assigned to a potentially druggable protein-family or category-based annotation.",
        "leakage_note": "Derived from protein-family/category annotations and should be interpreted as broad potential druggability evidence.",
    },
    "T8": {
        "column": "final_any_druggable_label",
        "name": "Broad Druggability Target",
        "description": "Broad evidence-union label; positive if a gene satisfied any druggability definition.",
        "leakage_note": "Union of all target definitions. Use cautiously because this is the broadest and most permissive label.",
    },
}


# =============================================================================
# FEATURE GROUP DEFINITIONS
# =============================================================================

INDIVIDUAL_GROUPS: Dict[str, List[int]] = {
    f"F{n:02d}": [n] for n in range(1, 23)
}

CUMULATIVE_GROUPS: Dict[str, List[int]] = {
    f"CUM{n:02d}": list(range(1, n + 1)) for n in range(1, 23)
}

THEMATIC_GROUPS: Dict[str, List[int]] = {
    "GRP_Structure": [4, 10],
    "GRP_Network": [2, 15],
    "GRP_Expression": [7, 16],
    "GRP_Constraint": [9, 19],
    "GRP_Pathway": [3, 14],
    "GRP_Functional": [5, 6, 11, 12],
    "GRP_Omics": [1, 18],
    "GRP_Annotation": [8, 13, 20, 21, 22],
    "GRP_Literature": [17],
    "GRP_NoStructure": [n for n in range(1, 23) if n not in (4, 10)],
    "GRP_Full": list(range(1, 23)),
}

ALL_GROUPS: Dict[str, List[int]] = {}
ALL_GROUPS.update(INDIVIDUAL_GROUPS)
ALL_GROUPS.update(CUMULATIVE_GROUPS)
ALL_GROUPS.update(THEMATIC_GROUPS)


# =============================================================================
# NON-FEATURE COLUMNS
# =============================================================================

NON_FEATURE_PATTERNS = [
    r"^gene_symbol$",
    r"^approved_symbol$",
    r"^symbol$",
    r"^hgnc_id$",
    r"^entrez_id$",
    r"^ensembl_gene_id$",
    r"_has_any_feature$",
    r"_n_nonmissing_features$",
    r"_has_any$",
    r"_ids$",
    r"_names$",
    r"_text$",
    r"_raw$",
    r"_audit_only$",
    r"_file$",
    r"_path$",
    r"_url$",
    r"_version$",
    r"_source$",
    r"locus_group$",
    r"locus_type$",
    r"gene_family",
    r"alias_symbol",
    r"prev_symbol",
    r"^name$",
]


# =============================================================================
# HELPERS
# =============================================================================

def log(msg: str) -> None:
    print(msg, flush=True)


def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def now_iso() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def normalize_symbol(x: Any) -> str:
    if x is None:
        return ""
    try:
        if isinstance(x, float) and np.isnan(x):
            return ""
    except Exception:
        pass
    s = str(x).strip().upper()
    s = re.sub(r"\s+", "", s)
    if s.lower() in {"nan", "none", "null", "na"}:
        return ""
    return s


def is_non_feature_col(col: str) -> bool:
    c = col.lower()
    for pat in NON_FEATURE_PATTERNS:
        if re.search(pat.lower(), c):
            return True
    return False


def get_tier(group_name: str) -> str:
    if group_name.startswith("F") and not group_name.startswith("Full"):
        return "Individual"
    if group_name.startswith("CUM"):
        return "Cumulative"
    return "Thematic"


def sort_group_key(group_name: str) -> Tuple[int, str]:
    if group_name.startswith("F") and not group_name.startswith("Full"):
        return (0, group_name)
    if group_name.startswith("CUM"):
        return (1, group_name)
    return (2, group_name)


def extract_feature_number(col: str) -> Optional[int]:
    m = re.match(r"^Feature(\d+)_", col)
    if not m:
        return None
    return int(m.group(1))


# =============================================================================
# LOAD DATA
# =============================================================================

def load_features(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Features file not found: {path}\n"
            "Run Merge_Features.py first to generate Dataset/Features_All.csv"
        )

    log(f"[LOAD FEATURES] {path}")
    df = pd.read_csv(path, low_memory=False)

    if "gene_symbol" not in df.columns:
        raise RuntimeError("Features file must contain a 'gene_symbol' column.")

    before = len(df)
    df["gene_symbol"] = df["gene_symbol"].apply(normalize_symbol)
    df = df[df["gene_symbol"] != ""].copy()

    duplicated = int(df["gene_symbol"].duplicated().sum())
    if duplicated > 0:
        log(f"  [WARNING] duplicate gene symbols in features: {duplicated}; keeping first occurrence")

    df = df.drop_duplicates("gene_symbol", keep="first").reset_index(drop=True)

    log(f"  rows before cleaning: {before:,}")
    log(f"  genes after cleaning: {df.shape[0]:,}")
    log(f"  feature columns:       {df.shape[1] - 1:,}")

    return df


def load_labels(path: Path, target_cols: List[str]) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Labels file not found: {path}\n"
            "Run the target-label construction pipeline first."
        )

    log(f"[LOAD LABELS]   {path}")
    df = pd.read_csv(path, low_memory=False)

    symbol_col = None
    for c in ["gene_symbol", "symbol", "Gene_Symbol", "approved_symbol"]:
        if c in df.columns:
            symbol_col = c
            break

    if symbol_col is None:
        raise RuntimeError(
            "Labels file must contain one of: gene_symbol, symbol, Gene_Symbol, approved_symbol"
        )

    before = len(df)
    df["gene_symbol"] = df[symbol_col].apply(normalize_symbol)
    df = df[df["gene_symbol"] != ""].copy()

    duplicated = int(df["gene_symbol"].duplicated().sum())
    if duplicated > 0:
        log(f"  [WARNING] duplicate gene symbols in labels: {duplicated}; keeping first occurrence")

    df = df.drop_duplicates("gene_symbol", keep="first").reset_index(drop=True)

    keep = ["gene_symbol"] + [c for c in target_cols if c in df.columns]
    df = df[keep].copy()

    for c in target_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)

    missing = [c for c in target_cols if c not in df.columns]
    if missing:
        log(f"  [WARNING] Target columns not found in labels file: {missing}")

    log(f"  rows before cleaning: {before:,}")
    log(f"  genes after cleaning: {df.shape[0]:,}")
    log(f"  target columns loaded: {len(keep) - 1}")

    return df


# =============================================================================
# FEATURE SELECTION
# =============================================================================

def get_all_feature_columns(features_df: pd.DataFrame) -> List[str]:
    cols = []
    for col in features_df.columns:
        if col == "gene_symbol":
            continue
        if is_non_feature_col(col):
            continue
        if extract_feature_number(col) is None:
            continue
        cols.append(col)
    return cols


def get_feature_columns(
    all_feature_cols: List[str],
    feature_numbers: List[int],
) -> List[str]:
    selected = []
    feature_numbers_set = set(feature_numbers)

    for col in all_feature_cols:
        n = extract_feature_number(col)
        if n in feature_numbers_set:
            selected.append(col)

    return selected


# =============================================================================
# AUDITS
# =============================================================================

def save_feature_label_overlap_summary(
    features_df: pd.DataFrame,
    labels_df: pd.DataFrame,
    outdir: Path,
) -> pd.DataFrame:
    feature_genes = set(features_df["gene_symbol"].dropna())
    label_genes = set(labels_df["gene_symbol"].dropna())
    overlap = feature_genes & label_genes

    rows = [
        {
            "n_genes_features": len(feature_genes),
            "n_genes_labels": len(label_genes),
            "n_genes_overlap": len(overlap),
            "n_features_only": len(feature_genes - label_genes),
            "n_labels_only": len(label_genes - feature_genes),
            "pct_features_retained_after_label_join": round(100 * len(overlap) / len(feature_genes), 3) if feature_genes else 0.0,
            "pct_labels_retained_after_feature_join": round(100 * len(overlap) / len(label_genes), 3) if label_genes else 0.0,
        }
    ]

    out = pd.DataFrame(rows)
    out.to_csv(outdir / "feature_label_overlap_summary.csv", index=False)

    pd.Series(sorted(feature_genes - label_genes), name="gene_symbol").to_csv(
        outdir / "genes_in_features_not_labels.csv",
        index=False,
    )
    pd.Series(sorted(label_genes - feature_genes), name="gene_symbol").to_csv(
        outdir / "genes_in_labels_not_features.csv",
        index=False,
    )

    return out


def build_missingness_summary(X: pd.DataFrame) -> pd.DataFrame:
    rows = []
    n = len(X)

    for col in X.columns:
        n_missing = int(X[col].isna().sum())
        rows.append(
            {
                "feature_column": col,
                "n_missing": n_missing,
                "pct_missing": round(100 * n_missing / n, 4) if n > 0 else 0.0,
                "n_nonmissing": int(X[col].notna().sum()),
            }
        )

    return pd.DataFrame(rows).sort_values(
        ["pct_missing", "feature_column"],
        ascending=[False, True],
    )


# =============================================================================
# BUILD ONE DATASET
# =============================================================================

def build_dataset(
    features_df: pd.DataFrame,
    labels_df: pd.DataFrame,
    feature_cols: List[str],
    target_col: str,
    group_name: str,
    target_key: str,
    dataset_id: int,
    outdir: Path,
    min_positives: int,
    dry_run: bool,
) -> Optional[Dict[str, Any]]:

    if not feature_cols:
        return None

    if target_col not in labels_df.columns:
        return None

    # Gene-level inner join.
    df = (
        labels_df[["gene_symbol", target_col]]
        .merge(
            features_df[["gene_symbol"] + feature_cols],
            on="gene_symbol",
            how="inner",
        )
    )

    df = df[df[target_col].notna()].copy()
    df[target_col] = df[target_col].astype(int)

    n_genes = int(len(df))
    if n_genes == 0:
        return None

    n_positives = int(df[target_col].sum())
    n_negatives = int((df[target_col] == 0).sum())
    prevalence = round(float(n_positives / n_genes), 6) if n_genes else 0.0

    if n_positives < min_positives:
        log(
            f"  [SKIP] Dataset{dataset_id:03d}  {group_name:<30s} {target_key}"
            f" — only {n_positives} positives (min={min_positives})"
        )
        return None

    X_all = df[feature_cols].copy()

    # Drop columns entirely missing in this dataset.
    all_nan_cols = [c for c in X_all.columns if X_all[c].isna().all()]
    X = X_all.drop(columns=all_nan_cols)

    if X.shape[1] == 0:
        log(
            f"  [SKIP] Dataset{dataset_id:03d}  {group_name:<30s} {target_key}"
            " — no usable feature columns after dropping all-NaN columns"
        )
        return None

    feature_cols_used = list(X.columns)
    feat_nums_used = sorted(
        {
            extract_feature_number(c)
            for c in feature_cols_used
            if extract_feature_number(c) is not None
        }
    )

    n_values = int(X.shape[0] * X.shape[1])
    n_missing_values = int(X.isna().sum().sum())
    pct_missing_values = round(100 * n_missing_values / n_values, 6) if n_values else 0.0

    dataset_name = f"Dataset{dataset_id:03d}_{group_name}_{target_key}"
    dataset_dir = mkdir(outdir / dataset_name)

    missingness_df = build_missingness_summary(X)

    meta = {
        "dataset_id": dataset_id,
        "dataset_name": dataset_name,
        "tier": get_tier(group_name),
        "group_name": group_name,
        "requested_feature_numbers": list(ALL_GROUPS[group_name]),
        "feature_numbers_used": feat_nums_used,
        "n_feature_blocks_used": len(feat_nums_used),
        "target_key": target_key,
        "target_name": TARGETS[target_key]["name"],
        "target_column": target_col,
        "target_description": TARGETS[target_key]["description"],
        "leakage_note": TARGETS[target_key]["leakage_note"],
        "n_genes": n_genes,
        "n_positives": n_positives,
        "n_negatives": n_negatives,
        "prevalence": prevalence,
        "n_features_requested": len(feature_cols),
        "n_features_used": len(feature_cols_used),
        "n_all_nan_features_dropped": len(all_nan_cols),
        "all_nan_features_dropped": all_nan_cols,
        "n_missing_values": n_missing_values,
        "pct_missing_values": pct_missing_values,
        "created_at": now_iso(),
        "note": "Missing values are retained and should be imputed within downstream training folds only.",
    }

    if not dry_run:
        X_out = pd.concat(
            [
                df[["gene_symbol"]].reset_index(drop=True),
                X.reset_index(drop=True),
            ],
            axis=1,
        )
        X_out.to_csv(dataset_dir / "X_train.csv", index=False)

        df[["gene_symbol", target_col]].to_csv(
            dataset_dir / "y_train.csv",
            index=False,
        )

        with open(dataset_dir / "meta.json", "w") as f:
            json.dump(meta, f, indent=2)

        with open(dataset_dir / "feature_columns.txt", "w") as f:
            for col in feature_cols_used:
                f.write(f"{col}\n")

        missingness_df.to_csv(dataset_dir / "missingness_summary.csv", index=False)

    log(
        f"  [OK]  Dataset{dataset_id:03d}  {group_name:<30s} {target_key}"
        f"  genes={n_genes:,} pos={n_positives:,}"
        f"  prev={prevalence:.3f} features={len(feature_cols_used):,}"
        f"  missing={pct_missing_values:.2f}%"
    )

    return meta


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate all gene-level (X, y) dataset pairs for druggability modelling."
    )
    parser.add_argument(
        "--features-file",
        default=str(DEFAULT_FEATURES_FILE),
        help=f"Path to merged features CSV. Default: {DEFAULT_FEATURES_FILE}",
    )
    parser.add_argument(
        "--labels-file",
        default=str(DEFAULT_LABELS_FILE),
        help=f"Path to druggability labels CSV. Default: {DEFAULT_LABELS_FILE}",
    )
    parser.add_argument(
        "--outdir",
        default=str(DEFAULT_OUTDIR),
        help=f"Output directory. Default: {DEFAULT_OUTDIR}",
    )
    parser.add_argument(
        "--targets",
        nargs="*",
        default=list(TARGETS.keys()),
        help="Subset of targets to generate, e.g. --targets T1 T3 T8",
    )
    parser.add_argument(
        "--groups",
        nargs="*",
        default=None,
        help="Subset of groups to generate, e.g. --groups F01 CUM10 GRP_Full",
    )
    parser.add_argument(
        "--min-positives",
        type=int,
        default=30,
        help="Skip datasets with fewer than this many positive examples. Default: 30",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be created without writing dataset files.",
    )

    args = parser.parse_args()

    outdir = mkdir(Path(args.outdir))

    log("=" * 100)
    log("ANALYSIS3 — MAKE DATASETS")
    log(f"Started:        {now_iso()}")
    log(f"Features file:  {args.features_file}")
    log(f"Labels file:    {args.labels_file}")
    log(f"Output dir:     {outdir.resolve()}")
    log(f"Targets:        {args.targets}")
    log(f"Min positives:  {args.min_positives}")
    log(f"Dry run:        {args.dry_run}")
    log("=" * 100)

    # Validate target keys.
    invalid_targets = [t for t in args.targets if t not in TARGETS]
    if invalid_targets:
        raise ValueError(f"Invalid target keys: {invalid_targets}. Valid keys: {list(TARGETS)}")

    target_cols = [TARGETS[k]["column"] for k in args.targets]

    # Load inputs.
    features_df = load_features(Path(args.features_file))
    labels_df = load_labels(Path(args.labels_file), target_cols)

    # Save feature-label overlap audit.
    overlap_df = save_feature_label_overlap_summary(features_df, labels_df, outdir)
    log("\n[FEATURE-LABEL OVERLAP]")
    print(overlap_df.to_string(index=False))

    all_feature_cols = get_all_feature_columns(features_df)

    log(f"\n[FEATURE COLUMNS AVAILABLE]  {len(all_feature_cols):,}")
    log(f"[GENES IN FEATURES]          {len(features_df):,}")
    log(f"[GENES IN LABELS]            {len(labels_df):,}")
    log(f"[TARGETS REQUESTED]          {args.targets}")

    # Select groups.
    groups_to_run = args.groups if args.groups else list(ALL_GROUPS.keys())
    invalid_groups = [g for g in groups_to_run if g not in ALL_GROUPS]
    if invalid_groups:
        raise ValueError(f"Invalid group names: {invalid_groups}")

    groups_to_run = sorted(groups_to_run, key=sort_group_key)

    n_individual = sum(1 for g in groups_to_run if get_tier(g) == "Individual")
    n_cumulative = sum(1 for g in groups_to_run if get_tier(g) == "Cumulative")
    n_thematic = sum(1 for g in groups_to_run if get_tier(g) == "Thematic")

    log(f"\n[GROUPS — INDIVIDUAL]        {n_individual}")
    log(f"[GROUPS — CUMULATIVE]        {n_cumulative}")
    log(f"[GROUPS — THEMATIC]          {n_thematic}")
    log(f"[GROUPS TOTAL]               {len(groups_to_run)}")
    log(f"[TARGETS]                    {len(args.targets)}")
    log(f"[MAX DATASETS POSSIBLE]      {len(groups_to_run) * len(args.targets)}")

    # Generate datasets.
    log("\n" + "=" * 100)
    log("GENERATING DATASETS")
    log("=" * 100)

    catalogue: List[Dict[str, Any]] = []
    dataset_id = 1

    for group_name in groups_to_run:
        feature_numbers = ALL_GROUPS[group_name]
        feature_cols = get_feature_columns(all_feature_cols, feature_numbers)

        if not feature_cols:
            log(f"\n  [SKIP GROUP] {group_name} — no matching columns in Features_All.csv")
            continue

        log(f"\n  {group_name:<30s} feature_blocks={feature_numbers} cols={len(feature_cols)}")

        for target_key in args.targets:
            target_col = TARGETS[target_key]["column"]

            meta = build_dataset(
                features_df=features_df,
                labels_df=labels_df,
                feature_cols=feature_cols,
                target_col=target_col,
                group_name=group_name,
                target_key=target_key,
                dataset_id=dataset_id,
                outdir=outdir,
                min_positives=args.min_positives,
                dry_run=args.dry_run,
            )

            if meta is not None:
                catalogue.append(meta)
                dataset_id += 1

    # Save summary files.
    log("\n" + "=" * 100)
    log("SUMMARY")
    log("=" * 100)

    cat_df = pd.DataFrame(catalogue)

    run_meta = {
        "created_at": now_iso(),
        "features_file": str(args.features_file),
        "labels_file": str(args.labels_file),
        "outdir": str(outdir.resolve()),
        "targets_requested": args.targets,
        "groups_requested": groups_to_run,
        "min_positives": args.min_positives,
        "dry_run": args.dry_run,
        "n_feature_columns_available": len(all_feature_cols),
        "n_genes_features": int(len(features_df)),
        "n_genes_labels": int(len(labels_df)),
        "n_datasets_generated": int(len(cat_df)),
    }

    if not args.dry_run:
        with open(outdir / "run_metadata.json", "w") as f:
            json.dump(run_meta, f, indent=2)

    if cat_df.empty:
        log("[WARNING] No datasets were generated. Check paths, targets, groups and --min-positives.")
        return

    if not args.dry_run:
        cat_path = outdir / "dataset_catalogue.csv"
        audit_path = outdir / "dataset_audit.csv"

        cat_df.to_csv(cat_path, index=False)

        audit_cols = [
            "dataset_id",
            "dataset_name",
            "tier",
            "group_name",
            "target_key",
            "target_name",
            "target_column",
            "n_genes",
            "n_positives",
            "n_negatives",
            "prevalence",
            "n_features_requested",
            "n_features_used",
            "n_all_nan_features_dropped",
            "n_missing_values",
            "pct_missing_values",
            "feature_numbers_used",
        ]
        cat_df[[c for c in audit_cols if c in cat_df.columns]].to_csv(audit_path, index=False)

        log(f"[SAVED CATALOGUE] {cat_path.resolve()}")
        log(f"[SAVED AUDIT]     {audit_path.resolve()}")

    log(f"[DATASETS GENERATED]   {len(cat_df)}")
    log(f"[TARGETS COVERED]      {cat_df['target_key'].nunique()}")
    log(f"[GROUPS COVERED]       {cat_df['group_name'].nunique()}")

    log("\nDatasets per tier:")
    print(cat_df.groupby("tier")["dataset_id"].count().rename("n_datasets").to_string())

    log("\nDatasets per target:")
    print(
        cat_df.groupby(["target_key", "target_column"])[
            ["n_genes", "n_positives", "prevalence", "n_features_used", "pct_missing_values"]
        ].mean().round(4).to_string()
    )

    log("\n[DONE]")
    log(f"Finished: {now_iso()}")
    log("=" * 100)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        log("\n[STOPPED] Interrupted by user.")
        sys.exit(130)
    except Exception as exc:
        log(f"\n[ERROR] {exc}")
        raise