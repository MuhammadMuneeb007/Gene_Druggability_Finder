#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Merge_Features.py

Merge all 22 feature files into one master matrix.

Output:
    Dataset/Features_All.csv

Every feature column is prefixed with its Feature number:
    Feature1_*
    Feature2_*
    ...
    Feature22_*

This version fixes Feature10 by accepting:
    feature10_*
    Feature10_*

and by checking both:
    feature10_pocket_geometry_gene_features.csv
    feature10_pocket_geometry_hgnc_merged.csv

Usage:
    python Merge_Features.py
    python Merge_Features.py --outdir Dataset
    python Merge_Features.py --verbose
"""

from __future__ import annotations

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, List, Optional, Tuple

import numpy as np
import pandas as pd


# =============================================================================
# SETTINGS
# =============================================================================

DEFAULT_OUTDIR = Path("Dataset")

# (feature_number, clean_label, list_of_possible_csv_paths, join_key_in_file)
FEATURE_FILES: List[Tuple[int, str, List[Path], str]] = [
    (
        1,
        "DepMap",
        [Path("feature1.dmapp.database") / "processed" / "feature1_depmap_gene_features.csv"],
        "gene_symbol",
    ),
    (
        2,
        "String",
        [Path("feature2.string.database") / "processed" / "feature2_string_gene_features.csv"],
        "gene_symbol",
    ),
    (
        3,
        "Reactome",
        [Path("feature3_pathway") / "processed" / "feature3_reactome_gene_features.csv"],
        "gene_symbol",
    ),
    (
        4,
        "Structure",
        [Path("feature4_structure") / "processed" / "feature4_structure_gene_features.csv"],
        "gene_symbol",
    ),
    (
        5,
        "InterProPfam",
        [Path("feature5_interpro_pfam") / "processed" / "feature5_interpro_pfam_gene_features.csv"],
        "gene_symbol",
    ),
    (
        6,
        "UniProt",
        [Path("feature6_uniprot") / "processed" / "feature6_uniprot_gene_features.csv"],
        "gene_symbol",
    ),
    (
        7,
        "GTEx",
        [Path("feature7_gtex") / "processed" / "feature7_gtex_gene_features.csv"],
        "gene_symbol",
    ),
    (
        8,
        "Ensembl",
        [Path("feature8_ensembl") / "processed" / "feature8_ensembl_gene_features.csv"],
        "gene_symbol",
    ),
    (
        9,
        "ConstraintDisorder",
        [Path("feature9_disorder_constraint") / "processed" / "feature9_disorder_constraint_gene_features.csv"],
        "gene_symbol",
    ),
    (
        10,
        "FPocket",
        [
            Path("feature10_pocket_geometry") / "processed" / "feature10_pocket_geometry_gene_features.csv",
            Path("feature10_pocket_geometry") / "processed" / "feature10_pocket_geometry_hgnc_merged.csv",
        ],
        "gene_symbol",
    ),
    (
        11,
        "ProteinSequence",
        [Path("feature11_protein_sequence") / "processed" / "feature11_protein_sequence_gene_features.csv"],
        "gene_symbol",
    ),
    (
        12,
        "Embeddings",
        [Path("feature12_protein_embeddings") / "processed" / "feature12_embedding_gene_pca_features.csv"],
        "gene_symbol",
    ),
    (
        13,
        "GWAS",
        [Path("feature13_gwas_catalog") / "processed" / "feature13_gwas_catalog_gene_features.csv"],
        "gene_symbol",
    ),
    (
        14,
        "GO",
        [Path("feature14_gene_ontology") / "processed" / "feature14_go_gene_features.csv"],
        "gene_symbol",
    ),
    (
        15,
        "BioGRID",
        [Path("feature15_biogrid") / "processed" / "feature15_biogrid_gene_features.csv"],
        "gene_symbol",
    ),
    (
        16,
        "HPA",
        [Path("feature16_hpa") / "processed" / "feature16_hpa_gene_features.csv"],
        "gene_symbol",
    ),
    (
        17,
        "CTD",
        [Path("feature17_ctd") / "processed" / "feature17_ctd_gene_features.csv"],
        "approved_symbol",
    ),
    (
        18,
        "MGI",
        [Path("feature18_mgi") / "processed" / "feature18_mgi_gene_features.csv"],
        "gene_symbol",
    ),
    (
        19,
        "gnomADFull",
        [Path("feature19_gnomad_full") / "processed" / "feature19_gnomad_full_gene_features.csv"],
        "gene_symbol",
    ),
    (
        20,
        "CORUM",
        [Path("feature20_corum") / "processed" / "feature20_corum_gene_features.csv"],
        "approved_symbol",
    ),
    (
        21,
        "PhosphoSitePlus",
        [Path("feature21_phosphositeplus") / "processed" / "feature21_psp_gene_features.csv"],
        "approved_symbol",
    ),
    (
        22,
        "Paralogues",
        [Path("feature22_paralogues") / "processed" / "feature22_paralogues_gene_features.csv"],
        "approved_symbol",
    ),
]


# Columns that are identifiers/metadata, not model features.
# Important: do NOT remove normal feature columns such as:
#     feature10_pocket_has_structure
#     feature10_pocket_structure_count
NON_FEATURE_PATTERNS = [
    r"^gene_symbol$",
    r"^approved_symbol$",
    r"^symbol$",
    r"^hgnc_id$",
    r"^entrez_id$",
    r"^ensembl_gene_id$",
    r"^uniprot_ids$",
    r"^uniprot_accession$",
    r"^pdb_id$",
    r"^structure_source$",
    r"^structure_id$",
    r"^alphafold_id$",
    r"^protein_id$",
    r"^transcript_id$",
    r"^mouse_gene_id$",
    r"^mgi_id$",
    r"^name$",
    r"^description$",
    r"^location$",
    r"^locus_group$",
    r"^locus_type$",
    r"^alias_symbol$",
    r"^prev_symbol$",
    r"_raw$",
    r"_audit_only$",
    r"_names$",
    r"_text$",
    r"_file$",
    r"_path$",
    r"_url$",
    r"_version$",
]


# =============================================================================
# HELPERS
# =============================================================================

def log(msg: str) -> None:
    print(msg, flush=True)


def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


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


def choose_existing_path(paths: List[Path]) -> Optional[Path]:
    for path in paths:
        if path.exists():
            return path
    return None


def strip_existing_feature_prefix(col: str, feat_num: int) -> str:
    """
    Remove existing feature prefix so we do not double-prefix too badly.

    Examples:
        feature10_pocket_has_structure -> pocket_has_structure
        Feature10_pocket_has_structure -> pocket_has_structure
        feature1_depmap_crispr_mean -> depmap_crispr_mean
    """
    col2 = re.sub(
        rf"^(?:Feature|feature){feat_num}_",
        "",
        col,
        count=1,
        flags=re.IGNORECASE,
    )
    return col2


def standardize_feature_column_name(col: str, feat_num: int, feat_label: str) -> str:
    """
    Final feature naming style:
        Feature10_FPocket_pocket_has_structure
    """
    stripped = strip_existing_feature_prefix(col, feat_num)
    stripped = stripped.strip("_")
    stripped = re.sub(r"[^A-Za-z0-9_]+", "_", stripped)
    stripped = re.sub(r"_+", "_", stripped)
    return f"Feature{feat_num}_{feat_label}_{stripped}"


def collapse_duplicate_genes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse duplicate gene rows.

    Numeric columns are averaged. Since this master file is numeric-only after
    filtering, this is safe.
    """
    if not df["gene_symbol"].duplicated().any():
        return df

    feature_cols = [c for c in df.columns if c != "gene_symbol"]
    agg_dict = {c: "mean" for c in feature_cols}
    return df.groupby("gene_symbol", as_index=False).agg(agg_dict)


# =============================================================================
# LOAD ONE FEATURE FILE
# =============================================================================

def load_feature_file(
    feat_num: int,
    feat_label: str,
    possible_paths: List[Path],
    join_key: str,
    verbose: bool,
) -> Optional[pd.DataFrame]:

    path = choose_existing_path(possible_paths)

    if path is None:
        log(f"  [MISSING]  Feature{feat_num} ({feat_label})")
        for p in possible_paths:
            log(f"             checked: {p}")
        return None

    log(f"  [LOADING]  Feature{feat_num} ({feat_label}): {path}")

    df = pd.read_csv(path, low_memory=False)
    df.columns = [str(c).strip() for c in df.columns]

    # Normalize the join key to gene_symbol.
    if join_key not in df.columns:
        fallback_found = None
        for fallback in ["gene_symbol", "approved_symbol", "symbol", "Gene_Symbol", "gene", "Gene"]:
            if fallback in df.columns:
                fallback_found = fallback
                break

        if fallback_found is None:
            log(f"  [WARNING]  No join key column found in {path.name} — skipping")
            return None

        join_key = fallback_found

    df["gene_symbol"] = df[join_key].apply(normalize_symbol)
    df = df[df["gene_symbol"] != ""].copy()

    # Select numeric feature columns.
    cols_to_keep = ["gene_symbol"]
    feature_cols = []

    for col in df.columns:
        if col == "gene_symbol":
            continue
        if col == join_key:
            continue
        if is_non_feature_col(col):
            continue

        converted = pd.to_numeric(df[col], errors="coerce")

        # Keep a column only if it contains at least one numeric value.
        if converted.notna().sum() == 0:
            continue

        df[col] = converted
        cols_to_keep.append(col)
        feature_cols.append(col)

    df = df[cols_to_keep].copy()

    if not feature_cols:
        log(f"  [WARNING]  No numeric feature columns in {path.name} — skipping")
        return None

    # Rename all feature columns to standard FeatureN_Label_* names.
    rename_map = {}
    for col in feature_cols:
        rename_map[col] = standardize_feature_column_name(col, feat_num, feat_label)

    df = df.rename(columns=rename_map)

    # Remove duplicate columns if any.
    df = df.loc[:, ~df.columns.duplicated()].copy()

    # Collapse duplicate genes.
    df = collapse_duplicate_genes(df)

    renamed_feature_cols = [c for c in df.columns if c != "gene_symbol"]

    log(f"  [LOADED]   Feature{feat_num}: {df.shape[0]:,} genes, {len(renamed_feature_cols):,} features")

    if verbose:
        log("             first columns:")
        for c in renamed_feature_cols[:10]:
            log(f"               - {c}")
        if len(renamed_feature_cols) > 10:
            log(f"               ... and {len(renamed_feature_cols) - 10} more")

    return df


# =============================================================================
# MAIN MERGE
# =============================================================================

def build_gene_universe(loaded: list[pd.DataFrame]) -> pd.DataFrame:
    """
    Union of all gene symbols across all loaded feature files.
    """
    all_genes = set()

    for df in loaded:
        all_genes.update(df["gene_symbol"].dropna().unique())

    all_genes.discard("")
    genes = sorted(all_genes)

    log(f"\n[GENE UNIVERSE]  {len(genes):,} unique genes across all loaded feature files")

    return pd.DataFrame({"gene_symbol": genes})


def merge_all(
    feature_files: List[Tuple[int, str, List[Path], str]],
    verbose: bool,
) -> pd.DataFrame:

    log("=" * 80)
    log("LOADING FEATURE FILES")
    log("=" * 80)

    loaded = []
    loaded_blocks = []
    missing_blocks = []

    for feat_num, feat_label, possible_paths, join_key in feature_files:
        df = load_feature_file(feat_num, feat_label, possible_paths, join_key, verbose)

        if df is not None:
            loaded.append(df)
            loaded_blocks.append(feat_num)
        else:
            missing_blocks.append(feat_num)

    log("")
    log(f"[FILES LOADED]  {len(loaded_blocks)}")
    log(f"[FILES MISSING] {len(missing_blocks)}")

    if missing_blocks:
        log("[MISSING BLOCKS] " + ", ".join(f"F{x:02d}" for x in missing_blocks))

    if not loaded:
        raise RuntimeError("No feature files could be loaded. Check paths above.")

    master = build_gene_universe(loaded)

    log("\n" + "=" * 80)
    log("MERGING ALL FEATURES")
    log("=" * 80)

    for df in loaded:
        before_cols = master.shape[1]
        master = master.merge(df, on="gene_symbol", how="left")
        after_cols = master.shape[1]
        log(f"  merged {after_cols - before_cols:>5} columns | master shape: {master.shape}")

    master = master.drop_duplicates(subset=["gene_symbol"], keep="first").copy()
    master = master.sort_values("gene_symbol").reset_index(drop=True)

    feature_cols = [c for c in master.columns if c != "gene_symbol"]

    log(f"\n[FINAL MATRIX]  {master.shape[0]:,} genes x {len(feature_cols):,} features")

    return master


# =============================================================================
# COVERAGE REPORT
# =============================================================================

def print_coverage_report(master: pd.DataFrame) -> None:
    log("\n" + "=" * 80)
    log("COVERAGE REPORT")
    log("=" * 80)

    blocks: dict[int, list[str]] = {}

    for col in master.columns:
        if col == "gene_symbol":
            continue

        m = re.match(r"^Feature(\d+)_", col)

        if not m:
            continue

        block_num = int(m.group(1))
        blocks.setdefault(block_num, []).append(col)

    for block_num in sorted(blocks):
        cols = blocks[block_num]
        mean_coverage = master[cols].notna().mean().mean() * 100
        gene_coverage = master[cols].notna().any(axis=1).mean() * 100

        log(
            f"  Feature{block_num:<2d} "
            f"{len(cols):>5d} cols  "
            f"gene coverage: {gene_coverage:6.2f}%  "
            f"value coverage: {mean_coverage:6.2f}%"
        )

    missing_blocks = [i for i in range(1, 23) if i not in blocks]

    if missing_blocks:
        log("")
        log("[WARNING] Missing from final matrix: " + ", ".join(f"F{x:02d}" for x in missing_blocks))

    total_missing = master.drop(columns=["gene_symbol"]).isna().mean().mean() * 100
    log(f"\n  Overall missing rate: {total_missing:.2f}%")


# =============================================================================
# FEATURE COUNT REPORT
# =============================================================================

def print_feature_count_report(master: pd.DataFrame) -> None:
    log("\n" + "=" * 80)
    log("FEATURE COLUMN COUNTS")
    log("=" * 80)

    total = 0

    print(f"{'Block':<8} {'N Features':>10}")
    print("=" * 80)

    for i in range(1, 23):
        prefix = f"Feature{i}_"
        n = len([c for c in master.columns if c.startswith(prefix)])
        total += n

        if n == 0:
            print(f"F{i:02d}     {n:>10}  <-- MISSING")
        else:
            print(f"F{i:02d}     {n:>10}")

    print("=" * 80)
    print(f"{'TOTAL':<8} {total:>10}")
    print("=" * 80)


# =============================================================================
# SAVE
# =============================================================================

def save(master: pd.DataFrame, outdir: Path) -> Path:
    mkdir(outdir)

    outpath = outdir / "Features_All.csv"
    master.to_csv(outpath, index=False)

    size_mb = outpath.stat().st_size / 1024 / 1024

    log(f"\n[SAVED]  {outpath.resolve()}")
    log(f"[SIZE]   {size_mb:.2f} MB")
    log(f"[SHAPE]  {master.shape[0]:,} rows x {master.shape[1]:,} columns")

    return outpath


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge all 22 feature files into Dataset/Features_All.csv"
    )
    parser.add_argument(
        "--outdir",
        default=str(DEFAULT_OUTDIR),
        help=f"Output directory. Default: {DEFAULT_OUTDIR}",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print per-file gene and feature column examples.",
    )

    args = parser.parse_args()
    outdir = Path(args.outdir)

    log("=" * 80)
    log("MERGE FEATURES — ALL 22 FEATURE BLOCKS")
    log(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"Output:  {(outdir / 'Features_All.csv').resolve()}")
    log("=" * 80)

    master = merge_all(FEATURE_FILES, verbose=args.verbose)

    print_feature_count_report(master)
    print_coverage_report(master)

    save(master, outdir)

    log("\n" + "=" * 80)
    log("DONE")
    log(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log("=" * 80)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        log("\n[STOPPED] Interrupted by user.")
        sys.exit(130)
    except Exception as exc:
        log(f"\n[ERROR] {exc}")
        raise