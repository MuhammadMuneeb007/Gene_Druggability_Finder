#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Analysis_FeatureSummary.py

Read Dataset/Features_All.csv and produce a summary of every feature column:
  - Which feature block it belongs to (Feature 1 – Feature 22)
  - The database / source name
  - The column name
  - Total gene count
  - Non-missing count
  - Missing count
  - Missingness (%)
  - Mean, Std, Min, Median, Max of non-missing values

Output
------
    Supplementary_Material_2.xlsx   (Sheet 1: "Feature Summary")

Also prints a per-block summary table to the screen.

Usage
-----
    python Analysis_FeatureSummary.py
    python Analysis_FeatureSummary.py --features Dataset/Features_All.csv
    python Analysis_FeatureSummary.py --out Supplementary_Material_2.xlsx
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import (
    Alignment,
    Border,
    Font,
    PatternFill,
    Side,
)
from openpyxl.utils import get_column_letter


# =============================================================================
# SETTINGS
# =============================================================================

DEFAULT_FEATURES = Path("Dataset") / "Features_All.csv"
DEFAULT_OUT      = Path("Supplementary_Material_2.xlsx")
SHEET_NAME       = "Feature Summary"

# Human-readable source database names per feature block number
SOURCE_NAMES: Dict[int, str] = {
    1:  "DepMap (CRISPR / Essentiality)",
    2:  "STRING v12.0 (PPI Network)",
    3:  "Reactome v88 (Pathways)",
    4:  "AlphaFold DB v4 + RCSB PDB (Structure)",
    5:  "InterPro 97.0 + Pfam 36.0 (Domains)",
    6:  "UniProt / Swiss-Prot (Protein Annotation)",
    7:  "GTEx v8 (Tissue Expression)",
    8:  "Ensembl release 111 (Gene Structure)",
    9:  "gnomAD v2.1.1 + MobiDB 5.0 + DisProt (Constraint & Disorder)",
    10: "fpocket v4.0 (Binding Pockets)",
    11: "UniProt FASTA + Biopython ProtParam (Physicochemical)",
    12: "ESM-2 Embeddings + PCA (Sequence)",
    13: "NHGRI-EBI GWAS Catalog (Disease Genetics)",
    14: "Gene Ontology Annotation / QuickGO (GO Terms)",
    15: "BioGRID v4.4 (Curated PPI)",
    16: "Human Protein Atlas v23.0 (Expression & Localisation)",
    17: "CTD 2024 (Chemical-Gene Interactions)",
    18: "Mouse Genome Informatics (Mouse Knockout)",
    19: "gnomAD v4.1 Full Constraint (LOEUF / pLI)",
    20: "CORUM 4.0 (Protein Complexes)",
    21: "PhosphoSitePlus v6.7 (PTM Sites)",
    22: "Ensembl BioMart / Compara (Paralogues)",
}


# =============================================================================
# HELPERS
# =============================================================================

def log(msg: str) -> None:
    print(msg, flush=True)


def section(title: str) -> None:
    print(flush=True)
    print("=" * 100, flush=True)
    print(title, flush=True)
    print("=" * 100, flush=True)


def parse_feature_number(col: str) -> int:
    """Extract the leading feature number from a column name like Feature1_... or Feature22_..."""
    m = re.match(r"^Feature(\d+)_", col)
    return int(m.group(1)) if m else -1


def strip_prefix(col: str) -> str:
    """
    Remove the double-prefix produced by Merge_Features.py.
    e.g. Feature1_Feature1_DepMap_depmap_crispr_... -> depmap_crispr_...
    """
    # Remove two prefix layers: Feature{N}_{Label}_Feature{N}_{Label}_  (double)
    # or one layer: Feature{N}_{Label}_
    stripped = re.sub(r"^Feature\d+_\w+?_Feature\d+_\w+?_", "", col, count=1)
    if stripped == col:
        stripped = re.sub(r"^Feature\d+_\w+?_", "", col, count=1)
    return stripped if stripped else col


# =============================================================================
# COMPUTE STATISTICS
# =============================================================================

def compute_feature_stats(features_path: Path) -> pd.DataFrame:
    section("LOADING FEATURES")
    log(f"  File: {features_path}")

    df = pd.read_csv(features_path, low_memory=False)
    n_genes = len(df)

    log(f"  Genes (rows)   : {n_genes:,}")
    log(f"  Columns total  : {df.shape[1]:,}")

    feature_cols = [c for c in df.columns if c != "gene_symbol" and re.match(r"^Feature\d+_", c)]
    log(f"  Feature columns: {len(feature_cols):,}")

    section("COMPUTING STATISTICS PER FEATURE COLUMN")

    rows: List[Dict] = []

    for col in feature_cols:
        feat_num    = parse_feature_number(col)
        source_name = SOURCE_NAMES.get(feat_num, f"Feature {feat_num}")
        short_name  = strip_prefix(col)

        series = pd.to_numeric(df[col], errors="coerce")
        non_missing = int(series.notna().sum())
        missing     = int(series.isna().sum())
        miss_pct    = round(100.0 * missing / n_genes, 2) if n_genes > 0 else 0.0

        valid = series.dropna()
        rows.append({
            "Feature Block Number": feat_num,
            "Feature Block Name":   f"Feature {feat_num}",
            "Source Database":      source_name,
            "Column Name (Full)":   col,
            "Column Name (Short)":  short_name,
            "Total Genes":          n_genes,
            "Non-Missing Count":    non_missing,
            "Missing Count":        missing,
            "Missingness (%)":      miss_pct,
            "Mean":                 round(float(valid.mean()), 6)   if len(valid) > 0 else np.nan,
            "Std":                  round(float(valid.std()),  6)   if len(valid) > 1 else np.nan,
            "Min":                  round(float(valid.min()),  6)   if len(valid) > 0 else np.nan,
            "Median":               round(float(valid.median()), 6) if len(valid) > 0 else np.nan,
            "Max":                  round(float(valid.max()),  6)   if len(valid) > 0 else np.nan,
        })

    stats = pd.DataFrame(rows).sort_values(
        ["Feature Block Number", "Column Name (Short)"]
    ).reset_index(drop=True)

    return stats, n_genes


# =============================================================================
# BLOCK-LEVEL SUMMARY
# =============================================================================

def block_summary(stats: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for feat_num in sorted(stats["Feature Block Number"].unique()):
        g = stats[stats["Feature Block Number"] == feat_num]
        n_cols        = len(g)
        n_complete    = int((g["Missing Count"] == 0).sum())
        n_partial     = int((g["Missing Count"] > 0).sum())
        mean_miss     = round(float(g["Missingness (%)"].mean()), 2)
        max_miss      = round(float(g["Missingness (%)"].max()), 2)
        min_miss      = round(float(g["Missingness (%)"].min()), 2)
        src           = g["Source Database"].iloc[0]

        rows.append({
            "Feature Block":             f"Feature {feat_num}",
            "Source Database":           src,
            "Number of Feature Columns": n_cols,
            "Columns with No Missing":   n_complete,
            "Columns with Some Missing": n_partial,
            "Mean Missingness (%)":      mean_miss,
            "Min Missingness (%)":       min_miss,
            "Max Missingness (%)":       max_miss,
        })

    return pd.DataFrame(rows)


# =============================================================================
# PRINT TO SCREEN
# =============================================================================

def print_block_summary(summary: pd.DataFrame, stats: pd.DataFrame, n_genes: int) -> None:
    section("FEATURE BLOCK SUMMARY")
    log(f"  Total genes : {n_genes:,}")
    log(f"  Total cols  : {len(stats):,}\n")
    print(summary.to_string(index=False))

    section("TOP 10 MOST MISSING FEATURE COLUMNS")
    top_miss = stats.nlargest(10, "Missingness (%)")[[
        "Feature Block Name", "Column Name (Short)",
        "Non-Missing Count", "Missing Count", "Missingness (%)"
    ]]
    print(top_miss.to_string(index=False))

    section("TOP 10 MOST COMPLETE FEATURE COLUMNS")
    top_complete = stats.nsmallest(10, "Missingness (%)")[[
        "Feature Block Name", "Column Name (Short)",
        "Non-Missing Count", "Missing Count", "Missingness (%)"
    ]]
    print(top_complete.to_string(index=False))


# =============================================================================
# EXCEL STYLING CONSTANTS
# =============================================================================

HEADER_FILL   = PatternFill("solid", start_color="1F4E79")   # dark blue
BLOCK_FILLS   = [
    PatternFill("solid", start_color="D6E4F0"),  # light blue  (odd blocks)
    PatternFill("solid", start_color="EAF4FB"),  # lighter blue (even blocks)
]
SUMMARY_FILL  = PatternFill("solid", start_color="E2EFDA")   # light green
HEADER_FONT   = Font(name="Arial", bold=True, color="FFFFFF", size=10)
BODY_FONT     = Font(name="Arial", size=9)
BOLD_FONT     = Font(name="Arial", bold=True, size=9)
CENTER        = Alignment(horizontal="center", vertical="center", wrap_text=False)
LEFT          = Alignment(horizontal="left",   vertical="center", wrap_text=False)
WRAP_LEFT     = Alignment(horizontal="left",   vertical="center", wrap_text=True)
THIN          = Side(style="thin", color="BFBFBF")
THIN_BORDER   = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)
MISS_HIGH     = PatternFill("solid", start_color="FFDDC1")   # orange — high missingness
MISS_MED      = PatternFill("solid", start_color="FFF2CC")   # yellow — medium


# =============================================================================
# WRITE EXCEL
# =============================================================================

def write_excel(stats: pd.DataFrame, summary: pd.DataFrame, out_path: Path, n_genes: int) -> None:
    section("WRITING EXCEL FILE")
    log(f"  Output: {out_path}")

    wb = Workbook()

    # ----------------------------------------------------------------
    # Sheet 1 — Feature Summary (per-column detail)
    # ----------------------------------------------------------------
    ws1 = wb.active
    ws1.title = "Feature Summary"

    # --- title row ---
    ws1.merge_cells("A1:N1")
    ws1["A1"] = "Supplementary Material 2 — Feature Construction and Statistics"
    ws1["A1"].font      = Font(name="Arial", bold=True, size=13, color="1F4E79")
    ws1["A1"].alignment = CENTER

    ws1.merge_cells("A2:N2")
    ws1["A2"] = (
        f"Total genes: {n_genes:,}    |    "
        f"Total feature columns: {len(stats):,}    |    "
        f"Feature blocks: 22    |    "
        f"Overall missingness: {stats['Missingness (%)'].mean():.1f}%"
    )
    ws1["A2"].font      = Font(name="Arial", italic=True, size=9, color="404040")
    ws1["A2"].alignment = CENTER

    ws1.row_dimensions[1].height = 22
    ws1.row_dimensions[2].height = 16

    # --- column headers (row 3) ---
    headers = [
        "Block #", "Block Name", "Source Database",
        "Column (Full Name)", "Column (Short Name)",
        "Total Genes", "Non-Missing", "Missing",
        "Missingness (%)", "Mean", "Std", "Min", "Median", "Max",
    ]

    col_widths = [9, 14, 38, 55, 38, 11, 11, 11, 14, 12, 12, 12, 12, 12]

    for ci, (h, w) in enumerate(zip(headers, col_widths), start=1):
        cell = ws1.cell(row=3, column=ci, value=h)
        cell.font      = HEADER_FONT
        cell.fill      = HEADER_FILL
        cell.alignment = CENTER
        cell.border    = THIN_BORDER
        ws1.column_dimensions[get_column_letter(ci)].width = w

    ws1.row_dimensions[3].height = 28
    ws1.freeze_panes = "A4"

    # --- data rows ---
    stat_cols = [
        "Feature Block Number", "Feature Block Name", "Source Database",
        "Column Name (Full)", "Column Name (Short)",
        "Total Genes", "Non-Missing Count", "Missing Count",
        "Missingness (%)", "Mean", "Std", "Min", "Median", "Max",
    ]

    prev_block = None
    fill_idx   = 0

    for ri, (_, row) in enumerate(stats.iterrows(), start=4):
        block = row["Feature Block Number"]

        # Alternate fill per block
        if block != prev_block:
            fill_idx = 1 - fill_idx
            prev_block = block

        row_fill = BLOCK_FILLS[fill_idx]

        for ci, col_key in enumerate(stat_cols, start=1):
            val  = row[col_key]
            cell = ws1.cell(row=ri, column=ci, value=val if not (isinstance(val, float) and np.isnan(val)) else "")
            cell.font   = BODY_FONT
            cell.fill   = row_fill
            cell.border = THIN_BORDER

            # Alignment
            if ci <= 5:
                cell.alignment = LEFT
            else:
                cell.alignment = CENTER

            # Highlight high missingness in Missingness (%) column (col 9)
            if ci == 9 and isinstance(val, (int, float)) and not np.isnan(val):
                if val >= 50:
                    cell.fill = MISS_HIGH
                elif val >= 20:
                    cell.fill = MISS_MED

            # Number format for stats
            if ci in (10, 11, 12, 13, 14):
                cell.number_format = "0.000000"
            if ci == 9:
                cell.number_format = "0.00"

        ws1.row_dimensions[ri].height = 14

    # ----------------------------------------------------------------
    # Sheet 2 — Block Summary
    # ----------------------------------------------------------------
    ws2 = wb.create_sheet("Block Summary")

    ws2.merge_cells("A1:H1")
    ws2["A1"] = "Feature Block Summary — Column Counts and Missingness per Feature Set"
    ws2["A1"].font      = Font(name="Arial", bold=True, size=13, color="1F4E79")
    ws2["A1"].alignment = CENTER
    ws2.row_dimensions[1].height = 22

    sum_headers = [
        "Feature Block", "Source Database",
        "Number of Feature Columns",
        "Columns with No Missing",
        "Columns with Some Missing",
        "Mean Missingness (%)",
        "Min Missingness (%)",
        "Max Missingness (%)",
    ]
    sum_widths = [16, 46, 24, 22, 24, 20, 18, 18]

    for ci, (h, w) in enumerate(zip(sum_headers, sum_widths), start=1):
        cell = ws2.cell(row=2, column=ci, value=h)
        cell.font      = HEADER_FONT
        cell.fill      = HEADER_FILL
        cell.alignment = CENTER
        cell.border    = THIN_BORDER
        ws2.column_dimensions[get_column_letter(ci)].width = w

    ws2.row_dimensions[2].height = 28
    ws2.freeze_panes = "A3"

    for ri, (_, row) in enumerate(summary.iterrows(), start=3):
        fill = BLOCK_FILLS[(ri - 3) % 2]
        for ci, col_key in enumerate(sum_headers, start=1):
            val  = row[col_key]
            cell = ws2.cell(row=ri, column=ci, value=val)
            cell.font      = BODY_FONT
            cell.fill      = fill
            cell.border    = THIN_BORDER
            cell.alignment = LEFT if ci <= 2 else CENTER
            if ci in (6, 7, 8):
                cell.number_format = "0.00"

        ws2.row_dimensions[ri].height = 16

    # --- Totals row ---
    tr = len(summary) + 3
    ws2.cell(row=tr, column=1, value="TOTAL").font  = BOLD_FONT
    ws2.cell(row=tr, column=1).fill   = SUMMARY_FILL
    ws2.cell(row=tr, column=1).border = THIN_BORDER
    ws2.cell(row=tr, column=1).alignment = LEFT

    ws2.cell(row=tr, column=2, value="All 22 feature blocks").font = BODY_FONT
    ws2.cell(row=tr, column=2).fill   = SUMMARY_FILL
    ws2.cell(row=tr, column=2).border = THIN_BORDER
    ws2.cell(row=tr, column=2).alignment = LEFT

    ws2.cell(row=tr, column=3, value=f"=SUM(C3:C{tr-1})").font  = BOLD_FONT
    ws2.cell(row=tr, column=3).fill   = SUMMARY_FILL
    ws2.cell(row=tr, column=3).border = THIN_BORDER
    ws2.cell(row=tr, column=3).alignment = CENTER

    ws2.cell(row=tr, column=4, value=f"=SUM(D3:D{tr-1})").font  = BOLD_FONT
    ws2.cell(row=tr, column=4).fill   = SUMMARY_FILL
    ws2.cell(row=tr, column=4).border = THIN_BORDER
    ws2.cell(row=tr, column=4).alignment = CENTER

    ws2.cell(row=tr, column=5, value=f"=SUM(E3:E{tr-1})").font  = BOLD_FONT
    ws2.cell(row=tr, column=5).fill   = SUMMARY_FILL
    ws2.cell(row=tr, column=5).border = THIN_BORDER
    ws2.cell(row=tr, column=5).alignment = CENTER

    ws2.cell(row=tr, column=6, value=f"=AVERAGE(F3:F{tr-1})").font  = BOLD_FONT
    ws2.cell(row=tr, column=6).fill   = SUMMARY_FILL
    ws2.cell(row=tr, column=6).border = THIN_BORDER
    ws2.cell(row=tr, column=6).alignment = CENTER
    ws2.cell(row=tr, column=6).number_format = "0.00"

    ws2.row_dimensions[tr].height = 16

    wb.save(out_path)
    log(f"  [SAVED] {out_path.resolve()}")
    log(f"  Sheet 1: '{SHEET_NAME}' — {len(stats):,} rows (one per feature column)")
    log(f"  Sheet 2: 'Block Summary'  — {len(summary):,} rows (one per feature block)")


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarise all feature columns from Features_All.csv and save to Excel."
    )
    parser.add_argument(
        "--features", default=str(DEFAULT_FEATURES),
        help=f"Path to Features_All.csv. Default: {DEFAULT_FEATURES}",
    )
    parser.add_argument(
        "--out", default=str(DEFAULT_OUT),
        help=f"Output Excel file. Default: {DEFAULT_OUT}",
    )
    args = parser.parse_args()

    features_path = Path(args.features)
    out_path      = Path(args.out)

    if not features_path.exists():
        log(f"[ERROR] Features file not found: {features_path}")
        log("        Run Merge_Features.py first.")
        sys.exit(1)

    stats, n_genes = compute_feature_stats(features_path)
    summary        = block_summary(stats)

    print_block_summary(summary, stats, n_genes)

    write_excel(stats, summary, out_path, n_genes)

    section("DONE")
    log(f"  Excel file : {out_path.resolve()}")
    log(f"  Genes      : {n_genes:,}")
    log(f"  Features   : {len(stats):,}")
    log(f"  Blocks     : {stats['Feature Block Number'].nunique()}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        log("\n[STOPPED] Interrupted by user.")
        sys.exit(130)
    except Exception as exc:
        log(f"\n[ERROR] {exc}")
        raise