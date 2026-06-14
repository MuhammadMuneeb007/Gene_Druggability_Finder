#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Analysis2.1-ListFeaturesStatistics.py

Print the number of feature columns in each feature block
from Dataset/Features_All.csv.

This version recognises both:
    Feature10_xxx
and:
    feature10_xxx

Usage:
    python Analysis2.1-ListFeaturesStatistics.py
    python Analysis2.1-ListFeaturesStatistics.py --features Dataset/Features_All.csv
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


DEFAULT_FEATURES = Path("Dataset") / "Features_All.csv"

SOURCE_NAMES = {
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


def detect_feature_block(col: str) -> int | None:
    """
    Detect feature block number from column names.

    Accepts:
        Feature10_xxx
        feature10_xxx

    Returns:
        block number as int, or None if not a feature column.
    """
    if col == "gene_symbol":
        return None

    m = re.match(r"^(?:Feature|feature)(\d+)_", col)
    if not m:
        return None

    return int(m.group(1))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Print feature column count per feature block."
    )
    parser.add_argument(
        "--features",
        default=str(DEFAULT_FEATURES),
        help=f"Path to Features_All.csv. Default: {DEFAULT_FEATURES}",
    )
    args = parser.parse_args()

    features_path = Path(args.features)

    if not features_path.exists():
        print(f"[ERROR] File not found: {features_path}")
        print("        Run Merge_Features.py first.")
        sys.exit(1)

    print(f"Loading: {features_path}")

    # Read headers only, so this is fast even for large files.
    df = pd.read_csv(features_path, nrows=0)

    counts: dict[int, int] = {}
    lowercase_feature_cols = []

    for col in df.columns:
        block = detect_feature_block(col)
        if block is None:
            continue

        counts[block] = counts.get(block, 0) + 1

        if col.startswith("feature"):
            lowercase_feature_cols.append(col)

    print()
    print("=" * 80)
    print(f"{'Block':<8} {'Feature Set':<42} {'N Features':>10}  Source")
    print("=" * 80)

    total = 0
    for n in sorted(counts):
        c = counts[n]
        src = SOURCE_NAMES.get(n, "Unknown")
        name = f"Feature {n}"
        print(f"F{n:02d}     {name:<42} {c:>10}  {src}")
        total += c

    print("=" * 80)
    print(f"{'TOTAL':<8} {str(len(counts)) + ' feature blocks':<42} {total:>10}")
    print("=" * 80)
    print()

    if lowercase_feature_cols:
        print("[NOTE]")
        print(
            f"Detected {len(lowercase_feature_cols)} lowercase feature columns, "
            "for example:"
        )
        for col in lowercase_feature_cols[:10]:
            print(f"  - {col}")
        if len(lowercase_feature_cols) > 10:
            print(f"  ... and {len(lowercase_feature_cols) - 10} more")
        print()
        print(
            "This is okay for counting, but for consistency you may want to "
            "rename these columns to uppercase FeatureXX_ during Merge_Features.py."
        )
        print()


if __name__ == "__main__":
    main()