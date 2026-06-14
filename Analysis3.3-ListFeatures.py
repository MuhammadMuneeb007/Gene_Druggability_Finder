#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Analysis3.2-ListDatasetFeatureCounts.py

Very simple script to list number of features in each generated dataset.

Reads:
    Datasets/DatasetXXX_*/X_train.csv
    Datasets/DatasetXXX_*/feature_columns.txt  if available

Outputs:
    Analysis3.2_DatasetFeatureCounts.csv

Run:
    python Analysis3.2-ListDatasetFeatureCounts.py
"""

from pathlib import Path
import re
import pandas as pd


DATASETS_DIR = Path("Datasets")
OUT_CSV = Path("Analysis3.2_DatasetFeatureCounts.csv")


def extract_dataset_id(name: str):
    m = re.match(r"Dataset(\d+)_", name)
    return int(m.group(1)) if m else 999999


def extract_group_target(name: str):
    m = re.match(r"Dataset\d+_(.+)_(T\d+)$", name)
    if not m:
        return "", ""
    return m.group(1), m.group(2)


def count_features(dataset_dir: Path):
    feature_txt = dataset_dir / "feature_columns.txt"
    x_train = dataset_dir / "X_train.csv"

    # Prefer feature_columns.txt because it is fast and direct.
    if feature_txt.exists():
        cols = [
            line.strip()
            for line in feature_txt.read_text().splitlines()
            if line.strip()
        ]
        return len(cols)

    # Fallback: read only header of X_train.csv.
    if x_train.exists():
        header = pd.read_csv(x_train, nrows=0).columns.tolist()
        return len([c for c in header if c != "gene_symbol"])

    return 0


def main():
    if not DATASETS_DIR.exists():
        raise FileNotFoundError(f"Datasets directory not found: {DATASETS_DIR}")

    rows = []

    dataset_dirs = sorted(
        [
            p for p in DATASETS_DIR.iterdir()
            if p.is_dir() and re.match(r"^Dataset\d+_", p.name)
        ],
        key=lambda p: extract_dataset_id(p.name),
    )

    for d in dataset_dirs:
        dataset_id = extract_dataset_id(d.name)
        group_name, target_key = extract_group_target(d.name)
        n_features = count_features(d)

        rows.append({
            "dataset_id": dataset_id,
            "dataset_name": d.name,
            "group_name": group_name,
            "target_key": target_key,
            "n_features": n_features,
        })

    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)

    print("=" * 80)
    print("DATASET FEATURE COUNTS")
    print("=" * 80)
    print(df.to_string(index=False))
    print("=" * 80)
    print(f"[SAVED] {OUT_CSV.resolve()}")
    print(f"[DATASETS] {len(df)}")
    print("=" * 80)


if __name__ == "__main__":
    main()