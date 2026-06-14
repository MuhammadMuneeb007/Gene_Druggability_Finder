#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import pandas as pd
import os

BASE_DIR = Path("/data/ascher02/uqmmune1/DrugableGeneFinder/Final")
DATASETS_DIR = BASE_DIR / "Datasets"
CATALOGUE = DATASETS_DIR / "dataset_catalogue.csv"
OUTFILE = BASE_DIR / "Dataset_Training_Audit.csv"

REQUIRED_DATASET_FILES = [
    "X_train.csv",
    "y_train.csv",
    "meta.json",
    "feature_columns.txt",
    "missingness_summary.csv",
]

REQUIRED_TRAINING_FILES = [
    "00_training_metadata.json",
    "05_model_comparison_summary.csv",
    "06_best_model_summary.csv",
    "12_permutation_test_results.csv",
]

def file_status(path: Path):
    if not path.exists():
        return "missing", 0
    size = path.stat().st_size
    if size == 0:
        return "empty", 0
    return "ok", size

def find_dataset_dir(dataset_id, dataset_name):
    direct = DATASETS_DIR / str(dataset_name)
    if direct.exists():
        return direct

    matches = sorted(DATASETS_DIR.glob(f"Dataset{int(dataset_id):03d}_*"))
    if matches:
        return matches[0]

    return direct

def main():
    print("=" * 100)
    print("DATASET + TRAINING AUDIT")
    print("=" * 100)
    print(f"Base dir     : {BASE_DIR}")
    print(f"Datasets dir : {DATASETS_DIR}")
    print(f"Catalogue    : {CATALOGUE}")
    print("=" * 100)

    if not CATALOGUE.exists():
        raise FileNotFoundError(f"Catalogue not found: {CATALOGUE}")

    cat = pd.read_csv(CATALOGUE, low_memory=False)
    print(f"Catalogue rows: {len(cat)}")

    rows = []

    for _, r in cat.iterrows():
        dataset_id = int(r["dataset_id"])
        dataset_name = str(r["dataset_name"])
        ddir = find_dataset_dir(dataset_id, dataset_name)

        row = {
            "dataset_id": dataset_id,
            "dataset_name": dataset_name,
            "dataset_dir": str(ddir),
            "dataset_folder_exists": ddir.exists(),
        }

        # Dataset generation files
        dataset_ok = True
        for fname in REQUIRED_DATASET_FILES:
            status, size = file_status(ddir / fname)
            row[f"{fname}_status"] = status
            row[f"{fname}_size"] = size
            if fname in ["X_train.csv", "y_train.csv", "meta.json"] and status != "ok":
                dataset_ok = False

        row["dataset_generated_ok"] = dataset_ok

        # Training folder and result files
        training_dir = ddir / "Training"
        row["training_folder_exists"] = training_dir.exists()

        training_results_ok = True
        for fname in REQUIRED_TRAINING_FILES:
            status, size = file_status(training_dir / fname)
            row[f"Training/{fname}_status"] = status
            row[f"Training/{fname}_size"] = size
            if fname in ["05_model_comparison_summary.csv", "06_best_model_summary.csv"] and status != "ok":
                training_results_ok = False

        row["training_results_ok"] = training_results_ok

        # Simple final status
        if not ddir.exists():
            final_status = "DATASET_FOLDER_MISSING"
        elif not dataset_ok:
            final_status = "DATASET_FILES_MISSING_OR_EMPTY"
        elif not training_dir.exists():
            final_status = "NO_TRAINING_FOLDER"
        elif not training_results_ok:
            final_status = "TRAINING_FOLDER_BUT_RESULTS_MISSING_OR_EMPTY"
        else:
            final_status = "OK"

        row["final_status"] = final_status
        rows.append(row)

    audit = pd.DataFrame(rows)
    audit.to_csv(OUTFILE, index=False)

    print("\nSUMMARY")
    print("=" * 100)
    print(audit["final_status"].value_counts(dropna=False).to_string())

    print("\nDATASET FILE COUNTS")
    print("=" * 100)
    print(f"Dataset folders existing : {audit['dataset_folder_exists'].sum()} / {len(audit)}")
    print(f"Datasets generated OK    : {audit['dataset_generated_ok'].sum()} / {len(audit)}")
    print(f"Training folders existing: {audit['training_folder_exists'].sum()} / {len(audit)}")
    print(f"Training results OK      : {audit['training_results_ok'].sum()} / {len(audit)}")

    print("\nPROBLEM DATASETS")
    print("=" * 100)
    bad = audit[audit["final_status"] != "OK"]
    if bad.empty:
        print("No problems found.")
    else:
        cols = ["dataset_id", "dataset_name", "final_status", "dataset_dir"]
        print(bad[cols].to_string(index=False))

    print("\nSaved full audit to:")
    print(OUTFILE)

if __name__ == "__main__":
    main()