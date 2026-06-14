#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Dataset1_CheckGroundTruth.py

Simple diagnostic script for the generated ground-truth table.

It prints:
1. Dataset shape
2. All columns
3. Missing values per column
4. Counts for important labels
5. Counts for evidence columns
6. Druggability class distribution
7. Label confidence distribution
8. Top genes by druggability score
9. Suspicious sanity checks
"""

from pathlib import Path
import pandas as pd


INPUT_FILE = Path("Step0_Output/03_HumanGene_DruggabilityLabels.csv")
SUMMARY_FILE = Path("Step0_Output/04_label_summary.csv")


def print_section(title):
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


def count_binary_column(df, col):
    if col not in df.columns:
        return None

    s = pd.to_numeric(df[col], errors="coerce").fillna(0)
    n1 = int((s == 1).sum())
    n0 = int((s == 0).sum())
    n_missing = int(df[col].isna().sum())
    percent_1 = round(100 * n1 / max(len(df), 1), 3)

    return {
        "column": col,
        "n_1": n1,
        "n_0": n0,
        "n_missing": n_missing,
        "percent_1": percent_1,
    }


def main():
    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Could not find file: {INPUT_FILE}")

    df = pd.read_csv(INPUT_FILE, low_memory=False)

    print_section("GROUND-TRUTH TABLE BASIC INFORMATION")
    print(f"Input file: {INPUT_FILE}")
    print(f"Rows/genes: {len(df)}")
    print(f"Columns: {len(df.columns)}")

    print_section("ALL COLUMNS")
    for i, col in enumerate(df.columns, start=1):
        print(f"{i:03d}. {col}")

    print_section("MISSING VALUES PER COLUMN")
    missing = (
        df.isna()
        .sum()
        .reset_index()
        .rename(columns={"index": "column", 0: "missing_count"})
    )
    missing["missing_percent"] = (100 * missing["missing_count"] / len(df)).round(3)
    print(missing.to_string(index=False))

    print_section("MAIN FINAL LABEL COUNTS")

    main_label_cols = [
        "final_any_druggable_label",
        "clinical_target_label",
        "clinical_investigation_label",
        "chemical_tractable_label",
        "small_molecule_druggable_label",
        "biologic_druggable_label",
        "dgidb_interaction_supported_label",
        "potentially_druggable_category_label",
        "unknown_dark_label",
    ]

    rows = []
    for col in main_label_cols:
        result = count_binary_column(df, col)
        if result is not None:
            rows.append(result)

    label_summary = pd.DataFrame(rows)
    print(label_summary.to_string(index=False))

    print_section("EVIDENCE COLUMN COUNTS")

    evidence_cols = [
        # ChEMBL
        "chembl_has_target",
        "chembl_has_mechanism",
        "chembl_has_phase4_or_approved",
        "chembl_has_potent_compound",

        # Open Targets
        "opentargets_has_target",
        "ot_has_clinical_indication",
        "ot_has_small_molecule_tractability",
        "ot_has_antibody_tractability",
        "ot_has_protac_tractability",
        "ot_has_other_modality_tractability",

        # DGIdb
        "dgidb_has_interaction",
        "dgidb_has_gene_record",
        "dgidb_has_any_category",

        # Pharos
        "pharos_api_available",
        "pharos_is_tclin",
        "pharos_is_tchem",
        "pharos_is_tbio",
        "pharos_is_tdark",
    ]

    rows = []
    for col in evidence_cols:
        result = count_binary_column(df, col)
        if result is not None:
            rows.append(result)

    evidence_summary = pd.DataFrame(rows)
    print(evidence_summary.to_string(index=False))

    print_section("DRUGGABILITY CLASS DISTRIBUTION")
    if "druggability_class" in df.columns:
        class_counts = (
            df["druggability_class"]
            .fillna("MISSING")
            .value_counts()
            .reset_index()
        )
        class_counts.columns = ["druggability_class", "count"]
        class_counts["percent"] = (100 * class_counts["count"] / len(df)).round(3)
        print(class_counts.to_string(index=False))
    else:
        print("Column not found: druggability_class")

    print_section("LABEL CONFIDENCE DISTRIBUTION")
    if "label_confidence" in df.columns:
        conf_counts = (
            df["label_confidence"]
            .fillna("MISSING")
            .value_counts()
            .reset_index()
        )
        conf_counts.columns = ["label_confidence", "count"]
        conf_counts["percent"] = (100 * conf_counts["count"] / len(df)).round(3)
        print(conf_counts.to_string(index=False))
    else:
        print("Column not found: label_confidence")

    print_section("PHAROS TDL DISTRIBUTION")
    if "pharos_tdl" in df.columns:
        pharos_counts = (
            df["pharos_tdl"]
            .fillna("MISSING")
            .value_counts()
            .reset_index()
        )
        pharos_counts.columns = ["pharos_tdl", "count"]
        pharos_counts["percent"] = (100 * pharos_counts["count"] / len(df)).round(3)
        print(pharos_counts.to_string(index=False))
    else:
        print("Column not found: pharos_tdl")

    print_section("DRUGGABILITY SCORE SUMMARY")
    if "druggability_score" in df.columns:
        score = pd.to_numeric(df["druggability_score"], errors="coerce")
        print(score.describe().to_string())

        print("\nTop 30 genes by druggability_score:")
        show_cols = [
            "gene_symbol",
            "gene_name",
            "druggability_score",
            "druggability_class",
            "label_confidence",
            "final_any_druggable_label",
            "clinical_target_label",
            "chemical_tractable_label",
            "biologic_druggable_label",
            "dgidb_interaction_supported_label",
            "potentially_druggable_category_label",
            "pharos_tdl",
            "evidence_sources",
        ]
        show_cols = [c for c in show_cols if c in df.columns]

        top = (
            df.assign(druggability_score_numeric=score)
            .sort_values("druggability_score_numeric", ascending=False)
            .head(30)
        )
        print(top[show_cols].to_string(index=False))
    else:
        print("Column not found: druggability_score")

    print_section("IMPORTANT SANITY CHECKS")

    checks = []

    def get_sum(col):
        if col not in df.columns:
            return None
        return int(pd.to_numeric(df[col], errors="coerce").fillna(0).sum())

    n_genes = len(df)

    final_any = get_sum("final_any_druggable_label")
    clinical = get_sum("clinical_target_label")
    chembl_phase4 = get_sum("chembl_has_phase4_or_approved")
    chembl_potent = get_sum("chembl_has_potent_compound")
    ot_sm = get_sum("ot_has_small_molecule_tractability")
    ot_ab = get_sum("ot_has_antibody_tractability")
    ot_protac = get_sum("ot_has_protac_tractability")
    ot_other = get_sum("ot_has_other_modality_tractability")
    biologic = get_sum("biologic_druggable_label")
    clinical_inv = get_sum("clinical_investigation_label")

    checks.append(("final_any_druggable_label", final_any))
    checks.append(("clinical_target_label", clinical))
    checks.append(("clinical_investigation_label", clinical_inv))
    checks.append(("chembl_has_phase4_or_approved", chembl_phase4))
    checks.append(("chembl_has_potent_compound", chembl_potent))
    checks.append(("ot_has_small_molecule_tractability", ot_sm))
    checks.append(("ot_has_antibody_tractability", ot_ab))
    checks.append(("ot_has_protac_tractability", ot_protac))
    checks.append(("ot_has_other_modality_tractability", ot_other))
    checks.append(("biologic_druggable_label", biologic))

    for name, value in checks:
        if value is not None:
            print(f"{name:45s}: {value}")

    print("\nPotential warnings:")

    if final_any is not None and final_any / n_genes > 0.70:
        print(
            f"[WARNING] final_any_druggable_label is very high: "
            f"{final_any}/{n_genes} = {100 * final_any / n_genes:.2f}%"
        )
        print("          This may be too broad for a binary druggable/not-druggable ground truth.")

    if biologic is not None and biologic / n_genes > 0.50:
        print(
            f"[WARNING] biologic_druggable_label is very high: "
            f"{biologic}/{n_genes} = {100 * biologic / n_genes:.2f}%"
        )
        print("          Open Targets biologic/other modality parsing may be too permissive.")

    if chembl_phase4 == 0:
        print("[WARNING] ChEMBL approved/phase4 evidence is zero.")
        print("          This is suspicious. Check ChEMBL mechanism/drug mechanism table parsing.")

    if chembl_potent == 0:
        print("[WARNING] ChEMBL potent compound evidence is zero.")
        print("          This is suspicious. Check activities.pchembl_value parsing and table joins.")

    if ot_sm == 0 and ot_ab == 0 and ot_protac == 0:
        print("[WARNING] Open Targets small molecule, antibody, and PROTAC tractability are all zero.")
        print("          This is suspicious. Check parsing of the Open Targets tractability column.")

    if clinical_inv == 0:
        print("[WARNING] clinical_investigation_label is zero.")
        print("          This may mean Open Targets clinical_indication was not mapped to target genes correctly.")

    print_section("RECOMMENDED GROUND-TRUTH COLUMNS TO USE")
    print("Strict clinical drug target label:")
    print("  clinical_target_label")

    print("\nBroad druggability label:")
    print("  final_any_druggable_label")

    print("\nRanking score:")
    print("  druggability_score")

    print("\nInterpretable class:")
    print("  druggability_class")

    print("\nConfidence:")
    print("  label_confidence")

    print("\nDo NOT use these as model input features:")
    print("  ChEMBL / OpenTargets / DGIdb / Pharos evidence columns")
    print("  They were used to construct the labels and would leak the answer.")

    print_section("DONE")


if __name__ == "__main__":
    main()