#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature10-FPocket-merge.py

Merge Feature 10 HPC chunk outputs into final gene-level tables.

Run after all SLURM array jobs finish:

    python Feature10-FPocket-merge.py

Final outputs:
    feature10_pocket_geometry/processed/feature10_structure_pocket_long.csv
    feature10_pocket_geometry/processed/feature10_pocket_geometry_gene_features.csv
    feature10_pocket_geometry/processed/feature10_pocket_geometry_hgnc_merged.csv
    feature10_pocket_geometry/feature10_pocket_geometry_summary.txt
    feature10_pocket_geometry/feature10_pocket_geometry_run_metadata.json
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd


DEFAULT_HGNC_PATH = Path("databases") / "HGNC" / "hgnc_complete_set.txt"
DEFAULT_OUTDIR = Path("feature10_pocket_geometry")


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
    return clean_text(x).upper()


def load_hgnc(hgnc_path: Path, protein_coding_only: bool = True) -> pd.DataFrame:
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

    hgnc["gene_symbol"] = hgnc["symbol"].map(normalize_symbol)

    if "uniprot_ids" not in hgnc.columns:
        hgnc["uniprot_ids"] = ""
    if "ensembl_gene_id" not in hgnc.columns:
        hgnc["ensembl_gene_id"] = ""
    if "entrez_id" not in hgnc.columns:
        hgnc["entrez_id"] = ""

    hgnc = hgnc.drop_duplicates(subset=["gene_symbol"], keep="first").copy()
    hgnc = hgnc.sort_values("gene_symbol").reset_index(drop=True)

    log(f"[HGNC ROWS] {before} -> {len(hgnc)}")
    log(f"[PROTEIN CODING ONLY] {protein_coding_only}")

    return hgnc


def aggregate_gene_features(hgnc: pd.DataFrame, long_df: pd.DataFrame) -> pd.DataFrame:
    genes = pd.DataFrame({"gene_symbol": sorted(hgnc["gene_symbol"].unique().tolist())})
    records = []

    by_gene = dict(tuple(long_df.groupby("gene_symbol"))) if not long_df.empty else {}

    metric_cols = [
        c for c in long_df.columns
        if c.startswith("feature10_")
        and c not in {"feature10_parse_error"}
    ] if not long_df.empty else []

    for gene in genes["gene_symbol"]:
        g = by_gene.get(gene, pd.DataFrame())

        if g.empty:
            records.append(
                {
                    "gene_symbol": gene,
                    "feature10_pocket_has_structure": 0,
                    "feature10_pocket_structure_count": 0,
                    "feature10_pocket_alphafold_structure_count": 0,
                    "feature10_pocket_pdb_structure_count": 0,
                    "feature10_pocket_unique_pdb_count": 0,
                    "feature10_pocket_unique_uniprot_count": 0,
                }
            )
            continue

        rec = {
            "gene_symbol": gene,
            "feature10_pocket_has_structure": 1,
            "feature10_pocket_structure_count": int(len(g)),
            "feature10_pocket_alphafold_structure_count": int(
                g["structure_source"].astype(str).str.contains("AlphaFold", case=False, na=False).sum()
            ),
            "feature10_pocket_pdb_structure_count": int(
                g["structure_source"].astype(str).str.contains("PDB", case=False, na=False).sum()
            ),
            "feature10_pocket_unique_pdb_count": int(
                g["pdb_id"].replace("", np.nan).dropna().nunique()
            ) if "pdb_id" in g.columns else 0,
            "feature10_pocket_unique_uniprot_count": int(
                g["uniprot_accession"].replace("", np.nan).dropna().nunique()
            ) if "uniprot_accession" in g.columns else 0,
        }

        for c in metric_cols:
            vals = pd.to_numeric(g[c], errors="coerce").dropna()
            if len(vals) == 0:
                continue

            base = c.replace("feature10_", "feature10_pocket_")

            rec[f"{base}_mean"] = float(vals.mean())
            rec[f"{base}_median"] = float(vals.median())
            rec[f"{base}_max"] = float(vals.max())
            rec[f"{base}_min"] = float(vals.min())

        focused = {
            "feature10_pocket_proxy_cluster_count_size_ge_3": "pocketlike_cluster_count",
            "feature10_pocket_proxy_largest_cluster_size": "largest_pocketlike_cluster_size",
            "feature10_surface_exposed_hydrophobic_fraction": "surface_hydrophobic_fraction",
            "feature10_surface_exposed_aromatic_fraction": "surface_aromatic_fraction",
            "feature10_surface_exposed_charged_fraction": "surface_charged_fraction",
            "feature10_surface_exposed_residue_fraction_proxy": "surface_exposed_fraction",
            "feature10_structure_compactness_residue_per_rg3": "compactness",
            "feature10_fpocket_pocket_count": "fpocket_pocket_count",
            "feature10_fpocket_max_druggability_score": "fpocket_druggability_score",
            "feature10_fpocket_best_volume": "fpocket_best_volume",
        }

        for original, short_name in focused.items():
            if original in g.columns:
                vals = pd.to_numeric(g[original], errors="coerce").dropna()
                if len(vals):
                    rec[f"feature10_summary_{short_name}_max"] = float(vals.max())
                    rec[f"feature10_summary_{short_name}_mean"] = float(vals.mean())

        records.append(rec)

    features = pd.DataFrame(records)

    for c in features.columns:
        if c.startswith("feature10_") and (
            "_count" in c
            or "_has_" in c
            or c.endswith("_ran")
        ):
            if features[c].dtype != object:
                features[c] = features[c].fillna(0)

    return features


def merge_with_hgnc(hgnc: pd.DataFrame, features: pd.DataFrame) -> pd.DataFrame:
    merged = hgnc.merge(features, on="gene_symbol", how="left")

    feature_cols = [c for c in features.columns if c != "gene_symbol"]

    merged["feature10_pocket_geometry_has_any_feature"] = (
        merged[feature_cols].notna().any(axis=1).astype(int)
        if feature_cols
        else 0
    )
    merged["feature10_pocket_geometry_n_nonmissing_features"] = (
        merged[feature_cols].notna().sum(axis=1)
        if feature_cols
        else 0
    )

    return merged


def write_summary(
    outdir: Path,
    hgnc: pd.DataFrame,
    long_df: pd.DataFrame,
    features: pd.DataFrame,
    merged: pd.DataFrame,
    chunk_files: List[Path],
    meta_files: List[Path],
    args: argparse.Namespace,
) -> None:
    path = outdir / "feature10_pocket_geometry_summary.txt"

    lines = []
    lines.append("Feature 10: Binding-pocket / surface-geometry features")
    lines.append("=" * 80)
    lines.append(f"Created at: {now_iso()}")
    lines.append("")
    lines.append("Inputs:")
    lines.append(f"HGNC: {args.hgnc}")
    lines.append(f"Outdir: {args.outdir}")
    lines.append(f"Chunks merged: {len(chunk_files)}")
    lines.append(f"Metadata files found: {len(meta_files)}")
    lines.append("")
    lines.append("Rows:")
    lines.append(f"HGNC rows: {hgnc.shape[0]}")
    lines.append(f"Structure long rows: {long_df.shape[0]}")
    lines.append(f"Genes with structure rows: {long_df['gene_symbol'].nunique() if not long_df.empty else 0}")
    lines.append(f"Gene feature rows: {features.shape[0]}")
    lines.append(f"Gene feature columns: {features.shape[1]}")
    lines.append(f"Merged rows: {merged.shape[0]}")
    lines.append(f"Merged columns: {merged.shape[1]}")
    lines.append("")
    lines.append("Coverage:")

    if "feature10_pocket_has_structure" in features.columns:
        lines.append(f"Structure coverage: {features['feature10_pocket_has_structure'].mean():.4f}")

    if "feature10_pocket_structure_count" in features.columns:
        lines.append(f"Median structure count per gene: {features['feature10_pocket_structure_count'].median():.2f}")

    if "feature10_summary_pocketlike_cluster_count_max" in features.columns:
        lines.append(
            "Median max pocket-like cluster count: "
            f"{features['feature10_summary_pocketlike_cluster_count_max'].median(skipna=True):.4f}"
        )

    if "feature10_summary_largest_pocketlike_cluster_size_max" in features.columns:
        lines.append(
            "Median largest pocket-like cluster size: "
            f"{features['feature10_summary_largest_pocketlike_cluster_size_max'].median(skipna=True):.4f}"
        )

    lines.append("")
    lines.append("Leakage policy:")
    lines.append("Included: structure geometry, surface exposure, residue chemistry, pocket-like cluster proxies, optional fpocket geometry scores.")
    lines.append("Excluded: ligand names, drug names, known drug-target labels, ChEMBL labels, DrugBank labels, DGIdb labels, Open Targets known-drug evidence.")
    lines.append("")
    lines.append("Interpretation:")
    lines.append("These are not known-druggability labels. They are geometry-based proxies for whether the structure has ligandable surface patches.")

    path.write_text("\n".join(lines) + "\n")
    log(f"[SUMMARY] {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge Feature 10 pocket geometry chunk outputs.")
    parser.add_argument("--hgnc", default=str(DEFAULT_HGNC_PATH), help="Path to HGNC complete set.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Feature 10 output directory.")
    parser.add_argument("--all-hgnc-genes", action="store_true", help="Use all HGNC genes, not only protein-coding.")
    args = parser.parse_args()

    outdir = mkdir(Path(args.outdir))
    processed_dir = mkdir(outdir / "processed")
    chunks_dir = mkdir(processed_dir / "chunks")

    log("=" * 100)
    log("FEATURE 10: MERGE POCKET GEOMETRY CHUNKS")
    log("=" * 100)
    log(f"[HGNC]       {args.hgnc}")
    log(f"[OUTDIR]     {outdir.resolve()}")
    log(f"[CHUNKS DIR] {chunks_dir.resolve()}")
    log("=" * 100)

    hgnc = load_hgnc(
        Path(args.hgnc),
        protein_coding_only=not args.all_hgnc_genes,
    )

    chunk_files = sorted(chunks_dir.glob("feature10_structure_pocket_long.chunk_*.csv"))
    meta_files = sorted(chunks_dir.glob("feature10_chunk_*.metadata.json")) + sorted(chunks_dir.glob("feature10_chunk_*.json")) + sorted(chunks_dir.glob("feature10_chunk_*_metadata.json"))

    if not chunk_files:
        raise FileNotFoundError(f"No chunk long files found in {chunks_dir}")

    log(f"[CHUNK LONG FILES FOUND] {len(chunk_files):,}")

    dfs = []
    empty_files = 0

    for i, path in enumerate(chunk_files, start=1):
        try:
            df = pd.read_csv(path, low_memory=False)
        except pd.errors.EmptyDataError:
            empty_files += 1
            continue

        if df.empty or len(df.columns) == 0:
            empty_files += 1
            continue

        dfs.append(df)

        if i % 25 == 0:
            log(f"[READ] {i:,}/{len(chunk_files):,}")

    if dfs:
        long_df = pd.concat(dfs, ignore_index=True)
        long_df = long_df.drop_duplicates()
    else:
        long_df = pd.DataFrame()

    long_out = processed_dir / "feature10_structure_pocket_long.csv"
    long_df.to_csv(long_out, index=False)

    log("=" * 100)
    log("[SAVED COMBINED LONG TABLE]")
    log(f"[PATH]        {long_out}")
    log(f"[SHAPE]       {long_df.shape}")
    log(f"[EMPTY FILES] {empty_files}")

    features = aggregate_gene_features(hgnc, long_df)

    gene_out = processed_dir / "feature10_pocket_geometry_gene_features.csv"
    features.to_csv(gene_out, index=False)

    log("=" * 100)
    log("[SAVED FINAL GENE FEATURES]")
    log(f"[PATH]  {gene_out}")
    log(f"[SHAPE] {features.shape}")

    merged = merge_with_hgnc(hgnc, features)

    merged_out = processed_dir / "feature10_pocket_geometry_hgnc_merged.csv"
    merged.to_csv(merged_out, index=False)

    log("=" * 100)
    log("[SAVED FINAL HGNC MERGED]")
    log(f"[PATH]  {merged_out}")
    log(f"[SHAPE] {merged.shape}")

    metadata = {
        "created_at": now_iso(),
        "hgnc": str(Path(args.hgnc)),
        "outdir": str(outdir.resolve()),
        "chunks_dir": str(chunks_dir.resolve()),
        "chunk_long_files_merged": len(chunk_files),
        "empty_chunk_files": empty_files,
        "metadata_files_found": len(meta_files),
        "outputs": {
            "structure_long": str(long_out),
            "gene_features": str(gene_out),
            "hgnc_merged": str(merged_out),
        },
        "rows": {
            "hgnc": int(hgnc.shape[0]),
            "structure_long": int(long_df.shape[0]),
            "genes_with_structure_rows": int(long_df["gene_symbol"].nunique()) if not long_df.empty else 0,
            "gene_features": int(features.shape[0]),
            "hgnc_merged": int(merged.shape[0]),
        },
        "leakage_policy": {
            "included": [
                "structure geometry",
                "surface exposure proxy",
                "surface residue chemistry",
                "pocket-like exposed hydrophobic/aromatic cluster features",
                "optional fpocket geometry scores",
            ],
            "excluded": [
                "ligand names",
                "drug names",
                "known drug-target labels",
                "ChEMBL labels",
                "DrugBank labels",
                "DGIdb labels",
                "Open Targets knownDrugs",
                "clinical target labels",
            ],
        },
    }

    metadata_out = outdir / "feature10_pocket_geometry_run_metadata.json"
    metadata_out.write_text(json.dumps(metadata, indent=2) + "\n")

    write_summary(
        outdir=outdir,
        hgnc=hgnc,
        long_df=long_df,
        features=features,
        merged=merged,
        chunk_files=chunk_files,
        meta_files=meta_files,
        args=args,
    )

    log("=" * 100)
    log("[DONE MERGE]")
    log(f"[STRUCTURE LONG] {long_out}")
    log(f"[GENE FEATURES]  {gene_out}")
    log(f"[HGNC MERGED]    {merged_out}")
    log(f"[METADATA]       {metadata_out}")
    log("=" * 100)


if __name__ == "__main__":
    main()