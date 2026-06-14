#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature12_ProteinEmbeddings.py

Feature 12: UniProt ProtT5 protein embeddings.

Purpose
-------
Download UniProt precomputed human per-protein embeddings and convert them into
gene-level features for ML.

Input embedding:
    feature_databases/ProteinEmbeddings/UniProt_ProtT5/UP000005640_9606/per-protein.h5

Embedding source:
    UniProt human reference proteome UP000005640_9606
    per-protein ProtT5 embeddings
    1024 dimensions per UniProt accession

HGNC:
    databases/HGNC/hgnc_complete_set.txt

Outputs
-------
feature12_protein_embeddings/
    processed/feature12_embedding_protein_long.csv
    processed/feature12_embedding_gene_raw_mean.csv.gz
    processed/feature12_embedding_gene_pca_features.csv
    processed/feature12_embedding_hgnc_merged.csv
    processed/feature12_embedding_pca_model.pkl
    feature12_embedding_summary.txt
    feature12_embedding_run_metadata.json

Run
---
    python Feature12_ProteinEmbeddings.py

Force download:
    python Feature12_ProteinEmbeddings.py --download

Use 128 PCA dimensions:
    python Feature12_ProteinEmbeddings.py --pca-components 128

Fast test:
    python Feature12_ProteinEmbeddings.py --limit-genes 1000

Leakage policy
--------------
Safe. These are sequence-derived protein language model embeddings.
They do NOT use:
    ChEMBL
    DrugBank
    DGIdb
    Open Targets
    Pharos
    clinical target labels
    known drug-target labels
"""

from __future__ import annotations

import argparse
import json
import pickle
import re
import sys
import time
import urllib.request
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple

import numpy as np
import pandas as pd


try:
    import h5py
except Exception:
    h5py = None

try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
except Exception:
    PCA = None
    StandardScaler = None


# =============================================================================
# DEFAULTS
# =============================================================================

DEFAULT_HGNC_PATH = Path("databases") / "HGNC" / "hgnc_complete_set.txt"
DEFAULT_DBDIR = Path("feature_databases") / "ProteinEmbeddings" / "UniProt_ProtT5" / "UP000005640_9606"
DEFAULT_OUTDIR = Path("feature12_protein_embeddings")

EMBEDDING_URL = (
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
    "knowledgebase/embeddings/UP000005640_9606/per-protein.h5"
)

METALINK_URL = (
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
    "knowledgebase/embeddings/UP000005640_9606/RELEASE.metalink"
)


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
    return clean_text(x).upper()


def clean_uniprot(x: Any) -> str:
    s = clean_text(x)
    if not s:
        return ""
    s = s.replace("UniProtKB:", "")
    s = s.replace("sp|", "")
    s = s.replace("tr|", "")
    s = s.split("|")[0]
    s = s.split("-")[0]
    s = s.split(".")[0]
    return s.strip()


def split_uniprot_ids(x: Any) -> List[str]:
    s = clean_text(x)
    if not s:
        return []

    out = []

    for part in re.split(r"[|;, ]+", s):
        acc = clean_uniprot(part)
        if acc and acc not in out:
            out.append(acc)

    return out


def clean_ensembl_gene_id(x: Any) -> str:
    s = clean_text(x)
    if not s:
        return ""
    return s.split(".")[0]


def download_file(url: str, outpath: Path) -> None:
    mkdir(outpath.parent)

    if outpath.exists() and outpath.stat().st_size > 0:
        log(f"[DOWNLOAD SKIP] Already exists: {outpath}")
        return

    log(f"[DOWNLOAD] {url}")
    log(f"[TO]       {outpath}")

    tmp = outpath.with_suffix(outpath.suffix + ".tmp")

    with urllib.request.urlopen(url) as response:
        total = response.length
        downloaded = 0
        last_print = time.time()

        with open(tmp, "wb") as f:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break

                f.write(chunk)
                downloaded += len(chunk)

                if time.time() - last_print > 5:
                    if total:
                        pct = 100 * downloaded / total
                        log(f"[DOWNLOAD] {downloaded / 1024 / 1024:.1f} MB / {total / 1024 / 1024:.1f} MB ({pct:.1f}%)")
                    else:
                        log(f"[DOWNLOAD] {downloaded / 1024 / 1024:.1f} MB")
                    last_print = time.time()

    tmp.rename(outpath)

    log(f"[DOWNLOAD DONE] {outpath} ({outpath.stat().st_size / 1024 / 1024:.2f} MB)")


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
        raise RuntimeError("HGNC file must contain a symbol column.")

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
    if "name" not in hgnc.columns:
        hgnc["name"] = ""

    hgnc["ensembl_gene_id_clean"] = hgnc["ensembl_gene_id"].map(clean_ensembl_gene_id)
    hgnc = hgnc.drop_duplicates(subset=["gene_symbol"], keep="first").copy()

    if limit_genes and limit_genes > 0:
        hgnc = hgnc.head(limit_genes).copy()

    log(f"[HGNC ROWS] {before} -> {len(hgnc)}")
    log(f"[PROTEIN CODING ONLY] {protein_coding_only}")
    log(f"[GENES WITH UNIPROT] {(hgnc['uniprot_ids'].astype(str).str.len() > 0).sum()}")

    return hgnc


def build_uniprot_to_genes(hgnc: pd.DataFrame) -> Dict[str, Set[str]]:
    mapping: Dict[str, Set[str]] = defaultdict(set)

    for _, row in hgnc.iterrows():
        gene = row["gene_symbol"]

        for acc in split_uniprot_ids(row.get("uniprot_ids", "")):
            mapping[acc].add(gene)

    log("=" * 100)
    log("[UNIPROT ? HGNC MAP]")
    log(f"[UNIPROT ACCESSIONS IN HGNC] {len(mapping)}")

    return mapping


# =============================================================================
# EMBEDDING LOADING
# =============================================================================

def inspect_h5_embedding_file(h5_path: Path) -> Tuple[int, int, str]:
    if h5py is None:
        raise ImportError("h5py is required. Install with: conda install -c conda-forge h5py -y")

    if not h5_path.exists():
        raise FileNotFoundError(f"Embedding H5 not found: {h5_path}")

    with h5py.File(h5_path, "r") as f:
        keys = list(f.keys())
        n = len(keys)

        if n == 0:
            raise RuntimeError("H5 file contains zero embeddings.")

        first_key = keys[0]
        dim = int(f[first_key].shape[0])

    log("=" * 100)
    log("[H5 INSPECTION]")
    log(f"[FILE] {h5_path}")
    log(f"[SIZE MB] {h5_path.stat().st_size / 1024 / 1024:.2f}")
    log(f"[N EMBEDDINGS] {n}")
    log(f"[FIRST ACCESSION] {first_key}")
    log(f"[EMBEDDING DIM] {dim}")

    return n, dim, first_key


def build_gene_embedding_arrays(
    h5_path: Path,
    hgnc: pd.DataFrame,
    uniprot_to_genes: Dict[str, Set[str]],
    processed_dir: Path,
) -> Tuple[pd.DataFrame, pd.DataFrame, np.ndarray, List[str], int]:
    """
    Read H5 embeddings, map UniProt accession ? HGNC gene,
    aggregate multiple protein embeddings by mean per gene.

    Returns:
        protein_long
        gene_meta
        gene_embedding_matrix
        gene_order
        embedding_dim
    """
    if h5py is None:
        raise ImportError("h5py is required. Install with: conda install -c conda-forge h5py -y")

    genes_all = sorted(hgnc["gene_symbol"].unique().tolist())

    gene_vectors: Dict[str, List[np.ndarray]] = defaultdict(list)
    protein_rows = []

    n_seen = 0
    n_mapped_accessions = 0
    n_gene_links = 0
    embedding_dim = None

    log("=" * 100)
    log("[READ EMBEDDINGS AND MAP TO GENES]")

    with h5py.File(h5_path, "r") as f:
        for accession in f.keys():
            n_seen += 1

            acc = clean_uniprot(accession)
            genes = uniprot_to_genes.get(acc, set())

            if not genes:
                continue

            vec = np.asarray(f[accession][:], dtype=np.float32)

            if embedding_dim is None:
                embedding_dim = int(vec.shape[0])

            if vec.ndim != 1:
                continue

            n_mapped_accessions += 1

            norm = float(np.linalg.norm(vec))

            for gene in genes:
                gene_vectors[gene].append(vec)
                n_gene_links += 1

                protein_rows.append(
                    {
                        "gene_symbol": gene,
                        "uniprot_accession": acc,
                        "embedding_dim": int(vec.shape[0]),
                        "embedding_norm": norm,
                        "embedding_mean_value": float(np.mean(vec)),
                        "embedding_std_value": float(np.std(vec)),
                        "embedding_min_value": float(np.min(vec)),
                        "embedding_max_value": float(np.max(vec)),
                    }
                )

            if n_seen % 5000 == 0:
                log(f"[H5] seen={n_seen:,} mapped_accessions={n_mapped_accessions:,} gene_links={n_gene_links:,}")

    if embedding_dim is None:
        raise RuntimeError("No embeddings mapped to HGNC genes. Check HGNC uniprot_ids mapping.")

    protein_long = pd.DataFrame(protein_rows)

    protein_long_path = processed_dir / "feature12_embedding_protein_long.csv"
    protein_long.to_csv(protein_long_path, index=False)

    log(f"[SAVED PROTEIN LONG] {protein_long_path}")
    log(f"[PROTEIN LONG SHAPE] {protein_long.shape}")

    gene_order = genes_all
    X = np.full((len(gene_order), embedding_dim), np.nan, dtype=np.float32)

    gene_meta_rows = []

    for i, gene in enumerate(gene_order):
        vectors = gene_vectors.get(gene, [])

        if vectors:
            arr = np.vstack(vectors).astype(np.float32)
            mean_vec = np.mean(arr, axis=0)
            X[i, :] = mean_vec

            norms = np.linalg.norm(arr, axis=1)

            gene_meta_rows.append(
                {
                    "gene_symbol": gene,
                    "feature12_embedding_has_embedding": 1,
                    "feature12_embedding_protein_count": int(arr.shape[0]),
                    "feature12_embedding_dim": int(embedding_dim),
                    "feature12_embedding_mean_norm": float(np.mean(norms)),
                    "feature12_embedding_max_norm": float(np.max(norms)),
                    "feature12_embedding_min_norm": float(np.min(norms)),
                    "feature12_embedding_std_norm": float(np.std(norms, ddof=1)) if len(norms) > 1 else 0.0,
                    "feature12_embedding_gene_mean_vector_norm": float(np.linalg.norm(mean_vec)),
                    "feature12_embedding_gene_mean_value": float(np.mean(mean_vec)),
                    "feature12_embedding_gene_std_value": float(np.std(mean_vec)),
                    "feature12_embedding_gene_min_value": float(np.min(mean_vec)),
                    "feature12_embedding_gene_max_value": float(np.max(mean_vec)),
                }
            )
        else:
            gene_meta_rows.append(
                {
                    "gene_symbol": gene,
                    "feature12_embedding_has_embedding": 0,
                    "feature12_embedding_protein_count": 0,
                    "feature12_embedding_dim": int(embedding_dim),
                    "feature12_embedding_mean_norm": np.nan,
                    "feature12_embedding_max_norm": np.nan,
                    "feature12_embedding_min_norm": np.nan,
                    "feature12_embedding_std_norm": np.nan,
                    "feature12_embedding_gene_mean_vector_norm": np.nan,
                    "feature12_embedding_gene_mean_value": np.nan,
                    "feature12_embedding_gene_std_value": np.nan,
                    "feature12_embedding_gene_min_value": np.nan,
                    "feature12_embedding_gene_max_value": np.nan,
                }
            )

    gene_meta = pd.DataFrame(gene_meta_rows)

    log("=" * 100)
    log("[GENE EMBEDDING MATRIX]")
    log(f"[GENES] {len(gene_order)}")
    log(f"[DIM] {embedding_dim}")
    log(f"[GENES WITH EMBEDDING] {int(gene_meta['feature12_embedding_has_embedding'].sum())}")
    log(f"[COVERAGE] {gene_meta['feature12_embedding_has_embedding'].mean():.4f}")

    return protein_long, gene_meta, X, gene_order, embedding_dim


# =============================================================================
# SAVE RAW MEAN EMBEDDINGS
# =============================================================================

def save_gene_raw_mean_embeddings(
    X: np.ndarray,
    gene_order: List[str],
    processed_dir: Path,
    save_raw: bool,
) -> Path:
    path = processed_dir / "feature12_embedding_gene_raw_mean.csv.gz"

    if not save_raw:
        log("[SKIP RAW EMBEDDINGS] --no-save-raw was used.")
        return path

    log("=" * 100)
    log("[SAVE RAW GENE MEAN EMBEDDINGS]")
    log("This can be a moderately large compressed CSV.")

    dim = X.shape[1]
    cols = [f"feature12_raw_embed_{i + 1:04d}" for i in range(dim)]

    df = pd.DataFrame(X, columns=cols)
    df.insert(0, "gene_symbol", gene_order)

    df.to_csv(path, index=False, compression="gzip")

    log(f"[SAVED RAW GENE MEAN EMBEDDINGS] {path}")
    log(f"[SHAPE] {df.shape}")

    return path


# =============================================================================
# PCA FEATURES
# =============================================================================

def build_pca_features(
    X: np.ndarray,
    gene_order: List[str],
    gene_meta: pd.DataFrame,
    pca_components: int,
    processed_dir: Path,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    if PCA is None or StandardScaler is None:
        raise ImportError("scikit-learn is required. Install with: conda install -c conda-forge scikit-learn -y")

    log("=" * 100)
    log("[BUILD PCA FEATURES]")

    has_embedding = gene_meta["feature12_embedding_has_embedding"].astype(int).to_numpy() == 1
    X_obs = X[has_embedding, :]

    if X_obs.shape[0] == 0:
        raise RuntimeError("No genes with embeddings for PCA.")

    n_components = min(pca_components, X_obs.shape[0] - 1, X_obs.shape[1])

    if n_components < 1:
        raise RuntimeError("Not enough genes with embeddings for PCA.")

    log(f"[PCA INPUT GENES] {X_obs.shape[0]}")
    log(f"[PCA INPUT DIM] {X_obs.shape[1]}")
    log(f"[PCA COMPONENTS] {n_components}")

    scaler = StandardScaler(with_mean=True, with_std=True)
    X_scaled = scaler.fit_transform(X_obs)

    pca = PCA(n_components=n_components, random_state=42)
    X_pca_obs = pca.fit_transform(X_scaled)

    X_pca_all = np.full((len(gene_order), n_components), np.nan, dtype=np.float32)
    X_pca_all[has_embedding, :] = X_pca_obs.astype(np.float32)

    pca_cols = [f"feature12_embed_pc{i + 1:03d}" for i in range(n_components)]

    pca_df = pd.DataFrame(X_pca_all, columns=pca_cols)
    pca_df.insert(0, "gene_symbol", gene_order)

    # Add metadata and norm summaries.
    pca_df = gene_meta.merge(pca_df, on="gene_symbol", how="left")

    # PCA summary features.
    if n_components >= 1:
        pca_values = pca_df[pca_cols].to_numpy(dtype=float)
        pca_df["feature12_embedding_pca_vector_norm"] = np.sqrt(np.nansum(pca_values * pca_values, axis=1))
        pca_df.loc[pca_df["feature12_embedding_has_embedding"] == 0, "feature12_embedding_pca_vector_norm"] = np.nan

    outpath = processed_dir / "feature12_embedding_gene_pca_features.csv"
    pca_df.to_csv(outpath, index=False)

    model_path = processed_dir / "feature12_embedding_pca_model.pkl"
    with open(model_path, "wb") as f:
        pickle.dump(
            {
                "scaler": scaler,
                "pca": pca,
                "pca_columns": pca_cols,
                "n_components": n_components,
                "explained_variance_ratio": pca.explained_variance_ratio_,
            },
            f,
        )

    explained = float(np.sum(pca.explained_variance_ratio_))

    log(f"[SAVED PCA FEATURES] {outpath}")
    log(f"[SAVED PCA MODEL] {model_path}")
    log(f"[PCA EXPLAINED VARIANCE TOTAL] {explained:.4f}")

    pca_info = {
        "n_components": int(n_components),
        "explained_variance_ratio_total": explained,
        "explained_variance_ratio_first_10": [float(x) for x in pca.explained_variance_ratio_[:10]],
        "pca_model_path": str(model_path),
        "pca_features_path": str(outpath),
    }

    return pca_df, pca_info


# =============================================================================
# MERGE WITH HGNC
# =============================================================================

def merge_with_hgnc(
    hgnc: pd.DataFrame,
    pca_features: pd.DataFrame,
    processed_dir: Path,
) -> pd.DataFrame:
    merged = hgnc.merge(pca_features, on="gene_symbol", how="left")

    feature_cols = [c for c in pca_features.columns if c != "gene_symbol"]

    merged["feature12_embedding_has_any_feature"] = merged[feature_cols].notna().any(axis=1).astype(int)
    merged["feature12_embedding_n_nonmissing_features"] = merged[feature_cols].notna().sum(axis=1)

    outpath = processed_dir / "feature12_embedding_hgnc_merged.csv"
    merged.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED HGNC MERGED]")
    log(f"[PATH]  {outpath}")
    log(f"[SHAPE] {merged.shape}")

    return merged


# =============================================================================
# SUMMARY
# =============================================================================

def write_summary(
    outdir: Path,
    hgnc: pd.DataFrame,
    protein_long: pd.DataFrame,
    gene_meta: pd.DataFrame,
    pca_features: pd.DataFrame,
    merged: pd.DataFrame,
    h5_path: Path,
    raw_path: Path,
    pca_info: Dict[str, Any],
    args: argparse.Namespace,
) -> None:
    path = outdir / "feature12_embedding_summary.txt"

    lines = []
    lines.append("Feature 12: UniProt ProtT5 protein embeddings")
    lines.append("=" * 80)
    lines.append(f"Created at: {now_iso()}")
    lines.append("")
    lines.append("Inputs:")
    lines.append(f"HGNC: {args.hgnc}")
    lines.append(f"Embedding H5: {h5_path}")
    lines.append(f"Embedding URL: {EMBEDDING_URL}")
    lines.append("")
    lines.append("Rows:")
    lines.append(f"HGNC rows: {hgnc.shape[0]}")
    lines.append(f"Protein embedding mapped rows: {protein_long.shape[0]}")
    lines.append(f"Genes with embeddings: {int(gene_meta['feature12_embedding_has_embedding'].sum())}")
    lines.append(f"Gene embedding coverage: {gene_meta['feature12_embedding_has_embedding'].mean():.4f}")
    lines.append(f"PCA feature rows: {pca_features.shape[0]}")
    lines.append(f"PCA feature columns: {pca_features.shape[1]}")
    lines.append(f"Merged rows: {merged.shape[0]}")
    lines.append(f"Merged columns: {merged.shape[1]}")
    lines.append("")
    lines.append("Embedding:")
    lines.append(f"Raw embedding dim: {int(gene_meta['feature12_embedding_dim'].max())}")
    lines.append(f"PCA components: {pca_info['n_components']}")
    lines.append(f"PCA total explained variance: {pca_info['explained_variance_ratio_total']:.4f}")
    lines.append(f"First 10 PCA explained variance ratios: {pca_info['explained_variance_ratio_first_10']}")
    lines.append("")
    lines.append("Outputs:")
    lines.append(f"Protein long table: feature12_protein_embeddings/processed/feature12_embedding_protein_long.csv")
    lines.append(f"Raw gene mean embeddings: {raw_path}")
    lines.append(f"PCA features: feature12_protein_embeddings/processed/feature12_embedding_gene_pca_features.csv")
    lines.append(f"HGNC merged: feature12_protein_embeddings/processed/feature12_embedding_hgnc_merged.csv")
    lines.append(f"PCA model: {pca_info['pca_model_path']}")
    lines.append("")
    lines.append("Leakage policy:")
    lines.append("Included: UniProt ProtT5 sequence embeddings only.")
    lines.append("Excluded: clinical labels, known drug-target labels, ChEMBL, DrugBank, DGIdb, Open Targets, Pharos.")

    path.write_text("\n".join(lines) + "\n")

    log(f"[SUMMARY] {path}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description="Build Feature 12 UniProt ProtT5 protein embedding features.")
    parser.add_argument("--hgnc", default=str(DEFAULT_HGNC_PATH), help="Path to HGNC complete set.")
    parser.add_argument("--dbdir", default=str(DEFAULT_DBDIR), help="Embedding database directory.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory.")
    parser.add_argument("--embedding-file", default="", help="Manual path to per-protein.h5.")
    parser.add_argument("--download", action="store_true", help="Download UniProt human per-protein.h5 if missing.")
    parser.add_argument("--force-download", action="store_true", help="Re-download even if file exists.")
    parser.add_argument("--pca-components", type=int, default=64, help="Number of PCA components to save.")
    parser.add_argument("--no-save-raw", action="store_true", help="Do not save raw 1024-dim mean gene embeddings CSV.")
    parser.add_argument("--all-hgnc-genes", action="store_true", help="Use all HGNC genes, not only protein-coding.")
    parser.add_argument("--limit-genes", type=int, default=0, help="Debug only: first N HGNC genes.")
    args = parser.parse_args()

    if h5py is None:
        raise ImportError(
            "h5py is missing. Install it with:\n"
            "  conda install -c conda-forge h5py -y\n"
            "or:\n"
            "  pip install h5py"
        )

    if PCA is None or StandardScaler is None:
        raise ImportError(
            "scikit-learn is missing. Install it with:\n"
            "  conda install -c conda-forge scikit-learn -y\n"
            "or:\n"
            "  pip install scikit-learn"
        )

    dbdir = mkdir(Path(args.dbdir))
    outdir = mkdir(Path(args.outdir))
    processed_dir = mkdir(outdir / "processed")

    h5_path = Path(args.embedding_file) if args.embedding_file else dbdir / "per-protein.h5"
    metalink_path = dbdir / "RELEASE.metalink"

    log("=" * 100)
    log("FEATURE 12: UNIPROT PROTT5 PROTEIN EMBEDDINGS")
    log("=" * 100)
    log(f"[HGNC]           {args.hgnc}")
    log(f"[DBDIR]          {dbdir.resolve()}")
    log(f"[OUTDIR]         {outdir.resolve()}")
    log(f"[EMBEDDING H5]   {h5_path}")
    log(f"[DOWNLOAD]       {args.download}")
    log(f"[FORCE DOWNLOAD] {args.force_download}")
    log(f"[PCA COMPONENTS] {args.pca_components}")
    log(f"[SAVE RAW]       {not args.no_save_raw}")
    log(f"[LIMIT GENES]    {args.limit_genes if args.limit_genes else 'none'}")
    log("=" * 100)

    if args.force_download and h5_path.exists():
        log(f"[REMOVE EXISTING] {h5_path}")
        h5_path.unlink()

    if args.download or not h5_path.exists():
        download_file(EMBEDDING_URL, h5_path)
        try:
            download_file(METALINK_URL, metalink_path)
        except Exception as exc:
            log(f"[WARNING] Could not download RELEASE.metalink: {exc}")

    inspect_h5_embedding_file(h5_path)

    hgnc = load_hgnc(
        Path(args.hgnc),
        protein_coding_only=not args.all_hgnc_genes,
        limit_genes=args.limit_genes,
    )

    uniprot_to_genes = build_uniprot_to_genes(hgnc)

    protein_long, gene_meta, X, gene_order, embedding_dim = build_gene_embedding_arrays(
        h5_path=h5_path,
        hgnc=hgnc,
        uniprot_to_genes=uniprot_to_genes,
        processed_dir=processed_dir,
    )

    raw_path = save_gene_raw_mean_embeddings(
        X=X,
        gene_order=gene_order,
        processed_dir=processed_dir,
        save_raw=not args.no_save_raw,
    )

    pca_features, pca_info = build_pca_features(
        X=X,
        gene_order=gene_order,
        gene_meta=gene_meta,
        pca_components=args.pca_components,
        processed_dir=processed_dir,
    )

    merged = merge_with_hgnc(
        hgnc=hgnc,
        pca_features=pca_features,
        processed_dir=processed_dir,
    )

    metadata = {
        "created_at": now_iso(),
        "hgnc": str(Path(args.hgnc)),
        "dbdir": str(dbdir.resolve()),
        "outdir": str(outdir.resolve()),
        "embedding_file": str(h5_path),
        "embedding_url": EMBEDDING_URL,
        "protein_coding_only": not args.all_hgnc_genes,
        "limit_genes": args.limit_genes,
        "embedding_dim": int(embedding_dim),
        "pca_components": int(pca_info["n_components"]),
        "pca_explained_variance_total": float(pca_info["explained_variance_ratio_total"]),
        "outputs": {
            "protein_long": str(processed_dir / "feature12_embedding_protein_long.csv"),
            "gene_raw_mean": str(raw_path),
            "gene_pca_features": str(processed_dir / "feature12_embedding_gene_pca_features.csv"),
            "hgnc_merged": str(processed_dir / "feature12_embedding_hgnc_merged.csv"),
            "pca_model": str(processed_dir / "feature12_embedding_pca_model.pkl"),
        },
        "leakage_policy": {
            "included": [
                "UniProt ProtT5 per-protein embeddings",
                "HGNC UniProt accession mapping",
                "gene-level mean embedding",
                "PCA-compressed embedding features",
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

    with open(outdir / "feature12_embedding_run_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    write_summary(
        outdir=outdir,
        hgnc=hgnc,
        protein_long=protein_long,
        gene_meta=gene_meta,
        pca_features=pca_features,
        merged=merged,
        h5_path=h5_path,
        raw_path=raw_path,
        pca_info=pca_info,
        args=args,
    )

    log("=" * 100)
    log("[DONE]")
    log(f"[PROTEIN LONG] {processed_dir / 'feature12_embedding_protein_long.csv'}")
    log(f"[RAW MEAN]     {raw_path}")
    log(f"[PCA FEATURES] {processed_dir / 'feature12_embedding_gene_pca_features.csv'}")
    log(f"[HGNC MERGED]  {processed_dir / 'feature12_embedding_hgnc_merged.csv'}")
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