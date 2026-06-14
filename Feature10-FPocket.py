#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature10-FPocket.py

Feature 10: Binding-pocket / surface-geometry features.

HPC array version
-----------------
Each SLURM array task processes one chunk of HGNC protein-coding genes.

Example:
    python Feature10-FPocket.py 1 --chunk-size 100

Chunk logic:
    job 1 -> genes[0:100]
    job 2 -> genes[100:200]
    job 3 -> genes[200:300]

Outputs per job:
    feature10_pocket_geometry/processed/chunks/
        feature10_structure_pocket_long.chunk_0001.csv
        feature10_pocket_geometry_gene_features.chunk_0001.csv
        feature10_chunk_0001_metadata.json

After all jobs finish:
    python Feature10-FPocket-merge.py
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict, deque
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pandas as pd


# =============================================================================
# DEFAULTS
# =============================================================================

DEFAULT_HGNC_PATH = Path("databases") / "HGNC" / "hgnc_complete_set.txt"
DEFAULT_FEATURE4_DIR = Path("feature4_structure")
DEFAULT_OUTDIR = Path("feature10_pocket_geometry")

HYDROPHOBIC = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO"}
AROMATIC = {"PHE", "TRP", "TYR", "HIS"}
POSITIVE = {"ARG", "LYS", "HIS"}
NEGATIVE = {"ASP", "GLU"}
POLAR = {"SER", "THR", "ASN", "GLN", "TYR", "CYS"}
CHARGED = POSITIVE | NEGATIVE
SMALL = {"GLY", "ALA", "SER", "CYS"}

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


# =============================================================================
# BASIC HELPERS
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
    return s.split("-")[0].split(".")[0]


def clean_pdb_id(x: Any) -> str:
    s = clean_text(x).upper()
    s = re.sub(r"[^A-Z0-9]", "", s)
    if len(s) >= 4:
        return s[:4]
    return s


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


def open_text_maybe_gzip(path: Path):
    name = str(path).lower()
    if name.endswith(".gz"):
        return gzip.open(path, "rt", errors="ignore")
    return open(path, "rt", errors="ignore")


def fpocket_available() -> bool:
    return shutil.which("fpocket") is not None


def infer_source_from_path(path: Path) -> str:
    name = path.name.lower()
    full = str(path).lower()

    if "alphafold" in full or name.startswith("af-") or "af_" in name:
        return "AlphaFold"

    if re.match(r"^[0-9a-z]{4}\.(pdb|cif|ent|pdb\.gz|ent\.gz)$", name):
        return "PDB"

    if name.endswith(".pdb") or name.endswith(".pdb.gz") or name.endswith(".cif") or name.endswith(".ent"):
        if "pdb" in full:
            return "PDB"

    return "Unknown"


# =============================================================================
# HGNC / MAPPING
# =============================================================================

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
    log(f"[GENES WITH UNIPROT] {(hgnc['uniprot_ids'].astype(str).str.len() > 0).sum()}")

    return hgnc


def build_uniprot_to_gene(hgnc: pd.DataFrame) -> Dict[str, Set[str]]:
    out: Dict[str, Set[str]] = defaultdict(set)

    for _, row in hgnc.iterrows():
        gene = row["gene_symbol"]
        for acc in split_uniprot_ids(row.get("uniprot_ids", "")):
            out[acc].add(gene)

    return out


def load_feature4_pdb_map(feature4_dir: Path) -> Dict[str, Set[str]]:
    possible = [
        feature4_dir / "processed" / "feature4_structure_pdb_long.csv",
        feature4_dir / "processed" / "feature4_pdb_long.csv",
        feature4_dir / "feature4_structure_pdb_long.csv",
    ]

    mapping: Dict[str, Set[str]] = defaultdict(set)

    path = None
    for p in possible:
        if p.exists():
            path = p
            break

    if path is None:
        log("[FEATURE4 PDB MAP] Not found. Will infer from file names only.")
        return mapping

    log(f"[READ FEATURE4 PDB MAP] {path}")

    df = pd.read_csv(path, low_memory=False)
    df.columns = [c.strip() for c in df.columns]

    gene_col = "gene_symbol" if "gene_symbol" in df.columns else None

    pdb_col = None
    for c in df.columns:
        if c.lower() in {"pdb_id", "pdb", "structure_id"}:
            pdb_col = c
            break

    if gene_col is None or pdb_col is None:
        log("[FEATURE4 PDB MAP] Could not find gene_symbol/pdb_id columns.")
        return mapping

    for _, row in df.iterrows():
        gene = normalize_symbol(row.get(gene_col, ""))
        pdb = clean_pdb_id(row.get(pdb_col, ""))
        if gene and pdb:
            mapping[pdb].add(gene)

    log(f"[FEATURE4 PDB MAP] PDB IDs mapped: {len(mapping)}")
    return mapping


def infer_genes_from_filename(
    path: Path,
    uniprot_to_genes: Dict[str, Set[str]],
    pdb_to_genes: Dict[str, Set[str]],
) -> Tuple[Set[str], str, str]:
    name = path.name
    genes: Set[str] = set()
    uniprot_acc = ""
    pdb_id = ""

    m = re.search(r"AF-([A-Z0-9]+)-F\d+", name, flags=re.IGNORECASE)
    if m:
        uniprot_acc = clean_uniprot(m.group(1))
        if uniprot_acc in uniprot_to_genes:
            genes.update(uniprot_to_genes[uniprot_acc])

    if not genes:
        for acc, gset in uniprot_to_genes.items():
            if acc and acc.lower() in name.lower():
                uniprot_acc = acc
                genes.update(gset)
                break

    m2 = re.match(r"^([0-9A-Za-z]{4})", name)
    if m2:
        pdb_id = clean_pdb_id(m2.group(1))
        if pdb_id in pdb_to_genes:
            genes.update(pdb_to_genes[pdb_id])

    return genes, uniprot_acc, pdb_id


# =============================================================================
# STRUCTURE DISCOVERY
# =============================================================================

def find_structure_files(feature4_dir: Path, extra_dirs: List[Path]) -> List[Path]:
    roots = [feature4_dir] + extra_dirs
    files: List[Path] = []

    for root in roots:
        if not root.exists():
            continue

        for pattern in ["*.pdb", "*.pdb.gz", "*.ent", "*.ent.gz", "*.cif", "*.cif.gz"]:
            files.extend(list(root.rglob(pattern)))

    files = [p for p in files if p.exists() and p.stat().st_size > 0]
    files = [p for p in files if "_out" not in str(p)]
    files = sorted(set(files))
    return files


# =============================================================================
# PDB PARSING
# =============================================================================

def parse_pdb_lines(lines: Iterable[str]) -> Dict[str, Any]:
    atom_coords = []
    ca_coords = []
    ca_resnames = []
    ca_reskeys = []

    seen_res = set()

    for line in lines:
        if not line.startswith("ATOM"):
            continue

        if len(line) < 54:
            continue

        atom_name = line[12:16].strip()
        resname = line[17:20].strip().upper()
        chain = line[21].strip()
        resseq = line[22:26].strip()
        icode = line[26].strip()

        if resname not in AA3_TO_1:
            continue

        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except Exception:
            continue

        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if element == "H":
            continue

        atom_coords.append((x, y, z))

        reskey = (chain, resseq, icode)
        if atom_name == "CA" and reskey not in seen_res:
            ca_coords.append((x, y, z))
            ca_resnames.append(resname)
            ca_reskeys.append(reskey)
            seen_res.add(reskey)

    return {
        "atom_coords": np.asarray(atom_coords, dtype=float) if atom_coords else np.empty((0, 3)),
        "ca_coords": np.asarray(ca_coords, dtype=float) if ca_coords else np.empty((0, 3)),
        "ca_resnames": ca_resnames,
        "ca_reskeys": ca_reskeys,
    }


def parse_pdb_file(path: Path) -> Dict[str, Any]:
    with open_text_maybe_gzip(path) as f:
        return parse_pdb_lines(f)


# =============================================================================
# GEOMETRY FEATURES
# =============================================================================

def pairwise_neighbor_counts(coords: np.ndarray, radius: float = 10.0, chunk: int = 1000) -> np.ndarray:
    n = coords.shape[0]
    if n == 0:
        return np.asarray([], dtype=int)

    counts = np.zeros(n, dtype=int)
    r2 = radius * radius

    for i in range(0, n, chunk):
        block = coords[i:i + chunk]
        diff = block[:, None, :] - coords[None, :, :]
        dist2 = np.sum(diff * diff, axis=2)
        c = np.sum(dist2 <= r2, axis=1) - 1
        counts[i:i + chunk] = c

    return counts


def cluster_points(coords: np.ndarray, distance_cutoff: float = 12.0) -> List[List[int]]:
    n = coords.shape[0]
    if n == 0:
        return []

    r2 = distance_cutoff * distance_cutoff
    visited = np.zeros(n, dtype=bool)
    clusters = []

    for start in range(n):
        if visited[start]:
            continue

        cluster = []
        queue = deque([start])
        visited[start] = True

        while queue:
            i = queue.popleft()
            cluster.append(i)

            diff = coords - coords[i]
            dist2 = np.sum(diff * diff, axis=1)
            neigh = np.where((dist2 <= r2) & (~visited))[0]

            for j in neigh:
                visited[j] = True
                queue.append(int(j))

        clusters.append(cluster)

    return clusters


def residue_fraction(resnames: List[str], residue_set: Set[str]) -> float:
    if not resnames:
        return np.nan
    return float(sum(1 for r in resnames if r in residue_set) / len(resnames))


def empty_geometry_features() -> Dict[str, Any]:
    keys = [
        "structure_radius_gyration",
        "structure_mean_distance_to_centroid",
        "structure_max_distance_to_centroid",
        "structure_bounding_box_volume",
        "structure_compactness_residue_per_rg3",
        "structure_residue_density_bbox",
        "surface_mean_neighbor_count_10a",
        "surface_median_neighbor_count_10a",
        "surface_min_neighbor_count_10a",
        "surface_max_neighbor_count_10a",
        "surface_mean_neighbor_count_12a",
        "surface_exposed_residue_count_proxy",
        "surface_exposed_residue_fraction_proxy",
        "surface_exposed_hydrophobic_fraction",
        "surface_exposed_aromatic_fraction",
        "surface_exposed_positive_fraction",
        "surface_exposed_negative_fraction",
        "surface_exposed_charged_fraction",
        "surface_exposed_polar_fraction",
        "surface_exposed_small_fraction",
        "surface_buried_hydrophobic_fraction",
        "surface_buried_aromatic_fraction",
        "pocket_proxy_candidate_residue_count",
        "pocket_proxy_candidate_residue_fraction",
        "pocket_proxy_cluster_count_all",
        "pocket_proxy_cluster_count_size_ge_3",
        "pocket_proxy_largest_cluster_size",
        "pocket_proxy_largest_cluster_fraction",
        "pocket_proxy_candidate_hydrophobic_fraction",
        "pocket_proxy_candidate_aromatic_fraction",
    ]
    return {k: np.nan for k in keys}


def calculate_structure_features(parsed: Dict[str, Any]) -> Dict[str, Any]:
    atom_coords = parsed["atom_coords"]
    ca_coords = parsed["ca_coords"]
    ca_resnames = parsed["ca_resnames"]

    n_atoms = int(atom_coords.shape[0])
    n_residues = int(ca_coords.shape[0])

    out = {
        "structure_atom_count": n_atoms,
        "structure_residue_count": n_residues,
    }

    if n_atoms == 0 or n_residues == 0:
        out.update(empty_geometry_features())
        return out

    centroid = np.mean(ca_coords, axis=0)
    centered = ca_coords - centroid
    dist_centroid = np.sqrt(np.sum(centered * centered, axis=1))

    radius_gyration = float(np.sqrt(np.mean(np.sum(centered * centered, axis=1))))
    max_distance_centroid = float(np.max(dist_centroid))
    mean_distance_centroid = float(np.mean(dist_centroid))

    min_xyz = np.min(ca_coords, axis=0)
    max_xyz = np.max(ca_coords, axis=0)
    bbox = max_xyz - min_xyz
    bbox_volume = float(np.prod(np.maximum(bbox, 1.0)))

    compactness = float(n_residues / (radius_gyration ** 3)) if radius_gyration > 0 else np.nan
    residue_density_bbox = float(n_residues / bbox_volume) if bbox_volume > 0 else np.nan

    neigh10 = pairwise_neighbor_counts(ca_coords, radius=10.0)
    neigh12 = pairwise_neighbor_counts(ca_coords, radius=12.0)

    exposed_mask = neigh10 <= 8

    exposed_resnames = [r for r, m in zip(ca_resnames, exposed_mask) if m]
    buried_resnames = [r for r, m in zip(ca_resnames, exposed_mask) if not m]

    exposed_count = int(np.sum(exposed_mask))
    exposed_fraction = float(exposed_count / n_residues) if n_residues else np.nan

    candidate_mask = []
    for res, exposed in zip(ca_resnames, exposed_mask):
        is_candidate = bool(exposed and (res in HYDROPHOBIC or res in AROMATIC))
        candidate_mask.append(is_candidate)

    candidate_mask = np.asarray(candidate_mask, dtype=bool)
    candidate_coords = ca_coords[candidate_mask]
    candidate_resnames = [r for r, m in zip(ca_resnames, candidate_mask) if m]

    clusters = cluster_points(candidate_coords, distance_cutoff=12.0)
    cluster_sizes = [len(c) for c in clusters]
    pocketlike_clusters = [c for c in clusters if len(c) >= 3]

    largest_cluster_size = max(cluster_sizes) if cluster_sizes else 0
    pocketlike_cluster_count = len(pocketlike_clusters)

    out.update(
        {
            "structure_radius_gyration": radius_gyration,
            "structure_mean_distance_to_centroid": mean_distance_centroid,
            "structure_max_distance_to_centroid": max_distance_centroid,
            "structure_bounding_box_volume": bbox_volume,
            "structure_compactness_residue_per_rg3": compactness,
            "structure_residue_density_bbox": residue_density_bbox,

            "surface_mean_neighbor_count_10a": float(np.mean(neigh10)) if len(neigh10) else np.nan,
            "surface_median_neighbor_count_10a": float(np.median(neigh10)) if len(neigh10) else np.nan,
            "surface_min_neighbor_count_10a": float(np.min(neigh10)) if len(neigh10) else np.nan,
            "surface_max_neighbor_count_10a": float(np.max(neigh10)) if len(neigh10) else np.nan,
            "surface_mean_neighbor_count_12a": float(np.mean(neigh12)) if len(neigh12) else np.nan,

            "surface_exposed_residue_count_proxy": exposed_count,
            "surface_exposed_residue_fraction_proxy": exposed_fraction,

            "surface_exposed_hydrophobic_fraction": residue_fraction(exposed_resnames, HYDROPHOBIC),
            "surface_exposed_aromatic_fraction": residue_fraction(exposed_resnames, AROMATIC),
            "surface_exposed_positive_fraction": residue_fraction(exposed_resnames, POSITIVE),
            "surface_exposed_negative_fraction": residue_fraction(exposed_resnames, NEGATIVE),
            "surface_exposed_charged_fraction": residue_fraction(exposed_resnames, CHARGED),
            "surface_exposed_polar_fraction": residue_fraction(exposed_resnames, POLAR),
            "surface_exposed_small_fraction": residue_fraction(exposed_resnames, SMALL),

            "surface_buried_hydrophobic_fraction": residue_fraction(buried_resnames, HYDROPHOBIC),
            "surface_buried_aromatic_fraction": residue_fraction(buried_resnames, AROMATIC),

            "pocket_proxy_candidate_residue_count": int(np.sum(candidate_mask)),
            "pocket_proxy_candidate_residue_fraction": float(np.mean(candidate_mask)) if len(candidate_mask) else np.nan,
            "pocket_proxy_cluster_count_all": int(len(clusters)),
            "pocket_proxy_cluster_count_size_ge_3": int(pocketlike_cluster_count),
            "pocket_proxy_largest_cluster_size": int(largest_cluster_size),
            "pocket_proxy_largest_cluster_fraction": float(largest_cluster_size / n_residues) if n_residues else np.nan,
            "pocket_proxy_candidate_hydrophobic_fraction": residue_fraction(candidate_resnames, HYDROPHOBIC),
            "pocket_proxy_candidate_aromatic_fraction": residue_fraction(candidate_resnames, AROMATIC),
        }
    )

    return out


# =============================================================================
# OPTIONAL FPOCKET
# =============================================================================

def run_fpocket_and_parse(path: Path, tmpdir: Path, timeout: int = 180) -> Dict[str, Any]:
    out = {
        "fpocket_ran": 0,
        "fpocket_pocket_count": np.nan,
        "fpocket_best_score": np.nan,
        "fpocket_mean_score": np.nan,
        "fpocket_max_druggability_score": np.nan,
        "fpocket_mean_druggability_score": np.nan,
        "fpocket_best_volume": np.nan,
        "fpocket_mean_volume": np.nan,
    }

    if not fpocket_available():
        return out

    work = tmpdir / re.sub(r"[^A-Za-z0-9_.-]", "_", path.stem)
    mkdir(work)

    local_path = work / path.name

    try:
        if str(path).endswith(".gz"):
            ungz_path = work / path.name.replace(".gz", "")
            with gzip.open(path, "rb") as fin, open(ungz_path, "wb") as fout:
                fout.write(fin.read())
            local_path = ungz_path
        else:
            shutil.copy2(path, local_path)
    except Exception:
        return out

    if not str(local_path).lower().endswith(".pdb"):
        return out

    try:
        cmd = ["fpocket", "-f", str(local_path)]
        subprocess.run(
            cmd,
            cwd=str(work),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=timeout,
            check=False,
        )
    except Exception:
        return out

    outdir = work / f"{local_path.stem}_out"
    info_file = outdir / f"{local_path.stem}_info.txt"

    if not info_file.exists():
        return out

    scores = []
    drug_scores = []
    volumes = []
    current = {}

    with open(info_file, "r", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.lower().startswith("pocket"):
                if current:
                    if "score" in current:
                        scores.append(current["score"])
                    if "druggability_score" in current:
                        drug_scores.append(current["druggability_score"])
                    if "volume" in current:
                        volumes.append(current["volume"])
                current = {}

            low = line.lower()

            if "score" in low and "druggability" not in low:
                nums = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)
                if nums:
                    current["score"] = float(nums[-1])

            if "druggability score" in low:
                nums = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)
                if nums:
                    current["druggability_score"] = float(nums[-1])

            if "volume" in low:
                nums = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)
                if nums:
                    current["volume"] = float(nums[-1])

    if current:
        if "score" in current:
            scores.append(current["score"])
        if "druggability_score" in current:
            drug_scores.append(current["druggability_score"])
        if "volume" in current:
            volumes.append(current["volume"])

    pocket_count = max(len(scores), len(drug_scores), len(volumes))

    out["fpocket_ran"] = 1
    out["fpocket_pocket_count"] = pocket_count

    if scores:
        out["fpocket_best_score"] = float(np.max(scores))
        out["fpocket_mean_score"] = float(np.mean(scores))

    if drug_scores:
        out["fpocket_max_druggability_score"] = float(np.max(drug_scores))
        out["fpocket_mean_druggability_score"] = float(np.mean(drug_scores))

    if volumes:
        out["fpocket_best_volume"] = float(np.max(volumes))
        out["fpocket_mean_volume"] = float(np.mean(volumes))

    return out


# =============================================================================
# PROCESSING
# =============================================================================

def get_chunk_genes(hgnc: pd.DataFrame, job_index: int, chunk_size: int) -> List[str]:
    if job_index < 1:
        raise ValueError("--job-index must be 1-based, e.g. 1, 2, 3.")

    genes = sorted(hgnc["gene_symbol"].dropna().map(normalize_symbol).unique().tolist())

    start = (job_index - 1) * chunk_size
    end = job_index * chunk_size

    return genes[start:end]


def process_structure_file(
    path: Path,
    genes: Set[str],
    uniprot_acc: str,
    pdb_id: str,
    run_fpocket: bool,
    tmpdir: Path,
) -> List[Dict[str, Any]]:
    source = infer_source_from_path(path)

    try:
        parsed = parse_pdb_file(path)
        feats = calculate_structure_features(parsed)
        parse_error = ""
    except Exception as exc:
        parse_error = str(exc)
        feats = empty_geometry_features()
        feats["structure_atom_count"] = 0
        feats["structure_residue_count"] = 0

    if run_fpocket:
        fp = run_fpocket_and_parse(path, tmpdir)
    else:
        fp = {
            "fpocket_ran": 0,
            "fpocket_pocket_count": np.nan,
            "fpocket_best_score": np.nan,
            "fpocket_mean_score": np.nan,
            "fpocket_max_druggability_score": np.nan,
            "fpocket_mean_druggability_score": np.nan,
            "fpocket_best_volume": np.nan,
            "fpocket_mean_volume": np.nan,
        }

    rows = []
    for gene in sorted(genes):
        rec = {
            "gene_symbol": gene,
            "structure_file": str(path),
            "structure_source": source,
            "uniprot_accession": uniprot_acc,
            "pdb_id": pdb_id,
            "feature10_parse_error": parse_error,
        }

        for k, v in feats.items():
            rec[f"feature10_{k}"] = v

        for k, v in fp.items():
            rec[f"feature10_{k}"] = v

        rows.append(rec)

    return rows


def aggregate_gene_features_for_chunk(
    chunk_genes: List[str],
    long_df: pd.DataFrame,
) -> pd.DataFrame:
    records = []

    by_gene = dict(tuple(long_df.groupby("gene_symbol"))) if not long_df.empty else {}

    metric_cols = [
        c for c in long_df.columns
        if c.startswith("feature10_")
        and c not in {"feature10_parse_error"}
    ] if not long_df.empty else []

    for gene in chunk_genes:
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


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description="Feature 10 pocket geometry array-job runner.")

    parser.add_argument(
        "job_index_positional",
        nargs="?",
        type=int,
        default=None,
        help="1-based chunk index. Can also use --job-index or SLURM_ARRAY_TASK_ID.",
    )

    parser.add_argument("--job-index", type=int, default=None, help="1-based chunk index.")
    parser.add_argument("--chunk-size", type=int, default=100, help="Genes per array job.")
    parser.add_argument("--hgnc", default=str(DEFAULT_HGNC_PATH), help="Path to HGNC complete set.")
    parser.add_argument("--feature4-dir", default=str(DEFAULT_FEATURE4_DIR), help="Feature4 structure output directory.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory.")
    parser.add_argument("--extra-structure-dirs", nargs="*", default=[], help="Extra directories to scan for PDB/CIF files.")
    parser.add_argument("--run-fpocket", action="store_true", help="Run fpocket if installed.")
    parser.add_argument("--all-hgnc-genes", action="store_true", help="Use all HGNC genes, not only protein-coding.")
    parser.add_argument("--limit-structures", type=int, default=0, help="Debug only: process first N structures.")
    parser.add_argument("--force", action="store_true", help="Overwrite existing chunk outputs.")

    args = parser.parse_args()

    job_index = args.job_index
    if job_index is None:
        job_index = args.job_index_positional
    if job_index is None:
        slurm_id = os.environ.get("SLURM_ARRAY_TASK_ID", "").strip()
        if slurm_id:
            job_index = int(slurm_id)

    if job_index is None:
        raise RuntimeError("No job index supplied. Use python Feature10-FPocket.py 1 or submit with SLURM_ARRAY_TASK_ID.")

    outdir = mkdir(Path(args.outdir))
    processed_dir = mkdir(outdir / "processed")
    chunks_dir = mkdir(processed_dir / "chunks")

    chunk_tag = f"chunk_{job_index:04d}"

    long_out = chunks_dir / f"feature10_structure_pocket_long.{chunk_tag}.csv"
    gene_out = chunks_dir / f"feature10_pocket_geometry_gene_features.{chunk_tag}.csv"
    meta_out = chunks_dir / f"feature10_{chunk_tag}_metadata.json"

    if long_out.exists() and gene_out.exists() and not args.force:
        log("=" * 100)
        log(f"[SKIP] Existing chunk outputs found for {chunk_tag}. Use --force to rerun.")
        log(f"[LONG] {long_out}")
        log(f"[GENE] {gene_out}")
        log("=" * 100)
        return

    feature4_dir = Path(args.feature4_dir)
    extra_dirs = [Path(x) for x in args.extra_structure_dirs]

    log("=" * 100)
    log("FEATURE 10: BINDING-POCKET / SURFACE-GEOMETRY FEATURES - ARRAY CHUNK")
    log("=" * 100)
    log(f"[JOB INDEX]             {job_index}")
    log(f"[CHUNK SIZE]            {args.chunk_size}")
    log(f"[HGNC]                  {args.hgnc}")
    log(f"[FEATURE4 DIR]          {feature4_dir}")
    log(f"[OUTDIR]                {outdir.resolve()}")
    log(f"[CHUNKS DIR]            {chunks_dir.resolve()}")
    log(f"[EXTRA STRUCTURE DIRS]  {extra_dirs}")
    log(f"[RUN FPOCKET]           {args.run_fpocket}")
    log(f"[FPOCKET AVAILABLE]     {fpocket_available()}")
    log(f"[LIMIT STRUCTURES]      {args.limit_structures if args.limit_structures else 'none'}")
    log("=" * 100)

    hgnc = load_hgnc(
        Path(args.hgnc),
        protein_coding_only=not args.all_hgnc_genes,
    )

    all_genes = sorted(hgnc["gene_symbol"].dropna().map(normalize_symbol).unique().tolist())
    chunk_genes = get_chunk_genes(hgnc, job_index=job_index, chunk_size=args.chunk_size)
    chunk_gene_set = set(chunk_genes)

    start = (job_index - 1) * args.chunk_size
    end = job_index * args.chunk_size

    log("=" * 100)
    log(f"[TOTAL GENES]           {len(all_genes):,}")
    log(f"[CHUNK RANGE]           genes[{start}:{end}]")
    log(f"[CHUNK GENES]           {len(chunk_genes):,}")

    if not chunk_genes:
        log("[EMPTY CHUNK] No genes assigned to this job. Writing empty outputs and exiting.")

        pd.DataFrame().to_csv(long_out, index=False)
        pd.DataFrame({"gene_symbol": []}).to_csv(gene_out, index=False)

        meta = {
            "created_at": now_iso(),
            "job_index": job_index,
            "chunk_size": args.chunk_size,
            "start": start,
            "end": end,
            "n_total_genes": len(all_genes),
            "n_chunk_genes": 0,
            "status": "empty_chunk",
        }
        meta_out.write_text(json.dumps(meta, indent=2) + "\n")
        return

    log(f"[FIRST GENE]            {chunk_genes[0]}")
    log(f"[LAST GENE]             {chunk_genes[-1]}")

    uniprot_to_genes_all = build_uniprot_to_gene(hgnc)
    pdb_to_genes_all = load_feature4_pdb_map(feature4_dir)

    structure_files = find_structure_files(feature4_dir, extra_dirs)

    log("=" * 100)
    log(f"[STRUCTURE FILES FOUND] {len(structure_files):,}")

    if args.limit_structures:
        structure_files = structure_files[: args.limit_structures]

    log(f"[STRUCTURE FILES TO SCAN] {len(structure_files):,}")
    log("=" * 100)

    rows = []
    scanned = 0
    matched_structures = 0
    matched_rows = 0

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)

        for i, path in enumerate(structure_files, start=1):
            scanned += 1

            genes_all, acc, pdb_id = infer_genes_from_filename(
                path,
                uniprot_to_genes=uniprot_to_genes_all,
                pdb_to_genes=pdb_to_genes_all,
            )

            genes_this_chunk = genes_all.intersection(chunk_gene_set)

            if not genes_this_chunk:
                continue

            matched_structures += 1

            rows_new = process_structure_file(
                path=path,
                genes=genes_this_chunk,
                uniprot_acc=acc,
                pdb_id=pdb_id,
                run_fpocket=args.run_fpocket,
                tmpdir=tmpdir,
            )

            rows.extend(rows_new)
            matched_rows += len(rows_new)

            if matched_structures % 50 == 0:
                log(
                    f"[PROCESS] scanned={scanned:,}/{len(structure_files):,} "
                    f"matched_structures={matched_structures:,} rows={matched_rows:,}"
                )

    long_df = pd.DataFrame(rows)

    if not long_df.empty:
        long_df = long_df.drop_duplicates()

    long_df.to_csv(long_out, index=False)

    features = aggregate_gene_features_for_chunk(
        chunk_genes=chunk_genes,
        long_df=long_df,
    )
    features.to_csv(gene_out, index=False)

    meta = {
        "created_at": now_iso(),
        "job_index": job_index,
        "chunk_size": args.chunk_size,
        "start": start,
        "end": end,
        "n_total_genes": len(all_genes),
        "n_chunk_genes": len(chunk_genes),
        "first_gene": chunk_genes[0],
        "last_gene": chunk_genes[-1],
        "hgnc": str(Path(args.hgnc)),
        "feature4_dir": str(feature4_dir),
        "outdir": str(outdir.resolve()),
        "extra_structure_dirs": [str(x) for x in extra_dirs],
        "run_fpocket": args.run_fpocket,
        "fpocket_available": fpocket_available(),
        "limit_structures": args.limit_structures,
        "structure_files_scanned": scanned,
        "matched_structures": matched_structures,
        "long_rows": int(long_df.shape[0]),
        "gene_feature_rows": int(features.shape[0]),
        "outputs": {
            "long_chunk": str(long_out),
            "gene_features_chunk": str(gene_out),
        },
        "status": "done",
    }

    meta_out.write_text(json.dumps(meta, indent=2) + "\n")

    log("=" * 100)
    log("[DONE CHUNK]")
    log(f"[JOB INDEX]          {job_index}")
    log(f"[GENES]              {len(chunk_genes):,}")
    log(f"[MATCHED STRUCTURES] {matched_structures:,}")
    log(f"[LONG ROWS]          {long_df.shape[0]:,}")
    log(f"[GENE ROWS]          {features.shape[0]:,}")
    log(f"[SAVED LONG]         {long_out}")
    log(f"[SAVED GENE]         {gene_out}")
    log(f"[SAVED META]         {meta_out}")
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