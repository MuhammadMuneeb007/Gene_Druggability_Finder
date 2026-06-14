 

"""
Feature3_Pathway.py

Feature 3: Reactome pathway-membership features for gene-level druggability modelling.

Purpose
-------
Download Reactome pathway mapping files, map them to the HGNC protein-coding gene
universe, and generate one gene-level pathway feature table.

Why download files instead of one API request per gene?
------------------------------------------------------
Reactome pathway membership is available as downloadable mapping files. Downloading
the files once is faster, more reproducible, easier to cache, and less likely to
hit API limits than sending thousands of per-gene REST requests.

Leakage policy
--------------
This script uses only curated pathway membership and pathway hierarchy information.
It does NOT use known drug targets, compound-target annotations, drug response,
ChEMBL labels, DrugBank labels, approved-drug status, or any target-development
status.

Default input
-------------
    databases/HGNC/hgnc_complete_set.txt

Default output directory
------------------------
    feature3_pathway/

Main outputs
------------
    feature3_pathway/downloads/
    feature3_pathway/processed/feature3_reactome_gene_features.csv
    feature3_pathway/processed/feature3_reactome_hgnc_merged.csv
    feature3_pathway/processed/feature3_reactome_gene_pathway_long.csv
    feature3_pathway/feature3_reactome_summary.txt
    feature3_pathway/feature3_reactome_run_metadata.json

Run
---
    python Feature3_Pathway.py

Use existing downloads only:
    python Feature3_Pathway.py --no-download

Use all HGNC genes, not only protein-coding:
    python Feature3_Pathway.py --all-hgnc-genes

Force re-download:
    python Feature3_Pathway.py --force-download
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
import time
from collections import defaultdict, deque
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import requests


# =============================================================================
# DEFAULTS
# =============================================================================

DEFAULT_HGNC_PATH = Path("databases") / "HGNC" / "hgnc_complete_set.txt"
DEFAULT_OUTDIR = Path("feature3_pathway")

REACTOME_BASE = "https://reactome.org/download/current"

# Core files. UniProt2Reactome is the main mapping file because HGNC has UniProt IDs.
# Ensembl2Reactome is used as a backup/second mapping route because HGNC has Ensembl IDs.
REACTOME_FILES = {
    "UniProt2Reactome.txt": f"{REACTOME_BASE}/UniProt2Reactome.txt",
    "Ensembl2Reactome.txt": f"{REACTOME_BASE}/Ensembl2Reactome.txt",
    "ReactomePathways.txt": f"{REACTOME_BASE}/ReactomePathways.txt",
    "ReactomePathwaysRelation.txt": f"{REACTOME_BASE}/ReactomePathwaysRelation.txt",
    "HumanDiseasePathways.txt": f"{REACTOME_BASE}/HumanDiseasePathways.txt",
}

# Optional but useful for pathway-type interpretation.
OPTIONAL_REACTOME_FILES = {
    "Pathways2GoTerms_human.txt": f"{REACTOME_BASE}/Pathways2GoTerms_human.txt",
}

HUMAN_SPECIES_NAME = "Homo sapiens"


# =============================================================================
# HELPERS
# =============================================================================

def now_iso() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str) -> None:
    print(msg, flush=True)


def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def normalize_gene_symbol(x: object) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


def split_pipe_values(x: object) -> List[str]:
    if pd.isna(x):
        return []
    s = str(x).strip()
    if not s or s.lower() in {"nan", "none"}:
        return []
    return [v.strip() for v in s.split("|") if v.strip()]


def clean_id(x: object) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip()


def safe_filename_key(pathway_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", pathway_id)


def download_file(url: str, outpath: Path, force: bool = False, retries: int = 5) -> bool:
    if outpath.exists() and outpath.stat().st_size > 0 and not force:
        log(f"[SKIP] Existing: {outpath} ({outpath.stat().st_size / 1024 / 1024:.2f} MB)")
        return True

    tmp = outpath.with_suffix(outpath.suffix + ".tmp")

    for attempt in range(1, retries + 1):
        try:
            log("=" * 100)
            log(f"[DOWNLOAD] {outpath.name}")
            log(f"[URL]      {url}")
            log(f"[SAVE TO]  {outpath}")

            with requests.get(url, stream=True, timeout=300) as r:
                if r.status_code != 200:
                    raise RuntimeError(f"HTTP {r.status_code}: {r.text[:500]}")

                total = int(r.headers.get("content-length", 0))
                downloaded = 0

                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if not chunk:
                            continue
                        f.write(chunk)
                        downloaded += len(chunk)

                        if total > 0:
                            pct = downloaded / total * 100
                            print(
                                f"\r       {downloaded / 1024 / 1024:.2f} MB / "
                                f"{total / 1024 / 1024:.2f} MB ({pct:.1f}%)",
                                end="",
                                flush=True,
                            )
                        else:
                            print(
                                f"\r       {downloaded / 1024 / 1024:.2f} MB",
                                end="",
                                flush=True,
                            )
                print()

            tmp.replace(outpath)

            if outpath.stat().st_size == 0:
                raise RuntimeError("Downloaded file is empty")

            log(f"[DONE] {outpath.name} ({outpath.stat().st_size / 1024 / 1024:.2f} MB)")
            return True

        except Exception as e:
            log(f"[FAILED] Attempt {attempt}/{retries}: {outpath.name}")
            log(f"[ERROR]  {e}")

            if tmp.exists():
                tmp.unlink()

            if attempt < retries:
                wait = 5 * attempt
                log(f"[RETRY] Waiting {wait} seconds...")
                time.sleep(wait)

    return False


# =============================================================================
# DOWNLOADS
# =============================================================================

def download_reactome_files(downloads_dir: Path, include_optional: bool, force: bool) -> Dict[str, str]:
    mkdir(downloads_dir)

    files = dict(REACTOME_FILES)
    if include_optional:
        files.update(OPTIONAL_REACTOME_FILES)

    metadata = {
        "created_at": now_iso(),
        "reactome_base": REACTOME_BASE,
        "downloads_dir": str(downloads_dir.resolve()),
        "files": {},
        "failed": [],
    }

    for fname, url in files.items():
        outpath = downloads_dir / fname
        ok = download_file(url, outpath, force=force)
        if ok:
            metadata["files"][fname] = {
                "url": url,
                "path": str(outpath),
                "size_bytes": outpath.stat().st_size,
            }
        else:
            metadata["failed"].append(fname)

    with open(downloads_dir / "feature3_reactome_download_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    if metadata["failed"]:
        log("=" * 100)
        log("[WARNING] Some Reactome files failed to download:")
        for f in metadata["failed"]:
            log(f"  - {f}")

    return {fname: str(downloads_dir / fname) for fname in files}


# =============================================================================
# HGNC
# =============================================================================

def load_hgnc(hgnc_path: Path, protein_coding_only: bool = True) -> pd.DataFrame:
    if not hgnc_path.exists():
        raise FileNotFoundError(
            f"HGNC file not found: {hgnc_path}\n"
            "Expected default path: databases/HGNC/hgnc_complete_set.txt"
        )

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
        else:
            log("[WARNING] Could not find locus_group/locus_type. Keeping all HGNC genes.")

    hgnc["gene_symbol"] = hgnc["symbol"].map(normalize_gene_symbol)

    if "entrez_id" in hgnc.columns:
        hgnc["entrez_id_str"] = hgnc["entrez_id"].astype(str).str.strip().replace({"nan": "", "None": ""})
    else:
        hgnc["entrez_id_str"] = ""

    if "ensembl_gene_id" in hgnc.columns:
        hgnc["ensembl_gene_id_str"] = hgnc["ensembl_gene_id"].astype(str).str.strip().replace({"nan": "", "None": ""})
    else:
        hgnc["ensembl_gene_id_str"] = ""

    if "uniprot_ids" in hgnc.columns:
        hgnc["uniprot_ids_str"] = hgnc["uniprot_ids"].astype(str).str.strip().replace({"nan": "", "None": ""})
    else:
        hgnc["uniprot_ids_str"] = ""

    hgnc = hgnc.drop_duplicates(subset=["gene_symbol"], keep="first").copy()

    log(f"[HGNC ROWS] {before} -> {len(hgnc)}")
    log(f"[PROTEIN CODING ONLY] {protein_coding_only}")

    return hgnc


def build_identifier_maps(hgnc: pd.DataFrame) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Build UniProt->gene_symbol and Ensembl->gene_symbol maps from HGNC.
    """
    uniprot_to_gene = {}
    ensembl_to_gene = {}

    for _, row in hgnc.iterrows():
        gene = row["gene_symbol"]

        for up in split_pipe_values(row.get("uniprot_ids_str", "")):
            if up and up not in uniprot_to_gene:
                uniprot_to_gene[up] = gene

        ens = clean_id(row.get("ensembl_gene_id_str", ""))
        if ens:
            # Remove version suffix if present.
            ens_base = ens.split(".")[0]
            if ens_base and ens_base not in ensembl_to_gene:
                ensembl_to_gene[ens_base] = gene

    log("=" * 100)
    log("[IDENTIFIER MAPS]")
    log(f"[UNIPROT IDS] {len(uniprot_to_gene)}")
    log(f"[ENSEMBL IDS] {len(ensembl_to_gene)}")

    return uniprot_to_gene, ensembl_to_gene


# =============================================================================
# REACTOME PARSING
# =============================================================================

def read_reactome_pathways(path: Path) -> pd.DataFrame:
    """
    ReactomePathways.txt:
      pathway_id    pathway_name    species
    """
    if not path.exists():
        raise FileNotFoundError(path)

    log(f"[READ] {path}")
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["pathway_id", "pathway_name", "species"],
        dtype=str,
        low_memory=False,
    )

    df = df[df["species"].astype(str).eq(HUMAN_SPECIES_NAME)].copy()
    df["pathway_id"] = df["pathway_id"].astype(str).str.strip()
    df["pathway_name"] = df["pathway_name"].astype(str).str.strip()
    df = df.drop_duplicates(subset=["pathway_id"], keep="first")

    log(f"[HUMAN PATHWAYS] {len(df)}")
    return df


def read_pathway_relations(path: Path, valid_pathway_ids: Set[str]) -> pd.DataFrame:
    """
    ReactomePathwaysRelation.txt:
      parent_pathway_id    child_pathway_id
    """
    if not path.exists():
        log(f"[WARNING] Missing relation file: {path}")
        return pd.DataFrame(columns=["parent_pathway_id", "child_pathway_id"])

    log(f"[READ] {path}")
    rel = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["parent_pathway_id", "child_pathway_id"],
        dtype=str,
        low_memory=False,
    )
    rel["parent_pathway_id"] = rel["parent_pathway_id"].astype(str).str.strip()
    rel["child_pathway_id"] = rel["child_pathway_id"].astype(str).str.strip()

    # Keep human pathway relationships by filtering to known human pathway IDs.
    rel = rel[
        rel["parent_pathway_id"].isin(valid_pathway_ids)
        & rel["child_pathway_id"].isin(valid_pathway_ids)
    ].drop_duplicates()

    log(f"[HUMAN RELATIONS] {len(rel)}")
    return rel


def build_pathway_hierarchy_features(
    pathways: pd.DataFrame,
    relations: pd.DataFrame,
) -> pd.DataFrame:
    """
    Compute pathway depth, number of children, and top-level pathway assignment.
    """
    ids = set(pathways["pathway_id"].tolist())

    children = defaultdict(list)
    parents = defaultdict(list)

    for _, row in relations.iterrows():
        p = row["parent_pathway_id"]
        c = row["child_pathway_id"]
        children[p].append(c)
        parents[c].append(p)

    top_level = sorted([pid for pid in ids if pid not in parents])
    top_level_set = set(top_level)

    # BFS from top-level nodes to assign minimum depth and root.
    depth = {pid: np.nan for pid in ids}
    root_id = {pid: "" for pid in ids}

    q = deque()
    for root in top_level:
        depth[root] = 0
        root_id[root] = root
        q.append(root)

    while q:
        node = q.popleft()
        for child in children.get(node, []):
            new_depth = int(depth[node]) + 1 if not pd.isna(depth[node]) else 1
            if pd.isna(depth.get(child, np.nan)) or new_depth < depth[child]:
                depth[child] = new_depth
                root_id[child] = root_id[node] or node
                q.append(child)

    out = pathways.copy()
    out["reactome_pathway_depth"] = out["pathway_id"].map(depth)
    out["reactome_pathway_child_count"] = out["pathway_id"].map(lambda x: len(children.get(x, [])))
    out["reactome_pathway_parent_count"] = out["pathway_id"].map(lambda x: len(parents.get(x, [])))
    out["reactome_is_top_level_pathway"] = out["pathway_id"].isin(top_level_set).astype(int)
    out["reactome_top_level_pathway_id"] = out["pathway_id"].map(root_id).fillna("")

    name_map = dict(zip(pathways["pathway_id"], pathways["pathway_name"]))
    out["reactome_top_level_pathway_name"] = out["reactome_top_level_pathway_id"].map(name_map).fillna("")

    return out


def read_human_disease_pathways(path: Path) -> Set[str]:
    """
    HumanDiseasePathways.txt format can vary slightly. We treat any R-HSA-* token
    in the file as a disease pathway ID.
    """
    if not path.exists():
        log(f"[WARNING] Missing disease pathway file: {path}")
        return set()

    text = path.read_text(errors="ignore")
    ids = set(re.findall(r"R-HSA-\d+", text))
    log(f"[DISEASE PATHWAY IDS] {len(ids)}")
    return ids


def read_identifier_to_reactome(
    path: Path,
    id_to_gene: Dict[str, str],
    source: str,
    valid_pathway_ids: Set[str],
) -> pd.DataFrame:
    """
    Generic parser for UniProt2Reactome.txt and Ensembl2Reactome.txt.

    Expected columns in current Reactome mapping files commonly include:
      identifier, pathway_id, url, pathway_name, evidence, species

    Some files have extra or fewer columns. This parser assigns the first six
    columns and ignores extras.
    """
    if not path.exists():
        log(f"[WARNING] Missing mapping file: {path}")
        return pd.DataFrame()

    log(f"[READ MAPPING] {path}")

    raw = pd.read_csv(
        path,
        sep="\t",
        header=None,
        dtype=str,
        low_memory=False,
    )

    if raw.shape[1] < 2:
        raise RuntimeError(f"Reactome mapping file has too few columns: {path}")

    # Assign first six columns by known Reactome convention.
    col_names = ["identifier", "pathway_id", "url", "pathway_name", "evidence", "species"]
    rename = {i: col_names[i] for i in range(min(raw.shape[1], len(col_names)))}
    raw = raw.rename(columns=rename)

    for c in col_names:
        if c not in raw.columns:
            raw[c] = ""

    raw["identifier"] = raw["identifier"].astype(str).str.strip()
    raw["pathway_id"] = raw["pathway_id"].astype(str).str.strip()
    raw["species"] = raw["species"].astype(str).str.strip()

    # Keep human records.
    raw = raw[raw["species"].eq(HUMAN_SPECIES_NAME)].copy()

    # Keep only known human pathway IDs.
    raw = raw[raw["pathway_id"].isin(valid_pathway_ids)].copy()

    raw["gene_symbol"] = raw["identifier"].map(id_to_gene)
    raw = raw[raw["gene_symbol"].notna()].copy()
    raw["gene_symbol"] = raw["gene_symbol"].map(normalize_gene_symbol)

    raw["mapping_source"] = source

    keep_cols = [
        "gene_symbol",
        "identifier",
        "pathway_id",
        "pathway_name",
        "evidence",
        "species",
        "mapping_source",
    ]

    out = raw[keep_cols].drop_duplicates().copy()

    log(f"[{source}] mapped gene-pathway rows: {len(out)}")
    log(f"[{source}] mapped genes: {out['gene_symbol'].nunique() if len(out) else 0}")

    return out


# =============================================================================
# FEATURE GENERATION
# =============================================================================

def classify_top_level_pathway_name(name: object) -> str:
    """
    Simplified broad category based on Reactome top-level pathway name.
    """
    if pd.isna(name):
        return "other"

    s = str(name).lower()

    if "immune" in s:
        return "immune_system"
    if "signal transduction" in s or "signaling" in s or "signalling" in s:
        return "signal_transduction"
    if "metabolism" in s:
        return "metabolism"
    if "gene expression" in s or "transcription" in s or "translation" in s:
        return "gene_expression"
    if "dna repair" in s or "cell cycle" in s or "cellular responses to stimuli" in s:
        return "cell_cycle_dna_response"
    if "transport" in s:
        return "transport"
    if "disease" in s:
        return "disease"
    if "developmental biology" in s:
        return "developmental_biology"
    if "hemostasis" in s:
        return "hemostasis"
    if "neuronal" in s or "nervous" in s:
        return "neuronal_system"
    if "extracellular matrix" in s or "cell-cell communication" in s:
        return "cell_communication_ecm"

    return "other"


def build_gene_pathway_long(
    hgnc: pd.DataFrame,
    downloads_dir: Path,
    processed_dir: Path,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Create long gene-pathway table using UniProt and Ensembl mappings.
    """
    pathways = read_reactome_pathways(downloads_dir / "ReactomePathways.txt")
    valid_ids = set(pathways["pathway_id"].tolist())

    relations = read_pathway_relations(downloads_dir / "ReactomePathwaysRelation.txt", valid_ids)
    pathway_meta = build_pathway_hierarchy_features(pathways, relations)

    disease_ids = read_human_disease_pathways(downloads_dir / "HumanDiseasePathways.txt")
    pathway_meta["reactome_is_disease_pathway"] = pathway_meta["pathway_id"].isin(disease_ids).astype(int)
    pathway_meta["reactome_broad_category"] = pathway_meta["reactome_top_level_pathway_name"].map(classify_top_level_pathway_name)

    uniprot_to_gene, ensembl_to_gene = build_identifier_maps(hgnc)

    up_map = read_identifier_to_reactome(
        downloads_dir / "UniProt2Reactome.txt",
        id_to_gene=uniprot_to_gene,
        source="UniProt",
        valid_pathway_ids=valid_ids,
    )

    ens_map = read_identifier_to_reactome(
        downloads_dir / "Ensembl2Reactome.txt",
        id_to_gene=ensembl_to_gene,
        source="Ensembl",
        valid_pathway_ids=valid_ids,
    )

    combined = pd.concat([up_map, ens_map], ignore_index=True) if len(up_map) or len(ens_map) else pd.DataFrame()

    if combined.empty:
        raise RuntimeError("No HGNC genes mapped to Reactome pathways.")

    # Merge metadata.
    combined = combined.merge(
        pathway_meta[
            [
                "pathway_id",
                "pathway_name",
                "reactome_pathway_depth",
                "reactome_pathway_child_count",
                "reactome_pathway_parent_count",
                "reactome_is_top_level_pathway",
                "reactome_top_level_pathway_id",
                "reactome_top_level_pathway_name",
                "reactome_is_disease_pathway",
                "reactome_broad_category",
            ]
        ],
        on="pathway_id",
        how="left",
        suffixes=("", "_meta"),
    )

    # Prefer pathway_name from metadata if available.
    if "pathway_name_meta" in combined.columns:
        combined["pathway_name"] = combined["pathway_name_meta"].fillna(combined["pathway_name"])
        combined = combined.drop(columns=["pathway_name_meta"])

    # One row per gene-pathway pair, but preserve source evidence as combined strings.
    grouped = (
        combined
        .groupby(["gene_symbol", "pathway_id"], as_index=False)
        .agg({
            "identifier": lambda x: ";".join(sorted(set([str(v) for v in x if pd.notna(v)]))[:20]),
            "pathway_name": "first",
            "evidence": lambda x: ";".join(sorted(set([str(v) for v in x if pd.notna(v) and str(v) != ""]))[:20]),
            "species": "first",
            "mapping_source": lambda x: ";".join(sorted(set([str(v) for v in x if pd.notna(v)]))),
            "reactome_pathway_depth": "first",
            "reactome_pathway_child_count": "first",
            "reactome_pathway_parent_count": "first",
            "reactome_is_top_level_pathway": "first",
            "reactome_top_level_pathway_id": "first",
            "reactome_top_level_pathway_name": "first",
            "reactome_is_disease_pathway": "first",
            "reactome_broad_category": "first",
        })
    )

    long_path = processed_dir / "feature3_reactome_gene_pathway_long.csv"
    grouped.to_csv(long_path, index=False)

    pathway_meta_path = processed_dir / "feature3_reactome_pathway_metadata.csv"
    pathway_meta.to_csv(pathway_meta_path, index=False)

    log("=" * 100)
    log("[SAVED LONG GENE-PATHWAY TABLE]")
    log(f"[PATH] {long_path}")
    log(f"[ROWS] {len(grouped)}")
    log(f"[GENES] {grouped['gene_symbol'].nunique()}")
    log(f"[PATHWAYS] {grouped['pathway_id'].nunique()}")
    log(f"[SAVED PATHWAY METADATA] {pathway_meta_path}")

    return grouped, pathway_meta


def build_gene_features(
    hgnc: pd.DataFrame,
    gene_pathway: pd.DataFrame,
    processed_dir: Path,
) -> pd.DataFrame:
    """
    Aggregate long gene-pathway table to one row per HGNC gene.
    """
    hgnc_genes = pd.DataFrame({"gene_symbol": sorted(hgnc["gene_symbol"].unique().tolist())})

    # Pathway size is number of genes mapped to each pathway.
    pathway_size = (
        gene_pathway
        .groupby("pathway_id")["gene_symbol"]
        .nunique()
        .rename("reactome_pathway_gene_count")
        .reset_index()
    )

    gp = gene_pathway.merge(pathway_size, on="pathway_id", how="left")
    gp["reactome_inverse_pathway_size"] = 1.0 / gp["reactome_pathway_gene_count"].replace(0, np.nan)

    # Broad category one-hot counts.
    category_counts = (
        gp
        .pivot_table(
            index="gene_symbol",
            columns="reactome_broad_category",
            values="pathway_id",
            aggfunc="nunique",
            fill_value=0,
        )
        .reset_index()
    )
    category_counts.columns = [
        "gene_symbol" if c == "gene_symbol" else f"feature3_reactome_category_count_{c}"
        for c in category_counts.columns
    ]

    records = []

    for gene, g in gp.groupby("gene_symbol"):
        pathway_ids = sorted(set(g["pathway_id"].dropna().astype(str).tolist()))
        pathway_names = sorted(set(g["pathway_name"].dropna().astype(str).tolist()))
        top_ids = sorted(set(g["reactome_top_level_pathway_id"].dropna().astype(str).tolist()))
        top_names = sorted(set(g["reactome_top_level_pathway_name"].dropna().astype(str).tolist()))

        sizes = pd.to_numeric(g["reactome_pathway_gene_count"], errors="coerce")
        inv_sizes = pd.to_numeric(g["reactome_inverse_pathway_size"], errors="coerce")
        depths = pd.to_numeric(g["reactome_pathway_depth"], errors="coerce")
        child_counts = pd.to_numeric(g["reactome_pathway_child_count"], errors="coerce")

        sources = set()
        for v in g["mapping_source"].dropna().astype(str):
            for part in v.split(";"):
                if part:
                    sources.add(part)

        rec = {
            "gene_symbol": gene,

            # Basic pathway membership.
            "feature3_reactome_has_pathway": 1,
            "feature3_reactome_pathway_count": len(pathway_ids),
            "feature3_reactome_top_level_pathway_count": len([x for x in top_ids if x]),
            "feature3_reactome_mapping_source_count": len(sources),
            "feature3_reactome_mapped_by_uniprot": int("UniProt" in sources),
            "feature3_reactome_mapped_by_ensembl": int("Ensembl" in sources),

            # Disease/pathway broad flags.
            "feature3_reactome_disease_pathway_count": int(g.loc[g["reactome_is_disease_pathway"].fillna(0).astype(int) == 1, "pathway_id"].nunique()),
            "feature3_reactome_has_disease_pathway": int((g["reactome_is_disease_pathway"].fillna(0).astype(int) == 1).any()),

            # Pathway hierarchy.
            "feature3_reactome_mean_pathway_depth": float(depths.mean(skipna=True)) if depths.notna().any() else np.nan,
            "feature3_reactome_max_pathway_depth": float(depths.max(skipna=True)) if depths.notna().any() else np.nan,
            "feature3_reactome_mean_pathway_child_count": float(child_counts.mean(skipna=True)) if child_counts.notna().any() else np.nan,
            "feature3_reactome_max_pathway_child_count": float(child_counts.max(skipna=True)) if child_counts.notna().any() else np.nan,

            # Pathway size/specificity.
            "feature3_reactome_mean_pathway_size": float(sizes.mean(skipna=True)) if sizes.notna().any() else np.nan,
            "feature3_reactome_median_pathway_size": float(sizes.median(skipna=True)) if sizes.notna().any() else np.nan,
            "feature3_reactome_min_pathway_size": float(sizes.min(skipna=True)) if sizes.notna().any() else np.nan,
            "feature3_reactome_max_pathway_size": float(sizes.max(skipna=True)) if sizes.notna().any() else np.nan,
            "feature3_reactome_mean_inverse_pathway_size": float(inv_sizes.mean(skipna=True)) if inv_sizes.notna().any() else np.nan,
            "feature3_reactome_sum_inverse_pathway_size": float(inv_sizes.sum(skipna=True)) if inv_sizes.notna().any() else np.nan,

            # Audit columns.
            "feature3_reactome_top_level_pathway_ids": ";".join([x for x in top_ids if x][:50]),
            "feature3_reactome_top_level_pathway_names": ";".join([x for x in top_names if x][:50]),
            "feature3_reactome_pathway_ids": ";".join(pathway_ids[:100]),
            "feature3_reactome_pathway_names": ";".join(pathway_names[:100]),
        }

        records.append(rec)

    features = pd.DataFrame(records)

    if features.empty:
        features = hgnc_genes.copy()

    features = hgnc_genes.merge(features, on="gene_symbol", how="left")

    # Add category counts.
    features = features.merge(category_counts, on="gene_symbol", how="left")

    # Fill zero for count/flag columns.
    for c in features.columns:
        if c.startswith("feature3_reactome_") and (
            "_count" in c
            or c.startswith("feature3_reactome_has_")
            or c.startswith("feature3_reactome_mapped_by_")
        ):
            if features[c].dtype == object:
                # Avoid filling audit text accidentally.
                continue
            features[c] = features[c].fillna(0)

    numeric_zero_cols = [
        "feature3_reactome_has_pathway",
        "feature3_reactome_pathway_count",
        "feature3_reactome_top_level_pathway_count",
        "feature3_reactome_mapping_source_count",
        "feature3_reactome_mapped_by_uniprot",
        "feature3_reactome_mapped_by_ensembl",
        "feature3_reactome_disease_pathway_count",
        "feature3_reactome_has_disease_pathway",
    ]
    for c in numeric_zero_cols:
        if c in features.columns:
            features[c] = features[c].fillna(0).astype(int)

    cat_cols = [c for c in features.columns if c.startswith("feature3_reactome_category_count_")]
    for c in cat_cols:
        features[c] = features[c].fillna(0).astype(int)

    text_cols = [
        "feature3_reactome_top_level_pathway_ids",
        "feature3_reactome_top_level_pathway_names",
        "feature3_reactome_pathway_ids",
        "feature3_reactome_pathway_names",
    ]
    for c in text_cols:
        if c in features.columns:
            features[c] = features[c].fillna("")

    outpath = processed_dir / "feature3_reactome_gene_features.csv"
    features.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED REACTOME GENE FEATURES]")
    log(f"[PATH]  {outpath}")
    log(f"[SHAPE] {features.shape[0]} genes x {features.shape[1]} columns")

    return features


def merge_with_hgnc(hgnc: pd.DataFrame, features: pd.DataFrame, processed_dir: Path) -> pd.DataFrame:
    log("=" * 100)
    log("[MERGE HGNC + REACTOME FEATURES]")

    merged = hgnc.merge(features, on="gene_symbol", how="left")

    feature_cols = [c for c in features.columns if c != "gene_symbol"]
    merged["feature3_reactome_has_any_feature"] = merged[feature_cols].notna().any(axis=1).astype(int)
    merged["feature3_reactome_n_nonmissing_features"] = merged[feature_cols].notna().sum(axis=1)

    outpath = processed_dir / "feature3_reactome_hgnc_merged.csv"
    merged.to_csv(outpath, index=False)

    coverage = merged["feature3_reactome_has_pathway"].mean() if "feature3_reactome_has_pathway" in merged.columns else np.nan

    log(f"[SAVED] {outpath}")
    log(f"[SHAPE] {merged.shape[0]} genes x {merged.shape[1]} columns")
    log(f"[PATHWAY COVERAGE] {coverage:.3f}")

    return merged


def write_summary_report(
    outdir: Path,
    hgnc: pd.DataFrame,
    gene_pathway: pd.DataFrame,
    features: pd.DataFrame,
    merged: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    report_path = outdir / "feature3_reactome_summary.txt"

    lines = []
    lines.append("Feature 3 Reactome pathway-membership features")
    lines.append("=" * 80)
    lines.append(f"Created at: {now_iso()}")
    lines.append(f"Reactome base URL: {REACTOME_BASE}")
    lines.append(f"Output directory: {outdir.resolve()}")
    lines.append("")
    lines.append("Leakage rule:")
    lines.append("  Included: Reactome curated pathway membership and pathway hierarchy.")
    lines.append("  Excluded: known drug targets, compound-target metadata, drug-response data, ChEMBL/DrugBank labels.")
    lines.append("")
    lines.append(f"HGNC rows: {hgnc.shape[0]}")
    lines.append(f"Gene-pathway rows: {gene_pathway.shape[0]}")
    lines.append(f"Mapped genes in gene-pathway table: {gene_pathway['gene_symbol'].nunique()}")
    lines.append(f"Unique Reactome pathways: {gene_pathway['pathway_id'].nunique()}")
    lines.append(f"Feature rows: {features.shape[0]}")
    lines.append(f"Feature columns: {features.shape[1]}")
    lines.append(f"Merged rows: {merged.shape[0]}")
    lines.append(f"Merged columns: {merged.shape[1]}")

    if "feature3_reactome_has_pathway" in merged.columns:
        lines.append(f"Pathway coverage: {merged['feature3_reactome_has_pathway'].mean():.4f}")
    if "feature3_reactome_pathway_count" in merged.columns:
        lines.append(f"Median pathway count: {merged['feature3_reactome_pathway_count'].median():.2f}")
        lines.append(f"Mean pathway count: {merged['feature3_reactome_pathway_count'].mean():.2f}")
        lines.append(f"Max pathway count: {merged['feature3_reactome_pathway_count'].max():.2f}")
    if "feature3_reactome_has_disease_pathway" in merged.columns:
        lines.append(f"Disease-pathway gene fraction: {merged['feature3_reactome_has_disease_pathway'].mean():.4f}")

    report_path.write_text("\n".join(lines) + "\n")
    log(f"[REPORT] {report_path}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build Reactome pathway-membership features for HGNC genes."
    )
    parser.add_argument("--hgnc", default=str(DEFAULT_HGNC_PATH), help="Path to hgnc_complete_set.txt.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory.")
    parser.add_argument("--all-hgnc-genes", action="store_true", help="Use all HGNC rows instead of protein-coding only.")
    parser.add_argument("--no-download", action="store_true", help="Use existing files in downloads folder.")
    parser.add_argument("--force-download", action="store_true", help="Re-download Reactome files.")
    parser.add_argument("--include-optional", action="store_true", help="Download optional Reactome helper files.")
    args = parser.parse_args()

    outdir = mkdir(Path(args.outdir))
    downloads_dir = mkdir(outdir / "downloads")
    processed_dir = mkdir(outdir / "processed")

    log("=" * 100)
    log("FEATURE 3: REACTOME PATHWAY MEMBERSHIP")
    log("=" * 100)
    log(f"[HGNC]             {args.hgnc}")
    log(f"[OUTDIR]           {outdir.resolve()}")
    log(f"[DOWNLOADS]        {downloads_dir.resolve()}")
    log(f"[PROCESSED]        {processed_dir.resolve()}")
    log(f"[NO DOWNLOAD]      {args.no_download}")
    log(f"[FORCE DOWNLOAD]   {args.force_download}")
    log("=" * 100)

    if not args.no_download:
        download_reactome_files(
            downloads_dir=downloads_dir,
            include_optional=args.include_optional,
            force=args.force_download,
        )

    hgnc = load_hgnc(
        hgnc_path=Path(args.hgnc),
        protein_coding_only=not args.all_hgnc_genes,
    )

    gene_pathway, pathway_meta = build_gene_pathway_long(
        hgnc=hgnc,
        downloads_dir=downloads_dir,
        processed_dir=processed_dir,
    )

    features = build_gene_features(
        hgnc=hgnc,
        gene_pathway=gene_pathway,
        processed_dir=processed_dir,
    )

    merged = merge_with_hgnc(
        hgnc=hgnc,
        features=features,
        processed_dir=processed_dir,
    )

    metadata = {
        "created_at": now_iso(),
        "reactome_base": REACTOME_BASE,
        "hgnc": str(Path(args.hgnc)),
        "outdir": str(outdir.resolve()),
        "protein_coding_only": not args.all_hgnc_genes,
        "leakage_policy": {
            "included": [
                "Reactome pathway membership",
                "Reactome pathway hierarchy",
                "Reactome broad top-level pathway categories",
                "Reactome disease pathway membership as pathway-category membership only",
            ],
            "excluded": [
                "known drug targets",
                "compound-target metadata",
                "drug-response data",
                "ChEMBL labels",
                "DrugBank labels",
                "approved-drug status",
            ],
        },
        "downloaded_files": list(REACTOME_FILES.keys()) + (list(OPTIONAL_REACTOME_FILES.keys()) if args.include_optional else []),
        "outputs": {
            "gene_pathway_long": str(processed_dir / "feature3_reactome_gene_pathway_long.csv"),
            "pathway_metadata": str(processed_dir / "feature3_reactome_pathway_metadata.csv"),
            "gene_features": str(processed_dir / "feature3_reactome_gene_features.csv"),
            "hgnc_merged": str(processed_dir / "feature3_reactome_hgnc_merged.csv"),
        },
    }

    with open(outdir / "feature3_reactome_run_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    write_summary_report(outdir, hgnc, gene_pathway, features, merged, args)

    log("=" * 100)
    log("[DONE]")
    log(f"[GENE-PATHWAY LONG] {processed_dir / 'feature3_reactome_gene_pathway_long.csv'}")
    log(f"[GENE FEATURES]     {processed_dir / 'feature3_reactome_gene_features.csv'}")
    log(f"[HGNC MERGED TABLE] {processed_dir / 'feature3_reactome_hgnc_merged.csv'}")
    log("=" * 100)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        log("\n[STOPPED] Interrupted by user.")
        sys.exit(130)
    except Exception as e:
        log("=" * 100)
        log("[ERROR]")
        log(str(e))
        log("=" * 100)
        raise
 