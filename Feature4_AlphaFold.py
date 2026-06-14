#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature4_Structure.py

Feature 4: PDB + AlphaFold structure downloader and gene-level structural features.

Design:
1. Start from local HGNC protein-coding genes.
2. Use HGNC UniProt IDs where available.
3. For each UniProt protein:
   - Query UniProt metadata.
   - Get linked experimental PDB IDs.
   - Download PDB structures first.
   - Always download AlphaFold structure and AlphaFold confidence JSON if available.
4. Save all structures locally for later pocket analysis.
5. Build:
   - protein/UniProt-level table
   - PDB long table
   - gene-level feature table
   - HGNC-merged table

Run:
    python Feature4_Structure.py

Fast test:
    python Feature4_Structure.py --limit-genes 100 --workers 8 --max-pdb-per-uniprot 3

Final recommended:
    python Feature4_Structure.py --workers 12 --max-pdb-per-uniprot 5

If you want all PDBs:
    python Feature4_Structure.py --workers 12 --max-pdb-per-uniprot 0

Outputs:
    feature4_structure/downloads/pdb/
    feature4_structure/downloads/alphafold_pdb/
    feature4_structure/downloads/alphafold_confidence/
    feature4_structure/downloads/metadata/
    feature4_structure/processed/feature4_structure_uniprot_table.csv
    feature4_structure/processed/feature4_structure_pdb_long.csv
    feature4_structure/processed/feature4_structure_gene_features.csv
    feature4_structure/processed/feature4_structure_hgnc_merged.csv
    feature4_structure/feature4_structure_summary.txt
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


# =============================================================================
# CONFIG
# =============================================================================

DEFAULT_HGNC_PATH = Path("databases") / "HGNC" / "hgnc_complete_set.txt"
DEFAULT_OUTDIR = Path("feature4_structure")

UNIPROT_ENTRY_API = "https://rest.uniprot.org/uniprotkb/{acc}.json"
UNIPROT_SEARCH_API = "https://rest.uniprot.org/uniprotkb/search"

RCSB_ENTRY_API = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_PDB_DOWNLOAD = "https://files.rcsb.org/download/{pdb_id}.pdb"

ALPHAFOLD_VERSIONS = ["v6", "v5", "v4", "v3", "v2"]
AF_PDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_{ver}.pdb"
AF_CONF_URL = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-confidence_{ver}.json"

HEADERS = {
    "User-Agent": "Feature4_Structure/1.0 academic structure feature downloader"
}

REQUEST_TIMEOUT = 45
MAX_RETRIES = 4
DEFAULT_WORKERS = 12


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


def safe_token(x: Any) -> str:
    s = "" if x is None else str(x).strip()
    if not s or s.lower() in {"nan", "none"}:
        return "NA"
    s = re.sub(r"[^A-Za-z0-9._+-]+", "_", s)
    return s.strip("._") or "NA"


def normalize_symbol(x: Any) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


def split_uniprot_ids(x: Any) -> List[str]:
    if pd.isna(x):
        return []
    s = str(x).strip()
    if not s or s.lower() in {"nan", "none"}:
        return []

    ids = []
    for part in re.split(r"[|;, ]+", s):
        part = part.strip()
        if part and part.lower() not in {"nan", "none"} and part not in ids:
            ids.append(part)
    return ids


def safe_float(x: Any) -> float:
    try:
        if x is None or pd.isna(x):
            return np.nan
        return float(x)
    except Exception:
        return np.nan


def make_session() -> requests.Session:
    session = requests.Session()
    retry = Retry(
        total=MAX_RETRIES,
        connect=MAX_RETRIES,
        read=MAX_RETRIES,
        backoff_factor=1.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "HEAD"],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=100, pool_maxsize=100)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update(HEADERS)
    return session


def read_json(path: Path) -> Optional[dict]:
    try:
        if path.exists() and path.stat().st_size > 0:
            return json.loads(path.read_text(errors="ignore"))
    except Exception:
        return None
    return None


def write_json(path: Path, data: dict) -> None:
    mkdir(path.parent)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2))
    tmp.replace(path)


def get_json_cached(
    session: requests.Session,
    url: str,
    cache_path: Path,
    force: bool = False,
) -> Optional[dict]:
    if cache_path.exists() and cache_path.stat().st_size > 0 and not force:
        return read_json(cache_path)

    try:
        r = session.get(url, timeout=REQUEST_TIMEOUT)
        if r.status_code != 200:
            return None
        data = r.json()
        write_json(cache_path, data)
        return data
    except Exception:
        return None


def url_exists(session: requests.Session, url: str) -> bool:
    try:
        r = session.head(url, timeout=15, allow_redirects=True)
        if r.status_code == 200:
            return True
        if r.status_code == 404:
            return False
    except Exception:
        pass

    try:
        r = session.get(url, timeout=15, stream=True, allow_redirects=True)
        ok = r.status_code == 200
        r.close()
        return ok
    except Exception:
        return False


def download_file(
    session: requests.Session,
    url: str,
    dest: Path,
    force: bool = False,
    min_bytes: int = 50,
) -> bool:
    if dest.exists() and dest.stat().st_size >= min_bytes and not force:
        return True

    mkdir(dest.parent)
    tmp = dest.with_suffix(dest.suffix + ".part")

    try:
        with session.get(url, timeout=REQUEST_TIMEOUT, stream=True, allow_redirects=True) as r:
            if r.status_code == 404:
                return False
            if r.status_code != 200:
                return False

            with open(tmp, "wb") as f:
                for chunk in r.iter_content(chunk_size=256 * 1024):
                    if chunk:
                        f.write(chunk)

        if tmp.exists() and tmp.stat().st_size >= min_bytes:
            tmp.replace(dest)
            return True

        if tmp.exists():
            tmp.unlink()
        return False

    except Exception:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass
        return False


# =============================================================================
# HGNC
# =============================================================================

def load_hgnc(hgnc_path: Path, protein_coding_only: bool = True) -> pd.DataFrame:
    if not hgnc_path.exists():
        raise FileNotFoundError(f"HGNC file not found: {hgnc_path}")

    log("=" * 100)
    log(f"[READ HGNC] {hgnc_path}")

    df = pd.read_csv(hgnc_path, sep="\t", dtype=str, low_memory=False)
    df.columns = [c.strip() for c in df.columns]

    if "symbol" not in df.columns:
        raise RuntimeError("HGNC file must contain a symbol column.")

    before = len(df)

    if protein_coding_only:
        if "locus_group" in df.columns:
            df = df[df["locus_group"].astype(str).str.lower().eq("protein-coding gene")].copy()
        elif "locus_type" in df.columns:
            df = df[df["locus_type"].astype(str).str.lower().str.contains("protein", na=False)].copy()

    df["gene_symbol"] = df["symbol"].map(normalize_symbol)

    if "uniprot_ids" not in df.columns:
        df["uniprot_ids"] = ""

    if "ensembl_gene_id" not in df.columns:
        df["ensembl_gene_id"] = ""

    if "entrez_id" not in df.columns:
        df["entrez_id"] = ""

    df = df.drop_duplicates(subset=["gene_symbol"], keep="first").copy()

    log(f"[HGNC ROWS] {before} -> {len(df)}")
    log(f"[PROTEIN CODING ONLY] {protein_coding_only}")
    log(f"[GENES WITH UNIPROT] {(df['uniprot_ids'].fillna('').astype(str).str.len() > 0).sum()}")

    return df


def build_tasks(hgnc: pd.DataFrame, limit_genes: int = 0) -> List[dict]:
    tasks = []

    sub = hgnc.copy()
    if limit_genes and limit_genes > 0:
        sub = sub.head(limit_genes).copy()

    for _, row in sub.iterrows():
        gene = row["gene_symbol"]
        accs = split_uniprot_ids(row.get("uniprot_ids", ""))

        if not accs:
            tasks.append({
                "gene_symbol": gene,
                "hgnc_id": row.get("hgnc_id", ""),
                "ensembl_gene_id": row.get("ensembl_gene_id", ""),
                "entrez_id": row.get("entrez_id", ""),
                "uniprot_accession": "",
                "is_primary_uniprot": 0,
            })
        else:
            for i, acc in enumerate(accs):
                tasks.append({
                    "gene_symbol": gene,
                    "hgnc_id": row.get("hgnc_id", ""),
                    "ensembl_gene_id": row.get("ensembl_gene_id", ""),
                    "entrez_id": row.get("entrez_id", ""),
                    "uniprot_accession": acc,
                    "is_primary_uniprot": 1 if i == 0 else 0,
                })

    return tasks


# =============================================================================
# UNIPROT
# =============================================================================

def search_uniprot_by_gene(
    session: requests.Session,
    gene_symbol: str,
    metadata_dir: Path,
    force: bool = False,
) -> Tuple[str, Optional[dict]]:
    cache = metadata_dir / "uniprot_search" / f"{safe_token(gene_symbol)}.json"

    if cache.exists() and cache.stat().st_size > 0 and not force:
        data = read_json(cache)
    else:
        params = {
            "query": f"gene_exact:{gene_symbol} AND organism_id:9606 AND reviewed:true",
            "format": "json",
            "fields": "accession,id,protein_name,gene_names,organism_name,organism_id,sequence,xref_pdb",
            "size": 1,
        }
        try:
            r = session.get(UNIPROT_SEARCH_API, params=params, timeout=REQUEST_TIMEOUT)
            if r.status_code != 200:
                return "", None
            data = r.json()
            write_json(cache, data)
        except Exception:
            return "", None

    results = data.get("results", []) if data else []
    if not results:
        return "", None

    acc = results[0].get("primaryAccession", "")
    return acc, results[0]


def fetch_uniprot_entry(
    session: requests.Session,
    acc: str,
    metadata_dir: Path,
    force: bool = False,
) -> Optional[dict]:
    if not acc:
        return None
    cache = metadata_dir / "uniprot" / f"{safe_token(acc)}.json"
    url = UNIPROT_ENTRY_API.format(acc=acc)
    return get_json_cached(session, url, cache, force=force)


def parse_uniprot(entry: Optional[dict]) -> Dict[str, Any]:
    if not entry:
        return {
            "uniprot_found": 0,
            "protein_name": "",
            "uniprot_id": "",
            "entry_type": "",
            "organism_name": "",
            "organism_id": "",
            "protein_length": np.nan,
            "pdb_ids": [],
        }

    protein_desc = entry.get("proteinDescription", {}) or {}
    rec_name = protein_desc.get("recommendedName", {}) or {}
    full_name = rec_name.get("fullName", {}) or {}

    protein_name = full_name.get("value", "")
    if not protein_name:
        subs = protein_desc.get("submissionNames", []) or []
        if subs:
            protein_name = (subs[0].get("fullName", {}) or {}).get("value", "")

    organism = entry.get("organism", {}) or {}
    seq = entry.get("sequence", {}) or {}

    pdb_ids = []
    for ref in entry.get("uniProtKBCrossReferences", []) or []:
        if ref.get("database") == "PDB" and ref.get("id"):
            pdb_ids.append(ref["id"].upper())

    return {
        "uniprot_found": 1,
        "protein_name": protein_name,
        "uniprot_id": entry.get("uniProtkbId", ""),
        "entry_type": entry.get("entryType", ""),
        "organism_name": organism.get("scientificName", ""),
        "organism_id": str(organism.get("taxonId", "")),
        "protein_length": seq.get("length", np.nan),
        "pdb_ids": sorted(set(pdb_ids)),
    }


# =============================================================================
# RCSB / PDB
# =============================================================================

def fetch_rcsb_entry(
    session: requests.Session,
    pdb_id: str,
    metadata_dir: Path,
    force: bool = False,
) -> Optional[dict]:
    pdb_id = pdb_id.upper()
    cache = metadata_dir / "rcsb_entry" / f"{pdb_id}.json"
    url = RCSB_ENTRY_API.format(pdb_id=pdb_id)
    return get_json_cached(session, url, cache, force=force)


def parse_rcsb_entry(pdb_id: str, data: Optional[dict]) -> Dict[str, Any]:
    if not data:
        return {
            "pdb_id": pdb_id,
            "pdb_metadata_found": 0,
            "pdb_title": "",
            "pdb_deposition_date": "",
            "pdb_release_date": "",
            "pdb_method": "",
            "pdb_resolution": np.nan,
            "pdb_chain_count": np.nan,
            "pdb_polymer_entity_count": np.nan,
            "pdb_nonpolymer_entity_count": np.nan,
            "pdb_atom_count": np.nan,
            "pdb_molecular_weight": np.nan,
        }

    methods = sorted({
        x.get("method", "")
        for x in data.get("exptl", []) or []
        if x.get("method")
    })

    entry_info = data.get("rcsb_entry_info", {}) or {}
    ids = data.get("rcsb_entry_container_identifiers", {}) or {}
    access = data.get("rcsb_accession_info", {}) or {}

    res = np.nan
    res_list = entry_info.get("resolution_combined") or []
    if res_list:
        res = safe_float(res_list[0])

    return {
        "pdb_id": pdb_id,
        "pdb_metadata_found": 1,
        "pdb_title": (data.get("struct", {}) or {}).get("title", ""),
        "pdb_deposition_date": access.get("deposit_date", ""),
        "pdb_release_date": access.get("initial_release_date", ""),
        "pdb_method": ";".join(methods),
        "pdb_resolution": res,
        "pdb_chain_count": entry_info.get("deposited_polymer_entity_instance_count", np.nan),
        "pdb_polymer_entity_count": len(ids.get("polymer_entity_ids", []) or []),
        "pdb_nonpolymer_entity_count": len(ids.get("non_polymer_entity_ids", []) or []),
        "pdb_atom_count": entry_info.get("deposited_atom_count", np.nan),
        "pdb_molecular_weight": entry_info.get("molecular_weight", np.nan),
    }


def download_pdb(
    session: requests.Session,
    pdb_id: str,
    gene_symbol: str,
    pdb_dir: Path,
    force: bool = False,
) -> str:
    pdb_id = pdb_id.upper()
    dest = pdb_dir / safe_token(gene_symbol) / f"{pdb_id}.pdb"
    url = RCSB_PDB_DOWNLOAD.format(pdb_id=pdb_id)
    ok = download_file(session, url, dest, force=force, min_bytes=100)
    return str(dest) if ok else ""


# =============================================================================
# ALPHAFOLD
# =============================================================================

def find_alphafold(
    session: requests.Session,
    acc: str,
) -> Tuple[str, str, str]:
    for ver in ALPHAFOLD_VERSIONS:
        pdb_url = AF_PDB_URL.format(acc=acc, ver=ver)
        if url_exists(session, pdb_url):
            conf_url = AF_CONF_URL.format(acc=acc, ver=ver)
            return ver, pdb_url, conf_url
    return "", "", ""


def download_alphafold(
    session: requests.Session,
    acc: str,
    gene_symbol: str,
    af_pdb_dir: Path,
    af_conf_dir: Path,
    force: bool = False,
) -> Dict[str, Any]:
    result = {
        "alphafold_available": 0,
        "alphafold_version": "",
        "alphafold_pdb_url": "",
        "alphafold_confidence_url": "",
        "alphafold_pdb_file": "",
        "alphafold_confidence_file": "",
        "alphafold_confidence_json_available": 0,
    }

    if not acc:
        return result

    ver, pdb_url, conf_url = find_alphafold(session, acc)
    if not ver:
        return result

    result["alphafold_available"] = 1
    result["alphafold_version"] = ver
    result["alphafold_pdb_url"] = pdb_url
    result["alphafold_confidence_url"] = conf_url

    pdb_dest = af_pdb_dir / safe_token(gene_symbol) / f"AF-{acc}-F1-model_{ver}.pdb"
    conf_dest = af_conf_dir / safe_token(gene_symbol) / f"AF-{acc}-F1-confidence_{ver}.json"

    if download_file(session, pdb_url, pdb_dest, force=force, min_bytes=100):
        result["alphafold_pdb_file"] = str(pdb_dest)

    if download_file(session, conf_url, conf_dest, force=force, min_bytes=10):
        result["alphafold_confidence_file"] = str(conf_dest)
        result["alphafold_confidence_json_available"] = 1

    return result


def parse_alphafold_confidence(conf_file: str) -> Dict[str, Any]:
    out = {
        "alphafold_residue_count": 0,
        "alphafold_mean_plddt": np.nan,
        "alphafold_median_plddt": np.nan,
        "alphafold_min_plddt": np.nan,
        "alphafold_max_plddt": np.nan,
        "alphafold_q05_plddt": np.nan,
        "alphafold_q25_plddt": np.nan,
        "alphafold_q75_plddt": np.nan,
        "alphafold_q95_plddt": np.nan,
        "alphafold_fraction_plddt_ge_90": np.nan,
        "alphafold_fraction_plddt_ge_70": np.nan,
        "alphafold_fraction_plddt_lt_50": np.nan,
        "alphafold_predicted_disorder_fraction": np.nan,
    }

    if not conf_file:
        return out

    path = Path(conf_file)
    if not path.exists() or path.stat().st_size == 0:
        return out

    try:
        data = json.loads(path.read_text(errors="ignore"))
    except Exception:
        return out

    scores = []
    if isinstance(data, dict):
        scores = data.get("confidenceScore") or data.get("plddt") or []
    elif isinstance(data, list):
        if data and isinstance(data[0], (int, float)):
            scores = data

    arr = np.asarray([safe_float(x) for x in scores], dtype=float)
    arr = arr[~np.isnan(arr)]

    if arr.size == 0:
        return out

    out.update({
        "alphafold_residue_count": int(arr.size),
        "alphafold_mean_plddt": float(np.mean(arr)),
        "alphafold_median_plddt": float(np.median(arr)),
        "alphafold_min_plddt": float(np.min(arr)),
        "alphafold_max_plddt": float(np.max(arr)),
        "alphafold_q05_plddt": float(np.quantile(arr, 0.05)),
        "alphafold_q25_plddt": float(np.quantile(arr, 0.25)),
        "alphafold_q75_plddt": float(np.quantile(arr, 0.75)),
        "alphafold_q95_plddt": float(np.quantile(arr, 0.95)),
        "alphafold_fraction_plddt_ge_90": float(np.mean(arr >= 90)),
        "alphafold_fraction_plddt_ge_70": float(np.mean(arr >= 70)),
        "alphafold_fraction_plddt_lt_50": float(np.mean(arr < 50)),
        "alphafold_predicted_disorder_fraction": float(np.mean(arr < 50)),
    })

    return out


def parse_pdb_geometry(pdb_file: str) -> Dict[str, Any]:
    out = {
        "alphafold_pdb_parse_ok": 0,
        "alphafold_atom_count": 0,
        "alphafold_ca_atom_count": 0,
        "alphafold_mean_bfactor_ca": np.nan,
        "alphafold_radius_of_gyration_ca": np.nan,
        "alphafold_mean_ca_distance_to_centroid": np.nan,
        "alphafold_max_ca_distance_to_centroid": np.nan,
        "alphafold_compactness_proxy_ca_count_over_rg": np.nan,
    }

    if not pdb_file:
        return out

    path = Path(pdb_file)
    if not path.exists() or path.stat().st_size == 0:
        return out

    coords = []
    ca_b = []
    atom_count = 0

    try:
        with open(path, "rt", errors="ignore") as f:
            for line in f:
                if not line.startswith("ATOM"):
                    continue

                atom_count += 1
                atom_name = line[12:16].strip()

                if atom_name != "CA":
                    continue

                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    b = float(line[60:66])
                except Exception:
                    continue

                coords.append((x, y, z))
                ca_b.append(b)

        out["alphafold_pdb_parse_ok"] = 1
        out["alphafold_atom_count"] = atom_count
        out["alphafold_ca_atom_count"] = len(coords)

        if ca_b:
            out["alphafold_mean_bfactor_ca"] = float(np.mean(ca_b))

        if coords:
            arr = np.asarray(coords, dtype=float)
            centroid = np.mean(arr, axis=0)
            d2 = np.sum((arr - centroid) ** 2, axis=1)
            d = np.sqrt(d2)
            rg = math.sqrt(float(np.mean(d2)))

            out["alphafold_radius_of_gyration_ca"] = rg
            out["alphafold_mean_ca_distance_to_centroid"] = float(np.mean(d))
            out["alphafold_max_ca_distance_to_centroid"] = float(np.max(d))
            out["alphafold_compactness_proxy_ca_count_over_rg"] = float(len(coords) / rg) if rg > 0 else np.nan

    except Exception:
        return out

    return out


# =============================================================================
# WORKER
# =============================================================================

def process_one(task: dict, cfg: dict) -> Tuple[dict, List[dict]]:
    session = make_session()

    outdir = Path(cfg["outdir"])
    force = bool(cfg["force"])
    no_download = bool(cfg["no_download"])
    max_pdb = int(cfg["max_pdb_per_uniprot"])

    downloads = outdir / "downloads"
    metadata_dir = downloads / "metadata"
    pdb_dir = downloads / "pdb"
    af_pdb_dir = downloads / "alphafold_pdb"
    af_conf_dir = downloads / "alphafold_confidence"

    gene = task["gene_symbol"]
    acc = task.get("uniprot_accession", "")

    protein_row = {
        "gene_symbol": gene,
        "hgnc_id": task.get("hgnc_id", ""),
        "ensembl_gene_id": task.get("ensembl_gene_id", ""),
        "entrez_id": task.get("entrez_id", ""),
        "uniprot_accession": acc,
        "is_primary_uniprot": task.get("is_primary_uniprot", 0),
        "uniprot_found": 0,
        "protein_name": "",
        "uniprot_id": "",
        "entry_type": "",
        "organism_name": "",
        "organism_id": "",
        "protein_length": np.nan,
        "pdb_available": 0,
        "pdb_count_uniprot": 0,
        "pdb_downloaded_count": 0,
        "pdb_ids_uniprot": "",
        "best_pdb_id": "",
        "best_pdb_resolution": np.nan,
        "best_pdb_method": "",
        "best_pdb_file": "",
        "alphafold_available": 0,
        "alphafold_version": "",
        "alphafold_pdb_file": "",
        "alphafold_confidence_file": "",
        "structure_source_priority": "none",
        "status": "",
    }

    pdb_rows = []

    if acc:
        uniprot_entry = fetch_uniprot_entry(session, acc, metadata_dir, force=force)
    else:
        found_acc, uniprot_entry = search_uniprot_by_gene(session, gene, metadata_dir, force=force)
        if found_acc:
            acc = found_acc
            protein_row["uniprot_accession"] = acc

    parsed = parse_uniprot(uniprot_entry)
    protein_row.update({k: v for k, v in parsed.items() if k != "pdb_ids"})

    pdb_ids = parsed.get("pdb_ids", []) or []
    protein_row["pdb_available"] = 1 if pdb_ids else 0
    protein_row["pdb_count_uniprot"] = len(pdb_ids)
    protein_row["pdb_ids_uniprot"] = ";".join(pdb_ids)

    selected_pdb_ids = pdb_ids if max_pdb == 0 else pdb_ids[:max_pdb]

    best_pdb_row = None
    best_res = np.inf
    pdb_downloaded = 0

    for pdb_id in selected_pdb_ids:
        meta_json = fetch_rcsb_entry(session, pdb_id, metadata_dir, force=force)
        meta = parse_rcsb_entry(pdb_id, meta_json)

        pdb_file = ""
        if not no_download:
            pdb_file = download_pdb(session, pdb_id, gene, pdb_dir, force=force)

        if pdb_file:
            pdb_downloaded += 1

        row = {
            "gene_symbol": gene,
            "uniprot_accession": acc,
            **meta,
            "pdb_file": pdb_file,
            "pdb_downloaded": 1 if pdb_file else 0,
            "source_pdb_api_url": RCSB_ENTRY_API.format(pdb_id=pdb_id),
            "source_pdb_download_url": RCSB_PDB_DOWNLOAD.format(pdb_id=pdb_id),
        }
        pdb_rows.append(row)

        res = safe_float(meta.get("pdb_resolution"))
        if not pd.isna(res) and res < best_res:
            best_res = res
            best_pdb_row = row
        elif best_pdb_row is None:
            best_pdb_row = row

    protein_row["pdb_downloaded_count"] = pdb_downloaded

    if best_pdb_row:
        protein_row["best_pdb_id"] = best_pdb_row.get("pdb_id", "")
        protein_row["best_pdb_resolution"] = best_pdb_row.get("pdb_resolution", np.nan)
        protein_row["best_pdb_method"] = best_pdb_row.get("pdb_method", "")
        protein_row["best_pdb_file"] = best_pdb_row.get("pdb_file", "")

    # Always AlphaFold, even if PDB exists.
    af = {
        "alphafold_available": 0,
        "alphafold_version": "",
        "alphafold_pdb_file": "",
        "alphafold_confidence_file": "",
        "alphafold_pdb_url": "",
        "alphafold_confidence_url": "",
    }

    if acc and not no_download:
        af = download_alphafold(
            session=session,
            acc=acc,
            gene_symbol=gene,
            af_pdb_dir=af_pdb_dir,
            af_conf_dir=af_conf_dir,
            force=force,
        )
    elif acc and no_download:
        gene_af_dir = af_pdb_dir / safe_token(gene)
        gene_conf_dir = af_conf_dir / safe_token(gene)
        for ver in ALPHAFOLD_VERSIONS:
            pdb_local = gene_af_dir / f"AF-{acc}-F1-model_{ver}.pdb"
            conf_local = gene_conf_dir / f"AF-{acc}-F1-confidence_{ver}.json"
            if pdb_local.exists() and pdb_local.stat().st_size > 100:
                af = {
                    "alphafold_available": 1,
                    "alphafold_version": ver,
                    "alphafold_pdb_file": str(pdb_local),
                    "alphafold_confidence_file": str(conf_local) if conf_local.exists() else "",
                    "alphafold_pdb_url": AF_PDB_URL.format(acc=acc, ver=ver),
                    "alphafold_confidence_url": AF_CONF_URL.format(acc=acc, ver=ver),
                }
                break

    protein_row["alphafold_available"] = af.get("alphafold_available", 0)
    protein_row["alphafold_version"] = af.get("alphafold_version", "")
    protein_row["alphafold_pdb_file"] = af.get("alphafold_pdb_file", "")
    protein_row["alphafold_confidence_file"] = af.get("alphafold_confidence_file", "")
    protein_row["source_alphafold_model_url"] = af.get("alphafold_pdb_url", "")
    protein_row["source_alphafold_confidence_url"] = af.get("alphafold_confidence_url", "")

    conf_features = parse_alphafold_confidence(af.get("alphafold_confidence_file", ""))
    geom_features = parse_pdb_geometry(af.get("alphafold_pdb_file", ""))

    protein_row.update(conf_features)
    protein_row.update(geom_features)

    # Fallback: AlphaFold PDB B-factor column is pLDDT-like.
    if pd.isna(protein_row.get("alphafold_mean_plddt", np.nan)):
        if not pd.isna(protein_row.get("alphafold_mean_bfactor_ca", np.nan)):
            protein_row["alphafold_mean_plddt"] = protein_row["alphafold_mean_bfactor_ca"]

    if protein_row["pdb_available"] and protein_row["alphafold_available"]:
        protein_row["structure_source_priority"] = "PDB_and_AlphaFold"
    elif protein_row["pdb_available"]:
        protein_row["structure_source_priority"] = "PDB_only"
    elif protein_row["alphafold_available"]:
        protein_row["structure_source_priority"] = "AlphaFold_only"
    else:
        protein_row["structure_source_priority"] = "none"

    if not acc:
        protein_row["status"] = "No UniProt accession"
    elif not protein_row["uniprot_found"]:
        protein_row["status"] = "UniProt not found"
    elif not protein_row["pdb_available"] and not protein_row["alphafold_available"]:
        protein_row["status"] = "No PDB or AlphaFold"
    else:
        protein_row["status"] = "OK"

    return protein_row, pdb_rows


# =============================================================================
# AGGREGATION
# =============================================================================

def summarise_values(values: List[float], prefix: str) -> Dict[str, Any]:
    arr = np.asarray([safe_float(v) for v in values], dtype=float)
    arr = arr[~np.isnan(arr)]

    if arr.size == 0:
        return {
            f"{prefix}_count": 0,
            f"{prefix}_mean": np.nan,
            f"{prefix}_median": np.nan,
            f"{prefix}_min": np.nan,
            f"{prefix}_max": np.nan,
        }

    return {
        f"{prefix}_count": int(arr.size),
        f"{prefix}_mean": float(np.mean(arr)),
        f"{prefix}_median": float(np.median(arr)),
        f"{prefix}_min": float(np.min(arr)),
        f"{prefix}_max": float(np.max(arr)),
    }


def aggregate_gene_features(hgnc: pd.DataFrame, protein_df: pd.DataFrame) -> pd.DataFrame:
    genes = pd.DataFrame({"gene_symbol": sorted(hgnc["gene_symbol"].unique().tolist())})
    records = []

    for gene, g in protein_df.groupby("gene_symbol", dropna=False):
        g = g.copy()

        pdb_available = int((pd.to_numeric(g["pdb_available"], errors="coerce").fillna(0) > 0).any())
        af_available = int((pd.to_numeric(g["alphafold_available"], errors="coerce").fillna(0) > 0).any())

        # Choose representative row:
        # primary UniProt > AlphaFold available > PDB available
        rep = (
            g.sort_values(
                ["is_primary_uniprot", "alphafold_available", "pdb_available"],
                ascending=False,
            )
            .iloc[0]
            .to_dict()
        )

        # Best PDB = lowest resolution among rows.
        best = g.copy()
        best["res_num"] = pd.to_numeric(best["best_pdb_resolution"], errors="coerce")
        best_with_res = best[best["res_num"].notna()]
        if len(best_with_res):
            best_row = best_with_res.sort_values("res_num").iloc[0].to_dict()
        else:
            best_row = best.iloc[0].to_dict()

        rec = {
            "gene_symbol": gene,

            "feature4_structure_uniprot_count": int(g["uniprot_accession"].replace("", np.nan).dropna().nunique()),
            "feature4_structure_primary_uniprot": rep.get("uniprot_accession", ""),
            "feature4_structure_primary_protein_name": rep.get("protein_name", ""),
            "feature4_structure_primary_protein_length": rep.get("protein_length", np.nan),

            "feature4_structure_pdb_available": pdb_available,
            "feature4_structure_pdb_count_total": int(pd.to_numeric(g["pdb_count_uniprot"], errors="coerce").fillna(0).sum()),
            "feature4_structure_pdb_downloaded_count": int(pd.to_numeric(g["pdb_downloaded_count"], errors="coerce").fillna(0).sum()),
            "feature4_structure_best_pdb_id": best_row.get("best_pdb_id", ""),
            "feature4_structure_best_pdb_resolution": best_row.get("best_pdb_resolution", np.nan),
            "feature4_structure_best_pdb_method": best_row.get("best_pdb_method", ""),
            "feature4_structure_best_pdb_file": best_row.get("best_pdb_file", ""),

            "feature4_structure_alphafold_available": af_available,
            "feature4_structure_alphafold_version": rep.get("alphafold_version", ""),
            "feature4_structure_alphafold_pdb_file": rep.get("alphafold_pdb_file", ""),
            "feature4_structure_alphafold_confidence_file": rep.get("alphafold_confidence_file", ""),

            "feature4_structure_alphafold_residue_count": rep.get("alphafold_residue_count", np.nan),
            "feature4_structure_alphafold_mean_plddt": rep.get("alphafold_mean_plddt", np.nan),
            "feature4_structure_alphafold_median_plddt": rep.get("alphafold_median_plddt", np.nan),
            "feature4_structure_alphafold_min_plddt": rep.get("alphafold_min_plddt", np.nan),
            "feature4_structure_alphafold_max_plddt": rep.get("alphafold_max_plddt", np.nan),
            "feature4_structure_alphafold_q05_plddt": rep.get("alphafold_q05_plddt", np.nan),
            "feature4_structure_alphafold_q25_plddt": rep.get("alphafold_q25_plddt", np.nan),
            "feature4_structure_alphafold_q75_plddt": rep.get("alphafold_q75_plddt", np.nan),
            "feature4_structure_alphafold_q95_plddt": rep.get("alphafold_q95_plddt", np.nan),
            "feature4_structure_alphafold_fraction_plddt_ge_90": rep.get("alphafold_fraction_plddt_ge_90", np.nan),
            "feature4_structure_alphafold_fraction_plddt_ge_70": rep.get("alphafold_fraction_plddt_ge_70", np.nan),
            "feature4_structure_alphafold_fraction_plddt_lt_50": rep.get("alphafold_fraction_plddt_lt_50", np.nan),
            "feature4_structure_alphafold_predicted_disorder_fraction": rep.get("alphafold_predicted_disorder_fraction", np.nan),

            "feature4_structure_alphafold_ca_atom_count": rep.get("alphafold_ca_atom_count", np.nan),
            "feature4_structure_alphafold_atom_count": rep.get("alphafold_atom_count", np.nan),
            "feature4_structure_alphafold_radius_of_gyration_ca": rep.get("alphafold_radius_of_gyration_ca", np.nan),
            "feature4_structure_alphafold_mean_ca_distance_to_centroid": rep.get("alphafold_mean_ca_distance_to_centroid", np.nan),
            "feature4_structure_alphafold_max_ca_distance_to_centroid": rep.get("alphafold_max_ca_distance_to_centroid", np.nan),
            "feature4_structure_alphafold_compactness_proxy_ca_count_over_rg": rep.get("alphafold_compactness_proxy_ca_count_over_rg", np.nan),

            "feature4_structure_source_priority": rep.get("structure_source_priority", ""),
            "feature4_structure_has_any_structure": int(pdb_available or af_available),
        }

        plddt_vals = pd.to_numeric(g["alphafold_mean_plddt"], errors="coerce").dropna().tolist()
        rec.update(summarise_values(plddt_vals, "feature4_structure_alphafold_mean_plddt_across_uniprot"))

        records.append(rec)

    out = pd.DataFrame(records)
    out = genes.merge(out, on="gene_symbol", how="left")

    zero_cols = [
        "feature4_structure_uniprot_count",
        "feature4_structure_pdb_available",
        "feature4_structure_pdb_count_total",
        "feature4_structure_pdb_downloaded_count",
        "feature4_structure_alphafold_available",
        "feature4_structure_has_any_structure",
    ]

    for c in zero_cols:
        if c in out.columns:
            out[c] = out[c].fillna(0).astype(int)

    text_cols = [
        "feature4_structure_primary_uniprot",
        "feature4_structure_primary_protein_name",
        "feature4_structure_best_pdb_id",
        "feature4_structure_best_pdb_method",
        "feature4_structure_best_pdb_file",
        "feature4_structure_alphafold_version",
        "feature4_structure_alphafold_pdb_file",
        "feature4_structure_alphafold_confidence_file",
        "feature4_structure_source_priority",
    ]

    for c in text_cols:
        if c in out.columns:
            out[c] = out[c].fillna("")

    return out


def merge_with_hgnc(hgnc: pd.DataFrame, gene_features: pd.DataFrame) -> pd.DataFrame:
    merged = hgnc.merge(gene_features, on="gene_symbol", how="left")
    feature_cols = [c for c in gene_features.columns if c != "gene_symbol"]
    merged["feature4_structure_has_any_feature"] = merged[feature_cols].notna().any(axis=1).astype(int)
    merged["feature4_structure_n_nonmissing_features"] = merged[feature_cols].notna().sum(axis=1)
    return merged


def write_summary(
    outdir: Path,
    hgnc: pd.DataFrame,
    protein_df: pd.DataFrame,
    pdb_df: pd.DataFrame,
    gene_df: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    path = outdir / "feature4_structure_summary.txt"

    lines = []
    lines.append("Feature 4: PDB + AlphaFold structure features")
    lines.append("=" * 80)
    lines.append(f"Created at: {now_iso()}")
    lines.append(f"HGNC genes: {hgnc.shape[0]}")
    lines.append(f"UniProt/protein rows: {protein_df.shape[0]}")
    lines.append(f"PDB long rows: {pdb_df.shape[0]}")
    lines.append(f"Gene feature rows: {gene_df.shape[0]}")
    lines.append(f"Gene feature columns: {gene_df.shape[1]}")
    lines.append("")
    lines.append("Run parameters")
    lines.append(f"workers: {args.workers}")
    lines.append(f"max_pdb_per_uniprot: {'all' if args.max_pdb_per_uniprot == 0 else args.max_pdb_per_uniprot}")
    lines.append(f"limit_genes: {args.limit_genes if args.limit_genes else 'none'}")
    lines.append(f"no_download: {args.no_download}")
    lines.append(f"force: {args.force}")
    lines.append("")
    lines.append("Coverage")
    if "feature4_structure_pdb_available" in gene_df.columns:
        lines.append(f"PDB available fraction: {gene_df['feature4_structure_pdb_available'].mean():.4f}")
    if "feature4_structure_alphafold_available" in gene_df.columns:
        lines.append(f"AlphaFold available fraction: {gene_df['feature4_structure_alphafold_available'].mean():.4f}")
    if "feature4_structure_has_any_structure" in gene_df.columns:
        lines.append(f"Any structure fraction: {gene_df['feature4_structure_has_any_structure'].mean():.4f}")
    if "feature4_structure_pdb_count_total" in gene_df.columns:
        lines.append(f"Total PDB crossrefs: {gene_df['feature4_structure_pdb_count_total'].sum():.0f}")
    if "feature4_structure_pdb_downloaded_count" in gene_df.columns:
        lines.append(f"Total PDB downloaded: {gene_df['feature4_structure_pdb_downloaded_count'].sum():.0f}")
    lines.append("")
    lines.append("Leakage policy")
    lines.append("Included: PDB structure metadata/files and AlphaFold structural confidence/geometric summaries.")
    lines.append("Excluded: drug-target labels, approved-drug annotations, ChEMBL/DrugBank/DGIdb labels, drug-response data.")

    path.write_text("\n".join(lines) + "\n")
    log(f"[SUMMARY] {path}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hgnc", default=str(DEFAULT_HGNC_PATH))
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR))
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--limit-genes", type=int, default=0)
    parser.add_argument("--all-hgnc-genes", action="store_true")
    parser.add_argument("--max-pdb-per-uniprot", type=int, default=5, help="0 = all PDBs; recommended 5 for speed/storage.")
    parser.add_argument("--no-download", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    outdir = mkdir(Path(args.outdir))
    downloads = mkdir(outdir / "downloads")
    mkdir(downloads / "pdb")
    mkdir(downloads / "alphafold_pdb")
    mkdir(downloads / "alphafold_confidence")
    mkdir(downloads / "metadata")
    mkdir(downloads / "metadata" / "uniprot")
    mkdir(downloads / "metadata" / "uniprot_search")
    mkdir(downloads / "metadata" / "rcsb_entry")
    processed = mkdir(outdir / "processed")

    log("=" * 100)
    log("FEATURE 4: STRUCTURE DOWNLOAD + STRUCTURE FEATURES")
    log("=" * 100)
    log(f"[HGNC]                {args.hgnc}")
    log(f"[OUTDIR]              {outdir.resolve()}")
    log(f"[WORKERS]             {args.workers}")
    log(f"[MAX PDB/UNIPROT]     {'all' if args.max_pdb_per_uniprot == 0 else args.max_pdb_per_uniprot}")
    log(f"[LIMIT GENES]         {args.limit_genes if args.limit_genes else 'none'}")
    log(f"[NO DOWNLOAD]         {args.no_download}")
    log(f"[FORCE]               {args.force}")
    log("=" * 100)

    hgnc = load_hgnc(Path(args.hgnc), protein_coding_only=not args.all_hgnc_genes)
    tasks = build_tasks(hgnc, limit_genes=args.limit_genes)

    task_manifest = processed / "feature4_structure_uniprot_task_manifest.csv"
    pd.DataFrame(tasks).to_csv(task_manifest, index=False)
    log(f"[TASKS] {len(tasks)}")
    log(f"[TASK MANIFEST] {task_manifest}")

    cfg = {
        "outdir": str(outdir),
        "force": args.force,
        "no_download": args.no_download,
        "max_pdb_per_uniprot": args.max_pdb_per_uniprot,
    }

    protein_rows = []
    pdb_rows = []

    partial_protein = processed / "feature4_structure_uniprot_table.partial.csv"
    partial_pdb = processed / "feature4_structure_pdb_long.partial.csv"

    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as ex:
        futures = {ex.submit(process_one, task, cfg): task for task in tasks}

        for i, future in enumerate(as_completed(futures), start=1):
            task = futures[future]

            try:
                protein_row, pdb_subrows = future.result()
            except Exception as exc:
                protein_row = {
                    "gene_symbol": task.get("gene_symbol", ""),
                    "uniprot_accession": task.get("uniprot_accession", ""),
                    "status": f"ERROR: {exc}",
                }
                pdb_subrows = []

            protein_rows.append(protein_row)
            pdb_rows.extend(pdb_subrows)

            if i % 100 == 0 or i == len(tasks):
                log(f"[PROGRESS] {i}/{len(tasks)}")
                pd.DataFrame(protein_rows).to_csv(partial_protein, index=False)
                if pdb_rows:
                    pd.DataFrame(pdb_rows).to_csv(partial_pdb, index=False)

    protein_df = pd.DataFrame(protein_rows)
    pdb_df = pd.DataFrame(pdb_rows)

    protein_out = processed / "feature4_structure_uniprot_table.csv"
    pdb_out = processed / "feature4_structure_pdb_long.csv"
    gene_out = processed / "feature4_structure_gene_features.csv"
    merged_out = processed / "feature4_structure_hgnc_merged.csv"

    protein_df.to_csv(protein_out, index=False)
    pdb_df.to_csv(pdb_out, index=False)

    gene_features = aggregate_gene_features(hgnc, protein_df)
    gene_features.to_csv(gene_out, index=False)

    merged = merge_with_hgnc(hgnc, gene_features)
    merged.to_csv(merged_out, index=False)

    metadata = {
        "created_at": now_iso(),
        "hgnc": str(Path(args.hgnc)),
        "outdir": str(outdir.resolve()),
        "workers": args.workers,
        "limit_genes": args.limit_genes,
        "protein_coding_only": not args.all_hgnc_genes,
        "max_pdb_per_uniprot": args.max_pdb_per_uniprot,
        "alphafold_versions": ALPHAFOLD_VERSIONS,
        "outputs": {
            "uniprot_table": str(protein_out),
            "pdb_long": str(pdb_out),
            "gene_features": str(gene_out),
            "hgnc_merged": str(merged_out),
        },
    }
    write_json(outdir / "feature4_structure_run_metadata.json", metadata)

    write_summary(outdir, hgnc, protein_df, pdb_df, gene_features, args)

    log("=" * 100)
    log("[DONE]")
    log(f"[UNIPROT TABLE] {protein_out}")
    log(f"[PDB LONG]      {pdb_out}")
    log(f"[GENE FEATURES] {gene_out}")
    log(f"[HGNC MERGED]   {merged_out}")
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