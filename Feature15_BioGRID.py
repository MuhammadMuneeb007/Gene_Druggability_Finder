#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature15_BioGRID.py

Feature 15: BioGRID curated interaction-network features.

Purpose
-------
Build leakage-aware gene-level BioGRID network features from curated human
protein/genetic interaction data.

This complements STRING:
    Feature 2  = STRING functional / predicted / integrated PPI features
    Feature 15 = BioGRID curated experimental interaction features

Inputs
------
HGNC:
    databases/HGNC/hgnc_complete_set.txt

BioGRID:
    feature_databases/BioGRID/BIOGRID-ORGANISM-LATEST.tab3.zip

or auto-downloaded from:
    https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip

Outputs
-------
feature15_biogrid/
    downloads/BIOGRID-ORGANISM-LATEST.tab3.zip
    processed/feature15_biogrid_interaction_long.csv
    processed/feature15_biogrid_gene_features.csv
    processed/feature15_biogrid_hgnc_merged.csv
    feature15_biogrid_summary.txt
    feature15_biogrid_run_metadata.json

Run
---
    python Feature15_BioGRID.py --download

Fast test:
    python Feature15_BioGRID.py --download --max-rows 100000

Leakage policy
--------------
Included:
    curated human physical/genetic interaction counts
    network degree
    publication count
    experimental-system diversity
    throughput diversity
    interaction-type diversity
    self-interaction flag

Excluded:
    BioGRID chemical interaction files
    drug/chemical names
    ChEMBL
    DrugBank
    DGIdb
    Open Targets
    Pharos
    known drug-target labels
    clinical target labels
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import shutil
import subprocess
import sys
import time
import urllib.request
import zipfile
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pandas as pd


# =============================================================================
# DEFAULTS
# =============================================================================

DEFAULT_HGNC_PATH = Path("databases") / "HGNC" / "hgnc_complete_set.txt"
DEFAULT_DBDIR = Path("feature_databases") / "BioGRID"
DEFAULT_OUTDIR = Path("feature15_biogrid")

BIOGRID_ORGANISM_TAB3_URL = (
    "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip"
)

# Fallback URLs for robust automated download.
# LATEST is preferred for routine use. The release-specific archive is useful
# when BioGRID blocks Latest-Release or when you want a reproducible frozen file.
BIOGRID_ORGANISM_TAB3_URLS = [
    "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip",
]

HUMAN_TAX_ID = "9606"

# Common BioGRID Tab 3.0 columns.
# We still auto-detect columns because exact capitalization may vary.
INTERACTION_TYPE_MAP = {
    "physical": "physical",
    "genetic": "genetic",
}


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
    s = clean_text(x).upper()
    s = re.sub(r"\s+", "", s)
    return s


def safe_float(x: Any) -> float:
    try:
        if x is None or pd.isna(x):
            return np.nan
        s = str(x).strip()
        if s in {"", ".", "NA", "NaN", "nan", "None", "-"}:
            return np.nan
        return float(s)
    except Exception:
        return np.nan


def shannon_entropy(items: Iterable[str]) -> float:
    items = [clean_text(x) for x in items if clean_text(x)]
    if not items:
        return 0.0

    c = Counter(items)
    n = sum(c.values())

    ent = 0.0
    for _, count in c.items():
        p = count / n
        if p > 0:
            ent -= p * math.log2(p)

    return float(ent)


def open_text_maybe_gzip(path: Path):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", errors="ignore")
    return open(path, "rt", errors="ignore")


def find_col(df: pd.DataFrame, candidates: List[str], contains: Optional[List[str]] = None) -> Optional[str]:
    lower_map = {str(c).strip().lower(): c for c in df.columns}

    for c in candidates:
        key = c.lower().strip()
        if key in lower_map:
            return lower_map[key]

    if contains:
        for c in df.columns:
            low = str(c).lower()
            if all(x.lower() in low for x in contains):
                return c

    return None


def is_valid_zip(path: Path) -> bool:
    """
    Return True only when path is a readable, non-empty ZIP file.

    This avoids a common BioGRID failure mode where an HTML 403/404 page is
    saved with a .zip filename.
    """
    if not path.exists() or path.stat().st_size == 0:
        return False

    try:
        with zipfile.ZipFile(path, "r") as z:
            names = z.namelist()
            if not names:
                log("[ZIP INVALID] Empty ZIP file.")
                return False
            bad_member = z.testzip()
            if bad_member is not None:
                log(f"[ZIP INVALID] First bad member: {bad_member}")
                return False
        return True
    except Exception as exc:
        log(f"[ZIP INVALID] {path}: {exc}")
        return False


def looks_like_html_error(path: Path) -> bool:
    """
    Detect downloaded HTML error pages saved as .zip files.
    """
    if not path.exists() or path.stat().st_size == 0:
        return True

    try:
        with open(path, "rb") as f:
            head = f.read(1000).lower()
        markers = [
            b"<html",
            b"<!doctype html",
            b"403 forbidden",
            b"404 not found",
            b"access denied",
            b"forbidden",
        ]
        return any(m in head for m in markers)
    except Exception:
        return False


def download_with_python(url: str, outpath: Path, timeout: int = 180) -> bool:
    """
    Download using urllib with browser-like headers.
    BioGRID may reject the default Python urllib user-agent with HTTP 403.
    """
    tmp = outpath.with_suffix(outpath.suffix + ".tmp")

    headers = {
        "User-Agent": (
            "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
            "(KHTML, like Gecko) Chrome/120.0 Safari/537.36"
        ),
        "Accept": "application/zip,application/octet-stream,*/*",
        "Accept-Language": "en-US,en;q=0.9",
        "Connection": "keep-alive",
    }

    req = urllib.request.Request(url, headers=headers)

    log(f"[DOWNLOAD PYTHON] {url}")
    log(f"[TO]              {outpath}")

    try:
        with urllib.request.urlopen(req, timeout=timeout) as response:
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
                            log(
                                f"[DOWNLOAD PYTHON] {downloaded / 1024 / 1024:.1f} MB / "
                                f"{total / 1024 / 1024:.1f} MB ({pct:.1f}%)"
                            )
                        else:
                            log(f"[DOWNLOAD PYTHON] {downloaded / 1024 / 1024:.1f} MB")
                        last_print = time.time()

        tmp.rename(outpath)

        if looks_like_html_error(outpath):
            log("[DOWNLOAD PYTHON FAILED] Downloaded file looks like an HTML/error page.")
            try:
                outpath.unlink()
            except Exception:
                pass
            return False

        if not is_valid_zip(outpath):
            log("[DOWNLOAD PYTHON FAILED] Downloaded file is not a valid ZIP.")
            try:
                outpath.unlink()
            except Exception:
                pass
            return False

        log(f"[DOWNLOAD PYTHON OK] {outpath} ({outpath.stat().st_size / 1024 / 1024:.2f} MB)")
        return True

    except Exception as exc:
        log(f"[DOWNLOAD PYTHON ERROR] {exc}")
        try:
            if tmp.exists():
                tmp.unlink()
        except Exception:
            pass
        return False


def download_with_command(url: str, outpath: Path, tool: str) -> bool:
    """
    Download using wget or curl as a fallback.
    """
    tmp = outpath.with_suffix(outpath.suffix + f".{tool}.tmp")

    if tool == "wget":
        cmd = [
            "wget",
            "--user-agent=Mozilla/5.0",
            "--tries=3",
            "--timeout=60",
            "-O",
            str(tmp),
            url,
        ]
    elif tool == "curl":
        cmd = [
            "curl",
            "-L",
            "--retry",
            "3",
            "--connect-timeout",
            "60",
            "-A",
            "Mozilla/5.0",
            "-o",
            str(tmp),
            url,
        ]
    else:
        return False

    log(f"[DOWNLOAD {tool.upper()}] {url}")
    log(f"[CMD] {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

        if result.returncode != 0:
            log(f"[DOWNLOAD {tool.upper()} ERROR] returncode={result.returncode}")
            if result.stderr:
                log(result.stderr[-2000:])
            try:
                if tmp.exists():
                    tmp.unlink()
            except Exception:
                pass
            return False

        tmp.rename(outpath)

        if looks_like_html_error(outpath):
            log(f"[DOWNLOAD {tool.upper()} FAILED] Downloaded file looks like an HTML/error page.")
            try:
                outpath.unlink()
            except Exception:
                pass
            return False

        if not is_valid_zip(outpath):
            log(f"[DOWNLOAD {tool.upper()} FAILED] Downloaded file is not a valid ZIP.")
            try:
                outpath.unlink()
            except Exception:
                pass
            return False

        log(f"[DOWNLOAD {tool.upper()} OK] {outpath} ({outpath.stat().st_size / 1024 / 1024:.2f} MB)")
        return True

    except FileNotFoundError:
        log(f"[DOWNLOAD {tool.upper()} SKIP] {tool} not found on this system.")
        return False
    except Exception as exc:
        log(f"[DOWNLOAD {tool.upper()} ERROR] {exc}")
        try:
            if tmp.exists():
                tmp.unlink()
        except Exception:
            pass
        return False


def download_file(url: str, outpath: Path, force: bool = False) -> None:
    """
    Robust BioGRID downloader.

    Strategy:
        1. Reuse an existing valid ZIP if present.
        2. Try Python urllib with browser-like headers.
        3. Try wget with browser-like user-agent.
        4. Try curl with browser-like user-agent.
        5. Try fallback BioGRID URLs.

    Raises RuntimeError only if all download attempts fail.
    """
    mkdir(outpath.parent)

    if outpath.exists() and outpath.stat().st_size > 0 and not force:
        if is_valid_zip(outpath):
            log(f"[DOWNLOAD SKIP] Existing valid ZIP: {outpath}")
            return
        log(f"[DOWNLOAD WARNING] Existing file is invalid. Removing: {outpath}")
        try:
            outpath.unlink()
        except Exception:
            pass

    if outpath.exists() and force:
        log(f"[REMOVE EXISTING] {outpath}")
        outpath.unlink()

    urls = [url]
    for u in BIOGRID_ORGANISM_TAB3_URLS:
        if u not in urls:
            urls.append(u)

    for candidate_url in urls:
        log("=" * 100)
        log(f"[BIOGRID DOWNLOAD ATTEMPT] {candidate_url}")

        if download_with_python(candidate_url, outpath):
            return

        if download_with_command(candidate_url, outpath, tool="wget"):
            return

        if download_with_command(candidate_url, outpath, tool="curl"):
            return

    raise RuntimeError(
        "Could not download BioGRID automatically. Manually download "
        "BIOGRID-ORGANISM-LATEST.tab3.zip or BIOGRID-ORGANISM-5.0.258.tab3.zip "
        "and provide --biogrid-file."
    )

def find_local_biogrid_file(dbdir: Path) -> Optional[Path]:
    if not dbdir.exists():
        return None

    patterns = [
        "BIOGRID-ORGANISM-LATEST.tab3.zip",
        "*BIOGRID*ORGANISM*tab3*.zip",
        "*BIOGRID*Homo_sapiens*tab3*.txt",
        "*BIOGRID*Homo_sapiens*tab3*.txt.gz",
        "*Homo_sapiens*tab3*.txt",
        "*9606*tab3*.txt",
        "*.tab3.txt",
        "*.tab3.txt.gz",
        "*.tab3.zip",
        "*.zip",
    ]

    files = []
    for pat in patterns:
        files.extend(list(dbdir.rglob(pat)))

    files = [p for p in files if p.exists() and p.stat().st_size > 0]

    if not files:
        return None

    ranked = []
    for p in files:
        name = p.name.lower()
        score = 0
        if "biogrid-organism-latest.tab3.zip" in name:
            score += 100
        if "organism" in name:
            score += 30
        if "homo_sapiens" in name or "homo-sapiens" in name or "9606" in name:
            score += 30
        if "tab3" in name:
            score += 20
        if name.endswith(".zip"):
            score += 5
        ranked.append((score, p.stat().st_size, p))

    ranked.sort(reverse=True)
    return ranked[0][2]


def extract_human_tab3_from_zip(zip_path: Path, extract_dir: Path, force: bool = False) -> Path:
    """
    Extract the Homo sapiens BioGRID Tab3 file from BIOGRID-ORGANISM-LATEST.tab3.zip.
    """
    mkdir(extract_dir)

    if not zip_path.exists():
        raise FileNotFoundError(zip_path)

    log("=" * 100)
    log(f"[INSPECT ZIP] {zip_path}")

    with zipfile.ZipFile(zip_path, "r") as z:
        names = z.namelist()

        # Prefer human organism file.
        human_candidates = []
        for name in names:
            low = name.lower()
            if (
                ("homo_sapiens" in low or "homo-sapiens" in low or "9606" in low)
                and low.endswith((".txt", ".tab3.txt"))
            ):
                human_candidates.append(name)

        if not human_candidates:
            # fallback: any tab3 txt that may be human
            for name in names:
                low = name.lower()
                if low.endswith(".txt") and "human" in low:
                    human_candidates.append(name)

        if not human_candidates:
            raise RuntimeError(
                "Could not find Homo sapiens file inside BioGRID organism zip. "
                f"First 20 zip entries: {names[:20]}"
            )

        # Usually only one Homo_sapiens file.
        member = sorted(human_candidates, key=len)[0]
        outpath = extract_dir / Path(member).name

        if outpath.exists() and outpath.stat().st_size > 0 and not force:
            log(f"[EXTRACT SKIP] Already exists: {outpath}")
            return outpath

        log(f"[EXTRACT HUMAN FILE] {member}")
        with z.open(member) as zin, open(outpath, "wb") as fout:
            shutil.copyfileobj(zin, fout)

    log(f"[EXTRACT DONE] {outpath} ({outpath.stat().st_size / 1024 / 1024:.2f} MB)")
    return outpath


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
        raise RuntimeError("HGNC file must contain symbol column.")

    before = len(hgnc)

    if protein_coding_only:
        if "locus_group" in hgnc.columns:
            hgnc = hgnc[hgnc["locus_group"].astype(str).str.lower().eq("protein-coding gene")].copy()
        elif "locus_type" in hgnc.columns:
            hgnc = hgnc[hgnc["locus_type"].astype(str).str.lower().str.contains("protein", na=False)].copy()

    hgnc["gene_symbol"] = hgnc["symbol"].map(normalize_symbol)

    for col in ["alias_symbol", "prev_symbol", "entrez_id", "ensembl_gene_id", "uniprot_ids", "name"]:
        if col not in hgnc.columns:
            hgnc[col] = ""

    hgnc = hgnc.drop_duplicates(subset=["gene_symbol"], keep="first").copy()

    if limit_genes and limit_genes > 0:
        hgnc = hgnc.head(limit_genes).copy()

    log(f"[HGNC ROWS] {before} -> {len(hgnc)}")
    log(f"[PROTEIN CODING ONLY] {protein_coding_only}")

    return hgnc


def build_hgnc_maps(hgnc: pd.DataFrame) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Returns:
        symbol_alias_to_current_symbol
        entrez_to_current_symbol
    """
    symbol_map: Dict[str, str] = {}
    entrez_map: Dict[str, str] = {}

    for _, row in hgnc.iterrows():
        gene = normalize_symbol(row.get("gene_symbol", ""))
        if not gene:
            continue

        symbol_map[gene] = gene

        entrez = clean_text(row.get("entrez_id", ""))
        if entrez:
            entrez_map[entrez] = gene

        for field in ["alias_symbol", "prev_symbol"]:
            raw = clean_text(row.get(field, ""))
            if not raw:
                continue

            for part in re.split(r"[|,;]+", raw):
                s = normalize_symbol(part)
                if s and s not in symbol_map:
                    symbol_map[s] = gene

    log("=" * 100)
    log("[HGNC MAPS]")
    log(f"[SYMBOLS + ALIASES] {len(symbol_map)}")
    log(f"[ENTREZ IDS]         {len(entrez_map)}")

    return symbol_map, entrez_map


# =============================================================================
# BIOGRID READING / MAPPING
# =============================================================================

def read_biogrid_tab3(path: Path, max_rows: int = 0) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)

    nrows = max_rows if max_rows and max_rows > 0 else None

    log("=" * 100)
    log(f"[READ BIOGRID TAB3] {path}")

    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        low_memory=False,
        compression="infer",
        nrows=nrows,
    )
    df.columns = [str(c).strip() for c in df.columns]

    log(f"[BIOGRID SHAPE] {df.shape}")
    log(f"[BIOGRID COLUMNS FIRST 20] {list(df.columns)[:20]}")

    return df


def detect_biogrid_columns(df: pd.DataFrame) -> Dict[str, Optional[str]]:
    cols = {
        "biogrid_id_a": find_col(df, ["BioGRID Interaction ID", "BIOGRID_INTERACTION_ID"], contains=["interaction", "id"]),

        "entrez_a": find_col(df, ["Entrez Gene Interactor A", "Entrez Gene Interactor A ", "ENTREZ_GENE_A"], contains=["entrez", "interactor", "a"]),
        "entrez_b": find_col(df, ["Entrez Gene Interactor B", "ENTREZ_GENE_B"], contains=["entrez", "interactor", "b"]),

        "symbol_a": find_col(df, ["Official Symbol Interactor A", "OFFICIAL_SYMBOL_A"], contains=["official", "symbol", "a"]),
        "symbol_b": find_col(df, ["Official Symbol Interactor B", "OFFICIAL_SYMBOL_B"], contains=["official", "symbol", "b"]),

        "aliases_a": find_col(df, ["Aliases Interactor A", "ALIASES_FOR_A"], contains=["aliases", "a"]),
        "aliases_b": find_col(df, ["Aliases Interactor B", "ALIASES_FOR_B"], contains=["aliases", "b"]),

        "organism_a": find_col(df, ["Organism Interactor A", "ORGANISM_A"], contains=["organism", "a"]),
        "organism_b": find_col(df, ["Organism Interactor B", "ORGANISM_B"], contains=["organism", "b"]),

        "experimental_system": find_col(df, ["Experimental System", "EXPERIMENTAL_SYSTEM"], contains=["experimental", "system"]),
        "experimental_system_type": find_col(df, ["Experimental System Type", "EXPERIMENTAL_SYSTEM_TYPE"], contains=["experimental", "system", "type"]),

        "author": find_col(df, ["Author", "AUTHOR"], contains=["author"]),
        "pubmed_id": find_col(df, ["Pubmed ID", "PUBMED_ID", "PUBMEDID"], contains=["pubmed"]),
        "throughput": find_col(df, ["Throughput", "THROUGHPUT"], contains=["throughput"]),
        "score": find_col(df, ["Score", "SCORE"], contains=["score"]),
        "modification": find_col(df, ["Modification", "MODIFICATION"], contains=["modification"]),
        "phenotypes": find_col(df, ["Phenotypes", "PHENOTYPES"], contains=["phenotypes"]),
        "qualifications": find_col(df, ["Qualifications", "QUALIFICATIONS"], contains=["qualifications"]),
        "tags": find_col(df, ["Tags", "TAGS"], contains=["tags"]),
        "source_database": find_col(df, ["Source Database", "SOURCE_DATABASE"], contains=["source", "database"]),
    }

    log("=" * 100)
    log("[DETECTED BIOGRID COLUMNS]")
    for k, v in cols.items():
        log(f"{k:28s} = {v}")

    return cols


def map_interactor_to_gene(
    symbol: Any,
    entrez: Any,
    aliases: Any,
    symbol_map: Dict[str, str],
    entrez_map: Dict[str, str],
) -> str:
    entrez_s = clean_text(entrez)
    if entrez_s and entrez_s in entrez_map:
        return entrez_map[entrez_s]

    sym = normalize_symbol(symbol)
    if sym and sym in symbol_map:
        return symbol_map[sym]

    alias_s = clean_text(aliases)
    if alias_s:
        for part in re.split(r"[|,;]+", alias_s):
            s = normalize_symbol(part)
            if s and s in symbol_map:
                return symbol_map[s]

    return ""


def is_human_interaction(row: pd.Series, cols: Dict[str, Optional[str]]) -> bool:
    oa = clean_text(row.get(cols.get("organism_a"), "")) if cols.get("organism_a") else ""
    ob = clean_text(row.get(cols.get("organism_b"), "")) if cols.get("organism_b") else ""

    # BioGRID organism fields are usually taxonomy IDs, e.g. 9606.
    if not oa and not ob:
        return True

    return oa == HUMAN_TAX_ID and ob == HUMAN_TAX_ID


def build_biogrid_long_table(
    biogrid: pd.DataFrame,
    symbol_map: Dict[str, str],
    entrez_map: Dict[str, str],
    processed_dir: Path,
    physical_only: bool = False,
    human_only: bool = True,
) -> pd.DataFrame:
    log("=" * 100)
    log("[BUILD BIOGRID LONG TABLE]")

    cols = detect_biogrid_columns(biogrid)

    required = ["symbol_a", "symbol_b"]
    if cols.get("symbol_a") is None or cols.get("symbol_b") is None:
        raise RuntimeError("Could not detect BioGRID Official Symbol Interactor A/B columns.")

    rows = []
    skipped_nonhuman = 0
    skipped_nonphysical = 0
    skipped_unmapped = 0
    skipped_self = 0

    for idx, row in biogrid.iterrows():
        if human_only and not is_human_interaction(row, cols):
            skipped_nonhuman += 1
            continue

        exp_type = clean_text(row.get(cols.get("experimental_system_type"), "")) if cols.get("experimental_system_type") else ""
        exp_type_norm = exp_type.lower()

        if physical_only and exp_type_norm != "physical":
            skipped_nonphysical += 1
            continue

        gene_a = map_interactor_to_gene(
            symbol=row.get(cols.get("symbol_a"), "") if cols.get("symbol_a") else "",
            entrez=row.get(cols.get("entrez_a"), "") if cols.get("entrez_a") else "",
            aliases=row.get(cols.get("aliases_a"), "") if cols.get("aliases_a") else "",
            symbol_map=symbol_map,
            entrez_map=entrez_map,
        )

        gene_b = map_interactor_to_gene(
            symbol=row.get(cols.get("symbol_b"), "") if cols.get("symbol_b") else "",
            entrez=row.get(cols.get("entrez_b"), "") if cols.get("entrez_b") else "",
            aliases=row.get(cols.get("aliases_b"), "") if cols.get("aliases_b") else "",
            symbol_map=symbol_map,
            entrez_map=entrez_map,
        )

        if not gene_a or not gene_b:
            skipped_unmapped += 1
            continue

        exp_system = clean_text(row.get(cols.get("experimental_system"), "")) if cols.get("experimental_system") else ""
        pubmed = clean_text(row.get(cols.get("pubmed_id"), "")) if cols.get("pubmed_id") else ""
        throughput = clean_text(row.get(cols.get("throughput"), "")) if cols.get("throughput") else ""
        score = clean_text(row.get(cols.get("score"), "")) if cols.get("score") else ""
        source_db = clean_text(row.get(cols.get("source_database"), "")) if cols.get("source_database") else ""
        author = clean_text(row.get(cols.get("author"), "")) if cols.get("author") else ""
        interaction_id = clean_text(row.get(cols.get("biogrid_id_a"), "")) if cols.get("biogrid_id_a") else ""

        is_self = int(gene_a == gene_b)

        # Keep self-interactions but mark them. They can be biologically meaningful.
        if gene_a == gene_b:
            skipped_self += 1

        pair = tuple(sorted([gene_a, gene_b]))

        base = {
            "biogrid_interaction_id": interaction_id,
            "gene_a": gene_a,
            "gene_b": gene_b,
            "gene_pair_sorted": f"{pair[0]}--{pair[1]}",
            "is_self_interaction": is_self,
            "experimental_system": exp_system,
            "experimental_system_type": exp_type,
            "throughput": throughput,
            "pubmed_id": pubmed,
            "author": author,
            "score": score,
            "source_database": source_db,
        }

        # Store two directional rows so aggregation per gene is easy.
        rows.append(
            {
                **base,
                "gene_symbol": gene_a,
                "partner_gene_symbol": gene_b,
                "direction": "A_to_B",
            }
        )
        rows.append(
            {
                **base,
                "gene_symbol": gene_b,
                "partner_gene_symbol": gene_a,
                "direction": "B_to_A",
            }
        )

        if (idx + 1) % 500000 == 0:
            log(f"[BIOGRID LONG] processed={idx + 1:,} long_rows={len(rows):,}")

    long_df = pd.DataFrame(rows)

    if not long_df.empty:
        long_df = long_df.drop_duplicates()

    outpath = processed_dir / "feature15_biogrid_interaction_long.csv"
    long_df.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED BIOGRID LONG]")
    log(f"[PATH] {outpath}")
    log(f"[SHAPE] {long_df.shape}")
    log(f"[GENES] {long_df['gene_symbol'].nunique() if not long_df.empty else 0}")
    log(f"[SKIPPED NONHUMAN] {skipped_nonhuman}")
    log(f"[SKIPPED NONPHYSICAL] {skipped_nonphysical}")
    log(f"[SKIPPED UNMAPPED] {skipped_unmapped}")
    log(f"[SELF INTERACTION ROWS BEFORE DIRECTIONAL DUP] {skipped_self}")

    return long_df


# =============================================================================
# GRAPH / AGGREGATION
# =============================================================================

def compute_local_clustering_for_gene(gene: str, partners: Set[str], undirected_edges: Set[Tuple[str, str]]) -> float:
    """
    Local clustering coefficient among known partners.
    For many genes this is cheap. For huge hubs, cap pair enumeration.
    """
    n = len(partners)
    if n < 2:
        return 0.0

    # Avoid insane O(k^2) for huge hubs.
    if n > 1000:
        return np.nan

    partners_list = sorted(partners)
    possible = n * (n - 1) / 2
    observed = 0

    for i in range(n):
        a = partners_list[i]
        for j in range(i + 1, n):
            b = partners_list[j]
            pair = tuple(sorted([a, b]))
            if pair in undirected_edges:
                observed += 1

    return float(observed / possible) if possible > 0 else 0.0


def aggregate_biogrid_features(
    hgnc: pd.DataFrame,
    long_df: pd.DataFrame,
    processed_dir: Path,
) -> pd.DataFrame:
    log("=" * 100)
    log("[AGGREGATE BIOGRID FEATURES PER GENE]")

    genes = sorted(hgnc["gene_symbol"].unique().tolist())

    if long_df.empty:
        features = pd.DataFrame({"gene_symbol": genes})
        features["feature15_biogrid_has_interaction"] = 0
        features["feature15_biogrid_interaction_count"] = 0
        features.to_csv(processed_dir / "feature15_biogrid_gene_features.csv", index=False)
        return features

    # Undirected edges for graph summaries.
    edges = set()
    for _, row in long_df[["gene_a", "gene_b"]].drop_duplicates().iterrows():
        a = normalize_symbol(row["gene_a"])
        b = normalize_symbol(row["gene_b"])
        if a and b and a != b:
            edges.add(tuple(sorted([a, b])))

    by_gene = dict(tuple(long_df.groupby("gene_symbol")))

    records = []

    for gene in genes:
        g = by_gene.get(gene, pd.DataFrame())

        if g.empty:
            records.append(
                {
                    "gene_symbol": gene,
                    "feature15_biogrid_has_interaction": 0,
                    "feature15_biogrid_interaction_count": 0,
                    "feature15_biogrid_degree": 0,
                    "feature15_biogrid_physical_interaction_count": 0,
                    "feature15_biogrid_genetic_interaction_count": 0,
                }
            )
            continue

        partners = set(g["partner_gene_symbol"].dropna().astype(str).map(normalize_symbol).tolist())
        partners.discard("")
        partners_no_self = {p for p in partners if p != gene}

        exp_types = g["experimental_system_type"].dropna().astype(str).tolist() if "experimental_system_type" in g.columns else []
        exp_systems = g["experimental_system"].dropna().astype(str).tolist() if "experimental_system" in g.columns else []
        throughputs = g["throughput"].dropna().astype(str).tolist() if "throughput" in g.columns else []
        pubmeds = g["pubmed_id"].replace("", np.nan).dropna().astype(str).tolist() if "pubmed_id" in g.columns else []
        source_dbs = g["source_database"].replace("", np.nan).dropna().astype(str).tolist() if "source_database" in g.columns else []

        exp_type_lower = [x.lower() for x in exp_types]
        throughput_lower = [x.lower() for x in throughputs]

        physical_count = sum(1 for x in exp_type_lower if x == "physical")
        genetic_count = sum(1 for x in exp_type_lower if x == "genetic")

        low_throughput_count = sum(1 for x in throughput_lower if "low" in x)
        high_throughput_count = sum(1 for x in throughput_lower if "high" in x)

        self_count = int(pd.to_numeric(g["is_self_interaction"], errors="coerce").fillna(0).sum()) if "is_self_interaction" in g.columns else 0

        rec = {
            "gene_symbol": gene,
            "feature15_biogrid_has_interaction": 1,
            "feature15_biogrid_interaction_count": int(len(g)),
            "feature15_biogrid_unique_interaction_pair_count": int(g["gene_pair_sorted"].nunique()) if "gene_pair_sorted" in g.columns else 0,
            "feature15_biogrid_degree": int(len(partners_no_self)),
            "feature15_biogrid_self_interaction_count": self_count,
            "feature15_biogrid_has_self_interaction": int(self_count > 0),

            "feature15_biogrid_physical_interaction_count": int(physical_count),
            "feature15_biogrid_genetic_interaction_count": int(genetic_count),
            "feature15_biogrid_physical_fraction": float(physical_count / len(g)) if len(g) else np.nan,
            "feature15_biogrid_genetic_fraction": float(genetic_count / len(g)) if len(g) else np.nan,

            "feature15_biogrid_low_throughput_count": int(low_throughput_count),
            "feature15_biogrid_high_throughput_count": int(high_throughput_count),
            "feature15_biogrid_low_throughput_fraction": float(low_throughput_count / len(g)) if len(g) else np.nan,
            "feature15_biogrid_high_throughput_fraction": float(high_throughput_count / len(g)) if len(g) else np.nan,

            "feature15_biogrid_unique_pubmed_count": int(pd.Series(pubmeds).nunique()) if pubmeds else 0,
            "feature15_biogrid_unique_experimental_system_count": int(pd.Series(exp_systems).nunique()) if exp_systems else 0,
            "feature15_biogrid_unique_experimental_system_type_count": int(pd.Series(exp_types).nunique()) if exp_types else 0,
            "feature15_biogrid_unique_throughput_count": int(pd.Series(throughputs).nunique()) if throughputs else 0,
            "feature15_biogrid_unique_source_database_count": int(pd.Series(source_dbs).nunique()) if source_dbs else 0,

            "feature15_biogrid_experimental_system_entropy": shannon_entropy(exp_systems),
            "feature15_biogrid_experimental_system_type_entropy": shannon_entropy(exp_types),
            "feature15_biogrid_throughput_entropy": shannon_entropy(throughputs),
            "feature15_biogrid_partner_entropy": shannon_entropy(list(partners_no_self)),

            "feature15_biogrid_local_clustering_coefficient": compute_local_clustering_for_gene(
                gene=gene,
                partners=partners_no_self,
                undirected_edges=edges,
            ),
        }

        # Common experimental-system keyword counts.
        exp_text = " | ".join(exp_systems).lower()

        keyword_groups = {
            "two_hybrid": ["two-hybrid", "two hybrid", "2-hybrid"],
            "affinity_capture": ["affinity capture", "pull down", "pulldown"],
            "co_fractionation": ["co-fractionation", "cofractionation"],
            "co_crystal_structure": ["co-crystal", "crystal structure", "x-ray"],
            "biochemical_activity": ["biochemical activity"],
            "synthetic_lethality": ["synthetic lethality"],
            "dosage_rescue": ["dosage rescue"],
            "phenotypic_enhancement": ["phenotypic enhancement"],
            "phenotypic_suppression": ["phenotypic suppression"],
            "reconstituted_complex": ["reconstituted complex"],
        }

        for group, kws in keyword_groups.items():
            count = 0
            for x in exp_systems:
                xl = x.lower()
                if any(kw in xl for kw in kws):
                    count += 1
            rec[f"feature15_biogrid_expsys_{group}_count"] = int(count)
            rec[f"feature15_biogrid_expsys_{group}_fraction"] = float(count / len(g)) if len(g) else np.nan

        records.append(rec)

    features = pd.DataFrame(records)

    # Derived transformed features.
    for col in [
        "feature15_biogrid_interaction_count",
        "feature15_biogrid_unique_interaction_pair_count",
        "feature15_biogrid_degree",
        "feature15_biogrid_unique_pubmed_count",
        "feature15_biogrid_unique_experimental_system_count",
    ]:
        if col in features.columns:
            features[f"{col}_log1p"] = np.log1p(pd.to_numeric(features[col], errors="coerce").fillna(0))

    if "feature15_biogrid_degree" in features.columns and "feature15_biogrid_unique_pubmed_count" in features.columns:
        features["feature15_biogrid_evidence_weighted_degree"] = (
            np.log1p(pd.to_numeric(features["feature15_biogrid_degree"], errors="coerce").fillna(0))
            * np.log1p(pd.to_numeric(features["feature15_biogrid_unique_pubmed_count"], errors="coerce").fillna(0))
        )

    # Fill count/flag columns.
    for c in features.columns:
        if c.startswith("feature15_") and (
            "_count" in c
            or "_has_" in c
            or c.endswith("_degree")
            or c.endswith("_log1p")
            or c.endswith("_weighted_degree")
        ):
            if features[c].dtype != object:
                features[c] = features[c].fillna(0)

    outpath = processed_dir / "feature15_biogrid_gene_features.csv"
    features.to_csv(outpath, index=False)

    log("[SAVED BIOGRID GENE FEATURES]")
    log(f"[PATH] {outpath}")
    log(f"[SHAPE] {features.shape}")
    log(f"[COVERAGE] {features['feature15_biogrid_has_interaction'].mean():.4f}")

    return features


# =============================================================================
# MERGE / SUMMARY
# =============================================================================

def merge_with_hgnc(hgnc: pd.DataFrame, features: pd.DataFrame, processed_dir: Path) -> pd.DataFrame:
    merged = hgnc.merge(features, on="gene_symbol", how="left")

    feature_cols = [c for c in features.columns if c != "gene_symbol"]

    merged["feature15_biogrid_has_any_feature"] = merged[feature_cols].notna().any(axis=1).astype(int)
    merged["feature15_biogrid_n_nonmissing_features"] = merged[feature_cols].notna().sum(axis=1)

    outpath = processed_dir / "feature15_biogrid_hgnc_merged.csv"
    merged.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED HGNC MERGED]")
    log(f"[PATH] {outpath}")
    log(f"[SHAPE] {merged.shape}")

    return merged


def write_summary(
    outdir: Path,
    hgnc: pd.DataFrame,
    biogrid_file: Path,
    biogrid: pd.DataFrame,
    long_df: pd.DataFrame,
    features: pd.DataFrame,
    merged: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    path = outdir / "feature15_biogrid_summary.txt"

    lines = []
    lines.append("Feature 15: BioGRID curated interaction-network features")
    lines.append("=" * 80)
    lines.append(f"Created at: {now_iso()}")
    lines.append("")
    lines.append("Inputs:")
    lines.append(f"HGNC: {args.hgnc}")
    lines.append(f"BioGRID file: {biogrid_file}")
    lines.append(f"BioGRID download URL: {BIOGRID_ORGANISM_TAB3_URL}")
    lines.append("")
    lines.append("Options:")
    lines.append(f"human_only: {not args.include_nonhuman}")
    lines.append(f"physical_only: {args.physical_only}")
    lines.append(f"max_rows: {args.max_rows}")
    lines.append("")
    lines.append("Rows:")
    lines.append(f"HGNC rows: {hgnc.shape[0]}")
    lines.append(f"BioGRID rows read: {biogrid.shape[0]}")
    lines.append(f"BioGRID directional long rows: {long_df.shape[0]}")
    lines.append(f"Genes with BioGRID interactions: {int(features['feature15_biogrid_has_interaction'].sum())}")
    lines.append(f"BioGRID coverage: {features['feature15_biogrid_has_interaction'].mean():.4f}")
    lines.append(f"Gene feature rows: {features.shape[0]}")
    lines.append(f"Gene feature columns: {features.shape[1]}")
    lines.append(f"Merged rows: {merged.shape[0]}")
    lines.append(f"Merged columns: {merged.shape[1]}")
    lines.append("")
    lines.append("Feature summaries:")

    for col in [
        "feature15_biogrid_interaction_count",
        "feature15_biogrid_degree",
        "feature15_biogrid_unique_pubmed_count",
        "feature15_biogrid_physical_interaction_count",
        "feature15_biogrid_genetic_interaction_count",
        "feature15_biogrid_evidence_weighted_degree",
        "feature15_biogrid_local_clustering_coefficient",
    ]:
        if col in features.columns:
            x = pd.to_numeric(features[col], errors="coerce")
            lines.append(
                f"{col}: median={x.median(skipna=True):.4f}, "
                f"max={x.max(skipna=True):.4f}, "
                f"nonmissing={int(x.notna().sum())}"
            )

    lines.append("")
    lines.append("Interpretation:")
    lines.append("These features capture curated experimental interaction burden and local network context.")
    lines.append("They are not known-drug or clinical-target labels.")
    lines.append("")
    lines.append("Leakage policy:")
    lines.append("Included: human physical/genetic interaction topology, publications, experimental systems, throughput.")
    lines.append("Excluded: chemical interaction files, drug names, ChEMBL, DrugBank, DGIdb, Open Targets, Pharos, clinical target labels.")

    path.write_text("\n".join(lines) + "\n")
    log(f"[SUMMARY] {path}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description="Build Feature 15 BioGRID interaction-network features.")
    parser.add_argument("--hgnc", default=str(DEFAULT_HGNC_PATH), help="Path to HGNC complete set.")
    parser.add_argument("--dbdir", default=str(DEFAULT_DBDIR), help="BioGRID database directory.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory.")
    parser.add_argument("--biogrid-file", default="", help="Manual BioGRID .tab3.txt or .zip file.")
    parser.add_argument("--download", action="store_true", help="Download BIOGRID-ORGANISM-LATEST.tab3.zip.")
    parser.add_argument("--force-download", action="store_true", help="Force re-download / re-extract.")
    parser.add_argument("--physical-only", action="store_true", help="Use only physical interactions. Default uses physical + genetic.")
    parser.add_argument("--include-nonhuman", action="store_true", help="Do not filter to taxon 9606. Not recommended.")
    parser.add_argument("--max-rows", type=int, default=0, help="Debug only: read first N BioGRID rows.")
    parser.add_argument("--limit-genes", type=int, default=0, help="Debug only: first N HGNC genes.")
    parser.add_argument("--all-hgnc-genes", action="store_true", help="Use all HGNC genes, not only protein-coding.")
    args = parser.parse_args()

    dbdir = mkdir(Path(args.dbdir))
    outdir = mkdir(Path(args.outdir))
    downloads_dir = mkdir(outdir / "downloads")
    processed_dir = mkdir(outdir / "processed")
    extract_dir = mkdir(dbdir / "extracted")

    default_zip = dbdir / "BIOGRID-ORGANISM-LATEST.tab3.zip"
    out_download_zip = downloads_dir / "BIOGRID-ORGANISM-LATEST.tab3.zip"

    log("=" * 100)
    log("FEATURE 15: BIOGRID CURATED INTERACTION FEATURES")
    log("=" * 100)
    log(f"[HGNC]             {args.hgnc}")
    log(f"[DBDIR]            {dbdir.resolve()}")
    log(f"[OUTDIR]           {outdir.resolve()}")
    log(f"[BIOGRID FILE]     {args.biogrid_file if args.biogrid_file else 'auto'}")
    log(f"[DOWNLOAD]         {args.download}")
    log(f"[FORCE DOWNLOAD]   {args.force_download}")
    log(f"[PHYSICAL ONLY]    {args.physical_only}")
    log(f"[HUMAN ONLY]       {not args.include_nonhuman}")
    log(f"[MAX ROWS]         {args.max_rows if args.max_rows else 'none'}")
    log(f"[LIMIT GENES]      {args.limit_genes if args.limit_genes else 'none'}")
    log("=" * 100)

    # -------------------------------------------------------------------------
    # BioGRID input detection / automatic download
    # -------------------------------------------------------------------------
    if args.biogrid_file:
        biogrid_input = Path(args.biogrid_file)
        if not biogrid_input.exists():
            raise FileNotFoundError(f"Provided --biogrid-file does not exist: {biogrid_input}")
    else:
        biogrid_input = find_local_biogrid_file(dbdir)

        # Fully automatic mode: download if no usable local BioGRID file exists.
        if biogrid_input is None or not biogrid_input.exists():
            log("=" * 100)
            log("[BIOGRID] No local BioGRID file found. Attempting automatic download.")
            download_file(BIOGRID_ORGANISM_TAB3_URL, default_zip, force=args.force_download)
            biogrid_input = default_zip

        # Explicit download mode: refresh/reuse the default ZIP.
        elif args.download or args.force_download:
            log("=" * 100)
            log("[BIOGRID] Local BioGRID file found, but --download/--force-download was requested.")
            download_file(BIOGRID_ORGANISM_TAB3_URL, default_zip, force=args.force_download)
            biogrid_input = default_zip

    if biogrid_input is None or not biogrid_input.exists():
        raise FileNotFoundError(
            "BioGRID file not found and automatic download failed. Provide --biogrid-file manually."
        )

    if biogrid_input.suffix.lower() == ".zip" and not is_valid_zip(biogrid_input):
        raise RuntimeError(f"BioGRID ZIP is invalid or incomplete: {biogrid_input}")

    # Copy raw zip/input for reproducibility.
    try:
        if biogrid_input.suffix.lower() == ".zip":
            if not out_download_zip.exists() or args.force_download:
                shutil.copy2(biogrid_input, out_download_zip)
    except Exception:
        pass

    # If ZIP, extract Homo sapiens tab3.
    if biogrid_input.suffix.lower() == ".zip":
        biogrid_file = extract_human_tab3_from_zip(
            zip_path=biogrid_input,
            extract_dir=extract_dir,
            force=args.force_download,
        )
    else:
        biogrid_file = biogrid_input

    hgnc = load_hgnc(
        Path(args.hgnc),
        protein_coding_only=not args.all_hgnc_genes,
        limit_genes=args.limit_genes,
    )

    symbol_map, entrez_map = build_hgnc_maps(hgnc)

    biogrid = read_biogrid_tab3(
        path=biogrid_file,
        max_rows=args.max_rows,
    )

    long_df = build_biogrid_long_table(
        biogrid=biogrid,
        symbol_map=symbol_map,
        entrez_map=entrez_map,
        processed_dir=processed_dir,
        physical_only=args.physical_only,
        human_only=not args.include_nonhuman,
    )

    features = aggregate_biogrid_features(
        hgnc=hgnc,
        long_df=long_df,
        processed_dir=processed_dir,
    )

    merged = merge_with_hgnc(
        hgnc=hgnc,
        features=features,
        processed_dir=processed_dir,
    )

    metadata = {
        "created_at": now_iso(),
        "hgnc": str(Path(args.hgnc)),
        "dbdir": str(dbdir.resolve()),
        "outdir": str(outdir.resolve()),
        "biogrid_input": str(biogrid_input),
        "biogrid_file_used": str(biogrid_file),
        "biogrid_download_url": BIOGRID_ORGANISM_TAB3_URL,
        "human_only": not args.include_nonhuman,
        "physical_only": args.physical_only,
        "max_rows": args.max_rows,
        "limit_genes": args.limit_genes,
        "protein_coding_only": not args.all_hgnc_genes,
        "outputs": {
            "interaction_long": str(processed_dir / "feature15_biogrid_interaction_long.csv"),
            "gene_features": str(processed_dir / "feature15_biogrid_gene_features.csv"),
            "hgnc_merged": str(processed_dir / "feature15_biogrid_hgnc_merged.csv"),
            "summary": str(outdir / "feature15_biogrid_summary.txt"),
        },
        "leakage_policy": {
            "included": [
                "human protein/genetic interaction counts",
                "degree",
                "publication counts",
                "experimental system counts",
                "throughput counts",
                "local clustering coefficient",
                "self-interaction flag",
            ],
            "excluded": [
                "BioGRID chemical interaction files",
                "drug names",
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

    with open(outdir / "feature15_biogrid_run_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    write_summary(
        outdir=outdir,
        hgnc=hgnc,
        biogrid_file=biogrid_file,
        biogrid=biogrid,
        long_df=long_df,
        features=features,
        merged=merged,
        args=args,
    )

    log("=" * 100)
    log("[DONE]")
    log(f"[BIOGRID LONG]  {processed_dir / 'feature15_biogrid_interaction_long.csv'}")
    log(f"[GENE FEATURES] {processed_dir / 'feature15_biogrid_gene_features.csv'}")
    log(f"[HGNC MERGED]   {processed_dir / 'feature15_biogrid_hgnc_merged.csv'}")
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