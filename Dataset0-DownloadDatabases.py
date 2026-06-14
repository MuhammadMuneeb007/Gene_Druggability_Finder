#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Dataset0-DownloadDatabases.py

Safe downloader for human gene druggability / tractability databases.

Downloads:
  1. HGNC
  2. ChEMBL SQLite
  3. Open Targets parquet datasets
  4. DGIdb TSV files
  5. TCRD / Pharos if discoverable, otherwise writes manual instructions

Important behaviour:
  - Safe to rerun.
  - Existing files are skipped unless --overwrite is used.
  - Existing extracted ChEMBL database is skipped unless --overwrite is used.
  - Existing Open Targets parquet folders are skipped unless --overwrite is used.
  - Open Targets directory recursion is restricted so it cannot climb to parent folders.

Recommended:
  conda create -n drugdb python=3.10 -y
  conda activate drugdb
  pip install requests tqdm

Run:
  python Dataset0-DownloadDatabases.py --all --opentargets-release 26.03
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import shutil
import sys
import tarfile
import time
import zipfile
from html.parser import HTMLParser
from pathlib import Path
from typing import Dict, List, Sequence, Tuple
from urllib.parse import urljoin, urlparse

import requests
from tqdm import tqdm


# =============================================================================
# URLs
# =============================================================================

HGNC_COMPLETE_TSV_URLS = [
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
    "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt",
]

CHEMBL_FTP_ROOT = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/"
OPENTARGETS_FTP_ROOT = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/"

TCRD_DOWNLOAD_ROOTS = [
    "http://juniper.health.unm.edu/tcrd/",
    "https://juniper.health.unm.edu/tcrd/",
]

DGIDB_DOWNLOAD_PAGE = "https://dgidb.org/downloads"

# For Open Targets 26.03, knownDrugs has effectively moved away from the older knownDrugs dataset.
# clinical_indication is now more appropriate for drug-disease clinical relationships.
DEFAULT_OPENTARGETS_DATASETS = [
    "target",
    "tractability",
    "clinical_indication",
]

OPENTARGETS_ALIASES = {
    "targets": "target",
    "target": "target",
    "tractability": "tractability",
    "knownDrugs": "clinical_indication",
    "known_drugs": "clinical_indication",
    "known_drug": "clinical_indication",
    "clinicalIndication": "clinical_indication",
    "clinical_indications": "clinical_indication",
    "clinical_indication": "clinical_indication",
    "drug": "drug",
    "drugs": "drug",
}

DGIDB_DEFAULT_FILES = [
    "interactions.tsv",
    "genes.tsv",
    "drugs.tsv",
    "categories.tsv",
]


# =============================================================================
# Generic helpers
# =============================================================================

class LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links: List[str] = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() != "a":
            return
        for k, v in attrs:
            if k.lower() == "href" and v:
                self.links.append(v)


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def now_string() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def save_json(path: Path, data: dict) -> None:
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)


def write_text(path: Path, text: str) -> None:
    ensure_dir(path.parent)
    path.write_text(text.strip() + "\n", encoding="utf-8")


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            block = f.read(block_size)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def get_html(url: str, timeout: int = 60) -> str:
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    return r.text


def get_links(url: str, timeout: int = 60) -> List[str]:
    html = get_html(url, timeout=timeout)
    parser = LinkParser()
    parser.feed(html)
    return parser.links


def file_exists_and_big_enough(path: Path, min_bytes: int = 1) -> bool:
    return path.exists() and path.is_file() and path.stat().st_size >= min_bytes


def any_file_exists(root: Path, patterns: Sequence[str], min_bytes: int = 1) -> bool:
    if not root.exists():
        return False

    for pattern in patterns:
        for p in root.rglob(pattern):
            if p.is_file() and p.stat().st_size >= min_bytes:
                return True

    return False


def is_probably_html_response(path: Path) -> bool:
    try:
        with open(path, "rb") as f:
            head = f.read(500).lower()
        return b"<html" in head or b"<!doctype html" in head
    except Exception:
        return False


def same_or_child_url(candidate_url: str, root_url: str) -> bool:
    """
    True only if candidate_url stays inside root_url.

    This prevents recursion from following Parent Directory links back to:
      https://ftp.ebi.ac.uk/pub/
    """
    candidate = urlparse(candidate_url)
    root = urlparse(root_url)

    if candidate.scheme != root.scheme:
        return False

    if candidate.netloc != root.netloc:
        return False

    root_path = root.path
    if not root_path.endswith("/"):
        root_path += "/"

    return candidate.path.startswith(root_path)


def download_file(
    url: str,
    out_path: Path,
    timeout: int = 180,
    chunk_size: int = 1024 * 1024,
    overwrite: bool = False,
    min_bytes: int = 1,
    allow_html: bool = False,
) -> Dict[str, object]:
    ensure_dir(out_path.parent)

    if file_exists_and_big_enough(out_path, min_bytes=min_bytes) and not overwrite:
        print(f"[SKIP] Existing file: {out_path}")
        return {
            "url": url,
            "path": str(out_path),
            "status": "exists",
            "bytes": out_path.stat().st_size,
            "sha256": sha256_file(out_path),
        }

    part_path = out_path.with_suffix(out_path.suffix + ".part")

    if overwrite and part_path.exists():
        part_path.unlink()

    headers = {}
    existing_size = 0

    if part_path.exists() and not overwrite:
        existing_size = part_path.stat().st_size
        if existing_size > 0:
            headers["Range"] = f"bytes={existing_size}-"

    with requests.get(url, stream=True, timeout=timeout, headers=headers) as r:
        if r.status_code == 200 and existing_size > 0:
            existing_size = 0
            part_path.unlink(missing_ok=True)

        r.raise_for_status()

        total = r.headers.get("content-length")
        total = int(total) + existing_size if total is not None else None

        mode = "ab" if existing_size > 0 else "wb"

        with open(part_path, mode) as f:
            with tqdm(
                total=total,
                initial=existing_size,
                unit="B",
                unit_scale=True,
                desc=out_path.name[:50],
            ) as pbar:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)
                        pbar.update(len(chunk))

    if part_path.stat().st_size < min_bytes:
        raise RuntimeError(f"Downloaded file is too small: {part_path}")

    part_path.replace(out_path)

    if not allow_html and is_probably_html_response(out_path):
        raise RuntimeError(f"Downloaded HTML instead of data file: {url} -> {out_path}")

    return {
        "url": url,
        "path": str(out_path),
        "status": "downloaded",
        "bytes": out_path.stat().st_size,
        "sha256": sha256_file(out_path),
    }


def extract_archive(path: Path, dest_dir: Path, overwrite: bool = False) -> Dict[str, object]:
    ensure_dir(dest_dir)
    marker = dest_dir / ".extracted.ok"

    if marker.exists() and not overwrite:
        print(f"[SKIP] Already extracted marker found: {dest_dir}")
        return {
            "archive": str(path),
            "dest_dir": str(dest_dir),
            "status": "already_extracted_marker",
        }

    if not overwrite:
        if any_file_exists(dest_dir, ["*.db", "*.sqlite", "*.sqlite3"], min_bytes=1024):
            print(f"[SKIP] Extracted database already exists: {dest_dir}")
            marker.write_text(now_string(), encoding="utf-8")
            return {
                "archive": str(path),
                "dest_dir": str(dest_dir),
                "status": "already_extracted_files_found",
            }

    print(f"[EXTRACT] {path} -> {dest_dir}")

    if path.name.endswith((".tar.gz", ".tgz")):
        with tarfile.open(path, "r:gz") as tar:
            tar.extractall(dest_dir)

    elif path.name.endswith(".tar"):
        with tarfile.open(path, "r:") as tar:
            tar.extractall(dest_dir)

    elif path.name.endswith(".zip"):
        with zipfile.ZipFile(path, "r") as z:
            z.extractall(dest_dir)

    elif path.name.endswith(".gz"):
        out_file = dest_dir / path.name[:-3]
        if out_file.exists() and not overwrite:
            print(f"[SKIP] Existing decompressed file: {out_file}")
        else:
            with gzip.open(path, "rb") as src, open(out_file, "wb") as dst:
                shutil.copyfileobj(src, dst)

    else:
        return {
            "archive": str(path),
            "dest_dir": str(dest_dir),
            "status": "not_archive",
        }

    marker.write_text(now_string(), encoding="utf-8")

    return {
        "archive": str(path),
        "dest_dir": str(dest_dir),
        "status": "extracted",
    }


# =============================================================================
# HGNC
# =============================================================================

def download_hgnc(base_dir: Path, overwrite: bool = False) -> Dict[str, object]:
    out_dir = base_dir / "HGNC"
    ensure_dir(out_dir)

    out_file = out_dir / "hgnc_complete_set.txt"

    if file_exists_and_big_enough(out_file, min_bytes=1000) and not overwrite:
        print(f"[HGNC] SKIP existing: {out_file}")
        meta = {
            "source": "HGNC",
            "status": "exists",
            "path": str(out_file),
            "bytes": out_file.stat().st_size,
            "sha256": sha256_file(out_file),
        }
        save_json(out_dir / "download_metadata.json", meta)
        return meta

    errors = []

    for url in HGNC_COMPLETE_TSV_URLS:
        try:
            meta = download_file(url, out_file, overwrite=overwrite, min_bytes=1000)
            meta["source"] = "HGNC"
            save_json(out_dir / "download_metadata.json", meta)
            return meta
        except Exception as e:
            errors.append({"url": url, "error": repr(e)})
            print(f"[HGNC] Failed: {url} | {e}")

    save_json(out_dir / "download_errors.json", {"errors": errors})
    raise RuntimeError("HGNC download failed.")


# =============================================================================
# ChEMBL
# =============================================================================

def discover_latest_chembl_release() -> Tuple[str, str]:
    links = get_links(CHEMBL_FTP_ROOT)
    releases = []

    for link in links:
        m = re.match(r"chembl_(\d+)/?$", link)
        if m:
            releases.append((int(m.group(1)), link))

    if not releases:
        raise RuntimeError("Could not discover ChEMBL releases.")

    version, link = sorted(releases)[-1]
    release_name = f"chembl_{version}"
    release_url = urljoin(CHEMBL_FTP_ROOT, link if link.endswith("/") else link + "/")

    return release_name, release_url


def discover_chembl_sqlite_url(release_url: str, release_name: str) -> str:
    links = get_links(release_url)

    candidates = []
    for link in links:
        lower = link.lower()
        if "sqlite" in lower and lower.endswith((".tar.gz", ".tgz", ".gz", ".zip")):
            candidates.append(link)

    exact = [c for c in candidates if release_name.lower() in c.lower()]

    if exact:
        return urljoin(release_url, exact[0])

    if candidates:
        return urljoin(release_url, candidates[0])

    raise RuntimeError(f"Could not find ChEMBL SQLite download in {release_url}")


def download_chembl(
    base_dir: Path,
    release: str = "latest",
    overwrite: bool = False,
    extract: bool = True,
) -> Dict[str, object]:
    out_dir = base_dir / "ChEMBL"
    ensure_dir(out_dir)

    if release == "latest":
        release_name, release_url = discover_latest_chembl_release()
    else:
        release_name = release if release.startswith("chembl_") else f"chembl_{release}"
        release_url = urljoin(CHEMBL_FTP_ROOT, release_name + "/")

    sqlite_url = discover_chembl_sqlite_url(release_url, release_name)
    file_name = sqlite_url.rstrip("/").split("/")[-1]

    out_file = out_dir / file_name
    extract_dir = out_dir / "extracted"

    print(f"[ChEMBL] Release: {release_name}")
    print(f"[ChEMBL] SQLite URL: {sqlite_url}")

    if (
        file_exists_and_big_enough(out_file, min_bytes=1024 * 1024)
        and any_file_exists(extract_dir, ["*.db", "*.sqlite", "*.sqlite3"], min_bytes=1024)
        and not overwrite
    ):
        print("[ChEMBL] SKIP existing archive and extracted SQLite database.")
        meta = {
            "source": "ChEMBL",
            "release": release_name,
            "status": "exists_downloaded_and_extracted",
            "archive": str(out_file),
            "extract_dir": str(extract_dir),
            "bytes": out_file.stat().st_size,
            "sha256": sha256_file(out_file),
        }
        save_json(out_dir / "download_metadata.json", meta)
        return meta

    meta = download_file(sqlite_url, out_file, overwrite=overwrite, min_bytes=1024 * 1024)

    meta.update(
        {
            "source": "ChEMBL",
            "release": release_name,
            "release_url": release_url,
            "sqlite_url": sqlite_url,
        }
    )

    if extract:
        meta["extract"] = extract_archive(out_file, extract_dir, overwrite=overwrite)
    else:
        meta["extract"] = {"status": "not_requested"}

    save_json(out_dir / "download_metadata.json", meta)
    return meta


# =============================================================================
# Open Targets
# =============================================================================

def normalise_opentargets_dataset_name(dataset: str) -> str:
    return OPENTARGETS_ALIASES.get(dataset, dataset)


def discover_latest_opentargets_release() -> Tuple[str, str]:
    links = get_links(OPENTARGETS_FTP_ROOT)
    releases = []

    for link in links:
        name = link.strip("/")
        if re.match(r"^\d{2}\.\d{2}$", name):
            year, month = name.split(".")
            releases.append((int(year), int(month), name, link))

    if not releases:
        raise RuntimeError("Could not discover Open Targets release directories.")

    _, _, name, link = sorted(releases)[-1]
    release_url = urljoin(OPENTARGETS_FTP_ROOT, link if link.endswith("/") else link + "/")

    return name, release_url


def opentargets_uses_new_layout(release_name: str) -> bool:
    try:
        y, m = release_name.split(".")
        return (int(y), int(m)) >= (25, 3)
    except Exception:
        return True


def local_dataset_has_parquet(dataset_out: Path) -> bool:
    if not dataset_out.exists():
        return False

    parquet_files = list(dataset_out.rglob("*.parquet"))
    parquet_files = [p for p in parquet_files if p.exists() and p.stat().st_size > 0]

    return len(parquet_files) > 0


def safe_list_directory(url: str) -> List[str]:
    """
    Return links from a remote FTP HTML page, excluding parent-directory links.
    """
    raw_links = get_links(url)
    clean_links = []

    for href in raw_links:
        href = href.strip()

        if not href:
            continue

        low = href.lower()

        if href in {"../", "./"}:
            continue

        if "parent directory" in low:
            continue

        if href.startswith("?"):
            continue

        if href.startswith("#"):
            continue

        # Avoid icons and Apache sorting links.
        if low.endswith((".gif", ".png", ".ico")):
            continue

        if "c=n;o=" in low or "c=m;o=" in low or "c=s;o=" in low or "c=d;o=" in low:
            continue

        full = urljoin(url, href)

        # Critical fix: never follow links outside this dataset URL.
        if not same_or_child_url(full, url):
            continue

        clean_links.append(href)

    return clean_links


def download_opentargets_directory(
    dataset_root_url: str,
    current_url: str,
    local_dir: Path,
    overwrite: bool = False,
    current_depth: int = 0,
    max_depth: int = 8,
) -> List[Dict[str, object]]:
    """
    Recursively download Open Targets files, but only inside dataset_root_url.
    """
    ensure_dir(local_dir)

    if current_depth > max_depth:
        raise RuntimeError(f"Maximum recursion depth reached at {current_url}")

    links = safe_list_directory(current_url)
    downloaded = []

    for href in links:
        full_url = urljoin(current_url, href)

        # Extra safety.
        if not same_or_child_url(full_url, dataset_root_url):
            continue

        clean_name = href.split("?")[0].strip("/")

        if not clean_name:
            continue

        if href.endswith("/"):
            subdir = local_dir / clean_name
            downloaded.extend(
                download_opentargets_directory(
                    dataset_root_url=dataset_root_url,
                    current_url=full_url,
                    local_dir=subdir,
                    overwrite=overwrite,
                    current_depth=current_depth + 1,
                    max_depth=max_depth,
                )
            )
            continue

        lower = clean_name.lower()

        keep_file = (
            lower.endswith(".parquet")
            or lower == "_success"
            or lower == "_metadata"
            or lower == "_common_metadata"
            or lower == "readme.txt"
            or lower.endswith(".json")
        )

        if not keep_file:
            continue

        out_file = local_dir / clean_name

        meta = download_file(
            url=full_url,
            out_path=out_file,
            overwrite=overwrite,
            min_bytes=0 if lower == "_success" else 1,
            allow_html=False,
        )
        downloaded.append(meta)

    return downloaded


def get_opentargets_candidate_urls(
    release_url: str,
    release_name: str,
    dataset: str,
) -> List[str]:
    d = normalise_opentargets_dataset_name(dataset)

    urls = []

    if d == "tractability":
        # Current useful tractability path is often here.
        urls.extend(
            [
                urljoin(release_url, "input/target/tractability/"),
                urljoin(release_url, "output/target/tractability/"),
                urljoin(release_url, "output/tractability/"),
                urljoin(release_url, "output/etl/parquet/tractability/"),
            ]
        )

    elif d == "target":
        urls.extend(
            [
                urljoin(release_url, "output/target/"),
                urljoin(release_url, "output/etl/parquet/targets/"),
                urljoin(release_url, "output/etl/parquet/target/"),
            ]
        )

    elif d == "clinical_indication":
        urls.extend(
            [
                urljoin(release_url, "output/clinical_indication/"),
                urljoin(release_url, "output/clinicalIndication/"),
                urljoin(release_url, "output/known_drug/"),
                urljoin(release_url, "output/knownDrugs/"),
                urljoin(release_url, "output/etl/parquet/knownDrugs/"),
                urljoin(release_url, "output/etl/parquet/known_drug/"),
            ]
        )

    else:
        if opentargets_uses_new_layout(release_name):
            urls.extend(
                [
                    urljoin(release_url, f"output/{d}/"),
                    urljoin(release_url, f"input/{d}/"),
                    urljoin(release_url, f"output/etl/parquet/{d}/"),
                ]
            )
        else:
            urls.extend(
                [
                    urljoin(release_url, f"output/etl/parquet/{d}/"),
                    urljoin(release_url, f"output/{d}/"),
                    urljoin(release_url, f"input/{d}/"),
                ]
            )

    # Deduplicate.
    seen = set()
    final = []
    for u in urls:
        if u not in seen:
            final.append(u)
            seen.add(u)

    return final


def download_opentargets_dataset(
    release_url: str,
    release_name: str,
    dataset: str,
    out_dir: Path,
    overwrite: bool = False,
) -> Dict[str, object]:
    d = normalise_opentargets_dataset_name(dataset)
    dataset_out = out_dir / d
    ensure_dir(dataset_out)

    print(f"[OpenTargets] Dataset requested: {dataset}")
    print(f"[OpenTargets] Dataset normalised: {d}")
    print(f"[OpenTargets] Local output: {dataset_out}")

    if local_dataset_has_parquet(dataset_out) and not overwrite:
        parquet_files = list(dataset_out.rglob("*.parquet"))
        print(f"[OpenTargets] SKIP existing dataset: {d} ({len(parquet_files)} parquet files)")
        meta = {
            "source": "OpenTargets",
            "dataset_requested": dataset,
            "dataset": d,
            "status": "exists",
            "local_dir": str(dataset_out),
            "num_existing_parquet_files": len(parquet_files),
        }
        save_json(dataset_out / "download_metadata.json", meta)
        return meta

    candidate_urls = get_opentargets_candidate_urls(release_url, release_name, dataset)
    errors = []

    for dataset_url in candidate_urls:
        print(f"[OpenTargets] Trying URL: {dataset_url}")

        try:
            # Test listing first.
            links = safe_list_directory(dataset_url)
            if not links:
                raise RuntimeError(f"No usable links found at {dataset_url}")

            downloaded = download_opentargets_directory(
                dataset_root_url=dataset_url,
                current_url=dataset_url,
                local_dir=dataset_out,
                overwrite=overwrite,
                max_depth=8,
            )

            parquet_files = list(dataset_out.rglob("*.parquet"))

            if len(parquet_files) == 0:
                save_json(
                    dataset_out / "listing_no_parquet.json",
                    {
                        "dataset_url": dataset_url,
                        "links": links,
                        "downloaded": downloaded,
                    },
                )
                raise RuntimeError(f"No parquet files downloaded from {dataset_url}")

            meta = {
                "source": "OpenTargets",
                "dataset_requested": dataset,
                "dataset": d,
                "status": "downloaded",
                "dataset_url": dataset_url,
                "local_dir": str(dataset_out),
                "num_local_parquet_files": len(parquet_files),
                "num_downloaded_or_existing_files": len(downloaded),
                "downloaded": downloaded,
                "errors": errors,
            }
            save_json(dataset_out / "download_metadata.json", meta)
            print(f"[OpenTargets] SUCCESS: {d} ({len(parquet_files)} parquet files)")
            return meta

        except Exception as e:
            errors.append({"dataset_url": dataset_url, "error": repr(e)})
            print(f"[OpenTargets] URL failed: {dataset_url} | {e}")

    meta = {
        "source": "OpenTargets",
        "dataset_requested": dataset,
        "dataset": d,
        "status": "failed",
        "candidate_urls": candidate_urls,
        "errors": errors,
    }
    save_json(dataset_out / "download_metadata_failed.json", meta)
    raise RuntimeError(f"Open Targets dataset failed: {dataset} / {d}")


def download_opentargets(
    base_dir: Path,
    release: str = "latest",
    datasets: Sequence[str] = DEFAULT_OPENTARGETS_DATASETS,
    overwrite: bool = False,
) -> Dict[str, object]:
    out_dir = base_dir / "OpenTargets"
    ensure_dir(out_dir)

    if release == "latest":
        release_name, release_url = discover_latest_opentargets_release()
    else:
        release_name = release
        release_url = urljoin(OPENTARGETS_FTP_ROOT, release_name + "/")

    print(f"[OpenTargets] Release: {release_name}")
    print(f"[OpenTargets] Release URL: {release_url}")

    release_out = out_dir / release_name
    ensure_dir(release_out)

    metas = []

    for dataset in datasets:
        try:
            metas.append(
                download_opentargets_dataset(
                    release_url=release_url,
                    release_name=release_name,
                    dataset=dataset,
                    out_dir=release_out,
                    overwrite=overwrite,
                )
            )
        except Exception as e:
            print(f"[OpenTargets] Dataset failed: {dataset} | {e}")
            metas.append(
                {
                    "source": "OpenTargets",
                    "dataset_requested": dataset,
                    "dataset": normalise_opentargets_dataset_name(dataset),
                    "status": "failed",
                    "error": repr(e),
                }
            )

    meta = {
        "source": "OpenTargets",
        "release": release_name,
        "release_url": release_url,
        "datasets_requested": list(datasets),
        "datasets": metas,
    }
    save_json(release_out / "download_metadata.json", meta)
    return meta


# =============================================================================
# DGIdb
# =============================================================================

def download_dgidb(base_dir: Path, overwrite: bool = False) -> Dict[str, object]:
    out_dir = base_dir / "DGIdb"
    ensure_dir(out_dir)

    if any_file_exists(out_dir, ["*.tsv", "*.tsv.gz", "*.csv", "*.csv.gz"], min_bytes=100) and not overwrite:
        print(f"[DGIdb] SKIP existing data file in: {out_dir}")
        meta = {
            "source": "DGIdb",
            "status": "exists",
            "local_dir": str(out_dir),
        }
        save_json(out_dir / "download_metadata.json", meta)
        return meta

    html = get_html(DGIDB_DOWNLOAD_PAGE)
    links = get_links(DGIDB_DOWNLOAD_PAGE)

    candidate_urls = []

    # First, collect real href links from the page.
    for href in links:
        low = href.lower()
        if any(name.lower() in low for name in DGIDB_DEFAULT_FILES):
            candidate_urls.append(urljoin(DGIDB_DOWNLOAD_PAGE, href))

    # If the page uses predictable relative routes but parser misses them,
    # try common DGIdb latest routes.
    fallback_urls = []
    for fname in DGIDB_DEFAULT_FILES:
        fallback_urls.extend(
            [
                f"https://dgidb.org/data/latest/{fname}",
                f"https://dgidb.org/downloads/latest/{fname}",
                f"https://dgidb.org/downloads/2024-Dec/{fname}",
                f"https://dgidb.org/data/2024-Dec/{fname}",
            ]
        )

    candidate_urls.extend(fallback_urls)

    # Deduplicate.
    seen = set()
    candidate_urls = [u for u in candidate_urls if not (u in seen or seen.add(u))]

    downloaded = []
    errors = []

    for url in candidate_urls:
        fname = url.rstrip("/").split("/")[-1]
        if fname not in DGIDB_DEFAULT_FILES:
            # Keep file names clean.
            continue

        out_file = out_dir / fname

        try:
            meta = download_file(
                url=url,
                out_path=out_file,
                overwrite=overwrite,
                min_bytes=100,
                allow_html=False,
            )
            downloaded.append(meta)
        except Exception as e:
            errors.append({"url": url, "error": repr(e)})

    if downloaded:
        meta = {
            "source": "DGIdb",
            "status": "downloaded",
            "local_dir": str(out_dir),
            "downloaded": downloaded,
            "num_downloaded": len(downloaded),
            "errors": errors,
        }
        save_json(out_dir / "download_metadata.json", meta)
        return meta

    instructions = """
DGIdb automatic download did not find working direct URLs.

Manual download:
  Open:
    https://dgidb.org/downloads

Download these latest files:
  interactions.tsv
  genes.tsv
  drugs.tsv
  categories.tsv

Place them here:
  databases/DGIdb/

Your pasted DGIdb page confirms these files are available for latest 2024-Dec.
"""
    write_text(out_dir / "MANUAL_DOWNLOAD_REQUIRED.txt", instructions)

    meta = {
        "source": "DGIdb",
        "status": "manual_download_required",
        "candidate_urls_tried": candidate_urls,
        "errors": errors,
    }
    save_json(out_dir / "download_metadata.json", meta)
    return meta


# =============================================================================
# TCRD / Pharos
# =============================================================================

def discover_tcrd_downloads(root_url: str) -> List[str]:
    links = get_links(root_url, timeout=60)
    candidates = []

    for link in links:
        low = link.lower()
        if any(x in low for x in ["tcrd", "pharos"]) and low.endswith(
            (".gz", ".zip", ".sql", ".dump", ".tar.gz", ".tgz", ".sql.gz", ".dump.gz")
        ):
            candidates.append(urljoin(root_url, link))

    return candidates


def download_tcrd(base_dir: Path, overwrite: bool = False) -> Dict[str, object]:
    out_dir = base_dir / "TCRD_Pharos"
    ensure_dir(out_dir)

    if any_file_exists(
        out_dir,
        ["*.sql", "*.sql.gz", "*.dump", "*.dump.gz", "*.db", "*.sqlite", "*.zip", "*.tar.gz", "*.tgz"],
        min_bytes=1024,
    ) and not overwrite:
        print(f"[TCRD/Pharos] SKIP existing file in: {out_dir}")
        meta = {
            "source": "TCRD_Pharos",
            "status": "exists",
            "local_dir": str(out_dir),
        }
        save_json(out_dir / "download_metadata.json", meta)
        return meta

    errors = []
    candidates = []

    for root in TCRD_DOWNLOAD_ROOTS:
        try:
            candidates.extend(discover_tcrd_downloads(root))
        except Exception as e:
            errors.append({"url": root, "error": repr(e)})
            print(f"[TCRD/Pharos] Discovery failed: {root} | {e}")

    candidates = sorted(set(candidates))

    if candidates:
        chosen = candidates[-1]
        fname = chosen.rstrip("/").split("/")[-1]
        out_file = out_dir / fname

        try:
            meta = download_file(chosen, out_file, overwrite=overwrite, min_bytes=1024)
            meta.update(
                {
                    "source": "TCRD_Pharos",
                    "status": "downloaded",
                    "chosen": chosen,
                    "candidates": candidates,
                    "errors": errors,
                }
            )
            save_json(out_dir / "download_metadata.json", meta)
            return meta
        except Exception as e:
            errors.append({"url": chosen, "error": repr(e)})

    instructions = """
TCRD / Pharos could not be downloaded automatically.

This is expected on many systems because the older public TCRD dump location:
  http://juniper.health.unm.edu/tcrd/

may be unavailable, slow, blocked, or no longer maintained as a stable bulk
download endpoint.

Recommended treatment:
  - Treat TCRD/Pharos as optional/manual.
  - Continue with HGNC + ChEMBL + Open Targets + DGIdb.
  - If TDL labels are required, manually download a TCRD dump and place it here:
      databases/TCRD_Pharos/

The merge script should not fail if TCRD is absent.
"""
    write_text(out_dir / "MANUAL_DOWNLOAD_REQUIRED.txt", instructions)

    meta = {
        "source": "TCRD_Pharos",
        "status": "manual_download_required",
        "errors": errors,
        "candidates": candidates,
    }
    save_json(out_dir / "download_metadata.json", meta)
    return meta


# =============================================================================
# Main
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Download druggability benchmark databases.")

    p.add_argument("--outdir", default="databases")
    p.add_argument(
        "--sources",
        nargs="+",
        default=[],
        choices=["hgnc", "chembl", "opentargets", "dgidb", "tcrd"],
    )
    p.add_argument("--all", action="store_true")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--no-extract", action="store_true")
    p.add_argument("--chembl-release", default="latest")
    p.add_argument("--opentargets-release", default="latest")
    p.add_argument(
        "--opentargets-datasets",
        nargs="+",
        default=DEFAULT_OPENTARGETS_DATASETS,
    )

    return p.parse_args()


def main() -> None:
    args = parse_args()

    base_dir = Path(args.outdir)
    ensure_dir(base_dir)

    sources = set(args.sources)

    if args.all:
        sources.update(["hgnc", "chembl", "opentargets", "dgidb", "tcrd"])

    if not sources:
        print("No sources selected. Use --all or --sources hgnc chembl opentargets dgidb tcrd")
        sys.exit(1)

    print("=" * 90)
    print("DRUGGABILITY DATABASE DOWNLOADER")
    print("=" * 90)
    print(f"Output directory: {base_dir.resolve()}")
    print(f"Sources: {', '.join(sorted(sources))}")
    print(f"Overwrite existing files: {args.overwrite}")
    print("=" * 90)

    run_meta = {
        "started": now_string(),
        "output_directory": str(base_dir.resolve()),
        "sources": sorted(sources),
        "overwrite": args.overwrite,
        "results": {},
    }

    if "hgnc" in sources:
        try:
            run_meta["results"]["hgnc"] = download_hgnc(base_dir, overwrite=args.overwrite)
        except Exception as e:
            run_meta["results"]["hgnc"] = {"status": "failed", "error": repr(e)}
            print(f"[ERROR] HGNC failed: {e}")

    if "chembl" in sources:
        try:
            run_meta["results"]["chembl"] = download_chembl(
                base_dir,
                release=args.chembl_release,
                overwrite=args.overwrite,
                extract=not args.no_extract,
            )
        except Exception as e:
            run_meta["results"]["chembl"] = {"status": "failed", "error": repr(e)}
            print(f"[ERROR] ChEMBL failed: {e}")

    if "opentargets" in sources:
        try:
            run_meta["results"]["opentargets"] = download_opentargets(
                base_dir,
                release=args.opentargets_release,
                datasets=args.opentargets_datasets,
                overwrite=args.overwrite,
            )
        except Exception as e:
            run_meta["results"]["opentargets"] = {"status": "failed", "error": repr(e)}
            print(f"[ERROR] Open Targets failed: {e}")

    if "dgidb" in sources:
        try:
            run_meta["results"]["dgidb"] = download_dgidb(base_dir, overwrite=args.overwrite)
        except Exception as e:
            run_meta["results"]["dgidb"] = {"status": "failed", "error": repr(e)}
            print(f"[ERROR] DGIdb failed: {e}")

    if "tcrd" in sources:
        try:
            run_meta["results"]["tcrd"] = download_tcrd(base_dir, overwrite=args.overwrite)
        except Exception as e:
            run_meta["results"]["tcrd"] = {"status": "failed", "error": repr(e)}
            print(f"[ERROR] TCRD/Pharos failed: {e}")

    run_meta["finished"] = now_string()
    save_json(base_dir / "download_run_metadata.json", run_meta)

    print("\n" + "=" * 90)
    print("SOURCE SUMMARY")
    print("=" * 90)

    for source, meta in run_meta["results"].items():
        status = meta.get("status", "completed") if isinstance(meta, dict) else "completed"
        print(f"{source:15s} {status}")

    print("=" * 90)
    print("DOWNLOAD STEP COMPLETE")
    print("=" * 90)
    print(f"Run metadata: {base_dir / 'download_run_metadata.json'}")
    print("Next step:")
    print("  Use the merge script to create:")
    print("  Step0_Output/03_HumanGene_DruggabilityLabels.csv")
    print("=" * 90)


if __name__ == "__main__":
    main()