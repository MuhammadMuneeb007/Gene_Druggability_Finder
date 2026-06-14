#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature1_DepMap.py

Download selected DepMap release files using the Figshare public API.

No merging.
No feature generation.
No processing.

Run:
    python Feature1_DepMap.py

Optional:
    python Feature1_DepMap.py --release "DepMap Public 26Q1"
    python Feature1_DepMap.py --include-optional
"""

import argparse
import time
from pathlib import Path

import requests


# =============================================================================
# SETTINGS
# =============================================================================

DEFAULT_RELEASE = "DepMap Public 26Q1"
OUTDIR = Path("Feature1_DepMap") / "downloads"

REQUIRED_FILES = [
    "README.txt",
    "Gene.csv",
    "Model.csv",
    "CRISPRGeneEffect.csv",
    "CRISPRGeneDependency.csv",
    "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
]

OPTIONAL_FILES = [
    "OmicsSomaticMutationsMatrixHotspot.csv",
    "OmicsSomaticMutationsMatrixDamaging.csv",
    "OmicsCNGeneWGS.csv",
    "OmicsFusionFiltered.csv",
    "SubtypeMatrix.csv",
    "SubtypeTree.csv",
]

DO_NOT_DOWNLOAD = [
    "PortalCompounds.csv",
]

FIGSHARE_API = "https://api.figshare.com/v2"


# =============================================================================
# API HELPERS
# =============================================================================

def make_session():
    s = requests.Session()
    s.headers.update({
        "User-Agent": "Mozilla/5.0 DepMapDownloader/1.0",
        "Accept": "application/json",
    })
    return s


def api_get(session, url, params=None, retries=5):
    for attempt in range(1, retries + 1):
        try:
            r = session.get(url, params=params, timeout=120)
            if r.status_code == 200:
                return r.json()

            print(f"[API WARNING] HTTP {r.status_code}: {url}")
            print(str(r.text)[:500])

        except Exception as e:
            print(f"[API WARNING] Attempt {attempt}/{retries} failed: {e}")

        if attempt < retries:
            time.sleep(3 * attempt)

    raise RuntimeError(f"API request failed: {url}")


def search_figshare_articles(session, release_name):
    """
    Search Figshare for the DepMap release article.
    """
    print("=" * 100)
    print("[API] Searching Figshare articles")
    print(f"[QUERY] {release_name}")

    queries = [
        release_name,
        release_name.replace("DepMap Public ", "DepMap "),
        release_name.replace("Public ", ""),
    ]

    all_hits = []

    for q in queries:
        url = f"{FIGSHARE_API}/articles"
        params = {
            "search_for": q,
            "page_size": 100,
            "order": "published_date",
            "order_direction": "desc",
        }

        hits = api_get(session, url, params=params)

        for h in hits:
            title = h.get("title", "")
            if "depmap" in title.lower():
                all_hits.append(h)

    # Deduplicate by article id
    dedup = {}
    for h in all_hits:
        dedup[h["id"]] = h

    hits = list(dedup.values())

    print(f"[API] Candidate articles found: {len(hits)}")
    for h in hits[:20]:
        print(f"  - ID={h.get('id')} | {h.get('title')}")

    return hits


def get_article_details(session, article_id):
    url = f"{FIGSHARE_API}/articles/{article_id}"
    return api_get(session, url)


def get_article_files(session, article_id):
    url = f"{FIGSHARE_API}/articles/{article_id}/files"
    return api_get(session, url)


def find_release_article_with_files(session, release_name, wanted_files):
    """
    Find the Figshare article that contains the requested DepMap files.
    """
    candidates = search_figshare_articles(session, release_name)

    if not candidates:
        raise RuntimeError(
            f"No Figshare article candidates found for release: {release_name}"
        )

    best = None
    best_score = -1
    best_files = None

    wanted_set = set(wanted_files)

    print("=" * 100)
    print("[API] Checking candidate article files")

    for c in candidates:
        article_id = c["id"]
        title = c.get("title", "")

        try:
            files = get_article_files(session, article_id)
        except Exception as e:
            print(f"[SKIP] Could not fetch files for article {article_id}: {e}")
            continue

        names = {f.get("name") for f in files}
        score = len(wanted_set.intersection(names))

        print(f"[CHECK] ID={article_id} | matched {score}/{len(wanted_set)} | {title}")

        if score > best_score:
            best = c
            best_score = score
            best_files = files

    if best is None or best_score == 0:
        raise RuntimeError(
            "Could not find a Figshare article containing the requested DepMap files. "
            "The release name may be different on Figshare."
        )

    print("=" * 100)
    print("[SELECTED ARTICLE]")
    print(f"[ID]    {best['id']}")
    print(f"[TITLE] {best.get('title')}")
    print(f"[MATCH] {best_score}/{len(wanted_set)} requested files found")

    return best, best_files


# =============================================================================
# DOWNLOAD
# =============================================================================

def download_url(session, url, outpath, expected_size=None, retries=5):
    if outpath.exists() and outpath.stat().st_size > 0:
        print(f"[SKIP] Already exists: {outpath}")
        print(f"       Size: {outpath.stat().st_size / 1024 / 1024:.2f} MB")
        return True

    tmp = outpath.with_suffix(outpath.suffix + ".tmp")

    for attempt in range(1, retries + 1):
        try:
            print("=" * 100)
            print(f"[DOWNLOAD] {outpath.name}")
            print(f"[SAVE TO]  {outpath}")

            with session.get(url, stream=True, timeout=300) as r:
                if r.status_code != 200:
                    print(f"[HTTP] {r.status_code}")
                    print(str(r.text)[:500])
                    raise RuntimeError(f"HTTP {r.status_code}")

                total = int(r.headers.get("content-length", 0))
                if total == 0 and expected_size:
                    total = int(expected_size)

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

            print(f"[DONE] {outpath.name}")
            print(f"[SIZE] {outpath.stat().st_size / 1024 / 1024:.2f} MB")
            return True

        except Exception as e:
            print()
            print(f"[FAILED] Attempt {attempt}/{retries}: {outpath.name}")
            print(f"[ERROR]  {e}")

            if tmp.exists():
                tmp.unlink()

            if attempt < retries:
                wait = 5 * attempt
                print(f"[RETRY] Waiting {wait} seconds...")
                time.sleep(wait)

    return False


def main():
    parser = argparse.ArgumentParser(
        description="Download selected DepMap files using Figshare API."
    )

    parser.add_argument(
        "--release",
        default=DEFAULT_RELEASE,
        help=f'DepMap release name. Default: "{DEFAULT_RELEASE}"',
    )

    parser.add_argument(
        "--include-optional",
        action="store_true",
        help="Also download optional mutation/CN/fusion/subtype files.",
    )

    parser.add_argument(
        "--outdir",
        default=str(OUTDIR),
        help=f"Output directory. Default: {OUTDIR}",
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    files_to_download = list(REQUIRED_FILES)

    if args.include_optional:
        files_to_download.extend(OPTIONAL_FILES)

    files_to_download = [f for f in files_to_download if f not in DO_NOT_DOWNLOAD]

    print("=" * 100)
    print("DEPMAP DOWNLOAD ONLY USING FIGSHARE API")
    print("=" * 100)
    print(f"[RELEASE] {args.release}")
    print(f"[OUTDIR]  {outdir.resolve()}")
    print("[FILES TO DOWNLOAD]")
    for f in files_to_download:
        print(f"  - {f}")

    print()
    print("[NOT DOWNLOADING - LEAKAGE RISK]")
    for f in DO_NOT_DOWNLOAD:
        print(f"  - {f}")

    session = make_session()

    article, article_files = find_release_article_with_files(
        session=session,
        release_name=args.release,
        wanted_files=files_to_download,
    )

    file_map = {f["name"]: f for f in article_files}

    print("=" * 100)
    print("[AVAILABLE MATCHED FILES]")
    for fname in files_to_download:
        if fname in file_map:
            f = file_map[fname]
            size_mb = f.get("size", 0) / 1024 / 1024
            print(f"  [FOUND]   {fname} ({size_mb:.2f} MB)")
        else:
            print(f"  [MISSING] {fname}")

    failed = []
    missing = []

    for fname in files_to_download:
        if fname not in file_map:
            missing.append(fname)
            continue

        f = file_map[fname]
        download_link = f.get("download_url")

        if not download_link:
            print(f"[NO DOWNLOAD URL] {fname}")
            failed.append(fname)
            continue

        ok = download_url(
            session=session,
            url=download_link,
            outpath=outdir / fname,
            expected_size=f.get("size"),
        )

        if not ok:
            failed.append(fname)

    print("=" * 100)
    print("[FINISHED]")
    print(f"[LOCATION] {outdir.resolve()}")

    if missing:
        print("[MISSING FROM FIGSHARE ARTICLE]")
        for f in missing:
            print(f"  - {f}")

    if failed:
        print("[FAILED DOWNLOADS]")
        for f in failed:
            print(f"  - {f}")

    if not missing and not failed:
        print("[ALL REQUESTED FILES DOWNLOADED]")

    print("=" * 100)


if __name__ == "__main__":
    main()