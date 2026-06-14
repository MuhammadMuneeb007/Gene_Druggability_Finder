 
"""
Step0B_Download_NoLeakage_Feature_Databases.py

Purpose
-------
Download LOCAL biological/protein feature databases for the next step of the
Gene Druggability Finder pipeline, WITHOUT downloading druggability-label
databases that would cause data leakage.

This script is for the feature-building step after you already have:
    Step0_Output/03_HumanGene_DruggabilityLabels.csv

The next pipeline step will read the genes from that labelled table and extract
independent biological/protein features locally.

What this script downloads
--------------------------
SAFE / NON-LABEL FEATURE SOURCES:
  1. HGNC complete gene set
     - official symbols, HGNC IDs, Ensembl IDs, Entrez IDs, UniProt IDs
  2. UniProt human proteome TSV
     - protein length, protein name, subcellular location, PDB IDs,
       InterPro/Pfam/PROSITE cross-references, RefSeq, Ensembl, STRING IDs
  3. UniProt human proteome FASTA
     - protein sequences for sequence-derived features
  4. Ensembl current human GTF and protein FASTA
     - genomic coordinates, transcripts, biotypes, protein sequences
  5. STRING human network files
     - protein links, physical links, protein aliases, protein info
  6. GTEx open expression files, if URL is available
     - tissue expression features
  7. AlphaFold human proteome structures, optional
     - predicted structures and pLDDT-derived structural confidence/disorder
  8. InterPro / Pfam bulk files, optional because some are very large
     - domain/family mapping if you want full offline domain annotation

What this script deliberately does NOT download
-----------------------------------------------
LEAKAGE / LABEL-EVIDENCE SOURCES:
  - ChEMBL
  - Open Targets tractability / knownDrugs
  - Pharos / TCRD target-development labels
  - DGIdb
  - DrugBank
  - Guide to Pharmacology drug-target evidence

Those sources should be used only for labels/evidence, not as predictive features
in a non-circular validation experiment.

Recommended installation
------------------------
    conda create -n drugfeatures python=3.10 -y
    conda activate drugfeatures
    pip install requests tqdm pandas

Quick test
----------
    python Step0B_Download_NoLeakage_Feature_Databases.py --minimal

Recommended feature download
----------------------------
    python Step0B_Download_NoLeakage_Feature_Databases.py --recommended

Full feature download, including large optional sources
------------------------------------------------------
    python Step0B_Download_NoLeakage_Feature_Databases.py --full

Custom source selection
-----------------------
    python Step0B_Download_NoLeakage_Feature_Databases.py --sources hgnc uniprot ensembl string gtex

Output
------
feature_databases/
    HGNC/
    UniProt/
    Ensembl/
    STRING/
    GTEx/
    AlphaFold/
    InterPro/
    MANIFEST_feature_sources.csv
    download_run_metadata.json
    DO_NOT_USE_AS_LABEL_FEATURES.txt
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import shutil
import sys
import tarfile
import time
import warnings
import zipfile
from dataclasses import dataclass, asdict
from html.parser import HTMLParser
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from urllib.parse import quote, urljoin

import pandas as pd
import requests
from requests.exceptions import RequestsDependencyWarning
from tqdm import tqdm

warnings.simplefilter("ignore", RequestsDependencyWarning)


# =============================================================================
# SOURCE URLS AND DEFAULTS
# =============================================================================

# HGNC
HGNC_COMPLETE_TSV_URLS = [
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
    "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt",
]

# UniProt
UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"
UNIPROT_HUMAN_PROTEOME = "UP000005640"

# Field list is intentionally feature-oriented and does not include druggability labels.
UNIPROT_FIELDS_FULL = [
    "accession",
    "id",
    "gene_names",
    "gene_primary",
    "protein_name",
    "organism_name",
    "length",
    "sequence",
    "cc_subcellular_location",
    "xref_pdb",
    "xref_alphafolddb",
    "xref_interpro",
    "xref_pfam",
    "xref_prosite",
    "xref_refseq",
    "xref_ensembl",
    "xref_string",
    "ec",
    "ft_domain",
]

# If UniProt rejects one of the above field names, use a smaller reliable set.
UNIPROT_FIELDS_FALLBACK = [
    "accession",
    "id",
    "gene_names",
    "gene_primary",
    "protein_name",
    "organism_name",
    "length",
    "cc_subcellular_location",
    "xref_pdb",
    "xref_interpro",
    "xref_pfam",
    "xref_refseq",
    "xref_ensembl",
    "xref_string",
    "ec",
]

# Ensembl
ENSEMBL_FTP_ROOT = "https://ftp.ensembl.org/pub/"
ENSEMBL_CURRENT_GTF_HUMAN = "https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/"
ENSEMBL_CURRENT_FASTA_HUMAN_PEP = "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/"
ENSEMBL_CURRENT_FASTA_HUMAN_CDNA = "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/"

# STRING
STRING_VERSION = "12.0"
STRING_BASE = "https://stringdb-downloads.org/download"
STRING_FILES = {
    "protein_info": f"{STRING_BASE}/protein.info.v{STRING_VERSION}/9606.protein.info.v{STRING_VERSION}.txt.gz",
    "protein_aliases": f"{STRING_BASE}/protein.aliases.v{STRING_VERSION}/9606.protein.aliases.v{STRING_VERSION}.txt.gz",
    "protein_links": f"{STRING_BASE}/protein.links.v{STRING_VERSION}/9606.protein.links.v{STRING_VERSION}.txt.gz",
    "protein_physical_links": f"{STRING_BASE}/protein.physical.links.v{STRING_VERSION}/9606.protein.physical.links.v{STRING_VERSION}.txt.gz",
    "protein_links_detailed": f"{STRING_BASE}/protein.links.detailed.v{STRING_VERSION}/9606.protein.links.detailed.v{STRING_VERSION}.txt.gz",
}

# AlphaFold
# Human proteome archive is large, so it is optional.
ALPHAFOLD_LATEST_DIR = "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/"
ALPHAFOLD_HUMAN_PROTEOME_PATTERNS = [
    r"^UP000005640_9606_HUMAN_v\d+\.tar$",
    r"^UP000005640_9606_HUMAN_v\d+\.tar\.gz$",
    r"^UP000005640_9606_HUMAN.*\.tar$",
    r"^UP000005640_9606_HUMAN.*\.tar\.gz$",
]

# InterPro / Pfam
INTERPRO_DOWNLOAD_ROOT = "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/"
INTERPRO_FILES = {
    # protein2ipr is huge. It maps proteins to InterPro entries.
    "protein2ipr": "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz",
    "entry_list": "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list",
    "interpro_xml": "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro.xml.gz",
    "parent_child_tree": "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/ParentChildTreeFile.txt",
}
PFAM_FILES = {
    # Pfam is maintained via InterPro/EBI. These files can be very large.
    "Pfam-A_hmm": "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
    "Pfam-A_full": "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz",
    "Pfam-A_regions": "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.regions.tsv.gz",
}

# GTEx: candidate open-access expression files. Some releases/URLs may change.
# The script tries candidates and writes instructions if not available.
GTEX_CANDIDATE_FILES = {
    "v8_gene_median_tpm": [
        "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
    ],
    "v8_gene_tpm": [
        "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz",
    ],
    "v10_gene_median_tpm": [
        "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/gene_tpm/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz",
        "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz",
    ],
}


LEAKAGE_DATABASES = [
    "ChEMBL",
    "Open Targets tractability",
    "Open Targets knownDrugs",
    "Pharos/TCRD TDL",
    "DGIdb",
    "DrugBank",
    "Guide to Pharmacology drug-target evidence",
]


# =============================================================================
# BASIC DOWNLOAD HELPERS
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


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            block = f.read(block_size)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def save_json(path: Path, data: dict) -> None:
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)


def save_text(path: Path, text: str) -> None:
    ensure_dir(path.parent)
    path.write_text(text.strip() + "\n", encoding="utf-8")


def get_links(url: str, timeout: int = 60) -> List[str]:
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    parser = LinkParser()
    parser.feed(r.text)
    return parser.links


def is_probably_html_response(path: Path) -> bool:
    try:
        with open(path, "rb") as f:
            head = f.read(500).lower()
        return b"<html" in head or b"<!doctype html" in head
    except Exception:
        return False


def find_existing_files(root: Path, patterns: Sequence[str], min_bytes: int = 1) -> List[Path]:
    if not root.exists():
        return []

    found: List[Path] = []
    seen = set()
    for pattern in patterns:
        for path in root.rglob(pattern):
            if not path.is_file():
                continue
            if path.stat().st_size < min_bytes:
                continue
            key = str(path.resolve()).lower()
            if key in seen:
                continue
            seen.add(key)
            found.append(path)

    return sorted(found)


def first_existing_file(root: Path, patterns: Sequence[str], min_bytes: int = 1) -> Optional[Path]:
    matches = find_existing_files(root, patterns, min_bytes=min_bytes)
    return matches[0] if matches else None


def existing_file_meta(
    path: Path,
    *,
    url: Optional[str] = None,
    source: Optional[str] = None,
    extra: Optional[Dict[str, object]] = None,
) -> Dict[str, object]:
    meta: Dict[str, object] = {
        "path": str(path),
        "status": "exists",
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }
    if url is not None:
        meta["url"] = url
    if source is not None:
        meta["source"] = source
    if extra:
        meta.update(extra)
    return meta


def http_head_ok(url: str, timeout: int = 30) -> bool:
    try:
        r = requests.head(url, allow_redirects=True, timeout=timeout)
        if r.status_code < 400:
            return True
        # Some servers do not support HEAD.
        r = requests.get(url, stream=True, timeout=timeout)
        ok = r.status_code < 400
        r.close()
        return ok
    except Exception:
        return False


def download_file(
    url: str,
    out_path: Path,
    timeout: int = 180,
    chunk_size: int = 1024 * 1024,
    overwrite: bool = False,
    min_bytes: int = 1,
    allow_html: bool = False,
) -> Dict[str, object]:
    """
    Streaming download with resume through a .part file.
    """
    ensure_dir(out_path.parent)

    if out_path.exists() and not overwrite and out_path.stat().st_size >= min_bytes:
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
                desc=out_path.name[:45],
            ) as pbar:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)
                        pbar.update(len(chunk))

    if part_path.stat().st_size < min_bytes:
        raise RuntimeError(f"Downloaded file too small: {part_path}")

    part_path.replace(out_path)

    if not allow_html and is_probably_html_response(out_path):
        raise RuntimeError(f"Downloaded HTML instead of data: {url}")

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
        return {"archive": str(path), "dest_dir": str(dest_dir), "status": "already_extracted"}

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
        with gzip.open(path, "rb") as src, open(out_file, "wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        return {"archive": str(path), "dest_dir": str(dest_dir), "status": "not_archive"}

    marker.write_text(now_string(), encoding="utf-8")
    return {"archive": str(path), "dest_dir": str(dest_dir), "status": "extracted"}


def write_leakage_warning(base_dir: Path) -> None:
    text = f"""
    DATA LEAKAGE WARNING

    This folder is intended for NON-LABEL biological/protein feature sources only.

    Do NOT place or use the following databases as model input features in
    non-circular validation:

    {chr(10).join("- " + x for x in LEAKAGE_DATABASES)}

    These sources can be used to create target labels/evidence, but they must be
    excluded from predictive feature matrices if you want to claim that the model
    recovers druggability from independent biological/protein features.

    Safe feature examples:
      - protein length
      - protein sequence
      - protein family/domain annotations
      - genomic coordinates
      - subcellular location
      - tissue expression
      - structure availability
      - AlphaFold confidence
      - binding-pocket features computed from structure
      - STRING network degree/connectivity

    Unsafe leakage examples:
      - approved drug count
      - known drug-target interactions
      - ChEMBL active compound counts
      - Open Targets small-molecule tractability labels
      - Pharos Tclin/Tchem/Tbio/Tdark as model features
    """
    save_text(base_dir / "DO_NOT_USE_AS_LABEL_FEATURES.txt", text)


# =============================================================================
# DOWNLOADERS
# =============================================================================

def download_hgnc(base_dir: Path, overwrite: bool = False) -> Dict[str, object]:
    out_dir = base_dir / "HGNC"
    ensure_dir(out_dir)
    out_file = out_dir / "hgnc_complete_set.txt"

    if out_file.exists() and not overwrite and out_file.stat().st_size >= 1000:
        meta = existing_file_meta(out_file, source="HGNC")
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
    return {"source": "HGNC", "status": "failed", "errors": errors}


def build_uniprot_stream_url(
    fields: Sequence[str],
    reviewed_only: bool = False,
    format_: str = "tsv",
    compressed: bool = True,
) -> str:
    query = f"proteome:{UNIPROT_HUMAN_PROTEOME}"
    if reviewed_only:
        query = f"({query}) AND (reviewed:true)"
    params = {
        "query": query,
        "format": format_,
        "fields": ",".join(fields) if format_ == "tsv" else None,
        "compressed": "true" if compressed else "false",
    }
    qs = []
    for k, v in params.items():
        if v is not None:
            qs.append(f"{quote(str(k))}={quote(str(v), safe=':,()')}")
    return UNIPROT_STREAM_URL + "?" + "&".join(qs)


def download_uniprot(base_dir: Path, overwrite: bool = False, reviewed_only: bool = False) -> Dict[str, object]:
    out_dir = base_dir / "UniProt"
    ensure_dir(out_dir)

    metas = []
    errors = []

    # 1. Human proteome TSV.
    existing_tsv = first_existing_file(
        out_dir,
        ["uniprot_human_proteome_full.tsv.gz", "uniprot_human_proteome_fallback.tsv.gz"],
        min_bytes=1000,
    )
    if existing_tsv and not overwrite:
        label = "full" if "full" in existing_tsv.name else "fallback"
        fields = UNIPROT_FIELDS_FULL if label == "full" else UNIPROT_FIELDS_FALLBACK
        metas.append(
            existing_file_meta(
                existing_tsv,
                source="UniProt",
                extra={
                    "type": "human_proteome_tsv",
                    "reviewed_only": reviewed_only,
                    "fields": list(fields),
                },
            )
        )
    else:
        for fields, label in [(UNIPROT_FIELDS_FULL, "full"), (UNIPROT_FIELDS_FALLBACK, "fallback")]:
            url = build_uniprot_stream_url(fields, reviewed_only=reviewed_only, format_="tsv", compressed=True)
            out_file = out_dir / f"uniprot_human_proteome_{label}.tsv.gz"
            try:
                meta = download_file(url, out_file, overwrite=overwrite, min_bytes=1000)
                meta.update({
                    "source": "UniProt",
                    "type": "human_proteome_tsv",
                    "reviewed_only": reviewed_only,
                    "fields": list(fields),
                })
                metas.append(meta)
                break
            except Exception as e:
                errors.append({"type": "tsv", "label": label, "url": url, "error": repr(e)})
                print(f"[UniProt] TSV {label} failed: {e}")

    # 2. Human proteome FASTA.
    fasta_out = out_dir / "uniprot_human_proteome.fasta.gz"
    if fasta_out.exists() and not overwrite and fasta_out.stat().st_size >= 1000:
        metas.append(
            existing_file_meta(
                fasta_out,
                source="UniProt",
                extra={
                    "type": "human_proteome_fasta",
                    "reviewed_only": reviewed_only,
                },
            )
        )
    else:
        fasta_url = build_uniprot_stream_url([], reviewed_only=reviewed_only, format_="fasta", compressed=True)
        try:
            meta = download_file(fasta_url, fasta_out, overwrite=overwrite, min_bytes=1000)
            meta.update({
                "source": "UniProt",
                "type": "human_proteome_fasta",
                "reviewed_only": reviewed_only,
            })
            metas.append(meta)
        except Exception as e:
            errors.append({"type": "fasta", "url": fasta_url, "error": repr(e)})
            print(f"[UniProt] FASTA failed: {e}")

    summary = {
        "source": "UniProt",
        "proteome": UNIPROT_HUMAN_PROTEOME,
        "reviewed_only": reviewed_only,
        "downloads": metas,
        "errors": errors,
    }
    save_json(out_dir / "download_metadata.json", summary)
    return summary


def discover_first_matching_file(directory_url: str, patterns: Sequence[str]) -> Optional[str]:
    links = get_links(directory_url)
    for pattern in patterns:
        regex = re.compile(pattern)
        for link in links:
            clean = link.split("/")[-1]
            if regex.search(clean):
                return urljoin(directory_url, link)
    return None


def discover_alphafold_human_archive_url() -> Optional[str]:
    links = get_links(ALPHAFOLD_LATEST_DIR)
    candidates: List[Tuple[int, int, str]] = []

    for link in links:
        clean = link.split("/")[-1]
        if not clean:
            continue

        if not any(re.match(pattern, clean) for pattern in ALPHAFOLD_HUMAN_PROTEOME_PATTERNS):
            continue

        version_match = re.search(r"_v(\d+)\.tar(?:\.gz)?$", clean)
        version = int(version_match.group(1)) if version_match else -1
        gz_score = 1 if clean.endswith(".tar.gz") else 0
        candidates.append((version, gz_score, clean))

    if not candidates:
        return None

    _, _, filename = sorted(candidates, reverse=True)[0]
    return urljoin(ALPHAFOLD_LATEST_DIR, filename)


def download_ensembl(base_dir: Path, overwrite: bool = False) -> Dict[str, object]:
    out_dir = base_dir / "Ensembl"
    ensure_dir(out_dir)

    tasks = [
        {
            "name": "gtf",
            "dir_url": ENSEMBL_CURRENT_GTF_HUMAN,
            "patterns": [r"Homo_sapiens\.GRCh38\.\d+\.gtf\.gz$"],
            "local_patterns": ["Homo_sapiens.GRCh38.*.gtf.gz"],
            "out_name": None,
        },
        {
            "name": "pep_all_fasta",
            "dir_url": ENSEMBL_CURRENT_FASTA_HUMAN_PEP,
            "patterns": [r"Homo_sapiens\.GRCh38\.pep\.all\.fa\.gz$"],
            "local_patterns": ["Homo_sapiens.GRCh38.pep.all.fa.gz"],
            "out_name": None,
        },
        {
            "name": "cdna_all_fasta",
            "dir_url": ENSEMBL_CURRENT_FASTA_HUMAN_CDNA,
            "patterns": [r"Homo_sapiens\.GRCh38\.cdna\.all\.fa\.gz$"],
            "local_patterns": ["Homo_sapiens.GRCh38.cdna.all.fa.gz"],
            "out_name": None,
        },
    ]

    metas = []
    errors = []
    for task in tasks:
        try:
            existing = None if overwrite else first_existing_file(
                out_dir,
                task["local_patterns"],
                min_bytes=1000,
            )
            if existing is not None:
                metas.append(
                    existing_file_meta(
                        existing,
                        source="Ensembl",
                        extra={
                            "type": task["name"],
                            "directory_url": task["dir_url"],
                        },
                    )
                )
                continue

            file_url = discover_first_matching_file(task["dir_url"], task["patterns"])
            if not file_url:
                raise RuntimeError(f"No matching file found at {task['dir_url']}")
            out_name = task["out_name"] or file_url.split("/")[-1]
            out_file = out_dir / out_name
            meta = download_file(file_url, out_file, overwrite=overwrite, min_bytes=1000)
            meta.update({
                "source": "Ensembl",
                "type": task["name"],
                "directory_url": task["dir_url"],
            })
            metas.append(meta)
        except Exception as e:
            errors.append({"task": task["name"], "error": repr(e), "dir_url": task["dir_url"]})
            print(f"[Ensembl] {task['name']} failed: {e}")

    summary = {"source": "Ensembl", "downloads": metas, "errors": errors}
    save_json(out_dir / "download_metadata.json", summary)
    return summary


def download_string(base_dir: Path, overwrite: bool = False, include_detailed: bool = False) -> Dict[str, object]:
    out_dir = base_dir / "STRING"
    ensure_dir(out_dir)

    selected = ["protein_info", "protein_aliases", "protein_links", "protein_physical_links"]
    if include_detailed:
        selected.append("protein_links_detailed")

    local_patterns = {
        "protein_info": ["9606.protein.info.v*.txt.gz"],
        "protein_aliases": ["9606.protein.aliases.v*.txt.gz"],
        "protein_links": ["9606.protein.links.v*.txt.gz"],
        "protein_physical_links": ["9606.protein.physical.links.v*.txt.gz"],
        "protein_links_detailed": ["9606.protein.links.detailed.v*.txt.gz"],
    }

    metas = []
    errors = []
    for name in selected:
        url = STRING_FILES[name]
        existing = None if overwrite else first_existing_file(out_dir, local_patterns[name], min_bytes=1000)
        if existing is not None:
            metas.append(
                existing_file_meta(
                    existing,
                    source="STRING",
                    extra={"type": name, "version": STRING_VERSION},
                )
            )
            continue

        out_file = out_dir / url.split("/")[-1]
        try:
            meta = download_file(url, out_file, overwrite=overwrite, min_bytes=1000)
            meta.update({"source": "STRING", "type": name, "version": STRING_VERSION})
            metas.append(meta)
        except Exception as e:
            errors.append({"type": name, "url": url, "error": repr(e)})
            print(f"[STRING] {name} failed: {e}")

    summary = {"source": "STRING", "version": STRING_VERSION, "downloads": metas, "errors": errors}
    save_json(out_dir / "download_metadata.json", summary)
    return summary


def download_gtex(base_dir: Path, overwrite: bool = False, include_full_tpm: bool = False) -> Dict[str, object]:
    out_dir = base_dir / "GTEx"
    ensure_dir(out_dir)

    selected_keys = ["v8_gene_median_tpm", "v10_gene_median_tpm"]
    if include_full_tpm:
        selected_keys.append("v8_gene_tpm")

    local_patterns = {
        "v8_gene_median_tpm": ["*v8*gene_median_tpm.gct.gz"],
        "v8_gene_tpm": ["*v8*gene_tpm.gct.gz"],
        "v10_gene_median_tpm": ["*v10*gene_median_tpm.gct.gz"],
    }

    metas = []
    errors = []
    for key in selected_keys:
        existing = None if overwrite else first_existing_file(
            out_dir,
            local_patterns.get(key, []),
            min_bytes=1000,
        )
        if existing is not None:
            metas.append(existing_file_meta(existing, source="GTEx", extra={"type": key}))
            continue

        downloaded = False
        for url in GTEX_CANDIDATE_FILES.get(key, []):
            out_file = out_dir / url.split("/")[-1]
            try:
                meta = download_file(url, out_file, overwrite=overwrite, min_bytes=1000)
                meta.update({"source": "GTEx", "type": key})
                metas.append(meta)
                downloaded = True
                break
            except Exception as e:
                errors.append({"type": key, "url": url, "error": repr(e)})
                print(f"[GTEx] candidate failed for {key}: {url} | {e}")
        if not downloaded:
            print(f"[GTEx] No candidate downloaded for {key}")

    if errors:
        save_text(
            out_dir / "GTEX_MANUAL_DOWNLOAD_NOTE.txt",
            """
            Some GTEx files could not be downloaded automatically.

            Go to the GTEx Portal open-access downloads page and download:
              - gene median TPM expression file
              - gene TPM matrix file, optional and large

            Put downloaded .gct.gz files in this folder:
              feature_databases/GTEx/

            For feature extraction, the gene median TPM file is usually enough.
            """
        )

    summary = {"source": "GTEx", "downloads": metas, "errors": errors}
    save_json(out_dir / "download_metadata.json", summary)
    return summary


def download_alphafold(base_dir: Path, overwrite: bool = False, extract: bool = False) -> Dict[str, object]:
    out_dir = base_dir / "AlphaFold"
    ensure_dir(out_dir)

    existing_archive = None if overwrite else first_existing_file(
        out_dir,
        ["UP000005640_9606_HUMAN*.tar", "UP000005640_9606_HUMAN*.tar.gz"],
        min_bytes=1024 * 1024,
    )
    if existing_archive is not None:
        meta = existing_file_meta(
            existing_archive,
            source="AlphaFold",
            extra={"type": "human_proteome_structures"},
        )
        if extract:
            meta["extract"] = extract_archive(existing_archive, out_dir / "extracted", overwrite=overwrite)
        save_json(out_dir / "download_metadata.json", meta)
        return meta

    errors = []
    archive_url = None
    try:
        archive_url = discover_alphafold_human_archive_url()
    except Exception as e:
        errors.append(
            {
                "stage": "discovery",
                "directory_url": ALPHAFOLD_LATEST_DIR,
                "error": repr(e),
            }
        )

    if archive_url:
        out_file = out_dir / archive_url.split("/")[-1]
        try:
            meta = download_file(archive_url, out_file, overwrite=overwrite, min_bytes=1024 * 1024)
            meta.update(
                {
                    "source": "AlphaFold",
                    "type": "human_proteome_structures",
                    "discovered_from": ALPHAFOLD_LATEST_DIR,
                }
            )
            if extract:
                meta["extract"] = extract_archive(out_file, out_dir / "extracted", overwrite=overwrite)
            save_json(out_dir / "download_metadata.json", meta)
            return meta
        except Exception as e:
            errors.append({"stage": "download", "url": archive_url, "error": repr(e)})
            print("[AlphaFold] Automatic download failed after archive discovery.")
    else:
        errors.append(
            {
                "stage": "discovery",
                "directory_url": ALPHAFOLD_LATEST_DIR,
                "error": "No matching human proteome archive found.",
            }
        )
        print("[AlphaFold] No matching human proteome archive was found in the latest AlphaFold directory.")

    save_text(
        out_dir / "ALPHAFOLD_MANUAL_DOWNLOAD_NOTE.txt",
        """
        Automatic AlphaFold human proteome download failed.

        Go to the official AlphaFold Protein Structure Database downloads page
        and download the Homo sapiens / Human / UP000005640 proteome archive.

        Put the archive in:
          feature_databases/AlphaFold/

        This source is optional at the download stage. You can also avoid downloading
        all AlphaFold structures and instead download structures only for selected
        genes later using UniProt accessions.
        """
    )
    summary = {"source": "AlphaFold", "status": "failed_or_manual_required", "errors": errors}
    save_json(out_dir / "download_metadata.json", summary)
    return summary


def download_interpro_pfam(
    base_dir: Path,
    overwrite: bool = False,
    include_huge: bool = False,
    include_pfam: bool = False,
) -> Dict[str, object]:
    out_dir = base_dir / "InterPro_Pfam"
    ensure_dir(out_dir)

    files = {
        "entry_list": INTERPRO_FILES["entry_list"],
        "interpro_xml": INTERPRO_FILES["interpro_xml"],
        "parent_child_tree": INTERPRO_FILES["parent_child_tree"],
    }
    if include_huge:
        files["protein2ipr"] = INTERPRO_FILES["protein2ipr"]
    if include_pfam:
        files.update(PFAM_FILES)

    metas = []
    errors = []
    for name, url in files.items():
        out_file = out_dir / url.split("/")[-1]
        try:
            meta = download_file(url, out_file, overwrite=overwrite, min_bytes=1000)
            meta.update({"source": "InterPro/Pfam", "type": name})
            metas.append(meta)
        except Exception as e:
            errors.append({"type": name, "url": url, "error": repr(e)})
            print(f"[InterPro/Pfam] {name} failed: {e}")

    if not include_huge:
        save_text(
            out_dir / "OPTIONAL_HUGE_FILES_NOTE.txt",
            """
            The full InterPro protein2ipr.dat.gz mapping is very large.
            It was skipped by default.

            To download it:
              python Step0B_Download_NoLeakage_Feature_Databases.py --sources interpro --include-interpro-huge

            Often you do not need the huge file because UniProt TSV already includes
            InterPro/Pfam cross-references for human proteins.
            """
        )

    summary = {"source": "InterPro/Pfam", "downloads": metas, "errors": errors}
    save_json(out_dir / "download_metadata.json", summary)
    return summary


# =============================================================================
# MANIFEST
# =============================================================================

def write_manifest(base_dir: Path) -> None:
    rows = [
        {
            "source": "HGNC",
            "purpose": "gene identifiers and official symbols",
            "safe_for_non_circular_features": "yes",
            "notes": "Identifier source; not a druggability label source.",
        },
        {
            "source": "UniProt",
            "purpose": "protein sequence, length, subcellular location, domain and database cross-references",
            "safe_for_non_circular_features": "yes",
            "notes": "Avoid fields that explicitly encode drug-target status if later added.",
        },
        {
            "source": "Ensembl",
            "purpose": "gene models, coordinates, transcripts, protein FASTA",
            "safe_for_non_circular_features": "yes",
            "notes": "Genomic annotation source.",
        },
        {
            "source": "STRING",
            "purpose": "protein network connectivity and interaction scores",
            "safe_for_non_circular_features": "yes_with_caution",
            "notes": "Use degree/connectivity features; do not use drug-target annotations.",
        },
        {
            "source": "GTEx",
            "purpose": "tissue expression features",
            "safe_for_non_circular_features": "yes",
            "notes": "Expression source, not direct druggability labels.",
        },
        {
            "source": "AlphaFold",
            "purpose": "predicted structures and pLDDT-derived structural confidence/disorder",
            "safe_for_non_circular_features": "yes",
            "notes": "Large optional source.",
        },
        {
            "source": "InterPro/Pfam",
            "purpose": "protein families and domains",
            "safe_for_non_circular_features": "yes",
            "notes": "Family/domain features may be strong predictors but are not direct drug evidence.",
        },
        {
            "source": "ChEMBL/OpenTargets/Pharos/DGIdb/DrugBank",
            "purpose": "LABELS/EVIDENCE ONLY",
            "safe_for_non_circular_features": "no",
            "notes": "Do not download or use here for model features; use separately for target labels only.",
        },
    ]
    pd.DataFrame(rows).to_csv(base_dir / "MANIFEST_feature_sources.csv", index=False)


# =============================================================================
# CLI
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Download no-leakage biological/protein feature databases."
    )

    p.add_argument("--outdir", default="feature_databases", help="Output directory.")
    p.add_argument(
        "--sources",
        nargs="+",
        default=[],
        choices=["hgnc", "uniprot", "ensembl", "string", "gtex", "alphafold", "interpro"],
        help="Sources to download.",
    )
    p.add_argument("--minimal", action="store_true", help="Download only HGNC + UniProt.")
    p.add_argument("--recommended", action="store_true", help="Download HGNC + UniProt + Ensembl + STRING + GTEx median candidates.")
    p.add_argument("--full", action="store_true", help="Download recommended + AlphaFold + InterPro/Pfam optional files.")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing files.")
    p.add_argument("--uniprot-reviewed-only", action="store_true", help="Use reviewed UniProt entries only.")
    p.add_argument("--include-string-detailed", action="store_true", help="Download detailed STRING links too.")
    p.add_argument("--include-gtex-full-tpm", action="store_true", help="Download full GTEx gene TPM matrix, if available. Large.")
    p.add_argument("--include-alphafold", action="store_true", help="Download AlphaFold human proteome archive. Very large.")
    p.add_argument("--extract-alphafold", action="store_true", help="Extract AlphaFold archive after download. Very large.")
    p.add_argument("--include-interpro-huge", action="store_true", help="Download huge InterPro protein2ipr.dat.gz.")
    p.add_argument("--include-pfam", action="store_true", help="Download Pfam bulk files. Can be large.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    base_dir = Path(args.outdir)
    ensure_dir(base_dir)

    sources = set(args.sources)

    if args.minimal:
        sources.update(["hgnc", "uniprot"])

    if args.recommended:
        sources.update(["hgnc", "uniprot", "ensembl", "string", "gtex"])

    if args.full:
        sources.update(["hgnc", "uniprot", "ensembl", "string", "gtex", "alphafold", "interpro"])
        args.include_alphafold = True
        args.include_interpro_huge = True
        args.include_pfam = True

    if args.include_alphafold:
        sources.add("alphafold")

    if not sources:
        print("No sources selected. Use --minimal, --recommended, --full, or --sources ...")
        sys.exit(1)

    write_leakage_warning(base_dir)
    write_manifest(base_dir)

    run_meta = {
        "started": now_string(),
        "output_directory": str(base_dir.resolve()),
        "sources": sorted(sources),
        "leakage_databases_excluded": LEAKAGE_DATABASES,
        "results": {},
    }

    print("=" * 100)
    print("NO-LEAKAGE FEATURE DATABASE DOWNLOADER")
    print("=" * 100)
    print(f"Output directory: {base_dir.resolve()}")
    print(f"Sources: {', '.join(sorted(sources))}")
    print("\nExcluded leakage sources:")
    for x in LEAKAGE_DATABASES:
        print(f"  - {x}")
    print("=" * 100)

    if "hgnc" in sources:
        run_meta["results"]["hgnc"] = download_hgnc(base_dir, overwrite=args.overwrite)

    if "uniprot" in sources:
        run_meta["results"]["uniprot"] = download_uniprot(
            base_dir,
            overwrite=args.overwrite,
            reviewed_only=args.uniprot_reviewed_only,
        )

    if "ensembl" in sources:
        run_meta["results"]["ensembl"] = download_ensembl(base_dir, overwrite=args.overwrite)

    if "string" in sources:
        run_meta["results"]["string"] = download_string(
            base_dir,
            overwrite=args.overwrite,
            include_detailed=args.include_string_detailed,
        )

    if "gtex" in sources:
        run_meta["results"]["gtex"] = download_gtex(
            base_dir,
            overwrite=args.overwrite,
            include_full_tpm=args.include_gtex_full_tpm,
        )

    if "alphafold" in sources:
        run_meta["results"]["alphafold"] = download_alphafold(
            base_dir,
            overwrite=args.overwrite,
            extract=args.extract_alphafold,
        )

    if "interpro" in sources:
        run_meta["results"]["interpro"] = download_interpro_pfam(
            base_dir,
            overwrite=args.overwrite,
            include_huge=args.include_interpro_huge,
            include_pfam=args.include_pfam,
        )

    run_meta["finished"] = now_string()
    save_json(base_dir / "download_run_metadata.json", run_meta)

    print("\n" + "=" * 100)
    print("DOWNLOAD COMPLETE")
    print("=" * 100)
    print(f"Run metadata: {base_dir / 'download_run_metadata.json'}")
    print(f"Feature-source manifest: {base_dir / 'MANIFEST_feature_sources.csv'}")
    print(f"Leakage warning: {base_dir / 'DO_NOT_USE_AS_LABEL_FEATURES.txt'}")
    print("\nNext step:")
    print("  Build a local Step1 master feature table by merging:")
    print("    Step0_Output/03_HumanGene_DruggabilityLabels.csv")
    print("    + feature_databases/HGNC")
    print("    + feature_databases/UniProt")
    print("    + feature_databases/Ensembl")
    print("    + feature_databases/STRING")
    print("    + feature_databases/GTEx")
    print("=" * 100)


if __name__ == "__main__":
    main()
 
