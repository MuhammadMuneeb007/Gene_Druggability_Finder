#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature5_InterProPfam.py

Feature 5: InterPro / Pfam protein-domain and family architecture features.

Purpose
-------
Build HGNC-level protein domain/family features from local InterPro and Pfam files.

Input files expected
--------------------
feature_databases/InterPro_Pfam/
    protein2ipr.dat.gz
    Pfam-A.regions.tsv.gz
    entry.list
    ParentChildTreeFile.txt
    interpro.xml.gz              optional, not required
    Pfam-A.full.gz               optional, not parsed by default because huge
    Pfam-A.hmm.gz                optional, not parsed by default

Default HGNC input:
    databases/HGNC/hgnc_complete_set.txt

Main outputs
------------
feature5_interpro_pfam/
    processed/feature5_interpro_long.csv
    processed/feature5_pfam_regions_long.csv
    processed/feature5_interpro_pfam_gene_features.csv
    processed/feature5_interpro_pfam_hgnc_merged.csv
    feature5_interpro_pfam_summary.txt
    feature5_interpro_pfam_run_metadata.json

Leakage policy
--------------
Included:
    InterPro IDs, InterPro names/types, Pfam accessions, Pfam names,
    domain counts, family counts, repeat counts, region counts,
    domain architecture complexity, Pfam coverage proxy,
    biological keyword counts from domain names.

Excluded:
    known drug-target labels,
    ChEMBL / DrugBank / DGIdb labels,
    approved-drug annotations,
    target development level,
    any "known druggable family" manual labels.

Run
---
    python Feature5_InterProPfam.py

Fast test:
    python Feature5_InterProPfam.py --max-lines 100000

If files are somewhere else:
    python Feature5_InterProPfam.py --dbdir feature_databases/InterPro_Pfam
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pandas as pd


# =============================================================================
# DEFAULTS
# =============================================================================

DEFAULT_HGNC_PATH = Path("databases") / "HGNC" / "hgnc_complete_set.txt"
DEFAULT_DBDIR = Path("feature_databases") / "InterPro_Pfam"
DEFAULT_OUTDIR = Path("feature5_interpro_pfam")

KEYWORD_GROUPS = {
    "kinase": ["kinase", "protein kinase", "tyrosine kinase", "serine/threonine kinase"],
    "phosphatase": ["phosphatase"],
    "protease": ["protease", "peptidase", "proteinase"],
    "receptor": ["receptor"],
    "gpcr": ["g protein-coupled receptor", "g-protein coupled receptor", "gpcr", "7tm", "seven transmembrane"],
    "ion_channel": ["ion channel", "channel", "voltage-gated", "ligand-gated"],
    "transporter": ["transporter", "transport", "solute carrier", "abc transporter"],
    "enzyme": ["enzyme", "oxidase", "reductase", "transferase", "hydrolase", "ligase", "isomerase", "lyase"],
    "dna_binding": ["dna-binding", "dna binding", "helix-turn-helix", "homeobox", "zinc finger"],
    "rna_binding": ["rna-binding", "rna binding", "rRM", "ribonucleoprotein"],
    "zinc_finger": ["zinc finger", "zf-", "zinc-finger"],
    "immunoglobulin": ["immunoglobulin", "ig-like", "ig domain"],
    "transmembrane": ["transmembrane", "membrane", "7tm", "tm domain"],
    "secreted_ecm": ["extracellular", "secreted", "collagen", "matrix", "ecm"],
    "repeat": ["repeat", "leucine-rich repeat", "ankyrin repeat", "tetratricopeptide"],
    "binding_site": ["binding site", "active site", "metal-binding", "calcium-binding", "nucleotide-binding"],
    "coiled_coil": ["coiled coil", "coiled-coil"],
    "disorder_related": ["low complexity", "intrinsically disordered", "disordered"],
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


def normalize_gene_symbol(x: Any) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


def clean_id(x: Any) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip()
    if not s or s.lower() in {"nan", "none"}:
        return ""
    return s


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


def uniprot_base(acc: str) -> str:
    """
    Normalise UniProt/Pfam accessions.

    Examples:
        P04637       -> P04637
        P04637-2     -> P04637
        P04637.1     -> P04637
    """
    acc = clean_id(acc)
    if not acc:
        return ""
    acc = acc.split(".")[0]
    acc = acc.split("-")[0]
    return acc


def safe_float(x: Any) -> float:
    try:
        if x is None or pd.isna(x):
            return np.nan
        return float(x)
    except Exception:
        return np.nan


def safe_int(x: Any) -> int:
    try:
        if x is None or pd.isna(x):
            return 0
        return int(float(x))
    except Exception:
        return 0


def open_text_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", errors="ignore")
    return open(path, "rt", errors="ignore")


def find_file(dbdir: Path, candidates: List[str], required: bool = True) -> Optional[Path]:
    for name in candidates:
        p = dbdir / name
        if p.exists() and p.stat().st_size > 0:
            return p

    if required:
        raise FileNotFoundError(
            f"Could not find any of these files in {dbdir}:\n"
            + "\n".join([f"  - {x}" for x in candidates])
        )
    return None


def keyword_counts(texts: Iterable[str]) -> Dict[str, int]:
    joined = " ; ".join([str(x).lower() for x in texts if pd.notna(x)])
    out = {}
    for group, keys in KEYWORD_GROUPS.items():
        out[f"feature5_domain_keyword_count_{group}"] = int(sum(joined.count(k.lower()) for k in keys))
        out[f"feature5_domain_keyword_has_{group}"] = int(any(k.lower() in joined for k in keys))
    return out


def interval_union_length(intervals: List[Tuple[int, int]]) -> int:
    clean = []
    for a, b in intervals:
        a = safe_int(a)
        b = safe_int(b)
        if a <= 0 or b <= 0:
            continue
        if b < a:
            a, b = b, a
        clean.append((a, b))

    if not clean:
        return 0

    clean.sort()
    total = 0
    cur_a, cur_b = clean[0]

    for a, b in clean[1:]:
        if a <= cur_b + 1:
            cur_b = max(cur_b, b)
        else:
            total += cur_b - cur_a + 1
            cur_a, cur_b = a, b

    total += cur_b - cur_a + 1
    return int(total)


# =============================================================================
# HGNC
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
        else:
            log("[WARNING] Could not find locus_group/locus_type. Keeping all HGNC genes.")

    hgnc["gene_symbol"] = hgnc["symbol"].map(normalize_gene_symbol)

    if "uniprot_ids" not in hgnc.columns:
        hgnc["uniprot_ids"] = ""

    if "ensembl_gene_id" not in hgnc.columns:
        hgnc["ensembl_gene_id"] = ""

    if "entrez_id" not in hgnc.columns:
        hgnc["entrez_id"] = ""

    hgnc = hgnc.drop_duplicates(subset=["gene_symbol"], keep="first").copy()

    log(f"[HGNC ROWS] {before} -> {len(hgnc)}")
    log(f"[PROTEIN CODING ONLY] {protein_coding_only}")
    log(f"[GENES WITH UNIPROT] {(hgnc['uniprot_ids'].fillna('').astype(str).str.len() > 0).sum()}")

    return hgnc


def build_uniprot_gene_map(hgnc: pd.DataFrame) -> Tuple[Dict[str, Set[str]], Dict[str, List[str]]]:
    """
    Returns:
        uniprot_to_genes: base UniProt accession -> gene symbols
        gene_to_uniprots: gene symbol -> list of base UniProt accessions
    """
    uniprot_to_genes: Dict[str, Set[str]] = defaultdict(set)
    gene_to_uniprots: Dict[str, List[str]] = defaultdict(list)

    for _, row in hgnc.iterrows():
        gene = row["gene_symbol"]
        accs = split_uniprot_ids(row.get("uniprot_ids", ""))

        for acc in accs:
            base = uniprot_base(acc)
            if not base:
                continue
            uniprot_to_genes[base].add(gene)
            if base not in gene_to_uniprots[gene]:
                gene_to_uniprots[gene].append(base)

    log("=" * 100)
    log("[UNIPROT MAP]")
    log(f"[UNIQUE UNIPROT IDS] {len(uniprot_to_genes)}")
    log(f"[GENES WITH UNIPROT IDS] {len(gene_to_uniprots)}")

    return uniprot_to_genes, gene_to_uniprots


# =============================================================================
# INTERPRO METADATA
# =============================================================================

def parse_entry_list(entry_list_path: Optional[Path]) -> Dict[str, Dict[str, str]]:
    """
    Parse InterPro entry.list if available.

    The file format can vary between releases, so this parser is permissive.
    It searches each tab-delimited line for an IPR accession and stores the
    remaining columns as type/name where possible.
    """
    meta: Dict[str, Dict[str, str]] = {}

    if not entry_list_path or not entry_list_path.exists():
        log("[WARNING] entry.list not found. InterPro type/name metadata will be limited.")
        return meta

    log("=" * 100)
    log(f"[READ ENTRY LIST] {entry_list_path}")

    with open_text_maybe_gzip(entry_list_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            ipr = ""
            ipr_idx = -1

            for i, part in enumerate(parts):
                if re.match(r"^IPR\d+$", part.strip()):
                    ipr = part.strip()
                    ipr_idx = i
                    break

            if not ipr:
                continue

            other = [p.strip() for j, p in enumerate(parts) if j != ipr_idx and p.strip()]

            entry_type = ""
            entry_name = ""

            # Common layouts are usually either:
            # IPRxxxx type name
            # type IPRxxxx name
            # IPRxxxx name type
            if len(other) >= 2:
                if other[0].lower() in {"domain", "family", "repeat", "site", "homologous_superfamily", "active_site", "binding_site"}:
                    entry_type = other[0]
                    entry_name = other[1]
                elif other[-1].lower() in {"domain", "family", "repeat", "site", "homologous_superfamily", "active_site", "binding_site"}:
                    entry_type = other[-1]
                    entry_name = other[0]
                else:
                    entry_type = other[0]
                    entry_name = other[1]
            elif len(other) == 1:
                entry_name = other[0]

            meta[ipr] = {
                "interpro_id": ipr,
                "interpro_type": entry_type,
                "interpro_name_from_entry_list": entry_name,
            }

    log(f"[INTERPRO METADATA ENTRIES] {len(meta)}")
    return meta


def parse_parent_child_tree(parent_child_path: Optional[Path]) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, int]]:
    """
    Parse ParentChildTreeFile.txt approximately.

    Returns:
        depth_map: IPR -> hierarchy depth proxy
        child_count: IPR -> number of direct children proxy
        parent_count: IPR -> number of parents proxy
    """
    depth_map: Dict[str, int] = {}
    child_count: Dict[str, int] = defaultdict(int)
    parent_count: Dict[str, int] = defaultdict(int)

    if not parent_child_path or not parent_child_path.exists():
        log("[WARNING] ParentChildTreeFile.txt not found. Hierarchy features will be zero.")
        return depth_map, child_count, parent_count

    log("=" * 100)
    log(f"[READ PARENT-CHILD TREE] {parent_child_path}")

    stack: List[str] = []

    with open_text_maybe_gzip(parent_child_path) as f:
        for line in f:
            if not line.strip():
                continue

            m = re.search(r"(IPR\d+)", line)
            if not m:
                continue

            ipr = m.group(1)

            # Depth proxy from leading non-word/tree characters.
            leading = len(line) - len(line.lstrip(" -+|\t."))
            depth = max(0, leading // 2)

            depth_map[ipr] = min(depth_map.get(ipr, depth), depth)

            while len(stack) > depth:
                stack.pop()

            if stack:
                parent = stack[-1]
                child_count[parent] += 1
                parent_count[ipr] += 1

            if len(stack) == depth:
                stack.append(ipr)
            elif len(stack) < depth:
                stack.append(ipr)

    log(f"[HIERARCHY IDS] {len(depth_map)}")
    return depth_map, child_count, parent_count


# =============================================================================
# PROTEIN2IPR PARSING
# =============================================================================

def parse_protein2ipr(
    protein2ipr_path: Path,
    uniprot_to_genes: Dict[str, Set[str]],
    interpro_meta: Dict[str, Dict[str, str]],
    depth_map: Dict[str, int],
    child_count: Dict[str, int],
    parent_count: Dict[str, int],
    processed_dir: Path,
    max_lines: int = 0,
) -> pd.DataFrame:
    """
    Stream protein2ipr.dat.gz and keep only rows mapping to our HGNC UniProt IDs.
    """
    out_path = processed_dir / "feature5_interpro_long.csv"

    log("=" * 100)
    log(f"[STREAM PROTEIN2IPR] {protein2ipr_path}")
    log(f"[OUTPUT] {out_path}")

    rows = []
    matched = 0
    total = 0
    t0 = time.time()

    with open_text_maybe_gzip(protein2ipr_path) as f:
        for line in f:
            total += 1

            if max_lines and total > max_lines:
                break

            if total % 5_000_000 == 0:
                elapsed = time.time() - t0
                log(f"[PROTEIN2IPR] scanned={total:,} matched={matched:,} elapsed={elapsed/60:.1f} min")

            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 2:
                continue

            protein_acc_raw = clean_id(parts[0])
            protein_acc = uniprot_base(protein_acc_raw)

            if protein_acc not in uniprot_to_genes:
                continue

            ipr = ""
            ipr_name_from_file = ""

            for p in parts[1:]:
                p2 = p.strip()
                if re.match(r"^IPR\d+$", p2):
                    ipr = p2
                    break

            if not ipr:
                continue

            # Name is often the field after InterPro ID or one of the text fields.
            for i, p in enumerate(parts):
                if p.strip() == ipr and i + 1 < len(parts):
                    ipr_name_from_file = parts[i + 1].strip()
                    break

            meta = interpro_meta.get(ipr, {})
            ipr_type = meta.get("interpro_type", "")
            ipr_name = meta.get("interpro_name_from_entry_list", "") or ipr_name_from_file

            for gene in uniprot_to_genes[protein_acc]:
                rows.append({
                    "gene_symbol": gene,
                    "uniprot_accession": protein_acc,
                    "uniprot_accession_raw": protein_acc_raw,
                    "interpro_id": ipr,
                    "interpro_name": ipr_name,
                    "interpro_type": ipr_type,
                    "interpro_hierarchy_depth": depth_map.get(ipr, np.nan),
                    "interpro_child_count": child_count.get(ipr, 0),
                    "interpro_parent_count": parent_count.get(ipr, 0),
                })
                matched += 1

    df = pd.DataFrame(rows)

    if not df.empty:
        df = df.drop_duplicates()

    df.to_csv(out_path, index=False)

    log(f"[PROTEIN2IPR DONE] scanned={total:,} matched_rows={len(df):,}")
    return df


# =============================================================================
# PFAM REGIONS PARSING
# =============================================================================

def detect_pfam_header(first_line: str) -> Tuple[bool, Dict[str, int]]:
    """
    Detect header for Pfam-A.regions.tsv.gz.

    Common columns include:
        pfamseq_acc, seq_version_crc64, seq_start, seq_end,
        hmm_acc, hmm_name, type, hmm_start, hmm_end, hmm_length,
        bit_score, evalue, significance, clan
    """
    parts = first_line.rstrip("\n").split("\t")
    lower = [p.strip().lower() for p in parts]

    if "pfamseq_acc" in lower or "hmm_acc" in lower or "seq_start" in lower:
        return True, {name: i for i, name in enumerate(lower)}

    # Fallback old/unknown layout.
    fallback = {
        "pfamseq_acc": 0,
        "seq_start": 2,
        "seq_end": 3,
        "hmm_acc": 4,
        "hmm_name": 5,
        "type": 6,
        "hmm_start": 7,
        "hmm_end": 8,
        "hmm_length": 9,
        "bit_score": 10,
        "evalue": 11,
        "significance": 12,
        "clan": 13,
    }
    return False, fallback


def get_part(parts: List[str], header: Dict[str, int], key: str, default: str = "") -> str:
    idx = header.get(key)
    if idx is None or idx < 0 or idx >= len(parts):
        return default
    return parts[idx].strip()


def parse_pfam_regions(
    pfam_regions_path: Path,
    uniprot_to_genes: Dict[str, Set[str]],
    processed_dir: Path,
    max_lines: int = 0,
) -> pd.DataFrame:
    """
    Stream Pfam-A.regions.tsv.gz and keep only regions mapping to our HGNC UniProt IDs.
    """
    out_path = processed_dir / "feature5_pfam_regions_long.csv"

    log("=" * 100)
    log(f"[STREAM PFAM REGIONS] {pfam_regions_path}")
    log(f"[OUTPUT] {out_path}")

    rows = []
    total = 0
    matched = 0
    t0 = time.time()

    with open_text_maybe_gzip(pfam_regions_path) as f:
        first = f.readline()
        if not first:
            df = pd.DataFrame()
            df.to_csv(out_path, index=False)
            return df

        has_header, header = detect_pfam_header(first)

        if not has_header:
            # Process first line as data.
            data_iter = [first]
        else:
            data_iter = []

        for line in data_iter:
            pass

        # Re-open easier, skipping header if needed.
    with open_text_maybe_gzip(pfam_regions_path) as f:
        if has_header:
            next(f, None)

        for line in f:
            total += 1

            if max_lines and total > max_lines:
                break

            if total % 5_000_000 == 0:
                elapsed = time.time() - t0
                log(f"[PFAM] scanned={total:,} matched={matched:,} elapsed={elapsed/60:.1f} min")

            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 5:
                continue

            acc_raw = get_part(parts, header, "pfamseq_acc", parts[0] if parts else "")
            acc = uniprot_base(acc_raw)

            if acc not in uniprot_to_genes:
                continue

            seq_start = safe_int(get_part(parts, header, "seq_start"))
            seq_end = safe_int(get_part(parts, header, "seq_end"))

            pfam_acc = get_part(parts, header, "hmm_acc")
            pfam_name = get_part(parts, header, "hmm_name")
            pfam_type = get_part(parts, header, "type")
            hmm_start = safe_int(get_part(parts, header, "hmm_start"))
            hmm_end = safe_int(get_part(parts, header, "hmm_end"))
            hmm_length = safe_int(get_part(parts, header, "hmm_length"))
            bit_score = safe_float(get_part(parts, header, "bit_score"))
            evalue = safe_float(get_part(parts, header, "evalue"))
            clan = get_part(parts, header, "clan")

            region_length = abs(seq_end - seq_start) + 1 if seq_start > 0 and seq_end > 0 else np.nan

            for gene in uniprot_to_genes[acc]:
                rows.append({
                    "gene_symbol": gene,
                    "uniprot_accession": acc,
                    "uniprot_accession_raw": acc_raw,
                    "pfam_acc": pfam_acc,
                    "pfam_name": pfam_name,
                    "pfam_type": pfam_type,
                    "seq_start": seq_start,
                    "seq_end": seq_end,
                    "region_length": region_length,
                    "hmm_start": hmm_start,
                    "hmm_end": hmm_end,
                    "hmm_length": hmm_length,
                    "bit_score": bit_score,
                    "evalue": evalue,
                    "clan": clan,
                })
                matched += 1

    df = pd.DataFrame(rows)

    if not df.empty:
        df = df.drop_duplicates()

    df.to_csv(out_path, index=False)

    log(f"[PFAM DONE] scanned={total:,} matched_rows={len(df):,}")
    return df


# =============================================================================
# FEATURE AGGREGATION
# =============================================================================

def count_by_category(df: pd.DataFrame, group_col: str, value_col: str, prefix: str) -> pd.DataFrame:
    if df.empty or value_col not in df.columns:
        return pd.DataFrame(columns=[group_col])

    tmp = df.copy()
    tmp[value_col] = tmp[value_col].fillna("").astype(str).str.strip()
    tmp = tmp[tmp[value_col] != ""].copy()

    if tmp.empty:
        return pd.DataFrame(columns=[group_col])

    wide = (
        tmp.groupby([group_col, value_col])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )

    rename = {}
    for c in wide.columns:
        if c == group_col:
            continue
        safe = re.sub(r"[^A-Za-z0-9]+", "_", str(c).lower()).strip("_")
        rename[c] = f"{prefix}_{safe}"

    wide = wide.rename(columns=rename)
    return wide


def build_gene_features(
    hgnc: pd.DataFrame,
    interpro_df: pd.DataFrame,
    pfam_df: pd.DataFrame,
    gene_to_uniprots: Dict[str, List[str]],
    processed_dir: Path,
) -> pd.DataFrame:
    genes = pd.DataFrame({"gene_symbol": sorted(hgnc["gene_symbol"].unique().tolist())})
    records = []

    interpro_by_gene = dict(tuple(interpro_df.groupby("gene_symbol"))) if not interpro_df.empty else {}
    pfam_by_gene = dict(tuple(pfam_df.groupby("gene_symbol"))) if not pfam_df.empty else {}

    for gene in genes["gene_symbol"]:
        ipr = interpro_by_gene.get(gene, pd.DataFrame())
        pf = pfam_by_gene.get(gene, pd.DataFrame())

        uniprots = gene_to_uniprots.get(gene, [])

        ipr_ids = sorted(set(ipr["interpro_id"].dropna().astype(str))) if not ipr.empty else []
        ipr_names = sorted(set(ipr["interpro_name"].dropna().astype(str))) if not ipr.empty and "interpro_name" in ipr.columns else []
        ipr_types = sorted(set(ipr["interpro_type"].dropna().astype(str))) if not ipr.empty and "interpro_type" in ipr.columns else []

        pfam_ids = sorted(set(pf["pfam_acc"].dropna().astype(str))) if not pf.empty else []
        pfam_names = sorted(set(pf["pfam_name"].dropna().astype(str))) if not pf.empty and "pfam_name" in pf.columns else []
        clans = sorted(set(pf["clan"].dropna().astype(str))) if not pf.empty and "clan" in pf.columns else []

        region_lengths = pd.to_numeric(pf["region_length"], errors="coerce").dropna().tolist() if not pf.empty else []
        bit_scores = pd.to_numeric(pf["bit_score"], errors="coerce").dropna().tolist() if not pf.empty else []
        evalues = pd.to_numeric(pf["evalue"], errors="coerce").dropna().tolist() if not pf.empty else []

        intervals = []
        max_end = 0
        if not pf.empty:
            for _, row in pf.iterrows():
                s = safe_int(row.get("seq_start"))
                e = safe_int(row.get("seq_end"))
                if s > 0 and e > 0:
                    intervals.append((s, e))
                    max_end = max(max_end, s, e)

        union_len = interval_union_length(intervals)
        coverage_proxy = union_len / max_end if max_end > 0 else np.nan

        hierarchy_depths = pd.to_numeric(ipr["interpro_hierarchy_depth"], errors="coerce").dropna().tolist() if not ipr.empty else []
        child_counts = pd.to_numeric(ipr["interpro_child_count"], errors="coerce").dropna().tolist() if not ipr.empty else []
        parent_counts = pd.to_numeric(ipr["interpro_parent_count"], errors="coerce").dropna().tolist() if not ipr.empty else []

        text_for_keywords = ipr_names + pfam_names + ipr_types
        kw = keyword_counts(text_for_keywords)

        rec = {
            "gene_symbol": gene,

            # Mapping
            "feature5_interpro_pfam_uniprot_count": len(uniprots),
            "feature5_interpro_pfam_uniprot_ids": ";".join(uniprots),

            # InterPro
            "feature5_interpro_has_annotation": int(len(ipr_ids) > 0),
            "feature5_interpro_row_count": int(len(ipr)) if not ipr.empty else 0,
            "feature5_interpro_unique_id_count": int(len(ipr_ids)),
            "feature5_interpro_unique_name_count": int(len(ipr_names)),
            "feature5_interpro_unique_type_count": int(len([x for x in ipr_types if x])),
            "feature5_interpro_ids": ";".join(ipr_ids[:100]),
            "feature5_interpro_names": ";".join(ipr_names[:100]),
            "feature5_interpro_types": ";".join([x for x in ipr_types if x][:50]),

            # InterPro hierarchy
            "feature5_interpro_mean_hierarchy_depth": float(np.mean(hierarchy_depths)) if hierarchy_depths else np.nan,
            "feature5_interpro_max_hierarchy_depth": float(np.max(hierarchy_depths)) if hierarchy_depths else np.nan,
            "feature5_interpro_mean_child_count": float(np.mean(child_counts)) if child_counts else np.nan,
            "feature5_interpro_max_child_count": float(np.max(child_counts)) if child_counts else np.nan,
            "feature5_interpro_mean_parent_count": float(np.mean(parent_counts)) if parent_counts else np.nan,
            "feature5_interpro_max_parent_count": float(np.max(parent_counts)) if parent_counts else np.nan,

            # Pfam
            "feature5_pfam_has_annotation": int(len(pfam_ids) > 0),
            "feature5_pfam_region_count": int(len(pf)) if not pf.empty else 0,
            "feature5_pfam_unique_domain_count": int(len(pfam_ids)),
            "feature5_pfam_unique_name_count": int(len(pfam_names)),
            "feature5_pfam_unique_clan_count": int(len([x for x in clans if x])),
            "feature5_pfam_repeat_region_count": int(max(0, len(pf) - len(pfam_ids))) if not pf.empty else 0,
            "feature5_pfam_domain_architecture_complexity": int(len(pfam_ids) + max(0, len(pf) - len(pfam_ids))) if not pf.empty else 0,
            "feature5_pfam_ids": ";".join(pfam_ids[:100]),
            "feature5_pfam_names": ";".join(pfam_names[:100]),
            "feature5_pfam_clans": ";".join([x for x in clans if x][:100]),

            # Pfam region lengths
            "feature5_pfam_mean_region_length": float(np.mean(region_lengths)) if region_lengths else np.nan,
            "feature5_pfam_median_region_length": float(np.median(region_lengths)) if region_lengths else np.nan,
            "feature5_pfam_min_region_length": float(np.min(region_lengths)) if region_lengths else np.nan,
            "feature5_pfam_max_region_length": float(np.max(region_lengths)) if region_lengths else np.nan,
            "feature5_pfam_total_union_region_length": int(union_len),
            "feature5_pfam_max_observed_residue_position": int(max_end),
            "feature5_pfam_coverage_proxy": float(coverage_proxy) if not pd.isna(coverage_proxy) else np.nan,

            # Pfam confidence/proxy scores
            "feature5_pfam_mean_bit_score": float(np.mean(bit_scores)) if bit_scores else np.nan,
            "feature5_pfam_max_bit_score": float(np.max(bit_scores)) if bit_scores else np.nan,
            "feature5_pfam_min_evalue": float(np.min(evalues)) if evalues else np.nan,
            "feature5_pfam_median_evalue": float(np.median(evalues)) if evalues else np.nan,

            # Overall
            "feature5_interpro_pfam_has_any_domain_feature": int(len(ipr_ids) > 0 or len(pfam_ids) > 0),
            "feature5_interpro_pfam_total_unique_annotation_count": int(len(set(ipr_ids + pfam_ids))),
        }

        rec.update(kw)
        records.append(rec)

    features = pd.DataFrame(records)

    # Add InterPro type count columns
    ipr_type_counts = count_by_category(
        interpro_df,
        group_col="gene_symbol",
        value_col="interpro_type",
        prefix="feature5_interpro_type_count",
    )
    if not ipr_type_counts.empty:
        features = features.merge(ipr_type_counts, on="gene_symbol", how="left")

    # Add Pfam type count columns if available
    pfam_type_counts = count_by_category(
        pfam_df,
        group_col="gene_symbol",
        value_col="pfam_type",
        prefix="feature5_pfam_type_count",
    )
    if not pfam_type_counts.empty:
        features = features.merge(pfam_type_counts, on="gene_symbol", how="left")

    # Fill count/flag columns
    for c in features.columns:
        if c.startswith("feature5_") and (
            "_count" in c
            or "_has_" in c
            or c.endswith("_complexity")
        ):
            if features[c].dtype != object:
                features[c] = features[c].fillna(0)

    text_cols = [
        "feature5_interpro_pfam_uniprot_ids",
        "feature5_interpro_ids",
        "feature5_interpro_names",
        "feature5_interpro_types",
        "feature5_pfam_ids",
        "feature5_pfam_names",
        "feature5_pfam_clans",
    ]
    for c in text_cols:
        if c in features.columns:
            features[c] = features[c].fillna("")

    out_path = processed_dir / "feature5_interpro_pfam_gene_features.csv"
    features.to_csv(out_path, index=False)

    log("=" * 100)
    log("[SAVED GENE FEATURES]")
    log(f"[PATH]  {out_path}")
    log(f"[SHAPE] {features.shape[0]} genes x {features.shape[1]} columns")

    return features


def merge_with_hgnc(hgnc: pd.DataFrame, features: pd.DataFrame, processed_dir: Path) -> pd.DataFrame:
    merged = hgnc.merge(features, on="gene_symbol", how="left")

    feature_cols = [c for c in features.columns if c != "gene_symbol"]
    merged["feature5_interpro_pfam_has_any_feature"] = merged[feature_cols].notna().any(axis=1).astype(int)
    merged["feature5_interpro_pfam_n_nonmissing_features"] = merged[feature_cols].notna().sum(axis=1)

    out_path = processed_dir / "feature5_interpro_pfam_hgnc_merged.csv"
    merged.to_csv(out_path, index=False)

    log("=" * 100)
    log("[SAVED HGNC MERGED TABLE]")
    log(f"[PATH]  {out_path}")
    log(f"[SHAPE] {merged.shape[0]} genes x {merged.shape[1]} columns")

    return merged


def write_summary(
    outdir: Path,
    hgnc: pd.DataFrame,
    interpro_df: pd.DataFrame,
    pfam_df: pd.DataFrame,
    features: pd.DataFrame,
    merged: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    path = outdir / "feature5_interpro_pfam_summary.txt"

    lines = []
    lines.append("Feature 5: InterPro/Pfam protein domain and family features")
    lines.append("=" * 80)
    lines.append(f"Created at: {now_iso()}")
    lines.append(f"HGNC rows: {hgnc.shape[0]}")
    lines.append(f"InterPro long rows: {interpro_df.shape[0]}")
    lines.append(f"Pfam long rows: {pfam_df.shape[0]}")
    lines.append(f"Gene feature rows: {features.shape[0]}")
    lines.append(f"Gene feature columns: {features.shape[1]}")
    lines.append(f"Merged rows: {merged.shape[0]}")
    lines.append(f"Merged columns: {merged.shape[1]}")
    lines.append("")
    lines.append("Input files:")
    lines.append(f"dbdir: {Path(args.dbdir).resolve()}")
    lines.append(f"hgnc: {Path(args.hgnc).resolve()}")
    lines.append("")
    lines.append("Coverage:")
    if "feature5_interpro_has_annotation" in features.columns:
        lines.append(f"InterPro annotation coverage: {features['feature5_interpro_has_annotation'].mean():.4f}")
    if "feature5_pfam_has_annotation" in features.columns:
        lines.append(f"Pfam annotation coverage: {features['feature5_pfam_has_annotation'].mean():.4f}")
    if "feature5_interpro_pfam_has_any_domain_feature" in features.columns:
        lines.append(f"Any InterPro/Pfam coverage: {features['feature5_interpro_pfam_has_any_domain_feature'].mean():.4f}")
    if "feature5_interpro_unique_id_count" in features.columns:
        lines.append(f"Mean InterPro IDs per gene: {features['feature5_interpro_unique_id_count'].mean():.2f}")
        lines.append(f"Median InterPro IDs per gene: {features['feature5_interpro_unique_id_count'].median():.2f}")
    if "feature5_pfam_unique_domain_count" in features.columns:
        lines.append(f"Mean Pfam domains per gene: {features['feature5_pfam_unique_domain_count'].mean():.2f}")
        lines.append(f"Median Pfam domains per gene: {features['feature5_pfam_unique_domain_count'].median():.2f}")
    lines.append("")
    lines.append("Leakage policy:")
    lines.append("Included: InterPro/Pfam domain, family, repeat, site, architecture and coverage features.")
    lines.append("Excluded: known drug-target labels, ChEMBL/DrugBank/DGIdb labels, approved drug annotations, target development level.")

    path.write_text("\n".join(lines) + "\n")
    log(f"[SUMMARY] {path}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build InterPro/Pfam domain architecture features for HGNC genes."
    )
    parser.add_argument("--hgnc", default=str(DEFAULT_HGNC_PATH), help="Path to HGNC complete set.")
    parser.add_argument("--dbdir", default=str(DEFAULT_DBDIR), help="Directory containing InterPro/Pfam files.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory.")
    parser.add_argument("--all-hgnc-genes", action="store_true", help="Use all HGNC rows, not only protein-coding.")
    parser.add_argument("--max-lines", type=int, default=0, help="Debug only: parse first N lines from each huge file.")
    parser.add_argument("--skip-pfam", action="store_true", help="Skip Pfam regions parsing.")
    parser.add_argument("--skip-interpro", action="store_true", help="Skip protein2ipr parsing.")
    args = parser.parse_args()

    hgnc_path = Path(args.hgnc)
    dbdir = Path(args.dbdir)
    outdir = mkdir(Path(args.outdir))
    processed_dir = mkdir(outdir / "processed")

    log("=" * 100)
    log("FEATURE 5: INTERPRO / PFAM DOMAIN ARCHITECTURE")
    log("=" * 100)
    log(f"[HGNC]        {hgnc_path}")
    log(f"[DBDIR]       {dbdir.resolve()}")
    log(f"[OUTDIR]      {outdir.resolve()}")
    log(f"[MAX LINES]   {args.max_lines if args.max_lines else 'none'}")
    log(f"[SKIP PFAM]   {args.skip_pfam}")
    log(f"[SKIP IPR]    {args.skip_interpro}")
    log("=" * 100)

    if not dbdir.exists():
        raise FileNotFoundError(f"Database directory not found: {dbdir}")

    protein2ipr_path = find_file(dbdir, ["protein2ipr.dat.gz", "protein2ipr.dat"], required=not args.skip_interpro)
    pfam_regions_path = find_file(dbdir, ["Pfam-A.regions.tsv.gz", "Pfam-A.regions.tsv"], required=not args.skip_pfam)
    entry_list_path = find_file(dbdir, ["entry.list", "entry.list.gz"], required=False)
    parent_child_path = find_file(dbdir, ["ParentChildTreeFile.txt", "ParentChildTreeFile.txt.gz"], required=False)

    hgnc = load_hgnc(hgnc_path, protein_coding_only=not args.all_hgnc_genes)
    uniprot_to_genes, gene_to_uniprots = build_uniprot_gene_map(hgnc)

    interpro_meta = parse_entry_list(entry_list_path)
    depth_map, child_count, parent_count = parse_parent_child_tree(parent_child_path)

    if args.skip_interpro:
        interpro_df = pd.DataFrame(columns=[
            "gene_symbol", "uniprot_accession", "interpro_id", "interpro_name", "interpro_type"
        ])
        interpro_df.to_csv(processed_dir / "feature5_interpro_long.csv", index=False)
    else:
        interpro_df = parse_protein2ipr(
            protein2ipr_path=protein2ipr_path,
            uniprot_to_genes=uniprot_to_genes,
            interpro_meta=interpro_meta,
            depth_map=depth_map,
            child_count=child_count,
            parent_count=parent_count,
            processed_dir=processed_dir,
            max_lines=args.max_lines,
        )

    if args.skip_pfam:
        pfam_df = pd.DataFrame(columns=[
            "gene_symbol", "uniprot_accession", "pfam_acc", "pfam_name", "seq_start", "seq_end"
        ])
        pfam_df.to_csv(processed_dir / "feature5_pfam_regions_long.csv", index=False)
    else:
        pfam_df = parse_pfam_regions(
            pfam_regions_path=pfam_regions_path,
            uniprot_to_genes=uniprot_to_genes,
            processed_dir=processed_dir,
            max_lines=args.max_lines,
        )

    features = build_gene_features(
        hgnc=hgnc,
        interpro_df=interpro_df,
        pfam_df=pfam_df,
        gene_to_uniprots=gene_to_uniprots,
        processed_dir=processed_dir,
    )

    merged = merge_with_hgnc(
        hgnc=hgnc,
        features=features,
        processed_dir=processed_dir,
    )

    metadata = {
        "created_at": now_iso(),
        "hgnc": str(hgnc_path),
        "dbdir": str(dbdir.resolve()),
        "outdir": str(outdir.resolve()),
        "protein_coding_only": not args.all_hgnc_genes,
        "max_lines": args.max_lines,
        "inputs": {
            "protein2ipr": str(protein2ipr_path) if protein2ipr_path else "",
            "pfam_regions": str(pfam_regions_path) if pfam_regions_path else "",
            "entry_list": str(entry_list_path) if entry_list_path else "",
            "parent_child_tree": str(parent_child_path) if parent_child_path else "",
        },
        "outputs": {
            "interpro_long": str(processed_dir / "feature5_interpro_long.csv"),
            "pfam_regions_long": str(processed_dir / "feature5_pfam_regions_long.csv"),
            "gene_features": str(processed_dir / "feature5_interpro_pfam_gene_features.csv"),
            "hgnc_merged": str(processed_dir / "feature5_interpro_pfam_hgnc_merged.csv"),
        },
        "leakage_policy": {
            "included": [
                "InterPro domain/family/repeat/site annotations",
                "Pfam domain regions",
                "domain architecture complexity",
                "domain keyword counts from annotation text",
                "Pfam coverage proxy",
            ],
            "excluded": [
                "known drug target labels",
                "approved drug annotations",
                "ChEMBL labels",
                "DrugBank labels",
                "DGIdb labels",
                "target development level",
                "manual druggable-family labels",
            ],
        },
    }

    with open(outdir / "feature5_interpro_pfam_run_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    write_summary(outdir, hgnc, interpro_df, pfam_df, features, merged, args)

    log("=" * 100)
    log("[DONE]")
    log(f"[INTERPRO LONG] {processed_dir / 'feature5_interpro_long.csv'}")
    log(f"[PFAM LONG]     {processed_dir / 'feature5_pfam_regions_long.csv'}")
    log(f"[GENE FEATURES] {processed_dir / 'feature5_interpro_pfam_gene_features.csv'}")
    log(f"[HGNC MERGED]   {processed_dir / 'feature5_interpro_pfam_hgnc_merged.csv'}")
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