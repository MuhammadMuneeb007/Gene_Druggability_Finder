#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature6_UniProt.py

Feature 6: UniProt protein annotation and sequence-level features.

Purpose
-------
Build HGNC-level UniProt protein annotation features using:
  1. local UniProt files if present, or
  2. UniProt REST batch download fallback.

This feature adds a protein-annotation layer after:
  Feature 1 = DepMap / cell-line functional genomics
  Feature 2 = STRING / PPI graph
  Feature 3 = Reactome / pathway membership
  Feature 4 = PDB + AlphaFold structure
  Feature 5 = InterPro/Pfam domains

What this feature captures
--------------------------
For each HGNC protein-coding gene:
  - UniProt reviewed/Swiss-Prot status
  - protein length and mass
  - sequence composition
  - subcellular location flags
  - membrane/transmembrane/signal peptide flags
  - GO annotation counts
  - UniProt keyword counts
  - active-site / binding-site / motif / region annotation counts
  - functional text keyword flags

Leakage policy
--------------
Included:
  - UniProt protein biology annotation
  - subcellular localization
  - sequence length/composition
  - GO/keyword counts
  - topology and feature-count summaries

Excluded:
  - known drug-target labels
  - approved drug annotations
  - ChEMBL / DrugBank / DGIdb labels
  - target development level
  - manually curated druggable-family labels

Run
---
    python Feature6_UniProt.py

Fast test:
    python Feature6_UniProt.py --limit-genes 500

Force UniProt REST re-download:
    python Feature6_UniProt.py --force-download

Use local only and do not query UniProt:
    python Feature6_UniProt.py --no-api

Outputs
-------
feature6_uniprot/
    downloads/feature6_uniprot_rest_download.tsv
    processed/feature6_uniprot_protein_table.csv
    processed/feature6_uniprot_gene_features.csv
    processed/feature6_uniprot_hgnc_merged.csv
    feature6_uniprot_summary.txt
    feature6_uniprot_run_metadata.json
"""

from __future__ import annotations

import argparse
import gzip
import io
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
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


# =============================================================================
# DEFAULTS
# =============================================================================

DEFAULT_HGNC_PATH = Path("databases") / "HGNC" / "hgnc_complete_set.txt"
DEFAULT_DBDIR = Path("feature_databases") / "UniProt"
DEFAULT_OUTDIR = Path("feature6_uniprot")

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

HEADERS = {
    "User-Agent": "Feature6_UniProt/1.0 academic gene-druggability feature builder"
}

REQUEST_TIMEOUT = 120
MAX_RETRIES = 4
BATCH_SIZE = 250
SLEEP_BETWEEN_BATCHES = 0.4

# UniProt REST fields. If the full list fails because of a release-specific field
# name, the script automatically falls back to BASIC_UNIPROT_FIELDS.
UNIPROT_FIELDS = [
    "accession",
    "id",
    "reviewed",
    "protein_name",
    "gene_names",
    "organism_id",
    "length",
    "mass",
    "cc_subcellular_location",
    "ft_transmem",
    "ft_intramem",
    "ft_topo_dom",
    "ft_signal",
    "ft_domain",
    "ft_region",
    "ft_motif",
    "ft_binding",
    "ft_act_site",
    "keyword",
    "go",
    "go_p",
    "go_c",
    "go_f",
    "cc_function",
    "cc_pathway",
    "cc_tissue_specificity",
    "cc_alternative_products",
    "sequence",
]

BASIC_UNIPROT_FIELDS = [
    "accession",
    "id",
    "reviewed",
    "protein_name",
    "gene_names",
    "organism_id",
    "length",
    "mass",
    "cc_subcellular_location",
    "keyword",
    "go",
    "go_p",
    "go_c",
    "go_f",
    "sequence",
]

AA_GROUPS = {
    "hydrophobic": set("AILMFWYV"),
    "polar": set("STNQCY"),
    "positive": set("KRH"),
    "negative": set("DE"),
    "charged": set("KRHDE"),
    "aromatic": set("FWY"),
    "small": set("AGSTC"),
    "proline": set("P"),
    "cysteine": set("C"),
    "glycine": set("G"),
}

SUBCELLULAR_KEYWORDS = {
    "membrane": ["membrane", "cell membrane", "plasma membrane"],
    "plasma_membrane": ["plasma membrane", "cell membrane"],
    "secreted": ["secreted", "extracellular", "extracellular space"],
    "extracellular": ["extracellular"],
    "nucleus": ["nucleus", "nuclear"],
    "cytoplasm": ["cytoplasm", "cytosol", "cytoplasmic"],
    "mitochondrion": ["mitochondrion", "mitochondrial"],
    "endoplasmic_reticulum": ["endoplasmic reticulum", "er membrane"],
    "golgi": ["golgi"],
    "lysosome": ["lysosome", "lysosomal"],
    "peroxisome": ["peroxisome", "peroxisomal"],
    "cell_surface": ["cell surface", "surface"],
}

FUNCTION_KEYWORDS = {
    "enzyme": ["enzyme", "catalytic", "oxidoreductase", "transferase", "hydrolase", "ligase", "isomerase", "lyase"],
    "kinase": ["kinase", "phosphorylation"],
    "phosphatase": ["phosphatase", "dephosphorylation"],
    "protease": ["protease", "peptidase", "proteinase"],
    "receptor": ["receptor"],
    "transporter": ["transporter", "transport", "solute carrier", "abc transporter"],
    "ion_channel": ["ion channel", "channel", "voltage-gated", "ligand-gated"],
    "transcription": ["transcription", "transcription factor", "dna-binding"],
    "dna_binding": ["dna-binding", "dna binding", "chromatin", "histone"],
    "rna_binding": ["rna-binding", "rna binding", "ribonucleoprotein"],
    "immune": ["immune", "immunity", "cytokine", "interleukin", "interferon", "complement"],
    "signaling": ["signal", "signaling", "signalling", "signal transduction"],
    "metabolism": ["metabolism", "metabolic", "biosynthesis", "catabolism"],
    "cell_cycle": ["cell cycle", "mitosis", "meiotic", "checkpoint"],
    "apoptosis": ["apoptosis", "cell death"],
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


def clean_text(x: Any) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip()
    if s.lower() in {"nan", "none"}:
        return ""
    return s


def clean_accession(x: Any) -> str:
    s = clean_text(x)
    if not s:
        return ""
    return s.split("-")[0].split(".")[0]


def split_values(x: Any) -> List[str]:
    s = clean_text(x)
    if not s:
        return []
    vals = []
    for part in re.split(r"[|;,]+", s):
        part = part.strip()
        if part and part.lower() not in {"nan", "none"} and part not in vals:
            vals.append(part)
    return vals


def split_uniprot_ids(x: Any) -> List[str]:
    s = clean_text(x)
    if not s:
        return []
    vals = []
    for part in re.split(r"[|;, ]+", s):
        part = clean_accession(part)
        if part and part not in vals:
            vals.append(part)
    return vals


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


def make_session() -> requests.Session:
    s = requests.Session()
    retry = Retry(
        total=MAX_RETRIES,
        connect=MAX_RETRIES,
        read=MAX_RETRIES,
        backoff_factor=1.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=30, pool_maxsize=30)
    s.mount("http://", adapter)
    s.mount("https://", adapter)
    s.headers.update(HEADERS)
    return s


def chunks(items: List[str], size: int) -> Iterable[List[str]]:
    for i in range(0, len(items), size):
        yield items[i:i + size]


def count_tokens(text: str, separators: str = r"[;|]") -> int:
    text = clean_text(text)
    if not text:
        return 0
    parts = [p.strip() for p in re.split(separators, text) if p.strip()]
    return len(parts)


def count_feature_blocks(text: str) -> int:
    """
    Count UniProt feature blocks in REST TSV text fields.
    These often appear as semicolon-separated annotations.
    """
    text = clean_text(text)
    if not text:
        return 0

    # Count semicolon chunks, but if no semicolon just return 1.
    parts = [p.strip() for p in re.split(r";\s*", text) if p.strip()]
    return len(parts) if parts else 1


def has_any_keyword(text: str, words: List[str]) -> int:
    low = clean_text(text).lower()
    if not low:
        return 0
    return int(any(w.lower() in low for w in words))


def keyword_count(text: str, words: List[str]) -> int:
    low = clean_text(text).lower()
    if not low:
        return 0
    return int(sum(low.count(w.lower()) for w in words))


def sequence_features(seq: str) -> Dict[str, Any]:
    seq = clean_text(seq).upper().replace(" ", "").replace("\n", "")
    seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq)

    out = {
        "seq_length_from_sequence": len(seq),
        "seq_hydrophobic_fraction": np.nan,
        "seq_polar_fraction": np.nan,
        "seq_positive_fraction": np.nan,
        "seq_negative_fraction": np.nan,
        "seq_charged_fraction": np.nan,
        "seq_aromatic_fraction": np.nan,
        "seq_small_fraction": np.nan,
        "seq_proline_fraction": np.nan,
        "seq_cysteine_fraction": np.nan,
        "seq_glycine_fraction": np.nan,
        "seq_instability_proxy_low_complexity_fraction": np.nan,
    }

    n = len(seq)
    if n == 0:
        return out

    for group, letters in AA_GROUPS.items():
        out[f"seq_{group}_fraction"] = sum(1 for aa in seq if aa in letters) / n

    # Simple low-complexity proxy: most frequent amino-acid fraction.
    aa_counts = {}
    for aa in seq:
        aa_counts[aa] = aa_counts.get(aa, 0) + 1
    out["seq_instability_proxy_low_complexity_fraction"] = max(aa_counts.values()) / n if aa_counts else np.nan

    return out


def find_local_uniprot_tsv(dbdir: Path) -> Optional[Path]:
    if not dbdir.exists():
        return None

    candidates = []
    for pattern in ["*.tsv", "*.tsv.gz", "*uniprot*.tab", "*uniprot*.tab.gz", "*.txt", "*.txt.gz"]:
        candidates.extend(list(dbdir.glob(pattern)))

    ranked = []
    for p in candidates:
        name = p.name.lower()
        score = 0
        if "uniprot" in name:
            score += 5
        if "human" in name or "9606" in name:
            score += 5
        if p.suffix in {".tsv", ".gz", ".txt"}:
            score += 1
        ranked.append((score, p.stat().st_size, p))

    if not ranked:
        return None

    ranked.sort(reverse=True)
    return ranked[0][2]


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


def build_uniprot_maps(hgnc: pd.DataFrame, limit_genes: int = 0) -> Tuple[Dict[str, Set[str]], Dict[str, List[str]], List[str]]:
    sub = hgnc.copy()
    if limit_genes and limit_genes > 0:
        sub = sub.head(limit_genes).copy()

    uniprot_to_genes: Dict[str, Set[str]] = defaultdict(set)
    gene_to_uniprots: Dict[str, List[str]] = defaultdict(list)

    for _, row in sub.iterrows():
        gene = row["gene_symbol"]
        accs = split_uniprot_ids(row.get("uniprot_ids", ""))

        for acc in accs:
            uniprot_to_genes[acc].add(gene)
            if acc not in gene_to_uniprots[gene]:
                gene_to_uniprots[gene].append(acc)

    accessions = sorted(uniprot_to_genes.keys())

    log("=" * 100)
    log("[UNIPROT MAP]")
    log(f"[ACCESSIONS] {len(accessions)}")
    log(f"[GENES WITH ACCESSIONS] {len(gene_to_uniprots)}")

    return uniprot_to_genes, gene_to_uniprots, accessions


# =============================================================================
# UNIPROT DOWNLOAD / LOAD
# =============================================================================

def download_uniprot_rest(accessions: List[str], downloads_dir: Path, force: bool = False) -> pd.DataFrame:
    mkdir(downloads_dir)
    outpath = downloads_dir / "feature6_uniprot_rest_download.tsv"

    if outpath.exists() and outpath.stat().st_size > 0 and not force:
        log(f"[SKIP] Existing UniProt REST table: {outpath}")
        return pd.read_csv(outpath, sep="\t", dtype=str, low_memory=False)

    if not accessions:
        return pd.DataFrame()

    session = make_session()
    rows = []
    fields = ",".join(UNIPROT_FIELDS)
    basic_fields = ",".join(BASIC_UNIPROT_FIELDS)

    batches = list(chunks(accessions, BATCH_SIZE))

    log("=" * 100)
    log("[DOWNLOAD UNIPROT REST]")
    log(f"[ACCESSIONS] {len(accessions)}")
    log(f"[BATCH SIZE] {BATCH_SIZE}")
    log(f"[BATCHES] {len(batches)}")

    for i, batch in enumerate(batches, start=1):
        query = "(" + " OR ".join([f"accession:{acc}" for acc in batch]) + ") AND organism_id:9606"

        params = {
            "query": query,
            "format": "tsv",
            "fields": fields,
            "size": len(batch),
        }

        log(f"[UNIPROT] batch {i}/{len(batches)} n={len(batch)}")

        try:
            r = session.get(UNIPROT_SEARCH_URL, params=params, timeout=REQUEST_TIMEOUT)
            if r.status_code != 200:
                # Fallback to basic field list.
                params["fields"] = basic_fields
                r = session.get(UNIPROT_SEARCH_URL, params=params, timeout=REQUEST_TIMEOUT)

            if r.status_code == 200 and r.text.strip():
                df = pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str)
                rows.append(df)
            else:
                log(f"[WARNING] UniProt batch failed HTTP {r.status_code}: {r.text[:200]}")
        except Exception as exc:
            log(f"[WARNING] UniProt batch failed: {exc}")

        time.sleep(SLEEP_BETWEEN_BATCHES)

    if rows:
        out = pd.concat(rows, ignore_index=True).drop_duplicates()
    else:
        out = pd.DataFrame()

    out.to_csv(outpath, sep="\t", index=False)

    log(f"[SAVED] {outpath}")
    log(f"[SHAPE] {out.shape}")

    return out


def load_local_or_download_uniprot(
    dbdir: Path,
    downloads_dir: Path,
    accessions: List[str],
    no_api: bool,
    force_download: bool,
) -> pd.DataFrame:
    local = find_local_uniprot_tsv(dbdir)

    if local is not None:
        log("=" * 100)
        log(f"[LOCAL UNIPROT TABLE FOUND] {local}")
        try:
            df = pd.read_csv(local, sep="\t", dtype=str, low_memory=False)
            log(f"[LOCAL SHAPE] {df.shape}")
            return df
        except Exception as exc:
            log(f"[WARNING] Could not read local UniProt table: {exc}")

    if no_api:
        raise RuntimeError(
            f"No usable local UniProt TSV found in {dbdir}, and --no-api was used."
        )

    return download_uniprot_rest(accessions, downloads_dir, force=force_download)


# =============================================================================
# COLUMN NORMALIZATION
# =============================================================================

def normalise_uniprot_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize UniProt TSV column names from REST/local files.
    """
    if df.empty:
        return df

    original_cols = list(df.columns)
    colmap = {}

    for c in original_cols:
        low = c.strip().lower()

        if low in {"entry", "accession", "primary accession", "primaryaccession"}:
            colmap[c] = "uniprot_accession"
        elif low in {"entry name", "id"}:
            colmap[c] = "uniprot_entry_name"
        elif low == "reviewed":
            colmap[c] = "reviewed"
        elif low in {"protein names", "protein_name", "protein name"}:
            colmap[c] = "protein_name"
        elif low in {"gene names", "gene_names"}:
            colmap[c] = "gene_names"
        elif low in {"organism (id)", "organism_id", "organism id"}:
            colmap[c] = "organism_id"
        elif low in {"length", "protein length"}:
            colmap[c] = "protein_length"
        elif low in {"mass", "molecular weight"}:
            colmap[c] = "protein_mass"
        elif "subcellular location" in low:
            colmap[c] = "subcellular_location"
        elif low in {"transmembrane", "ft_transmem"} or "transmembrane" in low:
            colmap[c] = "ft_transmem"
        elif low in {"intramembrane", "ft_intramem"} or "intramembrane" in low:
            colmap[c] = "ft_intramem"
        elif "topological domain" in low or "topo" in low:
            colmap[c] = "ft_topo_dom"
        elif "signal peptide" in low or low == "ft_signal":
            colmap[c] = "ft_signal"
        elif low in {"domain [ft]", "domain", "ft_domain"}:
            colmap[c] = "ft_domain"
        elif low in {"region", "region [ft]", "ft_region"}:
            colmap[c] = "ft_region"
        elif low in {"motif", "motif [ft]", "ft_motif"}:
            colmap[c] = "ft_motif"
        elif "binding site" in low or low == "ft_binding":
            colmap[c] = "ft_binding"
        elif "active site" in low or low == "ft_act_site":
            colmap[c] = "ft_act_site"
        elif low in {"keywords", "keyword"}:
            colmap[c] = "keywords"
        elif low in {"gene ontology ids", "go"}:
            colmap[c] = "go"
        elif low in {"gene ontology (biological process)", "go_p", "go biological process"}:
            colmap[c] = "go_p"
        elif low in {"gene ontology (cellular component)", "go_c", "go cellular component"}:
            colmap[c] = "go_c"
        elif low in {"gene ontology (molecular function)", "go_f", "go molecular function"}:
            colmap[c] = "go_f"
        elif "function [cc]" in low or low == "cc_function":
            colmap[c] = "cc_function"
        elif "pathway" in low and "reactome" not in low:
            colmap[c] = "cc_pathway"
        elif "tissue specificity" in low:
            colmap[c] = "cc_tissue_specificity"
        elif "alternative products" in low or low == "cc_alternative_products":
            colmap[c] = "cc_alternative_products"
        elif low in {"sequence"}:
            colmap[c] = "sequence"

    df = df.rename(columns=colmap).copy()

    if "uniprot_accession" not in df.columns:
        raise RuntimeError(
            "Could not identify UniProt accession column in UniProt table. "
            f"Columns were: {original_cols[:30]}"
        )

    df["uniprot_accession"] = df["uniprot_accession"].map(clean_accession)

    for c in [
        "uniprot_entry_name",
        "reviewed",
        "protein_name",
        "gene_names",
        "organism_id",
        "protein_length",
        "protein_mass",
        "subcellular_location",
        "ft_transmem",
        "ft_intramem",
        "ft_topo_dom",
        "ft_signal",
        "ft_domain",
        "ft_region",
        "ft_motif",
        "ft_binding",
        "ft_act_site",
        "keywords",
        "go",
        "go_p",
        "go_c",
        "go_f",
        "cc_function",
        "cc_pathway",
        "cc_tissue_specificity",
        "cc_alternative_products",
        "sequence",
    ]:
        if c not in df.columns:
            df[c] = ""

    df = df.drop_duplicates(subset=["uniprot_accession"], keep="first").copy()

    return df


# =============================================================================
# FEATURE BUILDING
# =============================================================================

def build_protein_table(
    uniprot_df: pd.DataFrame,
    uniprot_to_genes: Dict[str, Set[str]],
    processed_dir: Path,
) -> pd.DataFrame:
    if uniprot_df.empty:
        out = pd.DataFrame()
        out.to_csv(processed_dir / "feature6_uniprot_protein_table.csv", index=False)
        return out

    df = normalise_uniprot_columns(uniprot_df)

    rows = []

    for _, row in df.iterrows():
        acc = clean_accession(row.get("uniprot_accession", ""))

        if acc not in uniprot_to_genes:
            continue

        text_all = " ; ".join([
            clean_text(row.get("protein_name", "")),
            clean_text(row.get("subcellular_location", "")),
            clean_text(row.get("keywords", "")),
            clean_text(row.get("go", "")),
            clean_text(row.get("go_p", "")),
            clean_text(row.get("go_c", "")),
            clean_text(row.get("go_f", "")),
            clean_text(row.get("cc_function", "")),
            clean_text(row.get("cc_pathway", "")),
            clean_text(row.get("cc_tissue_specificity", "")),
            clean_text(row.get("ft_domain", "")),
            clean_text(row.get("ft_region", "")),
            clean_text(row.get("ft_motif", "")),
        ])

        seq_feats = sequence_features(row.get("sequence", ""))

        for gene in uniprot_to_genes[acc]:
            rec = {
                "gene_symbol": gene,
                "uniprot_accession": acc,
                "uniprot_entry_name": clean_text(row.get("uniprot_entry_name", "")),
                "protein_name": clean_text(row.get("protein_name", "")),
                "gene_names": clean_text(row.get("gene_names", "")),
                "organism_id": clean_text(row.get("organism_id", "")),
                "reviewed": clean_text(row.get("reviewed", "")),
                "feature6_uniprot_is_reviewed": int("reviewed" in clean_text(row.get("reviewed", "")).lower() or "swiss" in clean_text(row.get("reviewed", "")).lower()),
                "feature6_uniprot_protein_length": safe_float(row.get("protein_length")),
                "feature6_uniprot_protein_mass": safe_float(row.get("protein_mass")),

                # Feature block counts.
                "feature6_uniprot_transmembrane_count": count_feature_blocks(row.get("ft_transmem", "")),
                "feature6_uniprot_intramembrane_count": count_feature_blocks(row.get("ft_intramem", "")),
                "feature6_uniprot_topological_domain_count": count_feature_blocks(row.get("ft_topo_dom", "")),
                "feature6_uniprot_signal_peptide_count": count_feature_blocks(row.get("ft_signal", "")),
                "feature6_uniprot_domain_feature_count": count_feature_blocks(row.get("ft_domain", "")),
                "feature6_uniprot_region_feature_count": count_feature_blocks(row.get("ft_region", "")),
                "feature6_uniprot_motif_feature_count": count_feature_blocks(row.get("ft_motif", "")),
                "feature6_uniprot_binding_site_count": count_feature_blocks(row.get("ft_binding", "")),
                "feature6_uniprot_active_site_count": count_feature_blocks(row.get("ft_act_site", "")),

                # GO/keyword counts.
                "feature6_uniprot_keyword_count": count_tokens(row.get("keywords", "")),
                "feature6_uniprot_go_total_count": count_tokens(row.get("go", "")),
                "feature6_uniprot_go_biological_process_count": count_tokens(row.get("go_p", "")),
                "feature6_uniprot_go_cellular_component_count": count_tokens(row.get("go_c", "")),
                "feature6_uniprot_go_molecular_function_count": count_tokens(row.get("go_f", "")),

                # Text lengths as annotation-richness proxies.
                "feature6_uniprot_function_text_length": len(clean_text(row.get("cc_function", ""))),
                "feature6_uniprot_pathway_text_length": len(clean_text(row.get("cc_pathway", ""))),
                "feature6_uniprot_tissue_specificity_text_length": len(clean_text(row.get("cc_tissue_specificity", ""))),
                "feature6_uniprot_has_alternative_products": int(bool(clean_text(row.get("cc_alternative_products", "")))),

                # Raw audit text.
                "feature6_uniprot_subcellular_location_text": clean_text(row.get("subcellular_location", "")),
                "feature6_uniprot_keywords_text": clean_text(row.get("keywords", "")),
                "feature6_uniprot_go_text": clean_text(row.get("go", "")),
            }

            # Subcellular flags/counts.
            for group, words in SUBCELLULAR_KEYWORDS.items():
                rec[f"feature6_uniprot_subcellular_has_{group}"] = has_any_keyword(
                    clean_text(row.get("subcellular_location", "")) + " ; " + text_all,
                    words,
                )
                rec[f"feature6_uniprot_subcellular_keyword_count_{group}"] = keyword_count(
                    clean_text(row.get("subcellular_location", "")) + " ; " + text_all,
                    words,
                )

            # Functional keyword flags/counts.
            for group, words in FUNCTION_KEYWORDS.items():
                rec[f"feature6_uniprot_function_has_{group}"] = has_any_keyword(text_all, words)
                rec[f"feature6_uniprot_function_keyword_count_{group}"] = keyword_count(text_all, words)

            # Derived membrane/secreted flags.
            rec["feature6_uniprot_has_signal_peptide"] = int(
                rec["feature6_uniprot_signal_peptide_count"] > 0
                or rec.get("feature6_uniprot_subcellular_has_secreted", 0) == 1
            )
            rec["feature6_uniprot_has_transmembrane"] = int(
                rec["feature6_uniprot_transmembrane_count"] > 0
                or rec.get("feature6_uniprot_subcellular_has_membrane", 0) == 1
            )

            for k, v in seq_feats.items():
                rec[f"feature6_uniprot_{k}"] = v

            rows.append(rec)

    out = pd.DataFrame(rows)

    if not out.empty:
        out = out.drop_duplicates()

    outpath = processed_dir / "feature6_uniprot_protein_table.csv"
    out.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED UNIPROT PROTEIN TABLE]")
    log(f"[PATH]  {outpath}")
    log(f"[SHAPE] {out.shape}")

    return out


def aggregate_gene_features(
    hgnc: pd.DataFrame,
    protein_table: pd.DataFrame,
    gene_to_uniprots: Dict[str, List[str]],
    processed_dir: Path,
) -> pd.DataFrame:
    genes = pd.DataFrame({"gene_symbol": sorted(hgnc["gene_symbol"].unique().tolist())})
    records = []

    by_gene = dict(tuple(protein_table.groupby("gene_symbol"))) if not protein_table.empty else {}

    for gene in genes["gene_symbol"]:
        g = by_gene.get(gene, pd.DataFrame())
        uniprots = gene_to_uniprots.get(gene, [])

        if g.empty:
            rec = {
                "gene_symbol": gene,
                "feature6_uniprot_has_annotation": 0,
                "feature6_uniprot_accession_count": len(uniprots),
                "feature6_uniprot_accessions": ";".join(uniprots),
            }
            records.append(rec)
            continue

        # Prefer reviewed and longest protein as representative.
        g2 = g.copy()
        g2["reviewed_sort"] = pd.to_numeric(g2["feature6_uniprot_is_reviewed"], errors="coerce").fillna(0)
        g2["length_sort"] = pd.to_numeric(g2["feature6_uniprot_protein_length"], errors="coerce").fillna(0)
        rep = g2.sort_values(["reviewed_sort", "length_sort"], ascending=False).iloc[0].to_dict()

        numeric_cols = [
            c for c in g.columns
            if c.startswith("feature6_uniprot_")
            and c not in {
                "feature6_uniprot_subcellular_location_text",
                "feature6_uniprot_keywords_text",
                "feature6_uniprot_go_text",
            }
            and pd.api.types.is_numeric_dtype(pd.to_numeric(g[c], errors="coerce"))
        ]

        rec = {
            "gene_symbol": gene,
            "feature6_uniprot_has_annotation": 1,
            "feature6_uniprot_accession_count": int(g["uniprot_accession"].nunique()),
            "feature6_uniprot_accessions": ";".join(sorted(set(g["uniprot_accession"].dropna().astype(str)))),
            "feature6_uniprot_representative_accession": rep.get("uniprot_accession", ""),
            "feature6_uniprot_representative_protein_name": rep.get("protein_name", ""),
            "feature6_uniprot_any_reviewed": int((pd.to_numeric(g["feature6_uniprot_is_reviewed"], errors="coerce").fillna(0) > 0).any()),
        }

        # Representative scalar features.
        rep_keep = [
            "feature6_uniprot_protein_length",
            "feature6_uniprot_protein_mass",
            "feature6_uniprot_transmembrane_count",
            "feature6_uniprot_intramembrane_count",
            "feature6_uniprot_topological_domain_count",
            "feature6_uniprot_signal_peptide_count",
            "feature6_uniprot_domain_feature_count",
            "feature6_uniprot_region_feature_count",
            "feature6_uniprot_motif_feature_count",
            "feature6_uniprot_binding_site_count",
            "feature6_uniprot_active_site_count",
            "feature6_uniprot_keyword_count",
            "feature6_uniprot_go_total_count",
            "feature6_uniprot_go_biological_process_count",
            "feature6_uniprot_go_cellular_component_count",
            "feature6_uniprot_go_molecular_function_count",
            "feature6_uniprot_function_text_length",
            "feature6_uniprot_pathway_text_length",
            "feature6_uniprot_tissue_specificity_text_length",
            "feature6_uniprot_has_alternative_products",
            "feature6_uniprot_has_signal_peptide",
            "feature6_uniprot_has_transmembrane",
            "feature6_uniprot_seq_length_from_sequence",
            "feature6_uniprot_seq_hydrophobic_fraction",
            "feature6_uniprot_seq_polar_fraction",
            "feature6_uniprot_seq_positive_fraction",
            "feature6_uniprot_seq_negative_fraction",
            "feature6_uniprot_seq_charged_fraction",
            "feature6_uniprot_seq_aromatic_fraction",
            "feature6_uniprot_seq_small_fraction",
            "feature6_uniprot_seq_proline_fraction",
            "feature6_uniprot_seq_cysteine_fraction",
            "feature6_uniprot_seq_glycine_fraction",
            "feature6_uniprot_seq_instability_proxy_low_complexity_fraction",
        ]

        for c in rep_keep:
            if c in g.columns:
                rec[c] = rep.get(c, np.nan)

        # For all has_* columns, use max across isoforms/accessions.
        has_cols = [c for c in g.columns if c.startswith("feature6_uniprot_") and "_has_" in c]
        for c in has_cols:
            vals = pd.to_numeric(g[c], errors="coerce").fillna(0)
            rec[c] = int(vals.max()) if len(vals) else 0

        # For count columns, use max and sum variants.
        count_cols = [c for c in g.columns if c.startswith("feature6_uniprot_") and c.endswith("_count")]
        for c in count_cols:
            vals = pd.to_numeric(g[c], errors="coerce").dropna()
            if len(vals):
                rec[f"{c}_max_across_uniprot"] = float(vals.max())
                rec[f"{c}_sum_across_uniprot"] = float(vals.sum())

        records.append(rec)

    features = pd.DataFrame(records)

    # Fill count/flag columns.
    for c in features.columns:
        if c.startswith("feature6_uniprot_") and (
            "_count" in c
            or "_has_" in c
            or c.endswith("_reviewed")
            or c.startswith("feature6_uniprot_any_")
        ):
            if features[c].dtype != object:
                features[c] = features[c].fillna(0)

    text_cols = [
        "feature6_uniprot_accessions",
        "feature6_uniprot_representative_accession",
        "feature6_uniprot_representative_protein_name",
    ]
    for c in text_cols:
        if c in features.columns:
            features[c] = features[c].fillna("")

    outpath = processed_dir / "feature6_uniprot_gene_features.csv"
    features.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED GENE FEATURES]")
    log(f"[PATH]  {outpath}")
    log(f"[SHAPE] {features.shape}")

    return features


def merge_with_hgnc(hgnc: pd.DataFrame, features: pd.DataFrame, processed_dir: Path) -> pd.DataFrame:
    merged = hgnc.merge(features, on="gene_symbol", how="left")

    feature_cols = [c for c in features.columns if c != "gene_symbol"]
    merged["feature6_uniprot_has_any_feature"] = merged[feature_cols].notna().any(axis=1).astype(int)
    merged["feature6_uniprot_n_nonmissing_features"] = merged[feature_cols].notna().sum(axis=1)

    outpath = processed_dir / "feature6_uniprot_hgnc_merged.csv"
    merged.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED HGNC MERGED]")
    log(f"[PATH]  {outpath}")
    log(f"[SHAPE] {merged.shape}")

    return merged


def write_summary(
    outdir: Path,
    hgnc: pd.DataFrame,
    protein_table: pd.DataFrame,
    features: pd.DataFrame,
    merged: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    path = outdir / "feature6_uniprot_summary.txt"

    lines = []
    lines.append("Feature 6: UniProt protein annotation features")
    lines.append("=" * 80)
    lines.append(f"Created at: {now_iso()}")
    lines.append(f"HGNC rows: {hgnc.shape[0]}")
    lines.append(f"Protein table rows: {protein_table.shape[0]}")
    lines.append(f"Gene feature rows: {features.shape[0]}")
    lines.append(f"Gene feature columns: {features.shape[1]}")
    lines.append(f"Merged rows: {merged.shape[0]}")
    lines.append(f"Merged columns: {merged.shape[1]}")
    lines.append("")
    lines.append("Input:")
    lines.append(f"HGNC: {Path(args.hgnc).resolve()}")
    lines.append(f"UniProt dbdir: {Path(args.dbdir).resolve()}")
    lines.append("")
    lines.append("Coverage:")
    if "feature6_uniprot_has_annotation" in features.columns:
        lines.append(f"UniProt annotation coverage: {features['feature6_uniprot_has_annotation'].mean():.4f}")
    if "feature6_uniprot_any_reviewed" in features.columns:
        lines.append(f"Reviewed/Swiss-Prot coverage: {features['feature6_uniprot_any_reviewed'].mean():.4f}")
    if "feature6_uniprot_has_transmembrane" in features.columns:
        lines.append(f"Transmembrane fraction: {features['feature6_uniprot_has_transmembrane'].mean():.4f}")
    if "feature6_uniprot_has_signal_peptide" in features.columns:
        lines.append(f"Signal/secreted fraction: {features['feature6_uniprot_has_signal_peptide'].mean():.4f}")
    lines.append("")
    lines.append("Leakage policy:")
    lines.append("Included: UniProt sequence, localization, GO, keyword, topology and annotation-richness features.")
    lines.append("Excluded: known drug-target labels, approved-drug annotations, ChEMBL/DrugBank/DGIdb labels, target development level.")

    path.write_text("\n".join(lines) + "\n")
    log(f"[SUMMARY] {path}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build UniProt protein annotation features for HGNC genes."
    )
    parser.add_argument("--hgnc", default=str(DEFAULT_HGNC_PATH), help="Path to HGNC complete set.")
    parser.add_argument("--dbdir", default=str(DEFAULT_DBDIR), help="Local UniProt database directory.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory.")
    parser.add_argument("--all-hgnc-genes", action="store_true", help="Use all HGNC rows, not only protein-coding.")
    parser.add_argument("--limit-genes", type=int, default=0, help="Debug only: first N genes.")
    parser.add_argument("--no-api", action="store_true", help="Do not query UniProt REST if local file not found.")
    parser.add_argument("--force-download", action="store_true", help="Force UniProt REST re-download.")
    args = parser.parse_args()

    outdir = mkdir(Path(args.outdir))
    downloads_dir = mkdir(outdir / "downloads")
    processed_dir = mkdir(outdir / "processed")
    dbdir = Path(args.dbdir)

    log("=" * 100)
    log("FEATURE 6: UNIPROT PROTEIN ANNOTATION")
    log("=" * 100)
    log(f"[HGNC]            {args.hgnc}")
    log(f"[DBDIR]           {dbdir.resolve()}")
    log(f"[OUTDIR]          {outdir.resolve()}")
    log(f"[LIMIT GENES]     {args.limit_genes if args.limit_genes else 'none'}")
    log(f"[NO API]          {args.no_api}")
    log(f"[FORCE DOWNLOAD]  {args.force_download}")
    log("=" * 100)

    hgnc = load_hgnc(
        Path(args.hgnc),
        protein_coding_only=not args.all_hgnc_genes,
    )

    uniprot_to_genes, gene_to_uniprots, accessions = build_uniprot_maps(
        hgnc,
        limit_genes=args.limit_genes,
    )

    uniprot_raw = load_local_or_download_uniprot(
        dbdir=dbdir,
        downloads_dir=downloads_dir,
        accessions=accessions,
        no_api=args.no_api,
        force_download=args.force_download,
    )

    protein_table = build_protein_table(
        uniprot_df=uniprot_raw,
        uniprot_to_genes=uniprot_to_genes,
        processed_dir=processed_dir,
    )

    features = aggregate_gene_features(
        hgnc=hgnc.head(args.limit_genes).copy() if args.limit_genes else hgnc,
        protein_table=protein_table,
        gene_to_uniprots=gene_to_uniprots,
        processed_dir=processed_dir,
    )

    merged = merge_with_hgnc(
        hgnc=hgnc.head(args.limit_genes).copy() if args.limit_genes else hgnc,
        features=features,
        processed_dir=processed_dir,
    )

    metadata = {
        "created_at": now_iso(),
        "hgnc": str(Path(args.hgnc)),
        "dbdir": str(dbdir.resolve()),
        "outdir": str(outdir.resolve()),
        "protein_coding_only": not args.all_hgnc_genes,
        "limit_genes": args.limit_genes,
        "no_api": args.no_api,
        "force_download": args.force_download,
        "outputs": {
            "protein_table": str(processed_dir / "feature6_uniprot_protein_table.csv"),
            "gene_features": str(processed_dir / "feature6_uniprot_gene_features.csv"),
            "hgnc_merged": str(processed_dir / "feature6_uniprot_hgnc_merged.csv"),
        },
        "leakage_policy": {
            "included": [
                "UniProt protein annotation",
                "subcellular location",
                "sequence composition",
                "GO annotation counts",
                "keyword counts",
                "transmembrane/signal peptide/feature counts",
            ],
            "excluded": [
                "known drug-target labels",
                "approved-drug annotations",
                "ChEMBL labels",
                "DrugBank labels",
                "DGIdb labels",
                "target development level",
            ],
        },
    }

    with open(outdir / "feature6_uniprot_run_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    write_summary(outdir, hgnc, protein_table, features, merged, args)

    log("=" * 100)
    log("[DONE]")
    log(f"[PROTEIN TABLE] {processed_dir / 'feature6_uniprot_protein_table.csv'}")
    log(f"[GENE FEATURES] {processed_dir / 'feature6_uniprot_gene_features.csv'}")
    log(f"[HGNC MERGED]   {processed_dir / 'feature6_uniprot_hgnc_merged.csv'}")
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