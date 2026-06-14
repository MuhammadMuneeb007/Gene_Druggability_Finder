#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature11_ProteinSequence.py

Feature 11: Protein sequence composition and physicochemical features.

Purpose
-------
Build leakage-safe gene-level protein-sequence features from local Ensembl peptide FASTA.

Main input:
    feature_databases/Ensembl/Homo_sapiens.GRCh38.pep.all.fa.gz

HGNC:
    databases/HGNC/hgnc_complete_set.txt

Outputs:
    feature11_protein_sequence/processed/feature11_protein_sequence_long.csv
    feature11_protein_sequence/processed/feature11_protein_sequence_gene_features.csv
    feature11_protein_sequence/processed/feature11_protein_sequence_hgnc_merged.csv
    feature11_protein_sequence/feature11_protein_sequence_summary.txt
    feature11_protein_sequence/feature11_protein_sequence_run_metadata.json

Run:
    python Feature11_ProteinSequence.py

Fast test:
    python Feature11_ProteinSequence.py --max-records 5000

Manual FASTA:
    python Feature11_ProteinSequence.py \
      --pep-file feature_databases/Ensembl/Homo_sapiens.GRCh38.pep.all.fa.gz

Leakage policy:
    Safe. Uses only protein amino-acid sequence.
    Does not use ChEMBL, DrugBank, DGIdb, Open Targets, Pharos, clinical labels, or known drug-target labels.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# =============================================================================
# DEFAULTS
# =============================================================================

DEFAULT_HGNC_PATH = Path("databases") / "HGNC" / "hgnc_complete_set.txt"
DEFAULT_DBDIR = Path("feature_databases") / "Ensembl"
DEFAULT_OUTDIR = Path("feature11_protein_sequence")

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")

AA_GROUPS = {
    "hydrophobic": set("AVILMFWY"),
    "aliphatic": set("AVIL"),
    "aromatic": set("FWYH"),
    "positive": set("KRH"),
    "negative": set("DE"),
    "charged": set("KRHDE"),
    "polar": set("STNQCY"),
    "small": set("GASCTP"),
    "sulfur": set("CM"),
    "proline": set("P"),
    "glycine": set("G"),
    "cysteine": set("C"),
}

# Kyte-Doolittle hydropathy scale.
KD_HYDROPATHY = {
    "I": 4.5,
    "V": 4.2,
    "L": 3.8,
    "F": 2.8,
    "C": 2.5,
    "M": 1.9,
    "A": 1.8,
    "G": -0.4,
    "T": -0.7,
    "S": -0.8,
    "W": -0.9,
    "Y": -1.3,
    "P": -1.6,
    "H": -3.2,
    "E": -3.5,
    "Q": -3.5,
    "D": -3.5,
    "N": -3.5,
    "K": -3.9,
    "R": -4.5,
}

# Approximate residue masses, average residue mass after peptide bond formation.
AA_MASS = {
    "A": 89.09,
    "R": 174.20,
    "N": 132.12,
    "D": 133.10,
    "C": 121.15,
    "Q": 146.15,
    "E": 147.13,
    "G": 75.07,
    "H": 155.16,
    "I": 131.17,
    "L": 131.17,
    "K": 146.19,
    "M": 149.21,
    "F": 165.19,
    "P": 115.13,
    "S": 105.09,
    "T": 119.12,
    "W": 204.23,
    "Y": 181.19,
    "V": 117.15,
}

# Very simple pKa values for approximate charge/pI proxy.
PKA_NTERM = 9.69
PKA_CTERM = 2.34
PKA_SIDECHAIN_POS = {
    "K": 10.5,
    "R": 12.4,
    "H": 6.0,
}
PKA_SIDECHAIN_NEG = {
    "D": 3.9,
    "E": 4.1,
    "C": 8.3,
    "Y": 10.1,
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
    return clean_text(x).upper()


def clean_ensembl_gene_id(x: Any) -> str:
    s = clean_text(x)
    if not s:
        return ""
    s = s.replace("gene:", "")
    return s.split(".")[0]


def clean_ensembl_transcript_id(x: Any) -> str:
    s = clean_text(x)
    if not s:
        return ""
    s = s.replace("transcript:", "")
    return s.split(".")[0]


def clean_ensembl_protein_id(x: Any) -> str:
    s = clean_text(x)
    if not s:
        return ""
    s = s.replace("protein:", "")
    return s.split(".")[0]


def safe_float(x: Any) -> float:
    try:
        if x is None or pd.isna(x):
            return np.nan
        s = str(x).strip()
        if s in {"", ".", "NA", "NaN", "nan", "None"}:
            return np.nan
        return float(s)
    except Exception:
        return np.nan


def open_text_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", errors="ignore")
    return open(path, "rt", errors="ignore")


def find_pep_fasta(dbdir: Path) -> Optional[Path]:
    if not dbdir.exists():
        return None

    files = []
    for pattern in ["*pep*.fa.gz", "*pep*.fasta.gz", "*pep*.fa", "*pep*.fasta"]:
        files.extend(list(dbdir.rglob(pattern)))

    files = [p for p in files if p.exists() and p.stat().st_size > 0]

    if not files:
        return None

    ranked = []
    for p in files:
        name = p.name.lower()
        score = 0

        if "homo_sapiens" in name:
            score += 20
        if "grch38" in name:
            score += 10
        if "pep.all" in name:
            score += 20
        if name.endswith(".fa.gz") or name.endswith(".fasta.gz"):
            score += 5

        ranked.append((score, p.stat().st_size, p))

    ranked.sort(reverse=True)

    return ranked[0][2]


# =============================================================================
# HGNC
# =============================================================================

def load_hgnc(hgnc_path: Path, protein_coding_only: bool = True, limit_genes: int = 0) -> pd.DataFrame:
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

    if "ensembl_gene_id" not in hgnc.columns:
        hgnc["ensembl_gene_id"] = ""
    if "entrez_id" not in hgnc.columns:
        hgnc["entrez_id"] = ""
    if "uniprot_ids" not in hgnc.columns:
        hgnc["uniprot_ids"] = ""
    if "name" not in hgnc.columns:
        hgnc["name"] = ""

    hgnc["ensembl_gene_id_clean"] = hgnc["ensembl_gene_id"].map(clean_ensembl_gene_id)
    hgnc = hgnc.drop_duplicates(subset=["gene_symbol"], keep="first").copy()

    if limit_genes and limit_genes > 0:
        hgnc = hgnc.head(limit_genes).copy()

    log(f"[HGNC ROWS] {before} -> {len(hgnc)}")
    log(f"[PROTEIN CODING ONLY] {protein_coding_only}")
    log(f"[GENES WITH ENSEMBL] {(hgnc['ensembl_gene_id_clean'].astype(str).str.len() > 0).sum()}")

    return hgnc


def build_ensembl_to_gene(hgnc: pd.DataFrame) -> Dict[str, str]:
    mapping = {}

    for _, row in hgnc.iterrows():
        gene = row["gene_symbol"]
        ens = clean_ensembl_gene_id(row.get("ensembl_gene_id_clean", ""))

        if ens and gene:
            mapping[ens] = gene

    return mapping


# =============================================================================
# FASTA PARSING
# =============================================================================

def parse_ensembl_fasta_header(header: str) -> Dict[str, str]:
    """
    Example Ensembl peptide FASTA header:

    >ENSP00000354587.4 pep chromosome:GRCh38:17:7661779:7687550:1 gene:ENSG00000141510.18 transcript:ENST00000269305.9 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:TP53 description:p53 tumor suppressor [Source:HGNC Symbol;Acc:HGNC:11998]

    We parse:
        protein_id
        gene_id
        transcript_id
        gene_symbol
        gene_biotype
        transcript_biotype
    """
    h = header.strip()
    if h.startswith(">"):
        h = h[1:]

    parts = h.split()
    protein_id = clean_ensembl_protein_id(parts[0]) if parts else ""

    def get_tag(tag: str) -> str:
        m = re.search(rf"{re.escape(tag)}:([^\s]+)", h)
        return clean_text(m.group(1)) if m else ""

    gene_id = clean_ensembl_gene_id(get_tag("gene"))
    transcript_id = clean_ensembl_transcript_id(get_tag("transcript"))
    gene_symbol = normalize_symbol(get_tag("gene_symbol"))
    gene_biotype = clean_text(get_tag("gene_biotype"))
    transcript_biotype = clean_text(get_tag("transcript_biotype"))

    return {
        "protein_id": protein_id,
        "ensembl_gene_id": gene_id,
        "transcript_id": transcript_id,
        "gene_symbol_header": gene_symbol,
        "gene_biotype": gene_biotype,
        "transcript_biotype": transcript_biotype,
        "raw_header": header.strip(),
    }


def iter_fasta_records(path: Path, max_records: int = 0):
    header = None
    seq_chunks = []
    n = 0

    with open_text_maybe_gzip(path) as f:
        for line in f:
            line = line.rstrip("\n")

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    n += 1
                    yield header, "".join(seq_chunks)

                    if max_records and n >= max_records:
                        return

                header = line
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())

        if header is not None:
            n += 1
            yield header, "".join(seq_chunks)


# =============================================================================
# SEQUENCE FEATURES
# =============================================================================

def aa_composition(seq: str) -> Dict[str, float]:
    seq = seq.upper()
    n = len(seq)
    counts = Counter(seq)

    out = {}

    for aa in AA_LIST:
        out[f"aa_frac_{aa}"] = counts.get(aa, 0) / n if n > 0 else np.nan
        out[f"aa_count_{aa}"] = counts.get(aa, 0)

    return out


def group_fraction(seq: str, group: set) -> float:
    if not seq:
        return np.nan
    return sum(1 for aa in seq if aa in group) / len(seq)


def shannon_entropy(seq: str) -> float:
    if not seq:
        return np.nan

    counts = Counter(seq)
    n = len(seq)
    ent = 0.0

    for aa, c in counts.items():
        p = c / n
        if p > 0:
            ent -= p * math.log2(p)

    return ent


def normalized_entropy(seq: str) -> float:
    ent = shannon_entropy(seq)
    if pd.isna(ent):
        return np.nan
    return ent / math.log2(20)


def longest_run_same_aa(seq: str) -> int:
    if not seq:
        return 0

    best = 1
    cur = 1

    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1

    return best


def count_homopolymer_runs(seq: str, min_len: int = 5) -> int:
    if not seq:
        return 0

    count = 0
    cur = 1

    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
        else:
            if cur >= min_len:
                count += 1
            cur = 1

    if cur >= min_len:
        count += 1

    return count


def low_complexity_fraction(seq: str, window: int = 15, unique_threshold: int = 5) -> float:
    """
    Simple low-complexity proxy:
    fraction of residues belonging to windows with <= unique_threshold amino-acid types.
    """
    n = len(seq)

    if n == 0:
        return np.nan

    if n < window:
        return float(len(set(seq)) <= unique_threshold)

    marked = np.zeros(n, dtype=bool)

    for i in range(0, n - window + 1):
        w = seq[i:i + window]
        if len(set(w)) <= unique_threshold:
            marked[i:i + window] = True

    return float(marked.mean())


def hydropathy_values(seq: str) -> List[float]:
    return [KD_HYDROPATHY[aa] for aa in seq if aa in KD_HYDROPATHY]


def sliding_mean(values: List[float], window: int) -> List[float]:
    if len(values) < window:
        return []

    arr = np.asarray(values, dtype=float)
    csum = np.cumsum(np.insert(arr, 0, 0.0))
    return ((csum[window:] - csum[:-window]) / window).tolist()


def hydrophobic_segment_features(seq: str) -> Dict[str, Any]:
    """
    Transmembrane-like proxy:
    long hydrophobic windows with high average Kyte-Doolittle score.
    This is NOT a true TM predictor, just a sequence feature.
    """
    vals = hydropathy_values(seq)

    if not vals:
        return {
            "hydropathy_mean": np.nan,
            "hydropathy_median": np.nan,
            "hydropathy_max_window_19": np.nan,
            "hydropathy_max_window_21": np.nan,
            "tm_like_window_count_19": 0,
            "tm_like_window_count_21": 0,
            "has_tm_like_segment_proxy": 0,
        }

    mean_h = float(np.mean(vals))
    med_h = float(np.median(vals))

    win19 = sliding_mean(vals, 19)
    win21 = sliding_mean(vals, 21)

    max19 = float(np.max(win19)) if win19 else np.nan
    max21 = float(np.max(win21)) if win21 else np.nan

    tm19 = int(sum(1 for v in win19 if v >= 1.6))
    tm21 = int(sum(1 for v in win21 if v >= 1.6))

    return {
        "hydropathy_mean": mean_h,
        "hydropathy_median": med_h,
        "hydropathy_max_window_19": max19,
        "hydropathy_max_window_21": max21,
        "tm_like_window_count_19": tm19,
        "tm_like_window_count_21": tm21,
        "has_tm_like_segment_proxy": int(tm19 > 0 or tm21 > 0),
    }


def molecular_weight(seq: str) -> float:
    if not seq:
        return np.nan

    total = 18.015  # add water for termini

    for aa in seq:
        total += AA_MASS.get(aa, 0.0) - 18.015

    return float(total)


def charge_at_ph(seq: str, ph: float) -> float:
    """
    Approximate net charge at a given pH.
    """
    if not seq:
        return np.nan

    counts = Counter(seq)

    positive = 0.0
    negative = 0.0

    # N terminus.
    positive += 1.0 / (1.0 + 10 ** (ph - PKA_NTERM))

    # C terminus.
    negative += 1.0 / (1.0 + 10 ** (PKA_CTERM - ph))

    for aa, pka in PKA_SIDECHAIN_POS.items():
        positive += counts.get(aa, 0) / (1.0 + 10 ** (ph - pka))

    for aa, pka in PKA_SIDECHAIN_NEG.items():
        negative += counts.get(aa, 0) / (1.0 + 10 ** (pka - ph))

    return float(positive - negative)


def approximate_pI(seq: str) -> float:
    if not seq:
        return np.nan

    low = 0.0
    high = 14.0

    for _ in range(60):
        mid = (low + high) / 2
        charge = charge_at_ph(seq, mid)

        if charge > 0:
            low = mid
        else:
            high = mid

    return float((low + high) / 2)


def kmer_diversity(seq: str, k: int) -> float:
    if len(seq) < k:
        return np.nan

    kmers = [seq[i:i + k] for i in range(0, len(seq) - k + 1)]
    return float(len(set(kmers)) / len(kmers)) if kmers else np.nan


def sequence_features(seq: str) -> Dict[str, Any]:
    seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq.upper())
    n = len(seq)

    out = {
        "protein_length": n,
        "log1p_protein_length": float(np.log1p(n)) if n > 0 else np.nan,
    }

    if n == 0:
        for aa in AA_LIST:
            out[f"aa_frac_{aa}"] = np.nan
            out[f"aa_count_{aa}"] = 0

        for group in AA_GROUPS:
            out[f"group_frac_{group}"] = np.nan

        out.update(
            {
                "sequence_entropy": np.nan,
                "sequence_entropy_normalized": np.nan,
                "longest_same_aa_run": 0,
                "homopolymer_run_count_ge_5": 0,
                "low_complexity_fraction_w15": np.nan,
                "molecular_weight": np.nan,
                "net_charge_ph7": np.nan,
                "absolute_net_charge_ph7": np.nan,
                "estimated_pI": np.nan,
                "kmer_diversity_2": np.nan,
                "kmer_diversity_3": np.nan,
            }
        )
        out.update(hydrophobic_segment_features(seq))
        return out

    out.update(aa_composition(seq))

    for group_name, group in AA_GROUPS.items():
        out[f"group_frac_{group_name}"] = group_fraction(seq, group)

    out["sequence_entropy"] = shannon_entropy(seq)
    out["sequence_entropy_normalized"] = normalized_entropy(seq)

    out["longest_same_aa_run"] = longest_run_same_aa(seq)
    out["homopolymer_run_count_ge_5"] = count_homopolymer_runs(seq, min_len=5)
    out["low_complexity_fraction_w15"] = low_complexity_fraction(seq, window=15, unique_threshold=5)

    out["molecular_weight"] = molecular_weight(seq)
    out["net_charge_ph7"] = charge_at_ph(seq, 7.0)
    out["absolute_net_charge_ph7"] = abs(out["net_charge_ph7"])
    out["estimated_pI"] = approximate_pI(seq)

    out["kmer_diversity_2"] = kmer_diversity(seq, 2)
    out["kmer_diversity_3"] = kmer_diversity(seq, 3)

    out.update(hydrophobic_segment_features(seq))

    # Some interpretable flags.
    out["is_short_protein_lt_100aa"] = int(n < 100)
    out["is_long_protein_gt_1000aa"] = int(n > 1000)
    out["is_very_long_protein_gt_2000aa"] = int(n > 2000)
    out["high_cysteine_fraction_ge_0_05"] = int(out.get("group_frac_cysteine", 0) >= 0.05)
    out["high_positive_fraction_ge_0_15"] = int(out.get("group_frac_positive", 0) >= 0.15)
    out["high_negative_fraction_ge_0_15"] = int(out.get("group_frac_negative", 0) >= 0.15)
    out["high_low_complexity_fraction_ge_0_20"] = int(out.get("low_complexity_fraction_w15", 0) >= 0.20)

    return out


# =============================================================================
# FEATURE BUILDING
# =============================================================================

def build_protein_long_table(
    pep_file: Path,
    ensembl_to_gene: Dict[str, str],
    processed_dir: Path,
    max_records: int = 0,
) -> pd.DataFrame:
    if not pep_file.exists():
        raise FileNotFoundError(f"Peptide FASTA not found: {pep_file}")

    log("=" * 100)
    log(f"[STREAM ENSEMBL PEPTIDE FASTA] {pep_file}")

    rows = []
    total = 0
    matched = 0
    t0 = time.time()

    for header, seq in iter_fasta_records(pep_file, max_records=max_records):
        total += 1

        if total % 100000 == 0:
            elapsed = (time.time() - t0) / 60
            log(f"[FASTA] total={total:,} matched={matched:,} elapsed={elapsed:.1f} min")

        meta = parse_ensembl_fasta_header(header)
        gene_id = meta["ensembl_gene_id"]
        gene_symbol = ensembl_to_gene.get(gene_id, "")

        if not gene_symbol:
            # fallback if HGNC symbol appears in header
            header_symbol = meta.get("gene_symbol_header", "")
            if header_symbol and header_symbol in set(ensembl_to_gene.values()):
                gene_symbol = header_symbol

        if not gene_symbol:
            continue

        feats = sequence_features(seq)

        rec = {
            "gene_symbol": gene_symbol,
            "ensembl_gene_id": gene_id,
            "ensembl_transcript_id": meta["transcript_id"],
            "ensembl_protein_id": meta["protein_id"],
            "gene_symbol_header": meta["gene_symbol_header"],
            "gene_biotype": meta["gene_biotype"],
            "transcript_biotype": meta["transcript_biotype"],
        }

        for k, v in feats.items():
            rec[f"feature11_{k}"] = v

        rows.append(rec)
        matched += 1

    long_df = pd.DataFrame(rows)

    if not long_df.empty:
        long_df = long_df.drop_duplicates()

    outpath = processed_dir / "feature11_protein_sequence_long.csv"
    long_df.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED PROTEIN LONG TABLE]")
    log(f"[TOTAL FASTA RECORDS SCANNED] {total:,}")
    log(f"[MATCHED PROTEIN RECORDS] {matched:,}")
    log(f"[PATH] {outpath}")
    log(f"[SHAPE] {long_df.shape}")

    return long_df


def aggregate_gene_features(
    hgnc: pd.DataFrame,
    protein_long: pd.DataFrame,
    processed_dir: Path,
) -> pd.DataFrame:
    genes = pd.DataFrame({"gene_symbol": sorted(hgnc["gene_symbol"].unique().tolist())})
    records = []

    by_gene = dict(tuple(protein_long.groupby("gene_symbol"))) if not protein_long.empty else {}

    numeric_feature_cols = [
        c for c in protein_long.columns
        if c.startswith("feature11_")
    ] if not protein_long.empty else []

    for gene in genes["gene_symbol"]:
        g = by_gene.get(gene, pd.DataFrame())

        if g.empty:
            records.append(
                {
                    "gene_symbol": gene,
                    "feature11_has_protein_sequence": 0,
                    "feature11_protein_isoform_count": 0,
                }
            )
            continue

        rec = {
            "gene_symbol": gene,
            "feature11_has_protein_sequence": 1,
            "feature11_protein_isoform_count": int(len(g)),
            "feature11_unique_transcript_count": int(g["ensembl_transcript_id"].replace("", np.nan).dropna().nunique()) if "ensembl_transcript_id" in g.columns else 0,
            "feature11_unique_protein_count": int(g["ensembl_protein_id"].replace("", np.nan).dropna().nunique()) if "ensembl_protein_id" in g.columns else 0,
            "feature11_has_multiple_protein_isoforms": int(len(g) > 1),
        }

        for c in numeric_feature_cols:
            vals = pd.to_numeric(g[c], errors="coerce").dropna()

            if len(vals) == 0:
                continue

            base = c

            rec[f"{base}_mean"] = float(vals.mean())
            rec[f"{base}_median"] = float(vals.median())
            rec[f"{base}_min"] = float(vals.min())
            rec[f"{base}_max"] = float(vals.max())

            if len(vals) > 1:
                rec[f"{base}_std"] = float(vals.std(ddof=1))
            else:
                rec[f"{base}_std"] = 0.0

        # Focused simpler aliases for most important model features.
        important = [
            "feature11_protein_length",
            "feature11_hydropathy_mean",
            "feature11_hydropathy_max_window_19",
            "feature11_hydropathy_max_window_21",
            "feature11_tm_like_window_count_19",
            "feature11_has_tm_like_segment_proxy",
            "feature11_group_frac_hydrophobic",
            "feature11_group_frac_aromatic",
            "feature11_group_frac_positive",
            "feature11_group_frac_negative",
            "feature11_group_frac_charged",
            "feature11_group_frac_polar",
            "feature11_group_frac_cysteine",
            "feature11_sequence_entropy_normalized",
            "feature11_low_complexity_fraction_w15",
            "feature11_molecular_weight",
            "feature11_net_charge_ph7",
            "feature11_absolute_net_charge_ph7",
            "feature11_estimated_pI",
            "feature11_kmer_diversity_2",
            "feature11_kmer_diversity_3",
        ]

        for c in important:
            if c in g.columns:
                vals = pd.to_numeric(g[c], errors="coerce").dropna()
                if len(vals):
                    alias = c.replace("feature11_", "feature11_summary_")
                    rec[f"{alias}_mean"] = float(vals.mean())
                    rec[f"{alias}_max"] = float(vals.max())
                    rec[f"{alias}_min"] = float(vals.min())

        records.append(rec)

    features = pd.DataFrame(records)

    for c in features.columns:
        if c.startswith("feature11_") and (
            "_count" in c
            or "_has_" in c
            or c.startswith("feature11_has_")
            or c.endswith("_proxy")
            or "_ge_" in c
            or "_gt_" in c
            or "_lt_" in c
        ):
            if features[c].dtype != object:
                features[c] = features[c].fillna(0)

    outpath = processed_dir / "feature11_protein_sequence_gene_features.csv"
    features.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED GENE FEATURES]")
    log(f"[PATH]  {outpath}")
    log(f"[SHAPE] {features.shape}")

    return features


def merge_with_hgnc(
    hgnc: pd.DataFrame,
    features: pd.DataFrame,
    processed_dir: Path,
) -> pd.DataFrame:
    merged = hgnc.merge(features, on="gene_symbol", how="left")

    feature_cols = [c for c in features.columns if c != "gene_symbol"]

    merged["feature11_protein_sequence_has_any_feature"] = merged[feature_cols].notna().any(axis=1).astype(int)
    merged["feature11_protein_sequence_n_nonmissing_features"] = merged[feature_cols].notna().sum(axis=1)

    outpath = processed_dir / "feature11_protein_sequence_hgnc_merged.csv"
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
    features: pd.DataFrame,
    merged: pd.DataFrame,
    pep_file: Path,
    args: argparse.Namespace,
) -> None:
    path = outdir / "feature11_protein_sequence_summary.txt"

    lines = []
    lines.append("Feature 11: Protein sequence composition and physicochemical features")
    lines.append("=" * 80)
    lines.append(f"Created at: {now_iso()}")
    lines.append(f"HGNC: {args.hgnc}")
    lines.append(f"Peptide FASTA: {pep_file}")
    lines.append(f"HGNC rows: {hgnc.shape[0]}")
    lines.append(f"Protein long rows: {protein_long.shape[0]}")
    lines.append(f"Genes with protein rows: {protein_long['gene_symbol'].nunique() if not protein_long.empty else 0}")
    lines.append(f"Gene feature rows: {features.shape[0]}")
    lines.append(f"Gene feature columns: {features.shape[1]}")
    lines.append(f"Merged rows: {merged.shape[0]}")
    lines.append(f"Merged columns: {merged.shape[1]}")
    lines.append("")
    lines.append("Coverage:")

    if "feature11_has_protein_sequence" in features.columns:
        lines.append(f"Protein sequence coverage: {features['feature11_has_protein_sequence'].mean():.4f}")

    if "feature11_protein_isoform_count" in features.columns:
        lines.append(f"Median protein isoform count: {features['feature11_protein_isoform_count'].median():.2f}")
        lines.append(f"Max protein isoform count: {features['feature11_protein_isoform_count'].max():.2f}")

    if "feature11_summary_protein_length_mean" in features.columns:
        lines.append(f"Median mean protein length: {features['feature11_summary_protein_length_mean'].median(skipna=True):.2f}")

    if "feature11_summary_hydropathy_mean_mean" in features.columns:
        lines.append(f"Median mean hydropathy: {features['feature11_summary_hydropathy_mean_mean'].median(skipna=True):.4f}")

    if "feature11_summary_low_complexity_fraction_w15_mean" in features.columns:
        lines.append(
            "Median low-complexity fraction: "
            f"{features['feature11_summary_low_complexity_fraction_w15_mean'].median(skipna=True):.4f}"
        )

    lines.append("")
    lines.append("Feature interpretation:")
    lines.append("Protein length, amino-acid composition, charge, hydrophobicity, entropy, low-complexity, repeat and TM-like segment proxies.")
    lines.append("")
    lines.append("Leakage policy:")
    lines.append("Included: amino-acid sequence-derived features only.")
    lines.append("Excluded: drug labels, known target labels, clinical labels, ChEMBL, DrugBank, DGIdb, Open Targets, Pharos.")

    path.write_text("\n".join(lines) + "\n")

    log(f"[SUMMARY] {path}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description="Build Feature 11 protein sequence features.")
    parser.add_argument("--hgnc", default=str(DEFAULT_HGNC_PATH), help="Path to HGNC complete set.")
    parser.add_argument("--dbdir", default=str(DEFAULT_DBDIR), help="Ensembl database directory.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory.")
    parser.add_argument("--pep-file", default="", help="Manual Ensembl peptide FASTA file.")
    parser.add_argument("--all-hgnc-genes", action="store_true", help="Use all HGNC genes, not only protein-coding.")
    parser.add_argument("--limit-genes", type=int, default=0, help="Debug only: first N HGNC genes.")
    parser.add_argument("--max-records", type=int, default=0, help="Debug only: parse first N protein FASTA records.")
    args = parser.parse_args()

    outdir = mkdir(Path(args.outdir))
    processed_dir = mkdir(outdir / "processed")
    dbdir = Path(args.dbdir)

    log("=" * 100)
    log("FEATURE 11: PROTEIN SEQUENCE COMPOSITION / PHYSICOCHEMICAL FEATURES")
    log("=" * 100)
    log(f"[HGNC]        {args.hgnc}")
    log(f"[DBDIR]       {dbdir.resolve()}")
    log(f"[OUTDIR]      {outdir.resolve()}")
    log(f"[PEP FILE]    {args.pep_file if args.pep_file else 'auto-detect'}")
    log(f"[LIMIT GENES] {args.limit_genes if args.limit_genes else 'none'}")
    log(f"[MAX RECORDS] {args.max_records if args.max_records else 'none'}")
    log("=" * 100)

    pep_file = Path(args.pep_file) if args.pep_file else find_pep_fasta(dbdir)

    if pep_file is None or not pep_file.exists():
        raise FileNotFoundError(
            "Could not auto-detect Ensembl peptide FASTA. Provide it using --pep-file."
        )

    log(f"[USING PEP FILE] {pep_file}")

    hgnc = load_hgnc(
        Path(args.hgnc),
        protein_coding_only=not args.all_hgnc_genes,
        limit_genes=args.limit_genes,
    )

    ensembl_to_gene = build_ensembl_to_gene(hgnc)

    protein_long = build_protein_long_table(
        pep_file=pep_file,
        ensembl_to_gene=ensembl_to_gene,
        processed_dir=processed_dir,
        max_records=args.max_records,
    )

    features = aggregate_gene_features(
        hgnc=hgnc,
        protein_long=protein_long,
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
        "pep_file": str(pep_file),
        "protein_coding_only": not args.all_hgnc_genes,
        "limit_genes": args.limit_genes,
        "max_records": args.max_records,
        "outputs": {
            "protein_long": str(processed_dir / "feature11_protein_sequence_long.csv"),
            "gene_features": str(processed_dir / "feature11_protein_sequence_gene_features.csv"),
            "hgnc_merged": str(processed_dir / "feature11_protein_sequence_hgnc_merged.csv"),
        },
        "leakage_policy": {
            "included": [
                "amino-acid sequence length",
                "amino-acid composition",
                "physicochemical composition",
                "hydrophobicity",
                "charge",
                "entropy",
                "low-complexity proxy",
                "TM-like segment proxy",
                "molecular-weight proxy",
                "isoelectric-point proxy",
            ],
            "excluded": [
                "known drug labels",
                "clinical target labels",
                "ChEMBL",
                "DrugBank",
                "DGIdb",
                "Open Targets",
                "Pharos",
            ],
        },
    }

    with open(outdir / "feature11_protein_sequence_run_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    write_summary(
        outdir=outdir,
        hgnc=hgnc,
        protein_long=protein_long,
        features=features,
        merged=merged,
        pep_file=pep_file,
        args=args,
    )

    log("=" * 100)
    log("[DONE]")
    log(f"[PROTEIN LONG]  {processed_dir / 'feature11_protein_sequence_long.csv'}")
    log(f"[GENE FEATURES] {processed_dir / 'feature11_protein_sequence_gene_features.csv'}")
    log(f"[HGNC MERGED]   {processed_dir / 'feature11_protein_sequence_hgnc_merged.csv'}")
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