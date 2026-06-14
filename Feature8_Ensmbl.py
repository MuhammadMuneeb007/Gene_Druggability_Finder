#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature8_Ensembl.py

Feature 8: Ensembl genomic annotation features.

This script uses the files you have:

    feature_databases/Ensembl/Homo_sapiens.GRCh38.115.gtf.gz
    feature_databases/Ensembl/Homo_sapiens.GRCh38.cdna.all.fa.gz
    feature_databases/Ensembl/Homo_sapiens.GRCh38.pep.all.fa.gz

Main features:
    gene span
    chromosome flags
    strand
    transcript count
    exon count
    CDS count
    UTR count
    exon/CDS/UTR union length
    exon/CDS/UTR density
    transcript sequence length summaries from cDNA FASTA
    protein sequence length summaries from peptide FASTA

Outputs:
    feature8_ensembl/processed/feature8_ensembl_long_features.csv
    feature8_ensembl/processed/feature8_ensembl_gene_features.csv
    feature8_ensembl/processed/feature8_ensembl_hgnc_merged.csv
    feature8_ensembl/feature8_ensembl_summary.txt
    feature8_ensembl/feature8_ensembl_run_metadata.json

Run:
    python Feature8_Ensembl.py

Fast test:
    python Feature8_Ensembl.py --max-lines 500000

Manual file:
    python Feature8_Ensembl.py \
      --annotation-file feature_databases/Ensembl/Homo_sapiens.GRCh38.115.gtf.gz
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


DEFAULT_HGNC_PATH = Path("databases") / "HGNC" / "hgnc_complete_set.txt"
DEFAULT_DBDIR = Path("feature_databases") / "Ensembl"
DEFAULT_OUTDIR = Path("feature8_ensembl")

CANONICAL_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT", "M"}

GTF_FEATURES_TO_KEEP = {
    "gene",
    "transcript",
    "exon",
    "CDS",
    "five_prime_utr",
    "three_prime_utr",
    "UTR",
    "start_codon",
    "stop_codon",
}


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
    if s.lower() in {"nan", "none"}:
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


def safe_int(x: Any) -> int:
    try:
        if x is None or pd.isna(x):
            return 0
        return int(float(x))
    except Exception:
        return 0


def safe_float(x: Any) -> float:
    try:
        if x is None or pd.isna(x):
            return np.nan
        return float(x)
    except Exception:
        return np.nan


def open_text_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", errors="ignore")
    return open(path, "rt", errors="ignore")


def normalize_chromosome(chrom: Any) -> str:
    c = clean_text(chrom)
    c = c.replace("chr", "").replace("CHR", "")
    if c == "M":
        c = "MT"
    return c


def is_canonical_chromosome(chrom: Any) -> int:
    return int(normalize_chromosome(chrom) in CANONICAL_CHROMS)


def interval_union_length(intervals: List[Tuple[int, int]]) -> int:
    cleaned = []

    for a, b in intervals:
        a = safe_int(a)
        b = safe_int(b)

        if a <= 0 or b <= 0:
            continue

        if b < a:
            a, b = b, a

        cleaned.append((a, b))

    if not cleaned:
        return 0

    cleaned.sort()

    total = 0
    cur_a, cur_b = cleaned[0]

    for a, b in cleaned[1:]:
        if a <= cur_b + 1:
            cur_b = max(cur_b, b)
        else:
            total += cur_b - cur_a + 1
            cur_a, cur_b = a, b

    total += cur_b - cur_a + 1

    return int(total)


def parse_gtf_attributes(attr: str) -> Dict[str, str]:
    """
    Parse GTF/GFF attributes.

    GTF example:
        gene_id "ENSG..."; gene_name "TP53"; gene_biotype "protein_coding";

    GFF3 example:
        ID=gene:ENSG...;Name=TP53;biotype=protein_coding;
    """
    out = {}

    attr = clean_text(attr)

    if not attr:
        return out

    # GTF key "value";
    for m in re.finditer(r'([A-Za-z0-9_.:-]+)\s+"([^"]*)"', attr):
        out[m.group(1)] = m.group(2)

    # GFF3 key=value;
    for part in attr.split(";"):
        part = part.strip()

        if not part or "=" not in part:
            continue

        k, v = part.split("=", 1)
        k = k.strip()
        v = v.strip()

        if k and v and k not in out:
            out[k] = v

    return out


def find_file_by_patterns(dbdir: Path, patterns: List[str]) -> Optional[Path]:
    files = []

    for pattern in patterns:
        files.extend(list(dbdir.rglob(pattern)))

    files = [p for p in files if p.exists() and p.stat().st_size > 0]

    if not files:
        return None

    files = sorted(files, key=lambda p: p.stat().st_size, reverse=True)

    return files[0]


def find_ensembl_gtf(dbdir: Path) -> Optional[Path]:
    candidates = []
    for pattern in ["*.gtf.gz", "*.gtf", "*.gff3.gz", "*.gff3", "*.gff.gz", "*.gff"]:
        candidates.extend(list(dbdir.rglob(pattern)))

    if not candidates:
        return None

    ranked = []

    for p in candidates:
        name = p.name.lower()
        score = 0

        if "homo_sapiens" in name:
            score += 20
        if "grch38" in name:
            score += 15
        if name.endswith(".gtf.gz") or name.endswith(".gtf"):
            score += 20
        if "115" in name:
            score += 5
        if "abinitio" in name:
            score -= 20

        ranked.append((score, p.stat().st_size, p))

    ranked.sort(reverse=True)

    return ranked[0][2]


def find_cdna_fasta(dbdir: Path) -> Optional[Path]:
    return find_file_by_patterns(
        dbdir,
        [
            "*cdna*.fa.gz",
            "*cdna*.fasta.gz",
            "*cdna*.fa",
            "*cdna*.fasta",
        ],
    )


def find_pep_fasta(dbdir: Path) -> Optional[Path]:
    return find_file_by_patterns(
        dbdir,
        [
            "*pep*.fa.gz",
            "*pep*.fasta.gz",
            "*pep*.fa",
            "*pep*.fasta",
        ],
    )


def load_hgnc(hgnc_path: Path, protein_coding_only: bool = True) -> pd.DataFrame:
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
    if "location" not in hgnc.columns:
        hgnc["location"] = ""

    hgnc["ensembl_gene_id_clean"] = hgnc["ensembl_gene_id"].map(clean_ensembl_gene_id)
    hgnc = hgnc.drop_duplicates(subset=["gene_symbol"], keep="first").copy()

    log(f"[HGNC ROWS] {before} -> {len(hgnc)}")
    log(f"[PROTEIN CODING ONLY] {protein_coding_only}")
    log(f"[GENES WITH ENSEMBL] {(hgnc['ensembl_gene_id_clean'].astype(str).str.len() > 0).sum()}")

    return hgnc


def build_hgnc_maps(hgnc: pd.DataFrame) -> Tuple[Dict[str, str], Dict[str, str]]:
    ens_to_gene = {}
    symbol_to_gene = {}

    for _, row in hgnc.iterrows():
        gene = row["gene_symbol"]
        ens = clean_ensembl_gene_id(row.get("ensembl_gene_id_clean", ""))

        if ens:
            ens_to_gene[ens] = gene

        if gene:
            symbol_to_gene[gene] = gene

    return ens_to_gene, symbol_to_gene


def parse_gtf_or_gff(
    annotation_file: Path,
    ens_to_gene: Dict[str, str],
    symbol_to_gene: Dict[str, str],
    processed_dir: Path,
    max_lines: int = 0,
) -> pd.DataFrame:
    log("=" * 100)
    log(f"[STREAM ENSEMBL GTF/GFF] {annotation_file}")

    rows = []
    total = 0
    matched = 0
    t0 = time.time()

    with open_text_maybe_gzip(annotation_file) as f:
        for line in f:
            total += 1

            if max_lines and total > max_lines:
                break

            if total % 5_000_000 == 0:
                elapsed = (time.time() - t0) / 60
                log(f"[GTF] scanned={total:,} matched={matched:,} elapsed={elapsed:.1f} min")

            if not line or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 9:
                continue

            chrom, source, feature_type, start, end, score, strand, frame, attrs = parts[:9]

            if feature_type not in GTF_FEATURES_TO_KEEP:
                continue

            attr = parse_gtf_attributes(attrs)

            gene_id = clean_ensembl_gene_id(
                attr.get("gene_id")
                or attr.get("gene")
                or attr.get("ID")
                or ""
            )

            gene_name = normalize_symbol(
                attr.get("gene_name")
                or attr.get("Name")
                or attr.get("gene_symbol")
                or ""
            )

            gene_symbol = ""

            if gene_id in ens_to_gene:
                gene_symbol = ens_to_gene[gene_id]
            elif gene_name in symbol_to_gene:
                gene_symbol = gene_name

            if not gene_symbol:
                continue

            transcript_id = clean_ensembl_transcript_id(
                attr.get("transcript_id")
                or attr.get("transcript")
                or attr.get("Parent")
                or ""
            )

            exon_id = clean_text(
                attr.get("exon_id")
                or attr.get("exon_number")
                or ""
            )

            biotype = clean_text(
                attr.get("gene_biotype")
                or attr.get("gene_type")
                or attr.get("biotype")
                or attr.get("transcript_biotype")
                or attr.get("transcript_type")
                or ""
            )

            s = safe_int(start)
            e = safe_int(end)

            if s <= 0 or e <= 0:
                continue

            length = abs(e - s) + 1

            rows.append(
                {
                    "gene_symbol": gene_symbol,
                    "ensembl_gene_id": gene_id,
                    "gene_name_from_ensembl": gene_name,
                    "chromosome": normalize_chromosome(chrom),
                    "source": source,
                    "feature_type": feature_type,
                    "start": s,
                    "end": e,
                    "length": length,
                    "strand": strand,
                    "frame": frame,
                    "transcript_id": transcript_id,
                    "exon_id": exon_id,
                    "biotype": biotype,
                    "is_canonical_chromosome": is_canonical_chromosome(chrom),
                }
            )

            matched += 1

    long_df = pd.DataFrame(rows)

    if not long_df.empty:
        long_df = long_df.drop_duplicates()

    outpath = processed_dir / "feature8_ensembl_long_features.csv"
    long_df.to_csv(outpath, index=False)

    log(f"[GTF DONE] scanned={total:,} matched_rows={len(long_df):,}")
    log(f"[SAVED LONG] {outpath}")

    return long_df


def parse_fasta_gene_lengths(
    fasta_file: Optional[Path],
    ens_to_gene: Dict[str, str],
    prefix: str,
    max_records: int = 0,
) -> pd.DataFrame:
    """
    Parse Ensembl cDNA or peptide FASTA and summarize sequence lengths per HGNC gene.

    Ensembl FASTA headers usually contain:
        gene:ENSG00000141510.18
        transcript:ENST...
    """
    if fasta_file is None or not fasta_file.exists():
        log(f"[SKIP FASTA] No {prefix} FASTA found.")
        return pd.DataFrame(columns=["gene_symbol"])

    log("=" * 100)
    log(f"[STREAM FASTA] {fasta_file}")
    log(f"[PREFIX] {prefix}")

    gene_lengths: Dict[str, List[int]] = defaultdict(list)

    current_gene = ""
    current_len = 0
    records_seen = 0

    def flush_record():
        nonlocal current_gene, current_len, records_seen

        if current_gene and current_len > 0:
            gene_lengths[current_gene].append(current_len)

        current_gene = ""
        current_len = 0
        records_seen += 1

    with open_text_maybe_gzip(fasta_file) as f:
        for line in f:
            line = line.rstrip("\n")

            if line.startswith(">"):
                if current_len > 0 or current_gene:
                    flush_record()

                    if max_records and records_seen >= max_records:
                        break

                header = line[1:]

                gene_id = ""
                m = re.search(r"gene:(ENSG[0-9]+(?:\.[0-9]+)?)", header)
                if m:
                    gene_id = clean_ensembl_gene_id(m.group(1))

                gene = ens_to_gene.get(gene_id, "")

                current_gene = gene
                current_len = 0

            else:
                if current_gene:
                    current_len += len(line.strip())

        if current_len > 0 or current_gene:
            flush_record()

    rows = []

    for gene, lengths in gene_lengths.items():
        arr = np.asarray(lengths, dtype=float)

        rows.append(
            {
                "gene_symbol": gene,
                f"feature8_ensembl_{prefix}_sequence_count": int(arr.size),
                f"feature8_ensembl_{prefix}_mean_length": float(np.mean(arr)) if arr.size else np.nan,
                f"feature8_ensembl_{prefix}_median_length": float(np.median(arr)) if arr.size else np.nan,
                f"feature8_ensembl_{prefix}_min_length": float(np.min(arr)) if arr.size else np.nan,
                f"feature8_ensembl_{prefix}_max_length": float(np.max(arr)) if arr.size else np.nan,
                f"feature8_ensembl_{prefix}_std_length": float(np.std(arr, ddof=1)) if arr.size > 1 else 0.0,
            }
        )

    out = pd.DataFrame(rows)

    log(f"[FASTA DONE] {prefix} genes={out.shape[0]} records={records_seen:,}")

    return out


def build_gene_features(
    hgnc: pd.DataFrame,
    long_df: pd.DataFrame,
    cdna_features: pd.DataFrame,
    pep_features: pd.DataFrame,
    processed_dir: Path,
) -> pd.DataFrame:
    genes = pd.DataFrame({"gene_symbol": sorted(hgnc["gene_symbol"].unique().tolist())})
    records = []

    by_gene = dict(tuple(long_df.groupby("gene_symbol"))) if not long_df.empty else {}

    autosomes = {str(i) for i in range(1, 23)}

    for gene in genes["gene_symbol"]:
        g = by_gene.get(gene, pd.DataFrame())

        if g.empty:
            records.append({"gene_symbol": gene, "feature8_ensembl_has_annotation": 0})
            continue

        g = g.copy()

        gene_rows = g[g["feature_type"] == "gene"].copy()
        transcript_rows = g[g["feature_type"] == "transcript"].copy()
        exon_rows = g[g["feature_type"] == "exon"].copy()
        cds_rows = g[g["feature_type"] == "CDS"].copy()
        utr_rows = g[g["feature_type"].isin(["UTR", "five_prime_utr", "three_prime_utr"])].copy()
        start_rows = g[g["feature_type"] == "start_codon"].copy()
        stop_rows = g[g["feature_type"] == "stop_codon"].copy()

        representative = gene_rows.iloc[0].to_dict() if len(gene_rows) else g.iloc[0].to_dict()

        gene_start = safe_int(g["start"].min())
        gene_end = safe_int(g["end"].max())
        genomic_span = abs(gene_end - gene_start) + 1 if gene_start > 0 and gene_end > 0 else np.nan

        exon_intervals = [(safe_int(r["start"]), safe_int(r["end"])) for _, r in exon_rows.iterrows()]
        cds_intervals = [(safe_int(r["start"]), safe_int(r["end"])) for _, r in cds_rows.iterrows()]
        utr_intervals = [(safe_int(r["start"]), safe_int(r["end"])) for _, r in utr_rows.iterrows()]

        exon_union = interval_union_length(exon_intervals)
        cds_union = interval_union_length(cds_intervals)
        utr_union = interval_union_length(utr_intervals)

        intron_proxy = genomic_span - exon_union if not pd.isna(genomic_span) else np.nan
        if not pd.isna(intron_proxy):
            intron_proxy = max(0, intron_proxy)

        exon_lengths = pd.to_numeric(exon_rows["length"], errors="coerce").dropna().tolist() if len(exon_rows) else []
        cds_lengths = pd.to_numeric(cds_rows["length"], errors="coerce").dropna().tolist() if len(cds_rows) else []

        transcript_ids = sorted(set([x for x in g["transcript_id"].dropna().astype(str) if x]))
        biotypes = sorted(set([x for x in g["biotype"].dropna().astype(str) if x]))

        strand = clean_text(representative.get("strand", ""))
        chrom = clean_text(representative.get("chromosome", ""))

        rec = {
            "gene_symbol": gene,
            "feature8_ensembl_has_annotation": 1,

            "feature8_ensembl_gene_id": clean_text(representative.get("ensembl_gene_id", "")),
            "feature8_ensembl_chromosome": chrom,

            "feature8_ensembl_is_canonical_chromosome": is_canonical_chromosome(chrom),
            "feature8_ensembl_is_autosomal": int(chrom in autosomes),
            "feature8_ensembl_is_chrX": int(chrom == "X"),
            "feature8_ensembl_is_chrY": int(chrom == "Y"),
            "feature8_ensembl_is_mitochondrial": int(chrom in {"MT", "M"}),

            "feature8_ensembl_gene_start": gene_start,
            "feature8_ensembl_gene_end": gene_end,
            "feature8_ensembl_genomic_span_bp": genomic_span,
            "feature8_ensembl_log1p_genomic_span_bp": float(np.log1p(genomic_span)) if not pd.isna(genomic_span) else np.nan,

            "feature8_ensembl_strand_plus": int(strand == "+"),
            "feature8_ensembl_strand_minus": int(strand == "-"),

            "feature8_ensembl_feature_row_count": int(len(g)),
            "feature8_ensembl_gene_row_count": int(len(gene_rows)),
            "feature8_ensembl_transcript_count": int(len(transcript_ids)),
            "feature8_ensembl_exon_row_count": int(len(exon_rows)),
            "feature8_ensembl_cds_row_count": int(len(cds_rows)),
            "feature8_ensembl_utr_row_count": int(len(utr_rows)),
            "feature8_ensembl_start_codon_count": int(len(start_rows)),
            "feature8_ensembl_stop_codon_count": int(len(stop_rows)),

            "feature8_ensembl_exon_union_length_bp": int(exon_union),
            "feature8_ensembl_cds_union_length_bp": int(cds_union),
            "feature8_ensembl_utr_union_length_bp": int(utr_union),
            "feature8_ensembl_intron_length_proxy_bp": intron_proxy,

            "feature8_ensembl_exon_density": float(exon_union / genomic_span) if genomic_span and genomic_span > 0 else np.nan,
            "feature8_ensembl_cds_density": float(cds_union / genomic_span) if genomic_span and genomic_span > 0 else np.nan,
            "feature8_ensembl_utr_density": float(utr_union / genomic_span) if genomic_span and genomic_span > 0 else np.nan,
            "feature8_ensembl_intron_density_proxy": float(intron_proxy / genomic_span) if genomic_span and genomic_span > 0 and not pd.isna(intron_proxy) else np.nan,

            "feature8_ensembl_mean_exon_length": float(np.mean(exon_lengths)) if exon_lengths else np.nan,
            "feature8_ensembl_median_exon_length": float(np.median(exon_lengths)) if exon_lengths else np.nan,
            "feature8_ensembl_max_exon_length": float(np.max(exon_lengths)) if exon_lengths else np.nan,
            "feature8_ensembl_min_exon_length": float(np.min(exon_lengths)) if exon_lengths else np.nan,

            "feature8_ensembl_mean_cds_length": float(np.mean(cds_lengths)) if cds_lengths else np.nan,
            "feature8_ensembl_median_cds_length": float(np.median(cds_lengths)) if cds_lengths else np.nan,
            "feature8_ensembl_max_cds_length": float(np.max(cds_lengths)) if cds_lengths else np.nan,
            "feature8_ensembl_min_cds_length": float(np.min(cds_lengths)) if cds_lengths else np.nan,

            "feature8_ensembl_biotype_count": int(len(biotypes)),
            "feature8_ensembl_is_protein_coding_biotype": int(any("protein_coding" in b for b in biotypes)),
            "feature8_ensembl_has_multiple_transcripts": int(len(transcript_ids) > 1),
        }

        records.append(rec)

    features = pd.DataFrame(records)

    if not cdna_features.empty:
        features = features.merge(cdna_features, on="gene_symbol", how="left")

    if not pep_features.empty:
        features = features.merge(pep_features, on="gene_symbol", how="left")

    for c in features.columns:
        if c.startswith("feature8_") and (
            "_count" in c
            or "_has_" in c
            or "_is_" in c
            or c.endswith("_plus")
            or c.endswith("_minus")
        ):
            if features[c].dtype != object:
                features[c] = features[c].fillna(0)

    outpath = processed_dir / "feature8_ensembl_gene_features.csv"
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

    merged["feature8_ensembl_has_any_feature"] = merged[feature_cols].notna().any(axis=1).astype(int)
    merged["feature8_ensembl_n_nonmissing_features"] = merged[feature_cols].notna().sum(axis=1)

    outpath = processed_dir / "feature8_ensembl_hgnc_merged.csv"
    merged.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED HGNC MERGED]")
    log(f"[PATH]  {outpath}")
    log(f"[SHAPE] {merged.shape}")

    return merged


def write_summary(
    outdir: Path,
    hgnc: pd.DataFrame,
    long_df: pd.DataFrame,
    features: pd.DataFrame,
    merged: pd.DataFrame,
    annotation_file: Path,
    cdna_file: Optional[Path],
    pep_file: Optional[Path],
    args: argparse.Namespace,
) -> None:
    path = outdir / "feature8_ensembl_summary.txt"

    lines = []
    lines.append("Feature 8: Ensembl genomic annotation features")
    lines.append("=" * 80)
    lines.append(f"Created at: {now_iso()}")
    lines.append(f"Annotation file: {annotation_file}")
    lines.append(f"cDNA FASTA: {cdna_file if cdna_file else 'not found / skipped'}")
    lines.append(f"Peptide FASTA: {pep_file if pep_file else 'not found / skipped'}")
    lines.append(f"HGNC rows: {hgnc.shape[0]}")
    lines.append(f"Long annotation rows: {long_df.shape[0]}")
    lines.append(f"Genes with annotation rows: {long_df['gene_symbol'].nunique() if not long_df.empty else 0}")
    lines.append(f"Gene feature rows: {features.shape[0]}")
    lines.append(f"Gene feature columns: {features.shape[1]}")
    lines.append(f"Merged rows: {merged.shape[0]}")
    lines.append(f"Merged columns: {merged.shape[1]}")
    lines.append("")
    lines.append("Coverage:")

    if "feature8_ensembl_has_annotation" in features.columns:
        lines.append(f"Ensembl annotation coverage: {features['feature8_ensembl_has_annotation'].mean():.4f}")

    if "feature8_ensembl_transcript_count" in features.columns:
        lines.append(f"Median transcript count: {features['feature8_ensembl_transcript_count'].median():.2f}")

    if "feature8_ensembl_exon_row_count" in features.columns:
        lines.append(f"Median exon row count: {features['feature8_ensembl_exon_row_count'].median():.2f}")

    if "feature8_ensembl_genomic_span_bp" in features.columns:
        lines.append(f"Median genomic span bp: {features['feature8_ensembl_genomic_span_bp'].median():.2f}")

    if "feature8_ensembl_pep_sequence_count" in features.columns:
        lines.append(f"Genes with peptide sequence lengths: {(features['feature8_ensembl_pep_sequence_count'].fillna(0) > 0).mean():.4f}")

    if "feature8_ensembl_cdna_sequence_count" in features.columns:
        lines.append(f"Genes with cDNA sequence lengths: {(features['feature8_ensembl_cdna_sequence_count'].fillna(0) > 0).mean():.4f}")

    lines.append("")
    lines.append("Leakage policy:")
    lines.append("Included: genomic coordinates, transcript/exon/CDS/UTR structure, strand, chromosome, biotype, cDNA and peptide sequence lengths.")
    lines.append("Excluded: drug-target labels, clinical labels, disease associations, tractability labels and known-drug evidence.")

    path.write_text("\n".join(lines) + "\n")
    log(f"[SUMMARY] {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build Ensembl genomic annotation features for HGNC genes.")
    parser.add_argument("--hgnc", default=str(DEFAULT_HGNC_PATH), help="Path to HGNC complete set.")
    parser.add_argument("--dbdir", default=str(DEFAULT_DBDIR), help="Local Ensembl database directory.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory.")
    parser.add_argument("--annotation-file", default="", help="Direct path to Ensembl GTF/GFF.")
    parser.add_argument("--cdna-file", default="", help="Direct path to Ensembl cDNA FASTA.")
    parser.add_argument("--pep-file", default="", help="Direct path to Ensembl peptide FASTA.")
    parser.add_argument("--all-hgnc-genes", action="store_true", help="Use all HGNC rows, not only protein-coding.")
    parser.add_argument("--max-lines", type=int, default=0, help="Debug only: parse first N GTF lines.")
    parser.add_argument("--skip-fasta", action="store_true", help="Skip cDNA/peptide FASTA length features.")
    args = parser.parse_args()

    outdir = mkdir(Path(args.outdir))
    processed_dir = mkdir(outdir / "processed")
    dbdir = Path(args.dbdir)

    log("=" * 100)
    log("FEATURE 8: ENSEMBL GENOMIC ANNOTATION")
    log("=" * 100)
    log(f"[HGNC]            {args.hgnc}")
    log(f"[DBDIR]           {dbdir.resolve()}")
    log(f"[OUTDIR]          {outdir.resolve()}")
    log(f"[ANNOTATION FILE] {args.annotation_file if args.annotation_file else 'auto-detect'}")
    log(f"[CDNA FILE]       {args.cdna_file if args.cdna_file else 'auto-detect'}")
    log(f"[PEP FILE]        {args.pep_file if args.pep_file else 'auto-detect'}")
    log(f"[MAX LINES]       {args.max_lines if args.max_lines else 'none'}")
    log(f"[SKIP FASTA]      {args.skip_fasta}")
    log("=" * 100)

    annotation_file = Path(args.annotation_file) if args.annotation_file else find_ensembl_gtf(dbdir)

    if annotation_file is None or not annotation_file.exists():
        raise FileNotFoundError(
            "Could not auto-detect Ensembl GTF/GFF file. Provide one with --annotation-file."
        )

    cdna_file = None
    pep_file = None

    if not args.skip_fasta:
        cdna_file = Path(args.cdna_file) if args.cdna_file else find_cdna_fasta(dbdir)
        pep_file = Path(args.pep_file) if args.pep_file else find_pep_fasta(dbdir)

    log(f"[USING ANNOTATION FILE] {annotation_file}")
    log(f"[USING CDNA FILE]       {cdna_file if cdna_file else 'none'}")
    log(f"[USING PEP FILE]        {pep_file if pep_file else 'none'}")

    hgnc = load_hgnc(
        Path(args.hgnc),
        protein_coding_only=not args.all_hgnc_genes,
    )

    ens_to_gene, symbol_to_gene = build_hgnc_maps(hgnc)

    long_df = parse_gtf_or_gff(
        annotation_file=annotation_file,
        ens_to_gene=ens_to_gene,
        symbol_to_gene=symbol_to_gene,
        processed_dir=processed_dir,
        max_lines=args.max_lines,
    )

    if args.skip_fasta:
        cdna_features = pd.DataFrame(columns=["gene_symbol"])
        pep_features = pd.DataFrame(columns=["gene_symbol"])
    else:
        cdna_features = parse_fasta_gene_lengths(
            fasta_file=cdna_file,
            ens_to_gene=ens_to_gene,
            prefix="cdna",
        )
        pep_features = parse_fasta_gene_lengths(
            fasta_file=pep_file,
            ens_to_gene=ens_to_gene,
            prefix="pep",
        )

    features = build_gene_features(
        hgnc=hgnc,
        long_df=long_df,
        cdna_features=cdna_features,
        pep_features=pep_features,
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
        "annotation_file": str(annotation_file),
        "cdna_file": str(cdna_file) if cdna_file else "",
        "pep_file": str(pep_file) if pep_file else "",
        "protein_coding_only": not args.all_hgnc_genes,
        "max_lines": args.max_lines,
        "skip_fasta": args.skip_fasta,
        "outputs": {
            "long_features": str(processed_dir / "feature8_ensembl_long_features.csv"),
            "gene_features": str(processed_dir / "feature8_ensembl_gene_features.csv"),
            "hgnc_merged": str(processed_dir / "feature8_ensembl_hgnc_merged.csv"),
        },
        "leakage_policy": {
            "included": [
                "gene coordinates",
                "gene span",
                "transcript count",
                "exon count",
                "CDS/UTR structure",
                "strand",
                "chromosome",
                "biotype",
                "cDNA sequence length",
                "peptide sequence length",
            ],
            "excluded": [
                "drug-target labels",
                "clinical labels",
                "known drug evidence",
                "disease association evidence",
                "tractability labels",
            ],
        },
    }

    with open(outdir / "feature8_ensembl_run_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    write_summary(
        outdir=outdir,
        hgnc=hgnc,
        long_df=long_df,
        features=features,
        merged=merged,
        annotation_file=annotation_file,
        cdna_file=cdna_file,
        pep_file=pep_file,
        args=args,
    )

    log("=" * 100)
    log("[DONE]")
    log(f"[LONG FEATURES]  {processed_dir / 'feature8_ensembl_long_features.csv'}")
    log(f"[GENE FEATURES]  {processed_dir / 'feature8_ensembl_gene_features.csv'}")
    log(f"[HGNC MERGED]    {processed_dir / 'feature8_ensembl_hgnc_merged.csv'}")
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