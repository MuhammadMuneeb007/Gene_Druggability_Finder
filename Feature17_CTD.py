#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature17_CTD.py

FEATURE 17: CTD chemical-gene interaction burden features

This script downloads and processes the CTD chemical-gene interaction file:
    CTD_chem_gene_ixns.tsv.gz

Important fix:
    CTD stores the real header line with a leading '#':
        # ChemicalName ChemicalID CasRN GeneSymbol GeneID ...

    If pandas uses comment="#", it skips the real header and treats the first
    data row as column names. This script manually detects the CTD header,
    strips the '#', and then reads the table correctly.

Outputs:
    feature17_ctd/processed/feature17_ctd_long_human_hgnc_mapped.csv.gz
    feature17_ctd/processed/feature17_ctd_gene_features.csv
    feature17_ctd/processed/feature17_ctd_gene_features_hgnc_merged.csv
    feature17_ctd/processed/feature17_ctd_summary.txt

Example:
    python Feature17_CTD.py --download

Optional:
    python Feature17_CTD.py --download --force-download
    python Feature17_CTD.py --max-rows 100000
    python Feature17_CTD.py --limit-genes 500
    python Feature17_CTD.py --include-nonhuman
"""

import argparse
import gzip
import os
import re
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd


# =============================================================================
# Constants
# =============================================================================

SCRIPT_NAME = "Feature17_CTD.py"

CTD_URLS = [
    "http://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz",
    "https://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz",
]

DEFAULT_HGNC = "databases/HGNC/hgnc_complete_set.txt"
DEFAULT_DBDIR = "feature_databases/CTD"
DEFAULT_OUTDIR = "feature17_ctd"

CTD_FILENAME = "CTD_chem_gene_ixns.tsv.gz"

HUMAN_TAX_ID = "9606"

ACTION_KEYWORDS = {
    "affects": "affects",
    "increases": "increases",
    "decreases": "decreases",
    "binds": "binds",
    "cotreatment": "cotreatment",
    "response": "response",
    "expression": "expression",
    "activity": "activity",
    "reaction": "reaction",
    "localization": "localization",
    "transport": "transport",
    "secretion": "secretion",
    "abundance": "abundance",
    "phosphorylation": "phosphorylation",
    "methylation": "methylation",
    "acetylation": "acetylation",
    "ubiquitination": "ubiquitination",
    "mutagenesis": "mutagenesis",
    "cleavage": "cleavage",
    "metabolic processing": "metabolic_processing",
    "stability": "stability",
    "folding": "folding",
}


# =============================================================================
# Pretty printing
# =============================================================================

def line(char="=", n=100):
    print(char * n, flush=True)


def log(msg):
    print(msg, flush=True)


def die(msg):
    line()
    log("[ERROR]")
    log(str(msg))
    line()
    raise RuntimeError(str(msg))


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)
    return Path(path)


def file_size_mb(path):
    path = Path(path)
    if not path.exists():
        return 0.0
    return path.stat().st_size / (1024 * 1024)


# =============================================================================
# Download
# =============================================================================

def download_file(urls, dest, force=False, timeout=120):
    dest = Path(dest)
    ensure_dir(dest.parent)

    if dest.exists() and dest.stat().st_size > 0 and not force:
        log(f"[DOWNLOAD SKIP] {dest} already exists ({file_size_mb(dest):.2f} MB)")
        return dest

    last_error = None

    for url in urls:
        try:
            line()
            log(f"[DOWNLOAD ATTEMPT] {url}")
            log(f"[TO]               {dest}")

            req = urllib.request.Request(
                url,
                headers={
                    "User-Agent": "Mozilla/5.0 Feature17_CTD.py"
                },
            )

            with urllib.request.urlopen(req, timeout=timeout) as response:
                with open(dest, "wb") as out:
                    block_size = 1024 * 1024
                    total = 0
                    while True:
                        chunk = response.read(block_size)
                        if not chunk:
                            break
                        out.write(chunk)
                        total += len(chunk)

            if dest.exists() and dest.stat().st_size > 0:
                log(f"[DOWNLOAD DONE] {dest} ({file_size_mb(dest):.2f} MB)")
                return dest

        except Exception as e:
            last_error = e
            log(f"[DOWNLOAD FAILED] {url}")
            log(f"[REASON]          {repr(e)}")

    die(f"Could not download CTD file. Last error: {last_error}")


# =============================================================================
# CTD reader
# =============================================================================

def find_ctd_header(path):
    """
    CTD files contain metadata/comment lines before the table.

    The actual header is usually:
        # ChemicalName ChemicalID CasRN GeneSymbol GeneID GeneForms ...

    We must not use pandas comment='#' with header inference because that skips
    the real header and makes the first data row into the dataframe columns.
    """
    path = Path(path)
    opener = gzip.open if str(path).endswith(".gz") else open

    header = None
    skiprows = 0

    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        for line_raw in f:
            skiprows += 1
            line_clean = line_raw.strip()

            if not line_clean:
                continue

            line_no_hash = line_clean.lstrip("#").strip()
            fields = line_no_hash.split("\t")

            if len(fields) >= 8:
                lower_fields = [x.lower() for x in fields]

                has_chemical = "chemicalname" in lower_fields or "chemicalid" in lower_fields
                has_gene = "genesymbol" in lower_fields or "geneid" in lower_fields
                has_organism = "organism" in lower_fields or "organismid" in lower_fields

                if has_chemical and has_gene and has_organism:
                    header = fields
                    break

    if header is None:
        die(f"Could not find CTD header line in {path}")

    header = [str(x).strip().lstrip("#").strip() for x in header]

    return header, skiprows


def read_ctd_chem_gene_file(path, max_rows=None):
    """
    Read CTD_chem_gene_ixns.tsv.gz correctly.
    """
    path = Path(path)

    header, skiprows = find_ctd_header(path)

    line()
    log("[CTD HEADER DETECTED]")
    log(f"[SKIPROWS] {skiprows}")
    log(f"[HEADER]   {header}")

    df = pd.read_csv(
        path,
        sep="\t",
        names=header,
        skiprows=skiprows,
        nrows=max_rows,
        dtype=str,
        low_memory=False,
        compression="gzip" if str(path).endswith(".gz") else None,
        keep_default_na=False,
    )

    df.columns = [str(c).strip().lstrip("#").strip() for c in df.columns]

    return df


# =============================================================================
# HGNC utilities
# =============================================================================

def split_hgnc_multi_value(x):
    if pd.isna(x):
        return []
    x = str(x).strip()
    if x == "" or x.lower() == "nan":
        return []

    # HGNC fields are often pipe-separated.
    parts = re.split(r"\|", x)
    out = []
    for p in parts:
        p = str(p).strip()
        if p and p.lower() != "nan":
            out.append(p)
    return out


def clean_entrez_id(x):
    if pd.isna(x):
        return None
    x = str(x).strip()
    if x == "" or x.lower() == "nan":
        return None

    # Sometimes read as "367.0"
    if re.match(r"^\d+\.0$", x):
        x = x[:-2]

    if not re.match(r"^\d+$", x):
        return None

    return x


def read_hgnc(hgnc_path, protein_coding_only=True):
    hgnc_path = Path(hgnc_path)

    if not hgnc_path.exists():
        die(f"HGNC file not found: {hgnc_path}")

    line()
    log(f"[READ HGNC] {hgnc_path}")

    hgnc = pd.read_csv(
        hgnc_path,
        sep="\t",
        dtype=str,
        low_memory=False,
        keep_default_na=False,
    )

    n0 = len(hgnc)

    required = ["symbol"]
    for col in required:
        if col not in hgnc.columns:
            die(f"HGNC file missing required column: {col}")

    if protein_coding_only:
        if "locus_group" in hgnc.columns:
            hgnc = hgnc[hgnc["locus_group"].astype(str).str.lower().eq("protein-coding gene")].copy()
        elif "locus_type" in hgnc.columns:
            hgnc = hgnc[hgnc["locus_type"].astype(str).str.lower().str.contains("protein-coding", na=False)].copy()
        else:
            log("[WARNING] Could not find locus_group/locus_type; not applying protein-coding filter.")

    n1 = len(hgnc)

    log(f"[HGNC ROWS] {n0} -> {n1}")
    log(f"[PROTEIN CODING ONLY] {protein_coding_only}")

    hgnc["approved_symbol"] = hgnc["symbol"].astype(str).str.strip()

    if "entrez_id" not in hgnc.columns:
        hgnc["entrez_id"] = ""

    if "hgnc_id" not in hgnc.columns:
        hgnc["hgnc_id"] = ""

    if "name" not in hgnc.columns:
        hgnc["name"] = ""

    if "alias_symbol" not in hgnc.columns:
        hgnc["alias_symbol"] = ""

    if "prev_symbol" not in hgnc.columns:
        hgnc["prev_symbol"] = ""

    if "uniprot_ids" not in hgnc.columns:
        hgnc["uniprot_ids"] = ""

    # Symbol and alias map.
    symbol_to_approved = {}

    for _, row in hgnc.iterrows():
        approved = str(row["approved_symbol"]).strip()

        if not approved:
            continue

        candidates = [approved]
        candidates.extend(split_hgnc_multi_value(row.get("alias_symbol", "")))
        candidates.extend(split_hgnc_multi_value(row.get("prev_symbol", "")))

        for s in candidates:
            s = str(s).strip()
            if not s:
                continue
            symbol_to_approved[s.upper()] = approved

    # Entrez map.
    entrez_to_approved = {}

    for _, row in hgnc.iterrows():
        approved = str(row["approved_symbol"]).strip()
        entrez = clean_entrez_id(row.get("entrez_id", ""))

        if approved and entrez:
            entrez_to_approved[entrez] = approved

    line()
    log("[HGNC MAPS]")
    log(f"[SYMBOLS + ALIASES] {len(symbol_to_approved)}")
    log(f"[ENTREZ IDS]         {len(entrez_to_approved)}")

    keep_cols = [
        "hgnc_id",
        "approved_symbol",
        "name",
        "entrez_id",
        "uniprot_ids",
        "locus_group",
        "locus_type",
    ]

    keep_cols = [c for c in keep_cols if c in hgnc.columns]

    hgnc_keep = hgnc[keep_cols].drop_duplicates("approved_symbol").copy()

    return hgnc_keep, symbol_to_approved, entrez_to_approved


def resolve_hgnc_symbol(gene_symbol, gene_id, symbol_to_approved, entrez_to_approved):
    gene_id_clean = clean_entrez_id(gene_id)

    if gene_id_clean and gene_id_clean in entrez_to_approved:
        return entrez_to_approved[gene_id_clean]

    if pd.notna(gene_symbol):
        s = str(gene_symbol).strip()
        if s:
            s_upper = s.upper()
            if s_upper in symbol_to_approved:
                return symbol_to_approved[s_upper]

    return None


# =============================================================================
# CTD column detection
# =============================================================================

def normalize_colname(c):
    return re.sub(r"[^a-z0-9]+", "", str(c).strip().lower())


def detect_column(columns, candidates):
    norm_to_original = {normalize_colname(c): c for c in columns}

    for cand in candidates:
        key = normalize_colname(cand)
        if key in norm_to_original:
            return norm_to_original[key]

    return None


def detect_ctd_columns(ctd):
    cols = list(ctd.columns)

    detected = {
        "chemical_name": detect_column(cols, ["ChemicalName", "chemical_name", "Chemical Name"]),
        "chemical_id": detect_column(cols, ["ChemicalID", "chemical_id", "Chemical ID"]),
        "cas_rn": detect_column(cols, ["CasRN", "CAS", "CASRN", "Cas RN"]),
        "gene_symbol": detect_column(cols, ["GeneSymbol", "gene_symbol", "Gene Symbol"]),
        "gene_id": detect_column(cols, ["GeneID", "gene_id", "Gene ID"]),
        "gene_forms": detect_column(cols, ["GeneForms", "gene_forms", "Gene Forms"]),
        "organism": detect_column(cols, ["Organism", "organism"]),
        "organism_id": detect_column(cols, ["OrganismID", "organism_id", "Organism ID"]),
        "interaction": detect_column(cols, ["Interaction", "interaction"]),
        "interaction_actions": detect_column(cols, ["InteractionActions", "interaction_actions", "Interaction Actions"]),
        "pubmed_ids": detect_column(cols, ["PubMedIDs", "PubMedID", "pubmed_ids", "PubMed IDs"]),
    }

    line()
    log("[DETECTED CTD COLUMNS]")
    for k, v in detected.items():
        log(f"{k:<24} = {v}")

    if detected["gene_symbol"] is None and detected["gene_id"] is None:
        die("Could not detect GeneSymbol or GeneID columns in CTD file.")

    return detected


# =============================================================================
# Long table builder
# =============================================================================

def count_pubmed_ids(x):
    if pd.isna(x):
        return 0

    x = str(x).strip()

    if x == "" or x.lower() == "nan":
        return 0

    # CTD PubMed IDs are commonly pipe-separated.
    parts = [p.strip() for p in re.split(r"[|;,]", x) if p.strip()]
    return len(set(parts))


def first_pubmed_id(x):
    if pd.isna(x):
        return ""

    x = str(x).strip()

    if x == "" or x.lower() == "nan":
        return ""

    parts = [p.strip() for p in re.split(r"[|;,]", x) if p.strip()]
    return parts[0] if parts else ""


def build_ctd_long_table(
    ctd,
    detected,
    symbol_to_approved,
    entrez_to_approved,
    human_only=True,
    limit_genes=None,
):
    line()
    log("[BUILD CTD LONG TABLE]")

    def getcol(key):
        col = detected.get(key)
        if col is None:
            return pd.Series([""] * len(ctd), index=ctd.index, dtype=str)
        return ctd[col].astype(str).fillna("")

    long_df = pd.DataFrame({
        "chemical_name": getcol("chemical_name"),
        "chemical_id": getcol("chemical_id"),
        "cas_rn": getcol("cas_rn"),
        "ctd_gene_symbol": getcol("gene_symbol"),
        "ctd_gene_id": getcol("gene_id"),
        "gene_forms": getcol("gene_forms"),
        "organism": getcol("organism"),
        "organism_id": getcol("organism_id"),
        "interaction": getcol("interaction"),
        "interaction_actions": getcol("interaction_actions"),
        "pubmed_ids": getcol("pubmed_ids"),
    })

    n0 = len(long_df)

    if human_only:
        if "organism_id" in long_df.columns:
            long_df = long_df[long_df["organism_id"].astype(str).str.strip().eq(HUMAN_TAX_ID)].copy()
        elif "organism" in long_df.columns:
            long_df = long_df[long_df["organism"].astype(str).str.lower().eq("homo sapiens")].copy()

    n_human = len(long_df)

    log(f"[ROWS BEFORE HUMAN FILTER] {n0}")
    log(f"[ROWS AFTER HUMAN FILTER]  {n_human}")

    # Resolve to approved HGNC symbol.
    approved_symbols = []

    for gs, gid in zip(long_df["ctd_gene_symbol"].values, long_df["ctd_gene_id"].values):
        approved_symbols.append(resolve_hgnc_symbol(gs, gid, symbol_to_approved, entrez_to_approved))

    long_df["approved_symbol"] = approved_symbols

    n_before_map = len(long_df)
    long_df = long_df[long_df["approved_symbol"].notna()].copy()
    n_after_map = len(long_df)

    log(f"[ROWS BEFORE HGNC MAP] {n_before_map}")
    log(f"[ROWS AFTER HGNC MAP]  {n_after_map}")
    log(f"[UNIQUE HGNC GENES]    {long_df['approved_symbol'].nunique()}")

    if limit_genes is not None and int(limit_genes) > 0:
        keep_genes = sorted(long_df["approved_symbol"].dropna().unique())[: int(limit_genes)]
        long_df = long_df[long_df["approved_symbol"].isin(keep_genes)].copy()
        log(f"[LIMIT GENES]          {limit_genes}")
        log(f"[ROWS AFTER LIMIT]     {len(long_df)}")
        log(f"[GENES AFTER LIMIT]    {long_df['approved_symbol'].nunique()}")

    long_df["n_pubmed_ids_row"] = long_df["pubmed_ids"].apply(count_pubmed_ids).astype(np.int32)
    long_df["first_pubmed_id"] = long_df["pubmed_ids"].apply(first_pubmed_id)

    long_df["interaction_actions_lower"] = long_df["interaction_actions"].astype(str).str.lower()
    long_df["interaction_lower"] = long_df["interaction"].astype(str).str.lower()

    # Add row-level binary action flags.
    for raw_kw, safe_kw in ACTION_KEYWORDS.items():
        pattern = re.escape(raw_kw.lower())
        colname = f"ctd_row_has_action_{safe_kw}"
        long_df[colname] = (
            long_df["interaction_actions_lower"].str.contains(pattern, na=False)
            | long_df["interaction_lower"].str.contains(pattern, na=False)
        ).astype(np.int8)

    return long_df


# =============================================================================
# Feature construction
# =============================================================================

def safe_nunique(series):
    return series.replace("", np.nan).dropna().nunique()


def build_gene_features(long_df):
    line()
    log("[BUILD GENE-LEVEL CTD FEATURES]")

    if long_df.empty:
        die("Long CTD table is empty after filtering/mapping.")

    grouped = long_df.groupby("approved_symbol", sort=True)

    features = pd.DataFrame(index=sorted(long_df["approved_symbol"].unique()))
    features.index.name = "approved_symbol"

    features["ctd_chem_gene_interaction_rows"] = grouped.size().astype(np.int64)

    features["ctd_unique_chemicals"] = grouped["chemical_id"].agg(safe_nunique).astype(np.int64)
    features["ctd_unique_chemical_names"] = grouped["chemical_name"].agg(safe_nunique).astype(np.int64)
    features["ctd_unique_cas_rn"] = grouped["cas_rn"].agg(safe_nunique).astype(np.int64)

    features["ctd_unique_interaction_text"] = grouped["interaction"].agg(safe_nunique).astype(np.int64)
    features["ctd_unique_interaction_actions"] = grouped["interaction_actions"].agg(safe_nunique).astype(np.int64)

    features["ctd_total_pubmed_mentions"] = grouped["n_pubmed_ids_row"].sum().astype(np.int64)
    features["ctd_unique_pubmed_ids_approx"] = grouped["first_pubmed_id"].agg(safe_nunique).astype(np.int64)

    # GeneForms can include protein, mRNA, etc.
    if "gene_forms" in long_df.columns:
        gf = long_df[["approved_symbol", "gene_forms"]].copy()
        gf["gene_forms_lower"] = gf["gene_forms"].astype(str).str.lower()

        for form in ["protein", "mrna", "rna", "gene", "promoter"]:
            tmp = gf.assign(flag=gf["gene_forms_lower"].str.contains(form, na=False).astype(np.int8))
            features[f"ctd_rows_geneform_{form}"] = tmp.groupby("approved_symbol")["flag"].sum().reindex(features.index).fillna(0).astype(np.int64)

    # Action keyword counts.
    for raw_kw, safe_kw in ACTION_KEYWORDS.items():
        row_col = f"ctd_row_has_action_{safe_kw}"
        if row_col in long_df.columns:
            features[f"ctd_rows_action_{safe_kw}"] = grouped[row_col].sum().astype(np.int64)

    # Directional summary.
    inc_col = "ctd_rows_action_increases"
    dec_col = "ctd_rows_action_decreases"

    if inc_col in features.columns and dec_col in features.columns:
        features["ctd_directional_rows_increase_minus_decrease"] = (
            features[inc_col] - features[dec_col]
        ).astype(np.int64)

        denom = features[inc_col] + features[dec_col]
        features["ctd_directional_increase_fraction"] = np.where(
            denom > 0,
            features[inc_col] / denom,
            0.0,
        )

        features["ctd_directional_decrease_fraction"] = np.where(
            denom > 0,
            features[dec_col] / denom,
            0.0,
        )

    # Transformations.
    count_cols = [
        c for c in features.columns
        if c.startswith("ctd_")
        and not c.endswith("_fraction")
        and features[c].dtype.kind in "iuf"
    ]

    for c in count_cols:
        features[f"{c}_log1p"] = np.log1p(features[c].astype(float))

    features = features.reset_index()

    log(f"[GENE FEATURES SHAPE] {features.shape}")

    return features


def merge_with_hgnc_universe(hgnc_keep, features):
    line()
    log("[MERGE WITH HGNC UNIVERSE]")

    merged = hgnc_keep.merge(
        features,
        on="approved_symbol",
        how="left",
    )

    feature_cols = [c for c in merged.columns if c.startswith("ctd_")]

    for c in feature_cols:
        merged[c] = merged[c].fillna(0)

    merged["ctd_has_any_chemical_gene_interaction"] = (
        merged["ctd_chem_gene_interaction_rows"].astype(float) > 0
    ).astype(int)

    log(f"[HGNC MERGED SHAPE] {merged.shape}")
    log(f"[GENES WITH CTD INTERACTIONS] {int(merged['ctd_has_any_chemical_gene_interaction'].sum())}")
    log(f"[TOTAL HGNC GENES]            {len(merged)}")

    return merged


# =============================================================================
# Summary
# =============================================================================

def write_summary(path, args, ctd_path, raw_shape, long_df, features, merged):
    path = Path(path)

    with open(path, "w") as f:
        f.write("FEATURE 17: CTD CHEMICAL-GENE INTERACTION BURDEN\n")
        f.write("=" * 80 + "\n")
        f.write(f"Script: {SCRIPT_NAME}\n")
        f.write(f"CTD file: {ctd_path}\n")
        f.write(f"CTD file size MB: {file_size_mb(ctd_path):.2f}\n")
        f.write(f"Raw CTD shape: {raw_shape}\n")
        f.write(f"Human only: {not args.include_nonhuman}\n")
        f.write(f"Max rows: {args.max_rows}\n")
        f.write(f"Limit genes: {args.limit_genes}\n")
        f.write("\n")
        f.write(f"Long mapped rows: {len(long_df)}\n")
        f.write(f"Long mapped genes: {long_df['approved_symbol'].nunique() if not long_df.empty else 0}\n")
        f.write(f"Gene feature shape: {features.shape}\n")
        f.write(f"HGNC merged shape: {merged.shape}\n")
        f.write("\n")

        if "ctd_has_any_chemical_gene_interaction" in merged.columns:
            f.write(f"Genes with any CTD chemical-gene interaction: {int(merged['ctd_has_any_chemical_gene_interaction'].sum())}\n")

        if "ctd_chem_gene_interaction_rows" in merged.columns:
            s = merged["ctd_chem_gene_interaction_rows"].astype(float)
            f.write("\nInteraction row burden summary over HGNC universe:\n")
            f.write(str(s.describe()) + "\n")


# =============================================================================
# Main
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Feature 17: CTD chemical-gene interaction burden features"
    )

    parser.add_argument(
        "--hgnc",
        default=DEFAULT_HGNC,
        help=f"Path to HGNC complete set TSV. Default: {DEFAULT_HGNC}",
    )

    parser.add_argument(
        "--dbdir",
        default=DEFAULT_DBDIR,
        help=f"Directory for CTD database file. Default: {DEFAULT_DBDIR}",
    )

    parser.add_argument(
        "--outdir",
        default=DEFAULT_OUTDIR,
        help=f"Output directory. Default: {DEFAULT_OUTDIR}",
    )

    parser.add_argument(
        "--ctd-file",
        default="auto",
        help="Path to CTD_chem_gene_ixns.tsv.gz, or 'auto'. Default: auto",
    )

    parser.add_argument(
        "--download",
        action="store_true",
        help="Download CTD chemical-gene interaction file if missing.",
    )

    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Force re-download CTD file.",
    )

    parser.add_argument(
        "--include-nonhuman",
        action="store_true",
        help="Include non-human CTD rows. Default is human only.",
    )

    parser.add_argument(
        "--max-rows",
        type=int,
        default=None,
        help="Maximum CTD rows to read. Useful for testing.",
    )

    parser.add_argument(
        "--limit-genes",
        type=int,
        default=None,
        help="Limit to first N mapped HGNC genes. Useful for testing.",
    )

    parser.add_argument(
        "--no-protein-coding-filter",
        action="store_true",
        help="Do not restrict HGNC universe to protein-coding genes.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    t0 = time.time()

    hgnc_path = Path(args.hgnc)
    dbdir = ensure_dir(args.dbdir)
    outdir = ensure_dir(args.outdir)
    processed_dir = ensure_dir(outdir / "processed")

    if args.ctd_file == "auto":
        ctd_path = dbdir / CTD_FILENAME
    else:
        ctd_path = Path(args.ctd_file)

    line()
    log("FEATURE 17: CTD CHEMICAL-GENE INTERACTION BURDEN")
    line()
    log(f"[HGNC]           {hgnc_path}")
    log(f"[DBDIR]          {dbdir.resolve()}")
    log(f"[OUTDIR]         {outdir.resolve()}")
    log(f"[CTD FILE]       {args.ctd_file}")
    log(f"[DOWNLOAD]       {args.download}")
    log(f"[FORCE DOWNLOAD] {args.force_download}")
    log(f"[HUMAN ONLY]     {not args.include_nonhuman}")
    log(f"[MAX ROWS]       {args.max_rows if args.max_rows is not None else 'none'}")
    log(f"[LIMIT GENES]    {args.limit_genes if args.limit_genes is not None else 'none'}")
    line()

    if args.download:
        ctd_path = download_file(
            CTD_URLS,
            ctd_path,
            force=args.force_download,
        )

    if not ctd_path.exists():
        die(
            f"CTD file not found: {ctd_path}\n"
            f"Run with --download or provide --ctd-file /path/to/{CTD_FILENAME}"
        )

    hgnc_keep, symbol_to_approved, entrez_to_approved = read_hgnc(
        hgnc_path,
        protein_coding_only=not args.no_protein_coding_filter,
    )

    line()
    log(f"[READ CTD CHEM-GENE] {ctd_path}")

    ctd_raw = read_ctd_chem_gene_file(
        ctd_path,
        max_rows=args.max_rows,
    )

    log(f"[CTD SHAPE] {ctd_raw.shape}")
    log(f"[CTD COLUMNS] {list(ctd_raw.columns)}")

    detected = detect_ctd_columns(ctd_raw)

    long_df = build_ctd_long_table(
        ctd=ctd_raw,
        detected=detected,
        symbol_to_approved=symbol_to_approved,
        entrez_to_approved=entrez_to_approved,
        human_only=not args.include_nonhuman,
        limit_genes=args.limit_genes,
    )

    features = build_gene_features(long_df)

    merged = merge_with_hgnc_universe(
        hgnc_keep=hgnc_keep,
        features=features,
    )

    long_out = processed_dir / "feature17_ctd_long_human_hgnc_mapped.csv.gz"
    features_out = processed_dir / "feature17_ctd_gene_features.csv"
    merged_out = processed_dir / "feature17_ctd_gene_features_hgnc_merged.csv"
    summary_out = processed_dir / "feature17_ctd_summary.txt"

    line()
    log("[WRITE OUTPUTS]")

    long_df.to_csv(long_out, index=False, compression="gzip")
    log(f"[WRITE] {long_out}")

    features.to_csv(features_out, index=False)
    log(f"[WRITE] {features_out}")

    merged.to_csv(merged_out, index=False)
    log(f"[WRITE] {merged_out}")

    write_summary(
        summary_out,
        args=args,
        ctd_path=ctd_path,
        raw_shape=ctd_raw.shape,
        long_df=long_df,
        features=features,
        merged=merged,
    )
    log(f"[WRITE] {summary_out}")

    elapsed = time.time() - t0

    line()
    log("[DONE]")
    log(f"[RUNTIME SEC] {elapsed:.2f}")
    line()


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        line()
        log("[ERROR]")
        log(str(e))
        line()
        raise