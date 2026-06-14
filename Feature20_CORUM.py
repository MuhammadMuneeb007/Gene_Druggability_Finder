#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature20_CORUM.py

FEATURE 20: CORUM multi-file protein complex features

Purpose
-------
Build HGNC gene-level protein-complex features from local CORUM files.

This version is local-file first and robust to the CORUM filenames you already have:

    feature_databases/CORUM/corum_humanComplexes.txt
    feature_databases/CORUM/corum_allComplexes.txt
    feature_databases/CORUM/corum_fcg.txt
    feature_databases/CORUM/corum_uniprotCorumMapping.txt

Optional files are supported if present:

    drug-target complexes
    drug-complex interactions
    splice-variant complexes
    partial complexes

Main output
-----------
feature20_corum/processed/feature20_corum_gene_features_hgnc_merged.csv

Run
---
python Feature20_CORUM.py

Quick test
----------
python Feature20_CORUM.py --limit-rows 1000
"""

import argparse
import gzip
import io
import os
import re
import sys
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


# =============================================================================
# Defaults
# =============================================================================

SCRIPT_NAME = "Feature20_CORUM.py"

DEFAULT_HGNC = "databases/HGNC/hgnc_complete_set.txt"
DEFAULT_DBDIR = "feature_databases/CORUM"
DEFAULT_OUTDIR = "feature20_corum"

NULL_STRINGS = {"", "none", "null", "nan", "na", "-", "n/a"}

HUMAN_PATTERNS = [
    "human",
    "homo sapiens",
]

DATASET_PATTERNS = {
    "human_complexes": [
        "corum_humanComplexes.txt",
        "humanComplexes.txt",
        "humanComplexes.txt.zip",
        "*human*complex*.txt",
        "*human*complex*.zip",
    ],
    "complete_complexes": [
        "corum_allComplexes.txt",
        "allComplexes.txt",
        "allComplexes.txt.zip",
        "*all*complex*.txt",
        "*complete*complex*.txt",
        "*all*complex*.zip",
        "*complete*complex*.zip",
    ],
    "drug_target_complexes": [
        "corum_drugTargetComplexes.txt",
        "drugTargetComplexes.txt",
        "drugTargetComplexes.txt.zip",
        "*drug*target*complex*.txt",
        "*drug*target*complex*.zip",
        "*complex*drug*target*.txt",
        "*complex*drug*target*.zip",
    ],
    "drug_complex_interactions": [
        "corum_drugComplexInteractions.txt",
        "drugComplexInteractions.txt",
        "drugComplexInteractions.txt.zip",
        "*drug*complex*interaction*.txt",
        "*drug*complex*interaction*.zip",
        "*formal*description*.txt",
        "*formal*description*.zip",
    ],
    "splice_variant_complexes": [
        "corum_spliceVariantComplexes.txt",
        "spliceVariantComplexes.txt",
        "spliceVariantComplexes.txt.zip",
        "*splice*variant*complex*.txt",
        "*splice*variant*complex*.zip",
    ],
    "partial_complexes": [
        "corum_partialComplexes.txt",
        "partialComplexes.txt",
        "partialComplexes.txt.zip",
        "*partial*complex*.txt",
        "*partial*complex*.zip",
    ],
    "functional_complex_groups": [
        "corum_fcg.txt",
        "functionalComplexGroups.txt",
        "functionalComplexGroups.txt.zip",
        "*functional*complex*group*.txt",
        "*functional*complex*group*.zip",
        "*fcg*.txt",
        "*fcg*.zip",
    ],
    "uniprot_corum_mapping": [
        "corum_uniprotCorumMapping.txt",
        "uniprotCorumMapping.txt",
        "uniprotCorumMapping.txt.zip",
        "*uniprot*corum*mapping*.txt",
        "*uniprot*corum*mapping*.zip",
        "*uniprot*corum*.txt",
        "*uniprot*corum*.zip",
    ],
}


# =============================================================================
# Logging
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
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def file_size_mb(path):
    path = Path(path)
    if not path.exists():
        return 0.0
    return path.stat().st_size / (1024 * 1024)


# =============================================================================
# Generic helpers
# =============================================================================

def clean_text(x):
    if pd.isna(x):
        return ""
    x = str(x).strip()
    if x.lower() in NULL_STRINGS:
        return ""
    return x


def clean_entrez_id(x):
    x = clean_text(x)
    if not x:
        return None

    if re.match(r"^\d+\.0$", x):
        x = x[:-2]

    if re.match(r"^\d+$", x):
        return x

    return None


def split_hgnc_multi_value(x):
    x = clean_text(x)
    if not x:
        return []
    return [p.strip() for p in x.split("|") if p.strip()]


def split_multi_value(x):
    x = clean_text(x)
    if not x:
        return []

    parts = re.split(r"\s*[;,|]\s*", x)
    out = []

    for p in parts:
        p = clean_text(p)
        if p:
            out.append(p)

    return out


def normalize_colname(c):
    return re.sub(r"[^a-z0-9]+", "", str(c).strip().lower())


def detect_column(columns, candidates=None, contains_all=None, contains_any=None):
    columns = list(columns)
    norm_to_original = {normalize_colname(c): c for c in columns}

    if candidates:
        for cand in candidates:
            key = normalize_colname(cand)
            if key in norm_to_original:
                return norm_to_original[key]

    if contains_all:
        tokens = [normalize_colname(t) for t in contains_all]
        for c in columns:
            cn = normalize_colname(c)
            if all(t in cn for t in tokens):
                return c

    if contains_any:
        tokens = [normalize_colname(t) for t in contains_any]
        for c in columns:
            cn = normalize_colname(c)
            if any(t in cn for t in tokens):
                return c

    return None


def safe_nunique(series):
    if series is None:
        return 0
    return series.replace("", np.nan).dropna().nunique()


def is_human_organism(x):
    x = clean_text(x).lower()
    if not x:
        return False
    return any(p in x for p in HUMAN_PATTERNS)


def first_nonempty(values):
    for v in values:
        v = clean_text(v)
        if v:
            return v
    return ""


# =============================================================================
# Find local CORUM files
# =============================================================================

def find_dataset_file(dataset_key, dbdir):
    dbdir = Path(dbdir)

    patterns = DATASET_PATTERNS.get(dataset_key, [])

    for pattern in patterns:
        p = dbdir / pattern
        if "*" not in pattern and p.exists() and p.stat().st_size > 0:
            return p

    for pattern in patterns:
        if "*" in pattern:
            hits = sorted(dbdir.glob(pattern))
            hits = [h for h in hits if h.is_file() and h.stat().st_size > 0]
            if hits:
                return hits[0]

    return None


# =============================================================================
# Read HGNC
# =============================================================================

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
        keep_default_na=False,
        low_memory=False,
    )

    n0 = len(hgnc)

    if "symbol" not in hgnc.columns:
        die("HGNC file missing required column: symbol")

    if protein_coding_only:
        if "locus_group" in hgnc.columns:
            hgnc = hgnc[
                hgnc["locus_group"].astype(str).str.lower().eq("protein-coding gene")
            ].copy()
        elif "locus_type" in hgnc.columns:
            hgnc = hgnc[
                hgnc["locus_type"].astype(str).str.lower().str.contains("protein-coding", na=False)
            ].copy()
        else:
            log("[WARNING] HGNC locus_group/locus_type not found; protein-coding filter not applied.")

    n1 = len(hgnc)

    log(f"[HGNC ROWS] {n0} -> {n1}")
    log(f"[PROTEIN CODING ONLY] {protein_coding_only}")

    hgnc["approved_symbol"] = hgnc["symbol"].astype(str).str.strip()

    for col in [
        "hgnc_id",
        "name",
        "entrez_id",
        "alias_symbol",
        "prev_symbol",
        "uniprot_ids",
        "locus_group",
        "locus_type",
    ]:
        if col not in hgnc.columns:
            hgnc[col] = ""

    symbol_to_approved = {}
    entrez_to_approved = {}
    uniprot_to_approved = {}

    for _, row in hgnc.iterrows():
        approved = clean_text(row.get("approved_symbol", ""))

        if not approved:
            continue

        symbols = [approved]
        symbols.extend(split_hgnc_multi_value(row.get("alias_symbol", "")))
        symbols.extend(split_hgnc_multi_value(row.get("prev_symbol", "")))

        for s in symbols:
            s = clean_text(s)
            if s:
                symbol_to_approved[s.upper()] = approved

        entrez = clean_entrez_id(row.get("entrez_id", ""))
        if entrez:
            entrez_to_approved[entrez] = approved

        for u in split_hgnc_multi_value(row.get("uniprot_ids", "")):
            u = clean_text(u).upper()
            if u:
                uniprot_to_approved[u] = approved

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

    line()
    log("[HGNC MAPS]")
    log(f"[SYMBOLS + ALIASES] {len(symbol_to_approved)}")
    log(f"[ENTREZ IDS]         {len(entrez_to_approved)}")
    log(f"[UNIPROT IDS]        {len(uniprot_to_approved)}")

    return hgnc_keep, symbol_to_approved, entrez_to_approved, uniprot_to_approved


def resolve_gene(symbol, entrez, uniprot, symbol_to_approved, entrez_to_approved, uniprot_to_approved):
    entrez = clean_entrez_id(entrez)

    if entrez and entrez in entrez_to_approved:
        return entrez_to_approved[entrez], "entrez"

    uniprot = clean_text(uniprot).upper()

    # Remove isoform suffix if exact isoform not in HGNC map.
    uniprot_base = uniprot.split("-")[0] if uniprot else ""

    if uniprot and uniprot in uniprot_to_approved:
        return uniprot_to_approved[uniprot], "uniprot"

    if uniprot_base and uniprot_base in uniprot_to_approved:
        return uniprot_to_approved[uniprot_base], "uniprot_base"

    symbol = clean_text(symbol)

    if symbol and symbol.upper() in symbol_to_approved:
        return symbol_to_approved[symbol.upper()], "symbol"

    return None, "unmapped"


# =============================================================================
# Read CORUM table
# =============================================================================

def read_corum_table(path, dataset_key, limit_rows=None):
    path = Path(path)

    if not path.exists():
        die(f"CORUM file not found: {path}")

    line()
    log(f"[READ CORUM TABLE] {dataset_key}")
    log(f"[PATH] {path}")

    if str(path).lower().endswith(".zip"):
        with zipfile.ZipFile(path, "r") as z:
            names = z.namelist()
            members = [
                n for n in names
                if n.lower().endswith(".txt")
                or n.lower().endswith(".tsv")
                or n.lower().endswith(".csv")
            ]

            if not members:
                die(f"No text file found inside ZIP: {path}. Members: {names}")

            member = members[0]
            log(f"[ZIP MEMBER] {member}")

            with z.open(member) as f:
                raw = f.read()

            text = raw.decode("utf-8", errors="replace")

            df = pd.read_csv(
                io.StringIO(text),
                sep="\t",
                dtype=str,
                keep_default_na=False,
                low_memory=False,
                nrows=limit_rows,
            )

    elif str(path).lower().endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
            df = pd.read_csv(
                f,
                sep="\t",
                dtype=str,
                keep_default_na=False,
                low_memory=False,
                nrows=limit_rows,
            )

    else:
        df = pd.read_csv(
            path,
            sep="\t",
            dtype=str,
            keep_default_na=False,
            low_memory=False,
            nrows=limit_rows,
            encoding="utf-8",
            encoding_errors="replace",
        )

    df.columns = [str(c).strip() for c in df.columns]

    bad_cols = [c for c in df.columns if c == "" or c.lower().startswith("unnamed")]
    if bad_cols:
        df = df.drop(columns=bad_cols)

    log(f"[SHAPE]   {df.shape}")
    log(f"[COLUMNS] {list(df.columns)}")

    return df


# =============================================================================
# CORUM column detection
# =============================================================================

def detect_complex_columns(df):
    cols = list(df.columns)

    detected = {
        "complex_id": detect_column(
            cols,
            candidates=[
                "ComplexID",
                "Complex ID",
                "complex_id",
                "CORUM ID",
                "CORUM-ID",
                "Complex id",
            ],
            contains_all=["complex", "id"],
        ),
        "complex_name": detect_column(
            cols,
            candidates=[
                "ComplexName",
                "Complex name",
                "complex_name",
                "Name",
                "Complex",
            ],
            contains_all=["complex", "name"],
        ),
        "organism": detect_column(
            cols,
            candidates=["Organism", "Species", "organism"],
        ),
        "subunits_gene_symbol": detect_column(
            cols,
            candidates=[
                "subunits(Gene name)",
                "Subunits(Gene name)",
                "subunits gene name",
                "subunits_gene_name",
                "Gene name",
                "Gene names",
                "Genes",
                "Subunits",
            ],
            contains_all=["subunits", "gene"],
        ),
        "subunits_uniprot": detect_column(
            cols,
            candidates=[
                "subunits(UniProt IDs)",
                "Subunits(UniProt IDs)",
                "subunits UniProt IDs",
                "UniProt IDs",
                "Uniprot",
                "UniProt",
            ],
            contains_any=["uniprot"],
        ),
        "subunits_entrez": detect_column(
            cols,
            candidates=[
                "subunits(Entrez IDs)",
                "Subunits(Entrez IDs)",
                "subunits Entrez IDs",
                "Entrez IDs",
                "Entrez",
            ],
            contains_any=["entrez"],
        ),
        "pubmed_id": detect_column(
            cols,
            candidates=["PubMed ID", "PubMedID", "PMID", "PubMed"],
            contains_any=["pubmed"],
        ),
        "cell_line": detect_column(
            cols,
            candidates=["Cell line", "CellLine", "cell_line"],
            contains_all=["cell", "line"],
        ),
        "purification_method": detect_column(
            cols,
            candidates=["Purification method", "purification_method"],
            contains_any=["purification"],
        ),
        "go_id": detect_column(
            cols,
            candidates=["GO ID", "GOID", "GO term ID", "Gene Ontology ID", "GO"],
            contains_any=["go"],
        ),
        "function": detect_column(
            cols,
            candidates=["Function", "Functional comment", "Complex function", "Description"],
            contains_any=["function"],
        ),
    }

    return detected


def print_detected(title, detected):
    line()
    log(title)
    for k, v in detected.items():
        log(f"{k:<24} = {v}")


# =============================================================================
# UniProt-CORUM mapping helper
# =============================================================================

def build_corum_uniprot_to_complex_map(mapping_df):
    """
    Optional helper. CORUM mapping files can vary by release.
    We use this only to create extra long rows if a mapping file clearly has
    UniProt and complex ID columns.
    """
    if mapping_df is None or mapping_df.empty:
        return pd.DataFrame()

    cols = list(mapping_df.columns)

    uniprot_col = detect_column(
        cols,
        candidates=["UniProt", "UniProt ID", "Uniprot", "uniprot"],
        contains_any=["uniprot"],
    )
    complex_col = detect_column(
        cols,
        candidates=["ComplexID", "Complex ID", "CORUM ID", "CORUM-ID"],
        contains_all=["complex", "id"],
    )
    complex_name_col = detect_column(
        cols,
        candidates=["ComplexName", "Complex name"],
        contains_all=["complex", "name"],
    )

    if uniprot_col is None or complex_col is None:
        log("[UNIPROT-CORUM MAP] Could not detect UniProt/ComplexID columns. Skipping mapping expansion.")
        return pd.DataFrame()

    out = pd.DataFrame({
        "complex_id": mapping_df[complex_col].astype(str),
        "raw_uniprot_id": mapping_df[uniprot_col].astype(str),
        "complex_name": mapping_df[complex_name_col].astype(str) if complex_name_col else "",
    })

    out["raw_uniprot_id"] = out["raw_uniprot_id"].map(clean_text)
    out["complex_id"] = out["complex_id"].map(clean_text)
    out["complex_name"] = out["complex_name"].map(clean_text)

    out = out[(out["complex_id"] != "") & (out["raw_uniprot_id"] != "")].copy()
    out = out.drop_duplicates()

    log(f"[UNIPROT-CORUM MAP] usable rows: {len(out)}")

    return out


# =============================================================================
# Build long membership table
# =============================================================================

def build_complex_membership_long(
    df,
    dataset_key,
    symbol_to_approved,
    entrez_to_approved,
    uniprot_to_approved,
    human_only=True,
):
    if df is None or df.empty:
        return pd.DataFrame()

    detected = detect_complex_columns(df)
    print_detected(f"[DETECTED COLUMNS: {dataset_key}]", detected)

    has_subunit_info = (
        detected["subunits_gene_symbol"] is not None
        or detected["subunits_uniprot"] is not None
        or detected["subunits_entrez"] is not None
    )

    if not has_subunit_info:
        log(f"[WARNING] No subunit columns detected for {dataset_key}. Skipping membership expansion.")
        return pd.DataFrame()

    work = df.copy()

    if human_only and detected["organism"] is not None:
        n0 = len(work)
        work = work[work[detected["organism"]].apply(is_human_organism)].copy()
        log(f"[{dataset_key}] HUMAN FILTER: {n0} -> {len(work)}")
    elif human_only and dataset_key != "human_complexes":
        log(f"[{dataset_key}] WARNING: organism column not detected; human filter not applied.")

    rows = []

    for idx, row in work.iterrows():
        complex_id = clean_text(row.get(detected["complex_id"], "")) if detected["complex_id"] else f"{dataset_key}_{idx}"
        complex_name = clean_text(row.get(detected["complex_name"], "")) if detected["complex_name"] else ""
        organism = clean_text(row.get(detected["organism"], "")) if detected["organism"] else ""
        pubmed_id = clean_text(row.get(detected["pubmed_id"], "")) if detected["pubmed_id"] else ""
        cell_line = clean_text(row.get(detected["cell_line"], "")) if detected["cell_line"] else ""
        purification_method = clean_text(row.get(detected["purification_method"], "")) if detected["purification_method"] else ""
        go_id = clean_text(row.get(detected["go_id"], "")) if detected["go_id"] else ""
        function_text = clean_text(row.get(detected["function"], "")) if detected["function"] else ""

        symbols = split_multi_value(row.get(detected["subunits_gene_symbol"], "")) if detected["subunits_gene_symbol"] else []
        uniprots = split_multi_value(row.get(detected["subunits_uniprot"], "")) if detected["subunits_uniprot"] else []
        entrezs = split_multi_value(row.get(detected["subunits_entrez"], "")) if detected["subunits_entrez"] else []

        max_len = max(len(symbols), len(uniprots), len(entrezs), 0)

        if max_len == 0:
            continue

        symbols = symbols + [""] * (max_len - len(symbols))
        uniprots = uniprots + [""] * (max_len - len(uniprots))
        entrezs = entrezs + [""] * (max_len - len(entrezs))

        for i in range(max_len):
            approved, map_source = resolve_gene(
                symbols[i],
                entrezs[i],
                uniprots[i],
                symbol_to_approved,
                entrez_to_approved,
                uniprot_to_approved,
            )

            rows.append({
                "dataset_key": dataset_key,
                "complex_id": complex_id,
                "complex_name": complex_name,
                "organism": organism,
                "pubmed_id": pubmed_id,
                "cell_line": cell_line,
                "purification_method": purification_method,
                "go_id": go_id,
                "function_text": function_text,
                "raw_subunit_index": i + 1,
                "raw_complex_subunit_count": max_len,
                "raw_gene_symbol": symbols[i],
                "raw_uniprot_id": uniprots[i],
                "raw_entrez_id": entrezs[i],
                "approved_symbol": approved,
                "map_source": map_source,
            })

    long_df = pd.DataFrame(rows)

    if long_df.empty:
        log(f"[{dataset_key}] No expanded rows.")
        return long_df

    n0 = len(long_df)
    long_df = long_df[long_df["approved_symbol"].notna()].copy()

    log(f"[{dataset_key}] SUBUNIT ROWS BEFORE HGNC MAP: {n0}")
    log(f"[{dataset_key}] SUBUNIT ROWS AFTER HGNC MAP:  {len(long_df)}")

    if long_df.empty:
        return long_df

    long_df = long_df.drop_duplicates(
        subset=["dataset_key", "complex_id", "approved_symbol"]
    ).copy()

    size_df = (
        long_df.groupby(["dataset_key", "complex_id"])["approved_symbol"]
        .nunique()
        .rename("mapped_complex_size")
        .reset_index()
    )

    long_df = long_df.merge(size_df, on=["dataset_key", "complex_id"], how="left")

    log(f"[{dataset_key}] UNIQUE GENES:     {long_df['approved_symbol'].nunique()}")
    log(f"[{dataset_key}] UNIQUE COMPLEXES: {long_df['complex_id'].nunique()}")

    return long_df


# =============================================================================
# Drug-complex interaction parser
# =============================================================================

def parse_formal_interaction_string(x):
    x = clean_text(x)
    if not x:
        return {}

    x = x.strip().strip("[]")
    parts = [clean_text(p) for p in x.split(";")]

    out = {}
    if len(parts) >= 1:
        out["target_gene"] = parts[0]
    if len(parts) >= 2:
        out["complex_name"] = parts[1]
    if len(parts) >= 3:
        out["complex_id"] = parts[2]
    if len(parts) >= 4:
        out["drug"] = parts[3]
    if len(parts) >= 5:
        out["drug_activity_1"] = parts[4]
    if len(parts) >= 6:
        out["drug_activity_2"] = parts[5]
    if len(parts) >= 7:
        out["pubmed_id"] = parts[6]

    return out


def build_drug_complex_interactions_long(
    df,
    symbol_to_approved,
    entrez_to_approved,
    uniprot_to_approved,
):
    if df is None or df.empty:
        return pd.DataFrame()

    cols = list(df.columns)

    formal_col = None
    for c in cols:
        cn = normalize_colname(c)
        if "formal" in cn or "description" in cn or "interaction" in cn:
            formal_col = c
            break

    gene_col = detect_column(
        cols,
        candidates=["Drug target gene/protein", "Gene", "Target", "Target gene"],
        contains_any=["gene"],
    )
    complex_name_col = detect_column(
        cols,
        candidates=["Complex name", "ComplexName"],
        contains_all=["complex", "name"],
    )
    complex_id_col = detect_column(
        cols,
        candidates=["CORUM-ID", "CORUM ID", "ComplexID", "Complex ID"],
        contains_any=["id"],
    )
    drug_col = detect_column(cols, candidates=["Drug", "drug"])
    activity1_col = detect_column(cols, candidates=["Drug activity 1", "activity 1", "activity1"], contains_all=["activity"])
    activity2_col = detect_column(cols, candidates=["Drug activity 2", "activity 2", "activity2"], contains_all=["activity"])
    pubmed_col = detect_column(cols, candidates=["PMID", "PubMed ID", "PubMedID"], contains_any=["pubmed"])

    detected = {
        "formal_col": formal_col,
        "gene_col": gene_col,
        "complex_name_col": complex_name_col,
        "complex_id_col": complex_id_col,
        "drug_col": drug_col,
        "activity1_col": activity1_col,
        "activity2_col": activity2_col,
        "pubmed_col": pubmed_col,
    }

    print_detected("[DETECTED COLUMNS: drug_complex_interactions]", detected)

    rows = []

    for _, row in df.iterrows():
        parsed = {}

        if formal_col is not None:
            parsed = parse_formal_interaction_string(row.get(formal_col, ""))

        target_gene = clean_text(row.get(gene_col, "")) if gene_col else parsed.get("target_gene", "")
        complex_name = clean_text(row.get(complex_name_col, "")) if complex_name_col else parsed.get("complex_name", "")
        complex_id = clean_text(row.get(complex_id_col, "")) if complex_id_col else parsed.get("complex_id", "")
        drug = clean_text(row.get(drug_col, "")) if drug_col else parsed.get("drug", "")
        activity1 = clean_text(row.get(activity1_col, "")) if activity1_col else parsed.get("drug_activity_1", "")
        activity2 = clean_text(row.get(activity2_col, "")) if activity2_col else parsed.get("drug_activity_2", "")
        pubmed_id = clean_text(row.get(pubmed_col, "")) if pubmed_col else parsed.get("pubmed_id", "")

        approved, map_source = resolve_gene(
            target_gene,
            "",
            "",
            symbol_to_approved,
            entrez_to_approved,
            uniprot_to_approved,
        )

        rows.append({
            "target_gene_raw": target_gene,
            "approved_symbol": approved,
            "map_source": map_source,
            "complex_name": complex_name,
            "complex_id": complex_id,
            "drug": drug,
            "drug_activity_1": activity1,
            "drug_activity_2": activity2,
            "pubmed_id": pubmed_id,
        })

    out = pd.DataFrame(rows)

    if out.empty:
        return out

    n0 = len(out)
    out = out[out["approved_symbol"].notna()].copy()

    log(f"[DRUG COMPLEX INTERACTIONS] ROWS BEFORE HGNC MAP: {n0}")
    log(f"[DRUG COMPLEX INTERACTIONS] ROWS AFTER HGNC MAP:  {len(out)}")

    return out


# =============================================================================
# Feature construction
# =============================================================================

def add_membership_features(feature_index, long_df, dataset_key, prefix):
    out = pd.DataFrame(index=feature_index)

    base_cols = [
        f"{prefix}_complex_membership_count",
        f"{prefix}_complex_name_count",
        f"{prefix}_pubmed_count",
        f"{prefix}_complex_size_mean",
        f"{prefix}_complex_size_median",
        f"{prefix}_complex_size_min",
        f"{prefix}_complex_size_max",
        f"{prefix}_small_2_3_complex_count",
        f"{prefix}_medium_4_10_complex_count",
        f"{prefix}_large_gt10_complex_count",
        f"{prefix}_has_any_membership",
    ]

    if long_df is None or long_df.empty or "dataset_key" not in long_df.columns:
        for c in base_cols:
            out[c] = 0
        return out

    sub = long_df[long_df["dataset_key"] == dataset_key].copy()

    if sub.empty:
        for c in base_cols:
            out[c] = 0
        return out

    grouped = sub.groupby("approved_symbol", sort=True)

    out[f"{prefix}_complex_membership_count"] = grouped["complex_id"].nunique()
    out[f"{prefix}_complex_name_count"] = grouped["complex_name"].agg(safe_nunique)
    out[f"{prefix}_pubmed_count"] = grouped["pubmed_id"].agg(safe_nunique)

    out[f"{prefix}_complex_size_mean"] = grouped["mapped_complex_size"].mean()
    out[f"{prefix}_complex_size_median"] = grouped["mapped_complex_size"].median()
    out[f"{prefix}_complex_size_min"] = grouped["mapped_complex_size"].min()
    out[f"{prefix}_complex_size_max"] = grouped["mapped_complex_size"].max()

    tmp = sub[["approved_symbol", "complex_id", "mapped_complex_size"]].drop_duplicates().copy()

    tmp["small_2_3"] = tmp["mapped_complex_size"].between(2, 3, inclusive="both").astype(int)
    tmp["medium_4_10"] = tmp["mapped_complex_size"].between(4, 10, inclusive="both").astype(int)
    tmp["large_gt10"] = (tmp["mapped_complex_size"] > 10).astype(int)

    out[f"{prefix}_small_2_3_complex_count"] = tmp.groupby("approved_symbol")["small_2_3"].sum()
    out[f"{prefix}_medium_4_10_complex_count"] = tmp.groupby("approved_symbol")["medium_4_10"].sum()
    out[f"{prefix}_large_gt10_complex_count"] = tmp.groupby("approved_symbol")["large_gt10"].sum()

    out[f"{prefix}_has_any_membership"] = (
        out[f"{prefix}_complex_membership_count"].fillna(0) > 0
    ).astype(int)

    return out.reindex(feature_index).fillna(0)


def build_functional_group_features(feature_index, long_df):
    out = pd.DataFrame(index=feature_index)

    if long_df is None or long_df.empty or "dataset_key" not in long_df.columns:
        out["corum_functional_complex_group_membership_count"] = 0
        out["corum_functional_complex_group_name_count"] = 0
        out["corum_functional_complex_group_go_count"] = 0
        out["corum_has_functional_complex_group"] = 0
        return out

    sub = long_df[long_df["dataset_key"] == "functional_complex_groups"].copy()

    if sub.empty:
        out["corum_functional_complex_group_membership_count"] = 0
        out["corum_functional_complex_group_name_count"] = 0
        out["corum_functional_complex_group_go_count"] = 0
        out["corum_has_functional_complex_group"] = 0
        return out

    grouped = sub.groupby("approved_symbol", sort=True)

    out["corum_functional_complex_group_membership_count"] = grouped["complex_id"].nunique()
    out["corum_functional_complex_group_name_count"] = grouped["complex_name"].agg(safe_nunique)
    out["corum_functional_complex_group_go_count"] = grouped["go_id"].agg(safe_nunique)

    out["corum_has_functional_complex_group"] = (
        out["corum_functional_complex_group_membership_count"].fillna(0) > 0
    ).astype(int)

    return out.reindex(feature_index).fillna(0)


def build_drug_interaction_features(feature_index, drug_long):
    out = pd.DataFrame(index=feature_index)

    base_cols = [
        "corum_drug_complex_interaction_count",
        "corum_drug_complex_unique_drug_count",
        "corum_drug_complex_unique_complex_count",
        "corum_drug_complex_pubmed_count",
        "corum_drug_complex_activation_count",
        "corum_drug_complex_inhibition_count",
        "corum_drug_complex_modification_count",
        "corum_drug_complex_formation_count",
        "corum_drug_complex_dissociation_count",
        "corum_drug_complex_function_count",
        "corum_has_drug_complex_interaction",
    ]

    if drug_long is None or drug_long.empty:
        for c in base_cols:
            out[c] = 0
        return out

    grouped = drug_long.groupby("approved_symbol", sort=True)

    out["corum_drug_complex_interaction_count"] = grouped.size()
    out["corum_drug_complex_unique_drug_count"] = grouped["drug"].agg(safe_nunique)
    out["corum_drug_complex_unique_complex_count"] = grouped["complex_id"].agg(safe_nunique)
    out["corum_drug_complex_pubmed_count"] = grouped["pubmed_id"].agg(safe_nunique)

    tmp = drug_long.copy()
    tmp["act1"] = tmp["drug_activity_1"].astype(str).str.lower()
    tmp["act2"] = tmp["drug_activity_2"].astype(str).str.lower()

    keyword_map = {
        "activation": ["activate", "activates", "activation"],
        "inhibition": ["inhibit", "inhibits", "inhibition"],
        "modification": ["modify", "modifies", "modification"],
        "formation": ["formation"],
        "dissociation": ["dissociation", "dissociate"],
        "function": ["function"],
    }

    for label, kws in keyword_map.items():
        flag_col = f"flag_{label}"
        pattern = "|".join(re.escape(k) for k in kws)

        tmp[flag_col] = (
            tmp["act1"].str.contains(pattern, na=False)
            | tmp["act2"].str.contains(pattern, na=False)
        ).astype(int)

        out[f"corum_drug_complex_{label}_count"] = tmp.groupby("approved_symbol")[flag_col].sum()

    out["corum_has_drug_complex_interaction"] = (
        out["corum_drug_complex_interaction_count"].fillna(0) > 0
    ).astype(int)

    return out.reindex(feature_index).fillna(0)


def build_all_gene_features(hgnc_keep, membership_long, drug_long):
    line()
    log("[BUILD FINAL GENE FEATURES]")

    feature_index = pd.Index(
        sorted(hgnc_keep["approved_symbol"].dropna().unique()),
        name="approved_symbol",
    )

    dataset_prefixes = {
        "human_complexes": "corum_human",
        "complete_complexes": "corum_complete",
        "drug_target_complexes": "corum_drug_target_complex",
        "splice_variant_complexes": "corum_splice_variant",
        "partial_complexes": "corum_partial",
    }

    parts = []

    for dataset_key, prefix in dataset_prefixes.items():
        parts.append(add_membership_features(feature_index, membership_long, dataset_key, prefix))

    parts.append(build_functional_group_features(feature_index, membership_long))
    parts.append(build_drug_interaction_features(feature_index, drug_long))

    features = pd.concat(parts, axis=1)
    features = features.reindex(feature_index).fillna(0)

    features["corum_any_curated_human_complex_evidence"] = (
        (features.get("corum_human_has_any_membership", 0).astype(float) > 0)
        | (features.get("corum_drug_target_complex_has_any_membership", 0).astype(float) > 0)
        | (features.get("corum_drug_complex_interaction_count", 0).astype(float) > 0)
        | (features.get("corum_functional_complex_group_membership_count", 0).astype(float) > 0)
    ).astype(int)

    features["corum_total_complex_membership_count_combined"] = (
        features.get("corum_human_complex_membership_count", 0).astype(float)
        + features.get("corum_drug_target_complex_complex_membership_count", 0).astype(float)
        + features.get("corum_splice_variant_complex_membership_count", 0).astype(float)
        + features.get("corum_partial_complex_membership_count", 0).astype(float)
    )

    for c in list(features.columns):
        if c.startswith("corum_") and pd.api.types.is_numeric_dtype(features[c]):
            if not c.startswith("corum_has_") and not c.endswith("_mean") and not c.endswith("_median") and not c.endswith("_min") and not c.endswith("_max"):
                features[f"{c}_log1p"] = np.log1p(features[c].astype(float))

    features = features.reset_index()

    log(f"[FEATURE SHAPE] {features.shape}")

    return features


def merge_with_hgnc(hgnc_keep, features):
    line()
    log("[MERGE WITH HGNC UNIVERSE]")

    merged = hgnc_keep.merge(features, on="approved_symbol", how="left")

    feature_cols = [c for c in merged.columns if c.startswith("corum_")]

    for c in feature_cols:
        merged[c] = merged[c].fillna(0)

    log(f"[MERGED SHAPE] {merged.shape}")

    if "corum_any_curated_human_complex_evidence" in merged.columns:
        log(f"[GENES WITH ANY CORUM EVIDENCE] {int(merged['corum_any_curated_human_complex_evidence'].sum())}")

    return merged


# =============================================================================
# Summary
# =============================================================================

def write_summary(path, args, found_files, tables, membership_long, drug_long, features, merged):
    path = Path(path)

    with open(path, "w", encoding="utf-8") as f:
        f.write("FEATURE 20: CORUM MULTI-FILE PROTEIN COMPLEX FEATURES\n")
        f.write("=" * 80 + "\n")
        f.write(f"Script: {SCRIPT_NAME}\n")
        f.write(f"Human only: {not args.include_nonhuman}\n")
        f.write(f"Protein coding only: {not args.no_protein_coding_filter}\n")
        f.write(f"Limit rows: {args.limit_rows}\n")
        f.write("\nFound files:\n")

        for k, p in found_files.items():
            f.write(f"  {k}: {p}\n")

        f.write("\nLoaded table shapes:\n")
        for k, df in tables.items():
            if df is None:
                f.write(f"  {k}: NOT LOADED\n")
            else:
                f.write(f"  {k}: {df.shape}\n")

        f.write("\nLong tables:\n")
        f.write(f"  membership_long: {membership_long.shape if membership_long is not None else None}\n")
        f.write(f"  drug_long: {drug_long.shape if drug_long is not None else None}\n")

        if membership_long is not None and not membership_long.empty and "dataset_key" in membership_long.columns:
            f.write("\nMembership rows by dataset:\n")
            f.write(str(membership_long["dataset_key"].value_counts()) + "\n")

        f.write("\nFeature matrices:\n")
        f.write(f"  features: {features.shape}\n")
        f.write(f"  merged: {merged.shape}\n")

        if "corum_any_curated_human_complex_evidence" in merged.columns:
            f.write(
                f"  genes with any CORUM evidence: "
                f"{int(merged['corum_any_curated_human_complex_evidence'].sum())}\n"
            )


# =============================================================================
# CLI
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Feature 20: CORUM multi-file protein complex features"
    )

    parser.add_argument("--hgnc", default=DEFAULT_HGNC)
    parser.add_argument("--dbdir", default=DEFAULT_DBDIR)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)

    parser.add_argument("--include-nonhuman", action="store_true")
    parser.add_argument("--no-protein-coding-filter", action="store_true")
    parser.add_argument("--limit-rows", type=int, default=None)

    return parser.parse_args()


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()

    t0 = time.time()

    hgnc_path = Path(args.hgnc)
    dbdir = ensure_dir(args.dbdir)
    outdir = ensure_dir(args.outdir)
    processed_dir = ensure_dir(outdir / "processed")

    line()
    log("FEATURE 20: CORUM MULTI-FILE PROTEIN COMPLEX FEATURES")
    line()
    log(f"[HGNC]           {hgnc_path}")
    log(f"[DBDIR]          {dbdir.resolve()}")
    log(f"[OUTDIR]         {outdir.resolve()}")
    log(f"[DOWNLOAD]       disabled - local files only")
    log(f"[HUMAN ONLY]     {not args.include_nonhuman}")
    log(f"[LIMIT ROWS]     {args.limit_rows if args.limit_rows is not None else 'none'}")
    line()

    found_files = {}

    for dataset_key in DATASET_PATTERNS.keys():
        p = find_dataset_file(dataset_key, dbdir)
        found_files[dataset_key] = p

        if p is not None:
            log(f"[FOUND] {dataset_key}: {p}")
        else:
            if dataset_key in {"human_complexes", "complete_complexes"}:
                log(f"[MISSING IMPORTANT] {dataset_key}: not found in {dbdir}")
            else:
                log(f"[MISSING OPTIONAL] {dataset_key}: not found in {dbdir}")

    if found_files.get("human_complexes") is None and found_files.get("complete_complexes") is None:
        die(
            "No main CORUM complex file found. Expected one of:\n"
            "  corum_humanComplexes.txt\n"
            "  corum_allComplexes.txt\n"
            f"in {dbdir}"
        )

    hgnc_keep, symbol_to_approved, entrez_to_approved, uniprot_to_approved = read_hgnc(
        hgnc_path,
        protein_coding_only=not args.no_protein_coding_filter,
    )

    tables = {}

    for dataset_key, path in found_files.items():
        if path is None:
            tables[dataset_key] = None
            continue

        try:
            tables[dataset_key] = read_corum_table(
                path=path,
                dataset_key=dataset_key,
                limit_rows=args.limit_rows,
            )
        except Exception as e:
            log(f"[WARNING] Failed to read {dataset_key}: {repr(e)}")
            tables[dataset_key] = None

    membership_dataset_keys = [
        "human_complexes",
        "complete_complexes",
        "drug_target_complexes",
        "splice_variant_complexes",
        "partial_complexes",
        "functional_complex_groups",
    ]

    membership_parts = []

    for dataset_key in membership_dataset_keys:
        df = tables.get(dataset_key)

        if df is None or df.empty:
            log(f"[SKIP MEMBERSHIP] {dataset_key}")
            continue

        part = build_complex_membership_long(
            df=df,
            dataset_key=dataset_key,
            symbol_to_approved=symbol_to_approved,
            entrez_to_approved=entrez_to_approved,
            uniprot_to_approved=uniprot_to_approved,
            human_only=not args.include_nonhuman,
        )

        if part is not None and not part.empty:
            membership_parts.append(part)

    if membership_parts:
        membership_long = pd.concat(membership_parts, ignore_index=True)
    else:
        membership_long = pd.DataFrame(
            columns=[
                "dataset_key",
                "complex_id",
                "complex_name",
                "organism",
                "pubmed_id",
                "cell_line",
                "purification_method",
                "go_id",
                "function_text",
                "approved_symbol",
                "mapped_complex_size",
            ]
        )

    drug_long = pd.DataFrame()

    if tables.get("drug_complex_interactions") is not None:
        drug_long = build_drug_complex_interactions_long(
            df=tables.get("drug_complex_interactions"),
            symbol_to_approved=symbol_to_approved,
            entrez_to_approved=entrez_to_approved,
            uniprot_to_approved=uniprot_to_approved,
        )
    else:
        log("[SKIP DRUG COMPLEX INTERACTIONS] missing")

    features = build_all_gene_features(
        hgnc_keep=hgnc_keep,
        membership_long=membership_long,
        drug_long=drug_long,
    )

    merged = merge_with_hgnc(
        hgnc_keep=hgnc_keep,
        features=features,
    )

    membership_out = processed_dir / "feature20_corum_long_memberships.csv.gz"
    drug_out = processed_dir / "feature20_corum_drug_complex_interactions_long.csv.gz"
    features_out = processed_dir / "feature20_corum_gene_features.csv"
    merged_out = processed_dir / "feature20_corum_gene_features_hgnc_merged.csv"
    summary_out = processed_dir / "feature20_corum_summary.txt"

    line()
    log("[WRITE OUTPUTS]")

    membership_long.to_csv(membership_out, index=False, compression="gzip")
    log(f"[WRITE] {membership_out}")

    drug_long.to_csv(drug_out, index=False, compression="gzip")
    log(f"[WRITE] {drug_out}")

    features.to_csv(features_out, index=False)
    log(f"[WRITE] {features_out}")

    merged.to_csv(merged_out, index=False)
    log(f"[WRITE] {merged_out}")

    write_summary(
        path=summary_out,
        args=args,
        found_files=found_files,
        tables=tables,
        membership_long=membership_long,
        drug_long=drug_long,
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