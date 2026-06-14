#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature21_PhosphoSitePlus.py

FEATURE 21: PhosphoSitePlus PTM site count features.

Fixed version:
- FIXED: import io
- FIXED: removed low_memory=False from pd.read_csv(..., engine="python")
- FIXED: avoids adding log1p columns one-by-one to reduce fragmentation warnings
- Reads local PSP .gz files from feature_databases/PhosphoSitePlus/
- Builds HGNC gene-level PTM + kinase-substrate features

Run:
    python Feature21_PhosphoSitePlus.py

Optional:
    python Feature21_PhosphoSitePlus.py --limit-rows 1000
    python Feature21_PhosphoSitePlus.py --download
"""

import argparse
import gzip
import io
import re
import ssl
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd


# =============================================================================
# Defaults
# =============================================================================

SCRIPT_NAME = "Feature21_PhosphoSitePlus.py"

DEFAULT_HGNC = "databases/HGNC/hgnc_complete_set.txt"
DEFAULT_DBDIR = "feature_databases/PhosphoSitePlus"
DEFAULT_OUTDIR = "feature21_phosphositeplus"

BASE_URLS = [
    "https://www.phosphosite.org/downloads",
    "http://www.phosphosite.org/downloads",
]

NULL_STRINGS = {"", "none", "null", "nan", "na", "-", "n/a"}

PTM_DATASETS = {
    "phosphorylation": {
        "output_name": "Phosphorylation_site_dataset.gz",
        "candidate_filenames": ["Phosphorylation_site_dataset.gz", "Phosphorylation_site_dataset"],
    },
    "acetylation": {
        "output_name": "Acetylation_site_dataset.gz",
        "candidate_filenames": ["Acetylation_site_dataset.gz", "Acetylation_site_dataset"],
    },
    "ubiquitination": {
        "output_name": "Ubiquitination_site_dataset.gz",
        "candidate_filenames": ["Ubiquitination_site_dataset.gz", "Ubiquitination_site_dataset"],
    },
    "methylation": {
        "output_name": "Methylation_site_dataset.gz",
        "candidate_filenames": ["Methylation_site_dataset.gz", "Methylation_site_dataset"],
    },
    "sumoylation": {
        "output_name": "Sumoylation_site_dataset.gz",
        "candidate_filenames": ["Sumoylation_site_dataset.gz", "Sumoylation_site_dataset"],
    },
    "ogalnac": {
        "output_name": "O-GalNAc_site_dataset.gz",
        "candidate_filenames": [
            "O-GalNAc_site_dataset.gz",
            "O-GalNAc_site_dataset",
            "O_GalNAc_site_dataset.gz",
            "O_GalNAc_site_dataset",
        ],
    },
    "oglcna": {
        "output_name": "O-GlcNAc_site_dataset.gz",
        "candidate_filenames": [
            "O-GlcNAc_site_dataset.gz",
            "O-GlcNAc_site_dataset",
            "O_GlcNAc_site_dataset.gz",
            "O_GlcNAc_site_dataset",
        ],
    },
    "regulatory": {
        "output_name": "Regulatory_sites.gz",
        "candidate_filenames": ["Regulatory_sites.gz", "Regulatory_sites"],
    },
    "disease_associated": {
        "output_name": "Disease-associated_sites.gz",
        "candidate_filenames": [
            "Disease-associated_sites.gz",
            "Disease-associated_sites",
            "Disease_associated_sites.gz",
            "Disease_associated_sites",
        ],
    },
}

KINASE_DATASET = {
    "kinase_substrate": {
        "output_name": "Kinase_Substrate_Dataset.gz",
        "candidate_filenames": ["Kinase_Substrate_Dataset.gz", "Kinase_Substrate_Dataset"],
    }
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
    Path(path).mkdir(parents=True, exist_ok=True)
    return Path(path)


def file_size_mb(path):
    path = Path(path)
    if not path.exists():
        return 0.0
    return path.stat().st_size / (1024 * 1024)


# =============================================================================
# Download helpers
# =============================================================================

def make_ssl_context(verify_ssl=False):
    if verify_ssl:
        return None
    return ssl._create_unverified_context()


def candidate_urls(dataset_info):
    urls = []
    for base in BASE_URLS:
        for fname in dataset_info["candidate_filenames"]:
            urls.append(f"{base}/{fname}")
    return urls


def download_one(dataset_key, dataset_info, dbdir, force=False, verify_ssl=False):
    dbdir = Path(dbdir)
    ensure_dir(dbdir)

    dest = dbdir / dataset_info["output_name"]

    if dest.exists() and dest.stat().st_size > 0 and not force:
        log(f"[DOWNLOAD SKIP] {dataset_key}: {dest} ({file_size_mb(dest):.2f} MB)")
        return dest

    ssl_context = make_ssl_context(verify_ssl=verify_ssl)
    last_error = None

    for url in candidate_urls(dataset_info):
        try:
            line()
            log(f"[DOWNLOAD ATTEMPT] {dataset_key}")
            log(f"[URL] {url}")
            log(f"[TO]  {dest}")

            req = urllib.request.Request(
                url,
                headers={"User-Agent": "Mozilla/5.0 Feature21_PhosphoSitePlus.py"},
            )

            with urllib.request.urlopen(req, timeout=180, context=ssl_context) as response:
                data = response.read()

            if not data:
                raise RuntimeError("Downloaded zero bytes.")

            head = data[:500].decode("latin-1", errors="replace").lower()

            if "<html" in head or "<!doctype html" in head:
                raise RuntimeError("Downloaded HTML/login page instead of data file.")

            with open(dest, "wb") as f:
                f.write(data)

            log(f"[DOWNLOAD DONE] {dataset_key}: {dest} ({file_size_mb(dest):.2f} MB)")
            return dest

        except Exception as e:
            last_error = e
            log(f"[DOWNLOAD FAILED] {url}")
            log(f"[REASON]          {repr(e)}")

    log(f"[WARNING] Could not download {dataset_key}.")
    log(f"[WARNING] Please manually place the file here: {dest}")
    log(f"[WARNING] Last error: {repr(last_error)}")
    return None


def find_existing_file(dataset_key, dataset_info, dbdir):
    dbdir = Path(dbdir)

    preferred = dbdir / dataset_info["output_name"]
    if preferred.exists() and preferred.stat().st_size > 0:
        return preferred

    for fname in dataset_info["candidate_filenames"]:
        p = dbdir / fname
        if p.exists() and p.stat().st_size > 0:
            return p

    key_clean = dataset_key.lower().replace("_", "").replace("-", "")

    for p in dbdir.glob("*"):
        if not p.is_file():
            continue
        name_clean = p.name.lower().replace("_", "").replace("-", "")
        if key_clean in name_clean and p.stat().st_size > 0:
            return p

    return None


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
    return series.astype(str).replace("", np.nan).dropna().nunique()


def count_numeric_nonzero(series):
    vals = pd.to_numeric(series.astype(str).replace("", np.nan), errors="coerce")
    return int((vals.fillna(0) > 0).sum())


def count_nonempty(series):
    vals = series.astype(str).str.strip()
    vals = vals[~vals.str.lower().isin(NULL_STRINGS)]
    return int((vals != "").sum())


def is_human_value(x):
    x = clean_text(x).lower()
    return x in {"human, homo sapiens", "human", "homo sapiens"} or "homo sapiens" in x


def parse_site_position(site):
    site = clean_text(site)
    if not site:
        return np.nan
    m = re.search(r"([A-Za-z])\s*([0-9]+)", site)
    if not m:
        return np.nan
    return int(m.group(2))


def parse_site_residue(site):
    site = clean_text(site)
    if not site:
        return ""
    m = re.search(r"([A-Za-z])\s*([0-9]+)", site)
    if not m:
        return ""
    return m.group(1).upper()


# =============================================================================
# HGNC
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
        low_memory=False,
        keep_default_na=False,
    )

    n0 = len(hgnc)

    if "symbol" not in hgnc.columns:
        die("HGNC file missing required column: symbol")

    if protein_coding_only:
        if "locus_group" in hgnc.columns:
            hgnc = hgnc[hgnc["locus_group"].astype(str).str.lower().eq("protein-coding gene")].copy()
        elif "locus_type" in hgnc.columns:
            hgnc = hgnc[hgnc["locus_type"].astype(str).str.lower().str.contains("protein-coding", na=False)].copy()
        else:
            log("[WARNING] Could not find locus_group/locus_type; protein-coding filter not applied.")

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
    if uniprot and uniprot in uniprot_to_approved:
        return uniprot_to_approved[uniprot], "uniprot"

    symbol = clean_text(symbol)
    if symbol and symbol.upper() in symbol_to_approved:
        return symbol_to_approved[symbol.upper()], "symbol"

    return None, "unmapped"


# =============================================================================
# Robust PSP reader
# =============================================================================

def read_text_lines_safely(path):
    path = Path(path)

    if str(path).lower().endswith(".gz"):
        with gzip.open(path, "rb") as f:
            raw = f.read()
    else:
        with open(path, "rb") as f:
            raw = f.read()

    text = raw.decode("utf-8", errors="replace")
    return text.splitlines()


def is_standard_site_header(fields):
    norm = [normalize_colname(x) for x in fields]
    return (
        "gene" in norm
        and "protein" in norm
        and "accid" in norm
        and "organism" in norm
        and "modrsd" in norm
        and "sitegrpid" in norm
    )


def is_kinase_header(fields):
    norm = [normalize_colname(x) for x in fields]
    return (
        "gene" in norm
        and "kinase" in norm
        and "kinaccid" in norm
        and "kinorganism" in norm
        and "substrate" in norm
        and "subaccid" in norm
        and "subgene" in norm
        and "suborganism" in norm
        and "submodrsd" in norm
        and "sitegrpid" in norm
    )


def find_header_line(lines, dataset_key):
    for i, line_raw in enumerate(lines):
        fields = line_raw.rstrip("\n\r").split("\t")

        if dataset_key == "kinase_substrate":
            if is_kinase_header(fields):
                return i, fields
        else:
            if is_standard_site_header(fields):
                return i, fields

    return None, None


def read_psp_table(path, dataset_key, limit_rows=None):
    path = Path(path)

    if not path.exists():
        die(f"Missing PhosphoSitePlus file: {path}")

    line()
    log(f"[READ PSP] {dataset_key}")
    log(f"[PATH] {path}")

    lines = read_text_lines_safely(path)
    header_idx, header_fields = find_header_line(lines, dataset_key)

    if header_idx is None:
        die(
            f"Could not detect header for {dataset_key} in {path}. "
            "The file may not be a valid PSP dataset."
        )

    log(f"[HEADER LINE] {header_idx}")
    log(f"[HEADER] {header_fields}")

    table_text = "\n".join(lines[header_idx:])

    # IMPORTANT:
    # Do NOT use low_memory=False with engine="python".
    # Your HPC pandas raised:
    # ValueError: The 'low_memory' option is not supported with the 'python' engine
    df = pd.read_csv(
        io.StringIO(table_text),
        sep="\t",
        dtype=str,
        keep_default_na=False,
        nrows=limit_rows,
        engine="python",
    )

    df.columns = [str(c).strip() for c in df.columns]

    bad_cols = [c for c in df.columns if c == "" or c.lower().startswith("unnamed")]
    if bad_cols:
        df = df.drop(columns=bad_cols)

    log(f"[SHAPE]   {df.shape}")
    log(f"[COLUMNS] {list(df.columns)}")

    return df


# =============================================================================
# PTM site parsing
# =============================================================================

def detect_site_columns(df):
    cols = list(df.columns)

    detected = {
        "protein": detect_column(cols, candidates=["PROTEIN", "Protein"]),
        "acc_id": detect_column(cols, candidates=["ACC_ID", "Accession", "UniProt", "UniProt ID"], contains_any=["acc"]),
        "gene": detect_column(cols, candidates=["GENE", "Gene", "Gene Symbol", "GeneSymbol"]),
        "gene_id": detect_column(cols, candidates=["GENE_ID", "Gene ID", "Entrez Gene ID"], contains_all=["gene", "id"]),
        "organism": detect_column(cols, candidates=["ORGANISM", "Organism", "Species"]),
        "mod_rsd": detect_column(cols, candidates=["MOD_RSD", "SITE", "Site", "Modified Residue"], contains_any=["mod"]),
        "site_grp_id": detect_column(cols, candidates=["SITE_GRP_ID", "Site Group ID"], contains_all=["site", "grp"]),
        "lt_lit": detect_column(cols, candidates=["LT_LIT", "LTP_LIT", "Low-throughput literature"]),
        "ms_lit": detect_column(cols, candidates=["MS_LIT", "MS literature"]),
        "ms_cst": detect_column(cols, candidates=["MS_CST", "CST MS"]),
        "cst_cat": detect_column(cols, candidates=["CST_CAT#", "CST_CAT", "CST catalog"]),
        "domain": detect_column(cols, candidates=["DOMAIN", "Domain"]),
        "hu_chr_loc": detect_column(cols, candidates=["HU_CHR_LOC", "Human chromosome location"]),
    }

    return detected


def print_detected(title, detected):
    line()
    log(title)
    for k, v in detected.items():
        log(f"{k:<22} = {v}")


def build_ptm_site_long(
    df,
    dataset_key,
    symbol_to_approved,
    entrez_to_approved,
    uniprot_to_approved,
    human_only=True,
):
    if df is None or df.empty:
        return pd.DataFrame()

    detected = detect_site_columns(df)
    print_detected(f"[DETECTED PSP SITE COLUMNS: {dataset_key}]", detected)

    if detected["gene"] is None and detected["acc_id"] is None and detected["gene_id"] is None:
        log(f"[WARNING] Could not detect gene/accession columns for {dataset_key}. Skipping.")
        return pd.DataFrame()

    work = df.copy()

    if human_only and detected["organism"] is not None:
        n0 = len(work)
        work = work[work[detected["organism"]].apply(is_human_value)].copy()
        log(f"[{dataset_key}] HUMAN FILTER: {n0} -> {len(work)}")
    elif human_only:
        log(f"[{dataset_key}] WARNING: organism column not detected; human filter not applied.")

    rows = []

    for _, row in work.iterrows():
        gene = clean_text(row.get(detected["gene"], "")) if detected["gene"] else ""
        gene_id = clean_text(row.get(detected["gene_id"], "")) if detected["gene_id"] else ""
        acc_id = clean_text(row.get(detected["acc_id"], "")) if detected["acc_id"] else ""
        protein = clean_text(row.get(detected["protein"], "")) if detected["protein"] else ""
        organism = clean_text(row.get(detected["organism"], "")) if detected["organism"] else ""
        mod_rsd = clean_text(row.get(detected["mod_rsd"], "")) if detected["mod_rsd"] else ""
        site_grp_id = clean_text(row.get(detected["site_grp_id"], "")) if detected["site_grp_id"] else ""
        lt_lit = clean_text(row.get(detected["lt_lit"], "")) if detected["lt_lit"] else ""
        ms_lit = clean_text(row.get(detected["ms_lit"], "")) if detected["ms_lit"] else ""
        ms_cst = clean_text(row.get(detected["ms_cst"], "")) if detected["ms_cst"] else ""
        cst_cat = clean_text(row.get(detected["cst_cat"], "")) if detected["cst_cat"] else ""
        domain = clean_text(row.get(detected["domain"], "")) if detected["domain"] else ""
        hu_chr_loc = clean_text(row.get(detected["hu_chr_loc"], "")) if detected["hu_chr_loc"] else ""

        approved, map_source = resolve_gene(
            gene,
            gene_id,
            acc_id,
            symbol_to_approved,
            entrez_to_approved,
            uniprot_to_approved,
        )

        rows.append({
            "dataset_key": dataset_key,
            "ptm_type": dataset_key,
            "protein": protein,
            "acc_id": acc_id,
            "raw_gene_symbol": gene,
            "raw_gene_id": gene_id,
            "organism": organism,
            "mod_rsd": mod_rsd,
            "site_grp_id": site_grp_id,
            "site_residue": parse_site_residue(mod_rsd),
            "site_position": parse_site_position(mod_rsd),
            "lt_lit": lt_lit,
            "ms_lit": ms_lit,
            "ms_cst": ms_cst,
            "cst_cat": cst_cat,
            "domain": domain,
            "hu_chr_loc": hu_chr_loc,
            "approved_symbol": approved,
            "map_source": map_source,
        })

    out = pd.DataFrame(rows)

    if out.empty:
        return out

    n0 = len(out)
    out = out[out["approved_symbol"].notna()].copy()

    log(f"[{dataset_key}] ROWS BEFORE HGNC MAP: {n0}")
    log(f"[{dataset_key}] ROWS AFTER HGNC MAP:  {len(out)}")

    if out.empty:
        return out

    log(f"[{dataset_key}] UNIQUE HGNC GENES:    {out['approved_symbol'].nunique()}")

    out = out.drop_duplicates(
        subset=["dataset_key", "approved_symbol", "mod_rsd", "site_grp_id"]
    ).copy()

    return out


# =============================================================================
# Kinase-substrate parsing
# =============================================================================

def detect_kinase_columns(df):
    cols = list(df.columns)

    detected = {
        "gene": detect_column(cols, candidates=["GENE"]),
        "kinase": detect_column(cols, candidates=["KINASE"]),
        "kinase_acc_id": detect_column(cols, candidates=["KIN_ACC_ID"]),
        "kinase_organism": detect_column(cols, candidates=["KIN_ORGANISM"]),
        "substrate": detect_column(cols, candidates=["SUBSTRATE"]),
        "sub_gene_id": detect_column(cols, candidates=["SUB_GENE_ID"]),
        "substrate_acc_id": detect_column(cols, candidates=["SUB_ACC_ID"]),
        "sub_gene": detect_column(cols, candidates=["SUB_GENE"]),
        "substrate_organism": detect_column(cols, candidates=["SUB_ORGANISM"]),
        "site": detect_column(cols, candidates=["SUB_MOD_RSD"]),
        "site_grp_id": detect_column(cols, candidates=["SITE_GRP_ID"]),
        "in_vivo": detect_column(cols, candidates=["IN_VIVO_RXN"]),
        "in_vitro": detect_column(cols, candidates=["IN_VITRO_RXN"]),
        "cst_cat": detect_column(cols, candidates=["CST_CAT#"]),
    }

    return detected


def build_kinase_substrate_long(
    df,
    symbol_to_approved,
    entrez_to_approved,
    uniprot_to_approved,
    human_only=True,
):
    if df is None or df.empty:
        return pd.DataFrame()

    detected = detect_kinase_columns(df)
    print_detected("[DETECTED PSP KINASE-SUBSTRATE COLUMNS]", detected)

    required = ["kinase", "kinase_acc_id", "substrate_acc_id", "sub_gene", "substrate_organism", "site"]
    missing = [x for x in required if detected.get(x) is None]

    if missing:
        log(f"[WARNING] Missing kinase-substrate columns: {missing}. Skipping kinase-substrate.")
        return pd.DataFrame()

    work = df.copy()

    if human_only:
        n0 = len(work)

        if detected["kinase_organism"] is not None:
            work = work[work[detected["kinase_organism"]].apply(is_human_value)].copy()

        if detected["substrate_organism"] is not None:
            work = work[work[detected["substrate_organism"]].apply(is_human_value)].copy()

        log(f"[KINASE-SUBSTRATE] HUMAN KINASE+SUBSTRATE FILTER: {n0} -> {len(work)}")

    rows = []

    for _, row in work.iterrows():
        kinase_raw = clean_text(row.get(detected["kinase"], ""))
        kinase_acc = clean_text(row.get(detected["kinase_acc_id"], ""))
        kinase_org = clean_text(row.get(detected["kinase_organism"], "")) if detected["kinase_organism"] else ""

        substrate_raw = clean_text(row.get(detected["sub_gene"], ""))
        substrate_acc = clean_text(row.get(detected["substrate_acc_id"], ""))
        substrate_gene_id = clean_text(row.get(detected["sub_gene_id"], "")) if detected["sub_gene_id"] else ""
        substrate_org = clean_text(row.get(detected["substrate_organism"], ""))

        site = clean_text(row.get(detected["site"], ""))
        site_grp_id = clean_text(row.get(detected["site_grp_id"], "")) if detected["site_grp_id"] else ""

        in_vivo = clean_text(row.get(detected["in_vivo"], "")) if detected["in_vivo"] else ""
        in_vitro = clean_text(row.get(detected["in_vitro"], "")) if detected["in_vitro"] else ""
        cst_cat = clean_text(row.get(detected["cst_cat"], "")) if detected["cst_cat"] else ""

        kinase_approved, kinase_map_source = resolve_gene(
            kinase_raw,
            "",
            kinase_acc,
            symbol_to_approved,
            entrez_to_approved,
            uniprot_to_approved,
        )

        substrate_approved, substrate_map_source = resolve_gene(
            substrate_raw,
            substrate_gene_id,
            substrate_acc,
            symbol_to_approved,
            entrez_to_approved,
            uniprot_to_approved,
        )

        rows.append({
            "kinase_raw": kinase_raw,
            "kinase_acc_id": kinase_acc,
            "kinase_organism": kinase_org,
            "kinase_approved_symbol": kinase_approved,
            "kinase_map_source": kinase_map_source,
            "substrate_raw": substrate_raw,
            "substrate_gene_id": substrate_gene_id,
            "substrate_acc_id": substrate_acc,
            "substrate_organism": substrate_org,
            "substrate_approved_symbol": substrate_approved,
            "substrate_map_source": substrate_map_source,
            "site": site,
            "site_grp_id": site_grp_id,
            "site_residue": parse_site_residue(site),
            "site_position": parse_site_position(site),
            "in_vivo_rxn": in_vivo,
            "in_vitro_rxn": in_vitro,
            "cst_cat": cst_cat,
        })

    out = pd.DataFrame(rows)

    if out.empty:
        return out

    n0 = len(out)
    out = out[
        out["kinase_approved_symbol"].notna()
        & out["substrate_approved_symbol"].notna()
    ].copy()

    log(f"[KINASE-SUBSTRATE] ROWS BEFORE HGNC MAP: {n0}")
    log(f"[KINASE-SUBSTRATE] ROWS AFTER HGNC MAP:  {len(out)}")

    if out.empty:
        return out

    out = out.drop_duplicates(
        subset=[
            "kinase_approved_symbol",
            "substrate_approved_symbol",
            "site",
            "site_grp_id",
        ]
    ).copy()

    log(f"[KINASE-SUBSTRATE] UNIQUE KINASES:    {out['kinase_approved_symbol'].nunique()}")
    log(f"[KINASE-SUBSTRATE] UNIQUE SUBSTRATES: {out['substrate_approved_symbol'].nunique()}")

    return out


# =============================================================================
# Feature construction
# =============================================================================

def build_site_features(feature_index, ptm_long):
    out = pd.DataFrame(index=feature_index)
    base_ptm_types = list(PTM_DATASETS.keys())

    if ptm_long is None or ptm_long.empty:
        zero_cols = {}

        for ptm in base_ptm_types:
            zero_cols[f"psp_{ptm}_site_count"] = 0
            zero_cols[f"psp_{ptm}_unique_site_group_count"] = 0
            zero_cols[f"psp_{ptm}_unique_residue_count"] = 0
            zero_cols[f"psp_{ptm}_lt_lit_site_count"] = 0
            zero_cols[f"psp_{ptm}_ms_lit_site_count"] = 0
            zero_cols[f"psp_{ptm}_ms_cst_site_count"] = 0
            zero_cols[f"psp_has_{ptm}_site"] = 0

        zero_cols["psp_total_ptm_site_count"] = 0
        zero_cols["psp_total_distinct_ptm_type_count"] = 0
        zero_cols["psp_total_unique_site_group_count"] = 0
        zero_cols["psp_total_unique_residue_count"] = 0
        zero_cols["psp_total_domain_annotated_site_count"] = 0
        zero_cols["psp_has_any_ptm_site"] = 0

        return pd.DataFrame(zero_cols, index=feature_index)

    grouped_all = ptm_long.groupby("approved_symbol", sort=True)

    pieces = []

    overall = pd.DataFrame(index=feature_index)
    overall["psp_total_ptm_site_count"] = grouped_all.size()
    overall["psp_total_distinct_ptm_type_count"] = grouped_all["ptm_type"].nunique()
    overall["psp_total_unique_site_group_count"] = grouped_all["site_grp_id"].agg(safe_nunique)
    overall["psp_total_unique_residue_count"] = grouped_all["site_residue"].agg(safe_nunique)
    overall["psp_total_domain_annotated_site_count"] = grouped_all["domain"].agg(count_nonempty)
    overall["psp_has_any_ptm_site"] = (overall["psp_total_ptm_site_count"].fillna(0) > 0).astype(int)
    pieces.append(overall)

    for ptm in base_ptm_types:
        sub = ptm_long[ptm_long["ptm_type"] == ptm].copy()
        tmp = pd.DataFrame(index=feature_index)

        if sub.empty:
            tmp[f"psp_{ptm}_site_count"] = 0
            tmp[f"psp_{ptm}_unique_site_group_count"] = 0
            tmp[f"psp_{ptm}_unique_residue_count"] = 0
            tmp[f"psp_{ptm}_lt_lit_site_count"] = 0
            tmp[f"psp_{ptm}_ms_lit_site_count"] = 0
            tmp[f"psp_{ptm}_ms_cst_site_count"] = 0
            tmp[f"psp_has_{ptm}_site"] = 0
            pieces.append(tmp)
            continue

        g = sub.groupby("approved_symbol", sort=True)

        tmp[f"psp_{ptm}_site_count"] = g.size()
        tmp[f"psp_{ptm}_unique_site_group_count"] = g["site_grp_id"].agg(safe_nunique)
        tmp[f"psp_{ptm}_unique_residue_count"] = g["site_residue"].agg(safe_nunique)
        tmp[f"psp_{ptm}_lt_lit_site_count"] = g["lt_lit"].agg(count_numeric_nonzero)
        tmp[f"psp_{ptm}_ms_lit_site_count"] = g["ms_lit"].agg(count_numeric_nonzero)
        tmp[f"psp_{ptm}_ms_cst_site_count"] = g["ms_cst"].agg(count_numeric_nonzero)
        tmp[f"psp_has_{ptm}_site"] = (tmp[f"psp_{ptm}_site_count"].fillna(0) > 0).astype(int)

        pieces.append(tmp)

    out = pd.concat(pieces, axis=1)
    return out.reindex(feature_index).fillna(0)


def build_kinase_features(feature_index, kinase_long):
    base_cols = [
        "psp_as_kinase_reaction_count",
        "psp_as_kinase_unique_substrate_count",
        "psp_as_kinase_unique_site_count",
        "psp_as_kinase_in_vivo_count",
        "psp_as_kinase_in_vitro_count",
        "psp_as_substrate_reaction_count",
        "psp_as_substrate_unique_kinase_count",
        "psp_as_substrate_unique_site_count",
        "psp_as_substrate_in_vivo_count",
        "psp_as_substrate_in_vitro_count",
        "psp_has_kinase_activity_evidence",
        "psp_has_substrate_evidence",
    ]

    if kinase_long is None or kinase_long.empty:
        return pd.DataFrame({c: 0 for c in base_cols}, index=feature_index)

    tmp = kinase_long.copy()
    tmp["flag_in_vivo"] = tmp["in_vivo_rxn"].astype(str).str.upper().eq("X").astype(int)
    tmp["flag_in_vitro"] = tmp["in_vitro_rxn"].astype(str).str.upper().eq("X").astype(int)

    out = pd.DataFrame(index=feature_index)

    kg = tmp.groupby("kinase_approved_symbol", sort=True)
    out["psp_as_kinase_reaction_count"] = kg.size()
    out["psp_as_kinase_unique_substrate_count"] = kg["substrate_approved_symbol"].nunique()
    out["psp_as_kinase_unique_site_count"] = kg["site_grp_id"].agg(safe_nunique)
    out["psp_as_kinase_in_vivo_count"] = kg["flag_in_vivo"].sum()
    out["psp_as_kinase_in_vitro_count"] = kg["flag_in_vitro"].sum()

    sg = tmp.groupby("substrate_approved_symbol", sort=True)
    out["psp_as_substrate_reaction_count"] = sg.size()
    out["psp_as_substrate_unique_kinase_count"] = sg["kinase_approved_symbol"].nunique()
    out["psp_as_substrate_unique_site_count"] = sg["site_grp_id"].agg(safe_nunique)
    out["psp_as_substrate_in_vivo_count"] = sg["flag_in_vivo"].sum()
    out["psp_as_substrate_in_vitro_count"] = sg["flag_in_vitro"].sum()

    out["psp_has_kinase_activity_evidence"] = (
        out["psp_as_kinase_reaction_count"].fillna(0) > 0
    ).astype(int)

    out["psp_has_substrate_evidence"] = (
        out["psp_as_substrate_reaction_count"].fillna(0) > 0
    ).astype(int)

    return out.reindex(feature_index).fillna(0)


def build_gene_features(hgnc_keep, ptm_long, kinase_long):
    line()
    log("[BUILD GENE FEATURES]")

    feature_index = pd.Index(
        sorted(hgnc_keep["approved_symbol"].dropna().unique()),
        name="approved_symbol",
    )

    site_features = build_site_features(feature_index, ptm_long)
    kinase_features = build_kinase_features(feature_index, kinase_long)

    features = pd.concat([site_features, kinase_features], axis=1)
    features = features.reindex(feature_index).fillna(0)

    features["psp_total_ptm_and_kinase_evidence_count"] = (
        features.get("psp_total_ptm_site_count", 0).astype(float)
        + features.get("psp_as_kinase_reaction_count", 0).astype(float)
        + features.get("psp_as_substrate_reaction_count", 0).astype(float)
    )

    features["psp_has_any_evidence"] = (
        features["psp_total_ptm_and_kinase_evidence_count"].astype(float) > 0
    ).astype(int)

    numeric_cols = [
        c for c in features.columns
        if c.startswith("psp_")
        and not c.startswith("psp_has_")
        and pd.api.types.is_numeric_dtype(features[c])
    ]

    log1p_df = pd.DataFrame(
        {f"{c}_log1p": np.log1p(features[c].astype(float)) for c in numeric_cols},
        index=features.index,
    )

    features = pd.concat([features, log1p_df], axis=1).copy()
    features = features.reset_index()

    log(f"[FEATURE SHAPE] {features.shape}")

    return features


def merge_with_hgnc(hgnc_keep, features):
    line()
    log("[MERGE WITH HGNC UNIVERSE]")

    merged = hgnc_keep.merge(features, on="approved_symbol", how="left")

    feature_cols = [c for c in merged.columns if c.startswith("psp_")]
    for c in feature_cols:
        merged[c] = merged[c].fillna(0)

    log(f"[MERGED SHAPE] {merged.shape}")

    if "psp_has_any_evidence" in merged.columns:
        log(f"[GENES WITH ANY PSP EVIDENCE] {int(merged['psp_has_any_evidence'].sum())}")

    return merged


# =============================================================================
# Summary
# =============================================================================

def write_summary(path, args, files, ptm_tables, kinase_table, ptm_long, kinase_long, features, merged):
    path = Path(path)

    with open(path, "w", encoding="utf-8") as f:
        f.write("FEATURE 21: PHOSPHOSITEPLUS PTM SITE COUNT FEATURES\n")
        f.write("=" * 80 + "\n")
        f.write(f"Script: {SCRIPT_NAME}\n")
        f.write(f"Human only: {not args.include_nonhuman}\n")
        f.write(f"Protein coding only: {not args.no_protein_coding_filter}\n")
        f.write(f"Limit rows: {args.limit_rows}\n")
        f.write("\nFiles:\n")

        for k, p in files.items():
            f.write(f"  {k}: {p}\n")

        f.write("\nPTM table shapes:\n")
        for k, df in ptm_tables.items():
            if df is None:
                f.write(f"  {k}: NOT LOADED\n")
            else:
                f.write(f"  {k}: {df.shape}\n")

        if kinase_table is None:
            f.write("\nKinase-substrate table: NOT LOADED\n")
        else:
            f.write(f"\nKinase-substrate table: {kinase_table.shape}\n")

        f.write("\nLong tables:\n")
        f.write(f"  ptm_long: {ptm_long.shape if ptm_long is not None else None}\n")
        f.write(f"  kinase_long: {kinase_long.shape if kinase_long is not None else None}\n")

        if ptm_long is not None and not ptm_long.empty:
            f.write("\nPTM rows by type:\n")
            f.write(str(ptm_long["ptm_type"].value_counts()) + "\n")

        f.write("\nFeature matrices:\n")
        f.write(f"  features: {features.shape}\n")
        f.write(f"  merged: {merged.shape}\n")

        if "psp_has_any_evidence" in merged.columns:
            f.write(f"  genes with any PSP evidence: {int(merged['psp_has_any_evidence'].sum())}\n")


# =============================================================================
# CLI
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Feature 21: PhosphoSitePlus PTM site count features"
    )

    parser.add_argument("--hgnc", default=DEFAULT_HGNC)
    parser.add_argument("--dbdir", default=DEFAULT_DBDIR)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)

    parser.add_argument("--download", action="store_true")
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--verify-ssl", action="store_true")

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
    log("FEATURE 21: PHOSPHOSITEPLUS PTM SITE COUNT FEATURES")
    line()
    log(f"[HGNC]           {hgnc_path}")
    log(f"[DBDIR]          {dbdir.resolve()}")
    log(f"[OUTDIR]         {outdir.resolve()}")
    log(f"[DOWNLOAD]       {args.download}")
    log(f"[FORCE DOWNLOAD] {args.force_download}")
    log(f"[VERIFY SSL]     {args.verify_ssl}")
    log(f"[HUMAN ONLY]     {not args.include_nonhuman}")
    log(f"[LIMIT ROWS]     {args.limit_rows if args.limit_rows is not None else 'none'}")
    line()

    all_dataset_info = {}
    all_dataset_info.update(PTM_DATASETS)
    all_dataset_info.update(KINASE_DATASET)

    files = {}

    for dataset_key, dataset_info in all_dataset_info.items():
        existing = find_existing_file(dataset_key, dataset_info, dbdir)

        if args.download:
            p = download_one(
                dataset_key=dataset_key,
                dataset_info=dataset_info,
                dbdir=dbdir,
                force=args.force_download,
                verify_ssl=args.verify_ssl,
            )
            if p is None:
                p = existing
        else:
            p = existing

        files[dataset_key] = p

    if not any(p is not None for p in files.values()):
        die(
            "No PhosphoSitePlus files were found or downloaded.\n"
            "Run with --download or manually place PSP files in:\n"
            f"{dbdir}"
        )

    hgnc_keep, symbol_to_approved, entrez_to_approved, uniprot_to_approved = read_hgnc(
        hgnc_path=hgnc_path,
        protein_coding_only=not args.no_protein_coding_filter,
    )

    ptm_tables = {}
    ptm_long_parts = []

    for ptm_key in PTM_DATASETS.keys():
        p = files.get(ptm_key)

        if p is None:
            log(f"[SKIP] {ptm_key}: file missing")
            ptm_tables[ptm_key] = None
            continue

        try:
            df = read_psp_table(p, ptm_key, limit_rows=args.limit_rows)
            ptm_tables[ptm_key] = df

            long_part = build_ptm_site_long(
                df=df,
                dataset_key=ptm_key,
                symbol_to_approved=symbol_to_approved,
                entrez_to_approved=entrez_to_approved,
                uniprot_to_approved=uniprot_to_approved,
                human_only=not args.include_nonhuman,
            )

            if long_part is not None and not long_part.empty:
                ptm_long_parts.append(long_part)

        except Exception as e:
            log(f"[WARNING] Failed to process {ptm_key}: {repr(e)}")
            ptm_tables[ptm_key] = None

    if ptm_long_parts:
        ptm_long = pd.concat(ptm_long_parts, ignore_index=True)
    else:
        ptm_long = pd.DataFrame(
            columns=[
                "dataset_key",
                "ptm_type",
                "approved_symbol",
                "mod_rsd",
                "site_grp_id",
                "domain",
                "site_residue",
                "lt_lit",
                "ms_lit",
                "ms_cst",
            ]
        )

    kinase_table = None
    kinase_long = pd.DataFrame()

    kinase_path = files.get("kinase_substrate")

    if kinase_path is not None:
        try:
            kinase_table = read_psp_table(
                kinase_path,
                "kinase_substrate",
                limit_rows=args.limit_rows,
            )

            kinase_long = build_kinase_substrate_long(
                df=kinase_table,
                symbol_to_approved=symbol_to_approved,
                entrez_to_approved=entrez_to_approved,
                uniprot_to_approved=uniprot_to_approved,
                human_only=not args.include_nonhuman,
            )

        except Exception as e:
            log(f"[WARNING] Failed to process kinase_substrate: {repr(e)}")
            kinase_table = None
            kinase_long = pd.DataFrame()
    else:
        log("[SKIP] kinase_substrate: file missing")

    features = build_gene_features(
        hgnc_keep=hgnc_keep,
        ptm_long=ptm_long,
        kinase_long=kinase_long,
    )

    merged = merge_with_hgnc(
        hgnc_keep=hgnc_keep,
        features=features,
    )

    ptm_long_out = processed_dir / "feature21_psp_ptm_long_sites.csv.gz"
    kinase_long_out = processed_dir / "feature21_psp_kinase_substrate_long.csv.gz"
    features_out = processed_dir / "feature21_psp_gene_features.csv"
    merged_out = processed_dir / "feature21_psp_gene_features_hgnc_merged.csv"
    summary_out = processed_dir / "feature21_psp_summary.txt"

    line()
    log("[WRITE OUTPUTS]")

    ptm_long.to_csv(ptm_long_out, index=False, compression="gzip")
    log(f"[WRITE] {ptm_long_out}")

    kinase_long.to_csv(kinase_long_out, index=False, compression="gzip")
    log(f"[WRITE] {kinase_long_out}")

    features.to_csv(features_out, index=False)
    log(f"[WRITE] {features_out}")

    merged.to_csv(merged_out, index=False)
    log(f"[WRITE] {merged_out}")

    write_summary(
        path=summary_out,
        args=args,
        files=files,
        ptm_tables=ptm_tables,
        kinase_table=kinase_table,
        ptm_long=ptm_long,
        kinase_long=kinase_long,
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