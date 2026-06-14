#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature22_Paralogues.py

FEATURE 22: Ensembl BioMart human paralogue features.

Fixed version
-------------
This version fixes the BioMart error:

    Attributes from multiple attribute pages are not allowed

Cause:
    The previous query mixed normal gene attributes such as hgnc_symbol,
    external_gene_name, and entrezgene_id with homology/paralogue attributes.

Fix:
    Query only BioMart-compatible paralogue attributes:
        ensembl_gene_id
        hsapiens_paralog_ensembl_gene
        hsapiens_paralog_associated_gene_name
        hsapiens_paralog_perc_id
        hsapiens_paralog_perc_id_r1

    Then map query genes and paralogue genes back to HGNC locally using:
        databases/HGNC/hgnc_complete_set.txt

Default behavior:
    Batch BioMart queries using 100 Ensembl gene IDs per request.
    This avoids huge BioMart requests and avoids mixed attribute-page errors.

Run:
    python Feature22_Paralogues.py

Quick test:
    python Feature22_Paralogues.py --limit-genes 1000

Force refresh:
    python Feature22_Paralogues.py --force-download

Use a mirror:
    python Feature22_Paralogues.py --host https://www.ensembl.org
    python Feature22_Paralogues.py --host https://asia.ensembl.org

Outputs:
    feature22_paralogues/processed/feature22_paralogues_gene_features.csv
    feature22_paralogues/processed/feature22_paralogues_gene_features_hgnc_merged.csv
    feature22_paralogues/processed/feature22_paralogues_long.csv.gz
    feature22_paralogues/processed/feature22_paralogues_summary.txt
"""

import argparse
import gzip
import io
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np
import pandas as pd


# =============================================================================
# Defaults
# =============================================================================

SCRIPT_NAME = "Feature22_Paralogues.py"

DEFAULT_HGNC = "databases/HGNC/hgnc_complete_set.txt"
DEFAULT_DBDIR = "feature_databases/Ensembl_Paralogues"
DEFAULT_OUTDIR = "feature22_paralogues"

DEFAULT_HOSTS = [
    "https://www.ensembl.org",
    "https://asia.ensembl.org",
    "https://useast.ensembl.org",
    "https://uswest.ensembl.org",
]

RAW_CACHE_NAME = "ensembl_human_paralogues_biomart.tsv.gz"
RAW_META_NAME = "ensembl_human_paralogues_biomart.meta.txt"

NULL_STRINGS = {"", "none", "null", "nan", "na", "-", "n/a"}

# IMPORTANT:
# Do not add hgnc_symbol, external_gene_name, or entrezgene_id here.
# Those can trigger:
#   Attributes from multiple attribute pages are not allowed
BIOMART_PARALOGUE_ATTRIBUTES = [
    "ensembl_gene_id",
    "hsapiens_paralog_ensembl_gene",
    "hsapiens_paralog_associated_gene_name",
    "hsapiens_paralog_perc_id",
    "hsapiens_paralog_perc_id_r1",
]

RAW_COLUMNS = [
    "query_ensembl_gene_id",
    "paralogue_ensembl_gene_id",
    "paralogue_associated_gene_name",
    "paralogue_perc_id_target_identical_to_query",
    "paralogue_perc_id_query_identical_to_target",
]


# =============================================================================
# Logging and helpers
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
        return ""
    if re.match(r"^\d+\.0$", x):
        x = x[:-2]
    if re.match(r"^\d+$", x):
        return x
    return ""


def split_pipe_values(x):
    x = clean_text(x)
    if not x:
        return []
    return [p.strip() for p in x.split("|") if p.strip()]


def split_hgnc_multi_value(x):
    return split_pipe_values(x)


def to_float_series(s):
    return pd.to_numeric(s, errors="coerce")


def safe_nunique(series):
    if series is None:
        return 0
    return series.astype(str).replace("", np.nan).dropna().nunique()


def chunked(items, size):
    for i in range(0, len(items), size):
        yield items[i:i + size]


# =============================================================================
# HGNC parsing and mapping
# =============================================================================

def read_hgnc(hgnc_path, protein_coding_only=True, limit_genes=None):
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

    if "symbol" not in hgnc.columns:
        die("HGNC file missing required column: symbol")

    n0 = len(hgnc)

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
        "ensembl_gene_id",
        "alias_symbol",
        "prev_symbol",
        "locus_group",
        "locus_type",
    ]:
        if col not in hgnc.columns:
            hgnc[col] = ""

    if limit_genes is not None:
        hgnc = hgnc.head(int(limit_genes)).copy()
        log(f"[LIMIT HGNC GENES] {limit_genes}")

    symbol_to_approved = {}
    ensembl_to_approved = {}
    entrez_to_approved = {}
    approved_to_ensembl = {}

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

        ens_values = split_pipe_values(row.get("ensembl_gene_id", ""))
        if not ens_values:
            ens_single = clean_text(row.get("ensembl_gene_id", ""))
            if ens_single:
                ens_values = [ens_single]

        for ens in ens_values:
            ens = clean_text(ens)
            if ens:
                ensembl_to_approved[ens] = approved
                approved_to_ensembl.setdefault(approved, set()).add(ens)

    keep_cols = [
        "hgnc_id",
        "approved_symbol",
        "name",
        "entrez_id",
        "ensembl_gene_id",
        "locus_group",
        "locus_type",
    ]
    keep_cols = [c for c in keep_cols if c in hgnc.columns]
    hgnc_keep = hgnc[keep_cols].drop_duplicates("approved_symbol").copy()

    query_ensembl_ids = sorted(set(ensembl_to_approved.keys()))

    line()
    log("[HGNC MAPS]")
    log(f"[SYMBOLS + ALIASES] {len(symbol_to_approved)}")
    log(f"[ENSEMBL IDS]       {len(ensembl_to_approved)}")
    log(f"[ENTREZ IDS]        {len(entrez_to_approved)}")
    log(f"[QUERY ENSEMBL IDS]  {len(query_ensembl_ids)}")

    return (
        hgnc_keep,
        symbol_to_approved,
        ensembl_to_approved,
        entrez_to_approved,
        approved_to_ensembl,
        query_ensembl_ids,
    )


def resolve_gene(symbol, ensembl_id, entrez_id, symbol_to_approved, ensembl_to_approved, entrez_to_approved):
    ensembl_id = clean_text(ensembl_id)
    if ensembl_id and ensembl_id in ensembl_to_approved:
        return ensembl_to_approved[ensembl_id], "ensembl"

    entrez_id = clean_entrez_id(entrez_id)
    if entrez_id and entrez_id in entrez_to_approved:
        return entrez_to_approved[entrez_id], "entrez"

    symbol = clean_text(symbol)
    if symbol and symbol.upper() in symbol_to_approved:
        return symbol_to_approved[symbol.upper()], "symbol"

    return None, "unmapped"


# =============================================================================
# BioMart querying
# =============================================================================

def biomart_xml_for_batch(attributes, ensembl_gene_ids):
    attr_xml = "\n".join([f'    <Attribute name="{escape(a)}" />' for a in attributes])
    ids = ",".join([escape(x) for x in ensembl_gene_ids])

    xml = f'''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Filter name="ensembl_gene_id" value="{ids}" />
{attr_xml}
  </Dataset>
</Query>
'''
    return xml


def biomart_xml_full(attributes):
    attr_xml = "\n".join([f'    <Attribute name="{escape(a)}" />' for a in attributes])

    xml = f'''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Filter name="with_hsapiens_paralog" excluded="0" />
{attr_xml}
  </Dataset>
</Query>
'''
    return xml


def post_biomart_query(host, xml, timeout=300, retries=2, sleep_seconds=5):
    url = host.rstrip("/") + "/biomart/martservice"
    data = urllib.parse.urlencode({"query": xml}).encode("utf-8")

    last_error = None

    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(
                url,
                data=data,
                headers={
                    "User-Agent": f"{SCRIPT_NAME} Python urllib/3",
                    "Content-Type": "application/x-www-form-urlencoded",
                    "Accept": "text/plain,*/*",
                },
                method="POST",
            )

            with urllib.request.urlopen(req, timeout=timeout) as response:
                raw = response.read()

            if not raw:
                raise RuntimeError("BioMart returned zero bytes")

            text_head = raw[:1500].decode("utf-8", errors="replace")
            lower_head = text_head.lower()

            if (
                "query error" in lower_head
                or "exception" in lower_head
                or "<html" in lower_head
                or "<!doctype" in lower_head
            ):
                raise RuntimeError(f"BioMart returned error-like response: {text_head[:500]}")

            return raw

        except Exception as e:
            last_error = e
            if attempt < retries:
                time.sleep(sleep_seconds)

    raise RuntimeError(f"BioMart failed for host={host}: {repr(last_error)}")


def parse_biomart_raw_bytes(raw):
    text = raw.decode("utf-8", errors="replace").strip()

    if not text:
        return pd.DataFrame(columns=RAW_COLUMNS)

    df = pd.read_csv(
        io.StringIO(text),
        sep="\t",
        header=None,
        names=RAW_COLUMNS,
        dtype=str,
        keep_default_na=False,
    )

    return df


def query_batch_against_hosts(batch_ids, hosts, timeout, retries):
    xml = biomart_xml_for_batch(BIOMART_PARALOGUE_ATTRIBUTES, batch_ids)

    last_error = None

    for host in hosts:
        try:
            raw = post_biomart_query(
                host=host,
                xml=xml,
                timeout=timeout,
                retries=retries,
            )

            df = parse_biomart_raw_bytes(raw)
            return df, host

        except Exception as e:
            last_error = e
            log(f"[BATCH HOST FAILED] host={host} reason={repr(e)}")

    raise RuntimeError(f"All hosts failed for batch. Last error: {repr(last_error)}")


def query_full_against_hosts(hosts, timeout, retries):
    xml = biomart_xml_full(BIOMART_PARALOGUE_ATTRIBUTES)

    last_error = None

    for host in hosts:
        try:
            line("-")
            log(f"[BIOMART FULL QUERY] host={host}")
            raw = post_biomart_query(
                host=host,
                xml=xml,
                timeout=timeout,
                retries=retries,
            )
            df = parse_biomart_raw_bytes(raw)
            return df, host

        except Exception as e:
            last_error = e
            log(f"[FULL HOST FAILED] host={host} reason={repr(e)}")

    raise RuntimeError(f"All hosts failed for full query. Last error: {repr(last_error)}")


def download_biomart_paralogues(
    raw_cache_path,
    meta_path,
    query_ensembl_ids,
    host=None,
    timeout=300,
    retries=2,
    batch_size=100,
    use_full_query=False,
):
    raw_cache_path = Path(raw_cache_path)
    meta_path = Path(meta_path)
    ensure_dir(raw_cache_path.parent)

    hosts = [host] if host else DEFAULT_HOSTS

    line()
    log("[BIOMART DOWNLOAD]")
    log(f"[MODE]       {'full' if use_full_query else 'batch'}")
    log(f"[HOSTS]      {hosts}")
    log(f"[BATCH SIZE] {batch_size}")
    log(f"[GENE IDS]   {len(query_ensembl_ids)}")

    all_parts = []
    used_hosts = []
    failed_batches = []

    if use_full_query:
        try:
            df, used_host = query_full_against_hosts(
                hosts=hosts,
                timeout=timeout,
                retries=retries,
            )
            all_parts.append(df)
            used_hosts.append(used_host)
        except Exception as e:
            log(f"[WARNING] Full BioMart query failed. Falling back to batch mode. Reason: {repr(e)}")
            use_full_query = False

    if not use_full_query:
        batches = list(chunked(query_ensembl_ids, int(batch_size)))
        n_batches = len(batches)

        for i, batch_ids in enumerate(batches, start=1):
            line("-")
            log(f"[BIOMART BATCH] {i}/{n_batches} genes={len(batch_ids)} first={batch_ids[0]}")

            try:
                df, used_host = query_batch_against_hosts(
                    batch_ids=batch_ids,
                    hosts=hosts,
                    timeout=timeout,
                    retries=retries,
                )

                log(f"[BATCH DONE] rows={len(df)} host={used_host}")

                if df is not None and not df.empty:
                    all_parts.append(df)

                used_hosts.append(used_host)

                # Be polite to BioMart, but keep it fast.
                time.sleep(0.15)

            except Exception as e:
                log(f"[BATCH FAILED] {i}/{n_batches}: {repr(e)}")
                failed_batches.append((i, batch_ids[0], batch_ids[-1], repr(e)))

    if all_parts:
        raw_df = pd.concat(all_parts, ignore_index=True)
    else:
        raw_df = pd.DataFrame(columns=RAW_COLUMNS)

    for c in RAW_COLUMNS:
        if c not in raw_df.columns:
            raw_df[c] = ""

    raw_df = raw_df[RAW_COLUMNS].copy()
    raw_df = raw_df.drop_duplicates().copy()

    raw_df.to_csv(
        raw_cache_path,
        sep="\t",
        index=False,
        header=False,
        compression="gzip",
    )

    with open(meta_path, "w", encoding="utf-8") as f:
        f.write("FEATURE 22 BioMart paralogue cache\n")
        f.write("=" * 80 + "\n")
        f.write(f"download_time_unix\t{int(time.time())}\n")
        f.write(f"mode\t{'full' if use_full_query else 'batch'}\n")
        f.write(f"batch_size\t{batch_size}\n")
        f.write(f"query_gene_ids\t{len(query_ensembl_ids)}\n")
        f.write(f"raw_rows\t{len(raw_df)}\n")
        f.write(f"used_hosts\t{','.join(sorted(set(used_hosts)))}\n")
        f.write(f"attributes\t{','.join(BIOMART_PARALOGUE_ATTRIBUTES)}\n")
        f.write(f"columns\t{','.join(RAW_COLUMNS)}\n")
        f.write(f"failed_batch_count\t{len(failed_batches)}\n")
        if failed_batches:
            f.write("\nFAILED_BATCHES\n")
            for item in failed_batches:
                f.write("\t".join(map(str, item)) + "\n")

    log(f"[BIOMART CACHE WRITTEN] {raw_cache_path} ({file_size_mb(raw_cache_path):.2f} MB)")
    log(f"[BIOMART META WRITTEN]  {meta_path}")
    log(f"[RAW ROWS]             {len(raw_df)}")
    log(f"[FAILED BATCHES]       {len(failed_batches)}")

    return raw_cache_path


# =============================================================================
# Read and clean paralogue table
# =============================================================================

def read_raw_paralogues(raw_cache_path):
    raw_cache_path = Path(raw_cache_path)

    if not raw_cache_path.exists() or raw_cache_path.stat().st_size == 0:
        die(f"Raw BioMart cache not found: {raw_cache_path}")

    line()
    log(f"[READ RAW PARALOGUES] {raw_cache_path}")

    df = pd.read_csv(
        raw_cache_path,
        sep="\t",
        header=None,
        names=RAW_COLUMNS,
        dtype=str,
        keep_default_na=False,
        low_memory=False,
        compression="gzip" if str(raw_cache_path).endswith(".gz") else None,
    )

    for c in RAW_COLUMNS:
        if c not in df.columns:
            df[c] = ""

    df = df[RAW_COLUMNS].copy()

    log(f"[RAW SHAPE] {df.shape}")

    return df


def clean_paralogue_table(
    raw_df,
    symbol_to_approved,
    ensembl_to_approved,
    entrez_to_approved,
    limit_genes=None,
):
    line()
    log("[CLEAN PARALOGUE TABLE]")

    df = raw_df.copy()

    for c in RAW_COLUMNS:
        if c not in df.columns:
            df[c] = ""
        df[c] = df[c].map(clean_text)

    n0 = len(df)
    df = df[
        df["query_ensembl_gene_id"].ne("")
        & (
            df["paralogue_ensembl_gene_id"].ne("")
            | df["paralogue_associated_gene_name"].ne("")
        )
    ].copy()
    log(f"[DROP EMPTY PARALOGUE ROWS] {n0} -> {len(df)}")

    df["approved_symbol"] = df["query_ensembl_gene_id"].map(ensembl_to_approved)
    df["query_map_source"] = np.where(df["approved_symbol"].notna(), "ensembl", "unmapped")

    resolved_para = df.apply(
        lambda r: resolve_gene(
            r.get("paralogue_associated_gene_name", ""),
            r.get("paralogue_ensembl_gene_id", ""),
            "",
            symbol_to_approved,
            ensembl_to_approved,
            entrez_to_approved,
        ),
        axis=1,
        result_type="expand",
    )

    df["paralogue_approved_symbol"] = resolved_para[0]
    df["paralogue_map_source"] = resolved_para[1]

    n1 = len(df)
    df = df[df["approved_symbol"].notna()].copy()
    log(f"[QUERY GENE HGNC MAP] {n1} -> {len(df)}")

    n2 = len(df)
    df = df[
        (df["approved_symbol"].fillna("") != df["paralogue_approved_symbol"].fillna(""))
        | df["paralogue_approved_symbol"].isna()
    ].copy()
    log(f"[DROP SELF PAIRS] {n2} -> {len(df)}")

    df["paralogue_perc_id_target_identical_to_query"] = to_float_series(
        df["paralogue_perc_id_target_identical_to_query"]
    )
    df["paralogue_perc_id_query_identical_to_target"] = to_float_series(
        df["paralogue_perc_id_query_identical_to_target"]
    )

    df["paralogue_identity_mean"] = df[
        [
            "paralogue_perc_id_target_identical_to_query",
            "paralogue_perc_id_query_identical_to_target",
        ]
    ].mean(axis=1, skipna=True)

    df["paralogue_identity_min_symmetric"] = df[
        [
            "paralogue_perc_id_target_identical_to_query",
            "paralogue_perc_id_query_identical_to_target",
        ]
    ].min(axis=1, skipna=True)

    df["paralogue_unique_id"] = df["paralogue_approved_symbol"].fillna("").astype(str)
    mask = df["paralogue_unique_id"].eq("")
    df.loc[mask, "paralogue_unique_id"] = df.loc[mask, "paralogue_ensembl_gene_id"]
    mask = df["paralogue_unique_id"].eq("")
    df.loc[mask, "paralogue_unique_id"] = df.loc[mask, "paralogue_associated_gene_name"]

    n3 = len(df)
    df = df.sort_values(
        ["approved_symbol", "paralogue_unique_id", "paralogue_identity_mean"],
        ascending=[True, True, False],
        na_position="last",
    )
    df = df.drop_duplicates(
        subset=["approved_symbol", "paralogue_unique_id"],
        keep="first",
    ).copy()
    log(f"[DEDUP QUERY-PARALOGUE PAIRS] {n3} -> {len(df)}")

    if limit_genes is not None:
        genes = sorted(df["approved_symbol"].dropna().unique())[: int(limit_genes)]
        n4 = len(df)
        df = df[df["approved_symbol"].isin(genes)].copy()
        log(f"[LIMIT GENES AFTER CLEAN] {n4} -> {len(df)} rows for {len(genes)} genes")

    log(f"[UNIQUE QUERY GENES] {df['approved_symbol'].nunique()}")
    log(f"[UNIQUE PARALOGUES]  {df['paralogue_unique_id'].nunique()}")

    return df


# =============================================================================
# Feature construction
# =============================================================================

def build_gene_features(hgnc_keep, long_df):
    line()
    log("[BUILD GENE FEATURES]")

    feature_index = pd.Index(
        sorted(hgnc_keep["approved_symbol"].dropna().unique()),
        name="approved_symbol",
    )

    base_cols = [
        "ensembl_paralogue_count",
        "ensembl_paralogue_hgnc_mapped_count",
        "ensembl_paralogue_identity_mean",
        "ensembl_paralogue_identity_median",
        "ensembl_paralogue_identity_max",
        "ensembl_paralogue_identity_min",
        "ensembl_paralogue_identity_sd",
        "ensembl_paralogue_identity_min_symmetric_mean",
        "ensembl_paralogue_identity_min_symmetric_max",
        "ensembl_has_any_paralogue",
        "ensembl_has_high_identity_paralogue_ge30",
        "ensembl_has_high_identity_paralogue_ge50",
        "ensembl_has_high_identity_paralogue_ge70",
        "ensembl_has_high_identity_paralogue_ge90",
        "ensembl_high_identity_paralogue_count_ge30",
        "ensembl_high_identity_paralogue_count_ge50",
        "ensembl_high_identity_paralogue_count_ge70",
        "ensembl_high_identity_paralogue_count_ge90",
        "ensembl_paralogue_query_target_identity_mean",
        "ensembl_paralogue_query_target_identity_max",
        "ensembl_paralogue_target_query_identity_mean",
        "ensembl_paralogue_target_query_identity_max",
    ]

    if long_df is None or long_df.empty:
        features = pd.DataFrame({c: 0 for c in base_cols}, index=feature_index)
        return features.reset_index()

    grouped = long_df.groupby("approved_symbol", sort=True)

    features = pd.DataFrame(index=feature_index)

    features["ensembl_paralogue_count"] = grouped["paralogue_unique_id"].nunique()
    features["ensembl_paralogue_hgnc_mapped_count"] = grouped["paralogue_approved_symbol"].agg(safe_nunique)

    features["ensembl_paralogue_identity_mean"] = grouped["paralogue_identity_mean"].mean()
    features["ensembl_paralogue_identity_median"] = grouped["paralogue_identity_mean"].median()
    features["ensembl_paralogue_identity_max"] = grouped["paralogue_identity_mean"].max()
    features["ensembl_paralogue_identity_min"] = grouped["paralogue_identity_mean"].min()
    features["ensembl_paralogue_identity_sd"] = grouped["paralogue_identity_mean"].std(ddof=0)

    features["ensembl_paralogue_identity_min_symmetric_mean"] = grouped[
        "paralogue_identity_min_symmetric"
    ].mean()
    features["ensembl_paralogue_identity_min_symmetric_max"] = grouped[
        "paralogue_identity_min_symmetric"
    ].max()

    features["ensembl_paralogue_query_target_identity_mean"] = grouped[
        "paralogue_perc_id_query_identical_to_target"
    ].mean()
    features["ensembl_paralogue_query_target_identity_max"] = grouped[
        "paralogue_perc_id_query_identical_to_target"
    ].max()

    features["ensembl_paralogue_target_query_identity_mean"] = grouped[
        "paralogue_perc_id_target_identical_to_query"
    ].mean()
    features["ensembl_paralogue_target_query_identity_max"] = grouped[
        "paralogue_perc_id_target_identical_to_query"
    ].max()

    tmp = long_df[
        ["approved_symbol", "paralogue_unique_id", "paralogue_identity_mean"]
    ].drop_duplicates().copy()

    for threshold in [30, 50, 70, 90]:
        flag = f"ge{threshold}"
        tmp[f"flag_{flag}"] = (tmp["paralogue_identity_mean"] >= threshold).astype(int)

        features[f"ensembl_high_identity_paralogue_count_{flag}"] = (
            tmp.groupby("approved_symbol")[f"flag_{flag}"].sum()
        )

        features[f"ensembl_has_high_identity_paralogue_{flag}"] = (
            features[f"ensembl_high_identity_paralogue_count_{flag}"].fillna(0) > 0
        ).astype(int)

    features["ensembl_has_any_paralogue"] = (
        features["ensembl_paralogue_count"].fillna(0) > 0
    ).astype(int)

    features = features.reindex(feature_index).fillna(0)

    ordered = [c for c in base_cols if c in features.columns]
    extra = [c for c in features.columns if c not in ordered]
    features = features[ordered + extra]

    return features.reset_index()


# =============================================================================
# Output helpers
# =============================================================================

def write_summary(summary_path, args, raw_path, long_df, features, query_ensembl_count):
    summary_path = Path(summary_path)
    ensure_dir(summary_path.parent)

    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("FEATURE 22: Ensembl BioMart paralogue features\n")
        f.write("=" * 80 + "\n")
        f.write(f"Script\t{SCRIPT_NAME}\n")
        f.write(f"HGNC\t{args.hgnc}\n")
        f.write(f"DBDIR\t{args.dbdir}\n")
        f.write(f"OUTDIR\t{args.outdir}\n")
        f.write(f"RAW_CACHE\t{raw_path}\n")
        f.write(f"HOST\t{args.host or 'auto mirrors'}\n")
        f.write(f"DOWNLOAD_ALLOWED\t{not args.no_download}\n")
        f.write(f"FORCE_DOWNLOAD\t{args.force_download}\n")
        f.write(f"PROTEIN_CODING_ONLY\t{args.protein_coding_only}\n")
        f.write(f"LIMIT_GENES\t{args.limit_genes if args.limit_genes is not None else 'none'}\n")
        f.write(f"BATCH_SIZE\t{args.batch_size}\n")
        f.write(f"QUERY_ENSEMBL_IDS\t{query_ensembl_count}\n")
        f.write("\n")
        f.write(f"LONG_ROWS\t{0 if long_df is None else len(long_df)}\n")
        f.write(
            f"LONG_QUERY_GENES\t"
            f"{0 if long_df is None or long_df.empty else long_df['approved_symbol'].nunique()}\n"
        )
        f.write(
            f"LONG_PARALOGUES\t"
            f"{0 if long_df is None or long_df.empty else long_df['paralogue_unique_id'].nunique()}\n"
        )
        f.write(f"FEATURE_ROWS\t{len(features)}\n")
        f.write(f"FEATURE_COLUMNS\t{features.shape[1]}\n")
        f.write(
            f"GENES_WITH_PARALOGUE\t"
            f"{int(features['ensembl_has_any_paralogue'].sum()) if 'ensembl_has_any_paralogue' in features.columns else 0}\n"
        )
        f.write(
            f"MEAN_PARALOGUE_COUNT\t"
            f"{features['ensembl_paralogue_count'].mean() if 'ensembl_paralogue_count' in features.columns else 0}\n"
        )
        f.write("\nFEATURE_COLUMNS_LIST\n")
        for c in features.columns:
            f.write(f"{c}\n")

    log(f"[WRITE] {summary_path}")


# =============================================================================
# CLI and main
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Feature 22: Ensembl BioMart paralogue features"
    )

    p.add_argument("--hgnc", default=DEFAULT_HGNC, help="HGNC complete set TSV")
    p.add_argument("--dbdir", default=DEFAULT_DBDIR, help="Database/cache directory")
    p.add_argument("--outdir", default=DEFAULT_OUTDIR, help="Output directory")

    p.add_argument("--host", default=None, help="Specific Ensembl host, e.g. https://www.ensembl.org")
    p.add_argument("--force-download", action="store_true", help="Refresh cached BioMart paralogue table")
    p.add_argument("--no-download", action="store_true", help="Do not query BioMart; require cached raw file")
    p.add_argument("--timeout", type=int, default=300, help="BioMart request timeout in seconds")
    p.add_argument("--retries", type=int, default=2, help="Retries per BioMart host")

    p.add_argument("--batch-size", type=int, default=100, help="Number of Ensembl IDs per BioMart batch")
    p.add_argument("--full-query-first", action="store_true", help="Try one full BioMart query first, then fall back to batch")
    p.add_argument("--limit-genes", type=int, default=None, help="Quick test: first N HGNC genes")
    p.add_argument("--protein-coding-only", action=argparse.BooleanOptionalAction, default=True)

    return p.parse_args()


def main():
    args = parse_args()

    dbdir = ensure_dir(args.dbdir)
    outdir = ensure_dir(args.outdir)
    processed = ensure_dir(outdir / "processed")

    raw_path = dbdir / RAW_CACHE_NAME
    meta_path = dbdir / RAW_META_NAME

    line()
    log("FEATURE 22: ENSEMBL BIOMART HUMAN PARALOGUE FEATURES")
    line()
    log(f"[HGNC]              {args.hgnc}")
    log(f"[DBDIR]             {dbdir.resolve()}")
    log(f"[OUTDIR]            {outdir.resolve()}")
    log(f"[RAW CACHE]         {raw_path}")
    log(f"[HOST]              {args.host or 'auto mirrors'}")
    log(f"[DOWNLOAD ALLOWED]  {not args.no_download}")
    log(f"[FORCE DOWNLOAD]    {args.force_download}")
    log(f"[LIMIT GENES]       {args.limit_genes if args.limit_genes is not None else 'none'}")
    log(f"[PROTEIN CODING]    {args.protein_coding_only}")
    log(f"[BATCH SIZE]        {args.batch_size}")
    log(f"[FULL QUERY FIRST]  {args.full_query_first}")
    line()

    (
        hgnc_keep,
        symbol_to_approved,
        ensembl_to_approved,
        entrez_to_approved,
        approved_to_ensembl,
        query_ensembl_ids,
    ) = read_hgnc(
        args.hgnc,
        protein_coding_only=args.protein_coding_only,
        limit_genes=args.limit_genes,
    )

    if not query_ensembl_ids:
        die("No Ensembl gene IDs found in HGNC file. Cannot query BioMart.")

    if args.force_download or not raw_path.exists() or raw_path.stat().st_size == 0:
        if args.no_download:
            die(
                f"Raw paralogue cache is missing but --no-download was supplied. Expected file: {raw_path}"
            )

        download_biomart_paralogues(
            raw_cache_path=raw_path,
            meta_path=meta_path,
            query_ensembl_ids=query_ensembl_ids,
            host=args.host,
            timeout=args.timeout,
            retries=args.retries,
            batch_size=args.batch_size,
            use_full_query=args.full_query_first,
        )

    else:
        log(f"[CACHE HIT] {raw_path} ({file_size_mb(raw_path):.2f} MB)")

    raw_df = read_raw_paralogues(raw_path)

    long_df = clean_paralogue_table(
        raw_df,
        symbol_to_approved=symbol_to_approved,
        ensembl_to_approved=ensembl_to_approved,
        entrez_to_approved=entrez_to_approved,
        limit_genes=None,
    )

    features = build_gene_features(hgnc_keep, long_df)

    merged = hgnc_keep.merge(features, on="approved_symbol", how="left")

    feature_cols = [c for c in features.columns if c != "approved_symbol"]
    for c in feature_cols:
        merged[c] = pd.to_numeric(merged[c], errors="coerce").fillna(0)

    long_out = processed / "feature22_paralogues_long.csv.gz"
    features_out = processed / "feature22_paralogues_gene_features.csv"
    merged_out = processed / "feature22_paralogues_gene_features_hgnc_merged.csv"
    summary_out = processed / "feature22_paralogues_summary.txt"

    line()
    log("[WRITE OUTPUTS]")

    long_df.to_csv(long_out, index=False, compression="gzip")
    log(f"[WRITE] {long_out} ({file_size_mb(long_out):.2f} MB)")

    features.to_csv(features_out, index=False)
    log(f"[WRITE] {features_out} ({file_size_mb(features_out):.2f} MB)")

    merged.to_csv(merged_out, index=False)
    log(f"[WRITE] {merged_out} ({file_size_mb(merged_out):.2f} MB)")

    write_summary(
        summary_path=summary_out,
        args=args,
        raw_path=raw_path,
        long_df=long_df,
        features=features,
        query_ensembl_count=len(query_ensembl_ids),
    )

    line()
    log("[DONE] Feature 22 complete.")
    log(f"[GENES IN FEATURE MATRIX] {len(features)}")
    log(f"[GENES WITH PARALOGUE]    {int(features['ensembl_has_any_paralogue'].sum())}")
    log(f"[OUTPUT] {features_out}")
    line()


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        die(e)