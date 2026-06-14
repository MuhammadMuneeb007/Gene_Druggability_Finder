# -*- coding: utf-8 -*-

"""
Feature2_String.py

Feature 2: STRING PPI / functional-association graph features for gene-level
druggability modelling.

Why this feature is included
----------------------------
PPI/network topology is useful for druggability prediction because proteins that
occupy central, highly connected, dense, or bridge-like network regions may have
different therapeutic feasibility and biological impact profiles.

Leakage policy
--------------
This script uses only STRING protein-protein / functional-association topology:
  - STRING mapping
  - interaction counts
  - STRING confidence/evidence-channel summaries
  - graph degree/strength
  - local clustering/neighbourhood density
  - k-core
  - connected component features
  - PageRank
  - eigenvector centrality
  - approximate betweenness centrality
  - approximate closeness centrality

It deliberately DOES NOT use:
  - known drug-target neighbours
  - ChEMBL / DrugBank / DGIdb target labels
  - compound-target metadata
  - drug response
  - approved-drug annotations
  - manually curated drug-target families

STRING API endpoints used
-------------------------
  - get_string_ids
  - interaction_partners

The STRING interaction_partners endpoint returns interactions involving the query
proteins and their returned partners. Therefore graph features are computed on
the observed STRING query-partner graph produced by this run, not on the entire
downloaded STRING human interactome.

Recommended run
---------------
    python Feature2_String.py --required-score 700 --partner-limit 100

Faster/default run:
    python Feature2_String.py

Resume without API calls:
    python Feature2_String.py --no-api

If networkx is missing:
    pip install networkx

Inputs
------
Default HGNC input:
    databases/HGNC/hgnc_complete_set.txt

Outputs
-------
    feature2.string.database/raw/string_mapping.tsv
    feature2.string.database/raw/string_interaction_partners.tsv
    feature2.string.database/processed/feature2_string_gene_features.csv
    feature2.string.database/processed/feature2_string_hgnc_merged.csv
    feature2.string.database/feature2_string_summary.txt
    feature2.string.database/feature2_string_run_metadata.json
"""

from __future__ import annotations

import argparse
import io
import json
import math
import random
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import requests

try:
    import networkx as nx
except ImportError as e:
    nx = None


# =============================================================================
# DEFAULTS
# =============================================================================

DEFAULT_HGNC_PATH = Path("databases") / "HGNC" / "hgnc_complete_set.txt"
DEFAULT_OUTDIR = Path("feature2.string.database")

STRING_API_BASE = "https://string-db.org/api"
STRING_SPECIES = 9606
CALLER_IDENTITY = "Feature2_String_gene_druggability_model"

DEFAULT_BATCH_SIZE = 100
DEFAULT_REQUIRED_SCORE = 400
DEFAULT_PARTNER_LIMIT = 50
DEFAULT_SLEEP_SECONDS = 1.0

# For approximate global centrality. Increase for more stable estimates.
DEFAULT_BETWEENNESS_K = 1000
DEFAULT_CLOSENESS_LANDMARKS = 1000
DEFAULT_RANDOM_SEED = 42

SCORE_COLUMNS = [
    "score",
    "nscore",
    "fscore",
    "pscore",
    "ascore",
    "escore",
    "dscore",
    "tscore",
]


# =============================================================================
# HELPERS
# =============================================================================

def now_iso() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str) -> None:
    print(msg, flush=True)


def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def normalize_gene_symbol(x: object) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


def chunks(items: List[str], size: int) -> Iterable[List[str]]:
    for i in range(0, len(items), size):
        yield items[i:i + size]


def scale_score_series(s: pd.Series) -> pd.Series:
    """
    STRING API TSV usually returns scores in 0-1 scale, but many papers/docs
    discuss 0-1000 score thresholds. Convert all score columns to 0-1000.
    """
    s = pd.to_numeric(s, errors="coerce")
    mx = s.max(skipna=True)
    if pd.notna(mx) and mx <= 1.0:
        return s * 1000.0
    return s


def clean_string_id(x: object) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip()


def safe_float(x: object) -> float:
    try:
        if pd.isna(x):
            return np.nan
        return float(x)
    except Exception:
        return np.nan


def require_networkx() -> None:
    if nx is None:
        raise ImportError(
            "networkx is required for graph features.\n"
            "Install it with:\n"
            "    pip install networkx\n"
            "or in conda:\n"
            "    conda install -c conda-forge networkx"
        )


def post_string_tsv(
    endpoint: str,
    params: Dict[str, object],
    retries: int = 5,
    timeout: int = 180,
    sleep_seconds: float = 1.0,
) -> pd.DataFrame:
    """
    POST to a STRING TSV endpoint and return a DataFrame.
    """
    url = f"{STRING_API_BASE}/tsv/{endpoint}"

    last_error = None
    for attempt in range(1, retries + 1):
        try:
            r = requests.post(url, data=params, timeout=timeout)
            if r.status_code == 200:
                text = r.text.strip()
                if not text:
                    return pd.DataFrame()
                return pd.read_csv(io.StringIO(text), sep="\t", dtype=str)

            last_error = RuntimeError(f"HTTP {r.status_code}: {r.text[:500]}")
            log(f"[STRING WARNING] {endpoint} attempt {attempt}/{retries}: {last_error}")

        except Exception as e:
            last_error = e
            log(f"[STRING WARNING] {endpoint} attempt {attempt}/{retries}: {e}")

        if attempt < retries:
            time.sleep(sleep_seconds * attempt)

    raise RuntimeError(f"STRING API failed for endpoint={endpoint}: {last_error}")


# =============================================================================
# HGNC INPUT
# =============================================================================

def load_hgnc(hgnc_path: Path, protein_coding_only: bool = True) -> pd.DataFrame:
    if not hgnc_path.exists():
        raise FileNotFoundError(
            f"HGNC file not found: {hgnc_path}\n"
            "Expected default: databases/HGNC/hgnc_complete_set.txt"
        )

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

    if "entrez_id" in hgnc.columns:
        hgnc["entrez_id_str"] = hgnc["entrez_id"].astype(str).str.strip().replace({"nan": "", "None": ""})
    else:
        hgnc["entrez_id_str"] = ""

    if "ensembl_gene_id" in hgnc.columns:
        hgnc["ensembl_gene_id_str"] = hgnc["ensembl_gene_id"].astype(str).str.strip().replace({"nan": "", "None": ""})
    else:
        hgnc["ensembl_gene_id_str"] = ""

    hgnc = hgnc.drop_duplicates(subset=["gene_symbol"], keep="first").copy()

    log(f"[HGNC ROWS] {before} -> {len(hgnc)}")
    log(f"[PROTEIN CODING ONLY] {protein_coding_only}")

    return hgnc


# =============================================================================
# STRING API
# =============================================================================

def map_genes_to_string(
    genes: List[str],
    raw_dir: Path,
    batch_size: int,
    sleep_seconds: float,
    force: bool = False,
) -> pd.DataFrame:
    """
    Map HGNC symbols to STRING identifiers using get_string_ids.
    """
    outpath = raw_dir / "string_mapping.tsv"

    if outpath.exists() and outpath.stat().st_size > 0 and not force:
        log(f"[SKIP] Existing STRING mapping: {outpath}")
        return pd.read_csv(outpath, sep="\t", dtype=str)

    all_rows = []
    batches = list(chunks(genes, batch_size))

    log("=" * 100)
    log("[STRING MAPPING]")
    log(f"[GENES] {len(genes)}")
    log(f"[BATCH SIZE] {batch_size}")
    log(f"[BATCHES] {len(batches)}")

    for idx, batch in enumerate(batches, start=1):
        log(f"[MAP] Batch {idx}/{len(batches)} genes={len(batch)}")

        params = {
            "identifiers": "\r".join(batch),
            "species": STRING_SPECIES,
            "limit": 1,
            "echo_query": 1,
            "caller_identity": CALLER_IDENTITY,
        }

        df = post_string_tsv(
            endpoint="get_string_ids",
            params=params,
            retries=5,
            timeout=120,
            sleep_seconds=sleep_seconds,
        )

        if not df.empty:
            all_rows.append(df)

        time.sleep(sleep_seconds)

    if all_rows:
        mapping = pd.concat(all_rows, ignore_index=True)
    else:
        mapping = pd.DataFrame()

    if "queryItem" in mapping.columns:
        mapping["gene_symbol"] = mapping["queryItem"].map(normalize_gene_symbol)
    elif "preferredName" in mapping.columns:
        mapping["gene_symbol"] = mapping["preferredName"].map(normalize_gene_symbol)
    else:
        mapping["gene_symbol"] = ""

    mapping = mapping.drop_duplicates(subset=["gene_symbol"], keep="first").copy()
    mapping.to_csv(outpath, sep="\t", index=False)

    log(f"[SAVED MAPPING] {outpath}")
    log(f"[MAPPED GENES] {mapping['gene_symbol'].nunique() if 'gene_symbol' in mapping.columns else len(mapping)}")

    return mapping


def fetch_interaction_partners(
    mapping: pd.DataFrame,
    raw_dir: Path,
    batch_size: int,
    required_score: int,
    partner_limit: int,
    sleep_seconds: float,
    force: bool = False,
) -> pd.DataFrame:
    """
    Retrieve interaction partners for mapped STRING proteins.
    """
    outpath = raw_dir / "string_interaction_partners.tsv"

    if outpath.exists() and outpath.stat().st_size > 0 and not force:
        log(f"[SKIP] Existing STRING interaction partners: {outpath}")
        return pd.read_csv(outpath, sep="\t", dtype=str)

    if mapping.empty:
        raise RuntimeError("STRING mapping is empty; cannot fetch interaction partners.")

    string_col = None
    for c in ["stringId", "string_id", "STRING_ID"]:
        if c in mapping.columns:
            string_col = c
            break

    if string_col is None:
        raise RuntimeError(f"Could not find STRING ID column in mapping. Columns: {mapping.columns.tolist()}")

    ids = mapping[string_col].dropna().astype(str).unique().tolist()
    batches = list(chunks(ids, batch_size))

    all_rows = []

    log("=" * 100)
    log("[STRING INTERACTION PARTNERS]")
    log(f"[STRING IDS] {len(ids)}")
    log(f"[BATCH SIZE] {batch_size}")
    log(f"[BATCHES] {len(batches)}")
    log(f"[REQUIRED SCORE] {required_score}")
    log(f"[PARTNER LIMIT] {partner_limit}")

    for idx, batch in enumerate(batches, start=1):
        log(f"[INTERACTIONS] Batch {idx}/{len(batches)} proteins={len(batch)}")

        params = {
            "identifiers": "\r".join(batch),
            "species": STRING_SPECIES,
            "required_score": required_score,
            "limit": partner_limit,
            "caller_identity": CALLER_IDENTITY,
        }

        df = post_string_tsv(
            endpoint="interaction_partners",
            params=params,
            retries=5,
            timeout=240,
            sleep_seconds=sleep_seconds,
        )

        if not df.empty:
            df["feature2_query_batch"] = idx
            all_rows.append(df)

        time.sleep(sleep_seconds)

    if all_rows:
        interactions = pd.concat(all_rows, ignore_index=True)
    else:
        interactions = pd.DataFrame()

    interactions.to_csv(outpath, sep="\t", index=False)

    log(f"[SAVED INTERACTIONS] {outpath}")
    log(f"[INTERACTION ROWS] {len(interactions)}")

    return interactions


# =============================================================================
# INTERACTION CLEANING
# =============================================================================

def choose_string_columns(interactions: pd.DataFrame) -> Tuple[str, str, str, str]:
    """
    STRING interaction_partners usually returns:
      stringId_A, preferredName_A, stringId_B, preferredName_B
    """
    needed = ["stringId_A", "preferredName_A", "stringId_B", "preferredName_B"]
    if all(c in interactions.columns for c in needed):
        return "stringId_A", "preferredName_A", "stringId_B", "preferredName_B"

    cols = {c.lower(): c for c in interactions.columns}
    a_id = cols.get("stringid_a") or cols.get("string_id_a")
    a_name = cols.get("preferredname_a") or cols.get("preferred_name_a")
    b_id = cols.get("stringid_b") or cols.get("string_id_b")
    b_name = cols.get("preferredname_b") or cols.get("preferred_name_b")

    if not all([a_id, a_name, b_id, b_name]):
        raise RuntimeError(f"Could not identify STRING A/B columns. Columns: {interactions.columns.tolist()}")

    return a_id, a_name, b_id, b_name


def prepare_mapping_features(hgnc: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    genes = hgnc[["gene_symbol"]].drop_duplicates().copy()

    base = genes.copy()
    base["feature2_string_mapped"] = 0
    base["feature2_string_id"] = ""
    base["feature2_string_preferred_name"] = ""
    base["feature2_string_annotation"] = ""

    if mapping.empty:
        return base

    m = mapping.copy()

    if "queryItem" in m.columns:
        m["gene_symbol"] = m["queryItem"].map(normalize_gene_symbol)
    elif "preferredName" in m.columns:
        m["gene_symbol"] = m["preferredName"].map(normalize_gene_symbol)
    else:
        m["gene_symbol"] = ""

    m = m.drop_duplicates(subset=["gene_symbol"], keep="first").copy()

    keep = ["gene_symbol"]
    rename = {}

    if "stringId" in m.columns:
        keep.append("stringId")
        rename["stringId"] = "feature2_string_id"
    if "preferredName" in m.columns:
        keep.append("preferredName")
        rename["preferredName"] = "feature2_string_preferred_name"
    if "annotation" in m.columns:
        keep.append("annotation")
        rename["annotation"] = "feature2_string_annotation"

    m_small = m[keep].rename(columns=rename)
    m_small["feature2_string_mapped"] = 1

    out = genes.merge(m_small, on="gene_symbol", how="left")
    out["feature2_string_mapped"] = out["feature2_string_mapped"].fillna(0).astype(int)

    for c in ["feature2_string_id", "feature2_string_preferred_name", "feature2_string_annotation"]:
        if c not in out.columns:
            out[c] = ""
        out[c] = out[c].fillna("")

    return out


def clean_interactions(
    interactions: pd.DataFrame,
    mapping_features: pd.DataFrame,
) -> pd.DataFrame:
    """
    Clean STRING interactions and attach HGNC query gene symbols.
    """
    if interactions.empty:
        return pd.DataFrame()

    a_id, a_name, b_id, b_name = choose_string_columns(interactions)
    x = interactions.copy()

    x["string_id_a"] = x[a_id].map(clean_string_id)
    x["string_id_b"] = x[b_id].map(clean_string_id)
    x["gene_symbol_a_raw"] = x[a_name].map(normalize_gene_symbol)
    x["gene_symbol_b_raw"] = x[b_name].map(normalize_gene_symbol)

    # Recover original query HGNC gene symbol from STRING ID.
    id_to_hgnc = dict(
        zip(
            mapping_features["feature2_string_id"].astype(str),
            mapping_features["gene_symbol"].astype(str),
        )
    )

    x["gene_symbol_a"] = x["string_id_a"].map(id_to_hgnc)
    x["gene_symbol_a"] = x["gene_symbol_a"].fillna(x["gene_symbol_a_raw"]).map(normalize_gene_symbol)
    x["partner_symbol"] = x["gene_symbol_b_raw"].map(normalize_gene_symbol)

    # Score conversion.
    for c in SCORE_COLUMNS:
        if c in x.columns:
            x[f"{c}_1000"] = scale_score_series(x[c])
        else:
            x[f"{c}_1000"] = np.nan

    x["score_1000"] = pd.to_numeric(x["score_1000"], errors="coerce")

    # Remove malformed rows.
    x = x[(x["gene_symbol_a"] != "") & (x["partner_symbol"] != "")].copy()

    return x


# =============================================================================
# LOCAL STRING SUMMARY FEATURES
# =============================================================================

def build_local_interaction_features(
    hgnc: pd.DataFrame,
    mapping_features: pd.DataFrame,
    clean_edges: pd.DataFrame,
    required_score: int,
    partner_limit: int,
) -> pd.DataFrame:
    """
    Direct per-query STRING summary features.
    """
    out = mapping_features.copy()

    if clean_edges.empty:
        out["feature2_string_interaction_count"] = 0
        out["feature2_string_unique_partner_count"] = 0
        out["feature2_string_high_confidence_count_ge_700"] = 0
        out["feature2_string_very_high_confidence_count_ge_900"] = 0
        out["feature2_string_mean_score"] = np.nan
        out["feature2_string_median_score"] = np.nan
        out["feature2_string_max_score"] = np.nan
        out["feature2_string_q75_score"] = np.nan
        out["feature2_string_top_partners"] = ""
        return out

    records = []

    for gene, g in clean_edges.groupby("gene_symbol_a"):
        score = pd.to_numeric(g["score_1000"], errors="coerce")

        rec = {
            "gene_symbol": gene,
            "feature2_string_interaction_count": int(len(g)),
            "feature2_string_unique_partner_count": int(g["partner_symbol"].nunique()),
            "feature2_string_high_confidence_count_ge_700": int((score >= 700).sum(skipna=True)),
            "feature2_string_very_high_confidence_count_ge_900": int((score >= 900).sum(skipna=True)),
            "feature2_string_mean_score": float(score.mean(skipna=True)) if score.notna().any() else np.nan,
            "feature2_string_median_score": float(score.median(skipna=True)) if score.notna().any() else np.nan,
            "feature2_string_max_score": float(score.max(skipna=True)) if score.notna().any() else np.nan,
            "feature2_string_q75_score": float(score.quantile(0.75)) if score.notna().any() else np.nan,
            "feature2_string_required_score": required_score,
            "feature2_string_partner_limit": partner_limit,
        }

        for c in ["nscore", "fscore", "pscore", "ascore", "escore", "dscore", "tscore"]:
            sc = f"{c}_1000"
            vals = pd.to_numeric(g[sc], errors="coerce")
            rec[f"feature2_string_mean_{c}"] = float(vals.mean(skipna=True)) if vals.notna().any() else np.nan
            rec[f"feature2_string_max_{c}"] = float(vals.max(skipna=True)) if vals.notna().any() else np.nan

        top = g.copy()
        top["score_sort"] = pd.to_numeric(top["score_1000"], errors="coerce")
        top = top.sort_values("score_sort", ascending=False)

        partners = top["partner_symbol"].dropna().astype(str).tolist()
        partners = [p for p in partners if p and p != gene]
        rec["feature2_string_top_partners"] = ";".join(partners[:10])

        records.append(rec)

    local = pd.DataFrame(records)
    out = out.merge(local, on="gene_symbol", how="left")

    zero_cols = [
        "feature2_string_interaction_count",
        "feature2_string_unique_partner_count",
        "feature2_string_high_confidence_count_ge_700",
        "feature2_string_very_high_confidence_count_ge_900",
    ]
    for c in zero_cols:
        out[c] = out[c].fillna(0).astype(int)

    if "feature2_string_top_partners" in out.columns:
        out["feature2_string_top_partners"] = out["feature2_string_top_partners"].fillna("")

    return out


# =============================================================================
# GRAPH CONSTRUCTION AND FEATURES
# =============================================================================

def build_observed_string_graph(clean_edges: pd.DataFrame, hgnc_genes: Set[str]) -> "nx.Graph":
    """
    Build an undirected graph from observed STRING query-partner edges.

    Nodes are gene symbols. Edges are weighted by STRING combined score.
    If duplicate edges occur, the highest observed score is retained.
    """
    require_networkx()

    G = nx.Graph()

    # Include all HGNC genes as nodes, even if isolated.
    for gene in hgnc_genes:
        G.add_node(gene, is_hgnc_query=True)

    if clean_edges.empty:
        return G

    for _, row in clean_edges.iterrows():
        a = normalize_gene_symbol(row.get("gene_symbol_a"))
        b = normalize_gene_symbol(row.get("partner_symbol"))
        if not a or not b or a == b:
            continue

        score = safe_float(row.get("score_1000"))
        if pd.isna(score):
            score = 0.0

        # NetworkX shortest paths treat smaller weights as shorter. Convert
        # confidence score into a distance. Higher STRING score = shorter edge.
        distance = 1.0 / max(score, 1.0)

        if G.has_edge(a, b):
            old_score = G[a][b].get("weight", 0.0)
            if score > old_score:
                G[a][b]["weight"] = float(score)
                G[a][b]["distance"] = float(distance)
        else:
            G.add_edge(a, b, weight=float(score), distance=float(distance))

    return G


def approximate_closeness_from_landmarks(
    G: "nx.Graph",
    hgnc_nodes: List[str],
    landmarks: int,
    seed: int,
) -> Dict[str, float]:
    """
    Approximate closeness centrality for many nodes using sampled landmarks.

    Exact closeness for all nodes can be slow on large graphs. This approximation
    samples landmark source nodes, accumulates shortest path lengths from those
    landmarks, and estimates inverse average distance.

    Returned values are comparable within this run but should be labelled as
    approximate.
    """
    if G.number_of_nodes() == 0:
        return {n: np.nan for n in hgnc_nodes}

    rng = random.Random(seed)
    all_nodes = list(G.nodes())
    if not all_nodes:
        return {n: np.nan for n in hgnc_nodes}

    k = min(max(1, landmarks), len(all_nodes))
    sampled = rng.sample(all_nodes, k)

    dist_sum = defaultdict(float)
    reach_count = defaultdict(int)

    log(f"[GRAPH] Approx closeness using {k} landmarks")

    for i, src in enumerate(sampled, start=1):
        lengths = nx.single_source_shortest_path_length(G, src)
        for node, d in lengths.items():
            if d > 0:
                dist_sum[node] += float(d)
                reach_count[node] += 1

        if i % 100 == 0 or i == k:
            log(f"  closeness landmarks processed {i}/{k}")

    out = {}
    for node in hgnc_nodes:
        rc = reach_count.get(node, 0)
        if rc == 0:
            out[node] = 0.0
        else:
            avg_d = dist_sum[node] / rc
            reachable_fraction = rc / k
            out[node] = float(reachable_fraction / avg_d) if avg_d > 0 else 0.0

    return out


def calculate_graph_features(
    G: "nx.Graph",
    hgnc_genes: List[str],
    betweenness_k: int,
    closeness_landmarks: int,
    seed: int,
    skip_expensive: bool = False,
) -> pd.DataFrame:
    """
    Calculate graph-theory features for HGNC genes.
    """
    require_networkx()

    log("=" * 100)
    log("[GRAPH FEATURES]")
    log(f"[NODES] {G.number_of_nodes()}")
    log(f"[EDGES] {G.number_of_edges()}")
    log(f"[HGNC QUERY GENES] {len(hgnc_genes)}")

    hgnc_nodes = [g for g in hgnc_genes if g in G]
    n_total = G.number_of_nodes()
    m_total = G.number_of_edges()

    # Basic degree and strength.
    log("[GRAPH] Degree and weighted degree")
    degree = dict(G.degree())
    weighted_degree = dict(G.degree(weight="weight"))

    # Average neighbour degree.
    log("[GRAPH] Average neighbour degree")
    avg_neighbor_degree = nx.average_neighbor_degree(G) if G.number_of_edges() else {}

    log("[GRAPH] Weighted average neighbour degree")
    try:
        avg_neighbor_degree_weighted = nx.average_neighbor_degree(G, weight="weight") if G.number_of_edges() else {}
    except Exception:
        avg_neighbor_degree_weighted = {}

    # Clustering and triangles.
    log("[GRAPH] Clustering coefficients")
    clustering = nx.clustering(G) if G.number_of_edges() else {}
    try:
        weighted_clustering = nx.clustering(G, weight="weight") if G.number_of_edges() else {}
    except Exception:
        weighted_clustering = {}

    log("[GRAPH] Triangle counts")
    triangles = nx.triangles(G) if G.number_of_edges() else {}

    # k-core.
    log("[GRAPH] Core number")
    try:
        core_number = nx.core_number(G) if G.number_of_edges() else {}
    except Exception as e:
        log(f"[WARNING] core_number failed: {e}")
        core_number = {}

    # Connected components.
    log("[GRAPH] Connected components")
    component_size = {}
    component_rank = {}
    component_id = {}
    is_largest_component = {}

    components = sorted(nx.connected_components(G), key=len, reverse=True)
    for rank, comp in enumerate(components, start=1):
        size = len(comp)
        cid = rank
        for node in comp:
            component_size[node] = size
            component_rank[node] = rank
            component_id[node] = cid
            is_largest_component[node] = 1 if rank == 1 else 0

    # PageRank is global but scalable.
    log("[GRAPH] PageRank")
    try:
        pagerank = nx.pagerank(G, weight="weight", max_iter=100, tol=1e-06) if G.number_of_edges() else {}
    except Exception as e:
        log(f"[WARNING] PageRank failed: {e}")
        pagerank = {}

    # Eigenvector centrality can fail if convergence is difficult.
    log("[GRAPH] Eigenvector centrality")
    try:
        eigenvector = nx.eigenvector_centrality(G, max_iter=500, tol=1e-06, weight="weight") if G.number_of_edges() else {}
    except Exception as e:
        log(f"[WARNING] Weighted eigenvector failed: {e}")
        try:
            eigenvector = nx.eigenvector_centrality(G, max_iter=500, tol=1e-06) if G.number_of_edges() else {}
        except Exception as e2:
            log(f"[WARNING] Unweighted eigenvector failed: {e2}")
            eigenvector = {}

    if skip_expensive:
        log("[GRAPH] Skipping approximate betweenness/closeness because --skip-expensive-graph was used")
        betweenness = {}
        closeness_approx = {}
    else:
        # Approximate betweenness.
        log("[GRAPH] Approximate betweenness centrality")
        k = min(max(1, betweenness_k), G.number_of_nodes()) if G.number_of_nodes() else 0
        try:
            betweenness = nx.betweenness_centrality(
                G,
                k=k,
                normalized=True,
                weight="distance",
                seed=seed,
            ) if G.number_of_edges() and k > 0 else {}
        except Exception as e:
            log(f"[WARNING] Approximate betweenness failed: {e}")
            betweenness = {}

        # Approximate closeness.
        log("[GRAPH] Approximate closeness centrality")
        try:
            closeness_approx = approximate_closeness_from_landmarks(
                G=G,
                hgnc_nodes=hgnc_nodes,
                landmarks=closeness_landmarks,
                seed=seed,
            )
        except Exception as e:
            log(f"[WARNING] Approximate closeness failed: {e}")
            closeness_approx = {}

    # Local neighbourhood features: number of edges among neighbours and density.
    log("[GRAPH] Neighbourhood density features")
    records = []

    for i, gene in enumerate(hgnc_nodes, start=1):
        neighbors = list(G.neighbors(gene))
        k_deg = len(neighbors)

        if k_deg >= 2:
            sub = G.subgraph(neighbors)
            neigh_edges = sub.number_of_edges()
            neigh_possible = k_deg * (k_deg - 1) / 2.0
            neigh_density = neigh_edges / neigh_possible if neigh_possible > 0 else 0.0
        else:
            neigh_edges = 0
            neigh_density = 0.0

        deg = degree.get(gene, 0)
        wdeg = weighted_degree.get(gene, 0.0)
        comp_size = component_size.get(gene, 1)

        rec = {
            "gene_symbol": gene,

            # Graph size context, same value for all rows but useful for provenance.
            "feature2_string_graph_total_nodes": n_total,
            "feature2_string_graph_total_edges": m_total,

            # Degree / connectivity.
            "feature2_string_graph_degree": int(deg),
            "feature2_string_graph_weighted_degree": float(wdeg),
            "feature2_string_graph_degree_fraction": float(deg / (n_total - 1)) if n_total > 1 else 0.0,
            "feature2_string_graph_mean_edge_weight": float(wdeg / deg) if deg > 0 else np.nan,

            # Neighbourhood structure.
            "feature2_string_graph_average_neighbor_degree": float(avg_neighbor_degree.get(gene, 0.0)),
            "feature2_string_graph_weighted_average_neighbor_degree": float(avg_neighbor_degree_weighted.get(gene, 0.0)),
            "feature2_string_graph_clustering_coefficient": float(clustering.get(gene, 0.0)),
            "feature2_string_graph_weighted_clustering_coefficient": float(weighted_clustering.get(gene, 0.0)),
            "feature2_string_graph_triangle_count": int(triangles.get(gene, 0)),
            "feature2_string_graph_neighborhood_edge_count": int(neigh_edges),
            "feature2_string_graph_neighborhood_density": float(neigh_density),

            # Core / component.
            "feature2_string_graph_k_core_number": int(core_number.get(gene, 0)),
            "feature2_string_graph_component_size": int(comp_size),
            "feature2_string_graph_component_fraction": float(comp_size / n_total) if n_total > 0 else 0.0,
            "feature2_string_graph_component_rank": int(component_rank.get(gene, 0)),
            "feature2_string_graph_is_largest_component": int(is_largest_component.get(gene, 0)),

            # Global centrality.
            "feature2_string_graph_pagerank": float(pagerank.get(gene, 0.0)),
            "feature2_string_graph_eigenvector_centrality": float(eigenvector.get(gene, 0.0)),
            "feature2_string_graph_betweenness_centrality_approx": float(betweenness.get(gene, 0.0)),
            "feature2_string_graph_closeness_centrality_approx": float(closeness_approx.get(gene, 0.0)),
        }

        records.append(rec)

        if i % 5000 == 0:
            log(f"  graph features assembled {i}/{len(hgnc_nodes)}")

    out = pd.DataFrame(records)

    # Ensure all HGNC genes are represented, even if isolated/missing.
    all_genes_df = pd.DataFrame({"gene_symbol": hgnc_genes})
    out = all_genes_df.merge(out, on="gene_symbol", how="left")

    zero_cols = [
        "feature2_string_graph_total_nodes",
        "feature2_string_graph_total_edges",
        "feature2_string_graph_degree",
        "feature2_string_graph_weighted_degree",
        "feature2_string_graph_degree_fraction",
        "feature2_string_graph_average_neighbor_degree",
        "feature2_string_graph_weighted_average_neighbor_degree",
        "feature2_string_graph_clustering_coefficient",
        "feature2_string_graph_weighted_clustering_coefficient",
        "feature2_string_graph_triangle_count",
        "feature2_string_graph_neighborhood_edge_count",
        "feature2_string_graph_neighborhood_density",
        "feature2_string_graph_k_core_number",
        "feature2_string_graph_component_size",
        "feature2_string_graph_component_fraction",
        "feature2_string_graph_component_rank",
        "feature2_string_graph_is_largest_component",
        "feature2_string_graph_pagerank",
        "feature2_string_graph_eigenvector_centrality",
        "feature2_string_graph_betweenness_centrality_approx",
        "feature2_string_graph_closeness_centrality_approx",
    ]

    for c in zero_cols:
        if c in out.columns:
            out[c] = out[c].fillna(0)

    return out


# =============================================================================
# FINAL FEATURE TABLE
# =============================================================================

def build_all_features(
    hgnc: pd.DataFrame,
    mapping: pd.DataFrame,
    interactions: pd.DataFrame,
    processed_dir: Path,
    required_score: int,
    partner_limit: int,
    betweenness_k: int,
    closeness_landmarks: int,
    seed: int,
    skip_expensive_graph: bool,
) -> pd.DataFrame:
    mkdir(processed_dir)

    hgnc_genes = sorted(hgnc["gene_symbol"].dropna().astype(str).unique().tolist())
    hgnc_gene_set = set(hgnc_genes)

    mapping_features = prepare_mapping_features(hgnc, mapping)
    clean_edges = clean_interactions(interactions, mapping_features)

    clean_edges_path = processed_dir / "feature2_string_clean_edges.csv"
    clean_edges.to_csv(clean_edges_path, index=False)
    log(f"[SAVED CLEAN EDGES] {clean_edges_path}")
    log(f"[CLEAN EDGE ROWS] {len(clean_edges)}")

    local_features = build_local_interaction_features(
        hgnc=hgnc,
        mapping_features=mapping_features,
        clean_edges=clean_edges,
        required_score=required_score,
        partner_limit=partner_limit,
    )

    G = build_observed_string_graph(clean_edges, hgnc_gene_set)

    graph_features = calculate_graph_features(
        G=G,
        hgnc_genes=hgnc_genes,
        betweenness_k=betweenness_k,
        closeness_landmarks=closeness_landmarks,
        seed=seed,
        skip_expensive=skip_expensive_graph,
    )

    # Merge local + graph.
    features = local_features.merge(graph_features, on="gene_symbol", how="left")

    # Flag.
    features["feature2_string_has_ppi_feature"] = (
        (features["feature2_string_mapped"].fillna(0).astype(int) == 1) |
        (features["feature2_string_interaction_count"].fillna(0).astype(float) > 0) |
        (features["feature2_string_graph_degree"].fillna(0).astype(float) > 0)
    ).astype(int)

    outpath = processed_dir / "feature2_string_gene_features.csv"
    features.to_csv(outpath, index=False)

    log("=" * 100)
    log("[SAVED STRING GENE FEATURES]")
    log(f"[PATH]  {outpath}")
    log(f"[SHAPE] {features.shape[0]} genes x {features.shape[1]} columns")

    return features


def merge_with_hgnc(hgnc: pd.DataFrame, features: pd.DataFrame, processed_dir: Path) -> pd.DataFrame:
    log("=" * 100)
    log("[MERGE HGNC + STRING FEATURES]")

    merged = hgnc.merge(features, on="gene_symbol", how="left")

    feature_cols = [c for c in features.columns if c != "gene_symbol"]
    merged["feature2_string_has_any_feature"] = merged[feature_cols].notna().any(axis=1).astype(int)
    merged["feature2_string_n_nonmissing_features"] = merged[feature_cols].notna().sum(axis=1)

    outpath = processed_dir / "feature2_string_hgnc_merged.csv"
    merged.to_csv(outpath, index=False)

    coverage = merged["feature2_string_has_any_feature"].mean() if len(merged) else np.nan
    mapped = merged["feature2_string_mapped"].mean() if "feature2_string_mapped" in merged.columns else np.nan
    has_ppi = merged["feature2_string_has_ppi_feature"].mean() if "feature2_string_has_ppi_feature" in merged.columns else np.nan

    log(f"[SAVED] {outpath}")
    log(f"[SHAPE] {merged.shape[0]} genes x {merged.shape[1]} columns")
    log(f"[COVERAGE ANY FEATURE] {coverage:.3f}")
    log(f"[STRING MAPPED] {mapped:.3f}")
    log(f"[HAS PPI FEATURE] {has_ppi:.3f}")

    return merged


def write_summary_report(
    outdir: Path,
    hgnc: pd.DataFrame,
    features: pd.DataFrame,
    merged: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    report_path = outdir / "feature2_string_summary.txt"

    lines = []
    lines.append("Feature 2 STRING PPI / graph network features")
    lines.append("=" * 80)
    lines.append(f"Created at: {now_iso()}")
    lines.append(f"STRING API base: {STRING_API_BASE}")
    lines.append(f"STRING species: {STRING_SPECIES} human")
    lines.append(f"Required score: {args.required_score}")
    lines.append(f"Partner limit per query protein: {args.partner_limit}")
    lines.append(f"Batch size: {args.batch_size}")
    lines.append(f"Betweenness k: {args.betweenness_k}")
    lines.append(f"Closeness landmarks: {args.closeness_landmarks}")
    lines.append(f"Skip expensive graph features: {args.skip_expensive_graph}")
    lines.append(f"Output directory: {outdir.resolve()}")
    lines.append("")
    lines.append("Leakage rule:")
    lines.append("  Included: STRING PPI / functional-association topology and evidence-channel summaries.")
    lines.append("  Excluded: known drug-target proximity, compound-target metadata, ChEMBL/DrugBank labels, drug-response data.")
    lines.append("")
    lines.append("Important interpretation note:")
    lines.append("  Graph features are computed on the observed STRING query-partner graph generated by this run,")
    lines.append("  not on the full downloaded human STRING interactome.")
    lines.append("")
    lines.append(f"HGNC rows: {hgnc.shape[0]}")
    lines.append(f"Feature rows: {features.shape[0]}")
    lines.append(f"Feature columns: {features.shape[1]}")
    lines.append(f"Merged rows: {merged.shape[0]}")
    lines.append(f"Merged columns: {merged.shape[1]}")

    if "feature2_string_mapped" in merged.columns:
        lines.append(f"STRING mapped fraction: {merged['feature2_string_mapped'].mean():.4f}")
    if "feature2_string_interaction_count" in merged.columns:
        lines.append(f"Genes with >=1 returned partner: {(merged['feature2_string_interaction_count'] > 0).mean():.4f}")
        lines.append(f"Median returned partner count: {merged['feature2_string_interaction_count'].median():.2f}")
        lines.append(f"Mean returned partner count: {merged['feature2_string_interaction_count'].mean():.2f}")
    if "feature2_string_graph_degree" in merged.columns:
        lines.append(f"Genes with graph degree > 0: {(merged['feature2_string_graph_degree'] > 0).mean():.4f}")
        lines.append(f"Median graph degree: {merged['feature2_string_graph_degree'].median():.2f}")
        lines.append(f"Mean graph degree: {merged['feature2_string_graph_degree'].mean():.2f}")
    if "feature2_string_has_any_feature" in merged.columns:
        lines.append(f"Any feature coverage: {merged['feature2_string_has_any_feature'].mean():.4f}")

    report_path.write_text("\n".join(lines) + "\n")
    log(f"[REPORT] {report_path}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build STRING PPI graph features for HGNC protein-coding genes."
    )
    parser.add_argument("--hgnc", default=str(DEFAULT_HGNC_PATH), help="Path to hgnc_complete_set.txt.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory.")
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE, help="Genes/STRING IDs per API batch.")
    parser.add_argument("--required-score", type=int, default=DEFAULT_REQUIRED_SCORE, help="STRING required_score, 0-1000.")
    parser.add_argument("--partner-limit", type=int, default=DEFAULT_PARTNER_LIMIT, help="Max partners returned per query protein.")
    parser.add_argument("--sleep", type=float, default=DEFAULT_SLEEP_SECONDS, help="Seconds to sleep between API calls.")
    parser.add_argument("--all-hgnc-genes", action="store_true", help="Use all HGNC rows instead of protein-coding only.")
    parser.add_argument("--no-api", action="store_true", help="Do not call STRING API; use existing raw TSV files.")
    parser.add_argument("--force-api", action="store_true", help="Re-run STRING API calls even if raw files exist.")
    parser.add_argument("--mapping-only", action="store_true", help="Only map genes to STRING IDs; do not fetch interactions.")
    parser.add_argument("--betweenness-k", type=int, default=DEFAULT_BETWEENNESS_K, help="Sample size for approximate betweenness.")
    parser.add_argument("--closeness-landmarks", type=int, default=DEFAULT_CLOSENESS_LANDMARKS, help="Landmarks for approximate closeness.")
    parser.add_argument("--seed", type=int, default=DEFAULT_RANDOM_SEED, help="Random seed for approximate centrality.")
    parser.add_argument("--skip-expensive-graph", action="store_true", help="Skip approximate betweenness and closeness.")
    args = parser.parse_args()

    require_networkx()

    outdir = mkdir(Path(args.outdir))
    raw_dir = mkdir(outdir / "raw")
    processed_dir = mkdir(outdir / "processed")

    log("=" * 100)
    log("FEATURE 2: STRING PPI / GRAPH NETWORK FEATURES")
    log("=" * 100)
    log(f"[HGNC]                 {args.hgnc}")
    log(f"[OUTDIR]               {outdir.resolve()}")
    log(f"[RAW]                  {raw_dir.resolve()}")
    log(f"[PROCESSED]            {processed_dir.resolve()}")
    log(f"[BATCH SIZE]           {args.batch_size}")
    log(f"[REQUIRED SCORE]       {args.required_score}")
    log(f"[PARTNER LIMIT]        {args.partner_limit}")
    log(f"[BETWEENNESS K]        {args.betweenness_k}")
    log(f"[CLOSENESS LANDMARKS]  {args.closeness_landmarks}")
    log(f"[NO API]               {args.no_api}")
    log(f"[FORCE API]            {args.force_api}")
    log("=" * 100)

    hgnc = load_hgnc(
        hgnc_path=Path(args.hgnc),
        protein_coding_only=not args.all_hgnc_genes,
    )

    genes = sorted(hgnc["gene_symbol"].dropna().astype(str).unique().tolist())

    mapping_path = raw_dir / "string_mapping.tsv"
    interactions_path = raw_dir / "string_interaction_partners.tsv"

    if args.no_api:
        if not mapping_path.exists():
            raise FileNotFoundError(f"--no-api requested but missing: {mapping_path}")
        mapping = pd.read_csv(mapping_path, sep="\t", dtype=str)

        if interactions_path.exists():
            interactions = pd.read_csv(interactions_path, sep="\t", dtype=str)
        else:
            interactions = pd.DataFrame()
    else:
        mapping = map_genes_to_string(
            genes=genes,
            raw_dir=raw_dir,
            batch_size=args.batch_size,
            sleep_seconds=args.sleep,
            force=args.force_api,
        )

        if args.mapping_only:
            log("[DONE] Mapping-only mode finished.")
            return

        interactions = fetch_interaction_partners(
            mapping=mapping,
            raw_dir=raw_dir,
            batch_size=args.batch_size,
            required_score=args.required_score,
            partner_limit=args.partner_limit,
            sleep_seconds=args.sleep,
            force=args.force_api,
        )

    features = build_all_features(
        hgnc=hgnc,
        mapping=mapping,
        interactions=interactions,
        processed_dir=processed_dir,
        required_score=args.required_score,
        partner_limit=args.partner_limit,
        betweenness_k=args.betweenness_k,
        closeness_landmarks=args.closeness_landmarks,
        seed=args.seed,
        skip_expensive_graph=args.skip_expensive_graph,
    )

    merged = merge_with_hgnc(
        hgnc=hgnc,
        features=features,
        processed_dir=processed_dir,
    )

    metadata = {
        "created_at": now_iso(),
        "string_api_base": STRING_API_BASE,
        "species": STRING_SPECIES,
        "caller_identity": CALLER_IDENTITY,
        "hgnc": str(Path(args.hgnc)),
        "outdir": str(outdir.resolve()),
        "batch_size": args.batch_size,
        "required_score": args.required_score,
        "partner_limit": args.partner_limit,
        "sleep": args.sleep,
        "betweenness_k": args.betweenness_k,
        "closeness_landmarks": args.closeness_landmarks,
        "seed": args.seed,
        "skip_expensive_graph": args.skip_expensive_graph,
        "leakage_policy": {
            "included": [
                "STRING PPI topology",
                "STRING functional-association topology",
                "STRING confidence scores",
                "STRING evidence-channel summaries",
                "graph centrality",
                "local clustering",
                "component structure",
                "k-core",
            ],
            "excluded": [
                "known drug target proximity",
                "compound-target metadata",
                "drug-response data",
                "ChEMBL labels",
                "DrugBank labels",
                "approved drug annotations",
            ],
        },
        "outputs": {
            "mapping": str(mapping_path),
            "interactions": str(interactions_path),
            "clean_edges": str(processed_dir / "feature2_string_clean_edges.csv"),
            "gene_features": str(processed_dir / "feature2_string_gene_features.csv"),
            "hgnc_merged": str(processed_dir / "feature2_string_hgnc_merged.csv"),
        },
    }

    with open(outdir / "feature2_string_run_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    write_summary_report(outdir, hgnc, features, merged, args)

    log("=" * 100)
    log("[DONE]")
    log(f"[CLEAN EDGES]       {processed_dir / 'feature2_string_clean_edges.csv'}")
    log(f"[GENE FEATURES]     {processed_dir / 'feature2_string_gene_features.csv'}")
    log(f"[HGNC MERGED TABLE] {processed_dir / 'feature2_string_hgnc_merged.csv'}")
    log("=" * 100)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        log("\n[STOPPED] Interrupted by user.")
        sys.exit(130)
    except Exception as e:
        log("=" * 100)
        log("[ERROR]")
        log(str(e))
        log("=" * 100)
        raise