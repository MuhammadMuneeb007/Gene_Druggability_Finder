#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Dataset1-GenerateData.py

Build human gene-level druggability ground-truth/evidence table.

Inputs
------
databases/HGNC/hgnc_complete_set.txt
databases/ChEMBL/extracted/*/*.db
databases/OpenTargets/26.03/target/*.parquet
databases/OpenTargets/26.03/clinical_indication/*.parquet
databases/DGIdb/interactions.tsv
databases/DGIdb/genes.tsv
databases/DGIdb/categories.tsv

Optional:
Pharos GraphQL API for TDL labels:
https://pharos-api.ncats.io/graphql

Outputs
-------
Step0_Output/03_HumanGene_DruggabilityLabels.csv
Step0_Output/04_label_summary.csv
Step0_Output/raw_evidence/*.csv
Step0_Output/cache/pharos_api/*.json

Important
---------
This is a LABEL / GROUND-TRUTH / EVIDENCE table.
Do NOT use ChEMBL/OpenTargets/DGIdb/Pharos evidence columns as model input features,
otherwise the model will leak the answer.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import re
import sqlite3
import threading
import time
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import pandas as pd
import requests


PHAROS_GRAPHQL_URL = "https://pharos-api.ncats.io/graphql"
PHAROS_USER_AGENT = "DrugableGeneFinder/1.0 academic pipeline"
_PHAROS_SESSION_LOCAL = threading.local()


# =============================================================================
# Helpers
# =============================================================================

def safe_mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def normalise_symbol(x: Any) -> Optional[str]:
    if x is None or pd.isna(x):
        return None
    s = str(x).strip()
    if not s or s.lower() == "nan":
        return None
    return s.upper()


def normalise_id(x: Any) -> Optional[str]:
    if x is None or pd.isna(x):
        return None
    s = str(x).strip()
    if not s or s.lower() == "nan":
        return None
    return s


def split_multi_value(x: Any) -> List[str]:
    if x is None or pd.isna(x):
        return []
    s = str(x).strip()
    if not s or s.lower() == "nan":
        return []
    parts = re.split(r"[|,;]\s*", s)
    return [p.strip() for p in parts if p.strip()]


def clean_column_names(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    return df


def first_existing_column(df: pd.DataFrame, candidates: Sequence[str]) -> Optional[str]:
    lower_map = {str(c).lower(): c for c in df.columns}
    for c in candidates:
        if c in df.columns:
            return c
        if c.lower() in lower_map:
            return lower_map[c.lower()]
    return None


def numeric_series(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def save_json(path: Path, obj: Dict[str, Any]) -> None:
    safe_mkdir(path.parent)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, sort_keys=True)


def list_sqlite_tables(conn: sqlite3.Connection) -> List[str]:
    return pd.read_sql_query(
        "SELECT name FROM sqlite_master WHERE type='table'",
        conn,
    )["name"].tolist()


def sqlite_table_columns(conn: sqlite3.Connection, table: str) -> List[str]:
    return pd.read_sql_query(f'PRAGMA table_info("{table}")', conn)["name"].tolist()


def read_sql_table(conn: sqlite3.Connection, table: str, columns: Optional[List[str]] = None) -> pd.DataFrame:
    if columns:
        cols = ", ".join([f'"{c}"' for c in columns])
        q = f'SELECT {cols} FROM "{table}"'
    else:
        q = f'SELECT * FROM "{table}"'
    return pd.read_sql_query(q, conn)


def find_chembl_sqlite(databases_dir: Path) -> Path:
    chembl_dir = databases_dir / "ChEMBL" / "extracted"
    dbs = (
        list(chembl_dir.rglob("*.db"))
        + list(chembl_dir.rglob("*.sqlite"))
        + list(chembl_dir.rglob("*.sqlite3"))
    )

    if not dbs:
        candidates = [
            p for p in chembl_dir.rglob("*")
            if p.is_file() and "chembl" in p.name.lower() and p.stat().st_size > 1024 * 1024
        ]
        if candidates:
            return sorted(candidates, key=lambda p: p.stat().st_size, reverse=True)[0]

    if not dbs:
        raise FileNotFoundError(f"No ChEMBL SQLite DB found under {chembl_dir}")

    return sorted(dbs, key=lambda p: p.stat().st_size, reverse=True)[0]


def read_parquet_folder(folder: Path) -> pd.DataFrame:
    files = list(folder.rglob("*.parquet"))
    if not files:
        return pd.DataFrame()

    dfs = []
    for f in files:
        try:
            dfs.append(pd.read_parquet(f))
        except Exception as e:
            print(f"[WARN] Could not read parquet: {f} | {e}")

    if not dfs:
        return pd.DataFrame()

    return pd.concat(dfs, ignore_index=True)


def compact_unique_text(values: pd.Series, limit: int = 100) -> str:
    vals = []
    for v in values.dropna().astype(str):
        if not v or v.lower() in {"nan", "none"}:
            continue
        vals.extend([x for x in str(v).split("|") if x and x.lower() not in {"nan", "none"}])
    vals = sorted(set(vals))
    return "|".join(vals[:limit])


def is_scalar_missing(v: Any) -> bool:
    if v is None:
        return True
    if isinstance(v, (str, bytes, dict, list, tuple, set)):
        return False
    try:
        missing = pd.isna(v)
    except Exception:
        return False

    if hasattr(missing, "all") and not isinstance(missing, (bool, int)):
        try:
            return bool(missing.all())
        except Exception:
            return False
    return bool(missing)


def ensure_list_like(v: Any) -> List[Any]:
    if v is None:
        return []
    if isinstance(v, dict):
        return [v]
    if isinstance(v, (list, tuple, set)):
        return list(v)
    if hasattr(v, "tolist") and not isinstance(v, (str, bytes)):
        try:
            converted = v.tolist()
            if isinstance(converted, list):
                return converted
            if isinstance(converted, tuple):
                return list(converted)
            if converted is None:
                return []
            return [converted]
        except Exception:
            pass
    if is_scalar_missing(v):
        return []
    return [v]


def truthy_flag(value: Any, default_if_missing: bool = True) -> bool:
    if value is None:
        return default_if_missing
    if isinstance(value, (list, tuple, set)):
        return any(truthy_flag(v, default_if_missing=default_if_missing) for v in value)
    if hasattr(value, "tolist") and not isinstance(value, (str, bytes)):
        try:
            return truthy_flag(value.tolist(), default_if_missing=default_if_missing)
        except Exception:
            pass
    if isinstance(value, bool):
        return value
    if is_scalar_missing(value):
        return default_if_missing
    if isinstance(value, (int, float)):
        return value != 0
    text = str(value).strip().lower()
    if text in {"", "none", "nan"}:
        return default_if_missing
    if text in {"false", "f", "0", "no", "n"}:
        return False
    return True


def extract_ensembl_ids_from_value(value: Any) -> List[str]:
    """
    Recursively extract ENSG identifiers from nested Open Targets values.
    """
    out: List[str] = []

    if value is None:
        return out

    if isinstance(value, dict):
        for v in value.values():
            out.extend(extract_ensembl_ids_from_value(v))
        return out

    if isinstance(value, (list, tuple, set)):
        for item in value:
            out.extend(extract_ensembl_ids_from_value(item))
        return out

    if hasattr(value, "tolist") and not isinstance(value, (str, bytes)):
        try:
            return extract_ensembl_ids_from_value(value.tolist())
        except Exception:
            pass

    text = str(value)
    return re.findall(r"ENSG\d+", text)


def get_pharos_session() -> requests.Session:
    session = getattr(_PHAROS_SESSION_LOCAL, "session", None)
    if session is None:
        session = requests.Session()
        session.headers.update({"User-Agent": PHAROS_USER_AGENT})
        _PHAROS_SESSION_LOCAL.session = session
    return session


# =============================================================================
# HGNC
# =============================================================================

def load_hgnc(databases_dir: Path, protein_coding_only: bool = True) -> pd.DataFrame:
    path = databases_dir / "HGNC" / "hgnc_complete_set.txt"
    if not path.exists():
        raise FileNotFoundError(f"HGNC file not found: {path}")

    print(f"[HGNC] Reading: {path}")
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    df = clean_column_names(df)

    rename = {
        "symbol": "gene_symbol",
        "name": "gene_name",
        "hgnc_id": "hgnc_id",
        "ensembl_gene_id": "ensembl_gene_id",
        "entrez_id": "entrez_gene_id",
        "uniprot_ids": "uniprot_ids",
        "locus_type": "locus_type",
        "locus_group": "locus_group",
        "alias_symbol": "alias_symbol",
        "prev_symbol": "previous_symbol",
        "status": "status",
        "location": "location",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})

    if "gene_symbol" not in df.columns:
        raise ValueError("HGNC table does not contain symbol/gene_symbol.")

    if protein_coding_only and "locus_type" in df.columns:
        before = len(df)
        df = df[df["locus_type"].fillna("").str.lower().eq("gene with protein product")].copy()
        print(f"[HGNC] Protein-coding filter: {before} -> {len(df)}")

    keep = [
        "gene_symbol", "gene_name", "hgnc_id", "ensembl_gene_id", "entrez_gene_id",
        "uniprot_ids", "locus_type", "locus_group", "alias_symbol", "previous_symbol",
        "status", "location",
    ]
    keep = [c for c in keep if c in df.columns]
    df = df[keep].copy()

    for c in ["ensembl_gene_id", "entrez_gene_id", "uniprot_ids"]:
        if c not in df.columns:
            df[c] = None

    df["gene_symbol"] = df["gene_symbol"].astype(str).str.strip()
    df["gene_symbol_norm"] = df["gene_symbol"].apply(normalise_symbol)
    df["ensembl_gene_id"] = df["ensembl_gene_id"].apply(normalise_id)
    df["entrez_gene_id"] = df["entrez_gene_id"].apply(normalise_id)
    df["uniprot_id_primary"] = df["uniprot_ids"].apply(
        lambda x: split_multi_value(x)[0] if split_multi_value(x) else None
    )

    df = df[df["gene_symbol_norm"].notna()].copy()
    df = df.drop_duplicates(subset=["gene_symbol_norm"]).reset_index(drop=True)

    print(f"[HGNC] Gene universe rows: {len(df)}")
    return df


def build_hgnc_uniprot_bridge(hgnc: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in hgnc.iterrows():
        ids = split_multi_value(r.get("uniprot_ids"))
        for uid in ids:
            rows.append({
                "gene_symbol_norm": r["gene_symbol_norm"],
                "uniprot_id": uid,
            })
    bridge = pd.DataFrame(rows)
    if bridge.empty:
        return pd.DataFrame(columns=["gene_symbol_norm", "uniprot_id"])
    bridge = bridge.drop_duplicates()
    return bridge


# =============================================================================
# ChEMBL
# =============================================================================

def build_chembl_gene_evidence(databases_dir: Path, hgnc: pd.DataFrame) -> pd.DataFrame:
    db_path = find_chembl_sqlite(databases_dir)
    print(f"[ChEMBL] SQLite DB: {db_path}")

    conn = sqlite3.connect(str(db_path))
    tables = set(list_sqlite_tables(conn))
    print(f"[ChEMBL] Tables detected: {len(tables)}")

    needed = {"target_dictionary", "target_components", "component_sequences"}
    missing = [t for t in needed if t not in tables]
    if missing:
        print(f"[ChEMBL] Missing tables: {missing}")
        conn.close()
        return pd.DataFrame()

    target_map_query = """
    SELECT
        td.tid,
        td.chembl_id AS target_chembl_id,
        td.pref_name AS chembl_target_name,
        td.target_type AS chembl_target_type,
        td.organism AS chembl_organism,
        cs.accession AS uniprot_id,
        cs.component_id,
        cs.description AS component_description,
        cs.component_type
    FROM target_dictionary td
    JOIN target_components tc ON td.tid = tc.tid
    JOIN component_sequences cs ON tc.component_id = cs.component_id
    WHERE LOWER(COALESCE(td.organism, '')) LIKE '%homo sapiens%'
       OR LOWER(COALESCE(td.organism, '')) LIKE '%human%'
    """
    target_map = pd.read_sql_query(target_map_query, conn)
    target_map = clean_column_names(target_map)
    target_map["uniprot_id"] = target_map["uniprot_id"].apply(normalise_id)
    print(f"[ChEMBL] Human target-component rows: {len(target_map)}")

    mechanism_table = None
    for candidate in ["drug_mechanism", "mechanism"]:
        if candidate in tables:
            mechanism_table = candidate
            break
    print(f"[ChEMBL] Mechanism table: {mechanism_table or 'not found'}")

    # Mechanism evidence by tid.
    mech_by_target = pd.DataFrame()
    if mechanism_table is not None:
        mech_cols = sqlite_table_columns(conn, mechanism_table)
        cols = [c for c in ["tid", "molregno", "action_type", "mechanism_of_action", "max_phase"] if c in mech_cols]
        mechanism = read_sql_table(conn, mechanism_table, cols if cols else None)
        mechanism = clean_column_names(mechanism)

        if "tid" in mechanism.columns:
            if "max_phase" in mechanism.columns:
                mechanism["max_phase_num"] = numeric_series(mechanism["max_phase"])
            else:
                mechanism["max_phase_num"] = float("nan")

            mech_by_target = mechanism.groupby("tid").agg(
                chembl_num_mechanisms=("tid", "size"),
                chembl_max_mechanism_phase=("max_phase_num", "max"),
                chembl_num_phase4_mechanisms=("max_phase_num", lambda x: int((pd.to_numeric(x, errors="coerce") >= 4).sum())),
            ).reset_index()

    # Drug evidence by tid.
    drug_by_target = pd.DataFrame()
    if mechanism_table is not None and "molecule_dictionary" in tables:
        mech_cols = sqlite_table_columns(conn, mechanism_table)
        mol_cols = sqlite_table_columns(conn, "molecule_dictionary")

        if "tid" in mech_cols and "molregno" in mech_cols and "molregno" in mol_cols:
            select_cols = ["m.tid", "md.molregno"]
            if "chembl_id" in mol_cols:
                select_cols.append("md.chembl_id AS molecule_chembl_id")
            if "pref_name" in mol_cols:
                select_cols.append("md.pref_name AS molecule_name")
            if "max_phase" in mol_cols:
                select_cols.append("md.max_phase AS molecule_max_phase")
            if "first_approval" in mol_cols:
                select_cols.append("md.first_approval AS first_approval")

            q = f"""
            SELECT {", ".join(select_cols)}
            FROM "{mechanism_table}" m
            JOIN molecule_dictionary md ON m.molregno = md.molregno
            """
            drug_mech = pd.read_sql_query(q, conn)
            drug_mech = clean_column_names(drug_mech)

            if "molecule_max_phase" in drug_mech.columns:
                drug_mech["molecule_max_phase_num"] = numeric_series(drug_mech["molecule_max_phase"])
            else:
                drug_mech["molecule_max_phase_num"] = float("nan")

            if "first_approval" in drug_mech.columns:
                drug_mech["first_approval_num"] = numeric_series(drug_mech["first_approval"])
            else:
                drug_mech["first_approval_num"] = float("nan")

            drug_by_target = drug_mech.groupby("tid").agg(
                chembl_num_drugs=("molregno", pd.Series.nunique),
                chembl_max_drug_phase=("molecule_max_phase_num", "max"),
                chembl_num_phase4_drugs=("molecule_max_phase_num", lambda x: int((pd.to_numeric(x, errors="coerce") >= 4).sum())),
                chembl_num_first_approved_drugs=("first_approval_num", lambda x: int(pd.to_numeric(x, errors="coerce").notna().sum())),
            ).reset_index()

            if "molecule_name" in drug_mech.columns:
                names = drug_mech.groupby("tid")["molecule_name"].apply(
                    lambda x: compact_unique_text(x, limit=100)
                ).reset_index(name="chembl_drug_names")
                drug_by_target = drug_by_target.merge(names, on="tid", how="left")

    # Activity evidence by tid using SQL aggregation to avoid huge memory.
    activity_by_target = pd.DataFrame()
    if "activities" in tables and "assays" in tables:
        act_cols = sqlite_table_columns(conn, "activities")
        assay_cols = sqlite_table_columns(conn, "assays")
        if "assay_id" in act_cols and "pchembl_value" in act_cols and "assay_id" in assay_cols and "tid" in assay_cols:
            print("[ChEMBL] Aggregating activity evidence from SQLite.")
            standard_filter = ""
            if "standard_type" in act_cols:
                standard_filter = """
                AND UPPER(COALESCE(a.standard_type, '')) IN ('IC50', 'KI', 'KD', 'EC50', 'AC50')
                """

            q = f"""
            SELECT
                ass.tid AS tid,
                COUNT(*) AS chembl_num_activity_records,
                SUM(
                    CASE
                        WHEN CAST(a.pchembl_value AS REAL) >= 6.0
                        {standard_filter}
                        THEN 1 ELSE 0
                    END
                ) AS chembl_num_potent_activity_records,
                MAX(CAST(a.pchembl_value AS REAL)) AS chembl_max_pchembl
            FROM activities a
            JOIN assays ass ON a.assay_id = ass.assay_id
            WHERE ass.tid IS NOT NULL
            GROUP BY ass.tid
            """
            activity_by_target = pd.read_sql_query(q, conn)
            activity_by_target = clean_column_names(activity_by_target)

    conn.close()

    target_level = target_map.drop_duplicates(
        subset=["tid", "target_chembl_id", "uniprot_id"]
    ).copy()

    for extra in [mech_by_target, drug_by_target, activity_by_target]:
        if extra is not None and not extra.empty and "tid" in extra.columns:
            target_level = target_level.merge(extra, on="tid", how="left")

    zero_cols = [
        "chembl_num_mechanisms",
        "chembl_num_phase4_mechanisms",
        "chembl_num_drugs",
        "chembl_num_phase4_drugs",
        "chembl_num_first_approved_drugs",
        "chembl_num_activity_records",
        "chembl_num_potent_activity_records",
    ]
    for c in zero_cols:
        if c not in target_level.columns:
            target_level[c] = 0
        target_level[c] = pd.to_numeric(target_level[c], errors="coerce").fillna(0).astype(int)

    for c in ["chembl_max_mechanism_phase", "chembl_max_drug_phase", "chembl_max_pchembl"]:
        if c not in target_level.columns:
            target_level[c] = float("nan")
        target_level[c] = pd.to_numeric(target_level[c], errors="coerce")

    if "chembl_drug_names" not in target_level.columns:
        target_level["chembl_drug_names"] = None

    target_level["chembl_max_phase"] = target_level[
        ["chembl_max_mechanism_phase", "chembl_max_drug_phase"]
    ].max(axis=1, skipna=True)

    target_level["chembl_has_target"] = 1
    target_level["chembl_has_mechanism"] = (target_level["chembl_num_mechanisms"] > 0).astype(int)
    target_level["chembl_has_phase4_or_approved"] = (
        (target_level["chembl_max_phase"] >= 4)
        | (target_level["chembl_num_phase4_mechanisms"] > 0)
        | (target_level["chembl_num_phase4_drugs"] > 0)
        | (target_level["chembl_num_first_approved_drugs"] > 0)
    ).astype(int)
    target_level["chembl_has_potent_compound"] = (
        target_level["chembl_num_potent_activity_records"] > 0
    ).astype(int)

    # Aggregate by UniProt.
    by_uniprot = target_level[target_level["uniprot_id"].notna()].copy()
    if by_uniprot.empty:
        print("[ChEMBL] No UniProt evidence built.")
        return pd.DataFrame()

    chembl_uniprot = by_uniprot.groupby("uniprot_id").agg(
        chembl_has_target=("chembl_has_target", "max"),
        chembl_target_chembl_ids=("target_chembl_id", lambda x: "|".join(sorted(set(map(str, x.dropna())))[:100])),
        chembl_target_names=("chembl_target_name", lambda x: "|".join(sorted(set(map(str, x.dropna())))[:100])),
        chembl_num_targets=("target_chembl_id", pd.Series.nunique),
        chembl_num_mechanisms=("chembl_num_mechanisms", "sum"),
        chembl_num_drugs=("chembl_num_drugs", "sum"),
        chembl_max_phase=("chembl_max_phase", "max"),
        chembl_has_mechanism=("chembl_has_mechanism", "max"),
        chembl_has_phase4_or_approved=("chembl_has_phase4_or_approved", "max"),
        chembl_num_activity_records=("chembl_num_activity_records", "sum"),
        chembl_num_potent_activity_records=("chembl_num_potent_activity_records", "sum"),
        chembl_max_pchembl=("chembl_max_pchembl", "max"),
        chembl_has_potent_compound=("chembl_has_potent_compound", "max"),
        chembl_drug_names=("chembl_drug_names", lambda x: compact_unique_text(x, limit=100)),
    ).reset_index()

    # Map ChEMBL UniProt evidence to HGNC genes.
    bridge = build_hgnc_uniprot_bridge(hgnc)
    mapped = bridge.merge(chembl_uniprot, on="uniprot_id", how="inner")

    if mapped.empty:
        print("[ChEMBL] No HGNC genes matched through UniProt.")
        return pd.DataFrame()

    chembl_gene = mapped.groupby("gene_symbol_norm").agg(
        chembl_has_target=("chembl_has_target", "max"),
        chembl_target_chembl_ids=("chembl_target_chembl_ids", lambda x: compact_unique_text(x, limit=100)),
        chembl_target_names=("chembl_target_names", lambda x: compact_unique_text(x, limit=100)),
        chembl_num_targets=("chembl_num_targets", "sum"),
        chembl_num_mechanisms=("chembl_num_mechanisms", "sum"),
        chembl_num_drugs=("chembl_num_drugs", "sum"),
        chembl_max_phase=("chembl_max_phase", "max"),
        chembl_has_mechanism=("chembl_has_mechanism", "max"),
        chembl_has_phase4_or_approved=("chembl_has_phase4_or_approved", "max"),
        chembl_num_activity_records=("chembl_num_activity_records", "sum"),
        chembl_num_potent_activity_records=("chembl_num_potent_activity_records", "sum"),
        chembl_max_pchembl=("chembl_max_pchembl", "max"),
        chembl_has_potent_compound=("chembl_has_potent_compound", "max"),
        chembl_drug_names=("chembl_drug_names", lambda x: compact_unique_text(x, limit=100)),
    ).reset_index()

    print(f"[ChEMBL] Gene-level evidence rows: {len(chembl_gene)}")
    print(f"[ChEMBL] Genes with approved/phase4 evidence: {int(chembl_gene['chembl_has_phase4_or_approved'].sum())}")
    print(f"[ChEMBL] Genes with potent compound evidence: {int(chembl_gene['chembl_has_potent_compound'].sum())}")

    return chembl_gene


# =============================================================================
# Open Targets
# =============================================================================

def parse_ot_tractability_cell(v: Any) -> Dict[str, int]:
    """
    Conservative Open Targets tractability parser.

    Open Targets modality codes:
      SM = small molecule
      AB = antibody
      PR = targeted protein degradation / PROTAC-related evidence
      OC = other clinical modality

    We only count entries where value is truthy. PR and OC are handled
    conservatively to avoid inflating druggability labels from generic
    degradation or broad modality evidence.
    """
    out = {
        "ot_has_small_molecule_tractability": 0,
        "ot_has_antibody_tractability": 0,
        "ot_has_protac_tractability": 0,
        "ot_has_other_modality_tractability": 0,
    }

    strong_pr_ids = {
        "Approved Drug",
        "Advanced Clinical",
        "Phase 1 Clinical",
        "Literature",
        "Small Molecule Binder",
    }

    strong_oc_ids = {
        "Approved Drug",
        "Advanced Clinical",
        "Phase 1 Clinical",
    }

    items = ensure_list_like(v)
    if not items:
        return out

    for item in items:
        if not isinstance(item, dict):
            continue

        modality = str(item.get("modality", "")).strip().upper()
        criterion = str(item.get("id", "")).strip()
        value = item.get("value", False)

        if not truthy_flag(value, default_if_missing=False):
            continue

        if modality == "SM":
            out["ot_has_small_molecule_tractability"] = 1
        elif modality == "AB":
            out["ot_has_antibody_tractability"] = 1
        elif modality == "PR":
            if criterion in strong_pr_ids:
                out["ot_has_protac_tractability"] = 1
        elif modality == "OC":
            if criterion in strong_oc_ids:
                out["ot_has_other_modality_tractability"] = 1

    return out


def build_opentargets_evidence(databases_dir: Path, release: str) -> pd.DataFrame:
    ot_dir = databases_dir / "OpenTargets" / release
    target_dir = ot_dir / "target"
    clinical_indication_dir = ot_dir / "clinical_indication"
    clinical_target_dir = ot_dir / "clinical_target"
    drug_moa_dir = ot_dir / "drug_mechanism_of_action"

    print(f"[OpenTargets] Reading target dataset: {target_dir}")
    target_df = read_parquet_folder(target_dir)

    if target_df.empty:
        print("[OpenTargets] target dataset empty/missing.")
        target_evidence = pd.DataFrame()
    else:
        target_df = clean_column_names(target_df)
        print(f"[OpenTargets] target rows: {len(target_df)}, columns: {list(target_df.columns)[:35]}")

        ensembl_col = first_existing_column(target_df, ["id", "ensemblId", "ensembl_gene_id", "targetId"])
        symbol_col = first_existing_column(target_df, ["approvedSymbol", "approved_symbol", "symbol"])
        name_col = first_existing_column(target_df, ["approvedName", "approved_name", "name"])
        tract_col = first_existing_column(target_df, ["tractability"])

        if ensembl_col is None:
            target_evidence = pd.DataFrame()
        else:
            target_evidence = pd.DataFrame()
            target_evidence["ensembl_gene_id"] = target_df[ensembl_col].astype(str)
            target_evidence["opentargets_has_target"] = 1

            if symbol_col:
                target_evidence["opentargets_approved_symbol"] = target_df[symbol_col].astype(str)
                target_evidence["gene_symbol_norm"] = target_evidence["opentargets_approved_symbol"].apply(normalise_symbol)
            else:
                target_evidence["opentargets_approved_symbol"] = None
                target_evidence["gene_symbol_norm"] = None

            if name_col:
                target_evidence["opentargets_approved_name"] = target_df[name_col].astype(str)
            else:
                target_evidence["opentargets_approved_name"] = None

            if tract_col:
                tract_rows = target_df[tract_col].apply(parse_ot_tractability_cell).apply(pd.Series)
                target_evidence = pd.concat([target_evidence.reset_index(drop=True), tract_rows.reset_index(drop=True)], axis=1)
            else:
                target_evidence["ot_has_small_molecule_tractability"] = 0
                target_evidence["ot_has_antibody_tractability"] = 0
                target_evidence["ot_has_protac_tractability"] = 0
                target_evidence["ot_has_other_modality_tractability"] = 0

            target_evidence = target_evidence.drop_duplicates(subset=["ensembl_gene_id"])

    print(f"[OpenTargets] Reading clinical indication dataset: {clinical_indication_dir}")
    clinical_df = read_parquet_folder(clinical_indication_dir)
    clinical_by_target = pd.DataFrame()
    drug_target_map = pd.DataFrame()

    print(f"[OpenTargets] Reading drug mechanism dataset: {drug_moa_dir}")
    drug_moa_df = read_parquet_folder(drug_moa_dir)
    if not drug_moa_df.empty:
        drug_moa_df = clean_column_names(drug_moa_df)
        print(f"[OpenTargets] drug mechanism rows: {len(drug_moa_df)}, columns: {list(drug_moa_df.columns)[:30]}")

        drug_col = first_existing_column(drug_moa_df, ["drugId", "drug_id", "drug", "chemblId"])
        target_candidate_cols = [
            c for c in [
                first_existing_column(drug_moa_df, ["targetId", "target_id"]),
                first_existing_column(drug_moa_df, ["targets", "linkedTargets", "targetsInfo"]),
                first_existing_column(drug_moa_df, ["approvedTargets", "mechanismsOfAction"]),
            ] if c is not None
        ]

        if drug_col and target_candidate_cols:
            map_rows = []
            for _, row in drug_moa_df.iterrows():
                drug_id = normalise_id(row.get(drug_col))
                if not drug_id:
                    continue

                ensembl_ids: List[str] = []
                for col in target_candidate_cols:
                    ensembl_ids.extend(extract_ensembl_ids_from_value(row.get(col)))

                for ensembl_id in sorted(set(ensembl_ids)):
                    map_rows.append({
                        "drugId": drug_id,
                        "ensembl_gene_id": ensembl_id,
                    })

            if map_rows:
                drug_target_map = pd.DataFrame(map_rows).drop_duplicates()
                print(f"[OpenTargets] drug->target mappings: {len(drug_target_map)}")

    if clinical_df.empty and clinical_target_dir.exists():
        print(f"[OpenTargets] Reading fallback clinical target dataset: {clinical_target_dir}")
        clinical_df = read_parquet_folder(clinical_target_dir)

    if not clinical_df.empty:
        clinical_df = clean_column_names(clinical_df)
        print(f"[OpenTargets] clinical rows: {len(clinical_df)}, columns: {list(clinical_df.columns)[:30]}")

        target_col = first_existing_column(clinical_df, ["targetId", "target_id", "target", "targets"])
        phase_col = first_existing_column(clinical_df, ["maxClinicalStage", "clinicalPhase", "clinical_phase", "maxPhase", "max_phase", "phase"])
        drug_col = first_existing_column(clinical_df, ["drugId", "drug_id", "chemblId", "drug"])
        disease_col = first_existing_column(clinical_df, ["diseaseId", "disease_id", "disease", "efoId"])

        if target_col:
            tmp = clinical_df.copy()
            tmp["ensembl_gene_id"] = tmp[target_col].apply(
                lambda x: (lambda ids: ids[0] if ids else None)(extract_ensembl_ids_from_value(x))
            )
            tmp = tmp[tmp["ensembl_gene_id"].notna()].copy()
        elif drug_col and not drug_target_map.empty:
            tmp = clinical_df.copy()
            tmp["drugId"] = tmp[drug_col].apply(normalise_id)
            tmp = tmp.merge(drug_target_map, on="drugId", how="inner")
        else:
            tmp = pd.DataFrame()

        if not tmp.empty:
            tmp["_phase_num"] = numeric_series(tmp[phase_col]) if phase_col else float("nan")
            tmp["_drug"] = tmp[drug_col].astype(str) if drug_col else None
            tmp["_disease"] = tmp[disease_col].astype(str) if disease_col else None

            clinical_by_target = tmp.groupby("ensembl_gene_id").agg(
                ot_has_clinical_indication=("ensembl_gene_id", lambda x: 1),
                ot_num_clinical_indications=("ensembl_gene_id", "size"),
                ot_num_clinical_drugs=("_drug", lambda x: len(set([v for v in x.dropna().astype(str) if v and v != "None"]))),
                ot_num_clinical_diseases=("_disease", lambda x: len(set([v for v in x.dropna().astype(str) if v and v != "None"]))),
                ot_max_clinical_phase=("_phase_num", "max"),
            ).reset_index()

    if target_evidence.empty and clinical_by_target.empty:
        return pd.DataFrame()

    if target_evidence.empty:
        out = clinical_by_target.copy()
    elif clinical_by_target.empty:
        out = target_evidence.copy()
    else:
        out = target_evidence.merge(clinical_by_target, on="ensembl_gene_id", how="outer")

    zero_cols = [
        "opentargets_has_target",
        "ot_has_clinical_indication",
        "ot_num_clinical_indications",
        "ot_num_clinical_drugs",
        "ot_num_clinical_diseases",
        "ot_has_small_molecule_tractability",
        "ot_has_antibody_tractability",
        "ot_has_protac_tractability",
        "ot_has_other_modality_tractability",
    ]
    for c in zero_cols:
        if c not in out.columns:
            out[c] = 0
        out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0).astype(int)

    if "ot_max_clinical_phase" not in out.columns:
        out["ot_max_clinical_phase"] = float("nan")
    out["ot_max_clinical_phase"] = pd.to_numeric(out["ot_max_clinical_phase"], errors="coerce")

    print(f"[OpenTargets] Evidence rows: {len(out)}")
    print(f"[OpenTargets] Small molecule tractability genes: {int(out['ot_has_small_molecule_tractability'].sum())}")
    print(f"[OpenTargets] Antibody tractability genes: {int(out['ot_has_antibody_tractability'].sum())}")
    print(f"[OpenTargets] PROTAC tractability genes: {int(out['ot_has_protac_tractability'].sum())}")

    return out


# =============================================================================
# DGIdb
# =============================================================================

def load_dgidb_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    print(f"[DGIdb] Reading: {path}")
    return clean_column_names(pd.read_csv(path, sep="\t", dtype=str, low_memory=False))


def build_dgidb_evidence(databases_dir: Path) -> pd.DataFrame:
    dgidb_dir = databases_dir / "DGIdb"

    interactions = load_dgidb_tsv(dgidb_dir / "interactions.tsv")
    genes = load_dgidb_tsv(dgidb_dir / "genes.tsv")
    categories = load_dgidb_tsv(dgidb_dir / "categories.tsv")

    parts = []

    if not interactions.empty:
        print(f"[DGIdb] interactions rows: {len(interactions)}, columns: {list(interactions.columns)}")

        gene_col = first_existing_column(interactions, ["gene_name", "gene_claim_name", "gene", "symbol"])
        drug_col = first_existing_column(interactions, ["drug_name", "drug_claim_name", "drug"])
        interaction_col = first_existing_column(interactions, ["interaction_type", "interaction_types"])
        source_col = first_existing_column(interactions, ["interaction_source_db_name", "source_db_name", "source"])

        if gene_col:
            tmp = interactions.copy()
            tmp["gene_symbol_norm"] = tmp[gene_col].apply(normalise_symbol)
            tmp["_drug"] = tmp[drug_col].astype(str) if drug_col else None
            tmp["_interaction_type"] = tmp[interaction_col].astype(str) if interaction_col else None
            tmp["_source"] = tmp[source_col].astype(str) if source_col else None
            tmp = tmp[tmp["gene_symbol_norm"].notna()].copy()

            inter = tmp.groupby("gene_symbol_norm").agg(
                dgidb_has_interaction=("gene_symbol_norm", lambda x: 1),
                dgidb_num_interactions=("gene_symbol_norm", "size"),
                dgidb_num_drugs=("_drug", lambda x: len(set([v for v in x.dropna().astype(str) if v and v != "None"]))),
                dgidb_drug_names=("_drug", lambda x: compact_unique_text(x, limit=100)),
                dgidb_interaction_types=("_interaction_type", lambda x: compact_unique_text(x, limit=100)),
                dgidb_sources=("_source", lambda x: compact_unique_text(x, limit=100)),
            ).reset_index()

            parts.append(inter)

    if not genes.empty:
        print(f"[DGIdb] genes rows: {len(genes)}, columns: {list(genes.columns)}")

        gene_col = first_existing_column(genes, ["gene_name", "gene_claim_name", "gene", "symbol", "name"])
        source_col = first_existing_column(genes, ["source_db_name", "source"])

        if gene_col:
            tmp = genes.copy()
            tmp["gene_symbol_norm"] = tmp[gene_col].apply(normalise_symbol)
            tmp["_source"] = tmp[source_col].astype(str) if source_col else None
            tmp = tmp[tmp["gene_symbol_norm"].notna()].copy()

            gene_rec = tmp.groupby("gene_symbol_norm").agg(
                dgidb_has_gene_record=("gene_symbol_norm", lambda x: 1),
                dgidb_gene_sources=("_source", lambda x: compact_unique_text(x, limit=100)),
            ).reset_index()
            parts.append(gene_rec)

    if not categories.empty:
        print(f"[DGIdb] categories rows: {len(categories)}, columns: {list(categories.columns)}")

        gene_col = first_existing_column(categories, ["name", "gene_name", "gene", "symbol"])
        cat_col = first_existing_column(categories, ["name-2", "category", "categories", "source_db_name"])

        if gene_col:
            tmp = categories.copy()
            tmp["gene_symbol_norm"] = tmp[gene_col].apply(normalise_symbol)
            tmp["_category"] = tmp[cat_col].astype(str) if cat_col else None
            tmp = tmp[tmp["gene_symbol_norm"].notna()].copy()

            cat = tmp.groupby("gene_symbol_norm").agg(
                dgidb_categories_from_categories_file=("_category", lambda x: compact_unique_text(x, limit=100)),
            ).reset_index()
            cat["dgidb_has_druggable_category_from_categories_file"] = (
                cat["dgidb_categories_from_categories_file"].fillna("").ne("").astype(int)
            )
            parts.append(cat)

    if not parts:
        return pd.DataFrame()

    out = parts[0]
    for p in parts[1:]:
        out = out.merge(p, on="gene_symbol_norm", how="outer")

    for c in [
        "dgidb_has_interaction",
        "dgidb_num_interactions",
        "dgidb_num_drugs",
        "dgidb_has_gene_record",
        "dgidb_has_druggable_category_from_categories_file",
    ]:
        if c not in out.columns:
            out[c] = 0
        out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0).astype(int)

    out["dgidb_has_any_category"] = out["dgidb_has_druggable_category_from_categories_file"].astype(int)

    print(f"[DGIdb] Evidence rows: {len(out)}")
    return out


# =============================================================================
# Pharos API
# =============================================================================

def pharos_cache_key(symbol: str) -> str:
    symbol = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(symbol).strip())
    return f"{symbol}.json"


def query_pharos_api_for_symbol(
    gene_symbol: str,
    cache_dir: Path,
    sleep_seconds: float = 0.2,
    timeout: int = 60,
    force: bool = False,
) -> Dict[str, Any]:
    safe_mkdir(cache_dir)
    gene_symbol = str(gene_symbol).strip()
    cache_file = cache_dir / pharos_cache_key(gene_symbol)

    if cache_file.exists() and not force:
        try:
            with open(cache_file, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            pass

    query = """
    query targetDetails($sym: String!) {
      target(q: {sym: $sym}) {
        name
        tdl
        fam
        sym
        description
        novelty
      }
    }
    """

    payload = {
        "query": query,
        "variables": {"sym": gene_symbol},
    }

    time.sleep(sleep_seconds)

    try:
        r = get_pharos_session().post(
            PHAROS_GRAPHQL_URL,
            json=payload,
            timeout=timeout,
        )
        r.raise_for_status()
        data = r.json()

        out = {
            "gene_symbol": gene_symbol,
            "status": "ok",
            "response": data,
        }

    except Exception as e:
        out = {
            "gene_symbol": gene_symbol,
            "status": f"error:{type(e).__name__}",
            "error": str(e),
            "response": None,
        }

    with open(cache_file, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, sort_keys=True)

    return out


def parse_pharos_response(gene_symbol: str, result: Dict[str, Any]) -> Dict[str, Any]:
    out = {
        "gene_symbol_norm": normalise_symbol(gene_symbol),
        "pharos_api_available": 1,
        "pharos_api_status": result.get("status"),
        "pharos_name": None,
        "pharos_sym": None,
        "pharos_tdl": None,
        "pharos_fam": None,
        "pharos_novelty": None,
        "pharos_description": None,
        "pharos_is_tclin": 0,
        "pharos_is_tchem": 0,
        "pharos_is_tbio": 0,
        "pharos_is_tdark": 0,
    }

    if result.get("status") != "ok":
        out["pharos_api_available"] = 0
        return out

    response = result.get("response") or {}

    if response.get("errors"):
        out["pharos_api_status"] = "graphql_error"
        out["pharos_api_available"] = 0
        return out

    target = (response.get("data") or {}).get("target")
    if not target:
        out["pharos_api_status"] = "not_found"
        return out

    tdl = str(target.get("tdl") or "").upper().strip()

    out.update({
        "pharos_api_status": "ok",
        "pharos_name": target.get("name"),
        "pharos_sym": target.get("sym"),
        "pharos_tdl": tdl if tdl else None,
        "pharos_fam": target.get("fam"),
        "pharos_novelty": target.get("novelty"),
        "pharos_description": target.get("description"),
        "pharos_is_tclin": 1 if tdl == "TCLIN" else 0,
        "pharos_is_tchem": 1 if tdl == "TCHEM" else 0,
        "pharos_is_tbio": 1 if tdl == "TBIO" else 0,
        "pharos_is_tdark": 1 if tdl == "TDARK" else 0,
    })

    return out


def build_pharos_api_evidence(
    hgnc: pd.DataFrame,
    output_dir: Path,
    sleep_seconds: float = 0.0,
    timeout: int = 60,
    force: bool = False,
    max_genes: Optional[int] = None,
    workers: int = 16,
    save_every: int = 500,
) -> pd.DataFrame:
    cache_dir = output_dir / "cache" / "pharos_api"
    raw_dir = output_dir / "raw_evidence"
    safe_mkdir(cache_dir)
    safe_mkdir(raw_dir)

    genes = hgnc[["gene_symbol", "gene_symbol_norm"]].drop_duplicates().copy()

    if max_genes is not None and max_genes > 0:
        genes = genes.head(max_genes).copy()

    total = len(genes)
    workers = max(1, int(workers))
    save_every = max(1, int(save_every))
    progress_every = max(50, min(save_every, 200))
    rows_by_symbol: Dict[str, Dict[str, Any]] = {}

    print(f"[Pharos API] Querying {total} genes")
    print(f"[Pharos API] Cache: {cache_dir}")
    print(f"[Pharos API] Workers: {workers}")

    gene_symbols = [str(x).strip() for x in genes["gene_symbol"].tolist() if str(x).strip()]
    if workers == 1:
        completed = 0
        for gene_symbol in gene_symbols:
            result = query_pharos_api_for_symbol(
                gene_symbol=gene_symbol,
                cache_dir=cache_dir,
                sleep_seconds=sleep_seconds,
                timeout=timeout,
                force=force,
            )
            rows_by_symbol[gene_symbol] = parse_pharos_response(gene_symbol, result)
            completed += 1
            if completed % progress_every == 0 or completed == total:
                print(f"[Pharos API] Progress: {completed}/{total}")
            if completed % save_every == 0:
                partial = [rows_by_symbol[s] for s in gene_symbols if s in rows_by_symbol]
                pd.DataFrame(partial).to_csv(raw_dir / "pharos_api_tdl_evidence.partial.csv", index=False)
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            future_to_symbol = {
                executor.submit(
                    query_pharos_api_for_symbol,
                    gene_symbol,
                    cache_dir,
                    sleep_seconds,
                    timeout,
                    force,
                ): gene_symbol
                for gene_symbol in gene_symbols
            }

            completed = 0
            for future in as_completed(future_to_symbol):
                gene_symbol = future_to_symbol[future]
                try:
                    result = future.result()
                except Exception as e:
                    result = {
                        "gene_symbol": gene_symbol,
                        "status": f"error:{type(e).__name__}",
                        "error": str(e),
                        "response": None,
                    }
                rows_by_symbol[gene_symbol] = parse_pharos_response(gene_symbol, result)
                completed += 1

                if completed % progress_every == 0 or completed == total:
                    print(f"[Pharos API] Progress: {completed}/{total}")

                if completed % save_every == 0:
                    partial = [rows_by_symbol[s] for s in gene_symbols if s in rows_by_symbol]
                    pd.DataFrame(partial).to_csv(raw_dir / "pharos_api_tdl_evidence.partial.csv", index=False)

    ordered_rows = [rows_by_symbol[s] for s in gene_symbols if s in rows_by_symbol]
    out = pd.DataFrame(ordered_rows)
    out_path = raw_dir / "pharos_api_tdl_evidence.csv"
    out.to_csv(out_path, index=False)

    print(f"[Pharos API] Saved: {out_path}")
    print(f"[Pharos API] Rows: {len(out)}")
    if "pharos_tdl" in out.columns:
        print(out["pharos_tdl"].value_counts(dropna=False).to_string())

    return out


# =============================================================================
# Merge and labels
# =============================================================================

def merge_evidence(
    hgnc: pd.DataFrame,
    chembl: pd.DataFrame,
    opentargets: pd.DataFrame,
    dgidb: pd.DataFrame,
    pharos: pd.DataFrame,
) -> pd.DataFrame:
    df = hgnc.copy()

    if chembl is not None and not chembl.empty:
        df = df.merge(chembl, on="gene_symbol_norm", how="left")
    else:
        print("[MERGE] No ChEMBL evidence.")

    if opentargets is not None and not opentargets.empty:
        df = df.merge(opentargets, on="ensembl_gene_id", how="left", suffixes=("", "_ot"))
    else:
        print("[MERGE] No Open Targets evidence.")

    if dgidb is not None and not dgidb.empty:
        df = df.merge(dgidb, on="gene_symbol_norm", how="left")
    else:
        print("[MERGE] No DGIdb evidence.")

    if pharos is not None and not pharos.empty:
        df = df.merge(pharos, on="gene_symbol_norm", how="left")
    else:
        print("[MERGE] No Pharos API evidence.")

    return add_final_labels(df)


def add_final_labels(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    defaults = {
        "chembl_has_target": 0,
        "chembl_has_mechanism": 0,
        "chembl_has_phase4_or_approved": 0,
        "chembl_has_potent_compound": 0,
        "chembl_num_mechanisms": 0,
        "chembl_num_drugs": 0,
        "chembl_num_potent_activity_records": 0,
        "chembl_max_phase": float("nan"),

        "opentargets_has_target": 0,
        "ot_has_clinical_indication": 0,
        "ot_num_clinical_indications": 0,
        "ot_num_clinical_drugs": 0,
        "ot_max_clinical_phase": float("nan"),
        "ot_has_small_molecule_tractability": 0,
        "ot_has_antibody_tractability": 0,
        "ot_has_protac_tractability": 0,
        "ot_has_other_modality_tractability": 0,

        "dgidb_has_interaction": 0,
        "dgidb_num_interactions": 0,
        "dgidb_num_drugs": 0,
        "dgidb_has_any_category": 0,

        "pharos_api_available": 0,
        "pharos_is_tclin": 0,
        "pharos_is_tchem": 0,
        "pharos_is_tbio": 0,
        "pharos_is_tdark": 0,
    }

    for c, default in defaults.items():
        if c not in df.columns:
            df[c] = default
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(default)

    df["clinical_target_label"] = (
        (df["chembl_has_phase4_or_approved"] == 1)
        | (df["pharos_is_tclin"] == 1)
    ).astype(int)

    df["clinical_investigation_label"] = (
        (df["clinical_target_label"] == 0)
        & (
            (df["ot_has_clinical_indication"] == 1)
            | (df["ot_max_clinical_phase"] >= 1)
            | ((df["chembl_max_phase"] >= 1) & (df["chembl_max_phase"] < 4))
        )
    ).astype(int)

    df["chemical_tractable_label"] = (
        (df["chembl_has_potent_compound"] == 1)
        | (df["pharos_is_tchem"] == 1)
        | (df["chembl_num_potent_activity_records"] > 0)
        | (df["ot_has_small_molecule_tractability"] == 1)
        | (df["ot_has_protac_tractability"] == 1)
    ).astype(int)

    df["small_molecule_druggable_label"] = (
        (df["clinical_target_label"] == 1)
        | (df["chemical_tractable_label"] == 1)
        | (df["chembl_has_mechanism"] == 1)
        | (df["ot_has_small_molecule_tractability"] == 1)
        | (df["ot_has_protac_tractability"] == 1)
    ).astype(int)

    df["biologic_druggable_label"] = (
        (df["ot_has_antibody_tractability"] == 1)
        | (df["ot_has_other_modality_tractability"] == 1)
    ).astype(int)

    df["dgidb_interaction_supported_label"] = (
        (df["dgidb_has_interaction"] == 1)
        | (df["dgidb_num_interactions"] > 0)
    ).astype(int)

    df["potentially_druggable_category_label"] = (
        df["dgidb_has_any_category"] == 1
    ).astype(int)

    df["final_any_druggable_label"] = (
        (df["clinical_target_label"] == 1)
        | (df["clinical_investigation_label"] == 1)
        | (df["chemical_tractable_label"] == 1)
        | (df["small_molecule_druggable_label"] == 1)
        | (df["biologic_druggable_label"] == 1)
        | (df["dgidb_interaction_supported_label"] == 1)
        | (df["potentially_druggable_category_label"] == 1)
    ).astype(int)

    df["unknown_dark_label"] = (
        (df["final_any_druggable_label"] == 0)
        | (df["pharos_is_tdark"] == 1)
    ).astype(int)

    score = pd.Series(0.0, index=df.index)

    score += 5.0 * df["chembl_has_phase4_or_approved"]
    score += 5.0 * df["pharos_is_tclin"]

    score += 4.0 * (((df["chembl_max_phase"] >= 1) & (df["chembl_max_phase"] < 4)).astype(int))

    score += 3.0 * df["ot_has_clinical_indication"]
    score += 3.0 * df["chembl_has_mechanism"]

    score += 2.0 * df["chembl_has_potent_compound"]
    score += 2.0 * df["pharos_is_tchem"]
    score += 2.0 * df["dgidb_interaction_supported_label"]
    score += 2.0 * df["ot_has_small_molecule_tractability"]
    score += 2.0 * df["ot_has_protac_tractability"]

    score += 1.5 * df["ot_has_antibody_tractability"]
    score += 1.0 * df["potentially_druggable_category_label"]

    score += 0.5 * df["pharos_is_tbio"]
    score -= 1.0 * df["pharos_is_tdark"]

    df["druggability_score"] = score.round(3)

    def assign_class(row: pd.Series) -> str:
        if int(row["clinical_target_label"]) == 1:
            return "Clinically_drugged"
        if int(row["clinical_investigation_label"]) == 1:
            return "Clinically_investigated"
        if int(row["chemical_tractable_label"]) == 1:
            return "Chemically_tractable"
        if int(row["biologic_druggable_label"]) == 1:
            return "Biologic_or_modality_tractable"
        if int(row["dgidb_interaction_supported_label"]) == 1:
            return "Drug_gene_interaction_supported"
        if int(row["potentially_druggable_category_label"]) == 1:
            return "Potentially_druggable_category"
        if int(row["pharos_is_tbio"]) == 1:
            return "Biologically_supported_only"
        return "Unknown_or_dark"

    df["druggability_class"] = df.apply(assign_class, axis=1)

    def evidence_sources(row: pd.Series) -> str:
        sources = []

        if int(row.get("chembl_has_target", 0)) == 1:
            sources.append("ChEMBL_target")
        if int(row.get("chembl_has_phase4_or_approved", 0)) == 1:
            sources.append("ChEMBL_approved_phase4")
        if int(row.get("chembl_has_mechanism", 0)) == 1:
            sources.append("ChEMBL_mechanism")
        if int(row.get("chembl_has_potent_compound", 0)) == 1:
            sources.append("ChEMBL_potent_activity")

        if int(row.get("opentargets_has_target", 0)) == 1:
            sources.append("OpenTargets_target")
        if int(row.get("ot_has_clinical_indication", 0)) == 1:
            sources.append("OpenTargets_clinical_indication")
        if int(row.get("ot_has_small_molecule_tractability", 0)) == 1:
            sources.append("OpenTargets_small_molecule_tractability")
        if int(row.get("ot_has_antibody_tractability", 0)) == 1:
            sources.append("OpenTargets_antibody_tractability")
        if int(row.get("ot_has_protac_tractability", 0)) == 1:
            sources.append("OpenTargets_PROTAC_tractability")

        if int(row.get("dgidb_has_interaction", 0)) == 1:
            sources.append("DGIdb_interaction")
        if int(row.get("dgidb_has_any_category", 0)) == 1:
            sources.append("DGIdb_category")

        if int(row.get("pharos_api_available", 0)) == 1 and str(row.get("pharos_tdl", "")).strip() not in ["", "nan", "None"]:
            sources.append("Pharos_API_TDL")

        if not sources:
            sources.append("No_positive_evidence")

        return "|".join(sources)

    df["evidence_sources"] = df.apply(evidence_sources, axis=1)
    df["evidence_source_count"] = df["evidence_sources"].apply(
        lambda x: 0 if x == "No_positive_evidence" else len(str(x).split("|"))
    )

    def confidence(row: pd.Series) -> str:
        if int(row["clinical_target_label"]) == 1:
            return "high"
        if int(row["clinical_investigation_label"]) == 1:
            return "medium_high"
        if int(row["chemical_tractable_label"]) == 1 and int(row["evidence_source_count"]) >= 2:
            return "medium_high"
        if int(row["chemical_tractable_label"]) == 1:
            return "medium"
        if int(row["biologic_druggable_label"]) == 1:
            return "medium"
        if int(row["dgidb_interaction_supported_label"]) == 1:
            return "medium"
        if int(row["potentially_druggable_category_label"]) == 1:
            return "low_medium"
        return "low_or_unknown"

    df["label_confidence"] = df.apply(confidence, axis=1)
    df["do_not_use_evidence_columns_as_model_features"] = 1

    return df


# =============================================================================
# Output
# =============================================================================

def build_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []

    cols = [
        "clinical_target_label",
        "clinical_investigation_label",
        "chemical_tractable_label",
        "small_molecule_druggable_label",
        "biologic_druggable_label",
        "dgidb_interaction_supported_label",
        "potentially_druggable_category_label",
        "unknown_dark_label",
        "final_any_druggable_label",
    ]

    for c in cols:
        if c in df.columns:
            n = int(pd.to_numeric(df[c], errors="coerce").fillna(0).sum())
            rows.append({
                "metric": c,
                "count": n,
                "total": len(df),
                "percent": round(100 * n / max(len(df), 1), 3),
            })

    for col in ["druggability_class", "label_confidence", "pharos_tdl"]:
        if col in df.columns:
            for cls, n in df[col].value_counts(dropna=False).items():
                rows.append({
                    "metric": f"{col}={cls}",
                    "count": int(n),
                    "total": len(df),
                    "percent": round(100 * int(n) / max(len(df), 1), 3),
                })

    return pd.DataFrame(rows)


def reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    preferred = [
        "gene_symbol", "gene_name", "hgnc_id", "ensembl_gene_id", "entrez_gene_id",
        "uniprot_ids", "uniprot_id_primary", "locus_type", "locus_group", "status", "location",

        "druggability_class", "druggability_score", "final_any_druggable_label",
        "clinical_target_label", "clinical_investigation_label", "chemical_tractable_label",
        "small_molecule_druggable_label", "biologic_druggable_label",
        "dgidb_interaction_supported_label", "potentially_druggable_category_label",
        "unknown_dark_label", "label_confidence", "evidence_sources", "evidence_source_count",

        "chembl_has_target", "chembl_target_chembl_ids", "chembl_target_names",
        "chembl_num_targets", "chembl_num_mechanisms", "chembl_num_drugs",
        "chembl_max_phase", "chembl_has_mechanism", "chembl_has_phase4_or_approved",
        "chembl_num_activity_records", "chembl_num_potent_activity_records",
        "chembl_max_pchembl", "chembl_has_potent_compound", "chembl_drug_names",

        "opentargets_has_target", "opentargets_approved_symbol", "opentargets_approved_name",
        "ot_has_clinical_indication", "ot_num_clinical_indications", "ot_num_clinical_drugs",
        "ot_num_clinical_diseases", "ot_max_clinical_phase",
        "ot_has_small_molecule_tractability", "ot_has_antibody_tractability",
        "ot_has_protac_tractability", "ot_has_other_modality_tractability",

        "dgidb_has_interaction", "dgidb_num_interactions", "dgidb_num_drugs",
        "dgidb_drug_names", "dgidb_interaction_types", "dgidb_sources",
        "dgidb_has_gene_record", "dgidb_has_any_category",
        "dgidb_categories_from_categories_file",

        "pharos_api_available", "pharos_api_status", "pharos_tdl",
        "pharos_name", "pharos_sym", "pharos_fam", "pharos_novelty",
        "pharos_description", "pharos_is_tclin", "pharos_is_tchem",
        "pharos_is_tbio", "pharos_is_tdark",

        "do_not_use_evidence_columns_as_model_features",
        "gene_symbol_norm",
    ]

    existing = [c for c in preferred if c in df.columns]
    remaining = [c for c in df.columns if c not in existing]
    return df[existing + remaining].copy()


# =============================================================================
# CLI
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build human gene druggability ground-truth table.")

    p.add_argument("--databases-dir", default="databases")
    p.add_argument("--output-dir", default="Step0_Output")
    p.add_argument("--opentargets-release", default="26.03")
    p.add_argument("--all-locus-types", action="store_true")
    p.add_argument("--save-raw-evidence", action="store_true")

    # Pharos API is ON by default.
    p.add_argument("--skip-pharos-api", action="store_true", help="Disable Pharos API.")
    p.add_argument("--pharos-sleep", type=float, default=0.0)
    p.add_argument("--pharos-timeout", type=int, default=60)
    p.add_argument("--pharos-workers", type=int, default=16)
    p.add_argument("--pharos-save-every", type=int, default=500)
    p.add_argument("--force-pharos-api", action="store_true")
    p.add_argument("--max-pharos-genes", type=int, default=None)

    return p.parse_args()


def main() -> None:
    args = parse_args()

    databases_dir = Path(args.databases_dir)
    output_dir = Path(args.output_dir)
    raw_dir = output_dir / "raw_evidence"

    safe_mkdir(output_dir)
    safe_mkdir(raw_dir)

    use_pharos_api = not args.skip_pharos_api

    print("=" * 100)
    print("BUILD HUMAN GENE DRUGGABILITY GROUND-TRUTH TABLE")
    print("=" * 100)
    print(f"[DATABASES DIR]        {databases_dir.resolve()}")
    print(f"[OUTPUT DIR]           {output_dir.resolve()}")
    print(f"[OPEN TARGETS RELEASE] {args.opentargets_release}")
    print(f"[PROTEIN CODING ONLY]  {not args.all_locus_types}")
    print(f"[USE PHAROS API]       {use_pharos_api}")
    print(f"[PHAROS WORKERS]       {args.pharos_workers}")
    print("=" * 100)

    hgnc = load_hgnc(databases_dir, protein_coding_only=not args.all_locus_types)
    hgnc.to_csv(output_dir / "01_gene_universe_protein_coding.csv", index=False)

    try:
        chembl = build_chembl_gene_evidence(databases_dir, hgnc=hgnc)
    except Exception as e:
        print(f"[ERROR] ChEMBL evidence failed: {e}")
        chembl = pd.DataFrame()

    if args.save_raw_evidence and not chembl.empty:
        chembl.to_csv(raw_dir / "chembl_gene_evidence.csv", index=False)

    try:
        opentargets = build_opentargets_evidence(databases_dir, release=args.opentargets_release)
    except Exception as e:
        print(f"[ERROR] Open Targets evidence failed: {e}")
        opentargets = pd.DataFrame()

    if args.save_raw_evidence and not opentargets.empty:
        opentargets.to_csv(raw_dir / "opentargets_gene_evidence.csv", index=False)

    try:
        dgidb = build_dgidb_evidence(databases_dir)
    except Exception as e:
        print(f"[ERROR] DGIdb evidence failed: {e}")
        dgidb = pd.DataFrame()

    if args.save_raw_evidence and not dgidb.empty:
        dgidb.to_csv(raw_dir / "dgidb_gene_evidence.csv", index=False)

    if use_pharos_api:
        try:
            pharos = build_pharos_api_evidence(
                hgnc=hgnc,
                output_dir=output_dir,
                sleep_seconds=args.pharos_sleep,
                timeout=args.pharos_timeout,
                force=args.force_pharos_api,
                max_genes=args.max_pharos_genes,
                workers=args.pharos_workers,
                save_every=args.pharos_save_every,
            )
        except Exception as e:
            print(f"[ERROR] Pharos API evidence failed: {e}")
            pharos = pd.DataFrame()
    else:
        pharos = pd.DataFrame()

    if args.save_raw_evidence and not pharos.empty:
        pharos.to_csv(raw_dir / "pharos_api_tdl_evidence.csv", index=False)

    final_df = merge_evidence(
        hgnc=hgnc,
        chembl=chembl,
        opentargets=opentargets,
        dgidb=dgidb,
        pharos=pharos,
    )

    final_df = reorder_columns(final_df)

    final_path = output_dir / "03_HumanGene_DruggabilityLabels.csv"
    summary_path = output_dir / "04_label_summary.csv"

    final_df.to_csv(final_path, index=False)

    summary = build_summary(final_df)
    summary.to_csv(summary_path, index=False)

    meta = {
        "databases_dir": str(databases_dir.resolve()),
        "output_dir": str(output_dir.resolve()),
        "opentargets_release": args.opentargets_release,
        "protein_coding_only": not args.all_locus_types,
        "use_pharos_api": use_pharos_api,
        "pharos_sleep": args.pharos_sleep,
        "pharos_timeout": args.pharos_timeout,
        "max_pharos_genes": args.max_pharos_genes,
        "n_genes": int(len(final_df)),
        "n_columns": int(len(final_df.columns)),
        "final_table": str(final_path),
        "summary": str(summary_path),
        "important_note": "This is a ground-truth/evidence table. Do not use evidence columns as model features.",
    }
    save_json(output_dir / "00_build_metadata.json", meta)

    print("\n" + "=" * 100)
    print("DONE")
    print("=" * 100)
    print(f"[FINAL TABLE] {final_path}")
    print(f"[SUMMARY]     {summary_path}")
    print(f"[ROWS]        {len(final_df)}")
    print(f"[COLUMNS]     {len(final_df.columns)}")
    print("\nLABEL SUMMARY")
    print("-" * 100)
    print(summary.to_string(index=False))
    print("=" * 100)


if __name__ == "__main__":
    warnings.simplefilter(action="ignore", category=FutureWarning)
    main()
