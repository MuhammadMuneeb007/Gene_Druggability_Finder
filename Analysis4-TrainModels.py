#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Analysis4-TrainModels.py

Train classification models for gene druggability prediction.

Fixed version:
- Removes non-model identifier columns such as PDB ID.
- Removes known broken CTD -inf log feature.
- Converts all selected features to numeric.
- Replaces inf / -inf / too-large float32 values with NaN.
- Imputes missing values inside CV folds using Pipeline.
- Writes clean output files expected by Analysis5-ListPerformance.py.
"""

from __future__ import annotations

import argparse
import json
import platform
import re
import sys
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from scipy.stats import ttest_1samp

from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    cohen_kappa_score,
    confusion_matrix,
    f1_score,
    log_loss,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline

warnings.filterwarnings("ignore")

try:
    import joblib
    HAS_JOBLIB = True
except ImportError:
    HAS_JOBLIB = False

try:
    from xgboost import XGBClassifier
    HAS_XGBOOST = True
except ImportError:
    HAS_XGBOOST = False


# =============================================================================
# SETTINGS
# =============================================================================

DATASETS_DIR    = Path("Datasets")
RANDOM_SEED     = 42
N_CV_FOLDS      = 5
N_BOOTSTRAP     = 2000
N_PERMUTATIONS  = 100
MIN_POSITIVES   = 10
SAVE_MODELS     = True
TOP_N_FEATURES  = 50

FLOAT32_MAX = np.finfo(np.float32).max


# =============================================================================
# FEATURE CLEANING SETTINGS
# =============================================================================

# IMPORTANT:
# Keep the original Analysis4 pipeline exactly the same.
# Only remove unsafe / leakage / identifier / metadata columns before model training.
# Do NOT remove whole biological feature blocks here.

DROP_FEATURE_EXACT = {
    # -------------------------------------------------------------------------
    # Confirmed identifier leakage: these were ranked at the top and must not
    # be used as numeric model inputs.
    # -------------------------------------------------------------------------
    "Feature1_DepMap_depmap_genecsv_entrez_id",
    "Feature1_DepMap_depmap_entrez_id",

    # -------------------------------------------------------------------------
    # Non-biological structural identifier.
    # -------------------------------------------------------------------------
    "Feature4_Structure_structure_best_pdb_id",

    # -------------------------------------------------------------------------
    # Broken CTD feature with inf / -inf risk.
    # -------------------------------------------------------------------------
    "Feature17_CTD_ctd_directional_rows_increase_minus_decrease_log1p",

    # -------------------------------------------------------------------------
    # Ensembl genomic coordinates.
    # -------------------------------------------------------------------------
    "Feature8_Ensembl_ensembl_gene_start",
    "Feature8_Ensembl_ensembl_gene_end",

    # -------------------------------------------------------------------------
    # STRING graph-level / mapping artifacts, not gene-specific biology.
    # -------------------------------------------------------------------------
    "Feature2_String_string_graph_component_rank",
    "Feature2_String_string_graph_component_size",
    "Feature2_String_string_graph_component_fraction",
    "Feature2_String_string_graph_is_largest_component",
    "Feature2_String_string_graph_total_edges",
    "Feature2_String_string_graph_total_nodes",
    "Feature2_String_string_partner_limit",
    "Feature2_String_string_required_score",
    "Feature2_String_string_mapped",
    "Feature2_String_string_has_ppi_feature",

    # -------------------------------------------------------------------------
    # Embedding metadata / missingness flags, not biology.
    # -------------------------------------------------------------------------
    "Feature12_Embeddings_embedding_dim",
    "Feature12_Embeddings_embedding_protein_count",
    "Feature12_Embeddings_embedding_has_embedding",

    # -------------------------------------------------------------------------
    # Annotation history / literature time proxies.
    # -------------------------------------------------------------------------
    "Feature14_GO_go_first_annotation_year",
    "Feature14_GO_go_last_annotation_year",
    "Feature14_GO_go_annotation_year_span",
    "Feature13_GWAS_gwas_first_association_year",
    "Feature13_GWAS_gwas_last_association_year",
    "Feature13_GWAS_gwas_association_year_span",

    # -------------------------------------------------------------------------
    # Publication / database attention proxies.
    # -------------------------------------------------------------------------
    "Feature18_MGI_mgi_unique_pubmed_count",
    "Feature18_MGI_mgi_pubmed_mention_count",
    "Feature14_GO_go_annotation_count",
    "Feature14_GO_go_log1p_annotation_count",
    "Feature14_GO_go_unique_assigned_by_count",
    "Feature15_BioGRID_biogrid_unique_source_database_count",

    # -------------------------------------------------------------------------
    # gnomAD row/transcript count proxies.
    # -------------------------------------------------------------------------
    "Feature19_gnomADFull_gnomad_constraint_row_count",
    "Feature19_gnomADFull_gnomad_unique_transcript_count",

    # -------------------------------------------------------------------------
    # Direct druggability proxy.
    # -------------------------------------------------------------------------
    "Feature16_HPA_hpa_classical_druggability_proxy_score",

    # -------------------------------------------------------------------------
    # CORUM drug-target complex leakage features.
    # These encode drug-target-complex evidence, so keep generic CORUM complex
    # biology features but remove these drug-target-specific columns.
    # -------------------------------------------------------------------------
    "Feature20_CORUM_corum_drug_target_complex_complex_membership_count",
    "Feature20_CORUM_corum_drug_target_complex_complex_membership_count_log1p",
    "Feature20_CORUM_corum_drug_target_complex_complex_name_count",
    "Feature20_CORUM_corum_drug_target_complex_complex_name_count_log1p",
    "Feature20_CORUM_corum_drug_target_complex_complex_size_max",
    "Feature20_CORUM_corum_drug_target_complex_complex_size_mean",
    "Feature20_CORUM_corum_drug_target_complex_complex_size_median",
    "Feature20_CORUM_corum_drug_target_complex_complex_size_min",
    "Feature20_CORUM_corum_drug_target_complex_has_any_membership",
    "Feature20_CORUM_corum_drug_target_complex_has_any_membership_log1p",
    "Feature20_CORUM_corum_drug_target_complex_large_gt10_complex_count",
    "Feature20_CORUM_corum_drug_target_complex_large_gt10_complex_count_log1p",
    "Feature20_CORUM_corum_drug_target_complex_medium_4_10_complex_count",
    "Feature20_CORUM_corum_drug_target_complex_medium_4_10_complex_count_log1p",
    "Feature20_CORUM_corum_drug_target_complex_pubmed_count",
    "Feature20_CORUM_corum_drug_target_complex_pubmed_count_log1p",
    "Feature20_CORUM_corum_drug_target_complex_small_2_3_complex_count",
    "Feature20_CORUM_corum_drug_target_complex_small_2_3_complex_count_log1p",
}

# Keep this pattern list conservative.
# Do NOT use broad patterns like "_id" because that would incorrectly remove
# valid biological features such as paralogue_identity_mean.
DROP_FEATURE_PATTERNS = [
    "_pdb_id",
    "entrez_id",
    "gene_start",
    "gene_end",
    "graph_total_nodes",
    "graph_total_edges",
    "graph_component_rank",
    "graph_component_size",
    "graph_component_fraction",
    "graph_is_largest",
    "partner_limit",
    "required_score",
    "embedding_dim",
    "_raw",
    "_text",
    "_url",
    "_file",
    "_path",
    "identifier",
    "accession",
]

LEAKAGE_TERMS = [
    "chembl", "opentarget", "open_target", "dgidb", "drugbank",
    "pharos", "tcrd", "druggability", "drug_target", "tractability",
    "clinical_target", "label", "evidence", "known_drug",
]


# =============================================================================
# HELPERS
# =============================================================================

def log(msg: str) -> None:
    print(msg, flush=True)


def section(title: str) -> None:
    print(flush=True)
    print("=" * 100, flush=True)
    print(title, flush=True)
    print("=" * 100, flush=True)


def subsection(title: str) -> None:
    print(flush=True)
    print("-" * 80, flush=True)
    print(title, flush=True)
    print("-" * 80, flush=True)


def now_iso() -> str:
    return datetime.now().isoformat()


def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def safe_float(x: Any, default: float = np.nan) -> float:
    try:
        v = float(x)
        return v if np.isfinite(v) else default
    except Exception:
        return default


def normalize_symbol(x: Any) -> str:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return ""
    return re.sub(r"\s+", "", str(x).strip().upper())


def stars_from_p(p: Any) -> str:
    p = safe_float(p)
    if np.isnan(p):
        return "NA"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def save_csv(df: pd.DataFrame, path: Path) -> None:
    if df is None or df.empty:
        log(f"  [SKIP empty] {path}")
        return
    df.to_csv(path, index=False)
    log(f"  [SAVED] {path}")


# =============================================================================
# FIND DATASET FOLDER
# =============================================================================

def find_dataset_folder(dataset_id: int) -> Path:
    catalogue_path = DATASETS_DIR / "dataset_catalogue.csv"

    if catalogue_path.exists():
        cat = pd.read_csv(catalogue_path, low_memory=False)
        row = cat[cat["dataset_id"] == dataset_id]
        if not row.empty:
            name = str(row.iloc[0]["dataset_name"])
            folder = DATASETS_DIR / name
            if folder.exists():
                return folder

    pattern = f"Dataset{dataset_id:03d}_*"
    matches = sorted(DATASETS_DIR.glob(pattern))
    if matches:
        return matches[0]

    raise FileNotFoundError(
        f"Dataset {dataset_id} not found in {DATASETS_DIR}. "
        f"Expected folder matching {pattern}."
    )


# =============================================================================
# LOAD DATA WITH HARD CLEANING
# =============================================================================

def load_dataset(dataset_dir: Path) -> Tuple[
    pd.DataFrame, pd.Series, pd.Series, str, List[str], Dict
]:
    x_path    = dataset_dir / "X_train.csv"
    y_path    = dataset_dir / "y_train.csv"
    meta_path = dataset_dir / "meta.json"

    if not x_path.exists():
        raise FileNotFoundError(f"X_train.csv not found: {x_path}")
    if not y_path.exists():
        raise FileNotFoundError(f"y_train.csv not found: {y_path}")

    X_df = pd.read_csv(x_path, low_memory=False)
    y_df = pd.read_csv(y_path, low_memory=False)

    meta: Dict[str, Any] = {}
    if meta_path.exists():
        with open(meta_path, "r") as f:
            meta = json.load(f)

    if "gene_symbol" not in X_df.columns:
        raise RuntimeError(f"gene_symbol column missing from {x_path}")
    if "gene_symbol" not in y_df.columns:
        raise RuntimeError(f"gene_symbol column missing from {y_path}")

    X_df["gene_symbol"] = X_df["gene_symbol"].apply(normalize_symbol)
    y_df["gene_symbol"] = y_df["gene_symbol"].apply(normalize_symbol)

    merged = X_df.merge(y_df, on="gene_symbol", how="inner")
    merged = merged[merged["gene_symbol"] != ""].drop_duplicates("gene_symbol")

    y_cols = [c for c in y_df.columns if c != "gene_symbol"]
    if not y_cols:
        raise RuntimeError("y_train.csv has no target column.")
    target_col = y_cols[0]

    raw_feature_cols = [
        c for c in X_df.columns
        if c != "gene_symbol" and re.match(r"^Feature\d+_", c)
    ]

    if not raw_feature_cols:
        raise RuntimeError("No Feature{N}_* columns found in X_train.csv.")

    dropped_exact = [c for c in raw_feature_cols if c in DROP_FEATURE_EXACT]

    dropped_pattern = [
        c for c in raw_feature_cols
        if c not in DROP_FEATURE_EXACT
        and any(pat in c.lower() for pat in DROP_FEATURE_PATTERNS)
    ]

    feature_cols = [
        c for c in raw_feature_cols
        if c not in DROP_FEATURE_EXACT
        and not any(pat in c.lower() for pat in DROP_FEATURE_PATTERNS)
    ]

    if not feature_cols:
        raise RuntimeError("No usable feature columns after removing bad/non-numeric columns.")

    for c in feature_cols:
        merged[c] = pd.to_numeric(merged[c], errors="coerce")

    bad_before_rows = []
    for c in feature_cols:
        s = merged[c]
        n_inf = int(np.isinf(s).sum())
        n_big = int((s.abs() > FLOAT32_MAX).sum())
        if n_inf > 0 or n_big > 0:
            bad_before_rows.append({
                "feature": c,
                "n_inf_before_cleaning": n_inf,
                "n_too_large_float32_before_cleaning": n_big,
                "max_abs_before_cleaning": safe_float(s.abs().replace([np.inf, -np.inf], np.nan).max()),
            })

    merged[feature_cols] = merged[feature_cols].replace([np.inf, -np.inf], np.nan)
    merged[feature_cols] = merged[feature_cols].mask(
        merged[feature_cols].abs() > FLOAT32_MAX,
        np.nan,
    )

    merged[target_col] = (
        pd.to_numeric(merged[target_col], errors="coerce")
        .fillna(0)
        .astype(int)
    )

    all_nan_cols = [c for c in feature_cols if merged[c].isna().all()]
    valid_cols = [c for c in feature_cols if c not in all_nan_cols]

    if not valid_cols:
        raise RuntimeError("No usable feature columns after cleaning.")

    X     = merged[valid_cols].copy()
    y     = merged[target_col].copy()
    genes = merged["gene_symbol"].copy()

    arr = X.to_numpy(dtype=float)
    n_inf_final = int(np.isinf(arr).sum())
    n_big_final = int((np.abs(arr) > FLOAT32_MAX).sum())

    if n_inf_final > 0 or n_big_final > 0:
        raise RuntimeError(
            f"Bad values remain after cleaning: "
            f"n_inf={n_inf_final}, n_too_large_float32={n_big_final}"
        )

    cleaning_report = {
        "n_raw_feature_cols": len(raw_feature_cols),
        "n_dropped_exact_bad_features": len(dropped_exact),
        "dropped_exact_bad_features": dropped_exact,
        "n_dropped_pattern_features": len(dropped_pattern),
        "dropped_pattern_features": dropped_pattern,
        "n_bad_value_features_before_cleaning": len(bad_before_rows),
        "bad_value_features_before_cleaning": bad_before_rows,
        "n_all_nan_cols_after_cleaning": len(all_nan_cols),
        "all_nan_cols_after_cleaning": all_nan_cols,
        "n_final_feature_cols": len(valid_cols),
        "float32_max_threshold": float(FLOAT32_MAX),
    }

    meta["training_feature_cleaning"] = cleaning_report

    log(f"  Genes:              {len(X)}")
    log(f"  Raw feature columns:{len(raw_feature_cols)}")
    log(f"  Dropped exact bad:  {len(dropped_exact)}")
    for c in dropped_exact:
        log(f"    [DROP exact] {c}")
    log(f"  Dropped pattern:    {len(dropped_pattern)}")
    for c in dropped_pattern[:20]:
        log(f"    [DROP pattern] {c}")
    if len(dropped_pattern) > 20:
        log(f"    ... plus {len(dropped_pattern) - 20} more")
    log(f"  Bad-value features: {len(bad_before_rows)}")
    for r in bad_before_rows[:20]:
        log(
            f"    [CLEAN values] {r['feature']} "
            f"inf={r['n_inf_before_cleaning']} "
            f"big={r['n_too_large_float32_before_cleaning']}"
        )
    log(f"  All-NaN dropped:    {len(all_nan_cols)}")
    log(f"  Feature columns:    {len(valid_cols)}")
    log(f"  Target column:      {target_col}")
    log(f"  Positives (y=1):    {int(y.sum())}")
    log(f"  Negatives (y=0):    {int((y == 0).sum())}")
    log(f"  Prevalence:         {y.mean():.4f}")

    return X, y, genes, target_col, valid_cols, meta


# =============================================================================
# LEAKAGE AUDIT
# =============================================================================

def audit_features(feature_cols: List[str]) -> pd.DataFrame:
    rows = []
    for col in feature_cols:
        c = col.lower()
        matched = [term for term in LEAKAGE_TERMS if term in c]
        rows.append({
            "feature": col,
            "possible_leakage": bool(matched),
            "reason": ";".join(matched),
        })

    df = pd.DataFrame(rows)
    n_leaky = int(df["possible_leakage"].sum()) if not df.empty else 0

    if n_leaky > 0:
        log(f"  [LEAKAGE WARNING] {n_leaky} suspicious columns detected:")
        for _, r in df[df["possible_leakage"]].iterrows():
            log(f"    {r['feature']}  matched={r['reason']}")
    else:
        log("  [LEAKAGE AUDIT] Passed — no suspicious columns found.")

    return df


# =============================================================================
# MODEL DEFINITIONS
# =============================================================================

def build_models() -> Dict[str, Pipeline]:
    models: Dict[str, Pipeline] = {}

    models["RandomForest"] = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("model", RandomForestClassifier(
            n_estimators=500,
            max_depth=None,
            min_samples_leaf=3,
            max_features="sqrt",
            class_weight="balanced_subsample",
            n_jobs=-1,
            random_state=RANDOM_SEED,
        )),
    ])

    if HAS_XGBOOST:
        models["XGBoost"] = Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("model", XGBClassifier(
                n_estimators=500,
                max_depth=4,
                learning_rate=0.03,
                subsample=0.85,
                colsample_bytree=0.85,
                reg_lambda=2.0,
                reg_alpha=0.0,
                objective="binary:logistic",
                eval_metric="logloss",
                n_jobs=-1,
                random_state=RANDOM_SEED,
                verbosity=0,
            )),
        ])
    else:
        log("  [INFO] xgboost not installed — XGBoost skipped.")

    return models


# =============================================================================
# METRICS
# =============================================================================

def compute_classification_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    y_prob: np.ndarray,
) -> Dict[str, Any]:
    y_true = np.asarray(y_true, dtype=int)
    y_pred = np.asarray(y_pred, dtype=int)
    y_prob = np.asarray(y_prob, dtype=float)

    n = len(y_true)
    n_pos = int(y_true.sum())
    n_neg = int((y_true == 0).sum())

    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel() if cm.shape == (2, 2) else (0, 0, 0, 0)

    out = {
        "n_test": n,
        "n_pos_test": n_pos,
        "n_neg_test": n_neg,
        "tp": int(tp),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
    }

    out["auroc"] = safe_float(roc_auc_score(y_true, y_prob)) if n_pos > 0 and n_neg > 0 else np.nan
    out["auprc"] = safe_float(average_precision_score(y_true, y_prob)) if n_pos > 0 else np.nan
    out["accuracy"] = safe_float(accuracy_score(y_true, y_pred))
    out["balanced_accuracy"] = safe_float(balanced_accuracy_score(y_true, y_pred))
    out["f1_binary"] = safe_float(f1_score(y_true, y_pred, average="binary", zero_division=0))
    out["f1_macro"] = safe_float(f1_score(y_true, y_pred, average="macro", zero_division=0))
    out["f1_weighted"] = safe_float(f1_score(y_true, y_pred, average="weighted", zero_division=0))
    out["precision"] = safe_float(precision_score(y_true, y_pred, zero_division=0))
    out["recall"] = safe_float(recall_score(y_true, y_pred, zero_division=0))
    out["sensitivity"] = out["recall"]
    out["specificity"] = safe_float(tn / (tn + fp)) if (tn + fp) > 0 else np.nan
    out["ppv"] = out["precision"]
    out["npv"] = safe_float(tn / (tn + fn)) if (tn + fn) > 0 else np.nan
    out["mcc"] = safe_float(matthews_corrcoef(y_true, y_pred)) if len(np.unique(y_pred)) > 1 else np.nan
    out["cohen_kappa"] = safe_float(cohen_kappa_score(y_true, y_pred)) if len(np.unique(y_pred)) > 1 else np.nan
    out["brier_score"] = safe_float(brier_score_loss(y_true, y_prob))
    out["log_loss"] = safe_float(log_loss(y_true, y_prob)) if n_pos > 0 and n_neg > 0 else np.nan

    return out


def bootstrap_ci(
    values: List[float],
    n_boot: int = N_BOOTSTRAP,
    seed: int = RANDOM_SEED,
    ci: float = 0.95,
) -> Tuple[float, float, float]:
    vals = pd.Series(values).dropna().to_numpy(dtype=float)

    if len(vals) == 0:
        return np.nan, np.nan, np.nan
    if len(vals) == 1:
        return float(vals[0]), float(vals[0]), float(vals[0])

    rng = np.random.default_rng(seed)
    boot = [
        np.mean(rng.choice(vals, size=len(vals), replace=True))
        for _ in range(n_boot)
    ]

    lo = float(np.percentile(boot, (1 - ci) / 2 * 100))
    hi = float(np.percentile(boot, (1 + ci) / 2 * 100))
    return float(np.mean(vals)), lo, hi


def summarise_fold_metrics(fold_rows: List[Dict[str, Any]]) -> Dict[str, float]:
    df = pd.DataFrame(fold_rows)
    out: Dict[str, float] = {}

    if df.empty:
        return out

    numeric_metrics = [
        "auroc", "auprc", "accuracy", "balanced_accuracy",
        "f1_binary", "f1_macro", "f1_weighted",
        "precision", "recall", "sensitivity", "specificity",
        "ppv", "npv", "mcc", "cohen_kappa",
        "brier_score", "log_loss",
    ]

    for m in numeric_metrics:
        if m not in df.columns:
            continue

        vals = df[m].dropna()
        if len(vals) == 0:
            continue

        mean_v, ci_lo, ci_hi = bootstrap_ci(vals.tolist())

        out[f"{m}_mean"] = mean_v
        out[f"{m}_std"] = float(vals.std(ddof=1)) if len(vals) > 1 else 0.0
        out[f"{m}_se"] = float(vals.std(ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else 0.0
        out[f"{m}_median"] = float(vals.median())
        out[f"{m}_min"] = float(vals.min())
        out[f"{m}_max"] = float(vals.max())
        out[f"{m}_ci95_lo"] = ci_lo
        out[f"{m}_ci95_hi"] = ci_hi

        null = 0.5 if m == "auroc" else 0.0
        if len(vals) >= 2 and float(vals.std(ddof=1)) > 0:
            try:
                _, p = ttest_1samp(vals.to_numpy(), popmean=null, alternative="greater")
                out[f"{m}_ttest_p_vs_null"] = float(p)
            except Exception:
                out[f"{m}_ttest_p_vs_null"] = np.nan
        else:
            out[f"{m}_ttest_p_vs_null"] = np.nan

    return out


# =============================================================================
# PROBABILITY EXTRACTION
# =============================================================================

def _get_prob(pipeline: Pipeline, X: np.ndarray) -> np.ndarray:
    if hasattr(pipeline, "predict_proba"):
        p = pipeline.predict_proba(X)
        return p[:, 1] if p.ndim == 2 else p.ravel()

    if hasattr(pipeline, "decision_function"):
        d = pipeline.decision_function(X).ravel()
        mn, mx = d.min(), d.max()
        return (d - mn) / (mx - mn + 1e-9)

    return pipeline.predict(X).astype(float)


# =============================================================================
# CROSS VALIDATION
# =============================================================================

def get_importance(
    pipeline: Pipeline,
    feature_cols: List[str],
) -> Optional[pd.DataFrame]:
    fitted = pipeline.named_steps["model"]

    if not hasattr(fitted, "feature_importances_"):
        return None

    vals = np.asarray(fitted.feature_importances_, dtype=float)

    if len(vals) != len(feature_cols):
        return None

    return pd.DataFrame({
        "feature": feature_cols,
        "importance": vals,
        "abs_importance": np.abs(vals),
        "importance_type": "native_importance",
    }).sort_values("abs_importance", ascending=False).reset_index(drop=True)


def run_cv(
    X: pd.DataFrame,
    y: pd.Series,
    genes: pd.Series,
    feature_cols: List[str],
    models: Dict[str, Pipeline],
    n_folds: int,
    dataset_name: str,
    target_col: str,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:

    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=RANDOM_SEED)

    X_np = X.to_numpy(dtype=float)
    y_np = y.to_numpy(dtype=int)
    genes_np = genes.to_numpy()

    all_preds: List[Dict[str, Any]] = []
    all_fold_met: List[Dict[str, Any]] = []
    all_importance: List[pd.DataFrame] = []
    all_cm: List[Dict[str, Any]] = []

    for model_name, model_template in models.items():
        subsection(f"Model: {model_name}")

        model_fold_rows: List[Dict[str, Any]] = []

        for fold_i, (tr_idx, te_idx) in enumerate(skf.split(X_np, y_np), start=1):
            X_tr, X_te = X_np[tr_idx], X_np[te_idx]
            y_tr, y_te = y_np[tr_idx], y_np[te_idx]

            clf = clone(model_template)

            if model_name == "XGBoost" and HAS_XGBOOST:
                n_neg_tr = int((y_tr == 0).sum())
                n_pos_tr = max(1, int(y_tr.sum()))
                clf.named_steps["model"].set_params(
                    scale_pos_weight=n_neg_tr / n_pos_tr
                )

            try:
                clf.fit(X_tr, y_tr)
                prob = _get_prob(clf, X_te)
                pred = (prob >= 0.5).astype(int)
            except Exception as e:
                log(f"  [FAILED] fold={fold_i} model={model_name}: {e}")
                continue

            metrics = compute_classification_metrics(y_te, pred, prob)
            metrics.update({
                "dataset_name": dataset_name,
                "target_col": target_col,
                "model": model_name,
                "fold": fold_i,
                "n_train": len(tr_idx),
                "n_test": len(te_idx),
            })

            model_fold_rows.append(metrics)
            all_fold_met.append(metrics)

            for gene, yt, yp, pr in zip(genes_np[te_idx], y_te, pred, prob):
                all_preds.append({
                    "dataset_name": dataset_name,
                    "model": model_name,
                    "fold": fold_i,
                    "gene_symbol": gene,
                    "y_true": int(yt),
                    "y_pred": int(yp),
                    "y_prob": float(pr),
                    "correct": int(yt == yp),
                })

            all_cm.append({
                "dataset_name": dataset_name,
                "model": model_name,
                "fold": fold_i,
                "tn": int(metrics["tn"]),
                "fp": int(metrics["fp"]),
                "fn": int(metrics["fn"]),
                "tp": int(metrics["tp"]),
            })

            imp_df = get_importance(clf, feature_cols)
            if imp_df is not None:
                imp_df["dataset_name"] = dataset_name
                imp_df["model"] = model_name
                imp_df["fold"] = fold_i
                all_importance.append(imp_df)

            log(
                f"  fold={fold_i}/{n_folds}"
                f"  AUROC={metrics['auroc']:.4f}"
                f"  AUPRC={metrics['auprc']:.4f}"
                f"  BalAcc={metrics['balanced_accuracy']:.4f}"
                f"  F1={metrics['f1_binary']:.4f}"
                f"  MCC={metrics['mcc']:.4f}"
                f"  tp={metrics['tp']} fp={metrics['fp']}"
                f"  tn={metrics['tn']} fn={metrics['fn']}"
            )

        summary = summarise_fold_metrics(model_fold_rows)
        log(
            f"\n  [{model_name}] MEAN"
            f"  AUROC={summary.get('auroc_mean', np.nan):.4f}"
            f"  ±{summary.get('auroc_std', np.nan):.4f}"
            f"  CI95=[{summary.get('auroc_ci95_lo', np.nan):.4f},"
            f"{summary.get('auroc_ci95_hi', np.nan):.4f}]"
            f"  AUPRC={summary.get('auprc_mean', np.nan):.4f}"
            f"  BalAcc={summary.get('balanced_accuracy_mean', np.nan):.4f}"
            f"  F1={summary.get('f1_binary_mean', np.nan):.4f}"
            f"  MCC={summary.get('mcc_mean', np.nan):.4f}"
        )

    preds_df = pd.DataFrame(all_preds)
    folds_df = pd.DataFrame(all_fold_met)
    imp_df = pd.concat(all_importance, ignore_index=True) if all_importance else pd.DataFrame()
    cm_df = pd.DataFrame(all_cm)

    return preds_df, folds_df, imp_df, cm_df


# =============================================================================
# PERMUTATION TEST
# =============================================================================

def permutation_test_auroc(
    X: pd.DataFrame,
    y: pd.Series,
    model_pipeline: Pipeline,
    observed_auroc: float,
    n_permutations: int = N_PERMUTATIONS,
    n_folds: int = N_CV_FOLDS,
    seed: int = RANDOM_SEED,
) -> Dict[str, Any]:

    rng = np.random.default_rng(seed)
    X_np = X.to_numpy(dtype=float)
    y_np = y.to_numpy(dtype=int)

    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    perm_aurocs: List[float] = []

    for _ in range(n_permutations):
        y_perm = rng.permutation(y_np)

        fold_aurocs = []
        for tr_idx, te_idx in skf.split(X_np, y_perm):
            try:
                clf = clone(model_pipeline)

                if HAS_XGBOOST and "model" in clf.named_steps:
                    if isinstance(clf.named_steps["model"], XGBClassifier):
                        n_neg_tr = int((y_perm[tr_idx] == 0).sum())
                        n_pos_tr = max(1, int(y_perm[tr_idx].sum()))
                        clf.named_steps["model"].set_params(
                            scale_pos_weight=n_neg_tr / n_pos_tr
                        )

                clf.fit(X_np[tr_idx], y_perm[tr_idx])
                prob = _get_prob(clf, X_np[te_idx])

                if len(np.unique(y_perm[te_idx])) < 2:
                    continue

                fold_aurocs.append(roc_auc_score(y_perm[te_idx], prob))
            except Exception:
                continue

        if fold_aurocs:
            perm_aurocs.append(float(np.mean(fold_aurocs)))

    if not perm_aurocs:
        return {
            "permutation_auroc_mean": np.nan,
            "permutation_auroc_std": np.nan,
            "permutation_p_value": np.nan,
            "n_permutations_run": 0,
            "observed_auroc": observed_auroc,
            "auroc_above_random": False,
            "statistically_significant": False,
        }

    perm_arr = np.asarray(perm_aurocs, dtype=float)
    p_val = float((np.sum(perm_arr >= observed_auroc) + 1) / (len(perm_arr) + 1))

    return {
        "permutation_auroc_mean": float(np.mean(perm_arr)),
        "permutation_auroc_std": float(np.std(perm_arr, ddof=1)) if len(perm_arr) > 1 else 0.0,
        "permutation_p_value": p_val,
        "n_permutations_run": len(perm_arr),
        "observed_auroc": observed_auroc,
        "auroc_above_random": bool(observed_auroc > np.mean(perm_arr)),
        "statistically_significant": bool(p_val < 0.05),
    }


# =============================================================================
# SUMMARIES
# =============================================================================

def build_comparison_summary(
    folds_df: pd.DataFrame,
    perm_results: Dict[str, Dict[str, Any]],
) -> pd.DataFrame:

    rows = []

    if folds_df.empty or "model" not in folds_df.columns:
        return pd.DataFrame()

    for model_name, g in folds_df.groupby("model"):
        row: Dict[str, Any] = {
            "model": model_name,
            "n_folds": int(g["fold"].nunique()) if "fold" in g.columns else len(g),
        }

        row.update(summarise_fold_metrics(g.to_dict("records")))

        perm = perm_results.get(model_name, {})
        for k, v in perm.items():
            row[f"perm_{k}"] = v

        p = perm.get("permutation_p_value", np.nan)
        row["permutation_p_value"] = safe_float(p)
        row["significance"] = stars_from_p(p)

        rows.append(row)

    out = pd.DataFrame(rows)

    if not out.empty and "auroc_mean" in out.columns:
        out = out.sort_values("auroc_mean", ascending=False).reset_index(drop=True)
        out["rank"] = range(1, len(out) + 1)

    return out


def summarise_importance(all_imp: pd.DataFrame) -> pd.DataFrame:
    if all_imp.empty:
        return pd.DataFrame()

    rows = []

    for (model, feat), g in all_imp.groupby(["model", "feature"]):
        vals = g["importance"].dropna()
        abs_vals = g["abs_importance"].dropna()

        mean_v, ci_lo, ci_hi = bootstrap_ci(vals.tolist())

        rows.append({
            "model": model,
            "feature": feat,
            "importance_type": g["importance_type"].iloc[0],
            "mean_importance": mean_v,
            "std_importance": float(vals.std(ddof=1)) if len(vals) > 1 else 0.0,
            "median_importance": float(vals.median()) if len(vals) else np.nan,
            "ci95_lo": ci_lo,
            "ci95_hi": ci_hi,
            "mean_abs_importance": float(abs_vals.mean()) if len(abs_vals) else np.nan,
            "max_abs_importance": float(abs_vals.max()) if len(abs_vals) else np.nan,
            "n_folds": int(g["fold"].nunique()),
        })

    out = pd.DataFrame(rows)
    out = out.sort_values(["model", "mean_abs_importance"], ascending=[True, False])
    out["rank"] = out.groupby("model")["mean_abs_importance"].rank(
        method="first",
        ascending=False,
    ).astype(int)

    return out.reset_index(drop=True)


def aggregate_cm(cm_df: pd.DataFrame) -> pd.DataFrame:
    if cm_df.empty:
        return pd.DataFrame()

    rows = []

    for model, g in cm_df.groupby("model"):
        tp = int(g["tp"].sum())
        tn = int(g["tn"].sum())
        fp = int(g["fp"].sum())
        fn = int(g["fn"].sum())
        total = tp + tn + fp + fn

        rows.append({
            "model": model,
            "tp": tp,
            "tn": tn,
            "fp": fp,
            "fn": fn,
            "total": total,
            "sensitivity": tp / (tp + fn) if (tp + fn) > 0 else np.nan,
            "specificity": tn / (tn + fp) if (tn + fp) > 0 else np.nan,
            "precision": tp / (tp + fp) if (tp + fp) > 0 else np.nan,
            "npv": tn / (tn + fn) if (tn + fn) > 0 else np.nan,
            "accuracy": (tp + tn) / total if total > 0 else np.nan,
            "f1": 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) > 0 else np.nan,
            "mcc": (
                (tp * tn - fp * fn) /
                np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
                if (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) > 0
                else np.nan
            ),
        })

    return pd.DataFrame(rows)


# =============================================================================
# FINAL MODEL SAVE
# =============================================================================

def save_final_models(
    X: pd.DataFrame,
    y: pd.Series,
    models: Dict[str, Pipeline],
    outdir: Path,
) -> List[str]:

    if not SAVE_MODELS:
        return []

    if not HAS_JOBLIB:
        log("  [INFO] joblib unavailable — models not saved.")
        return []

    model_dir = mkdir(outdir / "models")
    saved = []

    X_np = X.to_numpy(dtype=float)
    y_np = y.to_numpy(dtype=int)

    for model_name, model_template in models.items():
        clf = clone(model_template)

        if model_name == "XGBoost" and HAS_XGBOOST:
            n_neg = int((y_np == 0).sum())
            n_pos = max(1, int(y_np.sum()))
            clf.named_steps["model"].set_params(scale_pos_weight=n_neg / n_pos)

        try:
            clf.fit(X_np, y_np)
            path = model_dir / f"{model_name}.joblib"
            joblib.dump(clf, path)
            saved.append(str(path))
            log(f"  [SAVED MODEL] {path}")
        except Exception as e:
            log(f"  [MODEL SAVE FAILED] {model_name}: {e}")

    return saved


# =============================================================================
# PUBLICATION SUMMARY
# =============================================================================

def write_publication_summary(
    outdir: Path,
    dataset_name: str,
    target_col: str,
    meta: Dict[str, Any],
    comparison: pd.DataFrame,
    best_row: pd.Series,
    top_feat: pd.DataFrame,
) -> None:

    md_path = outdir / "14_publication_summary.md"
    tex_path = outdir / "15_publication_summary.tex"

    best_model = best_row.get("model", "NA")
    auroc = best_row.get("auroc_mean", np.nan)
    auprc = best_row.get("auprc_mean", np.nan)
    mcc = best_row.get("mcc_mean", np.nan)
    f1 = best_row.get("f1_binary_mean", np.nan)
    p = best_row.get("permutation_p_value", np.nan)

    with open(md_path, "w") as f:
        f.write(f"# {dataset_name} training summary\n\n")
        f.write(f"Generated: {now_iso()}\n\n")
        f.write("## Dataset\n\n")
        f.write(f"- Dataset: `{dataset_name}`\n")
        f.write(f"- Target: `{target_col}`\n")
        f.write(f"- Group: `{meta.get('group_name', '')}`\n")
        f.write(f"- Tier: `{meta.get('tier', '')}`\n\n")

        f.write("## Best model\n\n")
        f.write(f"- Best model: **{best_model}**\n")
        f.write(f"- AUROC: {safe_float(auroc):.4f}\n")
        f.write(f"- AUPRC: {safe_float(auprc):.4f}\n")
        f.write(f"- F1: {safe_float(f1):.4f}\n")
        f.write(f"- MCC: {safe_float(mcc):.4f}\n")
        f.write(f"- Permutation p-value: {safe_float(p):.4g}\n\n")

        f.write("## Model comparison\n\n")
        if not comparison.empty:
            try:
                f.write(comparison.to_markdown(index=False))
            except Exception:
                f.write(comparison.to_string(index=False))
        else:
            f.write("No comparison table available.\n")

        f.write("\n\n## Top features\n\n")
        if not top_feat.empty:
            try:
                f.write(top_feat.head(50).to_markdown(index=False))
            except Exception:
                f.write(top_feat.head(50).to_string(index=False))
        else:
            f.write("No feature importance available.\n")

    with open(tex_path, "w") as f:
        f.write(f"% {dataset_name} training summary\n")
        f.write(f"% Generated: {now_iso()}\n\n")

        if not comparison.empty:
            f.write("% Model comparison\n")
            try:
                f.write(comparison.to_latex(index=False, escape=True, longtable=True))
            except Exception:
                f.write("% Could not export model comparison to LaTeX.\n")

        if not top_feat.empty:
            f.write("\n% Top features\n")
            try:
                f.write(top_feat.head(50).to_latex(index=False, escape=True, longtable=True))
            except Exception:
                f.write("% Could not export top features to LaTeX.\n")

    log(f"  [SAVED] {md_path}")
    log(f"  [SAVED] {tex_path}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    global DATASETS_DIR, SAVE_MODELS

    parser = argparse.ArgumentParser(
        description="Train RandomForest and XGBoost on one druggability dataset."
    )
    parser.add_argument("dataset_id", type=int)
    parser.add_argument("--datasets-dir", default=str(DATASETS_DIR))
    parser.add_argument("--cv", type=int, default=N_CV_FOLDS)
    parser.add_argument("--no-permutation", action="store_true")
    parser.add_argument("--no-save-models", action="store_true")

    args = parser.parse_args()

    DATASETS_DIR = Path(args.datasets_dir)
    SAVE_MODELS = not args.no_save_models

    section(f"ANALYSIS4 — TRAIN DATASET {args.dataset_id}")
    log(f"Started:         {now_iso()}")
    log(f"Dataset ID:      {args.dataset_id}")
    log(f"Datasets dir:    {DATASETS_DIR}")
    log(f"CV folds:        {args.cv}")
    log(f"Permutations:    {'skipped' if args.no_permutation else N_PERMUTATIONS}")
    log(f"Save models:     {SAVE_MODELS}")
    log(f"XGBoost:         {'available' if HAS_XGBOOST else 'NOT installed'}")

    section("LOADING DATASET")
    dataset_dir = find_dataset_folder(args.dataset_id)
    log(f"  Dataset folder: {dataset_dir}")

    X, y, genes, target_col, feature_cols, meta = load_dataset(dataset_dir)

    n_pos = int(y.sum())
    n_neg = int((y == 0).sum())
    n_genes = len(y)
    prevalence = float(y.mean())
    dataset_name = dataset_dir.name

    if n_pos < MIN_POSITIVES:
        log(f"[ABORT] Only {n_pos} positives — minimum is {MIN_POSITIVES}.")
        sys.exit(1)

    outdir = mkdir(dataset_dir / "Training")

    section("LEAKAGE AUDIT")
    audit_df = audit_features(feature_cols)

    section("MODELS")
    models = build_models()
    for name in models:
        log(f"  {name}")

    section("CROSS-VALIDATION")
    preds_df, folds_df, imp_df, cm_df = run_cv(
        X=X,
        y=y,
        genes=genes,
        feature_cols=feature_cols,
        models=models,
        n_folds=args.cv,
        dataset_name=dataset_name,
        target_col=target_col,
    )

    perm_results: Dict[str, Dict[str, Any]] = {}

    if not args.no_permutation:
        section("PERMUTATION TESTS")

        if folds_df.empty or "model" not in folds_df.columns:
            log("  [SKIP] No successful CV folds — skipping permutation tests.")
        else:
            for model_name, model_template in models.items():
                model_folds = folds_df[folds_df["model"] == model_name]["auroc"].dropna()

                if len(model_folds) == 0:
                    log(f"  [SKIP] {model_name}: no successful folds.")
                    continue

                observed_auroc = float(model_folds.mean())

                log(
                    f"  Running permutation test for {model_name} "
                    f"observed AUROC={observed_auroc:.4f}"
                )

                perm_res = permutation_test_auroc(
                    X=X,
                    y=y,
                    model_pipeline=clone(model_template),
                    observed_auroc=observed_auroc,
                    n_permutations=N_PERMUTATIONS,
                    n_folds=args.cv,
                    seed=RANDOM_SEED,
                )

                perm_results[model_name] = perm_res

                log(
                    f"  {model_name}: random_mean="
                    f"{perm_res.get('permutation_auroc_mean', np.nan):.4f} "
                    f"p={perm_res.get('permutation_p_value', np.nan):.4g}"
                )

    section("MODEL COMPARISON")

    if folds_df.empty or "model" not in folds_df.columns:
        log("[ERROR] No successful CV folds for any model.")
        log("[ERROR] Writing failure marker and stopping cleanly.")

        failure_info = {
            "dataset_id": args.dataset_id,
            "dataset_name": dataset_name,
            "dataset_dir": str(dataset_dir),
            "target_col": target_col,
            "n_genes": n_genes,
            "n_positives": n_pos,
            "n_negatives": n_neg,
            "prevalence": prevalence,
            "n_features": len(feature_cols),
            "reason": "No successful CV folds for any model.",
            "run_datetime": now_iso(),
        }

        with open(outdir / "FAILED_no_successful_cv_folds.json", "w") as f:
            json.dump(failure_info, f, indent=2)

        pd.DataFrame([failure_info]).to_csv(
            outdir / "FAILED_no_successful_cv_folds.csv",
            index=False,
        )

        log(f"  [SAVED] {outdir / 'FAILED_no_successful_cv_folds.json'}")
        log(f"  [SAVED] {outdir / 'FAILED_no_successful_cv_folds.csv'}")
        return

    comparison = build_comparison_summary(folds_df, perm_results)

    if comparison.empty:
        raise RuntimeError("Comparison summary is empty despite successful folds.")

    best_row = comparison.iloc[0]
    best_model_name = str(best_row["model"])

    display_cols = [
        "rank", "model", "n_folds",
        "auroc_mean", "auroc_std", "auroc_ci95_lo", "auroc_ci95_hi",
        "auprc_mean", "balanced_accuracy_mean",
        "f1_binary_mean", "mcc_mean",
        "permutation_p_value", "significance",
    ]
    dc = [c for c in display_cols if c in comparison.columns]
    print(comparison[dc].to_string(index=False))

    section("FEATURE IMPORTANCE")
    imp_summary = summarise_importance(imp_df)

    if not imp_summary.empty:
        top_feat = (
            imp_summary[imp_summary["model"] == best_model_name]
            .sort_values("mean_abs_importance", ascending=False)
            .head(TOP_N_FEATURES)
            .copy()
        )
    else:
        top_feat = pd.DataFrame()

    if not top_feat.empty:
        print(
            top_feat[
                ["rank", "feature", "mean_importance", "std_importance", "ci95_lo", "ci95_hi"]
            ].head(20).to_string(index=False)
        )
    else:
        log("  [INFO] No feature importance available.")

    section("CONFUSION MATRIX")
    agg_cm = aggregate_cm(cm_df)
    if not agg_cm.empty:
        print(agg_cm.to_string(index=False))

    section("SAVING FINAL MODELS")
    saved_models = save_final_models(X, y, models, outdir)

    section("WRITING OUTPUT FILES")

    cleaning_report = meta.get("training_feature_cleaning", {})

    metadata = {
        "script": Path(__file__).name,
        "dataset_id": args.dataset_id,
        "dataset_name": dataset_name,
        "dataset_dir": str(dataset_dir),
        "output_dir": str(outdir),
        "target_col": target_col,
        "target_description": meta.get("target_description", ""),
        "group_name": meta.get("group_name", ""),
        "tier": meta.get("tier", ""),
        "feature_numbers": meta.get("feature_numbers_used", meta.get("requested_feature_numbers", [])),
        "n_genes": n_genes,
        "n_positives": n_pos,
        "n_negatives": n_neg,
        "prevalence": prevalence,
        "n_features": len(feature_cols),
        "cv_folds": args.cv,
        "n_permutations": 0 if args.no_permutation else N_PERMUTATIONS,
        "n_bootstrap": N_BOOTSTRAP,
        "models": list(models.keys()),
        "best_model": best_model_name,
        "best_auroc_mean": safe_float(best_row.get("auroc_mean")),
        "best_auprc_mean": safe_float(best_row.get("auprc_mean")),
        "best_mcc_mean": safe_float(best_row.get("mcc_mean")),
        "random_seed": RANDOM_SEED,
        "xgboost_available": HAS_XGBOOST,
        "run_datetime": now_iso(),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "saved_models": saved_models,
        "training_feature_cleaning": cleaning_report,
    }

    with open(outdir / "00_training_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)
    log(f"  [SAVED] {outdir / '00_training_metadata.json'}")

    pd.DataFrame({"feature": feature_cols}).to_csv(
        outdir / "01_feature_columns_used.csv",
        index=False,
    )
    log(f"  [SAVED] {outdir / '01_feature_columns_used.csv'}")

    pd.DataFrame([{
        "target_col": target_col,
        "target_description": meta.get("target_description", ""),
        "n_positives": n_pos,
        "n_negatives": n_neg,
        "prevalence": prevalence,
    }]).to_csv(outdir / "02_target_info.csv", index=False)
    log(f"  [SAVED] {outdir / '02_target_info.csv'}")

    pd.DataFrame([cleaning_report]).to_csv(
        outdir / "00_feature_cleaning_report.csv",
        index=False,
    )
    log(f"  [SAVED] {outdir / '00_feature_cleaning_report.csv'}")

    save_csv(preds_df, outdir / "03_cv_all_predictions.csv")
    save_csv(folds_df, outdir / "04_cv_fold_metrics.csv")
    save_csv(comparison, outdir / "05_model_comparison_summary.csv")

    pd.DataFrame([best_row]).to_csv(
        outdir / "06_best_model_summary.csv",
        index=False,
    )
    log(f"  [SAVED] {outdir / '06_best_model_summary.csv'}")

    save_csv(imp_df, outdir / "07_feature_importance_all_folds.csv")
    save_csv(imp_summary, outdir / "08_feature_importance_summary.csv")
    save_csv(top_feat, outdir / "09_top_features_best_model.csv")
    save_csv(cm_df, outdir / "10_confusion_matrix_per_fold.csv")
    save_csv(agg_cm, outdir / "11_confusion_matrix_aggregate.csv")

    if perm_results:
        perm_rows = [{"model": k, **v} for k, v in perm_results.items()]
        save_csv(pd.DataFrame(perm_rows), outdir / "12_permutation_test_results.csv")
    else:
        empty_perm = pd.DataFrame(columns=[
            "model",
            "permutation_auroc_mean",
            "permutation_auroc_std",
            "permutation_p_value",
            "n_permutations_run",
            "observed_auroc",
            "auroc_above_random",
            "statistically_significant",
        ])
        empty_perm.to_csv(outdir / "12_permutation_test_results.csv", index=False)
        log(f"  [SAVED empty] {outdir / '12_permutation_test_results.csv'}")

    save_csv(audit_df, outdir / "13_leakage_audit.csv")

    write_publication_summary(
        outdir=outdir,
        dataset_name=dataset_name,
        target_col=target_col,
        meta=meta,
        comparison=comparison,
        best_row=best_row,
        top_feat=top_feat,
    )

    section("DONE")
    log(f"Dataset:      {dataset_name}")
    log(f"Best model:   {best_model_name}")
    log(f"AUROC:        {safe_float(best_row.get('auroc_mean')):.4f}")
    log(f"AUPRC:        {safe_float(best_row.get('auprc_mean')):.4f}")
    log(f"F1 binary:    {safe_float(best_row.get('f1_binary_mean')):.4f}")
    log(f"MCC:          {safe_float(best_row.get('mcc_mean')):.4f}")
    log(f"Output dir:   {outdir}")
    log("=" * 100)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        log("\n[STOPPED] Interrupted by user.")
        sys.exit(130)
    except Exception as exc:
        log(f"\n[ERROR] {exc}")
        raise