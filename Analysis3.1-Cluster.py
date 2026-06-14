#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Analysis3.1-Cluster.py

PCA and clustering visualisation of generated gene-level feature datasets.

This script reads the actual generated dataset folders:

    Datasets/Dataset177_CUM01_T1/
        X_train.csv
        y_train.csv
        meta.json
        feature_columns.txt
        missingness_summary.csv

It does not require GRP_* folders. If thematic GRP folders are absent, it
automatically selects available cumulative datasets for the requested target.

Major fixes
-----------
1. Reads existing DatasetXXX_<GROUP>_<TARGET> folders.
2. Handles missing GRP_* datasets.
3. Replaces inf/-inf with NaN before PCA.
4. Drops columns with extreme numeric values.
5. Clips remaining values before scaling.
6. Avoids matplotlib tight_layout recursion bug.
7. Saves PCA coordinates, panel summaries, KMeans summaries, PDF and PNG.

Run
---
python Analysis3.1-Cluster.py

Recommended
-----------
python Analysis3.1-Cluster.py --datasets-dir Datasets --target T1 --auto

Force safe groups
-----------------
python Analysis3.1-Cluster.py --datasets-dir Datasets --target T1 --groups CUM01 CUM04 CUM07 CUM11 CUM14

Outputs
-------
Analysis3.1_Cluster_Output/
    Analysis3.1_available_datasets_T1.csv
    Analysis3.1_PCA_T1.pdf
    Analysis3.1_PCA_T1.png
    Analysis3.1_PCA_coordinates_T1.csv
    Analysis3.1_PCA_panel_summary_T1.csv
    Analysis3.1_KMeans_summary_T1.csv
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import warnings
from pathlib import Path
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


# =============================================================================
# SETTINGS
# =============================================================================

DEFAULT_DATASETS_DIR = Path("Datasets")
DEFAULT_OUTDIR = Path("Analysis3.1_Cluster_Output")
DEFAULT_TARGET = "T1"

TARGET_NAMES = {
    "T1": "Clinical target",
    "T2": "Clinical investigation target",
    "T3": "Small-molecule target",
    "T4": "Chemical tractability target",
    "T5": "Biologic/modality target",
    "T6": "Drug-gene interaction target",
    "T7": "Potentially druggable family target",
    "T8": "Broad druggability target",
}

PREFERRED_THEMATIC_GROUPS = [
    "GRP_Full",
    "GRP_Functional",
    "GRP_Network",
    "GRP_Expression",
    "GRP_Constraint",
    "GRP_Annotation",
]

PREFERRED_CUM_GROUPS = [
    "CUM01", "CUM02", "CUM03", "CUM04", "CUM05", "CUM06",
    "CUM07", "CUM08", "CUM09", "CUM10", "CUM11", "CUM12",
    "CUM13", "CUM14", "CUM15", "CUM16", "CUM17", "CUM18",
    "CUM19", "CUM20", "CUM21", "CUM22",
]

PREFERRED_SINGLE_FEATURE_GROUPS = [
    "F01", "F02", "F03", "F04", "F05", "F06", "F07", "F08",
    "F09", "F10", "F11", "F12", "F13", "F14", "F15", "F16",
    "F17", "F18", "F19", "F20", "F21", "F22",
]


# =============================================================================
# BASIC HELPERS
# =============================================================================

def log(msg: str) -> None:
    print(msg, flush=True)


def now_iso() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except Exception:
        return {}


def extract_dataset_id(dataset_name: str) -> int:
    m = re.match(r"^Dataset(\d+)_", dataset_name)
    if not m:
        return 10**9
    return int(m.group(1))


def extract_group_target_from_folder(name: str) -> Tuple[str, str]:
    """
    Dataset177_CUM01_T1      -> CUM01, T1
    Dataset001_F01_T1        -> F01, T1
    Dataset353_GRP_Full_T1   -> GRP_Full, T1
    """
    m = re.match(r"^Dataset\d+_(.+)_(T\d+)$", name)
    if not m:
        return "", ""
    return m.group(1), m.group(2)


def classify_group(group: str) -> str:
    if group.startswith("GRP_"):
        return "Thematic"
    if group.startswith("CUM"):
        return "Cumulative"
    if re.match(r"^F\d+$", group):
        return "Individual"
    return "Unknown"


def group_sort_key(group: str) -> Tuple[int, int, str]:
    if group.startswith("GRP_"):
        return (0, 0, group)

    if group.startswith("CUM"):
        m = re.match(r"^CUM(\d+)$", group)
        return (1, int(m.group(1)) if m else 999, group)

    if re.match(r"^F\d+$", group):
        m = re.match(r"^F(\d+)$", group)
        return (2, int(m.group(1)) if m else 999, group)

    return (9, 999, group)


# =============================================================================
# DATA CLEANING
# =============================================================================

def clean_numeric_matrix(
    X: pd.DataFrame,
    extreme_threshold: float = 1e100,
    verbose: bool = True,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Convert to numeric, remove inf values, drop unusable columns, and clip values.

    This is deliberately conservative because PCA and StandardScaler cannot
    handle inf/-inf or absurd numeric values.
    """
    report: Dict[str, Any] = {
        "n_rows_input": int(X.shape[0]),
        "n_cols_input": int(X.shape[1]),
        "n_inf_values_replaced": 0,
        "n_cols_all_nan_dropped": 0,
        "n_cols_constant_dropped": 0,
        "n_cols_extreme_dropped": 0,
        "n_cols_output": 0,
    }

    X2 = X.copy()

    # Convert every column to numeric.
    for col in X2.columns:
        X2[col] = pd.to_numeric(X2[col], errors="coerce")

    # Replace inf/-inf.
    inf_mask = np.isinf(X2.to_numpy(dtype=np.float64, copy=True))
    n_inf = int(inf_mask.sum())
    report["n_inf_values_replaced"] = n_inf

    if n_inf > 0:
        X2 = X2.replace([np.inf, -np.inf], np.nan)

    # Drop columns with absurd values.
    extreme_cols = []
    for col in X2.columns:
        vals = X2[col].dropna()
        if vals.empty:
            continue

        max_abs = vals.abs().max()

        if pd.notna(max_abs) and max_abs > extreme_threshold:
            extreme_cols.append(col)

    if extreme_cols:
        X2 = X2.drop(columns=extreme_cols)

    report["n_cols_extreme_dropped"] = len(extreme_cols)

    # Drop all-NaN columns.
    all_nan_cols = X2.columns[X2.isna().all()].tolist()

    if all_nan_cols:
        X2 = X2.drop(columns=all_nan_cols)

    report["n_cols_all_nan_dropped"] = len(all_nan_cols)

    if X2.shape[1] == 0:
        report["n_cols_output"] = 0
        return X2, report

    # Drop constant columns.
    nunique = X2.nunique(dropna=True)
    constant_cols = nunique[nunique <= 1].index.tolist()

    if constant_cols:
        X2 = X2.drop(columns=constant_cols)

    report["n_cols_constant_dropped"] = len(constant_cols)

    if X2.shape[1] == 0:
        report["n_cols_output"] = 0
        return X2, report

    # Final safety clipping.
    X2 = X2.clip(lower=-extreme_threshold, upper=extreme_threshold)

    # Replace any remaining inf just in case.
    X2 = X2.replace([np.inf, -np.inf], np.nan)

    report["n_cols_output"] = int(X2.shape[1])

    if verbose:
        if (
            report["n_inf_values_replaced"] > 0
            or report["n_cols_extreme_dropped"] > 0
            or report["n_cols_all_nan_dropped"] > 0
            or report["n_cols_constant_dropped"] > 0
        ):
            log(
                "  [CLEAN] "
                f"inf_values={report['n_inf_values_replaced']:,}; "
                f"extreme_cols={report['n_cols_extreme_dropped']:,}; "
                f"all_nan_cols={report['n_cols_all_nan_dropped']:,}; "
                f"constant_cols={report['n_cols_constant_dropped']:,}; "
                f"features_kept={report['n_cols_output']:,}"
            )

    return X2, report


# =============================================================================
# DISCOVERY
# =============================================================================

def discover_datasets(datasets_dir: Path) -> pd.DataFrame:
    """
    Discover dataset folders by scanning Datasets/.

    This is more robust than depending on dataset_catalogue.csv because your
    generated folder set may be partial.
    """
    if not datasets_dir.exists():
        raise FileNotFoundError(f"Datasets directory not found: {datasets_dir}")

    rows = []

    for p in sorted(datasets_dir.iterdir(), key=lambda x: x.name):
        if not p.is_dir():
            continue

        group, target = extract_group_target_from_folder(p.name)

        if not group or not target:
            continue

        xfile = p / "X_train.csv"
        yfile = p / "y_train.csv"
        meta_file = p / "meta.json"
        meta = read_json(meta_file)

        rows.append({
            "dataset_id": meta.get("dataset_id", extract_dataset_id(p.name)),
            "dataset_name": p.name,
            "dataset_dir": str(p),
            "group_name": meta.get("group_name", group),
            "target_key": meta.get("target_key", target),
            "target_name": meta.get("target_name", TARGET_NAMES.get(target, target)),
            "tier": meta.get("tier", classify_group(group)),
            "has_X_train": int(xfile.exists()),
            "has_y_train": int(yfile.exists()),
            "has_meta_json": int(meta_file.exists()),
            "n_genes_meta": meta.get("n_genes", np.nan),
            "n_positives_meta": meta.get("n_positives", np.nan),
            "n_features_used_meta": meta.get("n_features_used", np.nan),
        })

    out = pd.DataFrame(rows)

    if out.empty:
        return out

    out = out.sort_values(["target_key", "dataset_id", "dataset_name"]).reset_index(drop=True)

    return out


def choose_groups_for_target(
    available: pd.DataFrame,
    target: str,
    requested_groups: Optional[List[str]],
    max_panels: int,
) -> List[str]:
    """
    Select groups for PCA.

    Priority:
    1. User-specified groups.
    2. Thematic groups if present.
    3. Cumulative groups if present.
    4. Individual groups if present.
    """
    sub = available[
        (available["target_key"].astype(str) == str(target))
        & (available["has_X_train"] == 1)
        & (available["has_y_train"] == 1)
    ].copy()

    groups_present = set(sub["group_name"].astype(str).tolist())

    if requested_groups:
        return requested_groups

    thematic = [g for g in PREFERRED_THEMATIC_GROUPS if g in groups_present]

    if thematic:
        return thematic[:max_panels]

    cumulative = [g for g in PREFERRED_CUM_GROUPS if g in groups_present]

    if cumulative:
        if len(cumulative) <= max_panels:
            return cumulative

        idx = np.linspace(0, len(cumulative) - 1, max_panels).round().astype(int)
        return [cumulative[i] for i in idx]

    individual = [g for g in PREFERRED_SINGLE_FEATURE_GROUPS if g in groups_present]

    if individual:
        if len(individual) <= max_panels:
            return individual

        idx = np.linspace(0, len(individual) - 1, max_panels).round().astype(int)
        return [individual[i] for i in idx]

    return sorted(groups_present, key=group_sort_key)[:max_panels]


def find_dataset_dir(
    available: pd.DataFrame,
    group: str,
    target: str,
) -> Optional[Path]:
    sub = available[
        (available["group_name"].astype(str) == str(group))
        & (available["target_key"].astype(str) == str(target))
        & (available["has_X_train"] == 1)
        & (available["has_y_train"] == 1)
    ].copy()

    if sub.empty:
        return None

    sub = sub.sort_values(["dataset_id", "dataset_name"])

    return Path(str(sub.iloc[0]["dataset_dir"]))


# =============================================================================
# LOAD DATA
# =============================================================================

def load_dataset(dataset_dir: Path) -> Tuple[pd.DataFrame, pd.Series, str]:
    xfile = dataset_dir / "X_train.csv"
    yfile = dataset_dir / "y_train.csv"

    if not xfile.exists():
        raise FileNotFoundError(f"Missing X_train.csv: {xfile}")

    if not yfile.exists():
        raise FileNotFoundError(f"Missing y_train.csv: {yfile}")

    X = pd.read_csv(xfile, low_memory=False)
    y = pd.read_csv(yfile, low_memory=False)

    if "gene_symbol" not in X.columns:
        raise RuntimeError(f"{xfile} must contain gene_symbol")

    if "gene_symbol" not in y.columns:
        raise RuntimeError(f"{yfile} must contain gene_symbol")

    target_cols = [c for c in y.columns if c != "gene_symbol"]

    if not target_cols:
        raise RuntimeError(f"{yfile} has no target column")

    target_col = target_cols[0]

    df = X.merge(
        y[["gene_symbol", target_col]],
        on="gene_symbol",
        how="inner",
    )

    genes = df["gene_symbol"].astype(str)
    yy = pd.to_numeric(df[target_col], errors="coerce").fillna(0).astype(int)

    Xmat = df.drop(columns=["gene_symbol", target_col], errors="ignore")
    Xmat.index = genes
    yy.index = genes

    Xmat, clean_report = clean_numeric_matrix(Xmat)

    if Xmat.shape[1] == 0:
        raise RuntimeError("No usable numeric features after cleaning")

    return Xmat, yy, target_col


def downsample_for_plot(
    X: pd.DataFrame,
    y: pd.Series,
    max_genes: int,
    random_state: int,
) -> Tuple[pd.DataFrame, pd.Series]:
    """
    Stratified downsample for speed and readability.
    """
    if max_genes <= 0 or len(X) <= max_genes:
        return X, y

    rng = np.random.default_rng(random_state)

    pos_idx = np.where(y.values == 1)[0]
    neg_idx = np.where(y.values == 0)[0]

    n_pos = len(pos_idx)
    n_neg = len(neg_idx)

    if n_pos == 0 or n_neg == 0:
        chosen = rng.choice(np.arange(len(X)), size=max_genes, replace=False)
    else:
        pos_frac = n_pos / len(X)

        keep_pos = max(1, int(round(max_genes * pos_frac)))
        keep_neg = max_genes - keep_pos

        keep_pos = min(keep_pos, n_pos)
        keep_neg = min(keep_neg, n_neg)

        chosen_pos = rng.choice(pos_idx, size=keep_pos, replace=False)
        chosen_neg = rng.choice(neg_idx, size=keep_neg, replace=False)

        chosen = np.concatenate([chosen_pos, chosen_neg])
        rng.shuffle(chosen)

    return X.iloc[chosen].copy(), y.iloc[chosen].copy()


# =============================================================================
# PCA / CLUSTERING
# =============================================================================

def run_pca(X: pd.DataFrame) -> Tuple[np.ndarray, PCA]:
    """
    Median-impute, scale and run PCA.

    This function performs another final safety clean before sklearn.
    """
    X2 = X.copy()
    X2 = X2.replace([np.inf, -np.inf], np.nan)
    X2 = X2.clip(lower=-1e100, upper=1e100)

    # Guarantee finite or NaN only.
    arr = X2.to_numpy(dtype=np.float64, copy=True)
    arr[~np.isfinite(arr)] = np.nan
    X2 = pd.DataFrame(arr, index=X2.index, columns=X2.columns)

    pipe = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("pca", PCA(n_components=2, random_state=42)),
        ]
    )

    coords = pipe.fit_transform(X2)
    pca = pipe.named_steps["pca"]

    return coords, pca


def run_kmeans(
    coords: np.ndarray,
    n_clusters: int,
    random_state: int,
) -> Tuple[np.ndarray, float]:
    if n_clusters <= 1 or coords.shape[0] <= n_clusters:
        return np.zeros(coords.shape[0], dtype=int), np.nan

    km = KMeans(n_clusters=n_clusters, n_init=20, random_state=random_state)
    labels = km.fit_predict(coords)

    try:
        sil = float(silhouette_score(coords, labels))
    except Exception:
        sil = np.nan

    return labels, sil


# =============================================================================
# PLOTTING
# =============================================================================

def make_safe_pca_plot(
    coords_df: pd.DataFrame,
    panel_results: List[Dict[str, Any]],
    target: str,
    target_name: str,
    outdir: Path,
    dpi: int,
) -> Tuple[Path, Path]:
    """
    Make PCA figure safely.

    Important:
    - Do NOT use fig.tight_layout()
    - Do NOT use matplotlib PDF backend
    - Save PNG with Agg backend
    - Convert PNG to PDF using PIL

    This avoids the matplotlib recursion crash on Python 3.14.
    """
    n_panels = len(panel_results)

    if n_panels == 0:
        raise RuntimeError("No panels to plot")

    ncols = min(3, n_panels)
    nrows = int(np.ceil(n_panels / ncols))

    fig_w = 5.4 * ncols
    fig_h = 4.8 * nrows

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(fig_w, fig_h),
        squeeze=False,
    )

    axes_flat = axes.ravel()

    for ax in axes_flat:
        ax.set_axis_off()

    for i, row in enumerate(panel_results):
        ax = axes_flat[i]
        group = row["group"]

        sub = coords_df[coords_df["group"] == group].copy()

        neg = sub[sub["target_label"] == 0]
        pos = sub[sub["target_label"] == 1]

        ax.set_axis_on()

        # Critical: remove tick machinery to avoid matplotlib recursion bug.
        ax.set_xticks([])
        ax.set_yticks([])
        ax.tick_params(
            left=False,
            bottom=False,
            labelleft=False,
            labelbottom=False,
        )

        ax.scatter(
            neg["PC1"],
            neg["PC2"],
            s=8,
            alpha=0.30,
            linewidths=0,
            marker="o",
            label="Negative",
        )

        ax.scatter(
            pos["PC1"],
            pos["PC2"],
            s=16,
            alpha=0.80,
            linewidths=0,
            marker="o",
            label="Positive",
        )

        ax.axhline(0, linewidth=0.5, alpha=0.30)
        ax.axvline(0, linewidth=0.5, alpha=0.30)

        # Avoid ax.set_title because it triggers title-position calculations.
        ax.text(
            0.5,
            1.04,
            (
                f"{group}\n"
                f"{int(row['n_features']):,} features | "
                f"PC1 {float(row['pc1_variance_pct']):.1f}% | "
                f"PC2 {float(row['pc2_variance_pct']):.1f}%"
            ),
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )

        ax.text(
            0.5,
            -0.06,
            "PC1",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=9,
        )

        ax.text(
            -0.06,
            0.5,
            "PC2",
            transform=ax.transAxes,
            ha="right",
            va="center",
            rotation=90,
            fontsize=9,
        )

        if i == 0:
            ax.legend(frameon=False, fontsize=9, loc="best")

    # Avoid fig.suptitle because it can trigger layout calculations.
    fig.text(
        0.5,
        0.985,
        f"PCA visualisation of gene-level feature spaces — {target} ({target_name})",
        ha="center",
        va="top",
        fontsize=14,
        fontweight="bold",
    )

    # Manual spacing only. Never tight_layout / constrained_layout.
    fig.subplots_adjust(
        left=0.06,
        right=0.98,
        bottom=0.08,
        top=0.88,
        wspace=0.28,
        hspace=0.55,
    )

    png_path = outdir / f"Analysis3.1_PCA_{target}.png"
    pdf_path = outdir / f"Analysis3.1_PCA_{target}.pdf"

    # Save PNG using Agg backend.
    fig.savefig(png_path, dpi=dpi)
    plt.close(fig)

    # Convert PNG to PDF using PIL, avoiding matplotlib PDF backend completely.
    try:
        from PIL import Image

        img = Image.open(png_path).convert("RGB")
        img.save(pdf_path, "PDF", resolution=float(dpi))

    except Exception as exc:
        log(f"[WARNING] Could not convert PNG to PDF using PIL: {exc}")
        log("[WARNING] PNG was still saved successfully.")
        pdf_path = Path("")

    return pdf_path, png_path


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def make_figure(
    datasets_dir: Path,
    outdir: Path,
    target: str,
    requested_groups: Optional[List[str]],
    max_panels: int,
    max_genes: int,
    n_clusters: int,
    random_state: int,
    dpi: int,
) -> None:

    mkdir(outdir)

    target_name = TARGET_NAMES.get(target, target)

    log("=" * 100)
    log("PCA VISUALISATION OF GENE-LEVEL FEATURE SPACES")
    log(f"Started:  {now_iso()}")
    log(f"Datasets: {datasets_dir.resolve()}")
    log(f"Target:   {target} ({target_name})")
    log("=" * 100)

    available = discover_datasets(datasets_dir)

    if available.empty:
        raise RuntimeError(f"No DatasetXXX_* folders found in {datasets_dir}")

    available_target = available[available["target_key"].astype(str) == str(target)].copy()

    available_path = outdir / f"Analysis3.1_available_datasets_{target}.csv"
    available_target.to_csv(available_path, index=False)

    log(f"[AVAILABLE DATASETS TOTAL]  {len(available):,}")
    log(f"[AVAILABLE FOR {target}]          {len(available_target):,}")
    log(f"[SAVED AVAILABLE TABLE]     {available_path}")

    if available_target.empty:
        log("")
        log("[ERROR] No datasets found for this target.")
        log("Available target counts:")
        print(available["target_key"].value_counts().sort_index().to_string())
        raise RuntimeError(f"No datasets found for target {target}")

    log("")
    log("[AVAILABLE GROUPS FOR TARGET]")

    groups_by_tier = (
        available_target.groupby(["tier", "group_name"])
        .size()
        .reset_index(name="n")
        .sort_values(["tier", "group_name"])
    )

    print(groups_by_tier.to_string(index=False))

    groups = choose_groups_for_target(
        available=available,
        target=target,
        requested_groups=requested_groups,
        max_panels=max_panels,
    )

    log("")
    log("[GROUPS SELECTED FOR PCA]")

    for g in groups:
        log(f"  - {g}")

    panel_results: List[Dict[str, Any]] = []
    coord_tables: List[pd.DataFrame] = []
    kmeans_rows: List[Dict[str, Any]] = []

    for group in groups:
        ddir = find_dataset_dir(
            available=available,
            group=group,
            target=target,
        )

        if ddir is None:
            log(f"[MISSING] {group}_{target} — not present in generated folders")
            continue

        try:
            X, y, target_col = load_dataset(ddir)
        except Exception as exc:
            log(f"[SKIP] {ddir.name} — load failed: {exc}")
            continue

        if X.shape[0] < 3:
            log(f"[SKIP] {ddir.name} — too few genes: {X.shape[0]}")
            continue

        if X.shape[1] < 2:
            log(f"[SKIP] {ddir.name} — too few usable numeric features: {X.shape[1]}")
            continue

        X_plot, y_plot = downsample_for_plot(
            X=X,
            y=y,
            max_genes=max_genes,
            random_state=random_state,
        )

        try:
            coords, pca = run_pca(X_plot)
        except Exception as exc:
            log(f"[SKIP] {ddir.name} — PCA failed after cleaning: {exc}")
            continue

        cluster_labels, sil = run_kmeans(
            coords=coords,
            n_clusters=n_clusters,
            random_state=random_state,
        )

        evr = pca.explained_variance_ratio_ * 100.0

        panel_result = {
            "group": group,
            "tier": classify_group(group),
            "dataset_name": ddir.name,
            "dataset_dir": str(ddir),
            "target": target,
            "target_column": target_col,
            "n_genes_total": int(X.shape[0]),
            "n_genes_plotted": int(X_plot.shape[0]),
            "n_features": int(X.shape[1]),
            "n_positive_plotted": int((y_plot == 1).sum()),
            "n_negative_plotted": int((y_plot == 0).sum()),
            "pc1_variance_pct": round(float(evr[0]), 6),
            "pc2_variance_pct": round(float(evr[1]), 6),
            "kmeans_n_clusters": int(n_clusters),
            "kmeans_silhouette_on_pc12": round(float(sil), 6) if pd.notna(sil) else np.nan,
        }

        panel_results.append(panel_result)

        coord_df = pd.DataFrame(
            {
                "gene_symbol": X_plot.index.astype(str),
                "group": group,
                "tier": classify_group(group),
                "dataset_name": ddir.name,
                "target": target,
                "target_label": y_plot.values.astype(int),
                "PC1": coords[:, 0],
                "PC2": coords[:, 1],
                "kmeans_cluster": cluster_labels.astype(int),
            }
        )

        coord_tables.append(coord_df)

        for cl in sorted(np.unique(cluster_labels)):
            idx = cluster_labels == cl
            yy = y_plot.values[idx]

            kmeans_rows.append(
                {
                    "group": group,
                    "dataset_name": ddir.name,
                    "target": target,
                    "cluster": int(cl),
                    "n_genes": int(idx.sum()),
                    "n_positive": int((yy == 1).sum()),
                    "n_negative": int((yy == 0).sum()),
                    "positive_fraction": round(float((yy == 1).mean()), 6) if idx.sum() else np.nan,
                }
            )

        sil_txt = f"{sil:.3f}" if pd.notna(sil) else "NA"

        log(
            f"[OK] {ddir.name:<35s} "
            f"genes={X.shape[0]:,} plotted={X_plot.shape[0]:,} "
            f"features={X.shape[1]:,} "
            f"PC1={evr[0]:.2f}% PC2={evr[1]:.2f}% "
            f"silhouette={sil_txt}"
        )

    if not panel_results:
        log("")
        log("[ERROR] No PCA panels could be generated.")
        log("")
        log("Try forcing known-good groups:")
        log(
            "  python Analysis3.1-Cluster.py "
            "--datasets-dir Datasets "
            "--target T1 "
            "--groups CUM01 CUM04 CUM07 CUM11 CUM14"
        )
        raise RuntimeError("No PCA panels could be generated.")

    summary_df = pd.DataFrame(panel_results)
    coords_df = pd.concat(coord_tables, ignore_index=True)
    kmeans_df = pd.DataFrame(kmeans_rows)

    summary_path = outdir / f"Analysis3.1_PCA_panel_summary_{target}.csv"
    coords_path = outdir / f"Analysis3.1_PCA_coordinates_{target}.csv"
    kmeans_path = outdir / f"Analysis3.1_KMeans_summary_{target}.csv"

    summary_df.to_csv(summary_path, index=False)
    coords_df.to_csv(coords_path, index=False)
    kmeans_df.to_csv(kmeans_path, index=False)

    pdf_path, png_path = make_safe_pca_plot(
        coords_df=coords_df,
        panel_results=panel_results,
        target=target,
        target_name=target_name,
        outdir=outdir,
        dpi=dpi,
    )

    log("")
    log("=" * 100)
    log("[DONE]")
    log(f"[SAVED FIGURE PDF]     {pdf_path.resolve()}")
    log(f"[SAVED FIGURE PNG]     {png_path.resolve()}")
    log(f"[SAVED PCA SUMMARY]    {summary_path.resolve()}")
    log(f"[SAVED PCA COORDS]     {coords_path.resolve()}")
    log(f"[SAVED KMEANS SUMMARY] {kmeans_path.resolve()}")
    log(f"Finished: {now_iso()}")
    log("=" * 100)


# =============================================================================
# CLI
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="PCA / clustering visualisation of generated gene-level datasets."
    )

    parser.add_argument(
        "--datasets-dir",
        default=str(DEFAULT_DATASETS_DIR),
        help=f"Directory containing DatasetXXX_* folders. Default: {DEFAULT_DATASETS_DIR}",
    )

    parser.add_argument(
        "--outdir",
        default=str(DEFAULT_OUTDIR),
        help=f"Output directory. Default: {DEFAULT_OUTDIR}",
    )

    parser.add_argument(
        "--target",
        default=DEFAULT_TARGET,
        help="Target key: T1, T2, ..., T8. Default: T1",
    )

    parser.add_argument(
        "--groups",
        nargs="*",
        default=None,
        help="Optional exact groups to plot, e.g. --groups CUM01 CUM05 CUM10 CUM15",
    )

    parser.add_argument(
        "--auto",
        action="store_true",
        help="Accepted for compatibility. Auto-selection is used whenever --groups is not provided.",
    )

    parser.add_argument(
        "--max-panels",
        type=int,
        default=6,
        help="Maximum number of panels to plot when selecting automatically. Default: 6",
    )

    parser.add_argument(
        "--max-genes",
        type=int,
        default=8000,
        help="Maximum genes per panel. Use 0 for all genes. Default: 8000",
    )

    parser.add_argument(
        "--n-clusters",
        type=int,
        default=4,
        help="KMeans clusters on PC1/PC2 for summary CSV. Default: 4",
    )

    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random seed. Default: 42",
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="PNG DPI. Default: 300",
    )

    args = parser.parse_args()

    warnings.filterwarnings("ignore", category=RuntimeWarning)

    make_figure(
        datasets_dir=Path(args.datasets_dir),
        outdir=Path(args.outdir),
        target=args.target,
        requested_groups=args.groups,
        max_panels=args.max_panels,
        max_genes=args.max_genes,
        n_clusters=args.n_clusters,
        random_state=args.random_state,
        dpi=args.dpi,
    )


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        log("\n[STOPPED] Interrupted by user.")
        sys.exit(130)
    except Exception as exc:
        log(f"\n[ERROR] {exc}")
        raise