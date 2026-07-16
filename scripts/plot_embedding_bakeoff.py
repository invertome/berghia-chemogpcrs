#!/usr/bin/env python3
"""Publication-quality model-ranking figure for the PLM embedding bake-off.

Reads the tidy bakeoff_metrics.tsv (from embedding_bakeoff_metrics.py) and draws
a paired dumbbell plot of LOO-famAcc per model, FULL-1094 vs VERIFIED-723.

Why LOO-famAcc and not LOFO-AUROC: on the 830 characterized class-A anchors the
novelty-detection AUROC saturates at ~1.0 for every competent model, so it does
not rank them; within-class family assignment accuracy does. The FULL-vs-VERIFIED
pairing simultaneously shows the reviewed-inclusion effect (bead cw3.9): a short
bar = the 371 reviewed anchors barely change discrimination.

Colours are Okabe-Ito (colourblind-safe); redundant marker shapes back up colour.

Usage:
  plot_embedding_bakeoff.py results/ranking/diagnostics/bakeoff_metrics.tsv \\
      --out figures/embedding/model_ranking
"""
from __future__ import annotations

import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# production scorer config: tied-cov Mahalanobis, raw calibration, multi-prototype
PROD_CONFIG = ("maha", "raw", "multi")
# Okabe-Ito
_FULL_C, _VER_C, _LINE_C = "#0072B2", "#D55E00", "#999999"


def load_metrics(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def ranking_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Rows for the production config only, models ordered by the ranking metric
    (VERIFIED famAcc where present, else FULL famAcc), descending."""
    s, c, ce = PROD_CONFIG
    rf = df[(df["score"] == s) & (df["calib"] == c) & (df["cent"] == ce)].copy()
    if rf.empty:
        return rf
    piv = rf.pivot_table(index="model", columns="refset", values="loo_famacc")
    rank_key = piv["VERIFIED"] if "VERIFIED" in piv else piv.get("FULL")
    if "VERIFIED" in piv and "FULL" in piv:
        rank_key = piv["VERIFIED"].fillna(piv["FULL"])
    order = rank_key.sort_values(ascending=False).index.tolist()
    rf["model"] = pd.Categorical(rf["model"], categories=order, ordered=True)
    return rf.sort_values("model")


def ranking_figure(df: pd.DataFrame) -> plt.Figure:
    rf = ranking_frame(df)
    if rf.empty:
        raise ValueError("no production-config rows to plot")
    piv = rf.pivot_table(index="model", columns="refset", values="loo_famacc", observed=True)
    models = list(piv.index)                       # already famAcc-ordered
    y = range(len(models))

    fig, ax = plt.subplots(figsize=(5.2, 0.32 * len(models) + 1.0))
    for yi, m in zip(y, models):
        full = piv.loc[m].get("FULL"); ver = piv.loc[m].get("VERIFIED")
        if pd.notna(full) and pd.notna(ver):       # dumbbell connector
            ax.plot([full, ver], [yi, yi], color=_LINE_C, lw=1.2, zorder=1)
        if pd.notna(full):
            ax.scatter(full, yi, s=34, color=_FULL_C, marker="o", zorder=3)
        if pd.notna(ver):
            ax.scatter(ver, yi, s=34, color=_VER_C, marker="D", zorder=3)

    # explicit legend proxies -- the top-ranked model may lack a VERIFIED point
    # (its harness timed out), so a positional label would silently drop it.
    n_full = int(piv.get("FULL", pd.Series(dtype=float)).notna().sum())
    n_ver = int(piv.get("VERIFIED", pd.Series(dtype=float)).notna().sum())
    handles = [
        plt.Line2D([], [], marker="o", ls="", color=_FULL_C, label=f"FULL (1094), n={n_full}"),
        plt.Line2D([], [], marker="D", ls="", color=_VER_C, label=f"VERIFIED (723), n={n_ver}"),
    ]
    ax.set_yticks(list(y)); ax.set_yticklabels(models)
    ax.invert_yaxis()                              # best model at top
    ax.set_xlabel("LOO family-assignment accuracy")
    ax.set_title("PLM embedding bake-off — family discrimination", fontsize=10)
    ax.axvline(0.37, color="0.7", ls=":", lw=0.9)  # majority baseline
    ax.text(0.37, len(models) - 0.3, " majority\n baseline", fontsize=6,
            color="0.5", va="bottom", ha="left")
    ax.set_xlim(0.3, min(1.0, float(piv.max().max()) + 0.06))
    ax.legend(handles=handles, frameon=False, fontsize=7, loc="lower right")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("metrics_tsv")
    p.add_argument("--out", required=True, help="output path stem (no extension)")
    a = p.parse_args(argv)

    df = load_metrics(a.metrics_tsv)
    try:
        fig = ranking_figure(df)
    except ValueError as e:
        sys.exit(f"[plot_embedding_bakeoff] {e}: {a.metrics_tsv}")
    import os
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[plot_embedding_bakeoff] wrote {a.out}.pdf / .png")


if __name__ == "__main__":
    main()
