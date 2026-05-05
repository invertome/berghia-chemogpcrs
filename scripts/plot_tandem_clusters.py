#!/usr/bin/env python3
"""plot_tandem_clusters.py — Visualize tandem-cluster detection output.

Consumes the CSV produced by ``compute_tandem_clusters.py`` and writes a
two-panel figure: (A) histogram of tandem-cluster sizes (singletons through
largest), and (B) bar chart of the top-N largest distinct clusters.

CSV schema (from compute_tandem_clusters.py):
    candidate_id, tandem_cluster_size, tandem_cluster_id

Usage:
    python plot_tandem_clusters.py <tandem_clusters.csv> <output_prefix>

Outputs ``<output_prefix>.{png,svg,pdf}``.

Design decisions:
    - Header-aware (uses pandas), tolerates blank tandem_cluster_id.
    - Empty / all-singleton input still produces a figure (with a clear
      "no clusters detected" annotation) so stage 09 can always link to it.
    - No external state; deterministic given the CSV.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def _setup_axes(ax: plt.Axes, title: str, xlabel: str, ylabel: str) -> None:
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot(df: pd.DataFrame, output_prefix: str, top_n: int = 15) -> None:
    """Render the two-panel figure.

    Parameters
    ----------
    df : pd.DataFrame
        Must have columns ``candidate_id``, ``tandem_cluster_size``,
        ``tandem_cluster_id``. Sizes may be NaN for candidates missing
        from the GFF.
    output_prefix : str
        Output paths will be ``<output_prefix>.{png,svg,pdf}``.
    top_n : int
        Number of largest distinct clusters to plot in panel B.
    """
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11, 4.5))

    sizes = pd.to_numeric(df["tandem_cluster_size"], errors="coerce").dropna().astype(int)

    # ---- Panel A: size histogram ----
    if len(sizes) == 0:
        ax_a.text(0.5, 0.5, "No tandem-cluster data",
                  ha="center", va="center", transform=ax_a.transAxes,
                  color="gray", fontsize=12)
        _setup_axes(ax_a, "Tandem-cluster size distribution", "size", "candidates")
    else:
        max_size = int(sizes.max())
        # Bin edges 1, 2, 3, ..., max+1; show "1" (singletons) explicitly so
        # the reader sees what's NOT in a cluster.
        bins = list(range(1, max_size + 2))
        ax_a.hist(sizes, bins=bins, color="#3a7ca5", edgecolor="white", align="left")
        ax_a.set_xticks(range(1, max_size + 1))
        _setup_axes(ax_a, "Tandem-cluster size distribution",
                    "cluster size (genes)", "candidates")

    # ---- Panel B: top-N largest distinct clusters ----
    has_cid = (df["tandem_cluster_id"].astype(str).str.len() > 0)
    clustered = df[has_cid].copy()
    if clustered.empty:
        ax_b.text(0.5, 0.5, "No tandem clusters detected",
                  ha="center", va="center", transform=ax_b.transAxes,
                  color="gray", fontsize=12)
        _setup_axes(ax_b, f"Top {top_n} largest tandem clusters", "cluster", "size")
    else:
        per_cluster = (clustered
                       .groupby("tandem_cluster_id")["tandem_cluster_size"]
                       .max()
                       .sort_values(ascending=False)
                       .head(top_n))
        ax_b.barh(range(len(per_cluster)), per_cluster.values,
                  color="#c14b51", edgecolor="white")
        ax_b.set_yticks(range(len(per_cluster)))
        ax_b.set_yticklabels(per_cluster.index, fontsize=8)
        ax_b.invert_yaxis()
        _setup_axes(ax_b, f"Top {top_n} largest tandem clusters",
                    "size (genes)", "")

    fig.suptitle("Berghia GPCR tandem-cluster detection", fontsize=13, y=1.02)
    fig.tight_layout()

    out = Path(output_prefix)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg", "pdf"):
        target = out.with_suffix(f".{ext}") if out.suffix else Path(f"{output_prefix}.{ext}")
        fig.savefig(target, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("csv", help="Path to tandem_clusters.csv")
    ap.add_argument("output_prefix",
                    help="Output prefix (no extension); will write .png/.svg/.pdf")
    ap.add_argument("--top-n", type=int, default=15,
                    help="Number of largest clusters to show in panel B (default 15)")
    args = ap.parse_args()

    csv_path = Path(args.csv)
    if not csv_path.exists():
        print(f"WARN: {csv_path} not found; emitting empty figure.", file=sys.stderr)
        df = pd.DataFrame(columns=["candidate_id", "tandem_cluster_size", "tandem_cluster_id"])
    else:
        df = pd.read_csv(csv_path)
        for col in ("candidate_id", "tandem_cluster_size", "tandem_cluster_id"):
            if col not in df.columns:
                print(f"ERROR: missing column '{col}' in {csv_path}", file=sys.stderr)
                return 2

    plot(df, args.output_prefix, top_n=args.top_n)
    print(f"Wrote {args.output_prefix}.{{png,svg,pdf}}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
