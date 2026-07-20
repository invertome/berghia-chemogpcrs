#!/usr/bin/env python3
"""embedding_channel_report.py — summarize a rebuilt embedding channel.

Bead berghia-chemogpcrs-nsew, deliverable 3. Reads the two channel TSVs the
rebuild produces and prints the distributions needed to judge whether the fix
worked:

  * emb_nonchemo_family distribution over the candidates (consensus channel —
    the PRODUCTION path, what stage 07 actually writes).
  * emb_nonchemo_sim distribution (cosine channel — the only path that emits a
    similarity; the stale 888-row calibration numbers came from here, so this
    is what makes "still saturated?" answerable).
  * The per-candidate joint view, so the strongest non-chemoreceptor
    assignments can be named.

CALIBRATION (the stale 888-row cosine channel, for comparison): 351
class-B-secretin, 273 lipid, 147 class-A-other, 74 peptide, 18 class-C, 14
aminergic, 8 opsin, 3 class-F-frizzled; 883/888 at sim >= 0.90, top ~0.99. A
rebuilt channel that is still saturated, or that still puts hundreds of class-A
candidates onto a class-B secretin prototype, is still broken.
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd


def _dist(series: pd.Series, label: str) -> None:
    s = pd.to_numeric(series, errors="coerce").dropna()
    if s.empty:
        print(f"  {label}: (no numeric values)")
        return
    qs = [0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100]
    print(f"  {label}: n={len(s)} mean={s.mean():.4f} sd={s.std():.4f}")
    print("    percentiles: " + "  ".join(
        f"p{q}={np.percentile(s, q):.4f}" for q in qs))
    for thr in (0.99, 0.95, 0.90, 0.80, 0.70, 0.50):
        n = int((s >= thr).sum())
        print(f"    >= {thr:.2f}: {n:4d} / {len(s)}  ({100*n/len(s):5.1f}%)")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--consensus-tsv", required=True)
    p.add_argument("--cosine-tsv")
    p.add_argument("--top", type=int, default=25)
    a = p.parse_args()

    print("=" * 72)
    print("CONSENSUS CHANNEL (production path — fusion_consensus.py)")
    print("=" * 72)
    con = pd.read_csv(a.consensus_tsv, sep="\t")
    print(f"rows: {len(con)}   columns: {list(con.columns)}")
    if len(con) == 0:
        print("!! HEADER-ONLY TSV — the channel is EMPTY. This is the trap-2 "
              "signature (or a scoring failure). Do not report this as a result.")
        return 1
    if "has_emb_data" in con:
        n_true = int(con["has_emb_data"].astype(str).str.lower().isin(
            ["true", "1"]).sum())
        print(f"has_emb_data True: {n_true}/{len(con)}")
        if n_true == 0:
            print("!! has_emb_data is False for EVERY candidate — trap-2 signature.")
            return 1

    if "emb_nonchemo_family" in con:
        print("\nemb_nonchemo_family distribution (nearest non-chemoreceptor "
              "family prototype):")
        vc = con["emb_nonchemo_family"].fillna("(none)").value_counts()
        for fam, n in vc.items():
            print(f"  {fam:<28} {n:5d}  ({100*n/len(con):5.1f}%)")
    if "emb_novelty" in con:
        print("\nemb_novelty distribution (positive ranking axis):")
        _dist(con["emb_novelty"], "emb_novelty")

    if not a.cosine_tsv:
        return 0

    print()
    print("=" * 72)
    print("COSINE EXCLUSION CHANNEL (emb_nonchemo_sim — calibration-comparable)")
    print("=" * 72)
    cos = pd.read_csv(a.cosine_tsv, sep="\t")
    print(f"rows: {len(cos)}   columns: {list(cos.columns)}")
    if len(cos) == 0:
        print("!! HEADER-ONLY TSV — cosine channel empty.")
        return 1

    if "emb_nonchemo_family" in cos:
        print("\nemb_nonchemo_family distribution:")
        vc = cos["emb_nonchemo_family"].fillna("(none)").value_counts()
        for fam, n in vc.items():
            print(f"  {fam:<28} {n:5d}  ({100*n/len(cos):5.1f}%)")
    if "emb_nonchemo_sim" in cos:
        print("\nemb_nonchemo_sim distribution (SATURATION CHECK — the stale "
              "channel had 883/888 >= 0.90):")
        _dist(cos["emb_nonchemo_sim"], "emb_nonchemo_sim")
    if "emb_classA_sim" in cos:
        print("\nemb_classA_sim distribution (class-A recall centroid):")
        _dist(cos["emb_classA_sim"], "emb_classA_sim")

    if {"emb_nonchemo_sim", "emb_nonchemo_family"} <= set(cos.columns):
        print(f"\nTop {a.top} strongest non-chemoreceptor assignments "
              "(highest emb_nonchemo_sim):")
        top = cos.sort_values("emb_nonchemo_sim", ascending=False).head(a.top)
        merged = top.merge(
            con[[c for c in ("id", "emb_novelty") if c in con.columns]],
            on="id", how="left", suffixes=("", "_con"))
        print(f"  {'id':<34} {'family':<24} {'sim':>7} {'novelty':>9}")
        for _, r in merged.iterrows():
            nov = r.get("emb_novelty")
            nov_s = f"{nov:9.4f}" if pd.notna(nov) else "       NA"
            print(f"  {str(r['id']):<34} {str(r['emb_nonchemo_family']):<24} "
                  f"{r['emb_nonchemo_sim']:7.4f} {nov_s}")

        print("\nPer-family strongest assignment (max emb_nonchemo_sim per family):")
        for fam, g in cos.groupby("emb_nonchemo_family"):
            best = g.loc[pd.to_numeric(g["emb_nonchemo_sim"],
                                       errors="coerce").idxmax()]
            print(f"  {str(fam):<24} n={len(g):<5d} best={best['id']} "
                  f"sim={best['emb_nonchemo_sim']:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
