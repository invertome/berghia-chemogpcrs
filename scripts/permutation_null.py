#!/usr/bin/env python3
"""permutation_null.py — whole-pipeline permutation null for the ranking.

Breaks the candidate<->signal association (shuffles each signal's values across
candidates independently), re-runs the SAME rank aggregation, and asks whether
the observed top-k separation exceeds the shuffled null -- the honest "is the
ranking distinguishable from noise?" test for a label-free prioritizer
(North et al. 2002 + Phipson & Smyth 2010 for the never-zero (b+1)/(n+1) rule).
"""
from __future__ import annotations

from typing import Dict, Mapping

import numpy as np

from rank_aggregation import rra_score


def topk_separation(order, agg, k):
    """Mean-|score| gap between the top-k of ``order`` and the rest, measured on
    the RRA significance scale (-log10 of the aggregator score).

    RRA's ``rra_score`` returns a p-value-like quantity that is heavily
    concentrated near 1.0 -- most candidates are unremarkable by construction
    -- so the discrimination between a real ranking and a shuffled null lives
    almost entirely in the near-zero tail. A raw mean gap flattens that tail
    (real and null both differ from 1.0 by essentially the same amount) and
    cannot tell signal from noise; the -log10 significance scale -- the scale
    on which Kolde et al. (2012) report and plot RRA scores -- amplifies it.
    ``abs`` keeps the gap orientation-agnostic (RRA scores are lower=better, so
    top-k sits at higher -log10 than the rest).

    Bead 8k8e made this strictly better-behaved: the score is now the EXACT
    null of rho rather than the ``min(rho*m, 1.0)`` bound, which used to CLIP
    46.5% of the real cohort to exactly 1.0 -- contributing an identical
    -log10 of 0.0 for both the observed and the shuffled runs across nearly
    half the candidates."""
    if len(order) <= k:
        return 0.0
    sig = -np.log10(np.clip([agg[i] for i in order], 1e-300, None))
    return float(abs(sig[k:].mean() - sig[:k].mean()))


def _order_and_sep(per_signal, k):
    agg = rra_score(per_signal)
    order = sorted(agg, key=lambda i: (agg[i], i))     # RRA: ascending rho
    return topk_separation(order, agg, k)


def permutation_null_pvalue(per_signal: Mapping[str, Mapping[str, float]],
                            k: int, n_perm: int = 1000, seed: int = 0) -> float:
    """Empirical p that observed top-k separation exceeds a candidate<->signal-
    shuffled null. Each permutation shuffles every signal's values across
    candidates independently. Returns (b+1)/(n_perm+1)."""
    observed = _order_and_sep(per_signal, k)
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(n_perm):
        permuted: Dict[str, Dict[str, float]] = {}
        for name, smap in per_signal.items():
            ids = list(smap)
            vals = rng.permutation(np.array([smap[i] for i in ids]))
            permuted[name] = {i: float(v) for i, v in zip(ids, vals)}
        if _order_and_sep(permuted, k) >= observed:
            ge += 1
    return (ge + 1) / (n_perm + 1)


def main(argv=None):
    import argparse, json
    import pandas as pd
    from rank_aggregation import build_ranklists_from_df, production_excluded_signals
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--ranked-csv", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--k", type=int, default=20)
    ap.add_argument("--n-perm", type=int, default=1000)
    a = ap.parse_args(argv)
    df = pd.read_csv(a.ranked_csv)
    # Mirror the production ranker: a zero-weight or orthology-quarantined signal
    # must not vote in the null either, or the null measures a signal set the
    # shortlist never used (bead od2f).
    per_signal = build_ranklists_from_df(df, excluded=production_excluded_signals())
    p = permutation_null_pvalue(per_signal, k=a.k, n_perm=a.n_perm, seed=0)
    with open(a.out, "w") as fh:
        json.dump({"observed_separation": _order_and_sep(per_signal, a.k),
                   "p_value": p, "k": a.k, "n_perm": a.n_perm,
                   "n_signals": len(per_signal)}, fh, indent=2)
    print(f"[permutation_null] p={p:.4g} (k={a.k}, n_perm={a.n_perm})")


if __name__ == "__main__":
    main()
