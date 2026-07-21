#!/usr/bin/env python3
"""Measure what SIGNAL CORRELATION actually does to the production RRA order.

Companion to scripts/audit_signal_ranking_independence.py, which answers "are
the per-signal ranklists correlated?" This script answers the follow-up an
expert council raised and the correlation matrix alone cannot settle: "and does
that correlation change the shortlist?"

Two blind spots in the raw correlation audit motivate this tool:

1. WRONG COLUMNS. The correlation audit reads the raw ``<signal>_score``
   columns. Rank aggregation votes on ``<signal>_score_norm`` where present.
   rank_candidates.py applies the dN/dS ``reliability_shrink`` to
   ``positive_score_norm`` / ``purifying_score_norm`` ONLY -- shrinking both
   selection signals toward their cohort median by the SAME per-candidate
   weight ``min(1, n_ref_cds / DNDS_RELIABILITY_FULL)``. That shared multiplier
   is a genuine correlation SOURCE between the two selection axes, and it is
   invisible in the raw columns. This script builds its matrix from
   rank_aggregation.build_ranklists_from_df(), i.e. the exact values RRA votes
   on -- _norm preference, has_*_data gating and exclusion-signal inversion all
   resolved in one place, so the two audits cannot drift apart.

2. NO CONSEQUENCE. |rho| >= 0.7 is an unfalsifiable grouping threshold on its
   own. What matters is whether fusing the flagged group reorders the top of
   the list. `fusion_impact` reports that directly (order Spearman, top-K
   Jaccard, max displacement), so the threshold can be argued from effect size.

It also reports TIE SATURATION. Every candidate sharing an RRA score is exactly
TIED, and aggregate() breaks ties by ascending id -- alphabetically, not by
evidence. A large tied block means that stretch of the "ranking" is not a
ranking at all, and correlated inputs change which candidates escape it.

This is the metric bead 8k8e was opened on. rra_score used to return the
Bonferroni bound ``min(rho * m, 1.0)``, which clipped 204 of the real 439
candidates (46.5%) into ONE tie block flattening 112 distinct rho values. The
exact null replaced it, and this report is now the standing check that no such
block returns: ``n_saturated`` (score == 1.0) should be ~0, and the residual
tied blocks should each hold a single shared rho -- candidates that agree on
the one order statistic RRA looks at, which no exactness can or should
separate.

Usage (once a REAL ranked CSV exists -- see the module note below):

    python3 scripts/audit_rra_correlation_sensitivity.py \\
        --ranked-csv results/ranking/ranked_candidates_sorted.csv \\
        --out-prefix results/ranking/rra \\
        [--groups-json results/ranking/signal_independence_groups.json] \\
        [--threshold 0.7] [--top-k 20]

With no --groups-json, the flagged pairs are grouped by single-linkage at
--threshold and that grouping is what fusion_impact evaluates.

NOTE ON DATA: as of 2026-07 no trustworthy ranked CSV exists (stage 05 never
ran, LSE classification never worked, several axes dormant). Numbers produced
from the preliminary CSV are MEANINGLESS and must not be interpreted; this tool
exists to be run the moment a real stage-07 run completes.
"""
from __future__ import annotations

import argparse
import json
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from rank_aggregation import aggregate, build_ranklists_from_df, rra_score


def signal_matrix(df, id_col="id"):
    """Per-candidate matrix of the values RRA actually votes on.

    Columns are SIGNAL_SPEC signal keys; rows are ``df``'s ids in order. NaN
    marks "this signal does not vote for this candidate" (absent column,
    has_*_data False, or a missing value). Delegates all column resolution to
    build_ranklists_from_df so this audit and the ranker can never disagree
    about which column, polarity or gate a signal uses.
    """
    ranklists = build_ranklists_from_df(df, id_col=id_col)
    ids = list(df[id_col])
    return pd.DataFrame(
        {name: [scores.get(i, np.nan) for i in ids]
         for name, scores in ranklists.items()},
        index=ids,
    )


def rank_correlation_matrix(matrix):
    """Spearman rank correlation between every pair of signal columns.

    Each pair uses its own pairwise-complete rows (a signal only votes where it
    has data, so dropping rows globally would discard most of the matrix).
    Pairs with fewer than 3 shared rows, or with no variance, yield 0.0 --
    "no measurable dependence", never NaN, so downstream thresholding is total.
    """
    cols = list(matrix.columns)
    out = pd.DataFrame(np.eye(len(cols)), index=cols, columns=cols)
    for a, b in combinations(cols, 2):
        pair = matrix[[a, b]].dropna()
        rho = 0.0
        if len(pair) >= 3:
            value = spearmanr(pair[a], pair[b]).correlation
            if value is not None and not np.isnan(value):
                rho = float(value)
        out.loc[a, b] = out.loc[b, a] = rho
    return out


def flag_correlated_pairs(corr, threshold=0.7):
    """Pairs with ``|rho| >= threshold``, strongest first.

    Returns ``[(signal_a, signal_b, rho), ...]``. Keys on the ABSOLUTE value:
    an exclusion signal that mirrors another axis is as redundant as one that
    tracks it.
    """
    pairs = [(a, b, float(corr.loc[a, b]))
             for a, b in combinations(list(corr.columns), 2)
             if abs(corr.loc[a, b]) >= threshold]
    return sorted(pairs, key=lambda p: -abs(p[2]))


def groups_from_pairs(pairs, columns):
    """Single-linkage grouping of the flagged pairs (union-find).

    Mirrors audit_signal_ranking_independence.group_correlated_signals so the
    grouping this tool evaluates is the grouping the ranker would receive.
    Signals in no flagged pair come back as singletons.
    """
    parent = {c: c for c in columns}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b, _ in pairs:
        if a in parent and b in parent:
            parent[find(a)] = find(b)
    groups = {}
    for c in columns:
        groups.setdefault(find(c), []).append(c)
    return [sorted(g) for g in sorted(groups.values(), key=lambda g: sorted(g)[0])]


def rra_with(ranklists, groups=None):
    """RRA score per id under an explicit grouping (thin, named wrapper)."""
    return rra_score(ranklists, groups=groups)


def saturation_report(ranklists, groups=None):
    """How much of the order is decided by TIES rather than by evidence.

    ``aggregate`` breaks equal scores by ascending id, i.e. alphabetically, so
    every exactly-tied block is a stretch of the "ranking" that carries no
    information. Reports the largest such block plus how many candidates sit in
    any block of size > 1.

    ``n_saturated`` (score == 1.0) is retained as a REGRESSION detector. Under
    the Bonferroni bound this module was written against it counted 46.5% of
    the real cohort; under the exact null (bead 8k8e) a score of exactly 1.0
    requires rho == 1, i.e. ranked dead last in every list it appears in, so a
    large value here now means the clip has come back.
    """
    scores = rra_score(ranklists, groups=groups)
    n = len(scores)
    values = list(scores.values())
    n_saturated = sum(1 for v in values if v >= 1.0)
    counts = pd.Series(values).value_counts() if values else pd.Series(dtype=int)
    n_tied = int(counts[counts > 1].sum()) if len(counts) else 0
    return {
        "n_candidates": n,
        "n_saturated": n_saturated,
        "fraction_saturated": (n_saturated / n) if n else 0.0,
        "largest_tie_block": int(counts.iloc[0]) if len(counts) else 0,
        "n_tied": n_tied,
        "fraction_tied": (n_tied / n) if n else 0.0,
        "n_distinct_scores": int(len(counts)),
    }


def fusion_impact(ranklists, groups, top_k=20):
    """Effect of fusing ``groups`` on the RRA order: ungrouped vs grouped.

    All three metrics are computed over the ids ranked by BOTH orders:
      order_spearman   -- rank correlation of the two position vectors
                          (1.0 = fusing changed nothing)
      top_k_jaccard    -- overlap of the two top-K sets
      n_top_k_changed  -- ids entering or leaving the top K
      max_displacement -- largest absolute position shift of any candidate
    """
    ungrouped = aggregate(ranklists, method="rra", groups=None)
    grouped = aggregate(ranklists, method="rra", groups=groups)
    pos_u = {i: n for n, i in enumerate(ungrouped)}
    pos_g = {i: n for n, i in enumerate(grouped)}
    shared = [i for i in ungrouped if i in pos_g]

    displacements = [abs(pos_u[i] - pos_g[i]) for i in shared]
    top_u, top_g = set(ungrouped[:top_k]), set(grouped[:top_k])
    union = top_u | top_g

    order_spearman = 1.0
    if len(shared) >= 3:
        value = spearmanr([pos_u[i] for i in shared],
                          [pos_g[i] for i in shared]).correlation
        if value is not None and not np.isnan(value):
            order_spearman = float(value)

    return {
        "top_k": top_k,
        "order_spearman": order_spearman,
        "top_k_jaccard": (len(top_u & top_g) / len(union)) if union else 1.0,
        "n_top_k_changed": len(top_u ^ top_g),
        "max_displacement": max(displacements) if displacements else 0,
    }


def analyze(df, threshold=0.7, top_k=20, groups=None, id_col="id"):
    """Full audit of one ranked dataframe. ``groups`` overrides the derived grouping."""
    matrix = signal_matrix(df, id_col=id_col)
    corr = rank_correlation_matrix(matrix)
    pairs = flag_correlated_pairs(corr, threshold)
    derived = groups_from_pairs(pairs, list(matrix.columns))
    effective = groups if groups is not None else derived
    ranklists = build_ranklists_from_df(df, id_col=id_col)
    return {
        "threshold": threshold,
        "signals": list(matrix.columns),
        "correlation": corr,
        "flagged_pairs": pairs,
        "groups": effective,
        "saturation": saturation_report(ranklists),
        "saturation_grouped": saturation_report(ranklists, groups=effective),
        "fusion_impact": fusion_impact(ranklists, effective, top_k=top_k),
    }


def _interpret_fusion_impact(imp) -> str:
    """State what the measurement shows, derived from the measurement.

    This sentence used to be a hardcoded claim that the correlation was inert
    -- printed unconditionally, four lines below the numbers, and on the real
    cohort those numbers were Spearman 0.6122 with 380 positions of
    displacement. A conclusion that cannot be false is not a conclusion, so it
    is now computed from the same values the block above prints.
    """
    spearman = float(imp["order_spearman"])
    displacement = int(imp["max_displacement"])
    changed = int(imp["n_top_k_changed"])
    inert = spearman >= 0.9999 and displacement == 0

    if inert:
        return (
            f"An order Spearman of {spearman:.4f} with {displacement} displacement "
            "means the measured correlation does NOT change the ranking, so the "
            "double-counting objection, while structurally real, is inert here."
        )
    return (
        f"Fusing the flagged groups DOES change the ranking: order Spearman "
        f"{spearman:.4f}, max displacement {displacement} positions, and "
        f"{changed} of the top-{imp['top_k']} ids change. The double-counting "
        "objection is NOT inert on this data -- the correlated signals are "
        "moving the order, so the grouped and ungrouped rankings are different "
        "answers and one of them has to be chosen deliberately."
    )


def write_report(result, out_prefix):
    """Emit ``<prefix>_correlation.tsv``, ``_sensitivity.json`` and ``_sensitivity.md``."""
    result["correlation"].to_csv(out_prefix + "_correlation.tsv", sep="\t")

    payload = {k: v for k, v in result.items() if k != "correlation"}
    payload["flagged_pairs"] = [[a, b, rho] for a, b, rho in result["flagged_pairs"]]
    with open(out_prefix + "_sensitivity.json", "w") as fh:
        json.dump(payload, fh, indent=2)

    sat, imp = result["saturation"], result["fusion_impact"]
    lines = [
        "# RRA correlation-sensitivity audit",
        "",
        f"Signals audited: {len(result['signals'])}",
        f"Flag threshold: |rho| >= {result['threshold']}",
        "",
        "## Flagged pairs",
    ]
    lines += ([f"- {a} <-> {b}: rho = {rho:+.3f}"
               for a, b, rho in result["flagged_pairs"]]
              or ["- none above threshold"])
    lines += [
        "",
        "## Tie saturation (ungrouped)",
        f"- candidates: {sat['n_candidates']}",
        f"- distinct RRA scores: {sat['n_distinct_scores']}",
        f"- in an exactly-tied block: {sat['n_tied']} "
        f"({sat['fraction_tied']:.1%}) "
        "-- ties break alphabetically by id, NOT by evidence",
        f"- largest exactly-tied block: {sat['largest_tie_block']}",
        f"- score == 1.0 (regression check, expect ~0 under the exact null): "
        f"{sat['n_saturated']} ({sat['fraction_saturated']:.1%})",
        "",
        "## Effect of fusing the flagged groups",
        f"- order Spearman (grouped vs ungrouped): {imp['order_spearman']:.4f}",
        f"- top-{imp['top_k']} Jaccard: {imp['top_k_jaccard']:.3f}",
        f"- ids entering/leaving top-{imp['top_k']}: {imp['n_top_k_changed']}",
        f"- max position displacement: {imp['max_displacement']}",
        "",
        "",
        _interpret_fusion_impact(imp),
    ]
    with open(out_prefix + "_sensitivity.md", "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ranked-csv", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--groups-json", default=None,
                    help="signal_independence_groups.json; default derives "
                         "groups from the flagged pairs at --threshold")
    ap.add_argument("--threshold", type=float, default=0.7)
    ap.add_argument("--top-k", type=int, default=20)
    args = ap.parse_args(argv)

    groups = None
    if args.groups_json:
        with open(args.groups_json) as fh:
            data = json.load(fh)
        raw = data.get("groups") if isinstance(data, dict) else data
        if raw:
            from rank_aggregation import normalize_group_names
            groups = normalize_group_names(raw)

    df = pd.read_csv(args.ranked_csv)
    result = analyze(df, threshold=args.threshold, top_k=args.top_k,
                     groups=groups)
    write_report(result, args.out_prefix)

    sat, imp = result["saturation"], result["fusion_impact"]
    print(f"[rra-audit] {len(result['signals'])} signals, "
          f"{len(result['flagged_pairs'])} flagged pair(s); "
          f"saturated {sat['fraction_saturated']:.1%}; "
          f"fusion max displacement {imp['max_displacement']}")


if __name__ == "__main__":
    main()
