#!/usr/bin/env python3
# compare_ranking_methods.py
# Purpose: Honestly compare the production weighted-sum ranking against the
#          label-free rank-aggregation ranking, WITHOUT changing the shortlist.
# Author: Katz Lab, University of Massachusetts, Amherst

"""Weighted-sum vs rank-aggregation candidate ranking -- a non-disruptive audit.

The production shortlist still comes from the hand-weighted sum
(``RANK_METHOD=weighted``). This module recomputes the label-free Robust Rank
Aggregation ordering over the SAME 12 per-signal ranklists (via the shared
``rank_aggregation.build_ranklists_from_df``) and reports how far the two
orderings agree: Spearman correlation, top-k overlap, and the biggest movers
each direction.

Honest uncertainty. There are ~0 in-species validated positives, so any
"enrichment" statement here is DESCRIPTIVE only: we report the raw rank
positions of whatever positive controls map into the candidate set plus a
permutation-null p-value, and we deliberately REFUSE to print a precision@k
estimate when fewer than 5 positives are available (it would be sampling noise
dressed up as a metric). ``references/hcr_positive_controls.csv`` is ~empty
today, so the mapper is best-effort and tolerates 0 matches.
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rank_aggregation import (  # noqa: E402
    aggregate,
    build_ranklists_from_df,
    normalize_group_names,
    production_excluded_signals,
)

MIN_POSITIVES_FOR_PRECISION = 5


# --------------------------------------------------------------------------- #
# Orderings
# --------------------------------------------------------------------------- #
def weighted_order(df, id_col="id", score_col="rank_score"):
    """Reconstruct the weighted ordering (best first) from the signal df.

    Mirrors exactly what rank_candidates.py does:
    ``df.sort_values('rank_score', ascending=False)``. Reconstructing it from
    the score column (rather than trusting the CSV's row order) keeps this
    correct even when the CSV was written in rankagg order.
    """
    return df.sort_values(score_col, ascending=False)[id_col].tolist()


def rankagg_order(df, groups=None, method="rra", id_col="id", excluded=None):
    """Label-free rank-aggregation ordering (best first) over the 12 signals.

    Any id covered by no signal (none should exist -- the 4 base signals are
    always present) is appended at the tail in df order. ``excluded`` names
    signals barred from voting (a zero-weight or orthology-quarantined axis,
    mirroring the production ranker -- bead od2f); the default ``None`` excludes
    nothing.
    """
    order = aggregate(build_ranklists_from_df(df, id_col=id_col, excluded=excluded),
                      method=method, groups=groups)
    covered = set(order)
    tail = [i for i in df[id_col].tolist() if i not in covered]
    return order + tail


# --------------------------------------------------------------------------- #
# Pure ordering primitives (testable on explicit toy orderings)
# --------------------------------------------------------------------------- #
def spearman_of_orders(order_a, order_b):
    """Spearman rank correlation between two orderings over their shared ids."""
    b_set = set(order_b)
    common = [i for i in order_a if i in b_set]
    if len(common) < 3:
        return float("nan")
    pa = {v: i for i, v in enumerate(order_a)}
    pb = {v: i for i, v in enumerate(order_b)}
    rho = spearmanr([pa[i] for i in common], [pb[i] for i in common]).correlation
    return float(rho)


def topk_overlap(order_a, order_b, ks=(20, 50)):
    """Fraction of the top-k set shared between two orderings, per k."""
    out = {}
    for k in ks:
        a = set(order_a[:k])
        b = set(order_b[:k])
        denom = min(k, len(order_a), len(order_b))
        out[k] = (len(a & b) / denom) if denom else float("nan")
    return out


def biggest_movers(order_a, order_b, n=10):
    """The n ids that move most between order_a (weighted) and order_b (rankagg).

    Returns ``(up, down)`` lists of ``(id, delta, pos_a, pos_b)`` (1-based
    positions). ``delta = pos_a - pos_b`` > 0 means the id moved UP (nearer the
    top) under order_b.
    """
    pa = {v: i for i, v in enumerate(order_a)}
    pb = {v: i for i, v in enumerate(order_b)}
    common = [i for i in order_a if i in pb]
    movers = [(i, pa[i] - pb[i], pa[i] + 1, pb[i] + 1) for i in common]
    up = sorted(movers, key=lambda t: (-t[1], t[0]))[:n]
    down = sorted(movers, key=lambda t: (t[1], t[0]))[:n]
    return up, down


def compare(df, groups=None, method="rra", id_col="id", top_ks=(20, 50),
            excluded=None):
    """Compare the weighted vs rank-aggregation orderings for a signal df.

    ``excluded`` is threaded into the rank-aggregation ordering only (a
    zero-weight/quarantined axis must not vote there -- bead od2f); the weighted
    ordering is read straight from ``rank_score`` and is unaffected.
    """
    w_order = weighted_order(df, id_col=id_col)
    r_order = rankagg_order(df, groups=groups, method=method, id_col=id_col,
                            excluded=excluded)
    up, down = biggest_movers(w_order, r_order, n=10)
    return {
        "weighted_order": w_order,
        "rankagg_order": r_order,
        "spearman": spearman_of_orders(w_order, r_order),
        "topk_overlap": topk_overlap(w_order, r_order, ks=top_ks),
        "movers_up": up,
        "movers_down": down,
        "method": method,
        "grouped": bool(groups),
    }


# --------------------------------------------------------------------------- #
# Honest uncertainty
# --------------------------------------------------------------------------- #
def enrichment_vs_null(order, positive_ids, n_perm=10000, seed=0):
    """Descriptive permutation-null enrichment of positives near the top.

    Reports each present positive's 1-based rank position and a permutation
    p-value for the mean position: ``p = P(mean rank of |positives| uniformly
    random items <= observed mean rank)`` under add-one smoothing (never 0).
    LOWER p = positives cluster nearer the top than chance. DESCRIPTIVE only --
    with ~0 validated positives this is not a trained-model metric.
    """
    n = len(order)
    pos_index = {idv: i + 1 for i, idv in enumerate(order)}  # 1-based
    present = [p for p in positive_ids if p in pos_index]
    positions = {p: pos_index[p] for p in present}
    k = len(present)
    result = {
        "n_items": n, "n_positives": k, "positions": positions,
        "observed_mean_rank": None, "p_value": None, "n_perm": n_perm,
    }
    if n == 0 or k == 0:
        return result
    observed = sum(positions.values()) / k
    result["observed_mean_rank"] = observed
    if k >= n:
        result["p_value"] = 1.0
        return result
    rng = np.random.default_rng(seed)
    rand = rng.random((n_perm, n))
    # the k smallest random draws per row = a uniform random k-subset of ranks
    idx = np.argpartition(rand, k, axis=1)[:, :k]
    means = (idx + 1).mean(axis=1)
    count = int((means <= observed).sum())
    result["p_value"] = (count + 1) / (n_perm + 1)
    return result


def _precision_at_k(order, positive_ids, k):
    pos = set(positive_ids)
    topk = order[:k]
    denom = min(k, len(order))
    return (sum(1 for i in topk if i in pos) / denom) if denom else float("nan")


def precision_at_k_vs_null(order, positive_ids, k, n_perm=10000, seed=0):
    """Precision@k of positives + permutation-null p (random control-sets of the
    same size). Returns (precision, p) or (None, None) if no positives are in
    ``order``."""
    pos = [i for i in positive_ids if i in order]
    if not pos:
        return (None, None)
    topk = set(order[:k])
    obs = sum(1 for i in pos if i in topk) / k
    rng = np.random.default_rng(seed)
    n, m = len(order), len(pos)
    ge = 0
    for _ in range(n_perm):
        draw = set(rng.choice(n, m, replace=False))
        prec = sum(1 for idx in range(k) if idx in draw) / k
        if prec >= obs:
            ge += 1
    return (obs, (ge + 1) / (n_perm + 1))


def loo_recovery(positive_ids, recompute_fn):
    """For each positive, the 1-based rank it recovers to when the ranking is
    recomputed with it held out. recompute_fn(id) -> best-first order list."""
    out = {}
    for pid in positive_ids:
        order = recompute_fn(pid)
        out[pid] = order.index(pid) + 1 if pid in order else None
    return out


def map_positive_controls(controls_csv, candidate_ids):
    """Best-effort map of hcr_positive_controls.csv rows to candidate ids.

    Matches (case-insensitively) each control's refseq_protein / gene_name /
    aliases against the candidate id set: exact match first, then a
    conservative substring match for longer (>=6 char, accession-like) tokens.
    Tolerates 0 matches -- the controls file is ~empty today. Returns a
    de-duplicated list of matched candidate ids.
    """
    ids = list(candidate_ids)
    exact = {str(c).lower(): c for c in ids}
    matched = []
    try:
        with open(controls_csv, newline="") as fh:
            for row in csv.DictReader(fh):
                tokens = []
                for field in ("refseq_protein", "gene_name", "aliases"):
                    val = (row.get(field) or "").strip()
                    if val:
                        tokens += [t.strip() for t in val.replace(";", ",").split(",")
                                   if t.strip()]
                for tok in tokens:
                    tl = tok.lower()
                    if tl in exact and exact[tl] not in matched:
                        matched.append(exact[tl])
                        continue
                    if len(tl) >= 6:
                        hit = next((c for c in ids if tl in str(c).lower()), None)
                        if hit and hit not in matched:
                            matched.append(hit)
    except (OSError, ValueError):
        return []
    return matched


def _fmt_p(p):
    return "n/a" if p is None else f"{p:.4g}"


def _fmt_prec(x):
    return "n/a" if x is None else f"{x:.3f}"


def write_report(df, out_path, positive_ids=None, groups=None, method="rra",
                 n_perm=10000, seed=0, id_col="id", recompute_fn=None,
                 excluded=None):
    """Write the weighted-vs-rankagg comparison markdown; return its text.

    ``excluded`` (default ``None`` = exclude nothing) is threaded into the
    rank-aggregation ordering so a zero-weight/quarantined axis does not vote
    there -- bead od2f.
    """
    positive_ids = list(positive_ids or [])
    res = compare(df, groups=groups, method=method, id_col=id_col,
                  excluded=excluded)
    n_pos = len(positive_ids)

    L = []
    L.append("# Ranking-method comparison: weighted sum vs rank aggregation")
    L.append("")
    L.append("**Non-disruptive audit.** The production shortlist is still produced")
    L.append("by the hand-weighted sum (`RANK_METHOD=weighted`); switching it is a")
    L.append("separate, gated decision. Below, the label-free reranker is Robust")
    L.append(f"Rank Aggregation (`method={method}`) over the 12 per-signal ranklists, "
             + ("WITH" if groups else "without") + " confound grouping.")
    L.append("")
    L.append("## Agreement")
    L.append("")
    L.append(f"- Spearman rank correlation (weighted vs rankagg): {res['spearman']:.4f}")
    for k, ov in sorted(res["topk_overlap"].items()):
        L.append(f"- Top-{k} overlap: {ov:.3f}")
    L.append("")
    L.append("## Biggest movers (weighted rank -> rankagg rank)")
    L.append("")
    L.append("Promoted under rank aggregation:")
    L.append("")
    L.append("| id | weighted rank | rankagg rank | delta |")
    L.append("|----|--------------:|-------------:|------:|")
    for idv, delta, wpos, rpos in res["movers_up"]:
        L.append(f"| {idv} | {wpos} | {rpos} | +{delta} |")
    L.append("")
    L.append("Demoted under rank aggregation:")
    L.append("")
    L.append("| id | weighted rank | rankagg rank | delta |")
    L.append("|----|--------------:|-------------:|------:|")
    for idv, delta, wpos, rpos in res["movers_down"]:
        L.append(f"| {idv} | {wpos} | {rpos} | {delta} |")
    L.append("")
    L.append("## Positive-control enrichment (descriptive only)")
    L.append("")
    if n_pos == 0:
        L.append("No positive-control ids were provided or matched")
        L.append("(`references/hcr_positive_controls.csv` is ~empty today), so no")
        L.append("enrichment is computed. This is expected, not an error.")
    else:
        w_enr = enrichment_vs_null(res["weighted_order"], positive_ids,
                                   n_perm=n_perm, seed=seed)
        r_enr = enrichment_vs_null(res["rankagg_order"], positive_ids,
                                   n_perm=n_perm, seed=seed)
        L.append(f"{n_pos} positive-control id(s) supplied; {w_enr['n_positives']} "
                 "present in the ranking. Raw 1-based rank positions:")
        L.append("")
        L.append("| positive id | weighted rank | rankagg rank |")
        L.append("|-------------|--------------:|-------------:|")
        for pid in positive_ids:
            wr = w_enr["positions"].get(pid, "absent")
            rr = r_enr["positions"].get(pid, "absent")
            L.append(f"| {pid} | {wr} | {rr} |")
        L.append("")
        L.append(f"- weighted mean-rank permutation-null p = {_fmt_p(w_enr['p_value'])} "
                 f"(n_perm={w_enr['n_perm']})")
        L.append(f"- rankagg  mean-rank permutation-null p = {_fmt_p(r_enr['p_value'])} "
                 f"(n_perm={r_enr['n_perm']})")
        L.append("")
        if n_pos >= MIN_POSITIVES_FOR_PRECISION:
            L.append("### Precision@k vs permutation null")
            L.append("")
            L.append("Null = random control-sets of the same size; "
                     "p via (b+1)/(n+1). LOWER p = positives concentrate in the "
                     "top-k more than chance placement.")
            L.append("")
            L.append("| k | weighted precision@k | null p | rankagg precision@k | null p |")
            L.append("|--:|---------------------:|-------:|--------------------:|-------:|")
            for k in (20, 50):
                wp, wpp = precision_at_k_vs_null(res["weighted_order"], positive_ids,
                                                 k, n_perm=n_perm, seed=seed)
                rp, rpp = precision_at_k_vs_null(res["rankagg_order"], positive_ids,
                                                 k, n_perm=n_perm, seed=seed)
                L.append(f"| {k} | {_fmt_prec(wp)} | {_fmt_p(wpp)} | "
                         f"{_fmt_prec(rp)} | {_fmt_p(rpp)} |")
        else:
            L.append(f"Fewer than {MIN_POSITIVES_FOR_PRECISION} positive controls are "
                     "available, so only the raw positions above are reported:")
            L.append("**insufficient positives for precision estimate.** Any k-wise")
            L.append("precision would be dominated by sampling noise and is omitted.")
        if recompute_fn is not None:
            present = [p for p in positive_ids if p in w_enr["positions"]]
            rec = loo_recovery(present, recompute_fn)
            ranks = [r for r in rec.values() if r is not None]
            L.append("")
            L.append("### Leave-one-control-out recovery")
            L.append("")
            L.append("Each control is held out, the ranking is recomputed, and we")
            L.append("record the 1-based rank the held-out control recovers to.")
            L.append("")
            if ranks:
                L.append("| positive id | recovered rank |")
                L.append("|-------------|---------------:|")
                for pid in present:
                    rr = rec.get(pid)
                    L.append(f"| {pid} | {rr if rr is not None else 'absent'} |")
                L.append("")
                med = float(np.median(ranks))
                L.append(f"- median recovered rank: {med:.1f} (n={len(ranks)})")
                for k in (20, 50):
                    frac = sum(1 for r in ranks if r <= k) / len(ranks)
                    L.append(f"- fraction recovered within top-{k}: {frac:.3f}")
            else:
                L.append("No present positive controls recovered a rank "
                         "(recompute produced no placement); nothing to report.")
    L.append("")

    text = "\n".join(L) + "\n"
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    return text


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _load_groups(groups_json):
    if not groups_json or not os.path.exists(groups_json):
        return None
    try:
        with open(groups_json) as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        return None
    g = data.get("groups") if isinstance(data, dict) else data
    return normalize_group_names(g) if g else None


def main(argv=None):
    import argparse
    import pandas as pd

    ap = argparse.ArgumentParser(description="Compare weighted vs rankagg ranking.")
    ap.add_argument("--ranked-csv", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--groups-json", default=None)
    ap.add_argument("--controls-csv", default=None)
    ap.add_argument("--method", default="rra", choices=["rra", "rrf"])
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)

    df = pd.read_csv(a.ranked_csv)
    groups = _load_groups(a.groups_json)
    positive_ids = []
    if a.controls_csv and os.path.exists(a.controls_csv):
        positive_ids = map_positive_controls(a.controls_csv, df["id"].tolist())
    write_report(df, a.out, positive_ids=positive_ids, groups=groups,
                 method=a.method, n_perm=a.n_perm, seed=a.seed,
                 excluded=production_excluded_signals())
    print(f"[compare] wrote {a.out} (positive controls matched: {len(positive_ids)})")


if __name__ == "__main__":
    main()
