#!/usr/bin/env python3
# shortlist_impact.py
# Purpose: Read-only diagnostic -- does the embedding-novelty voter actually
#          change the wet-lab top-k shortlist under Robust Rank Aggregation?
# Author: Katz Lab, University of Massachusetts, Amherst

"""Does the embedding-novelty axis move the shortlist?

The production candidate order is a label-free Robust Rank Aggregation (RRA)
fusion of ~12 per-signal ranklists (``rank_aggregation.aggregate``). One of the
voters is ``emb_novelty`` (embedding novelty -- how divergent a candidate is
from known class-A families). This tool answers, for a given ranked CSV,
whether that voter is load-bearing for the shortlist: it recomputes the RRA
order WITH every signal vs WITHOUT the toggled voter and reports how much the
top-k set changes -- Jaccard overlap, how many candidates the voter pulls into
/ pushes out of the top-k, and the biggest rank movers.

It never re-ranks the production output; it only compares two orders. Swapping
raw ``emb_novelty`` for a future length-deconfounded ``emb_novelty_residual``
is the same question with a different column, so ``--signal`` is parameterised.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rank_aggregation import aggregate, build_ranklists_from_df  # noqa: E402


def _positions(order):
    """1-based position lookup for an ordering (best = position 1)."""
    return {idv: i + 1 for i, idv in enumerate(order)}


def topk_jaccard(order_a, order_b, k):
    """Jaccard overlap of the two top-k id sets.

    ``|A ∩ B| / |A ∪ B|`` over the first ``k`` ids of each order. Identical
    top-k sets give 1.0, disjoint give 0.0. An empty union (``k <= 0`` or both
    orders empty) gives 1.0 -- there is no shortlist to disturb.
    """
    a = set(order_a[:k])
    b = set(order_b[:k])
    union = a | b
    if not union:
        return 1.0
    return len(a & b) / len(union)


def topk_movement(order_a, order_b, k):
    """How the top-k changes going from ``order_a`` to ``order_b``.

    Returns a dict:
      - ``n_entered`` / ``entered``: count + sorted ids in ``order_b``'s top-k
        but not ``order_a``'s (what moving a -> b *pulled in*).
      - ``n_left`` / ``left``: count + sorted ids in ``order_a``'s top-k but not
        ``order_b``'s (what it *pushed out*).
      - ``movers``: every id in the union of the two top-k sets as
        ``(id, rank_a, rank_b)`` (1-based positions in each order), sorted by
        absolute rank movement descending, then id ascending. An id absent from
        an order is scored at ``len(order)+1`` (worst) for that side.
    """
    a_set = set(order_a[:k])
    b_set = set(order_b[:k])
    entered = sorted(b_set - a_set)
    left = sorted(a_set - b_set)

    pa, pb = _positions(order_a), _positions(order_b)
    worst_a, worst_b = len(order_a) + 1, len(order_b) + 1
    movers = [
        (idv, pa.get(idv, worst_a), pb.get(idv, worst_b))
        for idv in (a_set | b_set)
    ]
    movers.sort(key=lambda t: (-abs(t[1] - t[2]), t[0]))
    return {
        "n_entered": len(entered),
        "n_left": len(left),
        "entered": entered,
        "left": left,
        "movers": movers,
    }


def shortlist_impact(df, signal_to_toggle="emb_novelty", k=20):
    """Impact of one voter on the top-k RRA shortlist.

    Builds the per-signal ranklists from ``df``
    (``rank_aggregation.build_ranklists_from_df``), computes the RRA order with
    every signal (the production order) and again with ``signal_to_toggle``
    removed (the counterfactual), and reports the top-k Jaccard plus the
    entered/left/mover breakdown attributable to the toggled voter.

    ``order_a`` is the counterfactual (voter removed) and ``order_b`` is the
    full order, so ``n_entered`` counts candidates the voter pulls INTO the
    top-k and ``n_left`` counts those it pushes OUT.

    If ``signal_to_toggle`` is not among the ranklists (its column/flag is
    absent, or every value was missing), it cannot move anything: Jaccard is
    1.0 and a ``note`` explains why.
    """
    ranklists = build_ranklists_from_df(df)
    order_full = aggregate(ranklists, method="rra")

    if signal_to_toggle not in ranklists:
        return {
            "jaccard": 1.0,
            "movement": {"n_entered": 0, "n_left": 0,
                         "entered": [], "left": [], "movers": []},
            "signal": signal_to_toggle,
            "signal_present": False,
            "note": (f"signal {signal_to_toggle!r} is absent from the ranklists "
                     "(its column/flag is missing or every value was empty); it "
                     "cannot change the shortlist, so Jaccard = 1.0."),
        }

    ablated = {name: scores for name, scores in ranklists.items()
               if name != signal_to_toggle}
    order_ablated = aggregate(ablated, method="rra")

    return {
        "jaccard": topk_jaccard(order_full, order_ablated, k),
        "movement": topk_movement(order_ablated, order_full, k),
        "signal": signal_to_toggle,
        "signal_present": True,
    }


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _format_report(res, k, signal):
    lines = [f"# shortlist-impact: does {signal!r} change the top-{k}?",
             f"top-{k} Jaccard (with vs without {signal}): {res['jaccard']:.4f}"]
    if not res.get("signal_present", True):
        lines.append(res.get("note", ""))
        return "\n".join(lines)
    mv = res["movement"]
    lines.append(f"pulled INTO the top-{k} by {signal}: "
                 f"{mv['n_entered']} {mv['entered']}")
    lines.append(f"pushed OUT of the top-{k} by {signal}: "
                 f"{mv['n_left']} {mv['left']}")
    lines.append("biggest rank movers, union of the two top-k "
                 "(id: without-rank -> with-rank):")
    for idv, rank_without, rank_with in mv["movers"][:10]:
        lines.append(f"  {idv}: {rank_without} -> {rank_with} "
                     f"(|delta|={abs(rank_without - rank_with)})")
    return "\n".join(lines)


def main(argv=None):
    import argparse

    import pandas as pd

    ap = argparse.ArgumentParser(
        description="Does the embedding-novelty voter change the top-k RRA "
                    "shortlist?")
    ap.add_argument("--ranked-csv", required=True,
                    help="ranked candidate CSV with the per-signal columns")
    ap.add_argument("--signal", default="emb_novelty",
                    help="signal/voter to toggle off (default: emb_novelty)")
    ap.add_argument("--k", type=int, default=20, help="shortlist size")
    a = ap.parse_args(argv)

    df = pd.read_csv(a.ranked_csv)
    res = shortlist_impact(df, signal_to_toggle=a.signal, k=a.k)
    print(_format_report(res, a.k, a.signal))
    return res


if __name__ == "__main__":
    main()
