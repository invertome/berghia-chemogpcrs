#!/usr/bin/env python3
# select_probe_batch.py
# Purpose: Active-learning probe-batch selector for the ML/PLM chemoreceptor
#   ranking reranker (Task 7, docs/plans/2026-07-01-ml-plm-chemoreceptor-
#   ranking.md) -- the bench-ready deliverable: choose the FIRST HCR probe
#   batch to maximize information gain, not pure predicted-positive rank.
"""Select the first HCR probe batch as an ACTIVE-LEARNING round, not a top-N pick.

Council finding (bead berghia-chemogpcrs-875): the real bottleneck blocking a
trained/validated chemoreceptor classifier is the absence of in-species
ground truth (~0 confirmed Berghia positives). The highest-value first
wet-lab round is therefore not "the N candidates our current signals like
best" -- it is the batch that will teach us the most: candidates where the
per-signal evidence DISAGREES (phylogeny says one thing, dN/dS another,
tandem-cluster another...), spread across distinct tandem/orthogroup
clusters so one lucky or unlucky cluster doesn't dominate the round.
Resolving those disagreements at the bench manufactures the first real
in-species positives/negatives -- exactly the labels a future calibrated
model (deferred Phase C of the plan) would need to be trained and validated
on.

THIS BATCH IS NOT "the N most likely chemoreceptors." A candidate can be
selected here specifically because it looks mediocre-but-contested, not
because it looks great across the board. For the conventional top-score
comparison baseline, run with `--mode topscore` (or `mode="topscore"`): it
selects on rank-aggregation standing alone -- still gated through the
non-chemoreceptor classifier and the diversity cap, just with no
disagreement re-ranking. Diffing the two modes' outputs is itself part of
the argument for running the active round.

Pipeline (`main`):
    1. read a ranked CSV (rank_candidates.py output)
    2. `rank_aggregation.build_ranklists_from_df(df)` -> raw per-signal scores
    3. `rank_aggregation.normalized_ranks(...)` on each -> `per_signal_ranks`,
       the ``{signal: {id: normalized_rank_in_(0,1]}}`` structure
       `disagreement_score` and `select_batch` consume (smaller = that
       signal likes this id more, per rank_aggregation's own convention)
    4. `select_batch(...)` -- classifier gate, mode ordering, diversity cap
    5. `write_batch(...)` -- bench-ready TSV with a per-row rationale

Dispersion metric: `disagreement_score` uses the population VARIANCE of a
candidate's per-signal normalized ranks, not entropy. Entropy would require
binning continuous (0, 1] values into discrete buckets, which is unstable
with as few signals as are typically present per candidate (as few as 2,
rarely more than ~17 with every optional evidence channel present) -- there
usually aren't enough draws for a bucket histogram to be a meaningful
density estimate. Variance needs no binning and is exactly zero for a
perfectly consistent candidate: "ranked #1 by one signal and last by another
must score HIGH; ranked consistently must score LOW."

`mode="active"` blend: rather than a hand-tuned weighted sum of standing and
disagreement -- which would reintroduce exactly the kind of hand-picked
weight the rank_aggregation module was built to avoid in the first place
(its own header: "no weights, no labels") -- the blend is a two-stage
SHORTLIST-then-REORDER: take the top ``ACTIVE_TIER_MULTIPLIER * batch_size``
classifier-gated candidates by label-free rank-aggregation standing (a
generous shortlist, not just batch_size-wide, so a genuinely-disagreeing-
but-not-#1 candidate has a chance), then sort *that shortlist* by descending
disagreement. Candidates beyond the shortlist keep their standing-order
position and are only reached if the diversity cap forces the selector past
the shortlist. No new weight parameter is introduced; the only knob is the
shortlist depth.

Diversity cap: applied in BOTH modes (not just "active"), so the
active-vs-topscore comparison isolates the disagreement-reordering effect
rather than incidental cluster clumping. No more than
``ceil(batch_size / n_clusters)`` selected candidates come from any one
`cluster_col` value, where `n_clusters` counts distinct clusters among the
classifier-gated eligible pool. `cluster_col` defaults to
``tandem_cluster_id`` -- the exact tandem-array column rank_candidates.py
writes (its output has no bare ``cluster`` column) -- so the cap spreads
the batch across distinct tandem arrays, the field's signature
chemoreceptor grouping. An empty/NaN `cluster_col` value (an isolated gene,
not part of any tandem array) is never grouped with other isolated genes --
each is its own singleton cluster, so isolated loci are never throttled
against one another; only genuine multi-member tandem clusters are capped.
If `cluster_col` is absent from the input, the cap is skipped entirely
(documented, not an error -- tandem-cluster data is optional evidence
throughout this pipeline).

Classifier gate: candidates whose `classifier_col` (default "classification")
value is "non-chemoreceptor" or "likely-non-chemoreceptor" (the existing
06c/07 consensus classifier's confident calls) are excluded outright, in
both modes. A missing/NaN classification is treated as "not excluded" --
consistent with the classifier's own default of "chemoreceptor-candidate"
for unclassified rows (see add_classification_columns.py).

This module never computes a precision@k or a calibrated P(chemoreceptor)
-- honest-uncertainty rule, per the plan's notes for the implementer: there
are ~0 in-species positives, so no such number would mean anything yet.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rank_aggregation import aggregate, build_ranklists_from_df, normalized_ranks  # noqa: E402

# The consensus classifier's two confident non-chemoreceptor calls (see
# scripts/classify_consensus.py); "chemoreceptor-candidate" (its default) and
# any missing/NaN value are both treated as NOT excluded here.
NON_CHEMORECEPTOR_LABELS = {"non-chemoreceptor", "likely-non-chemoreceptor"}

# Shortlist depth for mode="active": how many classifier-gated candidates
# (by rank-aggregation standing) are eligible to be pulled forward by
# disagreement. See the module docstring's "mode=active blend" section.
ACTIVE_TIER_MULTIPLIER = 3

_ID_COL = "id"


def disagreement_score(per_signal_ranks: Dict[str, Dict[str, float]], cand_id: str) -> float:
    """Dispersion of one candidate's normalized rank across signals.

    ``per_signal_ranks`` is ``{signal: {id: normalized_rank_in_(0,1]}}`` --
    each inner dict is the output of ``rank_aggregation.normalized_ranks`` on
    one signal's raw scores. Returns the population variance (``ddof=0``) of
    the candidate's normalized rank over every signal it appears in. HIGH =
    the signals conflict about this candidate (e.g. #1 in one, last in
    another); LOW = the signals agree.

    A candidate present in fewer than two signals has nothing to disagree
    with -- returns 0.0 rather than an undefined/NaN variance.
    """
    values = [ranks[cand_id] for ranks in per_signal_ranks.values() if cand_id in ranks]
    if len(values) < 2:
        return 0.0
    return float(np.var(values))


def _signal_extremes(per_signal_ranks: Dict[str, Dict[str, float]], cand_id: str):
    """The signal(s) that rank ``cand_id`` best / worst (smallest / largest
    normalized rank). Ties (equal normalized rank) include every tied signal,
    sorted alphabetically for determinism. ``([], [])`` if the candidate is
    in no signal at all.
    """
    vals = {sig: ranks[cand_id] for sig, ranks in per_signal_ranks.items() if cand_id in ranks}
    if not vals:
        return [], []
    best_val = min(vals.values())
    worst_val = max(vals.values())
    best = sorted(s for s, v in vals.items() if v == best_val)
    worst = sorted(s for s, v in vals.items() if v == worst_val)
    return best, worst


def _excluded_ids(df: pd.DataFrame, classifier_col: str, id_col: str = _ID_COL) -> set:
    """Ids whose ``classifier_col`` value is a confident non-chemoreceptor call.

    A ``classifier_col`` absent from ``df`` excludes nothing (the
    classification step is optional evidence, like every other gated
    channel in this pipeline); a missing/NaN value for a present column
    also excludes nothing (``Series.isin`` never matches NaN).
    """
    if classifier_col not in df.columns:
        return set()
    mask = df[classifier_col].isin(NON_CHEMORECEPTOR_LABELS)
    return set(df.loc[mask, id_col])


def _rank_agg_order(df: pd.DataFrame, id_col: str = _ID_COL, method: str = "rra",
                     groups=None) -> List[str]:
    """Label-free rank-aggregation ordering (best first) over every id in ``df``.

    Any id covered by no signal at all (none should exist in production --
    the base signals are always present when their columns exist -- but a
    minimal/test df may omit them) is appended at the tail in sorted id
    order, so every id in ``df`` gets a position.
    """
    raw_ranklists = build_ranklists_from_df(df, id_col=id_col)
    order = aggregate(raw_ranklists, method=method, groups=groups)
    covered = set(order)
    tail = sorted(i for i in df[id_col].tolist() if i not in covered)
    return order + tail


def _apply_diversity_cap(ordered_ids: List[str], df: pd.DataFrame, batch_size: int,
                          cluster_col: str, id_col: str = _ID_COL) -> List[str]:
    """Walk ``ordered_ids`` picking up to ``batch_size``, capping any one
    ``cluster_col`` value at ``ceil(batch_size / n_clusters)`` (n_clusters =
    distinct cluster values among ``ordered_ids``). A missing/NaN cluster
    value is NOT a shared bucket: an isolated gene (not part of a tandem
    array) is its own distinct locus, so each gets its own singleton
    cluster and is never throttled against other isolated genes (a cap is
    always >= 1, and a singleton cluster has only one member to begin
    with). Real (non-empty) tandem clusters are still capped as before. If
    ``cluster_col`` isn't in ``df``, the cap is skipped and this is just a
    truncation to ``batch_size``.
    """
    if cluster_col not in df.columns:
        return ordered_ids[:batch_size]

    cluster_series = df.set_index(id_col)[cluster_col]

    def _cluster_key(cand_id: str) -> object:
        raw = cluster_series.get(cand_id)
        return f"__isolated__{cand_id}" if pd.isna(raw) else raw

    n_clusters = len({_cluster_key(i) for i in ordered_ids}) or 1
    cap = math.ceil(batch_size / n_clusters)

    counts: Dict[object, int] = defaultdict(int)
    selected: List[str] = []
    for cand_id in ordered_ids:
        if len(selected) >= batch_size:
            break
        c = _cluster_key(cand_id)
        if counts[c] < cap:
            selected.append(cand_id)
            counts[c] += 1
    return selected


def _active_order(eligible_order: List[str], per_signal_ranks: Dict[str, Dict[str, float]],
                   batch_size: int) -> List[str]:
    """Shortlist-then-reorder: top tier by standing, sorted by descending
    disagreement within the tier; the remainder keeps its standing order.
    See the module docstring's "mode=active blend" section.
    """
    pos = {cand_id: idx for idx, cand_id in enumerate(eligible_order)}
    tier_n = min(len(eligible_order), max(batch_size * ACTIVE_TIER_MULTIPLIER, batch_size))
    tier = eligible_order[:tier_n]
    remainder = eligible_order[tier_n:]
    tier_sorted = sorted(
        tier, key=lambda i: (-disagreement_score(per_signal_ranks, i), pos[i])
    )
    return tier_sorted + remainder


def select_batch(df: pd.DataFrame, per_signal_ranks: Dict[str, Dict[str, float]],
                  batch_size: int, classifier_col: str = "classification",
                  cluster_col: str = "tandem_cluster_id", mode: str = "active") -> List[str]:
    """An ordered list of <= ``batch_size`` candidate ids for the probe batch.

    (a) Excludes ids whose ``classifier_col`` is a confident non-chemoreceptor
        call (see `_excluded_ids`).
    (b) ``mode="topscore"``: pure rank-aggregation standing (the comparison
        baseline). ``mode="active"``: standing shortlist reordered by
        disagreement (see `_active_order` / the module docstring).
    (c) Diversity: no more than ``ceil(batch_size / n_clusters)`` from any one
        ``cluster_col`` value, in both modes (see `_apply_diversity_cap`).
    """
    excluded = _excluded_ids(df, classifier_col, id_col=_ID_COL)
    full_order = _rank_agg_order(df, id_col=_ID_COL)
    eligible = [i for i in full_order if i not in excluded]

    if mode == "active":
        candidate_order = _active_order(eligible, per_signal_ranks, batch_size)
    elif mode == "topscore":
        candidate_order = eligible
    else:
        raise ValueError(f"unknown mode: {mode!r} (expected 'active' or 'topscore')")

    return _apply_diversity_cap(candidate_order, df, batch_size, cluster_col, id_col=_ID_COL)


def write_batch(selected: List[str], df: pd.DataFrame,
                 per_signal_ranks: Dict[str, Dict[str, float]], path: str,
                 method: str = "rra", groups=None,
                 cluster_col: str = "tandem_cluster_id") -> pd.DataFrame:
    """Write the bench-ready probe-batch TSV; return the written DataFrame.

    Columns: ``id, rank_agg_position, disagreement, cluster,
    top_supporting_signals, conflicting_signals, rationale``.
    ``rank_agg_position`` is the id's 1-based position in the FULL (i.e. not
    classifier-gated) rank-aggregation order over every id in ``df`` -- the
    single, mode-independent "how did the label-free signals rank this
    candidate overall" number. ``top_supporting_signals`` /
    ``conflicting_signals`` are the signal(s) that rank the candidate best /
    worst (see `_signal_extremes`). Creates parent directories if needed.
    """
    full_order = _rank_agg_order(df, id_col=_ID_COL, method=method, groups=groups)
    position = {cand_id: p for p, cand_id in enumerate(full_order, start=1)}
    n_total = len(df)

    cluster_series = df.set_index(_ID_COL)[cluster_col] if cluster_col in df.columns else None

    rows = []
    for cand_id in selected:
        disagreement = disagreement_score(per_signal_ranks, cand_id)
        best, worst = _signal_extremes(per_signal_ranks, cand_id)
        top_sig = ",".join(best)
        conf_sig = ",".join(worst)

        cluster = ""
        if cluster_series is not None and cand_id in cluster_series.index:
            raw_cluster = cluster_series.loc[cand_id]
            cluster = "" if pd.isna(raw_cluster) else raw_cluster

        pos = position.get(cand_id)
        rationale = (
            f"Rank-aggregation position {pos} of {n_total} candidates; "
            f"disagreement={disagreement:.3f} "
            f"(best-supporting: {top_sig or 'none'}, "
            f"most-conflicting: {conf_sig or 'none'}); "
            f"cluster={cluster or 'none'}. Selected for information gain "
            "(cross-signal disagreement + cluster diversity), not as a pure "
            "top-score pick -- see --mode topscore for that baseline."
        )
        rows.append(
            {
                "id": cand_id,
                "rank_agg_position": pos,
                "disagreement": disagreement,
                "cluster": cluster,
                "top_supporting_signals": top_sig,
                "conflicting_signals": conf_sig,
                "rationale": rationale,
            }
        )

    out_df = pd.DataFrame(
        rows,
        columns=[
            "id",
            "rank_agg_position",
            "disagreement",
            "cluster",
            "top_supporting_signals",
            "conflicting_signals",
            "rationale",
        ],
    )
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, sep="\t", index=False)
    return out_df


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0].strip())
    ap.add_argument("--ranked-csv", required=True,
                     help="Ranked candidate CSV (rank_candidates.py output)")
    ap.add_argument("--batch-size", type=int, default=24)
    ap.add_argument("--mode", default="active", choices=["active", "topscore"])
    ap.add_argument("--classifier-col", default="classification")
    ap.add_argument("--cluster-col", default="tandem_cluster_id")
    ap.add_argument("--out", default=None,
                     help="Output TSV (default: $RESULTS_DIR/ranking/"
                          "probe_batch_active_learning.tsv)")
    args = ap.parse_args(argv)

    df = pd.read_csv(args.ranked_csv)
    raw_ranklists = build_ranklists_from_df(df)
    per_signal_ranks = {sig: normalized_ranks(scores) for sig, scores in raw_ranklists.items()}

    selected = select_batch(
        df, per_signal_ranks, args.batch_size,
        classifier_col=args.classifier_col, cluster_col=args.cluster_col, mode=args.mode,
    )

    out_path = args.out or os.path.join(
        os.environ.get("RESULTS_DIR", "results"), "ranking", "probe_batch_active_learning.tsv"
    )
    write_batch(selected, df, per_signal_ranks, out_path, cluster_col=args.cluster_col)
    print(
        f"[select_probe_batch] mode={args.mode}: wrote {len(selected)} candidates "
        f"to {out_path}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
