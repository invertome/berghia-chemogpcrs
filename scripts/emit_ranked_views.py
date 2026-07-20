#!/usr/bin/env python3
"""emit_ranked_views.py — split the ranked candidates CSV into two views.

Bead 1nr. The pipeline emits ONE composite-sorted ranked CSV
(results/ranking/ranked_candidates_sorted.csv). Sorting everything by a
single rank_score conflates two distinct wet-lab priorities, so this script
re-projects the SAME scored rows into two complementary views (it never
re-scores anything):

  CONFIDENCE view  (results/ranking/ranked_candidates_confidence.csv)
      The safe-bet HCR shortlist. Keep rows that are:
        * classification == 'chemoreceptor-candidate'  (not a known
          non-chemoreceptor / not flagged likely-non-chemoreceptor), AND
        * needs_manual_review != 'yes'                 (the automated
          classifiers reached a confident call), AND
        * evidence_completeness >= CONFIDENCE_MIN_COMPLETENESS (default 0.7)
      Sorted by rank_score descending — the production composite order,
      restricted to well-supported chemoreceptor candidates.

  DISCOVERY view   (results/ranking/ranked_candidates_discovery.csv)
      High-novelty divergent lineage-specific-expansion (LSE) chemoreceptors
      that a single completeness-rewarding composite would bury. Keep rows
      that are chemoreceptor-candidate AND trip at least one divergence flag:
        * og_dnds_reliability == 'low'   (reference-poor OG => the paralog is
          divergent enough that few references have recoverable CDS), OR
        * needs_manual_review == 'yes'   (HMM + OG classifiers both went dark
          — the unclassifiable / divergent case), OR
        * the tandem-cluster signal is high
          (tandem_cluster_score_norm >= DISCOVERY_TANDEM_HIGH, default 0.5;
          intra-genome tandem arrays are the field's signature LSE-
          chemoreceptor evidence).
      Sorted by a DISCOVERY SCORE that rewards divergence/novelty rather than
      evidence completeness (see formula below).

DISCOVERY SCORE (documented, exact — cw-lkhu Task 2):

      discovery_score = -log10( rho )

  where rho is the Robust Rank Aggregation (Kolde et al. 2012) score from
  scripts/rank_aggregation.py's rra_score() over the discovery signals
      {tandem, positive, novelty, lse_divergence, lse_nesting_depth}
  read from the ranked CSV as tandem_cluster_score_norm / positive_score_norm /
  emb_novelty / lse_divergence_score_norm / lse_nesting_depth_score_norm (each
  falling back to its raw <signal>_score column). RRA ranks each candidate within every signal it has data for and
  takes the most extreme order statistic across signals, so the score is
  WEIGHT-FREE: no DISCOVERY_*_WEIGHT knob enters, and rank aggregation is
  invariant to each signal's monotone normalization. rho is LOWER = better, so
  -log10(rho) is HIGHER = more consistently novel; the view is sorted by
  discovery_score descending.

  The signal set deliberately excludes phylogeny, synteny, expression and
  evidence_completeness — a divergent LSE paralog is reference-poor and
  evidence-poor by nature, so those axes penalize exactly the candidates this
  view exists to surface. BOTH LSE depth axes are included -- lse_divergence
  (cumulative branch length) and lse_nesting_depth (duplication depth); each is
  a positive divergence signal, they measure different things inside the top
  quartile that shapes the order, and the earlier hand-weighted mean omitted
  both. A candidate votes in a signal only
  where it has data: a missing signal column, or a blank value, simply drops
  that candidate's vote for that signal (RRA over whatever remains); if NO
  discovery signal column is present, discovery_score is 0 for every row.
  Embedding novelty (emb_novelty) enters the discovery SCORE only — never the
  confidence composite, and never as a membership disjunct (the locked cw3.6
  decision).

Robustness: every filter column and every signal column is OPTIONAL. A missing
column is skipped (its filter is dropped / its signal is treated as absent)
with a warning to stderr rather than crashing — the ranked CSV's exact column
set drifts as upstream augmenters are added. An empty result set writes a
header-only CSV (all columns preserved, zero data rows).

Config (read from the environment in main(), with defaults):
    CONFIDENCE_MIN_COMPLETENESS   confidence evidence-completeness floor (0.7)
    DISCOVERY_TANDEM_HIGH         tandem_cluster_score_norm "high" cutoff (0.5)

Usage:
    python3 emit_ranked_views.py \\
        --ranked-csv results/ranking/ranked_candidates_sorted.csv \\
        --confidence-out results/ranking/ranked_candidates_confidence.csv \\
        --discovery-out  results/ranking/ranked_candidates_discovery.csv
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

CONFIDENCE_CLASS = "chemoreceptor-candidate"


def _warn(msg: str) -> None:
    print(f"emit_ranked_views: {msg}", file=sys.stderr)


def _numeric(series: pd.Series) -> pd.Series:
    """Coerce a (possibly string) column to numeric; unparseable -> NaN."""
    return pd.to_numeric(series, errors="coerce")


def _resolve_norm_signal(df: pd.DataFrame, primary: str,
                         fallback: str | None = None):
    """Return (Series, present) for a normalized [0,1] signal.

    Prefer ``primary``; fall back to ``fallback`` ONLY when the fallback is
    guaranteed to be on the same [0,1] scale (e.g. tandem_cluster_score, which
    equals tandem_cluster_score_norm by construction). ``present`` is False
    when neither column exists — the caller then drops the signal entirely.
    """
    if primary in df.columns:
        return _numeric(df[primary]).fillna(0.0), True
    if fallback and fallback in df.columns:
        _warn(f"column {primary!r} absent; using {fallback!r} "
              f"(identical [0,1] value by construction)")
        return _numeric(df[fallback]).fillna(0.0), True
    return pd.Series(0.0, index=df.index), False


def build_confidence_view(df: pd.DataFrame,
                          min_completeness: float = 0.7) -> pd.DataFrame:
    """The safe-bet shortlist: confident chemoreceptor candidates with
    complete evidence, sorted by rank_score descending. Every filter column is
    optional — a missing one is skipped (warned) instead of crashing."""
    mask = pd.Series(True, index=df.index)

    if "classification" in df.columns:
        mask &= df["classification"].astype(str) == CONFIDENCE_CLASS
    else:
        _warn("no 'classification' column; not restricting confidence view "
              "to chemoreceptor-candidate rows")

    if "needs_manual_review" in df.columns:
        mask &= df["needs_manual_review"].astype(str) != "yes"
    else:
        _warn("no 'needs_manual_review' column; not excluding "
              "manual-review rows from confidence view")

    if "evidence_completeness" in df.columns:
        mask &= _numeric(df["evidence_completeness"]) >= float(min_completeness)
    else:
        _warn("no 'evidence_completeness' column; completeness filter skipped")

    out = df[mask].copy()

    # Honor the PRODUCTION order. rank_candidates.py emits final_rank = the
    # 1-based order of the active RANK_METHOD (weighted OR rankagg); sort by it
    # so promoting rankagg actually reorders the safe-bet shortlist. rank_score
    # stays the weighted value even under rankagg, so it is only the fallback for
    # older ranked CSVs that predate the final_rank column.
    if "final_rank" in out.columns:
        order = _numeric(out["final_rank"]).sort_values(
            ascending=True, kind="mergesort").index
        out = out.loc[order]
    elif "rank_score" in out.columns:
        order = _numeric(out["rank_score"]).sort_values(
            ascending=False, kind="mergesort").index
        out = out.loc[order]
    else:
        _warn("no 'final_rank'/'rank_score' column; confidence view left unsorted")

    return out.reset_index(drop=True)


def build_discovery_view(df: pd.DataFrame,
                         tandem_high: float = 0.5) -> pd.DataFrame:
    """The divergent-LSE view: chemoreceptor candidates that trip at least one
    divergence flag (reference-poor OG / manual-review / high tandem signal),
    sorted by a WEIGHT-FREE discovery score = -log10 of the Robust Rank
    Aggregation rho over {tandem, positive, novelty, lse_divergence,
    lse_nesting_depth}. See the module
    header for the exact formula. Every filter and signal column is optional.

    The consensus embedding-novelty signal (``emb_novelty``, cw3.6) enters the
    discovery SCORE only (as one of the RRA voters); it is deliberately NOT a
    membership disjunct and NOT a term in the confidence composite (novelty is
    orthogonal to chemoreceptor-likelihood — the locked cw3.6 routing
    decision)."""
    # Resolve the tandem signal ONCE for the membership filter below (the
    # discovery SCORE reads its own signals independently via RRA).
    tandem, tandem_present = _resolve_norm_signal(
        df, "tandem_cluster_score_norm", "tandem_cluster_score")
    if not tandem_present:
        _warn("no tandem-cluster score column "
              "('tandem_cluster_score_norm'/'tandem_cluster_score'); "
              "tandem membership filter treated as absent")

    # Hard AND gate: only chemoreceptor candidates are eligible.
    if "classification" in df.columns:
        chemo = df["classification"].astype(str) == CONFIDENCE_CLASS
    else:
        _warn("no 'classification' column; not restricting discovery view "
              "to chemoreceptor-candidate rows")
        chemo = pd.Series(True, index=df.index)

    # OR over the divergence disjuncts (each optional).
    disj = pd.Series(False, index=df.index)
    any_disjunct = False
    if "og_dnds_reliability" in df.columns:
        disj |= df["og_dnds_reliability"].astype(str) == "low"
        any_disjunct = True
    else:
        _warn("no 'og_dnds_reliability' column; low-reliability disjunct skipped")
    if "needs_manual_review" in df.columns:
        disj |= df["needs_manual_review"].astype(str) == "yes"
        any_disjunct = True
    else:
        _warn("no 'needs_manual_review' column; manual-review disjunct skipped")
    if tandem_present:
        disj |= tandem >= float(tandem_high)
        any_disjunct = True
    if not any_disjunct:
        _warn("no discovery disjunct columns present; discovery view is empty")

    out = df[chemo & disj].copy()

    # cw-lkhu Task 2: order by a weight-free Robust Rank Aggregation (Kolde et
    # al. 2012) over the discovery signals, replacing the hand-picked weighted
    # mean. RRA ranks each candidate within every signal it has data for and
    # scores the most extreme order statistic across signals; a candidate
    # consistently top-ranked (or extreme on a signal with few competitors)
    # rises. lse_divergence is included (the weighted mean omitted it).
    from rank_aggregation import rra_score
    import numpy as _np
    _signal_cols = {
        "tandem": ("tandem_cluster_score_norm", "tandem_cluster_score"),
        "positive": ("positive_score_norm", "positive_score"),
        "novelty": ("emb_novelty", None),
        "lse_divergence": ("lse_divergence_score_norm", "lse_divergence_score"),
        # Bead hf3u: the TOPOLOGICAL companion to lse_divergence. This view is
        # explicitly about divergent LSE paralogs, and duplication depth (how
        # many nodes deep inside the expansion a candidate sits) is a divergence
        # signal in its own right -- it was simply missing here while the
        # patristic axis was included. The two measure different things: on this
        # repo's real trees they agree population-wide (spearman +0.897/+0.604)
        # but barely inside the top quartile that actually shapes an order
        # (+0.078/+0.034), so listing only one of them silently picked a
        # different scored set. RRA is weight-free, so adding a voter changes
        # no weighting decision; a candidate with no nesting measurement simply
        # does not vote in this signal.
        "lse_nesting_depth": ("lse_nesting_depth_score_norm",
                              "lse_nesting_depth_score"),
    }
    per_signal = {}
    for name, (primary, fallback) in _signal_cols.items():
        col = primary if primary in out.columns else (
            fallback if fallback and fallback in out.columns else None)
        if col is None:
            continue
        series = _numeric(out[col])
        vals = {cid: float(v) for cid, v in zip(out["id"], series) if pd.notna(v)}
        if vals:
            per_signal[name] = vals
    if per_signal:
        rho = rra_score(per_signal)                       # lower = better
        out["discovery_score"] = out["id"].map(
            lambda i: float(-_np.log10(max(rho.get(i, 1.0), 1e-300))))
    else:
        _warn("no discovery signal columns present; discovery_score = 0")
        out["discovery_score"] = 0.0

    out = out.sort_values("discovery_score", ascending=False, kind="mergesort")
    return out.reset_index(drop=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--ranked-csv", required=True,
                    help="Augmented ranked candidates CSV "
                         "(rank_candidates.py output after the column adders)")
    ap.add_argument("--confidence-out", required=True,
                    help="Output CSV for the confidence (safe-bet) view")
    ap.add_argument("--discovery-out", required=True,
                    help="Output CSV for the discovery (divergent-LSE) view")
    args = ap.parse_args(argv)

    min_completeness = float(os.getenv("CONFIDENCE_MIN_COMPLETENESS", "0.7"))
    tandem_high = float(os.getenv("DISCOVERY_TANDEM_HIGH", "0.5"))

    if not os.path.exists(args.ranked_csv):
        _warn(f"ranked CSV not found: {args.ranked_csv}")
        return 1

    df = pd.read_csv(args.ranked_csv, dtype=str, keep_default_na=False)
    # Accept the pre-rename `lse_depth_*` schema under its canonical
    # `lse_divergence_*` names (announced, never silent). Without this a ranked
    # CSV written before the rename would drop the divergence signal out of the
    # discovery RRA entirely -- indistinguishable, at the call site, from a
    # candidate that genuinely has no divergence measurement.
    from _rank_candidates_lib import apply_legacy_column_aliases
    df = apply_legacy_column_aliases(df, source=f"ranked CSV {args.ranked_csv}")

    conf = build_confidence_view(df, min_completeness=min_completeness)
    disc = build_discovery_view(df, tandem_high=tandem_high)

    for out_path, view in ((args.confidence_out, conf),
                           (args.discovery_out, disc)):
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        view.to_csv(out_path, index=False)

    _warn(f"confidence view: {len(conf)} rows -> {args.confidence_out}")
    _warn(f"discovery view:  {len(disc)} rows -> {args.discovery_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
