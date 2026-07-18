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

DISCOVERY SCORE (documented, exact):

      discovery_score = ( Wt * tandem_cluster_score_norm
                        +  Wp * positive_score_norm
                        +  Wn * emb_novelty_minmax )
                        / ( Wt + Wp + Wn )     # over signals actually present

  where Wt = DISCOVERY_TANDEM_WEIGHT (default 2.0),
        Wp = DISCOVERY_POSITIVE_WEIGHT (default 1.0), and
        Wn = DISCOVERY_NOVELTY_WEIGHT (default 1.0).

  tandem_cluster_score_norm and positive_score_norm are the pipeline's
  already-[0,1]-normalized signal columns (rank_candidates.py:
  tandem_cluster_score_norm is log1p(cluster size) rescaled to [0,1];
  positive_score_norm is the min-max-normalized aBSREL positive-selection
  score). emb_novelty is the consensus embedding-novelty channel (cw3.6,
  unbounded -log10 RRA p) min-max-normalized to [0,1] here. The score is their
  weight-normalized mean, so it stays in [0,1]. It deliberately excludes
  phylogeny, synteny, expression and evidence_completeness — a divergent LSE
  paralog is reference-poor and evidence-poor by nature, so those axes penalize
  exactly the candidates this view exists to surface. Tandem is up-weighted over
  positive selection and embedding novelty by default because an intra-genome
  tandem array is the strongest lineage-specific-expansion signal (it mirrors
  the composite's TANDEM_CLUSTER_WEIGHT 2.5 vs POSITIVE_WEIGHT 1.0 emphasis).
  Embedding novelty enters the discovery SCORE only (never the confidence
  composite, and not as a membership disjunct) — the locked cw3.6 decision.

  If a *_norm signal column is absent, that term is dropped from BOTH the
  numerator and the denominator (weight-normalized mean over whatever is
  present); if neither is present, discovery_score is 0 for every row. Because
  tandem_cluster_score_norm equals tandem_cluster_score by construction in
  rank_candidates.py, the tandem signal falls back to tandem_cluster_score
  when the _norm column is missing (identical value, different name).
  positive_score has no [0,1] fallback (it is not normalized), so it is simply
  dropped when positive_score_norm is absent.

Robustness: every filter column and every signal column is OPTIONAL. A missing
column is skipped (its filter is dropped / its signal is treated as absent)
with a warning to stderr rather than crashing — the ranked CSV's exact column
set drifts as upstream augmenters are added. An empty result set writes a
header-only CSV (all columns preserved, zero data rows).

Config (read from the environment in main(), with defaults):
    CONFIDENCE_MIN_COMPLETENESS   confidence evidence-completeness floor (0.7)
    DISCOVERY_TANDEM_HIGH         tandem_cluster_score_norm "high" cutoff (0.5)
    DISCOVERY_TANDEM_WEIGHT       discovery-score tandem weight Wt (2.0)
    DISCOVERY_POSITIVE_WEIGHT     discovery-score positive-selection weight Wp (1.0)
    DISCOVERY_NOVELTY_WEIGHT      discovery-score embedding-novelty weight Wn (1.0)

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


def _minmax_signal(df: pd.DataFrame, col: str):
    """Return (Series, present) min-max-normalizing an UNBOUNDED signal to [0,1].

    Unlike ``_resolve_norm_signal`` (which expects an already-[0,1] column), the
    consensus embedding novelty ``emb_novelty`` is unbounded (``-log10`` RRA p),
    so it is rescaled to [0,1] here. The range is taken over candidates that
    actually HAVE the signal (non-null), so a no-emb-data candidate (blank) does
    not define the floor/ceiling; blank rows map to 0 (no novelty evidence).
    A missing column, or one that is entirely blank, returns present=False so the
    caller drops the term from both numerator and denominator; a constant column
    normalizes to 0 for every row (present, contributes nothing).
    """
    if col not in df.columns:
        return pd.Series(0.0, index=df.index), False
    v = _numeric(df[col])
    valid = v.dropna()
    if len(valid) == 0:                      # present but entirely blank -> no signal
        return pd.Series(0.0, index=df.index), False
    lo, hi = float(valid.min()), float(valid.max())
    if hi > lo:
        return ((v - lo) / (hi - lo)).fillna(0.0), True
    return pd.Series(0.0, index=df.index), True


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
                         tandem_high: float = 0.5,
                         tandem_weight: float = 2.0,
                         positive_weight: float = 1.0,
                         novelty_weight: float = 1.0) -> pd.DataFrame:
    """The divergent-LSE view: chemoreceptor candidates that trip at least one
    divergence flag (reference-poor OG / manual-review / high tandem signal),
    sorted by a divergence-rewarding discovery score. See the module header for
    the exact score formula. Every filter and signal column is optional.

    The consensus embedding-novelty signal (``emb_novelty``, cw3.6) enters the
    discovery SCORE as a min-max-normalized [0,1] term weighted by
    ``novelty_weight``; it is deliberately NOT a membership disjunct and NOT a
    term in the confidence composite (novelty is orthogonal to chemoreceptor-
    likelihood — the locked cw3.6 routing decision)."""
    # Resolve the normalized signals ONCE (shared by the tandem filter and the
    # discovery score) so we warn at most once per missing column.
    tandem, tandem_present = _resolve_norm_signal(
        df, "tandem_cluster_score_norm", "tandem_cluster_score")
    if not tandem_present:
        _warn("no tandem-cluster score column "
              "('tandem_cluster_score_norm'/'tandem_cluster_score'); "
              "tandem filter + signal treated as absent")
    positive, positive_present = _resolve_norm_signal(df, "positive_score_norm")
    if not positive_present:
        _warn("no 'positive_score_norm' column; "
              "positive-selection signal dropped from discovery score")
    novelty, novelty_present = _minmax_signal(df, "emb_novelty")
    if not novelty_present:
        _warn("no 'emb_novelty' column; consensus embedding-novelty signal "
              "dropped from discovery score")

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

    # discovery_score = weight-normalized mean of the signals actually present.
    denom = (tandem_weight if tandem_present else 0.0) \
        + (positive_weight if positive_present else 0.0) \
        + (novelty_weight if novelty_present else 0.0)
    numer = pd.Series(0.0, index=out.index)
    if tandem_present:
        numer = numer + float(tandem_weight) * tandem.loc[out.index]
    if positive_present:
        numer = numer + float(positive_weight) * positive.loc[out.index]
    if novelty_present:
        numer = numer + float(novelty_weight) * novelty.loc[out.index]
    if denom > 0:
        out["discovery_score"] = numer / denom
    else:
        _warn("no discovery-score signal columns present; discovery_score = 0")
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
    tandem_weight = float(os.getenv("DISCOVERY_TANDEM_WEIGHT", "2.0"))
    positive_weight = float(os.getenv("DISCOVERY_POSITIVE_WEIGHT", "1.0"))
    novelty_weight = float(os.getenv("DISCOVERY_NOVELTY_WEIGHT", "1.0"))

    if not os.path.exists(args.ranked_csv):
        _warn(f"ranked CSV not found: {args.ranked_csv}")
        return 1

    df = pd.read_csv(args.ranked_csv, dtype=str, keep_default_na=False)

    conf = build_confidence_view(df, min_completeness=min_completeness)
    disc = build_discovery_view(df, tandem_high=tandem_high,
                                tandem_weight=tandem_weight,
                                positive_weight=positive_weight,
                                novelty_weight=novelty_weight)

    for out_path, view in ((args.confidence_out, conf),
                           (args.discovery_out, disc)):
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        view.to_csv(out_path, index=False)

    _warn(f"confidence view: {len(conf)} rows -> {args.confidence_out}")
    _warn(f"discovery view:  {len(disc)} rows -> {args.discovery_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
