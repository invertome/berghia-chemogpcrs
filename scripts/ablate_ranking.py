#!/usr/bin/env python3
"""ablate_ranking.py — Ranking ablation harness for rank_candidates.py.

A council review of the hand-weighted chemoreceptor ranking raised two
empirical questions about what today's ranking is actually driven by:

    (a) the ``has_*_data`` / ``evidence_completeness`` missingness machinery
        may be encoding "this gene sits in a detectable tandem array"
        (i.e. how MUCH evidence exists for a candidate) rather than
        biology (what the evidence actually says); and
    (b) a few correlated signal groups (e.g. phylo + og_confidence, both
        derived from the same gene trees) may dominate the composite
        score regardless of the other axes.

This module answers both by re-ranking the same candidate set under
targeted ablations — zeroing a signal group's weight, or neutralizing the
missingness flags so every axis looks "present" — and comparing the result
to the real (unablated) ranking via Spearman rank correlation and top-k
Jaccard overlap. It is a pure analysis module: it takes an already-scored
DataFrame (the *_score_norm / has_*_data columns rank_candidates.py
produces) and re-ranks in memory. It does not read or write any pipeline
files itself.

It reuses the EXACT production scorer that produces the ranking actually
written to ranked_candidates_sorted.csv: ``rank_candidates.py`` computes
its ``rank_score`` column (line ~2002, then sorts by it) via the row-wise
``calculate_fair_rank_score`` (rank_candidates.py:1953), which bridges each
candidate's per-axis normalized scores and the ``*_WEIGHT`` env constants
into ``_rank_candidates_lib.calculate_fair_rank_score`` with
``completeness_floor=0.4`` (bead -ce4). That lib function is explicitly
side-effect-free and importable, so this module imports it directly and
reconstructs the same 12-signal scores/weights bridge (INCLUDING the
tandem_cluster axis and the evidence-completeness multiplier). The earlier
``calculate_rank_score`` (rank_candidates.py:1382) is NOT the ranking
scorer — it is used only by the Monte-Carlo sensitivity analysis and lacks
both the completeness multiplier and the tandem_cluster signal, so it is
deliberately not used here.

Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts
"""
from __future__ import annotations

import os
import sys
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Union

import pandas as pd
from scipy import stats

# Sibling-import the pure scoring lib the production ranker uses. Matches the
# import-path setup rank_candidates.py and the other scripts in this dir do
# (this only mutates sys.path so the sibling import resolves — no data
# loading, file I/O, or env reads happen at import time).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _rank_candidates_lib import calculate_fair_rank_score  # noqa: E402

# The 12 signal axes of the production fair scorer, mapping each weight key
# to the DataFrame columns rank_candidates.py's calculate_fair_rank_score
# reads (rank_candidates.py:1970-1998). ``gated`` axes contribute their
# score only when the paired ``has_*_data`` flag is truthy (otherwise the
# score is None = "missing", which the fair scorer excludes from the
# available-weight numerator but still counts in the total-weight
# denominator — that is the missingness penalty). The four base axes
# (phylo/purifying/positive/lse_divergence) are always contributed. NOTE the
# 'expression' weight key is fed by the EXPR_WEIGHT env var, and
# tandem_cluster (the field's signature chemoreceptor signal) is the 12th.
_SIGNAL_SPEC = [
    # (weight_key, score_norm_column, has_flag_column_or_None, env_var, default)
    ("phylo", "phylo_score_norm", None, "PHYLO_WEIGHT", 2.0),
    # Gated on has_dnds_data, mirroring rank_candidates.py's production scorer.
    # These were ungated, so a candidate aBSREL never reported on contributed a
    # full-weight present 0.0 to the ablation. That matters most under the
    # orthology quarantine (ORTHOLOGY_SOURCE_TRUSTED=0), where the flag is False
    # for every candidate: ungated, the ablation would keep "measuring" the
    # contribution of an axis that is not voting in production at all.
    ("purifying", "purifying_score_norm", "has_dnds_data", "PURIFYING_WEIGHT", 1.0),
    ("positive", "positive_score_norm", "has_dnds_data", "POSITIVE_WEIGHT", 1.0),
    ("lse_divergence", "lse_divergence_score_norm", None, "LSE_DIVERGENCE_WEIGHT", 1.0),
    # Bead hf3u: the topological companion to lse_divergence (node count vs
    # cumulative branch length). Unlike the four base axes it is GATED on its
    # own has_*_data flag, matching how rank_candidates.py's production scorer
    # contributes it -- the axis must drop out where it could not be measured
    # rather than contribute a full-weight present-zero.
    ("lse_nesting_depth", "lse_nesting_depth_score_norm",
     "has_lse_nesting_depth_data", "LSE_NESTING_DEPTH_WEIGHT", 1.0),
    ("synteny", "synteny_score_norm", "has_synteny_data", "SYNTENY_WEIGHT", 3.0),
    ("expression", "expression_score_norm", "has_expression_data", "EXPR_WEIGHT", 1.0),
    ("chemosensory_expr", "chemosensory_expr_score_norm", "has_chemosensory_expr_data",
     "CHEMOSENSORY_EXPR_WEIGHT", 3.0),
    ("gprotein_coexpr", "gprotein_coexpr_score_norm", "has_gprotein_data",
     "GPROTEIN_COEXPR_WEIGHT", 2.0),
    ("ecl_divergence", "ecl_divergence_score_norm", "has_ecl_data",
     "ECL_DIVERGENCE_WEIGHT", 1.5),
    ("expansion", "expansion_score_norm", "has_expansion_data", "EXPANSION_WEIGHT", 1.5),
    ("og_confidence", "og_confidence_score_norm", "has_og_confidence_data",
     "OG_CONFIDENCE_WEIGHT", 1.0),
    ("tandem_cluster", "tandem_cluster_score_norm", "has_tandem_cluster_data",
     "TANDEM_CLUSTER_WEIGHT", 2.5),
]

_COMPLETENESS_FLOOR = 0.4  # matches rank_candidates.py:1999


# ---------------------------------------------------------------------------
# Top-k Jaccard
# ---------------------------------------------------------------------------


def topk_jaccard(rank_a: Sequence[Any], rank_b: Sequence[Any], k: int) -> float:
    """Jaccard index of the top-k sets of two ranked id lists.

    ``|top-k(a) ∩ top-k(b)| / |top-k(a) ∪ top-k(b)|``. If both top-k sets
    are empty (e.g. k=0, or both lists empty), the union is empty and this
    returns 1.0 (vacuously identical) rather than dividing by zero.
    """
    top_a = set(rank_a[:k])
    top_b = set(rank_b[:k])
    union = top_a | top_b
    if not union:
        return 1.0
    return len(top_a & top_b) / len(union)


# ---------------------------------------------------------------------------
# Production weights + per-row scoring bridge
# ---------------------------------------------------------------------------


def _production_weights() -> Dict[str, float]:
    """Construct the 12-signal weights dict exactly as rank_candidates.py's
    ``calculate_fair_rank_score`` does (rank_candidates.py:1985-1998): the
    module-level ``*_WEIGHT`` env vars (same names, same string defaults)
    mapped to the fair scorer's weight keys. Note the key is 'expression'
    (fed by EXPR_WEIGHT), not 'expr', and the tandem_cluster axis is
    included with default 2.5.
    """
    # getenv_renamed honours the pre-rename LSE_DEPTH_WEIGHT (announced) for
    # keys that were renamed, and behaves exactly like os.getenv for the rest.
    from _rank_candidates_lib import getenv_renamed
    return {key: getenv_renamed(env, default, cast=float)
            for key, _, _, env, default in _SIGNAL_SPEC}


def _row_scores(row: Mapping[str, Any]) -> Dict[str, Optional[float]]:
    """Build the per-candidate scores dict exactly as rank_candidates.py's
    ``calculate_fair_rank_score`` does (rank_candidates.py:1970-1984): base
    axes always contribute their normalized score; gated axes contribute
    their score only when the paired ``has_*_data`` flag is truthy,
    otherwise None (missing).
    """
    scores: Dict[str, Optional[float]] = {}
    for key, score_col, has_col, _, _ in _SIGNAL_SPEC:
        if has_col is None:
            scores[key] = row.get(score_col)
        else:
            scores[key] = row.get(score_col) if row.get(has_col) else None
    return scores


def _fair_scores(df: pd.DataFrame, weights: Mapping[str, float]) -> pd.Series:
    """Compute the production fair rank_score for every row of df under the
    given weights. Mirrors ``df.apply(calculate_fair_rank_score, axis=1)``
    from rank_candidates.py:2002, feeding the same lib function.
    """
    return df.apply(
        lambda row: calculate_fair_rank_score(
            _row_scores(row), weights, completeness_floor=_COMPLETENESS_FLOOR
        ),
        axis=1,
    )


# ---------------------------------------------------------------------------
# Ablation primitives
# ---------------------------------------------------------------------------


def _neutralize_missingness(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of df with every has_*_data column forced True and
    evidence_completeness forced to 1.0 (when present). This makes every
    gated axis contribute (through the fair scorer's completeness
    multiplier), isolating whether the ranking is driven by which signals
    a candidate happens to have data for rather than by what the signals
    themselves say (hypothesis a in the module docstring).
    """
    out = df.copy(deep=True)
    for col in out.columns:
        if col.startswith("has_"):
            out[col] = True
    if "evidence_completeness" in out.columns:
        out["evidence_completeness"] = 1.0
    return out


def _apply_drop_signals(
    weights: Mapping[str, float], drop_signals: Optional[Iterable[str]]
) -> Dict[str, float]:
    """Return a copy of weights with each key in drop_signals zeroed out.

    In the fair scorer a zeroed weight removes that axis from BOTH the
    available-weight numerator and the total-weight denominator, so the
    signal group is fully ablated (as if it did not exist).
    """
    w = dict(weights)
    if drop_signals:
        for key in drop_signals:
            w[key] = 0.0
    return w


def ablate(
    df: pd.DataFrame,
    drop_signals: Optional[Iterable[str]] = None,
    neutralize_missingness: bool = False,
    weights: Optional[Mapping[str, float]] = None,
) -> List[Any]:
    """Re-rank candidates under an ablation and return the re-ranked id list.

    Reuses the real production fair scorer (``_rank_candidates_lib.
    calculate_fair_rank_score`` via `_fair_scores`) so the result reflects
    the actual ranking math — completeness multiplier and all 12 signals
    included — not an approximation.

    Args:
        df: DataFrame with the *_score_norm / has_*_data columns
            rank_candidates.py's calculate_fair_rank_score consumes, plus
            an 'id' column. Not mutated.
        drop_signals: weight keys to zero out (e.g. ``['phylo', 'synteny']``;
            any of the 12 keys in `_production_weights`, incl.
            'tandem_cluster').
        neutralize_missingness: if True, force every has_*_data column to
            True and evidence_completeness to 1.0 before scoring, so the
            completeness penalty no longer distinguishes candidates.
        weights: full weights dict to use instead of the production
            defaults (`_production_weights()`). `drop_signals` is applied
            on top of whichever weights dict is in effect.

    Returns:
        List of ids from df['id'], sorted by the re-computed fair score in
        descending order (stable sort — ties keep their original relative
        order).
    """
    working = _neutralize_missingness(df) if neutralize_missingness else df.copy(deep=True)
    w = _apply_drop_signals(
        weights if weights is not None else _production_weights(), drop_signals
    )

    scores = _fair_scores(working, w)
    order = scores.sort_values(ascending=False, kind="mergesort").index
    return working.loc[order, "id"].tolist()


# ---------------------------------------------------------------------------
# Ablation report
# ---------------------------------------------------------------------------

GroupsType = Union[Mapping[str, Sequence[str]], Iterable[Sequence[str]]]


def _named_groups(groups: GroupsType) -> List[Any]:
    """Normalize `groups` into a list of (name, [signal keys]) pairs.

    Accepts either a mapping of group-name -> signal-key list, or a plain
    iterable of signal-key lists (name auto-derived by joining the keys
    with '+', e.g. ``["phylo", "og_confidence"]`` -> "phylo+og_confidence").
    """
    if isinstance(groups, Mapping):
        return list(groups.items())
    return [("+".join(g), list(g)) for g in groups]


def _compare_to_baseline(
    df: pd.DataFrame,
    baseline_rank: Sequence[Any],
    k: int,
    drop_signals: Optional[Iterable[str]] = None,
    neutralize_missingness: bool = False,
) -> Dict[str, float]:
    """Ablate df and compare the result to an already-computed baseline
    ranking. Shared by every entry `ablation_report` produces."""
    ablated_rank = ablate(
        df, drop_signals=drop_signals, neutralize_missingness=neutralize_missingness
    )

    baseline_pos = {cid: i for i, cid in enumerate(baseline_rank)}
    ablated_pos = {cid: i for i, cid in enumerate(ablated_rank)}
    common = [cid for cid in baseline_rank if cid in ablated_pos]
    if len(common) >= 2:
        rho, _ = stats.spearmanr(
            [baseline_pos[c] for c in common],
            [ablated_pos[c] for c in common],
        )
    else:
        rho = float("nan")

    baseline_topk = set(baseline_rank[:k])
    ablated_topk = set(ablated_rank[:k])

    return {
        "spearman": float(rho),
        "jaccard_at_k": topk_jaccard(baseline_rank, ablated_rank, k),
        "n_moved_into_topk": len(ablated_topk - baseline_topk),
    }


def ablation_report(
    df: pd.DataFrame, groups: GroupsType, k: int = 20
) -> Dict[str, Dict[str, float]]:
    """Compare the real ranking to a missingness-neutralized ranking, and to
    one ranking per signal group with that group's weights zeroed out.

    Args:
        df: same contract as `ablate`.
        groups: either a mapping of group-name -> list of weight keys to
            drop together (e.g. ``{"phylo_axis": ["phylo", "og_confidence"]}``)
            or an iterable of such lists (group name auto-derived by
            joining the keys with '+', e.g. "phylo+og_confidence").
        k: top-k size for the Jaccard / moved-into-topk metrics.

    Returns:
        dict with key "missingness" plus one "group:<name>" key per group,
        each mapping to ``{"spearman", "jaccard_at_k", "n_moved_into_topk"}``
        versus the real (unablated) ranking.
    """
    baseline_rank = ablate(df)

    report: Dict[str, Dict[str, float]] = {
        "missingness": _compare_to_baseline(df, baseline_rank, k, neutralize_missingness=True)
    }
    for name, signals in _named_groups(groups):
        report[f"group:{name}"] = _compare_to_baseline(df, baseline_rank, k, drop_signals=signals)
    return report


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------


def write_markdown(report: Mapping[str, Mapping[str, float]], path: Union[str, os.PathLike]) -> None:
    """Write `ablation_report`'s output as a readable markdown table."""
    lines = [
        "# Ranking ablation report",
        "",
        "Top-k churn of the hand-weighted chemoreceptor ranking under",
        "missingness neutralization and per-signal-group ablation, versus",
        "the real (unablated) ranking.",
        "",
        "| Ablation | Spearman rho | Jaccard@k | N moved into top-k |",
        "|---|---:|---:|---:|",
    ]
    for name, metrics in report.items():
        spearman = metrics.get("spearman", float("nan"))
        jaccard = metrics.get("jaccard_at_k", float("nan"))
        n_moved = metrics.get("n_moved_into_topk", 0)
        lines.append(f"| {name} | {spearman:.3f} | {jaccard:.3f} | {n_moved} |")
    lines.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
