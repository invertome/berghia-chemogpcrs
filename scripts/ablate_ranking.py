#!/usr/bin/env python3
"""ablate_ranking.py — Ranking ablation harness for rank_candidates.py.

A council review of the hand-weighted chemoreceptor ranking raised two
empirical questions about what today's ranking is actually driven by:

    (a) the `has_*_data` / `evidence_completeness` missingness machinery
        may be encoding "this gene sits in a detectable tandem array"
        (i.e. how MUCH evidence exists for a candidate) rather than
        biology (what the evidence actually says); and
    (b) a few correlated signal groups (e.g. phylo + og_confidence, both
        derived from the same gene trees) may dominate the composite
        score regardless of the other axes.

This module answers both by re-ranking the same candidate set under
targeted ablations — dropping a signal group's weight to zero, or
neutralizing the missingness flags so every axis looks "present" — and
comparing the result to the real (unablated) ranking via Spearman rank
correlation and top-k Jaccard overlap. It is a pure analysis module: it
takes an already-scored DataFrame (the *_score_norm / has_*_data columns
rank_candidates.py produces) and re-ranks in memory. It does not read or
write any pipeline files itself.

Reuses the REAL production scorer, ``rank_candidates.calculate_rank_score``,
so ablation results reflect the actual ranking math rather than an
approximation of it — see `_resolve_calculate_rank_score` for how that
function is recovered. rank_candidates.py is a one-shot CLI script (it
parses ``sys.argv[1:7]`` and loads pipeline files at module scope, all
before ``calculate_rank_score`` is even defined), so a plain
``from rank_candidates import calculate_rank_score`` reliably fails
under a test process — and, worse, if it ever got far enough with an
unrelated argv it could reach the script's own file-writing code. Instead
this module recovers the function by slicing it out of rank_candidates.py's
source text and exec'ing just that fragment: the same technique already
used by tests/unit/test_rank_candidates_dnds_reliability.py and
tests/unit/test_ranking_lib.py for the same reason. That gives the ACTUAL
current function body (no hand-copy drift risk) without executing any of
the script's argv-parsing / file-loading / file-writing side effects. A
hand-maintained mirror (`_local_calculate_rank_score`) is the last-resort
fallback if even that extraction fails.

Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts
"""
from __future__ import annotations

import os
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Union

import pandas as pd
from scipy import stats

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
# Production scorer resolution (see module docstring)
# ---------------------------------------------------------------------------

_PROD_SCORER: Optional[Any] = None
_PROD_RESOLVED = False


def _extract_calculate_rank_score_from_source() -> Any:
    """Recover the real `calculate_rank_score` by slicing it out of
    scripts/rank_candidates.py's source and exec'ing just that fragment.

    rank_candidates.py cannot be imported directly — it parses
    ``sys.argv[1:7]`` and loads pipeline files at module scope, long
    before ``calculate_rank_score`` is even defined, so a plain import
    raises during that module-level execution (and could, with a
    different argv, run all the way to the script's own CSV-writing
    code). This is the same technique
    tests/unit/test_rank_candidates_dnds_reliability.py and
    tests/unit/test_ranking_lib.py already use for the same script — it
    recovers the ACTUAL current function body (no hand-copy drift risk)
    without running any of the script's side effects.
    """
    src_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "rank_candidates.py")
    with open(src_path) as fh:
        src = fh.read()
    start = src.find("def calculate_rank_score(")
    if start == -1:
        raise RuntimeError("calculate_rank_score not found in rank_candidates.py")
    end = src.find("\n\ndef ", start + 1)
    fragment = src[start:end] if end != -1 else src[start:]
    namespace: Dict[str, Any] = {"pd": pd, "np": __import__("numpy")}
    exec(fragment, namespace)
    return namespace["calculate_rank_score"]


def _local_calculate_rank_score(df: pd.DataFrame, weights: Mapping[str, float]) -> pd.Series:
    """Hand-maintained mirror of rank_candidates.calculate_rank_score, used
    only if the source-extraction above fails (e.g. the source file
    becomes unreadable or the function is renamed/restructured).

    Mirrors the exact per-row logic as of the 2026-05 review-fix sweep:
    phylo / purifying / positive / lse_depth are always included (the
    latter two scaled by the per-OG ``dnds_reliability_weight``, default
    1.0 when absent); synteny / expr / chemosensory_expr /
    gprotein_coexpr / ecl_divergence / expansion / og_confidence are
    included only when their ``has_*_data`` flag is True. The composite
    is ``(available_score / available_weight) * sum(all weights)`` — a
    candidate scores as if its available axes were the WHOLE signal,
    then is rescaled to the full-weight scale (this asymmetry is exactly
    what this harness's ``neutralize_missingness`` option probes). Keep
    this in sync with rank_candidates.calculate_rank_score if that
    function changes.
    """
    max_possible_weight = sum(weights.values())

    def calc_row(row):
        dnds_rw = row.get("dnds_reliability_weight", 1.0)
        try:
            dnds_rw = float(dnds_rw)
        except (TypeError, ValueError):
            dnds_rw = 1.0

        score = (
            row["phylo_score_norm"] * weights.get("phylo", 0)
            + row["purifying_score_norm"] * weights.get("purifying", 0) * dnds_rw
            + row["positive_score_norm"] * weights.get("positive", 0) * dnds_rw
            + row["lse_depth_score_norm"] * weights.get("lse_depth", 0)
        )
        total_weight = (
            weights.get("phylo", 0)
            + weights.get("purifying", 0) * dnds_rw
            + weights.get("positive", 0) * dnds_rw
            + weights.get("lse_depth", 0)
        )

        if row.get("has_synteny_data", False):
            score += row["synteny_score_norm"] * weights.get("synteny", 0)
            total_weight += weights.get("synteny", 0)
        if row.get("has_expression_data", False):
            score += row["expression_score_norm"] * weights.get("expr", 0)
            total_weight += weights.get("expr", 0)
        if row.get("has_chemosensory_expr_data", False):
            score += row.get("chemosensory_expr_score_norm", 0) * weights.get("chemosensory_expr", 0)
            total_weight += weights.get("chemosensory_expr", 0)
        if row.get("has_gprotein_data", False):
            score += row.get("gprotein_coexpr_score_norm", 0) * weights.get("gprotein_coexpr", 0)
            total_weight += weights.get("gprotein_coexpr", 0)
        if row.get("has_ecl_data", False):
            score += row.get("ecl_divergence_score_norm", 0) * weights.get("ecl_divergence", 0)
            total_weight += weights.get("ecl_divergence", 0)
        if row.get("has_expansion_data", False):
            score += row.get("expansion_score_norm", 0) * weights.get("expansion", 0)
            total_weight += weights.get("expansion", 0)
        if row.get("has_og_confidence_data", False):
            score += row.get("og_confidence_score_norm", 0) * weights.get("og_confidence", 0)
            total_weight += weights.get("og_confidence", 0)

        if total_weight > 0:
            return (score / total_weight) * max_possible_weight
        return 0.0

    return df.apply(calc_row, axis=1)


def _resolve_calculate_rank_score() -> Any:
    """Return the production ``calculate_rank_score``, preferring the real
    function recovered from rank_candidates.py's source (see module
    docstring). Falls back to `_local_calculate_rank_score` only if that
    extraction fails. Memoized after the first call.
    """
    global _PROD_SCORER, _PROD_RESOLVED
    if _PROD_RESOLVED:
        return _PROD_SCORER
    _PROD_RESOLVED = True

    try:
        _PROD_SCORER = _extract_calculate_rank_score_from_source()
    except Exception:
        _PROD_SCORER = _local_calculate_rank_score
    return _PROD_SCORER


def _production_weights() -> Dict[str, float]:
    """Construct the weights dict exactly as rank_candidates.py does before
    calling calculate_rank_score: the module-level ``*_WEIGHT`` env vars
    (same names, same string defaults) mapped to the dict keys calc_row
    actually reads via ``weights.get(...)`` — note the key is 'expr', not
    'expression' (that's what calculate_rank_score looks up; the
    'expression' key used when *building* rank_candidates.py's
    sensitivity-analysis base_weights dict is a separate, pre-existing
    naming mismatch in that call site, not something to replicate here).
    """
    return {
        "phylo": float(os.getenv("PHYLO_WEIGHT", 2)),
        "purifying": float(os.getenv("PURIFYING_WEIGHT", 1)),
        "positive": float(os.getenv("POSITIVE_WEIGHT", 1)),
        "synteny": float(os.getenv("SYNTENY_WEIGHT", 3)),
        "expr": float(os.getenv("EXPR_WEIGHT", 1)),
        "lse_depth": float(os.getenv("LSE_DEPTH_WEIGHT", 1)),
        "chemosensory_expr": float(os.getenv("CHEMOSENSORY_EXPR_WEIGHT", 3)),
        "gprotein_coexpr": float(os.getenv("GPROTEIN_COEXPR_WEIGHT", 2)),
        "ecl_divergence": float(os.getenv("ECL_DIVERGENCE_WEIGHT", 1.5)),
        "expansion": float(os.getenv("EXPANSION_WEIGHT", 1.5)),
        "og_confidence": float(os.getenv("OG_CONFIDENCE_WEIGHT", 1)),
    }


# ---------------------------------------------------------------------------
# Ablation primitives
# ---------------------------------------------------------------------------


def _neutralize_missingness(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of df with every has_*_data column forced True and
    evidence_completeness forced to 1.0 (when present). Isolates whether
    the ranking is driven by which signals a candidate happens to have
    data for, rather than by what the signals themselves say (hypothesis a
    in the module docstring).
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
    """Return a copy of weights with each key in drop_signals zeroed out."""
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

    Reuses the real production scorer (see `_resolve_calculate_rank_score`)
    so the result reflects the actual ranking math, not an approximation.

    Args:
        df: DataFrame with the *_score_norm / has_*_data columns
            rank_candidates.calculate_rank_score consumes, plus an 'id'
            column. Not mutated.
        drop_signals: weight keys to zero out (e.g. ``['phylo', 'synteny']``).
        neutralize_missingness: if True, force every has_*_data column to
            True and evidence_completeness to 1.0 before scoring.
        weights: full weights dict to use instead of the production
            defaults (`_production_weights()`). `drop_signals` is applied
            on top of whichever weights dict is in effect.

    Returns:
        List of ids from df['id'], sorted by the re-computed score in
        descending order (stable sort — ties keep their original relative
        order).
    """
    working = _neutralize_missingness(df) if neutralize_missingness else df.copy(deep=True)
    w = _apply_drop_signals(
        weights if weights is not None else _production_weights(), drop_signals
    )

    scorer = _resolve_calculate_rank_score()
    scores = scorer(working, w)
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
