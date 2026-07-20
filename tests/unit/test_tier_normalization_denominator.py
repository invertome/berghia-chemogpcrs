"""Confidence-tier normalization denominator in rank_candidates.py.

Two coupled defects, one root cause (a hand-maintained second copy of the
weight list drifting away from the scoring path):

1. SCALE. Bead -ce4 (commit db30962) changed the composite scorer from
   ``(weighted_sum / avail_weight) * max_possible_weight`` (weight-sum scale)
   to ``(weighted_sum / avail_weight) * completeness`` (a weighted MEAN in
   [0, 1]). ``assign_confidence_tier`` kept normalizing by a SUM OF WEIGHTS
   (17.0 at production defaults), which was the correct denominator only for
   the OLD numerator. The result is a ~17x scale mismatch: a theoretically
   perfect candidate scores score_pct = 1/17 = 0.059, so the ``High`` gate
   (score_pct >= 0.5) is structurally unreachable and the ``Medium`` score
   clause (>= 0.3) is dead code. The tier column silently degenerates into a
   pure evidence-count threshold.

2. DRIFT. TANDEM_CLUSTER_WEIGHT is in the scoring ``weights`` dict (so it
   reaches the numerator) but was never added to the hand-written denominator
   sum. This is the same class of bug that will recur for the next axis added.

The fix derives the denominator from a single source of truth (SCORING_WEIGHTS)
by asking the scorer itself what an all-axes-perfect candidate scores, so the
denominator tracks BOTH the weight list and the scorer's scale automatically.

rank_candidates.py is not import-safe (top-level ``sys.argv``), so the
structural guards here read the source, following the convention in
tests/unit/test_rank_candidates_config_defaults.py.
"""
from __future__ import annotations

import ast
import os
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RANK = PROJECT_ROOT / "scripts" / "rank_candidates.py"
SCRIPTS = PROJECT_ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import _rank_candidates_lib as lib  # noqa: E402

# The wrapper pins these; the denominator derivation must use the same values.
COMPLETENESS_FLOOR = 0.4
DNDS_RELIABILITY = 1.0


# --------------------------------------------------------------------------
# source-reading helpers
# --------------------------------------------------------------------------

def _module_tree() -> ast.Module:
    return ast.parse(RANK.read_text())


def _find_function(name: str) -> ast.FunctionDef:
    for node in ast.walk(_module_tree()):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise AssertionError(f"function {name}() not found in rank_candidates.py")


def _dict_literal_keys(fn: ast.FunctionDef, var: str) -> set[str] | None:
    """Keys of a ``var = {...}`` dict literal assigned inside ``fn``."""
    for node in ast.walk(fn):
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Dict):
            for tgt in node.targets:
                if isinstance(tgt, ast.Name) and tgt.id == var:
                    return {
                        k.value for k in node.value.keys
                        if isinstance(k, ast.Constant)
                    }
    return None


def _module_namespace() -> dict:
    """Execute only the module-level constant definitions we need.

    Binds every ``*_WEIGHT`` constant at its production getenv default plus,
    once the fix lands, SCORING_WEIGHTS / MAX_POSSIBLE_RANK_SCORE. Everything
    else in the module (argv parsing, I/O) is skipped.
    """
    ns: dict = {
        "os": os,
        "_calculate_fair_rank_score_corrected": lib.calculate_fair_rank_score,
    }
    wanted = {"SCORING_WEIGHTS", "MAX_POSSIBLE_RANK_SCORE", "COMPLETENESS_FLOOR"}
    for node in _module_tree().body:
        if not isinstance(node, ast.Assign):
            continue
        tgt = node.targets[0]
        if not isinstance(tgt, ast.Name):
            continue
        if not (tgt.id.endswith("_WEIGHT") or tgt.id in wanted):
            continue
        mod = ast.Module(body=[node], type_ignores=[])
        try:
            exec(compile(mod, "<rank_candidates_const>", "exec"), ns)
        except Exception:  # pragma: no cover - dict/str weights we don't need
            pass
    return ns


def _scoring_weight_keys() -> set[str]:
    """The axis names that can reach the rank_score numerator."""
    wrapper = _find_function("calculate_fair_rank_score")
    keys = _dict_literal_keys(wrapper, "weights")
    if keys is not None:
        return keys  # pre-fix: literal dict inside the wrapper
    ns = _module_namespace()
    assert "SCORING_WEIGHTS" in ns, (
        "neither a literal weights dict in calculate_fair_rank_score() nor a "
        "module-level SCORING_WEIGHTS was found"
    )
    return set(ns["SCORING_WEIGHTS"])


def _tier_denominator() -> float:
    """Evaluate the denominator ``assign_confidence_tier`` divides by."""
    fn = _find_function("assign_confidence_tier")
    for node in ast.walk(fn):
        if isinstance(node, ast.Assign):
            for tgt in node.targets:
                if isinstance(tgt, ast.Name) and tgt.id == "max_possible":
                    expr = ast.Expression(body=node.value)
                    return float(eval(  # noqa: S307 - our own source
                        compile(expr, "<max_possible>", "eval"),
                        _module_namespace(),
                    ))
    raise AssertionError("max_possible assignment not found in assign_confidence_tier()")


def _attainable_max(weight_keys: set[str]) -> float:
    """What the production scorer returns for an all-axes-perfect candidate."""
    weights = {k: 1.0 for k in weight_keys}
    scores = {k: 1.0 for k in weight_keys}
    return lib.calculate_fair_rank_score(
        scores, weights,
        completeness_floor=COMPLETENESS_FLOOR,
        dnds_reliability=DNDS_RELIABILITY,
    )


# --------------------------------------------------------------------------
# fact 1: tandem_cluster really does reach the numerator
# --------------------------------------------------------------------------

def test_tandem_cluster_is_in_the_scoring_path() -> None:
    """TANDEM_CLUSTER_WEIGHT must be in the weights dict the scorer receives."""
    assert "tandem_cluster" in _scoring_weight_keys()


def test_tandem_cluster_score_changes_rank_score() -> None:
    """Behavioral proof that the tandem axis moves the numerator."""
    keys = _scoring_weight_keys()
    weights = {k: 1.0 for k in keys}
    base = {k: 0.5 for k in keys}
    hot = dict(base, tandem_cluster=1.0)
    assert lib.calculate_fair_rank_score(hot, weights) > \
        lib.calculate_fair_rank_score(base, weights)


# --------------------------------------------------------------------------
# fact 2: the denominator is on the wrong scale entirely
# --------------------------------------------------------------------------

def test_rank_score_numerator_is_bounded_by_one() -> None:
    """The scorer returns a weighted MEAN, not a weight-sum."""
    keys = _scoring_weight_keys()
    assert _attainable_max(keys) == pytest.approx(1.0)


def test_tier_denominator_equals_attainable_max_score() -> None:
    """score_pct must be able to reach 1.0 for a perfect candidate.

    Pre-fix this is 17.0 vs an attainable max of 1.0, so score_pct tops out at
    0.059 and the High tier is unreachable.
    """
    keys = _scoring_weight_keys()
    assert _tier_denominator() == pytest.approx(_attainable_max(keys))


def test_perfect_candidate_clears_the_high_tier_gate() -> None:
    """A candidate maxing every axis must be eligible for High.

    Reproduces the user-visible symptom rather than the arithmetic.
    """
    keys = _scoring_weight_keys()
    score_pct = _attainable_max(keys) / _tier_denominator()
    assert score_pct >= 0.5, (
        f"perfect candidate reaches only score_pct={score_pct:.4f}; the High "
        f"gate (>= 0.5) is unreachable"
    )


def test_medium_tier_score_clause_is_reachable() -> None:
    """The `score_pct >= 0.3` clause must not be dead code."""
    keys = _scoring_weight_keys()
    score_pct = _attainable_max(keys) / _tier_denominator()
    assert score_pct >= 0.3


# --------------------------------------------------------------------------
# anti-drift guards
# --------------------------------------------------------------------------

def test_denominator_covers_every_axis_that_can_score() -> None:
    """ANTI-DRIFT: every axis reaching the numerator must reach the denominator.

    Fails if a future weight is added to the scoring path but not to the
    denominator's source of truth -- the exact drift that lost
    TANDEM_CLUSTER_WEIGHT.
    """
    wrapper = _find_function("calculate_fair_rank_score")
    score_keys = _dict_literal_keys(wrapper, "scores")
    assert score_keys, "scores dict literal not found in calculate_fair_rank_score()"
    assert score_keys == _scoring_weight_keys(), (
        "axes in the scores dict do not match the weights source of truth; "
        f"only in scores: {sorted(score_keys - _scoring_weight_keys())}; "
        f"only in weights: {sorted(_scoring_weight_keys() - score_keys)}"
    )


def test_denominator_is_derived_not_hand_written() -> None:
    """ANTI-DRIFT: assign_confidence_tier must not re-grow a hand-summed list.

    A second hand-maintained sum of *_WEIGHT constants is what rotted the
    first time; forbid it structurally.
    """
    fn = _find_function("assign_confidence_tier")
    named_weights = {
        n.id for n in ast.walk(fn)
        if isinstance(n, ast.Name) and n.id.endswith("_WEIGHT")
    }
    assert not named_weights, (
        "assign_confidence_tier() references weight constants directly "
        f"({sorted(named_weights)}); derive the denominator from the shared "
        "SCORING_WEIGHTS instead of re-summing them by hand"
    )


def test_single_source_of_truth_is_shared_with_sensitivity_analysis() -> None:
    """ANTI-DRIFT: the sensitivity base_weights must not be a third copy."""
    src = RANK.read_text()
    assert "SCORING_WEIGHTS" in src, "no shared SCORING_WEIGHTS constant defined"
    ns = _module_namespace()
    assert set(ns["SCORING_WEIGHTS"]) == _scoring_weight_keys()
    # base_weights for run_sensitivity_analysis must reuse it, not re-literal it.
    assert "base_weights = dict(SCORING_WEIGHTS)" in src or \
        "base_weights = SCORING_WEIGHTS" in src, (
            "sensitivity analysis still builds its own literal weight dict"
        )
