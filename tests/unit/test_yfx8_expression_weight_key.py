"""Bead -yfx8: the expression axis' weight key must be the same string in every
scorer in rank_candidates.py.

THE DEFECT
----------
``SCORING_WEIGHTS`` keys the expression axis ``'expression'`` (and so do the
production scorer ``calculate_fair_rank_score`` and the Latin-Hypercube sweep),
but ``calculate_rank_score``'s ``_axes`` table keyed it ``'expr'`` and looked it
up with ``weights.get(wkey, 0)``. The miss returned 0.0 silently.

``run_sensitivity_analysis`` is fed ``dict(SCORING_WEIGHTS)``, so in every Monte
Carlo iteration the expression axis contributed nothing to ``weighted_sum`` AND
nothing to ``avail_weight``, while ``total_weight`` still counted EXPR_WEIGHT.
Everything downstream therefore measured the wrong function:
``sensitivity_analysis.csv``, the ``rank_stability`` / ``std_rank`` /
``rank_range`` columns merged into the shipped CSV, and
``weight_importance.json`` -- in which ``weight_importance['expression']`` was
guaranteed to report expression as unimportant regardless of the data.

``run_leave_one_out_crossval`` read ``base_weights['expr']`` from a locally
built dict -- a third vocabulary in the same file.

The surviving key is ``'expression'``: it is what the PRODUCTION scorer and
``SCORING_WEIGHTS`` already use, so converging on it changes no production
arithmetic, only repairs the two consumers that had drifted.
"""
from __future__ import annotations

import ast
import os
from typing import Any

import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(REPO, "scripts", "rank_candidates.py")


def _source() -> str:
    with open(SRC) as fh:
        return fh.read()


def _exec_fragment(fn_name: str) -> Any:
    # Slice by AST line range rather than scanning for the next `def`, so the
    # fragment can never swallow module-level code that follows the function.
    src = _source()
    lines = src.splitlines()
    bounds = {n.name: (n.lineno, n.end_lineno) for n in ast.parse(src).body
              if isinstance(n, ast.FunctionDef)}
    assert fn_name in bounds, "%s not found" % fn_name
    lo, hi = bounds[fn_name]
    fragment = "\n".join(lines[lo - 1:hi])
    ns: dict[str, Any] = {"pd": pd, "np": __import__("numpy")}
    exec(compile(fragment, "rank_candidates.py:%s" % fn_name, "exec"), ns)
    return ns[fn_name]


@pytest.fixture(scope="module")
def calc():
    return _exec_fragment("calculate_rank_score")


def _row(expression: float = 0.9) -> dict[str, Any]:
    """A candidate whose ONLY available axis is expression.

    This is the row shape ``run_sensitivity_analysis`` sees: the ``has_*_data``
    flags are exactly the columns rank_candidates.py sets before calling it.
    """
    return {
        "id": "cand_1",
        "phylo_score_norm": 0.0,
        "lse_divergence_score_norm": 0.0,
        "purifying_score_norm": 0.0,
        "positive_score_norm": 0.0,
        "synteny_score_norm": 0.0,
        "expression_score_norm": expression,
        "chemosensory_expr_score_norm": 0.0,
        "gprotein_coexpr_score_norm": 0.0,
        "ecl_divergence_score_norm": 0.0,
        "expansion_score_norm": 0.0,
        "og_confidence_score_norm": 0.0,
        "tandem_cluster_score_norm": 0.0,
        "has_phylo_data": False,
        "has_synteny_data": False,
        "has_expression_data": True,
        "has_chemosensory_expr_data": False,
        "has_gprotein_data": False,
        "has_ecl_data": False,
        "has_expansion_data": False,
        "has_og_confidence_data": False,
        "has_tandem_cluster_data": False,
        "dnds_reliability_weight": 1.0,
    }


# The production weight dict, verbatim in shape.
PROD_WEIGHTS = {
    "phylo": 2.0, "purifying": 0.0, "positive": 2.0, "lse_divergence": 1.0,
    "synteny": 1.0, "expression": 1.0, "chemosensory_expr": 3.0,
    "gprotein_coexpr": 2.0, "ecl_divergence": 2.0, "expansion": 1.0,
    "og_confidence": 1.0, "tandem_cluster": 2.5,
}


def test_expression_axis_reaches_the_sensitivity_scorer(calc) -> None:
    """Pre-fix this is 0.0: `_axes` asks for 'expr', SCORING_WEIGHTS has
    'expression', `weights.get('expr', 0)` misses, weight 0 -> avail_weight 0
    -> the scorer's `avail_weight <= 0` guard returns 0.0."""
    out = calc(pd.DataFrame([_row(0.9)]), PROD_WEIGHTS).iloc[0]
    assert out > 0.0, "expression must contribute to the Monte Carlo score"


def test_sensitivity_scorer_matches_production_scorer_on_expression() -> None:
    """The Monte Carlo scorer and the production scorer must agree when the
    expression axis is the one carrying the signal.

    The pre-existing anti-drift guard
    (test_rank_candidates_dnds_reliability.test_sensitivity_matches_production_lib_formula)
    could NOT see this drift: it passed a dict keyed 'expr' to one scorer and a
    hand-written dict keyed 'expression' to the other, and agreed only because
    no fixture row set has_expression_data=True. This row does.
    """
    import _rank_candidates_lib as lib

    calc_fn = _exec_fragment("calculate_rank_score")
    s_sens = calc_fn(pd.DataFrame([_row(0.9)]), PROD_WEIGHTS).iloc[0]

    # Mirror the production scores dict exactly (rank_candidates.py
    # calculate_fair_rank_score): the dN/dS axes are NOT has_*_data gated, so
    # they are always passed as present values, and only the gated axes go None.
    scores = {k: None for k in PROD_WEIGHTS}
    scores["purifying"] = 0.0
    scores["positive"] = 0.0
    scores["expression"] = 0.9
    s_prod = lib.calculate_fair_rank_score(
        scores, PROD_WEIGHTS, completeness_floor=0.4, dnds_reliability=1.0)
    assert s_sens == pytest.approx(s_prod)


def test_expression_weight_actually_changes_the_score(calc) -> None:
    """weight_importance['expression'] can only be meaningful if perturbing the
    expression weight moves the score."""
    lo = dict(PROD_WEIGHTS, expression=0.5)
    hi = dict(PROD_WEIGHTS, expression=3.0)
    s_lo = calc(pd.DataFrame([_row(0.9)]), lo).iloc[0]
    s_hi = calc(pd.DataFrame([_row(0.9)]), hi).iloc[0]
    assert s_lo != pytest.approx(s_hi)


def test_one_vocabulary_for_the_expression_weight_key() -> None:
    """Static guard: no scorer in rank_candidates.py may key expression 'expr'.

    A value test alone cannot pin this down -- the LOO cross-validation path
    (``base_weights['expr']``) is fed its own locally built dict, so it is
    internally consistent while still disagreeing with SCORING_WEIGHTS. This
    catches any of the three sites drifting again.
    """
    tree = ast.parse(_source())
    offenders = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and node.value == "expr":
            offenders.append(node.lineno)
    assert not offenders, (
        "rank_candidates.py still uses the 'expr' weight key at line(s) %s; "
        "the single vocabulary is 'expression'" % offenders)


def test_crossval_and_scoring_weights_share_their_keys() -> None:
    """The LOO dict built at module scope must be a strict subset of
    SCORING_WEIGHTS' keys, so ``base_weights[k]`` can never KeyError and can
    never silently address a different axis."""
    src = _source()
    tree = ast.parse(src)

    def _dict_keys_of(assign_name):
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign) and isinstance(node.value, ast.Dict):
                for tgt in node.targets:
                    if isinstance(tgt, ast.Name) and tgt.id == assign_name:
                        return {k.value for k in node.value.keys
                                if isinstance(k, ast.Constant)}
        return None

    scoring = _dict_keys_of("SCORING_WEIGHTS")
    base = _dict_keys_of("base_weights")
    assert scoring, "SCORING_WEIGHTS literal not found"
    assert base, "base_weights literal not found"
    assert base <= scoring, "cross-validation weight keys drifted: %s" % (base - scoring)
