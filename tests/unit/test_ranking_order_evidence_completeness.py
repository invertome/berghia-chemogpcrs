"""`evidence_completeness` must measure AVAILABILITY, not value.

rank_candidates.py's docstring for the block says the metric is "based on
data availability, not just score > 0" -- then uses `score > 0` for three of
its eleven sources:

    phylo_score > 0
    purifying_score > 0 or positive_score > 0
    lse_score > 0

The `lse_divergence` case is structural, not incidental. `lse_threshold` is
``np.percentile(all_depths, LSE_DIVERGENCE_PERCENTILE)`` with the percentile at
75, and get_lse_divergence_score() returns 0.0 for any candidate at or below that
threshold. So ~75% of candidates can NEVER score > 0 on lse_divergence BY
CONSTRUCTION, yet each is recorded as missing that evidence -- their depth
was measured, it just wasn't extreme.

Reproduced: a candidate with 7 genuine evidence sources plus a
measured-but-subthreshold depth reports 7/11 = 0.636 when the true
availability is 8/11 = 0.727. CONFIDENCE_MIN_COMPLETENESS defaults to 0.7,
so the understatement flips the confidence-view gate from True to False and
drops the candidate off the safe-bet shortlist.

`phylo` already has an availability flag (`has_phylo_data`, bead o98 -- tree
membership). `purifying`/`positive` are the only axes with no `has_*_data`
gate anywhere in the codebase; the fix adds `has_dnds_data`.

rank_candidates.py is not import-safe (top-level ``sys.argv``); these tests
exec a single top-level function fragment, the pattern established in
tests/unit/test_rank_candidates_dnds_reliability.py.
"""
from __future__ import annotations

import ast
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RANK = PROJECT_ROOT / "scripts" / "rank_candidates.py"


def _exec_function(name: str) -> Any:
    src = RANK.read_text()
    start = src.find(f"def {name}(")
    assert start != -1, f"{name} not found in rank_candidates.py"
    end = src.find("\n\ndef ", start + 1)
    fragment = src[start:end] if end != -1 else src[start:]
    ns: dict[str, Any] = {"pd": pd, "np": np, "os": __import__("os")}
    exec(fragment, ns)
    return ns[name]


@pytest.fixture(scope="module")
def completeness():
    return _exec_function("evidence_completeness")


ALL_ABSENT = dict(
    has_phylo_data=False, has_dnds_data=False, has_synteny=False,
    has_expression=False, has_lse_data=False, has_chemo_expr=False,
    has_gprotein=False, has_ecl=False, has_expansion=False,
    has_og_data=False, has_tandem=False,
)


def test_no_evidence_is_zero(completeness) -> None:
    assert completeness(**ALL_ABSENT) == 0.0


def test_all_evidence_is_one(completeness) -> None:
    assert completeness(**{k: True for k in ALL_ABSENT}) == 1.0


def test_eleven_sources_are_counted(completeness) -> None:
    """The denominator must stay 11; a 12th axis without a flag would skew it."""
    one = dict(ALL_ABSENT, has_tandem=True)
    assert completeness(**one) == pytest.approx(1 / 11)


def test_measured_but_subthreshold_depth_counts_as_present(completeness) -> None:
    """The reproduction: 7 genuine sources + a measured subthreshold depth.

    Because lse_threshold is the 75th percentile of all depths, a candidate
    scoring 0.0 on lse_divergence is the COMMON case, not the missing-data case.
    Counting it as absent understates completeness by 1/11 and crosses the
    0.7 CONFIDENCE_MIN_COMPLETENESS gate.
    """
    args = dict(
        ALL_ABSENT,
        has_phylo_data=True, has_dnds_data=True, has_synteny=True,
        has_expression=True, has_chemo_expr=True, has_og_data=True,
        has_tandem=True,
        has_lse_data=True,  # in the tree; depth measured, just not extreme
    )
    value = completeness(**args)
    assert value == pytest.approx(8 / 11), (
        f"got {value:.3f} (= {value * 11:.0f}/11); a measured-but-subthreshold "
        f"depth was scored as MISSING evidence"
    )
    assert value >= 0.7, "must clear the CONFIDENCE_MIN_COMPLETENESS gate"


def test_absent_from_the_tree_really_is_missing(completeness) -> None:
    """The genuine missing case must still be missing (no blanket credit)."""
    args = dict(
        ALL_ABSENT,
        has_dnds_data=True, has_synteny=True, has_expression=True,
        has_chemo_expr=True, has_og_data=True, has_tandem=True,
        has_phylo_data=False, has_lse_data=False,
    )
    assert completeness(**args) == pytest.approx(6 / 11)


def test_dnds_availability_is_independent_of_omega_value(completeness) -> None:
    """A candidate WITH aBSREL results scores as having that evidence.

    Under PURIFYING_WEIGHT=0 a purely-purifying candidate has
    purifying_score==0 and positive_score==0, so the old `score > 0` test
    called its dN/dS evidence missing even though aBSREL ran and reported.
    """
    with_dnds = completeness(**dict(ALL_ABSENT, has_dnds_data=True))
    without = completeness(**ALL_ABSENT)
    assert with_dnds > without


# --------------------------------------------------------------------------
# structural: the availability flags must reach the caller and the CSV
# --------------------------------------------------------------------------

def test_completeness_is_computed_by_the_shared_function() -> None:
    """The main loop must CALL evidence_completeness, not inline a copy.

    A hand-maintained second copy is the drift that produced this bug class
    (see tests/unit/test_tier_normalization_denominator.py).
    """
    src = RANK.read_text()
    tree = ast.parse(src)
    calls = [n for n in ast.walk(tree)
             if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
             and n.func.id == "evidence_completeness"]
    assert calls, "evidence_completeness() is never called"
    assert "evidence_sources = [" not in src, (
        "the inline evidence_sources list is still present; it will drift "
        "from evidence_completeness()"
    )


def test_has_dnds_data_is_emitted() -> None:
    """`has_dnds_data` must exist and reach the CSV like every other flag."""
    src = RANK.read_text()
    assert "'has_dnds_data'" in src or '"has_dnds_data"' in src, (
        "no has_dnds_data anywhere; purifying/positive remain the only axes "
        "with no availability gate"
    )
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and node.col_offset == 0:
            for tgt in node.targets:
                if (isinstance(tgt, ast.Name) and tgt.id == "output_cols"
                        and isinstance(node.value, ast.List)):
                    cols = [e.value for e in node.value.elts
                            if isinstance(e, ast.Constant)]
                    assert "has_dnds_data" in cols
                    return
    pytest.fail("no module-level output_cols list literal found")
