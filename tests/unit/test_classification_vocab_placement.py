"""Pin NaN handling in the LWR gate of the placement classifier.

Audit finding #17 (partial). `classify_placement_for_id` gates on
`p["lwr"] < lwr_threshold`. IEEE-754 makes every comparison with NaN
False, so `NaN < 0.80` is False and a degenerate placement -- one whose
likelihood weight ratio could not be computed at all -- falls straight
through to the CONFIDENT branch and is emitted as a real family call.

A missing LWR is the least confident outcome possible, so it must land in
the same bucket as an LWR below threshold.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

import classify_via_placement as cvp


EDGE_MAP = {7: ("aminergic", "5HT")}


def _placement(lwr, edge_num=7):
    return {"q1": {"lwr": lwr, "edge_num": edge_num}}


# ---------------------------------------------------------------------------
# 1. NaN must not yield a confident call
# ---------------------------------------------------------------------------

def test_nan_lwr_is_unclassified():
    """The bug: NaN < 0.80 is False, so this used to return 'aminergic'."""
    result = cvp.classify_placement_for_id(
        "q1", _placement(float("nan")), EDGE_MAP, lwr_threshold=0.80)
    assert result["family"] == "unclassified-placement"
    assert result["subfamily"] == ""


def test_nan_lwr_is_not_propagated_as_a_score():
    """A NaN must not leak into the output where it would compare False
    against every downstream threshold too."""
    result = cvp.classify_placement_for_id(
        "q1", _placement(float("nan")), EDGE_MAP, lwr_threshold=0.80)
    assert not math.isnan(result["lwr"])
    assert result["lwr"] == 0.0


def test_nan_lwr_via_the_batch_entry_point():
    results = cvp.classify_placement(
        _placement(float("nan")), EDGE_MAP, lwr_threshold=0.80)
    assert results["q1"]["family"] == "unclassified-placement"


# ---------------------------------------------------------------------------
# 2. Other non-finite and out-of-range values
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bad", [float("inf"), float("-inf")])
def test_infinite_lwr_is_unclassified(bad):
    """An LWR is a ratio in [0, 1]; +/-inf is degenerate, not confident."""
    result = cvp.classify_placement_for_id(
        "q1", _placement(bad), EDGE_MAP, lwr_threshold=0.80)
    assert result["family"] == "unclassified-placement"


def test_missing_lwr_key_is_unclassified():
    placements = {"q1": {"edge_num": 7}}
    result = cvp.classify_placement_for_id(
        "q1", placements, EDGE_MAP, lwr_threshold=0.80)
    assert result["family"] == "unclassified-placement"


def test_non_numeric_lwr_is_unclassified():
    result = cvp.classify_placement_for_id(
        "q1", _placement("n/a"), EDGE_MAP, lwr_threshold=0.80)
    assert result["family"] == "unclassified-placement"


# ---------------------------------------------------------------------------
# 3. Guard against over-correction
# ---------------------------------------------------------------------------

def test_confident_placement_still_classifies():
    result = cvp.classify_placement_for_id(
        "q1", _placement(0.95), EDGE_MAP, lwr_threshold=0.80)
    assert result["family"] == "aminergic"
    assert result["subfamily"] == "5HT"
    assert result["lwr"] == 0.95


def test_lwr_exactly_at_threshold_still_classifies():
    """The gate is `< threshold`; equality must remain a pass."""
    result = cvp.classify_placement_for_id(
        "q1", _placement(0.80), EDGE_MAP, lwr_threshold=0.80)
    assert result["family"] == "aminergic"


def test_low_lwr_is_still_unclassified_and_keeps_its_value():
    result = cvp.classify_placement_for_id(
        "q1", _placement(0.42), EDGE_MAP, lwr_threshold=0.80)
    assert result["family"] == "unclassified-placement"
    assert result["lwr"] == 0.42
