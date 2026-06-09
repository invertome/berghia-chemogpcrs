"""Unit tests for evaluate_anchor_divergence.py (C3).

The trees themselves are built on Unity (stage-04 align -> filter -> IQ-TREE);
this module is the eval that, given a without-anchor / with-anchor tree pair per
class, decides whether the out-group anchors (tiers 2-3) are safe to keep:
  - placement: an out-group anchor must NOT nest inside a well-supported
    in-group (Berghia/mollusc) clade,
  - in-group topology: adding anchors must not perturb the in-group tree
    (restricted Robinson-Foulds),
  - support: median in-group clade support must not drop materially.
Tier-1 (in-group) anchors are never gated.

Tests run on tiny synthetic ete3 trees — no tree inference here.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest
from ete3 import Tree

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import evaluate_anchor_divergence as ead  # noqa: E402


# ---------------------------------------------------------------------------
# anchor label parsing
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("label,tier", [
    ("ANCHOR_A_1_P31356", "1"),
    ("ANCHOR_A_2_PLAT1", "2"),
    ("ANCHOR_F_3_Q9VXD9", "3"),
    ("berghia_A1", None),
    ("ref_123", None),
])
def test_parse_anchor_tier(label, tier):
    assert ead.parse_anchor_tier(label) == tier


@pytest.mark.parametrize("label,expected", [
    ("ANCHOR_A_2_X", True),    # tier 2 = out-group
    ("ANCHOR_F_3_Y", True),    # tier 3 = out-group
    ("ANCHOR_A_1_Z", False),   # tier 1 = in-group, never gated
    ("berghia_1", False),
])
def test_is_outgroup_anchor(label, expected):
    assert ead.is_outgroup_anchor(label) == expected


# ---------------------------------------------------------------------------
# placement: out-group anchor nested in a supported in-group clade
# ---------------------------------------------------------------------------

def test_no_infiltration_when_anchors_form_their_own_clade():
    # Anchors are sisters to each other, outside the Berghia clade -> coherent.
    t = Tree("((berg1,berg2)95,(ANCHOR_A_2_X,ANCHOR_A_3_Y)90);", format=0)
    ingroup = {"berg1", "berg2"}
    assert ead.anchor_infiltrations(t, ingroup, support_threshold=80) == []


def test_infiltration_when_anchor_nested_in_ingroup_clade():
    # ANCHOR_A_3_Y sits inside a support-95 clade whose only other member is an
    # in-group leaf -> infiltration.
    t = Tree("((berg1,(berg2,ANCHOR_A_3_Y)95)90,ref1);", format=0)
    ingroup = {"berg1", "berg2", "ref1"}
    assert ead.anchor_infiltrations(t, ingroup, support_threshold=80) == ["ANCHOR_A_3_Y"]


def test_low_support_clade_is_not_infiltration():
    # Same nesting but the enclosing clade is weakly supported (<80) -> not gated.
    t = Tree("((berg1,(berg2,ANCHOR_A_3_Y)40)90,ref1);", format=0)
    ingroup = {"berg1", "berg2", "ref1"}
    assert ead.anchor_infiltrations(t, ingroup, support_threshold=80) == []


def test_tier1_anchor_not_gated():
    # A tier-1 (in-group) anchor nested among Berghia tips is NOT flagged.
    t = Tree("((berg1,(berg2,ANCHOR_A_1_Z)95)90,ref1);", format=0)
    ingroup = {"berg1", "berg2", "ref1"}
    assert ead.anchor_infiltrations(t, ingroup, support_threshold=80) == []


# ---------------------------------------------------------------------------
# in-group topology: restricted Robinson-Foulds
# ---------------------------------------------------------------------------

def test_restricted_rf_zero_when_ingroup_topology_unchanged():
    without = Tree("((berg1,berg2),(berg3,berg4));")
    with_anchor = Tree("(((berg1,berg2),(berg3,berg4)),ANCHOR_A_2_X);")
    shared = {"berg1", "berg2", "berg3", "berg4"}
    assert ead.restricted_rf(without, with_anchor, shared) == 0.0


def test_restricted_rf_positive_when_ingroup_topology_changes():
    without = Tree("((berg1,berg2),(berg3,berg4));")
    with_anchor = Tree("(((berg1,berg3),(berg2,berg4)),ANCHOR_A_2_X);")
    shared = {"berg1", "berg2", "berg3", "berg4"}
    assert ead.restricted_rf(without, with_anchor, shared) > 0.0


# ---------------------------------------------------------------------------
# support: median in-group clade support drop
# ---------------------------------------------------------------------------

def test_ingroup_support_drop():
    without = Tree("((berg1,berg2)95,(berg3,berg4)90);", format=0)
    with_anchor = Tree("(((berg1,berg2)80,(berg3,berg4)78),ANCHOR_A_2_X);", format=0)
    ingroup = {"berg1", "berg2", "berg3", "berg4"}
    drop = ead.ingroup_support_drop(without, with_anchor, ingroup)
    # median 92.5 -> 79 => drop ~13.5
    assert drop == pytest.approx(13.5, abs=0.6)


# ---------------------------------------------------------------------------
# verdict
# ---------------------------------------------------------------------------

def test_evaluate_class_includes_clean_anchors():
    without = Tree("((berg1,berg2)95,(berg3,berg4)90);", format=0)
    with_anchor = Tree("(((berg1,berg2)95,(berg3,berg4)90),(ANCHOR_A_2_X,ANCHOR_A_3_Y)88);", format=0)
    ingroup = {"berg1", "berg2", "berg3", "berg4"}
    res = ead.evaluate_class(without, with_anchor, ingroup)
    assert res["verdict"] == "include"
    assert res["n_infiltrations"] == 0
    assert res["rf"] == 0.0


def test_evaluate_class_excludes_on_infiltration():
    without = Tree("((berg1,berg2)95,(berg3,berg4)90);", format=0)
    with_anchor = Tree("((berg1,(berg2,ANCHOR_A_3_Y)95)90,(berg3,berg4)90);", format=0)
    ingroup = {"berg1", "berg2", "berg3", "berg4"}
    res = ead.evaluate_class(without, with_anchor, ingroup)
    assert res["verdict"] == "exclude"
    assert res["n_infiltrations"] >= 1


def test_evaluate_class_excludes_on_support_drop():
    without = Tree("((berg1,berg2)95,(berg3,berg4)92);", format=0)
    with_anchor = Tree("(((berg1,berg2)70,(berg3,berg4)68),ANCHOR_A_2_X);", format=0)
    ingroup = {"berg1", "berg2", "berg3", "berg4"}
    res = ead.evaluate_class(without, with_anchor, ingroup, support_drop_threshold=5)
    assert res["verdict"] == "exclude"
    assert res["support_drop"] > 5
