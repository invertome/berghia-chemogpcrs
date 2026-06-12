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

def test_load_tree_parses_iqtree_slash_support(tmp_path):
    # IQ-TREE writes node support as SH-aLRT/UFBoot (e.g. 57.1/73); ete3 cannot
    # parse that directly. load_tree must keep the UFBoot (second) value.
    p = tmp_path / "t.treefile"
    p.write_text("((a:0.1,b:0.1)57.1/73:0.05,(c:0.1,d:0.1)83.6/100:0.04);")
    t = ead.load_tree(str(p))
    sup = sorted(n.support for n in t.traverse() if not n.is_leaf() and not n.is_root())
    assert sup == [73.0, 100.0]


def test_load_tree_parses_iqtree_three_part_support(tmp_path):
    # With -abayes added to -alrt + -B, IQ-TREE writes a THREE-part support
    # aBayes/SH-aLRT/UFBoot (e.g. 0.072/0/86). load_tree must keep the UFBoot
    # (last) value, not the middle one.
    p = tmp_path / "t3.treefile"
    p.write_text("((a:0.1,b:0.1)0.072/0/86:0.05,(c:0.1,d:0.1)0.965/82.4/98:0.04);")
    t = ead.load_tree(str(p))
    sup = sorted(n.support for n in t.traverse() if not n.is_leaf() and not n.is_root())
    assert sup == [86.0, 98.0]


def test_load_tree_single_number_support(tmp_path):
    p = tmp_path / "t.nwk"
    p.write_text("((a,b)95,(c,d)88);")
    t = ead.load_tree(str(p))
    sup = sorted(n.support for n in t.traverse() if not n.is_leaf() and not n.is_root())
    assert sup == [88.0, 95.0]


def test_anchor_placements_reports_context():
    t = Tree("((berg1,(berg2,ANCHOR_A_3_Y)95)90,ref1);", format=0)
    berghia = {"berg1", "berg2"}
    ingroup = {"berg1", "berg2", "ref1"}
    placements = {p["anchor"]: p for p in
                  ead.anchor_placements(t, ingroup, berghia, support_threshold=80)}
    p = placements["ANCHOR_A_3_Y"]
    assert p["parent_support"] == 95
    assert p["sister_berghia"] == 1          # berg2 is the sister
    assert p["sister_size"] == 1
    assert p["in_berghia_clade"] is True


def test_evaluate_class_includes_placements():
    t = Tree("((berg1,berg2)95,(ANCHOR_A_2_X,ANCHOR_A_3_Y)88);", format=0)
    ingroup = {"berg1", "berg2"}
    res = ead.evaluate_class(t, t, ingroup, berghia_labels=ingroup)
    assert "anchor_placements" in res
    assert {p["anchor"] for p in res["anchor_placements"]} == {"ANCHOR_A_2_X", "ANCHOR_A_3_Y"}


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


# ---------------------------------------------------------------------------
# in-group classification from the pool-membership manifest
# ---------------------------------------------------------------------------

def _mock_is_mollusc(taxid):
    # Aplysia (6500) and Berghia (1287507) are molluscs; Platynereis (6359) is not.
    return taxid in {6447, 6500, 1287507}


def test_load_pool_labels(tmp_path):
    m = tmp_path / "pool_members_class_A.tsv"
    m.write_text(
        "seq_id\ttaxid\tsource\n"
        "berghia_A1\t1287507\tberghia\n"
        "ref_aplysia\t6500\tref\n"          # mollusc ref -> in-group
        "ref_platy\t6359\tref\n"            # non-mollusc ref -> NOT in-group
        "ANCHOR_A_1_P1\t6500\tanchor\n"     # anchor -> never in-group set
    )
    ingroup, berghia = ead.load_pool_labels(str(m), is_mollusc_fn=_mock_is_mollusc)
    assert berghia == {"berghia_A1"}
    assert ingroup == {"berghia_A1", "ref_aplysia"}


def test_infiltration_uses_berghia_labels_not_mollusc_refs():
    # An out-group anchor nested among MOLLUSC REFS (not Berghia) is NOT a
    # Berghia-clade infiltration (spec gate is Berghia clades specifically).
    t = Tree("((berg1,(ref_apl,ANCHOR_A_3_Y)95)90,berg2);", format=0)
    ingroup = {"berg1", "berg2", "ref_apl"}
    berghia = {"berg1", "berg2"}
    assert ead.anchor_infiltrations(t, berghia, support_threshold=80) == []
    res = ead.evaluate_class(t, t, ingroup, berghia_labels=berghia)
    assert res["n_infiltrations"] == 0
