"""Tests for the EPA-ng placement -> distance-to-nearest-reference producer (bead amiu).

Generalizes the A1 tree-distance confound: instead of rebuilding a de-novo tree
per species, candidates are PLACED on a fixed anchor backbone and each one's
distance to the nearest reference leaf is read from the placement.

Reference backbone used throughout::

    ((A:1{0},B:1{1}):1{2},C:2{3}):0{4};

so, writing `below[x]` = nearest-leaf distance below x and `above[x]` = nearest
leaf reached by going up/around from x:
    below: A=0 B=0 C=0 AB=1 root=2
    above: A=2 B=2 AB=3 C=4 root=inf
"""
import math

import pytest

from placement_distance import (
    load_backbone,
    nearest_ref_distance_from_jplace,
)

_TREE = "((A:1{0},B:1{1}):1{2},C:2{3}):0{4};"
_FIELDS = ["edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"]


def _jplace(placements, tree=_TREE, fields=None):
    return {"tree": tree, "fields": fields or _FIELDS,
            "placements": placements, "version": 3, "metadata": {}}


def test_edge_labels_map_to_the_right_nodes():
    _tree, edge_map = load_backbone(_TREE)
    assert edge_map[0].name == "A"
    assert edge_map[1].name == "B"
    assert edge_map[3].name == "C"
    # edge 2 is the internal (A,B) node's parent edge -> its leaves are A and B
    assert sorted(l.name for l in edge_map[2].get_leaves()) == ["A", "B"]


def test_distance_on_a_terminal_edge():
    # placed on edge 0 (the edge leading to leaf A, length 1), 0.25 from A,
    # pendant 0.5  ->  0.25 (down to A) + 0.5 = 0.75
    out = nearest_ref_distance_from_jplace(
        _jplace([{"n": ["q1"], "p": [[0, -100.0, 1.0, 0.25, 0.5]]}]))
    assert out["q1"] == pytest.approx(0.75)


def test_distance_on_an_internal_edge_takes_the_nearer_side():
    # placed on edge 2 (edge to the (A,B) node, length 1), 0.5 from that node,
    # pendant 0.1. down: 0.5 + below[AB]=1 -> 1.5 ; up: above[AB]=3 - 0.5 -> 2.5
    out = nearest_ref_distance_from_jplace(
        _jplace([{"n": ["q2"], "p": [[2, -100.0, 1.0, 0.5, 0.1]]}]))
    assert out["q2"] == pytest.approx(1.6)


def test_best_like_weight_ratio_placement_wins():
    # two candidate placements for one query -> the higher LWR one is used
    out = nearest_ref_distance_from_jplace(_jplace([{
        "n": ["q3"],
        "p": [[0, -100.0, 0.1, 0.25, 0.5],     # would give 0.75
              [2, -90.0, 0.9, 0.5, 0.1]],      # higher LWR -> 1.6
    }]))
    assert out["q3"] == pytest.approx(1.6)


def test_field_order_is_read_from_the_fields_header():
    # jplace permits any field order; the parser must not assume positions
    fields = ["distal_length", "edge_num", "pendant_length", "like_weight_ratio"]
    out = nearest_ref_distance_from_jplace(
        _jplace([{"n": ["q4"], "p": [[0.25, 0, 0.5, 1.0]]}], fields=fields))
    assert out["q4"] == pytest.approx(0.75)


def test_ref_ids_restricts_which_leaves_count_as_references():
    # only C is a reference: from the edge-0 attachment (0.25 above A) the route
    # is 0.75 up to AB, 1 to root, 2 down to C = 3.75, plus pendant 0.5
    out = nearest_ref_distance_from_jplace(
        _jplace([{"n": ["q5"], "p": [[0, -100.0, 1.0, 0.25, 0.5]]}]),
        ref_ids={"C"})
    assert out["q5"] == pytest.approx(4.25)


def test_multiplicity_name_form_is_supported():
    # jplace allows "nm": [[name, multiplicity], ...] instead of "n"
    out = nearest_ref_distance_from_jplace(
        _jplace([{"nm": [["q6", 1]], "p": [[0, -100.0, 1.0, 0.25, 0.5]]}]))
    assert out["q6"] == pytest.approx(0.75)


def test_placement_without_any_reference_leaf_is_omitted():
    # no leaf qualifies as a reference -> no finite distance, candidate dropped
    out = nearest_ref_distance_from_jplace(
        _jplace([{"n": ["q7"], "p": [[0, -100.0, 1.0, 0.25, 0.5]]}]),
        ref_ids={"NOT_A_LEAF"})
    assert "q7" not in out


def test_empty_placement_list_is_skipped_not_an_error():
    out = nearest_ref_distance_from_jplace(_jplace([{"n": ["q8"], "p": []}]))
    assert out == {}


def test_real_epa_ng_placement_pins_the_distal_length_convention():
    """Regression test built from ACTUAL EPA-ng output, not a hand-made fixture.

    A query that is an exact copy of leaf ``A`` was placed on A's terminal edge
    (length 1.0); EPA-ng reported distal_length=5e-5 and pendant_length=1e-4.
    Because the query is identical to A it must sit AT the tip, so distal_length
    is measured from the distal/child node. The distance is therefore
    pendant + distal = 1.5e-4. If the convention were read from the parent end
    instead, this would come out ~1.0001 — so this test is what keeps the
    orientation from silently flipping.
    """
    tree = ("((A:1.0000000000{0},B:0.0500000000{1}):0.1000000000{2},"
            "(C:0.0500000000{3},D:0.0500000000{4}):0.1000000000{5});")
    out = nearest_ref_distance_from_jplace(
        _jplace([{"n": ["q_copy_of_A"], "p": [[0, -100.0, 1.0, 0.00005, 0.0001]]}],
                tree=tree))
    assert out["q_copy_of_A"] == pytest.approx(0.00015, abs=1e-9)
    assert out["q_copy_of_A"] < 0.5, "distal_length must not be read from the parent end"


def test_distances_are_finite_and_non_negative():
    out = nearest_ref_distance_from_jplace(
        _jplace([{"n": ["q9"], "p": [[3, -100.0, 1.0, 1.0, 0.2]]}]))
    assert math.isfinite(out["q9"]) and out["q9"] >= 0
