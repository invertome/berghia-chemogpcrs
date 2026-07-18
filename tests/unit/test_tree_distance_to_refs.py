"""Tests for the class-A tree-distance-to-nearest-reference producer (A1)."""
import math

import pytest

from tree_distance_to_refs import nearest_ref_distance

# ((A:1,B:1):1,(C:1,R:1):1);  — R is the only reference.
#   C..R = 1 + 1                 = 2   (siblings under (C,R))
#   A..R = 1 + 1 + 1 + 1         = 4   (up to root, down to R)
#   B..R = 4
_NEWICK = "((A:1,B:1):1,(C:1,R:1):1);"


def test_min_patristic_distance_to_single_reference():
    d = nearest_ref_distance(_NEWICK, ref_ids={"R"})
    assert math.isclose(d["C"], 2.0)
    assert math.isclose(d["A"], 4.0)
    assert math.isclose(d["B"], 4.0)
    assert "R" not in d  # a reference leaf is not scored as a candidate


def test_takes_the_nearest_of_multiple_references():
    # ((A:1,B:1):1,(C:1,R:1):1); with BOTH B and R as references:
    #   A..B = 2, A..R = 4  -> nearest 2
    #   C..R = 2, C..B = 4  -> nearest 2
    d = nearest_ref_distance(_NEWICK, ref_ids={"B", "R"})
    assert math.isclose(d["A"], 2.0)
    assert math.isclose(d["C"], 2.0)
    assert set(d) == {"A", "C"}  # only non-reference leaves scored


def test_candidate_ids_restricts_the_scored_set():
    d = nearest_ref_distance(_NEWICK, ref_ids={"R"}, candidate_ids=["A"])
    assert set(d) == {"A"}
    assert math.isclose(d["A"], 4.0)


def test_candidate_absent_from_tree_is_omitted():
    d = nearest_ref_distance(_NEWICK, ref_ids={"R"}, candidate_ids=["A", "GHOST"])
    assert set(d) == {"A"}


def test_parses_iqtree_style_support_labels():
    # IQ-TREE .treefile carries internal support labels; distances use only
    # branch lengths, so the producer must still resolve leaf distances.
    nwk = "((A:1,B:1)95:1,(C:1,R:1)80:1);"
    d = nearest_ref_distance(nwk, ref_ids={"R"})
    assert math.isclose(d["C"], 2.0)
