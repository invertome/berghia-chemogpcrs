"""Tests for the amino-acid composition-distance confound producer (A2)."""
import math

import numpy as np

from composition_distance import aa_composition, composition_distance

_AA = "ACDEFGHIKLMNPQRSTVWY"


def test_aa_composition_is_a_normalized_frequency_vector():
    v = aa_composition("AAAA")
    assert v.shape == (20,)
    assert math.isclose(v.sum(), 1.0)
    assert math.isclose(v[_AA.index("A")], 1.0)


def test_composition_is_length_normalized():
    # frequency, not count: "AA" and "AAAA" have identical composition
    assert np.allclose(aa_composition("AA"), aa_composition("AAAA"))


def test_non_standard_residues_are_ignored():
    # X / gaps do not contribute; "AAXX" is pure-A composition
    assert np.allclose(aa_composition("AAXX"), aa_composition("AA"))


def test_distance_to_reference_centroid_is_zero_for_typical_composition():
    # a candidate whose composition equals the reference centroid -> distance 0
    refs = {"r1": "ACDE", "r2": "ACDE"}
    cands = {"typical": "ACDE", "odd": "WWWW"}
    d = composition_distance(cands, reference_seqs=refs)
    assert math.isclose(d["typical"], 0.0, abs_tol=1e-9)
    assert d["odd"] > d["typical"]


def test_self_centered_when_no_reference_given():
    # centroid falls back to the mean composition of the candidates themselves
    cands = {"a": "AAAA", "b": "AAAA", "c": "YYYY"}
    d = composition_distance(cands)
    # 'c' (the lone outlier) is farther from the self-centroid than the 'a'/'b' pair
    assert d["c"] > d["a"]


def test_deterministic():
    cands = {"a": "ACDEFG", "b": "WWWWYY"}
    refs = {"r": "ACDEFG"}
    assert composition_distance(cands, refs) == composition_distance(cands, refs)
