import math

import pytest
from scipy.special import betainc

import rank_aggregation as ra


def test_normalized_ranks_drop_missing():
    nr = ra.normalized_ranks({"a": 0.9, "b": float("nan"), "c": 0.1})
    assert set(nr) == {"a", "c"}
    assert nr["a"] < nr["c"]  # a (highest raw score) -> smallest normalized rank


def test_normalized_ranks_lower_is_better():
    nr = ra.normalized_ranks({"a": 0.9, "b": 0.1}, higher_is_better=False)
    assert nr["b"] < nr["a"]  # smaller raw value wins when higher_is_better=False


def test_normalized_ranks_empty_and_singleton():
    assert ra.normalized_ranks({}) == {}
    assert ra.normalized_ranks({"a": 0.5}) == {"a": 1.0}


def test_rra_promotes_consistent_top():
    # 'a' is #1 in every signal; 'z' is last everywhere -> a must rank first.
    lists = {s: {"a": 0.99, "m": 0.5, "z": 0.01} for s in ["s1", "s2", "s3"]}
    order = ra.aggregate(lists, method="rra", groups=None)
    assert order[0] == "a" and order[-1] == "z"


def test_rra_matches_hand_computed_beta():
    # 'a' is ranked #1-of-3 in three identical signals -> its normalized
    # rank is 1/3 in all three; the minimum order-statistic beta occurs at
    # the top order statistic (i=m=3): betainc(3, 1, 1/3) = (1/3)**3 = 1/27,
    # Bonferroni correction 3 * 1/27 = 1/9.
    lists = {s: {"a": 0.99, "m": 0.5, "z": 0.01} for s in ["s1", "s2", "s3"]}
    scores = ra.rra_score(lists)
    m, r = 3, 1.0 / 3.0
    rho = min(betainc(k, m - k + 1, r) for k in range(1, m + 1))
    assert math.isclose(rho, 1.0 / 27.0, rel_tol=1e-9)
    assert math.isclose(scores["a"], min(rho * m, 1.0), rel_tol=1e-9)
    assert math.isclose(scores["a"], 1.0 / 9.0, rel_tol=1e-9)


def test_rrf_matches_formula():
    lists = {"s1": {"a": 0.9, "b": 0.1}, "s2": {"a": 0.1, "b": 0.9}}
    sc = ra.rrf_score(lists, k=60)
    # a: rank1 in s1 (1/61) + rank2 in s2 (1/62); b symmetric -> tie
    assert math.isclose(sc["a"], 1.0 / 61.0 + 1.0 / 62.0, rel_tol=1e-9)
    assert math.isclose(sc["a"], sc["b"], rel_tol=1e-9)


def test_group_fusion_counts_confound_once():
    # 3 identical clone signals favor a>b>c; 1 independent signal favors
    # b>c>a. Ungrouped, the 3 clone votes outvote the 1 independent vote,
    # so 'a' wins. Grouped, the clones fuse into a SINGLE vote, making it
    # a fair 1-vs-1 fight -- 'b' wins because the independent signal ranks
    # 'a' dead last, while the clone group only ranks 'b' 2nd (not last).
    clones = {f"c{i}": {"a": 0.9, "b": 0.5, "c": 0.1} for i in range(3)}
    indep = {"ind": {"a": 0.1, "b": 0.9, "c": 0.5}}
    lists = {**clones, **indep}

    grouped = ra.aggregate(lists, method="rrf", groups=[["c0", "c1", "c2"], ["ind"]])
    ungrouped = ra.aggregate(lists, method="rrf", groups=None)

    assert grouped[0] == "b"     # 1 fused vote each -> independent's pick wins
    assert ungrouped[0] == "a"   # 3 raw clone votes outvote 1 independent vote


def test_partial_presence_rrf_matches_closed_form():
    # 'b' is absent from s2 -- it must not be treated as "worst" there,
    # it simply doesn't receive a vote from s2.
    lists = {"s1": {"a": 0.9, "b": 0.1}, "s2": {"a": 0.8}}
    rrf = ra.rrf_score(lists, k=60)
    assert math.isclose(rrf["a"], 2.0 / 61.0, rel_tol=1e-9)   # rank 1 in both signals
    assert math.isclose(rrf["b"], 1.0 / 62.0, rel_tol=1e-9)   # rank 2 in s1 only


def test_partial_presence_rra_invariant_to_unrelated_signal():
    # adding a signal that never mentions 'b' must not change b's score.
    only_a_and_b = {"s1": {"a": 0.9, "b": 0.1}}
    with_extra_signal = {"s1": {"a": 0.9, "b": 0.1}, "s2": {"a": 0.8, "x": 0.2}}
    base = ra.rra_score(only_a_and_b)
    extended = ra.rra_score(with_extra_signal)
    assert math.isclose(base["b"], extended["b"], rel_tol=1e-9)


def test_aggregate_dispatches_rra_ascending_and_rrf_descending():
    lists = {s: {"a": 0.99, "m": 0.5, "z": 0.01} for s in ["s1", "s2", "s3"]}
    assert ra.aggregate(lists, method="rra", groups=None) == ["a", "m", "z"]
    assert ra.aggregate(lists, method="rrf", groups=None) == ["a", "m", "z"]


def test_aggregate_tie_break_by_id():
    lists = {"s1": {"b": 0.5, "a": 0.5}}
    assert ra.aggregate(lists, method="rra", groups=None) == ["a", "b"]
    assert ra.aggregate(lists, method="rrf", groups=None) == ["a", "b"]


def test_aggregate_invalid_method_raises():
    with pytest.raises(ValueError):
        ra.aggregate({"s1": {"a": 0.5}}, method="bogus", groups=None)
