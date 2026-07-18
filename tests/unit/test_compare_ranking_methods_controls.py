"""Null-anchored positive-control recovery tests.

precision@k with a permutation null, and a leave-one-control-out recovery
curve. These make "known receptors rank high" a null-anchored statement,
without headlining a precision metric when the controls file is ~empty.
"""
from compare_ranking_methods import precision_at_k_vs_null, loo_recovery


def test_precision_at_k_enriched_controls_beats_null():
    order = [f"c{i}" for i in range(100)]          # best-first
    positives = {"c0", "c1", "c2"}                 # all in the top-3
    prec, p = precision_at_k_vs_null(order, positives, k=5, n_perm=2000, seed=0)
    assert prec == 3 / 5
    assert p < 0.01                                # enriched vs random placement


def test_precision_at_k_empty_controls_is_none():
    assert precision_at_k_vs_null([f"c{i}" for i in range(10)], set(), k=5) == (None, None)


def test_loo_recovery_reports_ranks():
    def recompute_fn(held_out):
        base = [f"c{i}" for i in range(10)]
        return [held_out] + [b for b in base if b != held_out]   # always rank 1
    rec = loo_recovery({"c3", "c7"}, recompute_fn)
    assert rec["c3"] == 1 and rec["c7"] == 1
