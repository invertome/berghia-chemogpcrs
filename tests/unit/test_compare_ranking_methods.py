"""Tests for the weighted-vs-rank-aggregation comparison + honest uncertainty.

Covers the pure ordering primitives (Spearman, top-k overlap, biggest movers),
the descriptive permutation-null enrichment, the <5-positives precision guard,
and the weighted-path regression guarantee (the RANK_METHOD toggle leaves the
default ordering untouched).
"""
import os

import numpy as np
import pandas as pd
import pytest

import compare_ranking_methods as cmp
import rank_aggregation as ra


# --------------------------------------------------------------------------- #
# Fixtures / helpers
# --------------------------------------------------------------------------- #
def _toy_df(n=8):
    """A ranked df with DISTINCT rank_score (no tie ambiguity) plus the four
    always-present base signals and one gated signal."""
    ids = [f"c{i}" for i in range(n)]
    rank_score = list(range(n, 0, -1))  # c0 best ... c{n-1} worst, all distinct
    return pd.DataFrame(
        {
            "id": ids,
            "rank_score": rank_score,
            "phylo_score_norm": np.linspace(1.0, 0.0, n),
            "purifying_score_norm": np.linspace(0.0, 1.0, n),
            "positive_score_norm": np.linspace(1.0, 0.0, n),
            "lse_divergence_score_norm": np.linspace(1.0, 0.0, n),
            "synteny_score_norm": np.linspace(0.5, 0.5, n),
            "has_synteny_data": [True] * n,
        }
    )


# --------------------------------------------------------------------------- #
# Ordering primitives on explicit toy orderings
# --------------------------------------------------------------------------- #
def test_spearman_of_orders_perfect_and_reversed():
    a = ["a", "b", "c", "d"]
    assert cmp.spearman_of_orders(a, a) == pytest.approx(1.0)
    assert cmp.spearman_of_orders(a, list(reversed(a))) == pytest.approx(-1.0)


def test_topk_overlap_counts_shared_prefix():
    a = ["a", "b", "c", "d"]
    b = ["a", "c", "x", "y"]
    ov = cmp.topk_overlap(a, b, ks=(2, 4))
    assert ov[2] == pytest.approx(0.5)   # {a,b} vs {a,c} -> share {a}
    # top-4: {a,b,c,d} vs {a,c,x,y} -> share {a,c} = 2/4
    assert ov[4] == pytest.approx(0.5)


def test_biggest_movers_direction_and_delta():
    # weighted: e last; rankagg: e first -> e is the biggest UP mover.
    weighted = ["a", "b", "c", "d", "e"]
    rankagg = ["e", "a", "b", "c", "d"]
    up, down = cmp.biggest_movers(weighted, rankagg, n=3)
    assert up[0][0] == "e"
    # e: weighted pos 5 (0-based 4) -> rankagg pos 1 (0-based 0); delta = 4
    assert up[0][1] == 4
    # a..d each slipped one spot -> deltas -1; they are the DOWN movers
    assert all(m[1] == -1 for m in down)


# --------------------------------------------------------------------------- #
# compare(df) end-to-end structure
# --------------------------------------------------------------------------- #
def test_compare_returns_expected_shape_and_weighted_order():
    df = _toy_df(8)
    res = cmp.compare(df)
    assert res["weighted_order"] == [f"c{i}" for i in range(8)]  # by rank_score
    assert set(res.keys()) >= {
        "weighted_order", "rankagg_order", "spearman",
        "topk_overlap", "movers_up", "movers_down",
    }
    assert -1.0 <= res["spearman"] <= 1.0
    assert set(res["rankagg_order"]) == set(res["weighted_order"])


# --------------------------------------------------------------------------- #
# Honest uncertainty: descriptive permutation-null enrichment
# --------------------------------------------------------------------------- #
def test_enrichment_planted_top_positive_gives_low_p():
    order = [f"g{i}" for i in range(100)]
    planted = order[0]  # ranked first
    res = cmp.enrichment_vs_null(order, [planted], n_perm=2000, seed=0)
    assert res["positions"][planted] == 1
    assert res["observed_mean_rank"] == pytest.approx(1.0)
    assert res["p_value"] < 0.05          # enriched at the top vs chance


def test_enrichment_planted_bottom_positive_gives_high_p():
    order = [f"g{i}" for i in range(100)]
    planted = order[-1]  # ranked last
    res = cmp.enrichment_vs_null(order, [planted], n_perm=2000, seed=0)
    assert res["positions"][planted] == 100
    assert res["p_value"] > 0.5


def test_enrichment_no_positives_is_graceful():
    order = [f"g{i}" for i in range(10)]
    res = cmp.enrichment_vs_null(order, [], n_perm=100, seed=0)
    assert res["n_positives"] == 0
    assert res["p_value"] is None
    assert res["positions"] == {}


# --------------------------------------------------------------------------- #
# <5 positives => NO precision@k (the guard), but the "insufficient" note
# --------------------------------------------------------------------------- #
def test_report_under_five_positives_omits_precision_at_k(tmp_path):
    df = _toy_df(8)
    out = tmp_path / "ranking_method_comparison.md"
    # one matched positive -> below the 5-positive floor
    text = cmp.write_report(df, str(out), positive_ids=["c0"], n_perm=200, seed=0)
    assert "precision@" not in text
    assert "insufficient positives for precision estimate" in text
    # file written and matches returned text
    assert out.read_text() == text


def test_report_five_or_more_positives_may_report_precision(tmp_path):
    df = _toy_df(12)
    out = tmp_path / "ranking_method_comparison.md"
    positives = ["c0", "c1", "c2", "c3", "c4"]  # exactly 5 -> allowed
    text = cmp.write_report(df, str(out), positive_ids=positives, n_perm=200, seed=0)
    assert "precision@" in text
    assert "insufficient positives for precision estimate" not in text


def test_report_zero_positives_states_expected_empty(tmp_path):
    df = _toy_df(8)
    out = tmp_path / "ranking_method_comparison.md"
    text = cmp.write_report(df, str(out), positive_ids=[], n_perm=100, seed=0)
    assert "precision@" not in text
    assert out.exists()


# --------------------------------------------------------------------------- #
# Best-effort positive-control mapping (tolerates 0 matches)
# --------------------------------------------------------------------------- #
def test_map_positive_controls_tolerates_no_matches(tmp_path):
    controls = tmp_path / "hcr_positive_controls.csv"
    controls.write_text(
        "gene_name,aliases,refseq_protein,expected_tissue\n"
        "Galpha_olf,GNAL;Gaolf,,rhinophore\n"
    )
    matched = cmp.map_positive_controls(str(controls), ["c0", "c1", "BersteEVm001"])
    assert matched == []  # empty refseq + non-matching aliases


def test_map_positive_controls_exact_refseq_match(tmp_path):
    controls = tmp_path / "hcr_positive_controls.csv"
    controls.write_text(
        "gene_name,aliases,refseq_protein,expected_tissue\n"
        "Foo,FOO1;FOO2,XP_012345.1,rhinophore\n"
    )
    matched = cmp.map_positive_controls(str(controls), ["XP_012345.1", "c1"])
    assert matched == ["XP_012345.1"]


# --------------------------------------------------------------------------- #
# WEIGHTED-PATH REGRESSION: the RANK_METHOD toggle must not alter the default
# --------------------------------------------------------------------------- #
def test_rerank_output_weighted_is_identity():
    df = _toy_df(8)
    df_sorted = df.sort_values("rank_score", ascending=False).reset_index(drop=True)
    # explicit "weighted"
    out = ra.rerank_output(df_sorted, "weighted")
    assert list(out["id"]) == list(df_sorted["id"])
    # unset env resolves to the weighted default -> still identity
    out2 = ra.rerank_output(df_sorted, os.getenv("RANK_METHOD", "weighted"))
    assert list(out2["id"]) == list(df_sorted["id"])


def test_weighted_order_equals_sort_by_rank_score():
    # compare's reconstruction of the weighted order must equal exactly what
    # rank_candidates.py does: df.sort_values('rank_score', ascending=False).
    df = _toy_df(8)
    expected = df.sort_values("rank_score", ascending=False)["id"].tolist()
    assert cmp.weighted_order(df) == expected


def test_rerank_output_rankagg_actually_reorders():
    # A df whose rankagg order differs from the weighted order, to prove the
    # rankagg branch is live (not a silent no-op).
    df = pd.DataFrame(
        {
            "id": ["a", "b", "c"],
            "rank_score": [3.0, 2.0, 1.0],          # weighted: a, b, c
            "phylo_score_norm": [0.0, 0.0, 1.0],    # rankagg favors c
            "purifying_score_norm": [0.0, 0.0, 1.0],
            "positive_score_norm": [0.0, 0.0, 1.0],
            "lse_divergence_score_norm": [0.0, 0.0, 1.0],
        }
    )
    weighted = df.sort_values("rank_score", ascending=False)["id"].tolist()
    out = ra.rerank_output(df, "rankagg")
    assert weighted == ["a", "b", "c"]
    assert list(out["id"])[0] == "c"                 # c promoted to the top
    assert list(out["id"]) != weighted
