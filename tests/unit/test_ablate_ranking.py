"""Tests for scripts/ablate_ranking.py — the ranking ablation harness.

A council review of the hand-weighted chemoreceptor ranking asked two
empirical questions: (a) is the ranking actually driven by which signals a
candidate happens to have data for (the has_*_data / evidence_completeness
machinery) rather than by what those signals say, and (b) do a few
correlated signal groups dominate the composite regardless of the rest?
These tests pin the harness's primitives (topk_jaccard, ablate,
ablation_report, write_markdown) against hand-computed expectations using
rank_candidates.calculate_rank_score's actual per-row formula, so a
regression in either the harness or its reuse of the production scorer
would be caught here.
"""
from __future__ import annotations

import pandas as pd
import pytest

import ablate_ranking as ar

# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------

# Explicit weights matching rank_candidates.py's literal *_WEIGHT env-var
# defaults (see ar._production_weights). Tests that need exact,
# hand-computable numbers pass this explicitly rather than relying on
# _production_weights()'s environment, which could vary by shell.
DEFAULT_WEIGHTS = {
    "phylo": 2.0,
    "purifying": 1.0,
    "positive": 1.0,
    "synteny": 3.0,
    "expr": 1.0,
    "lse_depth": 1.0,
    "chemosensory_expr": 3.0,
    "gprotein_coexpr": 2.0,
    "ecl_divergence": 1.5,
    "expansion": 1.5,
    "og_confidence": 1.0,
}  # sum = 18.0

WEIGHT_ENV_VARS = [
    "PHYLO_WEIGHT", "PURIFYING_WEIGHT", "POSITIVE_WEIGHT", "SYNTENY_WEIGHT",
    "EXPR_WEIGHT", "LSE_DEPTH_WEIGHT", "CHEMOSENSORY_EXPR_WEIGHT",
    "GPROTEIN_COEXPR_WEIGHT", "ECL_DIVERGENCE_WEIGHT", "EXPANSION_WEIGHT",
    "OG_CONFIDENCE_WEIGHT",
]


def _full_row(id_, **overrides):
    """A candidate row with every column calculate_rank_score touches,
    defaulting to "no data" (has_*=False, scores=0.0). Mirrors how
    rank_candidates.py always populates these columns for every row
    (missing raw scores normalize to 0.0; has_*_data is what actually
    gates inclusion) — see rank_candidates.py's normalize_cols loop.
    """
    row = {
        "id": id_,
        "phylo_score_norm": 0.0,
        "purifying_score_norm": 0.0,
        "positive_score_norm": 0.0,
        "lse_depth_score_norm": 0.0,
        "dnds_reliability_weight": 1.0,
        "synteny_score_norm": 0.0,
        "has_synteny_data": False,
        "expression_score_norm": 0.0,
        "has_expression_data": False,
        "chemosensory_expr_score_norm": 0.0,
        "has_chemosensory_expr_data": False,
        "gprotein_coexpr_score_norm": 0.0,
        "has_gprotein_data": False,
        "ecl_divergence_score_norm": 0.0,
        "has_ecl_data": False,
        "expansion_score_norm": 0.0,
        "has_expansion_data": False,
        "og_confidence_score_norm": 0.0,
        "has_og_confidence_data": False,
        "evidence_completeness": 0.0,
    }
    row.update(overrides)
    return row


@pytest.fixture
def sparse_vs_dense_df():
    """Two candidates with the SAME hand-computed baseline score gap used
    throughout this file:

    cand_sparse: only phylo_score_norm=1.0, no other data at all.
      score = 1.0*2 = 2.0; total_weight = 2+1+1+1 = 5 (base axes only)
      -> (2.0/5) * 18 = 7.2

    cand_dense: moderate (0.2-0.3) score on every axis, ALL has_*=True.
      score sums to 6.5 with total_weight == max_possible_weight (18)
      -> 6.5

    So the DEFAULT ranking (7.2 > 6.5) puts the sparse-but-strong-on-one-axis
    candidate FIRST — exactly the "fair scoring" behavior hypothesis (a)
    is concerned about. Dropping 'phylo' (removing cand_sparse's only
    signal) or neutralizing missingness (diluting cand_sparse's score by
    counting axes it has no data for) both flip the order.
    """
    sparse = _full_row("cand_sparse", phylo_score_norm=1.0)
    dense = _full_row(
        "cand_dense",
        phylo_score_norm=0.3, purifying_score_norm=0.2, positive_score_norm=0.2,
        lse_depth_score_norm=0.3,
        synteny_score_norm=0.4, has_synteny_data=True,
        expression_score_norm=0.4, has_expression_data=True,
        chemosensory_expr_score_norm=0.4, has_chemosensory_expr_data=True,
        gprotein_coexpr_score_norm=0.4, has_gprotein_data=True,
        ecl_divergence_score_norm=0.4, has_ecl_data=True,
        expansion_score_norm=0.4, has_expansion_data=True,
        og_confidence_score_norm=0.4, has_og_confidence_data=True,
        evidence_completeness=1.0,
    )
    return pd.DataFrame([sparse, dense])


# ---------------------------------------------------------------------------
# topk_jaccard
# ---------------------------------------------------------------------------


def test_topk_jaccard_identical_lists_is_one():
    a = ["g1", "g2", "g3", "g4"]
    assert ar.topk_jaccard(a, a, 4) == 1.0
    assert ar.topk_jaccard(a, a, 2) == 1.0


def test_topk_jaccard_disjoint_lists_is_zero():
    a = ["a", "b", "c"]
    b = ["x", "y", "z"]
    assert ar.topk_jaccard(a, b, 3) == 0.0


def test_topk_jaccard_partial_overlap_known_value():
    a = ["A", "B", "C", "D"]
    b = ["A", "B", "E", "F"]
    # top-3: {A,B,C} vs {A,B,E} -> intersection {A,B}=2, union {A,B,C,E}=4
    assert ar.topk_jaccard(a, b, 3) == pytest.approx(0.5)


def test_topk_jaccard_empty_union_returns_one():
    assert ar.topk_jaccard([], [], 5) == 1.0
    assert ar.topk_jaccard(["a", "b"], ["c", "d"], 0) == 1.0


def test_topk_jaccard_k_larger_than_lists_still_works():
    a = ["a", "b"]
    b = ["b", "a"]
    assert ar.topk_jaccard(a, b, 100) == 1.0


# ---------------------------------------------------------------------------
# Production scorer resolution
# ---------------------------------------------------------------------------


def test_resolve_calculate_rank_score_returns_memoized_callable():
    scorer1 = ar._resolve_calculate_rank_score()
    scorer2 = ar._resolve_calculate_rank_score()
    assert callable(scorer1)
    assert scorer1 is scorer2


def test_extract_calculate_rank_score_from_source_matches_hand_computed_value():
    fn = ar._extract_calculate_rank_score_from_source()
    df = pd.DataFrame([_full_row("g1", phylo_score_norm=1.0)])
    result = fn(df, DEFAULT_WEIGHTS)
    assert result.tolist() == pytest.approx([7.2])


def test_local_mirror_matches_hand_computed_values(sparse_vs_dense_df):
    result = ar._local_calculate_rank_score(sparse_vs_dense_df, DEFAULT_WEIGHTS)
    assert result.tolist() == pytest.approx([7.2, 6.5])


def test_resolved_scorer_matches_hand_computed_values_for_all_three_scenarios(sparse_vs_dense_df):
    """Pins the harness's actual reused scorer (whichever tier resolves)
    against hand-computed values for baseline, drop-phylo, and
    neutralize-missingness — the three scenarios exercised elsewhere in
    this file via the public `ablate` API."""
    scorer = ar._resolve_calculate_rank_score()

    baseline = scorer(sparse_vs_dense_df, DEFAULT_WEIGHTS)
    assert baseline.tolist() == pytest.approx([7.2, 6.5])

    dropped_weights = ar._apply_drop_signals(DEFAULT_WEIGHTS, ["phylo"])
    dropped = scorer(sparse_vs_dense_df, dropped_weights)
    assert dropped.tolist() == pytest.approx([0.0, 5.9])

    neutralized_df = ar._neutralize_missingness(sparse_vs_dense_df)
    neutralized = scorer(neutralized_df, DEFAULT_WEIGHTS)
    assert neutralized.tolist() == pytest.approx([2.0, 6.5])


# ---------------------------------------------------------------------------
# _production_weights
# ---------------------------------------------------------------------------


def test_production_weights_match_rank_candidates_defaults(monkeypatch):
    for var in WEIGHT_ENV_VARS:
        monkeypatch.delenv(var, raising=False)
    assert ar._production_weights() == DEFAULT_WEIGHTS


# ---------------------------------------------------------------------------
# _apply_drop_signals / _neutralize_missingness (isolated units)
# ---------------------------------------------------------------------------


def test_apply_drop_signals_zeroes_only_named_keys():
    weights = {"phylo": 2.0, "synteny": 3.0}
    dropped = ar._apply_drop_signals(weights, ["phylo"])
    assert dropped == {"phylo": 0.0, "synteny": 3.0}
    assert weights == {"phylo": 2.0, "synteny": 3.0}  # input untouched


def test_apply_drop_signals_none_returns_a_copy():
    weights = {"phylo": 2.0}
    result = ar._apply_drop_signals(weights, None)
    assert result == weights
    assert result is not weights


def test_neutralize_missingness_forces_has_columns_true_and_completeness_one(sparse_vs_dense_df):
    out = ar._neutralize_missingness(sparse_vs_dense_df)
    has_cols = [c for c in out.columns if c.startswith("has_")]
    assert has_cols, "fixture should have has_*_data columns"
    for c in has_cols:
        assert out[c].all()
    assert (out["evidence_completeness"] == 1.0).all()


def test_neutralize_missingness_does_not_mutate_input(sparse_vs_dense_df):
    before = sparse_vs_dense_df.copy(deep=True)
    ar._neutralize_missingness(sparse_vs_dense_df)
    pd.testing.assert_frame_equal(sparse_vs_dense_df, before)


# ---------------------------------------------------------------------------
# ablate()
# ---------------------------------------------------------------------------


def test_ablate_baseline_ranks_sparse_strong_signal_first(sparse_vs_dense_df):
    ranked = ar.ablate(sparse_vs_dense_df, weights=DEFAULT_WEIGHTS)
    assert ranked == ["cand_sparse", "cand_dense"]


def test_ablate_drop_dominant_signal_reorders_topk(sparse_vs_dense_df):
    baseline = ar.ablate(sparse_vs_dense_df, weights=DEFAULT_WEIGHTS)
    dropped = ar.ablate(sparse_vs_dense_df, drop_signals=["phylo"], weights=DEFAULT_WEIGHTS)
    assert dropped == ["cand_dense", "cand_sparse"]
    assert ar.topk_jaccard(baseline, dropped, k=1) == 0.0


def test_ablate_neutralize_missingness_changes_order_when_completeness_differs(sparse_vs_dense_df):
    baseline = ar.ablate(sparse_vs_dense_df, weights=DEFAULT_WEIGHTS)
    neutralized = ar.ablate(sparse_vs_dense_df, neutralize_missingness=True, weights=DEFAULT_WEIGHTS)
    assert baseline == ["cand_sparse", "cand_dense"]
    assert neutralized == ["cand_dense", "cand_sparse"]
    assert ar.topk_jaccard(baseline, neutralized, k=1) == 0.0


def test_ablate_does_not_mutate_input(sparse_vs_dense_df):
    before = sparse_vs_dense_df.copy(deep=True)
    ar.ablate(
        sparse_vs_dense_df, drop_signals=["phylo"], neutralize_missingness=True,
        weights=DEFAULT_WEIGHTS,
    )
    pd.testing.assert_frame_equal(sparse_vs_dense_df, before)


def test_ablate_uses_production_defaults_when_weights_not_given(sparse_vs_dense_df, monkeypatch):
    for var in WEIGHT_ENV_VARS:
        monkeypatch.delenv(var, raising=False)
    ranked = ar.ablate(sparse_vs_dense_df)
    assert ranked == ["cand_sparse", "cand_dense"]


def test_ablate_drop_signals_applies_on_top_of_custom_weights():
    df = pd.DataFrame([
        _full_row("g1", phylo_score_norm=1.0, lse_depth_score_norm=0.0),
        _full_row("g2", phylo_score_norm=0.0, lse_depth_score_norm=1.0),
    ])
    custom_weights = {"phylo": 5.0, "lse_depth": 1.0}
    # g1 wins on phylo under custom weights...
    assert ar.ablate(df, weights=custom_weights) == ["g1", "g2"]
    # ...but dropping 'phylo' flips it to g2.
    assert ar.ablate(df, drop_signals=["phylo"], weights=custom_weights) == ["g2", "g1"]


# ---------------------------------------------------------------------------
# ablation_report()
# ---------------------------------------------------------------------------


def test_ablation_report_has_missingness_and_group_keys(sparse_vs_dense_df):
    report = ar.ablation_report(sparse_vs_dense_df, groups=[["phylo"], ["synteny", "expr"]], k=1)
    assert set(report.keys()) == {"missingness", "group:phylo", "group:synteny+expr"}
    for metrics in report.values():
        assert set(metrics.keys()) == {"spearman", "jaccard_at_k", "n_moved_into_topk"}


def test_ablation_report_supports_dict_groups_input(sparse_vs_dense_df):
    report = ar.ablation_report(
        sparse_vs_dense_df, groups={"phylo_axis": ["phylo", "og_confidence"]}, k=1
    )
    assert "group:phylo_axis" in report
    assert "missingness" in report


def test_ablation_report_missingness_reflects_hand_computed_flip(sparse_vs_dense_df):
    report = ar.ablation_report(sparse_vs_dense_df, groups=[], k=1)
    assert report["missingness"]["jaccard_at_k"] == 0.0
    assert report["missingness"]["n_moved_into_topk"] == 1
    assert report["missingness"]["spearman"] == pytest.approx(-1.0)


def test_ablation_report_no_op_group_matches_baseline_exactly():
    # Dropping an axis neither candidate has data for is a true no-op:
    # ranking, top-k membership, and rank order must all be unchanged.
    df = pd.DataFrame([
        _full_row("g1", phylo_score_norm=0.9),
        _full_row("g2", phylo_score_norm=0.1),
    ])
    report = ar.ablation_report(df, groups=[["og_confidence"]], k=2)
    assert report["group:og_confidence"]["spearman"] == pytest.approx(1.0)
    assert report["group:og_confidence"]["jaccard_at_k"] == 1.0
    assert report["group:og_confidence"]["n_moved_into_topk"] == 0


# ---------------------------------------------------------------------------
# write_markdown()
# ---------------------------------------------------------------------------


def test_write_markdown_creates_readable_table(tmp_path, sparse_vs_dense_df):
    report = ar.ablation_report(sparse_vs_dense_df, groups=[["phylo"]], k=1)
    out_path = tmp_path / "ablation_report.md"
    ar.write_markdown(report, out_path)

    assert out_path.exists()
    content = out_path.read_text()
    assert "missingness" in content
    assert "group:phylo" in content
    assert "Jaccard" in content
    assert "Spearman" in content
    # Formatted numeric cells (3 decimal places) actually landed in the file.
    assert "0.000" in content or "1.000" in content or "-1.000" in content


def test_write_markdown_handles_nan_spearman_gracefully(tmp_path):
    report = {
        "missingness": {"spearman": float("nan"), "jaccard_at_k": 1.0, "n_moved_into_topk": 0},
    }
    out_path = tmp_path / "nan_report.md"
    ar.write_markdown(report, out_path)
    content = out_path.read_text()
    assert "nan" in content.lower()
