"""Tests for scripts/ablate_ranking.py — the ranking ablation harness.

A council review of the hand-weighted chemoreceptor ranking asked two
empirical questions: (a) is the ranking actually driven by which signals a
candidate happens to have data for (the has_*_data / evidence_completeness
machinery) rather than by what those signals say, and (b) do a few
correlated signal groups dominate the composite regardless of the rest?

The harness must measure the REAL ranking — the ``rank_score`` column
rank_candidates.py writes and sorts on (line ~2002/2047), which is produced
by the row-wise ``calculate_fair_rank_score`` (rank_candidates.py:1953) →
``_rank_candidates_lib.calculate_fair_rank_score(scores, weights,
completeness_floor=0.4)``. That scorer has the evidence-completeness
multiplier (floored at 0.4 — literally the missingness effect the ablation
exists to measure) and a 12th ``tandem_cluster`` signal. These tests pin
the harness against that scorer, with hand-computed literals derived from
the fair formula ``score = (weighted_sum / avail_weight) * max(floor,
avail_weight / total_weight)``.
"""
from __future__ import annotations

import pandas as pd
import pytest

import _rank_candidates_lib
import ablate_ranking as ar

# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------

# The full 12-signal production weights (rank_candidates.py's *_WEIGHT env
# defaults). NOTE the fair scorer's key is 'expression' (fed by EXPR_WEIGHT),
# not 'expr', and includes tandem_cluster (2.5). Sum = 20.5. Passed
# explicitly wherever a test needs exact, hand-computable numbers.
DEFAULT_WEIGHTS = {
    "phylo": 2.0,
    "purifying": 1.0,
    "positive": 1.0,
    "lse_depth": 1.0,
    "synteny": 3.0,
    "expression": 1.0,
    "chemosensory_expr": 3.0,
    "gprotein_coexpr": 2.0,
    "ecl_divergence": 1.5,
    "expansion": 1.5,
    "og_confidence": 1.0,
    "tandem_cluster": 2.5,
}  # sum = 20.5
TOTAL_WEIGHT = 20.5

WEIGHT_ENV_VARS = [
    "PHYLO_WEIGHT", "PURIFYING_WEIGHT", "POSITIVE_WEIGHT", "SYNTENY_WEIGHT",
    "EXPR_WEIGHT", "LSE_DEPTH_WEIGHT", "CHEMOSENSORY_EXPR_WEIGHT",
    "GPROTEIN_COEXPR_WEIGHT", "ECL_DIVERGENCE_WEIGHT", "EXPANSION_WEIGHT",
    "OG_CONFIDENCE_WEIGHT", "TANDEM_CLUSTER_WEIGHT",
]


def _full_row(id_, **overrides):
    """A candidate row with every column the fair scorer's bridge touches,
    defaulting to "no data" (has_*=False, scores=0.0). Mirrors how
    rank_candidates.py always populates these columns for every row
    (missing raw scores normalize to 0.0; has_*_data is what actually
    gates inclusion) — see rank_candidates.py's normalize_cols loop and
    calculate_fair_rank_score at lines 1970-1998.
    """
    row = {
        "id": id_,
        "phylo_score_norm": 0.0,
        "purifying_score_norm": 0.0,
        "positive_score_norm": 0.0,
        "lse_depth_score_norm": 0.0,
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
        "tandem_cluster_score_norm": 0.0,
        "has_tandem_cluster_data": False,
        "evidence_completeness": 0.0,
    }
    row.update(overrides)
    return row


def _direct_fair_score(row, weights):
    """Independent ground-truth: build the scores dict inline (an explicit
    replication of rank_candidates.py:1970-1984) and call the production
    lib scorer directly. Deliberately does NOT reuse ar._row_scores, so a
    bug in the harness's own bridge is caught here.
    """
    scores = {
        "phylo": row.get("phylo_score_norm"),
        "purifying": row.get("purifying_score_norm"),
        "positive": row.get("positive_score_norm"),
        "lse_depth": row.get("lse_depth_score_norm"),
        "synteny": row.get("synteny_score_norm") if row.get("has_synteny_data") else None,
        "expression": row.get("expression_score_norm") if row.get("has_expression_data") else None,
        "chemosensory_expr": row.get("chemosensory_expr_score_norm") if row.get("has_chemosensory_expr_data") else None,
        "gprotein_coexpr": row.get("gprotein_coexpr_score_norm") if row.get("has_gprotein_data") else None,
        "ecl_divergence": row.get("ecl_divergence_score_norm") if row.get("has_ecl_data") else None,
        "expansion": row.get("expansion_score_norm") if row.get("has_expansion_data") else None,
        "og_confidence": row.get("og_confidence_score_norm") if row.get("has_og_confidence_data") else None,
        "tandem_cluster": row.get("tandem_cluster_score_norm") if row.get("has_tandem_cluster_data") else None,
    }
    return _rank_candidates_lib.calculate_fair_rank_score(scores, weights, completeness_floor=0.4)


@pytest.fixture
def sparse_vs_dense_df():
    """Two candidates whose baseline order is set by the completeness floor:

    cand_sparse: only phylo_score_norm=1.0, no other data.
      avail_weight = 5 (the 4 always-present base axes); weighted_sum = 2.0
      (phylo weight 2 * 1.0). completeness_raw = 5/20.5 = 0.244, BELOW the
      0.4 floor -> completeness = 0.4. score = (2.0/5)*0.4 = 0.16.

    cand_mid: moderate signal on synteny/chemo/gprotein (has_*=True).
      avail_weight = 5 + 3 + 3 + 2 = 13; weighted_sum = 3*0.25 + 3*0.25 +
      2*0.4 = 2.3. completeness_raw = 13/20.5 = 0.634 (ABOVE floor).
      score = (2.3/13)*(13/20.5) = 2.3/20.5 = 0.11220.

    Baseline: sparse (0.16) > mid (0.1122). The sparse candidate wins ONLY
    because the floor protects its low completeness — exactly the
    missingness effect hypothesis (a) probes. Neutralizing missingness
    (removing the floor's protection: sparse -> 2.0/20.5 = 0.0976) or
    dropping 'phylo' (removing sparse's only signal -> 0.0) both flip it.
    """
    sparse = _full_row("cand_sparse", phylo_score_norm=1.0)
    mid = _full_row(
        "cand_mid",
        synteny_score_norm=0.25, has_synteny_data=True,
        chemosensory_expr_score_norm=0.25, has_chemosensory_expr_data=True,
        gprotein_coexpr_score_norm=0.4, has_gprotein_data=True,
    )
    return pd.DataFrame([sparse, mid])


@pytest.fixture
def tandem_vs_plain_df():
    """Isolates the 12th signal (tandem_cluster). cand_tandem leans on
    tandem_cluster; cand_plain has none. Dropping tandem_cluster's weight
    must flip their order — proving the tandem axis is wired into the
    harness.

    cand_tandem: phylo=0.5, tandem_cluster=1.0 (has_tandem_cluster_data=True).
      avail = 5 + 2.5 = 7.5; weighted_sum = 2*0.5 + 2.5*1.0 = 3.5.
      completeness_raw = 7.5/20.5 = 0.366 (below floor) -> 0.4.
      score = (3.5/7.5)*0.4 = 0.18667.
    cand_plain: phylo=0.6 only. avail = 5; weighted_sum = 1.2.
      completeness -> 0.4. score = (1.2/5)*0.4 = 0.096.
    Baseline: [tandem, plain]. Drop tandem_cluster -> tandem 0.08 < plain
    0.096 -> [plain, tandem].
    """
    tandem = _full_row(
        "cand_tandem", phylo_score_norm=0.5,
        tandem_cluster_score_norm=1.0, has_tandem_cluster_data=True,
    )
    plain = _full_row("cand_plain", phylo_score_norm=0.6)
    return pd.DataFrame([tandem, plain])


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
# Targets the REAL production fair scorer (not calculate_rank_score)
# ---------------------------------------------------------------------------


def test_harness_reuses_lib_fair_rank_score_object():
    """The harness must import the production lib scorer itself — not a copy,
    mirror, or the sensitivity-only calculate_rank_score."""
    assert ar.calculate_fair_rank_score is _rank_candidates_lib.calculate_fair_rank_score


def test_harness_dropped_the_wrong_scorer_apparatus():
    """The earlier revision targeted calculate_rank_score via source-slicing
    + a hand mirror. That must be gone: no extract/exec or local-mirror
    fallback remains, and the module does not carry its own
    calculate_rank_score."""
    assert not hasattr(ar, "_extract_calculate_rank_score_from_source")
    assert not hasattr(ar, "_local_calculate_rank_score")
    assert not hasattr(ar, "calculate_rank_score")


# ---------------------------------------------------------------------------
# _production_weights
# ---------------------------------------------------------------------------


def test_production_weights_are_the_full_12_signal_set(monkeypatch):
    for var in WEIGHT_ENV_VARS:
        monkeypatch.delenv(var, raising=False)
    weights = ar._production_weights()
    assert weights == DEFAULT_WEIGHTS
    # Explicitly pin the two properties the rework hinged on.
    assert "tandem_cluster" in weights and weights["tandem_cluster"] == 2.5
    assert "expression" in weights and "expr" not in weights


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
# Ground truth: ablate() reproduces the real fair rank_score ordering
# ---------------------------------------------------------------------------


def test_ablate_reproduces_direct_fair_score_ordering(sparse_vs_dense_df):
    """ablate(df) with no ablation must equal df sorted by the real fair
    rank_score obtained by calling _rank_candidates_lib.calculate_fair_rank_score
    directly per row."""
    direct = [
        (r["id"], _direct_fair_score(r, DEFAULT_WEIGHTS))
        for _, r in sparse_vs_dense_df.iterrows()
    ]
    expected_order = [cid for cid, _ in sorted(direct, key=lambda t: t[1], reverse=True)]
    assert ar.ablate(sparse_vs_dense_df, weights=DEFAULT_WEIGHTS) == expected_order


def test_ablate_baseline_matches_hand_computed_fair_scores(sparse_vs_dense_df):
    scores = {
        r["id"]: _direct_fair_score(r, DEFAULT_WEIGHTS)
        for _, r in sparse_vs_dense_df.iterrows()
    }
    assert scores["cand_sparse"] == pytest.approx(0.16)          # floor-protected
    assert scores["cand_mid"] == pytest.approx(2.3 / TOTAL_WEIGHT)  # 0.11220
    assert ar.ablate(sparse_vs_dense_df, weights=DEFAULT_WEIGHTS) == ["cand_sparse", "cand_mid"]


# ---------------------------------------------------------------------------
# ablate() — ablations reorder as the fair formula predicts
# ---------------------------------------------------------------------------


def test_ablate_drop_dominant_signal_reorders_topk(sparse_vs_dense_df):
    baseline = ar.ablate(sparse_vs_dense_df, weights=DEFAULT_WEIGHTS)
    dropped = ar.ablate(sparse_vs_dense_df, drop_signals=["phylo"], weights=DEFAULT_WEIGHTS)
    assert baseline == ["cand_sparse", "cand_mid"]
    assert dropped == ["cand_mid", "cand_sparse"]
    assert ar.topk_jaccard(baseline, dropped, k=1) == 0.0


def test_ablate_neutralize_missingness_changes_order_when_completeness_differs(sparse_vs_dense_df):
    """Proves the evidence-completeness multiplier is in the loop: forcing
    every axis 'present' removes the floor bonus the sparse candidate
    relied on, flipping the order."""
    baseline = ar.ablate(sparse_vs_dense_df, weights=DEFAULT_WEIGHTS)
    neutralized = ar.ablate(
        sparse_vs_dense_df, neutralize_missingness=True, weights=DEFAULT_WEIGHTS
    )
    assert baseline == ["cand_sparse", "cand_mid"]
    assert neutralized == ["cand_mid", "cand_sparse"]
    assert ar.topk_jaccard(baseline, neutralized, k=1) == 0.0


def test_ablate_drop_tandem_cluster_reorders(tandem_vs_plain_df):
    """The 12th signal must be wired: a tandem-led candidate loses to a
    plain one once tandem_cluster's weight is zeroed."""
    baseline = ar.ablate(tandem_vs_plain_df, weights=DEFAULT_WEIGHTS)
    dropped = ar.ablate(tandem_vs_plain_df, drop_signals=["tandem_cluster"], weights=DEFAULT_WEIGHTS)
    assert baseline == ["cand_tandem", "cand_plain"]
    assert dropped == ["cand_plain", "cand_tandem"]


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
    assert ranked == ["cand_sparse", "cand_mid"]


def test_ablate_drop_signals_applies_on_top_of_custom_weights():
    df = pd.DataFrame([
        _full_row("g1", phylo_score_norm=1.0),
        _full_row("g2", lse_depth_score_norm=1.0),
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
    report = ar.ablation_report(sparse_vs_dense_df, groups=[["phylo"], ["synteny", "expression"]], k=1)
    assert set(report.keys()) == {"missingness", "group:phylo", "group:synteny+expression"}
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


def test_ablation_report_order_stable_group_shows_no_churn():
    # Two candidates ordered by phylo alone; dropping og_confidence (which
    # neither has and which is below both floors) leaves the ORDER unchanged,
    # so the churn metrics report no movement.
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
