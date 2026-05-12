"""Unit tests for scripts/_rank_candidates_lib.py.

These cover the corrected versions of:
- categorize_reference  (bead -wux part 2)
- benjamini_hochberg    (bead -wux part 1)
- extract_branch_omega  (bead -ea9 part 2)
- get_selection_scores  (bead -ea9 part 1)
- calculate_fair_rank_score (bead -ce4)
- normalize_synteny_counts (bead -ce4 / -mqt)
"""
import math

import numpy as np
import pytest

from _rank_candidates_lib import (
    benjamini_hochberg,
    calculate_fair_rank_score,
    categorize_reference,
    extract_branch_omega,
    get_selection_scores,
    load_meme_concordance,
    normalize_synteny_counts,
)


# -----------------------------------------------------------------------------
# categorize_reference (bead -wux part 2)
# -----------------------------------------------------------------------------

class TestCategorizeReference:
    def test_olfactory_receptor_named_canonical(self):
        # Mammalian-style olfactory receptor names
        for name in ["Olfactory receptor 7A1", "OLFR1234", "Olfr12-ps", "OR1A1_HUMAN", "Or22a"]:
            assert categorize_reference(name, chemoreceptor_weight=2.0) == 2.0, name

    def test_vomeronasal_taste_taar(self):
        for name in ["Vomeronasal receptor V1R", "V2R-1", "Taste receptor type 2",
                     "TASR1", "T1R3", "T2R10", "TAAR5", "Trace amine-associated receptor 1"]:
            assert categorize_reference(name) == 2.0, name

    def test_non_chemoreceptor_or_substring_no_match(self):
        # The 'or' substring used to match these. Word-boundary regex must not.
        for name in [
            "Orexin receptor type 1",
            "Orphanin FQ receptor",
            "Origin recognition complex protein",
            "Adrenergic receptor alpha-2A",
            "Histamine receptor H1",
            "Dopamine receptor D2",
            "Cannabinoid receptor 1",
            "Glycoprotein hormone receptor",
        ]:
            assert categorize_reference(name) == 1.0, name

    def test_explicit_category_map_wins(self):
        m = {"weird_name_42": 7.5}
        assert categorize_reference("weird_name_42", explicit_category_map=m) == 7.5

    def test_explicit_chemoreceptor_set_second(self):
        s = {"some_lse_paralog_xyz"}
        # Has no keyword match but is explicitly in the set -> chemoreceptor weight
        assert categorize_reference("some_lse_paralog_xyz",
                                    explicit_chemoreceptor_set=s) == 2.0

    def test_empty_or_none_name(self):
        assert categorize_reference("") == 1.0
        assert categorize_reference(None) == 1.0  # type: ignore[arg-type]


# -----------------------------------------------------------------------------
# benjamini_hochberg (bead -wux part 1)
# -----------------------------------------------------------------------------

class TestBenjaminiHochberg:
    def test_empty_input(self):
        assert benjamini_hochberg([]) == []

    def test_against_statsmodels_reference(self):
        # Compare against statsmodels directly on a random-but-fixed input
        np.random.seed(0)
        pvals = np.random.uniform(0, 1, size=100).tolist()
        # Inject some small p-values so FDR has signal
        pvals[5] = 1e-5
        pvals[42] = 1e-3
        pvals[77] = 0.012
        from statsmodels.stats.multitest import multipletests
        _, expected_q, _, _ = multipletests(pvals, method="fdr_bh")
        out = benjamini_hochberg(pvals)
        np.testing.assert_allclose(out, expected_q, atol=1e-12)

    def test_monotone_q_in_sorted_p_order(self):
        # q-values, when sorted by p-value, must be non-decreasing
        pvals = [0.5, 0.01, 0.3, 0.001, 0.04, 0.9]
        q = benjamini_hochberg(pvals)
        sorted_pairs = sorted(zip(pvals, q), key=lambda x: x[0])
        qs = [x[1] for x in sorted_pairs]
        for a, b in zip(qs, qs[1:]):
            assert b >= a - 1e-9

    def test_nan_propagation(self):
        pvals = [0.01, float("nan"), 0.5, 0.0001, float("nan")]
        q = benjamini_hochberg(pvals)
        assert math.isnan(q[1])
        assert math.isnan(q[4])
        # The non-NaN entries should be FDR-corrected among themselves
        assert all(0 <= q[i] <= 1 for i in (0, 2, 3))

    def test_all_nan(self):
        q = benjamini_hochberg([float("nan"), float("nan")])
        assert all(math.isnan(x) for x in q)


# -----------------------------------------------------------------------------
# extract_branch_omega (bead -ea9 part 2)
# -----------------------------------------------------------------------------

class TestExtractBranchOmega:
    @staticmethod
    def _branch(rate_classes):
        return {"Rate Distributions": [list(rc) for rc in rate_classes]}

    def test_episodic_positive_selection_preserved(self):
        # 90% sites at omega=0.1, 10% sites at omega=20
        out = extract_branch_omega(self._branch([(0.1, 0.9), (20.0, 0.1)]))
        assert out["omega_max"] == pytest.approx(20.0)
        assert out["omega_mean"] == pytest.approx(0.1 * 0.9 + 20.0 * 0.1)
        assert out["weight_at_max"] == pytest.approx(0.1)
        assert out["n_rate_classes"] == 2

    def test_pure_purifying(self):
        out = extract_branch_omega(self._branch([(0.05, 1.0)]))
        assert out["omega_max"] == pytest.approx(0.05)
        assert out["omega_mean"] == pytest.approx(0.05)
        assert out["weight_at_max"] == pytest.approx(1.0)

    def test_empty_rate_classes(self):
        out = extract_branch_omega({"Rate Distributions": []})
        assert math.isnan(out["omega_max"])
        assert out["n_rate_classes"] == 0

    def test_missing_rate_distributions_key(self):
        out = extract_branch_omega({})
        assert math.isnan(out["omega_max"])

    def test_three_rate_classes(self):
        out = extract_branch_omega(self._branch([(0.01, 0.5), (1.0, 0.4), (15.0, 0.1)]))
        assert out["omega_max"] == 15.0
        assert out["weight_at_max"] == pytest.approx(0.1)
        assert out["omega_mean"] == pytest.approx(0.01 * 0.5 + 1.0 * 0.4 + 15.0 * 0.1)


# -----------------------------------------------------------------------------
# get_selection_scores (bead -ea9 part 1)
# -----------------------------------------------------------------------------

class TestGetSelectionScores:
    def test_purifying_weight_zero_does_not_score_purifying(self):
        # Default purifying_weight=0 (chemoreceptor-discovery mode)
        s = get_selection_scores(omega=0.05, p_corrected=0.01,
                                 purifying_weight=0.0, positive_weight=1.0)
        assert s["purifying_score"] == 0.0
        assert s["positive_score"] == 0.0  # omega < 1 -> no positive score

    def test_strong_positive_selection_rewarded_with_significance_boost(self):
        # omega=20, log10=1.301; p<0.05 -> 1.5x boost
        s = get_selection_scores(omega=20.0, p_corrected=0.001,
                                 purifying_weight=0.0, positive_weight=1.0)
        assert s["positive_score"] == pytest.approx(1.30103 * 1.5, rel=1e-3)
        assert s["is_significant"] is True

    def test_positive_without_significance(self):
        s = get_selection_scores(omega=20.0, p_corrected=0.5,
                                 purifying_weight=0.0, positive_weight=1.0)
        assert s["positive_score"] == pytest.approx(1.30103, rel=1e-3)
        assert s["is_significant"] is False

    def test_purifying_when_explicitly_enabled(self):
        # If user sets purifying_weight > 0, conserved-function GPCRs are scored
        s = get_selection_scores(omega=0.05, p_corrected=0.01,
                                 purifying_weight=1.0, positive_weight=1.0)
        # |log10(0.05)| = 1.301; with 1.5x significance boost
        assert s["purifying_score"] == pytest.approx(1.30103 * 1.5, rel=1e-3)
        assert s["positive_score"] == 0.0

    def test_neutral_omega_one(self):
        s = get_selection_scores(omega=1.0, p_corrected=0.5,
                                 purifying_weight=1.0, positive_weight=1.0)
        assert s["purifying_score"] == 0.0
        assert s["positive_score"] == 0.0

    def test_nan_omega(self):
        s = get_selection_scores(omega=float("nan"), p_corrected=float("nan"),
                                 purifying_weight=0.0, positive_weight=1.0)
        assert s["purifying_score"] == 0.0
        assert s["positive_score"] == 0.0
        assert s["is_significant"] is False


# -----------------------------------------------------------------------------
# calculate_fair_rank_score (bead -ce4)
# -----------------------------------------------------------------------------

class TestCalculateFairRankScore:
    def test_full_evidence_outranks_sparse_at_same_per_axis(self):
        weights = {"phylo": 2, "expr": 1, "dnds": 1}
        sparse = calculate_fair_rank_score({"phylo": 1.0}, weights)
        full = calculate_fair_rank_score({"phylo": 1.0, "expr": 1.0, "dnds": 1.0}, weights)
        assert full > sparse

    def test_evidence_completeness_multiplier_floor(self):
        # 1 of 3 axes available -> raw completeness 0.5 (phylo weight 2 / total 4)
        # but never below floor 0.4
        weights = {"phylo": 2, "expr": 1, "dnds": 1}
        out = calculate_fair_rank_score({"phylo": 1.0}, weights, return_diagnostics=True)
        assert out["evidence_completeness"] >= 0.4
        assert out["evidence_completeness_raw"] == pytest.approx(0.5)

    def test_no_data_floors_at_zero(self):
        weights = {"phylo": 2}
        assert calculate_fair_rank_score({}, weights) == 0.0

    def test_nan_score_treated_as_missing(self):
        weights = {"phylo": 1, "expr": 1}
        out = calculate_fair_rank_score(
            {"phylo": 1.0, "expr": float("nan")}, weights, return_diagnostics=True
        )
        # Only phylo contributed; weight_avail = 1, total = 2 -> raw 0.5
        assert out["available_weight"] == pytest.approx(1.0)
        assert out["evidence_completeness_raw"] == pytest.approx(0.5)

    def test_full_data_completeness_is_one(self):
        weights = {"a": 1, "b": 2, "c": 3}
        scores = {"a": 0.5, "b": 0.5, "c": 0.5}
        out = calculate_fair_rank_score(scores, weights, return_diagnostics=True)
        assert out["evidence_completeness"] == pytest.approx(1.0)
        assert out["score"] == pytest.approx(0.5)

    def test_zero_weights_safe(self):
        assert calculate_fair_rank_score({"phylo": 1.0}, {"phylo": 0}) == 0.0
        assert calculate_fair_rank_score({}, {}) == 0.0


# -----------------------------------------------------------------------------
# normalize_synteny_counts (bead -ce4 part 2 / -mqt)
# -----------------------------------------------------------------------------

class TestNormalizeSyntenyCounts:
    def test_skip_when_max_below_threshold(self):
        out = normalize_synteny_counts({"a": 0, "b": 1, "c": 2}, min_max_anchors=5)
        assert all(v is None for v in out.values())

    def test_log_scale_when_real_data(self):
        out = normalize_synteny_counts({"a": 1, "b": 5, "c": 50}, min_max_anchors=5)
        assert out["c"] == pytest.approx(1.0)
        assert 0.0 < out["b"] < 1.0
        assert 0.0 < out["a"] < out["b"]

    def test_zero_counts_score_zero(self):
        out = normalize_synteny_counts({"a": 0, "b": 10}, min_max_anchors=5)
        assert out["a"] == 0.0
        assert out["b"] == pytest.approx(1.0)

    def test_empty_input(self):
        assert normalize_synteny_counts({}) == {}


# -----------------------------------------------------------------------------
# load_meme_concordance (bead -7cy step 2)
# -----------------------------------------------------------------------------

class TestLoadMemeConcordance:
    """parse_meme.py's dual-mode writes a per-OG concordance CSV with
    high_confidence/lenient_only/strict_only counts + robustness index.
    load_meme_concordance reads it for rank_candidates.py."""

    def _write_csv(self, path, rows):
        import csv
        fields = ["og_name", "n_strict_positive_sites",
                  "n_lenient_positive_sites", "high_confidence_sites_n",
                  "lenient_only_sites_n", "strict_only_sites_n",
                  "alignment_robustness_index"]
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)

    def test_reads_per_og_concordance_fields(self, tmp_path):
        f = tmp_path / "meme_concordance.csv"
        self._write_csv(f, [
            {"og_name": "OG001", "n_strict_positive_sites": 2,
             "n_lenient_positive_sites": 3, "high_confidence_sites_n": 2,
             "lenient_only_sites_n": 1, "strict_only_sites_n": 0,
             "alignment_robustness_index": 0.6667},
            {"og_name": "OG002", "n_strict_positive_sites": 0,
             "n_lenient_positive_sites": 4, "high_confidence_sites_n": 0,
             "lenient_only_sites_n": 4, "strict_only_sites_n": 0,
             "alignment_robustness_index": 0.0},
        ])
        out = load_meme_concordance(str(f))
        assert set(out.keys()) == {"OG001", "OG002"}
        assert out["OG001"]["high_confidence_sites_n"] == 2
        assert out["OG001"]["lenient_only_sites_n"] == 1
        assert out["OG001"]["alignment_robustness_index"] == pytest.approx(0.6667)
        # The all-alignment-sensitive case: lenient finds signal, strict misses entirely
        assert out["OG002"]["high_confidence_sites_n"] == 0
        assert out["OG002"]["alignment_robustness_index"] == 0.0

    def test_missing_file_returns_empty_dict(self, tmp_path):
        out = load_meme_concordance(str(tmp_path / "absent.csv"))
        assert out == {}

    def test_skips_rows_with_blank_og_name(self, tmp_path):
        f = tmp_path / "c.csv"
        self._write_csv(f, [
            {"og_name": "", "n_strict_positive_sites": 1,
             "n_lenient_positive_sites": 1, "high_confidence_sites_n": 1,
             "lenient_only_sites_n": 0, "strict_only_sites_n": 0,
             "alignment_robustness_index": 1.0},
            {"og_name": "OGreal", "n_strict_positive_sites": 1,
             "n_lenient_positive_sites": 1, "high_confidence_sites_n": 1,
             "lenient_only_sites_n": 0, "strict_only_sites_n": 0,
             "alignment_robustness_index": 1.0},
        ])
        out = load_meme_concordance(str(f))
        assert set(out.keys()) == {"OGreal"}

    def test_robustness_index_coerced_to_float(self, tmp_path):
        f = tmp_path / "c.csv"
        self._write_csv(f, [{
            "og_name": "OG", "n_strict_positive_sites": 1,
            "n_lenient_positive_sites": 2, "high_confidence_sites_n": 1,
            "lenient_only_sites_n": 1, "strict_only_sites_n": 0,
            "alignment_robustness_index": "0.5",   # stringified
        }])
        out = load_meme_concordance(str(f))
        assert isinstance(out["OG"]["alignment_robustness_index"], float)
        assert out["OG"]["alignment_robustness_index"] == 0.5
