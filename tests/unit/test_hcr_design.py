"""Unit tests for scripts/_hcr_design_lib.py (bead -xqz)."""
import pytest

from _hcr_design_lib import (
    closest_paralog_identity,
    compute_pairwise_identity,
    is_hcr_probe_friendly,
    windowed_max_identity,
)


class TestComputePairwiseIdentity:
    def test_identical(self):
        assert compute_pairwise_identity("ACGT", "ACGT") == 1.0

    def test_completely_different(self):
        assert compute_pairwise_identity("AAAA", "TTTT") == 0.0

    def test_half_match(self):
        assert compute_pairwise_identity("ACGT", "ACAA") == 0.5

    def test_gaps_excluded(self):
        # Aligned columns where neither side is gap: position 0 (A/A, match)
        # and position 2 (C/C, match). 2/2 = 1.0
        assert compute_pairwise_identity("A-C-", "ACCT") == 1.0

    def test_length_mismatch_raises(self):
        with pytest.raises(ValueError):
            compute_pairwise_identity("ACG", "ACGT")

    def test_all_gaps_returns_zero(self):
        assert compute_pairwise_identity("----", "ACGT") == 0.0


class TestWindowedMaxIdentity:
    def test_short_sequence_uses_full_overlap(self):
        assert windowed_max_identity("ACGT", "ACGT", window_aa=300) == 1.0

    def test_long_sequences_picks_best_window(self):
        a = "A" * 300 + "G" * 300
        b = "T" * 300 + "G" * 300
        assert windowed_max_identity(a, b, window_aa=300, step_aa=50) == 1.0

    def test_uniformly_divergent(self):
        a = "A" * 600
        b = "T" * 600
        assert windowed_max_identity(a, b, window_aa=300) == 0.0


class TestClosestParalogIdentity:
    def test_finds_closest(self):
        seqs = {
            "target": "AAAA" * 100,
            "very_similar": "AAAA" * 100,
            "different": "TTTT" * 100,
        }
        closest, ident = closest_paralog_identity("target", seqs, window_aa=200)
        assert closest == "very_similar"
        assert ident == pytest.approx(1.0)

    def test_target_missing(self):
        closest, ident = closest_paralog_identity("missing", {"a": "A" * 50})
        assert closest is None
        assert ident == 0.0

    def test_no_paralogs(self):
        closest, ident = closest_paralog_identity("only", {"only": "A" * 50})
        assert closest is None
        assert ident == 0.0


class TestIsHcrProbeFriendly:
    def test_fully_friendly(self):
        assert is_hcr_probe_friendly(cds_length_bp=900,
                                      paralog_min_identity=0.5,
                                      tandem_cluster_size=1)

    def test_too_short_unfriendly(self):
        assert not is_hcr_probe_friendly(cds_length_bp=200,
                                         paralog_min_identity=0.3,
                                         tandem_cluster_size=1)

    def test_paralog_too_similar(self):
        assert not is_hcr_probe_friendly(cds_length_bp=900,
                                         paralog_min_identity=0.95,
                                         tandem_cluster_size=1)

    def test_tandem_cluster_warning(self):
        assert not is_hcr_probe_friendly(cds_length_bp=900,
                                         paralog_min_identity=0.5,
                                         tandem_cluster_size=8)

    def test_unknown_paralog_identity_treated_as_friendly(self):
        assert is_hcr_probe_friendly(cds_length_bp=900,
                                      paralog_min_identity=None,
                                      tandem_cluster_size=1)

    def test_missing_cds_length_unfriendly(self):
        assert not is_hcr_probe_friendly(cds_length_bp=None,
                                         paralog_min_identity=0.3)
