"""Unit tests for scripts/_validation_lib.py (bead -bdu)."""
import pytest

from _validation_lib import (
    average_rank,
    precision_at_n,
    recall_at_n,
    recall_curve,
)


class TestRecallAtN:
    def test_perfect_recall(self):
        ranked = ["a", "b", "c", "d", "e"]
        truth = {"a", "b"}
        assert recall_at_n(ranked, truth, 2) == 1.0

    def test_partial_recall(self):
        ranked = ["a", "b", "c", "d"]
        truth = {"a", "c"}
        assert recall_at_n(ranked, truth, 2) == 0.5
        assert recall_at_n(ranked, truth, 3) == 1.0

    def test_zero_recall(self):
        assert recall_at_n(["a", "b"], {"x"}, 2) == 0.0

    def test_empty_truth(self):
        assert recall_at_n(["a"], [], 1) == 0.0

    def test_n_capped_at_ranked_length(self):
        # N larger than ranked length -> capped
        assert recall_at_n(["a"], {"a"}, 100) == 1.0


class TestPrecisionAtN:
    def test_basic(self):
        ranked = ["a", "b", "c"]
        truth = {"a", "b"}
        assert precision_at_n(ranked, truth, 3) == pytest.approx(2/3)

    def test_zero_n(self):
        assert precision_at_n(["a"], {"a"}, 0) == 0.0


class TestAverageRank:
    def test_basic(self):
        ranked = ["a", "b", "c"]
        # Average of ranks of a (1) and c (3) = 2.0
        assert average_rank(ranked, {"a", "c"}) == 2.0

    def test_no_truth_found(self):
        assert average_rank(["a"], {"x"}) is None


class TestRecallCurve:
    def test_columns_present(self):
        df = recall_curve(["a", "b", "c", "d"], {"a", "c"}, ns=(2, 4))
        assert list(df.columns) == ["N", "recall", "precision",
                                    "n_truth_total", "n_truth_found_in_topN"]
        assert len(df) == 2
