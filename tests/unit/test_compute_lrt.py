"""Unit tests for scripts/compute_lrt.py (bead -cio)."""
import math

import pytest

from compute_lrt import lrt_pvalue, extract_lnL


class TestExtractLnL:
    def test_paml_format_with_ntime_np(self, tmp_path):
        f = tmp_path / "rst.txt"
        f.write_text("Time used: 12.3\nlnL(ntime:  5  np: 10):  -1234.567 +0.0\n")
        assert extract_lnL(str(f)) == pytest.approx(-1234.567)

    def test_simple_lnL_equals(self, tmp_path):
        f = tmp_path / "rst.txt"
        f.write_text("foo bar\nlnL = -987.654\n")
        assert extract_lnL(str(f)) == pytest.approx(-987.654)

    def test_positive_lnL(self, tmp_path):
        # Edge case: very short alignments can produce positive lnL
        f = tmp_path / "rst.txt"
        f.write_text("lnL(ntime: 2 np: 3): 12.345 +0.0\n")
        assert extract_lnL(str(f)) == pytest.approx(12.345)

    def test_no_lnL_raises(self, tmp_path):
        f = tmp_path / "rst.txt"
        f.write_text("nothing useful here\n")
        with pytest.raises(ValueError, match="No lnL"):
            extract_lnL(str(f))


class TestLrtPvalue:
    def test_zero_or_negative_stat_returns_one(self):
        assert lrt_pvalue(0, df=1) == 1.0
        assert lrt_pvalue(-1.5, df=2) == 1.0

    def test_df_argument_changes_p_value(self):
        chi2_stat = 5.0
        p1 = lrt_pvalue(chi2_stat, df=1)
        p2 = lrt_pvalue(chi2_stat, df=2)
        assert p1 < p2

    def test_boundary_halves_p_value(self):
        chi2_stat = 4.0
        p_normal = lrt_pvalue(chi2_stat, df=1, boundary=False)
        p_boundary = lrt_pvalue(chi2_stat, df=1, boundary=True)
        assert p_boundary == pytest.approx(0.5 * p_normal)

    def test_known_value_df1(self):
        # chi2(1).sf(3.84) ≈ 0.05
        assert lrt_pvalue(3.84, df=1) == pytest.approx(0.05, rel=1e-2)

    def test_known_value_df2(self):
        # chi2(2).sf(5.99) ≈ 0.05
        assert lrt_pvalue(5.99, df=2) == pytest.approx(0.05, rel=1e-2)
