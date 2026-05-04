"""Unit tests for scripts/parse_busted.py (bead -urk)."""
import json
from pathlib import Path

import pytest

from parse_busted import FIELDNAMES, extract_busted, parse_one


def _busted_json(p_value, omegas_props):
    """Build a minimal BUSTED-style JSON dict."""
    classes = [{"omega": o, "proportion": w} for o, w in omegas_props]
    return {
        "test results": {"p-value": p_value, "LRT": 12.34},
        "fits": {
            "Unconstrained model": {
                "Rate Distributions": {"Test": classes}
            }
        }
    }


class TestExtractBusted:
    def test_three_class_distribution(self):
        d = _busted_json(0.001, [(0.05, 0.7), (1.0, 0.25), (15.0, 0.05)])
        out = extract_busted(d)
        assert out["p_value"] == pytest.approx(0.001)
        assert out["lrt"] == pytest.approx(12.34)
        assert out["omega_low"] == pytest.approx(0.05)
        assert out["omega_mid"] == pytest.approx(1.0)
        assert out["omega_high"] == pytest.approx(15.0)
        assert out["prop_low"] == pytest.approx(0.7)
        assert out["prop_high"] == pytest.approx(0.05)

    def test_two_class_padded_to_three(self):
        d = _busted_json(0.5, [(0.1, 0.9), (5.0, 0.1)])
        out = extract_busted(d)
        # mid class = NaN when only 2 fitted classes
        import math
        assert math.isnan(out["omega_mid"])
        assert out["omega_high"] == pytest.approx(5.0)

    def test_missing_test_results(self):
        d = {"fits": {}}
        out = extract_busted(d)
        import math
        assert math.isnan(out["p_value"])

    def test_list_format_rate_classes(self):
        d = {
            "test results": {"p-value": 0.04, "LRT": 5.0},
            "fits": {
                "Unconstrained model": {
                    "Rate Distributions": {"Test": [[0.1, 0.8], [10.0, 0.2]]}
                }
            }
        }
        out = extract_busted(d)
        assert out["omega_high"] == pytest.approx(10.0)


class TestParseOne:
    def test_significance_flag(self, tmp_path):
        f = tmp_path / "x_busted_s.json"
        f.write_text(json.dumps(_busted_json(0.01, [(0.1, 0.95), (10.0, 0.05)])))
        row = parse_one(str(f), "OGxx", "S", alpha=0.05)
        assert row["og_name"] == "OGxx"
        assert row["model_variant"] == "S"
        assert row["is_significant"] == 1

    def test_below_alpha_threshold_not_significant(self, tmp_path):
        f = tmp_path / "x_busted_s.json"
        f.write_text(json.dumps(_busted_json(0.5, [(0.1, 0.95), (10.0, 0.05)])))
        row = parse_one(str(f), "OGyy", "S", alpha=0.05)
        assert row["is_significant"] == 0

    def test_schema_fields_present(self, tmp_path):
        f = tmp_path / "x.json"
        f.write_text(json.dumps(_busted_json(0.001, [(0.1, 1.0)])))
        row = parse_one(str(f), "OGzz", "MH", alpha=0.05)
        for k in FIELDNAMES:
            assert k in row, k
