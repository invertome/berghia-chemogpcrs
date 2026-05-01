"""Unit tests for scripts/parse_absrel.py.

Verifies that the JSON-to-CSV parser correctly extracts:
- omega_max (the chemoreceptor-relevant episodic signal)
- omega_mean (legacy weighted mean)
- weight_at_max (site fraction at max omega)
- proper CSV schema
"""
import csv
import json
import os
from pathlib import Path

import pytest

from parse_absrel import parse_absrel_json, FIELDNAMES


def _make_absrel_json(tmp_path: Path, branches: list[dict]) -> Path:
    """Helper to build a minimal aBSREL JSON file with given branch attrs."""
    data = {
        "branch attributes": {
            "0": {  # partition 0
                b["name"]: {
                    "original name": b["name"],
                    "Rate Distributions": [list(rc) for rc in b.get("rate_classes", [])],
                    "Uncorrected P-value": b.get("p", 0.5),
                    "Corrected P-value": b.get("p_corr", b.get("p", 0.5)),
                }
                for b in branches
            }
        }
    }
    f = tmp_path / "absrel.json"
    f.write_text(json.dumps(data))
    return f


def test_parses_episodic_positive_selection_branch(tmp_path):
    j = _make_absrel_json(tmp_path, [
        {"name": "berghia_TRINITY_DN1", "rate_classes": [(0.1, 0.9), (20.0, 0.1)],
         "p": 0.001, "p_corr": 0.01},
    ])
    out = tmp_path / "out.csv"
    n = parse_absrel_json(str(j), str(out))
    assert n == 1
    rows = list(csv.DictReader(open(out)))
    assert rows[0]["branch_id"] == "berghia_TRINITY_DN1"
    assert float(rows[0]["omega_max"]) == pytest.approx(20.0)
    assert float(rows[0]["omega_mean"]) == pytest.approx(0.1 * 0.9 + 20.0 * 0.1)
    assert float(rows[0]["weight_at_max"]) == pytest.approx(0.1)
    assert int(rows[0]["n_rate_classes"]) == 2
    assert int(rows[0]["is_selected"]) == 1   # p_corr<0.05 AND omega_max>1


def test_purifying_only_branch_marked_unselected(tmp_path):
    j = _make_absrel_json(tmp_path, [
        {"name": "ref_taxid_1234_5", "rate_classes": [(0.05, 1.0)],
         "p": 0.5, "p_corr": 0.99},
    ])
    out = tmp_path / "out.csv"
    parse_absrel_json(str(j), str(out))
    rows = list(csv.DictReader(open(out)))
    assert float(rows[0]["omega_max"]) == pytest.approx(0.05)
    assert int(rows[0]["is_selected"]) == 0


def test_csv_header_matches_schema(tmp_path):
    j = _make_absrel_json(tmp_path, [
        {"name": "b1", "rate_classes": [(1.0, 1.0)]},
    ])
    out = tmp_path / "out.csv"
    parse_absrel_json(str(j), str(out))
    with open(out) as f:
        header = next(csv.reader(f))
    assert header == FIELDNAMES


def test_atomic_mode_overwrites_existing(tmp_path):
    j1 = _make_absrel_json(tmp_path, [{"name": "b1", "rate_classes": [(0.5, 1.0)]}])
    j2 = _make_absrel_json(tmp_path, [{"name": "b2", "rate_classes": [(2.0, 1.0)]}])
    out = tmp_path / "out.csv"
    parse_absrel_json(str(j1), str(out))
    parse_absrel_json(str(j2), str(out))  # atomic mode: should overwrite
    rows = list(csv.DictReader(open(out)))
    assert len(rows) == 1
    assert rows[0]["branch_id"] == "b2"


def test_append_mode_validates_schema(tmp_path):
    """Pre-existing file with wrong schema must raise (bead -mqt)."""
    out = tmp_path / "out.csv"
    # Write a CSV with an OLD schema (missing the new columns)
    out.write_text("branch_id,omega,p_value\nfoo,1.0,0.5\n")
    j = _make_absrel_json(tmp_path, [{"name": "b1", "rate_classes": [(0.5, 1.0)]}])
    with pytest.raises(ValueError, match="schema"):
        parse_absrel_json(str(j), str(out), append_mode='append')


def test_handles_missing_rate_distributions(tmp_path):
    """Branch with no Rate Distributions falls back to scalar omega field."""
    data = {
        "branch attributes": {
            "0": {
                "branchx": {
                    "original name": "branchx",
                    "omega": 0.7,
                    "Uncorrected P-value": 0.3,
                    "Corrected P-value": 0.6,
                }
            }
        }
    }
    j = tmp_path / "absrel.json"
    j.write_text(json.dumps(data))
    out = tmp_path / "out.csv"
    parse_absrel_json(str(j), str(out))
    rows = list(csv.DictReader(open(out)))
    assert float(rows[0]["omega"]) == pytest.approx(0.7)
    assert int(rows[0]["n_rate_classes"]) == 0
