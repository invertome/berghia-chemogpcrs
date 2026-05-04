"""Unit tests for scripts/parse_meme.py (bead -urk)."""
import json
from pathlib import Path

import pytest

from parse_meme import OG_FIELDS, SITE_FIELDS, aggregate, extract_sites, parse_one


def _meme_json(rows):
    """Build a minimal MEME-style JSON dict.
    rows is a list of [alpha, beta_plus, prop_plus, p_value, n_branches] sublists.
    """
    return {
        "MLE": {
            "headers": [
                ["alpha"], ["beta-"], ["p-"], ["beta+"], ["p+"],
                ["LRT"], ["p-value"], ["# branches"], ["post-mean"]
            ],
            "content": {
                "0": [
                    [r[0], 0.0, 0.0, r[1], r[2], 0.0, r[3], r[4], 0.0]
                    for r in rows
                ]
            }
        }
    }


class TestExtractSites:
    def test_basic_three_sites(self):
        d = _meme_json([
            [0.5, 5.0, 0.1, 0.001, 3],
            [0.5, 1.0, 0.0, 0.5,   0],
            [0.5, 8.0, 0.2, 0.04,  2],
        ])
        sites = extract_sites(d)
        assert len(sites) == 3
        assert sites[0]["beta_plus"] == 5.0
        assert sites[0]["p_value"] == pytest.approx(0.001)
        assert sites[2]["n_branches_episodic"] == 2

    def test_empty(self):
        assert extract_sites({}) == []


class TestAggregate:
    def test_counts_episodic_sites_correctly(self):
        sites = [
            {"p_value": 0.001, "beta_plus": 10.0},
            {"p_value": 0.04,  "beta_plus": 5.0},
            {"p_value": 0.5,   "beta_plus": 1.0},
            {"p_value": 0.02,  "beta_plus": 0.5},
        ]
        out = aggregate(sites, alpha=0.05)
        assert out["n_sites_total"] == 4
        assert out["n_sites_episodic"] == 2
        assert out["fraction_episodic"] == pytest.approx(0.5)

    def test_no_episodic(self):
        sites = [{"p_value": 0.5, "beta_plus": 0.5}]
        out = aggregate(sites, alpha=0.05)
        assert out["n_sites_episodic"] == 0
        assert out["fraction_episodic"] == 0.0


class TestParseOne:
    def test_schema_fields(self, tmp_path):
        f = tmp_path / "OGxx_meme.json"
        f.write_text(json.dumps(_meme_json([[0.5, 8.0, 0.1, 0.01, 2]])))
        sites, og = parse_one(str(f), "OGxx", alpha=0.05)
        assert og["og_name"] == "OGxx"
        for k in OG_FIELDS:
            assert k in og, k
        for k in ("og_name", "site", "is_episodic"):
            assert k in sites[0], k

    def test_episodic_flag_per_site(self, tmp_path):
        f = tmp_path / "x_meme.json"
        f.write_text(json.dumps(_meme_json([
            [0.5, 8.0, 0.1, 0.001, 2],
            [0.5, 0.5, 0.0, 0.001, 0],
        ])))
        sites, og = parse_one(str(f), "OGyy", alpha=0.05)
        assert sites[0]["is_episodic"] == 1
        assert sites[1]["is_episodic"] == 0
        assert og["n_sites_episodic"] == 1
