"""Unit tests for scripts/calibrate_depth_threshold.py (bead -m6k)."""
import json

import pytest

from calibrate_depth_threshold import simulate_null_depths, summarize


def test_small_simulation_runs():
    d = simulate_null_depths(n_taxa=20, n_sims=5, birth_rate=1.0,
                             death_rate=0.0, seed=42)
    assert len(d["per_tree_max"]) == 5
    assert len(d["all_depths"]) > 0
    assert all(x >= 0 for x in d["all_depths"])


def test_summarize_returns_expected_keys():
    d = simulate_null_depths(n_taxa=20, n_sims=5, seed=7)
    s = summarize(d)
    for k in ("n_simulations", "n_internal_nodes", "depth_p95", "depth_p99",
             "per_tree_max_p95"):
        assert k in s, k
    assert s["depth_p95"] >= s["depth_median"] >= 0
    assert s["depth_p99"] >= s["depth_p95"]


def test_higher_birth_rate_compresses_depths():
    """Higher birth rate -> tree builds faster -> shorter root-to-tip distances."""
    slow = simulate_null_depths(n_taxa=30, n_sims=20, birth_rate=0.5,
                                death_rate=0.0, seed=1)
    fast = simulate_null_depths(n_taxa=30, n_sims=20, birth_rate=2.0,
                                death_rate=0.0, seed=1)
    s_slow = summarize(slow)
    s_fast = summarize(fast)
    assert s_slow["depth_mean"] > s_fast["depth_mean"]


def test_seed_reproducibility():
    d1 = simulate_null_depths(n_taxa=20, n_sims=3, seed=99)
    d2 = simulate_null_depths(n_taxa=20, n_sims=3, seed=99)
    # dendropy's birth_death_tree isn't bit-reproducible across calls even with an
    # identical seeded rng: the draw COUNT is deterministic but draw-to-lineage order
    # is object-identity (memory-address) based, giving ULP-level (~1e-15) depth
    # differences (not fixed by PYTHONHASHSEED — it's id-based, not str-hash). The
    # percentile thresholds this module computes are robust to that, so assert
    # numerical, not bit-exact, reproducibility; a real break (e.g. a different seed)
    # differs by O(0.1) and still fails this.
    assert d1["per_tree_max"] == pytest.approx(d2["per_tree_max"], rel=1e-9)
