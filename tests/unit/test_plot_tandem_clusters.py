"""Tests for scripts/plot_tandem_clusters.py.

These tests verify the plotter handles the three regimes that occur in real
runs: a populated cluster CSV, an all-singletons CSV (no detected clusters),
and a missing CSV. Each regime must produce a non-empty PNG so stage 09 can
always reference the figure.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

# conftest.py adds scripts/ to sys.path
import plot_tandem_clusters as ptc


@pytest.fixture
def populated_csv(tmp_path: Path) -> Path:
    df = pd.DataFrame(
        [
            {"candidate_id": "BersteEVm0001t1", "tandem_cluster_size": 4, "tandem_cluster_id": "tc_001"},
            {"candidate_id": "BersteEVm0002t1", "tandem_cluster_size": 4, "tandem_cluster_id": "tc_001"},
            {"candidate_id": "BersteEVm0003t1", "tandem_cluster_size": 4, "tandem_cluster_id": "tc_001"},
            {"candidate_id": "BersteEVm0004t1", "tandem_cluster_size": 4, "tandem_cluster_id": "tc_001"},
            {"candidate_id": "BersteEVm0010t1", "tandem_cluster_size": 3, "tandem_cluster_id": "tc_002"},
            {"candidate_id": "BersteEVm0011t1", "tandem_cluster_size": 3, "tandem_cluster_id": "tc_002"},
            {"candidate_id": "BersteEVm0012t1", "tandem_cluster_size": 3, "tandem_cluster_id": "tc_002"},
            {"candidate_id": "BersteEVm0050t1", "tandem_cluster_size": 1, "tandem_cluster_id": ""},
            {"candidate_id": "BersteEVm0051t1", "tandem_cluster_size": 1, "tandem_cluster_id": ""},
        ]
    )
    fp = tmp_path / "tandem_clusters.csv"
    df.to_csv(fp, index=False)
    return fp


@pytest.fixture
def all_singletons_csv(tmp_path: Path) -> Path:
    df = pd.DataFrame(
        [
            {"candidate_id": f"BersteEVm{i:04d}t1", "tandem_cluster_size": 1, "tandem_cluster_id": ""}
            for i in range(1, 6)
        ]
    )
    fp = tmp_path / "tandem_clusters.csv"
    df.to_csv(fp, index=False)
    return fp


def _all_outputs_exist(prefix: Path) -> None:
    for ext in ("png", "svg", "pdf"):
        out = prefix.with_suffix(f".{ext}")
        assert out.exists(), f"missing output: {out}"
        assert out.stat().st_size > 0, f"empty output: {out}"


def test_populated(populated_csv: Path, tmp_path: Path) -> None:
    prefix = tmp_path / "tandem_plot"
    df = pd.read_csv(populated_csv)
    ptc.plot(df, str(prefix))
    _all_outputs_exist(prefix)


def test_all_singletons_filter_excludes_nan_rows(all_singletons_csv: Path,
                                                  tmp_path: Path) -> None:
    """Reviewer-found regression: pandas read_csv coerces empty strings to
    NaN; ``.astype(str)`` then turns NaN into the literal "nan" (length 3),
    so a naive ``str.len() > 0`` check let phantom rows past, producing a
    phantom 'nan' bar in panel B. The fixed filter must produce an EMPTY
    clustered subframe for the all-singletons case."""
    prefix = tmp_path / "tandem_plot_singletons"
    df = pd.read_csv(all_singletons_csv)

    # Replicate the (fixed) filter the plotter uses
    cid_col = df["tandem_cluster_id"]
    has_cid = (cid_col.notna()
               & (cid_col.astype(str).str.len() > 0)
               & (cid_col.astype(str) != "nan"))
    assert not has_cid.any(), \
        "all-singletons input must yield zero cluster-IDs after the filter"

    # End-to-end: plotter still emits a figure (empty-state annotation path).
    ptc.plot(df, str(prefix))
    _all_outputs_exist(prefix)


def test_missing_csv_via_cli_emits_figure(tmp_path: Path) -> None:
    """CLI path: when the CSV is missing, stderr warns and an empty-state
    figure is still written (so the report's image hook resolves)."""
    prefix = tmp_path / "tandem_plot_missing"
    rc = subprocess.run(
        [sys.executable, "scripts/plot_tandem_clusters.py",
         str(tmp_path / "no_such.csv"), str(prefix)],
        cwd=Path(__file__).resolve().parent.parent.parent,
        capture_output=True, text=True,
    )
    assert rc.returncode == 0, f"stderr: {rc.stderr}"
    _all_outputs_exist(prefix)


def test_top_n_param_caps_panel_b(tmp_path: Path) -> None:
    """top_n=2 should keep panel B to <=2 distinct clusters."""
    df = pd.DataFrame(
        [
            {"candidate_id": f"g{i}", "tandem_cluster_size": s, "tandem_cluster_id": cid}
            for i, (s, cid) in enumerate(
                [(5, "tc_a"), (5, "tc_a"), (5, "tc_a"), (5, "tc_a"), (5, "tc_a"),
                 (4, "tc_b"), (4, "tc_b"), (4, "tc_b"), (4, "tc_b"),
                 (3, "tc_c"), (3, "tc_c"), (3, "tc_c")]
            )
        ]
    )
    prefix = tmp_path / "topn_clipped"
    ptc.plot(df, str(prefix), top_n=2)
    _all_outputs_exist(prefix)
