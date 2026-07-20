"""Tests for scripts/check_positive_controls.py.

The script reads a ranked-candidates CSV and a controls CSV, finds each
control's rank/percentile in the ranking, and emits a TSV carrying a `status`
(found_healthy / found_below_percentile / not_found) plus an `alert` flag.
Only a found-but-below-percentile control raises `alert`; absence is
informational (bead 444). The TSV is consumed by stage 09 (positive-control
sanity-check section).

Tests cover the four match strategies, the percentile arithmetic, and the
missing-control case. The status/alert split itself is pinned in more detail
by test_positive_controls_status.py.
"""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent.parent


def _make_ranked_csv(path: Path, ids: list[str], scores: list[float] | None = None) -> None:
    """Ranked CSV with descending scores. The script sorts by rank_score so
    order in the input doesn't matter, but we keep it descending for clarity."""
    if scores is None:
        scores = [10.0 - i for i in range(len(ids))]
    pd.DataFrame({"id": ids, "rank_score": scores}).to_csv(path, index=False)


def _make_controls_csv(path: Path, rows: list[dict]) -> None:
    """Controls CSV with at least gene_name; aliases / refseq_protein optional."""
    fieldnames = ["gene_name", "aliases", "refseq_protein", "notes"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})


def _run_cli(ranked: Path, controls: Path, out: Path,
             alert_percentile: float | None = None) -> subprocess.CompletedProcess:
    args = [sys.executable, "scripts/check_positive_controls.py",
            "--ranked-csv", str(ranked),
            "--controls-csv", str(controls),
            "--out", str(out)]
    if alert_percentile is not None:
        args += ["--alert-percentile", str(alert_percentile)]
    return subprocess.run(args, cwd=REPO_ROOT, capture_output=True, text=True)


def _read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_exact_id_match_top_ranked_no_alert(tmp_path: Path) -> None:
    """Control whose gene_name matches a top-ranked id exactly: found, top
    percentile, no alert."""
    ranked = tmp_path / "r.csv"
    controls = tmp_path / "c.csv"
    out = tmp_path / "out.tsv"
    _make_ranked_csv(ranked, ["Galpha_olf", "x", "y", "z"])
    _make_controls_csv(controls, [{"gene_name": "Galpha_olf",
                                   "notes": "HCR-validated"}])

    rc = _run_cli(ranked, controls, out)
    assert rc.returncode == 0, f"stderr: {rc.stderr}"
    assert "OK:" in rc.stderr  # success-summary line

    rows = _read_tsv(out)
    assert len(rows) == 1
    r = rows[0]
    assert r["gene_name"] == "Galpha_olf"
    assert r["found"] == "True"
    assert r["candidate_id"] == "Galpha_olf"
    assert int(r["rank"]) == 1
    assert float(r["percentile"]) == 100.0
    assert r["matched_strategy"] == "exact_id"
    assert r["alert"] == "False"


def test_not_found_is_reported_without_alert(tmp_path: Path) -> None:
    """Control absent from the ranking: blank rank/percentile,
    matched_strategy=none, status=not_found.

    Absence is informational, NOT an alert (bead 444) — the shipped control
    Ga_olf is a G-protein alpha subunit and can never appear in a class-A
    ranked CSV, so alerting on it fired red every run.
    """
    ranked = tmp_path / "r.csv"
    controls = tmp_path / "c.csv"
    out = tmp_path / "out.tsv"
    _make_ranked_csv(ranked, ["a", "b", "c"])
    _make_controls_csv(controls, [{"gene_name": "MissingGene"}])

    rc = _run_cli(ranked, controls, out)
    assert rc.returncode == 0, f"stderr: {rc.stderr}"
    # no rank regression -> quiet branch, but the absence is still reported
    assert "OK:" in rc.stderr
    assert "not_found=1" in rc.stderr

    rows = _read_tsv(out)
    assert len(rows) == 1
    r = rows[0]
    assert r["found"] == "False"
    assert r["candidate_id"] == ""
    assert r["rank"] == ""
    assert r["percentile"] == ""
    assert r["matched_strategy"] == "none"
    assert r["status"] == "not_found"
    assert r["alert"] == "False"


def test_below_percentile_threshold_emits_alert(tmp_path: Path) -> None:
    """Control found but below the alert percentile (default 50): alert=True
    even though it's present in the ranking."""
    ranked = tmp_path / "r.csv"
    controls = tmp_path / "c.csv"
    out = tmp_path / "out.tsv"
    # 10 candidates; sneaky_ctrl is rank 8 of 10 -> percentile ≈ 30
    ids = [f"id_{i}" for i in range(10)]
    ids[7] = "sneaky_ctrl"
    _make_ranked_csv(ranked, ids)
    _make_controls_csv(controls, [{"gene_name": "sneaky_ctrl"}])

    rc = _run_cli(ranked, controls, out)
    assert rc.returncode == 0, f"stderr: {rc.stderr}"
    rows = _read_tsv(out)
    r = rows[0]
    assert r["found"] == "True"
    assert int(r["rank"]) == 8
    assert float(r["percentile"]) < 50.0
    assert r["alert"] == "True"
    assert "WARN:" in rc.stderr


def test_custom_alert_percentile_overrides_default(tmp_path: Path) -> None:
    """Lowering the threshold lets a borderline control pass without alert."""
    ranked = tmp_path / "r.csv"
    controls = tmp_path / "c.csv"
    out = tmp_path / "out.tsv"
    ids = [f"id_{i}" for i in range(10)]
    ids[7] = "borderline"  # percentile ≈ 30
    _make_ranked_csv(ranked, ids)
    _make_controls_csv(controls, [{"gene_name": "borderline"}])

    rc = _run_cli(ranked, controls, out, alert_percentile=20.0)
    assert rc.returncode == 0, f"stderr: {rc.stderr}"
    rows = _read_tsv(out)
    assert rows[0]["alert"] == "False"
    assert "OK:" in rc.stderr


def test_multiple_controls_mixed_results(tmp_path: Path) -> None:
    """Three controls — one top-ranked, one below threshold, one missing —
    each gets its own row with the right semantics. The summary line counts
    alerts."""
    ranked = tmp_path / "r.csv"
    controls = tmp_path / "c.csv"
    out = tmp_path / "out.tsv"
    ids = [f"id_{i}" for i in range(20)]
    ids[0] = "top_dog"
    ids[18] = "tail_lurker"
    _make_ranked_csv(ranked, ids)
    _make_controls_csv(controls, [
        {"gene_name": "top_dog"},
        {"gene_name": "tail_lurker"},
        {"gene_name": "ghost"},
    ])

    rc = _run_cli(ranked, controls, out)
    assert rc.returncode == 0, f"stderr: {rc.stderr}"
    rows = _read_tsv(out)
    assert len(rows) == 3

    by_name = {r["gene_name"]: r for r in rows}
    assert by_name["top_dog"]["alert"] == "False"
    assert by_name["tail_lurker"]["alert"] == "True"
    # absent control is reported but does not alert (bead 444)
    assert by_name["ghost"]["alert"] == "False"
    assert by_name["ghost"]["status"] == "not_found"
    # only the genuine rank regression counts as an alert
    assert "WARN: 1/3" in rc.stderr
    assert "not_found=1" in rc.stderr


def test_empty_controls_csv_is_valid_no_rows(tmp_path: Path) -> None:
    """A controls CSV with only headers produces a TSV with only headers."""
    ranked = tmp_path / "r.csv"
    controls = tmp_path / "c.csv"
    out = tmp_path / "out.tsv"
    _make_ranked_csv(ranked, ["a", "b", "c"])
    _make_controls_csv(controls, [])

    rc = _run_cli(ranked, controls, out)
    assert rc.returncode == 0, f"stderr: {rc.stderr}"
    rows = _read_tsv(out)
    assert rows == []
    # 0 of 0 alerts: summary is the OK branch
    assert "OK: all 0" in rc.stderr


def test_output_columns_are_exactly_the_documented_set(tmp_path: Path) -> None:
    """Stage 09 reads gene_name, found, rank, percentile, alert by name —
    schema must be stable."""
    ranked = tmp_path / "r.csv"
    controls = tmp_path / "c.csv"
    out = tmp_path / "out.tsv"
    _make_ranked_csv(ranked, ["x"])
    _make_controls_csv(controls, [{"gene_name": "x"}])

    rc = _run_cli(ranked, controls, out)
    assert rc.returncode == 0, f"stderr: {rc.stderr}"
    with open(out) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    expected = {"gene_name", "found", "candidate_id", "rank", "percentile",
                "rank_score", "matched_strategy", "status", "alert", "notes"}
    assert set(header) == expected, header
