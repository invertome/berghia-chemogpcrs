"""Tests for the positive-control status split (bead 444, change 2).

``scripts/check_positive_controls.py`` used to set ``alert=True`` whenever a
control was not found in the ranked CSV. The only control shipped in
``references/hcr_positive_controls.csv`` is ``Galpha_olf`` — a G-protein alpha
subunit, NOT a class-A GPCR — so it can never appear in a class-A
chemoreceptor ranked CSV. The check therefore rendered a red ALERT on every
single run, which trains the reader to ignore it and destroys its value as a
drift detector.

The fix separates the two conditions:

  * ``not_found``             -> informational, ``alert=False``
  * ``found_below_percentile`` -> the real drift signal, ``alert=True``
  * ``found_healthy``          -> ``alert=False``

Absence must stay fully visible (its own row, its own ``status``, its own
count in the stderr summary) — the point is to stop crying wolf, not to
suppress information.

These tests complement ``test_check_positive_controls.py`` (match strategies,
percentile arithmetic); here we pin the status/alert semantics and the summary
counts.
"""
from __future__ import annotations

import csv
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

STATUS_NOT_FOUND = "not_found"
STATUS_BELOW = "found_below_percentile"
STATUS_HEALTHY = "found_healthy"


def _make_ranked_csv(path: Path, ids: list[str]) -> None:
    scores = [10.0 - i for i in range(len(ids))]
    pd.DataFrame({"id": ids, "rank_score": scores}).to_csv(path, index=False)


def _make_controls_csv(path: Path, names: list[str]) -> None:
    fieldnames = ["gene_name", "aliases", "refseq_protein", "notes"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for n in names:
            w.writerow({"gene_name": n, "aliases": "", "refseq_protein": "",
                        "notes": ""})


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


def _summary_counts(stderr: str) -> dict[str, int]:
    """Parse the machine-readable count line the script emits.

    Format: ``positive controls: total=3 healthy=1 below_percentile=1 not_found=1``
    """
    m = re.search(
        r"positive controls: total=(\d+) healthy=(\d+) "
        r"below_percentile=(\d+) not_found=(\d+)", stderr)
    assert m, f"no parseable summary count line in stderr:\n{stderr}"
    return {"total": int(m.group(1)), "healthy": int(m.group(2)),
            "below_percentile": int(m.group(3)), "not_found": int(m.group(4))}


# --- absence is informational, never an alert --------------------------------

def test_absent_control_does_not_raise_alert(tmp_path: Path) -> None:
    """The Galpha_olf case: a control that cannot appear in a class-A ranked
    CSV must NOT be flagged as a ranking regression."""
    ranked, controls, out = (tmp_path / "r.csv", tmp_path / "c.csv",
                             tmp_path / "o.tsv")
    _make_ranked_csv(ranked, ["a", "b", "c"])
    _make_controls_csv(controls, ["Galpha_olf"])

    rc = _run_cli(ranked, controls, out)
    assert rc.returncode == 0, rc.stderr

    row = _read_tsv(out)[0]
    assert row["found"] == "False"
    assert row["status"] == STATUS_NOT_FOUND
    assert row["alert"] == "False", "absence must not set alert"


def test_absent_control_is_still_fully_reported(tmp_path: Path) -> None:
    """Not hidden, not dropped: its own row in the TSV and its own count in
    the stderr summary."""
    ranked, controls, out = (tmp_path / "r.csv", tmp_path / "c.csv",
                             tmp_path / "o.tsv")
    _make_ranked_csv(ranked, ["a", "b", "c"])
    _make_controls_csv(controls, ["Galpha_olf", "AlsoMissing"])

    rc = _run_cli(ranked, controls, out)
    rows = _read_tsv(out)
    assert len(rows) == 2, "absent controls must still get rows"
    assert {r["gene_name"] for r in rows} == {"Galpha_olf", "AlsoMissing"}
    assert all(r["status"] == STATUS_NOT_FOUND for r in rows)
    assert all(r["matched_strategy"] == "none" for r in rows)

    assert _summary_counts(rc.stderr)["not_found"] == 2
    assert "not found" in rc.stderr.lower(), "absence must be named in the summary"


def test_all_absent_uses_ok_branch_not_warn(tmp_path: Path) -> None:
    """No rank regression anywhere -> the verdict line is OK, not WARN."""
    ranked, controls, out = (tmp_path / "r.csv", tmp_path / "c.csv",
                             tmp_path / "o.tsv")
    _make_ranked_csv(ranked, ["a", "b", "c"])
    _make_controls_csv(controls, ["Galpha_olf"])

    rc = _run_cli(ranked, controls, out)
    assert "OK:" in rc.stderr
    assert "WARN:" not in rc.stderr


# --- found-but-low is the real drift signal ----------------------------------

def test_found_below_percentile_raises_alert(tmp_path: Path) -> None:
    ranked, controls, out = (tmp_path / "r.csv", tmp_path / "c.csv",
                             tmp_path / "o.tsv")
    ids = [f"id_{i}" for i in range(10)]
    ids[7] = "sinker"  # rank 8/10 -> percentile 30
    _make_ranked_csv(ranked, ids)
    _make_controls_csv(controls, ["sinker"])

    rc = _run_cli(ranked, controls, out)
    row = _read_tsv(out)[0]
    assert row["found"] == "True"
    assert float(row["percentile"]) < 50.0
    assert row["status"] == STATUS_BELOW
    assert row["alert"] == "True"
    assert "WARN:" in rc.stderr


def test_found_and_healthy_raises_no_alert(tmp_path: Path) -> None:
    ranked, controls, out = (tmp_path / "r.csv", tmp_path / "c.csv",
                             tmp_path / "o.tsv")
    _make_ranked_csv(ranked, ["winner", "b", "c", "d"])
    _make_controls_csv(controls, ["winner"])

    rc = _run_cli(ranked, controls, out)
    row = _read_tsv(out)[0]
    assert row["found"] == "True"
    assert row["status"] == STATUS_HEALTHY
    assert row["alert"] == "False"
    assert "OK:" in rc.stderr


def test_alert_percentile_boundary_moves_status(tmp_path: Path) -> None:
    """The same control flips status/alert purely on the threshold."""
    ranked, controls, out = (tmp_path / "r.csv", tmp_path / "c.csv",
                             tmp_path / "o.tsv")
    ids = [f"id_{i}" for i in range(10)]
    ids[7] = "borderline"  # percentile 30
    _make_ranked_csv(ranked, ids)
    _make_controls_csv(controls, ["borderline"])

    strict = _run_cli(ranked, controls, out, alert_percentile=50.0)
    assert _read_tsv(out)[0]["status"] == STATUS_BELOW
    assert _read_tsv(out)[0]["alert"] == "True"
    assert "WARN:" in strict.stderr

    lax = _run_cli(ranked, controls, out, alert_percentile=20.0)
    assert _read_tsv(out)[0]["status"] == STATUS_HEALTHY
    assert _read_tsv(out)[0]["alert"] == "False"
    assert "OK:" in lax.stderr


# --- summary counts each category separately ---------------------------------

def test_summary_counts_each_category(tmp_path: Path) -> None:
    """One healthy, one below, one absent -> three distinct counts, and only
    the below-percentile one drives the alert total."""
    ranked, controls, out = (tmp_path / "r.csv", tmp_path / "c.csv",
                             tmp_path / "o.tsv")
    ids = [f"id_{i}" for i in range(20)]
    ids[0] = "top_dog"
    ids[18] = "tail_lurker"
    _make_ranked_csv(ranked, ids)
    _make_controls_csv(controls, ["top_dog", "tail_lurker", "ghost"])

    rc = _run_cli(ranked, controls, out)
    assert rc.returncode == 0, rc.stderr

    by_name = {r["gene_name"]: r for r in _read_tsv(out)}
    assert by_name["top_dog"]["status"] == STATUS_HEALTHY
    assert by_name["top_dog"]["alert"] == "False"
    assert by_name["tail_lurker"]["status"] == STATUS_BELOW
    assert by_name["tail_lurker"]["alert"] == "True"
    assert by_name["ghost"]["status"] == STATUS_NOT_FOUND
    assert by_name["ghost"]["alert"] == "False"

    counts = _summary_counts(rc.stderr)
    assert counts == {"total": 3, "healthy": 1,
                      "below_percentile": 1, "not_found": 1}
    # only the genuine regression is counted as an alert
    assert "WARN: 1/3" in rc.stderr


def test_summary_counts_with_no_controls(tmp_path: Path) -> None:
    ranked, controls, out = (tmp_path / "r.csv", tmp_path / "c.csv",
                             tmp_path / "o.tsv")
    _make_ranked_csv(ranked, ["a", "b"])
    _make_controls_csv(controls, [])

    rc = _run_cli(ranked, controls, out)
    assert rc.returncode == 0, rc.stderr
    assert _read_tsv(out) == []
    assert _summary_counts(rc.stderr) == {
        "total": 0, "healthy": 0, "below_percentile": 0, "not_found": 0}


# --- schema / contract -------------------------------------------------------

def test_existing_columns_are_preserved_and_status_added(tmp_path: Path) -> None:
    """Stage 09 reads gene_name/found/rank/percentile/alert by name — those
    must survive. `status` is ADDED, nothing renamed."""
    ranked, controls, out = (tmp_path / "r.csv", tmp_path / "c.csv",
                             tmp_path / "o.tsv")
    _make_ranked_csv(ranked, ["x"])
    _make_controls_csv(controls, ["x"])

    rc = _run_cli(ranked, controls, out)
    assert rc.returncode == 0, rc.stderr
    with open(out) as fh:
        header = fh.readline().rstrip("\n").split("\t")

    legacy = {"gene_name", "found", "candidate_id", "rank", "percentile",
              "rank_score", "matched_strategy", "alert", "notes"}
    assert legacy.issubset(set(header)), f"dropped legacy columns: {legacy - set(header)}"
    assert "status" in header


def test_exit_code_is_zero_even_when_alerting(tmp_path: Path) -> None:
    """Non-fatal contract: this is a sanity check, not a gate."""
    ranked, controls, out = (tmp_path / "r.csv", tmp_path / "c.csv",
                             tmp_path / "o.tsv")
    ids = [f"id_{i}" for i in range(10)]
    ids[9] = "bottom"
    _make_ranked_csv(ranked, ids)
    _make_controls_csv(controls, ["bottom", "ghost"])

    rc = _run_cli(ranked, controls, out)
    assert rc.returncode == 0, rc.stderr


def test_module_docstring_describes_split_not_always_alert() -> None:
    """The docstring must not still claim absence raises an alert."""
    text = (REPO_ROOT / "scripts" / "check_positive_controls.py").read_text()
    doc = text.split('"""')[1]
    assert "status" in doc, "docstring must document the new status column"
    lowered = doc.lower()
    assert "not found" in lowered or "not_found" in lowered
