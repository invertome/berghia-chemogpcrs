"""Tests for ``scripts/plot_ranking.py`` displaying ``rank_tier``.

Companion to ``test_rank_tier_display_stage09.py``. The ranking figure's tier
panel, its per-candidate tier markers, and the sidecar summary CSV all keyed on
the legacy ``confidence_tier`` (an uncalibrated evidence-completeness heuristic
derived from the demoted hand-weighted composite). They now key on
``rank_tier`` -- the signal-bootstrap P(in top-k) tier from
``scripts/rank_confidence.py``.

``plot_ranking.py`` is a top-level script (it reads ``sys.argv`` at import), so
these tests EXECUTE it against synthetic ranked CSVs and assert on the artefacts
it writes, rather than importing it.
"""
from __future__ import annotations

import csv
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
PLOT_RANKING = REPO_ROOT / "scripts" / "plot_ranking.py"

TIERS = ("high", "plausible", "tail")

# Component columns gate plot_ranking's multi-panel path (has_components).
COMPONENT_COLUMNS = [
    "phylo_score", "purifying_score", "positive_score",
    "synteny_score", "expression_score", "lse_divergence_score",
]
RANKED_COLUMNS = (
    ["id", "rank_score", "confidence_tier", "rank_tier", "p_top_k",
     "evidence_completeness", "selection_significant"] + COMPONENT_COLUMNS
)

pytest.importorskip("matplotlib")
pytest.importorskip("seaborn")


def _row(idx: int, tier: str) -> dict:
    row = {
        "id": f"cand{idx}",
        "rank_score": f"{1.0 - idx * 0.01:.3f}",
        "confidence_tier": "High",          # legacy column: present, ignored
        "rank_tier": tier,
        "p_top_k": "0.5",
        "evidence_completeness": "0.8",
        "selection_significant": "False",
    }
    for c in COMPONENT_COLUMNS:
        row[c] = "0.5"
    return row


def _run_plot(tmp_path: Path, rows: list[dict], columns=None) -> dict:
    """Run plot_ranking.py; return its summary CSV as a Metric->Value dict."""
    cols = list(columns if columns is not None else RANKED_COLUMNS)
    src = tmp_path / "ranked.csv"
    with src.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({c: row.get(c, "") for c in cols})

    out = tmp_path / "ranking_plot"
    proc = subprocess.run(
        [sys.executable, str(PLOT_RANKING), str(src), str(out)],
        capture_output=True, text=True, timeout=300,
        env={"MPLBACKEND": "Agg", "PATH": "/usr/bin:/bin", "HOME": str(tmp_path)},
    )
    assert proc.returncode == 0, f"plot_ranking failed:\n{proc.stderr[-4000:]}"
    assert (tmp_path / "ranking_plot.png").exists(), "no PNG produced"

    summary = {}
    with (tmp_path / "ranking_plot_summary.csv").open() as fh:
        for rec in csv.DictReader(fh):
            summary[rec["Metric"]] = rec["Value"]
    return summary


def _mixed_rows() -> list[dict]:
    rows = [_row(i, "high") for i in range(3)]
    rows += [_row(10 + i, "plausible") for i in range(4)]
    rows += [_row(20 + i, "tail") for i in range(5)]
    rows += [_row(30 + i, "") for i in range(2)]
    return rows


# --- Summary CSV ------------------------------------------------------------

def test_summary_counts_each_rank_tier(tmp_path: Path) -> None:
    summary = _run_plot(tmp_path, _mixed_rows())
    assert summary["Rank tier: high"] == "3"
    assert summary["Rank tier: plausible"] == "4"
    assert summary["Rank tier: tail"] == "5"
    assert summary["Rank tier: not scored"] == "2"


def test_summary_tiers_partition_the_candidate_set(tmp_path: Path) -> None:
    rows = _mixed_rows()
    summary = _run_plot(tmp_path, rows)
    total = sum(int(summary[f"Rank tier: {t}"]) for t in TIERS)
    total += int(summary["Rank tier: not scored"])
    assert total == int(summary["Total candidates"]) == len(rows)


def test_summary_no_longer_reports_legacy_confidence_tiers(tmp_path: Path) -> None:
    """The legacy column is still in the CSV but must not be summarised."""
    summary = _run_plot(tmp_path, _mixed_rows())
    for legacy in ("High confidence", "Medium confidence", "Low confidence"):
        assert legacy not in summary, f"legacy metric still emitted: {legacy}"


def test_summary_reports_not_computed_when_column_absent(tmp_path: Path) -> None:
    """rank_confidence.py is non-fatal; a run where it failed has no column.

    The summary must say the statistic is unavailable rather than emit zeros,
    which would read as a genuine "nothing ranks stably" finding.
    """
    cols = [c for c in RANKED_COLUMNS if c != "rank_tier"]
    summary = _run_plot(tmp_path, [_row(i, "high") for i in range(5)], columns=cols)
    joined = " ".join(f"{k}={v}" for k, v in summary.items()).lower()
    assert "not computed" in joined
    for tier in TIERS:
        assert f"Rank tier: {tier}" not in summary, (
            f"absent column still reported a '{tier}' count")


# --- Robustness -------------------------------------------------------------

def test_absent_column_still_produces_all_figure_formats(tmp_path: Path) -> None:
    cols = [c for c in RANKED_COLUMNS if c != "rank_tier"]
    _run_plot(tmp_path, [_row(i, "high") for i in range(5)], columns=cols)
    for ext in ("png", "svg", "pdf"):
        assert (tmp_path / f"ranking_plot.{ext}").exists(), f"missing .{ext}"


def test_blank_and_unknown_tiers_do_not_crash(tmp_path: Path) -> None:
    rows = [_row(0, "high"), _row(1, ""), _row(2, "some_future_tier"),
            _row(3, "tail")]
    summary = _run_plot(tmp_path, rows)
    # Blank and unrecognised both fall into "not scored" -- never binned as tail.
    assert summary["Rank tier: not scored"] == "2"
    assert summary["Rank tier: tail"] == "1"


def test_all_one_tier_does_not_crash(tmp_path: Path) -> None:
    summary = _run_plot(tmp_path, [_row(i, "tail") for i in range(6)])
    assert summary["Rank tier: tail"] == "6"
    assert summary["Rank tier: high"] == "0"


# --- Source-level: palette + labels -----------------------------------------

def _source() -> str:
    return PLOT_RANKING.read_text()


def test_tier_palette_keys_are_the_rank_tier_vocabulary() -> None:
    """The tier colour map must key on the lowercase rank_tier values."""
    src = _source()
    m = re.search(r"^TIER_COLORS\s*=\s*\{(.*?)^\}", src, flags=re.S | re.M)
    assert m, "TIER_COLORS palette not found"
    keys = set(re.findall(r"'([^']+)'\s*:", m.group(1)))
    assert keys == set(TIERS), f"palette keys are {keys}, expected {set(TIERS)}"


def test_tail_is_not_coloured_red() -> None:
    """``tail`` is a rank statement, not a failure; red would misrepresent it.

    Most candidates are necessarily tail (only k fit in the top k), so the
    colour must read as neutral rather than as an error state.
    """
    src = _source()
    m = re.search(r"^TIER_COLORS\s*=\s*\{(.*?)^\}", src, flags=re.S | re.M)
    tail = re.search(r"'tail'\s*:\s*'(#[0-9a-fA-F]{6})'", m.group(1))
    assert tail, "no colour assigned to the tail tier"
    assert tail.group(1).lower() != "#d62728", "tail is still coloured alert red"


def test_plot_no_longer_reads_confidence_tier() -> None:
    """Display switched to rank_tier; the legacy column must not be plotted."""
    assert "confidence_tier" not in _source()


def test_tier_panel_labels_describe_the_topk_probability() -> None:
    """The axis label / title must not still call this a confidence grade."""
    src = _source()
    # The stale labels: a bare "Confidence Tier" axis and the old panel title,
    # both of which read as an evidence grade. "Rank-Confidence Tier" is fine.
    assert not re.search(r"(?<!Rank-)Confidence Tier", src)
    assert "Candidate Confidence Distribution" not in src
    assert "No confidence tier data" not in src
    assert re.search(r"P\(top|top-k|Rank tier|rank confidence", src, flags=re.I), (
        "no panel label mentions the top-k rank-confidence semantics")
