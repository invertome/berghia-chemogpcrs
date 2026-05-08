"""Tests for scripts/drop_near_all_gap_rows.py.

The helper drops alignment rows whose gap fraction is at or above a
threshold. Used as a defensive filter between ClipKit and IQ-TREE 3 to
prevent "Unknown sequence type" errors caused by near-all-gap rows
confusing IQ-TREE's alphabet auto-detection.
"""
from __future__ import annotations

from pathlib import Path

import pytest

import drop_near_all_gap_rows as dnagr


# ---- gap_fraction --------------------------------------------------------

def test_gap_fraction_all_residues() -> None:
    assert dnagr.gap_fraction("ACDEFG") == 0.0


def test_gap_fraction_all_gaps() -> None:
    assert dnagr.gap_fraction("------") == 1.0


def test_gap_fraction_mixed() -> None:
    assert dnagr.gap_fraction("AC--EF") == pytest.approx(2 / 6)


def test_gap_fraction_empty_seq_treated_as_fully_gappy() -> None:
    """An empty sequence has no information — treat as fully gappy so
    it gets dropped by any sensible threshold."""
    assert dnagr.gap_fraction("") == 1.0


# ---- filter_alignment ----------------------------------------------------

def test_filter_keeps_clean_rows(tmp_path: Path) -> None:
    """Rows below the threshold are kept verbatim."""
    inp = tmp_path / "in.fa"
    inp.write_text(">a\nACDEFG\n>b\nACDEFG\n")
    out = tmp_path / "out.fa"
    n_kept, n_dropped = dnagr.filter_alignment(str(inp), str(out), 0.95)
    assert n_kept == 2
    assert n_dropped == 0
    assert out.read_text() == ">a\nACDEFG\n>b\nACDEFG\n"


def test_filter_drops_near_all_gap_rows(tmp_path: Path) -> None:
    """Rows at or above the threshold are dropped. Default threshold
    is 0.95, so a 19/20-gap row counts."""
    inp = tmp_path / "in.fa"
    inp.write_text(">good\n" + "A" * 20 + "\n>gappy\n" + "-" * 19 + "A\n")
    out = tmp_path / "out.fa"
    n_kept, n_dropped = dnagr.filter_alignment(str(inp), str(out), 0.95)
    assert n_kept == 1
    assert n_dropped == 1
    text = out.read_text()
    assert ">good" in text
    assert ">gappy" not in text


def test_filter_handles_multiline_sequences(tmp_path: Path) -> None:
    """FASTA bodies wrapped across multiple lines must be re-joined
    before computing gap fraction."""
    inp = tmp_path / "in.fa"
    inp.write_text(">a\nAC\nDE\nFG\n>b\n--\n--\n--\n")
    out = tmp_path / "out.fa"
    n_kept, n_dropped = dnagr.filter_alignment(str(inp), str(out), 0.95)
    assert n_kept == 1
    assert n_dropped == 1
    assert out.read_text().startswith(">a\n")


def test_filter_threshold_strict_inequality(tmp_path: Path) -> None:
    """Drop when fraction >= threshold, keep when fraction < threshold.
    A row at exactly the threshold is dropped (defensive: the IQ-TREE
    auto-detect failure happens at the boundary)."""
    inp = tmp_path / "in.fa"
    # Exactly 95 % gap: 19 gaps + 1 residue
    inp.write_text(">borderline\n" + "-" * 19 + "A\n")
    out = tmp_path / "out.fa"
    n_kept, n_dropped = dnagr.filter_alignment(str(inp), str(out), 0.95)
    assert n_kept == 0
    assert n_dropped == 1


def test_filter_drops_all_writes_empty_output(tmp_path: Path) -> None:
    """If every sequence is near-all-gap, the output is an empty (but
    valid) FASTA. Downstream code can decide what to do."""
    inp = tmp_path / "in.fa"
    inp.write_text(">x\n----\n>y\n----\n")
    out = tmp_path / "out.fa"
    n_kept, n_dropped = dnagr.filter_alignment(str(inp), str(out), 0.95)
    assert n_kept == 0
    assert n_dropped == 2
    assert out.read_text() == ""


def test_filter_creates_parent_directory(tmp_path: Path) -> None:
    """Output parent dir is created if missing — matches the Path
    conventions used by the rest of the pipeline."""
    inp = tmp_path / "in.fa"
    inp.write_text(">a\nACDE\n")
    out = tmp_path / "deep" / "nested" / "out.fa"
    dnagr.filter_alignment(str(inp), str(out), 0.95)
    assert out.exists()
