"""Tests for the duplicate-FASTA-id guard in functions.sh.

``assert_no_duplicate_fasta_ids <fasta> [label]`` must fail loud (exit 1)
when any sequence ID (the first whitespace-delimited header token) occurs
more than once. This is the guard that prevents the stage-04 per-class
combined pool from carrying duplicate Berghia records into IQ-TREE, which
aborts on duplicate taxon names. See bead berghia-chemogpcrs-3sd.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
FUNCTIONS_SH = PROJECT_ROOT / "functions.sh"


def _run_guard(fasta: Path, tmp_path: Path, label: str = "") -> subprocess.CompletedProcess:
    env = os.environ.copy()
    env["LOGS_DIR"] = str(tmp_path)  # log() appends to ${LOGS_DIR}/pipeline.log
    label_arg = f' "{label}"' if label else ""
    script = f'source "{FUNCTIONS_SH}"; assert_no_duplicate_fasta_ids "{fasta}"{label_arg}'
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=env)


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> Path:
    with path.open("w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")
    return path


def test_unique_ids_pass(tmp_path):
    fa = _write_fasta(
        tmp_path / "clean.fa",
        [("ref_1", "ACDEFG"), ("ref_2", "GHIKLM"), ("1287507_bs_1", "MNPQRS")],
    )
    result = _run_guard(fa, tmp_path)
    assert result.returncode == 0, result.stderr


def test_duplicate_id_fails_loud(tmp_path):
    fa = _write_fasta(
        tmp_path / "dup.fa",
        [("ref_1", "ACDEFG"), ("1287507_bs_1", "GHIKLM"), ("1287507_bs_1", "MNPQRS")],
    )
    result = _run_guard(fa, tmp_path)
    assert result.returncode == 1
    assert "1287507_bs_1" in result.stderr
    assert "uplicate" in result.stderr  # "Duplicate" / "duplicate"


def test_duplicate_id_ignores_header_description(tmp_path):
    # Same ID token, different free-text description after whitespace → still
    # a duplicate ID (IQ-TREE keys on the ID token, not the full header).
    fa = _write_fasta(
        tmp_path / "dup_desc.fa",
        [("ref_1 class A", "ACDEFG"), ("ref_1 class B", "GHIKLM")],
    )
    result = _run_guard(fa, tmp_path)
    assert result.returncode == 1
    assert "ref_1" in result.stderr


def test_missing_file_fails(tmp_path):
    result = _run_guard(tmp_path / "nope.fa", tmp_path)
    assert result.returncode == 1
