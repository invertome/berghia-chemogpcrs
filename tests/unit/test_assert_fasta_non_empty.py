"""Tests for assert_fasta_non_empty in functions.sh (bead 444, item 3).

Stage 04 guards the per-class P2 pool with a bare existence test::

    REFS_FA="${PER_CLASS_POOL_DIR}/refs_class_${CLASS}.fa"
    if [ ! -f "${REFS_FA}" ]; then
        log "ERROR: missing P2 pool ${REFS_FA} — run P2 first"
        exit 1
    fi

`[ -f ]` is satisfied by a zero-record file. build_per_class_reference_pools.py
creates one pool file per class, so a class with no members yields a present-but-
empty pool. That file then flows into `cat REFS_FA OUTGROUP_FA > COMBINED`, the
length filter, MAFFT, ClipKit, FastTree and IQ-TREE — every one of which
"succeeds" on empty or near-empty input — and the stage ends up emitting a
degenerate tree for that class with nothing in the log but a soft
`0 Berghia` warning.

`assert_fasta_non_empty` makes the empty-pool case fail loudly at the top of the
per-class loop, matching the existing `assert_no_duplicate_fasta_ids` /
`assert_array_covers_manifest` guard pattern.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
FUNCTIONS_SH = PROJECT_ROOT / "functions.sh"
STAGE04 = PROJECT_ROOT / "04_phylogenetic_analysis.sh"


def _assert_non_empty(fasta: Path, tmp_path: Path, label: str = "test pool"):
    env = os.environ.copy()
    env["LOGS_DIR"] = str(tmp_path)
    script = (
        f'source "{FUNCTIONS_SH}"; trap - EXIT; '
        f'assert_fasta_non_empty "{fasta}" "{label}"'
    )
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=env)


def test_pool_with_records_passes(tmp_path):
    fa = tmp_path / "refs_class_A.fa"
    fa.write_text(">seq1\nMKF\n>seq2\nMKV\n")
    result = _assert_non_empty(fa, tmp_path)
    assert result.returncode == 0, result.stderr


def test_single_record_pool_passes(tmp_path):
    fa = tmp_path / "refs_class_A.fa"
    fa.write_text(">only\nMKF\n")
    result = _assert_non_empty(fa, tmp_path)
    assert result.returncode == 0, result.stderr


def test_zero_byte_pool_fails_loud(tmp_path):
    # The exact silent-pass case: `[ -f ]` is true, but there are no records.
    fa = tmp_path / "refs_class_C.fa"
    fa.touch()
    result = _assert_non_empty(fa, tmp_path, "class C pool")
    assert result.returncode == 1
    assert "0 fasta records" in result.stderr.lower()
    assert "class C pool" in result.stderr


def test_pool_with_content_but_no_headers_fails_loud(tmp_path):
    # Truncated/garbled pool: bytes present, no '>' records.
    fa = tmp_path / "refs_class_F.fa"
    fa.write_text("MKFMKVMKL\n")
    result = _assert_non_empty(fa, tmp_path, "class F pool")
    assert result.returncode == 1
    assert "0 fasta records" in result.stderr.lower()


def test_missing_pool_fails_loud(tmp_path):
    result = _assert_non_empty(tmp_path / "nope.fa", tmp_path, "class B pool")
    assert result.returncode == 1
    assert "class B pool" in result.stderr


def test_stage04_guards_the_per_class_pool():
    """Regression guard: stage 04 must assert the P2 pool is non-empty."""
    text = STAGE04.read_text()
    assert "assert_fasta_non_empty" in text, \
        "stage 04 must call assert_fasta_non_empty on the per-class P2 pool"
    # The guard has to cover REFS_FA (the pool), not some other FASTA.
    assert 'assert_fasta_non_empty "${REFS_FA}"' in text, \
        "the empty-pool guard must be applied to ${REFS_FA}"


def test_stage04_still_parses():
    result = subprocess.run(["bash", "-n", str(STAGE04)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
