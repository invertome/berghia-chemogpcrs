"""Tests for count_sequences in functions.sh (the bead-444 `grep -c` idiom).

``count_sequences`` was written as::

    count_sequences() {
        local fasta="$1"
        grep -c "^>" "$fasta" 2>/dev/null || echo 0
    }

``grep -c`` prints ``0`` *and* exits 1 when nothing matches, so the ``|| echo 0``
fallback appends a SECOND ``0``. The helper then emits the two-line string
``"0\\n0"`` instead of one integer.

This is NOT a cosmetic helper. It has six callers inside functions.sh, and three
of them feed the value straight into shell arithmetic, where ``"0\\n0"`` is a
hard failure rather than a wrong number::

    functions.sh:1296  total_seqs=$((total_seqs + $(count_sequences "$f")))
    functions.sh:1325  local num_branches=$((2 * num_seqs - 2))
    functions.sh:1440  total_seqs=$((total_seqs + $(count_sequences "$f")))

  bash: t + 0
  0: arithmetic syntax error in expression (error token is "0")

The remaining callers (estimate_memory_for_alignment, estimate_memory_for_tree,
estimate_memory_for_hyphy) pass it to ``awk -v n=...``, which silently truncates
at the newline -- a wrong memory estimate rather than a crash.

The no-match path is reachable in normal operation: any zero-length or
header-less FASTA produced by an upstream filter (PREQUAL/CLOAK/TAPER dropping
every sequence, an empty orthogroup, a failed extraction) hits it.

The fix is the established in-repo pattern from ``count_taxid_occurrences``:
capture with ``|| true`` so the count is kept without the exit status, then
emit ``${n:-0}`` so a missing file still yields exactly one integer.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
FUNCTIONS_SH = PROJECT_ROOT / "functions.sh"


def _run_bash(body: str, tmp_path: Path) -> subprocess.CompletedProcess:
    """Source functions.sh and run *body*.

    Sourcing functions.sh installs an EXIT trap (finalize_pipeline) that logs to
    stdout; `trap - EXIT` disarms it so only the helper's output is captured
    (same technique as tests/unit/test_count_taxid_occurrences.py).
    """
    env = os.environ.copy()
    env["LOGS_DIR"] = str(tmp_path)
    script = f'source "{FUNCTIONS_SH}"; trap - EXIT; {body}'
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=env)


def _fasta(tmp_path: Path, name: str, contents: str) -> Path:
    path = tmp_path / name
    path.write_text(contents)
    return path


# --- the counted value itself -------------------------------------------------


def test_counts_sequences_in_a_normal_fasta(tmp_path):
    fa = _fasta(tmp_path, "two.fa", ">a\nACGT\n>b\nACGT\n")
    result = _run_bash(f'count_sequences "{fa}"', tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "2"


def test_headerless_fasta_emits_a_single_zero_line(tmp_path):
    """The defect: `grep -c` prints 0 AND exits 1, so `|| echo 0` doubled it."""
    fa = _fasta(tmp_path, "headerless.fa", "ACGT\nTTTT\n")
    result = _run_bash(f'count_sequences "{fa}"', tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout.split() == ["0"], f"expected one '0' line, got {result.stdout!r}"


def test_empty_file_emits_a_single_zero_line(tmp_path):
    fa = _fasta(tmp_path, "empty.fa", "")
    result = _run_bash(f'count_sequences "{fa}"', tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout.split() == ["0"], f"expected one '0' line, got {result.stdout!r}"


def test_missing_file_emits_a_single_zero_line(tmp_path):
    missing = tmp_path / "does_not_exist.fa"
    result = _run_bash(f'count_sequences "{missing}"', tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout.split() == ["0"], f"expected one '0' line, got {result.stdout!r}"


# --- what the callers actually do with it -------------------------------------


def test_result_survives_an_integer_test(tmp_path):
    fa = _fasta(tmp_path, "headerless.fa", "ACGT\n")
    result = _run_bash(
        f'n=$(count_sequences "{fa}"); if [ "$n" -gt 0 ]; then echo YES; else echo NO; fi',
        tmp_path,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "NO"
    assert "integer expression expected" not in result.stderr
    assert "integer expected" not in result.stderr


def test_result_survives_shell_arithmetic(tmp_path):
    """functions.sh:1296 and :1440 do `$((total + $(count_sequences "$f")))`."""
    fa = _fasta(tmp_path, "headerless.fa", "ACGT\n")
    result = _run_bash(
        f'total=0; total=$((total + $(count_sequences "{fa}"))); echo "total=$total"',
        tmp_path,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "total=0"
    assert "arithmetic syntax error" not in result.stderr


def test_get_dataset_stats_handles_a_directory_of_headerless_fastas(tmp_path):
    """Real caller: get_dataset_stats sums count_sequences over a directory."""
    d = tmp_path / "seqs"
    d.mkdir()
    (d / "a.fa").write_text("ACGT\n")  # no headers -> the no-match path
    (d / "b.fa").write_text(">x\nACGT\n")
    result = _run_bash(f'get_dataset_stats "{d}"', tmp_path)
    assert result.returncode == 0, result.stderr
    assert "arithmetic syntax error" not in result.stderr
    assert "Total sequences: 1" in result.stdout, result.stdout


def test_estimate_memory_for_hyphy_handles_a_headerless_alignment(tmp_path):
    """Real caller: `$((2 * num_seqs - 2))` on functions.sh:1325."""
    aln = _fasta(tmp_path, "aln.fa", "ACGT\n")
    result = _run_bash(f'estimate_memory_for_hyphy "{aln}"', tmp_path)
    assert result.returncode == 0, result.stderr
    assert "arithmetic syntax error" not in result.stderr
    assert result.stdout.strip().endswith("G"), result.stdout


# --- regression guard ---------------------------------------------------------


def test_count_sequences_no_longer_uses_the_doubling_idiom():
    text = FUNCTIONS_SH.read_text()
    start = text.index("count_sequences() {")
    body = text[start : text.index("}", start)]
    assert "|| echo 0" not in body, \
        "count_sequences still uses the `|| echo 0` idiom that doubles the zero"


def test_functions_sh_still_parses():
    result = subprocess.run(["bash", "-n", str(FUNCTIONS_SH)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
