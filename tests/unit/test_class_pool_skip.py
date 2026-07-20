"""Tests for the empty-per-class-pool SKIP behaviour in stage 04 (bead 444).

Stage 04 guarded the per-class P2 pool with::

    assert_fasta_non_empty "${REFS_FA}" "class ${CLASS} P2 pool" || exit 1

inside ``for CLASS in ${GPCR_CLASSES:-A B C F}``. The guard itself is right (a
zero-record pool would otherwise flow through cat/length-filter/MAFFT/ClipKit/
FastTree/IQ-TREE, all of which "succeed" on empty input, and emit a degenerate
tree). The *reaction* was wrong: ``build_per_class_reference_pools.py`` writes
``refs_class_<C>.fa`` for EVERY class unconditionally, so a legitimately empty
class — e.g. no class-F members in this dataset — aborted the whole stage. That
costs the class-A tree, which is the scientific point of this pipeline.

USER DECISION (bead 444): skip the empty class loudly and continue to the next
one. The skip must be unmistakable and greppable in the log, never silent.

``assert_fasta_non_empty`` itself is unchanged — other call sites rely on its
return semantics — so these tests only pin stage 04's reaction to it.
"""
from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
FUNCTIONS_SH = PROJECT_ROOT / "functions.sh"
STAGE04 = PROJECT_ROOT / "04_phylogenetic_analysis.sh"


def _extract_guard_block() -> str:
    """Pull the live empty-pool guard out of stage 04 so we test the real code.

    The block starts at the `if ! assert_fasta_non_empty "${REFS_FA}"` line and
    ends at the first `fi` at the same (4-space, inside the class loop) indent.
    """
    lines = STAGE04.read_text().splitlines()
    start = None
    for i, line in enumerate(lines):
        if line.strip().startswith("if ! assert_fasta_non_empty") and "${REFS_FA}" in line:
            start = i
            break
    if start is None:
        pytest.fail(
            "stage 04 has no `if ! assert_fasta_non_empty \"${REFS_FA}\"` guard block; "
            "the empty per-class pool must be handled by a skip, not `|| exit 1`"
        )
    indent = re.match(r"\s*", lines[start]).group(0)
    for j in range(start + 1, len(lines)):
        if lines[j] == f"{indent}fi":
            return "\n".join(lines[start : j + 1])
    pytest.fail("could not find the closing `fi` of the empty-pool guard block")


def _run_loop_with_guard(tmp_path: Path, fa_by_class: dict[str, Path]) -> subprocess.CompletedProcess:
    """Run the real guard block inside a stand-in for stage 04's class loop."""
    env = os.environ.copy()
    env["LOGS_DIR"] = str(tmp_path)
    assignments = "\n".join(f'FA_{c}="{p}"' for c, p in fa_by_class.items())
    classes = " ".join(fa_by_class)
    script = f"""
source "{FUNCTIONS_SH}"; trap - EXIT
{assignments}
for CLASS in {classes}; do
    eval REFS_FA="\\$FA_${{CLASS}}"
{_extract_guard_block()}
    echo "REACHED_BODY:${{CLASS}}"
done
echo "LOOP_FINISHED"
"""
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=env)


def _pool(tmp_path: Path, cls: str, records: int) -> Path:
    fa = tmp_path / f"refs_class_{cls}.fa"
    fa.write_text("".join(f">{cls}_seq{i}\nMKFLV\n" for i in range(records)))
    return fa


def test_empty_class_is_skipped_and_the_loop_continues(tmp_path):
    """The core user decision: empty F must not cost the A tree."""
    fas = {"F": _pool(tmp_path, "F", 0), "A": _pool(tmp_path, "A", 3)}
    result = _run_loop_with_guard(tmp_path, fas)
    assert result.returncode == 0, result.stderr
    assert "REACHED_BODY:F" not in result.stdout, "empty class F must be skipped"
    assert "REACHED_BODY:A" in result.stdout, "class A must still be built after F is skipped"
    assert "LOOP_FINISHED" in result.stdout, "the stage must not abort on an empty class"


def test_non_empty_classes_all_run(tmp_path):
    fas = {"A": _pool(tmp_path, "A", 2), "C": _pool(tmp_path, "C", 5)}
    result = _run_loop_with_guard(tmp_path, fas)
    assert result.returncode == 0, result.stderr
    assert "REACHED_BODY:A" in result.stdout
    assert "REACHED_BODY:C" in result.stdout


def test_missing_pool_file_is_also_skipped_not_fatal(tmp_path):
    """A pool that was never written is the same non-fatal situation."""
    fas = {"F": tmp_path / "refs_class_F.fa", "A": _pool(tmp_path, "A", 2)}
    result = _run_loop_with_guard(tmp_path, fas)
    assert result.returncode == 0, result.stderr
    assert "REACHED_BODY:F" not in result.stdout
    assert "REACHED_BODY:A" in result.stdout


def test_skip_is_logged_loudly_and_greppably(tmp_path):
    """Not silent: the log must name the class and say no tree is produced."""
    fas = {"F": _pool(tmp_path, "F", 0)}
    result = _run_loop_with_guard(tmp_path, fas)
    combined = result.stdout + result.stderr
    assert "SKIPPING CLASS F" in combined, (
        "the skip must emit an unmistakable, greppable marker naming the class; "
        f"got:\n{combined}"
    )
    assert "no tree" in combined.lower(), \
        "the skip message must state that no tree will be produced for the class"
    # ERROR level goes to stderr (see log() in functions.sh).
    assert "SKIPPING CLASS F" in result.stderr
    # assert_fasta_non_empty's own ERROR must still be there for context.
    assert "0 fasta records" in result.stderr.lower()


def test_skip_reaches_the_pipeline_log_file(tmp_path):
    fas = {"F": _pool(tmp_path, "F", 0)}
    _run_loop_with_guard(tmp_path, fas)
    log_text = (tmp_path / "pipeline.log").read_text()
    assert "SKIPPING CLASS F" in log_text


def test_stage04_no_longer_aborts_on_an_empty_pool():
    """Regression guard on the source: the `|| exit 1` reaction must be gone."""
    text = STAGE04.read_text()
    assert 'assert_fasta_non_empty "${REFS_FA}" "class ${CLASS} P2 pool" || exit 1' not in text, \
        "an empty class pool must skip the class, not abort the stage"
    block = _extract_guard_block()
    assert "continue" in block, "the empty-pool guard must `continue` to the next class"
    assert "exit 1" not in block, "the empty-pool guard must not exit the stage"


def test_stage04_still_parses():
    result = subprocess.run(["bash", "-n", str(STAGE04)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
