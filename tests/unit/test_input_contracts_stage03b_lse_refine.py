"""03b discards lse_refine.py's diagnostics entirely (audit finding F12.4).

``03b:86-91`` ran the classifier as::

    python3 "${SCRIPTS_DIR}/lse_refine.py" "$og" ... --exclude-refs 2>/dev/null || true

``2>/dev/null`` throws away the only channel that says WHY a classification
failed, and ``|| true`` throws away the exit status. Every orthogroup can fail
-- a missing NCBI taxonomy cache, an unparseable id_map, an import error -- and
the stage still reaches ``touch step_completed_lse_classification.txt`` at
:180 and reports "LSE classification completed."

Scope of the fix, deliberately narrow: stderr now goes to a per-orthogroup log
under ``${LOGS_DIR}`` (matching 03a:161's existing
``2>>"${LOGS_DIR}/busco_${base}_drop_gappy.log" || true`` idiom), and failures
are counted and reported. Control flow is UNCHANGED -- a failing orthogroup is
still skipped rather than aborting the stage, and the completion flag is still
written. Making the completion flag conditional on the failure count would
change pipeline behaviour and is reported for the user's decision instead.

These tests execute the REAL orthogroup loop lifted out of the shipping stage
03b script, with a stubbed ``python3`` on PATH standing in for lse_refine.py.
"""
from __future__ import annotations

import os
import stat
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
STAGE03B = PROJECT_ROOT / "03b_lse_classification.sh"


def _dedent(lines: list[str], indent: int) -> str:
    return "\n".join(l[indent:] if len(l) > indent else l for l in lines)


def _find_loop_bounds() -> tuple[list[str], int, int, int]:
    """Locate the `for og in "$OG_DIR"/OG*.fa; do ... done` loop in stage 03b."""
    lines = STAGE03B.read_text().splitlines()
    start = next((i for i, l in enumerate(lines) if l.lstrip().startswith("for og in")), None)
    assert start is not None, "03b no longer loops over orthogroups"
    indent = len(lines[start]) - len(lines[start].lstrip())
    for j in range(start + 1, len(lines)):
        if lines[j].strip() == "done" and (len(lines[j]) - len(lines[j].lstrip())) == indent:
            return lines, start, j, indent
    pytest.fail("could not find the closing `done` of 03b's orthogroup loop")


def _extract_og_loop() -> str:
    """Just the loop body -- used by the static stderr-redirect guards."""
    lines, start, end, indent = _find_loop_bounds()
    return _dedent(lines[start : end + 1], indent)


def _extract_og_loop_with_summary() -> str:
    """The loop PLUS the failure-summary block that follows it, so the executed
    snippet covers everything the fix touches."""
    lines, start, end, indent = _find_loop_bounds()
    tail = end
    for j in range(end + 1, min(end + 20, len(lines))):
        stripped = lines[j].strip()
        if stripped.startswith("if [") and "lse_refine_failures" in stripped:
            for k in range(j + 1, len(lines)):
                if lines[k].strip() == "fi" and (len(lines[k]) - len(lines[k].lstrip())) == indent:
                    tail = k
                    break
            break
        if stripped and not stripped.startswith("#"):
            break
    return _dedent(lines[start : tail + 1], indent)


def _stub_python3(bin_dir: Path, *, exit_code: int, stderr_text: str = "") -> None:
    """Put a fake `python3` on PATH standing in for lse_refine.py."""
    bin_dir.mkdir(parents=True, exist_ok=True)
    stub = bin_dir / "python3"
    stub.write_text(f'#!/bin/bash\n>&2 echo "{stderr_text}"\nexit {exit_code}\n')
    stub.chmod(stub.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)


def _run_loop(tmp_path: Path, *, n_ogs: int, exit_code: int, stderr_text: str = "boom"):
    og_dir = tmp_path / "ogs"
    og_dir.mkdir(parents=True, exist_ok=True)
    for i in range(n_ogs):
        (og_dir / f"OG{i:07d}.fa").write_text(">s1\nMKFLV\n")

    results = tmp_path / "results"
    (results / "lse_classification").mkdir(parents=True, exist_ok=True)
    logs = tmp_path / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    bin_dir = tmp_path / "bin"
    _stub_python3(bin_dir, exit_code=exit_code, stderr_text=stderr_text)

    snippet = "\n".join(
        [
            'log() { echo "LOG $*"; }',
            f'OG_DIR="{og_dir}"',
            f'RESULTS_DIR="{results}"',
            f'LOGS_DIR="{logs}"',
            f'SCRIPTS_DIR="{PROJECT_ROOT / "scripts"}"',
            'ID_MAP_ARG=""',
            'GENE_TREE_ARG=""',
            'SYNTENY_ARG=""',
            'GENE_TREE_DIR=""',
            _extract_og_loop_with_summary(),
        ]
    )
    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}:{env['PATH']}"
    return subprocess.run(["bash", "-c", snippet], capture_output=True, text=True, env=env), logs


def _log_text(logs: Path) -> str:
    return "".join(p.read_text() for p in logs.rglob("*") if p.is_file())


# --------------------------------------------------------------------------
# the regression: diagnostics must survive
# --------------------------------------------------------------------------
def test_lse_refine_stderr_is_captured_to_a_log(tmp_path):
    """The core bug: `2>/dev/null` deleted the only explanation of a failure."""
    proc, logs = _run_loop(tmp_path, n_ogs=2, exit_code=1, stderr_text="TaxIDLookupError")

    assert "TaxIDLookupError" in _log_text(logs), (
        "lse_refine.py's stderr must be recoverable, not discarded:\n" + proc.stdout + proc.stderr
    )


def test_stderr_of_a_successful_run_is_also_kept(tmp_path):
    """Warnings on a zero-exit run are diagnostics too."""
    _, logs = _run_loop(tmp_path, n_ogs=1, exit_code=0, stderr_text="WARN partial lineage")

    assert "WARN partial lineage" in _log_text(logs)


def test_failures_are_counted_and_reported(tmp_path):
    """A total wipe-out must be visible in the stage log, not silent."""
    proc, _ = _run_loop(tmp_path, n_ogs=3, exit_code=1)

    combined = proc.stdout + proc.stderr
    assert "3" in combined, "the number of failed orthogroups must be reported: " + combined
    assert "LOG" in combined


def test_all_success_reports_no_failures(tmp_path):
    """The failure report must discriminate, not cry wolf on every run."""
    proc, _ = _run_loop(tmp_path, n_ogs=3, exit_code=0)

    combined = proc.stdout + proc.stderr
    assert "fail" not in combined.lower() or "0" in combined, combined


# --------------------------------------------------------------------------
# control flow must be UNCHANGED (widening it needs the user's approval)
# --------------------------------------------------------------------------
def test_a_failing_orthogroup_does_not_abort_the_loop(tmp_path):
    """Skip-and-continue is existing behaviour and is preserved by this fix."""
    proc, _ = _run_loop(tmp_path, n_ogs=4, exit_code=1)

    assert proc.returncode == 0, (
        "the fix must not turn a per-orthogroup failure into a stage abort "
        "(that is a behaviour change requiring approval):\n" + proc.stdout + proc.stderr
    )


def test_completion_flag_line_is_untouched():
    """Gating the flag on the failure count is a behaviour change -- reported, not made."""
    src = STAGE03B.read_text()
    assert 'touch "${RESULTS_DIR}/step_completed_lse_classification.txt"' in src


# --------------------------------------------------------------------------
# static guards
# --------------------------------------------------------------------------
def test_no_stderr_discard_remains_in_the_orthogroup_loop():
    loop = _extract_og_loop()
    assert "2>/dev/null" not in loop, "stderr suppression re-hides every lse_refine.py failure"


def test_stderr_is_redirected_into_logs_dir():
    loop = _extract_og_loop()
    assert "LOGS_DIR" in loop, "lse_refine.py's stderr must land under ${LOGS_DIR}"
