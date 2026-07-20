"""Tests for the ECL-divergence producer wiring (bead 444, change 1).

Historically there were TWO producers writing
``${RESULTS_DIR}/ecl_analysis/ecl_divergence.csv``:

  * a dead "Phase 3" block in ``07_candidate_ranking.sh`` that guarded on a
    DeepTMHMM directory stage 02 never writes and globbed an alignment path
    stage 04 never writes (so it never ran), and
  * ``04b_ecl_analysis.sh``, which resolves the right DeepTMHMM directory and
    the canonical un-CLOAK-masked alignment — but was never registered in
    ``run_pipeline.sh``, so it never ran either.

Net effect: ``ecl_divergence.csv`` was never produced and the
``ecl_divergence`` ranking signal (``ECL_DIVERGENCE_WEIGHT``) was permanently
dormant.

These tests pin the fix: 04b is the SOLE producer and is registered in the
pipeline between its stage-04 input and its stage-07 consumer; stage 07 is a
pure CONSUMER that no longer tries to produce the file itself.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

STAGE_07 = REPO_ROOT / "07_candidate_ranking.sh"
STAGE_04B = REPO_ROOT / "04b_ecl_analysis.sh"
RUN_PIPELINE = REPO_ROOT / "run_pipeline.sh"


def _dry_run_steps(*args: str) -> list[str]:
    """Run run_pipeline.sh in dry-run mode and return the submitted step ids.

    The submission lines look like ``[04b] Submitting: ECL Divergence Analysis``.

    The exit code is deliberately NOT asserted: run_pipeline.sh's `log` helper
    appends to ${LOGS_DIR}/pipeline.log, which does not exist in a fresh
    checkout (results/ is generated), and under `set -e` that makes the script
    exit 1 even on a fully successful dry run. That is pre-existing and
    orthogonal to step registration, which is what these tests cover — so we
    assert on the submission lines the dry run emits instead.
    """
    proc = subprocess.run(
        ["./run_pipeline.sh", "--dry-run", "--no-validate", *args],
        cwd=REPO_ROOT, capture_output=True, text=True, timeout=120,
    )
    combined = proc.stdout + proc.stderr
    steps = re.findall(r"^\[([0-9a-z]+)\] Submitting:", combined, flags=re.M)
    assert steps, f"dry run submitted nothing:\n{combined}"
    return steps


# --- Change 1a: stage 07 no longer produces ECL divergence -------------------

def test_stage_07_does_not_invoke_analyze_ecl() -> None:
    """The dead Phase 3 producer block is gone; 07 only CONSUMES the CSV."""
    text = STAGE_07.read_text()
    assert "analyze_ecl.py" not in text, (
        "07_candidate_ranking.sh still invokes analyze_ecl.py — 04b is the "
        "sole producer of ecl_divergence.csv"
    )


def test_stage_07_has_no_stale_ecl_paths() -> None:
    """The two stale paths that made the old block unreachable must be gone.

    ``${RESULTS_DIR}/deeptmhmm`` (stage 02 writes chemogpcrs/deeptmhmm_berghia)
    and the ``aligned_*.fasta`` glob (stage 04 writes class_A/class_A_aligned.fa).
    """
    text = STAGE_07.read_text()
    assert '${RESULTS_DIR}/deeptmhmm' not in text, "stale DeepTMHMM path remains"
    assert "aligned_*.fasta" not in text, "stale alignment glob remains"


def test_stage_07_still_creates_ecl_output_dir_and_neighbours_intact() -> None:
    """Surgical deletion: the surrounding G-protein block and the
    candidate_ids.txt extraction must be untouched."""
    text = STAGE_07.read_text()
    # neighbour above
    assert "classify_gproteins.py" in text
    assert "coexpression_analysis.py" in text
    # neighbour below
    assert "ranking/candidate_ids.txt" in text
    # the mkdir that provisions ecl_analysis/ is shared with other dirs
    assert 'mkdir -p "${RESULTS_DIR}/expression"' in text


def test_stage_07_is_syntactically_valid() -> None:
    proc = subprocess.run(["bash", "-n", str(STAGE_07)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr


# --- Change 1b: 04b registered in run_pipeline.sh ----------------------------

def test_04b_script_exists_and_is_the_producer() -> None:
    assert STAGE_04B.is_file()
    text = STAGE_04B.read_text()
    assert "analyze_ecl.py" in text
    assert "ecl_analysis/ecl_divergence.csv" in text


def test_run_pipeline_registers_04b_with_expected_format() -> None:
    """PIPELINE_STEPS entry must follow script:description:is_array exactly."""
    text = RUN_PIPELINE.read_text()
    m = re.search(r'^\s*\["04b"\]="([^"]+)"\s*$', text, flags=re.M)
    assert m, "no PIPELINE_STEPS entry for 04b"
    parts = m.group(1).split(":")
    assert len(parts) == 3, f"expected script:description:is_array, got {parts}"
    script, description, is_array = parts
    assert script == "04b_ecl_analysis.sh"
    assert description.strip(), "description must be non-empty"
    assert is_array == "no", "04b is not a SLURM array job"


def test_run_pipeline_step_order_places_04b_after_04_and_before_07() -> None:
    """04b consumes stage 04's alignment and must produce before 07 consumes."""
    text = RUN_PIPELINE.read_text()
    m = re.search(r"^STEP_ORDER=\((.*?)\)\s*$", text, flags=re.M | re.S)
    assert m, "STEP_ORDER not found"
    order = re.findall(r'"([^"]+)"', m.group(1))
    assert "04b" in order, f"04b missing from STEP_ORDER: {order}"
    assert order.index("04") < order.index("04b") < order.index("07")


def test_dry_run_all_submits_04b_between_04_and_07() -> None:
    """Behavioural check: `run_pipeline.sh all` actually schedules 04b."""
    steps = _dry_run_steps("all")
    assert "04b" in steps, f"04b never submitted by `all`: {steps}"
    assert steps.index("04") < steps.index("04b") < steps.index("07")


def test_dry_run_accepts_04b_as_an_explicit_step() -> None:
    steps = _dry_run_steps("04b")
    assert steps == ["04b"], steps


def test_dry_run_range_04_to_07_includes_04b() -> None:
    """Ranges expand via STEP_ORDER, so 04-07 must pick up 04b."""
    steps = _dry_run_steps("04-07")
    assert "04b" in steps, f"range 04-07 skipped 04b: {steps}"


def test_run_pipeline_is_syntactically_valid() -> None:
    proc = subprocess.run(["bash", "-n", str(RUN_PIPELINE)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
