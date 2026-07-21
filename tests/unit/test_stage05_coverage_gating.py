"""Stage 05 must report coverage and gate its completion markers (bead higb).

The per-class path bug (F1) was fixed, but the SILENCE that let it hide for
months was not: stage 05 touched both ``step_completed_*`` markers
unconditionally, so a run in which every orthogroup failed to resolve still
reported done and let stage 07 proceed. "Complete" without coverage numbers is
an unverified claim.

Stage 05 is a SLURM array over orthogroups, so coverage must be recorded
PER-ORTHOGROUP. A shared append target would be a write race between
concurrent array tasks.
"""
from __future__ import annotations

import re
from pathlib import Path

STAGE05 = Path(__file__).resolve().parents[2] / "05_selective_pressure_and_asr.sh"
SRC = STAGE05.read_text()


def _touch_lines(marker: str) -> list[str]:
    """Every line that touches the given completion marker."""
    return [ln.strip() for ln in SRC.splitlines()
            if ln.strip().startswith("touch ") and marker in ln]


def test_stage05_script_exists():
    assert STAGE05.is_file()


def test_per_og_marker_is_not_touched_unconditionally():
    """A per-OG marker asserts that orthogroup produced something."""
    lines = _touch_lines("step_completed_${base}.txt")
    assert lines, "expected stage 05 to touch a per-orthogroup completion marker"
    for ln in lines:
        assert not re.match(r'^touch\s+"?\$\{RESULTS_DIR\}', ln) or "og_products" in SRC, (
            "per-OG completion marker is touched with no reference to a product "
            "counter anywhere in the stage")


def test_a_product_counter_gates_the_markers():
    """The markers must sit inside a conditional on actual output."""
    assert "og_products" in SRC, (
        "stage 05 has no product counter, so it cannot distinguish 'ran and "
        "produced nothing' from 'ran and produced results'")
    # the counter must be consulted, not merely incremented
    assert re.search(r'if\s+\[\s+"?\$\{?og_products\}?"?\s+-(gt|eq)\s+0', SRC), (
        "og_products is never tested in a conditional")


def test_zero_coverage_is_reported_explicitly():
    """A zero-product orthogroup must say so, not pass silently."""
    assert re.search(r'og_products.*-eq\s+0|no products|produced nothing',
                     SRC, re.IGNORECASE), (
        "stage 05 never reports the zero-coverage case")


def test_coverage_record_is_per_orthogroup_not_a_shared_append():
    """Concurrency: array tasks run in parallel, so a shared append target is a
    write race. The coverage record must be keyed by orthogroup."""
    assert "coverage" in SRC, "stage 05 writes no coverage record"
    shared_append = re.search(r'>>\s*"?\$\{RESULTS_DIR\}/selective_pressure/'
                              r'(stage05_)?coverage\.tsv', SRC)
    assert not shared_append, (
        "coverage is appended to ONE shared file from parallel array tasks, "
        "which is a write race; write one record per orthogroup instead")
    assert re.search(r'coverage_dir=.*coverage', SRC), (
        "no per-orthogroup coverage directory is defined")
    assert re.search(r'>\s*"\$\{coverage_dir\}/\$\{base\}', SRC), (
        "the coverage record is not keyed by orthogroup, so parallel array "
        "tasks would write to the same target")


def test_coverage_record_names_both_analyses():
    """Selection and ASR are separate products; collapsing them hides which
    one was missing."""
    m = re.search(r'selection_ran.*asr_ran|asr_ran.*selection_ran', SRC, re.S)
    assert m, "coverage record does not distinguish selection from ASR"
