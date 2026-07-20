"""Tests for three stage-registration / producer-wiring defects (F2, F5, F7).

All three share the ECL bug's signature: a correct component exists but is
never reached, so a scored axis or an entire data track is silently dormant.

F2 -- the CAFE expansion axis is dormant.
    ``scripts/rank_candidates.py::load_cafe_expansion`` reads
    ``results/cafe/expansion_interpretation.csv``. ``scripts/interpret_cafe.py``
    is the intended producer (its output columns match the reader's keys
    exactly) but was invoked by no stage, so 03c ran the full CAFE5
    computation only for the result to be discarded.

F5 -- ``02c_genome_reconcile.sh`` was not registered in ``run_pipeline.sh`` at
    all. ``RUN_GENOME_TRACK`` defaults to 1, but the only producer of
    ``reconciliation/reconciled_candidates.faa`` never ran, so
    ``berghia_candidate_fasta()`` silently fell back to the transcriptome-only
    set and the entire genome track was absent from orthology onward.

F7 -- ``STEP_ORDER`` placed 03d BEFORE 04, but 03d hard-requires
    ``step_completed_04.txt`` via ``check_file`` (which exits). So
    ``./run_pipeline.sh all`` always aborted at 03d and the registry order
    could not execute end to end.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

RUN_PIPELINE = REPO_ROOT / "run_pipeline.sh"
STAGE_02C = REPO_ROOT / "02c_genome_reconcile.sh"
STAGE_03C = REPO_ROOT / "03c_cafe_analysis.sh"
INTERPRET_CAFE = REPO_ROOT / "scripts" / "interpret_cafe.py"
RANK_CANDIDATES = REPO_ROOT / "scripts" / "rank_candidates.py"


# --- helpers ----------------------------------------------------------------


def _dry_run_steps(*args: str) -> list[str]:
    """Run run_pipeline.sh in dry-run mode and return the submitted step ids.

    Submission lines look like ``[02c] Submitting: Genome-Track Reconciliation``.

    The exit code is deliberately NOT asserted: run_pipeline.sh's `log` helper
    appends to ${LOGS_DIR}/pipeline.log, which does not exist in a fresh
    checkout, and under `set -e` that makes the script exit 1 even on a fully
    successful dry run. That is pre-existing and orthogonal to registration.
    """
    proc = subprocess.run(
        ["./run_pipeline.sh", "--dry-run", "--no-validate", *args],
        cwd=REPO_ROOT, capture_output=True, text=True, timeout=120,
    )
    combined = proc.stdout + proc.stderr
    steps = re.findall(r"^\[([0-9a-z]+)\] Submitting:", combined, flags=re.M)
    assert steps, f"dry run submitted nothing:\n{combined}"
    return steps


def _step_order() -> list[str]:
    m = re.search(r"^STEP_ORDER=\((.*?)\)\s*$", RUN_PIPELINE.read_text(),
                  flags=re.M | re.S)
    assert m, "STEP_ORDER not found"
    return re.findall(r'"([^"]+)"', m.group(1))


def _pipeline_step_keys() -> list[str]:
    text = RUN_PIPELINE.read_text()
    m = re.search(r"declare -A PIPELINE_STEPS=\((.*?)^\)\s*$", text,
                  flags=re.M | re.S)
    assert m, "PIPELINE_STEPS not found"
    return re.findall(r'^\s*\["([^"]+)"\]=', m.group(1), flags=re.M)


def _pipeline_step_entry(step: str) -> tuple[str, str, str]:
    """Return (script, description, is_array) for a PIPELINE_STEPS entry."""
    m = re.search(rf'^\s*\["{re.escape(step)}"\]="([^"]+)"\s*$',
                  RUN_PIPELINE.read_text(), flags=re.M)
    assert m, f"no PIPELINE_STEPS entry for {step}"
    parts = m.group(1).split(":")
    assert len(parts) == 3, f"expected script:description:is_array, got {parts}"
    return parts[0], parts[1], parts[2]


def _interpret_cafe_block() -> str:
    """The 03c shell lines that invoke interpret_cafe.py (incl. continuations)."""
    text = STAGE_03C.read_text()
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.lstrip().startswith("#"):
            continue  # doc comments mention the script but do not invoke it
        if "interpret_cafe.py" in line:
            # walk back over backslash continuations to the start of the command
            start = i
            while start > 0 and lines[start - 1].rstrip().endswith("\\"):
                start -= 1
            end = i
            while end < len(lines) - 1 and lines[end].rstrip().endswith("\\"):
                end += 1
            return "\n".join(lines[start:end + 1])
    return ""


def _interpret_cafe_required_flags() -> list[str]:
    """Required flags straight from interpret_cafe.py's own argparse parser.

    Asserting against the parser itself (rather than a hand-copied list) means
    the test keeps tracking the real contract if the parser changes.
    """
    proc = subprocess.run(
        ["python3", str(INTERPRET_CAFE)],
        cwd=REPO_ROOT, capture_output=True, text=True, timeout=60,
    )
    stderr = proc.stderr
    if "the following arguments are required" not in stderr:
        pytest.skip(f"could not introspect interpret_cafe.py parser: {stderr}")
    tail = stderr.split("the following arguments are required:", 1)[1]
    flags = re.findall(r"--[a-z0-9-]+", tail)
    assert flags, f"no required flags parsed from: {stderr}"
    return flags


# === F2: 03c must invoke interpret_cafe.py ==================================


def test_interpret_cafe_columns_match_the_ranking_reader() -> None:
    """The producer's columns are exactly the keys load_cafe_expansion reads."""
    producer = INTERPRET_CAFE.read_text()
    reader = RANK_CANDIDATES.read_text()
    for column in ("orthogroup", "expansion_pvalue", "taxonomic_level",
                   "expansion_fold", "berghia_copies", "ancestral_copies"):
        assert f"'{column}'" in producer, f"{column} missing from producer"
        assert f"'{column}'" in reader, f"{column} missing from reader"


def test_stage_03c_invokes_interpret_cafe() -> None:
    """F2 regression: the producer must actually be reached by a stage."""
    assert "interpret_cafe.py" in STAGE_03C.read_text(), (
        "03c_cafe_analysis.sh does not invoke scripts/interpret_cafe.py — the "
        "expansion ranking axis stays dormant"
    )


def test_stage_03c_passes_every_required_interpret_cafe_argument() -> None:
    """Every flag interpret_cafe.py's parser marks required must be passed."""
    block = _interpret_cafe_block()
    assert block, "no interpret_cafe.py invocation found in 03c"
    for flag in _interpret_cafe_required_flags():
        assert flag in block, (
            f"03c's interpret_cafe.py call omits required {flag}:\n{block}"
        )


def test_stage_03c_writes_the_csv_the_ranking_stage_reads() -> None:
    """Producer output path must resolve to the path load_cafe_expansion opens.

    The ``--output`` value is expanded one level against 03c's own shell
    assignments (``EXPANSION_CSV=...``, ``CAFE_DIR=...``) so the test pins the
    resolved path, not whichever variable happens to be used.
    """
    block = _interpret_cafe_block()
    m = re.search(r'--output\s+"?([^"\s\\]+)"?', block)
    assert m, f"no --output argument in the interpret_cafe.py call:\n{block}"

    stage = STAGE_03C.read_text()
    value = m.group(1)
    for _ in range(5):  # expand ${VAR} references defined in the stage itself
        ref = re.search(r"\$\{([A-Z_]+)\}", value)
        if not ref:
            break
        assign = re.search(rf'^{ref.group(1)}="([^"]+)"', stage, flags=re.M)
        if not assign:
            break
        value = value.replace(ref.group(0), assign.group(1))

    assert value.endswith("/cafe/expansion_interpretation.csv"), (
        f"--output resolves to {value!r}, not results/cafe/"
        "expansion_interpretation.csv"
    )
    # and that is exactly what the ranking stage opens
    assert "'cafe', 'expansion_interpretation.csv'" in RANK_CANDIDATES.read_text()


def test_stage_03c_interpret_cafe_uses_cafe_results_and_ultrametric_tree() -> None:
    """It must be fed CAFE5's own output dir and the tree CAFE5 actually ran on."""
    block = _interpret_cafe_block()
    assert "CAFE_OUTPUT" in block or "cafe_results" in block, (
        "--cafe-output must point at CAFE5's output directory"
    )
    assert "ULTRAMETRIC_TREE" in block, (
        "--species-tree must be the ultrametric tree CAFE5 was run against, so "
        "branch names in the CAFE output resolve"
    )


def test_stage_03c_is_syntactically_valid() -> None:
    proc = subprocess.run(["bash", "-n", str(STAGE_03C)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr


# === F5: 02c must be registered =============================================


def test_02c_script_exists_and_is_the_genome_track_producer() -> None:
    assert STAGE_02C.is_file()
    assert "reconciled_candidates.faa" in STAGE_02C.read_text()


def test_run_pipeline_registers_02c_with_expected_format() -> None:
    script, description, is_array = _pipeline_step_entry("02c")
    assert script == "02c_genome_reconcile.sh"
    assert description.strip(), "description must be non-empty"
    assert is_array == "no", "02c is not a SLURM array job"


def test_run_pipeline_step_order_places_02c_after_02b_and_before_03() -> None:
    """02c consumes stage-02 outputs; stage 03 consumes its reconciled FASTA."""
    order = _step_order()
    assert "02c" in order, f"02c missing from STEP_ORDER: {order}"
    assert order.index("02") < order.index("02c") < order.index("03")
    assert order.index("02b") < order.index("02c")


def test_dry_run_all_submits_02c_between_02b_and_03() -> None:
    steps = _dry_run_steps("all")
    assert "02c" in steps, f"02c never submitted by `all`: {steps}"
    assert steps.index("02b") < steps.index("02c") < steps.index("03")


def test_dry_run_accepts_02c_as_an_explicit_step() -> None:
    assert _dry_run_steps("02c") == ["02c"]


def test_dry_run_range_02_to_03_includes_02c() -> None:
    steps = _dry_run_steps("02-03")
    assert "02c" in steps, f"range 02-03 skipped 02c: {steps}"


# === F7: 03d must run after 04 ==============================================


def test_run_pipeline_step_order_places_03d_after_04() -> None:
    """03d hard-requires step_completed_04.txt via an exiting check_file."""
    order = _step_order()
    assert order.index("04") < order.index("03d"), (
        f"03d still ordered before 04, so `all` aborts at 03d: {order}"
    )


def test_stage_03d_hard_requires_stage_04_completion() -> None:
    """Pin the requirement this ordering is derived from."""
    text = (REPO_ROOT / "03d_notung_reconciliation.sh").read_text()
    assert "step_completed_04.txt" in text
    assert "step_completed_busco_species_tree.txt" in text


def test_dry_run_all_submits_03d_after_04_and_03a() -> None:
    steps = _dry_run_steps("all")
    assert "03d" in steps
    assert steps.index("04") < steps.index("03d")
    assert steps.index("03a") < steps.index("03d")


# === registry-wide internal consistency =====================================


def test_step_order_and_pipeline_steps_cover_the_same_steps() -> None:
    order, keys = _step_order(), _pipeline_step_keys()
    assert sorted(order) == sorted(keys), (
        f"STEP_ORDER {sorted(order)} != PIPELINE_STEPS {sorted(keys)}"
    )
    assert len(order) == len(set(order)), f"duplicate step in STEP_ORDER: {order}"


def test_every_registered_script_exists() -> None:
    for step in _pipeline_step_keys():
        script, _, _ = _pipeline_step_entry(step)
        assert (REPO_ROOT / script).is_file(), f"{step} -> missing {script}"


def test_dry_run_all_submits_exactly_step_order() -> None:
    """Behavioural: `all` schedules every registered step, in registry order."""
    assert _dry_run_steps("all") == _step_order()


def test_usage_text_enumerates_every_registered_step() -> None:
    """The help text is a step list; a step missing from it is undiscoverable."""
    text = RUN_PIPELINE.read_text()
    m = re.search(r"^STEPS:$(.*?)^OPTIONS:$", text, flags=re.M | re.S)
    assert m, "STEPS block not found in usage text"
    listed = re.findall(r"^\s{2}([0-9]{2}[a-z]?)\s{2,}\S", m.group(1), flags=re.M)
    missing = [s for s in _step_order() if s not in listed]
    assert not missing, f"usage text omits steps: {missing}"


def test_run_pipeline_is_syntactically_valid() -> None:
    proc = subprocess.run(["bash", "-n", str(RUN_PIPELINE)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
