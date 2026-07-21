"""run_pipeline.sh must submit through a path sbatch will actually accept (bead mban).

Eleven numbered stages carry ``#SBATCH --time=${DEFAULT_TIME}`` /
``#SBATCH $(scale_resources)``. sbatch parses directives BEFORE the body runs,
so it sees those tokens verbatim and rejects the job -- verified on Unity:

    sbatch --test-only 05_selective_pressure_and_asr.sh
    -> sbatch: error: Invalid --time specification

A command-line ``--time`` override does NOT rescue it, because the in-file
directive is still parsed. ``run_pipeline.sh`` submitted with plain ``sbatch``
and, when no job id came back, printed a warning and CONTINUED -- so the
orchestrator would walk the whole pipeline, submit nothing for 11 of 17 stages,
and exit reporting warnings. An orchestrator that completes while submitting
nothing is the same silent-success shape as the stage-05 no-op, one level up.

``scripts/unity/submit_stage.sh`` already exists for exactly this and was
hardened under bead mqme: substitution touches only ``#SBATCH`` lines, each
variable must resolve non-empty or it aborts, and it refuses to hand sbatch any
directive still containing ``$``.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
RUN_PIPELINE = REPO / "run_pipeline.sh"
SUBMIT_STAGE = REPO / "scripts/unity/submit_stage.sh"
SRC = RUN_PIPELINE.read_text()

# The stages whose #SBATCH headers cannot survive plain sbatch.
UNEXPANDED = re.compile(r'^#SBATCH\s+(--\w[\w-]*=)?[^\n]*(\$\{|\$\()', re.M)


def _stages_with_unexpanded_directives() -> list[str]:
    out = []
    for path in sorted(REPO.glob("[0-9]*.sh")):
        header = [ln for ln in path.read_text().splitlines()
                  if ln.startswith("#SBATCH")]
        header = [ln for ln in header
                  if not re.search(r'--(output|error|mail-user)=', ln)]
        if any(("${" in ln) or ("$(" in ln) for ln in header):
            out.append(path.name)
    return out


def test_the_defect_is_real_and_widespread():
    """Guard the premise: if this ever returns empty the fix can be reverted."""
    affected = _stages_with_unexpanded_directives()
    assert len(affected) >= 5, (
        f"expected many stages with unexpanded #SBATCH directives, got {affected}")
    assert "05_selective_pressure_and_asr.sh" in affected


def test_run_pipeline_is_syntactically_valid():
    subprocess.run(["bash", "-n", str(RUN_PIPELINE)], check=True)


def test_the_hardened_submitter_exists():
    assert SUBMIT_STAGE.is_file(), "the documented submission route is missing"


def test_run_pipeline_submits_through_the_hardened_submitter():
    """Plain `sbatch "$script_path"` cannot submit 11 of the 17 stages."""
    assert "submit_stage.sh" in SRC, (
        "run_pipeline.sh does not route submission through submit_stage.sh, so "
        "every stage with an unexpanded #SBATCH directive is rejected at parse time")


def test_no_bare_sbatch_submission_of_a_stage_script_remains():
    """A surviving bare `sbatch ... "$script_path"` would reintroduce the defect."""
    bare = [ln.strip() for ln in SRC.splitlines()
            if re.search(r'(?<!_)\bsbatch\$\{?slurm_opts\}?\s+"\$script_path"', ln)]
    assert not bare, f"bare sbatch submission still present: {bare}"


def test_a_failed_submission_is_fatal_not_a_warning():
    """The load-bearing gate.

    Detection that does not gate has done nothing: previously the empty-job-id
    branch warned and returned "", and the orchestrator moved to the next stage.
    """
    assert not re.search(r'Warning:\s*Could not parse job ID', SRC), (
        "a submission that produced no job id is still only a warning")
    assert re.search(r'(Could not parse job ID|submission failed)[^\n]*', SRC), (
        "the failed-submission branch disappeared entirely; it must still be handled")
    # the branch must terminate, not fall through
    m = re.search(r'(Could not parse job ID|submission failed).{0,400}', SRC, re.S)
    assert m, "the failed-submission branch could not be located"
    assert re.search(r'\b(exit|return)\s+[1-9]', m.group(0)), (
        "the failed-submission branch does not exit non-zero")


def test_the_call_site_aborts_when_a_stage_cannot_be_submitted():
    """submit_step exits inside a command substitution, where the exit status
    is easy to lose. The loop must check it explicitly."""
    m = re.search(r'for step in "\$\{steps\[@\]\}".*?\n    done', SRC, re.S)
    assert m, "the submission loop could not be located"
    body = m.group(0)
    assert "submit_step" in body
    assert re.search(r'if\s+!\s+prev_job_id=\$\(submit_step|\|\|\s*\{', body), (
        "the loop assigns submit_step's output without checking whether the "
        "submission actually succeeded")
    assert re.search(r'exit\s+[1-9]', body), "the loop does not abort on failure"
