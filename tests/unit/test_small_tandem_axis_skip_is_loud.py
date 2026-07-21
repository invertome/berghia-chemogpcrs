"""Stage 06 must never skip the tandem-cluster axis silently.

The defect
----------
The tandem-cluster detection block was gated on a two-clause conditional with
no ``else``::

    if [ "${RUN_TANDEM_DETECTION:-1}" = "1" ] && [ -f "${GENOME_GFF:-}" ]; then
        ...
    fi                      # <-- nothing said when the gate fails

So an unset/missing ``GENOME_GFF`` produced no message, no marker, and no
non-zero status. Stage 06 then printed "Synteny and mapping analysis
completed." and touched its completion marker, and stage 07 went on to rank
candidates with the tandem axis absent.

Why that is expensive: intra-genome tandem arrays are the field's signature
chemoreceptor signal, and the axis carries ``TANDEM_CLUSTER_WEIGHT=2.5`` --
the largest single weight in the ranking. Losing it is not a rounding error,
and losing it without a trace makes the resulting shortlist unfalsifiable.

The two failure modes are also NOT the same thing and must not print the same
message:

  * ``RUN_TANDEM_DETECTION=0``  -- a deliberate opt-out.
  * ``GENOME_GFF`` missing      -- a required input is unavailable.

The cross-cutting invariant applies to the stale-file case too: if a PREVIOUS
run left ``tandem_clusters.csv`` behind and this run skips detection, stage 07
reads that file and treats last run's measurements as this run's. Detection
being skipped must be announced loudly enough that the stale artifact is not
mistaken for a fresh measurement.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
FUNCTIONS_SH = PROJECT_ROOT / "functions.sh"
STAGE06 = PROJECT_ROOT / "06_synteny_and_mapping.sh"


def _extract_tandem_block() -> str:
    """Lift the live tandem-detection block out of stage 06."""
    lines = STAGE06.read_text().splitlines()
    start = None
    for i, line in enumerate(lines):
        if "intra-genome tandem-cluster detection" in line and line.lstrip().startswith("#"):
            start = i
            break
    if start is None:
        pytest.fail("stage 06 has no tandem-cluster detection block")
    for j in range(start, len(lines)):
        if lines[j].startswith('touch "${RESULTS_DIR}/step_completed_synteny.txt"'):
            return "\n".join(lines[start:j])
    pytest.fail("could not find the end of the tandem block in stage 06")


def _run_tandem_block(tmp_path: Path, *, run_flag=None, gff=True,
                      candidates=True, stale_csv=False) -> subprocess.CompletedProcess:
    """Run the block with a stubbed detector, under the requested conditions."""
    results = tmp_path / "results"
    (results / "chemogpcrs").mkdir(parents=True, exist_ok=True)
    (results / "synteny").mkdir(parents=True, exist_ok=True)
    logs = tmp_path / "logs"
    logs.mkdir(exist_ok=True)
    scripts = tmp_path / "scripts"
    scripts.mkdir(exist_ok=True)

    # Stub the detector so the test never invokes real gffutils work.
    (scripts / "run_tandem_detection.sh").write_text(
        '#!/bin/bash\necho "STUB_DETECTOR_RAN $1"\n'
    )

    if candidates:
        # Real candidate FASTA shape: short IDs as update_headers.py mints them.
        (results / "chemogpcrs" / "chemogpcrs_berghia.fa").write_text(
            ">1287507_1\nMKTLLVLAV\n>1287507_2\nMKTLLVLAC\n"
        )

    tandem_csv = results / "synteny" / "tandem_clusters.csv"
    if stale_csv:
        # Real producer shape: compute_tandem_clusters.py writes these columns.
        tandem_csv.write_text(
            "candidate_id,tandem_cluster_size,tandem_cluster_id\n"
            "1287507_1,4,cluster_0001\n"
        )

    gff_path = tmp_path / "berghia.gff"
    if gff:
        gff_path.write_text(
            "##gff-version 3\n"
            "NC_088360.1\tRefSeq\tgene\t4611\t8065\t.\t+\t.\tID=gene-LOC144963509\n"
        )

    env = os.environ.copy()
    env["LOGS_DIR"] = str(logs)
    env.pop("RUN_TANDEM_DETECTION", None)
    if run_flag is not None:
        env["RUN_TANDEM_DETECTION"] = run_flag

    script = f"""
source "{FUNCTIONS_SH}"; trap - EXIT
RESULTS_DIR="{results}"
SCRIPTS_DIR="{scripts}"
LOGS_DIR="{logs}"
GENOME_GFF="{gff_path if gff else ''}"
TANDEM_CLUSTERS_FILE="{tandem_csv}"
{_extract_tandem_block()}
"""
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=env)


# --------------------------------------------------------------------------
# the happy path still works
# --------------------------------------------------------------------------
def test_detection_runs_when_inputs_are_present(tmp_path):
    proc = _run_tandem_block(tmp_path)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "STUB_DETECTOR_RAN" in proc.stdout + proc.stderr


# --------------------------------------------------------------------------
# the defect: each skip path must announce itself
# --------------------------------------------------------------------------
def test_missing_genome_gff_is_announced(tmp_path):
    """A missing GFF is an UNAVAILABLE input and must be logged as such."""
    proc = _run_tandem_block(tmp_path, gff=False)
    combined = proc.stdout + proc.stderr

    assert "STUB_DETECTOR_RAN" not in combined
    assert "WARN" in combined, (
        "the tandem axis (weight 2.5) was skipped with no warning:\n" + combined
    )
    assert "GENOME_GFF" in combined, (
        "the warning does not name the missing input:\n" + combined
    )


def test_opt_out_is_announced(tmp_path):
    """RUN_TANDEM_DETECTION=0 is deliberate, but still must leave a record."""
    proc = _run_tandem_block(tmp_path, run_flag="0")
    combined = proc.stdout + proc.stderr

    assert "STUB_DETECTOR_RAN" not in combined
    assert "RUN_TANDEM_DETECTION" in combined, (
        "an explicit opt-out left no trace in the log:\n" + combined
    )


def test_opt_out_and_missing_input_are_distinguishable(tmp_path):
    """Deliberate opt-out and unavailable input are different facts.

    Collapsing them would make 'we chose not to measure' indistinguishable from
    'we could not measure', which is exactly the distinction the ranking needs.
    """
    opt_out = _run_tandem_block(tmp_path / "a", run_flag="0")
    no_gff = _run_tandem_block(tmp_path / "b", gff=False)
    assert (opt_out.stdout + opt_out.stderr) != (no_gff.stdout + no_gff.stderr)


def test_skip_warns_about_a_stale_csv_from_a_previous_run(tmp_path):
    """A leftover tandem CSV must not be silently reused as a fresh measurement."""
    proc = _run_tandem_block(tmp_path, gff=False, stale_csv=True)
    combined = proc.stdout + proc.stderr

    assert "stale" in combined.lower(), (
        "detection was skipped while a previous run's tandem_clusters.csv is "
        "still on disk; stage 07 will read it as current data:\n" + combined
    )


def test_no_stale_warning_when_there_is_no_leftover(tmp_path):
    """The stale warning must be about a real file, not fired unconditionally."""
    proc = _run_tandem_block(tmp_path, gff=False, stale_csv=False)
    assert "stale" not in (proc.stdout + proc.stderr).lower()


def test_missing_candidate_fasta_still_announced(tmp_path):
    """The pre-existing inner else must keep working."""
    proc = _run_tandem_block(tmp_path, candidates=False)
    combined = proc.stdout + proc.stderr
    assert "STUB_DETECTOR_RAN" not in combined
    assert "WARN" in combined


# --------------------------------------------------------------------------
# the consumer contract: absent data must read as absent, never as zero
# --------------------------------------------------------------------------
def test_rank_candidates_reports_the_axis_unavailable_rather_than_zero():
    """Stage 07 must gate the tandem vote on presence, not emit a present 0.

    Read-only assertion on rank_candidates.py (owned elsewhere): it already
    gates on ``has_tandem_cluster_data``. Pinning it here means a future change
    that turns 'no data' into a scoring 0.0 -- which would rank an unmeasured
    candidate identically to one measured as having no tandem neighbours --
    fails loudly instead of quietly degrading the shortlist.
    """
    text = (PROJECT_ROOT / "scripts" / "rank_candidates.py").read_text()
    assert "'tandem_cluster': row.get('tandem_cluster_score_norm') if row.get('has_tandem_cluster_data') else None" in text, (
        "the tandem axis is no longer gated on has_tandem_cluster_data; an "
        "unmeasured candidate would vote as a measured zero"
    )
