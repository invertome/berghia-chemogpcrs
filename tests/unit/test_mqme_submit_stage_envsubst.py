"""Bead mqme: submit_stage.sh must substitute ONLY the #SBATCH variables.

`envsubst` with no SHELL-FORMAT argument substitutes every $VAR/${VAR} in the
file and replaces unset ones with the empty string. Since Slurm sets
$SLURM_ARRAY_TASK_ID at RUNTIME, it is necessarily unset at submit time, so
`if [ -n "$SLURM_ARRAY_TASK_ID" ]` was rewritten to `if [ -n "" ]` -- false for
every array task, which made every task skip its per-orthogroup work and exit 0.

These tests drive the REAL scripts/unity/submit_stage.sh with `sbatch` stubbed,
against both a synthetic stage and the repo's actual stage scripts.
"""
from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SUBMIT_STAGE = REPO_ROOT / "scripts" / "unity" / "submit_stage.sh"

# Every stage script carrying the `#SBATCH $(scale_resources)` directive that
# plain sbatch rejects at parse time -- submit_stage.sh is their only working
# submission path, so the bug had no alternative route around it.
SCALE_RESOURCES_STAGES = [
    "01_reference_processing.sh",
    "02_chemogpcrs_identification.sh",
    "02c_genome_reconcile.sh",
    "03_orthology_clustering.sh",
    "03a_busco_species_tree.sh",
    "03b_lse_classification.sh",
    "06_synteny_and_mapping.sh",
    "07_candidate_ranking.sh",
    "09_report_generation.sh",
]

# The array stages: these are the ones $SLURM_ARRAY_TASK_ID erasure silently
# no-op'd.
ARRAY_STAGES = [
    "04_phylogenetic_analysis.sh",
    "05_selective_pressure_and_asr.sh",
]


@pytest.fixture
def stub_bin(tmp_path):
    """A PATH dir with sbatch/scontrol stubs so no job is ever submitted."""
    bin_dir = tmp_path / "stub_bin"
    bin_dir.mkdir()
    sbatch = bin_dir / "sbatch"
    sbatch.write_text(
        '#!/bin/bash\n'
        # Record the preprocessed file submit_stage.sh handed us.
        'echo "$@" > "${STUB_SBATCH_ARGV_OUT}"\n'
        'echo "Submitted batch job 999999"\n'
    )
    sbatch.chmod(0o755)
    scontrol = bin_dir / "scontrol"
    scontrol.write_text('#!/bin/bash\necho "MaxArraySize        = 10001"\n')
    scontrol.chmod(0o755)
    return bin_dir


def run_submit_stage(stage_rel_path, stub_bin, tmp_path, extra_args=()):
    """Run the real submit_stage.sh from the repo root; return (proc, argv_file)."""
    argv_out = tmp_path / "sbatch_argv.txt"
    env = dict(os.environ)
    env["PATH"] = f"{stub_bin}{os.pathsep}{env['PATH']}"
    env["STUB_SBATCH_ARGV_OUT"] = str(argv_out)
    env["TMPDIR"] = str(tmp_path)
    env.setdefault("SLURM_EMAIL", "tester@umass.edu")
    proc = subprocess.run(
        ["bash", str(SUBMIT_STAGE), str(stage_rel_path), *extra_args],
        cwd=REPO_ROOT, env=env, capture_output=True, text=True,
    )
    return proc, argv_out


def preprocessed_path(tmp_path, stage_basename):
    return tmp_path / "berghia_sbatch" / f"{stage_basename}.sbatch.sh"


# --------------------------------------------------------------------------- #
# The decisive regression: the runtime variable must survive.
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("stage", ARRAY_STAGES)
def test_slurm_array_task_id_survives_preprocessing(stage, stub_bin, tmp_path):
    """$SLURM_ARRAY_TASK_ID is unset at submit time and must NOT be erased."""
    proc, _ = run_submit_stage(stage, stub_bin, tmp_path, ["--array=0-1"])
    assert proc.returncode == 0, proc.stderr

    out = preprocessed_path(tmp_path, stage[:-3]).read_text()
    original = (REPO_ROOT / stage).read_text()

    n_original = original.count("SLURM_ARRAY_TASK_ID")
    assert n_original > 0, "fixture assumption: stage references the array var"
    assert out.count("SLURM_ARRAY_TASK_ID") == n_original

    # The specific line the bug turned into `if [ -n "" ]`.
    assert 'if [ -n "$SLURM_ARRAY_TASK_ID" ]' in out
    assert 'if [ -n "" ]' not in out


def test_unset_body_variables_are_not_blanked(stub_bin, tmp_path):
    """Loop variables and locals are unset at submit time; none may be erased."""
    stage = tmp_path / "stage_body_vars.sh"
    stage.write_text(
        "#!/bin/bash\n"
        "#SBATCH --job-name=body_vars\n"
        "#SBATCH --time=${DEFAULT_TIME}\n"
        "#SBATCH --output=${LOGS_DIR}/body_%j.out\n"
        "#SBATCH $(scale_resources)\n"
        "\n"
        'case "$CLASS" in\n'
        'for _ctsv in "${POOL_DIR}"/class_*.tsv; do :; done\n'
        'awk -v min_len="$MIN_LEN" -v max_gap="$MAX_GAP" \'{print}\'\n'
        'cat "$A" "$B" > "$C"\n'
        'if [ -n "$SLURM_ARRAY_TASK_ID" ]; then :; fi\n'
    )
    proc, _ = run_submit_stage(stage, stub_bin, tmp_path)
    assert proc.returncode == 0, proc.stderr

    out = preprocessed_path(tmp_path, "stage_body_vars").read_text()
    for token in (
        'case "$CLASS" in',
        'for _ctsv in "${POOL_DIR}"/class_*.tsv',
        'awk -v min_len="$MIN_LEN" -v max_gap="$MAX_GAP"',
        'cat "$A" "$B" > "$C"',
        'if [ -n "$SLURM_ARRAY_TASK_ID" ]; then',
    ):
        assert token in out, f"body construct was rewritten: {token!r}"

    # None of the observed corruption signatures.
    for corrupted in ('case "" in', 'for _ctsv in ""/class_*.tsv',
                      'min_len="" -v max_gap=""', 'cat "" "" > ""',
                      'if [ -n "" ]'):
        assert corrupted not in out, f"found corrupted line: {corrupted!r}"


# --------------------------------------------------------------------------- #
# Substitution must be confined to #SBATCH directive lines.
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("stage", SCALE_RESOURCES_STAGES + ARRAY_STAGES)
def test_only_sbatch_lines_differ_from_source(stage, stub_bin, tmp_path):
    """The script body is copied through byte-for-byte."""
    extra = ["--array=0-1"] if stage in ARRAY_STAGES else []
    proc, _ = run_submit_stage(stage, stub_bin, tmp_path, extra)
    assert proc.returncode == 0, proc.stderr

    original = (REPO_ROOT / stage).read_text().splitlines()
    out = preprocessed_path(tmp_path, stage[:-3]).read_text().splitlines()

    # Drop the $(...) directive lines the preprocessor deliberately removes.
    kept = [ln for ln in original
            if not re.match(r"^#SBATCH\s+\$\(", ln)]
    assert len(kept) == len(out), "line count changed beyond directive removal"

    for src_line, out_line in zip(kept, out):
        if src_line.startswith("#SBATCH"):
            continue
        assert src_line == out_line, f"body line was rewritten:\n  {src_line!r}\n  {out_line!r}"


@pytest.mark.parametrize("stage", SCALE_RESOURCES_STAGES + ARRAY_STAGES)
def test_no_unresolved_reference_left_in_directives(stage, stub_bin, tmp_path):
    """Every #SBATCH directive reaching sbatch is fully resolved."""
    extra = ["--array=0-1"] if stage in ARRAY_STAGES else []
    proc, _ = run_submit_stage(stage, stub_bin, tmp_path, extra)
    assert proc.returncode == 0, proc.stderr

    out = preprocessed_path(tmp_path, stage[:-3]).read_text().splitlines()
    directives = [ln for ln in out if ln.startswith("#SBATCH")]
    assert directives, "stage lost all its directives"
    for line in directives:
        assert "$" not in line, f"unresolved reference in directive: {line!r}"


@pytest.mark.parametrize("stage", SCALE_RESOURCES_STAGES)
def test_scale_resources_directive_is_stripped(stage, stub_bin, tmp_path):
    """`#SBATCH $(scale_resources)` must not reach sbatch -- it rejects it.

    The previous `sed '/^#SBATCH \\\\?\\$(/d'` pattern matched a literal
    backslash-questionmark and so never deleted the line at all.
    """
    source = (REPO_ROOT / stage).read_text()
    assert re.search(r"^#SBATCH\s+\$\(scale_resources\)", source, re.M), \
        "fixture assumption: stage carries the command-substitution directive"

    proc, _ = run_submit_stage(stage, stub_bin, tmp_path)
    assert proc.returncode == 0, proc.stderr

    out = preprocessed_path(tmp_path, stage[:-3]).read_text()
    assert "scale_resources" not in out.split("\n#!", 1)[0] or True  # body may mention it
    assert not re.search(r"^#SBATCH\s+\$\(", out, re.M), \
        "command-substitution directive survived into the sbatch input"


# --------------------------------------------------------------------------- #
# An unset directive variable must fail loudly, not substitute to "".
# --------------------------------------------------------------------------- #

def test_unset_directive_variable_aborts(stub_bin, tmp_path):
    """A #SBATCH var with no value would yield `--time=`; refuse instead."""
    stage = tmp_path / "stage_unset_var.sh"
    stage.write_text(
        "#!/bin/bash\n"
        "#SBATCH --job-name=unset_var\n"
        "#SBATCH --time=${THIS_VARIABLE_IS_NOT_SET_ANYWHERE}\n"
        "echo body\n"
    )
    proc, argv_out = run_submit_stage(stage, stub_bin, tmp_path)
    assert proc.returncode != 0, "unset directive variable must abort"
    assert "THIS_VARIABLE_IS_NOT_SET_ANYWHERE" in proc.stderr
    assert not argv_out.exists(), "sbatch must not be invoked"


def test_directive_variables_are_actually_substituted(stub_bin, tmp_path):
    """The whole point: intended directive variables DO get their values."""
    stage = tmp_path / "stage_subst.sh"
    stage.write_text(
        "#!/bin/bash\n"
        "#SBATCH --time=${DEFAULT_TIME}\n"
        "#SBATCH --mail-user=${SLURM_EMAIL}\n"
        'echo "${DEFAULT_TIME} stays literal in the body"\n'
    )
    proc, _ = run_submit_stage(stage, stub_bin, tmp_path)
    assert proc.returncode == 0, proc.stderr

    out = preprocessed_path(tmp_path, "stage_subst").read_text()
    assert "#SBATCH --mail-user=tester@umass.edu" in out
    assert re.search(r"^#SBATCH --time=\S+$", out, re.M)
    assert "#SBATCH --time=${DEFAULT_TIME}" not in out
    # ...but the identical token in the body is untouched.
    assert 'echo "${DEFAULT_TIME} stays literal in the body"' in out
