"""Bead raky: every numbered stage must self-activate the project conda env.

0 of the numbered stages activated conda, and `set -u` is armed nowhere, so a
clean-shell `sbatch` submit ran on the base PATH and the stage's tools (TMbed,
HMMER, IQ-TREE, ...) were simply missing -- this is what killed stage 02c. The
fix adds a shared helper `activate_conda_env` (functions.sh) that each stage
calls right after sourcing config/functions, using the single-source-of-truth
`CONDA_ENV` from config.sh (default "berghia-gpcr").

These tests pin the helper's contract:
  1. it invokes `conda activate <CONDA_ENV>` with the value config.sh exports;
  2. it does not abort a caller running under `set -e`;
  3. with no conda present (a unit-test host) it warns and returns 0 rather than
     tripping `set -e`;
and a source scan that every in-scope stage actually calls it.

The helper is exercised through a stub `conda` on PATH that records its args, so
no real conda/env is needed.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
CONFIG = PROJECT_ROOT / "config.sh"
FUNCTIONS = PROJECT_ROOT / "functions.sh"

# Scope of bead raky: the 9 primary stages (0N_*.sh, underscore right after the
# digits) plus the two `set -eo pipefail` specials 02c and 06c. The lettered
# sub-stages (02a/02b/03a-d/04b) carry the identical defect (independent Slurm
# jobs that source config+functions and activate nothing) and were brought into
# the same fix, so they are covered here too.
PRIMARY_STAGES = sorted(PROJECT_ROOT.glob("0[1-9]_*.sh"))
SPECIAL_STAGES = [
    PROJECT_ROOT / "02c_genome_reconcile.sh",
    PROJECT_ROOT / "06c_classify_non_chemoreceptors.sh",
]
LETTERED_STAGES = [
    PROJECT_ROOT / "02a_cluster_sequences.sh",
    PROJECT_ROOT / "02b_classify_gpcrs.sh",
    PROJECT_ROOT / "03a_busco_species_tree.sh",
    PROJECT_ROOT / "03b_lse_classification.sh",
    PROJECT_ROOT / "03c_cafe_analysis.sh",
    PROJECT_ROOT / "03d_notung_reconciliation.sh",
    PROJECT_ROOT / "04b_ecl_analysis.sh",
]
IN_SCOPE_STAGES = PRIMARY_STAGES + SPECIAL_STAGES + LETTERED_STAGES


def _write_recording_conda_stub(bindir: Path, record: Path) -> None:
    """A `conda` on PATH that appends its args to `record` and exits 0."""
    bindir.mkdir(parents=True, exist_ok=True)
    stub = bindir / "conda"
    stub.write_text(
        "#!/usr/bin/env bash\n"
        f'printf "%s\\n" "$*" >> "{record}"\n'
        "exit 0\n"
    )
    stub.chmod(0o755)


def _prelude(tmp_path: Path) -> str:
    """Source config + functions with side effects redirected to tmp_path."""
    logs = tmp_path / "results" / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    return (
        f'source "{CONFIG}"\n'
        # Redirect the EXIT-trap / cleanup side effects away from the real repo.
        f'export RESULTS_DIR="{tmp_path / "results"}"\n'
        f'export LOGS_DIR="{logs}"\n'
        f'source "{FUNCTIONS}"\n'
        # Ensure the "already active" short-circuit does not fire in CI/dev shells
        # that happen to run inside the target env.
        'unset CONDA_DEFAULT_ENV\n'
    )


def test_activate_invokes_conda_activate_with_config_env(tmp_path: Path) -> None:
    """activate_conda_env must call `conda activate <CONDA_ENV>` where CONDA_ENV
    is the value config.sh exports (default berghia-gpcr)."""
    record = tmp_path / "conda_args.txt"
    bindir = tmp_path / "bin"
    _write_recording_conda_stub(bindir, record)

    snippet = (
        _prelude(tmp_path)
        + f'export PATH="{bindir}:$PATH"\n'
        + 'echo "CONDA_ENV=${CONDA_ENV}"\n'
        + 'activate_conda_env\n'
    )
    proc = subprocess.run(["bash", "-c", snippet], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr

    conda_env = [ln.split("=", 1)[1] for ln in proc.stdout.splitlines()
                 if ln.startswith("CONDA_ENV=")]
    assert conda_env, f"config.sh did not export CONDA_ENV:\n{proc.stdout}"
    env = conda_env[-1]
    assert env == "berghia-gpcr", f"unexpected CONDA_ENV default: {env!r}"

    assert record.exists(), f"stub conda was never called\nstderr:\n{proc.stderr}"
    calls = [ln for ln in record.read_text().splitlines() if ln]
    assert calls == [f"activate {env}"], (
        f"helper must invoke `conda activate {env}`; got {calls}")


def test_activate_does_not_error_under_set_e(tmp_path: Path) -> None:
    """Called under `set -e` with a working (stub) conda, the helper returns 0
    and execution continues past it."""
    record = tmp_path / "conda_args.txt"
    bindir = tmp_path / "bin"
    _write_recording_conda_stub(bindir, record)

    snippet = (
        _prelude(tmp_path)
        + f'export PATH="{bindir}:$PATH"\n'
        + 'set -e\n'
        + 'activate_conda_env\n'
        + 'echo "REACHED_END=1"\n'
    )
    proc = subprocess.run(["bash", "-c", snippet], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert "REACHED_END=1" in proc.stdout, (
        f"set -e aborted after activate_conda_env\nstdout:\n{proc.stdout}"
        f"\nstderr:\n{proc.stderr}")


def test_activate_no_conda_returns_zero_and_warns(tmp_path: Path) -> None:
    """With no conda on PATH and no profile to source, the helper must WARN and
    return 0 -- never abort a caller's `set -e` on a host without conda."""
    missing_profile = tmp_path / "nonexistent-conda.sh"
    snippet = (
        _prelude(tmp_path)
        + f'export CONDA_PROFILE="{missing_profile}"\n'
        # A minimal PATH with coreutils but no conda binary.
        + 'export PATH="/usr/bin:/bin"\n'
        + 'set -e\n'
        + 'activate_conda_env\n'
        + 'echo "REACHED_END=1"\n'
    )
    proc = subprocess.run(["bash", "-c", snippet], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert "REACHED_END=1" in proc.stdout, (
        f"helper aborted set -e when conda was absent\nstderr:\n{proc.stderr}")
    assert "conda not found" in proc.stderr, (
        f"expected a WARN about missing conda; stderr:\n{proc.stderr}")


def test_all_in_scope_stages_call_activate_conda_env() -> None:
    """Every in-scope numbered stage (9 primaries + 02c + 06c) invokes the
    shared helper. Guards against a stage being added or reverted without it."""
    assert len(PRIMARY_STAGES) == 9, (
        f"expected 9 primary stages 0N_*.sh, found {len(PRIMARY_STAGES)}: "
        f"{[p.name for p in PRIMARY_STAGES]}")

    missing = []
    for stage in IN_SCOPE_STAGES:
        assert stage.exists(), f"stage not found: {stage}"
        calls_it = any(
            ln.lstrip().split("#", 1)[0].strip() == "activate_conda_env"
            for ln in stage.read_text().splitlines()
        )
        if not calls_it:
            missing.append(stage.name)
    assert not missing, f"stages missing an activate_conda_env call: {missing}"
