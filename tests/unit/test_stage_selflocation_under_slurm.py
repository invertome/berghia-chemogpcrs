"""A stage must not locate itself with ``dirname $0`` under Slurm (bead rian).

Slurm COPIES the batch script into a per-job spool directory and executes the
copy, so under ``sbatch`` ``$0`` is ``/var/spool/slurm/slurmd/job<id>/slurm_script``
and ``dirname "$0"`` is that spool directory -- which contains no ``config.sh``.
Observed on Unity job 62037984: FAILED in 1 second, stdout empty, stderr::

    /var/spool/slurm/slurmd/job62037984/slurm_script: line 43: config.sh: No such file or directory

The job's WorkDir was correct; the script cd'd AWAY from it before sourcing
anything. Interactively ``$0`` is the real path, so this works when run by hand
and fails only under sbatch -- which is why it survived.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
STAGES = sorted(REPO.glob("[0-9]*.sh"))

# `SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"` and friends
SELF_LOCATE = re.compile(r'\$\(\s*cd\s+"\$\(dirname\s+"\$(?:0|\{BASH_SOURCE\[0\]\})"')


def _stages_that_self_locate() -> list[Path]:
    return [p for p in STAGES if SELF_LOCATE.search(p.read_text())]


def test_every_self_locating_stage_prefers_the_submit_dir():
    """dirname $0 is the spool copy under sbatch; SLURM_SUBMIT_DIR is the truth."""
    offenders = []
    for path in _stages_that_self_locate():
        if "SLURM_SUBMIT_DIR" not in path.read_text():
            offenders.append(path.name)
    assert not offenders, (
        "these stages resolve their own directory from $0 with no "
        f"SLURM_SUBMIT_DIR fallback, so they cd into the Slurm spool dir: {offenders}")


def test_a_self_locating_stage_fails_loud_when_config_is_absent():
    """A wrong directory must report itself, not emit a bare 'No such file'."""
    for path in _stages_that_self_locate():
        src = path.read_text()
        assert re.search(r'if\s+\[\s+!\s+-f\s+"?config\.sh|\[\s+-f\s+"?config\.sh', src), (
            f"{path.name} sources config.sh without first asserting it is present, "
            "so a wrong working directory produces an unattributed shell error")


def test_stage_06c_is_syntactically_valid():
    subprocess.run(["bash", "-n", str(REPO / "06c_classify_non_chemoreceptors.sh")],
                   check=True)


def test_the_defect_is_the_one_we_think_it_is():
    """Premise guard: 06c is the stage that failed, so it must self-locate."""
    names = {p.name for p in _stages_that_self_locate()}
    assert "06c_classify_non_chemoreceptors.sh" in names, (
        "06c no longer self-locates; if that pattern was removed entirely this "
        "test file and its guards can be revisited")
