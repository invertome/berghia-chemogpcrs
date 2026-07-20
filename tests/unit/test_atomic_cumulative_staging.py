"""REPAIR 1 — cumulative-CSV staging must be unique CLUSTER-wide, not per-node.

``05_selective_pressure_and_asr.sh`` runs as ``#SBATCH --array=0-999%50``:
up to 50 array tasks run concurrently on DIFFERENT nodes, all writing into
the same shared ``${RESULTS_DIR}/selective_pressure/``.

``concat_per_og_csv`` rebuilds the cumulative results CSV by globbing every
per-OG CSV and publishing the result with ``mv``. It staged through
``${cumulative}.tmp.$$``. That was treated as the project's canonical fix,
but ``$$`` is a PID, and PIDs are unique only WITHIN a node -- the staging
directory is on shared storage, so two array tasks on two nodes can hold
the same PID and therefore the same staging path. Task A's ``: > "$tmp"``
then truncates the file task B is mid-append to, and B publishes the
truncated result over the cumulative CSV.

``scripts/hpc/run_selection_stack.sh`` fixed the same class of bug with
``TASK_TAG="${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID}.$$"``, which is unique
cluster-wide. This test pins that pattern into stage 05 too, and pins the
general rule that no staging path in either stage script is namespaced by
the PID alone.
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
STAGE05 = PROJECT_ROOT / "05_selective_pressure_and_asr.sh"
SELECTION_STACK = PROJECT_ROOT / "scripts" / "hpc" / "run_selection_stack.sh"

TASK_TAG_PATTERN = '"${SLURM_JOB_ID:-nojob}.${SLURM_ARRAY_TASK_ID:-0}.$$"'


def _extract_shell_function(script: Path, name: str) -> str:
    lines = script.read_text().splitlines()
    start = None
    for i, line in enumerate(lines):
        if re.match(rf"^{re.escape(name)}\(\)\s*\{{\s*$", line):
            start = i
            break
    if start is None:
        pytest.fail(f"{script.name} defines no `{name}()` function")
    for j in range(start + 1, len(lines)):
        if lines[j] == "}":
            return "\n".join(lines[start: j + 1])
    pytest.fail(f"could not find the closing brace of `{name}()` in {script.name}")


def _run_concat(tmp_path: Path, *, job: str, task: str) -> subprocess.CompletedProcess:
    """Run the real concat function with a simulated SLURM array identity.

    ``bash -x`` traces the staging path so the test can assert on the name
    the function actually used, not on the source text.
    """
    env = os.environ.copy()
    env["LOGS_DIR"] = str(tmp_path / "logs")
    (tmp_path / "logs").mkdir(exist_ok=True)
    script = f"""
set -x
SLURM_JOB_ID="{job}"
SLURM_ARRAY_TASK_ID="{task}"
{_extract_shell_function(STAGE05, "concat_per_og_csv")}
concat_per_og_csv "{tmp_path / 'sp'}" "_absrel.csv" "{tmp_path / 'sp' / 'absrel_results.csv'}"
"""
    return subprocess.run(["bash", "-c", script], capture_output=True,
                          text=True, env=env)


def _staging_names(proc: subprocess.CompletedProcess) -> set[str]:
    return set(re.findall(r"absrel_results\.csv\.tmp\.\S+", proc.stderr))


@pytest.fixture
def sp_dir(tmp_path: Path) -> Path:
    d = tmp_path / "sp"
    d.mkdir()
    (d / "OG0000001_absrel.csv").write_text("og,omega\nOG0000001,1.5\n")
    (d / "OG0000002_absrel.csv").write_text("og,omega\nOG0000002,0.2\n")
    return d


# --------------------------------------------------------------------------
# the staging name
# --------------------------------------------------------------------------
def test_concat_staging_is_unique_across_nodes(tmp_path, sp_dir):
    """Two array tasks that share a PID must still stage through distinct paths.

    Same job, different array index -- the exact case ``$$`` alone cannot
    separate when the two tasks land on different nodes.
    """
    a = _run_concat(tmp_path, job="900100", task="7")
    b = _run_concat(tmp_path, job="900100", task="8")
    assert a.returncode == 0 and b.returncode == 0, (a.stderr, b.stderr)

    names_a, names_b = _staging_names(a), _staging_names(b)
    assert names_a and names_b, "no staging path appeared in the bash trace"
    assert not (names_a & names_b), (
        f"array tasks 7 and 8 staged through the SAME path {names_a & names_b}; "
        f"on shared storage that is a truncate-under-append race"
    )


def test_concat_staging_carries_the_slurm_identity(tmp_path, sp_dir):
    proc = _run_concat(tmp_path, job="900100", task="7")
    names = _staging_names(proc)
    assert names, "no staging path appeared in the bash trace"
    assert all("900100" in n and ".7." in n for n in names), (
        f"staging path {names} does not carry the SLURM job/array identity; "
        f"a PID-only tag is unique per node, not per cluster"
    )


def test_concat_still_publishes_the_complete_cumulative(tmp_path, sp_dir):
    """The staging fix must not change what gets published."""
    proc = _run_concat(tmp_path, job="900100", task="7")
    assert proc.returncode == 0, proc.stderr
    published = (sp_dir / "absrel_results.csv").read_text()
    assert published.splitlines()[0] == "og,omega", "header lost"
    assert "OG0000001" in published and "OG0000002" in published
    assert published.count("og,omega") == 1, "header duplicated across per-OG files"


def test_concat_leaves_no_staging_debris(tmp_path, sp_dir):
    _run_concat(tmp_path, job="900100", task="7")
    debris = sorted(p.name for p in sp_dir.iterdir() if ".tmp." in p.name)
    assert debris == [], f"staging debris left in a globbed directory: {debris}"


# --------------------------------------------------------------------------
# the pattern must not regress, and must stay consistent with its sibling
# --------------------------------------------------------------------------
def test_task_tag_matches_the_sibling_script():
    """Stage 05 and run_selection_stack.sh must define TASK_TAG identically."""
    for script in (STAGE05, SELECTION_STACK):
        text = script.read_text()
        assert f"TASK_TAG={TASK_TAG_PATTERN}" in text, (
            f"{script.name} does not define the cluster-unique TASK_TAG; the "
            f"correct fix must be propagated, not applied in one place only"
        )


@pytest.mark.parametrize("script", [STAGE04, STAGE05], ids=lambda p: p.name)
def test_no_pid_only_staging_paths_remain(script):
    """No staging path in either array stage may be namespaced by $$ alone."""
    offenders = [
        f"{script.name}:{i}: {line.strip()}"
        for i, line in enumerate(script.read_text().splitlines(), 1)
        # A trailing `.$$` is only safe when the SLURM array identity is
        # composed into the same expression.
        if re.search(r'\.\$\$(?:\}|")', line)
        and "SLURM_ARRAY_TASK_ID" not in line
    ]
    assert offenders == [], (
        "PID-only staging paths remain; PIDs are unique per node but the "
        "staging directory is shared storage:\n" + "\n".join(offenders)
    )
