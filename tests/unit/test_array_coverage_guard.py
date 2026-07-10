"""Tests for assert_array_covers_manifest in functions.sh.

Stages 04 and 05 submit the per-orthogroup work as a SLURM array with a
hardcoded `#SBATCH --array=0-999`. If the orthogroup manifest has more than
1000 entries, every OG at index > the array's top index is silently never
scheduled — no tree, no dN/dS, no error. This guard fails loud when the
submitted array range does not cover the manifest. See bead
berghia-chemogpcrs-fxx.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
FUNCTIONS_SH = PROJECT_ROOT / "functions.sh"


def _run(og_count: str, array_max: str | None, tmp_path: Path) -> subprocess.CompletedProcess:
    env = os.environ.copy()
    env["LOGS_DIR"] = str(tmp_path)
    env.pop("SLURM_ARRAY_TASK_MAX", None)
    arg2 = "" if array_max is None else f' "{array_max}"'
    script = f'source "{FUNCTIONS_SH}"; assert_array_covers_manifest "{og_count}"{arg2}'
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=env)


def test_array_covers_manifest_passes(tmp_path):
    # array 0-1499 (max index 1499) covers a 1500-OG manifest.
    result = _run("1500", "1499", tmp_path)
    assert result.returncode == 0, result.stderr


def test_exact_boundary_passes(tmp_path):
    # array 0-999 (max 999) exactly covers a 1000-OG manifest.
    result = _run("1000", "999", tmp_path)
    assert result.returncode == 0, result.stderr


def test_too_small_array_fails_loud(tmp_path):
    # array 0-999 (max 999) cannot cover a 1500-OG manifest → fail loud.
    result = _run("1500", "999", tmp_path)
    assert result.returncode == 1
    assert "too small" in result.stderr.lower()
    assert "500" in result.stderr  # 1500 - 999 - 1 = 500 orthogroups skipped


def test_not_an_array_job_passes(tmp_path):
    # No array max (local / non-array run) → nothing to check.
    result = _run("1500", None, tmp_path)
    assert result.returncode == 0, result.stderr


def test_unknown_og_count_does_not_block(tmp_path):
    # Count unknown (0) → warn, do not block.
    result = _run("0", "999", tmp_path)
    assert result.returncode == 0, result.stderr
