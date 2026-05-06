"""Tests for the alignment-filter wrapper scripts.

The four wrappers — run_prequal.sh, run_alignment_ensemble.sh,
run_cloak.sh, run_taper.sh — are thin shell-out layers over upstream
binaries. The tests here exercise the CLI argument validation and
error-handling paths *without* requiring the upstream binaries to be
installed. Integration tests with the real binaries are run on Unity
once `scripts/unity/install_alignment_filters.sh` has completed.

Coverage:
    1. --help exits 0 and prints usage.
    2. Missing required args fails with exit 1 and an error message.
    3. Missing input file fails with exit 1.
    4. Missing upstream binary fails with exit 2 (not 1, so the caller
       can distinguish "user error" from "tool not installed").
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
SCRIPTS_DIR = PROJECT_ROOT / "scripts"


WRAPPERS = {
    "prequal": {
        "script": SCRIPTS_DIR / "run_prequal.sh",
        "required": ["--input=", "--output="],
        "missing_tool_env": {"PREQUAL": "/nonexistent/prequal_xyz_xyz"},
    },
    "alignment_ensemble": {
        "script": SCRIPTS_DIR / "run_alignment_ensemble.sh",
        "required": ["--input=", "--output-dir="],
        # MAFFT is a required precondition — point env at a missing path
        "missing_tool_env": {"MAFFT": "/nonexistent/mafft_xyz_xyz"},
    },
    "cloak": {
        "script": SCRIPTS_DIR / "run_cloak.sh",
        "required": ["--ensemble-dir=", "--output="],
        "missing_tool_env": {"CLOAK": "/nonexistent/cloak.py"},
    },
    "taper": {
        "script": SCRIPTS_DIR / "run_taper.sh",
        "required": ["--input=", "--output="],
        "missing_tool_env": {"TAPER": "/nonexistent/taper.jl"},
    },
}


def _run(script: Path, args: list[str], env_extra: dict | None = None,
         input_file: Path | None = None) -> subprocess.CompletedProcess:
    env = os.environ.copy()
    if env_extra:
        env.update(env_extra)
    cmd = ["bash", str(script), *args]
    return subprocess.run(cmd, capture_output=True, text=True, env=env)


# ---- 1. --help exits 0 and prints usage --------------------------------

@pytest.mark.parametrize("name", list(WRAPPERS.keys()))
def test_help_exits_zero(name: str) -> None:
    """`bash run_X.sh --help` must exit 0 and print a usage block to stderr."""
    result = _run(WRAPPERS[name]["script"], ["--help"])
    assert result.returncode == 0, (
        f"{name} --help exited {result.returncode}\nstderr={result.stderr}")
    # Usage should mention the script name or 'Required'/'Usage'
    combined = result.stdout + result.stderr
    assert "Required" in combined or "Usage" in combined, (
        f"{name} --help output missing usage block:\n{combined}")


# ---- 2. Missing required args fails ------------------------------------

@pytest.mark.parametrize("name", list(WRAPPERS.keys()))
def test_no_args_fails(name: str) -> None:
    """Running with no args must exit non-zero with an ERROR message."""
    result = _run(WRAPPERS[name]["script"], [])
    assert result.returncode != 0, (
        f"{name} with no args unexpectedly succeeded")
    assert "ERROR" in result.stderr or "required" in result.stderr.lower(), (
        f"{name} stderr lacks ERROR message:\n{result.stderr}")


@pytest.mark.parametrize("name", list(WRAPPERS.keys()))
def test_partial_args_fails(name: str) -> None:
    """Providing only the first required arg must still fail."""
    first = WRAPPERS[name]["required"][0]
    # Use a path that exists so the failure is from the *missing other arg*
    result = _run(WRAPPERS[name]["script"], [f"{first}/dev/null"])
    assert result.returncode != 0


# ---- 3. Missing input file fails ---------------------------------------

@pytest.mark.parametrize("name", ["prequal", "alignment_ensemble", "taper"])
def test_missing_input_file_fails(name: str, tmp_path: Path) -> None:
    """If --input points at a nonexistent file, fail with exit 1."""
    out = tmp_path / "out.fa"
    nonexistent = tmp_path / "does_not_exist.fa"
    args = [f"--input={nonexistent}"]
    if name == "alignment_ensemble":
        args.append(f"--output-dir={tmp_path / 'ens'}")
    else:
        args.append(f"--output={out}")
    result = _run(WRAPPERS[name]["script"], args)
    assert result.returncode == 1, (
        f"{name} with missing input expected exit 1, got {result.returncode}\n"
        f"stderr={result.stderr}")
    assert "not found" in result.stderr.lower()


def test_cloak_missing_ensemble_dir_fails(tmp_path: Path) -> None:
    """CLOAK with a nonexistent --ensemble-dir must fail with exit 1."""
    out = tmp_path / "out.fa"
    nonexistent = tmp_path / "no_ensemble"
    result = _run(WRAPPERS["cloak"]["script"],
                  [f"--ensemble-dir={nonexistent}", f"--output={out}"])
    assert result.returncode == 1
    assert "not found" in result.stderr.lower()


# ---- 4. Missing upstream binary fails with exit 2 ----------------------

def _make_min_fasta(path: Path) -> None:
    path.write_text(">a\nMKLLR\n>b\nMKLLR\n>c\nMKLLR\n")


def test_prequal_missing_binary_fails_with_exit_2(tmp_path: Path) -> None:
    """If $PREQUAL points at a missing binary, exit 2 (not 1) so the caller
    can distinguish 'tool not installed' from user-input errors."""
    inp = tmp_path / "in.fa"
    _make_min_fasta(inp)
    out = tmp_path / "out.fa"
    result = _run(WRAPPERS["prequal"]["script"],
                  [f"--input={inp}", f"--output={out}",
                   "--prequal-bin=/nonexistent/prequal_xyz_xyz"])
    assert result.returncode == 2, (
        f"prequal missing-bin: expected 2, got {result.returncode}\n{result.stderr}")
    assert "not found" in result.stderr.lower() or "install" in result.stderr.lower()


def test_taper_missing_julia_fails_with_exit_2(tmp_path: Path) -> None:
    """TAPER without julia binary must exit 2."""
    inp = tmp_path / "in.fa"
    _make_min_fasta(inp)
    out = tmp_path / "out.fa"
    fake_jl = tmp_path / "fake.jl"
    fake_jl.write_text("# placeholder\n")
    result = _run(WRAPPERS["taper"]["script"],
                  [f"--input={inp}", f"--output={out}",
                   f"--taper-jl={fake_jl}",
                   "--julia-bin=/nonexistent/julia_xyz_xyz"])
    assert result.returncode == 2


def test_taper_missing_jl_fails_with_exit_2(tmp_path: Path) -> None:
    """TAPER without correction_multi.jl must exit 2."""
    inp = tmp_path / "in.fa"
    _make_min_fasta(inp)
    out = tmp_path / "out.fa"
    result = _run(WRAPPERS["taper"]["script"],
                  [f"--input={inp}", f"--output={out}",
                   "--taper-jl=/nonexistent/correction_multi.jl"])
    assert result.returncode == 2


def test_cloak_missing_py_fails_with_exit_2(tmp_path: Path) -> None:
    """CLOAK without cloak.py must exit 2."""
    ens = tmp_path / "ens"
    ens.mkdir()
    (ens / "a.fa").write_text(">x\nMK\n")
    (ens / "b.fa").write_text(">x\nMK\n")
    out = tmp_path / "out.fa"
    result = _run(WRAPPERS["cloak"]["script"],
                  [f"--ensemble-dir={ens}", f"--output={out}",
                   "--cloak-py=/nonexistent/cloak.py"])
    assert result.returncode == 2


# ---- 5. CLOAK requires >= 2 ensemble members ---------------------------

def test_cloak_fails_with_single_member_ensemble(tmp_path: Path) -> None:
    """CLOAK with only 1 alignment in the ensemble dir is an error: the
    consensus algorithm requires at least 2 alignments to compare."""
    ens = tmp_path / "ens"
    ens.mkdir()
    (ens / "only.fa").write_text(">x\nMK\n")  # just one
    out = tmp_path / "out.fa"
    # Provide a stub cloak.py path so we get past the binary-check and
    # hit the >=2-member check.
    fake_cloak = tmp_path / "cloak.py"
    fake_cloak.write_text("import sys\n")
    result = _run(WRAPPERS["cloak"]["script"],
                  [f"--ensemble-dir={ens}", f"--output={out}",
                   f"--cloak-py={fake_cloak}"])
    assert result.returncode == 3, (
        f"expected exit 3 for too-few ensemble members, got {result.returncode}\n"
        f"{result.stderr}")
    assert "2" in result.stderr  # the threshold message
