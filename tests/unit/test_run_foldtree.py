"""Tests for scripts/run_foldtree.sh — the stage-08 -> DessimozLab/fold_tree adapter.

Stage 08 calls `${FOLDTREE} --input_dir <dir> --output <tre> --method <m>`, but
the real foldtree CLI is `foldtree --folder <dir w/ structs/> --custom-structs`.
This adapter bridges them: maps --method to a fold_tree metric
(mattypes = foldtree|alntmscore|lddt), stages PDBs into <run>/structs/, and
(on a real run) returns the selected metric tree.

These tests exercise the conda-free paths via --dry-run with .pdb-only input
(CIF->PDB needs gemmi in the foldtree env, validated on Unity).
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent.parent / "scripts" / "run_foldtree.sh"


def _run(args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", str(SCRIPT), *args],
                          capture_output=True, text=True)


def _mk_inputs(d: Path, n_pdb: int = 2) -> Path:
    d.mkdir(parents=True, exist_ok=True)
    for i in range(n_pdb):
        (d / f"prot{i}.pdb").write_text("ATOM  1 N\n")
    return d


class TestMethodMapping:
    @pytest.mark.parametrize("method,metric", [
        ("upgma", "foldtree"),       # stage 08's current value
        ("foldtree", "foldtree"),
        ("tm", "alntmscore"),
        ("alntmscore", "alntmscore"),
        ("lddt", "lddt"),
        ("whatever", "foldtree"),    # unknown -> default foldtree
    ])
    def test_method_maps_to_metric(self, tmp_path: Path, method: str, metric: str) -> None:
        inp = _mk_inputs(tmp_path / "in")
        out = tmp_path / "foldtree.tre"
        r = _run(["--input_dir", str(inp), "--output", str(out),
                  "--method", method, "--dry-run"])
        assert r.returncode == 0, r.stderr
        assert f"metric={metric}" in r.stdout


class TestStructsStaging:
    def test_pdb_files_staged_into_structs(self, tmp_path: Path) -> None:
        inp = _mk_inputs(tmp_path / "in", n_pdb=3)
        out = tmp_path / "foldtree.tre"
        r = _run(["--input_dir", str(inp), "--output", str(out), "--dry-run"])
        assert r.returncode == 0, r.stderr
        assert "structs=3" in r.stdout
        run_dir = tmp_path / "foldtree_run"
        staged = sorted(p.name for p in (run_dir / "structs").glob("*.pdb"))
        assert staged == ["prot0.pdb", "prot1.pdb", "prot2.pdb"]
        # custom-structs requires an (empty) identifiers.txt
        assert (run_dir / "identifiers.txt").exists()

    def test_dry_run_prints_real_cli(self, tmp_path: Path) -> None:
        inp = _mk_inputs(tmp_path / "in")
        out = tmp_path / "foldtree.tre"
        r = _run(["--input_dir", str(inp), "--output", str(out), "--dry-run"])
        assert "foldtree --folder" in r.stdout
        assert "--custom-structs" in r.stdout


class TestArgValidation:
    def test_missing_output_errors(self, tmp_path: Path) -> None:
        inp = _mk_inputs(tmp_path / "in")
        r = _run(["--input_dir", str(inp), "--dry-run"])
        assert r.returncode != 0

    def test_missing_input_dir_errors(self, tmp_path: Path) -> None:
        r = _run(["--input_dir", str(tmp_path / "nope"),
                  "--output", str(tmp_path / "o.tre"), "--dry-run"])
        assert r.returncode != 0
