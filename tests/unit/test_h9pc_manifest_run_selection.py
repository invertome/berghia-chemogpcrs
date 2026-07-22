"""Bead h9pc: create_orthogroup_manifest must pick the OrthoFinder run
DETERMINISTICALLY, and pick the SAME run every other stage does.

The function built the manifest (which stages 04/05 index into for OG NAMES)
by resolving Orthogroups.tsv with an unsorted ``find ... | head -1``. `find`
order is filesystem-dependent and not even stable between invocations, so with
more than one ``Results_*`` run present the manifest could bind OG names to a
DIFFERENT run than the SEQUENCES those stages read (they resolve via
resolve_orthogroup_sequences_dir -> resolve_orthofinder_run). The same
orthogroup id then denotes different gene sets while every join still succeeds.

The fix routes the manifest through the shared authoritative-run rule
(scripts/orthofinder_paths.sh: newest Orthogroups/Orthogroups.tsv by mtime, ties
broken by reverse path order). This test builds two runs whose mtime order and
lexicographic (name) order DISAGREE, and asserts the manifest follows mtime --
i.e. the same run resolve_orthofinder_run returns -- not the unsorted find and
not a name sort.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
OGPATHS = PROJECT_ROOT / "scripts" / "orthofinder_paths.sh"
FUNCTIONS = PROJECT_ROOT / "functions.sh"

HEADER = "Orthogroup\tberghia_1287507\tref_1\n"


def _make_run(root: Path, stamp: str, mtime: float, og_name: str,
              berghia_cells: str, ref_cells: str) -> Path:
    """One OrthoFinder Results_<stamp> tree with a one-row Orthogroups.tsv."""
    run = root / "orthogroups" / "input" / "OrthoFinder" / stamp
    (run / "Orthogroups").mkdir(parents=True, exist_ok=True)
    (run / "Orthogroup_Sequences").mkdir(parents=True, exist_ok=True)
    tsv = run / "Orthogroups" / "Orthogroups.tsv"
    tsv.write_text(HEADER + f"{og_name}\t{berghia_cells}\t{ref_cells}\n")
    os.utime(tsv, (mtime, mtime))
    os.utime(run, (mtime, mtime))
    return run


def _run_manifest(root_results: Path, search_root: str, manifest: str):
    logs = root_results / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    snippet = (
        f'set -o pipefail\n'
        f'export RESULTS_DIR="{root_results}"\n'
        f'export LOGS_DIR="{logs}"\n'
        f'export BERGHIA_TAXID=1287507\n'
        f'source "{OGPATHS}"\n'
        f'source "{FUNCTIONS}"\n'
        f'create_orthogroup_manifest "{search_root}" "{manifest}" >/dev/null\n'
        # Marker-prefixed so an EXIT-trap log line in functions.sh cannot be
        # mistaken for the resolver's output.
        f'echo "RESOLVED_RUN=$(resolve_orthofinder_run "{search_root}")"\n'
    )
    return subprocess.run(["bash", "-c", snippet],
                          capture_output=True, text=True)


def _manifest_orthogroups(path: Path):
    rows = [ln.split("\t") for ln in path.read_text().splitlines()[1:] if ln]
    return {r[1]: r[2] for r in rows}   # {orthogroup: num_seqs}


@pytest.fixture
def results(tmp_path: Path) -> Path:
    d = tmp_path / "results"
    (d / "orthogroups").mkdir(parents=True)
    return d


def test_manifest_follows_newest_mtime_not_name_order(results: Path) -> None:
    """mtime and lexicographic order DISAGREE: 'Results_Sep05' sorts above
    'Results_Dec01' (S > D) but Dec01 has the newer mtime. The manifest must
    reflect the newer-mtime run (OG_NEWRUN), proving it is not a name sort and
    not the old unsorted find|head."""
    root = results
    # older mtime, lexicographically HIGHER stamp -> a name sort would pick this
    _make_run(root, "Results_Sep05", mtime=1000.0, og_name="OG_OLDRUN",
              berghia_cells="bste_old", ref_cells="ref_1")
    # newer mtime, lexicographically LOWER stamp -> the correct (mtime) winner
    _make_run(root, "Results_Dec01", mtime=2000.0, og_name="OG_NEWRUN",
              berghia_cells="bste_a, bste_b", ref_cells="ref_1, ref_2")

    search_root = str(root / "orthogroups")
    manifest = root / "orthogroup_manifest.tsv"
    proc = _run_manifest(root, search_root, str(manifest))
    assert proc.returncode == 0, proc.stderr

    ogs = _manifest_orthogroups(manifest)
    assert "OG_NEWRUN" in ogs, (
        f"manifest bound to the wrong run (a name sort or unsorted find would "
        f"pick OG_OLDRUN); got {ogs}\nstderr:\n{proc.stderr}")
    assert "OG_OLDRUN" not in ogs
    assert ogs["OG_NEWRUN"] == "4"   # bste_a,bste_b,ref_1,ref_2

    # Consistency: the manifest's run is exactly the run every other stage
    # resolves via resolve_orthofinder_run.
    resolved = [ln for ln in proc.stdout.splitlines()
                if ln.startswith("RESOLVED_RUN=")]
    assert resolved, f"no resolver marker in stdout:\n{proc.stdout}"
    chosen_run = resolved[-1].split("=", 1)[1]
    assert chosen_run.endswith("Results_Dec01"), (
        f"manifest run must match resolve_orthofinder_run; got {chosen_run!r}")


def test_manifest_selection_is_stable_across_invocations(results: Path) -> None:
    """Two runs written in the same second (mtime tie) must resolve to the same
    run every time -- the reverse-path tie-break, not filesystem order."""
    root = results
    _make_run(root, "Results_run_a", mtime=1500.0, og_name="OG_RUNA",
              berghia_cells="bste_a", ref_cells="ref_1")
    _make_run(root, "Results_run_b", mtime=1500.0, og_name="OG_RUNB",
              berghia_cells="bste_b", ref_cells="ref_1")
    search_root = str(root / "orthogroups")

    seen = set()
    for _ in range(3):
        manifest = root / "orthogroup_manifest.tsv"
        proc = _run_manifest(root, search_root, str(manifest))
        assert proc.returncode == 0, proc.stderr
        seen.add(frozenset(_manifest_orthogroups(manifest)))
    assert len(seen) == 1, f"selection was not deterministic across runs: {seen}"
