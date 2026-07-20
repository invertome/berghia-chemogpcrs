"""REPAIR 3 (bead myz9.18) — the OrthoFinder fallback enumeration must be
deterministic and must fail LOUDLY when it matches nothing.

Stages 04 and 05 pick their orthogroup from ``orthogroup_manifest.tsv`` and
fall back to enumerating OrthoFinder's per-OG FASTAs when the manifest is
missing. The original fallback globbed::

    orthogroups/OrthoFinder/Results*/Orthogroups/OG*.fa

which is wrong on BOTH components. The verified layout on the cluster is::

    results/orthogroups/input/OrthoFinder/Results_May15/Orthogroup_Sequences/

-- the root carries an ``input/`` level, and the per-OG FASTAs live in
``Orthogroup_Sequences/``; ``Orthogroups/`` holds TSV summaries. The root
and subdirectory were repaired by routing through
``resolve_orthogroup_sequences_dir``. This module pins the two properties
that repair did NOT establish:

1. **A zero-match enumeration must abort.** Bash leaves an unmatched glob as
   its literal pattern (nullglob is off), so ``ORTHOGROUPS=(dir/OG*.fa)``
   yields a one-element array holding ``dir/OG*.fa`` and ``basename`` turns
   it into the orthogroup name ``OG*``. The stage then runs a full analysis
   against a nonexistent orthogroup and exits 0. A fallback that cannot
   match is worse than no fallback: it converts a clear failure into a
   silent empty result.

2. **Selection among candidates must be deterministically ordered.** The
   SLURM array index selects positionally, and the 50 concurrent tasks run
   on different nodes. If ordering came from filesystem traversal order or
   a locale-dependent collation, index N could denote different orthogroups
   on different nodes -- an identity corruption that every downstream join
   would accept.
"""
from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
HELPER = PROJECT_ROOT / "scripts" / "orthofinder_paths.sh"
STAGE04 = PROJECT_ROOT / "04_phylogenetic_analysis.sh"
STAGE05 = PROJECT_ROOT / "05_selective_pressure_and_asr.sh"
STAGES = [STAGE04, STAGE05]
IDS = [p.name for p in STAGES]

# The layout verified on the cluster (results/orthogroups/...).
RUN_DIR = "input/OrthoFinder/Results_May15"

OG_TSV_TEXT = (
    "Orthogroup\tberghia\tref_species\n"
    "OG0000000\tbste_cand_a\tref_1\n"
)


def _make_run(root: Path, ogs: list[str]) -> Path:
    """Build the real OrthoFinder layout under a synthetic orthogroups root."""
    run = root / RUN_DIR
    (run / "Orthogroups").mkdir(parents=True, exist_ok=True)
    (run / "Orthogroups" / "Orthogroups.tsv").write_text(OG_TSV_TEXT)
    seqs = run / "Orthogroup_Sequences"
    seqs.mkdir(parents=True, exist_ok=True)
    for og in ogs:
        (seqs / f"{og}.fa").write_text(f">{og}_seq1\nMAAA\n")
    return run


def _extract_fallback(script: Path) -> str:
    """Lift the array-mode fallback branch out of a real stage script."""
    text = script.read_text()
    m = re.search(
        r"_og_seq_dir=\$\(resolve_orthogroup_sequences_dir.*?"
        r"base=\$\(basename \"\$og\" \.fa\)",
        text, re.S,
    )
    if not m:
        pytest.fail(
            f"{script.name}: could not locate the array-mode OrthoFinder "
            f"fallback; it must resolve through "
            f"resolve_orthogroup_sequences_dir and set `base`"
        )
    return m.group(0)


def _run_fallback(script: Path, tmp_path: Path, *, task_id: str) -> subprocess.CompletedProcess:
    env = os.environ.copy()
    env["LOGS_DIR"] = str(tmp_path / "logs")
    (tmp_path / "logs").mkdir(exist_ok=True)
    body = f"""
source "{HELPER}"
log() {{ echo "LOG: $*" >&2; }}
RESULTS_DIR="{tmp_path / 'results'}"
SLURM_ARRAY_TASK_ID="{task_id}"
{_extract_fallback(script)}
echo "BASE=${{base}}"
"""
    return subprocess.run(["bash", "-c", body], capture_output=True,
                          text=True, env=env)


@pytest.fixture
def og_root(tmp_path: Path) -> Path:
    return tmp_path / "results" / "orthogroups"


# --------------------------------------------------------------------------
# 1. a zero-match enumeration must be LOUD
# --------------------------------------------------------------------------
@pytest.mark.parametrize("script", STAGES, ids=IDS)
def test_empty_sequences_dir_aborts_instead_of_inventing_an_orthogroup(
    script, tmp_path, og_root
):
    """An Orthogroup_Sequences/ with no OG*.fa must abort, not yield "OG*"."""
    _make_run(og_root, ogs=[])  # directory exists, but contains no OG*.fa
    proc = _run_fallback(script, tmp_path, task_id="0")

    assert proc.returncode != 0, (
        f"{script.name}: a zero-match enumeration exited 0; the stage would "
        f"analyse a nonexistent orthogroup and report success"
    )
    assert "BASE=OG*" not in proc.stdout, (
        f"{script.name}: the unmatched glob pattern leaked through basename "
        f"as the orthogroup name"
    )
    combined = proc.stdout + proc.stderr
    assert re.search(r"no OG\*?\.fa|no orthogroup", combined, re.I), (
        f"{script.name}: the abort must say WHAT was empty. Got:\n{combined}"
    )


@pytest.mark.parametrize("script", STAGES, ids=IDS)
def test_out_of_range_array_index_aborts(script, tmp_path, og_root):
    """An index past the end of the enumeration must abort, not run empty."""
    _make_run(og_root, ogs=["OG0000000", "OG0000001"])
    proc = _run_fallback(script, tmp_path, task_id="99")

    assert proc.returncode != 0, (
        f"{script.name}: array index 99 against 2 orthogroups exited 0; "
        f"an out-of-range index silently yields an empty orthogroup name"
    )
    m = re.search(r"BASE=(\S*)", proc.stdout)
    assert m is None or m.group(1) == "", (
        f"{script.name}: an out-of-range index produced an orthogroup name "
        f"{m.group(1)!r} that does not exist"
    )


# --------------------------------------------------------------------------
# 2. the enumeration must be deterministically ordered
# --------------------------------------------------------------------------
@pytest.mark.parametrize("script", STAGES, ids=IDS)
def test_index_selects_the_same_orthogroup_every_time(script, tmp_path, og_root):
    """Index N must denote one orthogroup, not whatever the FS returns first."""
    ogs = ["OG0000003", "OG0000000", "OG0000002", "OG0000001"]
    _make_run(og_root, ogs=ogs)

    bases = []
    for _ in range(3):
        proc = _run_fallback(script, tmp_path, task_id="1")
        assert proc.returncode == 0, proc.stderr
        m = re.search(r"BASE=(\S+)", proc.stdout)
        assert m, proc.stdout
        bases.append(m.group(1))

    assert len(set(bases)) == 1, f"{script.name}: index 1 was unstable: {bases}"
    assert bases[0] == "OG0000001", (
        f"{script.name}: index 1 gave {bases[0]}, expected the second entry in "
        f"sorted order (OG0000001); the enumeration is not sorted"
    )


def test_both_stages_enumerate_identically(tmp_path, og_root):
    """Stage 04 and stage 05 must agree on which OG an index denotes.

    They are separate array jobs over the same orthogroup set; if their
    enumerations diverged, stage 05 would analyse a different orthogroup
    than stage 04 built the tree for, and every join would still succeed.
    """
    _make_run(og_root, ogs=["OG0000003", "OG0000000", "OG0000002"])
    got = {}
    for script in STAGES:
        proc = _run_fallback(script, tmp_path, task_id="2")
        assert proc.returncode == 0, proc.stderr
        got[script.name] = re.search(r"BASE=(\S+)", proc.stdout).group(1)
    assert len(set(got.values())) == 1, f"stages disagree on index 2: {got}"


# --------------------------------------------------------------------------
# 3. the dead glob must never come back
# --------------------------------------------------------------------------
@pytest.mark.parametrize("script", STAGES, ids=IDS)
def test_the_dead_glob_is_gone(script):
    text = script.read_text()
    for dead in ("orthogroups/OrthoFinder/Results*/Orthogroups/OG",
                 "Results*/Orthogroups/OG*.fa"):
        offending = [
            line for line in text.splitlines()
            if dead in line and not line.lstrip().startswith("#")
        ]
        assert offending == [], (
            f"{script.name}: the dead fallback glob is back (it has both the "
            f"wrong root and the wrong subdirectory, so it can never match): "
            f"{offending}"
        )
