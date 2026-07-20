"""OrthoFinder output discovery in rank_candidates.py (beads -mgci / F3, F4).

Two loaders in ``scripts/rank_candidates.py`` each hardcoded their own guess at
where OrthoFinder writes ``Orthogroups.tsv``, and both guessed a root that the
pipeline never creates:

  F3  ``load_gene_to_orthogroup`` probed ``results/orthofinder/...`` (and a
      ``results/orthology/gene_orthogroups.csv`` that has no producer anywhere
      in the repo).  ``gene_to_og`` was therefore ``{}`` on every run, which in
      turn killed the ``expansion`` axis (weight 1.5), the ``og_confidence``
      axis (weight 1.0), the ``load_og_trees`` call (gated on the dict being
      non-empty) and the BUSTED/MEME corroboration join.

  F4  ``load_og_dnds_reliability`` probed ``results/orthogroups/OrthoFinder``
      and ``results/orthofinder``.  This one is worse than inert: an empty
      weight map means ``dnds_reliability_weight`` is 0.0 for every candidate,
      so ``reliability_shrink`` pulls ``positive_score_norm`` /
      ``purifying_score_norm`` toward the COHORT MEDIAN for the whole cohort —
      the dN/dS axis keeps its weight in the denominator but stops
      discriminating.

The real layout is set by ``03_orthology_clustering.sh``, which runs
``orthofinder -f ${RESULTS_DIR}/orthogroups/input``:

    results/orthogroups/input/OrthoFinder/Results_<stamp>/Orthogroups/Orthogroups.tsv

``Results_<stamp>`` is not knowable ahead of time, which is why the shell never
hardcodes it.  Stages 04, 06c and 07 all use the same recursive idiom:

    find "${RESULTS_DIR}/orthogroups" -name "Orthogroups.tsv" \\
         -path "*/Orthogroups/*" 2>/dev/null | head -1

These tests pin a single Python helper with those semantics, and assert that
BOTH loaders route through it (so the two can no longer drift apart again).

``rank_candidates.py`` is not import-safe (top-level ``sys.argv``), so — as in
tests/unit/test_rank_candidates_dnds_reliability.py and
tests/unit/test_rank_candidates_emits_orthogroup.py — the functions under test
are extracted from the source by AST and executed in a controlled namespace.
"""
from __future__ import annotations

import ast
import fnmatch
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RANK = PROJECT_ROOT / "scripts" / "rank_candidates.py"
SCRIPTS = PROJECT_ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

HELPER = "find_orthogroups_tsv"


# ---------------------------------------------------------------------------
# source extraction (rank_candidates.py is a script, not an importable module)
# ---------------------------------------------------------------------------

def _source() -> str:
    return RANK.read_text()


def _func_node(name: str) -> ast.FunctionDef:
    for node in ast.walk(ast.parse(_source())):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise AssertionError(f"function {name}() not found in scripts/rank_candidates.py")


def _func_src(name: str) -> str:
    seg = ast.get_source_segment(_source(), _func_node(name))
    assert seg, f"could not extract source for {name}()"
    return seg


def _load(*names: str, extra: dict[str, Any] | None = None) -> dict[str, Any]:
    """Exec the named top-level functions into one shared namespace.

    Sharing the namespace is deliberate: it is what lets the loaders call the
    discovery helper by name, which is exactly the coupling under test.
    """
    ns: dict[str, Any] = {
        "os": os,
        "sys": sys,
        "fnmatch": fnmatch,
        "Path": Path,
        "pd": pd,
        "__file__": str(RANK),
    }
    if extra:
        ns.update(extra)
    for name in names:
        exec(_func_src(name), ns)
    return ns


@pytest.fixture(scope="module")
def discover():
    return _load(HELPER)[HELPER]


# ---------------------------------------------------------------------------
# layout builders
# ---------------------------------------------------------------------------

OG_TSV_TEXT = (
    "Orthogroup\tberghia\tref_species\n"
    "OG0000001\tbste_cand_a, bste_cand_b\tref_1, ref_2, ref_3\n"
    "OG0000002\tbste_cand_c\tref_4\n"
)


def _make_results_layout(root: Path, stamp: str = "Results_Jul19") -> Path:
    """The real OrthoFinder layout produced by 03_orthology_clustering.sh:
    ``orthofinder -f ${RESULTS_DIR}/orthogroups/input``.
    """
    og_dir = root / "orthogroups" / "input" / "OrthoFinder" / stamp / "Orthogroups"
    og_dir.mkdir(parents=True, exist_ok=True)
    tsv = og_dir / "Orthogroups.tsv"
    tsv.write_text(OG_TSV_TEXT)
    # Sibling files OrthoFinder really writes — none of them may be picked up.
    (og_dir / "Orthogroups_UnassignedGenes.tsv").write_text("Orthogroup\n")
    (og_dir / "Orthogroups.txt").write_text("OG0000001: bste_cand_a\n")
    wd = root / "orthogroups" / "input" / "OrthoFinder" / stamp / "WorkingDirectory"
    wd.mkdir(parents=True, exist_ok=True)
    (wd / "SpeciesIDs.txt").write_text("0: berghia\n")
    return tsv


def _shell_find(results_dir: Path) -> list[str]:
    """The exact shell idiom stages 04/06c/07 use, for parity checking."""
    proc = subprocess.run(
        ["find", str(results_dir / "orthogroups"),
         "-name", "Orthogroups.tsv", "-path", "*/Orthogroups/*"],
        capture_output=True, text=True,
    )
    return [ln for ln in proc.stdout.splitlines() if ln]


# ---------------------------------------------------------------------------
# 1. discovery finds a realistically nested Results_* layout
# ---------------------------------------------------------------------------

def test_discovery_finds_realistic_results_layout(discover, tmp_path: Path) -> None:
    """The whole point of F3/F4: the path is 5 levels below results/ and
    contains an unpredictable ``Results_<stamp>`` component."""
    results = tmp_path / "results"
    expected = _make_results_layout(results)

    found = discover(str(results))

    assert found is not None, "discovery must locate the real OrthoFinder layout"
    assert os.path.realpath(found) == os.path.realpath(str(expected))


def test_discovery_does_not_pick_sibling_orthofinder_files(discover, tmp_path: Path) -> None:
    """Orthogroups_UnassignedGenes.tsv / Orthogroups.txt sit in the same
    directory and must never be returned."""
    results = tmp_path / "results"
    _make_results_layout(results)

    assert os.path.basename(discover(str(results))) == "Orthogroups.tsv"


@pytest.mark.skipif(shutil.which("find") is None, reason="find(1) not available")
def test_discovery_matches_shell_find_semantics(discover, tmp_path: Path) -> None:
    """Parity with 07_candidate_ranking.sh:373 — Python must select from the
    SAME match set the shell's find produces, so stage 07's OG-coverage columns
    and rank_candidates.py's internal weights describe the same OrthoFinder run.
    """
    results = tmp_path / "results"
    _make_results_layout(results)
    # A decoy under the correct root but NOT inside an Orthogroups/ directory:
    # the shell's -path "*/Orthogroups/*" filter excludes it, so must we.
    decoy = results / "orthogroups" / "input" / "OrthoFinder" / "Results_Jul19"
    (decoy / "Comparative_Genomics_Statistics").mkdir(parents=True, exist_ok=True)
    (decoy / "Comparative_Genomics_Statistics" / "Orthogroups.tsv").write_text("x\n")

    shell_matches = {os.path.realpath(p) for p in _shell_find(results)}
    found = os.path.realpath(discover(str(results)))

    assert shell_matches, "test setup: shell find should match at least one file"
    assert found in shell_matches
    assert os.path.realpath(
        str(decoy / "Comparative_Genomics_Statistics" / "Orthogroups.tsv")
    ) not in shell_matches


# ---------------------------------------------------------------------------
# 2. the multi-match case is deterministic
# ---------------------------------------------------------------------------

def test_discovery_multi_match_is_deterministic(discover, tmp_path: Path) -> None:
    """Several Results_* dirs accumulate when OrthoFinder is re-run (and
    03_orthology_clustering.sh's ``-fg`` resume path creates more).  The shell's
    ``head -1`` takes find's readdir order, which is not reproducible; the
    Python helper must return the same answer every call.
    """
    results = tmp_path / "results"
    for stamp in ("Results_Jul19", "Results_Jul19_1", "Results_Jul19_2"):
        _make_results_layout(results, stamp=stamp)

    answers = {discover(str(results)) for _ in range(10)}

    assert len(answers) == 1, f"discovery is nondeterministic: {answers}"


def test_discovery_multi_match_picks_newest_results_dir(discover, tmp_path: Path) -> None:
    """Tie-break rule, matching ``sort -r | head -1`` in
    03_orthology_clustering.sh:108 — a same-day OrthoFinder re-run appends
    ``_1``, ``_2``, ... so the newest run sorts LAST lexicographically.
    """
    results = tmp_path / "results"
    for stamp in ("Results_Jul19", "Results_Jul19_1", "Results_Jul19_2"):
        _make_results_layout(results, stamp=stamp)

    assert "Results_Jul19_2" in discover(str(results))


# ---------------------------------------------------------------------------
# 3. absent OrthoFinder output degrades gracefully (no crash)
# ---------------------------------------------------------------------------

def test_discovery_returns_none_when_results_dir_missing(discover, tmp_path: Path) -> None:
    assert discover(str(tmp_path / "nope")) is None


def test_discovery_returns_none_when_orthogroups_root_absent(discover, tmp_path: Path) -> None:
    """A results/ that exists but where stage 03 never ran."""
    results = tmp_path / "results"
    (results / "ranking").mkdir(parents=True)

    assert discover(str(results)) is None


def test_discovery_returns_none_when_orthofinder_never_finished(discover, tmp_path: Path) -> None:
    """Stage 03 created the input dir but OrthoFinder produced no Orthogroups."""
    results = tmp_path / "results"
    (results / "orthogroups" / "input").mkdir(parents=True)
    (results / "orthogroups" / "input" / "proteome.fa").write_text(">a\nMA\n")

    assert discover(str(results)) is None


# ---------------------------------------------------------------------------
# 4. F3 — load_gene_to_orthogroup finds the real layout
# ---------------------------------------------------------------------------

def test_load_gene_to_orthogroup_finds_real_orthofinder_output(tmp_path: Path) -> None:
    """The F3 regression itself.  Pre-fix this returned {}, silently zeroing the
    expansion (1.5) and og_confidence (1.0) axes and the BUSTED/MEME join.
    """
    results = tmp_path / "results"
    _make_results_layout(results)

    ns = _load(HELPER, "load_gene_to_orthogroup")
    gene_to_og = ns["load_gene_to_orthogroup"](str(results))

    assert gene_to_og, "gene_to_og must not be empty for a real OrthoFinder run"
    assert gene_to_og["bste_cand_a"] == "OG0000001"
    assert gene_to_og["bste_cand_b"] == "OG0000001"
    assert gene_to_og["bste_cand_c"] == "OG0000002"
    assert gene_to_og["ref_3"] == "OG0000001"


def test_load_gene_to_orthogroup_empty_without_orthofinder(tmp_path: Path) -> None:
    """Degrade gracefully: no OrthoFinder output must yield {} , not a crash."""
    results = tmp_path / "results"
    results.mkdir()

    ns = _load(HELPER, "load_gene_to_orthogroup")

    assert ns["load_gene_to_orthogroup"](str(results)) == {}


def test_load_gene_to_orthogroup_no_longer_probes_orthofinder_root(tmp_path: Path) -> None:
    """``results/orthofinder/`` is a root nothing in the pipeline ever writes.
    Guard against the wrong-root guess creeping back in.
    """
    src = _func_src("load_gene_to_orthogroup")

    assert "'orthofinder'" not in src and '"orthofinder"' not in src


# ---------------------------------------------------------------------------
# 5. F4 — load_og_dnds_reliability finds the real layout
# ---------------------------------------------------------------------------

def _dnds_ns() -> dict[str, Any]:
    return _load(HELPER, "load_og_dnds_reliability",
                 extra={"DNDS_RELIABILITY_FULL": 5.0})


def test_load_og_dnds_reliability_finds_real_orthofinder_output(tmp_path: Path) -> None:
    """The F4 regression.  Pre-fix every OG was absent from the weight map, so
    every candidate got dnds_reliability_weight = 0.0 and the whole cohort's
    selection scores were shrunk to the cohort median.
    """
    results = tmp_path / "results"
    _make_results_layout(results)
    cds_dir = results / "reference_sequences" / "cds"
    cds_dir.mkdir(parents=True)
    # 5 reference CDS, all in OG0000001 -> full credit at DNDS_RELIABILITY_FULL=5.
    (cds_dir / "all_references_cds.fna").write_text(
        "".join(f">ref_{i}\nATGAAA\n" for i in range(1, 6))
    )

    weights = _dnds_ns()["load_og_dnds_reliability"](str(results))

    assert weights, "reliability weights must not be empty for a real run"
    assert weights["OG0000001"] == pytest.approx(3 / 5)   # ref_1..ref_3 in OG1
    assert weights["OG0000002"] == pytest.approx(1 / 5)   # ref_4 in OG2


def test_load_og_dnds_reliability_empty_without_orthofinder(tmp_path: Path) -> None:
    """No OrthoFinder output -> {} (callers already default to 0.0), no crash."""
    results = tmp_path / "results"
    cds_dir = results / "reference_sequences" / "cds"
    cds_dir.mkdir(parents=True)
    (cds_dir / "all_references_cds.fna").write_text(">ref_1\nATGAAA\n")

    assert _dnds_ns()["load_og_dnds_reliability"](str(results)) == {}


# ---------------------------------------------------------------------------
# 6. both loaders share ONE discovery path (anti-drift)
# ---------------------------------------------------------------------------

def _calls(fn_name: str) -> set[str]:
    out: set[str] = set()
    for node in ast.walk(_func_node(fn_name)):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            out.add(node.func.id)
    return out


@pytest.mark.parametrize("loader", ["load_gene_to_orthogroup", "load_og_dnds_reliability"])
def test_both_loaders_call_the_shared_discovery_helper(loader: str) -> None:
    """F3 and F4 were the SAME bug written twice.  Requiring both loaders to
    call one helper is what stops a third copy from drifting in.
    """
    assert HELPER in _calls(loader), (
        f"{loader}() must resolve Orthogroups.tsv via {HELPER}(), not its own probe list"
    )


@pytest.mark.parametrize("loader", ["load_gene_to_orthogroup", "load_og_dnds_reliability"])
def test_loaders_do_not_walk_the_filesystem_themselves(loader: str) -> None:
    """No second, private discovery implementation inside either loader."""
    assert "os.walk" not in _func_src(loader)


def test_discovery_helper_documents_the_shell_idiom() -> None:
    """The helper is the Python twin of the shell's find; keep the pointer so a
    future edit to one side is obviously an edit to a contract.
    """
    doc = ast.get_docstring(_func_node(HELPER)) or ""

    assert "find" in doc and "Orthogroups" in doc
