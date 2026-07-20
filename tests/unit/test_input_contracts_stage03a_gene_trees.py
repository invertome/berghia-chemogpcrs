"""03a's ASTRAL input list cannot detect "no gene trees" (audit finding F12.2).

``03a:175`` was::

    find "${RESULTS_DIR}/busco/gene_trees/" -name "*.treefile" \\
        > "${RESULTS_DIR}/busco/gene_trees_list.txt" || { log "Error: ..."; exit 1; }

``find`` exits 0 when it matches nothing -- a zero-match search is a successful
search, not an error -- so the ``|| exit 1`` guard is unreachable for the case
it was written to catch. An EMPTY ``gene_trees_list.txt`` then flowed straight
into ASTRAL on the next line.

The guard only fires when ``find`` itself fails (e.g. the gene_trees/ directory
does not exist), which is a different and rarer condition. Note that 03a's
earlier degenerate-case branch exits before this point only when there are no
trimmed ALIGNMENTS or <3 samples; alignments existing but every IQ-TREE run
failing lands here with an empty list.

These tests execute the REAL find-and-guard block lifted out of the shipping
stage 03a script.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
STAGE03A = PROJECT_ROOT / "03a_busco_species_tree.sh"

ASTRAL_MARKER = 'run_command "astral_species_tree"'


def _extract_gene_tree_list_block() -> str:
    """Lift everything between the ASTRAL section header and the ASTRAL call."""
    lines = STAGE03A.read_text().splitlines()
    start = next(
        (i for i, l in enumerate(lines) if l.startswith("# --- Generate species tree with ASTRAL")),
        None,
    )
    assert start is not None, "03a's ASTRAL section header moved"
    end = next((i for i, l in enumerate(lines) if ASTRAL_MARKER in l), None)
    assert end is not None, "03a no longer invokes ASTRAL"
    assert end > start
    return "\n".join(lines[start + 1 : end])


def _run_block(tmp_path: Path, n_trees: int, *, make_dir: bool = True) -> subprocess.CompletedProcess:
    """Run the real block against a synthetic busco/gene_trees/ tree."""
    results = tmp_path / "results"
    gene_trees = results / "busco" / "gene_trees"
    if make_dir:
        gene_trees.mkdir(parents=True)
        for i in range(n_trees):
            (gene_trees / f"{i}at6447.treefile").write_text("(a:0.1,b:0.2);\n")
    else:
        results.mkdir(parents=True)

    snippet = "\n".join(
        [
            'log() { echo "LOG $*"; }',
            f'RESULTS_DIR="{results}"',
            _extract_gene_tree_list_block(),
            'echo "REACHED_ASTRAL"',
        ]
    )
    return subprocess.run(["bash", "-c", snippet], capture_output=True, text=True)


def _list_lines(tmp_path: Path) -> list[str]:
    p = tmp_path / "results" / "busco" / "gene_trees_list.txt"
    return [l for l in p.read_text().splitlines() if l.strip()] if p.exists() else []


# --------------------------------------------------------------------------
# the regression: zero gene trees must NOT reach ASTRAL
# --------------------------------------------------------------------------
def test_zero_gene_trees_aborts_before_astral(tmp_path):
    """The core bug: `find` exits 0 on no matches, so an empty list flowed on."""
    proc = _run_block(tmp_path, n_trees=0)

    assert "REACHED_ASTRAL" not in proc.stdout, (
        "an empty gene_trees_list.txt must not be handed to ASTRAL:\n" + proc.stdout + proc.stderr
    )
    assert proc.returncode != 0
    assert _list_lines(tmp_path) == []


def test_zero_gene_trees_logs_an_explanatory_error(tmp_path):
    proc = _run_block(tmp_path, n_trees=0)

    combined = proc.stdout + proc.stderr
    assert "LOG" in combined and "rror" in combined, combined


def test_missing_gene_trees_directory_also_aborts(tmp_path):
    """The condition the original guard DID catch must keep working."""
    proc = _run_block(tmp_path, n_trees=0, make_dir=False)

    assert "REACHED_ASTRAL" not in proc.stdout, proc.stdout + proc.stderr
    assert proc.returncode != 0


# --------------------------------------------------------------------------
# the happy paths must be untouched
# --------------------------------------------------------------------------
def test_one_gene_tree_reaches_astral(tmp_path):
    """A single tree is degenerate for ASTRAL but is NOT this guard's job.

    Widening the guard to a minimum tree count would be a behaviour change;
    the fix only closes the zero-match hole.
    """
    proc = _run_block(tmp_path, n_trees=1)

    assert "REACHED_ASTRAL" in proc.stdout, proc.stdout + proc.stderr
    assert proc.returncode == 0
    assert len(_list_lines(tmp_path)) == 1


def test_many_gene_trees_reach_astral_with_a_complete_list(tmp_path):
    proc = _run_block(tmp_path, n_trees=7)

    assert "REACHED_ASTRAL" in proc.stdout, proc.stdout + proc.stderr
    assert proc.returncode == 0
    lines = _list_lines(tmp_path)
    assert len(lines) == 7
    assert all(l.endswith(".treefile") for l in lines)


def test_only_treefiles_are_listed(tmp_path):
    """Non-treefile IQ-TREE by-products must not pad the list into looking non-empty."""
    results = tmp_path / "results"
    gene_trees = results / "busco" / "gene_trees"
    gene_trees.mkdir(parents=True)
    for name in ("g1.log", "g1.iqtree", "g1.ckp.gz"):
        (gene_trees / name).write_text("x\n")

    snippet = "\n".join(
        [
            'log() { echo "LOG $*"; }',
            f'RESULTS_DIR="{results}"',
            _extract_gene_tree_list_block(),
            'echo "REACHED_ASTRAL"',
        ]
    )
    proc = subprocess.run(["bash", "-c", snippet], capture_output=True, text=True)

    assert "REACHED_ASTRAL" not in proc.stdout, (
        "IQ-TREE by-products are not gene trees; the list is still empty:\n" + proc.stdout
    )
    assert proc.returncode != 0


# --------------------------------------------------------------------------
# static guard: the emptiness check must exist and must not be `find`'s status
# --------------------------------------------------------------------------
def test_emptiness_is_checked_separately_from_find_exit_status():
    block = _extract_gene_tree_list_block()
    assert "gene_trees_list.txt" in block
    assert "-s " in block or "-s\"" in block or "wc -l" in block, (
        "03a must test the LIST for content; `find`'s exit status cannot detect zero matches"
    )
