"""Stage 03d must find per-class per-OG alignments, and both backends must
agree about what a gene tree is (bead F6).

Two defects, same root cause as the stage 05 no-op:

1. The GeneRax ``families.txt`` builder looked for each family's alignment at
   the pre-refactor FLAT paths ``phylogenies/protein/<base>_trimmed.fa`` /
   ``_aligned.fa``. Nothing writes those, so ``[ -z "$aln" ] && continue``
   dropped every family and ``families.txt`` ended up as a bare
   ``[FAMILIES]`` header -- GeneRax reconciled nothing.

2. The bash ``find`` that collects gene trees is RECURSIVE, but the Python
   fallback globbed ``Path(phylo_dir).glob("*.treefile")`` -- non-recursive, so
   it saw zero trees. The two reconciliation backends disagreed about what a
   gene tree even is.

The alignment is now located relative to the tree's OWN directory, which is
correct for all three layouts stage 04 produces (per-OG ``class_<C>/<base>/``,
per-class global ``class_<C>/``, and the flat ``all_berghia_refs.*``
back-compat aliases) without stage 03d needing to re-derive the class.
"""
from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
FUNCTIONS_SH = PROJECT_ROOT / "functions.sh"
STAGE03D = PROJECT_ROOT / "03d_notung_reconciliation.sh"


# --------------------------------------------------------------------------
# harness
# --------------------------------------------------------------------------
def _extract_families_block() -> str:
    """Lift the live families.txt builder out of stage 03d."""
    lines = STAGE03D.read_text().splitlines()
    start = None
    for i, line in enumerate(lines):
        if line.strip() == "{" and i + 1 < len(lines) and '[FAMILIES]' in lines[i + 1]:
            start = i
            break
    if start is None:
        pytest.fail("stage 03d has no `{ echo \"[FAMILIES]\" ...` block")
    for j in range(start + 1, len(lines)):
        if lines[j].strip().startswith("} >"):
            return "\n".join(lines[start : j + 1])
    pytest.fail("could not find the end of the families.txt block in stage 03d")


def _extract_find_command() -> str:
    """Lift the live gene-tree `find` out of stage 03d (following continuations)."""
    lines = STAGE03D.read_text().splitlines()
    for i, line in enumerate(lines):
        if line.startswith("GENE_TREES=$(find"):
            block = [line]
            while block[-1].rstrip().endswith("\\") and i + len(block) < len(lines):
                block.append(lines[i + len(block)])
            return "\n".join(block)
    pytest.fail("stage 03d has no `GENE_TREES=$(find ...)` line")


def _make_perog_tree(results: Path, cls: str, base: str, *, align=True) -> Path:
    d = results / "phylogenies" / "protein" / f"class_{cls}" / base
    d.mkdir(parents=True, exist_ok=True)
    tree = d / f"{base}.treefile"
    tree.write_text("(s1:0.1,s2:0.2);\n")
    if align:
        (d / f"{base}_trimmed.fa").write_text(">s1\nMKFLV\n>s2\nMKFLA\n")
    return tree


def _run_families_builder(tmp_path: Path) -> tuple[subprocess.CompletedProcess, Path]:
    results = tmp_path / "results"
    env = os.environ.copy()
    env["LOGS_DIR"] = str(tmp_path / "logs")
    (tmp_path / "logs").mkdir(exist_ok=True)
    out = tmp_path / "families.txt"
    generax_out = tmp_path / "generax"
    generax_out.mkdir(exist_ok=True)
    script = f"""
source "{FUNCTIONS_SH}"; trap - EXIT
RESULTS_DIR="{results}"
PHYLO_DIR="{results}/phylogenies"
GENERAX_OUT="{generax_out}"
FAMILIES_TXT="{out}"
{_extract_find_command()}
{_extract_families_block()}
"""
    proc = subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=env)
    return proc, out


# --------------------------------------------------------------------------
# defect 1: families.txt was empty
# --------------------------------------------------------------------------
def test_families_txt_is_non_empty_for_perclass_perog_trees(tmp_path):
    """The layout stage 04 actually writes must yield real GeneRax families."""
    results = tmp_path / "results"
    _make_perog_tree(results, "A", "OG0000123")
    _make_perog_tree(results, "A", "OG0000124")

    proc, out = _run_families_builder(tmp_path)
    assert proc.returncode == 0, proc.stdout + proc.stderr

    text = out.read_text()
    assert text.splitlines()[0] == "[FAMILIES]"
    assert "- OG0000123" in text, f"families.txt has no entries:\n{text}"
    assert "- OG0000124" in text
    assert text.count("starting_gene_tree = ") == 2
    assert text.count("alignment = ") == 2


def test_families_entry_points_at_the_perog_alignment(tmp_path):
    """The alignment path recorded must be the per-OG one that exists on disk."""
    results = tmp_path / "results"
    tree = _make_perog_tree(results, "C", "OG0000900")

    proc, out = _run_families_builder(tmp_path)
    assert proc.returncode == 0, proc.stdout + proc.stderr

    expected = tree.parent / "OG0000900_trimmed.fa"
    assert f"alignment = {expected}" in out.read_text()


def test_family_without_an_alignment_is_still_skipped(tmp_path):
    """A tree with no alignment beside it remains excluded (unchanged contract)."""
    results = tmp_path / "results"
    _make_perog_tree(results, "A", "OG0000123")
    _make_perog_tree(results, "A", "OG0000555", align=False)

    proc, out = _run_families_builder(tmp_path)
    assert proc.returncode == 0, proc.stdout + proc.stderr

    text = out.read_text()
    assert "- OG0000123" in text
    assert "OG0000555" not in text


def test_aligned_fa_is_accepted_when_trimmed_is_absent(tmp_path):
    """_aligned.fa remains the documented second choice after _trimmed.fa."""
    results = tmp_path / "results"
    d = results / "phylogenies" / "protein" / "class_A" / "OG0000123"
    d.mkdir(parents=True)
    (d / "OG0000123.treefile").write_text("(s1:0.1,s2:0.2);\n")
    (d / "OG0000123_aligned.fa").write_text(">s1\nMKFLV\n>s2\nMKFLA\n")

    proc, out = _run_families_builder(tmp_path)
    assert proc.returncode == 0, proc.stdout + proc.stderr

    assert f"alignment = {d / 'OG0000123_aligned.fa'}" in out.read_text()


def test_no_flat_perog_alignment_paths_remain(tmp_path):
    """No writer exists for protein/<base>_trimmed.fa -- no reader may remain."""
    src = STAGE03D.read_text()
    for dead in (
        '"${RESULTS_DIR}/phylogenies/protein/${base}_trimmed.fa"',
        '"${RESULTS_DIR}/phylogenies/protein/${base}_aligned.fa"',
    ):
        assert dead not in src, f"stage 03d still reads the unwritten flat path {dead}"


# --------------------------------------------------------------------------
# defect 2: the two backends disagreed about what a gene tree is
# --------------------------------------------------------------------------
def _extract_python_tree_discovery() -> str:
    """Lift the Python fallback's gene-tree discovery line out of stage 03d."""
    for line in STAGE03D.read_text().splitlines():
        if line.startswith("gene_tree_files = "):
            return line
    pytest.fail("stage 03d's Python fallback has no `gene_tree_files = ...` line")


def test_python_fallback_sees_the_same_trees_as_bash_find(tmp_path):
    """The Python glob and the bash find must return the same set of trees."""
    results = tmp_path / "results"
    phylo = results / "phylogenies"
    _make_perog_tree(results, "A", "OG0000123")
    _make_perog_tree(results, "C", "OG0000900")
    # a per-class global tree, one level up from the per-OG dirs
    (phylo / "protein" / "class_A" / "class_A.treefile").write_text("(s1:0.1,s2:0.2);\n")
    # a legacy *_tree.tre, which the bash find also matches
    lse = phylo / "lse"
    lse.mkdir(parents=True)
    (lse / "level1_tree.tre").write_text("(s1:0.1,s2:0.2);\n")

    env = os.environ.copy()
    env["LOGS_DIR"] = str(tmp_path / "logs")
    (tmp_path / "logs").mkdir(exist_ok=True)
    bash_proc = subprocess.run(
        ["bash", "-c", f'PHYLO_DIR="{phylo}"\n{_extract_find_command()}\necho "$GENE_TREES"'],
        capture_output=True, text=True, env=env,
    )
    assert bash_proc.returncode == 0, bash_proc.stderr
    bash_trees = {Path(p).resolve() for p in bash_proc.stdout.split() if p}

    py_proc = subprocess.run(
        ["python3", "-c",
         "from pathlib import Path\n"
         f'phylo_dir = "{phylo}"\n'
         f"{_extract_python_tree_discovery()}\n"
         "print('\\n'.join(str(p) for p in gene_tree_files))"],
        capture_output=True, text=True,
    )
    assert py_proc.returncode == 0, py_proc.stderr
    py_trees = {Path(p).resolve() for p in py_proc.stdout.splitlines() if p}

    assert py_trees, "Python fallback still finds zero gene trees"
    assert py_trees == bash_trees, (
        f"backends disagree.\n  bash only: {sorted(bash_trees - py_trees)}"
        f"\n  python only: {sorted(py_trees - bash_trees)}"
    )


def test_python_fallback_is_recursive(tmp_path):
    """Pins the direction: a non-recursive glob is what made the fallback blind."""
    src = _extract_python_tree_discovery()
    assert "rglob" in src, "stage 03d's Python fallback must recurse like the bash find"


def test_treeshrink_rollback_copies_are_not_counted_as_gene_trees(tmp_path):
    """`<tree>.original.treefile` is the SAME tree, pre-TreeShrink.

    run_treeshrink.sh drops a verbatim rollback copy beside every tree it
    edits. It matches ``*.treefile``, so both backends would otherwise
    reconcile each shrunk orthogroup twice and double-count its duplications.
    """
    results = tmp_path / "results"
    phylo = results / "phylogenies"
    tree = _make_perog_tree(results, "A", "OG0000123")
    (tree.parent / "OG0000123.treefile.original.treefile").write_text("(s1:0.1,s2:0.2);\n")

    env = os.environ.copy()
    env["LOGS_DIR"] = str(tmp_path / "logs")
    (tmp_path / "logs").mkdir(exist_ok=True)
    bash_proc = subprocess.run(
        ["bash", "-c", f'PHYLO_DIR="{phylo}"\n{_extract_find_command()}\necho "$GENE_TREES"'],
        capture_output=True, text=True, env=env,
    )
    assert bash_proc.returncode == 0, bash_proc.stderr
    bash_trees = {Path(p).name for p in bash_proc.stdout.split() if p}
    assert bash_trees == {"OG0000123.treefile"}, bash_trees

    py_proc = subprocess.run(
        ["python3", "-c",
         "from pathlib import Path\n"
         f'phylo_dir = "{phylo}"\n'
         f"{_extract_python_tree_discovery()}\n"
         "print('\\n'.join(p.name for p in gene_tree_files))"],
        capture_output=True, text=True,
    )
    assert py_proc.returncode == 0, py_proc.stderr
    assert {n for n in py_proc.stdout.splitlines() if n} == {"OG0000123.treefile"}
