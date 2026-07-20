"""Bead -i3w9: branch-support values must be parsed out of the Newick internal
LABEL, not read off ``node.support``.

WHY THIS TEST EXISTS
--------------------
``rank_candidates.py`` loads trees with ``Tree(tree_file, format=1)``. In ete3,
``format=1`` puts the internal label into ``node.name`` and leaves
``node.support`` at the *constructor default of 1.0*. Every consumer then read
support via ``getattr(node, 'support', None)`` and compared it to
``BOOTSTRAP_THRESHOLD = 70``, so every internal node scored 1.0 and every node
failed the threshold. Net effect: ``phylo_score`` was 0.0 for every candidate
(a weight-2 axis contributing zero discriminative power), the LSE depth
threshold silently fell back to a hardcoded 0.5, and ``og_confidence`` was a
constant 0.0.

WHERE THE FIXTURE FORMATS WERE VERIFIED
---------------------------------------
Every label form asserted below was measured from files that are actually in
this repository, on 2026-07-20, not assumed:

  * ``preliminary/results/phylogenies/protein/v1/all_berghia_refs.treefile``
    -> composite ``0.777/97``, ``1.000/100``, ``0.229/100``  (SH-aLRT / UFBoot).
    This is what production emits: ``04_phylogenetic_analysis.sh`` line 63 sets
    ``IQTREE_BOOT_FLAGS="-B ${IQTREE_BOOTSTRAP} -alrt 1000"``.
  * ``.../all_berghia_refs.contree``      -> bare integer ``100`` (UFBoot only, percent).
  * ``.../all_berghia_refs_fasttree.tre`` -> bare float ``0.768``, ``1.000`` (proportion).
  * ``.../orthogroups/input/OrthoFinder/Results_Feb24_1/Gene_Trees/*_tree.txt``
    -> bare ``1``, ``0.999``, ``0.998`` (proportion) and empty labels.
  * ``.../Results_Feb24_1/Resolved_Gene_Trees/*_tree.txt``
    -> ``n1``, ``n2``, ``n5`` ... these are NODE NAMES, NOT support values.
    Measured over 150 files: 6409 distinct labels, all of the form ``n<int>``.
    ``get_og_confidence_score`` globs *Resolved_Gene_Trees*, so the OG-confidence
    axis has no support data at all and must report that, not fabricate a score.
  * ``.../OG0000339.treefile`` -> no internal labels whatsoever.

The bare-``1`` case is the reason scale is inferred per tree rather than per
label: in an OrthoFinder Gene_Tree ``1`` means 1.0 == 100% (its siblings are
0.999, 0.998), while in an IQ-TREE ``.contree`` a bare integer is already a
percent. Composite IQ-TREE labels are exempt from the inference because the
trailing UFBoot field is always an integer percent.
"""
from __future__ import annotations

import os
from typing import Any

import pytest

from _rank_candidates_lib import (          # conftest puts scripts/ on sys.path
    parse_support_label,
    resolve_support_scale,
)

ete3 = pytest.importorskip("ete3")
Tree = ete3.Tree

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
PHYLO_V1 = os.path.join(REPO, "preliminary", "results", "phylogenies", "protein", "v1")
OF_RESULTS = os.path.join(
    REPO, "preliminary", "results", "orthogroups", "input", "OrthoFinder", "Results_Feb24_1"
)

IQTREE_TREEFILE = os.path.join(PHYLO_V1, "all_berghia_refs.treefile")
IQTREE_CONTREE = os.path.join(PHYLO_V1, "all_berghia_refs.contree")
FASTTREE_TRE = os.path.join(PHYLO_V1, "all_berghia_refs_fasttree.tre")


def _extract(*names: str) -> dict[str, Any]:
    """Exec the named top-level functions of rank_candidates.py in isolation.

    rank_candidates.py is a script with heavy import-time side effects; this is
    the same fragment-extraction trick test_rank_candidates_dnds_reliability.py
    uses, generalized to several interdependent functions.
    """
    src_path = os.path.join(REPO, "scripts", "rank_candidates.py")
    with open(src_path) as fh:
        src = fh.read()

    import numpy as np

    import _rank_candidates_lib as lib

    ns: dict[str, Any] = {
        "np": np,
        "Tree": Tree,
        "categorize_reference": lib.categorize_reference,
        "parse_support_label": lib.parse_support_label,
        "resolve_support_scale": lib.resolve_support_scale,
        "BOOTSTRAP_THRESHOLD": 70,
        "sys": __import__("sys"),
    }
    # Slice by AST line range, not by scanning for the next `def`: these helpers
    # sit directly above module-level code, so a text scan would swallow it.
    import ast

    lines = src.splitlines()
    bounds = {n.name: (n.lineno, n.end_lineno) for n in ast.parse(src).body
              if isinstance(n, ast.FunctionDef)}
    for name in names:
        assert name in bounds, "%s not found in rank_candidates.py" % name
        lo, hi = bounds[name]
        fragment = "\n".join(lines[lo - 1:hi])
        exec(compile(fragment, "rank_candidates.py:%s" % name, "exec"), ns)
    return ns


# --------------------------------------------------------------------------
# 0a. Module-level definition order.
# --------------------------------------------------------------------------

def test_module_level_helpers_are_defined_before_they_are_called() -> None:
    """rank_candidates.py is a SCRIPT: its module body runs top to bottom.

    ``annotate_tree_support`` is called from the module-level tree-loading block,
    so it must be defined above that block or the script dies with a NameError
    at import time. Every test in this file execs isolated function fragments
    and so is structurally blind to ordering -- this check is what covers it.
    """
    import ast

    src_path = os.path.join(REPO, "scripts", "rank_candidates.py")
    with open(src_path) as fh:
        tree = ast.parse(fh.read())

    defined_at = {n.name: n.lineno for n in tree.body
                  if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))}

    first_module_level_call: dict[str, int] = {}
    for stmt in tree.body:
        if isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            continue                      # bodies run later, not at import
        for node in ast.walk(stmt):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                name = node.func.id
                if name in defined_at:
                    first_module_level_call.setdefault(name, node.lineno)

    late = {name: (defined_at[name], line)
            for name, line in first_module_level_call.items()
            if defined_at[name] > line}
    assert not late, (
        "these functions are called at module level before they are defined "
        "(NameError at import): %s" % late)


# --------------------------------------------------------------------------
# 0. The defect itself, asserted against ete3 on the repo's own tree.
# --------------------------------------------------------------------------

def test_ete3_format1_does_not_populate_node_support() -> None:
    """Root-cause proof. ete3 format=1 leaves .support at 1.0 for EVERY node.

    This is what made ``getattr(node, 'support', None) >= 70`` unsatisfiable.
    If a future ete3 changes this, the fix below is still correct (it reads the
    label) but this test documents why the fix was needed.
    """
    t = Tree(IQTREE_TREEFILE, format=1)
    internals = [n for n in t.traverse() if not n.is_leaf()]
    assert len(internals) == 437
    assert {n.support for n in internals} == {1.0}
    assert sum(1 for n in internals if n.support >= 70) == 0
    # ...while the real support IS present, in .name:
    assert any("/" in n.name for n in internals)


# --------------------------------------------------------------------------
# 1. parse_support_label — every real label form.
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "label,expected",
    [
        # IQ-TREE `-B N -alrt N` composite: SH-aLRT / UFBoot. Use UFBoot, which
        # is what BOOTSTRAP_THRESHOLD=70 conventionally means, and which is
        # always an integer percent (so scale is known).
        ("0.777/97", (97.0, True)),
        ("1.000/100", (100.0, True)),
        ("0.229/100", (100.0, True)),
        # `-alrt N -abayes -B N` writes SH-aLRT/aBayes/UFBoot; UFBoot stays last.
        ("88.3/0.99/95", (95.0, True)),
        # Bare values: real number, scale not decidable from one label.
        ("97", (97.0, False)),
        ("100", (100.0, False)),
        ("1.000", (1.0, False)),
        ("0.768", (0.768, False)),
        ("1", (1.0, False)),
        # No support information at all.
        ("", None),
        ("   ", None),
        (None, None),
        # OrthoFinder Resolved_Gene_Trees node NAMES are not support values.
        ("n5", None),
        ("n137", None),
        # Genuinely unparseable -> None (degrade), never a silent default.
        ("0.9/", None),
        ("//", None),
        ("abc", None),
    ],
)
def test_parse_support_label_real_forms(label, expected) -> None:
    assert parse_support_label(label) == expected


def test_resolve_support_scale() -> None:
    # All <= 1.0 -> proportion scale, multiply by 100.
    assert resolve_support_scale([0.768, 1.0, 0.999]) == 100.0
    # Anything above 1.0 -> already percent.
    assert resolve_support_scale([97.0, 100.0, 1.0]) == 1.0
    assert resolve_support_scale([]) == 1.0


# --------------------------------------------------------------------------
# 2. annotate_tree_support on the repo's REAL trees.
# --------------------------------------------------------------------------

def test_annotate_real_iqtree_treefile_recovers_ufboot() -> None:
    ns = _extract("annotate_tree_support")
    t = Tree(IQTREE_TREEFILE, format=1)
    n_internal, n_with, n_labelled = ns["annotate_tree_support"](t, source="class_A")

    assert n_internal == 437
    assert n_with == 437, "every internal label in this file carries a UFBoot value"
    vals = [n._parsed_support for n in t.traverse() if not n.is_leaf()]
    assert all(0.0 <= v <= 100.0 for v in vals)
    # The whole point: nodes now actually clear BOOTSTRAP_THRESHOLD=70.
    n_pass = sum(1 for v in vals if v >= 70)
    assert n_pass > 0
    # ...and it is not a degenerate "everything passes" either.
    assert len(set(vals)) > 1


def test_annotate_real_contree_bare_integer_percent() -> None:
    """.contree carries bare integer UFBoot percentages (verified: `)100:`)."""
    ns = _extract("annotate_tree_support")
    t = Tree(IQTREE_CONTREE, format=1)
    _, n_with, _ = ns["annotate_tree_support"](t, source="contree")
    vals = [n._parsed_support for n in t.traverse()
            if not n.is_leaf() and n._parsed_support is not None]
    assert n_with > 0
    assert max(vals) == 100.0, "a bare 100 must stay 100, not become 10000"
    assert sum(1 for v in vals if v >= 70) > 0


def test_annotate_real_fasttree_proportions_are_rescaled() -> None:
    """FastTree writes local support as a 0-1 proportion (verified: `)0.768:`).

    Without rescaling, a perfectly supported FastTree node (1.000) would read as
    1% and fail the 70 threshold -- the same class of bug as the ete3 default.
    """
    ns = _extract("annotate_tree_support")
    t = Tree(FASTTREE_TRE, format=1)
    _, n_with, _ = ns["annotate_tree_support"](t, source="fasttree")
    vals = [n._parsed_support for n in t.traverse()
            if not n.is_leaf() and n._parsed_support is not None]
    assert n_with > 0
    assert max(vals) == pytest.approx(100.0)
    assert sum(1 for v in vals if v >= 70) > 0


def test_annotate_orthofinder_resolved_tree_reports_no_support() -> None:
    """Resolved_Gene_Trees label internal nodes 'n1','n2',... -- NOT support.

    Measured over 150 real files: 6409 distinct internal labels, all `n<int>`.
    The annotator must report zero support values (so the OG-confidence axis is
    marked unavailable) rather than coercing node names into numbers.
    """
    ns = _extract("annotate_tree_support")
    d = os.path.join(OF_RESULTS, "Resolved_Gene_Trees")
    if not os.path.isdir(d):
        pytest.skip("OrthoFinder Resolved_Gene_Trees not present")
    files = sorted(f for f in os.listdir(d) if f.endswith("_tree.txt"))
    if not files:
        pytest.skip("no resolved gene trees")
    t = Tree(os.path.join(d, files[0]), format=1)
    n_internal, n_with, n_labelled = ns["annotate_tree_support"](t, source=files[0])
    assert n_internal > 0
    assert n_labelled == n_internal, "labels ARE present, they are just node names"
    assert n_with == 0, "node names must not be read as support values"


def test_annotate_orthofinder_gene_tree_proportions() -> None:
    """Unresolved Gene_Trees DO carry support, as 0-1 proportions (1, 0.999...)."""
    ns = _extract("annotate_tree_support")
    d = os.path.join(OF_RESULTS, "Gene_Trees")
    if not os.path.isdir(d):
        pytest.skip("OrthoFinder Gene_Trees not present")
    files = sorted(f for f in os.listdir(d) if f.endswith("_tree.txt"))
    if not files:
        pytest.skip("no gene trees")
    t = Tree(os.path.join(d, files[0]), format=1)
    _, n_with, _ = ns["annotate_tree_support"](t, source=files[0])
    vals = [n._parsed_support for n in t.traverse()
            if not n.is_leaf() and n._parsed_support is not None]
    if n_with == 0:
        pytest.skip("this particular OG tree has no internal labels")
    assert max(vals) <= 100.0
    # A bare `1` here means 1.0 == 100%, not 1%.
    assert max(vals) == pytest.approx(100.0)


def test_annotate_unlabelled_tree_yields_no_support() -> None:
    """OG0000339.treefile has no internal labels at all (verified)."""
    ns = _extract("annotate_tree_support")
    t = Tree("((A:0.1,B:0.1):0.1,C:0.1);", format=1)
    n_internal, n_with, n_labelled = ns["annotate_tree_support"](t, source="bare")
    assert n_internal > 0
    assert n_labelled == 0
    assert n_with == 0


def test_annotate_format0_uses_native_support() -> None:
    """The OG loader falls back to ``Tree(..., format=0)``; there ete3 DOES put
    real support in .support and leaves .name empty. The annotator must honour
    that instead of reporting 'no support'."""
    ns = _extract("annotate_tree_support")
    t = Tree("((A:0.1,B:0.1)95:0.1,C:0.1);", format=0)
    _, n_with, _ = ns["annotate_tree_support"](t, source="og", newick_format=0)
    assert n_with > 0
    vals = [n._parsed_support for n in t.traverse()
            if not n.is_leaf() and n._parsed_support is not None]
    assert 95.0 in vals


# --------------------------------------------------------------------------
# 3. The scoring consumers.
# --------------------------------------------------------------------------

# Synthetic tree in the VERIFIED production label format (SH-aLRT/UFBoot, as
# emitted by `-B 1000 -alrt 1000`). Leaf names follow the ref_/candidate
# convention rank_candidates.py expects; the repo's own preliminary tree is
# Berghia-only (0 ref_ leaves) so it cannot exercise the ref-distance path.
SUPPORTED_NEWICK = (
    "(((cand_1:0.1,ref_or1a1:0.1)0.99/98:0.2,"
    "(cand_2:0.1,ref_gpcr_x:0.1)0.98/95:0.2)1.000/100:0.3,ref_out:0.5);"
)
UNSUPPORTED_NEWICK = (
    "(((cand_1:0.1,ref_or1a1:0.1)0.50/12:0.2,"
    "(cand_2:0.1,ref_gpcr_x:0.1)0.40/9:0.2)0.30/11:0.3,ref_out:0.5);"
)


def _phylo_ns():
    return _extract(
        "annotate_tree_support",
        "precompute_phylo_data",
        "fast_path_bootstrap_confidence",
        "weighted_distance_to_refs_fast",
    )


def test_phylo_score_is_not_identically_zero() -> None:
    """The headline defect: phylo_score was 0.0 for EVERY candidate.

    Pre-fix this asserts False -- ete3's 1.0 fails `>= 70`, confidence is 0.0,
    and weighted_distance_to_refs_fast multiplies by it.
    """
    ns = _phylo_ns()
    t = Tree(SUPPORTED_NEWICK, format=1)
    ns["annotate_tree_support"](t)
    refs = [l.name for l in t if l.name.startswith("ref_")]
    dist, supports, lookup = ns["precompute_phylo_data"](t, ["cand_1", "cand_2"], refs)
    scores = [ns["weighted_distance_to_refs_fast"](c, dist, lookup, supports)
              for c in ("cand_1", "cand_2")]
    assert all(s > 0 for s in scores)
    assert len(set(scores)) > 1, "the axis must discriminate between candidates"


def test_low_support_tree_scores_below_high_support_tree() -> None:
    """Support must actually modulate the score, in the right direction."""
    ns = _phylo_ns()
    out = {}
    for tag, nwk in (("hi", SUPPORTED_NEWICK), ("lo", UNSUPPORTED_NEWICK)):
        t = Tree(nwk, format=1)
        ns["annotate_tree_support"](t)
        refs = [l.name for l in t if l.name.startswith("ref_")]
        dist, sup, lookup = ns["precompute_phylo_data"](t, ["cand_1"], refs)
        out[tag] = ns["weighted_distance_to_refs_fast"]("cand_1", dist, lookup, sup)
    assert out["hi"] > out["lo"]


def test_bootstrap_confidence_reports_no_penalty_only_when_flagged() -> None:
    """A tree with NO support annotation must degrade explicitly, not silently.

    ``tree_has_support=False`` (set by the caller after the annotator reports
    zero parsed values, and warned about) means distance-only scoring. With
    ``tree_has_support=True`` an unparseable path node is a real evidence gap.
    """
    ns = _phylo_ns()
    t = Tree("((cand_1:0.1,ref_a:0.1):0.2,ref_b:0.5);", format=1)
    ns["annotate_tree_support"](t)
    dist, sup, lookup = ns["precompute_phylo_data"](t, ["cand_1"], ["ref_a", "ref_b"])
    fp = ns["fast_path_bootstrap_confidence"]
    assert fp("cand_1", "ref_a", lookup, sup, 70, tree_has_support=False) == 1.0
    assert fp("cand_1", "ref_a", lookup, sup, 70, tree_has_support=True) == 0.0


# --------------------------------------------------------------------------
# 4. The related trap: "could not measure" must not be recorded as "perfect".
# --------------------------------------------------------------------------

def test_og_confidence_returns_no_data_when_tree_has_no_support() -> None:
    """RELATED TRAP (bead -i3w9). ``get_og_confidence_score`` returned
    ``1.0, True`` when ``supports`` was empty -- "could not measure" recorded as
    "scored perfectly, WITH evidence". Today that is masked because ete3 always
    supplies 1.0 so the `>= 70` branch fires instead and it returns 0.0.
    Repairing the parse without repairing this converts a constant 0 into a
    constant 1, which is worse: it would push every OG-tree candidate to a
    full-weight perfect score on a fabricated axis.

    Real trigger: OrthoFinder Resolved_Gene_Trees (what the loader globs) label
    internal nodes 'n1','n2',... so there is never any support to read.
    """
    ns = _extract("annotate_tree_support", "get_og_confidence_score")
    t = Tree("((cand_1:0.1,ref_lse_1:0.1)n1:0.2,ref_lse_2:0.5)n2;", format=1)
    _, n_with, _ = ns["annotate_tree_support"](t, source="OG0000001")
    assert n_with == 0
    score, has_data = ns["get_og_confidence_score"](
        "cand_1", {"cand_1": "OG0000001"}, {"OG0000001": t}, 70)
    assert (score, has_data) == (0.0, False), (
        "no support data must be reported as unavailable, not as a perfect score")


def test_og_confidence_scores_when_support_is_present() -> None:
    ns = _extract("annotate_tree_support", "get_og_confidence_score")
    t = Tree("((cand_1:0.1,ref_lse_1:0.1)0.99/98:0.2,ref_lse_2:0.5)1.000/100;",
             format=1)
    _, n_with, _ = ns["annotate_tree_support"](t, source="OG0000002")
    assert n_with > 0
    score, has_data = ns["get_og_confidence_score"](
        "cand_1", {"cand_1": "OG0000002"}, {"OG0000002": t}, 70)
    assert has_data is True
    assert score == pytest.approx(1.0)

    lo = Tree("((cand_1:0.1,ref_lse_1:0.1)0.50/10:0.2,ref_lse_2:0.5)0.40/12;",
              format=1)
    ns["annotate_tree_support"](lo, source="OG0000003")
    score_lo, has_lo = ns["get_og_confidence_score"](
        "cand_1", {"cand_1": "OG0000003"}, {"OG0000003": lo}, 70)
    assert has_lo is True
    assert score_lo == pytest.approx(0.0)
    assert score_lo < score


# --------------------------------------------------------------------------
# 5. LSE depth threshold: LSE_DEPTH_PERCENTILE must stop being dead config.
# --------------------------------------------------------------------------

def test_lse_depth_collects_depths_on_a_supported_tree() -> None:
    """The ``all_depths`` loop broke on every candidate (ete3's 1.0 < 70), so
    ``np.percentile`` was never reached and the code fell back to a hardcoded
    ``lse_threshold = 0.5`` -- making LSE_DEPTH_PERCENTILE dead config.

    SUPERSEDED IN PART by bead hf3u. Parsing support correctly was necessary
    but not sufficient: at the real per-node pass rate (0.8238) over a median
    27-node path, requiring EVERY node to clear 70 still returned 0 of 439
    candidates. ``collect_candidate_depths`` no longer takes a support
    threshold at all -- the percentile population is every candidate in the
    tree, because support anti-correlates with depth (spearman -0.59) and any
    support filter biases the threshold toward shallow candidates. See
    tests/unit/test_lsedepth_percentile_population.py.

    What this test still guards is the i3w9 half: depths are collected from a
    tree whose labels parse, and they are real root distances.
    """
    ns = _extract("annotate_tree_support", "collect_candidate_depths")
    t = Tree(SUPPORTED_NEWICK, format=1)
    ns["annotate_tree_support"](t)
    depths = ns["collect_candidate_depths"](t, ["cand_1", "cand_2"])
    assert len(depths) == 2
    assert all(d > 0 for d in depths)

    # hf3u: a WEAKLY supported tree is still measured. Low support is a
    # statement about confidence in the topology, not about the depth of a
    # leaf within it, and it must no longer empty the percentile population.
    lo = Tree(UNSUPPORTED_NEWICK, format=1)
    ns["annotate_tree_support"](lo)
    assert len(ns["collect_candidate_depths"](lo, ["cand_1", "cand_2"])) == 2
