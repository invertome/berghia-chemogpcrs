"""``lse_divergence``'s threshold must be a percentile of the WHOLE candidate set.

Bead hf3u. ``collect_candidate_depths`` walked each candidate's path to the
root and ``break``-ed at the first node below ``BOOTSTRAP_THRESHOLD`` (70),
so a candidate was admitted to the percentile population only if EVERY node
on its path cleared 70. Measured on this repo's own 439-leaf class-A-style
tree (``preliminary/results/phylogenies/protein/v1/all_berghia_refs.treefile``)
on 2026-07-20:

    internal nodes with parsed support   437/437
    per-node pass rate at >= 70          0.8238
    median path length                   27 nodes
    candidates admitted                  0 of 439

P(all 27 nodes pass) ~ 0.824**27 ~ 0.5%, so the population was empty on every
run, ``np.percentile`` was never reached, and ``LSE_DIVERGENCE_PERCENTILE`` was
dead config. The hardcoded ``lse_threshold = 0.5`` that stood in for it is
not a neutral default -- it is INVERTED: 98.6% of candidates sit above 0.5,
where the design intent (75th percentile) is 25%. A weight-1 axis was
near-saturated, and backwards, in every ranking the project has produced.

WHY THE GATE IS DROPPED RATHER THAN SOFTENED TO A FRACTION
----------------------------------------------------------
The obvious repair is to converge on the fraction-of-path criterion the two
sibling functions already use (``get_og_confidence_score`` and
``path_bootstrap_confidence`` both return "fraction of path nodes with
bootstrap >= threshold"). Measured on the same real tree, that repair does
not work here, because on this tree support ANTI-correlates with depth:

    support fraction vs root-path node count   spearman = -0.6792
    support fraction vs patristic depth        spearman = -0.5903

Deep candidates have longer paths and therefore systematically weaker
whole-path support. Filtering the percentile population on support thus
preferentially discards DEEP candidates, drags the 75th percentile DOWN, and
re-inflates the pass fraction -- the very defect being fixed:

    filter            kept   p75      % of ALL candidates above
    none (correct)     439   6.2516        25.1%   <- design intent
    frac >= 0.80       261   5.0269        52.8%
    frac >= 0.85       223   4.9342        56.7%
    frac >= 0.90       151   4.7428        61.7%
    frac >= 1.00         0    --            --     <- the shipped behaviour

The all-or-nothing gate is not a special bug; it is the ``frac >= 1.00``
endpoint of a systematically biased family. Only the unfiltered population
reproduces the intended 25%. Branch support is a STATISTICAL quantity and
depth is a TOPOLOGICAL one; support belongs on the axes that already carry it
(``phylo``, ``og_confidence``), not fused into this one.

rank_candidates.py is not import-safe (top-level ``sys.argv``); these tests
exec single top-level function fragments, the pattern established in
tests/unit/test_rank_candidates_dnds_reliability.py.
"""
from __future__ import annotations

import ast
import re
from pathlib import Path
from typing import Any

import numpy as np
import pytest

ete3 = pytest.importorskip("ete3")
from ete3 import Tree  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RANK = PROJECT_ROOT / "scripts" / "rank_candidates.py"


def _exec_function(name: str) -> Any:
    """Exec one top-level function, bounded by its AST extent.

    Not a ``src.find("\\n\\ndef ")`` scan: ``collect_candidate_depths`` is
    followed by the module-level tree-loading block, which that scan would
    swallow and then fail on an undefined ``phylo_dir``.
    """
    src = RANK.read_text()
    tree = ast.parse(src)
    node = next((n for n in tree.body
                 if isinstance(n, ast.FunctionDef) and n.name == name), None)
    assert node is not None, f"{name} not found at module level in rank_candidates.py"
    fragment = ast.get_source_segment(src, node)
    ns: dict[str, Any] = {"np": np, "os": __import__("os")}
    exec(fragment, ns)
    return ns[name]


@pytest.fixture(scope="module")
def collect():
    return _exec_function("collect_candidate_depths")


@pytest.fixture(scope="module")
def depth_score():
    return _exec_function("get_lse_divergence_score")


def _caterpillar(n_leaves: int, weak_from: int = 0, weak_support: float = 10.0):
    """A ladder tree: leaf *i* hangs off the *i*-th internal node.

    Depth increases monotonically with *i*, so the leaf name encodes its own
    rank. Internal nodes at or beyond ``weak_from`` get ``weak_support``,
    reproducing the real tree's anti-correlation between depth and support.
    """
    root = Tree()
    root.name = "root"
    spine = root
    leaves = []
    for i in range(n_leaves):
        internal = spine.add_child(name=f"n{i}")
        internal.dist = 1.0
        internal._parsed_support = weak_support if i >= weak_from else 99.0
        leaf = internal.add_child(name=f"cand{i:03d}")
        leaf.dist = 0.0
        leaves.append(leaf.name)
        spine = internal
    return root, leaves


def test_one_weak_node_does_not_erase_the_candidate(collect):
    """A single sub-threshold node on a long path must not drop the candidate.

    This is the defect in miniature: the ``break`` made admission
    all-or-nothing over a median-27-node path.
    """
    tree, leaves = _caterpillar(30, weak_from=29)
    depths = collect(tree, leaves)
    assert len(depths) == 30, (
        "one weak node near the root erased %d of 30 candidates" % (30 - len(depths))
    )


def test_every_candidate_in_the_tree_enters_the_population(collect):
    """The real regression: uniformly imperfect support emptied the list.

    At the measured 82.4% per-node pass rate over 27-node paths, essentially
    no candidate survives; the shipped code returned 0 of 439.
    """
    tree, leaves = _caterpillar(40, weak_from=5)
    depths = collect(tree, leaves)
    assert len(depths) == 40, (
        "support-gated population kept %d of 40; the percentile is then taken "
        "over a biased subset (or, when empty, not at all)" % len(depths)
    )


def test_percentile_leaves_the_intended_quarter_above_threshold(collect):
    """The design intent: LSE_DIVERGENCE_PERCENTILE=75 => ~25% above threshold.

    Guards the end-to-end property the hardcoded 0.5 inverted (98.6% above).
    """
    tree, leaves = _caterpillar(400, weak_from=40)
    depths = collect(tree, leaves)
    threshold = np.percentile(depths, 75)
    above = sum(1 for d in depths if d > threshold) / len(depths)
    assert 0.20 <= above <= 0.30, (
        "%.1f%% of candidates above the 75th-percentile threshold; expected "
        "~25%%" % (100 * above)
    )


def test_deep_candidates_are_not_preferentially_discarded(collect):
    """Support-based filtering biases the population toward SHALLOW candidates.

    Measured on the real tree: mean depth of kept-at-frac>=0.85 is 4.262 vs
    5.834 for the discarded. A population that loses its deep tail cannot
    yield a meaningful upper-quartile depth threshold.
    """
    tree, leaves = _caterpillar(60, weak_from=30)
    depths = collect(tree, leaves)
    deepest = max(leaf_depth for leaf_depth in depths)
    # the deepest leaf sits 60 nodes down, every one of them weak-supported
    assert deepest >= 59.0, (
        "deepest candidate (weakest path) was discarded; max retained depth "
        "%.1f" % deepest
    )


def test_missing_candidates_do_not_enter_the_population(collect):
    """Only tree members are measured; unknown ids are absent, not zero."""
    tree, leaves = _caterpillar(10)
    depths = collect(tree, leaves + ["not_in_tree_1", "not_in_tree_2"])
    assert len(depths) == 10


def test_raw_depth_is_reported_even_below_threshold(depth_score):
    """``raw_tree_depth`` is a measurement; a sub-threshold depth is not 0.

    The scored value is correctly gated at the threshold, but zeroing the raw
    measurement too makes "shallow" indistinguishable from "never measured" in
    the emitted column -- the same present-but-unmeasured defect class.
    """
    tree, leaves = _caterpillar(10)
    score, raw = depth_score("cand002", tree, threshold=100.0)
    assert score == 0.0, "a depth below threshold must not score"
    assert raw == pytest.approx(3.0), (
        "raw depth reported as %r; the candidate sits 3.0 from the root and "
        "that was measured" % raw
    )


def test_threshold_has_no_hardcoded_fallback():
    """An underivable threshold must disable the axis, not invent a constant.

    The shipped fallback ``lse_threshold = 0.5`` fired on every run and put
    98.6% of candidates above threshold.
    """
    src = RANK.read_text()
    assert not re.search(r"^\s*lse_threshold\s*=\s*0\.5\s*$", src, re.M), (
        "rank_candidates.py still assigns a hardcoded lse_threshold fallback"
    )


def test_axis_reports_itself_unavailable_when_threshold_underivable():
    """There must be an availability flag distinct from tree membership.

    ``has_lse_data`` was aliased to ``has_phylo`` (tree membership), so a
    candidate in the tree claimed lse_divergence evidence even when no threshold
    could be derived at all.
    """
    src = RANK.read_text()
    assert "LSE_DIVERGENCE_AVAILABLE" in src, (
        "no flag records whether an lse_divergence threshold could be derived"
    )
    assert "has_lse_divergence_data" in src, (
        "lse_divergence availability is still aliased to has_phylo_data"
    )
    assert re.search(r"'lse_divergence',\s*'lse_divergence_score_norm',\s*'has_lse_divergence_data'", src), (
        "the composite's lse_divergence axis still gates on has_phylo_data"
    )
