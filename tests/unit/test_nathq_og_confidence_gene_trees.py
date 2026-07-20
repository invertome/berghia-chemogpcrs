"""Bead g416: og_confidence must read Gene_Trees/, on the right SUPPORT SCALE.

THE ORIGINAL DEFECT
-------------------
``load_og_trees`` globbed OrthoFinder's ``Resolved_Gene_Trees/``, whose internal
labels are NODE NAMES ('n1', 'n10', ...) carrying no support at all. So
``get_og_confidence_score`` never found a support value and the axis reported
itself unavailable for essentially every candidate -- a weight-1 axis
contributing nothing, silently.

Measured on Unity: ``Gene_Trees/`` carries real support in 368 of 430 trees,
with 897 distinct values (1, 0.999, 0.998, ...). It is the only one of the two
directories that can answer the question this axis asks.

THE TRAP IN THE FIX
-------------------
Those values are PROPORTIONS (0-1). ``BOOTSTRAP_THRESHOLD`` is 70, on the 0-100
percent scale. A naive repoint therefore compares ``0.999 >= 70``, which fails
at EVERY node, and produces a constant 0.0 for every candidate -- the exact
inverse of the bug being fixed, equally silent, and now with a plausible-looking
"we measured it and support was poor" story attached. ``resolve_support_scale``
exists for precisely this: it decides the scale ONCE PER TREE (a single value is
genuinely ambiguous -- a bare '1' is 100% in an OrthoFinder gene tree, where
siblings read 0.999, but 1% in an IQ-TREE .contree).

The tests below drive the real parsing helpers over real Newick strings in both
scale conventions, rather than asserting on source text.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RANK = PROJECT_ROOT / "scripts" / "rank_candidates.py"
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from _rank_candidates_lib import (  # noqa: E402
    parse_support_label, resolve_support_scale,
)

BOOTSTRAP_THRESHOLD = 70.0

# Real OrthoFinder Gene_Trees support values (proportions).
ORTHOFINDER_PROPORTIONS = ["1", "0.999", "0.998", "0.87", "0.42"]
# Real IQ-TREE .contree support values (already percent).
IQTREE_PERCENTS = ["100", "97", "88", "42"]


def _resolve(labels):
    """Parse a whole tree's labels and return support on the 0-100 scale.

    ``None`` entries mean the label carried no support at all and must stay
    None -- never be defaulted.
    """
    parsed = [parse_support_label(x) for x in labels]
    undetermined = [v for p in parsed if p is not None and not p[1] for v in (p[0],)]
    mult = resolve_support_scale(undetermined)
    out = []
    for p in parsed:
        if p is None:
            out.append(None)
            continue
        value, is_percent = p
        out.append(float(value) if is_percent else float(value) * mult)
    return out


# ---------------------------------------------------------------------------
# 1. The repoint itself
# ---------------------------------------------------------------------------

def test_load_og_trees_globs_gene_trees_not_resolved_gene_trees() -> None:
    """Would have caught: the original glob, which pointed at the one directory
    of the two that provably carries no support."""
    src = RANK.read_text()
    m = re.search(r"Path\(results_dir\)\.glob\(\s*\n?\s*'([^']+)'", src)
    assert m, "load_og_trees' glob pattern not found"
    pattern = m.group(1)
    assert pattern.endswith("/Gene_Trees"), (
        f"og_confidence still reads {pattern!r}; Resolved_Gene_Trees labels "
        f"internal nodes 'n1','n2',... and carries no support")
    assert "Resolved_Gene_Trees" not in pattern


# ---------------------------------------------------------------------------
# 2. The scale trap -- the requirement that makes or breaks this change
# ---------------------------------------------------------------------------

def test_orthofinder_proportions_are_rescaled_before_the_threshold() -> None:
    """0.999 is 99.9% support and MUST clear BOOTSTRAP_THRESHOLD=70.

    Would have caught: a naive repoint that compares the raw proportion against
    70. Every node fails, og_confidence becomes a constant 0.0 for the whole
    cohort, and -- worse than the bug it replaced -- it looks like a measurement.
    """
    supports = _resolve(ORTHOFINDER_PROPORTIONS)
    assert supports[:3] == [pytest.approx(100.0), pytest.approx(99.9),
                            pytest.approx(99.8)]
    strong = [s for s in supports if s >= BOOTSTRAP_THRESHOLD]
    assert len(strong) == 4, (
        "OrthoFinder proportions were not rescaled; nodes at 0.999 support are "
        "being read as 0.999% and failing a threshold of 70")
    # ...and a genuinely weak node must still fail, or the rescale would have
    # simply made everything pass.
    assert supports[-1] == pytest.approx(42.0)
    assert supports[-1] < BOOTSTRAP_THRESHOLD


def test_percent_scaled_trees_are_left_alone() -> None:
    """The rescale must not fire on a tree that is already in percent.

    Would have caught: an unconditional x100, which turns 97% into 9700% and
    makes every node pass -- the mirror-image silent failure.
    """
    supports = _resolve(IQTREE_PERCENTS)
    assert supports == [pytest.approx(v) for v in (100.0, 97.0, 88.0, 42.0)]


def test_scale_is_decided_per_tree_not_globally() -> None:
    """A single bare value is genuinely ambiguous, so the decision has to be
    made from the whole tree -- and separately for each tree, because these two
    producers appear in the same run.

    Would have caught: resolving the scale once across all OG trees, where one
    percent-scaled tree would suppress the rescale for every proportion-scaled
    one.
    """
    assert resolve_support_scale([1.0, 0.999, 0.87]) == 100.0
    assert resolve_support_scale([100.0, 97.0, 1.0]) == 1.0
    # A lone '1' resolves differently depending on its tree -- which is exactly
    # why it cannot be decided label-by-label.
    assert _resolve(["1"]) == [pytest.approx(100.0)]
    assert _resolve(["1", "95"]) == [pytest.approx(1.0), pytest.approx(95.0)]


# ---------------------------------------------------------------------------
# 3. The 62 trees with no support
# ---------------------------------------------------------------------------

def test_unsupported_labels_stay_None_and_are_never_defaulted() -> None:
    """The ~62 of 430 Gene_Trees that carry no support must yield None, so their
    candidates report has_og_confidence_data = False.

    Would have caught: substituting a default (the historical ``return 1.0,
    True`` bug -- "no bootstrap data, don't penalize" -- which records "could
    not measure" as "scored perfectly, WITH evidence" and hands a weight-1 axis
    full marks on nothing).
    """
    assert parse_support_label("") is None
    assert parse_support_label(None) is None
    assert parse_support_label("n1") is None
    assert parse_support_label("n10") is None
    assert _resolve(["n1", "n2", ""]) == [None, None, None]


def test_get_og_confidence_score_reports_unavailable_on_an_unsupported_tree() -> None:
    """End-to-end on the real helper: no supports on the path -> (0.0, False).

    Would have caught: returning (0.0, True), which is a present zero -- it
    would drag the candidate's weighted mean down and count toward its evidence
    completeness on evidence that does not exist.
    """
    ete3 = pytest.importorskip("ete3")
    src = RANK.read_text()
    start = src.index("def get_og_confidence_score")
    end = src.index("\ndef ", start + 1)
    ns: dict = {}
    exec(compile(src[start:end], "<og_conf>", "exec"), ns)
    get_og_confidence_score = ns["get_og_confidence_score"]

    # OrthoFinder Resolved-style labels: node NAMES, no support anywhere.
    tree = ete3.Tree("((cand_1:0.1,ref_lse_A:0.1)n2:0.2,ref_lse_B:0.3)n1;",
                     format=1)
    for node in tree.traverse():
        node._parsed_support = None

    score, has_data = get_og_confidence_score(
        "cand_1", {"cand_1": "OG1"}, {"OG1": tree}, BOOTSTRAP_THRESHOLD)
    assert (score, has_data) == (0.0, False), (
        "a tree with no support must report the axis UNAVAILABLE, not 0.0")


def test_get_og_confidence_score_scores_a_supported_tree() -> None:
    """The positive control for the test above: with support present on the
    0-100 scale, the axis reports available and scores the fraction of
    well-supported nodes on the path."""
    ete3 = pytest.importorskip("ete3")
    src = RANK.read_text()
    start = src.index("def get_og_confidence_score")
    end = src.index("\ndef ", start + 1)
    ns: dict = {}
    exec(compile(src[start:end], "<og_conf>", "exec"), ns)
    get_og_confidence_score = ns["get_og_confidence_score"]

    tree = ete3.Tree("((cand_1:0.1,ref_lse_A:0.1)n2:0.2,ref_lse_B:0.3)n1;",
                     format=1)
    # 99.9% support, i.e. what an OrthoFinder '0.999' becomes after rescaling.
    for node in tree.traverse():
        node._parsed_support = 99.9 if not node.is_leaf() else None

    score, has_data = get_og_confidence_score(
        "cand_1", {"cand_1": "OG1"}, {"OG1": tree}, BOOTSTRAP_THRESHOLD)
    assert has_data is True
    assert score == pytest.approx(1.0), (
        "rescaled OrthoFinder support must clear BOOTSTRAP_THRESHOLD=70")


# ---------------------------------------------------------------------------
# 4. The axis is still quarantined
# ---------------------------------------------------------------------------

def test_the_repoint_does_not_reactivate_a_quarantined_axis() -> None:
    """og_confidence is orthogroup-derived, so it stays dormant under
    ORTHOLOGY_SOURCE_TRUSTED=0 regardless of this fix. The docstring must say
    so, or a reader will expect the axis to start voting.

    Would have caught: quietly restoring a Nath-derived axis to the ranking
    while fixing an unrelated correctness bug in it.
    """
    src = RANK.read_text()
    start = src.index("def load_og_trees")
    end = src.index("\ndef ", start + 1)
    body = src[start:end]
    assert "ORTHOLOGY_SOURCE_TRUSTED" in body, (
        "load_og_trees does not record that this axis stays dormant under the "
        "orthology quarantine")
