"""Bead hf3u: topological nesting depth is a SEPARATE axis from patristic depth.

WHY THIS AXIS EXISTS
--------------------
The ``lse_divergence`` axis is named for nesting depth inside a lineage-specific
expansion but measures ``node.get_distance(tree)`` -- cumulative branch length,
which is evolutionary RATE multiplied by TIME. For chemoreceptors specifically
those two quantities come apart, because chemoreceptor radiations are fast
evolving under diversifying selection: a fast-evolving shallow paralog can
outscore a slow-evolving deeply-nested one, which inverts the signal the axis
is supposed to contribute.

The pipeline is SUBTRACTIVE -- every other molluscan GPCR family is classified
and the unclassified residual is the enriched chemoreceptor set -- so there is
no positive control against which to arbitrate the two quantities directly.
Both are therefore emitted as independent voters and
``audit_signal_ranking_independence.py`` decides empirically (via the Spearman
grouping consumed by ``_load_signal_groups``) whether they are redundant.

MEASURED, NOT ASSUMED (2026-07-20, this repo's own trees)
---------------------------------------------------------
Population = every leaf; depth = ``get_distance(tree)`` vs nesting =
``len(get_ancestors())``; selection = above the 75th percentile of each.

  preliminary/results/phylogenies/protein/v1/all_berghia_refs.treefile (439 leaves)
      population spearman   +0.8973
      within the branch-length top quartile   +0.0782   (n=110)
      selected sets: 110 vs 89, 69 shared -> 62.73% of the scored set retained

  preliminary/results/phylogenies/protein/v2/gpcrs.treefile (2371 leaves)
      population spearman   +0.6038
      within the branch-length top quartile   +0.0339   (n=593)
      selected sets: 593 vs 593, 252 shared -> 42.50% of the scored set retained

Population-wide the two look nearly interchangeable; inside the only range that
can affect a ranking they are all but uncorrelated. That is the entire case for
carrying both.

WHY ``len(get_ancestors())`` AND NOT NODES-TO-THE-EXPANSION-MRCA
----------------------------------------------------------------
The alternative -- counting nodes up to the root of the maximal conspecific
(candidate-only) clade -- was measured on the same two trees and rejected:

  * v1 has no ``ref_`` leaves at all, so every leaf is a candidate and the
    measure degenerates to exactly ``len(get_ancestors())``.
  * v2 is 20.1% zero-inflated (a candidate whose sister is a reference gets 0)
    and ANTI-correlates with patristic depth (spearman -0.5016), because the
    measure is a function of which references happened to be sampled into the
    tree rather than of the candidate's own duplication history.

Reference sampling is exactly the dataset-specific dependence the
generalize-to-other-species constraint forbids. ``len(get_ancestors())`` is a
pure node count, defined for every leaf, needing no clade delineation, no
species labels and no external data. It inherits a dependence on the tree's
rooting, but that is the SAME dependence ``get_distance(tree)`` already has, so
it introduces no new one.
"""
from __future__ import annotations

import ast
import os
import re
from pathlib import Path
from typing import Any

import pytest

ete3 = pytest.importorskip("ete3")
Tree = ete3.Tree

REPO = Path(__file__).resolve().parent.parent.parent
RANK = REPO / "scripts" / "rank_candidates.py"
PHYLO_V1 = REPO / "preliminary" / "results" / "phylogenies" / "protein" / "v1"
PHYLO_V2 = REPO / "preliminary" / "results" / "phylogenies" / "protein" / "v2"
TREE_V1 = PHYLO_V1 / "all_berghia_refs.treefile"
TREE_V2 = PHYLO_V2 / "gpcrs.treefile"

# ((a,b),(c,d)) newick: cand_1 sits 2 nodes deep on a SHORT total branch length,
# cand_far sits 2 nodes deep on a LONG one, cand_shallow sits 1 node deep on a
# long one. Any implementation that returns branch length rather than a node
# count fails the ordering assertions below.
MIXED_NEWICK = (
    "(((cand_1:0.01,cand_2:0.01)n1:0.01,(cand_far:3.0,cand_deep2:3.0)n2:3.0)n3:0.1,"
    "cand_shallow:9.0)root;"
)


def _extract(*names: str) -> dict[str, Any]:
    """Exec named top-level functions of rank_candidates.py in isolation.

    rank_candidates.py is a script with heavy import-time side effects (it reads
    sys.argv at module level), so functions are sliced out by AST line range --
    the same trick tests/unit/test_i3w9_support_parsing.py uses.
    """
    import numpy as np

    src = RANK.read_text()
    ns: dict[str, Any] = {"np": np, "Tree": Tree, "sys": __import__("sys")}
    lines = src.splitlines()
    bounds = {n.name: (n.lineno, n.end_lineno) for n in ast.parse(src).body
              if isinstance(n, ast.FunctionDef)}
    for name in names:
        assert name in bounds, f"{name} not found in rank_candidates.py"
        lo, hi = bounds[name]
        exec(compile("\n".join(lines[lo - 1:hi]),
                     f"rank_candidates.py:{name}", "exec"), ns)
    return ns


def _rank_src() -> str:
    return RANK.read_text()


# --------------------------------------------------------------------------
# 1. The measurement itself, pinned on the repo's REAL trees.
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "tree_path,n_leaves,pop_rho,quartile_rho,retained",
    [
        (TREE_V1, 439, 0.8973, 0.0782, 0.6273),
        (TREE_V2, 2371, 0.6038, 0.0339, 0.4250),
    ],
)
def test_nesting_and_patristic_depth_diverge_inside_the_scored_quartile(
    tree_path: Path, n_leaves: int, pop_rho: float, quartile_rho: float,
    retained: float,
) -> None:
    """The two quantities agree population-wide and NOT where ranking happens.

    This is the measurement that justifies carrying both axes, so it is pinned
    rather than described: if a future change makes the topological measure
    track branch length again (for instance by reintroducing a length term),
    the within-quartile correlation rises and this fails.
    """
    import numpy as np
    from scipy.stats import spearmanr

    if not tree_path.exists():
        pytest.skip(f"real tree not present: {tree_path}")

    t = Tree(str(tree_path), format=1)
    leaves = list(t)
    assert len(leaves) == n_leaves

    ids = np.array([leaf.name for leaf in leaves])
    dist = np.array([leaf.get_distance(t) for leaf in leaves])
    nest = np.array([len(leaf.get_ancestors()) for leaf in leaves], dtype=float)

    assert spearmanr(dist, nest).correlation == pytest.approx(pop_rho, abs=5e-4)

    thr_d = np.percentile(dist, 75)
    thr_n = np.percentile(nest, 75)
    mask = dist > thr_d
    assert spearmanr(dist[mask], nest[mask]).correlation == pytest.approx(
        quartile_rho, abs=5e-4)

    sel_d, sel_n = set(ids[mask]), set(ids[nest > thr_n])
    assert len(sel_d & sel_n) / len(sel_d) == pytest.approx(retained, abs=5e-4)
    # ...and the disagreement is two-sided: candidates enter as well as leave.
    assert sel_n - sel_d, "no candidate is promoted by the topological measure"
    assert sel_d - sel_n, "no candidate is demoted by the topological measure"


# --------------------------------------------------------------------------
# 2. collect_candidate_nesting_depths: node counts, and ABSENT != zero.
# --------------------------------------------------------------------------

def test_nesting_depths_are_integer_node_counts_not_branch_lengths() -> None:
    ns = _extract("collect_candidate_nesting_depths")
    t = Tree(MIXED_NEWICK, format=1)
    depths = ns["collect_candidate_nesting_depths"](
        t, ["cand_1", "cand_far", "cand_shallow"])

    assert depths == [3, 3, 1], (
        "must be len(get_ancestors()) node counts; branch lengths would give "
        "floats and would rank cand_shallow (9.0) deepest")
    assert all(isinstance(d, int) for d in depths)


def test_absent_candidate_is_omitted_from_the_population_never_scored_zero() -> None:
    """The percentile population must contain only MEASURED candidates.

    Admitting an absent candidate as 0 would drag the threshold down and let
    unmeasured genes vote -- the defect class that produced the hardcoded
    ``lse_threshold = 0.5``.
    """
    ns = _extract("collect_candidate_nesting_depths")
    t = Tree(MIXED_NEWICK, format=1)
    depths = ns["collect_candidate_nesting_depths"](
        t, ["cand_1", "not_in_this_tree", "cand_shallow"])
    assert depths == [3, 1]


def test_nesting_population_is_every_candidate_in_the_tree() -> None:
    if not TREE_V1.exists():
        pytest.skip("real tree not present")
    ns = _extract("collect_candidate_nesting_depths")
    t = Tree(str(TREE_V1), format=1)
    names = [leaf.name for leaf in t]
    depths = ns["collect_candidate_nesting_depths"](t, names)
    assert len(depths) == 439, "no support filter may thin the population"
    assert min(depths) >= 1 and max(depths) == 39


# --------------------------------------------------------------------------
# 3. get_lse_nesting_depth_score: gated score, but raw ALWAYS reported.
# --------------------------------------------------------------------------

def test_subthreshold_candidate_reports_its_raw_nesting_depth() -> None:
    """Score 0.0 means "measured, shallow"; raw 0 would mean "never measured".

    Collapsing the two is exactly the bug the patristic axis had before hf3u.
    """
    ns = _extract("get_lse_nesting_depth_score")
    t = Tree(MIXED_NEWICK, format=1)

    score, raw = ns["get_lse_nesting_depth_score"]("cand_shallow", t, threshold=2)
    assert score == 0.0
    assert raw == 1, "raw nesting depth must survive the threshold gate"

    score, raw = ns["get_lse_nesting_depth_score"]("cand_1", t, threshold=2)
    assert score == 3 and raw == 3


def test_absent_candidate_scores_zero_zero() -> None:
    ns = _extract("get_lse_nesting_depth_score")
    t = Tree(MIXED_NEWICK, format=1)
    assert ns["get_lse_nesting_depth_score"]("nope", t, threshold=2) == (0.0, 0.0)


def test_infinite_threshold_scores_nobody_but_still_measures() -> None:
    """The unavailable path substitutes NO constant.

    When no threshold can be derived the threshold is +inf, so every candidate
    scores 0.0 while its raw depth is still reported, and the has_*_data flag
    (asserted separately below) is what tells the aggregator to drop the axis.
    """
    ns = _extract("get_lse_nesting_depth_score")
    t = Tree(MIXED_NEWICK, format=1)
    score, raw = ns["get_lse_nesting_depth_score"](
        "cand_1", t, threshold=float("inf"))
    assert score == 0.0 and raw == 3


# --------------------------------------------------------------------------
# 4. The two axes are SEPARATE voters, each with its own availability flag.
# --------------------------------------------------------------------------

def test_rank_aggregation_registers_both_depth_signals_separately() -> None:
    import rank_aggregation as ra

    keys = [spec[0] for spec in ra.SIGNAL_SPEC]
    assert "lse_divergence" in keys
    assert "lse_nesting_depth" in keys
    assert len(keys) == len(set(keys)), "duplicate signal key in SIGNAL_SPEC"

    flags = {spec[0]: spec[1] for spec in ra.SIGNAL_SPEC}
    assert flags["lse_nesting_depth"] == "has_lse_nesting_depth_data", (
        "the nesting axis must gate on its OWN measurement, not on the "
        "patristic axis's flag or on bare tree membership")


def test_both_depth_signals_produce_independent_ranklists() -> None:
    import pandas as pd

    import rank_aggregation as ra

    # cand_a is patristically deep but topologically shallow; cand_c is the
    # reverse. If the two axes were fused (or aliased) their ranklists would be
    # identical, and this asserts they are not.
    df = pd.DataFrame({
        "id": ["cand_a", "cand_b", "cand_c"],
        "lse_divergence_score": [9.0, 5.0, 1.0],
        "has_lse_divergence_data": [True, True, True],
        "lse_nesting_depth_score": [1.0, 5.0, 9.0],
        "has_lse_nesting_depth_data": [True, True, True],
    })
    rl = ra.build_ranklists_from_df(df)
    assert "lse_divergence" in rl and "lse_nesting_depth" in rl
    assert rl["lse_divergence"] != rl["lse_nesting_depth"]
    assert rl["lse_nesting_depth"]["cand_c"] > rl["lse_nesting_depth"]["cand_a"]


def test_nesting_axis_drops_out_when_it_could_not_be_measured() -> None:
    """Unavailable must mean ABSENT from the vote, never a substituted value."""
    import pandas as pd

    import rank_aggregation as ra

    df = pd.DataFrame({
        "id": ["cand_a", "cand_b"],
        "lse_divergence_score": [9.0, 1.0],
        "has_lse_divergence_data": [True, True],
        "lse_nesting_depth_score": [1.0, 9.0],
        "has_lse_nesting_depth_data": [False, False],
    })
    rl = ra.build_ranklists_from_df(df)
    assert "lse_divergence" in rl, "the patristic axis must be unaffected"
    assert "lse_nesting_depth" not in rl

    # A CSV with no nesting columns at all (every ranked CSV written before this
    # axis existed) must also skip it rather than fall back to voting.
    legacy = df.drop(columns=["lse_nesting_depth_score",
                              "has_lse_nesting_depth_data"])
    assert "lse_nesting_depth" not in ra.build_ranklists_from_df(legacy)
    assert "lse_nesting_depth" not in ra._OPTIONAL_FLAG_SIGNALS, (
        "a brand-new axis has no legacy CSVs to be lenient toward; a missing "
        "flag column means unmeasured, not 'predates the flag'")


# --------------------------------------------------------------------------
# 5. The independence audit must actually SEE both signals.
# --------------------------------------------------------------------------

def test_independence_audit_enumerates_both_depth_signals() -> None:
    """The audit is what decides empirically whether to fuse the two axes.

    It can only do that if both reach its correlation matrix, so its column
    list -- and the has_*_data flag it derives by suffix-swap -- are asserted
    here rather than assumed.
    """
    import audit_signal_ranking_independence as audit

    assert "lse_divergence_score" in audit.SIGNAL_COLUMNS
    assert "lse_nesting_depth_score" in audit.SIGNAL_COLUMNS
    # suffix-swap must land on the real flag names; no FLAG_OVERRIDES needed
    for col, expected in (("lse_divergence_score", "has_lse_divergence_data"),
                          ("lse_nesting_depth_score",
                           "has_lse_nesting_depth_data")):
        flag = audit.FLAG_OVERRIDES.get(col, "has_" + col.replace("_score", "_data"))
        assert flag == expected


def test_audit_masks_each_depth_signal_by_its_own_flag(tmp_path: Path) -> None:
    import audit_signal_ranking_independence as audit
    import numpy as np

    csv = tmp_path / "ranked.csv"
    csv.write_text(
        "id,lse_divergence_score,has_lse_divergence_data,"
        "lse_nesting_depth_score,has_lse_nesting_depth_data\n"
        "a,9.0,True,1.0,True\n"
        "b,5.0,True,5.0,False\n"
        "c,1.0,False,9.0,True\n"
    )
    m = audit.load_signal_matrix(str(csv))
    assert np.isnan(m.loc["c", "lse_divergence_score"])
    assert np.isnan(m.loc["b", "lse_nesting_depth_score"])
    assert m.loc["a", "lse_divergence_score"] == 9.0
    assert m.loc["a", "lse_nesting_depth_score"] == 1.0


# --------------------------------------------------------------------------
# 6. Registration at every enumeration site inside rank_candidates.py.
#    (the script is not import-safe, so these are source-level guards)
# --------------------------------------------------------------------------

def test_nesting_axis_has_its_own_weight_and_percentile_knobs() -> None:
    src = _rank_src()
    m = re.search(
        r"LSE_NESTING_DEPTH_WEIGHT\s*=\s*float\(\s*os\.getenv\("
        r"\s*['\"]LSE_NESTING_DEPTH_WEIGHT['\"]\s*,\s*([0-9.]+)\s*\)", src)
    assert m is not None, "LSE_NESTING_DEPTH_WEIGHT getenv default not found"
    assert float(m.group(1)) > 0.0, (
        "a weight of exactly 0 is read as an EXCLUSION by "
        "excluded_signals_from_weights, which would stop the new axis voting "
        "under the production RANK_METHOD=rankagg default")
    assert "LSE_NESTING_DEPTH_PERCENTILE" in src, (
        "the threshold must be a percentile (scale-free), not a constant")


def test_nesting_axis_registered_at_every_enumeration_site() -> None:
    src = _rank_src()
    # normalization, weights table, written columns, per-signal correlation
    for needle in (
        "'lse_nesting_depth_score',",
        "'lse_nesting_depth': LSE_NESTING_DEPTH_WEIGHT",
        "'raw_nesting_depth'",
        "'has_lse_nesting_depth_data'",
    ):
        assert needle in src, f"nesting axis not registered: {needle}"


def test_nesting_axis_gated_by_its_own_flag_in_the_production_scorer() -> None:
    src = _rank_src()
    start = src.find("def calculate_fair_rank_score(row):")
    assert start != -1
    end = src.find("\n\ndef ", start + 1)
    body = src[start:end] if end != -1 else src[start:]
    assert ("'lse_nesting_depth': row.get('lse_nesting_depth_score_norm') "
            "if row.get('has_lse_nesting_depth_data') else None") in body
    # the patristic axis must be untouched by this change
    assert ("'lse_divergence': row.get('lse_divergence_score_norm') "
            "if row.get('has_lse_divergence_data') else None") in body


def test_availability_is_derived_not_aliased() -> None:
    """``has_lse_nesting_depth_data`` must come from the nesting measurement.

    Aliasing it to ``has_lse_divergence_data`` would make the audit's two columns
    co-missing by construction and bias the very redundancy estimate this axis
    exists to let the audit make.
    """
    src = _rank_src()
    assert "LSE_NESTING_DEPTH_AVAILABLE" in src
    assert re.search(r"has_lse_nesting_depth\s*=\s*bool\(\s*has_phylo\s+and\s+"
                     r"LSE_NESTING_DEPTH_AVAILABLE\s*\)", src), \
        "has_lse_nesting_depth must be derived from LSE_NESTING_DEPTH_AVAILABLE"
    assert "has_lse_nesting_depth = has_lse_divergence" not in src
