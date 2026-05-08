"""Tests for scripts/classify_via_placement.py.

Phase 4 Task 4.3 — phylogenetic placement wrapper (EPA-ng).

Tests cover the local (no-subprocess, no-EPA-ng-binary-required) parts:
  - jplace JSON parsing (placements, edge_num, LWR per query)
  - edge-to-family mapping via reference tree + leaf-annotation TSV
  - best-placement selection per query
  - LWR threshold gating ('unclassified-placement' if best LWR < 0.80)
  - two-stage logic: backbone family call -> conditional subtree placement

The actual EPA-ng + mafft --add invocation is a thin subprocess wrapper
tested via integration on Unity (where EPA-ng is installed).
"""
from __future__ import annotations

import json
from pathlib import Path

import classify_via_placement as cvp


# ---- jplace parsing -----------------------------------------------------

def _make_jplace(tree: str, placements: list[dict]) -> dict:
    return {
        "tree": tree,
        "placements": placements,
        "metadata": {"epa-ng-version": "0.3.8"},
        "version": 3,
        "fields": ["edge_num", "likelihood", "like_weight_ratio",
                   "distal_length", "pendant_length"],
    }


def test_parse_jplace_extracts_best_placement_per_query(tmp_path: Path) -> None:
    """For each query in the jplace, return its best placement
    (highest LWR)."""
    jp = _make_jplace(
        tree="((A:0.1{0},B:0.2{1}):0.05{2},C:0.3{3});",
        placements=[
            {"p": [[0, -100.0, 0.85, 0.05, 0.02],
                   [1, -110.0, 0.10, 0.05, 0.02]],
             "nm": [["query_a", 1.0]]},
            {"p": [[3, -90.0, 0.95, 0.05, 0.02]],
             "nm": [["query_b", 1.0]]},
        ],
    )
    f = tmp_path / "out.jplace"
    f.write_text(json.dumps(jp))
    placements = cvp.parse_jplace(str(f))
    assert placements["query_a"]["edge_num"] == 0
    assert placements["query_a"]["lwr"] == 0.85
    assert placements["query_b"]["edge_num"] == 3
    assert placements["query_b"]["lwr"] == 0.95


def test_parse_jplace_handles_empty_placements(tmp_path: Path) -> None:
    jp = _make_jplace(tree="(A,B);", placements=[])
    f = tmp_path / "empty.jplace"
    f.write_text(json.dumps(jp))
    assert cvp.parse_jplace(str(f)) == {}


def test_parse_jplace_handles_missing_file(tmp_path: Path) -> None:
    """Gracefully handle missing jplace (e.g., EPA-ng failed)."""
    assert cvp.parse_jplace(str(tmp_path / "missing.jplace")) == {}


# ---- edge-to-family map -------------------------------------------------

def test_edge_to_family_leaf_edge() -> None:
    """An edge that connects to a leaf inherits that leaf's family."""
    leaf_annotations = {
        "P28223": ("aminergic", "5HT"),
        "P50406": ("aminergic", "5HT"),
        "P12345": ("opsin", ""),
    }
    # Tree: ((P28223:0.1{0},P50406:0.2{1}):0.05{2},P12345:0.3{3});
    # Edge 0 -> P28223 (leaf edge, aminergic/5HT)
    # Edge 3 -> P12345 (leaf edge, opsin/'')
    tree_newick = "((P28223:0.1{0},P50406:0.2{1}):0.05{2},P12345:0.3{3});"
    edge_map = cvp.build_edge_to_family_map(tree_newick, leaf_annotations)
    assert edge_map[0] == ("aminergic", "5HT")
    assert edge_map[1] == ("aminergic", "5HT")
    assert edge_map[3] == ("opsin", "")


def test_edge_to_family_internal_edge_monophyletic() -> None:
    """An internal edge subtending a monophyletic clade gets that clade's
    family. Subfamily empty if clade has multiple subfamilies."""
    leaf_annotations = {
        "X": ("aminergic", "5HT"),
        "Y": ("aminergic", "5HT"),
        "Z": ("opsin", ""),
    }
    tree_newick = "((X:0.1{0},Y:0.2{1}):0.05{2},Z:0.3{3});"
    edge_map = cvp.build_edge_to_family_map(tree_newick, leaf_annotations)
    # Edge 2 subtends (X, Y) — both aminergic/5HT; assign aminergic/5HT
    assert edge_map[2] == ("aminergic", "5HT")


def test_edge_to_family_pipe_format_leaves() -> None:
    """Tree leaves are written by select_backbone_reps.py as
    `accession|family|subfamily|species`; the leaf-annotation TSV is
    keyed by accession alone. The map builder must split the pipe
    suffix before looking up annotations, otherwise every edge gets
    `?unknown-leaf?` and no placement is ever classified.

    Regression for 2026-05-08 smoke postmortem (4th EPA-ng bug).
    """
    leaf_annotations = {
        "P28223": ("aminergic", "5HT"),
        "P50406": ("aminergic", "5HT"),
        "P12345": ("opsin", ""),
    }
    tree_newick = (
        "((P28223|aminergic|5HT|Homosapiens:0.1{0},"
        "P50406|aminergic|5HT|Homosapiens:0.2{1}):0.05{2},"
        "P12345|opsin||Homosapiens:0.3{3});"
    )
    edge_map = cvp.build_edge_to_family_map(tree_newick, leaf_annotations)
    assert edge_map[0] == ("aminergic", "5HT")
    assert edge_map[2] == ("aminergic", "5HT")
    assert edge_map[3] == ("opsin", "")


def test_edge_to_family_mixed_clade_returns_empty() -> None:
    """Internal edge subtending a polyphyletic clade returns empty
    family — placement on this edge is ambiguous."""
    leaf_annotations = {
        "X": ("aminergic", "5HT"),
        "Y": ("opsin", ""),
        "Z": ("opsin", ""),
    }
    tree_newick = "((X:0.1{0},Y:0.2{1}):0.05{2},Z:0.3{3});"
    edge_map = cvp.build_edge_to_family_map(tree_newick, leaf_annotations)
    # Edge 2 subtends (X aminergic, Y opsin) — mixed
    assert edge_map[2] == ("", "")


# ---- LWR threshold gating -----------------------------------------------

def test_classify_placement_passes_threshold() -> None:
    placements = {"q1": {"edge_num": 0, "lwr": 0.95}}
    edge_map = {0: ("aminergic", "5HT")}
    result = cvp.classify_placement(placements, edge_map, lwr_threshold=0.80)
    assert result["q1"]["family"] == "aminergic"
    assert result["q1"]["subfamily"] == "5HT"
    assert result["q1"]["lwr"] == 0.95


def test_classify_placement_below_threshold_unclassified() -> None:
    placements = {"q1": {"edge_num": 0, "lwr": 0.50}}
    edge_map = {0: ("aminergic", "5HT")}
    result = cvp.classify_placement(placements, edge_map, lwr_threshold=0.80)
    assert result["q1"]["family"] == "unclassified-placement"


def test_classify_placement_ambiguous_edge() -> None:
    """Even with high LWR, if the edge has no family (mixed clade),
    return unclassified-placement."""
    placements = {"q1": {"edge_num": 2, "lwr": 0.95}}
    edge_map = {2: ("", "")}
    result = cvp.classify_placement(placements, edge_map, lwr_threshold=0.80)
    assert result["q1"]["family"] == "unclassified-placement"


def test_classify_placement_query_not_placed() -> None:
    """Query not in placements (EPA-ng dropped it) -> unclassified."""
    placements: dict = {}
    edge_map = {0: ("aminergic", "5HT")}
    result = cvp.classify_placement_for_id("q1", placements, edge_map,
                                            lwr_threshold=0.80)
    assert result["family"] == "unclassified-placement"


# ---- leaf-annotation TSV loader -----------------------------------------

def test_load_leaf_annotations(tmp_path: Path) -> None:
    """Load per-leaf annotations from the tree's TSV (Phase 3 output)."""
    f = tmp_path / "backbone.tsv"
    f.write_text(
        "accession\tfamily\tsubfamily\tspecies\tgene\n"
        "P28223\taminergic\t5HT\tHomo sapiens\tHTR2A\n"
        "P12345\topsin\t\tHomo sapiens\tRHO\n"
    )
    out = cvp.load_leaf_annotations(str(f))
    assert out["P28223"] == ("aminergic", "5HT")
    assert out["P12345"] == ("opsin", "")


# ---- IQ-TREE model extraction (for EPA-ng --model) ----------------------

def test_read_iqtree_model_present(tmp_path: Path) -> None:
    """When sibling .iqtree exists, return the inferred model string."""
    (tmp_path / "tree.iqtree").write_text(
        "Some header\n"
        "Best-fit model according to BIC: LG+I+G4\n"
        "Model of substitution: LG+I+G4\n"
        "Other content\n"
    )
    (tmp_path / "tree.treefile").write_text("((A,B),C);")
    assert cvp._read_iqtree_model(str(tmp_path / "tree.treefile")) == "LG+I+G4"


def test_read_iqtree_model_missing_falls_back(tmp_path: Path) -> None:
    """No .iqtree file -> safe protein default (LG+G), not GTR+G."""
    (tmp_path / "tree.treefile").write_text("((A,B),C);")
    assert cvp._read_iqtree_model(str(tmp_path / "tree.treefile")) == "LG+G"


def test_read_iqtree_model_handles_contree(tmp_path: Path) -> None:
    """Works for .contree paths too (consensus tree variant)."""
    (tmp_path / "tree.iqtree").write_text("Model of substitution: VT+R8\n")
    (tmp_path / "tree.contree").write_text("((A,B),C);")
    assert cvp._read_iqtree_model(str(tmp_path / "tree.contree")) == "VT+R8"
