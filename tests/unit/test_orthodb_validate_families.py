"""Unit tests for the OrthoDB orthogroup-vs-curated-family agreement metrics.

These metrics decide whether OrthoDB is reported as usable, so they are tested
against hand-computable cases and cross-checked against scikit-learn's
reference implementations when it is installed.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "orthodb_validate_families",
    Path(__file__).resolve().parents[2] / "scripts" / "orthodb_validate_families.py",
)
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


# ------------------------------------------------------------- ARI: basics --

def test_ari_perfect_agreement_is_one():
    fams = ["opsin", "opsin", "peptide", "peptide"]
    ogs = ["OG1", "OG1", "OG2", "OG2"]
    assert mod.adjusted_rand_index(fams, ogs) == pytest.approx(1.0)


def test_ari_is_invariant_to_cluster_relabeling():
    fams = ["a", "a", "b", "b", "c", "c"]
    ogs1 = ["X", "X", "Y", "Y", "Z", "Z"]
    ogs2 = ["Z", "Z", "X", "X", "Y", "Y"]
    assert mod.adjusted_rand_index(fams, ogs1) == pytest.approx(
        mod.adjusted_rand_index(fams, ogs2)
    )


def test_ari_total_conflation_scores_zero():
    """Every family dumped into one orthogroup carries no information."""
    fams = ["opsin", "opsin", "peptide", "peptide", "aminergic", "aminergic"]
    ogs = ["OG_ALL"] * 6
    assert mod.adjusted_rand_index(fams, ogs) == pytest.approx(0.0)


def test_ari_total_fragmentation_scores_zero():
    """Every anchor in its own orthogroup also carries no information."""
    fams = ["opsin", "opsin", "peptide", "peptide", "aminergic", "aminergic"]
    ogs = [f"OG{i}" for i in range(6)]
    assert mod.adjusted_rand_index(fams, ogs) == pytest.approx(0.0)


def test_ari_penalises_partial_conflation():
    """Merging two families must score below keeping them apart."""
    fams = ["opsin"] * 3 + ["peptide"] * 3
    clean = ["A", "A", "A", "B", "B", "B"]
    merged = ["A"] * 6
    assert mod.adjusted_rand_index(fams, clean) > mod.adjusted_rand_index(fams, merged)


def test_ari_rejects_mismatched_lengths():
    with pytest.raises(ValueError):
        mod.adjusted_rand_index(["a", "b"], ["X"])


def test_ari_rejects_single_item():
    with pytest.raises(ValueError):
        mod.adjusted_rand_index(["a"], ["X"])


# ------------------------------ homogeneity / completeness map to the axes --

def test_homogeneity_is_one_when_no_orthogroup_mixes_families():
    fams = ["opsin", "opsin", "peptide", "peptide"]
    ogs = ["A", "B", "C", "D"]  # fragmented, but never mixed
    hom, comp, _ = mod.homogeneity_completeness_vmeasure(fams, ogs)
    assert hom == pytest.approx(1.0)
    assert comp < 1.0  # fragmentation is visible in completeness, not homogeneity


def test_completeness_is_one_when_each_family_is_contained():
    fams = ["opsin", "opsin", "opsin", "opsin"]
    ogs = ["A", "A", "A", "A"]
    hom, comp, _ = mod.homogeneity_completeness_vmeasure(fams, ogs)
    assert comp == pytest.approx(1.0)


def test_conflation_drives_homogeneity_down_not_completeness():
    """A level that conflates must be caught by homogeneity specifically."""
    fams = ["opsin", "opsin", "peptide", "peptide"]
    ogs = ["A", "A", "A", "A"]
    hom, comp, _ = mod.homogeneity_completeness_vmeasure(fams, ogs)
    assert hom == pytest.approx(0.0)
    assert comp == pytest.approx(1.0)


def test_vmeasure_is_harmonic_mean():
    fams = ["a", "a", "b", "b", "b", "c"]
    ogs = ["X", "Y", "Y", "Y", "Z", "Z"]
    hom, comp, v = mod.homogeneity_completeness_vmeasure(fams, ogs)
    assert v == pytest.approx(2 * hom * comp / (hom + comp))


def test_hcv_rejects_empty_input():
    with pytest.raises(ValueError):
        mod.homogeneity_completeness_vmeasure([], [])


# ------------------------------------------------ conflation / fragmentation --

def test_conflation_counts_mixed_orthogroups_and_the_anchors_in_them():
    pairs = [
        ("opsin", "OG1"),
        ("peptide", "OG1"),   # OG1 is mixed: 2 anchors
        ("aminergic", "OG2"),
        ("aminergic", "OG2"),  # OG2 is pure
        ("lipid", "OG3"),
    ]
    s = mod.conflation_stats(pairs)
    assert s["orthogroups_total"] == 3
    assert s["orthogroups_mixed_family"] == 1
    assert s["anchors_in_mixed_orthogroup"] == 2
    assert s["frac_anchors_in_mixed_orthogroup"] == pytest.approx(2 / 5)
    assert s["max_families_in_one_orthogroup"] == 2


def test_conflation_of_a_pure_partition_is_zero():
    pairs = [("opsin", "A"), ("opsin", "A"), ("peptide", "B")]
    s = mod.conflation_stats(pairs)
    assert s["orthogroups_mixed_family"] == 0
    assert s["frac_anchors_in_mixed_orthogroup"] == pytest.approx(0.0)


def test_fragmentation_normalises_by_family_size():
    """A big family in 3 groups is less fragmented than a small one in 3."""
    pairs = [("big", f"OG{i}") for i in range(3)]
    pairs += [("big", "OG0")] * 7          # 10 anchors, 3 orthogroups
    pairs += [("small", f"OG{i}") for i in range(10, 13)]  # 3 anchors, 3 groups
    s = mod.fragmentation_stats(pairs)
    assert s["per_family"]["big"]["orthogroups"] == 3
    assert s["per_family"]["big"]["orthogroups_per_anchor"] == pytest.approx(0.3)
    assert s["per_family"]["small"]["orthogroups_per_anchor"] == pytest.approx(1.0)


def test_fragmentation_perfectly_contained_family():
    pairs = [("opsin", "OG1")] * 5
    s = mod.fragmentation_stats(pairs)
    assert s["per_family"]["opsin"]["orthogroups"] == 1
    assert s["max_orthogroups_per_anchor"] == pytest.approx(0.2)


# ----------------------------------------------------------------- the join --

def _snap(acc, fam, primary=True):
    return {"accession": acc, "family": fam, "use_primary": str(primary)}


def test_assert_join_reports_mapping_rate():
    snapshot = [_snap("P1", "opsin"), _snap("P2", "peptide"), _snap("P3", "lipid")]
    gene_map = [
        {"accession": "P1", "odb_gene_id": "9606_0:001", "raw_xref": "P1"},
        {"accession": "P2", "odb_gene_id": "9606_0:002", "raw_xref": "P2"},
    ]
    rep = mod.assert_join(snapshot, gene_map)
    assert rep["anchors_mapped_any"] == 2
    assert rep["mapping_rate_class_a"] == pytest.approx(2 / 3)


def test_assert_join_flags_accessions_hitting_multiple_odb_genes():
    snapshot = [_snap("P1", "opsin")]
    gene_map = [
        {"accession": "P1", "odb_gene_id": "9606_0:001", "raw_xref": "P1"},
        {"accession": "P1", "odb_gene_id": "9606_0:002", "raw_xref": "P1"},
    ]
    rep = mod.assert_join(snapshot, gene_map)
    assert rep["anchors_mapping_to_multiple_odb_genes"] == 1


def test_assert_join_aborts_on_empty_intersection():
    """A join that silently means nothing is the defect this guards."""
    snapshot = [_snap("P1", "opsin")]
    with pytest.raises(SystemExit):
        mod.assert_join(snapshot, [])


def test_assert_join_aborts_when_gene_map_has_unknown_accessions():
    snapshot = [_snap("P1", "opsin")]
    gene_map = [{"accession": "ZZZ", "odb_gene_id": "9606_0:001", "raw_xref": "ZZZ"}]
    with pytest.raises(SystemExit):
        mod.assert_join(snapshot, gene_map)


def test_assert_join_excludes_non_primary_from_primary_rate():
    snapshot = [_snap("P1", "opsin"), _snap("P2", "orphan", primary=False)]
    gene_map = [{"accession": "P1", "odb_gene_id": "9606_0:001", "raw_xref": "P1"}]
    rep = mod.assert_join(snapshot, gene_map)
    assert rep["anchors_primary"] == 1
    assert rep["mapping_rate_primary"] == pytest.approx(1.0)


# ------------------------------------------------ cross-check with sklearn --

sklearn_metrics = pytest.importorskip(
    "sklearn.metrics", reason="scikit-learn not installed; skipping cross-check"
)


@pytest.mark.parametrize(
    "fams,ogs",
    [
        (["a", "a", "b", "b", "c", "c"], ["X", "X", "Y", "Z", "Z", "Z"]),
        (["a"] * 5 + ["b"] * 3 + ["c"] * 2, list("ABCDEFGHIJ")),
        (["a", "b", "a", "b", "a", "b"], ["X", "X", "X", "Y", "Y", "Y"]),
        (["opsin"] * 4 + ["peptide"] * 6, ["P"] * 10),
    ],
)
def test_ari_matches_sklearn(fams, ogs):
    assert mod.adjusted_rand_index(fams, ogs) == pytest.approx(
        sklearn_metrics.adjusted_rand_score(fams, ogs), abs=1e-9
    )


@pytest.mark.parametrize(
    "fams,ogs",
    [
        (["a", "a", "b", "b", "c", "c"], ["X", "X", "Y", "Z", "Z", "Z"]),
        (["a"] * 5 + ["b"] * 3 + ["c"] * 2, list("ABCDEFGHIJ")),
        (["opsin"] * 4 + ["peptide"] * 6, ["P"] * 10),
    ],
)
def test_hcv_matches_sklearn(fams, ogs):
    hom, comp, v = mod.homogeneity_completeness_vmeasure(fams, ogs)
    s_hom, s_comp, s_v = sklearn_metrics.homogeneity_completeness_v_measure(fams, ogs)
    assert hom == pytest.approx(s_hom, abs=1e-9)
    assert comp == pytest.approx(s_comp, abs=1e-9)
    assert v == pytest.approx(s_v, abs=1e-9)
