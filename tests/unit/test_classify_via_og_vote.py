"""Tests for scripts/classify_via_og_vote.py.

Phase 4 Task 4.2 — orthogroup-vote classifier. For each candidate:
look up its OG, get all OG members, vote family among annotated
references. Require ≥80% consensus AND ≥3 annotated members.
"""
from __future__ import annotations

from pathlib import Path

import classify_via_og_vote as ogv


# ---- load_reference_annotations -----------------------------------------

def test_load_reference_annotations(tmp_path: Path) -> None:
    """Read the curated reference TSV and return dict: accession -> (fam, sub)."""
    tsv = tmp_path / "refs.tsv"
    tsv.write_text(
        "accession\tfamily\tsubfamily\tspecies\tgene\tsource\tlength\n"
        "P28223\taminergic\t5HT\tHomo sapiens\tHTR2A\tswissprot\t471\n"
        "P12345\topsin\t\tHomo sapiens\tRHO\tswissprot\t348\n"
    )
    annotations = ogv.load_reference_annotations(str(tsv))
    assert annotations["P28223"] == ("aminergic", "5HT")
    assert annotations["P12345"] == ("opsin", "")


# ---- load_og_members (from add_og_coverage_columns; identical format) --

def test_load_og_members(tmp_path: Path) -> None:
    tsv = tmp_path / "Orthogroups.tsv"
    tsv.write_text(
        "Orthogroup\tspeciesA\tspeciesB\n"
        "OG0000001\tg1, g2\tg3\n"
        "OG0000002\tg4\t\n"
    )
    og = ogv.load_og_members(str(tsv))
    assert og["OG0000001"] == ["g1", "g2", "g3"]
    assert og["OG0000002"] == ["g4"]


# ---- vote_for_og --------------------------------------------------------

def test_vote_for_og_consensus_above_threshold() -> None:
    """If ≥80% of annotated OG members agree on a family, the OG votes
    that family."""
    members = ["g1", "g2", "g3", "g4", "g5"]
    annotations = {
        "g1": ("peptide", "tachykinin"),
        "g2": ("peptide", "tachykinin"),
        "g3": ("peptide", "opioid"),
        "g4": ("peptide", "opioid"),
        # g5 not annotated (no entry in the dict)
    }
    result = ogv.vote_for_og(members, annotations,
                              consensus_threshold=0.8, min_annotated=3)
    # 4 annotated, all 4 are 'peptide' → 100% peptide
    assert result["family"] == "peptide"
    assert result["consensus_fraction"] == 1.0
    assert result["n_annotated"] == 4


def test_vote_for_og_below_consensus_threshold() -> None:
    """Mixed-annotation OG should not be classified."""
    members = ["g1", "g2", "g3", "g4"]
    annotations = {
        "g1": ("aminergic", "5HT"),
        "g2": ("aminergic", "5HT"),
        "g3": ("peptide", "tachykinin"),
        "g4": ("opsin", ""),
    }
    # Top family (aminergic) has 2/4 = 50% < 80% threshold
    result = ogv.vote_for_og(members, annotations,
                              consensus_threshold=0.8, min_annotated=3)
    assert result["family"] == "unclassified-og"


def test_vote_for_og_below_min_annotated() -> None:
    """If <3 OG members are annotated, can't vote — even if all agree."""
    members = ["g1", "g2", "g3"]
    annotations = {
        "g1": ("aminergic", "5HT"),
        "g2": ("aminergic", "5HT"),
        # g3 not annotated → only 2 annotated members
    }
    result = ogv.vote_for_og(members, annotations,
                              consensus_threshold=0.8, min_annotated=3)
    assert result["family"] == "unclassified-og"
    assert result["n_annotated"] == 2


def test_vote_subfamily_when_consensus_clear() -> None:
    """If the family has clear consensus AND the subfamilies also agree,
    return the medium-granularity subfamily too."""
    members = ["g1", "g2", "g3", "g4"]
    annotations = {
        "g1": ("aminergic", "5HT"),
        "g2": ("aminergic", "5HT"),
        "g3": ("aminergic", "5HT"),
        "g4": ("aminergic", "5HT"),
    }
    result = ogv.vote_for_og(members, annotations,
                              consensus_threshold=0.8, min_annotated=3)
    assert result["family"] == "aminergic"
    assert result["subfamily"] == "5HT"


def test_vote_subfamily_falls_back_to_empty_on_disagreement() -> None:
    """If the family is consistent but subfamilies disagree, return
    the family with empty subfamily."""
    members = ["g1", "g2", "g3", "g4"]
    annotations = {
        "g1": ("aminergic", "5HT"),
        "g2": ("aminergic", "5HT"),
        "g3": ("aminergic", "dopamine"),
        "g4": ("aminergic", "dopamine"),
    }
    result = ogv.vote_for_og(members, annotations,
                              consensus_threshold=0.8, min_annotated=3)
    assert result["family"] == "aminergic"
    assert result["subfamily"] == ""  # subfamily disagreement


def test_vote_for_empty_og() -> None:
    """OG with no members at all (shouldn't happen in practice) returns
    unclassified."""
    result = ogv.vote_for_og([], {},
                              consensus_threshold=0.8, min_annotated=3)
    assert result["family"] == "unclassified-og"
    assert result["n_annotated"] == 0


# ---- classify_candidates_via_og -----------------------------------------

def test_classify_candidates_end_to_end(tmp_path: Path) -> None:
    """End-to-end: candidate IDs in -> classifications out."""
    refs_tsv = tmp_path / "refs.tsv"
    refs_tsv.write_text(
        "accession\tfamily\tsubfamily\tspecies\tgene\tsource\tlength\n"
        "P28223\taminergic\t5HT\tHomo sapiens\tHTR2A\tswissprot\t471\n"
        "P28335\taminergic\t5HT\tHomo sapiens\tHTR2C\tswissprot\t458\n"
        "P50406\taminergic\t5HT\tHomo sapiens\tHTR1A\tswissprot\t422\n"
        "P12345\topsin\t\tHomo sapiens\tRHO\tswissprot\t348\n"
    )
    og_tsv = tmp_path / "Orthogroups.tsv"
    og_tsv.write_text(
        "Orthogroup\tA\tB\tC\n"
        # OG1: 3 5HT receptors + 1 Berghia → strong aminergic/5HT vote
        "OG0000001\tP28223, P28335, P50406\tbste_gene_a\t\n"
        # OG2: just one Berghia, no annotated refs
        "OG0000002\tbste_gene_b\t\t\n"
        # OG3: opsin + Berghia (only 1 annotated member, below min)
        "OG0000003\tP12345\tbste_gene_c\t\n"
    )
    cands = ["bste_gene_a", "bste_gene_b", "bste_gene_c"]
    out_tsv = tmp_path / "out.tsv"

    ogv.classify_candidates_via_og(
        candidate_ids=cands,
        orthogroups_tsv=str(og_tsv),
        annotations_tsv=str(refs_tsv),
        output_tsv=str(out_tsv),
    )

    rows = list(open(out_tsv).read().strip().split("\n"))
    # header + 3 candidate rows
    assert len(rows) == 4
    # bste_gene_a -> aminergic/5HT (3/3 consensus)
    line_a = next(r for r in rows[1:] if "bste_gene_a" in r)
    assert "aminergic" in line_a and "5HT" in line_a
    # bste_gene_c -> unclassified-og (only 1 annotated member, below min=3)
    line_c = next(r for r in rows[1:] if "bste_gene_c" in r)
    assert "unclassified-og" in line_c
