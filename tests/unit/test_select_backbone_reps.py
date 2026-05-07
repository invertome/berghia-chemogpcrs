"""Tests for scripts/select_backbone_reps.py.

Phase 3 Task 3.1 — select N representative sequences per coarse family
for the backbone reference tree. Tests cover:
  - per-family quota (5-6 per family by default)
  - taxonomic diversity preference (1+ mammal AND 1+ invertebrate per family
    where the family has both available)
  - small-family fallback (use all available if family has < quota)
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

import select_backbone_reps as sbr


def _make_records(rows: list[dict]) -> list[dict]:
    return rows


# ---- per-family quota ---------------------------------------------------

def test_select_reps_respects_quota_per_family() -> None:
    records = (
        [{"accession": f"AMI{i}", "family": "aminergic", "subfamily": "5HT",
          "species": "Homo sapiens", "gene": f"HTR{i}"} for i in range(20)]
        + [{"accession": f"PEP{i}", "family": "peptide", "subfamily": "tachykinin",
            "species": "Homo sapiens", "gene": f"TACR{i}"} for i in range(15)]
    )
    selected = sbr.select_reps(records, quota_per_family=5)
    counts = {}
    for r in selected:
        counts[r["family"]] = counts.get(r["family"], 0) + 1
    assert counts["aminergic"] == 5
    assert counts["peptide"] == 5


def test_small_family_returns_all_available() -> None:
    """If a family has fewer than `quota` records, return all of them."""
    records = (
        [{"accession": "OPS1", "family": "opsin", "subfamily": "",
          "species": "Homo sapiens", "gene": "RHO"},
         {"accession": "OPS2", "family": "opsin", "subfamily": "",
          "species": "Homo sapiens", "gene": "OPN1"}]
    )
    selected = sbr.select_reps(records, quota_per_family=5)
    counts = {}
    for r in selected:
        counts[r["family"]] = counts.get(r["family"], 0) + 1
    assert counts["opsin"] == 2  # only 2 available


# ---- taxonomic diversity preference -------------------------------------

def test_prefers_mix_of_taxa_when_available() -> None:
    """Selection should include AT LEAST one invertebrate when family has
    both mammalian and invertebrate refs. The point is that the placement
    algorithm needs an invertebrate anchor in each family for our
    chemoreceptor candidates (mostly invertebrate) to land correctly."""
    records = [
        # 5 mammalian, 2 invertebrate aminergic refs
        {"accession": "M1", "family": "aminergic", "subfamily": "5HT",
         "species": "Homo sapiens", "gene": "HTR1"},
        {"accession": "M2", "family": "aminergic", "subfamily": "5HT",
         "species": "Mus musculus", "gene": "HTR1"},
        {"accession": "M3", "family": "aminergic", "subfamily": "dopamine",
         "species": "Rattus norvegicus", "gene": "DRD1"},
        {"accession": "M4", "family": "aminergic", "subfamily": "histamine",
         "species": "Homo sapiens", "gene": "HRH1"},
        {"accession": "M5", "family": "aminergic", "subfamily": "norepinephrine",
         "species": "Homo sapiens", "gene": "ADRB1"},
        {"accession": "I1", "family": "aminergic", "subfamily": "octopamine",
         "species": "Drosophila melanogaster", "gene": "OctR"},
        {"accession": "I2", "family": "aminergic", "subfamily": "tyramine",
         "species": "Caenorhabditis elegans", "gene": "tyra-2"},
    ]
    selected = sbr.select_reps(records, quota_per_family=5)
    # At least one invertebrate must be in the selection
    invertebrate_taxa = {"Drosophila melanogaster", "Caenorhabditis elegans"}
    has_invertebrate = any(
        r["species"] in invertebrate_taxa for r in selected)
    assert has_invertebrate, "Selection should include at least one invertebrate"


def test_prefers_diverse_subfamilies_within_family() -> None:
    """When picking from aminergic, prefer covering different subfamilies
    (5HT, dopamine, etc.) over picking 5 5HTs."""
    records = (
        [{"accession": f"HT{i}", "family": "aminergic", "subfamily": "5HT",
          "species": "Homo sapiens", "gene": "HTR"} for i in range(10)]
        + [{"accession": "DA1", "family": "aminergic", "subfamily": "dopamine",
            "species": "Homo sapiens", "gene": "DRD1"},
           {"accession": "NE1", "family": "aminergic", "subfamily": "norepinephrine",
            "species": "Homo sapiens", "gene": "ADRB1"},
           {"accession": "HI1", "family": "aminergic", "subfamily": "histamine",
            "species": "Homo sapiens", "gene": "HRH1"}]
    )
    selected = sbr.select_reps(records, quota_per_family=5)
    am = [r for r in selected if r["family"] == "aminergic"]
    sub_set = {r["subfamily"] for r in am}
    # We have 4 distinct subfamilies in input; selection should cover ≥3
    assert len(sub_set) >= 3, (
        f"Expected ≥3 distinct subfamilies in 5 picks, got {sub_set}")


# ---- TSV input parsing --------------------------------------------------

def test_load_records_from_tsv(tmp_path: Path) -> None:
    """Loader reads the curated reference TSV in the format produced by
    Task 1.1 / Task 1.6."""
    tsv = tmp_path / "refs.tsv"
    tsv.write_text(
        "accession\tfamily\tsubfamily\tspecies\tgene\tsource\tlength\n"
        "P28223\taminergic\t5HT\tHomo sapiens\tHTR2A\tswissprot\t471\n"
        "P12345\topsin\t\tHomo sapiens\tRHO\tswissprot\t348\n"
    )
    rows = sbr.load_records_from_tsv(str(tsv))
    assert len(rows) == 2
    assert rows[0]["accession"] == "P28223"
    assert rows[0]["family"] == "aminergic"
    assert rows[1]["subfamily"] == ""


# ---- write FASTA subset -------------------------------------------------

def test_write_subset_fasta_extracts_only_selected(tmp_path: Path) -> None:
    """write_subset_fasta reads the full reference FASTA and writes only
    the selected accessions, preserving the header convention."""
    full = tmp_path / "all.fasta"
    full.write_text(
        ">P28223|aminergic|5HT|Homo sapiens\nMDIL\n"
        ">P50406|aminergic|5HT|Mus musculus\nMSEK\n"
        ">P12345|opsin||Homo sapiens\nMNG\n"
    )
    selected = [
        {"accession": "P28223", "family": "aminergic", "subfamily": "5HT",
         "species": "Homo sapiens"},
        {"accession": "P12345", "family": "opsin", "subfamily": "",
         "species": "Homo sapiens"},
    ]
    out = tmp_path / "subset.fasta"
    sbr.write_subset_fasta(selected, str(full), str(out))
    text = out.read_text()
    assert "P28223" in text
    assert "P12345" in text
    assert "P50406" not in text
    assert "MDIL" in text  # sequence preserved
