"""Tests for scripts/build_classification_hmms.py.

Phase 2 Task 2.1 of non-chemoreceptor classification — build per-family
custom HMMs from the curated reference set, with leave-one-out
cross-validation in Task 2.2.

Tests cover the LOCAL (no-network) parts of the script:
  - parsing the reference TSV
  - grouping by family / subfamily with size thresholds
  - writing per-family input FASTAs
  - manifest TSV format

The actual mafft + hmmbuild invocations are tested via integration
(see test_build_classification_hmms_integration via subprocess in CI).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

import build_classification_hmms as bch


# ---- parse_reference_tsv ------------------------------------------------

def test_parse_reference_tsv_basic(tmp_path: Path) -> None:
    tsv = tmp_path / "refs.tsv"
    tsv.write_text(
        "accession\tfamily\tsubfamily\tspecies\tgene\tsource\tlength\n"
        "P28223\taminergic\t5HT\tHomo sapiens\tHTR2A\tswissprot\t471\n"
        "P28335\taminergic\t5HT\tHomo sapiens\tHTR2C\tswissprot\t458\n"
        "P12345\topsin\t\tHomo sapiens\tRHO\tswissprot\t348\n"
    )
    records = bch.parse_reference_tsv(str(tsv))
    assert len(records) == 3
    assert records[0]["accession"] == "P28223"
    assert records[0]["family"] == "aminergic"
    assert records[0]["subfamily"] == "5HT"
    assert records[2]["subfamily"] == ""


# ---- parse_reference_fasta ----------------------------------------------

def test_parse_reference_fasta_basic(tmp_path: Path) -> None:
    fa = tmp_path / "refs.fasta"
    fa.write_text(
        ">P28223|aminergic|5HT|Homo sapiens\nMDILCEENTSL\n"
        ">P12345|opsin||Homo sapiens\nMNGTEGPNF\n"
    )
    seqs = bch.parse_reference_fasta(str(fa))
    assert seqs["P28223"] == "MDILCEENTSL"
    assert seqs["P12345"] == "MNGTEGPNF"


# ---- group_by_family ----------------------------------------------------

def test_group_by_family_includes_all_records() -> None:
    records = [
        {"accession": "A", "family": "aminergic", "subfamily": "5HT"},
        {"accession": "B", "family": "aminergic", "subfamily": "dopamine"},
        {"accession": "C", "family": "opsin", "subfamily": ""},
        {"accession": "D", "family": "aminergic", "subfamily": "5HT"},
    ]
    groups = bch.group_by_family(records)
    assert sorted(groups["aminergic"]) == ["A", "B", "D"]
    assert groups["opsin"] == ["C"]


def test_group_by_subfamily_aminergic_peptide_only() -> None:
    """Subfamily groups are built ONLY for aminergic and peptide families
    (per design — coarse-only for other families)."""
    records = [
        {"accession": "A", "family": "aminergic", "subfamily": "5HT"},
        {"accession": "B", "family": "aminergic", "subfamily": "dopamine"},
        {"accession": "C", "family": "peptide", "subfamily": "tachykinin"},
        {"accession": "D", "family": "opsin", "subfamily": ""},
        {"accession": "E", "family": "lipid", "subfamily": ""},
    ]
    sub_groups = bch.group_by_subfamily(records)
    assert "aminergic/5HT" in sub_groups
    assert sub_groups["aminergic/5HT"] == ["A"]
    assert sub_groups["peptide/tachykinin"] == ["C"]
    # opsin/lipid should NOT be in subfamily groups (coarse-only)
    assert not any(k.startswith("opsin/") for k in sub_groups)
    assert not any(k.startswith("lipid/") for k in sub_groups)


def test_group_by_subfamily_skips_empty_subfamily() -> None:
    """Records with empty subfamily (e.g. aminergic/'' that didn't match
    a medium-granularity pattern) are NOT included in subfamily groups —
    they only contribute to the coarse-family HMM."""
    records = [
        {"accession": "A", "family": "aminergic", "subfamily": ""},
        {"accession": "B", "family": "aminergic", "subfamily": "5HT"},
    ]
    groups = bch.group_by_subfamily(records)
    assert groups == {"aminergic/5HT": ["B"]}


# ---- min_size filter ----------------------------------------------------

def test_filter_by_min_size_drops_small_groups() -> None:
    groups = {
        "aminergic": ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"],  # 10
        "peptide": ["X", "Y", "Z"],  # 3 -> dropped
        "opsin": [str(i) for i in range(15)],  # 15
    }
    out = bch.filter_by_min_size(groups, min_size=10)
    assert set(out.keys()) == {"aminergic", "opsin"}


# ---- write_family_fasta -------------------------------------------------

def test_write_family_fasta_preserves_input_order(tmp_path: Path) -> None:
    seqs = {"A": "MAAA", "B": "MBBB", "C": "MCCC"}
    accs = ["B", "A", "C"]  # NOT sorted
    out = tmp_path / "fam.fasta"
    bch.write_family_fasta(accs, seqs, str(out))
    text = out.read_text()
    # Verify accessions appear in input order
    pos_b = text.index(">B")
    pos_a = text.index(">A")
    pos_c = text.index(">C")
    assert pos_b < pos_a < pos_c


# ---- manifest TSV format -----------------------------------------------

def test_manifest_row_format() -> None:
    row = bch.manifest_row(
        family_label="aminergic", subfamily_label="",
        n_train=145, hmm_path="/path/aminergic.hmm",
        align_path="/path/aminergic.aln",
    )
    assert row["family"] == "aminergic"
    assert row["subfamily"] == ""
    assert row["n_train"] == 145
    assert row["hmm_path"] == "/path/aminergic.hmm"
    assert row["alignment_path"] == "/path/aminergic.aln"


def test_manifest_row_with_subfamily() -> None:
    row = bch.manifest_row(
        family_label="aminergic", subfamily_label="5HT",
        n_train=39, hmm_path="/p/5ht.hmm", align_path="/p/5ht.aln",
    )
    assert row["family"] == "aminergic"
    assert row["subfamily"] == "5HT"
