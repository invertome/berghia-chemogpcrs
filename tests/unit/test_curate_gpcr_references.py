"""Tests for scripts/curate_gpcr_references.py.

Reference-curation script that queries UniProt's REST API for reviewed
(Swiss-Prot) GPCR entries per family, restricting to high-confidence
species (human, mouse, rat, Drosophila, C. elegans). Output: per-source
FASTA + provenance TSV.

Tests stub the UniProt API (no network in unit tests) and verify:
  - TSV response parsing handles edge cases (empty, header-only, normal)
  - Per-family query construction
  - Record canonicalization (family/subfamily/source labels)
  - Output file format (FASTA header convention, TSV columns)
  - Deduplication across queries (a receptor appearing in two families
    should only appear once in the consolidated output)
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

import curate_gpcr_references as cgr


# ---- parse_uniprot_tsv --------------------------------------------------

EXAMPLE_TSV = (
    "Entry\tReviewed\tProtein names\tOrganism\tGene Names\tLength\tSequence\n"
    "P28223\treviewed\t5-hydroxytryptamine receptor 2A (5-HT-2A)\tHomo sapiens (Human)\tHTR2A\t471\tMDILCEENTSL\n"
    "P28335\treviewed\t5-hydroxytryptamine receptor 2C\tHomo sapiens (Human)\tHTR2C\t458\tMVNLRNAVHSF\n"
)


def test_parse_uniprot_tsv_normal_response() -> None:
    rows = cgr.parse_uniprot_tsv(EXAMPLE_TSV)
    assert len(rows) == 2
    assert rows[0]["Entry"] == "P28223"
    assert rows[0]["Reviewed"] == "reviewed"
    assert rows[0]["Gene Names"] == "HTR2A"
    assert rows[0]["Sequence"].startswith("MDILCEENTSL")


def test_parse_uniprot_tsv_empty_response() -> None:
    """An empty UniProt response (no hits) should return [] without error."""
    assert cgr.parse_uniprot_tsv("") == []


def test_parse_uniprot_tsv_header_only_response() -> None:
    """Header but no data rows = no hits."""
    header_only = "Entry\tReviewed\tProtein names\tOrganism\tGene Names\tLength\tSequence\n"
    assert cgr.parse_uniprot_tsv(header_only) == []


def test_parse_uniprot_tsv_filters_unreviewed() -> None:
    """parse_uniprot_tsv should skip rows where Reviewed != 'reviewed' even
    if the API returns mixed results (defensive — our query asks for reviewed
    only but UniProt sometimes leaks unreviewed entries through)."""
    mixed = (
        "Entry\tReviewed\tProtein names\tOrganism\tGene Names\tLength\tSequence\n"
        "P28223\treviewed\tFoo\tHomo sapiens\tHTR2A\t471\tMDILCEENTSL\n"
        "A0A123\tunreviewed\tBar\tFoo bar\tXX1\t300\tMVKL\n"
    )
    rows = cgr.parse_uniprot_tsv(mixed)
    assert len(rows) == 1
    assert rows[0]["Entry"] == "P28223"


# ---- family_queries ---------------------------------------------------

def test_family_queries_returns_all_ten_coarse_families() -> None:
    queries = cgr.family_queries()
    expected_families = {
        "aminergic", "peptide", "opsin", "lipid", "nucleotide",
        "metabotropic-neurotransmitter", "glycoprotein-hormone",
        "class-B-secretin", "class-C", "class-F-frizzled",
    }
    assert set(queries.keys()) == expected_families


def test_family_queries_restrict_to_curated_species() -> None:
    """Each family query must restrict to human, mouse, rat, Drosophila, C. elegans
    via taxon IDs to avoid the auto-annotated invertebrate Swiss-Prot leakage."""
    queries = cgr.family_queries()
    # Taxon IDs: 9606 (human), 10090 (mouse), 10116 (rat), 7227 (D. mel),
    # 6239 (C. elegans). All queries should reference these taxa.
    for family, q in queries.items():
        # Every query must have the taxonomy filter
        assert ("9606" in q or "10090" in q or "10116" in q
                or "7227" in q or "6239" in q), (
            f"Query for {family!r} missing taxon-ID filter: {q}")
        # All queries must be reviewed-only
        assert "reviewed:true" in q, (
            f"Query for {family!r} not reviewed-only: {q}")


# ---- make_curated_record ------------------------------------------------

def test_make_curated_record_basic() -> None:
    uniprot_row = {
        "Entry": "P28223",
        "Reviewed": "reviewed",
        "Protein names": "5-hydroxytryptamine receptor 2A (5-HT-2A)",
        "Organism": "Homo sapiens (Human)",
        "Gene Names": "HTR2A HTR2",
        "Length": "471",
        "Sequence": "MDILCEENTSL",
    }
    rec = cgr.make_curated_record(uniprot_row, family="aminergic",
                                   subfamily="5HT")
    assert rec["accession"] == "P28223"
    assert rec["family"] == "aminergic"
    assert rec["subfamily"] == "5HT"
    assert rec["species"] == "Homo sapiens"
    assert rec["gene"] == "HTR2A"  # primary gene name only
    assert rec["source"] == "swissprot"
    assert rec["sequence"] == "MDILCEENTSL"


def test_make_curated_record_strips_organism_common_name() -> None:
    """Strip ' (Human)' / ' (Mouse)' / etc. from organism field — keep
    binomial only for stable downstream matching."""
    row = {
        "Entry": "P0DMS8", "Reviewed": "reviewed", "Protein names": "x",
        "Organism": "Drosophila melanogaster (Fruit fly)",
        "Gene Names": "Octb1R", "Length": "100", "Sequence": "MMM",
    }
    rec = cgr.make_curated_record(row, family="aminergic", subfamily="octopamine")
    assert rec["species"] == "Drosophila melanogaster"


def test_make_curated_record_handles_missing_subfamily() -> None:
    """Coarse-only families (opsin, lipid, etc.) have no medium subfamily
    — empty subfamily must be stored as empty string, not None."""
    row = {
        "Entry": "P12345", "Reviewed": "reviewed", "Protein names": "Opsin",
        "Organism": "Homo sapiens", "Gene Names": "RHO", "Length": "300",
        "Sequence": "MNG",
    }
    rec = cgr.make_curated_record(row, family="opsin", subfamily="")
    assert rec["subfamily"] == ""


# ---- output formats --------------------------------------------------

def test_write_fasta(tmp_path: Path) -> None:
    records = [
        {"accession": "P28223", "family": "aminergic", "subfamily": "5HT",
         "species": "Homo sapiens", "gene": "HTR2A", "source": "swissprot",
         "sequence": "MDILCEENTSL"},
        {"accession": "P0DMS8", "family": "aminergic", "subfamily": "octopamine",
         "species": "Drosophila melanogaster", "gene": "Octb1R",
         "source": "swissprot", "sequence": "MMMM"},
    ]
    out = tmp_path / "out.fasta"
    cgr.write_fasta(records, str(out))
    text = out.read_text()
    # Header convention: >accession|family|subfamily|species
    assert ">P28223|aminergic|5HT|Homo sapiens" in text
    assert ">P0DMS8|aminergic|octopamine|Drosophila melanogaster" in text
    assert "MDILCEENTSL" in text
    assert "MMMM" in text


def test_write_tsv(tmp_path: Path) -> None:
    records = [
        {"accession": "P28223", "family": "aminergic", "subfamily": "5HT",
         "species": "Homo sapiens", "gene": "HTR2A", "source": "swissprot",
         "sequence": "MDILCEENTSL"},
    ]
    out = tmp_path / "out.tsv"
    cgr.write_tsv(records, str(out))
    text = out.read_text()
    lines = text.strip().split("\n")
    assert lines[0].split("\t") == [
        "accession", "family", "subfamily", "species", "gene",
        "source", "length",
    ]
    assert lines[1].split("\t")[0] == "P28223"
    assert lines[1].split("\t")[1] == "aminergic"


# ---- end-to-end with stubbed API --------------------------------------

def test_curate_one_family_with_stubbed_api(tmp_path: Path,
                                            monkeypatch: pytest.MonkeyPatch) -> None:
    """End-to-end: a single family's curation produces correct outputs when
    the UniProt API is stubbed."""
    fake_response = (
        "Entry\tReviewed\tProtein names\tOrganism\tGene Names\tLength\tSequence\n"
        "P28223\treviewed\t5-hydroxytryptamine receptor 2A\tHomo sapiens (Human)\tHTR2A\t471\tMDILCEENTSL\n"
    )

    def fake_query(query: str, fields: list[str] | None = None) -> str:
        return fake_response

    monkeypatch.setattr(cgr, "query_uniprot", fake_query)

    records = cgr.curate_family("aminergic", "5HT", "fake_query")
    assert len(records) == 1
    assert records[0]["accession"] == "P28223"
    assert records[0]["family"] == "aminergic"
    assert records[0]["subfamily"] == "5HT"


def test_dedup_across_subqueries() -> None:
    """If a receptor appears in two subfamily queries (e.g. some receptors
    overlap aminergic/5HT and aminergic/melatonin queries due to query
    looseness), the consolidated output should only contain it once.
    Family/subfamily of the FIRST match wins."""
    a = {"accession": "P12345", "family": "aminergic", "subfamily": "5HT",
         "species": "Homo sapiens", "gene": "X", "source": "swissprot",
         "sequence": "M"}
    b = {"accession": "P12345", "family": "aminergic", "subfamily": "melatonin",
         "species": "Homo sapiens", "gene": "X", "source": "swissprot",
         "sequence": "M"}
    out = cgr.dedup_records([a, b])
    assert len(out) == 1
    assert out[0]["subfamily"] == "5HT"  # first match wins
