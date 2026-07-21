"""Tests for resolving anchor sequences from the (uniprot TSV + additions FASTA) pair.

The composite-id store ``anchor_set_FINAL_clean.fasta`` holds only the anchors
that predate the 2026-07-21 reference widening. The widening's additions were
written to a sidecar FASTA, so no single FASTA holds the whole reference set.
The production sources are the PAIR: ``anchor_set_PROD_uniprot.tsv`` (keyed by
``queried_accession``) and ``anchor_set_PROD_additions.fasta`` (keyed by the
composite id), which is the resolution order
``test_chemoreceptor_probe_separation`` already treats as truth.
"""
from pathlib import Path

import pytest

from extract_prod_anchors import resolve_anchor_sequences

_PROD_TSV = "\t".join(["accession", "tier", "taxid", "species", "family", "class", "evidence"]) + "\n" + \
    "\n".join([
        "\t".join(["P31356", "1", "9606", "Homo sapiens", "ADRB2", "A", "reviewed"]),
        "\t".join(["O15973", "9", "9606", "Homo sapiens", "GPR1", "A", "orthodb-harvest"]),
        "\t".join(["Q99999", "2", "10090", "Mus musculus", "mGluR", "C", "reviewed"]),
    ]) + "\n"

_UNIPROT_TSV = "\t".join(["queried_accession", "Entry", "Length", "Sequence"]) + "\n" + \
    "\n".join([
        "\t".join(["P31356", "P31356", "4", "MAAA"]),
        "\t".join(["Q99999", "Q99999", "4", "MCCC"]),
    ]) + "\n"

_ADDITIONS_FA = ">ANCHOR_A_9_O15973 harvested ortholog\nMBBB\n"


def _write(tmp_path: Path, name: str, text: str) -> str:
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def _sources(tmp_path: Path, uniprot=_UNIPROT_TSV, additions=_ADDITIONS_FA):
    return (_write(tmp_path, "prod.tsv", _PROD_TSV),
            _write(tmp_path, "uniprot.tsv", uniprot),
            _write(tmp_path, "additions.fa", additions))


def test_resolves_class_a_from_both_sources(tmp_path: Path):
    """The two sources compose: one anchor from the dump, one from the sidecar."""
    prod, up, add = _sources(tmp_path)
    got = resolve_anchor_sequences(prod, up, add, gpcr_class="A")
    assert got == {"ANCHOR_A_1_P31356": "MAAA", "ANCHOR_A_9_O15973": "MBBB"}


def test_other_classes_are_excluded(tmp_path: Path):
    prod, up, add = _sources(tmp_path)
    assert set(resolve_anchor_sequences(prod, up, add, gpcr_class="C")) == {"ANCHOR_C_2_Q99999"}


def test_unresolved_anchor_raises_and_names_it(tmp_path: Path):
    """An anchor absent from BOTH sources is an error, never a silent drop."""
    prod, up, add = _sources(tmp_path, additions=">ANCHOR_A_9_SOMETHINGELSE x\nMZZZ\n")
    with pytest.raises(ValueError, match="ANCHOR_A_9_O15973"):
        resolve_anchor_sequences(prod, up, add, gpcr_class="A")


def test_empty_sequence_cell_counts_as_unresolved(tmp_path: Path):
    """A present-but-blank Sequence must not resolve to an empty string.

    An empty sequence writes a FASTA record with a header and no residues, which
    parses cleanly and silently corrupts every downstream alignment.
    """
    up = "\t".join(["queried_accession", "Entry", "Length", "Sequence"]) + "\n" + \
         "\t".join(["P31356", "P31356", "0", "   "]) + "\n"
    prod, up_path, add = _sources(tmp_path, uniprot=up)
    with pytest.raises(ValueError, match="ANCHOR_A_1_P31356"):
        resolve_anchor_sequences(prod, up_path, add, gpcr_class="A")


def test_uniprot_dump_takes_precedence_over_additions(tmp_path: Path):
    """Resolution order is dump-then-sidecar, matching the probe-separation test."""
    prod, up, add = _sources(tmp_path, additions=(">ANCHOR_A_1_P31356 shadow\nMSHADOW\n"
                                                  ">ANCHOR_A_9_O15973 x\nMBBB\n"))
    got = resolve_anchor_sequences(prod, up, add, gpcr_class="A")
    assert got["ANCHOR_A_1_P31356"] == "MAAA"


def test_missing_additions_file_is_tolerated_when_dump_covers_everything(tmp_path: Path):
    """The sidecar is optional; the dump alone is a valid complete source."""
    prod_only_dump = "\t".join(["accession", "tier", "taxid", "species", "family", "class", "evidence"]) + "\n" + \
        "\t".join(["P31356", "1", "9606", "Homo sapiens", "ADRB2", "A", "reviewed"]) + "\n"
    prod = _write(tmp_path, "prod_one.tsv", prod_only_dump)
    up = _write(tmp_path, "uniprot.tsv", _UNIPROT_TSV)
    got = resolve_anchor_sequences(prod, up, str(tmp_path / "absent.fa"), gpcr_class="A")
    assert got == {"ANCHOR_A_1_P31356": "MAAA"}
