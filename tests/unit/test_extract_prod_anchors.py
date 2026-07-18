"""Tests for the class-specific anchor_set_PROD FASTA subsetter (A1 tree build)."""
from pathlib import Path

import pytest

from extract_prod_anchors import composite_anchor_ids, subset_fasta_by_ids

_PROD_TSV = "\t".join(["accession", "tier", "taxid", "species", "family", "class", "evidence"]) + "\n" + \
    "\n".join([
        "\t".join(["P31356", "1", "9606", "Homo sapiens", "ADRB2", "A", "reviewed"]),
        "\t".join(["O15973", "1", "9606", "Homo sapiens", "GPR1", "A", "reviewed"]),
        "\t".join(["Q99999", "2", "10090", "Mus musculus", "mGluR", "C", "reviewed"]),
        "\t".join(["P11111", "1", "7227", "Drosophila", "SMO", "F", "reviewed"]),
    ]) + "\n"


def _write(tmp_path: Path, name: str, text: str) -> str:
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def test_composite_ids_are_class_filtered_and_reconstructed(tmp_path: Path):
    tsv = _write(tmp_path, "prod.tsv", _PROD_TSV)
    ids = composite_anchor_ids(tsv, gpcr_class="A")
    assert ids == {"ANCHOR_A_1_P31356", "ANCHOR_A_1_O15973"}


def test_composite_ids_other_class(tmp_path: Path):
    tsv = _write(tmp_path, "prod.tsv", _PROD_TSV)
    assert composite_anchor_ids(tsv, gpcr_class="C") == {"ANCHOR_C_2_Q99999"}


def test_subset_fasta_selects_only_requested_ids(tmp_path: Path):
    # source FASTA has extra records + descriptions after the id token
    fasta = _write(tmp_path, "src.fa",
                   ">ANCHOR_A_1_P31356 beta-2 adrenergic\nMAAA\n"
                   ">ANCHOR_A_1_O15973 gpr1\nMBBB\n"
                   ">ANCHOR_C_2_Q99999 mglur\nMCCC\n")
    out = tmp_path / "sub.fa"
    n = subset_fasta_by_ids(fasta, {"ANCHOR_A_1_P31356", "ANCHOR_A_1_O15973"}, str(out))
    assert n == 2
    got = [ln[1:].split()[0] for ln in out.read_text().splitlines() if ln.startswith(">")]
    assert set(got) == {"ANCHOR_A_1_P31356", "ANCHOR_A_1_O15973"}
    # sequence content preserved
    assert "MAAA" in out.read_text() and "MCCC" not in out.read_text()


def test_subset_raises_if_requested_id_missing(tmp_path: Path):
    # ID integrity: a requested anchor absent from the source is an ERROR,
    # never silently dropped (write-once id contract).
    fasta = _write(tmp_path, "src.fa", ">ANCHOR_A_1_P31356 x\nMAAA\n")
    with pytest.raises(ValueError, match="missing"):
        subset_fasta_by_ids(fasta, {"ANCHOR_A_1_P31356", "ANCHOR_A_1_GHOST"},
                            str(tmp_path / "o.fa"))
