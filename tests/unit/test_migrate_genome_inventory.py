"""Tests for scripts/migrate_genome_inventory.py — lossless M1+M2 merge."""
from __future__ import annotations

from pathlib import Path

import migrate_genome_inventory as mig
from build_braker4_samples_csv import read_targets


def _write_m1(path: Path) -> None:
    path.write_text(
        "taxid\tbinomial\tclade\tsource\taccession\tassembly_level\t"
        "annotation_status\test_protein_count\tsubmission_date\tdrop_reason\t"
        "contig_n50\ttotal_length_bp\n"
        "100\tFoo bar\tgastropoda\tGenBank\tGCA_100\tScaffold\t\t0\t2020-01-01\t\t5000\t1000000\n"
    )


def _write_m2(path: Path) -> None:
    path.write_text(
        "taxid\tbinomial\tpolicy_class\tclade_name\tsource\taccession\t"
        "assembly_level\tannotation_status\test_protein_count\tsubmission_date\t"
        "contig_n50\ttotal_length_bp\n"
        "200\tBaz qux\theterobranchia\tHeterobranchia\tGenBank\tGCA_200\t"
        "Chromosome\t\t0\t2021-02-02\t90000\t2000000\n"
    )


def test_superset_header_and_provenance(tmp_path: Path) -> None:
    m1, m2, out = tmp_path / "m1.tsv", tmp_path / "m2.tsv", tmp_path / "u.tsv"
    _write_m1(m1); _write_m2(m2)
    mig.write_unified(out, mig.merge_manifests(m1, m2))
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == list(mig.UNIFIED_COLUMNS)
    rows = {r.split("\t")[0]: r.split("\t") for r in lines[1:]}
    hdr = list(mig.UNIFIED_COLUMNS)
    r100 = dict(zip(hdr, rows["100"]))
    r200 = dict(zip(hdr, rows["200"]))
    assert r100["source_batch"] == "nath_phase1e"
    assert r100["clade"] == "gastropoda" and r100["policy_class"] == ""
    assert r200["source_batch"] == "extension_phase1d"
    assert r200["clade"] == "Heterobranchia" and r200["policy_class"] == "heterobranchia"
    assert r200["drop_reason"] == ""


def test_sorted_by_taxid(tmp_path: Path) -> None:
    m1, m2, out = tmp_path / "m1.tsv", tmp_path / "m2.tsv", tmp_path / "u.tsv"
    _write_m1(m1); _write_m2(m2)
    mig.write_unified(out, mig.merge_manifests(m1, m2))
    taxids = [int(r.split("\t")[0]) for r in out.read_text().splitlines()[1:]]
    assert taxids == sorted(taxids)


def test_lossless_vs_two_manifest_build(tmp_path: Path) -> None:
    """The unified manifest must yield the identical SpeciesTarget set."""
    m1, m2, out = tmp_path / "m1.tsv", tmp_path / "m2.tsv", tmp_path / "u.tsv"
    _write_m1(m1); _write_m2(m2)
    mig.write_unified(out, mig.merge_manifests(m1, m2))
    two = [(t.taxid, t.binomial, t.clade, t.accession) for t in read_targets(m1, m2)]
    one = [(t.taxid, t.binomial, t.clade, t.accession) for t in read_targets(out)]
    assert two == one
