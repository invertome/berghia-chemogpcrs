"""Tests for scripts/build_full_scan_manifest.py.

P6 prep (bead kib): build the full-557 scan manifest by enumerating the
consolidated proteome directory (orthogroups_genome_wide/input/<taxid>_<binomial>.fa)
into the same 3-column manifest the scan array consumes:

    taxid  binomial  proteome_path
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

import build_full_scan_manifest as bfsm


def _touch_fasta(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(">x\nACGT\n")


# ---------------------------------------------------------------------------
# parse_leaf_filename
# ---------------------------------------------------------------------------

class TestParseLeafFilename:
    def test_taxid_and_binomial(self) -> None:
        assert bfsm.parse_leaf_filename("6182_Schistosoma_japonicum.fa") == (
            6182, "Schistosoma japonicum")

    def test_multiword_and_high_taxid(self) -> None:
        assert bfsm.parse_leaf_filename("1093978_Elysia_marginata.fa") == (
            1093978, "Elysia marginata")

    def test_non_taxid_prefix_returns_none(self) -> None:
        assert bfsm.parse_leaf_filename("README.fa") is None

    def test_taxid_only_no_binomial_returns_none(self) -> None:
        assert bfsm.parse_leaf_filename("6182.fa") is None


# ---------------------------------------------------------------------------
# enumerate_proteomes
# ---------------------------------------------------------------------------

class TestEnumerateProteomes:
    def test_enumerates_matching_and_skips_others(self, tmp_path: Path) -> None:
        d = tmp_path / "input"
        _touch_fasta(d / "6182_Schistosoma_japonicum.fa")
        _touch_fasta(d / "100452_Candidula_unifasciata.fa")
        _touch_fasta(d / "README.fa")               # no taxid → skipped
        (d / "notes.txt").write_text("ignore")      # wrong ext → not globbed
        rows = bfsm.enumerate_proteomes(d, pattern="*.fa")
        taxids = sorted(r.taxid for r in rows)
        assert taxids == [6182, 100452]

    def test_paths_are_absolute(self, tmp_path: Path) -> None:
        d = tmp_path / "input"
        _touch_fasta(d / "6182_Schistosoma_japonicum.fa")
        rows = bfsm.enumerate_proteomes(d, pattern="*.fa")
        assert rows[0].proteome_path.is_absolute()


# ---------------------------------------------------------------------------
# write + CLI
# ---------------------------------------------------------------------------

class TestWriteAndCLI:
    def test_output_format_sorted_lf_header(self, tmp_path: Path) -> None:
        d = tmp_path / "input"
        _touch_fasta(d / "100452_Candidula_unifasciata.fa")
        _touch_fasta(d / "6182_Schistosoma_japonicum.fa")
        out = tmp_path / "p6_manifest.tsv"
        rc = bfsm.main(["--proteomes-dir", str(d), "--out", str(out)])
        assert rc == 0

        raw = out.read_bytes()
        assert b"\r" not in raw, "manifest must be LF (scan wrapper cut -f3 breaks on CR)"
        with out.open() as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
            body = list(reader)
        assert header == ["taxid", "binomial", "proteome_path"]
        # sorted by numeric taxid: 6182 before 100452
        assert [r[0] for r in body] == ["6182", "100452"]
        assert body[0][1] == "Schistosoma japonicum"
        assert body[0][2].endswith("6182_Schistosoma_japonicum.fa")

    def test_idempotent_skip_without_force(self, tmp_path: Path) -> None:
        d = tmp_path / "input"
        _touch_fasta(d / "6182_Schistosoma_japonicum.fa")
        out = tmp_path / "p6_manifest.tsv"
        out.write_text("taxid\tbinomial\tproteome_path\n999\tOld sp\t/old.fa\n")
        rc = bfsm.main(["--proteomes-dir", str(d), "--out", str(out)])
        assert rc == 0
        assert "999" in out.read_text()      # not overwritten

    def test_empty_dir_writes_header_only(self, tmp_path: Path) -> None:
        d = tmp_path / "input"
        d.mkdir()
        out = tmp_path / "p6_manifest.tsv"
        rc = bfsm.main(["--proteomes-dir", str(d), "--out", str(out)])
        assert rc == 0
        assert out.read_text().rstrip("\n") == "taxid\tbinomial\tproteome_path"
