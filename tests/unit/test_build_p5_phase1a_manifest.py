"""Tests for scripts/build_p5_phase1a_manifest.py.

P5 Phase 1a validation infrastructure — manifest builder.

Reads proteome_manifest.tsv (Phase 1a, 91 annotated species), resolves
proteome cache paths via sanitize_sample_name convention, emits a TSV
with columns: taxid, binomial, proteome_path.

Idempotent: skips if output already exists unless --force.
"""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

import build_p5_phase1a_manifest as bpm


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

MANIFEST_COLUMNS = (
    "taxid", "binomial", "clade", "source", "accession",
    "assembly_level", "annotation_status", "est_protein_count",
    "submission_date", "drop_reason",
)


def _write_manifest(path: Path, rows: list[dict]) -> None:
    """Write a Phase 1a manifest TSV."""
    with path.open("w") as f:
        f.write("\t".join(MANIFEST_COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(r.get(c, "") for c in MANIFEST_COLUMNS) + "\n")


def _write_fasta(path: Path, n_seqs: int = 3) -> None:
    """Write a minimal FASTA file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        for i in range(n_seqs):
            f.write(f">seq{i}\nACGT\n")


# ---------------------------------------------------------------------------
# sanitize_sample_name
# ---------------------------------------------------------------------------

class TestSanitizeSampleName:
    def test_spaces_become_underscores(self) -> None:
        assert bpm.sanitize_sample_name("Helix pomatia") == "Helix_pomatia"

    def test_special_chars_become_underscores(self) -> None:
        assert bpm.sanitize_sample_name("Aplysia californica") == "Aplysia_californica"

    def test_leading_trailing_stripped(self) -> None:
        result = bpm.sanitize_sample_name(" Mytilus edulis ")
        assert not result.startswith("_")
        assert not result.endswith("_")

    def test_multiple_spaces_collapsed(self) -> None:
        result = bpm.sanitize_sample_name("A  B")
        assert "__" not in result


# ---------------------------------------------------------------------------
# read_phase1a_manifest
# ---------------------------------------------------------------------------

class TestReadPhase1aManifest:
    def test_returns_only_rows_with_accession(self, tmp_path: Path) -> None:
        mf = tmp_path / "manifest.tsv"
        _write_manifest(mf, [
            {"taxid": "100", "binomial": "Aaa bbb", "accession": "GCA_001"},
            {"taxid": "101", "binomial": "Ccc ddd", "accession": ""},      # no accession
            {"taxid": "102", "binomial": "Eee fff", "accession": "GCA_003"},
        ])
        rows = bpm.read_phase1a_manifest(mf)
        taxids = [r.taxid for r in rows]
        assert 100 in taxids
        assert 102 in taxids
        assert 101 not in taxids

    def test_skips_rows_with_no_taxid(self, tmp_path: Path) -> None:
        mf = tmp_path / "manifest.tsv"
        _write_manifest(mf, [
            {"taxid": "", "binomial": "Aaa bbb", "accession": "GCA_001"},
        ])
        rows = bpm.read_phase1a_manifest(mf)
        assert rows == []

    def test_row_has_taxid_and_binomial(self, tmp_path: Path) -> None:
        mf = tmp_path / "manifest.tsv"
        _write_manifest(mf, [
            {"taxid": "6183", "binomial": "Schistosoma mansoni", "accession": "GCF_000237925.1"},
        ])
        rows = bpm.read_phase1a_manifest(mf)
        assert len(rows) == 1
        assert rows[0].taxid == 6183
        assert rows[0].binomial == "Schistosoma mansoni"


# ---------------------------------------------------------------------------
# resolve_proteome_path
# ---------------------------------------------------------------------------

class TestResolveProteomePath:
    def test_returns_expected_path_for_known_species(self, tmp_path: Path) -> None:
        proteomes_dir = tmp_path / "proteomes"
        path = bpm.resolve_proteome_path(
            taxid=6183,
            binomial="Schistosoma mansoni",
            proteomes_dir=proteomes_dir,
        )
        # Should be <proteomes_dir>/<taxid>_<sanitized>.faa
        assert path.parent == proteomes_dir
        assert path.name == "6183_Schistosoma_mansoni.faa"

    def test_sanitized_binomial_used_in_filename(self, tmp_path: Path) -> None:
        proteomes_dir = tmp_path / "proteomes"
        path = bpm.resolve_proteome_path(
            taxid=9999,
            binomial="A. californica",
            proteomes_dir=proteomes_dir,
        )
        # Non-word chars sanitized in the stem (not the .faa suffix)
        stem = path.stem  # e.g. "9999_A_californica"
        assert "." not in stem
        assert path.suffix == ".faa"


# ---------------------------------------------------------------------------
# build_manifest
# ---------------------------------------------------------------------------

class TestBuildManifest:
    def test_emits_rows_for_existing_files(self, tmp_path: Path) -> None:
        mf = tmp_path / "manifest.tsv"
        proteomes_dir = tmp_path / "proteomes"
        _write_manifest(mf, [
            {"taxid": "100", "binomial": "Aaa bbb", "accession": "GCA_001"},
            {"taxid": "101", "binomial": "Ccc ddd", "accession": "GCA_002"},
        ])
        # Create proteome file for taxid 100 only
        p100 = proteomes_dir / "100_Aaa_bbb.faa"
        _write_fasta(p100)

        rows = bpm.build_manifest(mf, proteomes_dir=proteomes_dir)
        assert len(rows) == 1
        assert rows[0].taxid == 100

    def test_warns_for_missing_proteome_and_skips(
        self, tmp_path: Path, capsys
    ) -> None:
        mf = tmp_path / "manifest.tsv"
        proteomes_dir = tmp_path / "proteomes"
        _write_manifest(mf, [
            {"taxid": "200", "binomial": "Xxx yyy", "accession": "GCA_999"},
        ])
        # proteome not created → missing

        rows = bpm.build_manifest(mf, proteomes_dir=proteomes_dir)
        assert rows == []
        captured = capsys.readouterr()
        assert "200" in captured.err or "Xxx" in captured.err

    def test_emitted_row_contains_proteome_path(self, tmp_path: Path) -> None:
        mf = tmp_path / "manifest.tsv"
        proteomes_dir = tmp_path / "proteomes"
        _write_manifest(mf, [
            {"taxid": "300", "binomial": "Ppp qqq", "accession": "GCA_300"},
        ])
        faa = proteomes_dir / "300_Ppp_qqq.faa"
        _write_fasta(faa)

        rows = bpm.build_manifest(mf, proteomes_dir=proteomes_dir)
        assert len(rows) == 1
        assert rows[0].proteome_path == faa


# ---------------------------------------------------------------------------
# write_manifest_tsv + idempotency
# ---------------------------------------------------------------------------

class TestWriteManifestTsv:
    def test_output_has_header_row(self, tmp_path: Path) -> None:
        out = tmp_path / "out.tsv"
        bpm.write_manifest_tsv([], out)
        lines = out.read_text().splitlines()
        assert lines[0] == "taxid\tbinomial\tproteome_path"

    def test_output_rows_match_input(self, tmp_path: Path) -> None:
        out = tmp_path / "out.tsv"
        mf = tmp_path / "manifest.tsv"
        proteomes_dir = tmp_path / "proteomes"
        _write_manifest(mf, [
            {"taxid": "400", "binomial": "Alpha beta", "accession": "GCA_400"},
        ])
        faa = proteomes_dir / "400_Alpha_beta.faa"
        _write_fasta(faa)

        rows = bpm.build_manifest(mf, proteomes_dir=proteomes_dir)
        bpm.write_manifest_tsv(rows, out)

        with out.open() as f:
            reader = csv.DictReader(f, delimiter="\t")
            data = list(reader)
        assert len(data) == 1
        assert data[0]["taxid"] == "400"
        assert data[0]["binomial"] == "Alpha beta"
        assert data[0]["proteome_path"].endswith(".faa")

    def test_idempotent_skip_without_force(self, tmp_path: Path) -> None:
        out = tmp_path / "out.tsv"
        out.write_text("taxid\tbinomial\tproteome_path\n999\tFoo bar\t/some/path.faa\n")

        # Calling write_manifest_tsv with empty rows should NOT overwrite
        # when the file already exists (tested via main() CLI idempotency)
        # Direct function test: file is already present, main should skip
        rc = bpm.main([
            "--manifest", str(tmp_path / "nonexistent.tsv"),  # won't be read
            "--proteomes-dir", str(tmp_path),
            "--out", str(out),
        ])
        # Exit should be 0 and file should still have old content
        assert rc == 0
        content = out.read_text()
        assert "999" in content

    def test_force_overwrites_existing(self, tmp_path: Path) -> None:
        mf = tmp_path / "manifest.tsv"
        proteomes_dir = tmp_path / "proteomes"
        out = tmp_path / "out.tsv"
        _write_manifest(mf, [
            {"taxid": "500", "binomial": "New species", "accession": "GCA_500"},
        ])
        faa = proteomes_dir / "500_New_species.faa"
        _write_fasta(faa)
        # Write stale content
        out.write_text("taxid\tbinomial\tproteome_path\n999\told\t/old.faa\n")

        rc = bpm.main([
            "--manifest", str(mf),
            "--proteomes-dir", str(proteomes_dir),
            "--out", str(out),
            "--force",
        ])
        assert rc == 0
        content = out.read_text()
        assert "500" in content
        assert "999" not in content


# ---------------------------------------------------------------------------
# CLI end-to-end
# ---------------------------------------------------------------------------

class TestCLIEndToEnd:
    def test_main_returns_0_on_success(self, tmp_path: Path) -> None:
        mf = tmp_path / "manifest.tsv"
        proteomes_dir = tmp_path / "proteomes"
        out = tmp_path / "out.tsv"
        _write_manifest(mf, [
            {"taxid": "600", "binomial": "Gamma delta", "accession": "GCA_600"},
        ])
        _write_fasta(proteomes_dir / "600_Gamma_delta.faa")

        rc = bpm.main([
            "--manifest", str(mf),
            "--proteomes-dir", str(proteomes_dir),
            "--out", str(out),
        ])
        assert rc == 0
        assert out.exists()

    def test_main_returns_0_even_with_all_missing_proteomes(
        self, tmp_path: Path
    ) -> None:
        """All proteomes missing → 0 rows but still succeeds (just warn)."""
        mf = tmp_path / "manifest.tsv"
        out = tmp_path / "out.tsv"
        _write_manifest(mf, [
            {"taxid": "700", "binomial": "No file", "accession": "GCA_700"},
        ])
        rc = bpm.main([
            "--manifest", str(mf),
            "--proteomes-dir", str(tmp_path / "proteomes"),
            "--out", str(out),
        ])
        assert rc == 0
        # Header should still be there
        assert out.read_text().startswith("taxid\t")
