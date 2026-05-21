"""Tests for scripts/build_braker4_samples_csv.py.

Bead -p49 (Phase 1f BRAKER4). Generates a BRAKER4-format samples.csv for
the 134 unannotated-genome species, in protein-only EP mode using
OrthoDB metazoa as evidence.

BRAKER4 samples.csv schema (14 columns) per the upstream README
(Gaius-Augustus/BRAKER4):
  sample_name, genome, genome_masked, protein_fasta, bam_files,
  fastq_r1, fastq_r2, sra_ids, varus_genus, varus_species,
  isoseq_bam, isoseq_fastq, busco_lineage, reference_gtf

For EP mode: only sample_name, genome, protein_fasta, busco_lineage are
populated. genome_masked stays empty so BRAKER4 triggers its built-in
RepeatMasker masking step.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

import build_braker4_samples_csv as gen


BRAKER4_COLUMNS = (
    "sample_name", "genome", "genome_masked", "protein_fasta",
    "bam_files", "fastq_r1", "fastq_r2", "sra_ids",
    "varus_genus", "varus_species", "isoseq_bam", "isoseq_fastq",
    "busco_lineage", "reference_gtf",
)


PHASE1E_COLUMNS = (
    "taxid", "binomial", "clade", "source", "accession",
    "assembly_level", "annotation_status", "est_protein_count",
    "submission_date", "drop_reason", "contig_n50", "total_length_bp",
)


def _write_phase1e_manifest(path: Path, rows: list[dict]) -> None:
    with path.open("w") as f:
        f.write("\t".join(PHASE1E_COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(r.get(c, "") for c in PHASE1E_COLUMNS) + "\n")


# ----------------------------------------------------------------------
# sanitize_sample_name
# ----------------------------------------------------------------------

class TestSanitizeSampleName:
    """BRAKER4 sample_name must be a unique identifier — restrict to
    [a-zA-Z0-9_] so any per-species path it constructs works. Spaces,
    hyphens, apostrophes, etc. get normalized to underscores.
    """

    def test_simple_binomial(self) -> None:
        assert gen.sanitize_sample_name("6161_Dugesia japonica") == "6161_Dugesia_japonica"

    def test_already_clean(self) -> None:
        assert gen.sanitize_sample_name("6161_Dugesia_japonica") == "6161_Dugesia_japonica"

    def test_hyphenated_binomial(self) -> None:
        assert gen.sanitize_sample_name("12345_Sub-genus example") == "12345_Sub_genus_example"

    def test_apostrophe(self) -> None:
        assert gen.sanitize_sample_name("999_Genus o'brien") == "999_Genus_o_brien"

    def test_collapse_multiple_underscores(self) -> None:
        # Defensive: don't end up with `Genus__species` ladders.
        assert gen.sanitize_sample_name("9_Genus   species") == "9_Genus_species"

    def test_strip_leading_trailing(self) -> None:
        assert gen.sanitize_sample_name(" 123_x ") == "123_x"


# ----------------------------------------------------------------------
# build_samples_row
# ----------------------------------------------------------------------

class TestBuildSamplesRow:
    """One species → one row with the 14 BRAKER4 columns. For EP mode:
    sample_name, genome, protein_fasta, busco_lineage populated; all
    other 10 cols empty.
    """

    def test_ep_mode_row(self, tmp_path: Path) -> None:
        target = gen.SpeciesTarget(
            taxid=6161, binomial="Dugesia japonica",
            clade="platyhelminthes", accession="GCA_001938525.1",
        )
        genome_cache = tmp_path / "genomes"
        protein_db = tmp_path / "orthodb_metazoa.fa"
        # The genome FASTA must exist for the row to be emitted.
        genome_cache.mkdir()
        genome_fa = genome_cache / "6161_Dugesia_japonica.fasta"
        genome_fa.write_text(">chr1\nACGT\n")
        protein_db.write_text(">orthodb_p1\nMAA\n")

        row = gen.build_samples_row(
            target, genome_cache, protein_db,
            busco_lineage="metazoa_odb12",
        )

        assert row["sample_name"] == "6161_Dugesia_japonica"
        assert row["genome"] == str(genome_fa.resolve())
        assert row["protein_fasta"] == str(protein_db.resolve())
        assert row["busco_lineage"] == "metazoa_odb12"
        # All other 10 columns must be empty string (not missing)
        for col in (
            "genome_masked", "bam_files", "fastq_r1", "fastq_r2",
            "sra_ids", "varus_genus", "varus_species",
            "isoseq_bam", "isoseq_fastq", "reference_gtf",
        ):
            assert row[col] == "", f"expected empty {col}, got {row[col]!r}"

    def test_returns_none_when_genome_missing(self, tmp_path: Path) -> None:
        target = gen.SpeciesTarget(
            taxid=6161, binomial="Dugesia japonica",
            clade="platyhelminthes", accession="GCA_001938525.1",
        )
        genome_cache = tmp_path / "genomes"
        genome_cache.mkdir()  # but no genome FASTA inside
        protein_db = tmp_path / "p.fa"
        protein_db.write_text(">x\nM\n")
        row = gen.build_samples_row(
            target, genome_cache, protein_db, busco_lineage="metazoa_odb12",
        )
        assert row is None


# ----------------------------------------------------------------------
# main + CLI smoke test
# ----------------------------------------------------------------------

class TestMain:
    """End-to-end: manifest + genomes-on-disk + protein db → samples.csv."""

    def test_writes_samples_csv_with_correct_schema(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        # Prepare fixture state
        manifest = tmp_path / "genome_inventory_unannotated.tsv"
        _write_phase1e_manifest(manifest, [
            {"taxid": "6161", "binomial": "Dugesia japonica",
             "clade": "platyhelminthes", "accession": "GCA_001938525.1"},
            {"taxid": "6162", "binomial": "Girardia tigrina",
             "clade": "platyhelminthes", "accession": "GCA_001938485.1"},
        ])
        genome_cache = tmp_path / "genomes"
        genome_cache.mkdir()
        (genome_cache / "6161_Dugesia_japonica.fasta").write_text(">c\nACGT\n")
        (genome_cache / "6162_Girardia_tigrina.fasta").write_text(">c\nACGT\n")
        protein_db = tmp_path / "orthodb.fa"
        protein_db.write_text(">x\nM\n")
        out = tmp_path / "samples.csv"

        rc = gen.main([
            "--manifest", str(manifest),
            "--genome-cache", str(genome_cache),
            "--protein-db", str(protein_db),
            "--out", str(out),
        ])

        assert rc == 0
        # Validate CSV structure
        with out.open() as f:
            reader = csv.DictReader(f)
            assert reader.fieldnames == list(BRAKER4_COLUMNS)
            rows = list(reader)
        assert len(rows) == 2
        assert {r["sample_name"] for r in rows} == {
            "6161_Dugesia_japonica", "6162_Girardia_tigrina",
        }
        for r in rows:
            assert r["protein_fasta"] == str(protein_db.resolve())
            assert r["busco_lineage"] == "metazoa_odb12"

    def test_skips_species_without_cached_genome(
        self, tmp_path: Path,
    ) -> None:
        manifest = tmp_path / "m.tsv"
        _write_phase1e_manifest(manifest, [
            {"taxid": "6161", "binomial": "Dugesia japonica",
             "clade": "x", "accession": "GCA_1.1"},
            {"taxid": "6162", "binomial": "Girardia tigrina",
             "clade": "x", "accession": "GCA_2.1"},
        ])
        genome_cache = tmp_path / "g"
        genome_cache.mkdir()
        # Only one of the two species has a genome FASTA present.
        (genome_cache / "6161_Dugesia_japonica.fasta").write_text(">c\nACGT\n")
        protein_db = tmp_path / "p.fa"
        protein_db.write_text(">x\nM\n")
        out = tmp_path / "samples.csv"

        rc = gen.main([
            "--manifest", str(manifest),
            "--genome-cache", str(genome_cache),
            "--protein-db", str(protein_db),
            "--out", str(out),
        ])
        # rc=0 — missing genomes are a warning, not an error.
        assert rc == 0
        with out.open() as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) == 1
        assert rows[0]["sample_name"] == "6161_Dugesia_japonica"

    def test_empty_manifest_yields_header_only(self, tmp_path: Path) -> None:
        manifest = tmp_path / "m.tsv"
        manifest.write_text("\t".join(PHASE1E_COLUMNS) + "\n")
        genome_cache = tmp_path / "g"
        genome_cache.mkdir()
        protein_db = tmp_path / "p.fa"
        protein_db.write_text(">x\nM\n")
        out = tmp_path / "samples.csv"

        rc = gen.main([
            "--manifest", str(manifest),
            "--genome-cache", str(genome_cache),
            "--protein-db", str(protein_db),
            "--out", str(out),
        ])
        assert rc == 0
        with out.open() as f:
            lines = f.read().splitlines()
        assert len(lines) == 1  # header only
        assert lines[0] == ",".join(BRAKER4_COLUMNS)


class TestSubprocessSmoke:
    """As with -9dn / -p49 download: catches function-ordering bugs that
    pytest module imports mask. Runs the script via subprocess.
    """

    def test_script_runs_without_NameError(self, tmp_path: Path) -> None:
        import subprocess
        manifest = tmp_path / "m.tsv"
        manifest.write_text("\t".join(PHASE1E_COLUMNS) + "\n")
        genome_cache = tmp_path / "g"
        genome_cache.mkdir()
        protein_db = tmp_path / "p.fa"
        protein_db.write_text(">x\nM\n")
        out = tmp_path / "samples.csv"

        script = (
            Path(__file__).resolve().parent.parent.parent
            / "scripts" / "build_braker4_samples_csv.py"
        )
        result = subprocess.run(
            [
                "python", str(script),
                "--manifest", str(manifest),
                "--genome-cache", str(genome_cache),
                "--protein-db", str(protein_db),
                "--out", str(out),
            ],
            capture_output=True, text=True, timeout=30,
        )
        assert "NameError" not in result.stderr, result.stderr
        assert "ImportError" not in result.stderr, result.stderr
        assert result.returncode == 0, (
            f"rc={result.returncode}, stderr={result.stderr!r}"
        )
