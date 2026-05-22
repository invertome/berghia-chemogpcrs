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


PHASE1D_COLUMNS = (
    "taxid", "binomial", "policy_class", "clade_name", "source",
    "accession", "assembly_level", "annotation_status",
    "est_protein_count", "submission_date", "contig_n50",
    "total_length_bp",
)


def _write_phase1d_manifest(path: Path, rows: list[dict]) -> None:
    with path.open("w") as f:
        f.write("\t".join(PHASE1D_COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(r.get(c, "") for c in PHASE1D_COLUMNS) + "\n")


# ----------------------------------------------------------------------
# read_targets — multi-manifest + clade_name column support
# ----------------------------------------------------------------------

class TestReadTargetsMultiManifest:
    """Bead -vqh: Phase 1d extension manifest must flow into the same
    samples.csv as Phase 1e, with per-taxid dedup (Phase 1e wins) and
    transparent handling of the `clade` vs `clade_name` column rename.
    """

    def test_single_phase1e_manifest_still_works(self, tmp_path: Path) -> None:
        m = tmp_path / "p1e.tsv"
        _write_phase1e_manifest(m, [
            {"taxid": "10", "binomial": "Aaa bbb", "clade": "gastropoda",
             "accession": "GCA_1"},
        ])
        targets = gen.read_targets(m)
        assert len(targets) == 1
        assert targets[0].taxid == 10
        assert targets[0].clade == "gastropoda"

    def test_phase1d_clade_name_column_recognized(self, tmp_path: Path) -> None:
        m = tmp_path / "p1d.tsv"
        _write_phase1d_manifest(m, [
            {"taxid": "20", "binomial": "Ccc ddd",
             "policy_class": "heterobranchia",
             "clade_name": "Heterobranchia", "accession": "GCA_2"},
        ])
        targets = gen.read_targets(m)
        assert len(targets) == 1
        assert targets[0].taxid == 20
        # The Phase 1d clade_name maps into SpeciesTarget.clade
        assert targets[0].clade == "Heterobranchia"

    def test_combines_phase1e_and_phase1d_dedup_first_wins(
        self, tmp_path: Path,
    ) -> None:
        m1 = tmp_path / "p1e.tsv"
        m2 = tmp_path / "p1d.tsv"
        _write_phase1e_manifest(m1, [
            {"taxid": "10", "binomial": "Phase1e sp", "clade": "p1e_clade",
             "accession": "GCA_E1"},
            {"taxid": "20", "binomial": "Other sp", "clade": "p1e_clade",
             "accession": "GCA_E2"},
        ])
        _write_phase1d_manifest(m2, [
            # taxid 10 overlaps with Phase 1e — Phase 1e must win
            {"taxid": "10", "binomial": "Phase1d sp",
             "policy_class": "heterobranchia",
             "clade_name": "Heterobranchia", "accession": "GCA_D1"},
            # taxid 30 is Phase 1d-only — must be picked up
            {"taxid": "30", "binomial": "New sp",
             "policy_class": "other_mollusca",
             "clade_name": "Bivalvia", "accession": "GCA_D3"},
        ])
        targets = gen.read_targets(m1, m2)
        by_taxid = {t.taxid: t for t in targets}
        assert set(by_taxid) == {10, 20, 30}
        # Phase 1e wins the dedup — clade and accession come from m1, not m2
        assert by_taxid[10].clade == "p1e_clade"
        assert by_taxid[10].accession == "GCA_E1"
        # New Phase 1d entry comes through cleanly
        assert by_taxid[30].clade == "Bivalvia"

    def test_empty_manifest_in_combination_does_no_harm(
        self, tmp_path: Path,
    ) -> None:
        m1 = tmp_path / "p1e.tsv"
        m2 = tmp_path / "p1d.tsv"
        _write_phase1e_manifest(m1, [
            {"taxid": "10", "binomial": "A b", "clade": "x",
             "accession": "GCA_1"},
        ])
        _write_phase1d_manifest(m2, [])  # header only
        targets = gen.read_targets(m1, m2)
        assert len(targets) == 1
        assert targets[0].taxid == 10

    def test_main_accepts_multiple_manifests(self, tmp_path: Path) -> None:
        m1 = tmp_path / "p1e.tsv"
        m2 = tmp_path / "p1d.tsv"
        _write_phase1e_manifest(m1, [
            {"taxid": "10", "binomial": "Aaa bbb", "clade": "p1e",
             "accession": "GCA_1"},
        ])
        _write_phase1d_manifest(m2, [
            {"taxid": "20", "binomial": "Ccc ddd",
             "policy_class": "heterobranchia",
             "clade_name": "Heterobranchia", "accession": "GCA_2"},
        ])
        genome_cache = tmp_path / "g"
        genome_cache.mkdir()
        (genome_cache / "10_Aaa_bbb.fasta").write_text(">x\nA\n")
        (genome_cache / "20_Ccc_ddd.fasta").write_text(">x\nA\n")
        protein_db = tmp_path / "p.fa"
        protein_db.write_text(">x\nM\n")
        out = tmp_path / "samples.csv"

        rc = gen.main([
            "--manifest", str(m1), str(m2),
            "--genome-cache", str(genome_cache),
            "--protein-db", str(protein_db),
            "--out", str(out),
        ])
        assert rc == 0
        with out.open() as f:
            rows = list(csv.DictReader(f))
        sample_names = sorted(r["sample_name"] for r in rows)
        assert sample_names == ["10_Aaa_bbb", "20_Ccc_ddd"]


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
