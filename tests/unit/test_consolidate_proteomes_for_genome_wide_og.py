"""Tests for scripts/consolidate_proteomes_for_genome_wide_og.py.

Covers the 5 required acceptance-criteria test groups:
  (a) sanitize_sample_name
  (b) rewrite_proteome — header rewriting + sequence passthrough
  (c) missing-proteome status path
  (d) idempotent re-run skips existing target
  (e) duplicate-taxid handling via consolidate_one + force flag

Import style matches test_build_braker4_samples_csv.py: conftest.py adds
scripts/ to sys.path so imports go through the bare module name.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

import consolidate_proteomes_for_genome_wide_og as cons


# ---------------------------------------------------------------------------
# (a) sanitize_sample_name
# ---------------------------------------------------------------------------

class TestSanitizeSampleName:
    """Rule must match build_braker4_samples_csv.sanitize_sample_name:
    collapse non-[A-Za-z0-9_] runs → '_', dedupe, strip leading/trailing '_'.
    """

    def test_space_replaced_by_underscore(self) -> None:
        assert cons.sanitize_sample_name("Berghia stephanieae") == "Berghia_stephanieae"

    def test_complex_name_with_punctuation(self) -> None:
        result = cons.sanitize_sample_name("Ctena cf. galapagana NHMUK 20210086")
        assert result == "Ctena_cf_galapagana_NHMUK_20210086"

    def test_leading_trailing_whitespace_stripped(self) -> None:
        assert cons.sanitize_sample_name("   leading_trailing   ") == "leading_trailing"

    def test_already_clean_unchanged(self) -> None:
        assert cons.sanitize_sample_name("1287507_Berghia_stephanieae") == "1287507_Berghia_stephanieae"

    def test_multiple_consecutive_spaces_collapsed(self) -> None:
        assert cons.sanitize_sample_name("Genus   species") == "Genus_species"

    def test_hyphen_becomes_underscore(self) -> None:
        assert cons.sanitize_sample_name("Sub-genus example") == "Sub_genus_example"


# ---------------------------------------------------------------------------
# (b) rewrite_proteome — header rewriting + sequence passthrough
# ---------------------------------------------------------------------------

class TestRewriteProteome:
    """Headers rewritten to >{leaf_name}|{original_id}; sequences untouched."""

    def test_two_sequence_fasta(self, tmp_path: Path) -> None:
        src = tmp_path / "in.faa"
        src.write_text(">XP_001 some description\nMSEQ\n>XP_002\nMOTH\n")
        dst = tmp_path / "out.fa"
        n = cons.rewrite_proteome(src, dst, "1287507_Berghia_stephanieae")
        assert n == 2
        content = dst.read_text()
        assert content.startswith(">1287507_Berghia_stephanieae|XP_001\n")
        assert ">1287507_Berghia_stephanieae|XP_002\n" in content
        assert "MSEQ" in content
        assert "MOTH" in content

    def test_description_stripped_from_header(self, tmp_path: Path) -> None:
        src = tmp_path / "in.faa"
        src.write_text(">PROT_1 [organism=Berghia stephanieae] extra\nMAA\n")
        dst = tmp_path / "out.fa"
        cons.rewrite_proteome(src, dst, "leaf")
        line = dst.read_text().splitlines()[0]
        # Only the first whitespace-delimited token of the original header survives
        assert line == ">leaf|PROT_1"

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        src = tmp_path / "in.faa"
        src.write_text(">P1\nM\n")
        dst = tmp_path / "deep" / "nested" / "out.fa"
        n = cons.rewrite_proteome(src, dst, "leaf")
        assert n == 1
        assert dst.exists()

    def test_returns_correct_count(self, tmp_path: Path) -> None:
        src = tmp_path / "in.faa"
        src.write_text(">A\nM\n>B\nM\n>C\nM\n")
        dst = tmp_path / "out.fa"
        assert cons.rewrite_proteome(src, dst, "leaf") == 3


# ---------------------------------------------------------------------------
# (c) missing-proteome status path
# ---------------------------------------------------------------------------

class TestMissingProteome:
    def test_missing_fasta_path_returns_missing_status(self, tmp_path: Path) -> None:
        src = cons.ProteomeSource(
            taxid=999,
            binomial="Imaginarius examplus",
            fasta_path=tmp_path / "nope.faa",
            phase="1a",
        )
        status = cons.consolidate_one(src, out_dir=tmp_path / "out")
        assert status.status == "missing_proteome"
        assert status.taxid == 999
        assert "nope.faa" in status.message

    def test_empty_fasta_returns_empty_status(self, tmp_path: Path) -> None:
        src_fa = tmp_path / "empty.faa"
        src_fa.write_text("")  # exists but has no sequences
        src = cons.ProteomeSource(
            taxid=42,
            binomial="Empty species",
            fasta_path=src_fa,
            phase="1a",
        )
        status = cons.consolidate_one(src, out_dir=tmp_path / "out")
        assert status.status == "empty_fasta"


# ---------------------------------------------------------------------------
# (d) idempotent re-run skips existing target
# ---------------------------------------------------------------------------

class TestIdempotentRerun:
    def test_second_run_skips_existing_file(self, tmp_path: Path) -> None:
        src_fa = tmp_path / "in.faa"
        src_fa.write_text(">a\nMS\n")
        src = cons.ProteomeSource(taxid=1, binomial="A b", fasta_path=src_fa, phase="1a")
        out_dir = tmp_path / "out"

        s1 = cons.consolidate_one(src, out_dir=out_dir)
        assert s1.status == "ok"
        assert s1.n_seqs == 1

        s2 = cons.consolidate_one(src, out_dir=out_dir)
        assert s2.status == "ok"
        assert "skipped" in s2.message.lower()

    def test_force_overwrites_existing(self, tmp_path: Path) -> None:
        src_fa = tmp_path / "in.faa"
        src_fa.write_text(">a\nMS\n")
        src = cons.ProteomeSource(taxid=1, binomial="A b", fasta_path=src_fa, phase="1a")
        out_dir = tmp_path / "out"

        s1 = cons.consolidate_one(src, out_dir=out_dir)
        assert s1.status == "ok"

        # Corrupt the target then force-overwrite
        target = out_dir / "1_A_b.fa"
        target.write_text("corrupted\n")

        s2 = cons.consolidate_one(src, out_dir=out_dir, force=True)
        assert s2.status == "ok"
        # Rewritten file should have proper header prefix
        assert target.read_text().startswith(">1_A_b|a")


# ---------------------------------------------------------------------------
# (e) duplicate-taxid handling
# ---------------------------------------------------------------------------

class TestDuplicateTaxid:
    """load_sources tracks duplicate taxids and returns them as duplicate_taxid
    status entries. First-seen taxid wins (early phases processed first).
    """

    def test_existing_target_without_force_returns_ok_skip(self, tmp_path: Path) -> None:
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        target = out_dir / "1_A_b.fa"
        target.write_text(">1_A_b|existing\nM\n")

        src_fa = tmp_path / "in.faa"
        src_fa.write_text(">x\nM\n")
        src = cons.ProteomeSource(taxid=1, binomial="A b", fasta_path=src_fa, phase="1f")
        status = cons.consolidate_one(src, out_dir=out_dir, force=False)
        # Idempotent skip — target already populated
        assert status.status == "ok"
        # File content unchanged (not overwritten)
        assert "existing" in target.read_text()

    def test_duplicate_taxid_emitted_by_load_sources(self, tmp_path: Path) -> None:
        """When same taxid appears in multiple manifests, second copy is marked
        as duplicate_taxid in the returned list of dropped duplicates."""
        phase1a_m = tmp_path / "proteome_manifest.tsv"
        phase1e_m = tmp_path / "genome_inventory_unannotated.tsv"
        proteome_dir = tmp_path / "cache" / "proteomes"
        proteome_dir.mkdir(parents=True)
        braker4_dir = tmp_path / "braker4_output"

        # Create proteome files
        (proteome_dir / "100_Shared_taxid.faa").write_text(">p1\nM\n")
        (braker4_dir / "100_Shared_taxid" / "results").mkdir(parents=True)
        (braker4_dir / "100_Shared_taxid" / "results" / "braker.aa").write_text(">p2\nM\n")

        # Phase 1a has taxid 100
        _write_tsv(phase1a_m, _PHASE1A_COLS, [
            {"taxid": "100", "binomial": "Shared taxid",
             "accession": "GCF_000001.1"},
        ])
        # Phase 1e ALSO has taxid 100 (duplicate)
        _write_tsv(phase1e_m, _PHASE1E_COLS, [
            {"taxid": "100", "binomial": "Shared taxid",
             "accession": "GCA_000002.1"},
        ])

        import argparse
        args = argparse.Namespace(
            phase1a_manifest=phase1a_m,
            phase1d_manifest=None,
            phase1e_manifest=phase1e_m,
            phase1g_manifest=None,
            berghia_proteome=None,
            base_dir=tmp_path,
            braker4_output_dir=braker4_dir,
        )
        sources, duplicates = cons.load_sources(args)

        # Should have 1 source (phase 1a wins, first-seen)
        assert len(sources) == 1
        assert sources[0].taxid == 100
        assert sources[0].phase == "1a"

        # Should have 1 duplicate marked as duplicate_taxid
        assert len(duplicates) == 1
        dup = duplicates[0]
        assert dup.taxid == 100
        assert dup.status == "duplicate_taxid"
        assert dup.phase == "1f"  # phase 1e uses "1f" label
        assert "already seen in phase 1a" in dup.message

    def test_write_consolidation_report_includes_all_statuses(
        self, tmp_path: Path
    ) -> None:
        statuses = [
            cons.ConsolidationStatus(
                taxid=1, binomial="A b", phase="1a", status="ok",
                n_seqs=5, target=tmp_path / "1_A_b.fa",
            ),
            cons.ConsolidationStatus(
                taxid=2, binomial="C d", phase="1f", status="missing_proteome",
                message="not found",
            ),
            cons.ConsolidationStatus(
                taxid=3, binomial="E f", phase="1f", status="duplicate_taxid",
                message="already seen in phase 1a",
            ),
        ]
        report_path = tmp_path / "report.tsv"
        cons.write_consolidation_report(statuses, report_path)

        with report_path.open() as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)

        assert len(rows) == 3
        assert rows[0]["status"] == "ok"
        assert rows[0]["n_seqs"] == "5"
        assert rows[1]["status"] == "missing_proteome"
        assert rows[1]["taxid"] == "2"
        assert rows[2]["status"] == "duplicate_taxid"
        assert rows[2]["taxid"] == "3"


# ---------------------------------------------------------------------------
# Integration: load_sources reads all three manifest schemas
# ---------------------------------------------------------------------------

# Manifest column sets for fixture helpers
_PHASE1A_COLS = (
    "taxid", "binomial", "clade", "source", "accession",
    "assembly_level", "annotation_status", "est_protein_count",
    "submission_date", "drop_reason",
)
_PHASE1D_COLS = (
    "taxid", "binomial", "policy_class", "clade_name", "source",
    "accession", "assembly_level", "annotation_status",
    "est_protein_count", "submission_date", "contig_n50", "total_length_bp",
)
_PHASE1E_COLS = (
    "taxid", "binomial", "clade", "source", "accession",
    "assembly_level", "annotation_status", "est_protein_count",
    "submission_date", "drop_reason", "contig_n50", "total_length_bp",
)


def _write_tsv(path: Path, cols: tuple, rows: list[dict]) -> None:
    with path.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(r.get(c, "") for c in cols) + "\n")


class TestLoadSources:
    """load_sources must handle all three manifest schemas + Berghia arg."""

    def test_phase1a_with_accession_yields_source(self, tmp_path: Path) -> None:
        m = tmp_path / "proteome_manifest.tsv"
        proteome_dir = tmp_path / "cache" / "proteomes"
        proteome_dir.mkdir(parents=True)
        (proteome_dir / "100_Testus_species.faa").write_text(">p1\nM\n")
        _write_tsv(m, _PHASE1A_COLS, [
            {"taxid": "100", "binomial": "Testus species",
             "accession": "GCF_000001.1"},
        ])
        import argparse
        args = argparse.Namespace(
            phase1a_manifest=m,
            phase1d_manifest=None,
            phase1e_manifest=None,
            phase1g_manifest=None,
            berghia_proteome=None,
            base_dir=tmp_path,
            braker4_output_dir=tmp_path / "braker4_output",
        )
        sources, duplicates = cons.load_sources(args)
        assert len(sources) == 1
        assert len(duplicates) == 0
        assert sources[0].taxid == 100
        assert sources[0].phase == "1a"

    def test_rows_with_empty_accession_skipped(self, tmp_path: Path) -> None:
        m = tmp_path / "proteome_manifest.tsv"
        proteome_dir = tmp_path / "cache" / "proteomes"
        proteome_dir.mkdir(parents=True)
        _write_tsv(m, _PHASE1A_COLS, [
            {"taxid": "200", "binomial": "No accession sp", "accession": ""},
            {"taxid": "201", "binomial": "Has accession sp", "accession": "GCF_1"},
        ])
        import argparse
        args = argparse.Namespace(
            phase1a_manifest=m,
            phase1d_manifest=None,
            phase1e_manifest=None,
            phase1g_manifest=None,
            berghia_proteome=None,
            base_dir=tmp_path,
            braker4_output_dir=tmp_path / "braker4_output",
        )
        sources, duplicates = cons.load_sources(args)
        taxids = {s.taxid for s in sources}
        assert 200 not in taxids
        assert 201 in taxids
        assert len(duplicates) == 0

    def test_berghia_proteome_arg_creates_source(self, tmp_path: Path) -> None:
        berghia_fa = tmp_path / "berghia.proteins.fa"
        berghia_fa.write_text(">prot1\nM\n")
        import argparse
        args = argparse.Namespace(
            phase1a_manifest=None,
            phase1d_manifest=None,
            phase1e_manifest=None,
            phase1g_manifest=None,
            berghia_proteome=berghia_fa,
            base_dir=tmp_path,
            braker4_output_dir=tmp_path / "braker4_output",
        )
        sources, duplicates = cons.load_sources(args)
        assert len(sources) == 1
        assert len(duplicates) == 0
        assert sources[0].phase == "berghia"
        assert sources[0].taxid == 1287507

    def test_phase1d_clade_name_column(self, tmp_path: Path) -> None:
        m = tmp_path / "extension_inventory.tsv"
        proteome_dir = tmp_path / "cache" / "proteomes_braker4"
        proteome_dir.mkdir(parents=True)
        (proteome_dir / "300_Testus_extensus.faa").write_text(">p\nM\n")
        _write_tsv(m, _PHASE1D_COLS, [
            {"taxid": "300", "binomial": "Testus extensus",
             "policy_class": "heterobranchia",
             "clade_name": "Heterobranchia",
             "accession": "GCA_000002.1"},
        ])
        import argparse
        args = argparse.Namespace(
            phase1a_manifest=None,
            phase1d_manifest=m,
            phase1e_manifest=None,
            phase1g_manifest=None,
            berghia_proteome=None,
            base_dir=tmp_path,
            braker4_output_dir=tmp_path / "braker4_output",
        )
        sources, duplicates = cons.load_sources(args)
        assert len(sources) == 1
        assert sources[0].taxid == 300
        assert sources[0].phase == "1f"
        assert len(duplicates) == 0

    def test_phase1f_braker_aa_fallback_path_used(self, tmp_path: Path) -> None:
        """When canonical proteome_braker4 cache doesn't exist, _locate_proteome
        returns the raw braker.aa fallback path. consolidate_one successfully
        processes the fallback."""
        m = tmp_path / "genome_inventory_unannotated.tsv"
        braker4_output = tmp_path / "braker4_output"

        # Do NOT create the canonical cache path
        # Instead, create the raw BRAKER4 output structure
        leaf = "400_Fallback_species"
        raw_braker_dir = braker4_output / leaf / "results"
        raw_braker_dir.mkdir(parents=True)
        (raw_braker_dir / "braker.aa").write_text(">prot1\nMSEQ\n>prot2\nMOTH\n")

        _write_tsv(m, _PHASE1E_COLS, [
            {"taxid": "400", "binomial": "Fallback species",
             "accession": "GCA_000003.1"},
        ])

        import argparse
        args = argparse.Namespace(
            phase1a_manifest=None,
            phase1d_manifest=None,
            phase1e_manifest=m,
            phase1g_manifest=None,
            berghia_proteome=None,
            base_dir=tmp_path,
            braker4_output_dir=braker4_output,
        )
        sources, duplicates = cons.load_sources(args)

        assert len(sources) == 1
        assert len(duplicates) == 0

        # Verify that the source's fasta_path points to the raw braker.aa
        assert sources[0].fasta_path.name == "braker.aa"
        assert "400_Fallback_species" in str(sources[0].fasta_path)

        # Now consolidate and verify it worked
        out_dir = tmp_path / "orthofinder_input"
        status = cons.consolidate_one(sources[0], out_dir=out_dir)

        assert status.status == "ok"
        assert status.n_seqs == 2
        target = out_dir / f"{leaf}.fa"
        assert target.exists()
        content = target.read_text()
        assert ">400_Fallback_species|prot1" in content
        assert "MSEQ" in content
