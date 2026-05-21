"""Tests for scripts/build_species_tree_phase1a_inventory.py.

Bead -c5d (Phase 1a of -dnk epic): proteome inventory only, no
downloads. Queries NCBI Datasets per taxid to build a download
manifest listing the best annotated assembly per species, so we can
see expected coverage (120-180/239 estimated) before committing to
the bulk-download step (Phase 1b).
"""
from __future__ import annotations

from pathlib import Path

import pytest

import build_species_tree_phase1a_inventory as inv


# ----------------------------------------------------------------------
# parse_reference_filename
# ----------------------------------------------------------------------

class TestParseReferenceFilename:
    """Filename schema: `<taxid>_<binomial>.faa`, where binomial uses
    underscores for spaces (e.g. `100452_Candidula_unifasciata.faa`).
    Real-world taxa we have to handle:
      - simple: `100452_Candidula_unifasciata.faa`
      - subspecies: `9606_Homo_sapiens.faa` (only two name tokens, but our refs include this style)
      - multi-word epithet: `216498_Chrysomallon_squamiferum.faa`
    Edge cases:
      - missing taxid prefix should raise
      - non-integer taxid should raise
      - missing .faa extension should raise
    """

    def test_simple_binomial(self) -> None:
        taxid, binomial = inv.parse_reference_filename("100452_Candidula_unifasciata.faa")
        assert taxid == 100452
        assert binomial == "Candidula unifasciata"

    def test_multi_token_genus(self) -> None:
        taxid, binomial = inv.parse_reference_filename("216498_Chrysomallon_squamiferum.faa")
        assert taxid == 216498
        assert binomial == "Chrysomallon squamiferum"

    def test_path_input_works(self) -> None:
        # Function should accept Path or str
        taxid, binomial = inv.parse_reference_filename(
            Path("/some/path/100452_Candidula_unifasciata.faa")
        )
        assert taxid == 100452
        assert binomial == "Candidula unifasciata"

    def test_strips_directory_and_extension(self) -> None:
        taxid, binomial = inv.parse_reference_filename(
            "references/nath_et_al/one_to_one_ortholog/gastropoda/6500_Aplysia_californica.faa"
        )
        assert taxid == 6500
        assert binomial == "Aplysia californica"

    def test_missing_taxid_raises(self) -> None:
        with pytest.raises(ValueError, match="taxid"):
            inv.parse_reference_filename("Aplysia_californica.faa")

    def test_non_integer_taxid_raises(self) -> None:
        with pytest.raises(ValueError, match="taxid"):
            inv.parse_reference_filename("ABC_Aplysia_californica.faa")

    def test_missing_extension_raises(self) -> None:
        with pytest.raises(ValueError, match="extension"):
            inv.parse_reference_filename("6500_Aplysia_californica")


# ----------------------------------------------------------------------
# collect_reference_taxa
# ----------------------------------------------------------------------

class TestCollectReferenceTaxa:
    """Scans references/nath_et_al/one_to_one_ortholog/<clade>/ for
    .faa files. Returns list of ReferenceTaxon (taxid, binomial, clade)
    sorted by taxid for determinism.
    """

    def _make_refs_tree(self, root: Path, layout: dict) -> Path:
        """Build a fake references dir under `root`, return it.
        `layout` is {clade: [filenames]}.
        """
        refs = root / "one_to_one_ortholog"
        for clade, files in layout.items():
            d = refs / clade
            d.mkdir(parents=True, exist_ok=True)
            for fn in files:
                (d / fn).write_text(">stub\nACDEFG\n")
        return refs

    def test_simple_two_clades(self, tmp_path: Path) -> None:
        refs = self._make_refs_tree(
            tmp_path,
            {
                "gastropoda": [
                    "6500_Aplysia_californica.faa",
                    "225164_Lottia_gigantea.faa",
                ],
                "bivalvia": ["29159_Crassostrea_gigas.faa"],
            },
        )
        taxa = inv.collect_reference_taxa(refs)
        assert len(taxa) == 3
        # Sorted by taxid: 6500, 29159, 225164
        taxids = [t.taxid for t in taxa]
        assert taxids == sorted(taxids)
        # Each entry has the right shape
        a = next(t for t in taxa if t.taxid == 6500)
        assert a.binomial == "Aplysia californica"
        assert a.clade == "gastropoda"

    def test_empty_dir_returns_empty_list(self, tmp_path: Path) -> None:
        refs = self._make_refs_tree(tmp_path, {"gastropoda": []})
        taxa = inv.collect_reference_taxa(refs)
        assert taxa == []

    def test_skips_non_faa_files(self, tmp_path: Path) -> None:
        refs = self._make_refs_tree(
            tmp_path,
            {"gastropoda": ["6500_Aplysia_californica.faa", "README.md", "junk.txt"]},
        )
        taxa = inv.collect_reference_taxa(refs)
        assert len(taxa) == 1
        assert taxa[0].taxid == 6500

    def test_missing_refs_dir_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            inv.collect_reference_taxa(tmp_path / "does_not_exist")

    def test_unparseable_filename_skipped_with_warning(
        self, tmp_path: Path, capsys: pytest.CaptureFixture
    ) -> None:
        # Filename without taxid prefix: skip + emit warning on stderr,
        # don't crash. Important because real refs may contain odd files
        # (test fixtures, README leftovers).
        refs = self._make_refs_tree(
            tmp_path,
            {
                "gastropoda": [
                    "6500_Aplysia_californica.faa",
                    "weird_no_taxid.faa",
                ]
            },
        )
        taxa = inv.collect_reference_taxa(refs)
        assert len(taxa) == 1
        assert taxa[0].taxid == 6500
        captured = capsys.readouterr()
        assert "weird_no_taxid.faa" in captured.err


# ----------------------------------------------------------------------
# pick_best_assembly
# ----------------------------------------------------------------------

class TestPickBestAssembly:
    """Selection priority from NCBI Datasets summary records:
        RefSeq + annotated  >  GenBank + annotated  >  none (drop)

    Tie-break within priority class: higher assembly level wins
    (Chromosome > Scaffold > Contig > "" / unknown), then most recent
    submission date (lex sort on submission_date string is fine for
    ISO-like YYYY-MM-DD format).

    Returns None if no candidate passes (no annotated assemblies).
    """

    @staticmethod
    def _refseq_record(acc: str, level: str = "Scaffold", date: str = "2024-01-01",
                       n_proteins: int = 20000) -> dict:
        return {
            "accession": acc,
            "source_database": "SOURCE_DATABASE_REFSEQ",
            "annotation_info": {
                "name": f"{acc}-RS_2024",
                "status": "Current",
                "stats": {"gene_counts": {"protein_coding": n_proteins}},
            },
            "assembly_info": {
                "assembly_level": level,
                "submission_date": date,
            },
        }

    @staticmethod
    def _genbank_record(acc: str, level: str = "Scaffold", date: str = "2024-01-01",
                        annotated: bool = True, n_proteins: int = 20000) -> dict:
        r = {
            "accession": acc,
            "source_database": "SOURCE_DATABASE_GENBANK",
            "assembly_info": {
                "assembly_level": level,
                "submission_date": date,
            },
        }
        if annotated:
            r["annotation_info"] = {
                "name": f"{acc}-annotation",
                "status": "Current",
                "stats": {"gene_counts": {"protein_coding": n_proteins}},
            }
        return r

    def test_refseq_annotated_beats_genbank(self) -> None:
        records = [
            self._genbank_record("GCA_000001234.1"),
            self._refseq_record("GCF_000002075.1"),
        ]
        best = inv.pick_best_assembly(records)
        assert best is not None
        assert best.accession == "GCF_000002075.1"
        assert best.source == "RefSeq"
        assert best.annotation_status == "Current"

    def test_genbank_annotated_when_no_refseq(self) -> None:
        records = [self._genbank_record("GCA_000001234.1")]
        best = inv.pick_best_assembly(records)
        assert best is not None
        assert best.accession == "GCA_000001234.1"
        assert best.source == "GenBank"

    def test_unannotated_genbank_is_dropped(self) -> None:
        records = [self._genbank_record("GCA_000001234.1", annotated=False)]
        best = inv.pick_best_assembly(records)
        assert best is None

    def test_empty_records_returns_none(self) -> None:
        assert inv.pick_best_assembly([]) is None

    def test_chromosome_level_beats_scaffold_within_refseq(self) -> None:
        records = [
            self._refseq_record("GCF_111111.1", level="Scaffold"),
            self._refseq_record("GCF_222222.1", level="Chromosome"),
        ]
        best = inv.pick_best_assembly(records)
        assert best.accession == "GCF_222222.1"

    def test_more_recent_date_breaks_tie_within_refseq_same_level(self) -> None:
        records = [
            self._refseq_record("GCF_OLD.1", level="Chromosome", date="2020-01-01"),
            self._refseq_record("GCF_NEW.1", level="Chromosome", date="2024-09-01"),
        ]
        best = inv.pick_best_assembly(records)
        assert best.accession == "GCF_NEW.1"

    def test_protein_count_carried_into_result(self) -> None:
        records = [self._refseq_record("GCF_000002075.1", n_proteins=17654)]
        best = inv.pick_best_assembly(records)
        assert best.est_protein_count == 17654

    def test_assembly_level_carried_into_result(self) -> None:
        records = [self._refseq_record("GCF_000002075.1", level="Chromosome")]
        best = inv.pick_best_assembly(records)
        assert best.assembly_level == "Chromosome"

    def test_missing_assembly_level_handled(self) -> None:
        # Real NCBI records sometimes omit assembly_level; should not crash
        rec = self._refseq_record("GCF_000002075.1")
        del rec["assembly_info"]["assembly_level"]
        best = inv.pick_best_assembly([rec])
        assert best is not None
        assert best.assembly_level == ""

    def test_missing_protein_count_handled(self) -> None:
        rec = self._refseq_record("GCF_000002075.1")
        rec["annotation_info"]["stats"] = {}
        best = inv.pick_best_assembly([rec])
        assert best is not None
        assert best.est_protein_count == 0


# ----------------------------------------------------------------------
# Manifest TSV writer
# ----------------------------------------------------------------------

class TestWriteManifestTsv:
    """The output schema (columns, exact names):
        taxid | binomial | clade | source | accession | assembly_level |
        annotation_status | est_protein_count | submission_date | drop_reason

    Dropped rows have empty values for the assembly-derived columns
    and a populated drop_reason. Successful rows have empty drop_reason.

    Rows are written sorted by taxid (matches collect_reference_taxa).
    """

    def test_columns_in_header(self, tmp_path: Path) -> None:
        out = tmp_path / "manifest.tsv"
        inv.write_manifest_tsv(out, [])
        header = out.read_text().splitlines()[0].split("\t")
        assert header == [
            "taxid",
            "binomial",
            "clade",
            "source",
            "accession",
            "assembly_level",
            "annotation_status",
            "est_protein_count",
            "submission_date",
            "drop_reason",
        ]

    def test_successful_entry_no_drop_reason(self, tmp_path: Path) -> None:
        choice = inv.AssemblyChoice(
            accession="GCF_000002075.1",
            source="RefSeq",
            assembly_level="Scaffold",
            annotation_status="Current",
            est_protein_count=17654,
            submission_date="2024-01-01",
        )
        entry = inv.ManifestEntry(
            taxon=inv.ReferenceTaxon(taxid=6500, binomial="Aplysia californica", clade="gastropoda"),
            choice=choice,
            drop_reason="",
        )
        out = tmp_path / "manifest.tsv"
        inv.write_manifest_tsv(out, [entry])
        lines = out.read_text().splitlines()
        assert len(lines) == 2  # header + 1 entry
        row = lines[1].split("\t")
        assert row[0] == "6500"
        assert row[1] == "Aplysia californica"
        assert row[2] == "gastropoda"
        assert row[3] == "RefSeq"
        assert row[4] == "GCF_000002075.1"
        assert row[5] == "Scaffold"
        assert row[6] == "Current"
        assert row[7] == "17654"
        assert row[8] == "2024-01-01"
        assert row[9] == ""  # drop_reason

    def test_dropped_entry_has_empty_assembly_fields(self, tmp_path: Path) -> None:
        entry = inv.ManifestEntry(
            taxon=inv.ReferenceTaxon(
                taxid=99999, binomial="Foo barensis", clade="annelida"
            ),
            choice=None,
            drop_reason="no_proteome_in_ncbi",
        )
        out = tmp_path / "manifest.tsv"
        inv.write_manifest_tsv(out, [entry])
        row = out.read_text().splitlines()[1].split("\t")
        assert row[0] == "99999"
        assert row[3:9] == ["", "", "", "", "", ""]  # all assembly cols empty
        assert row[9] == "no_proteome_in_ncbi"

    def test_rows_sorted_by_taxid(self, tmp_path: Path) -> None:
        entries = [
            inv.ManifestEntry(
                taxon=inv.ReferenceTaxon(taxid=225164, binomial="Lottia gigantea", clade="gastropoda"),
                choice=None,
                drop_reason="no_proteome_in_ncbi",
            ),
            inv.ManifestEntry(
                taxon=inv.ReferenceTaxon(taxid=6500, binomial="Aplysia californica", clade="gastropoda"),
                choice=None,
                drop_reason="no_proteome_in_ncbi",
            ),
        ]
        out = tmp_path / "manifest.tsv"
        inv.write_manifest_tsv(out, entries)
        rows = out.read_text().splitlines()[1:]
        taxids = [int(r.split("\t")[0]) for r in rows]
        assert taxids == [6500, 225164]


# ----------------------------------------------------------------------
# datasets-CLI wrapper (subprocess)
# ----------------------------------------------------------------------

class TestQueryDatasetsForTaxon:
    """Tests for the thin subprocess wrapper. Real network calls happen
    only on Unity; we test the call construction + JSON-lines parsing
    with subprocess.run mocked.
    """

    def test_passes_taxid_to_datasets_cli(self, monkeypatch) -> None:
        captured = {}

        def fake_run(cmd, **kwargs):
            captured["cmd"] = cmd
            captured["kwargs"] = kwargs

            class _Result:
                returncode = 0
                stdout = ""
                stderr = ""
            return _Result()

        monkeypatch.setattr(inv.subprocess, "run", fake_run)
        inv.query_datasets_for_taxon(6500)
        assert "datasets" in captured["cmd"][0] or captured["cmd"][0] == "datasets"
        assert "summary" in captured["cmd"]
        assert "genome" in captured["cmd"]
        assert "taxon" in captured["cmd"]
        assert "6500" in captured["cmd"]
        assert "--as-json-lines" in captured["cmd"]

    def test_parses_json_lines_stdout(self, monkeypatch) -> None:
        # Two JSON records, one per line
        sample_stdout = (
            '{"accession":"GCF_000002075.1","source_database":"SOURCE_DATABASE_REFSEQ"}\n'
            '{"accession":"GCA_000111222.1","source_database":"SOURCE_DATABASE_GENBANK"}\n'
        )

        class _Result:
            returncode = 0
            stdout = sample_stdout
            stderr = ""

        monkeypatch.setattr(inv.subprocess, "run", lambda *a, **k: _Result())
        records = inv.query_datasets_for_taxon(6500)
        assert len(records) == 2
        assert records[0]["accession"] == "GCF_000002075.1"
        assert records[1]["accession"] == "GCA_000111222.1"

    def test_empty_stdout_returns_empty_list(self, monkeypatch) -> None:
        class _Result:
            returncode = 0
            stdout = ""
            stderr = ""

        monkeypatch.setattr(inv.subprocess, "run", lambda *a, **k: _Result())
        records = inv.query_datasets_for_taxon(99999)
        assert records == []

    def test_no_results_message_returns_empty_list(self, monkeypatch) -> None:
        # NCBI datasets sometimes emits a "no results" JSON to stderr;
        # stdout stays empty. Wrapper should treat this as 'no hits',
        # not an error.
        class _Result:
            returncode = 0
            stdout = ""
            stderr = "Error: No assemblies found"

        monkeypatch.setattr(inv.subprocess, "run", lambda *a, **k: _Result())
        records = inv.query_datasets_for_taxon(99999)
        assert records == []

    def test_no_genome_data_with_rc1_returns_empty(self, monkeypatch) -> None:
        # NCBI Datasets returns rc=1 with stderr like:
        #   "Error: The taxonomy ID '6426' is valid for 'Riftia pachyptila',
        #    but no genome data is currently available for this taxon."
        # This is a "valid taxon, no proteome" outcome, NOT a query error.
        # Real-world: hit on Riftia pachyptila (taxid 6426) during the
        # 2026-05-21 Phase 1a run.
        class _Result:
            returncode = 1
            stdout = ""
            stderr = (
                "Error: The taxonomy ID '6426' is valid for "
                "'Riftia pachyptila', but no genome data is currently "
                "available for this taxon."
            )

        monkeypatch.setattr(inv.subprocess, "run", lambda *a, **k: _Result())
        records = inv.query_datasets_for_taxon(6426)
        assert records == []  # should map to no_proteome_in_ncbi downstream

    def test_non_zero_returncode_raises_with_stderr(self, monkeypatch) -> None:
        class _Result:
            returncode = 2
            stdout = ""
            stderr = "datasets: command not found"

        monkeypatch.setattr(inv.subprocess, "run", lambda *a, **k: _Result())
        with pytest.raises(RuntimeError, match="datasets"):
            inv.query_datasets_for_taxon(6500)

    def test_malformed_json_line_skipped(self, monkeypatch) -> None:
        # If one JSON line is corrupt, skip it (warn) rather than abort
        # the whole batch — robust against partial NCBI failures.
        sample_stdout = (
            '{"accession":"GCF_000002075.1","source_database":"SOURCE_DATABASE_REFSEQ"}\n'
            '{ this is not valid json\n'
            '{"accession":"GCA_000111222.1","source_database":"SOURCE_DATABASE_GENBANK"}\n'
        )

        class _Result:
            returncode = 0
            stdout = sample_stdout
            stderr = ""

        monkeypatch.setattr(inv.subprocess, "run", lambda *a, **k: _Result())
        records = inv.query_datasets_for_taxon(6500)
        assert len(records) == 2
        assert {r["accession"] for r in records} == {
            "GCF_000002075.1",
            "GCA_000111222.1",
        }


# ----------------------------------------------------------------------
# build_inventory_manifest orchestrator
# ----------------------------------------------------------------------

class TestBuildInventoryManifest:
    """End-to-end orchestrator: walk refs dir, query NCBI per taxon
    (mocked), pick best assembly, write manifest TSV, return summary.

    Uses tmp_path for refs root and output TSV. NCBI queries are
    mocked through a passed-in `query_fn` to avoid network in tests.
    """

    def _make_refs(self, tmp_path: Path) -> Path:
        refs = tmp_path / "refs"
        (refs / "gastropoda").mkdir(parents=True)
        (refs / "bivalvia").mkdir(parents=True)
        (refs / "gastropoda" / "6500_Aplysia_californica.faa").write_text(">s\nM\n")
        (refs / "gastropoda" / "225164_Lottia_gigantea.faa").write_text(">s\nM\n")
        (refs / "bivalvia" / "29159_Crassostrea_gigas.faa").write_text(">s\nM\n")
        return refs

    def test_happy_path_three_taxa(self, tmp_path: Path) -> None:
        refs = self._make_refs(tmp_path)
        out = tmp_path / "manifest.tsv"

        # Mock: all three taxa have a RefSeq proteome
        def fake_query(taxid: int) -> list[dict]:
            return [{
                "accession": f"GCF_{taxid:08d}.1",
                "source_database": "SOURCE_DATABASE_REFSEQ",
                "annotation_info": {
                    "status": "Current",
                    "stats": {"gene_counts": {"protein_coding": 20000}},
                },
                "assembly_info": {
                    "assembly_level": "Scaffold",
                    "submission_date": "2024-01-01",
                },
            }]

        summary = inv.build_inventory_manifest(
            refs_root=refs, out_path=out, query_fn=fake_query,
        )
        assert summary.n_total == 3
        assert summary.n_with_proteome == 3
        assert summary.n_dropped == 0
        # All rows have a RefSeq accession
        rows = out.read_text().splitlines()[1:]
        for row in rows:
            cols = row.split("\t")
            assert cols[3] == "RefSeq"
            assert cols[4].startswith("GCF_")
            assert cols[9] == ""  # no drop reason

    def test_mixed_outcomes(self, tmp_path: Path) -> None:
        refs = self._make_refs(tmp_path)
        out = tmp_path / "manifest.tsv"

        # Mock: 6500 = RefSeq, 29159 = GenBank only, 225164 = nothing
        def fake_query(taxid: int) -> list[dict]:
            if taxid == 6500:
                return [{
                    "accession": "GCF_000002075.1",
                    "source_database": "SOURCE_DATABASE_REFSEQ",
                    "annotation_info": {"status": "Current",
                        "stats": {"gene_counts": {"protein_coding": 17654}}},
                    "assembly_info": {"assembly_level": "Scaffold",
                                       "submission_date": "2024-01-01"},
                }]
            if taxid == 29159:
                return [{
                    "accession": "GCA_999.1",
                    "source_database": "SOURCE_DATABASE_GENBANK",
                    "annotation_info": {"status": "Current",
                        "stats": {"gene_counts": {"protein_coding": 9000}}},
                    "assembly_info": {"assembly_level": "Scaffold",
                                       "submission_date": "2023-06-01"},
                }]
            return []  # taxid=225164 -> no proteome

        summary = inv.build_inventory_manifest(
            refs_root=refs, out_path=out, query_fn=fake_query,
        )
        assert summary.n_total == 3
        assert summary.n_with_proteome == 2
        assert summary.n_dropped == 1
        # Per-source breakdown carried in summary
        assert summary.n_refseq == 1
        assert summary.n_genbank == 1

    def test_writes_summary_to_stderr(
        self, tmp_path: Path, capsys: pytest.CaptureFixture
    ) -> None:
        refs = self._make_refs(tmp_path)
        out = tmp_path / "manifest.tsv"
        inv.build_inventory_manifest(
            refs_root=refs, out_path=out, query_fn=lambda taxid: [],
        )
        captured = capsys.readouterr()
        assert "3" in captured.err  # n_total mentioned
        assert "drop" in captured.err.lower() or "no_proteome" in captured.err.lower()

    def test_query_failure_records_drop_reason(self, tmp_path: Path) -> None:
        refs = self._make_refs(tmp_path)
        out = tmp_path / "manifest.tsv"

        def flaky_query(taxid: int) -> list[dict]:
            if taxid == 29159:
                raise RuntimeError("simulated NCBI 503")
            return []  # others: just no proteome

        summary = inv.build_inventory_manifest(
            refs_root=refs, out_path=out, query_fn=flaky_query,
        )
        # Query failures don't crash the whole batch — they get a
        # distinct drop_reason so we can retry the specific taxa later.
        assert summary.n_total == 3
        rows = out.read_text().splitlines()[1:]
        taxid_to_drop = {r.split("\t")[0]: r.split("\t")[9] for r in rows}
        assert taxid_to_drop["29159"] == "query_error"
        assert taxid_to_drop["6500"] == "no_proteome_in_ncbi"


# ----------------------------------------------------------------------
# Phase 1e extensions: --allow-unannotated + --limit-taxids-from
# ----------------------------------------------------------------------

class TestPickBestAssemblyUnannotated:
    """When require_annotation=False, unannotated GenBank assemblies
    become a 3rd-priority candidate: RefSeq+annot > GenBank+annot >
    GenBank+unannotated > drop. This is Phase 1e behavior — Nath et al.
    annotated unannotated GenBank assemblies themselves with BRAKER3.
    """

    @staticmethod
    def _genbank_unannotated(acc: str, level: str = "Scaffold",
                              date: str = "2024-01-01",
                              n50: int = 50000,
                              total_len: int = 500_000_000) -> dict:
        # No annotation_info key → unannotated
        return {
            "accession": acc,
            "source_database": "SOURCE_DATABASE_GENBANK",
            "assembly_info": {
                "assembly_level": level,
                "submission_date": date,
            },
            "assembly_stats": {
                "contig_n50": n50,
                "total_sequence_length": total_len,
            },
        }

    def test_strict_mode_rejects_unannotated(self) -> None:
        # require_annotation default behavior is unchanged
        records = [self._genbank_unannotated("GCA_999.1")]
        assert inv.pick_best_assembly(records) is None

    def test_allow_unannotated_picks_unannotated_when_present(self) -> None:
        records = [self._genbank_unannotated("GCA_999.1")]
        best = inv.pick_best_assembly(records, require_annotation=False)
        assert best is not None
        assert best.accession == "GCA_999.1"
        assert best.source == "GenBank"
        assert best.annotation_status == ""  # unannotated

    def test_allow_unannotated_still_prefers_annotated(self) -> None:
        # When both annotated and unannotated are present, annotated wins
        annotated = {
            "accession": "GCF_AAA.1",
            "source_database": "SOURCE_DATABASE_REFSEQ",
            "annotation_info": {"status": "Current",
                "stats": {"gene_counts": {"protein_coding": 20000}}},
            "assembly_info": {"assembly_level": "Scaffold",
                              "submission_date": "2024-01-01"},
        }
        unannotated = self._genbank_unannotated("GCA_BBB.1")
        best = inv.pick_best_assembly(
            [unannotated, annotated], require_annotation=False,
        )
        assert best.accession == "GCF_AAA.1"
        assert best.annotation_status == "Current"

    def test_assembly_stats_carried_into_choice(self) -> None:
        records = [self._genbank_unannotated(
            "GCA_999.1", n50=125000, total_len=1_500_000_000,
        )]
        best = inv.pick_best_assembly(records, require_annotation=False)
        assert best.contig_n50 == 125000
        assert best.total_length_bp == 1_500_000_000


class TestManifestTsvUnannotatedColumns:
    """Manifest TSV gains contig_n50 + total_length_bp columns when
    --allow-unannotated mode is on (Phase 1e). These help Phase 1f
    prioritize species with higher-quality assemblies for BRAKER3.
    """

    def test_extra_columns_appended(self, tmp_path: Path) -> None:
        out = tmp_path / "manifest.tsv"
        inv.write_manifest_tsv(out, [], include_assembly_stats=True)
        header = out.read_text().splitlines()[0].split("\t")
        assert "contig_n50" in header
        assert "total_length_bp" in header

    def test_default_schema_unchanged(self, tmp_path: Path) -> None:
        # Phase 1a's existing output schema (without the new columns)
        out = tmp_path / "manifest.tsv"
        inv.write_manifest_tsv(out, [])
        header = out.read_text().splitlines()[0].split("\t")
        assert header == [
            "taxid", "binomial", "clade", "source", "accession",
            "assembly_level", "annotation_status", "est_protein_count",
            "submission_date", "drop_reason",
        ]


class TestLimitTaxidsFromManifest:
    """When provided a prior Phase 1a manifest, build_inventory_manifest
    can restrict to dropped taxa only — supports Phase 1e re-query of
    only the 148 species that had no proteome in Phase 1a, without
    re-querying the 91 successful species.
    """

    def test_dropped_subset_picked(self, tmp_path: Path) -> None:
        # Write a fake prior manifest with one success + one drop
        prior = tmp_path / "prior.tsv"
        cols = ("taxid","binomial","clade","source","accession",
                "assembly_level","annotation_status","est_protein_count",
                "submission_date","drop_reason")
        with prior.open("w") as f:
            f.write("\t".join(cols) + "\n")
            f.write("6500\tAplysia californica\tgastropoda\tRefSeq\t"
                    "GCF_000002075.1\tScaffold\tCurrent\t17654\t2024-01-01\t\n")
            f.write("6359\tPlatynereis dumerilii\tannelida\t\t\t\t\t\t\t"
                    "no_proteome_in_ncbi\n")
        skip = inv.read_taxids_to_skip(prior, mode="dropped_only")
        # In dropped_only mode, only 6359 is queriable; 6500 is skipped
        assert 6500 in skip
        assert 6359 not in skip


class TestUnannotatedQueryFlag:
    """The datasets-CLI wrapper supports an --assembly-source flag.
    For Phase 1e we may need a separate query path that doesn't
    filter to annotated assemblies. Tested via query construction.
    """

    def test_include_unannotated_passes_proper_args(self, monkeypatch) -> None:
        captured: dict = {}

        class _R:
            returncode = 0
            stdout = ""
            stderr = ""

        def fake_run(cmd, **kwargs):
            captured["cmd"] = cmd
            return _R()

        monkeypatch.setattr(inv.subprocess, "run", fake_run)
        inv.query_datasets_for_taxon(6359)
        # The default query doesn't pass --reference or filter flags;
        # NCBI returns both annotated and unannotated assemblies by
        # default. Our pick_best_assembly is what filters by annotation.
        # This test just confirms the CLI invocation stays minimal so
        # Phase 1e can rely on the same wrapper.
        assert "datasets" in captured["cmd"][0] or captured["cmd"][0] == "datasets"
        assert "--reference" not in captured["cmd"]
        # We DO want to make sure no flag silently filters out
        # unannotated assemblies on the CLI side.
        assert "--annotated" not in captured["cmd"]
