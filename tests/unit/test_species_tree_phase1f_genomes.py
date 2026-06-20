"""Tests for scripts/download_species_tree_phase1f_genomes.py.

Bead -p49 (Phase 1f BRAKER4, under -dnk species-tree epic). Companion to
the -9dn protein+CDS downloader; this one fetches only the genome FASTA
(`datasets download genome accession <ACC> --include genome`) for the
134 species in `references/species_tree/genome_inventory_unannotated.tsv`
that have a public GenBank assembly but no public annotation.

Output: `<taxid>_<binomial>.fasta` per species. BRAKER4 will annotate them
in Phase 1f.

The 'gotcha' classes from -9dn that this re-tests:
  - main() callable as a subprocess without NameError
  - paired-status / idempotency end-state semantics
  - failed accessions don't abort the batch
"""
from __future__ import annotations

from pathlib import Path

import pytest

import download_species_tree_phase1f_genomes as dl


# ----------------------------------------------------------------------
# read_download_targets — Phase 1e manifest schema
# ----------------------------------------------------------------------

class TestReadDownloadTargets:
    """Phase 1e manifest has the standard Phase-1a columns PLUS two extras
    (contig_n50, total_length_bp). The downloader only cares about the
    standard columns; the extras get ignored.
    """

    PHASE1E_COLUMNS = (
        "taxid", "binomial", "clade", "source", "accession",
        "assembly_level", "annotation_status", "est_protein_count",
        "submission_date", "drop_reason", "contig_n50", "total_length_bp",
    )

    def _write_manifest(self, path: Path, rows: list[dict]) -> None:
        with path.open("w") as f:
            f.write("\t".join(self.PHASE1E_COLUMNS) + "\n")
            for r in rows:
                f.write("\t".join(r.get(c, "") for c in self.PHASE1E_COLUMNS) + "\n")

    def test_returns_only_rows_with_accession(self, tmp_path: Path) -> None:
        manifest = tmp_path / "genome_inventory_unannotated.tsv"
        self._write_manifest(manifest, [
            {
                "taxid": "6161", "binomial": "Dugesia japonica",
                "clade": "platyhelminthes", "source": "GenBank",
                "accession": "GCA_001938525.1",
                "assembly_level": "Scaffold",
                "drop_reason": "",
                "contig_n50": "1166",
                "total_length_bp": "854187087",
            },
            {
                "taxid": "9999", "binomial": "Dropped example",
                "clade": "x", "accession": "",
                "drop_reason": "no_genome_in_ncbi",
            },
        ])

        targets = dl.read_download_targets(manifest)

        assert len(targets) == 1
        assert targets[0].taxid == 6161
        assert targets[0].binomial == "Dugesia japonica"
        assert targets[0].clade == "platyhelminthes"
        assert targets[0].accession == "GCA_001938525.1"

    def test_extra_columns_are_ignored(self, tmp_path: Path) -> None:
        # The contig_n50 + total_length_bp columns shouldn't break the parser.
        manifest = tmp_path / "m.tsv"
        self._write_manifest(manifest, [
            {"taxid": "100", "binomial": "Aaa example", "clade": "x",
             "accession": "GCA_100.1", "contig_n50": "999", "total_length_bp": "1234567"},
        ])
        targets = dl.read_download_targets(manifest)
        assert len(targets) == 1
        assert targets[0].accession == "GCA_100.1"

    def test_sorted_by_taxid(self, tmp_path: Path) -> None:
        manifest = tmp_path / "m.tsv"
        self._write_manifest(manifest, [
            {"taxid": "9999", "binomial": "Zzz", "clade": "x", "accession": "GCA_9.1"},
            {"taxid": "100",  "binomial": "Aaa", "clade": "x", "accession": "GCA_1.1"},
            {"taxid": "5000", "binomial": "Mmm", "clade": "x", "accession": "GCA_5.1"},
        ])
        targets = dl.read_download_targets(manifest)
        assert [t.taxid for t in targets] == [100, 5000, 9999]


# ----------------------------------------------------------------------
# read_download_targets — Phase 1d extension support + dedup
# (Bead -vqh: download script must also consume extension_inventory.tsv)
# ----------------------------------------------------------------------

class TestReadDownloadTargetsExtensionManifest:

    PHASE1E_COLUMNS = (
        "taxid", "binomial", "clade", "source", "accession",
        "assembly_level", "annotation_status", "est_protein_count",
        "submission_date", "drop_reason", "contig_n50", "total_length_bp",
    )
    PHASE1D_COLUMNS = (
        "taxid", "binomial", "policy_class", "clade_name", "source",
        "accession", "assembly_level", "annotation_status",
        "est_protein_count", "submission_date", "contig_n50",
        "total_length_bp",
    )

    def _write_p1e(self, path: Path, rows: list[dict]) -> None:
        with path.open("w") as f:
            f.write("\t".join(self.PHASE1E_COLUMNS) + "\n")
            for r in rows:
                f.write("\t".join(r.get(c, "") for c in self.PHASE1E_COLUMNS) + "\n")

    def _write_p1d(self, path: Path, rows: list[dict]) -> None:
        with path.open("w") as f:
            f.write("\t".join(self.PHASE1D_COLUMNS) + "\n")
            for r in rows:
                f.write("\t".join(r.get(c, "") for c in self.PHASE1D_COLUMNS) + "\n")

    def test_phase1d_clade_name_column_recognized(self, tmp_path: Path) -> None:
        m = tmp_path / "p1d.tsv"
        self._write_p1d(m, [
            {"taxid": "20", "binomial": "Ccc ddd",
             "policy_class": "heterobranchia",
             "clade_name": "Heterobranchia", "accession": "GCA_D2"},
        ])
        targets = dl.read_download_targets(m)
        assert len(targets) == 1
        assert targets[0].clade == "Heterobranchia"
        assert targets[0].accession == "GCA_D2"

    def test_combined_dedup_first_wins(self, tmp_path: Path) -> None:
        m1 = tmp_path / "p1e.tsv"
        m2 = tmp_path / "p1d.tsv"
        self._write_p1e(m1, [
            {"taxid": "10", "binomial": "Phase1e", "clade": "p1e_clade",
             "accession": "GCA_E1"},
            {"taxid": "20", "binomial": "Other", "clade": "p1e_clade",
             "accession": "GCA_E2"},
        ])
        self._write_p1d(m2, [
            # taxid 10 overlaps Phase 1e — Phase 1e wins
            {"taxid": "10", "binomial": "Phase1d",
             "policy_class": "heterobranchia",
             "clade_name": "Heterobranchia", "accession": "GCA_D1"},
            # taxid 30 is new
            {"taxid": "30", "binomial": "Brand new",
             "policy_class": "other_mollusca",
             "clade_name": "Bivalvia", "accession": "GCA_D3"},
        ])
        targets = dl.read_download_targets(m1, m2)
        by_taxid = {t.taxid: t for t in targets}
        assert set(by_taxid) == {10, 20, 30}
        # Phase 1e wins dedup for taxid 10
        assert by_taxid[10].accession == "GCA_E1"
        assert by_taxid[10].clade == "p1e_clade"
        # New Phase 1d entry came through
        assert by_taxid[30].clade == "Bivalvia"
        assert by_taxid[30].accession == "GCA_D3"

    def test_single_phase1d_only(self, tmp_path: Path) -> None:
        m = tmp_path / "only_p1d.tsv"
        self._write_p1d(m, [
            {"taxid": "100", "binomial": "Aa bb",
             "policy_class": "rare_basal",
             "clade_name": "Polyplacophora", "accession": "GCA_X"},
        ])
        targets = dl.read_download_targets(m)
        assert len(targets) == 1
        assert targets[0].clade == "Polyplacophora"


# ----------------------------------------------------------------------
# target_output_path — single FASTA, not paired
# ----------------------------------------------------------------------

class TestTargetOutputPath:
    """Phase 1f genomes are a single FASTA per species (.fasta), not paired
    like -9dn's protein + CDS. The naming follows the same
    `<taxid>_<binomial>` scheme with spaces -> underscores.
    """

    def test_single_path(self, tmp_path: Path) -> None:
        target = dl.DownloadTarget(
            taxid=6161, binomial="Dugesia japonica",
            clade="platyhelminthes", accession="GCA_001938525.1",
        )
        path = dl.target_output_path(target, tmp_path)
        assert path == tmp_path / "6161_Dugesia_japonica.fasta"

    def test_subspecies_name(self, tmp_path: Path) -> None:
        target = dl.DownloadTarget(
            taxid=12345, binomial="Apis mellifera carnica",
            clade="other", accession="GCA_x.1",
        )
        path = dl.target_output_path(target, tmp_path)
        assert path.name == "12345_Apis_mellifera_carnica.fasta"


# ----------------------------------------------------------------------
# is_already_downloaded — single-file end state
# ----------------------------------------------------------------------

class TestIsAlreadyDownloaded:
    """Single end state: the .fasta file exists and is non-empty."""

    def test_file_present_and_non_empty(self, tmp_path: Path) -> None:
        target = dl.DownloadTarget(
            taxid=1, binomial="Aaa bbb", clade="x", accession="GCA_1.1",
        )
        path = dl.target_output_path(target, tmp_path)
        path.write_text(">chr1\nACGT\n")
        assert dl.is_already_downloaded(target, tmp_path) is True

    def test_file_present_but_empty(self, tmp_path: Path) -> None:
        target = dl.DownloadTarget(
            taxid=1, binomial="Aaa bbb", clade="x", accession="GCA_1.1",
        )
        path = dl.target_output_path(target, tmp_path)
        path.write_text("")
        assert dl.is_already_downloaded(target, tmp_path) is False

    def test_file_missing(self, tmp_path: Path) -> None:
        target = dl.DownloadTarget(
            taxid=1, binomial="Aaa bbb", clade="x", accession="GCA_1.1",
        )
        assert dl.is_already_downloaded(target, tmp_path) is False


# ----------------------------------------------------------------------
# build_datasets_download_argv — --include genome
# ----------------------------------------------------------------------

class TestBuildDatasetsDownloadArgv:
    def test_includes_genome_only(self) -> None:
        argv = dl.build_datasets_download_argv(
            accession="GCA_001938525.1",
            output_zip=Path("/tmp/out.zip"),
            datasets_bin="datasets",
        )
        assert argv[0] == "datasets"
        assert "download" in argv and "genome" in argv and "accession" in argv
        assert "GCA_001938525.1" in argv
        include_idx = argv.index("--include")
        # Critical: GENOME only, not protein/cds.
        assert argv[include_idx + 1] == "genome"
        fn_idx = argv.index("--filename")
        assert argv[fn_idx + 1] == "/tmp/out.zip"


# ----------------------------------------------------------------------
# find_extracted_genome — different glob than protein.faa
# ----------------------------------------------------------------------

class TestFindExtractedGenome:
    """NCBI Datasets unzip layout for genome-only: the FASTA is named
    `<accession>_<asm_name>_genomic.fna` inside
    `ncbi_dataset/data/<accession>/`. The `cds_from_genomic.fna` file
    is NOT shipped because we requested only --include genome.
    """

    def _make_ncbi_layout(self, root: Path, accession: str,
                          with_genome: bool = True) -> Path:
        data_dir = root / "ncbi_dataset" / "data" / accession
        data_dir.mkdir(parents=True)
        if with_genome:
            (data_dir / f"{accession}_ASM193852v1_genomic.fna").write_text(
                ">chr1\nACGT\n",
            )
        return data_dir

    def test_finds_genome_fasta(self, tmp_path: Path) -> None:
        self._make_ncbi_layout(tmp_path, "GCA_001938525.1")
        path = dl.find_extracted_genome(tmp_path, "GCA_001938525.1")
        assert path is not None
        assert path.name.endswith("_genomic.fna")
        assert path.name.startswith("GCA_001938525.1_")

    def test_missing_genome_returns_none(self, tmp_path: Path) -> None:
        self._make_ncbi_layout(tmp_path, "GCA_xx.1", with_genome=False)
        path = dl.find_extracted_genome(tmp_path, "GCA_xx.1")
        assert path is None

    def test_wrong_accession_returns_none(self, tmp_path: Path) -> None:
        self._make_ncbi_layout(tmp_path, "GCA_AAA.1")
        path = dl.find_extracted_genome(tmp_path, "GCA_BBB.1")
        assert path is None

    def test_does_not_match_cds_from_genomic(self, tmp_path: Path) -> None:
        # Defensive: if a stray `cds_from_genomic.fna` ended up in the
        # directory (it shouldn't, since we requested only --include genome,
        # but a previous run could have polluted the workdir), we MUST NOT
        # return that as the genome — the pattern `GC[AF]_*_genomic.fna`
        # explicitly excludes the `cds_from_genomic.fna` form.
        data_dir = tmp_path / "ncbi_dataset" / "data" / "GCA_xx.1"
        data_dir.mkdir(parents=True)
        (data_dir / "cds_from_genomic.fna").write_text(">cds1\nATG\n")
        path = dl.find_extracted_genome(tmp_path, "GCA_xx.1")
        assert path is None


# ----------------------------------------------------------------------
# download_one — happy path + variants
# ----------------------------------------------------------------------

def _fake_fetcher_factory(genome_content: str | None = ">chr1\nACGTACGT\n",
                          error: str = ""):
    def fetcher(accession: str, work_dir: Path) -> dl.FetchResult:
        if error:
            return dl.FetchResult(genome_fa=None, error=error, ok=False)
        work_dir.mkdir(parents=True, exist_ok=True)
        if genome_content is None:
            return dl.FetchResult(genome_fa=None, error="", ok=True)
        genome_path = work_dir / f"{accession}_ASM_genomic.fna"
        genome_path.write_text(genome_content)
        return dl.FetchResult(genome_fa=genome_path, error="", ok=True)
    return fetcher


class TestDownloadOne:
    def test_happy_path(self, tmp_path: Path) -> None:
        cache = tmp_path / "cache"
        target = dl.DownloadTarget(
            taxid=6161, binomial="Dugesia japonica",
            clade="platyhelminthes", accession="GCA_001938525.1",
        )
        result = dl.download_one(
            target, cache, fetcher=_fake_fetcher_factory(
                genome_content=">chr1\nACGT\n>chr2\nGGGG\n",
            ),
        )
        assert result.status == "ok"
        expected_path = dl.target_output_path(target, cache)
        assert result.genome_path == expected_path
        assert expected_path.exists()
        assert result.n_seqs == 2

    def test_download_failure(self, tmp_path: Path) -> None:
        cache = tmp_path / "cache"
        target = dl.DownloadTarget(
            taxid=9999, binomial="Failz example",
            clade="x", accession="GCA_BAD.1",
        )
        result = dl.download_one(
            target, cache,
            fetcher=_fake_fetcher_factory(error="connection reset"),
        )
        assert result.status == "download_failed"
        assert "connection reset" in result.error
        assert result.genome_path is None

    def test_idempotent_skip(self, tmp_path: Path) -> None:
        cache = tmp_path / "cache"
        cache.mkdir()
        target = dl.DownloadTarget(
            taxid=6161, binomial="Dugesia japonica",
            clade="platyhelminthes", accession="GCA_001938525.1",
        )
        path = dl.target_output_path(target, cache)
        path.write_text(">chr1\nACGT\n")

        def exploding_fetcher(*args, **kwargs):
            raise AssertionError("fetcher should not have been called")

        result = dl.download_one(
            target, cache, fetcher=exploding_fetcher,
        )
        assert result.status == "skipped"
        assert result.genome_path == path

    def test_missing_genome_in_archive(self, tmp_path: Path) -> None:
        cache = tmp_path / "cache"
        target = dl.DownloadTarget(
            taxid=1, binomial="Aaa bbb", clade="x", accession="GCA_x.1",
        )
        result = dl.download_one(
            target, cache,
            fetcher=_fake_fetcher_factory(genome_content=None),
        )
        assert result.status == "download_failed"
        assert "genome" in result.error.lower()


# ----------------------------------------------------------------------
# Report TSV roundtrip
# ----------------------------------------------------------------------

class TestReportRoundtrip:
    REPORT_COLUMNS = (
        "taxid", "binomial", "clade", "accession",
        "status", "n_seqs", "size_bytes", "error",
    )

    def test_columns_match(self) -> None:
        assert dl.REPORT_COLUMNS == self.REPORT_COLUMNS

    def test_write_then_read(self, tmp_path: Path) -> None:
        results = [
            dl.DownloadResult(
                target=dl.DownloadTarget(
                    taxid=6161, binomial="Dugesia japonica",
                    clade="platyhelminthes", accession="GCA_001938525.1",
                ),
                status="ok",
                genome_path=tmp_path / "6161_Dugesia_japonica.fasta",
                n_seqs=42, size_bytes=854_187_087,
            ),
            dl.DownloadResult(
                target=dl.DownloadTarget(
                    taxid=9999, binomial="Failz example",
                    clade="x", accession="GCA_BAD.1",
                ),
                status="download_failed",
                genome_path=None,
                error="connection reset",
            ),
        ]
        report = tmp_path / "report.tsv"
        dl.write_report_tsv(report, results)
        prior = dl.read_prior_status(report)
        assert prior == {
            "GCA_001938525.1": "ok",
            "GCA_BAD.1": "download_failed",
        }


# ----------------------------------------------------------------------
# download_all
# ----------------------------------------------------------------------

class TestDownloadAll:
    def test_loops_in_order_with_failure_isolation(self, tmp_path: Path) -> None:
        cache = tmp_path / "cache"
        targets = [
            dl.DownloadTarget(taxid=1, binomial="Aaa", clade="x", accession="GCA_GOOD.1"),
            dl.DownloadTarget(taxid=2, binomial="Bbb", clade="x", accession="GCA_BAD.1"),
            dl.DownloadTarget(taxid=3, binomial="Ccc", clade="x", accession="GCA_GOOD2.1"),
        ]

        def selective_fetcher(accession: str, work_dir: Path) -> dl.FetchResult:
            if "BAD" in accession:
                return dl.FetchResult(genome_fa=None, error="simulated", ok=False)
            return _fake_fetcher_factory()(accession, work_dir)

        results = dl.download_all(
            targets, cache, fetcher=selective_fetcher,
        )
        assert [r.status for r in results] == ["ok", "download_failed", "ok"]
        assert [r.target.taxid for r in results] == [1, 2, 3]


# ----------------------------------------------------------------------
# make_ncbi_fetcher with mocked subprocess
# ----------------------------------------------------------------------

class TestMakeNcbiFetcher:
    def _make_zip(self, zip_path: Path, accession: str,
                  with_genome: bool = True) -> None:
        import zipfile
        zip_path.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zip_path, "w") as z:
            if with_genome:
                z.writestr(
                    f"ncbi_dataset/data/{accession}/{accession}_ASM1_genomic.fna",
                    ">chr1\nACGT\n",
                )

    def test_success(self, tmp_path: Path) -> None:
        import subprocess
        accession = "GCA_001938525.1"

        def fake_runner(argv, **kwargs):
            include_idx = argv.index("--include")
            assert argv[include_idx + 1] == "genome"
            zip_path = Path(argv[argv.index("--filename") + 1])
            self._make_zip(zip_path, accession)
            return subprocess.CompletedProcess(argv, 0, "", "")

        fetcher = dl.make_ncbi_fetcher("datasets", runner=fake_runner)
        result = fetcher(accession, tmp_path / "work")
        assert result.ok
        assert result.genome_fa is not None and result.genome_fa.exists()

    def test_subprocess_failure(self, tmp_path: Path) -> None:
        import subprocess

        def fake_runner(argv, **kwargs):
            return subprocess.CompletedProcess(argv, 1, "", "no such accession\n")

        fetcher = dl.make_ncbi_fetcher("datasets", runner=fake_runner)
        result = fetcher("GCA_BAD.1", tmp_path / "work")
        assert not result.ok
        assert "no such accession" in result.error


# ----------------------------------------------------------------------
# CLI + subprocess smoke test (catches name-ordering bugs)
# ----------------------------------------------------------------------

class TestArgparser:
    def test_required_args(self) -> None:
        parser = dl._build_argparser()
        args = parser.parse_args([
            "--manifest", "m.tsv", "--cache-dir", "out",
        ])
        # nargs='+' wraps the value in a list
        assert args.manifest == [Path("m.tsv")]
        assert args.cache_dir == Path("out")
        assert args.datasets_bin == "datasets"

    def test_default_manifest_is_unified_genome_inventory(self) -> None:
        parser = dl._build_argparser()
        args = parser.parse_args(["--cache-dir", "out"])
        # Default should point at the unified manifest, not the old split files.
        assert args.manifest == [Path("references/species_tree/genome_inventory.tsv")]


# ----------------------------------------------------------------------
# Unified genome_inventory.tsv schema (14-col superset)
# ----------------------------------------------------------------------

def test_reads_single_unified_manifest(tmp_path: Path) -> None:
    m = tmp_path / "genome_inventory.tsv"
    m.write_text(
        "taxid\tbinomial\tclade\tpolicy_class\tsource\taccession\t"
        "assembly_level\tannotation_status\test_protein_count\tsubmission_date\t"
        "contig_n50\ttotal_length_bp\tdrop_reason\tsource_batch\n"
        "100\tFoo bar\tgastropoda\t\tGenBank\tGCA_100\tScaffold\t\t0\t2020\t5000\t1000000\t\tnath_phase1e\n"
    )
    targets = dl.read_download_targets(m)
    assert [(t.taxid, t.accession) for t in targets] == [(100, "GCA_100")]


class TestSubprocessSmoke:
    def test_script_runs_without_name_error_on_empty_manifest(
        self, tmp_path: Path,
    ) -> None:
        """Repeats the -9dn lesson: pytest imports the module, missing the
        function-ordering bug. Only an actual subprocess invocation reproduces
        the runtime issue.
        """
        import subprocess
        columns = (
            "taxid", "binomial", "clade", "source", "accession",
            "assembly_level", "annotation_status", "est_protein_count",
            "submission_date", "drop_reason", "contig_n50", "total_length_bp",
        )
        manifest = tmp_path / "m.tsv"
        manifest.write_text("\t".join(columns) + "\n")
        script = (
            Path(__file__).resolve().parent.parent.parent
            / "scripts" / "download_species_tree_phase1f_genomes.py"
        )
        result = subprocess.run(
            [
                "python", str(script),
                "--manifest", str(manifest),
                "--cache-dir", str(tmp_path / "cache"),
            ],
            capture_output=True, text=True, timeout=30,
        )
        assert "NameError" not in result.stderr, result.stderr
        assert "ImportError" not in result.stderr, result.stderr
        assert result.returncode == 0, (
            f"rc={result.returncode}, stderr={result.stderr!r}"
        )
