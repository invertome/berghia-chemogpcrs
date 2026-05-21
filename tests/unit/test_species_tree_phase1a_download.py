"""Tests for scripts/download_species_tree_phase1a.py.

Bead -9dn (Phase 1a-CDS, under -dnk species-tree epic): download paired
protein + CDS FASTAs for the 91 species that Phase 1a inventory
identified as having a usable NCBI annotated assembly. Per-species
output: `<taxid>_<binomial>.faa` + `<taxid>_<binomial>.cds.fna` in a
single cache dir, ready for stage 05 dN/dS to consume directly (instead
of recovering CDS via miniprot).

The downloader sits on top of the canonical pattern from
`scripts/fetch_berghia_genome.sh`:

    datasets download genome accession <ACC> --include protein,cds
        --filename <ACC>.zip
    unzip <ACC>.zip
    symlink ncbi_dataset/data/<ACC>/protein.faa            -> <taxid>_<binomial>.faa
    symlink ncbi_dataset/data/<ACC>/cds_from_genomic.fna   -> <taxid>_<binomial>.cds.fna

Most edge cases come from CDS handling: GenBank-only assemblies often
ship `protein.faa` but no `cds_from_genomic.fna`, so per-accession
"no CDS in download" must be recorded (not crash the batch). Stage 05's
recover_cds remains the fallback for those.
"""
from __future__ import annotations

from pathlib import Path

import pytest

import download_species_tree_phase1a as dl


# ----------------------------------------------------------------------
# read_download_targets
# ----------------------------------------------------------------------

class TestReadDownloadTargets:
    """Read the Phase 1a inventory manifest TSV and return only the rows
    that have a non-empty accession (the species we can actually
    download). Dropped rows (no_proteome_in_ncbi, query_error) are
    skipped.
    """

    def _write_manifest(self, path: Path, rows: list[dict]) -> None:
        from build_species_tree_phase1a_inventory import MANIFEST_COLUMNS
        with path.open("w") as f:
            f.write("\t".join(MANIFEST_COLUMNS) + "\n")
            for r in rows:
                f.write("\t".join(r.get(c, "") for c in MANIFEST_COLUMNS) + "\n")

    def test_returns_only_rows_with_accession(self, tmp_path: Path) -> None:
        manifest = tmp_path / "proteome_manifest.tsv"
        self._write_manifest(manifest, [
            {
                "taxid": "6183", "binomial": "Schistosoma mansoni",
                "clade": "platyhelminthes", "source": "RefSeq",
                "accession": "GCF_000237925.1",
                "assembly_level": "Chromosome",
                "annotation_status": "Full annotation",
                "est_protein_count": "10711",
                "submission_date": "", "drop_reason": "",
            },
            {
                "taxid": "6161", "binomial": "Dugesia japonica",
                "clade": "platyhelminthes", "source": "",
                "accession": "",
                "assembly_level": "",
                "annotation_status": "",
                "est_protein_count": "",
                "submission_date": "", "drop_reason": "no_proteome_in_ncbi",
            },
        ])

        targets = dl.read_download_targets(manifest)

        assert len(targets) == 1
        assert targets[0].taxid == 6183
        assert targets[0].binomial == "Schistosoma mansoni"
        assert targets[0].clade == "platyhelminthes"
        assert targets[0].accession == "GCF_000237925.1"

    def test_sorted_by_taxid(self, tmp_path: Path) -> None:
        manifest = tmp_path / "m.tsv"
        self._write_manifest(manifest, [
            {"taxid": "9999", "binomial": "Zzz example", "clade": "x",
             "accession": "GCA_999.1", "drop_reason": ""},
            {"taxid": "100", "binomial": "Aaa example", "clade": "x",
             "accession": "GCF_100.1", "drop_reason": ""},
            {"taxid": "5000", "binomial": "Mmm example", "clade": "x",
             "accession": "GCF_5000.1", "drop_reason": ""},
        ])
        targets = dl.read_download_targets(manifest)
        assert [t.taxid for t in targets] == [100, 5000, 9999]

    def test_empty_manifest(self, tmp_path: Path) -> None:
        manifest = tmp_path / "m.tsv"
        self._write_manifest(manifest, [])
        assert dl.read_download_targets(manifest) == []


# ----------------------------------------------------------------------
# target_output_paths
# ----------------------------------------------------------------------

class TestTargetOutputPaths:
    """The paired output convention is:

        <cache_dir>/<taxid>_<Genus_species>.faa
        <cache_dir>/<taxid>_<Genus_species>.cds.fna

    matching the input filename schema used elsewhere in the
    species-tree pipeline (e.g. references/nath_et_al/one_to_one_ortholog/
    <clade>/<taxid>_<Genus_species>.faa, parsed by
    build_species_tree_phase1a_inventory.parse_reference_filename).

    Spaces in binomial are converted to underscores in the filename.
    """

    def test_paired_paths(self, tmp_path: Path) -> None:
        target = dl.DownloadTarget(
            taxid=6183, binomial="Schistosoma mansoni",
            clade="platyhelminthes", accession="GCF_000237925.1",
        )
        faa, cds = dl.target_output_paths(target, tmp_path)
        assert faa == tmp_path / "6183_Schistosoma_mansoni.faa"
        assert cds == tmp_path / "6183_Schistosoma_mansoni.cds.fna"

    def test_subspecies_with_three_tokens(self, tmp_path: Path) -> None:
        # Some refs have three-token names (e.g. "Genus species subsp"). Spaces
        # all become underscores; downstream consumers parse on the FIRST
        # underscore for taxid + everything else for binomial.
        target = dl.DownloadTarget(
            taxid=12345, binomial="Apis mellifera carnica",
            clade="other", accession="GCA_999.1",
        )
        faa, cds = dl.target_output_paths(target, tmp_path)
        assert faa.name == "12345_Apis_mellifera_carnica.faa"
        assert cds.name == "12345_Apis_mellifera_carnica.cds.fna"


# ----------------------------------------------------------------------
# is_already_downloaded
# ----------------------------------------------------------------------

class TestIsAlreadyDownloaded:
    """Idempotency check: re-running the downloader should skip species
    that already have BOTH paired outputs present and non-empty. A
    half-finished download (only .faa, no .cds.fna) should NOT be
    treated as complete — those need retrying because the protein-only
    state is indistinguishable from a CDS-less GenBank assembly until we
    actually retry.

    Distinguishing "successfully downloaded but no CDS available" from
    "interrupted mid-download" is handled by the report TSV: completed
    species without CDS get logged with status='ok_no_cds', so
    subsequent runs can use that status to skip the retry. Tested in a
    later case below.
    """

    def _touch(self, p: Path, content: str = ">seq\nATG\n") -> None:
        p.write_text(content)

    def test_both_files_present_means_done(self, tmp_path: Path) -> None:
        target = dl.DownloadTarget(
            taxid=1, binomial="Aaa bbb", clade="x", accession="GCF_1.1",
        )
        faa, cds = dl.target_output_paths(target, tmp_path)
        self._touch(faa)
        self._touch(cds)
        assert dl.is_already_downloaded(target, tmp_path, prior_status=None) is True

    def test_only_protein_means_retry(self, tmp_path: Path) -> None:
        target = dl.DownloadTarget(
            taxid=1, binomial="Aaa bbb", clade="x", accession="GCF_1.1",
        )
        faa, _cds = dl.target_output_paths(target, tmp_path)
        self._touch(faa)
        assert dl.is_already_downloaded(target, tmp_path, prior_status=None) is False

    def test_empty_files_dont_count(self, tmp_path: Path) -> None:
        target = dl.DownloadTarget(
            taxid=1, binomial="Aaa bbb", clade="x", accession="GCF_1.1",
        )
        faa, cds = dl.target_output_paths(target, tmp_path)
        faa.write_text("")
        cds.write_text("")
        assert dl.is_already_downloaded(target, tmp_path, prior_status=None) is False

    def test_protein_present_and_prior_status_ok_no_cds_means_done(self, tmp_path: Path) -> None:
        # When a prior run recorded 'ok_no_cds' for this species, having
        # just the .faa present is the legitimate end state — don't retry.
        target = dl.DownloadTarget(
            taxid=1, binomial="Aaa bbb", clade="x", accession="GCA_1.1",
        )
        faa, _cds = dl.target_output_paths(target, tmp_path)
        self._touch(faa)
        assert dl.is_already_downloaded(
            target, tmp_path, prior_status="ok_no_cds",
        ) is True

    def test_no_files_at_all(self, tmp_path: Path) -> None:
        target = dl.DownloadTarget(
            taxid=1, binomial="Aaa bbb", clade="x", accession="GCF_1.1",
        )
        assert dl.is_already_downloaded(target, tmp_path, prior_status=None) is False


# ----------------------------------------------------------------------
# download_one — happy path + variants
# ----------------------------------------------------------------------

def _fake_fetcher_factory(
    faa_content: str | None = ">p1\nMAA\n",
    cds_content: str | None = ">c1\nATGGCAGCA\n",
    error: str = "",
):
    """Build an injectable fetcher that simulates a successful NCBI Datasets
    'download + unzip' for tests. Setting faa_content=None simulates a
    download failure (RefSeq archive is missing protein.faa, which
    shouldn't happen in practice — kept only as a paranoid corner).
    Setting cds_content=None simulates a GenBank-only assembly with no
    cds_from_genomic.fna.
    """
    def fetcher(accession: str, work_dir: Path) -> dl.FetchResult:
        if error:
            return dl.FetchResult(protein_faa=None, cds_fna=None, error=error, ok=False)
        work_dir.mkdir(parents=True, exist_ok=True)
        faa_path = None
        cds_path = None
        if faa_content is not None:
            faa_path = work_dir / "protein.faa"
            faa_path.write_text(faa_content)
        if cds_content is not None:
            cds_path = work_dir / "cds_from_genomic.fna"
            cds_path.write_text(cds_content)
        return dl.FetchResult(
            protein_faa=faa_path, cds_fna=cds_path, error="", ok=True,
        )
    return fetcher


class TestDownloadOne:
    """Per-species orchestration: idempotent skip → fetch → place
    canonical paths → return a DownloadResult.
    """

    def test_happy_path_produces_paired_outputs(self, tmp_path: Path) -> None:
        cache = tmp_path / "cache"
        cache.mkdir()
        target = dl.DownloadTarget(
            taxid=6183, binomial="Schistosoma mansoni",
            clade="platyhelminthes", accession="GCF_000237925.1",
        )

        result = dl.download_one(
            target, cache,
            fetcher=_fake_fetcher_factory(
                faa_content=">p1\nMAA\n>p2\nMBB\n",
                cds_content=">c1\nATG\n>c2\nGGG\n",
            ),
        )

        assert result.status == "ok"
        assert result.error == ""
        faa_expected, cds_expected = dl.target_output_paths(target, cache)
        assert result.faa_path == faa_expected
        assert result.cds_path == cds_expected
        # Files were placed at canonical paths
        assert faa_expected.exists()
        assert cds_expected.exists()
        # FASTA sequence counts get reported
        assert result.n_proteins == 2
        assert result.n_cds == 2

    def test_no_cds_in_download_gets_ok_no_cds_status(self, tmp_path: Path) -> None:
        # GenBank-only assemblies sometimes ship protein.faa without
        # cds_from_genomic.fna. The result should still be ok-ish, but
        # status flagged so stage 05 can fall back to recover_cds.
        cache = tmp_path / "cache"
        cache.mkdir()
        target = dl.DownloadTarget(
            taxid=6182, binomial="Schistosoma japonicum",
            clade="platyhelminthes", accession="GCA_021461655.1",
        )

        result = dl.download_one(
            target, cache,
            fetcher=_fake_fetcher_factory(
                faa_content=">p1\nMAA\n",
                cds_content=None,
            ),
        )

        assert result.status == "ok_no_cds"
        assert result.error == ""
        faa_expected, cds_expected = dl.target_output_paths(target, cache)
        assert result.faa_path == faa_expected
        assert result.cds_path is None
        assert faa_expected.exists()
        assert not cds_expected.exists()  # NO cds file written
        assert result.n_proteins == 1
        assert result.n_cds == 0

    def test_download_failure_recorded(self, tmp_path: Path) -> None:
        cache = tmp_path / "cache"
        cache.mkdir()
        target = dl.DownloadTarget(
            taxid=9999, binomial="Failz example", clade="x",
            accession="GCA_BAD.1",
        )

        result = dl.download_one(
            target, cache,
            fetcher=_fake_fetcher_factory(
                error="datasets returncode=1: connection reset",
            ),
        )

        assert result.status == "download_failed"
        assert "connection reset" in result.error
        assert result.faa_path is None
        assert result.cds_path is None
        # No partial files left behind in the cache
        faa, cds = dl.target_output_paths(target, cache)
        assert not faa.exists()
        assert not cds.exists()

    def test_idempotent_skip(self, tmp_path: Path) -> None:
        cache = tmp_path / "cache"
        cache.mkdir()
        target = dl.DownloadTarget(
            taxid=6183, binomial="Schistosoma mansoni",
            clade="platyhelminthes", accession="GCF_000237925.1",
        )
        # Pre-populate the canonical paths so is_already_downloaded sees
        # them as the legitimate end state.
        faa, cds = dl.target_output_paths(target, cache)
        faa.write_text(">p1\nMAA\n")
        cds.write_text(">c1\nATG\n")

        # Use a fetcher that would explode if invoked — proves we skipped.
        def exploding_fetcher(*args, **kwargs):
            raise AssertionError("fetcher should not have been called")

        result = dl.download_one(
            target, cache, fetcher=exploding_fetcher,
        )

        assert result.status == "skipped"
        assert result.faa_path == faa
        assert result.cds_path == cds
        assert result.error == ""

    def test_idempotent_skip_honors_prior_ok_no_cds(self, tmp_path: Path) -> None:
        # GenBank species previously recorded as 'ok_no_cds' — only the
        # protein file is present, but that IS the end state.
        cache = tmp_path / "cache"
        cache.mkdir()
        target = dl.DownloadTarget(
            taxid=6182, binomial="Schistosoma japonicum",
            clade="platyhelminthes", accession="GCA_021461655.1",
        )
        faa, _cds = dl.target_output_paths(target, cache)
        faa.write_text(">p1\nMAA\n")

        def exploding_fetcher(*args, **kwargs):
            raise AssertionError("fetcher should not have been called")

        result = dl.download_one(
            target, cache,
            fetcher=exploding_fetcher,
            prior_status="ok_no_cds",
        )

        assert result.status == "skipped"
        assert result.faa_path == faa
        assert result.cds_path is None


# ----------------------------------------------------------------------
# Report TSV: write + read roundtrip
# ----------------------------------------------------------------------

class TestReportRoundtrip:
    """The per-species report TSV is how status persists across runs:
    write at end of one run, read at start of next, so 'ok_no_cds'
    species don't get retried forever.
    """

    REPORT_COLUMNS = (
        "taxid", "binomial", "clade", "accession",
        "status", "n_proteins", "n_cds", "error",
    )

    def test_columns_match(self) -> None:
        assert dl.REPORT_COLUMNS == self.REPORT_COLUMNS

    def test_write_then_read(self, tmp_path: Path) -> None:
        results = [
            dl.DownloadResult(
                target=dl.DownloadTarget(
                    taxid=6183, binomial="Schistosoma mansoni",
                    clade="platyhelminthes", accession="GCF_000237925.1",
                ),
                status="ok",
                faa_path=tmp_path / "6183_Schistosoma_mansoni.faa",
                cds_path=tmp_path / "6183_Schistosoma_mansoni.cds.fna",
                n_proteins=10711, n_cds=10711,
            ),
            dl.DownloadResult(
                target=dl.DownloadTarget(
                    taxid=6182, binomial="Schistosoma japonicum",
                    clade="platyhelminthes", accession="GCA_021461655.1",
                ),
                status="ok_no_cds",
                faa_path=tmp_path / "6182_Schistosoma_japonicum.faa",
                cds_path=None,
                n_proteins=9715, n_cds=0,
            ),
            dl.DownloadResult(
                target=dl.DownloadTarget(
                    taxid=9999, binomial="Failz example",
                    clade="x", accession="GCA_BAD.1",
                ),
                status="download_failed",
                faa_path=None, cds_path=None,
                error="connection reset",
            ),
        ]
        report = tmp_path / "report.tsv"
        dl.write_report_tsv(report, results)

        # Roundtrip: prior-status lookup returns the right map
        prior = dl.read_prior_status(report)
        assert prior == {
            "GCF_000237925.1": "ok",
            "GCA_021461655.1": "ok_no_cds",
            "GCA_BAD.1": "download_failed",
        }

    def test_read_prior_status_returns_empty_if_no_report(self, tmp_path: Path) -> None:
        nonexistent = tmp_path / "no_such_report.tsv"
        assert dl.read_prior_status(nonexistent) == {}


# ----------------------------------------------------------------------
# download_all
# ----------------------------------------------------------------------

class TestDownloadAll:
    """Top-level loop. Honors prior-status to skip ok_no_cds species
    without re-invoking fetcher, returns one DownloadResult per target,
    preserves order from the input list.
    """

    def test_loops_over_targets_in_order(self, tmp_path: Path) -> None:
        cache = tmp_path / "cache"
        targets = [
            dl.DownloadTarget(taxid=2, binomial="Bbb example",
                              clade="x", accession="GCF_2.1"),
            dl.DownloadTarget(taxid=1, binomial="Aaa example",
                              clade="x", accession="GCF_1.1"),
        ]
        # Both succeed with both files
        results = dl.download_all(
            targets, cache,
            fetcher=_fake_fetcher_factory(),
            prior_statuses=None,
        )
        assert [r.target.taxid for r in results] == [2, 1]
        assert all(r.status == "ok" for r in results)

    def test_prior_ok_no_cds_keeps_skipped(self, tmp_path: Path) -> None:
        cache = tmp_path / "cache"
        cache.mkdir()
        target = dl.DownloadTarget(
            taxid=6182, binomial="Schistosoma japonicum",
            clade="platyhelminthes", accession="GCA_021461655.1",
        )
        # Drop just the .faa to simulate a prior 'ok_no_cds' run.
        faa, _cds = dl.target_output_paths(target, cache)
        faa.write_text(">p1\nMAA\n")

        def exploding(*args, **kwargs):
            raise AssertionError("fetcher should not be called")

        results = dl.download_all(
            [target], cache,
            fetcher=exploding,
            prior_statuses={"GCA_021461655.1": "ok_no_cds"},
        )
        assert len(results) == 1
        assert results[0].status == "skipped"

    def test_failed_target_does_not_abort_batch(self, tmp_path: Path) -> None:
        cache = tmp_path / "cache"
        targets = [
            dl.DownloadTarget(taxid=1, binomial="Aaa example",
                              clade="x", accession="GCF_GOOD.1"),
            dl.DownloadTarget(taxid=2, binomial="Bbb example",
                              clade="x", accession="GCA_BAD.1"),
            dl.DownloadTarget(taxid=3, binomial="Ccc example",
                              clade="x", accession="GCF_GOOD2.1"),
        ]

        def selective_fetcher(accession: str, work_dir: Path) -> dl.FetchResult:
            if "BAD" in accession:
                return dl.FetchResult(
                    protein_faa=None, cds_fna=None,
                    error="simulated failure", ok=False,
                )
            return _fake_fetcher_factory()(accession, work_dir)

        results = dl.download_all(
            targets, cache, fetcher=selective_fetcher, prior_statuses=None,
        )
        assert [r.status for r in results] == ["ok", "download_failed", "ok"]


# ----------------------------------------------------------------------
# NCBI Datasets CLI argv + post-unzip file finder
# ----------------------------------------------------------------------

class TestBuildDatasetsDownloadArgv:
    """Builds the `datasets download genome accession ...` argv vector.
    The flag set is the contract — `--include protein,cds` is the
    central feature of -9dn (vs. inventory's `summary` calls).
    """

    def test_includes_protein_and_cds(self) -> None:
        argv = dl.build_datasets_download_argv(
            accession="GCF_000237925.1",
            output_zip=Path("/tmp/out.zip"),
            datasets_bin="datasets",
        )
        assert argv[0] == "datasets"
        assert "download" in argv and "genome" in argv and "accession" in argv
        assert "GCF_000237925.1" in argv
        # Critical: BOTH protein AND cds must be requested
        include_idx = argv.index("--include")
        assert argv[include_idx + 1] == "protein,cds"
        # The filename flag points at where the zip should land
        fn_idx = argv.index("--filename")
        assert argv[fn_idx + 1] == "/tmp/out.zip"

    def test_respects_custom_datasets_bin(self) -> None:
        argv = dl.build_datasets_download_argv(
            accession="GCA_999.1",
            output_zip=Path("/tmp/x.zip"),
            datasets_bin="/opt/datasets/bin/datasets",
        )
        assert argv[0] == "/opt/datasets/bin/datasets"


class TestFindExtractedFiles:
    """After unzipping, locate protein.faa + cds_from_genomic.fna inside
    the NCBI Datasets layout (`ncbi_dataset/data/<accession>/...`).
    Handles assemblies that ship without `cds_from_genomic.fna` (returns
    None for cds).
    """

    def _make_ncbi_layout(
        self, root: Path, accession: str,
        with_protein: bool = True, with_cds: bool = True,
    ) -> Path:
        data_dir = root / "ncbi_dataset" / "data" / accession
        data_dir.mkdir(parents=True)
        if with_protein:
            (data_dir / "protein.faa").write_text(">p1\nMAA\n")
        if with_cds:
            (data_dir / "cds_from_genomic.fna").write_text(">c1\nATG\n")
        return data_dir

    def test_finds_both_files(self, tmp_path: Path) -> None:
        self._make_ncbi_layout(tmp_path, "GCF_000237925.1")
        faa, cds = dl.find_extracted_files(tmp_path, "GCF_000237925.1")
        assert faa is not None
        assert faa.name == "protein.faa"
        assert faa.read_text() == ">p1\nMAA\n"
        assert cds is not None
        assert cds.name == "cds_from_genomic.fna"

    def test_finds_protein_only_when_no_cds_shipped(self, tmp_path: Path) -> None:
        self._make_ncbi_layout(
            tmp_path, "GCA_021461655.1", with_protein=True, with_cds=False,
        )
        faa, cds = dl.find_extracted_files(tmp_path, "GCA_021461655.1")
        assert faa is not None
        assert faa.name == "protein.faa"
        assert cds is None

    def test_missing_protein_returns_none(self, tmp_path: Path) -> None:
        self._make_ncbi_layout(
            tmp_path, "GCA_xx.1", with_protein=False, with_cds=False,
        )
        faa, cds = dl.find_extracted_files(tmp_path, "GCA_xx.1")
        assert faa is None
        assert cds is None

    def test_wrong_accession_subdir_returns_none(self, tmp_path: Path) -> None:
        # Defensive: datasets layout is keyed by accession; if the
        # caller's accession doesn't match the dir, we treat it as
        # not-found (rather than picking some other assembly's files).
        self._make_ncbi_layout(tmp_path, "GCF_AAA.1")
        faa, cds = dl.find_extracted_files(tmp_path, "GCF_BBB.1")
        assert faa is None
        assert cds is None


# ----------------------------------------------------------------------
# make_ncbi_fetcher (real fetcher, with mocked subprocess)
# ----------------------------------------------------------------------

class TestMakeNcbiFetcher:
    """End-to-end fetcher composition: argv → subprocess → unzip → find.
    We mock only the subprocess call (no real `datasets` CLI invocation),
    using a fake `runner` that drops a real zip file at the requested
    path. The fetcher's unzip step uses the real `zipfile` module so
    the layout-discovery code is exercised against actual archives.
    """

    def _make_zip(
        self, zip_path: Path, accession: str,
        with_cds: bool = True, with_protein: bool = True,
    ) -> None:
        import zipfile
        zip_path.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zip_path, "w") as z:
            if with_protein:
                z.writestr(
                    f"ncbi_dataset/data/{accession}/protein.faa",
                    ">p1\nMAA\n>p2\nMBB\n",
                )
            if with_cds:
                z.writestr(
                    f"ncbi_dataset/data/{accession}/cds_from_genomic.fna",
                    ">c1\nATG\n>c2\nGGG\n",
                )

    def test_success_with_both_files(self, tmp_path: Path) -> None:
        import subprocess
        accession = "GCF_000237925.1"

        def fake_runner(argv, **kwargs):
            # Verify the contract: protein+cds requested.
            include_idx = argv.index("--include")
            assert argv[include_idx + 1] == "protein,cds"
            zip_path = Path(argv[argv.index("--filename") + 1])
            self._make_zip(zip_path, accession)
            return subprocess.CompletedProcess(argv, 0, "", "")

        fetcher = dl.make_ncbi_fetcher("datasets", runner=fake_runner)
        result = fetcher(accession, tmp_path / "work")
        assert result.ok
        assert result.protein_faa is not None and result.protein_faa.exists()
        assert result.cds_fna is not None and result.cds_fna.exists()
        assert result.error == ""

    def test_success_with_no_cds_in_archive(self, tmp_path: Path) -> None:
        # GenBank-only assembly case: protein.faa but no cds_from_genomic.fna
        import subprocess
        accession = "GCA_021461655.1"

        def fake_runner(argv, **kwargs):
            zip_path = Path(argv[argv.index("--filename") + 1])
            self._make_zip(zip_path, accession, with_cds=False)
            return subprocess.CompletedProcess(argv, 0, "", "")

        fetcher = dl.make_ncbi_fetcher("datasets", runner=fake_runner)
        result = fetcher(accession, tmp_path / "work")
        assert result.ok
        assert result.protein_faa is not None and result.protein_faa.exists()
        assert result.cds_fna is None
        assert result.error == ""

    def test_subprocess_failure_returns_ok_false(self, tmp_path: Path) -> None:
        import subprocess
        accession = "GCA_BAD.1"

        def fake_runner(argv, **kwargs):
            return subprocess.CompletedProcess(
                argv, 1, "", "connection reset by NCBI\n",
            )

        fetcher = dl.make_ncbi_fetcher("datasets", runner=fake_runner)
        result = fetcher(accession, tmp_path / "work")
        assert not result.ok
        assert "connection reset" in result.error
        assert result.protein_faa is None
        assert result.cds_fna is None

    def test_archive_missing_protein_returns_ok_false(self, tmp_path: Path) -> None:
        # Bad archive with neither protein nor cds — treat as failure
        # (rather than silent "ok" with no files).
        import subprocess
        accession = "GCA_xx.1"

        def fake_runner(argv, **kwargs):
            zip_path = Path(argv[argv.index("--filename") + 1])
            self._make_zip(
                zip_path, accession, with_protein=False, with_cds=False,
            )
            return subprocess.CompletedProcess(argv, 0, "", "")

        fetcher = dl.make_ncbi_fetcher("datasets", runner=fake_runner)
        result = fetcher(accession, tmp_path / "work")
        assert not result.ok
        assert "protein" in result.error.lower()

    def test_unused_field_import_compat(self) -> None:
        # FetchResult uses defaults; this guards against accidental
        # signature breakage if dataclass field defaults are reordered.
        fr = dl.FetchResult(protein_faa=None, cds_fna=None)
        assert fr.ok is True
        assert fr.error == ""

    def test_clean_workdir_per_call(self, tmp_path: Path) -> None:
        # Stale work_dir contents shouldn't leak in: a stale protein.faa
        # from a prior accession must not be returned when the new
        # archive is empty.
        import subprocess
        stale_accession = "GCA_OLD.1"
        new_accession = "GCA_NEW.1"
        work_dir = tmp_path / "shared_work"
        # Pre-populate stale state under stale_accession
        stale_layout = work_dir / "ncbi_dataset" / "data" / stale_accession
        stale_layout.mkdir(parents=True)
        (stale_layout / "protein.faa").write_text(">stale\nMSTL\n")

        def fake_runner(argv, **kwargs):
            # Don't drop a real zip — simulate a failed download.
            return subprocess.CompletedProcess(
                argv, 1, "", "no such accession\n",
            )

        fetcher = dl.make_ncbi_fetcher("datasets", runner=fake_runner)
        result = fetcher(new_accession, work_dir)
        # Even though stale_layout has files, finder should return None
        # since it scopes the lookup to NEW accession.
        assert not result.ok
        assert result.protein_faa is None
        assert result.cds_fna is None


# ----------------------------------------------------------------------
# CLI argparser
# ----------------------------------------------------------------------

class TestArgparser:
    """The CLI is the entry point for the sbatch wrapper. Required args:
    --manifest (input), --cache-dir (output). Optional: --datasets-bin,
    --report.
    """

    def test_required_args(self) -> None:
        parser = dl._build_argparser()
        args = parser.parse_args([
            "--manifest", "ref/manifest.tsv",
            "--cache-dir", "scratch/proteomes_ncbi",
        ])
        assert args.manifest == Path("ref/manifest.tsv")
        assert args.cache_dir == Path("scratch/proteomes_ncbi")
        # Default datasets binary is "datasets" on PATH
        assert args.datasets_bin == "datasets"

    def test_report_defaults_relative_to_cache(self) -> None:
        # Default --report is implicit; main() resolves it relative to
        # --cache-dir as 'download_report.tsv'. The argparser leaves it
        # as None.
        parser = dl._build_argparser()
        args = parser.parse_args([
            "--manifest", "m.tsv", "--cache-dir", "out",
        ])
        assert args.report is None

    def test_custom_report_path(self) -> None:
        parser = dl._build_argparser()
        args = parser.parse_args([
            "--manifest", "m.tsv", "--cache-dir", "out",
            "--report", "logs/dl_report.tsv",
        ])
        assert args.report == Path("logs/dl_report.tsv")


# ----------------------------------------------------------------------
# main() end-to-end smoke test
# ----------------------------------------------------------------------

class TestMainSmoke:
    """Catches the kind of bug that unit tests on individual functions
    miss: name-resolution ordering, top-level imports, side effects of
    actually invoking the full pipeline. The 2026-05-21 first run hit a
    NameError because `if __name__ == '__main__'` was placed before
    function definitions it referenced.
    """

    def _write_minimal_manifest(self, path: Path) -> None:
        from build_species_tree_phase1a_inventory import MANIFEST_COLUMNS
        with path.open("w") as f:
            f.write("\t".join(MANIFEST_COLUMNS) + "\n")
            f.write("\t".join([
                "6183", "Schistosoma mansoni", "platyhelminthes",
                "RefSeq", "GCF_000237925.1", "Chromosome",
                "Full annotation", "10711", "", "",
            ]) + "\n")

    def test_main_runs_end_to_end_with_mocked_subprocess(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        # Mock subprocess.run to drop a real zip and return rc=0
        import subprocess
        import zipfile

        def fake_subprocess_run(argv, **kwargs):
            # --filename is the zip output path; --include must be protein,cds
            include_idx = argv.index("--include")
            assert argv[include_idx + 1] == "protein,cds"
            zip_path = Path(argv[argv.index("--filename") + 1])
            zip_path.parent.mkdir(parents=True, exist_ok=True)
            with zipfile.ZipFile(zip_path, "w") as z:
                z.writestr(
                    "ncbi_dataset/data/GCF_000237925.1/protein.faa",
                    ">p1\nMAA\n",
                )
                z.writestr(
                    "ncbi_dataset/data/GCF_000237925.1/cds_from_genomic.fna",
                    ">c1\nATG\n",
                )
            return subprocess.CompletedProcess(argv, 0, "", "")

        monkeypatch.setattr(subprocess, "run", fake_subprocess_run)

        manifest = tmp_path / "manifest.tsv"
        cache = tmp_path / "cache"
        self._write_minimal_manifest(manifest)

        rc = dl.main([
            "--manifest", str(manifest),
            "--cache-dir", str(cache),
        ])
        assert rc == 0
        # Paired outputs landed at canonical paths
        assert (cache / "6183_Schistosoma_mansoni.faa").exists()
        assert (cache / "6183_Schistosoma_mansoni.cds.fna").exists()
        # Report was written
        report = cache / "download_report.tsv"
        assert report.exists()
        # Round-trip the report and confirm status='ok' for the species
        prior = dl.read_prior_status(report)
        assert prior == {"GCF_000237925.1": "ok"}

    def test_script_as_subprocess_does_not_raise_NameError(
        self, tmp_path: Path,
    ) -> None:
        """Catch the order-of-definition trap: when the script runs as
        `python <script>.py`, the body executes top-to-bottom, so any
        `if __name__ == '__main__': main()` block placed BEFORE the
        function definitions it transitively references will hit a
        NameError. Pytest imports the module (running all definitions
        before any test) and so misses this — only an actual subprocess
        invocation reproduces it.

        2026-05-21: first sbatch run died at line 398 with
            NameError: name 'read_download_targets' is not defined
        because main() was wired up above the function definitions.
        """
        import subprocess
        manifest = tmp_path / "empty_manifest.tsv"
        # Write just a header so read_download_targets returns []. This
        # exercises the full main() path without any actual download work.
        from build_species_tree_phase1a_inventory import MANIFEST_COLUMNS
        manifest.write_text("\t".join(MANIFEST_COLUMNS) + "\n")

        script = Path(__file__).resolve().parent.parent.parent / "scripts" / "download_species_tree_phase1a.py"
        result = subprocess.run(
            [
                "python", str(script),
                "--manifest", str(manifest),
                "--cache-dir", str(tmp_path / "cache"),
            ],
            capture_output=True, text=True, timeout=30,
        )
        # NameError, ImportError, anything top-level-broken should NOT occur.
        assert "NameError" not in result.stderr, result.stderr
        assert "ImportError" not in result.stderr, result.stderr
        # rc=0 is expected (no targets, nothing to fail).
        assert result.returncode == 0, (
            f"rc={result.returncode}, stderr={result.stderr!r}"
        )
