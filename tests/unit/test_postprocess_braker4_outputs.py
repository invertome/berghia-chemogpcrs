"""Tests for scripts/postprocess_braker4_outputs.py.

Bead -p49 (Phase 1f BRAKER4). After BRAKER4 finishes, walk its per-species
output directories, gunzip braker.aa.gz + braker.codingseq.gz, and
copy/rename to canonical paths in our cache:

    references/species_tree/cache/proteomes_braker4/<sample_name>.aa.fna
    references/species_tree/cache/proteomes_braker4/<sample_name>.cds.fna

`sample_name` is what we put in samples.csv (sanitized
`<taxid>_<binomial>`). Per-species status is recorded in a TSV that
parallels Phase 1a-CDS's download_report.tsv.
"""
from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pytest

import postprocess_braker4_outputs as pp


BRAKER4_COLUMNS = (
    "sample_name", "genome", "genome_masked", "protein_fasta",
    "bam_files", "fastq_r1", "fastq_r2", "sra_ids",
    "varus_genus", "varus_species", "isoseq_bam", "isoseq_fastq",
    "busco_lineage", "reference_gtf",
)


def _write_samples_csv(path: Path, sample_names: list[str]) -> None:
    """Write a minimal samples.csv: just sample_name column populated."""
    import csv as _csv
    with path.open("w", newline="") as f:
        writer = _csv.DictWriter(f, fieldnames=list(BRAKER4_COLUMNS))
        writer.writeheader()
        for name in sample_names:
            row = {col: "" for col in BRAKER4_COLUMNS}
            row["sample_name"] = name
            writer.writerow(row)


def _make_braker4_output(
    output_root: Path, sample_name: str,
    aa_content: str | None = ">p1\nMAA\n>p2\nMBB\n",
    cds_content: str | None = ">c1\nATG\n>c2\nGGG\n",
) -> Path:
    """Create a synthetic BRAKER4 output dir with gzipped files."""
    results = output_root / sample_name / "results"
    results.mkdir(parents=True, exist_ok=True)
    if aa_content is not None:
        with gzip.open(results / "braker.aa.gz", "wt") as f:
            f.write(aa_content)
    if cds_content is not None:
        with gzip.open(results / "braker.codingseq.gz", "wt") as f:
            f.write(cds_content)
    return results


# ----------------------------------------------------------------------
# read_sample_names
# ----------------------------------------------------------------------

class TestReadSampleNames:
    def test_returns_list_in_order(self, tmp_path: Path) -> None:
        path = tmp_path / "samples.csv"
        _write_samples_csv(path, ["6161_Dugesia_japonica", "6162_Girardia_tigrina"])
        names = pp.read_sample_names(path)
        assert names == ["6161_Dugesia_japonica", "6162_Girardia_tigrina"]

    def test_skips_empty_rows(self, tmp_path: Path) -> None:
        path = tmp_path / "s.csv"
        _write_samples_csv(path, ["6161_x", "", "  ", "9999_y"])
        names = pp.read_sample_names(path)
        assert names == ["6161_x", "9999_y"]


# ----------------------------------------------------------------------
# postprocess_one — gunzip + rename + count
# ----------------------------------------------------------------------

class TestPostprocessOne:
    def test_happy_path_produces_paired_outputs(self, tmp_path: Path) -> None:
        output_root = tmp_path / "output"
        cache = tmp_path / "cache"
        sample = "6161_Dugesia_japonica"
        _make_braker4_output(output_root, sample)

        result = pp.postprocess_one(sample, output_root, cache)

        assert result.status == "ok"
        assert result.aa_path == cache / f"{sample}.aa.fna"
        assert result.cds_path == cache / f"{sample}.cds.fna"
        assert result.aa_path.exists()
        assert result.cds_path.exists()
        # Gunzipped content matches original
        assert result.aa_path.read_text().startswith(">p1\nMAA")
        assert result.cds_path.read_text().startswith(">c1\nATG")
        # Record counts
        assert result.n_proteins == 2
        assert result.n_cds == 2

    def test_missing_aa_gives_failed(self, tmp_path: Path) -> None:
        output_root = tmp_path / "output"
        cache = tmp_path / "cache"
        sample = "9999_failz"
        # BRAKER4 didn't produce anything (e.g. failed before training).
        (output_root / sample / "results").mkdir(parents=True)

        result = pp.postprocess_one(sample, output_root, cache)

        assert result.status == "failed"
        assert result.aa_path is None
        assert result.cds_path is None
        # No partial files left behind in cache
        assert not (cache / f"{sample}.aa.fna").exists()

    def test_missing_cds_gives_ok_no_cds(self, tmp_path: Path) -> None:
        output_root = tmp_path / "output"
        cache = tmp_path / "cache"
        sample = "6162_partial"
        _make_braker4_output(output_root, sample, cds_content=None)

        result = pp.postprocess_one(sample, output_root, cache)

        assert result.status == "ok_no_cds"
        assert result.aa_path == cache / f"{sample}.aa.fna"
        assert result.aa_path.exists()
        assert result.cds_path is None
        assert not (cache / f"{sample}.cds.fna").exists()

    def test_idempotent_skip(self, tmp_path: Path) -> None:
        output_root = tmp_path / "output"
        cache = tmp_path / "cache"
        sample = "1_aaa"
        cache.mkdir()
        (cache / f"{sample}.aa.fna").write_text(">p1\nMAA\n")
        (cache / f"{sample}.cds.fna").write_text(">c1\nATG\n")
        # No BRAKER4 output needed; the canonical files already exist.

        result = pp.postprocess_one(sample, output_root, cache)

        assert result.status == "skipped"
        # Canonical paths still report
        assert result.aa_path == cache / f"{sample}.aa.fna"
        assert result.cds_path == cache / f"{sample}.cds.fna"

    def test_force_overwrites_existing(self, tmp_path: Path) -> None:
        output_root = tmp_path / "output"
        cache = tmp_path / "cache"
        sample = "1_aaa"
        cache.mkdir()
        # Pre-existing stale content
        (cache / f"{sample}.aa.fna").write_text(">stale\nMSTALE\n")
        (cache / f"{sample}.cds.fna").write_text(">stale\nATG\n")
        # Fresh BRAKER4 output with different content
        _make_braker4_output(
            output_root, sample,
            aa_content=">fresh\nMFRESH\n", cds_content=">fresh\nGGG\n",
        )

        result = pp.postprocess_one(sample, output_root, cache, force=True)

        assert result.status == "ok"
        # Content matches fresh, not stale
        assert "MFRESH" in result.aa_path.read_text()
        assert "GGG" in result.cds_path.read_text()


# ----------------------------------------------------------------------
# Report TSV
# ----------------------------------------------------------------------

class TestReport:
    REPORT_COLUMNS = ("sample_name", "status", "n_proteins", "n_cds", "error")

    def test_columns_match(self) -> None:
        assert pp.REPORT_COLUMNS == self.REPORT_COLUMNS

    def test_write_then_read(self, tmp_path: Path) -> None:
        results = [
            pp.PostprocessResult(
                sample_name="6161_Dugesia_japonica",
                status="ok",
                aa_path=tmp_path / "6161_Dugesia_japonica.aa.fna",
                cds_path=tmp_path / "6161_Dugesia_japonica.cds.fna",
                n_proteins=15000, n_cds=15000,
            ),
            pp.PostprocessResult(
                sample_name="9999_failz",
                status="failed",
                aa_path=None, cds_path=None,
                error="braker.aa.gz missing",
            ),
        ]
        report = tmp_path / "report.tsv"
        pp.write_report_tsv(report, results)

        with report.open() as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
        assert [r["sample_name"] for r in rows] == [
            "6161_Dugesia_japonica", "9999_failz",
        ]
        assert rows[0]["status"] == "ok"
        assert rows[0]["n_proteins"] == "15000"
        assert rows[1]["status"] == "failed"
        assert "missing" in rows[1]["error"]


# ----------------------------------------------------------------------
# main + CLI smoke
# ----------------------------------------------------------------------

class TestMain:
    def test_end_to_end(self, tmp_path: Path) -> None:
        samples_csv = tmp_path / "samples.csv"
        _write_samples_csv(samples_csv, ["6161_a", "9999_failz"])
        output_root = tmp_path / "output"
        cache = tmp_path / "cache"
        _make_braker4_output(output_root, "6161_a")
        # 9999_failz: no BRAKER4 output dir at all
        report = tmp_path / "report.tsv"

        rc = pp.main([
            "--braker4-output", str(output_root),
            "--samples-csv", str(samples_csv),
            "--cache-dir", str(cache),
            "--report", str(report),
        ])

        # rc=2 because one species failed
        assert rc == 2
        assert (cache / "6161_a.aa.fna").exists()
        assert (cache / "6161_a.cds.fna").exists()
        assert not (cache / "9999_failz.aa.fna").exists()
        # Report contains both rows
        with report.open() as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
        assert len(rows) == 2

    def test_all_ok_returns_zero(self, tmp_path: Path) -> None:
        samples_csv = tmp_path / "samples.csv"
        _write_samples_csv(samples_csv, ["6161_a"])
        output_root = tmp_path / "output"
        cache = tmp_path / "cache"
        _make_braker4_output(output_root, "6161_a")
        report = tmp_path / "report.tsv"

        rc = pp.main([
            "--braker4-output", str(output_root),
            "--samples-csv", str(samples_csv),
            "--cache-dir", str(cache),
            "--report", str(report),
        ])
        assert rc == 0


class TestSubprocessSmoke:
    def test_script_runs_without_NameError(self, tmp_path: Path) -> None:
        import subprocess
        samples_csv = tmp_path / "samples.csv"
        _write_samples_csv(samples_csv, [])
        report = tmp_path / "r.tsv"

        script = (
            Path(__file__).resolve().parent.parent.parent
            / "scripts" / "postprocess_braker4_outputs.py"
        )
        result = subprocess.run(
            [
                "python", str(script),
                "--braker4-output", str(tmp_path / "out"),
                "--samples-csv", str(samples_csv),
                "--cache-dir", str(tmp_path / "cache"),
                "--report", str(report),
            ],
            capture_output=True, text=True, timeout=30,
        )
        assert "NameError" not in result.stderr, result.stderr
        assert "ImportError" not in result.stderr, result.stderr
        assert result.returncode == 0, (
            f"rc={result.returncode}, stderr={result.stderr!r}"
        )
