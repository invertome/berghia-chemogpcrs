"""Bead spnk: the Phase-1f producer/consumer filename seam.

`postprocess_braker4_outputs.postprocess_one` writes the per-species proteome
as `<cache>/<sample>.aa.fna`, while
`consolidate_proteomes_for_genome_wide_og._locate_proteome` looked for
`<cache>/<leaf>.faa` in the SAME directory. Same dir, different suffix, no
probing — so every Phase-1f species resolved to a non-existent path and was
reported `missing_proteome` while the job exited 0.

FIXTURE PROVENANCE — these fixtures encode what the producer actually emits,
not what the consumer assumed. Verified 2026-07-20 against the live Unity
checkout at
  /scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
  references/species_tree/cache/proteomes_braker4/   -> 264 files,
      132 x '<leaf>.aa.fna' + 132 x '<leaf>.cds.fna', ZERO '*.faa'
  species_tree_data/consolidation_report.tsv         -> 538 missing_proteome,
      16 duplicate_taxid, 0 ok
  species_tree_data/braker4_run/output/<leaf>/results/ -> 'braker.aa.gz'
      (gzipped; there is no plain 'braker.aa' on disk)

The pre-existing tests wrote '<leaf>.faa' fixtures — the ASSUMED name — which
is exactly why the suite stayed green while 100% of real species failed.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

import consolidate_proteomes_for_genome_wide_og as cons  # noqa: E402
import postprocess_braker4_outputs as post  # noqa: E402


_PHASE1E_COLS = ["taxid", "binomial", "clade", "accession"]


def _write_tsv(path: Path, cols: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _args(tmp_path: Path, **over):
    base = dict(
        phase1a_manifest=None,
        phase1d_manifest=None,
        phase1e_manifest=None,
        phase1g_manifest=None,
        manifest=None,
        berghia_proteome=None,
        base_dir=tmp_path,
        braker4_output_dir=tmp_path / "braker4_output",
    )
    base.update(over)
    return argparse.Namespace(**base)


def _braker4_cache(tmp_path: Path) -> Path:
    d = tmp_path / "references" / "species_tree" / "cache" / "proteomes_braker4"
    d.mkdir(parents=True, exist_ok=True)
    return d


class TestPhase1fSuffixProbe:
    """_locate_proteome must find what the producer actually writes."""

    def test_locates_real_producer_aa_fna_suffix(self, tmp_path: Path) -> None:
        """THE BUG: producer writes '<leaf>.aa.fna'; consumer must find it.

        Fixture name copied from the live Unity cache listing (see module
        docstring) rather than from the consumer's assumption.
        """
        cache = _braker4_cache(tmp_path)
        (cache / "102321_Hyotissa_hyotis.aa.fna").write_text(">p1\nMAA\n")

        got = cons._locate_proteome(
            102321, "Hyotissa hyotis", "102321_Hyotissa_hyotis", "1f",
            tmp_path, tmp_path / "braker4_output",
        )
        assert got.exists(), f"did not resolve producer output; probed {got}"
        assert got.name == "102321_Hyotissa_hyotis.aa.fna"

    def test_still_locates_legacy_faa(self, tmp_path: Path) -> None:
        """Back-compat: a '.faa' proteome must keep resolving."""
        cache = _braker4_cache(tmp_path)
        (cache / "300_Testus_extensus.faa").write_text(">p\nM\n")

        got = cons._locate_proteome(
            300, "Testus extensus", "300_Testus_extensus", "1f",
            tmp_path, tmp_path / "braker4_output",
        )
        assert got.exists()
        assert got.name == "300_Testus_extensus.faa"

    def test_faa_wins_when_both_present(self, tmp_path: Path) -> None:
        """Deterministic precedence: the conventional protein suffix wins.

        Keeps a future '.aa.fna' -> '.faa' migration forward-compatible: once
        a species has been migrated, the migrated file is preferred even if
        the legacy one is still on disk.
        """
        cache = _braker4_cache(tmp_path)
        (cache / "300_Testus_extensus.faa").write_text(">canonical\nM\n")
        (cache / "300_Testus_extensus.aa.fna").write_text(">legacy\nM\n")

        got = cons._locate_proteome(
            300, "Testus extensus", "300_Testus_extensus", "1f",
            tmp_path, tmp_path / "braker4_output",
        )
        assert got.name == "300_Testus_extensus.faa"

    def test_cds_fna_is_never_mistaken_for_a_proteome(self, tmp_path: Path) -> None:
        """The producer writes '<leaf>.cds.fna' NUCLEOTIDES beside the protein
        file. A loose '.fna' probe would feed CDS into OrthoFinder as protein.
        """
        cache = _braker4_cache(tmp_path)
        (cache / "102321_Hyotissa_hyotis.cds.fna").write_text(">c1\nATGATG\n")

        got = cons._locate_proteome(
            102321, "Hyotissa hyotis", "102321_Hyotissa_hyotis", "1f",
            tmp_path, tmp_path / "braker4_output",
        )
        assert not got.name.endswith(".cds.fna"), (
            "resolved the nucleotide CDS file as a proteome"
        )

    def test_missing_everything_still_reports_a_canonical_path(
        self, tmp_path: Path,
    ) -> None:
        """With nothing on disk the returned path must still be the canonical
        one, so the missing_proteome message names a useful location."""
        _braker4_cache(tmp_path)
        got = cons._locate_proteome(
            999, "Absent species", "999_Absent_species", "1f",
            tmp_path, tmp_path / "braker4_output",
        )
        assert not got.exists()
        assert "proteomes_braker4" in str(got)


class TestProducerConsumerSeam:
    """The canary: run the real producer, then the real consumer's locator.

    This is the test that would have caught bead spnk. It never hard-codes a
    filename — it lets postprocess_one name the file and asks the consumer to
    find it, so the two sides cannot drift apart again.
    """

    def test_postprocess_output_is_found_by_consolidate(self, tmp_path: Path) -> None:
        leaf = "102321_Hyotissa_hyotis"
        cache = _braker4_cache(tmp_path)
        results = tmp_path / "braker4_output" / leaf / "results"
        results.mkdir(parents=True)
        with gzip.open(results / "braker.aa.gz", "wt") as f:
            f.write(">g1.t1\nMAAWQ\n>g2.t1\nMTTYK\n")
        with gzip.open(results / "braker.codingseq.gz", "wt") as f:
            f.write(">g1.t1\nATGGCC\n")

        produced = post.postprocess_one(
            leaf, tmp_path / "braker4_output", cache,
        )
        assert produced.status == "ok"

        located = cons._locate_proteome(
            102321, "Hyotissa hyotis", leaf, "1f",
            tmp_path, tmp_path / "braker4_output",
        )
        assert located == produced.aa_path, (
            f"producer wrote {produced.aa_path}; consumer looked at {located}"
        )
        assert located.exists()

    def test_full_consolidation_of_producer_output(self, tmp_path: Path) -> None:
        """End to end: producer -> manifest -> load_sources -> consolidate_one."""
        leaf = "102321_Hyotissa_hyotis"
        cache = _braker4_cache(tmp_path)
        results = tmp_path / "braker4_output" / leaf / "results"
        results.mkdir(parents=True)
        with gzip.open(results / "braker.aa.gz", "wt") as f:
            f.write(">g1.t1\nMAAWQ\n>g2.t1\nMTTYK\n")

        post.postprocess_one(leaf, tmp_path / "braker4_output", cache)

        m = tmp_path / "genome_inventory.tsv"
        _write_tsv(m, _PHASE1E_COLS, [{
            "taxid": "102321", "binomial": "Hyotissa hyotis",
            "clade": "Bivalvia", "accession": "GCA_000001.1",
        }])

        sources, dups = cons.load_sources(_args(tmp_path, phase1e_manifest=m))
        assert len(sources) == 1 and not dups

        status = cons.consolidate_one(
            sources[0], out_dir=tmp_path / "of_input",
        )
        assert status.status == "ok", status.message
        assert status.n_seqs == 2
        assert (tmp_path / "of_input" / f"{leaf}.fa").exists()


class TestZeroMatchFailsLoudly:
    """A consolidation that resolves NOTHING must not exit 0.

    Real consequence on Unity: consolidation_report.tsv held 538
    missing_proteome / 0 ok and the job still returned 0, while the wrapper's
    own log told the operator that missing_proteome is the signal to act on.
    """

    def _run_main(self, tmp_path: Path, manifest: Path) -> int:
        return cons.main([
            "--out-dir", str(tmp_path / "of_input" / "input"),
            "--base-dir", str(tmp_path),
            "--braker4-output-dir", str(tmp_path / "braker4_output"),
            "--phase1e-manifest", str(manifest),
        ])

    def test_main_nonzero_when_nothing_consolidates(self, tmp_path: Path) -> None:
        _braker4_cache(tmp_path)
        m = tmp_path / "genome_inventory.tsv"
        _write_tsv(m, _PHASE1E_COLS, [
            {"taxid": "102321", "binomial": "Hyotissa hyotis",
             "clade": "Bivalvia", "accession": "GCA_000001.1"},
            {"taxid": "102329", "binomial": "Pinctada margaritifera",
             "clade": "Bivalvia", "accession": "GCA_000002.1"},
        ])
        assert self._run_main(tmp_path, m) != 0, (
            "zero-match consolidation exited 0 — the silent-success failure"
        )

    def test_main_zero_when_something_consolidates(self, tmp_path: Path) -> None:
        """One resolvable species is enough: partial gaps stay a WARN, and 03c
        keeps the hard cleanliness assertion."""
        cache = _braker4_cache(tmp_path)
        (cache / "102321_Hyotissa_hyotis.aa.fna").write_text(">p1\nMAA\n")
        m = tmp_path / "genome_inventory.tsv"
        _write_tsv(m, _PHASE1E_COLS, [
            {"taxid": "102321", "binomial": "Hyotissa hyotis",
             "clade": "Bivalvia", "accession": "GCA_000001.1"},
            {"taxid": "102329", "binomial": "Pinctada margaritifera",
             "clade": "Bivalvia", "accession": "GCA_000002.1"},
        ])
        assert self._run_main(tmp_path, m) == 0

    def test_main_zero_when_there_is_nothing_to_do(self, tmp_path: Path) -> None:
        """An empty manifest is not a zero-match failure — there were no
        sources to resolve, so there is nothing to be loud about."""
        _braker4_cache(tmp_path)
        m = tmp_path / "genome_inventory.tsv"
        _write_tsv(m, _PHASE1E_COLS, [])
        assert self._run_main(tmp_path, m) == 0
