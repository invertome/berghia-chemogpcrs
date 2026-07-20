"""Bead v91g: read_download_targets must dedup on assembly accession, not taxid alone.

`download_species_tree_phase1f_genomes.read_download_targets` deduplicated on
taxid only. (The bead text called it `read_targets`; that is the name of the
SAME-SHAPED reader in build_braker4_samples_csv.py, which carries an identical
taxid-only dedup and is reported separately as an unowned sibling defect.) The accession guard lives in `build_genome_inventory.append_entries` —
the WRITE seam — so it blocks new appends that collide on accession but does
nothing about rows ALREADY in the manifest. `find_duplicate_accessions`
reports such rows but never acts on them.

Consequence (real, not hypothetical): GCA_055670145.1 sat in the manifest
under both taxid 205083 (Dreissena rostriformis) and 427924 (D. r. bugensis)
and was downloaded twice — ~1.5 GB of byte-identical duplicate, then annotated
twice by BRAKER4 and treated as two species by OrthoFinder/CAFE.

FIXTURE PROVENANCE — verified 2026-07-20 against the tracked manifest
references/species_tree/genome_inventory.tsv (466 rows):
  * the Dreissena pair is already reconciled in the DATA by the orchestrator
    (commit 3f16f6b); this covers the CODE defect so it cannot recur.
  * 2 duplicate version-stripped accessions remain — GCA_034509925 and
    GCA_015776775 — but BOTH pairs share one taxid (synonym binomials, e.g.
    Crassostrea/Magallana hongkongensis), so taxid dedup already collapses
    them and this change is behaviour-neutral on today's manifest.
  * 0 GCA/GCF cross-namespace numeric collisions exist, so that gap is latent.
  * accession column values carry the '.N' version suffix, so version
    stripping is load-bearing, not cosmetic.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

import download_species_tree_phase1f_genomes as dl  # noqa: E402
from build_genome_inventory import base_accession  # noqa: E402


_COLS = ["taxid", "binomial", "clade", "accession", "drop_reason"]


def _write_manifest(path: Path, rows: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=_COLS, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in _COLS})
    return path


class TestAccessionDedup:
    def test_same_accession_two_taxids_downloads_once(self, tmp_path: Path) -> None:
        """THE BUG: the real Dreissena pair — one assembly, two taxonomic ranks.

        Taxid dedup alone lets both through and the genome is fetched twice.
        """
        m = _write_manifest(tmp_path / "genome_inventory.tsv", [
            {"taxid": "205083", "binomial": "Dreissena rostriformis",
             "clade": "Bivalvia", "accession": "GCA_055670145.1"},
            {"taxid": "427924", "binomial": "Dreissena rostriformis bugensis",
             "clade": "Bivalvia", "accession": "GCA_055670145.1"},
        ])
        targets = dl.read_download_targets(m)
        assert len(targets) == 1, (
            "one assembly enrolled under two taxids was downloaded twice: "
            f"{[(t.taxid, t.accession) for t in targets]}"
        )
        assert targets[0].taxid == 205083, "first-seen row must win"

    def test_dedup_is_version_insensitive(self, tmp_path: Path) -> None:
        """Versions are revisions of ONE physical assembly, per base_accession."""
        m = _write_manifest(tmp_path / "genome_inventory.tsv", [
            {"taxid": "205083", "binomial": "Dreissena rostriformis",
             "clade": "Bivalvia", "accession": "GCA_055670145.1"},
            {"taxid": "427924", "binomial": "Dreissena rostriformis bugensis",
             "clade": "Bivalvia", "accession": "GCA_055670145.2"},
        ])
        assert len(dl.read_download_targets(m)) == 1

    def test_dedup_spans_multiple_manifests(self, tmp_path: Path) -> None:
        """Phase 1d and Phase 1e are read as separate manifests; the accession
        guard must hold ACROSS them, not just within one file."""
        a = _write_manifest(tmp_path / "phase1e.tsv", [
            {"taxid": "205083", "binomial": "Dreissena rostriformis",
             "clade": "Bivalvia", "accession": "GCA_055670145.1"},
        ])
        b = _write_manifest(tmp_path / "phase1d.tsv", [
            {"taxid": "427924", "binomial": "Dreissena rostriformis bugensis",
             "clade": "Bivalvia", "accession": "GCA_055670145.1"},
        ])
        assert len(dl.read_download_targets(a, b)) == 1

    def test_distinct_accessions_are_both_kept(self, tmp_path: Path) -> None:
        """The guard must not over-collapse genuinely different assemblies."""
        m = _write_manifest(tmp_path / "genome_inventory.tsv", [
            {"taxid": "205083", "binomial": "Dreissena rostriformis",
             "clade": "Bivalvia", "accession": "GCA_055670145.1"},
            {"taxid": "6550", "binomial": "Mytilus edulis",
             "clade": "Bivalvia", "accession": "GCA_019925275.1"},
        ])
        assert len(dl.read_download_targets(m)) == 2

    def test_blank_accession_rows_are_skipped_not_collapsed(
        self, tmp_path: Path,
    ) -> None:
        """Blank is never a dedup key (base_accession contract). Rows without
        an accession are already skipped outright, so a blank must not become
        a bucket that swallows a second blank-accession species."""
        m = _write_manifest(tmp_path / "genome_inventory.tsv", [
            {"taxid": "1", "binomial": "Alpha one", "clade": "X",
             "accession": ""},
            {"taxid": "2", "binomial": "Beta two", "clade": "X",
             "accession": ""},
        ])
        assert dl.read_download_targets(m) == []

    def test_taxid_dedup_still_applies(self, tmp_path: Path) -> None:
        """Pre-existing behaviour: repeated taxid collapses to the first row."""
        m = _write_manifest(tmp_path / "genome_inventory.tsv", [
            {"taxid": "205083", "binomial": "Dreissena rostriformis",
             "clade": "Bivalvia", "accession": "GCA_055670145.1"},
            {"taxid": "205083", "binomial": "Dreissena rostriformis",
             "clade": "Bivalvia", "accession": "GCA_999999999.1"},
        ])
        targets = dl.read_download_targets(m)
        assert len(targets) == 1
        assert targets[0].accession == "GCA_055670145.1"

    def test_dropped_rows_do_not_reserve_their_accession(
        self, tmp_path: Path,
    ) -> None:
        """A row with a drop_reason is skipped before dedup, so it must not
        block a live row that legitimately carries the same assembly."""
        m = _write_manifest(tmp_path / "genome_inventory.tsv", [
            {"taxid": "205083", "binomial": "Dreissena rostriformis",
             "clade": "Bivalvia", "accession": "GCA_055670145.1",
             "drop_reason": "harvested_annotated"},
            {"taxid": "427924", "binomial": "Dreissena rostriformis bugensis",
             "clade": "Bivalvia", "accession": "GCA_055670145.1"},
        ])
        targets = dl.read_download_targets(m)
        assert len(targets) == 1
        assert targets[0].taxid == 427924


class TestSharedHelperReuse:
    """The accession key must come from ONE definition, not a second copy.

    Two independently-written normalisers drifting apart is how the write seam
    and the read seam disagreed in the first place.
    """

    def test_read_targets_uses_build_genome_inventory_base_accession(self) -> None:
        assert dl.base_accession is base_accession

    def test_version_stripping_contract(self) -> None:
        assert base_accession("GCA_055670145.1") == "GCA_055670145"
        assert base_accession("") == ""


class TestRealManifestInvariant:
    """Guard the live manifest against the defect this bead closes."""

    def test_tracked_manifest_has_no_cross_taxid_duplicate_accession(self) -> None:
        manifest = (Path(__file__).resolve().parents[2]
                    / "references" / "species_tree" / "genome_inventory.tsv")
        if not manifest.exists():
            pytest.skip("genome_inventory.tsv not present in this checkout")

        by_acc: dict[str, set[str]] = {}
        with manifest.open() as f:
            for row in csv.DictReader(f, delimiter="\t"):
                acc = base_accession(row.get("accession", ""))
                if acc:
                    by_acc.setdefault(acc, set()).add(str(row["taxid"]))

        offenders = {a: t for a, t in by_acc.items() if len(t) > 1}
        assert not offenders, (
            "one assembly enrolled under multiple taxids — the Dreissena "
            f"failure mode has recurred: {offenders}"
        )
