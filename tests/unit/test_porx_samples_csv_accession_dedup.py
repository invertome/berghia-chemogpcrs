"""porx: build_braker4_samples_csv.read_targets must dedup on assembly accession.

Taxid-only dedup lets ONE physical genome enrol as TWO organisms, and this
reader is the one that decides who gets BRAKER4-annotated — so a duplicate here
costs a multi-week annotation slot, not just disk.

The fixture below is the real record, copied from
``references/species_tree/genome_inventory.tsv`` on the Unity checkout
``/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05``
(verified 2026-07-20, lines 171 and 243): assembly ``GCA_055670145.1`` appears
under taxid 205083 (*Dreissena rostriformis*) and taxid 427924
(*D. r. bugensis*), both with an empty ``drop_reason``. Both rows reached the
live ``samples_frozen_fix.csv``; only ``205083_Dreissena_rostriformis`` has a
genome in ``species_tree_data/braker4_genomes/`` and an output directory in
``braker4_run/output/``.

An accession sweep of that manifest (408 eligible rows, 405 distinct base
accessions) found exactly three duplicate groups; in all three the row this
dedup drops has NO BRAKER4 output, so the guard invalidates no completed work.
"""
from __future__ import annotations

from pathlib import Path

import build_braker4_samples_csv as bs
import build_genome_inventory

# Header of the real genome_inventory.tsv (verified 2026-07-20).
_INVENTORY_COLS = (
    "taxid", "binomial", "clade", "policy_class", "source", "accession",
    "assembly_level", "annotation_status", "est_protein_count",
    "submission_date", "contig_n50", "total_length_bp", "drop_reason",
    "source_batch",
)

_DREISSENA_ROWS = [
    {"taxid": "205083", "binomial": "Dreissena rostriformis", "clade": "bivalvia",
     "source": "GenBank", "accession": "GCA_055670145.1",
     "assembly_level": "Chromosome", "source_batch": "nath_phase1e"},
    {"taxid": "427924", "binomial": "Dreissena rostriformis bugensis",
     "clade": "Bivalvia", "policy_class": "other_mollusca", "source": "GenBank",
     "accession": "GCA_055670145.1", "assembly_level": "Chromosome",
     "source_batch": "extension_phase1d"},
]


def _write_tsv(path: Path, rows) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("\t".join(_INVENTORY_COLS) + "\n")
        for r in rows:
            f.write("\t".join(r.get(c, "") for c in _INVENTORY_COLS) + "\n")
    return path


class TestAccessionDedup:

    def test_one_assembly_under_two_taxids_yields_one_target(self, tmp_path: Path) -> None:
        """The concrete GCA_055670145.1 case. Taxid dedup cannot see it."""
        m = _write_tsv(tmp_path / "genome_inventory.tsv", _DREISSENA_ROWS)
        targets = bs.read_targets(m)
        assert len(targets) == 1
        # First-seen wins, and the manifest is taxid-sorted, so the retained row
        # is 205083 — the one that actually has the cached genome and the
        # completed BRAKER4 output directory on Unity.
        assert targets[0].taxid == 205083

    def test_duplicate_across_two_manifests_is_caught(self, tmp_path: Path) -> None:
        """The legacy two-manifest form must dedup across files, not just within."""
        a = _write_tsv(tmp_path / "phase1e.tsv", _DREISSENA_ROWS[:1])
        b = _write_tsv(tmp_path / "phase1d.tsv", _DREISSENA_ROWS[1:])
        targets = bs.read_targets(a, b)
        assert [t.taxid for t in targets] == [205083]

    def test_accession_versions_are_the_same_assembly(self, tmp_path: Path) -> None:
        """`.1` and `.2` are revisions of one genome, not two genomes."""
        rows = [
            {"taxid": "111", "binomial": "Genus alpha", "accession": "GCA_012345678.1"},
            {"taxid": "222", "binomial": "Genus beta", "accession": "GCA_012345678.2"},
        ]
        m = _write_tsv(tmp_path / "genome_inventory.tsv", rows)
        assert [t.taxid for t in bs.read_targets(m)] == [111]

    def test_reuses_the_inventory_builders_normaliser(self) -> None:
        """One normaliser, not two.

        ``build_genome_inventory.append_entries`` guards the manifest WRITE seam
        with this key; a divergent copy on a READ seam is how the two seams
        disagreed in the first place.
        """
        assert bs.base_accession is build_genome_inventory.base_accession

    def test_blank_accession_never_becomes_a_dedup_key(self, tmp_path: Path) -> None:
        """Two blank-accession rows must not collapse into one another.

        They are already skipped for being blank; this pins that the accession
        guard did not turn "" into a shared key that swallows real species.
        """
        rows = [
            {"taxid": "111", "binomial": "Genus alpha", "accession": ""},
            {"taxid": "222", "binomial": "Genus beta", "accession": ""},
            {"taxid": "333", "binomial": "Genus gamma", "accession": "GCA_9.1"},
        ]
        m = _write_tsv(tmp_path / "genome_inventory.tsv", rows)
        assert [t.taxid for t in bs.read_targets(m)] == [333]

    def test_distinct_assemblies_all_survive(self, tmp_path: Path) -> None:
        """Regression: the guard must not over-collapse unrelated species."""
        rows = [
            {"taxid": "111", "binomial": "Genus alpha", "accession": "GCA_000000001.1"},
            {"taxid": "222", "binomial": "Genus beta", "accession": "GCA_000000002.1"},
            {"taxid": "333", "binomial": "Genus gamma", "accession": "GCF_000000003.1"},
        ]
        m = _write_tsv(tmp_path / "genome_inventory.tsv", rows)
        assert [t.taxid for t in bs.read_targets(m)] == [111, 222, 333]

    def test_taxid_dedup_still_applies(self, tmp_path: Path) -> None:
        """Regression: same taxid, different assemblies -> still one target.

        Real case from the manifest sweep: taxid 1348078 appears twice under
        synonym binomials (Hirudinaria / Poecilobdella manillensis).
        """
        rows = [
            {"taxid": "1348078", "binomial": "Hirudinaria manillensis",
             "accession": "GCA_034509925.1"},
            {"taxid": "1348078", "binomial": "Poecilobdella manillensis",
             "accession": "GCA_034509926.1"},
        ]
        m = _write_tsv(tmp_path / "genome_inventory.tsv", rows)
        targets = bs.read_targets(m)
        assert len(targets) == 1
        assert targets[0].binomial == "Hirudinaria manillensis"

    def test_dropped_rows_are_not_enrolled_and_do_not_claim_an_accession(
        self, tmp_path: Path
    ) -> None:
        """A drop_reason row must not reserve its accession against a live row.

        The harvest step marks annotated species drop_reason='harvested_annotated'
        in genome_inventory.tsv; if that row still claimed the accession key it
        would evict a legitimate de-novo target sharing it.
        """
        rows = [
            {"taxid": "111", "binomial": "Genus alpha", "accession": "GCA_000000001.1",
             "drop_reason": "harvested_annotated"},
            {"taxid": "222", "binomial": "Genus beta", "accession": "GCA_000000001.1"},
        ]
        m = _write_tsv(tmp_path / "genome_inventory.tsv", rows)
        assert [t.taxid for t in bs.read_targets(m)] == [222]
