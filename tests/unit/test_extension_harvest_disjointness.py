"""End-to-end regression for the harvest-annotation exclusion invariant.

Regression target: the original design assumed re-deriving the inventory
would drop harvested species, but the inventory builder is append-only and
the consumers filtered only on accession — so harvested taxids silently
stayed in the BRAKER de-novo annotation set.

This test chains the REAL mark output into BOTH BRAKER consumers and
asserts the single invariant that would have caught that defect:

  A taxid present in proteome_manifest.tsv, after mark, is absent from
  the output of both read_targets (braker4_samples_csv) and
  read_download_targets (phase1f genome downloader); a control taxid
  that is only in genome_inventory with an empty drop_reason remains
  in both.
"""
from __future__ import annotations

from pathlib import Path

import harvest_extension_proteomes as hep
import build_braker4_samples_csv as gen
import download_species_tree_phase1f_genomes as dl


# 14-column genome_inventory header (canonical CRLF production schema).
_GI_HEADER = (
    "taxid\tbinomial\tclade\tpolicy_class\tsource\taccession\t"
    "assembly_level\tannotation_status\test_protein_count\tsubmission_date\t"
    "contig_n50\ttotal_length_bp\tdrop_reason\tsource_batch"
)

_GI_COLS = _GI_HEADER.split("\t")


def _make_gi(tmp_path: Path, rows: list[dict], linesep: str = "\r\n") -> Path:
    r"""Write a minimal genome_inventory.tsv with the standard 14-col header.

    Bytes are written verbatim so `linesep` controls the on-disk terminator.
    Default is CRLF to match the canonical production genome_inventory.tsv.
    """
    p = tmp_path / "genome_inventory.tsv"
    lines = [_GI_HEADER]
    for r in rows:
        lines.append("\t".join(r.get(c, "") for c in _GI_COLS))
    p.write_bytes((linesep.join(lines) + linesep).encode("utf-8"))
    return p


def _make_pm(tmp_path: Path, rows: list[dict]) -> Path:
    """Write a minimal proteome_manifest.tsv (10-col, LF)."""
    cols = (
        "taxid", "binomial", "clade", "source", "accession",
        "assembly_level", "annotation_status", "est_protein_count",
        "submission_date", "drop_reason",
    )
    p = tmp_path / "proteome_manifest.tsv"
    lines = ["\t".join(cols)]
    for r in rows:
        lines.append("\t".join(r.get(c, "") for c in cols))
    p.write_text("\n".join(lines) + "\n")
    return p


class TestHarvestExclusionEndToEnd:
    """Regression: mark → read_targets + read_download_targets disjointness."""

    def test_harvested_taxid_absent_from_both_consumers_control_present(
        self, tmp_path: Path
    ) -> None:
        """After mark, the harvested taxid is excluded by both BRAKER consumers;
        the control taxid (empty drop_reason) remains in both outputs.

        Regression for the no-op-re-derive defect: append-only inventory +
        accession-only consumer filtering left harvested species in the
        BRAKER de-novo annotation set.
        """
        # ------------------------------------------------------------------ #
        # 1. genome_inventory: two rows, both with valid accessions and empty
        #    drop_reason.  111111 will be "harvested"; 222222 is the control.
        # ------------------------------------------------------------------ #
        gi = _make_gi(tmp_path, [
            {
                "taxid": "111111",
                "binomial": "Harvested species",
                "clade": "gastropoda",
                "accession": "GCA_111111.1",
            },
            {
                "taxid": "222222",
                "binomial": "Control species",
                "clade": "gastropoda",
                "accession": "GCA_222222.1",
            },
        ])

        # ------------------------------------------------------------------ #
        # 2. proteome_manifest: contains 111111 but NOT 222222.
        # ------------------------------------------------------------------ #
        pm = _make_pm(tmp_path, [
            {
                "taxid": "111111",
                "binomial": "Harvested species",
                "clade": "gastropoda",
                "source": "GenBank",
                "accession": "GCA_111111.1",
                "assembly_level": "Chromosome",
                "est_protein_count": "12345",
                "submission_date": "2024-01-01",
            },
        ])

        # ------------------------------------------------------------------ #
        # 3. Mark: 111111 must be stamped harvested_annotated; 222222 untouched.
        # ------------------------------------------------------------------ #
        n_marked = hep.mark_harvested_in_genome_inventory(str(gi), str(pm))
        assert n_marked == 1, (
            f"expected mark() to update 1 row, got {n_marked}"
        )

        # ------------------------------------------------------------------ #
        # 4. Both BRAKER consumers read the now-marked genome_inventory.
        # ------------------------------------------------------------------ #
        braker4_targets = gen.read_targets(gi)
        dl_targets = dl.read_download_targets(gi)

        # Extract taxids (both consumers return objects with .taxid int).
        braker4_taxids = {t.taxid for t in braker4_targets}
        dl_taxids = {t.taxid for t in dl_targets}

        # ------------------------------------------------------------------ #
        # 5. Invariant: harvested taxid absent; control taxid present.
        # ------------------------------------------------------------------ #
        assert 111111 not in braker4_taxids, (
            "read_targets() returned the harvested taxid 111111; "
            "consumers must filter on drop_reason='harvested_annotated'"
        )
        assert 111111 not in dl_taxids, (
            "read_download_targets() returned the harvested taxid 111111; "
            "consumers must filter on drop_reason='harvested_annotated'"
        )
        assert 222222 in braker4_taxids, (
            "read_targets() dropped control taxid 222222 (empty drop_reason); "
            "it must be included"
        )
        assert 222222 in dl_taxids, (
            "read_download_targets() dropped control taxid 222222 (empty drop_reason); "
            "it must be included"
        )
