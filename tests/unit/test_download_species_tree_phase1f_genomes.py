"""Tests for scripts/download_species_tree_phase1f_genomes.py — read_download_targets.

Covers the drop_reason exclusion behavior: rows with a non-empty drop_reason
must be skipped even when they carry a valid accession; rows with an empty
drop_reason and a valid accession must be included.
"""
from __future__ import annotations

from pathlib import Path

import pytest

import download_species_tree_phase1f_genomes as dl


# 14-column genome_inventory header (same schema as genome_inventory.tsv)
_GI_COLS = (
    "taxid", "binomial", "clade", "policy_class", "source",
    "accession", "assembly_level", "annotation_status",
    "est_protein_count", "submission_date", "contig_n50",
    "total_length_bp", "drop_reason", "source_batch",
)


def _write_gi(path: Path, rows: list[dict]) -> None:
    with path.open("w") as f:
        f.write("\t".join(_GI_COLS) + "\n")
        for r in rows:
            f.write("\t".join(r.get(c, "") for c in _GI_COLS) + "\n")


class TestReadDownloadTargetsDropReason:
    """read_download_targets must honor drop_reason in the manifest.

    A non-empty drop_reason means the row was intentionally excluded from
    the BRAKER target set (e.g. its proteome was already harvested from NCBI).
    Rows with an empty drop_reason and a valid accession are still included.
    """

    def test_drop_reason_row_excluded(self, tmp_path: Path) -> None:
        m = tmp_path / "gi.tsv"
        _write_gi(m, [
            {"taxid": "231223", "binomial": "Elysia crispata",
             "clade": "gastropoda", "accession": "GCA_033675545.1",
             "drop_reason": "harvested_annotated"},
        ])
        targets = dl.read_download_targets(m)
        assert targets == []

    def test_empty_drop_reason_row_included(self, tmp_path: Path) -> None:
        m = tmp_path / "gi.tsv"
        _write_gi(m, [
            {"taxid": "6161", "binomial": "Dugesia japonica",
             "clade": "platyhelminthes", "accession": "GCA_001938525.1",
             "drop_reason": ""},
        ])
        targets = dl.read_download_targets(m)
        assert len(targets) == 1
        assert targets[0].taxid == 6161

    def test_mixed_rows_only_empty_drop_reason_included(self, tmp_path: Path) -> None:
        m = tmp_path / "gi.tsv"
        _write_gi(m, [
            {"taxid": "231223", "binomial": "Elysia crispata",
             "clade": "gastropoda", "accession": "GCA_033675545.1",
             "drop_reason": "harvested_annotated"},
            {"taxid": "6161", "binomial": "Dugesia japonica",
             "clade": "platyhelminthes", "accession": "GCA_001938525.1",
             "drop_reason": ""},
        ])
        targets = dl.read_download_targets(m)
        assert len(targets) == 1
        assert targets[0].taxid == 6161

    def test_whitespace_only_drop_reason_treated_as_empty(
        self, tmp_path: Path
    ) -> None:
        """A drop_reason that is whitespace only after strip() is falsy — included."""
        m = tmp_path / "gi.tsv"
        _write_gi(m, [
            {"taxid": "231223", "binomial": "Elysia crispata",
             "clade": "gastropoda", "accession": "GCA_033675545.1",
             "drop_reason": "  "},
        ])
        targets = dl.read_download_targets(m)
        assert len(targets) == 1
        assert targets[0].taxid == 231223

    def test_absent_drop_reason_column_still_included(self, tmp_path: Path) -> None:
        """Manifests without a drop_reason column (old schema) still work."""
        p = tmp_path / "old.tsv"
        # Write without drop_reason column — simulates older manifest schema
        cols = ("taxid", "binomial", "clade", "accession")
        with p.open("w") as f:
            f.write("\t".join(cols) + "\n")
            f.write("6161\tDugesia japonica\tplatyhelminthes\tGCA_001938525.1\n")
        targets = dl.read_download_targets(p)
        assert len(targets) == 1
        assert targets[0].taxid == 6161
