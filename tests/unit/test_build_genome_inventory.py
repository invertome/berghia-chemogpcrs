"""Tests for scripts/build_genome_inventory.py — unified append-mode builder."""
from __future__ import annotations

from pathlib import Path

import build_genome_inventory as bgi
from build_species_tree_phase1d_extension_inventory import ExtensionEntry
from build_species_tree_phase1a_inventory import AssemblyChoice
from migrate_genome_inventory import UNIFIED_COLUMNS


def _entry(taxid: int, policy_class="heterobranchia", clade="Heterobranchia") -> ExtensionEntry:
    return ExtensionEntry(
        taxid=taxid, binomial=f"Genus sp{taxid}",
        policy_class=policy_class, clade_name=clade,
        choice=AssemblyChoice(
            accession=f"GCA_{taxid}", source="GenBank",
            assembly_level="Chromosome", annotation_status="",
            est_protein_count=0, submission_date="2020-01-01",
            contig_n50=90000, total_length_bp=2_000_000,
        ),
    )


def test_append_adds_new_rows_tagged(tmp_path: Path) -> None:
    existing = [{c: "" for c in UNIFIED_COLUMNS}]
    existing[0].update(taxid="100", binomial="Old one", clade="gastropoda",
                       accession="GCA_100", source_batch="nath_phase1e")
    out = bgi.append_entries(existing, [_entry(200)], source_batch="datasets_20260619")
    by_taxid = {r["taxid"]: r for r in out}
    assert set(by_taxid) == {"100", "200"}
    assert by_taxid["100"]["source_batch"] == "nath_phase1e"      # unchanged
    assert by_taxid["200"]["source_batch"] == "datasets_20260619"  # tagged
    assert by_taxid["200"]["clade"] == "Heterobranchia"
    assert by_taxid["200"]["policy_class"] == "heterobranchia"


def test_append_is_idempotent_on_existing_taxid(tmp_path: Path) -> None:
    existing = [{c: "" for c in UNIFIED_COLUMNS}]
    existing[0].update(taxid="200", binomial="Already here", source_batch="extension_phase1d")
    out = bgi.append_entries(existing, [_entry(200)], source_batch="datasets_20260619")
    assert len(out) == 1
    assert out[0]["source_batch"] == "extension_phase1d"  # not overwritten


def test_append_sorted_by_taxid() -> None:
    existing = [{c: "" for c in UNIFIED_COLUMNS}]
    existing[0].update(taxid="300", source_batch="x")
    out = bgi.append_entries(existing, [_entry(100), _entry(200)], source_batch="b")
    assert [int(r["taxid"]) for r in out] == [100, 200, 300]
