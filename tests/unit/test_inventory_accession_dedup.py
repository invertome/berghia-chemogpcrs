"""Accession-level deduplication in the unified genome-inventory builder.

Regression cover for the 2026-07 audit finding: ``append_entries`` deduped
by taxid ONLY, so the same physical assembly could enter the manifest twice
under two different taxids. Observed live in
``references/species_tree/genome_inventory.tsv``:

    line 171  taxid 205083  Dreissena rostriformis            GCA_055670145.1
    line 243  taxid 427924  Dreissena rostriformis bugensis   GCA_055670145.1

NCBI (datasets v2alpha, verified 2026-07-19) reports GCA_055670145.1 as
taxid 427924 / *Dreissena rostriformis bugensis*, so the 205083 row is a
mislabelled copy of the same genome. Consequence: one genome downloaded and
BRAKER4-annotated twice, then entering OrthoFinder as two species.

Two further same-accession pairs already exist in the manifest (both are
same-taxid synonym binomials, so the taxid guard did not catch them either):

    lines 335/336  taxid 1348078  GCA_034509925.1  Hirudinaria / Poecilobdella manillensis
    lines 401/402  taxid 2653900  GCA_015776775.1  Crassostrea / Magallana hongkongensis
"""
from __future__ import annotations

import build_genome_inventory as bgi
from build_species_tree_phase1a_inventory import AssemblyChoice
from build_species_tree_phase1d_extension_inventory import ExtensionEntry
from migrate_genome_inventory import UNIFIED_COLUMNS


def _entry(taxid: int, accession: str) -> ExtensionEntry:
    return ExtensionEntry(
        taxid=taxid, binomial=f"Genus sp{taxid}",
        policy_class="other_mollusca", clade_name="Bivalvia",
        choice=AssemblyChoice(
            accession=accession, source="GenBank",
            assembly_level="Chromosome", annotation_status="",
            est_protein_count=0, submission_date="",
            contig_n50=1_391_067, total_length_bp=1_613_145_916,
        ),
    )


def _row(taxid: str, accession: str, **kw) -> dict:
    row = {c: "" for c in UNIFIED_COLUMNS}
    row.update(taxid=taxid, accession=accession, **kw)
    return row


# ---------------------------------------------------------------- dedup

def test_append_skips_entry_whose_accession_already_present() -> None:
    """The live Dreissena case: same accession, different taxid."""
    existing = [_row("427924", "GCA_055670145.1", source_batch="extension_phase1d")]
    out = bgi.append_entries(existing, [_entry(205083, "GCA_055670145.1")],
                             source_batch="datasets_20260720")
    assert len(out) == 1, "same-accession entry must not be appended under a new taxid"
    assert out[0]["taxid"] == "427924"


def test_append_dedups_accession_ignoring_version_suffix() -> None:
    """GCA_X.1 and GCA_X.2 are revisions of one assembly, not two genomes.

    Downloading + annotating both burns the same compute twice, so the guard
    matches on the versionless base accession.
    """
    existing = [_row("427924", "GCA_055670145.1")]
    out = bgi.append_entries(existing, [_entry(205083, "GCA_055670145.2")],
                             source_batch="b")
    assert len(out) == 1


def test_append_dedups_accession_among_incoming_entries() -> None:
    """Two new entries sharing one accession: only the first is kept."""
    out = bgi.append_entries([], [_entry(205083, "GCA_055670145.1"),
                                  _entry(427924, "GCA_055670145.1")],
                             source_batch="b")
    assert len(out) == 1


def test_append_still_dedups_by_taxid() -> None:
    existing = [_row("200", "GCA_200.1", source_batch="extension_phase1d")]
    out = bgi.append_entries(existing, [_entry(200, "GCA_999.1")], source_batch="b")
    assert len(out) == 1
    assert out[0]["source_batch"] == "extension_phase1d"


def test_append_allows_distinct_accessions_and_taxids() -> None:
    existing = [_row("100", "GCA_100.1")]
    out = bgi.append_entries(existing, [_entry(200, "GCA_200.1")], source_batch="b")
    assert len(out) == 2
    assert {r["taxid"] for r in out} == {"100", "200"}


def test_append_ignores_blank_accessions_when_deduping() -> None:
    """Rows with no accession must not collapse into each other."""
    existing = [_row("100", ""), _row("101", "")]
    out = bgi.append_entries(existing, [_entry(200, "")], source_batch="b")
    assert len(out) == 3


# ------------------------------------------------- pre-existing dup audit

def test_find_duplicate_accessions_reports_existing_collisions() -> None:
    """The builder must be able to SURFACE dups already in the manifest.

    ``append_entries`` never mutates existing rows, so the three live
    duplicate pairs stay until the user reconciles them. The audit helper
    makes them visible instead of silent.
    """
    rows = [
        _row("205083", "GCA_055670145.1", binomial="Dreissena rostriformis"),
        _row("427924", "GCA_055670145.1", binomial="Dreissena rostriformis bugensis"),
        _row("100", "GCA_100.1", binomial="Unique one"),
    ]
    dups = bgi.find_duplicate_accessions(rows)
    assert set(dups) == {"GCA_055670145"}
    assert [r["taxid"] for r in dups["GCA_055670145"]] == ["205083", "427924"]


def test_find_duplicate_accessions_clean_manifest_is_empty() -> None:
    rows = [_row("100", "GCA_100.1"), _row("200", "GCA_200.1"), _row("300", "")]
    assert bgi.find_duplicate_accessions(rows) == {}
