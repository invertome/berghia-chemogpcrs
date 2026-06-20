"""Unified genome-inventory builder (replaces the phase1e + phase1d split).

Reads clade_policies.tsv + the existing unified genome_inventory.tsv, queries
NCBI datasets per clade (reusing the proven selection core), and APPENDS only
newly-discovered species tagged with --batch-tag. Never mutates existing rows;
re-run with no new deposits → no change. See
docs/plans/2026-06-19-unified-genome-inventory-sourcing-design.md.
"""
from __future__ import annotations

import argparse
import csv
import datetime
from pathlib import Path

from build_species_tree_phase1d_extension_inventory import (
    ExtensionEntry,
    build_extension_inventory,
    load_existing_taxids,
    query_datasets_for_clade,
)
from genome_inventory_policies import DEFAULT_POLICY_CONFIG, load_clade_policies
from migrate_genome_inventory import UNIFIED_COLUMNS, write_unified


def read_existing(path: str | Path) -> list[dict]:
    """Read an existing unified manifest into superset-keyed dicts. Returns []
    when the file is absent (first run)."""
    p = Path(path)
    if not p.exists():
        return []
    with p.open() as f:
        return [dict(row) for row in csv.DictReader(f, delimiter="\t")]


def _entry_to_unified(entry: ExtensionEntry, source_batch: str) -> dict:
    c = entry.choice
    row = {col: "" for col in UNIFIED_COLUMNS}
    row.update(
        taxid=str(entry.taxid), binomial=entry.binomial,
        clade=entry.clade_name, policy_class=entry.policy_class,
        source=c.source, accession=c.accession,
        assembly_level=c.assembly_level, annotation_status=c.annotation_status,
        est_protein_count=str(c.est_protein_count), submission_date=c.submission_date,
        contig_n50=str(c.contig_n50), total_length_bp=str(c.total_length_bp),
        drop_reason="", source_batch=source_batch,
    )
    return row


def append_entries(existing: list[dict], entries: list[ExtensionEntry],
                   source_batch: str) -> list[dict]:
    """Append new entries (by taxid) to existing rows; existing rows are never
    modified. Returns the combined list sorted by integer taxid."""
    present = {str(r["taxid"]) for r in existing}
    combined = list(existing)
    for e in entries:
        if str(e.taxid) in present:
            continue
        combined.append(_entry_to_unified(e, source_batch))
        present.add(str(e.taxid))
    combined.sort(key=lambda r: int(r["taxid"]))
    return combined


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Unified genome-inventory builder (append-mode)")
    p.add_argument("--manifest", type=Path,
                   default=Path("references/species_tree/genome_inventory.tsv"))
    p.add_argument("--policies", type=Path, default=DEFAULT_POLICY_CONFIG)
    p.add_argument("--refs-root", type=Path,
                   default=Path("references/nath_et_al/one_to_one_ortholog"))
    p.add_argument("--datasets-bin", default="datasets")
    p.add_argument("--timeout", type=int, default=600)
    p.add_argument("--batch-tag", default=None,
                   help="source_batch tag for new rows (default: datasets_YYYYMMDD)")
    args = p.parse_args(argv)

    batch_tag = args.batch_tag or f"datasets_{datetime.date.today():%Y%m%d}"
    existing = read_existing(args.manifest)
    exclude = load_existing_taxids(args.refs_root, args.manifest)
    policies = load_clade_policies(args.policies)

    def query_fn(clade_name: str) -> list[dict]:
        return query_datasets_for_clade(
            clade_name, datasets_bin=args.datasets_bin, timeout=args.timeout)

    entries = build_extension_inventory(
        policies=policies, exclude_taxids=exclude, query_fn=query_fn)
    combined = append_entries(existing, entries, source_batch=batch_tag)
    write_unified(args.manifest, combined)

    n_new = len(combined) - len(existing)
    print(f"[build_genome_inventory] {n_new} new species appended "
          f"(batch={batch_tag}); manifest now {len(combined)} rows -> {args.manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
