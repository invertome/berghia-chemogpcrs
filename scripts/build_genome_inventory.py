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
import sys
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


def base_accession(accession: str) -> str:
    """Version-stripped assembly accession ('GCA_055670145.1' -> 'GCA_055670145').

    Versions are revisions of ONE assembly, so two rows differing only by
    version denote the same physical genome — downloading and BRAKER4-annotating
    both burns the same compute twice. Blank stays blank (never a dedup key).
    """
    acc = (accession or "").strip()
    return acc.split(".", 1)[0] if acc else ""


def find_duplicate_accessions(rows: list[dict]) -> dict[str, list[dict]]:
    """Rows sharing one assembly accession, keyed by version-stripped accession.

    ``append_entries`` never mutates existing rows, so duplicates already in the
    manifest survive until the user reconciles them. This surfaces them instead
    of letting them stay silent. Rows without an accession are ignored.
    """
    by_acc: dict[str, list[dict]] = {}
    for r in rows:
        acc = base_accession(r.get("accession", ""))
        if acc:
            by_acc.setdefault(acc, []).append(r)
    return {acc: rs for acc, rs in by_acc.items() if len(rs) > 1}


def append_entries(existing: list[dict], entries: list[ExtensionEntry],
                   source_batch: str) -> list[dict]:
    """Append new entries to existing rows; existing rows are never modified.

    An entry is skipped when its taxid OR its assembly accession is already
    present. The accession guard is what stops one genome from entering the
    manifest twice under two taxids — the failure mode that put
    GCA_055670145.1 in as both taxid 205083 and taxid 427924, i.e. one genome
    downloaded and annotated twice and then treated as two species by
    OrthoFinder/CAFE and as a zero-length sister pair in the species tree.

    Returns the combined list sorted by integer taxid.
    """
    present_taxids = {str(r["taxid"]) for r in existing}
    present_accessions = {a for a in (base_accession(r.get("accession", ""))
                                      for r in existing) if a}
    combined = list(existing)
    for e in entries:
        if str(e.taxid) in present_taxids:
            continue
        acc = base_accession(e.choice.accession)
        if acc and acc in present_accessions:
            print(f"[build_genome_inventory] SKIP taxid {e.taxid} "
                  f"({e.binomial}): accession {e.choice.accession} is already "
                  f"in the manifest under a different taxid — one genome, one row.",
                  file=sys.stderr)
            continue
        combined.append(_entry_to_unified(e, source_batch))
        present_taxids.add(str(e.taxid))
        if acc:
            present_accessions.add(acc)
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

    # Duplicates already in the manifest are NOT auto-removed: rows are
    # write-once and other work indexes into this file. Report them so the
    # user can reconcile at the post-drain cleanup step.
    dups = find_duplicate_accessions(combined)
    for acc, rows in sorted(dups.items()):
        who = "; ".join(f"taxid {r['taxid']} ({r['binomial']}, {r['source_batch']})"
                        for r in rows)
        print(f"[build_genome_inventory] DUPLICATE ACCESSION {acc}: {who} "
              f"— same genome in {len(rows)} rows; reconcile before annotation.",
              file=sys.stderr)
    if dups:
        print(f"[build_genome_inventory] {len(dups)} duplicated accession(s) "
              f"present in the manifest.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
