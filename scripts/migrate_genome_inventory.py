"""One-time merge of the two genome-inventory manifests into one.

M1 = references/species_tree/genome_inventory_unannotated.tsv (Phase 1e, per-taxid)
M2 = references/species_tree/extension_inventory.tsv          (Phase 1d, per-clade)

Harmonizes to a 14-column superset with a `source_batch` provenance tag.
Does NOT re-query NCBI — frozen at the current union. See
docs/plans/2026-06-19-unified-genome-inventory-sourcing-design.md.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

UNIFIED_COLUMNS = (
    "taxid", "binomial", "clade", "policy_class", "source", "accession",
    "assembly_level", "annotation_status", "est_protein_count",
    "submission_date", "contig_n50", "total_length_bp", "drop_reason",
    "source_batch",
)


def _row_to_unified(row: dict, source_batch: str) -> dict:
    """Map an M1 (clade/drop_reason) or M2 (clade_name/policy_class) row onto
    the superset. Missing columns default to ''; M2's clade_name fills clade."""
    out = {col: "" for col in UNIFIED_COLUMNS}
    for col in UNIFIED_COLUMNS:
        val = row.get(col)
        if val is not None:
            out[col] = val
    if not out["clade"]:
        out["clade"] = (row.get("clade_name") or "").strip()
    out["source_batch"] = source_batch
    return out


def read_manifest_rows(path: str | Path, source_batch: str) -> list[dict]:
    with Path(path).open() as f:
        rows: list[dict] = []
        for row in csv.DictReader(f, delimiter="\t"):
            if not (row.get("taxid") or "").strip():
                continue  # skip blank lines / taxid-less rows (read_targets skips them too)
            rows.append(_row_to_unified(row, source_batch))
        return rows


def merge_manifests(m1_path: str | Path, m2_path: str | Path) -> list[dict]:
    rows = read_manifest_rows(m1_path, "nath_phase1e")
    rows += read_manifest_rows(m2_path, "extension_phase1d")
    rows.sort(key=lambda r: int(r["taxid"]))
    return rows


def write_unified(path: str | Path, rows: list[dict]) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(UNIFIED_COLUMNS), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Merge M1+M2 → unified genome_inventory.tsv")
    p.add_argument("--m1", type=Path,
                   default=Path("references/species_tree/genome_inventory_unannotated.tsv"))
    p.add_argument("--m2", type=Path,
                   default=Path("references/species_tree/extension_inventory.tsv"))
    p.add_argument("--out", type=Path,
                   default=Path("references/species_tree/genome_inventory.tsv"))
    args = p.parse_args(argv)
    rows = merge_manifests(args.m1, args.m2)
    write_unified(args.out, rows)
    print(f"[migrate_genome_inventory] {len(rows)} rows -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
