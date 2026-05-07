#!/usr/bin/env python3
"""select_backbone_reps.py — Pick representative sequences per coarse
GPCR family for the backbone reference tree.

Phase 3 Task 3.1 of the non-chemoreceptor classification pipeline.

For each coarse family in the curated reference TSV, select N
representatives chosen to:
  (a) hit the per-family quota (default 5; small families return all),
  (b) include at least one invertebrate (Drosophila / C. elegans) ref
      when the family has any (anchors phylogenetic placement of our
      mostly-invertebrate Berghia candidates),
  (c) cover as many distinct subfamilies as possible within the quota
      (so a 5-rep aminergic selection samples 5HT + dopamine + NE +
      histamine + tyramine rather than 5 5HTs).

Output: subset FASTA + TSV ready for MAFFT + IQ-TREE 3 backbone build.

Usage:
    python3 select_backbone_reps.py \\
        --reference-fasta references/non_chemo_gpcr/all_references.fasta \\
        --reference-tsv references/non_chemo_gpcr/all_references.tsv \\
        --output-fasta results/classification/trees/backbone.fasta \\
        --output-tsv results/classification/trees/backbone.tsv \\
        --quota-per-family 6
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

INVERTEBRATE_TAXA = {"Drosophila melanogaster", "Caenorhabditis elegans"}


def load_records_from_tsv(path: str) -> list[dict]:
    """Read the curated reference TSV. Returns list of dicts with the
    fields used by the selection logic."""
    rows: list[dict] = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append({
                "accession": row.get("accession", "").strip(),
                "family": row.get("family", "").strip(),
                "subfamily": row.get("subfamily", "").strip(),
                "species": row.get("species", "").strip(),
                "gene": row.get("gene", "").strip(),
            })
    return rows


def _select_within_family(records: list[dict], quota: int) -> list[dict]:
    """Pick `quota` reps from `records` (all in the same family). Strategy:

    1. Always include >=1 invertebrate ref if the family has one.
    2. Round-robin through subfamilies to maximize subfamily coverage.
    3. If still below quota, fill with whatever's left (mammalian preference
       just to be reproducible — alphabetical accession order).
    """
    if len(records) <= quota:
        return list(records)

    selected: list[dict] = []

    # Step 1: pick one invertebrate (preference: Drosophila over C. elegans
    # for tree-rooting consistency; alphabetical accession within taxon).
    invertebrates = sorted(
        [r for r in records if r["species"] in INVERTEBRATE_TAXA],
        key=lambda r: (r["species"], r["accession"]))
    if invertebrates:
        selected.append(invertebrates[0])

    # Step 2: round-robin through subfamilies for the remaining slots.
    remaining = [r for r in records if r not in selected]
    by_sub: dict[str, list[dict]] = defaultdict(list)
    for r in remaining:
        by_sub[r["subfamily"]].append(r)
    # Sort each subfamily list deterministically (mammalian first,
    # then alphabetical accession)
    for sub in by_sub:
        by_sub[sub].sort(
            key=lambda r: (r["species"] in INVERTEBRATE_TAXA, r["accession"]))
    sub_keys = sorted(by_sub.keys())

    while len(selected) < quota:
        added = False
        for sub in sub_keys:
            if not by_sub[sub]:
                continue
            selected.append(by_sub[sub].pop(0))
            added = True
            if len(selected) >= quota:
                break
        if not added:
            break

    return selected


def select_reps(records: list[dict], quota_per_family: int = 5
                ) -> list[dict]:
    """Apply per-family rep selection across all families."""
    by_family: dict[str, list[dict]] = defaultdict(list)
    for r in records:
        by_family[r["family"]].append(r)
    out: list[dict] = []
    for family in sorted(by_family.keys()):
        out.extend(_select_within_family(by_family[family], quota_per_family))
    return out


def write_subset_fasta(selected: list[dict], full_fasta_path: str,
                       out_path: str) -> None:
    """Read full FASTA, write only the records whose accession is in
    `selected`. Preserves header convention from the source."""
    accs_to_keep = {r["accession"] for r in selected}
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(full_fasta_path) as f, open(out_path, "w") as out:
        keep = False
        for line in f:
            if line.startswith(">"):
                # Header: >accession|family|subfamily|species
                acc = line[1:].split("|", 1)[0].strip()
                keep = acc in accs_to_keep
            if keep:
                out.write(line)


def write_subset_tsv(selected: list[dict], out_path: str) -> None:
    """Write a small per-leaf annotation TSV (used by tree visualization
    + EPA-ng placement post-processing)."""
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    cols = ["accession", "family", "subfamily", "species", "gene"]
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in selected:
            w.writerow([r.get(c, "") for c in cols])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--reference-fasta", required=True)
    ap.add_argument("--reference-tsv", required=True)
    ap.add_argument("--output-fasta", required=True)
    ap.add_argument("--output-tsv", required=True)
    ap.add_argument("--quota-per-family", type=int, default=5)
    args = ap.parse_args()

    records = load_records_from_tsv(args.reference_tsv)
    selected = select_reps(records, quota_per_family=args.quota_per_family)

    print(f"[backbone] Selected {len(selected)} reps "
          f"from {len(records)} input records:", file=sys.stderr)
    by_fam: dict[str, int] = defaultdict(int)
    for r in selected:
        by_fam[r["family"]] += 1
    for fam in sorted(by_fam):
        print(f"  {fam:<30} {by_fam[fam]:>3}", file=sys.stderr)

    write_subset_fasta(selected, args.reference_fasta, args.output_fasta)
    write_subset_tsv(selected, args.output_tsv)
    print(f"\n[backbone] Wrote {args.output_fasta} + {args.output_tsv}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
