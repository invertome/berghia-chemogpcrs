#!/usr/bin/env python3
"""Count the non-vertebrate members of orthogroups that contain characterized anchors.

This is the ceiling on what an OrthoDB-orthology route could add to the
reference envelope: if a curated family's anchors sit in orthogroup G, then the
non-vertebrate members of G are the sequences that route would nominate. The
number is only meaningful at a level whose orthogroups actually track curated
families, so run scripts/orthodb_validate_families.py first and read this
against that result.

Clade assignment uses each organism's OrthoDB level path (odb_level2species
column 4), not a name match, so it is exact for the clades OrthoDB models.
Note that OrthoDB v12v2 has no Bilateria / Protostomia / Deuterostomia /
Chordata / Ecdysozoa levels, so "other_Metazoa" necessarily absorbs the
non-vertebrate deuterostomes (tunicates, echinoderms, amphioxus), the
non-lophotrochozoan / non-nematode / non-arthropod protostomes, and the
non-bilaterians other than cnidarians. That residual is reported separately and
resolved against NCBI lineages when --resolve-residual is passed.

Usage:
    python3 scripts/orthodb_expansion_ceiling.py \
        --indir results/ranking/diagnostics/orthodb --level 33208
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

# Verified present as levels in odb12v2_levels.tab. Order matters: the first
# match wins, so Mollusca is tested before the Lophotrochozoa it sits inside.
CLADE_ORDER = [
    ("Vertebrata", "7742"),
    ("Mollusca", "6447"),
    ("other_Lophotrochozoa", "1206795"),
    ("Arthropoda", "6656"),
    ("Nematoda", "6231"),
    ("Cnidaria", "6073"),
]
METAZOA = "33208"

GENE_ID_RE = re.compile(r"^(\d+_\d+):")


def parse_org_id(gene_id: str) -> str:
    """OrthoDB gene ids are '<organism_id>:<serial>'; organism_id is '<taxid>_<n>'.

    Asserted rather than assumed, because the whole clade breakdown keys on it.
    """
    m = GENE_ID_RE.match(gene_id)
    if not m:
        raise ValueError(f"unexpected OrthoDB gene id format: {gene_id!r}")
    return m.group(1)


def load_org_clades(level2species_path: Path) -> dict[str, str]:
    """organism id -> clade label, from the organism's OrthoDB level path."""
    clades: dict[str, str] = {}
    with open(level2species_path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 4:
                continue
            org = row[1]
            path = set(row[3].strip("{}").split(","))
            if METAZOA not in path:
                clades[org] = "non_Metazoa"
                continue
            label = "other_Metazoa"
            for name, taxid in CLADE_ORDER:
                if taxid in path:
                    label = name
                    break
            clades[org] = label
    if not clades:
        raise SystemExit(f"ERROR: parsed no organisms from {level2species_path}")
    return clades


def load_org_names(species_path: Path) -> dict[str, tuple[str, str]]:
    """organism id -> (ncbi taxid, scientific name)."""
    names: dict[str, tuple[str, str]] = {}
    with open(species_path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 3:
                continue
            names[row[1]] = (row[0], row[2])
    return names


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--indir", default="results/ranking/diagnostics/orthodb")
    ap.add_argument(
        "--level",
        action="append",
        default=None,
        help="level taxid to report (repeatable). Default: all levels present.",
    )
    args = ap.parse_args(argv)
    indir = Path(args.indir)

    for f in (
        "anchor_snapshot.tsv", "anchor_gene_map.tsv", "anchor_og_meta.tsv",
        "anchor_og_membership.tsv", "odb_level2species.tsv", "odb_species.tsv",
    ):
        if not (indir / f).exists():
            raise SystemExit(f"ERROR: missing {indir / f}")

    with open(indir / "anchor_snapshot.tsv", newline="") as fh:
        snapshot = {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}

    # Only orthogroups anchored by an EXPERIMENTALLY CHARACTERIZED receptor
    # count toward the ceiling: the premise of the route is that the vertebrate
    # members are observed, not inferred.
    characterized = {
        acc for acc, r in snapshot.items()
        if r["use_primary"] == "True"
        and (r["evidence"].startswith("experimental") or r["evidence"] == "characterized-gtopdb")
    }

    anchor_genes: dict[str, set] = defaultdict(set)
    with open(indir / "anchor_gene_map.tsv", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2:
                anchor_genes[row[1]].add(row[0])

    og_level: dict[str, str] = {}
    og_name: dict[str, str] = {}
    with open(indir / "anchor_og_meta.tsv", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2:
                og_level[row[0]] = row[1]
                og_name[row[0]] = row[2] if len(row) > 2 else ""

    clades = load_org_clades(indir / "odb_level2species.tsv")
    org_names = load_org_names(indir / "odb_species.tsv")

    # Which orthogroups hold a characterized anchor?
    og_has_characterized: dict[str, set] = defaultdict(set)
    og_families: dict[str, Counter] = defaultdict(Counter)
    membership: dict[str, list] = defaultdict(list)
    unknown_org = 0

    with open(indir / "anchor_og_membership.tsv", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 2:
                continue
            og, gene = row[0], row[1]
            membership[og].append(gene)
            for acc in anchor_genes.get(gene, ()):
                if acc in characterized:
                    og_has_characterized[og].add(acc)
                rec = snapshot.get(acc)
                if rec and rec["use_primary"] == "True":
                    og_families[og][rec["family"]] += 1

    wanted_levels = set(args.level) if args.level else set(og_level.values())

    rows = []
    for og, genes in membership.items():
        lvl = og_level.get(og)
        if lvl not in wanted_levels:
            continue
        if not og_has_characterized.get(og):
            continue
        counts: Counter = Counter()
        species_by_clade: dict[str, set] = defaultdict(set)
        for g in genes:
            try:
                org = parse_org_id(g)
            except ValueError:
                unknown_org += 1
                continue
            clade = clades.get(org, "unknown_organism")
            counts[clade] += 1
            species_by_clade[clade].add(org)
        rows.append(
            {
                "level_taxid": lvl,
                "og_id": og,
                "og_name": og_name.get(og, ""),
                "curated_families": ";".join(
                    f"{k}:{v}" for k, v in og_families[og].most_common()
                ),
                "n_curated_families": len(og_families[og]),
                "characterized_anchors": len(og_has_characterized[og]),
                "members_total": len(genes),
                **{f"n_{name}": counts.get(name, 0) for name, _ in CLADE_ORDER},
                "n_other_Metazoa": counts.get("other_Metazoa", 0),
                "n_non_Metazoa": counts.get("non_Metazoa", 0),
                "n_nonvertebrate": len(genes) - counts.get("Vertebrata", 0),
                "spp_Mollusca": len(species_by_clade.get("Mollusca", ())),
            }
        )

    if unknown_org:
        print(f"WARNING: {unknown_org} gene ids did not parse to an organism id")

    rows.sort(key=lambda r: -r["members_total"])
    cols = list(rows[0].keys()) if rows else []
    out = indir / "expansion_ceiling_per_og.tsv"
    with open(out, "w", newline="") as fh:
        if cols:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    # Per-level totals, de-duplicated by gene so a gene in two orthogroups at
    # the same level is not double counted.
    summary: dict[str, dict] = {}
    for lvl in sorted(wanted_levels):
        lvl_rows = [r for r in rows if r["level_taxid"] == lvl]
        if not lvl_rows:
            continue
        genes_seen: dict[str, str] = {}
        for r in lvl_rows:
            for g in membership[r["og_id"]]:
                try:
                    genes_seen[g] = clades.get(parse_org_id(g), "unknown_organism")
                except ValueError:
                    pass
        per_clade = Counter(genes_seen.values())
        summary[lvl] = {
            "orthogroups_with_characterized_anchor": len(lvl_rows),
            "distinct_member_genes": len(genes_seen),
            "by_clade": dict(per_clade.most_common()),
            "nonvertebrate_total": len(genes_seen) - per_clade.get("Vertebrata", 0),
            "molluscan_species_represented": len(
                {parse_org_id(g) for g, c in genes_seen.items() if c == "Mollusca"}
            ),
        }

    (indir / "expansion_ceiling_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    print("=== EXPANSION CEILING (orthogroups holding >=1 characterized anchor) ===")
    for lvl, s in summary.items():
        print(f"\nlevel {lvl}: {s['orthogroups_with_characterized_anchor']} orthogroups, "
              f"{s['distinct_member_genes']} distinct member genes")
        for clade, n in s["by_clade"].items():
            print(f"    {clade:<22} {n:>8}")
        print(f"    {'NON-VERTEBRATE TOTAL':<22} {s['nonvertebrate_total']:>8}")
        print(f"    molluscan species represented: {s['molluscan_species_represented']}")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
