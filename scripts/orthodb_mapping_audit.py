#!/usr/bin/env python3
"""Characterise WHICH anchors fail to map into OrthoDB, by clade and by family.

The whole point of testing OrthoDB is to find structured evidence that can
judge a MOLLUSCAN query. The reference envelope is already ~89% vertebrate and
arthropod; if the anchors that fail to map into OrthoDB are disproportionately
the non-vertebrate ones, then every conflation/fragmentation statistic is
measured on the vertebrate-heavy subset and says nothing about the sequences we
actually care about. That is a negative result about the method, not a caveat
on the numbers, so it is computed and reported before any agreement statistic.

The failures are also split into two kinds, because they have different
remedies:

  COVERAGE gap    the anchor's species is absent from OrthoDB's organism
                  catalogue entirely. A hard ceiling: no join key fixes it.
  IDENTIFIER gap  the species IS in OrthoDB but the UniProt accession carries
                  no cross-reference to an OrthoDB gene. Potentially fixable
                  with a different join key (sequence search, NCBI/Ensembl id).

Clade assignment uses NCBI taxonomy lineages fetched from E-utilities for the
anchors' taxids, NOT genus string matching and NOT OrthoDB's level tree.
OrthoDB has no Deuterostomia or Ecdysozoa level, so its tree cannot express the
breakdown the question asks for; NCBI's can. Lineages are cached to disk with
provenance so the assignment is reproducible.

Usage:
    python3 scripts/orthodb_mapping_audit.py --indir results/ranking/diagnostics/orthodb
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BATCH = 150

# NCBI taxonomy ids. Order matters: first match in the lineage wins, so nested
# clades are listed before the clades that contain them.
CLADE_RULES = [
    ("Vertebrata", 7742),
    ("Mollusca", 6447),
    ("other_Lophotrochozoa", 1206795),
    ("Ecdysozoa", 1206794),
    ("non_vertebrate_Deuterostomia", 33511),
    ("Cnidaria", 6073),
    ("other_Metazoa", 33208),
]


def fetch_lineages(taxids: list[str], cache: Path) -> dict[str, list[int]]:
    """taxid -> list of ancestor taxids, cached on disk."""
    known: dict[str, list[int]] = {}
    if cache.exists():
        with open(cache, newline="") as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if len(row) >= 2:
                    known[row[0]] = [int(x) for x in row[1].split(",") if x]

    todo = sorted({t for t in taxids if t and t not in known})
    if todo:
        print(f"fetching NCBI lineages for {len(todo)} taxids ...", file=sys.stderr)
    for i in range(0, len(todo), BATCH):
        chunk = todo[i : i + BATCH]
        data = urllib.parse.urlencode(
            {"db": "taxonomy", "id": ",".join(chunk), "retmode": "xml"}
        ).encode()
        for attempt in range(4):
            try:
                with urllib.request.urlopen(EUTILS, data=data, timeout=120) as fh:
                    root = ET.parse(fh).getroot()
                break
            except Exception as exc:  # transient NCBI failures are routine
                if attempt == 3:
                    raise SystemExit(f"ERROR: NCBI taxonomy fetch failed: {exc}")
                time.sleep(3 * (attempt + 1))
        for taxon in root.findall("Taxon"):
            tid = taxon.findtext("TaxId")
            lineage = [
                int(t.findtext("TaxId"))
                for t in taxon.findall("./LineageEx/Taxon")
                if t.findtext("TaxId")
            ]
            if tid:
                lineage.append(int(tid))
                known[tid] = lineage
                # AkaTaxIds: the queried id may have been merged into another.
                for aka in taxon.findall("./AkaTaxIds/TaxId"):
                    if aka.text:
                        known[aka.text] = lineage
        time.sleep(0.4)

    with open(cache, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        for t, lin in sorted(known.items()):
            w.writerow([t, ",".join(str(x) for x in lin)])
    return known


def clade_of(lineage: list[int]) -> str:
    s = set(lineage)
    for name, taxid in CLADE_RULES:
        if taxid in s:
            return name
    return "non_Metazoa" if lineage else "unknown"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--indir", default="results/ranking/diagnostics/orthodb")
    args = ap.parse_args(argv)
    indir = Path(args.indir)

    with open(indir / "anchor_snapshot.tsv", newline="") as fh:
        snapshot = list(csv.DictReader(fh, delimiter="\t"))

    gene_map_path = indir / "anchor_gene_map.tsv"
    if not gene_map_path.exists():
        raise SystemExit(f"ERROR: {gene_map_path} missing; run the Unity join first.")
    mapped: set[str] = set()
    with open(gene_map_path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if row:
                mapped.add(row[0])

    # OrthoDB organism catalogue, for coverage-vs-identifier attribution.
    odb_taxids: set[int] = set()
    sp_path = indir / "odb_species.tsv"
    if sp_path.exists():
        with open(sp_path, newline="") as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if row and row[0].isdigit():
                    odb_taxids.add(int(row[0]))
    else:
        print(f"WARNING: {sp_path} missing; coverage attribution skipped", file=sys.stderr)

    lineages = fetch_lineages(
        [r["taxid"] for r in snapshot], indir / "anchor_taxid_lineage.tsv"
    )
    missing_lineage = [r["accession"] for r in snapshot if r["taxid"] not in lineages]
    if missing_lineage:
        print(
            f"WARNING: {len(missing_lineage)} anchors have no NCBI lineage "
            f"(e.g. {missing_lineage[:5]}); they are binned as 'unknown'",
            file=sys.stderr,
        )

    rows = []
    for r in snapshot:
        lin = lineages.get(r["taxid"], [])
        in_odb = bool(set(lin) & odb_taxids) if odb_taxids else None
        is_mapped = r["accession"] in mapped
        if is_mapped:
            gap = ""
        elif in_odb is None:
            gap = "unattributed"
        else:
            gap = "identifier_gap" if in_odb else "coverage_gap"
        rows.append(
            {
                **r,
                "clade": clade_of(lin),
                "species_in_orthodb": "" if in_odb is None else str(in_odb),
                "mapped": str(is_mapped),
                "gap_type": gap,
            }
        )

    out = indir / "mapping_audit.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    primary = [r for r in rows if r["use_primary"] == "True"]

    def table(recs, key):
        agg: dict[str, dict] = defaultdict(
            lambda: {"total": 0, "mapped": 0, "coverage_gap": 0, "identifier_gap": 0}
        )
        for r in recs:
            a = agg[r[key]]
            a["total"] += 1
            if r["mapped"] == "True":
                a["mapped"] += 1
            elif r["gap_type"] in a:
                a[r["gap_type"]] += 1
        for k, a in agg.items():
            a["rate"] = a["mapped"] / a["total"] if a["total"] else 0.0
        return dict(sorted(agg.items(), key=lambda kv: -kv[1]["total"]))

    by_clade = table(primary, "clade")
    by_family = table(primary, "family")
    overall = {
        "anchors_class_a": len(rows),
        "anchors_primary": len(primary),
        "primary_mapped": sum(1 for r in primary if r["mapped"] == "True"),
        "primary_unmapped": sum(1 for r in primary if r["mapped"] != "True"),
        "primary_mapping_rate": (
            sum(1 for r in primary if r["mapped"] == "True") / len(primary)
            if primary else 0.0
        ),
        "unmapped_coverage_gap": sum(1 for r in primary if r["gap_type"] == "coverage_gap"),
        "unmapped_identifier_gap": sum(
            1 for r in primary if r["gap_type"] == "identifier_gap"
        ),
    }

    # Skew: how much does mapping shift the vertebrate share of the set?
    vert_before = sum(1 for r in primary if r["clade"] == "Vertebrata") / len(primary)
    mapped_primary = [r for r in primary if r["mapped"] == "True"]
    vert_after = (
        sum(1 for r in mapped_primary if r["clade"] == "Vertebrata") / len(mapped_primary)
        if mapped_primary else 0.0
    )
    skew = {
        "vertebrate_share_before_mapping": round(vert_before, 4),
        "vertebrate_share_after_mapping": round(vert_after, 4),
        "vertebrate_share_shift_pp": round((vert_after - vert_before) * 100, 2),
        "nonvertebrate_before": sum(1 for r in primary if r["clade"] != "Vertebrata"),
        "nonvertebrate_after": sum(
            1 for r in mapped_primary if r["clade"] != "Vertebrata"
        ),
    }

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "overall": overall,
        "by_clade": by_clade,
        "by_family": by_family,
        "representativeness": skew,
    }
    (indir / "mapping_audit_summary.json").write_text(json.dumps(report, indent=2) + "\n")

    print("=== MAPPING RATE BY CLADE (primary anchors) ===")
    hdr = f"{'clade':<30}{'total':>7}{'mapped':>8}{'rate':>8}{'cover':>7}{'ident':>7}"
    print(hdr); print("-" * len(hdr))
    for k, a in by_clade.items():
        print(f"{k:<30}{a['total']:>7}{a['mapped']:>8}{a['rate']*100:>7.1f}%"
              f"{a['coverage_gap']:>7}{a['identifier_gap']:>7}")
    print(f"\n{'ALL':<30}{overall['anchors_primary']:>7}{overall['primary_mapped']:>8}"
          f"{overall['primary_mapping_rate']*100:>7.1f}%"
          f"{overall['unmapped_coverage_gap']:>7}{overall['unmapped_identifier_gap']:>7}")

    print("\n=== MAPPING RATE BY CURATED FAMILY (primary anchors) ===")
    print(hdr); print("-" * len(hdr))
    for k, a in by_family.items():
        print(f"{k:<30}{a['total']:>7}{a['mapped']:>8}{a['rate']*100:>7.1f}%"
              f"{a['coverage_gap']:>7}{a['identifier_gap']:>7}")

    print("\n=== REPRESENTATIVENESS OF THE MAPPED SUBSET ===")
    print(json.dumps(skew, indent=2))
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
