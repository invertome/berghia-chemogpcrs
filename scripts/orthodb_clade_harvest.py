#!/usr/bin/env python3
"""Per-clade harvest analysis of anchor-seeded OrthoDB orthogroups.

An earlier pass measured molluscan content and concluded the harvest was
unusable. That conclusion does not follow: the reference envelope's defect is
that it is VERTEBRATE-dominated, so its within-family spread encodes vertebrate
idiosyncrasy. Any well-characterized invertebrate widens it toward genuine
family variation across Metazoa. A cnidarian or an echinoderm does that job
without being a mollusc, and the anchor set currently holds zero of either.

So this script asks, per clade rather than for molluscs: which orthogroups are
tight enough for a harvest to be family-specific, and what do those tight
groups actually contain?

Clade assignment uses NCBI taxonomy lineages for OrthoDB's own organisms, NOT
OrthoDB's level tree. That is deliberate: OrthoDB has no Deuterostomia,
Ecdysozoa or Bilateria level, so its tree literally cannot express
"non-vertebrate deuterostome" - the single clade with the most headroom. NCBI's
taxonomy can. Lineages are cached with provenance.

Usage:
    python3 scripts/orthodb_clade_harvest.py --indir results/ranking/diagnostics/orthodb
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BATCH = 150
GENE_ID_RE = re.compile(r"^(\d+_\d+):")

# NCBI taxonomy ids, first match in the lineage wins. Finer than OrthoDB's
# level tree: Deuterostomia and Ecdysozoa exist here and do not exist there.
CLADE_RULES = [
    ("Vertebrata", 7742),
    ("Mollusca", 6447),
    ("other_Lophotrochozoa", 1206795),
    ("Arthropoda", 6656),
    ("Nematoda", 6231),
    ("other_Ecdysozoa", 1206794),
    ("non_vertebrate_Deuterostomia", 33511),
    ("Cnidaria", 6073),
    ("other_non_Bilateria", 33208),
]

# Clades with zero representation in the current anchor set carry the most
# headroom; adding one there changes the envelope more than another arthropod.
ANCHOR_CLADE_COUNTS_PRIMARY = {
    "Vertebrata": 399,
    "Arthropoda": 0,   # filled at runtime from the audit
}


def fetch_lineages(taxids: list[str], cache: Path) -> dict[str, list[int]]:
    known: dict[str, list[int]] = {}
    if cache.exists():
        with open(cache, newline="") as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if len(row) >= 2:
                    known[row[0]] = [int(x) for x in row[1].split(",") if x]
    todo = sorted({t for t in taxids if t and t not in known})
    if todo:
        print(f"fetching NCBI lineages for {len(todo)} organisms ...", file=sys.stderr)
    for i in range(0, len(todo), BATCH):
        chunk = todo[i : i + BATCH]
        data = urllib.parse.urlencode(
            {"db": "taxonomy", "id": ",".join(chunk), "retmode": "xml"}
        ).encode()
        for attempt in range(4):
            try:
                with urllib.request.urlopen(EUTILS, data=data, timeout=180) as fh:
                    root = ET.parse(fh).getroot()
                break
            except Exception as exc:
                if attempt == 3:
                    raise SystemExit(f"ERROR: NCBI taxonomy fetch failed: {exc}")
                time.sleep(3 * (attempt + 1))
        for taxon in root.findall("Taxon"):
            tid = taxon.findtext("TaxId")
            lin = [
                int(t.findtext("TaxId"))
                for t in taxon.findall("./LineageEx/Taxon")
                if t.findtext("TaxId")
            ]
            if tid:
                lin.append(int(tid))
                known[tid] = lin
                for aka in taxon.findall("./AkaTaxIds/TaxId"):
                    if aka.text:
                        known[aka.text] = lin
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
    ap.add_argument("--level", default="33208")
    ap.add_argument(
        "--tight-max", type=int, default=2500,
        help="orthogroups at or below this member count are treated as tight "
             "enough for a harvest to be plausibly family-specific",
    )
    args = ap.parse_args(argv)
    indir = Path(args.indir)

    # OrthoDB organism id -> ncbi taxid, name
    org_taxid: dict[str, str] = {}
    org_name: dict[str, str] = {}
    with open(indir / "odb_species.tsv", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 3:
                org_taxid[row[1]] = row[0]
                org_name[row[1]] = row[2]

    # anchor-seeded orthogroups at the chosen level
    og_level: dict[str, str] = {}
    og_name: dict[str, str] = {}
    with open(indir / "anchor_og_meta.tsv", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2:
                og_level[row[0]] = row[1]
                og_name[row[0]] = row[2] if len(row) > 2 else ""

    # Only organisms that actually occur in the orthogroups under test need a
    # lineage; fetching OrthoDB's full 32k-organism catalogue would hammer
    # NCBI for records we never look at.
    needed_orgs: set[str] = set()
    with open(indir / "anchor_og_membership.tsv", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2 and og_level.get(row[0]) == args.level:
                m = GENE_ID_RE.match(row[1])
                if m:
                    needed_orgs.add(m.group(1))
    missing_org = needed_orgs - set(org_taxid)
    if missing_org:
        raise SystemExit(
            f"ERROR: {len(missing_org)} organism ids in the membership table are "
            f"absent from the OrthoDB species table (e.g. {sorted(missing_org)[:5]}). "
            "Refusing to bin members into clades on an incomplete organism map."
        )
    lineages = fetch_lineages(
        [org_taxid[o] for o in sorted(needed_orgs)],
        indir / "orthodb_organism_lineage.tsv",
    )
    org_clade = {
        org: clade_of(lineages.get(org_taxid[org], [])) for org in needed_orgs
    }
    unresolved = [o for o, c in org_clade.items() if c == "unknown"]
    if unresolved:
        print(
            f"WARNING: {len(unresolved)} organisms have no NCBI lineage and are "
            f"binned as 'unknown' (e.g. {unresolved[:5]})",
            file=sys.stderr,
        )

    with open(indir / "anchor_snapshot.tsv", newline="") as fh:
        snapshot = {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}
    characterized = {
        a for a, r in snapshot.items()
        if r["use_primary"] == "True"
        and (r["evidence"].startswith("experimental")
             or r["evidence"] == "characterized-gtopdb")
    }
    gene_to_accs: dict[str, set] = defaultdict(set)
    with open(indir / "anchor_gene_map.tsv", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2:
                gene_to_accs[row[1]].add(row[0])

    og_members: dict[str, list[str]] = defaultdict(list)
    with open(indir / "anchor_og_membership.tsv", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2 and og_level.get(row[0]) == args.level:
                og_members[row[0]].append(row[1])

    og_seeded: dict[str, set] = defaultdict(set)
    og_fams: dict[str, Counter] = defaultdict(Counter)
    for og, genes in og_members.items():
        for g in genes:
            for acc in gene_to_accs.get(g, ()):
                if acc in characterized:
                    og_seeded[og].add(acc)
                r = snapshot.get(acc)
                if r and r["use_primary"] == "True":
                    og_fams[og][r["family"]] += 1

    clades = [c for c, _ in CLADE_RULES] + ["non_Metazoa", "unknown"]
    rows = []
    for og, genes in og_members.items():
        if not og_seeded.get(og):
            continue
        cc: Counter = Counter()
        spp: dict[str, set] = defaultdict(set)
        for g in genes:
            m = GENE_ID_RE.match(g)
            if not m:
                continue
            org = m.group(1)
            c = org_clade.get(org, "unknown")
            cc[c] += 1
            spp[c].add(org)
        rows.append({
            "og_id": og,
            "og_name": og_name.get(og, ""),
            "members_total": len(genes),
            "characterized_anchors": len(og_seeded[og]),
            "curated_families": ";".join(f"{k}:{v}" for k, v in og_fams[og].most_common()),
            "n_curated_families": len(og_fams[og]),
            **{f"n_{c}": cc.get(c, 0) for c in clades},
            **{f"spp_{c}": len(spp.get(c, ())) for c in clades},
        })

    rows.sort(key=lambda r: r["members_total"])
    out = indir / "clade_harvest_per_og.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    report_clades = [c for c, _ in CLADE_RULES if c != "Vertebrata"]
    totals = {c: sum(r[f"n_{c}"] for r in rows) for c in report_clades}

    print(f"=== ANCHOR-SEEDED ORTHOGROUPS AT LEVEL {args.level}: {len(rows)} ===\n")
    print("=== member counts per clade, by orthogroup size band ===")
    bands = [(0, 1000), (1000, 2500), (2500, 5000), (5000, 10000), (10000, 10**9)]
    hdr = f"{'size band':<14}{'OGs':>5}" + "".join(f"{c[:13]:>14}" for c in report_clades)
    print(hdr); print("-" * len(hdr))
    for lo, hi in bands:
        sel = [r for r in rows if lo <= r["members_total"] < hi]
        lbl = f"{lo}-{hi if hi < 10**9 else '+'}"
        print(f"{lbl:<14}{len(sel):>5}" +
              "".join(f"{sum(r[f'n_{c}'] for r in sel):>14}" for c in report_clades))
    print(f"{'TOTAL':<14}{len(rows):>5}" +
          "".join(f"{totals[c]:>14}" for c in report_clades))

    print(f"\n=== share of each clade's members sitting in TIGHT OGs (<= {args.tight_max}) ===")
    tight = [r for r in rows if r["members_total"] <= args.tight_max]
    print(f"{'clade':<32}{'total':>10}{'in tight OGs':>14}{'% tight':>10}{'tight OGs w/ clade':>20}")
    for c in report_clades:
        t = totals[c]
        n = sum(r[f"n_{c}"] for r in tight)
        k = sum(1 for r in tight if r[f"n_{c}"] > 0)
        print(f"{c:<32}{t:>10}{n:>14}{(n/t*100 if t else 0):>9.1f}%{k:>20}")

    print(f"\n=== the {len(tight)} tight orthogroups (<= {args.tight_max} members) ===")
    print(f"{'og_name':<38}{'memb':>6}" + "".join(f"{c[:9]:>10}" for c in report_clades))
    for r in tight:
        print(f"{r['og_name'][:37]:<38}{r['members_total']:>6}" +
              "".join(f"{r[f'n_{c}']:>10}" for c in report_clades))

    summary = {
        "level": args.level,
        "orthogroups_seeded": len(rows),
        "tight_max": args.tight_max,
        "tight_orthogroups": len(tight),
        "clade_totals": totals,
        "clade_totals_tight": {c: sum(r[f"n_{c}"] for r in tight) for c in report_clades},
        "species_represented": {
            c: len({r[f"spp_{c}"] for r in rows if r[f"spp_{c}"]}) and
               max(r[f"spp_{c}"] for r in rows) for c in report_clades
        },
    }
    (indir / "clade_harvest_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
