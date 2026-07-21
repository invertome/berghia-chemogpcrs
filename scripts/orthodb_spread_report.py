#!/usr/bin/env python3
"""Report an anchor set's taxonomic spread: distinct taxa, never sequence counts.

Sequence count is the wrong measure of a reference envelope. Thirty-one
chemokine sequences drawn from three species in one phylum describe vertebrate
chemokine idiosyncrasy, not the family's variation across Metazoa, and no
amount of additional vertebrate sequences fixes that. So every figure here
counts DISTINCT species, orders and phyla.

Ranked lineages are resolved from NCBI and cached; a taxid that will not
resolve is reported rather than dropped, because a silently missing taxon
inflates nothing but hides a broken join.

Usage:
    python3 scripts/orthodb_spread_report.py \
        --anchors references/anchors/anchor_set_PROD.tsv \
        --cache results/ranking/diagnostics/orthodb/ranked_lineage.tsv
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
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BATCH = 150
RANKS = ["phylum", "class", "order", "family", "genus"]


def fetch_ranked(taxids: list[str], cache: Path) -> dict[str, dict]:
    known: dict[str, dict] = {}
    if cache.exists():
        with open(cache, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                known[r["taxid"]] = r
    todo = sorted({t for t in taxids if t and t not in known})
    if todo:
        print(f"fetching ranked lineages for {len(todo)} taxids ...", file=sys.stderr)
    for i in range(0, len(todo), BATCH):
        chunk = todo[i:i + BATCH]
        data = urllib.parse.urlencode(
            {"db": "taxonomy", "id": ",".join(chunk), "retmode": "xml"}).encode()
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
            if not tid:
                continue
            rec = {"taxid": tid, "name": taxon.findtext("ScientificName") or ""}
            for rk in RANKS:
                rec[rk] = ""
            for t in taxon.findall("./LineageEx/Taxon"):
                rk = (t.findtext("Rank") or "").lower()
                if rk in RANKS:
                    rec[rk] = t.findtext("ScientificName") or ""
            own = (taxon.findtext("Rank") or "").lower()
            if own in RANKS and not rec[own]:
                rec[own] = rec["name"]
            known[tid] = rec
            for aka in taxon.findall("./AkaTaxIds/TaxId"):
                if aka.text:
                    known[aka.text] = dict(rec, taxid=aka.text)
        time.sleep(0.4)
    with open(cache, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["taxid", "name"] + RANKS, delimiter="\t")
        w.writeheader()
        for t in sorted(known):
            w.writerow({k: known[t].get(k, "") for k in ["taxid", "name"] + RANKS})
    return known


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--anchors", required=True)
    ap.add_argument("--cache",
                    default="results/ranking/diagnostics/orthodb/ranked_lineage.tsv")
    ap.add_argument("--label", default="")
    ap.add_argument("--out-json")
    args = ap.parse_args(argv)

    with open(args.anchors, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    tax = fetch_ranked([r["taxid"] for r in rows], Path(args.cache))

    unresolved = sorted({r["taxid"] for r in rows if r["taxid"] not in tax})
    if unresolved:
        print(f"WARNING: {len(unresolved)} taxid(s) did not resolve "
              f"(e.g. {unresolved[:5]}) -- they contribute no rank counts",
              file=sys.stderr)

    inv: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))
    for r in rows:
        t = tax.get(r["taxid"], {})
        for scope in (r["family"], "ALL"):
            for rk in ("phylum", "class", "order"):
                if t.get(rk):
                    inv[scope][rk].add(t[rk])
            inv[scope]["species"].add(r["taxid"])

    counts = Counter(r["family"] for r in rows)
    title = f"SPREAD INVENTORY{' -- ' + args.label if args.label else ''}"
    print(f"=== {title} ({args.anchors}) ===")
    hdr = f"{'family':<24}{'seqs':>6}{'species':>9}{'orders':>8}{'classes':>9}{'phyla':>7}"
    print(hdr)
    print("-" * len(hdr))
    for fam in sorted(counts, key=lambda f: -counts[f]):
        d = inv[fam]
        print(f"{fam:<24}{counts[fam]:>6}{len(d['species']):>9}"
              f"{len(d['order']):>8}{len(d['class']):>9}{len(d['phylum']):>7}")
    d = inv["ALL"]
    print(f"{'ALL':<24}{len(rows):>6}{len(d['species']):>9}"
          f"{len(d['order']):>8}{len(d['class']):>9}{len(d['phylum']):>7}")

    print(f"\nphyla present ({len(inv['ALL']['phylum'])}): "
          f"{sorted(inv['ALL']['phylum'])}")

    summary = {
        "anchors": args.anchors,
        "rows": len(rows),
        "species": len(inv["ALL"]["species"]),
        "orders": len(inv["ALL"]["order"]),
        "classes": len(inv["ALL"]["class"]),
        "phyla": len(inv["ALL"]["phylum"]),
        "phyla_list": sorted(inv["ALL"]["phylum"]),
        "per_family": {
            f: {"seqs": counts[f], "species": len(inv[f]["species"]),
                "orders": len(inv[f]["order"]), "classes": len(inv[f]["class"]),
                "phyla": len(inv[f]["phylum"])}
            for f in sorted(counts)
        },
        "unresolved_taxids": unresolved,
    }
    if args.out_json:
        Path(args.out_json).write_text(json.dumps(summary, indent=2) + "\n")
        print(f"\nwrote {args.out_json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
