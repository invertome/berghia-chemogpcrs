#!/usr/bin/env python3
"""Measure the anchor set's taxonomic spread, then harvest to widen it.

The objective is PHYLOGENETIC DIVERSITY, not sequence count. Thirty-eight
percent of anchors being Arthropoda sounds like coverage until you look at
which arthropods: a handful of model organisms (Drosophila, Aedes, Bombyx,
Rhodnius, Ixodes, Tribolium, Anopheles) spanning three dipterans, one
lepidopteran, one hemipteran, one beetle and a single tick - out of Crustacea,
Myriapoda, the rest of Chelicerata and a dozen unrepresented insect orders. For
widening a family's envelope that is nearly as narrow as having none. So the
inventory counts DISTINCT TAXA at phylum / class / order, never sequences, and
the harvest is selected to maximise distinct orders added per sequence.

Selection metric, stated explicitly: greedy maximisation of newly-added
distinct ORDERS, per receptor family. Order is the operational rank because it
is the finest rank at which NCBI's classification is populated consistently
across Metazoa (many invertebrate lineages have no assigned class, and family
rank is too sparse in non-model clades). Ties are broken by newly-added class,
then phylum, then by preferring an organism with fewer sequences already taken.
Per-taxon redundancy is capped so a clade OrthoDB happens to hold a lot of
cannot swamp the result - that would deepen a different bias rather than fix
the original one.

The granularity gate applied here is the one DERIVED by
scripts/orthodb_gate_loo.py, not a chosen one: strict single-family purity of
the characterized seed anchors, and a minimum seed count. Leave-one-out on the
198 invertebrate anchors showed a size cap destroys recovery without improving
precision, so size is reported as a metric but is NOT a gate.

Usage:
    python3 scripts/orthodb_spread_and_harvest.py \
        --indir results/ranking/diagnostics/orthodb --min-seeds 5
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


def is_characterized(rec: dict) -> bool:
    ev = rec.get("evidence", "")
    return ev.startswith("experimental") or ev == "characterized-gtopdb"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--indir", default="results/ranking/diagnostics/orthodb")
    ap.add_argument("--level", default="33208")
    ap.add_argument("--min-seeds", type=int, default=5,
                    help="derived gate: minimum characterized same-family seeds")
    ap.add_argument("--max-per-order", type=int, default=3,
                    help="redundancy cap: sequences taken per (family, order)")
    ap.add_argument("--max-per-species", type=int, default=1,
                    help="redundancy cap: sequences taken per (family, species)")
    args = ap.parse_args(argv)
    indir = Path(args.indir)

    with open(indir / "mapping_audit.tsv", newline="") as fh:
        audit = list(csv.DictReader(fh, delimiter="\t"))
    primary = [r for r in audit if r["use_primary"] == "True"]

    org_taxid, org_name = {}, {}
    with open(indir / "odb_species.tsv", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 3:
                org_taxid[row[1]] = row[0]
                org_name[row[1]] = row[2]

    og_level, og_name = {}, {}
    with open(indir / "anchor_og_meta.tsv", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2:
                og_level[row[0]] = row[1]
                og_name[row[0]] = row[2] if len(row) > 2 else ""

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

    by_acc = {r["accession"]: r for r in audit}
    needed = [r["taxid"] for r in primary]
    member_orgs = {
        GENE_ID_RE.match(g).group(1)
        for genes in og_members.values() for g in genes if GENE_ID_RE.match(g)
    }
    needed += [org_taxid[o] for o in member_orgs if o in org_taxid]
    tax = fetch_ranked(needed, indir / "ranked_lineage.tsv")

    # ---------------------------------------------------- spread inventory --
    print("=== SPREAD INVENTORY OF THE CURRENT ANCHOR SET (distinct taxa, not counts) ===")
    inv: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))
    for r in primary:
        t = tax.get(r["taxid"], {})
        for rk in ("phylum", "class", "order"):
            if t.get(rk):
                inv[r["family"]][rk].add(t[rk])
        inv[r["family"]]["species"].add(r["taxid"])
        for rk in ("phylum", "class", "order"):
            if t.get(rk):
                inv["ALL"][rk].add(t[rk])
        inv["ALL"]["species"].add(r["taxid"])
    hdr = f"{'family':<24}{'seqs':>6}{'species':>9}{'orders':>8}{'classes':>9}{'phyla':>7}"
    print(hdr); print("-" * len(hdr))
    fam_counts = Counter(r["family"] for r in primary)
    for fam in sorted(inv, key=lambda f: -fam_counts.get(f, 0)):
        if fam == "ALL":
            continue
        d = inv[fam]
        print(f"{fam:<24}{fam_counts[fam]:>6}{len(d['species']):>9}"
              f"{len(d['order']):>8}{len(d['class']):>9}{len(d['phylum']):>7}")
    d = inv["ALL"]
    print(f"{'ALL':<24}{len(primary):>6}{len(d['species']):>9}"
          f"{len(d['order']):>8}{len(d['class']):>9}{len(d['phylum']):>7}")

    print("\n=== arthropod orders present in the anchor set ===")
    arth = {tax[r['taxid']].get('order','') for r in primary
            if tax.get(r['taxid'],{}).get('phylum') == 'Arthropoda'}
    print(sorted(x for x in arth if x))
    print("\n=== phyla present in the anchor set ===")
    print(sorted({tax[r['taxid']].get('phylum','') for r in primary
                  if tax.get(r['taxid'],{}).get('phylum')}))

    # ------------------------------------------------------- apply the gate --
    og_seeds: dict[str, Counter] = defaultdict(Counter)
    for og, genes in og_members.items():
        for g in genes:
            for acc in gene_to_accs.get(g, ()):
                r = by_acc.get(acc)
                if r and r["use_primary"] == "True" and is_characterized(r):
                    og_seeds[og][r["family"]] += 1

    passed = {}
    for og, seeds in og_seeds.items():
        if len(seeds) != 1:
            continue
        fam, n = next(iter(seeds.items()))
        if n < args.min_seeds:
            continue
        passed[og] = fam
    print(f"\n=== GATE (strict purity, min_seeds={args.min_seeds}, no size cap) ===")
    print(f"orthogroups seeded: {len(og_seeds)}   passing the gate: {len(passed)}")
    for og, fam in sorted(passed.items(), key=lambda kv: -og_seeds[kv[0]][kv[1]]):
        print(f"  {og_name.get(og,'')[:44]:<46} {og_seeds[og][fam]:>3} seeds  "
              f"{len(og_members[og]):>6} members  family={fam}")

    # ----------------------------------------- greedy diversity-first harvest --
    anchor_orders = {rk: set(inv["ALL"][rk]) for rk in ("phylum", "class", "order")}
    fam_have: dict[str, dict[str, set]] = {
        f: {rk: set(inv[f][rk]) for rk in ("phylum", "class", "order")} for f in inv
    }

    cands = []
    for og, fam in passed.items():
        for g in og_members[og]:
            m = GENE_ID_RE.match(g)
            if not m:
                continue
            org = m.group(1)
            tx = org_taxid.get(org)
            t = tax.get(tx, {})
            if t.get("phylum") == "Chordata" and t.get("class") in (
                "Mammalia", "Aves", "Actinopteri", "Amphibia", "Lepidosauria"):
                continue  # vertebrates do not widen the envelope
            cands.append({
                "og_id": og, "family": fam, "gene_id": g, "org_id": org,
                "taxid": tx or "", "species": org_name.get(org, ""),
                "phylum": t.get("phylum", ""), "class": t.get("class", ""),
                "order": t.get("order", ""),
                "og_members": len(og_members[og]),
                "og_seeds": og_seeds[og][fam],
                "og_name": og_name.get(og, ""),
            })

    taken, per_order, per_species = [], Counter(), Counter()
    added: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))

    def gain(c):
        f = c["family"]
        have = fam_have.get(f, {"phylum": set(), "class": set(), "order": set()})
        g_ord = 1 if c["order"] and c["order"] not in have["order"] | added[f]["order"] else 0
        g_cls = 1 if c["class"] and c["class"] not in have["class"] | added[f]["class"] else 0
        g_phy = 1 if c["phylum"] and c["phylum"] not in have["phylum"] | added[f]["phylum"] else 0
        return (g_ord, g_cls, g_phy)

    pool = list(cands)
    while pool:
        pool = [c for c in pool
                if per_order[(c["family"], c["order"])] < args.max_per_order
                and per_species[(c["family"], c["taxid"])] < args.max_per_species]
        if not pool:
            break
        scored = [(gain(c), -per_species[(c["family"], c["taxid"])], c) for c in pool]
        best = max(scored, key=lambda x: (sum(x[0]), x[0][0], x[0][1], x[0][2], x[1]))
        if sum(best[0]) == 0:
            break
        c = best[2]
        taken.append({**c, "gain_order": best[0][0], "gain_class": best[0][1],
                      "gain_phylum": best[0][2]})
        f = c["family"]
        for rk in ("phylum", "class", "order"):
            if c[rk]:
                added[f][rk].add(c[rk])
        per_order[(f, c["order"])] += 1
        per_species[(f, c["taxid"])] += 1
        pool.remove(c)

    out = indir / "harvest_selection.tsv"
    if taken:
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(taken[0].keys()), delimiter="\t")
            w.writeheader(); w.writerows(taken)

    print(f"\n=== HARVEST SELECTION (greedy, maximise distinct ORDERS added) ===")
    print(f"candidates in gate-passing orthogroups: {len(cands)}   selected: {len(taken)}")
    print(f"\n{'family':<22}{'selected':>9}{'new orders':>12}{'new classes':>13}{'new phyla':>11}")
    for f in sorted({t['family'] for t in taken}):
        sel = [t for t in taken if t["family"] == f]
        print(f"{f:<22}{len(sel):>9}{len(added[f]['order']):>12}"
              f"{len(added[f]['class']):>13}{len(added[f]['phylum']):>11}")
    print(f"\n=== selected sequences per phylum ===")
    for p, n in Counter(t["phylum"] for t in taken).most_common():
        print(f"  {p or '(unranked)':<28}{n:>5}")
    print(f"\n=== new ORDERS added, by family ===")
    for f in sorted(added):
        if added[f]["order"]:
            print(f"  {f}: {sorted(added[f]['order'])}")

    (indir / "spread_inventory.json").write_text(json.dumps(
        {f: {rk: sorted(v) for rk, v in d.items() if rk != "species"}
         for f, d in inv.items()}, indent=2) + "\n")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
