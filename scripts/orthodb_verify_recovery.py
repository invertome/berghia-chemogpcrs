#!/usr/bin/env python3
"""Verify and merge the recovered identifier-gap joins.

A join key that maps an accession to the WRONG OrthoDB gene is worse than no
mapping, because it silently injects a foreign sequence into the reference
envelope under a curated family label. So nothing recovered by an alternate key
is accepted on the strength of the string match. Every candidate must clear:

  ORGANISM   the OrthoDB gene id encodes its organism as a taxid prefix. That
             taxid must equal the taxid UniProt reports for the accession. A
             mismatch means the key crossed species and the hit is discarded.
  SEQUENCE   where the exact-sequence route also fired for the same accession,
             the identifier route must agree with it. Disagreement is reported
             as a CONFLICT, never silently resolved in favour of either.

Recovered joins are tiered by how strongly they are verified, so downstream use
can choose its own risk tolerance rather than inheriting mine:

  tier1_sequence_exact   exact protein sequence match. No identifier needed.
  tier2_id_seq_agree     identifier route, corroborated by the sequence route.
  tier3_id_organism_only identifier route, organism checks out, no sequence
                         corroboration available.

Usage:
    python3 scripts/orthodb_verify_recovery.py --indir results/ranking/diagnostics/orthodb
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

GID = re.compile(r"^(\d+)_\d+:")


def read_rows(path: Path) -> list[list[str]]:
    if not path.exists():
        return []
    with open(path, newline="") as fh:
        return [r for r in csv.reader(fh, delimiter="\t") if r]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--indir", default="results/ranking/diagnostics/orthodb")
    args = ap.parse_args(argv)
    indir = Path(args.indir)

    # accession -> UniProt-reported taxid (the organism authority)
    alt: dict[str, dict] = {}
    with open(indir / "identifier_gap_alt_xrefs.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            alt[r["accession"]] = r

    # external id / gene name -> accessions that claimed it
    id_owner: dict[str, set] = defaultdict(set)
    name_owner: dict[str, set] = defaultdict(set)
    for acc, field, val in read_rows(indir / "identifier_gap_alt_ids.txt"):
        if field == "gene_names":
            name_owner[val.strip().lower()].add(acc)
        else:
            id_owner[val.strip()].add(acc)
            id_owner[re.sub(r"\.\d+$", "", val.strip())].add(acc)

    cand: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))

    for row in read_rows(indir / "recover_routeA_hits.tsv"):
        if len(row) >= 2:
            for acc in id_owner.get(row[0].strip(), ()):
                cand[acc]["A"].add(row[1])
    for row in read_rows(indir / "recover_routeB_hits.tsv"):
        if len(row) >= 2:
            for acc in name_owner.get(row[0].strip().lower(), ()):
                cand[acc]["B"].add(row[1])
    for row in read_rows(indir / "recover_routeC_hits.tsv"):
        if len(row) >= 2:
            cand[row[0].strip()]["C"].add(row[1])

    results, conflicts = [], []
    for acc in sorted(alt):
        want_tax = alt[acc]["taxid"]
        routes = cand.get(acc, {})
        seqhits = {g for g in routes.get("C", set())}
        idhits = {g for g in routes.get("A", set()) | routes.get("B", set())}

        def org_ok(g):
            m = GID.match(g)
            return bool(m) and m.group(1) == want_tax

        seq_ok = {g for g in seqhits if org_ok(g)}
        id_ok = {g for g in idhits if org_ok(g)}
        seq_badorg = seqhits - seq_ok
        id_badorg = idhits - id_ok

        if seq_ok and id_ok and not (seq_ok & id_ok):
            conflicts.append({
                "accession": acc, "sequence_route": ";".join(sorted(seq_ok)),
                "identifier_route": ";".join(sorted(id_ok)),
            })
            tier, chosen = "conflict", ""
        elif seq_ok:
            tier = "tier2_id_seq_agree" if (seq_ok & id_ok) else "tier1_sequence_exact"
            chosen = sorted(seq_ok)[0]
        elif id_ok:
            tier, chosen = "tier3_id_organism_only", sorted(id_ok)[0]
        else:
            tier, chosen = "unrecovered", ""

        results.append({
            "accession": acc,
            "uniprot_taxid": want_tax,
            "species": alt[acc].get("uniprot_id", ""),
            "routeA_hits": len(routes.get("A", ())),
            "routeB_hits": len(routes.get("B", ())),
            "routeC_hits": len(routes.get("C", ())),
            "rejected_wrong_organism": len(seq_badorg) + len(id_badorg),
            "tier": tier,
            "odb_gene_id": chosen,
        })

    # Many-to-one collisions. Two distinct accessions resolving to ONE OrthoDB
    # gene means at least one of them is wrong, or they are isoforms sharing a
    # gene-level identifier. Either way the weaker-tier claim must not be
    # silently kept: it would put two anchor labels on one sequence.
    TIER_RANK = {"tier1_sequence_exact": 3, "tier2_id_seq_agree": 3,
                 "tier3_id_organism_only": 1}
    by_gene: dict[str, list] = defaultdict(list)
    for r in results:
        if r["odb_gene_id"]:
            by_gene[r["odb_gene_id"]].append(r)
    collisions = []
    for gene, claimants in by_gene.items():
        if len(claimants) < 2:
            continue
        best = max(TIER_RANK.get(c["tier"], 0) for c in claimants)
        winners = [c for c in claimants if TIER_RANK.get(c["tier"], 0) == best]
        collisions.append({
            "odb_gene_id": gene,
            "claimants": ";".join(sorted(c["accession"] for c in claimants)),
            "tiers": ";".join(sorted(c["tier"] for c in claimants)),
            "resolution": ("kept_" + winners[0]["accession"]) if len(winners) == 1
                          else "all_demoted_ambiguous",
        })
        for c in claimants:
            if len(winners) != 1 or c is not winners[0]:
                c["tier"] = "collision_demoted"
                c["odb_gene_id"] = ""

    if collisions:
        with open(indir / "recovery_collisions.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(collisions[0].keys()), delimiter="\t")
            w.writeheader(); w.writerows(collisions)

    with open(indir / "recovery_verdicts.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(results)

    if conflicts:
        with open(indir / "recovery_conflicts.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(conflicts[0].keys()), delimiter="\t")
            w.writeheader(); w.writerows(conflicts)

    # Augmented gene map: original join plus everything verified here.
    recovered = [r for r in results if r["odb_gene_id"]]
    with open(indir / "anchor_gene_map.tsv", newline="") as fh:
        base = [r for r in csv.reader(fh, delimiter="\t") if r]
    with open(indir / "anchor_gene_map_recovered.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerows(base)
        for r in recovered:
            w.writerow([r["accession"], r["odb_gene_id"], f"RECOVERED:{r['tier']}"])

    tiers = Counter(r["tier"] for r in results)
    recovered = [r for r in results if r["odb_gene_id"]]
    summary = {
        "identifier_gap_anchors": len(results),
        "recovered_total": len(recovered),
        "by_tier": dict(tiers.most_common()),
        "conflicts": len(conflicts),
        "rejected_wrong_organism_total": sum(r["rejected_wrong_organism"] for r in results),
    }
    (indir / "recovery_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if conflicts:
        print("\nCONFLICTS (identifier route disagrees with sequence route):")
        for c in conflicts[:20]:
            print(f"  {c['accession']}: seq={c['sequence_route']} id={c['identifier_route']}")
    print(f"\nwrote {indir / 'recovery_verdicts.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
