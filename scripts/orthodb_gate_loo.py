#!/usr/bin/env python3
"""Derive and validate the harvest gate by leave-one-out on invertebrate anchors.

The molluscan harvest failed on granularity: a median 2,600-member orthogroup
with its family attested by ~2 anchors is 0.084% coverage, and the largest was
68,077 members named for the superfamily. So a harvest is only defensible from
orthogroups where the family assignment is densely enough attested. Three
properties define the gate, and none of their thresholds is chosen a priori -
they are swept and the curve is reported so the tradeoff is visible:

  SIZE                 members in the orthogroup.
  ATTESTATION DENSITY  characterized anchors of the modal curated family
                       seeding it, absolutely (min_seeds) and as a fraction of
                       members (min_density).
  PURITY               whether anchors of a DIFFERENT curated family are also
                       present. Mixed groups are disqualified or flagged.

The gate is then validated against ground truth we already hold. Every
invertebrate anchor with a curated family is held out in turn and we ask:
would the harvest have recovered it, and with the correct family? That gives
two numbers per gate setting:

  RECOVERY   fraction of held-out anchors whose orthogroup passes the gate,
             i.e. the harvest would have found them at all.
  PRECISION  among those recovered, fraction assigned the family curation
             actually gives them. This is the number that matters, because a
             confident wrong family assignment is worse than no assignment.

Two hold-out regimes, because they answer different questions:

  LOO          leave one anchor out. Other anchors of its clade remain, which
               matches production for clades we already have anchors in.
  LEAVE-CLADE  leave the anchor's entire clade out. This is the honest test for
               Cnidaria, non-vertebrate deuterostomes and non-bilaterians,
               which have LITERALLY ZERO anchors - for them, production has no
               same-clade anchor to lean on, so LOO would flatter the gate.

Usage:
    python3 scripts/orthodb_gate_loo.py --indir results/ranking/diagnostics/orthodb
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

SIZE_GRID = [500, 1000, 1500, 2000, 2500, 5000, 10000, 10**9]
SEED_GRID = [1, 2, 3, 5]
PURITY_MODES = ["strict", "modal80"]


def gate_verdict(
    members: int,
    seed_families: Counter,
    max_size: int,
    min_seeds: int,
    purity: str,
    min_density: float = 0.0,
) -> tuple[bool, str, str]:
    """-> (passes, assigned_family, reason_if_failed).

    seed_families counts CHARACTERIZED anchors per curated family, already
    excluding whatever is being held out.
    """
    if not seed_families:
        return False, "", "no_seed_anchors"
    if members > max_size:
        return False, "", "too_large"

    total = sum(seed_families.values())
    top = seed_families.most_common()
    modal_family, modal_n = top[0]
    if len(top) > 1 and top[1][1] == modal_n:
        return False, "", "tied_seed_families"

    if purity == "strict" and len(seed_families) > 1:
        return False, modal_family, "mixed_family_seeds"
    if purity == "modal80" and modal_n / total < 0.8:
        return False, modal_family, "impure_seeds"

    if modal_n < min_seeds:
        return False, modal_family, "too_few_seeds"
    if min_density and members and (modal_n / members) < min_density:
        return False, modal_family, "too_sparse"
    return True, modal_family, ""


def load(indir: Path, level: str) -> dict:
    with open(indir / "mapping_audit.tsv", newline="") as fh:
        audit = {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}

    gene_to_accs: dict[str, set] = defaultdict(set)
    with open(indir / "anchor_gene_map.tsv", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2:
                gene_to_accs[row[1]].add(row[0])

    og_level, og_name = {}, {}
    with open(indir / "anchor_og_meta.tsv", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2:
                og_level[row[0]] = row[1]
                og_name[row[0]] = row[2] if len(row) > 2 else ""

    anchor_ogs: dict[str, set] = defaultdict(set)
    with open(indir / "anchor_og2genes.tsv", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2 and og_level.get(row[0]) == level:
                for acc in gene_to_accs.get(row[1], ()):
                    anchor_ogs[acc].add(row[0])

    og_size: Counter = Counter()
    with open(indir / "anchor_og_membership.tsv", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2 and og_level.get(row[0]) == level:
                og_size[row[0]] += 1

    og_anchors: dict[str, list] = defaultdict(list)
    for acc, ogs in anchor_ogs.items():
        rec = audit.get(acc)
        if rec is None or rec["use_primary"] != "True":
            continue
        for og in ogs:
            og_anchors[og].append(rec)

    return {"audit": audit, "anchor_ogs": anchor_ogs, "og_size": og_size,
            "og_anchors": og_anchors, "og_name": og_name}


def is_characterized(rec: dict) -> bool:
    ev = rec.get("evidence", "")
    return ev.startswith("experimental") or ev == "characterized-gtopdb"


def evaluate(data: dict, regime: str, max_size: int, min_seeds: int,
             purity: str, min_density: float = 0.0) -> dict:
    """Hold out each invertebrate anchor and see if the gate recovers it."""
    audit = data["audit"]
    targets = [
        r for r in audit.values()
        if r["use_primary"] == "True"
        and r["mapped"] == "True"
        and r["clade"] not in ("Vertebrata", "non_Metazoa", "unknown")
        and data["anchor_ogs"].get(r["accession"])
    ]

    per_clade: dict[str, Counter] = defaultdict(Counter)
    rows = []
    for tgt in targets:
        acc = tgt["accession"]
        clade = tgt["clade"]
        held = {acc}
        if regime == "leave_clade":
            held = {a for a, r in audit.items() if r["clade"] == clade}

        best = None
        for og in data["anchor_ogs"][acc]:
            seeds: Counter = Counter()
            for m in data["og_anchors"].get(og, []):
                if m["accession"] in held or not is_characterized(m):
                    continue
                seeds[m["family"]] += 1
            ok, fam, why = gate_verdict(
                data["og_size"][og], seeds, max_size, min_seeds, purity, min_density
            )
            cand = (ok, fam, why, og, data["og_size"][og], sum(seeds.values()))
            if best is None or (ok and not best[0]) or (
                ok == best[0] and cand[4] < best[4]
            ):
                best = cand
        if best is None:
            continue
        ok, fam, why, og, size, nseed = best
        correct = ok and fam == tgt["family"]
        bucket = "recovered_correct" if correct else (
            "recovered_wrong" if ok else "not_recovered")
        per_clade[clade][bucket] += 1
        per_clade["ALL"][bucket] += 1
        rows.append({
            "accession": acc, "clade": clade, "true_family": tgt["family"],
            "og_id": og, "og_size": size, "seed_anchors": nseed,
            "assigned_family": fam, "passed_gate": str(ok),
            "outcome": bucket, "fail_reason": why,
        })

    def stats(c: Counter) -> dict:
        rec = c["recovered_correct"] + c["recovered_wrong"]
        tot = rec + c["not_recovered"]
        return {
            "n": tot,
            "recovered": rec,
            "recovery": rec / tot if tot else 0.0,
            "precision": c["recovered_correct"] / rec if rec else 0.0,
            "correct": c["recovered_correct"],
            "wrong": c["recovered_wrong"],
        }

    return {"per_clade": {k: stats(v) for k, v in per_clade.items()}, "rows": rows}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--indir", default="results/ranking/diagnostics/orthodb")
    ap.add_argument("--level", default="33208")
    args = ap.parse_args(argv)
    indir = Path(args.indir)
    data = load(indir, args.level)

    sweep = []
    for regime in ("loo", "leave_clade"):
        for purity in PURITY_MODES:
            for max_size in SIZE_GRID:
                for min_seeds in SEED_GRID:
                    res = evaluate(data, regime, max_size, min_seeds, purity)
                    for clade, s in res["per_clade"].items():
                        sweep.append({
                            "regime": regime, "purity": purity,
                            "max_size": max_size, "min_seeds": min_seeds,
                            "clade": clade, **s,
                        })

    cols = ["regime", "purity", "max_size", "min_seeds", "clade", "n",
            "recovered", "recovery", "precision", "correct", "wrong"]
    with open(indir / "gate_loo_sweep.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in sweep:
            w.writerow({**r, "recovery": round(r["recovery"], 4),
                        "precision": round(r["precision"], 4)})

    def show(regime, clade, purity="strict"):
        print(f"\n=== {regime.upper()} / clade={clade} / purity={purity} ===")
        print(f"{'max_size':>9}" + "".join(f"{'K='+str(k):>22}" for k in SEED_GRID))
        print(f"{'':>9}" + "".join(f"{'recov/prec (n)':>22}" for _ in SEED_GRID))
        for ms in SIZE_GRID:
            cells = []
            for k in SEED_GRID:
                m = [r for r in sweep if r["regime"] == regime and r["purity"] == purity
                     and r["max_size"] == ms and r["min_seeds"] == k
                     and r["clade"] == clade]
                if not m:
                    cells.append(f"{'-':>22}"); continue
                s = m[0]
                cells.append(f"{s['recovery']*100:>7.0f}%/{s['precision']*100:>3.0f}% ({s['recovered']:>3}/{s['n']:>3})")
            lbl = "inf" if ms >= 10**9 else str(ms)
            print(f"{lbl:>9}" + "".join(cells))

    for clade in ("ALL", "Ecdysozoa", "Mollusca"):
        show("loo", clade)
    show("leave_clade", "Ecdysozoa")
    show("leave_clade", "Mollusca")
    show("loo", "ALL", purity="modal80")

    counts = Counter(r["clade"] for r in evaluate(data, "loo", 10**9, 1, "strict")["rows"])
    print(f"\n=== observations available per clade ===\n{json.dumps(dict(counts), indent=2)}")
    print(f"\nwrote {indir / 'gate_loo_sweep.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
