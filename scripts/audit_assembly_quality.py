#!/usr/bin/env python3
"""Flag assemblies too poor to yield synteny-usable annotation — REVIEW ONLY.

This audit never edits `drop_reason`: it writes a flag table for a human to
adjudicate. Taxon sampling is a scientific decision, not a threshold's to make.

Why the rule is CLADE-RELATIVE AND two-signal (both conditions required):

  * Absolute size fails. `Gyrodactylus salaris` (67 Mb) and `Dipylidium caninum`
    (82 Mb) are legitimately compact parasite genomes that annotate fine
    (14,568 / 13,062 proteins). Any absolute Mb floor that catches a 34 Mb
    "nudibranch" also wrongly condemns them. Comparing each species to its OWN
    clade's median separates "compact genome" from "fragment of a genome".

  * N50 alone fails, and would have been a serious mistake here. Empirically,
    low-N50 assemblies DO annotate: of completed species, the N50 <3 kb band
    still averaged ~38k proteins and `Melibe leonina` (N50 2,931 bp, 565 Mb)
    yielded 59,607. An N50<10 kb gate alone would flag 107/453 species —
    including most of the nudibranch sampling — for no reason.

  * Size alone ALSO fails on the other side: `Hirudinaria manillensis` sits at
    19.4% of the annelid median but has N50 = 9.6 Mb — a chromosome-scale
    assembly of a genuinely compact leech genome. Excellent data.

So a species is flagged only when it is BOTH grossly incomplete for its clade
AND too fragmented to carry a gene: size < MIN_CLADE_FRACTION x clade median
AND contig_n50 < MIN_N50_BP. N50 floor is grounded in biology, not a percentile:
a metazoan gene spans ~10-30 kb, so below ~10 kb a gene cannot sit on one contig.

Clades with fewer than MIN_CLADE_N members have no stable median; their species
are reported as UNASSESSED rather than silently passed.
"""
from __future__ import annotations

import argparse
import csv
import statistics
import sys
from collections import defaultdict
from typing import Dict, List

MIN_CLADE_FRACTION = 0.20   # grossly incomplete relative to its own clade
MIN_N50_BP = 10_000         # a metazoan gene (~10-30 kb) cannot sit on one contig below this
MIN_CLADE_N = 3             # fewer members than this -> no stable median


def load(path: str) -> List[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def clade_key(clade: str) -> str:
    """Normalised clade label for grouping.

    The inventory mixes two naming conventions from different source batches --
    'Bivalvia'/'bivalvia', 'Annelida'/'annelida', 'Platyhelminthes'/
    'platyhelminthes' are the SAME clade. Grouping on the raw string splits each
    in two, shrinking both medians' support and stranding members of a
    well-sampled clade in the 'too small for a median' bucket.
    """
    return (clade or "").strip().lower()


def audit(rows: List[dict]) -> List[dict]:
    by_clade: Dict[str, List[int]] = defaultdict(list)
    for r in rows:
        try:
            L = int(r.get("total_length_bp") or 0)
        except ValueError:
            L = 0
        if L > 0:
            by_clade[clade_key(r.get("clade", ""))].append(L)
    medians = {c: statistics.median(v) for c, v in by_clade.items() if len(v) >= MIN_CLADE_N}

    out = []
    for r in rows:
        try:
            L = int(r.get("total_length_bp") or 0)
            n50 = int(r.get("contig_n50") or 0)
        except ValueError:
            L, n50 = 0, 0
        clade = r.get("clade", "")
        med = medians.get(clade_key(clade))
        if L <= 0:
            verdict, frac = "UNASSESSED_no_size", ""
        elif med is None:
            verdict, frac = "UNASSESSED_small_clade", ""
        else:
            frac = round(100.0 * L / med, 1)
            too_small = L < MIN_CLADE_FRACTION * med
            too_frag = 0 < n50 < MIN_N50_BP
            verdict = "FLAG_unusable" if (too_small and too_frag) else (
                "watch_small_but_contiguous" if too_small else
                "watch_fragmented_but_complete" if too_frag else "ok")
        out.append({
            "taxid": r.get("taxid", ""), "binomial": r.get("binomial", ""),
            "clade": clade, "accession": r.get("accession", ""),
            "total_length_bp": L, "pct_of_clade_median": frac,
            "contig_n50": n50, "verdict": verdict,
            "existing_drop_reason": r.get("drop_reason", ""),
        })
    return out


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inventory", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)

    rows = audit(load(a.inventory))
    cols = list(rows[0].keys())
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    counts = defaultdict(int)
    for r in rows:
        counts[r["verdict"]] += 1
    print(f"audited {len(rows)} assemblies -> {a.out}")
    for k in sorted(counts):
        print(f"  {k:32s} {counts[k]}")
    print(f"\n=== FLAG_unusable (size < {int(MIN_CLADE_FRACTION*100)}% of clade median AND n50 < {MIN_N50_BP}) ===")
    for r in sorted((x for x in rows if x["verdict"] == "FLAG_unusable"),
                    key=lambda x: x["pct_of_clade_median"]):
        print(f"  {r['taxid']:>8}  {r['binomial']:<32} {r['clade']:<18} "
              f"{r['total_length_bp']/1e6:8.1f} Mb = {r['pct_of_clade_median']:>5}%  n50={r['contig_n50']:>8}"
              f"{'  [already dropped: ' + r['existing_drop_reason'] + ']' if r['existing_drop_reason'] else ''}")
    print("\n=== watch_small_but_contiguous (compact genome, GOOD contiguity -> keep) ===")
    for r in (x for x in rows if x["verdict"] == "watch_small_but_contiguous"):
        print(f"  {r['binomial']:<32} {r['total_length_bp']/1e6:8.1f} Mb = {r['pct_of_clade_median']}%  n50={r['contig_n50']}")
    print("\n=== UNASSESSED (clade too small for a median) ===")
    for r in (x for x in rows if x["verdict"].startswith("UNASSESSED")):
        print(f"  {r['binomial']:<32} clade={r['clade']:<28} {r['verdict']}")


if __name__ == "__main__":
    main()
