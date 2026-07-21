#!/usr/bin/env python3
"""Assign each harvested sequence a taxonomic stratum.

Three strata, because they carry different weight as reference material and
lumping them would hide that:

  invertebrate           the harvest's purpose -- clades the anchor set is
                         thin or empty in.
  invertebrate_chordate  tunicates and cephalochordates. Chordates, but not
                         vertebrates, so they widen the envelope rather than
                         deepening its existing vertebrate bias.
  divergent_vertebrate   vertebrates that survived the selection stage's
                         exclusion of Mammalia/Aves/Actinopteri/Amphibia/
                         Lepidosauria -- chondrichthyans, lampreys, coelacanth,
                         lungfish. Legitimate, but they deepen a bias the
                         anchor set already has, so they are labelled and
                         counted separately.

Membership is decided from the NCBI LINEAGE (Vertebrata, taxid 7742), not from
a hand-maintained list of class names. A class-name whitelist silently
misfiles every organism NCBI leaves unranked -- and 24 of the previous
harvest's vertebrates carry no class rank at all.

Usage:
    python3 scripts/orthodb_harvest_stratify.py \
        --quality results/ranking/diagnostics/orthodb/harvest_quality.tsv \
        --lineage results/ranking/diagnostics/orthodb/orthodb_organism_lineage.tsv \
        --species results/ranking/diagnostics/orthodb/odb_species.tsv \
        --out     results/ranking/diagnostics/orthodb/harvest_final.tsv
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import Counter
from pathlib import Path

VERTEBRATA = 7742
CHORDATA = 7711
GENE_ID_RE = re.compile(r"^(\d+_\d+):")


def load_lineages(path: Path) -> dict[str, set[int]]:
    out: dict[str, set[int]] = {}
    with open(path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2:
                out[row[0]] = {int(x) for x in row[1].split(",") if x}
    return out


def load_org_taxid(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with open(path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2:
                out[row[1]] = row[0]
    return out


def stratum_of(lineage: set[int]) -> str:
    if VERTEBRATA in lineage:
        return "divergent_vertebrate"
    if CHORDATA in lineage:
        return "invertebrate_chordate"
    return "invertebrate"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quality", required=True)
    ap.add_argument("--lineage", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    lineages = load_lineages(Path(args.lineage))
    org_taxid = load_org_taxid(Path(args.species))

    with open(args.quality, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    unresolved: list[str] = []
    out_rows = []
    for r in rows:
        m = GENE_ID_RE.match(r["gene_id"])
        org = m.group(1) if m else ""
        tx = org_taxid.get(org, "")
        lin = lineages.get(tx)
        if lin is None:
            unresolved.append(r["gene_id"])
            stratum = "unresolved"
        else:
            stratum = stratum_of(lin)
        cols = {k: v for k, v in r.items() if k != "sequence"}
        cols["stratum"] = stratum
        cols["sequence"] = r.get("sequence", "")
        out_rows.append(cols)

    if unresolved:
        raise SystemExit(
            f"ERROR: {len(unresolved)} harvested gene(s) have no cached NCBI "
            f"lineage (e.g. {unresolved[:5]}). Refusing to stratify on an "
            "incomplete lineage map -- an unresolved organism would silently "
            "land in the 'invertebrate' stratum and be counted as breadth."
        )

    out = Path(args.out)
    tmp = out.with_suffix(out.suffix + ".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    os.replace(tmp, out)

    passing = [r for r in out_rows if r["quality_verdict"] == "PASS"]
    print(f"stratified {len(out_rows)} rows ({len(passing)} PASS)")
    print("all rows :", dict(Counter(r["stratum"] for r in out_rows)))
    print("PASS only:", dict(Counter(r["stratum"] for r in passing)))
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
