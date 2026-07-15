#!/usr/bin/env python3
"""build_tier4_anchors.py — Tier-4 deuterostome anchor layer for the reference
expansion (bead cw3.5, council-directed full sweep).

Sources reviewed human/mouse/zebrafish GPCRs (streamed via the UniProt skill to
tier4_raw.tsv), classifies each into family + subfamily using the SAME curated
pattern classifier as the molluscan tiers, HARD-excludes vertebrate chemosensory
(belt-and-suspenders on top of the query filter), then subsamples per subfamily
(cross-species round-robin, cap CAP) to keep phylogenetic spread without letting
the well-annotated deuterostome pool swamp the in-group molluscan landmarks.

Writes references/anchors/tier4_anchors.{tsv,fasta} in the anchor format
(composite header ANCHOR_<class>_4_<accession>; TSV cols match anchor_set.tsv
plus a `subfamily` column that load_ref_labels safely ignores).
"""
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(__file__))
from build_anchor_set import classify_family, family_to_class  # noqa: E402
from curate_gpcr_references import (  # noqa: E402
    _AMINERGIC_SUBFAMILY_PATTERNS,
    _PEPTIDE_SUBFAMILY_PATTERNS,
)

RAW = sys.argv[1] if len(sys.argv) > 1 else "tier4_raw.tsv"
TIER = os.environ.get("ANCHOR_TIER", "4")
OUT_TSV = os.environ.get("OUT_TSV", "references/anchors/tier4_anchors.tsv")
OUT_FA = os.environ.get("OUT_FA", "references/anchors/tier4_anchors.fasta")
CAP = int(os.environ.get("TIER4_SUBFAMILY_CAP", "12"))
# Residual chemosensory guard (never label these non-chemoreceptor) — covers both
# vertebrate (OR/TAAR/T2R/V1R/V2R/MRGPR) and invertebrate (insect OR/GR, nematode
# serpentine srh/str/sri) chemoreceptor families.
CHEMO = re.compile(
    r"olfactory receptor|odorant receptor|taste receptor|gustatory receptor|"
    r"vomeronasal|trace amine|mas-related|pheromone receptor|serpentine receptor",
    re.I)


def classify_subfamily(name, family):
    if family == "aminergic":
        for pat, sf in _AMINERGIC_SUBFAMILY_PATTERNS:
            if pat.search(name):
                return f"aminergic-{sf}"
        return "aminergic-other"
    if family == "peptide":
        for pat, sf in _PEPTIDE_SUBFAMILY_PATTERNS:
            if pat.search(name):
                return f"peptide-{sf}"
        return "peptide-other"
    return family  # coarse: family IS the subfamily for the rest


records = []
with open(RAW) as fh:
    header = fh.readline()
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 5:
            continue
        acc, name, taxid, organism, seq = p[0], p[1], p[2], p[3], p[4]
        if not seq or CHEMO.search(name):
            continue
        fam = classify_family(name)
        records.append({
            "acc": acc, "name": name, "taxid": taxid,
            "species": organism.split(" (")[0],
            "family": fam, "class": family_to_class(fam),
            "subfamily": classify_subfamily(name, fam), "seq": seq,
        })

# subsample per subfamily, cross-species round-robin
by_sf = defaultdict(list)
for r in records:
    by_sf[r["subfamily"]].append(r)

kept = []
for sf, recs in by_sf.items():
    if len(recs) <= CAP:
        kept.extend(recs)
        continue
    by_org = defaultdict(list)
    for r in recs:
        by_org[r["taxid"]].append(r)
    orgs = sorted(by_org)
    for o in orgs:
        by_org[o].sort(key=lambda r: r["acc"])
    picked, i = [], 0
    while len(picked) < CAP:
        progressed = False
        for o in orgs:
            if i < len(by_org[o]):
                picked.append(by_org[o][i])
                progressed = True
                if len(picked) >= CAP:
                    break
        if not progressed:
            break
        i += 1
    kept.extend(picked)

os.makedirs("references/anchors", exist_ok=True)
with open(OUT_TSV, "w") as t, open(OUT_FA, "w") as f:
    t.write("accession\ttier\ttaxid\tspecies\tfamily\tclass\tevidence\tsubfamily\n")
    for r in sorted(kept, key=lambda r: (r["family"], r["subfamily"], r["acc"])):
        t.write(f"{r['acc']}\t{TIER}\t{r['taxid']}\t{r['species']}\t{r['family']}\t"
                f"{r['class']}\treviewed\t{r['subfamily']}\n")
        f.write(f">ANCHOR_{r['class']}_{TIER}_{r['acc']}\n{r['seq']}\n")

fam_ct = defaultdict(int)
for r in kept:
    fam_ct[r["family"]] += 1
print(f"# streamed {len(records)} classifiable (post chemo-guard) -> kept {len(kept)} "
      f"(cap {CAP}/subfamily, {len(by_sf)} subfamilies)")
print(f"# tier-4 by family: {dict(sorted(fam_ct.items(), key=lambda x: -x[1]))}")
print(f"# wrote {OUT_TSV} + {OUT_FA}")
