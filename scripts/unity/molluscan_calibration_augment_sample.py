#!/usr/bin/env python3
"""molluscan_calibration_augment_sample.py — add HMM/Pfam detection evidence to
the null sample table, so the chemo-gate enrichment caveat can be bounded from
inside the sample as well as across the TM gate.

WHY
---
The 40,026 phase-1a candidates passed a 6TM *chemoreceptor* gate, so as a
calibration null they are enriched for chemoreceptor-like sequences rather than
being a neutral sample of molluscan class-A GPCRs. Comparing post-gate against
pre-gate-only bounds the TM gate's contribution, but it says nothing about the
HMM detection step that ran BEFORE it -- and that step used chemoreceptor-biased
profiles.

`class_phase1a.tsv` records which Pfam actually hit. Names AND clan membership
below were resolved through the InterPro API, not recalled:

    accession  name                                          clan
    PF00001    7 transmembrane receptor (rhodopsin family)    CL0192 GPCR_A
    PF10324    Serpentine 7TM GPCR chemoreceptor Srw          CL0192 GPCR_A
    PF05296    Taste receptor protein (TAS2R)                 CL0192 GPCR_A
    PF08395    7tm Chemosensory receptor  (7tm_7)             CL0176 Chemosens_recp
    PF02949    7tm Odorant receptor       (7tm_6)             CL0176 Chemosens_recp

The clan boundary is NOT the same as the chemoreceptor boundary, and conflating
them would be wrong in a project that gates strictly on class A. Srw and TAS2R
are chemoreceptor families that nonetheless sit INSIDE the rhodopsin clan, so
they are legitimately class A. 7tm_7 and 7tm_6 are a DIFFERENT clan altogether
-- yet classify_gpcr_by_class.py labels them class=A, so they enter a class-A
null as clan contaminants. Both facts are worth separating, so the stratum is
three-way rather than two-way:

    chemo_classA   PF10324 / PF05296 -- chemoreceptor AND rhodopsin clan
    chemo_offclan  PF08395 / PF02949 -- chemoreceptor, NOT rhodopsin clan
    generic_7tm1   PF00001           -- the generic rhodopsin profile

Contrasting `chemo_classA` against `generic_7tm1` bounds the enrichment WITHIN
class A, which is the comparison the caveat actually needs. `chemo_offclan` is
tracked separately so its size is visible and a clan-clean sensitivity run is
possible; it is a contamination estimate, not an enrichment estimate.

If the null shifts substantially between explicitly-chemoreceptor detections
and generic class-A detections, the enrichment matters; if it does not, the
null is dominated by "these are molluscs and the references are not", which is
the hypothesis under test.

Writes in place, atomically. The seq_id set is NEVER changed -- only columns are
added -- so the FASTAs and the embeddings already built from them stay valid.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

# clan membership resolved via the InterPro API (see module docstring)
CHEMO_CLASSA_PFAM = {"PF10324", "PF05296"}   # chemoreceptor, clan CL0192 GPCR_A
CHEMO_OFFCLAN_PFAM = {"PF08395", "PF02949"}  # chemoreceptor, clan CL0176
GENERIC_PFAM = {"PF00001"}                   # generic rhodopsin, clan CL0192


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample-tsv", required=True)
    p.add_argument("--class-tsv", required=True)
    a = p.parse_args()

    ev = {}
    with open(a.class_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ev[r["seq_id"]] = (r.get("evidence_pfam", ""),
                               r.get("evidence_family_hmm", ""))
    print(f"[augment] evidence rows: {len(ev)}")

    with open(a.sample_tsv) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        print("[augment] FATAL: empty sample table", file=sys.stderr)
        return 1

    post = [r for r in rows if r["stratum"] == "postgate_classA"]
    hit = sum(1 for r in post if r["seq_id"] in ev)
    print(f"[augment] key overlap on REAL data: {hit}/{len(post)} post-gate rows "
          f"resolve in the class table")
    if hit != len(post):
        print(f"[augment] FATAL: {len(post) - hit} post-gate rows unresolved; "
              "the detection-evidence stratification would be silently partial",
              file=sys.stderr)
        return 1

    counts = {"chemo_classA": 0, "chemo_offclan": 0, "generic_7tm1": 0,
              "other": 0, "pregate": 0}
    for r in rows:
        pf, fam = ev.get(r["seq_id"], ("", ""))
        r["evidence_pfam"] = pf
        r["evidence_family_hmm"] = fam
        if r["stratum"] != "postgate_classA":
            det = "pregate"
        elif pf in CHEMO_CLASSA_PFAM:
            det = "chemo_classA"
        elif pf in CHEMO_OFFCLAN_PFAM:
            det = "chemo_offclan"
        elif pf in GENERIC_PFAM:
            det = "generic_7tm1"
        else:
            det = "other"
        r["detection_class"] = det
        r["pfam_clan"] = ("CL0176" if pf in CHEMO_OFFCLAN_PFAM
                          else "CL0192" if pf in (CHEMO_CLASSA_PFAM | GENERIC_PFAM)
                          else "")
        counts[det] += 1
    print("[augment] detection_class: " +
          ", ".join(f"{k}={v}" for k, v in counts.items()))

    fields = list(rows[0])
    tmp = a.sample_tsv + ".tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, a.sample_tsv)
    print(f"[augment] wrote {a.sample_tsv} ({len(rows)} rows, seq_id set unchanged)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
