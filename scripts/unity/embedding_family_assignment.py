#!/usr/bin/env python3
"""embedding_family_assignment.py — is a candidate actually IN a known family?

Bead berghia-chemogpcrs-nsew, deliverable 4. `emb_nonchemo_family` names the
NEAREST non-chemoreceptor family prototype for every candidate unconditionally —
including candidates that are nowhere near any of them. Nearest-family alone
therefore cannot answer "are any Berghia candidates actually bioamine / opsin /
peptide receptors"; on its own it would label all 790 as something.

This script calibrates the assignment against the reference set's OWN spread.
For each family F it computes the squared-Mahalanobis distance from every
labelled reference member of F to F's prototypes, giving F's empirical
within-family distance distribution. A candidate's distance to F is then
expressed as a percentile of that distribution:

    pct <= 95   the candidate sits inside the envelope F's own members occupy
                -> a genuine "looks like an F" call
    pct  > 100  the candidate is further from F than EVERY reference F is
                -> nearest-family is an artifact of argmin, not membership

Specificity is reported alongside: the margin (in the same units) between the
best family and the runner-up, so a candidate that is ambiguously between two
families is visible as such.

The same tied precision / multi-prototype construction as
embedding_evidence.mahalanobis_channel is used, so these distances are exactly
the ones the production novelty axis is built from — this is a readout of the
production geometry, not a parallel scoring scheme.

CAVEAT (must survive into any report): embedding proximity is NOT orthology. A
high-confidence call here is one line of evidence that warrants checking with
the dedicated three-source classifier (HMM scan + orthogroup vote +
phylogenetic placement), which has never been run. It is not a verdict.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from build_embedding_channel import (  # noqa: E402
    load_ref_labels,
    novelty_reference_labels,
)
from embedding_evidence import (  # noqa: E402
    family_prototypes,
    load_embeddings,
    mahalanobis_sq,
    shrinkage_precision,
)


def dist_to_family(vec: np.ndarray, protos: np.ndarray, precision: np.ndarray) -> float:
    return min(mahalanobis_sq(vec, p, precision) for p in np.atleast_2d(protos))


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidate-npz", required=True)
    p.add_argument("--ref-npz", required=True)
    p.add_argument("--ref-labels", required=True)
    p.add_argument("--k", type=int, default=3)
    p.add_argument("--pct-threshold", type=float, default=95.0,
                   help="percentile of the family's own member-distance "
                        "distribution below which a candidate counts as inside "
                        "the family envelope (default 95)")
    p.add_argument("--out-tsv", required=True)
    a = p.parse_args()

    cand = load_embeddings(a.candidate_npz)
    ref = load_embeddings(a.ref_npz)
    labels = novelty_reference_labels(load_ref_labels(a.ref_labels))
    labels = {i: f for i, f in labels.items() if i in ref}
    print(f"[assign] candidates={len(cand)} references_resolved={len(labels)}")
    if not cand or not labels:
        print("[assign] FATAL: empty inputs", file=sys.stderr)
        return 1

    protos = family_prototypes(ref, labels, a.k)
    by_family: dict = {}
    for i, f in labels.items():
        by_family.setdefault(f, []).append(np.asarray(ref[i], dtype=float))
    centered = np.vstack([np.array(v) - np.array(v).mean(axis=0)
                          for v in by_family.values()])
    precision = shrinkage_precision(centered)

    # Empirical within-family distance distribution, per family.
    envelope: dict = {}
    for fam, protomat in protos.items():
        d = np.array([dist_to_family(v, protomat, precision) for v in by_family[fam]])
        envelope[fam] = d
        print(f"[assign] {fam:<24} n={len(d):<5d} "
              f"median={np.median(d):9.2f} p95={np.percentile(d,95):9.2f} "
              f"max={d.max():9.2f}")

    fams = sorted(protos)
    rows = []
    for cid, vec in cand.items():
        vec = np.asarray(vec, dtype=float)
        d = {f: dist_to_family(vec, protos[f], precision) for f in fams}
        order = sorted(fams, key=lambda f: d[f])
        best, second = order[0], (order[1] if len(order) > 1 else None)
        env = envelope[best]
        pct = float((env < d[best]).sum()) / len(env) * 100.0
        rows.append({
            "id": cid,
            "best_family": best,
            "best_dist": d[best],
            "pct_within_family": pct,
            "inside_envelope": pct <= a.pct_threshold,
            "runner_up_family": second,
            "runner_up_dist": d[second] if second else np.nan,
            "margin": (d[second] - d[best]) if second else np.nan,
            **{f"d_{f}": d[f] for f in fams},
        })

    df = pd.DataFrame(rows).sort_values(["inside_envelope", "pct_within_family"],
                                        ascending=[False, True])
    df.to_csv(a.out_tsv, sep="\t", index=False)
    print(f"\n[assign] wrote {a.out_tsv} ({len(df)} rows)")

    inside = df[df["inside_envelope"]]
    print(f"\n=== CANDIDATES INSIDE A KNOWN FAMILY ENVELOPE "
          f"(pct <= {a.pct_threshold}) : {len(inside)}/{len(df)} ===")
    if len(inside):
        print(inside["best_family"].value_counts().to_string())
        print(f"\n{'id':<34} {'family':<22} {'pct':>7} {'dist':>11} {'margin':>11}")
        for _, r in inside.sort_values("pct_within_family").head(40).iterrows():
            print(f"{r['id']:<34} {r['best_family']:<22} "
                  f"{r['pct_within_family']:7.2f} {r['best_dist']:11.2f} "
                  f"{r['margin']:11.2f}")
    else:
        print("(none — no candidate falls inside any family's own member spread)")

    print(f"\n=== nearest-family assignment for ALL {len(df)} candidates "
          "(argmin, uncalibrated) ===")
    print(df["best_family"].value_counts().to_string())
    print("\npct_within_family distribution (100 = further than every real member):")
    s = df["pct_within_family"]
    for q in (0, 1, 5, 10, 25, 50, 75, 90, 100):
        print(f"  p{q}={np.percentile(s, q):.2f}")
    print(f"  == 100.00 (beyond every reference member): "
          f"{int((s >= 100.0).sum())}/{len(s)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
