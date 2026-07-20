#!/usr/bin/env python3
"""molluscan_calibration_pathcheck.py — prove the molluscan null is embedded in
the SAME geometry as the candidates and references it calibrates.

WHY THIS EXISTS (the failure it is designed to catch)
-----------------------------------------------------
The production novelty channel is built on proteinclip3b + protrek.
proteinclip3b is not a model you can run directly: it is a THREE-step
derivation, and every step is a place the null could silently diverge.

    1. mean-pooled ESM-2 3B  (facebook/esm2_t36_3B_UR50D, MODE=esm in
       scratch_hf_auto_embed.py; mean over non-special tokens, MAXLEN=1024,
       NO L2 normalization at this stage)              -> 2560-d
    2. L2 normalization of that vector -- applied INSIDE
       ONNXModel.predict(apply_norm=True), the default; it is easy to miss
       because no caller passes it explicitly            -> 2560-d, unit norm
    3. the released ProteinCLIP contrastive head
       pretrained/proteinclip_esm2_36.onnx               -> 128-d

If the molluscan background were left in raw ESM-2 space at step 1 while the
candidates and references sit in ProteinCLIP space at step 3, every distance
the null produces would be computed in a different geometry from the thing it
calibrates, and every downstream number would be meaningless -- WITHOUT any
error being raised, because both are just float matrices of the right rank.

The existing `inside_vert_prodbasis` check cannot catch this: it re-derives the
vertebrate readout from the EXISTING production npz and never touches the newly
produced molluscan embeddings at all.

WHAT THIS CHECKS INSTEAD
------------------------
Round-trip reproduction. A deterministic subset of sequences whose production
embeddings ALREADY EXIST is pushed through the EXACT pipeline the null uses,
and the result is compared against the production npz. Agreement proves the new
path reproduces the production path end to end -- model, pooling, truncation,
normalization, checkpoint and projection alike.

Agreement is judged on two levels, because vector-level agreement alone is not
the quantity that matters:

  * VECTOR level: cosine similarity and max abs difference per id, plus the
    vector NORMS -- a normalization divergence shows up in the norms even when
    cosine similarity stays at 1.0.
  * DECISION level: the per-family squared Mahalanobis distances actually
    consumed downstream, and the family argmin they imply. A path that
    reproduces cosines but reorders any family call is still a failure.

Exit code is non-zero on failure, so a calling sbatch with `set -e` stops
before any null is computed from a divergent geometry.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

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


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--produced-npz", required=True,
                   help="npz built by the NULL pipeline for ids that also have "
                        "a production embedding")
    p.add_argument("--production-npz", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--ref-npz", default="", help="enables the decision-level check")
    p.add_argument("--ref-labels", default="")
    p.add_argument("--k", type=int, default=3)
    p.add_argument("--min-cosine", type=float, default=0.9999)
    p.add_argument("--out-json", default="")
    a = p.parse_args()

    new = load_embeddings(a.produced_npz)
    old = load_embeddings(a.production_npz)
    shared = sorted(set(new) & set(old))
    print(f"[pathcheck] {a.label}: produced={len(new)} production={len(old)} "
          f"shared={len(shared)}")
    if not shared:
        print("[pathcheck] FATAL: zero shared ids -- the round trip compares "
              "nothing. This is the exact 'assert key overlap on REAL data' "
              "failure a fixture cannot catch.", file=sys.stderr)
        return 1

    dims_new = {np.asarray(new[i]).shape[-1] for i in shared}
    dims_old = {np.asarray(old[i]).shape[-1] for i in shared}
    print(f"[pathcheck] dim produced={sorted(dims_new)} production={sorted(dims_old)}")
    if dims_new != dims_old:
        print(f"[pathcheck] FATAL: dimensionality differs -- the null is in a "
              f"DIFFERENT SPACE from what it calibrates. This is the "
              f"raw-ESM-2-vs-ProteinCLIP failure mode.", file=sys.stderr)
        return 1

    cos, mad, n_new, n_old = [], [], [], []
    for i in shared:
        u = np.asarray(new[i], dtype=np.float64)
        v = np.asarray(old[i], dtype=np.float64)
        nu, nv = np.linalg.norm(u), np.linalg.norm(v)
        n_new.append(nu)
        n_old.append(nv)
        cos.append(float(u @ v / (nu * nv)) if nu and nv else 0.0)
        mad.append(float(np.abs(u - v).max()))
    cos, mad = np.asarray(cos), np.asarray(mad)
    n_new, n_old = np.asarray(n_new), np.asarray(n_old)

    print(f"[pathcheck] cosine   min={cos.min():.8f} median={np.median(cos):.8f}")
    print(f"[pathcheck] max|diff| median={np.median(mad):.3e} max={mad.max():.3e}")
    print(f"[pathcheck] norms    produced median={np.median(n_new):.6f} "
          f"production median={np.median(n_old):.6f} "
          f"(ratio {np.median(n_new) / np.median(n_old):.6f})")

    ok = bool(cos.min() >= a.min_cosine)
    if not ok:
        print(f"[pathcheck] FATAL: {int((cos < a.min_cosine).sum())} of "
              f"{len(shared)} ids fall below cosine {a.min_cosine}; the null "
              "pipeline does NOT reproduce the production embedding path",
              file=sys.stderr)

    # ---- decision-level: do the family distances and calls agree? ---------
    dec = {}
    if a.ref_npz and a.ref_labels:
        ref = load_embeddings(a.ref_npz)
        labels = novelty_reference_labels(load_ref_labels(a.ref_labels))
        labels = {i: f for i, f in labels.items() if i in ref}
        protos = family_prototypes(ref, labels, a.k)
        by_family: dict = {}
        for i, f in labels.items():
            by_family.setdefault(f, []).append(np.asarray(ref[i], dtype=float))
        centered = np.vstack([np.array(v) - np.array(v).mean(axis=0)
                              for v in by_family.values()])
        precision = shrinkage_precision(centered)
        fams = sorted(protos)

        def dvec(store, i):
            x = np.asarray(store[i], dtype=float)
            return np.array([min(mahalanobis_sq(x, q, precision)
                                 for q in np.atleast_2d(protos[f])) for f in fams])

        D_new = np.vstack([dvec(new, i) for i in shared])
        D_old = np.vstack([dvec(old, i) for i in shared])
        rel = np.abs(D_new - D_old) / np.maximum(np.abs(D_old), 1e-12)
        arg_new = [fams[j] for j in D_new.argmin(axis=1)]
        arg_old = [fams[j] for j in D_old.argmin(axis=1)]
        n_disagree = sum(1 for x, y in zip(arg_new, arg_old) if x != y)
        dec = {"max_rel_dist_diff": float(rel.max()),
               "median_rel_dist_diff": float(np.median(rel)),
               "family_argmin_disagreements": n_disagree,
               "n_compared": len(shared)}
        print(f"[pathcheck] family distances: median rel diff "
              f"{np.median(rel):.3e}, max {rel.max():.3e}")
        print(f"[pathcheck] family argmin disagreements: {n_disagree}/{len(shared)}")
        if n_disagree:
            ok = False
            print("[pathcheck] FATAL: the reproduced path assigns a DIFFERENT "
                  "family to at least one sequence; vector-level agreement is "
                  "not sufficient when the call itself moves", file=sys.stderr)

    result = {"label": a.label, "n_shared": len(shared),
              "dim": sorted(dims_new)[0],
              "cosine_min": float(cos.min()),
              "cosine_median": float(np.median(cos)),
              "maxabs_median": float(np.median(mad)),
              "norm_ratio_median": float(np.median(n_new) / np.median(n_old)),
              "passed": ok, **dec}
    if a.out_json:
        with open(a.out_json, "w") as fh:
            json.dump(result, fh, indent=2)
    print(f"[pathcheck] {a.label}: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
