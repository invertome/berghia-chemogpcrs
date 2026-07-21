#!/usr/bin/env python3
"""molluscan_calibration_negcontrol_norm.py — build a SAME-DIMENSION negative
control by deliberately dropping the L2 normalization step.

WHY THIS SPECIFIC CONTROL
-------------------------
The first negative control (raw ESM-2 3B vs production proteinclip3b) failed on
DIMENSIONALITY, 2560 vs 128. That is the cheapest possible way to catch a space
mismatch, so it proved the checker catches a GROSS divergence and nothing more.

The realistic failure mode is subtler and same-dimension. In
`proteinclip.model_utils.ONNXModel.predict`, the input L2 normalization is a
DEFAULT ARGUMENT:

    def predict(self, x, apply_norm: bool = True):
        if apply_norm:
            x /= np.linalg.norm(x)

No caller anywhere passes it explicitly, so it is invisible at every call site
and is exactly the step a reimplementation would silently drop. Dropping it
still yields 128-d vectors, so dimensionality cannot rescue the check --
only the cosine / max|diff| / norm-ratio comparisons can, and whether THOSE
discriminate is what remains unproven.

This script reproduces the projection with `apply_norm=False` and nothing else
changed, so the resulting npz differs from production in exactly one step. The
caller then asserts the pathcheck FAILS on it. If it were to pass, the norms
instrument is not load-bearing and a same-dimension divergence could reach a
production null undetected.

Reuses the ESM-2 3B probe embeddings already produced by
molluscan_calibration_pathcheck.sh -- no GPU forward pass is repeated.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--base-npz", required=True,
                   help="ESM-2 3B probe embeddings (already built; not recomputed)")
    p.add_argument("--out-npz", required=True)
    p.add_argument("--proteinclip-repo", required=True)
    p.add_argument("--head", type=int, default=36)
    a = p.parse_args()

    sys.path.insert(0, a.proteinclip_repo)
    from proteinclip import model_utils  # noqa: E402

    m = model_utils.load_proteinclip("esm", a.head)
    base = np.load(a.base_npz)

    out, in_norms, out_norms = {}, [], []
    for k in base.files:
        v = np.asarray(base[k], dtype=np.float32)
        in_norms.append(float(np.linalg.norm(v)))
        # THE ONE DELIBERATE DEVIATION: apply_norm=False.
        # Production leaves this at its True default.
        w = m.predict(v, apply_norm=False).astype(np.float32)
        out[k] = w
        out_norms.append(float(np.linalg.norm(w)))

    os.makedirs(os.path.dirname(a.out_npz) or ".", exist_ok=True)
    np.savez(a.out_npz + ".tmp.npz", **out)
    os.replace(a.out_npz + ".tmp.npz", a.out_npz)

    dim = next(iter(out.values())).shape[-1]
    print(f"[negctl] projected {len(out)} probes with apply_norm=FALSE -> dim={dim}")
    print(f"[negctl] ESM-2 input norms : median {np.median(in_norms):.4f} "
          f"(range {min(in_norms):.4f}-{max(in_norms):.4f})")
    print(f"[negctl] output norms      : median {np.median(out_norms):.6f} "
          f"(range {min(out_norms):.6f}-{max(out_norms):.6f})")
    print("[negctl] If the output norms are ~1.0 anyway, the ONNX graph "
          "normalizes its OWN output, which would mean the norm-ratio "
          "comparison CANNOT catch a dropped input normalization and the "
          "cosine / max|diff| comparisons are carrying the check alone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
