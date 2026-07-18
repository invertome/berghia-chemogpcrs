#!/usr/bin/env python3
"""swe_pooling.py — Sliced-Wasserstein distributional pooling (A5).

A5 (epic v4bs): PLM per-residue embeddings are normally mean-pooled into one
per-protein vector, which discards the *distribution* of residue embeddings and
is length-sensitive. Sliced-Wasserstein pooling (Kolouri/Naderializadeh SWE;
Shaw et al.) pools a SET of L per-residue vectors into a fixed-length vector by
sorting their 1-D projections onto random directions and reading off evenly
spaced quantiles. The result is distribution-aware, permutation-invariant, and
length-robust — a candidate bake-off entry for the embedding channel.

Pure/numpy (no torch/GPU/esm); the per-residue embeddings come from HPC
elsewhere. This is a DORMANT producer: the pooling function plus a thin npz
CLI. Randomness flows only through a passed-in seed / ``np.random.Generator``
(never a global), so pooling is fully deterministic.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np


def random_unit_directions(
    n_slices: int, n_dim: int, rng: np.random.Generator
) -> np.ndarray:
    """``n_slices`` uniform random unit directions in R^``n_dim``.

    Gaussian samples normalized to unit length are uniform on the sphere.
    Returns an ``(n_slices, n_dim)`` array with each row of unit L2 norm.
    """
    dirs = rng.standard_normal((n_slices, n_dim))
    norms = np.linalg.norm(dirs, axis=1, keepdims=True)
    norms[norms == 0.0] = 1.0  # degenerate all-zero draw -> leave as-is
    return dirs / norms


def swe_pool(
    residue_vectors: np.ndarray,
    n_slices: int = 64,
    n_quantiles: int = 8,
    seed: int = 0,
) -> np.ndarray:
    """Sliced-Wasserstein pool an ``(L, D)`` set of residue vectors.

    Draws ``n_slices`` seeded random unit directions in R^D, projects the L
    vectors onto each (a dot product), sorts the L projections, and samples
    ``n_quantiles`` evenly-spaced quantiles (``np.linspace(0, 1, n_quantiles)``,
    endpoints inclusive; ``np.quantile`` linear interpolation). The per-slice
    quantile blocks are concatenated into a fixed-length
    ``(n_slices * n_quantiles,)`` vector.

    Sorting makes the pool invariant to the order of the L rows (a SET
    function) and robust to L (quantiles estimate the same 1-D distribution).
    """
    vecs = np.asarray(residue_vectors, dtype=float)
    n_dim = vecs.shape[1]
    rng = np.random.default_rng(seed)
    dirs = random_unit_directions(n_slices, n_dim, rng)

    # (L, n_slices) projections of every residue onto every slice direction
    proj = vecs @ dirs.T
    q_positions = np.linspace(0.0, 1.0, n_quantiles)
    # (n_quantiles, n_slices) -> (n_slices, n_quantiles) -> flat per-slice blocks
    quants = np.quantile(proj, q_positions, axis=0)
    return quants.T.reshape(-1)


def main(argv: Optional[Sequence[str]] = None) -> None:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument(
        "--in-npz", required=True,
        help="input .npz mapping protein_id -> its (L, D) per-residue array",
    )
    ap.add_argument(
        "--out-npz", required=True,
        help="output .npz mapping protein_id -> pooled (n_slices*n_quantiles,) vector",
    )
    ap.add_argument("--n-slices", type=int, default=64)
    ap.add_argument("--n-quantiles", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    data = np.load(args.in_npz)
    pooled = {
        pid: swe_pool(
            data[pid],
            n_slices=args.n_slices,
            n_quantiles=args.n_quantiles,
            seed=args.seed,
        )
        for pid in data.files
    }
    np.savez(args.out_npz, **pooled)
    print(f"[swe_pooling] pooled {len(pooled)} proteins -> {args.out_npz}")


if __name__ == "__main__":
    main()
