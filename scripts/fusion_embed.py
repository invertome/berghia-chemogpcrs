"""Embedding-space concat fusion for the bake-off (bead cw3.16). Pure numpy.

Builds fused reference/candidate npz; scoring is done by the UNCHANGED
scratch_lofo_bakeoff.py (never re-implemented here). See the design spec at
docs/plans/2026-07-16-fusion-consensus-harness.md (Module B).
"""
from __future__ import annotations

import argparse
import sys
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

sys.path.insert(0, "scripts")
from embedding_evidence import load_embeddings


def _l2(v: np.ndarray) -> np.ndarray:
    """L2-normalize a vector; a zero vector is returned unchanged (no div-by-0)."""
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    return v if n == 0 else v / n


def concat_fuse(
    a: Dict[str, np.ndarray],
    b: Dict[str, np.ndarray],
    pca_dim: Optional[int] = None,
) -> Tuple[Dict[str, np.ndarray], List[str]]:
    """Fuse two `{id: vector}` embedding dicts by id-intersection.

    For every shared id, L2-normalize each model's vector separately then
    concatenate (output dim = dim_a + dim_b, or `pca_dim` when given via a
    shared-basis PCA fit on the stacked fused matrix). `dropped` is the
    symmetric-difference of ids (present in only one model) so the caller can
    report — never silently shrink — the orphaned candidates.
    """
    shared = sorted(set(a) & set(b))
    dropped = sorted(set(a) ^ set(b))
    fused = {k: np.concatenate([_l2(a[k]), _l2(b[k])]) for k in shared}
    if pca_dim and fused:
        from sklearn.decomposition import PCA

        keys = list(fused)
        M = np.vstack([fused[k] for k in keys])
        M = PCA(n_components=pca_dim, random_state=0).fit_transform(M)
        fused = {k: M[i] for i, k in enumerate(keys)}
    return fused, dropped


def _parse_tagged(spec: str) -> Tuple[str, str]:
    """Split a `tag:path` CLI argument into (tag, path)."""
    tag, sep, path = spec.partition(":")
    if not sep:
        raise argparse.ArgumentTypeError(f"expected tag:path, got {spec!r}")
    return tag, path


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--pair",
        nargs=2,
        required=True,
        metavar=("tagA:npzA", "tagB:npzB"),
        help="two model embeddings to fuse, each as tag:path-to-npz",
    )
    ap.add_argument("--out-prefix", required=True, help="output written to <out-prefix>.npz")
    ap.add_argument("--pca-dim", type=int, default=None, help="optional shared-basis PCA dim")
    args = ap.parse_args(argv)

    (tag_a, npz_a), (tag_b, npz_b) = (_parse_tagged(s) for s in args.pair)
    emb_a = load_embeddings(npz_a)
    emb_b = load_embeddings(npz_b)
    print(f"loaded {tag_a}: {len(emb_a)} ids from {npz_a}")
    print(f"loaded {tag_b}: {len(emb_b)} ids from {npz_b}")

    fused, dropped = concat_fuse(emb_a, emb_b, pca_dim=args.pca_dim)
    out_npz = f"{args.out_prefix}.npz"
    np.savez(out_npz, **fused)
    print(f"wrote {out_npz}")
    print(f"fused {len(fused)}, dropped {len(dropped)}: {dropped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
