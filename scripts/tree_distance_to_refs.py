#!/usr/bin/env python3
"""tree_distance_to_refs.py — per-candidate patristic distance to nearest reference.

A1 (epic v4bs): the embedding-novelty axis conflates phylogeny with function
(distance-to-centroid ~ phylogenetic distance, Spearman ~0.87). To measure
novelty as *excess divergence beyond phylogenetic expectation*, we residualize
novelty on a phylogenetic confound. This producer computes that confound: for
each candidate leaf in the per-class (class-A) IQ-TREE tree, the minimum
cophenetic (patristic, branch-length) distance to any characterized-reference
leaf. The output TSV (id, tree_distance) is fed to
``fusion_consensus.py --residual-confound tree_distance``.

Pure/tree-only (no torch/GPU); the ete3 import is lazy so the module stays
light for unit collection. Real trees come from stage 04
(``results/phylogenies/protein/*.treefile``); tests exercise tiny Newicks.
"""
from __future__ import annotations

from typing import Dict, Iterable, Optional


def _load_tree(newick: str):
    """Parse a Newick string, tolerating IQ-TREE internal support labels.

    ete3 format 0 reads flexible internal support; format 1 reads internal
    names. Distances use only branch lengths, so either resolves leaves — try
    the permissive format first and fall back.
    """
    from ete3 import Tree  # lazy: keep module import torch/Qt-free

    for fmt in (0, 1):
        try:
            return Tree(newick, format=fmt)
        except Exception:  # noqa: BLE001 — try the next Newick dialect
            continue
    return Tree(newick)  # last resort: let the default raise a clear error


def nearest_ref_distance(
    newick: str,
    ref_ids: Iterable[str],
    candidate_ids: Optional[Iterable[str]] = None,
) -> Dict[str, float]:
    """Min patristic distance from each candidate leaf to any reference leaf.

    ``ref_ids`` are reference leaf names; ``candidate_ids`` restricts the scored
    set (default: every non-reference leaf in the tree). Candidates or
    references absent from the tree are skipped. Returns ``{candidate_id:
    distance}`` (reference leaves are never scored as candidates).
    """
    tree = _load_tree(newick)
    leaves = {leaf.name: leaf for leaf in tree.get_leaves()}

    ref_nodes = [leaves[r] for r in ref_ids if r in leaves]

    if candidate_ids is None:
        ref_name_set = {r for r in ref_ids}
        cand_names = [n for n in leaves if n not in ref_name_set]
    else:
        cand_names = [c for c in candidate_ids if c in leaves and c not in set(ref_ids)]

    out: Dict[str, float] = {}
    if not ref_nodes:
        return out
    for name in cand_names:
        node = leaves[name]
        out[name] = float(min(node.get_distance(r) for r in ref_nodes))
    return out


def main(argv=None) -> None:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--tree", required=True, help="Newick tree file (stage-04 .treefile)")
    ap.add_argument(
        "--ref-ids", required=True,
        help="text file of reference leaf names, one per line",
    )
    ap.add_argument(
        "--candidate-ids", default="",
        help="optional text file of candidate leaf names (default: all non-ref leaves)",
    )
    ap.add_argument("--out", required=True, help="output TSV (id\\ttree_distance)")
    args = ap.parse_args(argv)

    with open(args.tree) as fh:
        newick = fh.read()
    with open(args.ref_ids) as fh:
        ref_ids = [ln.strip() for ln in fh if ln.strip()]
    candidate_ids = None
    if args.candidate_ids:
        with open(args.candidate_ids) as fh:
            candidate_ids = [ln.strip() for ln in fh if ln.strip()]

    dist = nearest_ref_distance(newick, ref_ids, candidate_ids)
    with open(args.out, "w") as fh:
        fh.write("id\ttree_distance\n")
        for cid in sorted(dist):
            fh.write(f"{cid}\t{dist[cid]:.6g}\n")
    print(f"[tree_distance_to_refs] wrote {len(dist)} candidates -> {args.out}")


if __name__ == "__main__":
    main()
