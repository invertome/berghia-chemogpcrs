#!/usr/bin/env python3
"""placement_distance.py — distance-to-nearest-reference from EPA-ng placements.

Bead amiu: the GENERALIZABLE form of the A1 tree-distance confound. Building a
de-novo tree of {candidates + anchors} per species does not scale — only the
ANCHOR set is species-independent. So the anchor backbone is built ONCE
(scripts/unity/build_a1_backbone.sh) and each species' candidates are PLACED on
it with EPA-ng; this module reads that placement (.jplace) and emits each
candidate's patristic distance to the nearest reference leaf.

Output schema is deliberately identical to tree_distance_to_refs.py
(``id``/``tree_distance``), so it is a drop-in source for
``fusion_consensus.py --confound-tsv tree_distance:<path>``:

    de-novo tree  -> tree_distance_to_refs.py   (single-species one-off)
    placement     -> placement_distance.py      (any-species pipeline)

Geometry. A placement sits on the backbone edge ``edge_num`` at
``distal_length`` d from that edge's DISTAL (away-from-root) node v, with a
``pendant_length`` p branch to the query. With ``below[v]`` = nearest reference
leaf below v and ``above[v]`` = nearest reference leaf reached by going up and
around from v, the distance is::

    p + min( d + below[v],  above[v] - d )

``below``/``above`` come from one post-order and one pre-order pass, so each
placement is O(1). When several placements are reported for a query the one with
the highest ``like_weight_ratio`` is used.
"""
from __future__ import annotations

import re
from typing import Dict, Iterable, Optional, Tuple

INF = float("inf")

# `NAME:LEN{EDGE}` -> `NAME@@EDGE:LEN`, so the edge number survives Newick
# parsing as part of the node name (no reliance on traversal order).
_EDGE_LABEL_RE = re.compile(r"([^(),:;\[\]]*):([0-9.eE+\-]+)\{(\d+)\}")


def load_backbone(newick: str) -> Tuple[object, Dict[int, object]]:
    """Parse a jplace reference tree into ``(tree, {edge_num: node})``.

    Edge ``N`` is the edge leading to the returned node (its parent edge), which
    is the jplace convention.
    """
    from ete3 import Tree  # lazy: keep module import light

    clean = _EDGE_LABEL_RE.sub(r"\1@@\3:\2", newick)
    tree = Tree(clean, format=1)

    edge_map: Dict[int, object] = {}
    for node in tree.traverse():
        if "@@" in node.name:
            name, _, edge = node.name.rpartition("@@")
            node.name = name
            edge_map[int(edge)] = node
    return tree, edge_map


def _below_above(tree, ref_ids: Optional[Iterable[str]] = None):
    """Nearest-reference-leaf distance below each node, and up-and-around it."""
    refs = None if ref_ids is None else set(ref_ids)

    def _is_ref(node) -> bool:
        return node.is_leaf() and (refs is None or node.name in refs)

    below: Dict[object, float] = {}
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            below[node] = 0.0 if _is_ref(node) else INF
        else:
            below[node] = min((c.dist + below[c] for c in node.children), default=INF)

    above: Dict[object, float] = {tree: INF}
    for node in tree.traverse("preorder"):
        for child in node.children:
            sib_best = min(
                (s.dist + below[s] for s in node.children if s is not child),
                default=INF,
            )
            above[child] = child.dist + min(above[node], sib_best)
    return below, above


def _query_names(placement: dict):
    """jplace allows ``n: [names]`` or ``nm: [[name, multiplicity], ...]``."""
    if "n" in placement:
        n = placement["n"]
        return list(n) if isinstance(n, (list, tuple)) else [n]
    return [entry[0] for entry in placement.get("nm", [])]


def nearest_ref_distance_from_jplace(
    jplace: dict,
    ref_ids: Optional[Iterable[str]] = None,
) -> Dict[str, float]:
    """Per query: patristic distance to the nearest reference leaf.

    ``ref_ids`` restricts which backbone leaves count as references (default:
    every leaf, which is correct for an anchors-only backbone). Queries with no
    placement, or with no reachable reference, are omitted rather than given a
    sentinel value.
    """
    tree, edge_map = load_backbone(jplace["tree"])
    below, above = _below_above(tree, ref_ids)

    fields = list(jplace.get("fields", []))
    idx = {name: fields.index(name) for name in fields}
    i_edge = idx["edge_num"]
    i_distal = idx["distal_length"]
    i_pendant = idx["pendant_length"]
    i_lwr = idx.get("like_weight_ratio")

    out: Dict[str, float] = {}
    for placement in jplace.get("placements", []):
        rows = placement.get("p") or []
        if not rows:
            continue
        best = max(rows, key=lambda r: r[i_lwr]) if i_lwr is not None else rows[0]
        node = edge_map.get(int(best[i_edge]))
        if node is None:
            continue
        d = float(best[i_distal])
        pendant = float(best[i_pendant])
        reach = min(d + below[node], above[node] - d)
        if not (reach < INF):
            continue                      # no reference reachable -> omit
        for name in _query_names(placement):
            out[name] = float(pendant + reach)
    return out


def main(argv=None) -> None:
    import argparse
    import json

    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--jplace", required=True, help="EPA-ng placement file (.jplace)")
    ap.add_argument("--ref-ids", default="",
                    help="optional file of reference leaf names, one per line "
                         "(default: every backbone leaf counts as a reference)")
    ap.add_argument("--out", required=True, help="output TSV (id\\ttree_distance)")
    args = ap.parse_args(argv)

    with open(args.jplace) as fh:
        jplace = json.load(fh)
    ref_ids = None
    if args.ref_ids:
        with open(args.ref_ids) as fh:
            ref_ids = [ln.strip() for ln in fh if ln.strip()]

    dist = nearest_ref_distance_from_jplace(jplace, ref_ids)
    with open(args.out, "w") as fh:
        fh.write("id\ttree_distance\n")
        for cid in sorted(dist):
            fh.write(f"{cid}\t{dist[cid]:.6g}\n")
    print(f"[placement_distance] wrote {len(dist)} candidates -> {args.out}")


if __name__ == "__main__":
    main()
