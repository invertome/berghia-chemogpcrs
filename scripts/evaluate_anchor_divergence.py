#!/usr/bin/env python3
"""evaluate_anchor_divergence.py — C3 of the per-class anchor work
(bead berghia-chemogpcrs-521.4).

One-time calibration that decides, per class, whether the OUT-GROUP anchors
(tier 2 Platynereis, tier 3 fly/worm) are safe to keep in the per-class tree, or
whether long-branch attraction pulls them into the focal-species (Berghia)
clades and corrupts the lineage-specific-expansion reading. Tier-1 (in-group
mollusc) anchors are never gated.

Given a without-anchor / with-anchor tree pair per class (built on Unity by the
stage-04 align -> filter -> IQ-TREE machinery), it computes:
  - placement: does an out-group anchor nest inside a well-supported in-group
    clade? (an anchor's smallest clade is supported and made only of in-group
    tips),
  - in-group topology: restricted Robinson-Foulds between the two trees on the
    shared in-group tips (does adding anchors move the in-group around?),
  - support: median drop in in-group clade support (matched bipartitions).

Decision gate (tiers 2-3; thresholds are starting values, re-tune on the pilot):
exclude the out-group anchors for a class if ANY of
  - an anchor nests in a >=support_threshold in-group clade, OR
  - restricted RF > rf_threshold, OR
  - median in-group support drops by > support_drop_threshold.

Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from pathlib import Path
from typing import Optional

from ete3 import Tree

DEFAULT_SUPPORT_THRESHOLD = 80.0   # UFBoot >=80 (use 0.80 for TBE-scaled trees)
DEFAULT_RF_THRESHOLD = 0.10
DEFAULT_SUPPORT_DROP_THRESHOLD = 5.0

BERGHIA_TAXID = 1287507
MOLLUSCA_TAXID = 6447


# ---------------------------------------------------------------------------
# taxonomy + pool-membership manifest
# ---------------------------------------------------------------------------

_NCBI_TAXA = None


def get_lineage(taxid: int) -> list:
    """NCBI lineage (root→leaf) via ete3.NCBITaxa (lazy-init). Replaced in tests."""
    global _NCBI_TAXA
    if _NCBI_TAXA is None:
        try:
            from ete3 import NCBITaxa
            _NCBI_TAXA = NCBITaxa()
        except Exception:
            _NCBI_TAXA = False
    if not _NCBI_TAXA:
        return [taxid]
    try:
        return _NCBI_TAXA.get_lineage(taxid) or [taxid]
    except Exception:
        return [taxid]


def is_mollusc(taxid: int, lineage_fn=get_lineage) -> bool:
    """True if *taxid* is within Mollusca (6447)."""
    return MOLLUSCA_TAXID in set(lineage_fn(taxid))


def load_pool_labels(members_tsv: str, *, is_mollusc_fn=is_mollusc,
                     berghia_taxid: int = BERGHIA_TAXID) -> tuple[set, set]:
    """Read a pool-membership manifest (seq_id, taxid, source) and return
    (ingroup_labels, berghia_labels):
      - berghia_labels = leaves whose source is 'berghia',
      - ingroup_labels = Berghia + mollusc-reference leaves (anchors excluded).
    """
    ingroup: set[str] = set()
    berghia: set[str] = set()
    with open(members_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            seq_id = row["seq_id"]
            source = row.get("source", "")
            try:
                taxid = int(row["taxid"])
            except (KeyError, ValueError):
                taxid = 0
            if source == "berghia":
                berghia.add(seq_id)
                ingroup.add(seq_id)
            elif source == "ref" and is_mollusc_fn(taxid):
                ingroup.add(seq_id)
    return ingroup, berghia


# ---------------------------------------------------------------------------
# anchor labels
# ---------------------------------------------------------------------------

def parse_anchor_tier(label: str) -> Optional[str]:
    """Return the tier ('1'/'2'/'3') from an ANCHOR_<class>_<tier>_<acc> label,
    or None if the label is not an anchor."""
    parts = (label or "").split("_")
    if len(parts) >= 4 and parts[0] == "ANCHOR":
        return parts[2]
    return None


def is_outgroup_anchor(label: str) -> bool:
    """True for tier-2/3 (out-group) anchors only — the ones subject to the gate."""
    return parse_anchor_tier(label) in {"2", "3"}


# ---------------------------------------------------------------------------
# metrics
# ---------------------------------------------------------------------------

def anchor_infiltrations(tree: Tree, ingroup_labels: set,
                         support_threshold: float = DEFAULT_SUPPORT_THRESHOLD) -> list:
    """Out-group anchors whose smallest clade is well-supported and composed
    only of in-group tips — i.e. the anchor is confidently nested *inside* an
    in-group (Berghia/mollusc) clade rather than sitting outside it."""
    flagged: list[str] = []
    for leaf in tree.get_leaves():
        if not is_outgroup_anchor(leaf.name):
            continue
        parent = leaf.up
        if parent is None:
            continue
        others = [n for n in parent.get_leaf_names() if n != leaf.name]
        if (others and parent.support >= support_threshold
                and all(o in ingroup_labels for o in others)):
            flagged.append(leaf.name)
    return sorted(set(flagged))


def _prune_copy(tree: Tree, labels: set) -> Optional[Tree]:
    """Return a copy of *tree* pruned to the leaves in *labels* present in it."""
    present = [l for l in tree.get_leaf_names() if l in labels]
    if len(present) < 3:
        return None
    t = tree.copy()
    t.prune(present, preserve_branch_length=True)
    return t


def restricted_rf(tree_without: Tree, tree_with: Tree, shared_labels: set) -> float:
    """Normalized Robinson-Foulds between the two trees on their shared in-group
    tips (0 = identical in-group topology)."""
    a = _prune_copy(tree_without, shared_labels)
    b = _prune_copy(tree_with, shared_labels)
    if a is None or b is None:
        return 0.0
    res = a.robinson_foulds(b, unrooted_trees=True)
    rf, max_rf = res[0], res[1]
    return (rf / max_rf) if max_rf else 0.0


def _ingroup_bipartition_support(tree: Tree, ingroup_labels: set) -> dict:
    """Map frozenset(clade leaves) → support for each non-trivial in-group clade
    of *tree* pruned to the in-group tips."""
    t = _prune_copy(tree, ingroup_labels)
    if t is None:
        return {}
    all_leaves = frozenset(t.get_leaf_names())
    out: dict[frozenset, float] = {}
    for node in t.traverse():
        if node.is_leaf() or node.is_root():
            continue
        clade = frozenset(node.get_leaf_names())
        if 2 <= len(clade) < len(all_leaves):
            out[clade] = node.support
    return out


def ingroup_support_drop(tree_without: Tree, tree_with: Tree,
                         ingroup_labels: set) -> float:
    """Median support drop (without − with) over in-group clades (bipartitions)
    present in BOTH trees. Positive = anchors weakened in-group support."""
    a = _ingroup_bipartition_support(tree_without, ingroup_labels)
    b = _ingroup_bipartition_support(tree_with, ingroup_labels)
    shared = set(a) & set(b)
    if not shared:
        return 0.0
    deltas = [a[k] - b[k] for k in shared]
    return statistics.median(deltas)


# ---------------------------------------------------------------------------
# verdict
# ---------------------------------------------------------------------------

def evaluate_class(tree_without: Tree, tree_with: Tree, ingroup_labels: set, *,
                   berghia_labels: set = None,
                   support_threshold: float = DEFAULT_SUPPORT_THRESHOLD,
                   rf_threshold: float = DEFAULT_RF_THRESHOLD,
                   support_drop_threshold: float = DEFAULT_SUPPORT_DROP_THRESHOLD,
                   ) -> dict:
    """Decide whether out-group anchors are safe for this class. Returns a metrics
    dict with an include/exclude verdict. Tier-1 anchors are never gated.

    Placement infiltration is checked against *berghia_labels* (the gate is about
    out-group anchors nesting in Berghia clades, spec §7); RF and support use
    *ingroup_labels* (Berghia + mollusc, spec §6). berghia_labels defaults to
    ingroup_labels when not supplied.
    """
    if berghia_labels is None:
        berghia_labels = ingroup_labels
    infiltrations = anchor_infiltrations(tree_with, berghia_labels, support_threshold)
    rf = restricted_rf(tree_without, tree_with, ingroup_labels)
    drop = ingroup_support_drop(tree_without, tree_with, ingroup_labels)

    reasons = []
    if infiltrations:
        reasons.append(f"{len(infiltrations)} anchor(s) nest in a supported in-group clade")
    if rf > rf_threshold:
        reasons.append(f"in-group RF {rf:.3f} > {rf_threshold}")
    if drop > support_drop_threshold:
        reasons.append(f"in-group support drop {drop:.1f} > {support_drop_threshold}")

    return {
        "verdict": "exclude" if reasons else "include",
        "reasons": reasons,
        "n_infiltrations": len(infiltrations),
        "infiltrating_anchors": infiltrations,
        "rf": rf,
        "support_drop": drop,
        "thresholds": {
            "support": support_threshold,
            "rf": rf_threshold,
            "support_drop": support_drop_threshold,
        },
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _load_labels(path: str) -> set:
    with open(path) as fh:
        return {ln.strip() for ln in fh if ln.strip()}


def build_args_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Evaluate whether out-group anchors are safe to keep in a "
                    "per-class tree (with vs without anchors).")
    p.add_argument("--class", dest="klass", required=True, help="GPCR class label (A/B/C/F)")
    p.add_argument("--tree-without", required=True, help="Newick tree built WITHOUT out-group anchors")
    p.add_argument("--tree-with", required=True, help="Newick tree built WITH out-group anchors")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--pool-members",
                     help="pool_members_class_<X>.tsv manifest — derives in-group "
                          "(Berghia + mollusc) and Berghia leaf sets by taxid")
    src.add_argument("--ingroup-labels",
                     help="File of in-group leaf labels, one per line "
                          "(Berghia labels default to the same set)")
    p.add_argument("--out", required=True, help="Output verdict+metrics JSON path")
    p.add_argument("--support-threshold", type=float, default=DEFAULT_SUPPORT_THRESHOLD)
    p.add_argument("--rf-threshold", type=float, default=DEFAULT_RF_THRESHOLD)
    p.add_argument("--support-drop-threshold", type=float, default=DEFAULT_SUPPORT_DROP_THRESHOLD)
    return p


def main(argv=None) -> None:
    args = build_args_parser().parse_args(argv)
    without = Tree(args.tree_without, format=0)
    with_anchor = Tree(args.tree_with, format=0)
    if args.pool_members:
        ingroup, berghia = load_pool_labels(args.pool_members)
    else:
        ingroup = _load_labels(args.ingroup_labels)
        berghia = None
    res = evaluate_class(
        without, with_anchor, ingroup,
        berghia_labels=berghia,
        support_threshold=args.support_threshold,
        rf_threshold=args.rf_threshold,
        support_drop_threshold=args.support_drop_threshold,
    )
    res["class"] = args.klass
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(res, fh, indent=2)
    print(f"[evaluate_anchor_divergence] class {args.klass}: {res['verdict']} "
          f"(infiltrations={res['n_infiltrations']}, rf={res['rf']:.3f}, "
          f"support_drop={res['support_drop']:.1f}) -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
