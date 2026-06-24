#!/usr/bin/env python3
"""validate_landmark_placement.py — one-time pilot study: place dropped tier-2/3
landmarks on the clean per-class trees (EPA-ng) and run three directional
comparisons. See docs/plans/2026-06-24-landmark-placement-validation-design.md.

Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts
"""
from __future__ import annotations
import argparse
import csv
import json
import os
import sys
import warnings
from collections import Counter
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import classify_via_placement as cvp


def normalize_family(f) -> str:
    t = (f or "").strip().lower()
    return t if t else "?unknown?"


def _read_fasta(path: str) -> dict:
    seqs, cur = {}, None
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                cur = ln[1:].split()[0]
                seqs[cur] = ""
            elif cur is not None:
                seqs[cur] += ln.strip()
    return seqs


def load_landmarks(anchor_tsv: str, anchor_fasta: str, klass: str) -> list:
    """tier-2/3 anchors for `klass`: [{id, family, seq}]. id == ANCHOR_<class>_<tier>_<acc>."""
    seqs = _read_fasta(anchor_fasta)
    out = []
    with open(anchor_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("class") != klass or row.get("tier") not in {"2", "3"}:
                continue
            lid = f"ANCHOR_{klass}_{row['tier']}_{row['accession']}"
            rec = {"id": lid, "family": row.get("family", ""), "seq": seqs.get(lid, "")}
            if not rec["seq"]:
                warnings.warn(f"load_landmarks: no sequence for {lid}")
            out.append(rec)
    return out


def run_epang_placement(ref_msa: str, ref_tree: str, query_fasta: str,
                        outdir: str, threads: int = 4) -> str:
    """Align query_fasta onto ref_msa (mafft --add), split, run EPA-ng.
    Returns path to the epa_result.jplace file.
    Reuses helpers from scripts/classify_via_placement.py.
    Intermediate alignments (combined.aln, query.aln, ref.aln) are written
    into outdir and retained for inspection (not accidental leakage)."""
    Path(outdir).mkdir(parents=True, exist_ok=True)
    combined = os.path.join(outdir, "combined.aln")
    query_aln = os.path.join(outdir, "query.aln")
    ref_aln = os.path.join(outdir, "ref.aln")
    if not cvp._run_mafft_add(query_fasta, ref_msa, combined, threads):
        raise RuntimeError(f"mafft --add failed for {query_fasta}")
    ref_ids = set(_read_fasta(ref_msa).keys())
    cvp._split_query_ref(combined, ref_ids, query_aln, ref_aln)
    if not cvp._run_epa_ng(query_aln, ref_aln, ref_tree, outdir, threads):
        raise RuntimeError(f"epa-ng failed (query={query_fasta}, tree={ref_tree}, outdir={outdir})")
    return os.path.join(outdir, "epa_result.jplace")


def place_landmarks(ref_msa: str, ref_tree: str, landmarks: list,
                    outdir: str, threads: int = 4) -> tuple:
    """Place landmarks on the reference tree via EPA-ng.

    Returns (placements, edge_to_leaves) where:
      placements = {landmark_id: {"edge": int, "lwr": float}}
      edge_to_leaves = {edge_num: [leaf_name, ...]}
    """
    Path(outdir).mkdir(parents=True, exist_ok=True)
    qfa = os.path.join(outdir, "landmarks_query.fa")
    with open(qfa, "w") as fh:
        for lm in landmarks:
            fh.write(f">{lm['id']}\n{lm['seq']}\n")
    jplace = run_epang_placement(ref_msa, ref_tree, qfa, outdir, threads)
    raw = cvp.parse_jplace(jplace)
    placements = {k: {"edge": v["edge_num"], "lwr": v["lwr"]} for k, v in raw.items()}
    with open(jplace) as fh:
        tree = json.load(fh)["tree"]
    edge_to_leaves = cvp._parse_newick_with_edge_labels(tree)
    return placements, edge_to_leaves


def _placement_node(tree, edge_to_leaves, edge):
    """The reference-tree node a landmark sits on: MRCA of the leaves under its
    EPA-ng placement edge (or the single leaf if the edge subtends one)."""
    leaves = edge_to_leaves.get(edge, [])
    present = [leaf for leaf in leaves if leaf in set(tree.get_leaf_names())]
    present = list(dict.fromkeys(present))
    if not present:
        return None
    if len(present) == 1:
        return tree.search_nodes(name=present[0])[0]
    return tree.get_common_ancestor(present)


def nearest_landmark_per_candidate(tree, placements, landmark_family, berghia_ids,
                                   lwr_min=0.80, edge_to_leaves=None):
    """For each Berghia candidate, find the nearest LWR-qualifying landmark.

    Returns a list of {candidate, family, landmark, distance, lwr}, one per
    Berghia id that has at least one reachable qualifying landmark.
    Ties (equidistant landmarks) broken by landmark id (alphabetical), for
    reproducibility.
    """
    edge_to_leaves = edge_to_leaves or {}
    nodes = {lid: _placement_node(tree, edge_to_leaves, p["edge"])
             for lid, p in placements.items() if p["lwr"] >= lwr_min}
    nodes = {k: v for k, v in nodes.items() if v is not None}
    rows = []
    leafset = set(tree.get_leaf_names())
    for bid in berghia_ids:
        if bid not in leafset:
            continue
        bleaf = tree.search_nodes(name=bid)[0]
        best = min(((lid, tree.get_distance(bleaf, n, topology_only=True))
                    for lid, n in nodes.items()), key=lambda x: (x[1], x[0]), default=None)
        if best is None:
            continue
        lid, dist = best
        rows.append({"candidate": bid,
                     "family": normalize_family(landmark_family.get(lid, "")),
                     "landmark": lid, "distance": dist,
                     "lwr": placements[lid]["lwr"]})
    return rows


def axis1_vs_classifier(per_candidate_rows, class_berghia_tsv):
    """Per-candidate placement family vs classifier evidence_family_hmm.
    Candidates absent from the baseline are scored against '?unknown?'.
    confusion keys are (placement_family, classifier_family) TUPLES; callers
    must flatten (e.g. 'plc|clf') before JSON."""
    base = {}
    with open(class_berghia_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or "seq_id" not in reader.fieldnames:
            raise KeyError(
                f"axis1_vs_classifier: expected 'seq_id' column in "
                f"{class_berghia_tsv!r}; got {reader.fieldnames}")
        for row in reader:
            base[row["seq_id"]] = normalize_family(row.get("evidence_family_hmm", ""))
    n = 0
    n_agree = 0
    conf = Counter()
    disc = []
    for r in per_candidate_rows:
        plc = normalize_family(r["family"])
        clf = base.get(r["candidate"], "?unknown?")
        n += 1
        conf[(plc, clf)] += 1
        if plc == clf and plc != "?unknown?":
            n_agree += 1
        else:
            disc.append({"candidate": r["candidate"], "placement": plc, "classifier": clf})
    return {"n": n, "n_agree": n_agree,
            "concordance": (n_agree / n if n else 0.0),
            "confusion": dict(conf), "discordant": disc}


def _jaccard(a, b):
    sa, sb = set(a), set(b)
    union = sa | sb
    if not union:
        return 1.0
    return len(sa & sb) / len(union)


def axis2_position_concordance(with_tree, placements, edge_to_leaves,
                               landmark_ids, berghia_ids,
                               jaccard_min=0.5, supp_min=80):
    """Directional axis-2 comparison. For each tier-2/3 landmark, compare its
    Berghia neighborhood in the WITH tree (sisters of the landmark leaf) vs the
    Berghia in its placement edge's pendant subtree on the clean tree
    (edge_to_leaves). infiltrating_in_with_tree mirrors C3's anchor_infiltrations.
    No clean_tree arg: edge_to_leaves already encodes the clean-tree adjacency,
    and placed_berghia is that pendant-subtree proxy."""
    bset = set(berghia_ids)
    with_leaves = set(with_tree.get_leaf_names())
    out = []
    for lid in landmark_ids:
        with_berghia = set()
        infiltrating = False
        if lid in with_leaves:
            parent = with_tree.search_nodes(name=lid)[0].up
            if parent is not None:
                sisters = [n for n in parent.get_leaf_names() if n != lid]
                with_berghia = set(sisters) & bset
                infiltrating = bool(sisters and parent.support >= supp_min
                                    and all(s in bset for s in sisters))
        placed_berghia = set()
        # Limitation: placed_berghia is the pendant (descendant-only) subtree of the
        # placement edge; for edges near internal nodes with many Berghia on the
        # other side this undercounts Berghia adjacency and biases Jaccard downward
        # (can show reproduced=False for a landmark that truly lands in the Berghia region).
        if lid in placements:
            edge = placements[lid]["edge"]
            placed_berghia = set(edge_to_leaves.get(edge, [])) & bset
        jac = _jaccard(with_berghia, placed_berghia)
        out.append({"landmark": lid,
                    "with_berghia": sorted(with_berghia),
                    "placed_berghia": sorted(placed_berghia),
                    "jaccard": jac,
                    "reproduced": jac >= jaccard_min,
                    "infiltrating_in_with_tree": infiltrating})
    return out


def clade_consensus(tree, per_candidate_rows, berghia_ids, supp_min=80):
    """Consensus family for each supported all-Berghia internal clade.

    An LSE clade is an internal node with support >= supp_min whose leaf set
    is a subset of berghia_ids only (no reference taxa) and has >= 2 tips.
    Consensus = majority family among clade members present in per_candidate_rows;
    ties broken by lowest-sorted family name (deterministic).

    Returns a list of {clade_leaves, consensus_family, n_called}.
    Clades with no called members are skipped.
    """
    fam_by_cand = {r["candidate"]: normalize_family(r["family"]) for r in per_candidate_rows}
    bset = set(berghia_ids)
    out = []
    for node in tree.traverse():
        if node.is_leaf():
            continue
        leaves = set(node.get_leaf_names())
        berg = leaves & bset
        if len(berg) < 2 or leaves != berg or node.support < supp_min:
            continue
        fams = [fam_by_cand[b] for b in berg if b in fam_by_cand]
        if not fams:
            continue
        cnt = Counter(fams)
        top = sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        out.append({"clade_leaves": sorted(berg), "consensus_family": top, "n_called": len(fams)})
    return out
