#!/usr/bin/env python3
"""validate_landmark_placement.py — one-time pilot study: place dropped tier-2/3
landmarks on the clean per-class trees (EPA-ng) and run three directional
comparisons. See docs/plans/2026-06-24-landmark-placement-validation-design.md.

Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts
"""
from __future__ import annotations
import argparse
import contextlib
import csv
import json
import os
import sys
import warnings
from collections import Counter
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import classify_via_placement as cvp
import evaluate_anchor_divergence as ead
from _classification_labels import label_for_hmm


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
    must flatten (e.g. 'plc|clf') before JSON.

    Both sides are coarse-collapsed via label_for_hmm before normalize_family:
    medium HMM tokens such as 'aminergic_5HT' are reduced to their coarse
    family 'aminergic' (split on first underscore) before comparison.
    Anchor families never contain underscores so label_for_hmm is a no-op
    on the placement side; applying it to both keeps the comparison symmetric."""
    base = {}
    with open(class_berghia_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or "seq_id" not in reader.fieldnames:
            raise KeyError(
                f"axis1_vs_classifier: expected 'seq_id' column in "
                f"{class_berghia_tsv!r}; got {reader.fieldnames}")
        for row in reader:
            base[row["seq_id"]] = normalize_family(label_for_hmm(row.get("evidence_family_hmm", ""))[0])
    n = 0
    n_agree = 0
    conf = Counter()
    disc = []
    for r in per_candidate_rows:
        plc = normalize_family(label_for_hmm(r["family"])[0])
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


def _prune_to(tree, labels):
    present = [leaf for leaf in tree.get_leaf_names() if leaf in labels]
    if len(present) < 3:
        return None
    t = tree.copy()
    t.prune(present, preserve_branch_length=True)
    return t


def _supported_clades(tree, supp_min):
    all_leaves = frozenset(tree.get_leaf_names())
    out = set()
    for node in tree.traverse():
        if node.is_leaf() or node.is_root():
            continue
        clade = frozenset(node.get_leaf_names())
        if 2 <= len(clade) < len(all_leaves) and node.support >= supp_min:
            out.add(clade)
    return out


def axis3_lse_robustness(with_tree, clean_tree, berghia_ids, supp_min=80):
    """Berghia-restricted RF + shared supported-clade fraction between the
    with-anchor and clean trees. Low rf = anchors didn't perturb the LSE
    structure (drop harmless); high rf = they did (clean tree is LBA-free, so
    that supports the drop). Clean tree is the shared-fraction denominator."""
    neutral = {"rf": 0.0, "max_rf": 0, "n_shared_clades": 0,
               "n_clean_clades": 0, "shared_fraction": 1.0}
    bset = set(berghia_ids)
    a = _prune_to(with_tree, bset)
    b = _prune_to(clean_tree, bset)
    if a is None or b is None:
        return neutral
    shared = set(a.get_leaf_names()) & set(b.get_leaf_names())
    a = _prune_to(a, shared)
    b = _prune_to(b, shared)
    if a is None or b is None:
        return neutral
    rf, max_rf = a.robinson_foulds(b, unrooted_trees=True)[:2]
    norm_rf = (rf / max_rf) if max_rf else 0.0
    with_clades = _supported_clades(a, supp_min)
    clean_clades = _supported_clades(b, supp_min)
    n_shared = len(with_clades & clean_clades)
    n_clean = len(clean_clades)
    return {"rf": norm_rf, "max_rf": max_rf, "n_shared_clades": n_shared,
            "n_clean_clades": n_clean,
            "shared_fraction": (n_shared / n_clean) if n_clean else 1.0}


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


# ---------------------------------------------------------------------------
# Task 6: per-class orchestration + report writer + CLI
# ---------------------------------------------------------------------------


def _flatten_confusion(confusion: dict) -> dict:
    """Flatten (placement_family, classifier_family) tuple keys to 'plc|clf'
    strings so the confusion matrix is JSON/TSV/markdown safe (never a raw
    Python tuple)."""
    return {f"{k[0]}|{k[1]}": v for k, v in confusion.items()}


def _load_berghia_ids(pool_members_tsv: str) -> list:
    """Berghia seq_ids from a pool_members_class_<X>.tsv (columns seq_id, taxid,
    source); a Berghia tip is a row with source == 'berghia'."""
    out = []
    with open(pool_members_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("source") == "berghia":
                out.append(row["seq_id"])
    return out


def run_class(klass, calibration_dir, anchor_tsv, anchor_fasta, class_berghia_tsv,
              outdir, threads=4, lwr_min=0.80, supp_min=80, jaccard_min=0.5):
    """Wire Tasks 1-5 for one class and return its per-class result dict.

    This is the only function that calls real tools (mafft/epa-ng via
    place_landmarks); it is exercised for real in Task 8 (not unit-tested).
    Trees carry IQ-TREE SH-aLRT/UFBoot supports, so they are loaded through the
    support-normalizing loader (evaluate_anchor_divergence.load_tree).
    """
    trees_dir = os.path.join(calibration_dir, "trees_nocloak")
    clean_tree_path = os.path.join(trees_dir, f"without_class_{klass}.treefile")
    with_tree_path = os.path.join(trees_dir, f"with_class_{klass}.treefile")
    clean_msa_candidates = [
        os.path.join(trees_dir, f"without_class_{klass}_trimmed.fa"),
        os.path.join(trees_dir, f"without_class_{klass}_aligned.fa"),
        os.path.join(trees_dir, f"_fs_without_class_{klass}", f"without_class_{klass}_tapered.fa"),
    ]
    clean_msa = next((p for p in clean_msa_candidates if os.path.exists(p)), None)
    if clean_msa is None:
        raise FileNotFoundError(
            f"run_class: no clean-tree MSA for class {klass}; tried "
            + ", ".join(clean_msa_candidates))

    pool_members = os.path.join(calibration_dir, "pools", f"pool_members_class_{klass}.tsv")
    berghia_ids = _load_berghia_ids(pool_members)

    clean_tree = ead.load_tree(clean_tree_path)
    with_tree = ead.load_tree(with_tree_path)

    landmarks = load_landmarks(anchor_tsv, anchor_fasta, klass)
    placements, edge_to_leaves = place_landmarks(
        clean_msa, clean_tree_path, landmarks,
        os.path.join(outdir, f"class_{klass}"), threads)
    landmark_family = {lm["id"]: lm["family"] for lm in landmarks}
    landmark_ids = [lm["id"] for lm in landmarks]

    per_candidate = nearest_landmark_per_candidate(
        clean_tree, placements, landmark_family, berghia_ids, lwr_min, edge_to_leaves)
    clade = clade_consensus(clean_tree, per_candidate, berghia_ids, supp_min)
    axis1 = axis1_vs_classifier(per_candidate, class_berghia_tsv)
    axis2 = axis2_position_concordance(
        with_tree, placements, edge_to_leaves, landmark_ids, berghia_ids, jaccard_min, supp_min)
    axis3 = axis3_lse_robustness(with_tree, clean_tree, berghia_ids, supp_min)

    return {"class": klass, "per_candidate": per_candidate, "clade_consensus": clade,
            "axis1": axis1, "axis2": axis2, "axis3": axis3}


def write_report(outdir, per_class_results):
    """Emit the 6 output files (5 TSVs + markdown) summarising the study,
    aggregating rows across classes with a leading `class` column. Confusion
    tuple keys are flattened to 'plc|clf' before writing."""
    od = Path(outdir)
    od.mkdir(parents=True, exist_ok=True)

    with contextlib.ExitStack() as stack:
        def _writer(name, header):
            fh = stack.enter_context(open(od / name, "w", newline=""))
            w = csv.writer(fh, delimiter="\t")
            w.writerow(header)
            return w

        pc_w = _writer("per_candidate_calls.tsv",
                       ["class", "candidate", "family", "landmark", "distance", "lwr"])
        cc_w = _writer("clade_consensus.tsv",
                       ["class", "clade_leaves", "consensus_family", "n_called"])
        a1_w = _writer("axis1_concordance.tsv",
                       ["class", "n", "n_agree", "concordance"])
        a2_w = _writer("axis2_position.tsv",
                       ["class", "landmark", "with_berghia", "placed_berghia",
                        "jaccard", "reproduced", "infiltrating_in_with_tree"])
        a3_w = _writer("axis3_lse_robustness.tsv",
                       ["class", "rf", "max_rf", "n_shared_clades",
                        "n_clean_clades", "shared_fraction"])

        for res in per_class_results:
            klass = res["class"]
            for r in res["per_candidate"]:
                pc_w.writerow([klass, r["candidate"], r["family"], r["landmark"],
                               r["distance"], r["lwr"]])
            for c in res["clade_consensus"]:
                cc_w.writerow([klass, ";".join(c["clade_leaves"]),
                               c["consensus_family"], c["n_called"]])
            a1 = res["axis1"]
            a1_w.writerow([klass, a1["n"], a1["n_agree"], a1["concordance"]])
            for r in res["axis2"]:
                a2_w.writerow([klass, r["landmark"], ";".join(r["with_berghia"]),
                               ";".join(r["placed_berghia"]), r["jaccard"],
                               r["reproduced"], r["infiltrating_in_with_tree"]])
            a3 = res["axis3"]
            a3_w.writerow([klass, a3["rf"], a3["max_rf"], a3["n_shared_clades"],
                           a3["n_clean_clades"], a3["shared_fraction"]])

    _write_markdown(od / "landmark_placement_report.md", per_class_results)


def _write_markdown(path, per_class_results):
    """Markdown narrative: per-class summary + a section per axis with the
    directional-interpretation sentences."""
    lines = []
    lines.append("# Landmark-placement validation report")
    lines.append("")
    lines.append("One-time pilot: tier-2/3 out-group landmarks (dropped by the C3 "
                 "anchor-divergence gate) are placed back onto the clean per-class "
                 "trees with EPA-ng, then read along three directional axes.")
    lines.append("")

    lines.append("## Per-class summary")
    lines.append("")
    lines.append("| class | candidates called | LSE clades called | axis-1 concordance | axis-3 rf | axis-3 shared_fraction |")
    lines.append("|---|---|---|---|---|---|")
    for res in per_class_results:
        a1 = res["axis1"]
        a3 = res["axis3"]
        lines.append(f"| {res['class']} | {len(res['per_candidate'])} | "
                     f"{len(res['clade_consensus'])} | {a1['concordance']:.3f} | "
                     f"{a3['rf']:.3f} | {a3['shared_fraction']:.3f} |")
    lines.append("")

    lines.append("## Axis 1 — functional read vs classifier")
    lines.append("")
    lines.append("Per-candidate placement family vs the classifier's "
                 "`evidence_family_hmm`. High concordance means the placement-based "
                 "functional read agrees with the independent HMM classifier.")
    lines.append("")
    for res in per_class_results:
        a1 = res["axis1"]
        lines.append(f"### Class {res['class']} — concordance "
                     f"{a1['concordance']:.3f} ({a1['n_agree']}/{a1['n']})")
        lines.append("")
        lines.append("Confusion (placement|classifier -> count):")
        lines.append("")
        lines.append("| placement|classifier | count |")
        lines.append("|---|---|")
        for k, v in sorted(_flatten_confusion(a1["confusion"]).items()):
            lines.append(f"| {k} | {v} |")
        lines.append("")
        if a1["discordant"]:
            lines.append("Discordant candidates:")
            lines.append("")
            for d in a1["discordant"]:
                lines.append(f"- {d['candidate']}: placement={d['placement']}, "
                             f"classifier={d['classifier']}")
            lines.append("")
        else:
            lines.append("No discordant candidates.")
            lines.append("")

    lines.append("## Axis 2 — placement vs in-inference position")
    lines.append("")
    lines.append("Directional: where a landmark is soundly placed we expect the "
                 "EPA-ng position to *reproduce* its in-inference (with-anchor) "
                 "Berghia neighborhood; where the with-anchor tree shows the "
                 "landmark infiltrating an all-Berghia clade (long-branch "
                 "attraction), we expect placement to NOT reproduce that "
                 "infiltration — i.e. `reproduced=False` for an infiltrating "
                 "landmark is the desired outcome and supports the drop.")
    lines.append("")
    for res in per_class_results:
        rows = res["axis2"]
        n_repro = sum(1 for r in rows if r["reproduced"])
        n_infil = sum(1 for r in rows if r["infiltrating_in_with_tree"])
        lines.append(f"### Class {res['class']} — {len(rows)} landmarks; "
                     f"{n_repro} reproduced, {n_infil} infiltrating in with-tree")
        lines.append("")
        lines.append("| landmark | with_berghia | placed_berghia | jaccard | reproduced | infiltrating_in_with_tree |")
        lines.append("|---|---|---|---|---|---|")
        for r in rows:
            lines.append(f"| {r['landmark']} | {';'.join(r['with_berghia'])} | "
                         f"{';'.join(r['placed_berghia'])} | {r['jaccard']:.3f} | "
                         f"{r['reproduced']} | {r['infiltrating_in_with_tree']} |")
        lines.append("")

    lines.append("## Axis 3 — LSE robustness")
    lines.append("")
    lines.append("Berghia-restricted Robinson-Foulds between the with-anchor and "
                 "clean trees. A low rf means adding the anchors did not perturb "
                 "the lineage-specific-expansion structure, so dropping them is "
                 "harmless. A high rf means the anchors did perturb the Berghia "
                 "structure; since the clean tree is the LBA-free reference, a high "
                 "rf supports the drop.")
    lines.append("")
    lines.append("| class | rf | shared_fraction | n_shared_clades | n_clean_clades |")
    lines.append("|---|---|---|---|---|")
    for res in per_class_results:
        a3 = res["axis3"]
        lines.append(f"| {res['class']} | {a3['rf']:.3f} | "
                     f"{a3['shared_fraction']:.3f} | {a3['n_shared_clades']} | "
                     f"{a3['n_clean_clades']} |")
    lines.append("")

    Path(path).write_text("\n".join(lines))


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Place dropped tier-2/3 landmarks on the clean per-class "
                    "trees and run the three directional validation axes.")
    p.add_argument("--calibration-dir", required=True,
                   help="C3 anchor_calibration directory (contains trees_nocloak/, pools/)")
    p.add_argument("--anchor-tsv", required=True, help="anchor_set.tsv (landmark metadata)")
    p.add_argument("--anchor-fasta", required=True, help="anchor_set.fasta (landmark sequences)")
    p.add_argument("--class-berghia-tsv", required=True,
                   help="class_berghia.tsv (axis-1 classifier baseline)")
    p.add_argument("--out", required=True, help="output directory for the 6 report files")
    p.add_argument("--classes", default="A B C F",
                   help="whitespace-separated class list (default: 'A B C F')")
    p.add_argument("--threads", type=int, default=4, help="threads for mafft/epa-ng")
    args = p.parse_args(argv)

    missing = []
    if not os.path.isdir(args.calibration_dir):
        missing.append(f"--calibration-dir directory not found: {args.calibration_dir}")
    for label, pth in (("--anchor-tsv", args.anchor_tsv),
                       ("--anchor-fasta", args.anchor_fasta),
                       ("--class-berghia-tsv", args.class_berghia_tsv)):
        if not os.path.isfile(pth):
            missing.append(f"{label} file not found: {pth}")
    if missing:
        for m in missing:
            print(f"ERROR: {m}", file=sys.stderr)
        sys.exit(2)

    results = []
    for klass in args.classes.split():
        try:
            res = run_class(klass, args.calibration_dir, args.anchor_tsv,
                            args.anchor_fasta, args.class_berghia_tsv, args.out,
                            threads=args.threads)
        except FileNotFoundError as exc:
            print(f"WARNING: skipping class {klass}: {exc}", file=sys.stderr)
            continue
        results.append(res)

    write_report(args.out, results)
    print(f"Wrote landmark-placement report for {len(results)} class(es) to {args.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
