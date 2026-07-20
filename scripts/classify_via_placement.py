#!/usr/bin/env python3
"""classify_via_placement.py — phylogenetic placement classifier (EPA-ng).

Phase 4 Task 4.3 of the non-chemoreceptor classification pipeline.

Two-stage placement:
    1. Backbone placement: align candidates to the backbone reference
       alignment (mafft --add) and run EPA-ng against the backbone tree.
       Each candidate gets an edge + LWR. Edge -> family via the
       leaf-annotation TSV (Phase 3 backbone.tsv).
    2. Subtree placement: candidates classified at the backbone as
       'aminergic' or 'peptide' (the two coarse families with medium-
       drill-down trees) are re-aligned and placed on the corresponding
       subtree to refine the subfamily call.

Confidence gate: best placement must satisfy LWR >= 0.80
(default; tunable). Below that, the call is 'unclassified-placement'.

Output TSV: candidate_id, placement_family, placement_subfamily,
placement_lwr, placement_edge.

Usage:
    python3 classify_via_placement.py \\
        --candidate-fasta path/to/candidates.fa \\
        --tree-dir results/classification/trees/ \\
        --output-tsv results/classification/candidate_placement.tsv \\
        --threads 4
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Optional

LWR_THRESHOLD = 0.80


# ---- jplace parsing ----------------------------------------------------

def parse_jplace(path: str) -> dict[str, dict]:
    """Read EPA-ng jplace output. Returns dict: query_id -> dict with
    keys edge_num + lwr (best placement; max lwr)."""
    if not os.path.exists(path):
        return {}
    try:
        with open(path) as f:
            data = json.load(f)
    except (json.JSONDecodeError, OSError):
        return {}

    fields = data.get("fields", [])
    try:
        edge_idx = fields.index("edge_num")
        lwr_idx = fields.index("like_weight_ratio")
    except ValueError:
        return {}

    out: dict[str, dict] = {}
    for entry in data.get("placements", []):
        # Query identification can be in `nm` (named multiplicity) or `n`
        names: list[str] = []
        if "nm" in entry:
            names = [nm[0] for nm in entry["nm"]]
        elif "n" in entry:
            names = list(entry["n"])
        if not names:
            continue
        # Best placement = max LWR
        best = max(entry.get("p", []),
                    key=lambda p: p[lwr_idx], default=None)
        if best is None:
            continue
        for name in names:
            out[name] = {
                "edge_num": int(best[edge_idx]),
                "lwr": float(best[lwr_idx]),
            }
    return out


# ---- edge-to-family mapping --------------------------------------------

# Match an edge label `{N}` after each branch in EPA-ng-style newick.
# We use a simple recursive-descent parser to walk the tree and assign
# every edge_num to its descendant leaf set.
_EDGE_RE = re.compile(r"\{(\d+)\}")


def _strip_edge_labels(newick: str) -> tuple[str, dict[int, str]]:
    """Remove {N} edge labels from newick; return (clean_newick, edge_id_map).
    edge_id_map: edge_num -> the substring it labelled (kept for debugging)."""
    edge_map: dict[int, str] = {}
    for m in _EDGE_RE.finditer(newick):
        edge_map[int(m.group(1))] = newick[max(0, m.start() - 30):m.end()]
    clean = _EDGE_RE.sub("", newick)
    return clean, edge_map


def _parse_newick_with_edge_labels(newick: str) -> dict[int, list[str]]:
    """Parse a newick tree with EPA-ng edge labels, return dict mapping
    edge_num -> list of leaf names in that edge's descendant subtree.

    Lightweight parser (avoids ete3 dependency for this small tree).
    Handles: leaf names (alphanumeric + - _ . / : are name chars), edge
    labels {N}, branch lengths :X.Y, parens, commas, semicolon."""

    # Tokenize, preserving structure characters separately
    tokens: list[str] = []
    i = 0
    n = len(newick)
    while i < n:
        c = newick[i]
        if c in "(),;":
            tokens.append(c)
            i += 1
        elif c == ":":
            # Branch length: consume until next struct char or {
            j = i + 1
            while j < n and newick[j] not in "(),;{":
                j += 1
            tokens.append(("BRANCH", newick[i + 1:j]))
            i = j
        elif c == "{":
            j = newick.index("}", i)
            tokens.append(("EDGE", int(newick[i + 1:j])))
            i = j + 1
        elif c.isspace():
            i += 1
        else:
            # Leaf or internal node name
            j = i
            while j < n and newick[j] not in "(),;:{":
                j += 1
            tokens.append(("NAME", newick[i:j]))
            i = j

    # Recursive descent: build a tree of (children_or_name, edge_num)
    # and walk to map edge -> leaf set.
    pos = [0]
    edge_to_leaves: dict[int, list[str]] = {}

    def parse_clade() -> tuple[list[str], int | None]:
        """Returns (leaves_in_subtree, edge_num_above) for this clade."""
        tok = tokens[pos[0]]
        leaves: list[str] = []
        if tok == "(":
            pos[0] += 1  # consume '('
            while True:
                child_leaves, child_edge = parse_clade()
                leaves.extend(child_leaves)
                if child_edge is not None:
                    edge_to_leaves[child_edge] = list(child_leaves)
                if pos[0] < len(tokens) and tokens[pos[0]] == ",":
                    pos[0] += 1
                    continue
                break
            assert tokens[pos[0]] == ")"
            pos[0] += 1  # consume ')'
            # Optional internal node name
            if pos[0] < len(tokens) and isinstance(tokens[pos[0]], tuple) \
                    and tokens[pos[0]][0] == "NAME":
                pos[0] += 1
        elif isinstance(tok, tuple) and tok[0] == "NAME":
            leaves.append(tok[1])
            pos[0] += 1
        else:
            return leaves, None
        # Optional :branch_length and {edge_num}
        edge_num = None
        if pos[0] < len(tokens) and isinstance(tokens[pos[0]], tuple) \
                and tokens[pos[0]][0] == "BRANCH":
            pos[0] += 1
        if pos[0] < len(tokens) and isinstance(tokens[pos[0]], tuple) \
                and tokens[pos[0]][0] == "EDGE":
            edge_num = tokens[pos[0]][1]
            pos[0] += 1
        return leaves, edge_num

    leaves, root_edge = parse_clade()
    if root_edge is not None:
        edge_to_leaves[root_edge] = list(leaves)

    return edge_to_leaves


def build_edge_to_family_map(tree_newick: str,
                              leaf_annotations: dict[str, tuple[str, str]]
                              ) -> dict[int, tuple[str, str]]:
    """For each edge in the EPA-ng-labelled newick, determine the
    (family, subfamily) it represents:
      - leaf edge: that leaf's family/subfamily
      - internal monophyletic edge: clade's shared family; clade's
        shared subfamily if all members agree
      - mixed clade: ('', '')
    """
    edge_to_leaves = _parse_newick_with_edge_labels(tree_newick)
    edge_to_family: dict[int, tuple[str, str]] = {}
    for edge_num, leaves in edge_to_leaves.items():
        fams = set()
        subs = set()
        for leaf in leaves:
            # Tree leaves are written as `accession|family|subfamily|species`
            # by select_backbone_reps.py; the leaf-annotation TSV is keyed
            # by accession alone. Strip pipe-suffix before lookup, falling
            # back to the full name for backward-compat with simpler trees.
            leaf_acc = leaf.split("|", 1)[0]
            ann = leaf_annotations.get(leaf_acc) or leaf_annotations.get(leaf)
            if ann is None:
                fams.add("?unknown-leaf?")
            else:
                fams.add(ann[0])
                if ann[1]:
                    subs.add(ann[1])
        if len(fams) == 1 and "?unknown-leaf?" not in fams:
            family = next(iter(fams))
            subfamily = next(iter(subs)) if len(subs) == 1 else ""
            edge_to_family[edge_num] = (family, subfamily)
        else:
            edge_to_family[edge_num] = ("", "")
    return edge_to_family


# ---- LWR-gated classification ------------------------------------------

def _finite_lwr(value) -> float | None:
    """Coerce a likelihood weight ratio to a finite float, else None.

    IEEE-754 makes every comparison with NaN False, so a bare
    `lwr < threshold` gate lets a NaN fall through to the CONFIDENT branch —
    a placement whose LWR could not be computed at all was being emitted as
    a real family call. A missing or degenerate LWR is the least confident
    outcome possible, so it must land in the same bucket as a low one.
    """
    try:
        lwr = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(lwr):
        return None
    return lwr

def classify_placement_for_id(
    candidate_id: str,
    placements: dict[str, dict],
    edge_map: dict[int, tuple[str, str]],
    lwr_threshold: float = LWR_THRESHOLD,
) -> dict:
    if candidate_id not in placements:
        return {
            "family": "unclassified-placement",
            "subfamily": "",
            "lwr": 0.0,
            "edge": -1,
        }
    p = placements[candidate_id]
    lwr = _finite_lwr(p.get("lwr"))
    if lwr is None or lwr < lwr_threshold:
        return {
            "family": "unclassified-placement",
            "subfamily": "",
            # A degenerate LWR is reported as 0.0, never propagated: a NaN
            # would compare False against every downstream threshold too,
            # recreating this bug further along the pipeline.
            "lwr": 0.0 if lwr is None else lwr,
            "edge": p.get("edge_num", -1),
        }
    fam, sub = edge_map.get(p["edge_num"], ("", ""))
    if not fam:
        return {
            "family": "unclassified-placement",
            "subfamily": "",
            "lwr": p["lwr"],
            "edge": p["edge_num"],
        }
    return {"family": fam, "subfamily": sub,
            "lwr": p["lwr"], "edge": p["edge_num"]}


def classify_placement(placements: dict[str, dict],
                       edge_map: dict[int, tuple[str, str]],
                       lwr_threshold: float = LWR_THRESHOLD
                       ) -> dict[str, dict]:
    """Apply classify_placement_for_id over all queries in `placements`."""
    return {q: classify_placement_for_id(q, placements, edge_map, lwr_threshold)
            for q in placements}


# ---- leaf-annotation TSV loader -----------------------------------------

def load_leaf_annotations(tsv_path: str) -> dict[str, tuple[str, str]]:
    """Load (accession -> (family, subfamily)) from a Phase 3 tree-companion
    TSV (e.g. backbone.tsv produced by select_backbone_reps.py)."""
    out: dict[str, tuple[str, str]] = {}
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            acc = row.get("accession", "").strip()
            fam = row.get("family", "").strip()
            sub = row.get("subfamily", "").strip()
            if acc and fam:
                out[acc] = (fam, sub)
    return out


# ---- EPA-ng + mafft --add invocation -----------------------------------

def _run_mafft_add(query_fasta: str, ref_alignment: str,
                   out_path: str, threads: int) -> bool:
    """Add candidate sequences to an existing alignment via mafft --add.
    Returns True on success."""
    try:
        with open(out_path, "w") as out:
            subprocess.run(
                ["mafft", "--add", query_fasta, "--keeplength",
                 "--quiet", "--thread", str(threads), ref_alignment],
                stdout=out, stderr=subprocess.DEVNULL, check=True)
        return os.path.getsize(out_path) > 0
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def _split_query_ref(combined_aln: str, ref_ids: set[str],
                     query_path: str, ref_path: str) -> None:
    """Split a combined alignment into ref-only and query-only files."""
    cur_id: Optional[str] = None
    cur_lines: list[str] = []
    cur_is_ref: bool = False
    with open(combined_aln) as f, \
            open(query_path, "w") as q, \
            open(ref_path, "w") as r:
        def _flush() -> None:
            if cur_id is None:
                return
            target = r if cur_is_ref else q
            target.write(f">{cur_id}\n")
            target.write("".join(cur_lines))
        for line in f:
            line_stripped = line.rstrip("\n")
            if line_stripped.startswith(">"):
                _flush()
                cur_id = line_stripped[1:].split()[0]
                cur_lines = []
                cur_is_ref = cur_id in ref_ids
            else:
                cur_lines.append(line)
        _flush()


def _read_iqtree_model(ref_tree: str) -> str:
    """Read the IQ-TREE-inferred model from the sibling .iqtree file.
    EPA-ng defaults to GTR+G (DNA) which crashes on AA data; we must
    pass the protein model explicitly. Falls back to 'LG+G' if the
    .iqtree file is missing or unparsable (a safe default for GPCRs)."""
    iqtree_path = ref_tree
    for suffix in (".treefile", ".contree"):
        if iqtree_path.endswith(suffix):
            iqtree_path = iqtree_path[: -len(suffix)] + ".iqtree"
            break
    try:
        with open(iqtree_path) as f:
            for line in f:
                if line.startswith("Model of substitution:"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return "LG+G"


def _run_epa_ng(query_aln: str, ref_aln: str, ref_tree: str,
                out_dir: str, threads: int) -> bool:
    """Run EPA-ng against the reference tree. Output: out_dir/epa_result.jplace.

    Captures stderr so the failure is visible if EPA-ng crashes (the
    previous DEVNULL-everything approach hid a GTR+G default that
    silently broke placement for protein data — see 2026-05-07
    smoke-test postmortem).
    """
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    model = _read_iqtree_model(ref_tree)
    try:
        result = subprocess.run(
            ["epa-ng",
             "--redo",
             "--tree", ref_tree,
             "--ref-msa", ref_aln,
             "--query", query_aln,
             "--out-dir", out_dir,
             "--threads", str(threads),
             "--model", model],
            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, check=True)
        return os.path.exists(os.path.join(out_dir, "epa_result.jplace"))
    except subprocess.CalledProcessError as e:
        err = (e.stderr or b"").decode("utf-8", errors="replace")
        print(f"  EPA-ng stderr: {err.strip()[:500]}", file=sys.stderr)
        return False
    except FileNotFoundError:
        print(f"  EPA-ng not on PATH", file=sys.stderr)
        return False


# ---- Main orchestration -------------------------------------------------

def place_against_tree(query_fasta: str, ref_fasta: str, ref_aln: str,
                       ref_tree: str, ref_tsv: str, work_dir: str,
                       threads: int = 1
                       ) -> dict[str, dict]:
    """Run mafft --add + EPA-ng for one tree; return placement results
    per query (after LWR gating + edge->family mapping)."""
    Path(work_dir).mkdir(parents=True, exist_ok=True)
    combined = os.path.join(work_dir, "combined.aln")
    query_aln = os.path.join(work_dir, "query.aln")
    splitref_aln = os.path.join(work_dir, "splitref.aln")
    epa_dir = os.path.join(work_dir, "epa")

    if not _run_mafft_add(query_fasta, ref_aln, combined, threads):
        print(f"  WARN: mafft --add failed for {ref_aln}", file=sys.stderr)
        return {}

    # Split combined alignment back into query/ref portions
    ref_ids: set[str] = set()
    with open(ref_aln) as f:
        for line in f:
            if line.startswith(">"):
                ref_ids.add(line[1:].split()[0])
    _split_query_ref(combined, ref_ids, query_aln, splitref_aln)

    if not _run_epa_ng(query_aln, splitref_aln, ref_tree, epa_dir, threads):
        print(f"  WARN: EPA-ng failed for {ref_tree}", file=sys.stderr)
        return {}

    jplace_path = os.path.join(epa_dir, "epa_result.jplace")
    placements = parse_jplace(jplace_path)
    leaf_ann = load_leaf_annotations(ref_tsv)
    # Read the EPA-ng-labeled tree from inside the jplace file (it has the
    # `{N}` edge labels that placements reference). The IQ-TREE treefile on
    # disk has SH-aLRT/UFBoot labels but no `{N}` labels — using it gives
    # an empty edge_map and 'unclassified-placement' for every candidate.
    tree_text = ""
    try:
        with open(jplace_path) as f:
            jdata = json.load(f)
        tree_text = (jdata.get("tree") or "").strip()
    except (OSError, json.JSONDecodeError):
        pass
    if not tree_text:
        with open(ref_tree) as f:
            tree_text = f.read().strip()
    edge_map = build_edge_to_family_map(tree_text, leaf_ann)
    return classify_placement(placements, edge_map)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--candidate-fasta", required=True)
    ap.add_argument("--tree-dir", required=True,
                    help="Directory with backbone + subtree files "
                         "(from build_classification_reference_trees.sh)")
    ap.add_argument("--output-tsv", required=True)
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    tree_dir = Path(args.tree_dir)
    # Prefer the post-defensive-filter alignment that the tree was actually
    # built from. build_classification_reference_trees.sh drops ≥95%-gap rows
    # into <name>.cleaned.aln and feeds THAT to IQ-TREE, so the tree's leaf
    # set matches cleaned.aln (not trimmed.aln). Using trimmed.aln gives a
    # silent EPA-ng MSA-tree mismatch and 0 placements.
    backbone_aln = tree_dir / "backbone.cleaned.aln"
    if not backbone_aln.exists():
        backbone_aln = tree_dir / "backbone.trimmed.aln"
    backbone_tre = tree_dir / "backbone.treefile"
    backbone_tsv = tree_dir / "backbone.tsv"

    if not (backbone_aln.exists() and backbone_tre.exists()
            and backbone_tsv.exists()):
        print(f"ERROR: backbone tree artifacts missing in {tree_dir}",
              file=sys.stderr)
        return 1

    work_root = tempfile.mkdtemp(prefix="classify_placement_")
    try:
        # --- Stage 1: backbone placement ---
        print("[placement] Stage 1: backbone placement...", file=sys.stderr)
        backbone_results = place_against_tree(
            args.candidate_fasta, str(backbone_aln), str(backbone_aln),
            str(backbone_tre), str(backbone_tsv),
            os.path.join(work_root, "backbone"),
            threads=args.threads)
        n_classified = sum(1 for r in backbone_results.values()
                            if r["family"] != "unclassified-placement")
        print(f"  -> {n_classified}/{len(backbone_results)} candidates "
              f"placed on backbone with LWR >= {LWR_THRESHOLD}",
              file=sys.stderr)

        # --- Stage 2: subtree placement (aminergic / peptide only) ---
        final_results: dict[str, dict] = dict(backbone_results)
        for fam, name in [("aminergic", "aminergic_subtree"),
                          ("peptide", "peptide_subtree")]:
            sub_aln = tree_dir / f"{name}.cleaned.aln"
            if not sub_aln.exists():
                sub_aln = tree_dir / f"{name}.trimmed.aln"
            sub_tre = tree_dir / f"{name}.treefile"
            sub_tsv = tree_dir / f"{name}.tsv"
            if not (sub_aln.exists() and sub_tre.exists() and sub_tsv.exists()):
                print(f"  WARN: {name} artifacts missing, skipping subfamily "
                      f"refinement", file=sys.stderr)
                continue
            # Candidates that placed at this family on backbone -> drill down
            cand_ids = [cid for cid, r in backbone_results.items()
                         if r["family"] == fam]
            if not cand_ids:
                continue
            # Write a subset query FASTA
            sub_query = os.path.join(work_root, f"{name}_query.fa")
            keep = set(cand_ids)
            with open(args.candidate_fasta) as src, \
                    open(sub_query, "w") as dst:
                keep_block = False
                for line in src:
                    if line.startswith(">"):
                        cid = line[1:].split()[0]
                        keep_block = cid in keep
                    if keep_block:
                        dst.write(line)
            print(f"[placement] Stage 2: {name} ({len(cand_ids)} candidates)",
                  file=sys.stderr)
            sub_results = place_against_tree(
                sub_query, str(sub_aln), str(sub_aln), str(sub_tre),
                str(sub_tsv), os.path.join(work_root, name),
                threads=args.threads)
            for cid, sr in sub_results.items():
                if sr["family"] != "unclassified-placement":
                    # Override with finer-grained subfamily call
                    final_results[cid]["subfamily"] = sr["subfamily"]
                    final_results[cid]["edge"] = sr["edge"]
                    final_results[cid]["lwr"] = sr["lwr"]

        # --- Output TSV ---
        Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
        candidate_ids: list[str] = []
        with open(args.candidate_fasta) as f:
            for line in f:
                if line.startswith(">"):
                    candidate_ids.append(line[1:].split()[0])
        with open(args.output_tsv, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["candidate_id", "placement_family",
                        "placement_subfamily", "placement_lwr",
                        "placement_edge"])
            for cid in candidate_ids:
                r = final_results.get(cid, {
                    "family": "unclassified-placement",
                    "subfamily": "", "lwr": 0.0, "edge": -1,
                })
                w.writerow([cid, r["family"], r["subfamily"],
                            f"{r['lwr']:.3f}", r["edge"]])
        print(f"[placement] DONE -> {args.output_tsv}", file=sys.stderr)
    finally:
        shutil.rmtree(work_root, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
