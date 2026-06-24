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
