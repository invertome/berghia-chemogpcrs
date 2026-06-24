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
