#!/usr/bin/env python3
"""classify_via_hmm.py — HMM-scan-based classification for the
non-chemoreceptor classification pipeline.

Phase 4 Task 4.1.

For each candidate sequence: hmmscan against the consolidated HMM
library (custom HMMs from Phase 2 + Pfam fallback), apply per-HMM
E-value thresholds from the LOO benchmark (Phase 2 Task 2.2), and
assign the best PASSING hit's family + subfamily. Candidates whose
top hits all fail their thresholds get 'unclassified-hmm'.

Output: TSV with columns
    candidate_id, hmm_family, hmm_subfamily, evalue, score, evidence

Usage:
    python3 classify_via_hmm.py \\
        --candidate-fasta path/to/candidates.fa \\
        --hmm-dir results/classification/hmms/ \\
        --pfam-fallback-dir results/classification/hmms/pfam_fallback/ \\
        --loo-metrics results/classification/loo/loo_metrics.tsv \\
        --output-tsv results/classification/candidate_hmm_assignments.tsv \\
        --threads 4
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Optional


DEFAULT_THRESHOLD = 1e-10  # Conservative default when LOO metrics are absent
                            # or a family has no specific threshold (e.g. Pfam
                            # fallback HMMs). Looser than the LOO-derived
                            # thresholds so candidates can still classify.


# ---- HMMER tblout parsing -----------------------------------------------

def parse_hmmscan_tblout(path: str) -> dict[str, list[tuple[str, float, float]]]:
    """Parse HMMER's --tblout format.

    Returns dict: query_id -> list of (hmm_target, evalue, score)
    sorted ascending by E-value (best hit first).
    """
    if not os.path.exists(path):
        return {}
    by_query: dict[str, list[tuple[str, float, float]]] = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            target = parts[0]
            query = parts[2]
            try:
                evalue = float(parts[4])
                score = float(parts[5])
            except ValueError:
                continue
            by_query.setdefault(query, []).append((target, evalue, score))
    for q in by_query:
        by_query[q].sort(key=lambda t: t[1])
    return by_query


# ---- LOO threshold loading ----------------------------------------------

def load_thresholds(loo_metrics_path: str) -> dict[str, float]:
    """Read per-family E-value thresholds from validate_classification_hmms's
    output (Task 2.2). Returns {} if file missing — caller falls back to
    DEFAULT_THRESHOLD per family."""
    if not os.path.exists(loo_metrics_path):
        return {}
    out: dict[str, float] = {}
    with open(loo_metrics_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            fam = row.get("family", "").strip()
            thr = row.get("evalue_threshold", "").strip()
            if not fam or not thr:
                continue
            try:
                out[fam] = float(thr)
            except ValueError:
                continue
    return out


# ---- HMM name -> family/subfamily (shared helper) ----------------------

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent))
from _classification_labels import label_for_hmm  # noqa: E402, F401


# ---- Classification assignment -------------------------------------------

def assign_classification(hits: list[tuple[str, float, float]],
                          thresholds: dict[str, float]
                          ) -> dict:
    """For a single candidate's hit list (already sorted best E-value first),
    return the best-PASSING-threshold classification, or unclassified-hmm
    if all hits fail."""
    for hmm_target, evalue, score in hits:
        family, subfamily = label_for_hmm(hmm_target)
        threshold = thresholds.get(family, DEFAULT_THRESHOLD)
        if evalue <= threshold:
            return {
                "family": family,
                "subfamily": subfamily,
                "evalue": evalue,
                "score": score,
                "hmm_target": hmm_target,
            }
    return {
        "family": "unclassified-hmm",
        "subfamily": "",
        "evalue": float("inf"),
        "score": 0.0,
        "hmm_target": "",
    }


# ---- HMM library consolidation + scan -----------------------------------

def _load_pfam_manifest(pfam_dir: str) -> dict[str, str]:
    """Read the Pfam fallback manifest to map original HMM NAME (e.g.
    '7tm_1') to our classification label (e.g. 'class-A-7tm'). Returns
    empty dict if no manifest is present.

    Without this mapping, the HMM names '7tm_1' / '7tm_2' / '7tm_3'
    would be parsed by label_for_hmm as ('7tm', '1') etc — splitting
    Pfam's family numbering as if it were a subfamily separator. We
    rename the Pfam HMMs at consolidation time to avoid that.
    """
    manifest = os.path.join(pfam_dir, "manifest.tsv")
    if not os.path.exists(manifest):
        return {}
    out: dict[str, str] = {}
    with open(manifest) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pfam_name = row.get("pfam_name", "").strip()
            label = row.get("classification_label", "").strip()
            if pfam_name and label:
                out[pfam_name] = label
    return out


def _rewrite_hmm_name(hmm_text: str, new_name: str) -> str:
    """Replace the NAME line of an HMM block. HMMER format: line starts
    with 'NAME' followed by whitespace and the name. We only rewrite the
    FIRST NAME line per block (subsequent NAME-prefixed text in alignments
    can't legitimately exist at column 0 in HMMER format)."""
    out_lines: list[str] = []
    renamed = False
    for line in hmm_text.splitlines(keepends=True):
        if not renamed and line.startswith("NAME ") and len(line.split()) >= 2:
            out_lines.append(f"NAME  {new_name}\n")
            renamed = True
        else:
            out_lines.append(line)
    return "".join(out_lines)


def consolidate_hmm_library(hmm_dir: str,
                             pfam_dir: Optional[str],
                             out_path: str) -> None:
    """Concatenate all .hmm files from hmm_dir + pfam_dir into a single
    file and hmmpress it for hmmscan.

    Pfam HMMs are renamed at consolidation time using the Pfam manifest's
    classification_label, so 7tm_1 becomes class-A-7tm etc — avoids the
    label_for_hmm parser splitting Pfam family numbering on the
    underscore separator.
    """
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    pfam_label_map = _load_pfam_manifest(pfam_dir) if pfam_dir else {}

    with open(out_path, "w") as out:
        # Custom HMMs (already have meaningful NAMEs from build step).
        for fname in sorted(os.listdir(hmm_dir)):
            if not fname.endswith(".hmm"):
                continue
            full = os.path.join(hmm_dir, fname)
            if not os.path.isfile(full):
                continue
            with open(full) as f:
                shutil.copyfileobj(f, out)
        # Pfam fallback HMMs — rewrite NAME using the Pfam manifest mapping.
        if pfam_dir and os.path.isdir(pfam_dir):
            for fname in sorted(os.listdir(pfam_dir)):
                if not fname.endswith(".hmm"):
                    continue
                full = os.path.join(pfam_dir, fname)
                with open(full) as f:
                    text = f.read()
                # Find original NAME and remap
                orig_name = ""
                for line in text.splitlines():
                    if line.startswith("NAME "):
                        orig_name = line.split(None, 1)[1].strip()
                        break
                new_name = pfam_label_map.get(orig_name, orig_name)
                if new_name != orig_name:
                    text = _rewrite_hmm_name(text, new_name)
                out.write(text)
    # hmmpress for fast hmmscan; remove existing index files first.
    for ext in ("h3f", "h3i", "h3m", "h3p"):
        p = out_path + "." + ext
        if os.path.exists(p):
            os.remove(p)
    subprocess.run(["hmmpress", out_path],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                   check=True)


def run_hmmscan(query_fasta: str, hmm_db: str, tbl_out: str,
                threads: int = 1) -> None:
    """Run hmmscan; output goes to tbl_out (HMMER --tblout format)."""
    subprocess.run(
        ["hmmscan", "--cpu", str(threads), "--tblout", tbl_out,
         "-E", "1.0", hmm_db, query_fasta],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)


# ---- Main ---------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--candidate-fasta", required=True)
    ap.add_argument("--hmm-dir", required=True,
                    help="Directory with custom family HMMs (from Task 2.1)")
    ap.add_argument("--pfam-fallback-dir", default="",
                    help="Optional directory with Pfam fallback HMMs (Task 2.3)")
    ap.add_argument("--loo-metrics", default="",
                    help="Optional LOO metrics TSV with per-family E-value "
                         "thresholds (from Task 2.2). If absent, the default "
                         f"threshold {DEFAULT_THRESHOLD} is used per family.")
    ap.add_argument("--output-tsv", required=True)
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    thresholds = load_thresholds(args.loo_metrics) if args.loo_metrics else {}
    if thresholds:
        print(f"[hmm-classify] Loaded {len(thresholds)} per-family E-value "
              f"thresholds", file=sys.stderr)

    # Consolidate library + run hmmscan
    work_dir = tempfile.mkdtemp(prefix="hmm_classify_")
    try:
        lib_path = os.path.join(work_dir, "library.hmm")
        consolidate_hmm_library(
            args.hmm_dir,
            args.pfam_fallback_dir if args.pfam_fallback_dir else None,
            lib_path)
        print(f"[hmm-classify] Consolidated library at {lib_path}",
              file=sys.stderr)

        tbl_path = os.path.join(work_dir, "scan.tbl")
        run_hmmscan(args.candidate_fasta, lib_path, tbl_path,
                    threads=args.threads)
        hits = parse_hmmscan_tblout(tbl_path)
        print(f"[hmm-classify] hmmscan: {len(hits)} candidates with hits",
              file=sys.stderr)

        # Read candidate IDs from FASTA so we report unclassified ones too
        candidate_ids: list[str] = []
        with open(args.candidate_fasta) as f:
            for line in f:
                if line.startswith(">"):
                    candidate_ids.append(line[1:].split()[0])

        # Write output TSV
        Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output_tsv, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["candidate_id", "hmm_family", "hmm_subfamily",
                        "evalue", "score", "hmm_target"])
            n_classified = 0
            for cid in candidate_ids:
                cand_hits = hits.get(cid, [])
                result = assign_classification(cand_hits, thresholds)
                if result["family"] != "unclassified-hmm":
                    n_classified += 1
                evalue_str = (f"{result['evalue']:.3e}"
                              if result["evalue"] != float("inf") else "inf")
                w.writerow([
                    cid, result["family"], result["subfamily"],
                    evalue_str, f"{result['score']:.2f}",
                    result["hmm_target"],
                ])

        print(f"[hmm-classify] DONE: {n_classified}/{len(candidate_ids)} "
              f"candidates classified -> {args.output_tsv}",
              file=sys.stderr)
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
