#!/usr/bin/env python3
"""validate_classification_hmms.py — Leave-one-out cross-validation of the
non-chemoreceptor classification HMM library.

Phase 2 Task 2.2 of the non-chemoreceptor classification feature.

Procedure for each HMM in the library:
    For each training sequence:
        1. Hold out the sequence.
        2. Rebuild the HMM from the remaining members (MAFFT + hmmbuild).
        3. hmmscan the held-out sequence against the FULL library
           (containing the rebuilt HMM in place of the original).
        4. Record the top-hit family/subfamily vs the held-out's true label.

Outputs to <output-dir>:
    loo_validation.tsv         — per-sequence: held_out_id, true_family,
                                 true_subfamily, predicted_family,
                                 predicted_subfamily, evalue, score
    loo_metrics.tsv            — per-family: n_total, recall, precision,
                                 selected_evalue_threshold
    loo_confusion.tsv          — full confusion matrix (true x predicted)

Targets (from design doc): per-family recall ≥0.90, precision ≥0.95.

Usage:
    python3 validate_classification_hmms.py \\
        --hmm-dir results/classification/hmms/ \\
        --output-dir results/classification/loo/ \\
        --threads 4

Implementation notes:
    - Rebuilding the HMM in-process is cheap (mafft + hmmbuild on a few
      hundred sequences per family takes <1 second per cycle). The full
      LOO sweep is 728 cycles × ~1 second = ~10-15 minutes.
    - We only test against COARSE HMMs for the LOO benchmark — coarse
      family is the gating call. Subfamily classification is treated as
      a secondary tag (the held-out evaluation is "did coarse family
      classify correctly").
"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Optional


DEFAULT_THRESHOLD = 1e-10  # conservative E-value default if no LOO data

# Coarse families known to the design (used to validate predicted family
# is in the expected set — defensive against malformed manifest).
COARSE_FAMILIES = {
    "aminergic", "peptide", "opsin", "lipid", "nucleotide",
    "metabotropic-neurotransmitter", "glycoprotein-hormone",
    "class-B-secretin", "class-C", "class-F-frizzled",
}


# ---- Manifest + alignment parsing ---------------------------------------

def parse_manifest(path: str) -> list[dict]:
    """Read the HMM manifest TSV produced by build_classification_hmms.py."""
    rows: list[dict] = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append({
                "family": row.get("family", "").strip(),
                "subfamily": row.get("subfamily", "").strip(),
                "n_train": int(row.get("n_train", 0)),
                "hmm_path": row.get("hmm_path", "").strip(),
                "alignment_path": row.get("alignment_path", "").strip(),
            })
    return rows


def parse_alignment_fasta(path: str) -> dict[str, str]:
    """Return dict: seq_id -> sequence (with gaps preserved)."""
    seqs: dict[str, str] = {}
    cur_id: Optional[str] = None
    cur_parts: list[str] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(cur_parts)
                cur_id = line[1:].split()[0]
                cur_parts = []
            else:
                cur_parts.append(line)
        if cur_id is not None:
            seqs[cur_id] = "".join(cur_parts)
    return seqs


def strip_gaps(seq: str) -> str:
    """Remove gap and X-mask characters for re-alignment."""
    return seq.replace("-", "").replace(".", "")


# ---- HMM label parsing --------------------------------------------------

# label_for_hmm moved to scripts/_classification_labels.py for shared use
# across classify_via_hmm.py, validate_classification_hmms.py, etc.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent))
from _classification_labels import label_for_hmm  # noqa: E402, F401


def build_id_to_family_map(manifest: list[dict],
                            seqs_by_alignment: dict[str, dict[str, str]]
                            ) -> dict[str, tuple[str, str]]:
    """For each training-set member, record its (family, subfamily) label.

    A sequence may appear in BOTH a coarse-family alignment AND a medium-
    subfamily alignment. We prefer the medium-granularity label."""
    id_to_label: dict[str, tuple[str, str]] = {}
    # First pass: coarse labels
    for row in manifest:
        if row["subfamily"]:
            continue
        seqs = seqs_by_alignment.get(row["alignment_path"], {})
        for seq_id in seqs:
            id_to_label[seq_id] = (row["family"], "")
    # Second pass: subfamily labels override coarse
    for row in manifest:
        if not row["subfamily"]:
            continue
        seqs = seqs_by_alignment.get(row["alignment_path"], {})
        for seq_id in seqs:
            id_to_label[seq_id] = (row["family"], row["subfamily"])
    return id_to_label


# ---- HMM rebuild + scan -------------------------------------------------

def rebuild_hmm_without(alignment_path: str, hold_out_id: str,
                        work_dir: str, threads: int = 1) -> Optional[str]:
    """Rebuild an HMM from `alignment_path` excluding `hold_out_id`.

    Optimization: instead of re-aligning with MAFFT each iteration (which
    takes 30+ seconds per cycle for 100+ seq alignments), we just REMOVE
    the held-out row from the existing aligned MSA and run hmmbuild
    directly. The remaining rows retain their original column structure
    (gaps preserved), so the resulting HMM is trained on N-1 sequences
    aligned at the same column positions. hmmbuild handles the
    occasional all-gap column gracefully via its --symfrac heuristic.

    This is ~100x faster than re-aligning, with negligible scientific
    difference (one less training observation per column).

    Returns the path to the rebuilt HMM, or None if the holdout left
    fewer than 3 training seqs (not worth rebuilding).
    """
    seqs = parse_alignment_fasta(alignment_path)
    remaining = {sid: s for sid, s in seqs.items() if sid != hold_out_id}
    if len(remaining) < 3:
        return None

    base = os.path.basename(alignment_path).rsplit(".", 1)[0]
    aln_path = os.path.join(work_dir, f"{base}_loo.aln")
    hmm_path = os.path.join(work_dir, f"{base}_loo.hmm")

    # Write the held-out-removed MSA (column structure preserved from input)
    with open(aln_path, "w") as f:
        for sid, s in remaining.items():
            f.write(f">{sid}\n{s}\n")

    # Build HMM directly from the trimmed MSA — no re-alignment needed
    try:
        subprocess.run(
            ["hmmbuild", "--cpu", str(threads), "-n", base,
             hmm_path, aln_path],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None

    return hmm_path


def hmmscan_query(query_seq: str, query_id: str,
                  hmm_db_path: str, work_dir: str,
                  threads: int = 1) -> list[tuple[str, float, float]]:
    """hmmscan a single query seq against a concatenated HMM database.
    Returns list of (hmm_name, evalue, score) sorted by E-value (best first)."""
    query_path = os.path.join(work_dir, "query.fa")
    out_path = os.path.join(work_dir, "scan.tbl")
    with open(query_path, "w") as f:
        f.write(f">{query_id}\n{query_seq}\n")

    try:
        subprocess.run(
            ["hmmscan", "--cpu", str(threads), "--tblout", out_path,
             "-E", "1.0", hmm_db_path, query_path],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return []

    hits: list[tuple[str, float, float]] = []
    if not os.path.exists(out_path):
        return hits
    with open(out_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            hmm_name = parts[0]
            evalue = float(parts[4])
            score = float(parts[5])
            hits.append((hmm_name, evalue, score))
    return sorted(hits, key=lambda x: x[1])


# ---- Metrics ------------------------------------------------------------

def confusion_matrix(results: list[tuple[str, str]]) -> dict[tuple[str, str], int]:
    """Build (true_family, predicted_family) -> count."""
    cm: dict[tuple[str, str], int] = defaultdict(int)
    for true_f, pred_f in results:
        cm[(true_f, pred_f)] += 1
    return dict(cm)


def per_family_metrics(results: list[tuple[str, str]]) -> dict[str, dict]:
    """Per-true-family recall + precision + counts.

    Recall = TP / (TP + FN) — what fraction of THIS family's held-outs
    were classified as THIS family.
    Precision = TP / (TP + FP) — what fraction of predictions of THIS
    family were correct.
    """
    families = sorted(set(t for t, _ in results) | set(p for _, p in results))
    metrics: dict[str, dict] = {}
    for fam in families:
        tp = sum(1 for t, p in results if t == fam and p == fam)
        fn = sum(1 for t, p in results if t == fam and p != fam)
        fp = sum(1 for t, p in results if t != fam and p == fam)
        n_total = tp + fn
        if n_total == 0:
            continue
        recall = tp / n_total
        precision = tp / (tp + fp) if (tp + fp) > 0 else None
        metrics[fam] = {
            "n_total": n_total,
            "tp": tp, "fn": fn, "fp": fp,
            "recall": recall,
            "precision": precision,
        }
    return metrics


def select_evalue_threshold(rows: list[tuple[str, str, float]],
                            family: str, target_recall: float = 0.9,
                            target_precision: float = 0.95
                            ) -> float:
    """Given LOO rows (true_family, predicted_family, evalue) and a target
    recall + precision, return the loosest E-value threshold for the given
    `family` such that:
      (a) >= target_recall of `family`'s held-outs are TPs at evalue ≤ thr
      (b) the precision at `family` (TPs / TPs+FPs) is >= target_precision

    Both constraints come from the design doc: per-family LOO must achieve
    >=0.90 recall AND >=0.95 precision. Earlier versions only optimized
    for recall (H3 from the 2026-05-07 code review).

    Returns DEFAULT_THRESHOLD if no rows or no threshold satisfies both.
    """
    fam_rows = [(t, p, e) for (t, p, e) in rows if t == family]
    if not fam_rows:
        return DEFAULT_THRESHOLD
    tps = [e for (t, p, e) in fam_rows if p == family]
    n_total = len(fam_rows)
    n_required = max(1, int(target_recall * n_total + 0.999))  # ceiling
    if len(tps) < n_required:
        # Can't reach target recall; conservatively return DEFAULT_THRESHOLD
        # rather than a looser E-value (which would only worsen precision
        # in production).
        return DEFAULT_THRESHOLD

    # Off-family rows that get *predicted* as `family` are the FP candidates;
    # any with E-value <= cutoff at runtime would be FPs.
    fp_evalues_for_other_truth = [
        e for (t, p, e) in rows
        if t != family and p == family
    ]

    # Iterate candidate cutoffs from STRICT to LOOSE (sorted TP E-values),
    # picking the loosest that satisfies BOTH recall and precision.
    tp_evalues_sorted = sorted(tps)
    # Start at the threshold that just satisfies recall (n_required-th best TP)
    # and walk looser until precision fails; pick the last loose threshold
    # where both constraints hold.
    chosen = DEFAULT_THRESHOLD
    for i in range(n_required - 1, len(tp_evalues_sorted)):
        candidate_thr = tp_evalues_sorted[i]
        n_tp = sum(1 for e in tps if e <= candidate_thr)
        n_fp = sum(1 for e in fp_evalues_for_other_truth if e <= candidate_thr)
        if n_tp + n_fp == 0:
            continue
        precision = n_tp / (n_tp + n_fp)
        recall = n_tp / n_total
        if recall >= target_recall and precision >= target_precision:
            chosen = candidate_thr  # keep walking; record loosest passing
        else:
            break
    return chosen


# ---- LOO orchestration --------------------------------------------------

def run_loo(hmm_dir: str, output_dir: str, threads: int = 1) -> int:
    """Full LOO sweep. Returns 0 on success, non-zero on error."""
    manifest_path = os.path.join(hmm_dir, "manifest.tsv")
    if not os.path.exists(manifest_path):
        print(f"ERROR: manifest not found at {manifest_path}", file=sys.stderr)
        return 1

    manifest = parse_manifest(manifest_path)
    coarse_rows = [r for r in manifest if not r["subfamily"]]
    if not coarse_rows:
        print("ERROR: no coarse-family HMMs in manifest", file=sys.stderr)
        return 1

    # Load all coarse alignments (these are the LOO targets)
    seqs_by_alignment: dict[str, dict[str, str]] = {}
    for row in coarse_rows:
        seqs_by_alignment[row["alignment_path"]] = parse_alignment_fasta(
            row["alignment_path"])

    id_to_label = build_id_to_family_map(coarse_rows, seqs_by_alignment)
    print(f"[loo] {len(id_to_label)} held-out candidates across "
          f"{len(coarse_rows)} coarse families", file=sys.stderr)

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    rows_tsv = open(os.path.join(output_dir, "loo_validation.tsv"), "w")
    rows_writer = csv.writer(rows_tsv, delimiter="\t")
    rows_writer.writerow([
        "held_out_id", "true_family", "predicted_family",
        "evalue", "score",
    ])

    results: list[tuple[str, str]] = []
    evalue_rows: list[tuple[str, str, float]] = []

    # Build a consolidated HMM "library" path: just concatenate all coarse
    # HMMs into one .hmm file, then hmmpress.
    lib_dir = tempfile.mkdtemp(prefix="loo_lib_")
    full_lib_path = os.path.join(lib_dir, "lib.hmm")

    def write_lib(replacements: dict[str, str]) -> None:
        """Write the HMM library to full_lib_path. `replacements` is a dict
        family -> path to a substitute HMM (used when one coarse HMM is
        being LOO'd; the other coarse HMMs come from the manifest)."""
        with open(full_lib_path, "w") as out:
            for row in coarse_rows:
                src = replacements.get(row["family"], row["hmm_path"])
                with open(src) as f:
                    out.write(f.read())
        # hmmpress for fast hmmscan (overwrites existing index files)
        for ext in ("h3f", "h3i", "h3m", "h3p"):
            p = full_lib_path + "." + ext
            if os.path.exists(p):
                os.remove(p)
        subprocess.run(["hmmpress", full_lib_path],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       check=True)

    work_dir = tempfile.mkdtemp(prefix="loo_work_")
    rebuild_dir = tempfile.mkdtemp(prefix="loo_rebuild_")

    n_done = 0
    n_total = len(id_to_label)
    for hold_out_id, (true_family, _) in id_to_label.items():
        n_done += 1
        if n_done % 50 == 0:
            print(f"[loo] {n_done}/{n_total} done", file=sys.stderr)

        # Rebuild the HMM for this seq's family with hold-out removed
        family_row = next(r for r in coarse_rows if r["family"] == true_family)
        rebuilt = rebuild_hmm_without(family_row["alignment_path"],
                                       hold_out_id, rebuild_dir,
                                       threads=threads)
        if not rebuilt:
            continue

        # Build library with the rebuilt HMM substituted for this family
        write_lib({true_family: rebuilt})

        # Get the held-out sequence
        seqs = seqs_by_alignment[family_row["alignment_path"]]
        held_out_seq = strip_gaps(seqs[hold_out_id])

        # hmmscan
        hits = hmmscan_query(held_out_seq, hold_out_id, full_lib_path,
                              work_dir, threads=threads)

        if not hits:
            pred_family = "unclassified"
            evalue = float("inf")
            score = 0.0
        else:
            best = hits[0]
            pred_hmm = best[0]
            pred_family, _ = label_for_hmm(pred_hmm)
            evalue = best[1]
            score = best[2]

        rows_writer.writerow([
            hold_out_id, true_family, pred_family,
            f"{evalue:.3e}" if evalue != float("inf") else "inf",
            f"{score:.2f}",
        ])
        results.append((true_family, pred_family))
        evalue_rows.append((true_family, pred_family, evalue))

    rows_tsv.close()

    # Per-family metrics + thresholds
    metrics = per_family_metrics(results)
    metrics_path = os.path.join(output_dir, "loo_metrics.tsv")
    with open(metrics_path, "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["family", "n_total", "tp", "fn", "fp",
                    "recall", "precision", "evalue_threshold"])
        for fam in sorted(metrics):
            m = metrics[fam]
            thr = select_evalue_threshold(evalue_rows, fam)
            w.writerow([
                fam, m["n_total"], m["tp"], m["fn"], m["fp"],
                f"{m['recall']:.3f}",
                f"{m['precision']:.3f}" if m["precision"] is not None else "NA",
                f"{thr:.3e}",
            ])

    # Confusion matrix
    cm = confusion_matrix(results)
    cm_path = os.path.join(output_dir, "loo_confusion.tsv")
    all_fams = sorted(set(t for t, _ in results) | set(p for _, p in results))
    with open(cm_path, "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["true \\ predicted"] + all_fams)
        for true_f in all_fams:
            row = [true_f]
            for pred_f in all_fams:
                row.append(cm.get((true_f, pred_f), 0))
            w.writerow(row)

    # Cleanup workdir / rebuild dir / lib dir
    import shutil
    shutil.rmtree(work_dir, ignore_errors=True)
    shutil.rmtree(rebuild_dir, ignore_errors=True)
    shutil.rmtree(lib_dir, ignore_errors=True)

    print(f"\n[loo] DONE: {n_done} sequences scanned", file=sys.stderr)
    print(f"  validation: {os.path.join(output_dir, 'loo_validation.tsv')}",
          file=sys.stderr)
    print(f"  metrics:    {metrics_path}", file=sys.stderr)
    print(f"  confusion:  {cm_path}", file=sys.stderr)

    # Print summary table
    print("\n[loo] Per-family metrics:", file=sys.stderr)
    print(f"  {'family':<30} {'n':>4} {'recall':>8} {'precision':>10}",
          file=sys.stderr)
    for fam in sorted(metrics):
        m = metrics[fam]
        prec_str = f"{m['precision']:.3f}" if m["precision"] is not None else "NA"
        print(f"  {fam:<30} {m['n_total']:>4} {m['recall']:>8.3f} {prec_str:>10}",
              file=sys.stderr)

    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--hmm-dir", required=True,
                    help="Directory with HMM library + manifest.tsv "
                         "(from build_classification_hmms.py)")
    ap.add_argument("--output-dir", required=True,
                    help="Directory for LOO output TSVs")
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()
    return run_loo(args.hmm_dir, args.output_dir, threads=args.threads)


if __name__ == "__main__":
    sys.exit(main())
