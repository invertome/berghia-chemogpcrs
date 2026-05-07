#!/usr/bin/env python3
"""validate_classification_against_aplysia.py — Phase 6 Task 6.3
retrospective validation harness for the non-chemoreceptor classification.

Reads a ranked candidates CSV (with classification columns from stage 07)
and a curated Aplysia truth CSV; computes per-entry agreement and
overall metrics:
    - chemoreceptor_correctness: fraction of known Aplysia chemoreceptors
      (Cummins 2009) correctly LEFT as 'chemoreceptor-candidate' (TN rate).
      A 'false positive' here = pipeline polluting the wet-lab shortlist
      by mis-flagging a real chemoreceptor as non-chemo.
    - non_chemoreceptor_recall: fraction of known Aplysia non-chemoreceptors
      (5HT receptors, allatostatin receptor, etc.) correctly marked as
      'non-chemoreceptor' or 'likely-non-chemoreceptor'.
    - family_accuracy: among correctly-marked non-chemo, fraction with the
      right family label.
    - subfamily_accuracy: among correctly-family non-chemo, fraction with
      the right subfamily label (only meaningful for aminergic/peptide).

Usage:
    python3 validate_classification_against_aplysia.py \\
        --ranked-csv path/to/aplysia_ranked.csv \\
        --truth-csv references/aplysia_classification_validation.csv \\
        --out results/validation/aplysia_classification_report.md
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import pandas as pd


NON_CHEMO_LABELS = {"non-chemoreceptor", "likely-non-chemoreceptor"}


def load_truth(path: str) -> dict[str, dict]:
    """Read the curated truth CSV; return dict: aplysia_protein_id -> dict
    of expected_classification + expected_family + expected_subfamily."""
    out: dict[str, dict] = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            acc = row.get("aplysia_protein_id", "").strip()
            if not acc:
                continue
            out[acc] = {
                "expected_classification": row.get("expected_classification", "").strip(),
                "expected_family": row.get("expected_family", "").strip(),
                "expected_subfamily": row.get("expected_subfamily", "").strip(),
                "source_citation": row.get("source_citation", "").strip(),
                "notes": row.get("notes", "").strip(),
            }
    return out


def verdict_for_entry(truth: dict, obs: dict) -> dict:
    """Per-entry verdict.

    Returns dict with:
      correct: bool
      error_type: 'TP' (correctly marked non-chemo) | 'TN' (correctly
                  unmarked chemo) | 'FP' (chemoreceptor falsely marked
                  non-chemo) | 'FN' (non-chemo missed)
      family_correct: bool (only when error_type == TP)
      subfamily_correct: bool (only when error_type == TP)
    """
    expected = truth["expected_classification"]
    observed = obs.get("classification", "")
    expected_is_non_chemo = expected == "non-chemoreceptor"
    observed_is_non_chemo = observed in NON_CHEMO_LABELS

    if not expected_is_non_chemo and not observed_is_non_chemo:
        return {"correct": True, "error_type": "TN",
                "family_correct": None, "subfamily_correct": None}
    if expected_is_non_chemo and observed_is_non_chemo:
        fam_ok = (truth["expected_family"]
                   == obs.get("classification_family", ""))
        sub_ok = (truth["expected_subfamily"]
                   == obs.get("classification_subfamily", ""))
        return {"correct": True, "error_type": "TP",
                "family_correct": fam_ok, "subfamily_correct": sub_ok}
    if not expected_is_non_chemo and observed_is_non_chemo:
        return {"correct": False, "error_type": "FP",
                "family_correct": False, "subfamily_correct": False}
    return {"correct": False, "error_type": "FN",
            "family_correct": False, "subfamily_correct": False}


def run_validation(ranked_csv: str, truth_csv: str,
                    out_md_path: str, id_col: str = "id") -> dict:
    """End-to-end: read inputs, compute per-entry verdicts + summary,
    write markdown report. Returns summary metrics dict."""
    truth = load_truth(truth_csv)
    ranked = pd.read_csv(ranked_csv, keep_default_na=False, dtype=str)

    by_id: dict[str, dict] = {}
    if id_col in ranked.columns:
        for _, row in ranked.iterrows():
            by_id[str(row[id_col])] = row.to_dict()

    rows: list[dict] = []
    n_not_classified = 0
    for acc, t in truth.items():
        obs = by_id.get(acc)
        if obs is None:
            n_not_classified += 1
            rows.append({"id": acc, "verdict": "NOT_IN_RANKED",
                         "expected": t["expected_classification"],
                         "observed": "—",
                         "expected_family": t["expected_family"],
                         "observed_family": "—",
                         "notes": t["notes"]})
            continue
        v = verdict_for_entry(t, obs)
        rows.append({
            "id": acc,
            "verdict": v["error_type"],
            "expected": t["expected_classification"],
            "observed": obs.get("classification", ""),
            "expected_family": t["expected_family"],
            "observed_family": obs.get("classification_family", ""),
            "expected_subfamily": t["expected_subfamily"],
            "observed_subfamily": obs.get("classification_subfamily", ""),
            "notes": t["notes"],
            "family_correct": v["family_correct"],
            "subfamily_correct": v["subfamily_correct"],
        })

    n_total = len(truth)
    n_classified = n_total - n_not_classified
    n_correct = sum(1 for r in rows
                     if r["verdict"] in ("TP", "TN"))
    n_tp = sum(1 for r in rows if r["verdict"] == "TP")
    n_fp = sum(1 for r in rows if r["verdict"] == "FP")
    n_tn = sum(1 for r in rows if r["verdict"] == "TN")
    n_fn = sum(1 for r in rows if r["verdict"] == "FN")
    n_chemo_truth = sum(1 for t in truth.values()
                          if t["expected_classification"] == "chemoreceptor-candidate")
    n_non_chemo_truth = sum(1 for t in truth.values()
                              if t["expected_classification"] == "non-chemoreceptor")

    metrics = {
        "n_total": n_total,
        "n_classified": n_classified,
        "n_not_classified": n_not_classified,
        "n_correct": n_correct,
        "n_tp": n_tp, "n_fp": n_fp, "n_tn": n_tn, "n_fn": n_fn,
        "accuracy": n_correct / n_total if n_total else 0,
        "chemoreceptor_correctness": (n_tn / n_chemo_truth
                                       if n_chemo_truth else None),
        "non_chemoreceptor_recall": (n_tp / n_non_chemo_truth
                                      if n_non_chemo_truth else None),
    }

    # Markdown report
    Path(out_md_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_md_path, "w") as f:
        f.write("# Aplysia retrospective classification validation\n\n")
        f.write("Compares non-chemoreceptor classifier output against the "
                "curated Aplysia truth set.\n\n")
        f.write(f"## Summary\n\n")
        f.write(f"- Total truth entries: {n_total}\n")
        f.write(f"- Classified (in ranked CSV): {n_classified}\n")
        f.write(f"- Not classified (missing from ranked): {n_not_classified}\n")
        f.write(f"- Correct verdicts: {n_correct}/{n_total} "
                f"({metrics['accuracy']:.1%})\n")
        f.write(f"- True positives (non-chemo correctly flagged): {n_tp}\n")
        f.write(f"- True negatives (chemoreceptor correctly preserved): {n_tn}\n")
        f.write(f"- False positives (chemoreceptor falsely flagged): {n_fp}\n")
        f.write(f"- False negatives (non-chemo missed): {n_fn}\n\n")
        if metrics["chemoreceptor_correctness"] is not None:
            f.write(f"- Chemoreceptor correctness "
                    f"(TN / known chemoreceptors): "
                    f"{metrics['chemoreceptor_correctness']:.1%} "
                    f"({n_tn}/{n_chemo_truth})\n")
        if metrics["non_chemoreceptor_recall"] is not None:
            f.write(f"- Non-chemoreceptor recall "
                    f"(TP / known non-chemo): "
                    f"{metrics['non_chemoreceptor_recall']:.1%} "
                    f"({n_tp}/{n_non_chemo_truth})\n\n")
        f.write("## Per-entry verdicts\n\n")
        f.write("| ID | Verdict | Expected | Observed | "
                "Expected family | Observed family | Notes |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for r in rows:
            f.write(f"| {r['id']} | {r['verdict']} | {r['expected']} | "
                    f"{r['observed']} | {r.get('expected_family', '')} | "
                    f"{r.get('observed_family', '')} | "
                    f"{r['notes'][:60]} |\n")
    return metrics


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--ranked-csv", required=True,
                    help="Ranked candidates CSV with classification columns "
                         "(stage 07 output)")
    ap.add_argument("--truth-csv", required=True,
                    help="Curated Aplysia truth CSV "
                         "(references/aplysia_classification_validation.csv)")
    ap.add_argument("--out", required=True,
                    help="Output markdown report path")
    ap.add_argument("--id-col", default="id")
    args = ap.parse_args()
    metrics = run_validation(args.ranked_csv, args.truth_csv, args.out,
                              args.id_col)
    print(f"[validate-aplysia] {metrics['n_correct']}/{metrics['n_total']} "
          f"correct ({metrics['accuracy']:.1%})", file=sys.stderr)
    if metrics.get("chemoreceptor_correctness") is not None:
        print(f"  chemoreceptor correctness: "
              f"{metrics['chemoreceptor_correctness']:.1%}", file=sys.stderr)
    if metrics.get("non_chemoreceptor_recall") is not None:
        print(f"  non-chemoreceptor recall:  "
              f"{metrics['non_chemoreceptor_recall']:.1%}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
