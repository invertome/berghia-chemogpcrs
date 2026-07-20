#!/usr/bin/env python3
"""classify_consensus.py — Three-source consensus classifier.

Phase 4 Task 4.4 of the non-chemoreceptor classification pipeline.

Joins per-candidate output from three independent classifiers:
  - HMM scan          (Task 4.1, scripts/classify_via_hmm.py)
  - Orthogroup vote   (Task 4.2, scripts/classify_via_og_vote.py)
  - Phylogenetic placement (Task 4.3, optional)

For each candidate, applies the consensus logic settled in the design:

    3-of-3 family agreement
        -> classification = 'non-chemoreceptor'
           classification_confidence = 'high'
    2-of-3 family agreement (HMM + OG agree; placement disagrees or absent)
        -> classification = 'likely-non-chemoreceptor'
           classification_confidence = 'medium'
    <2 agreement
        -> classification = 'chemoreceptor-candidate' (default)
           classification_confidence = 'NA'

Subfamily is emitted only when ALL agreeing sources also agree on the
medium-granularity subfamily.

Output columns (joined into ranked CSV by add_classification_columns.py
in stage 07):
    candidate_id, classification, classification_confidence,
    classification_family, classification_subfamily, classification_evidence

Usage:
    python3 classify_consensus.py \\
        --hmm-tsv results/classification/candidate_hmm_assignments.tsv \\
        --og-tsv results/classification/candidate_og_assignments.tsv \\
        --placement-tsv results/classification/candidate_placement.tsv \\
        --out results/classification/candidate_classifications.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


# Sentinel labels emitted by upstream classifiers when they decline to call.
UNCLASSIFIED_LABELS = {"unclassified-hmm", "unclassified-og",
                        "unclassified-placement", ""}


# classify_via_hmm.py tags calls made against the unvalidated fallback
# threshold. Such a call is ~100 orders of magnitude looser than any
# leave-one-out-benchmarked family, and the Pfam library that produces it
# also annotates the reference set the OG vote reads — so HMM and OG agree
# by construction rather than by evidence. Those calls are reported but may
# not carry a candidate into a demotion tier.
UNVALIDATED_THRESHOLD_SOURCE = "unvalidated-default"


def _is_real_call(family: str) -> bool:
    """Return True if `family` is a non-empty, non-unclassified label."""
    return bool(family) and family not in UNCLASSIFIED_LABELS


def _consensus(hmm_fam: str, og_fam: str, placement_fam: str | None,
               hmm_validated: bool = True
               ) -> tuple[str, str, int]:
    """Determine consensus from up to three family calls.

    Returns (consensus_family, classification_label, n_agree) where:
        consensus_family: the family for a high/medium call; "" otherwise
        classification_label:
            'non-chemoreceptor'         — all three sources agree (high)
            'likely-non-chemoreceptor'  — HMM + OG agree (medium)
            'chemoreceptor-candidate'   — otherwise (default)
        n_agree: 3 (high), 2 (medium), 1 (some real call, no override), 0 (none)

    The medium tier requires the two ANNOTATION-based sources (HMM + OG) to
    agree. Phylogenetic placement alone cannot substitute for either: placement
    of divergent invertebrate LSE sequences is the least reliable of the three
    signals, so an HMM+placement or OG+placement pair (with the HMM/OG pair
    itself split) stays a conservative 'chemoreceptor-candidate'.

    `hmm_validated` is False when the HMM call came from the unvalidated
    fallback threshold. Such a call cannot contribute to either demotion
    tier: demoting a real chemoreceptor on unbenchmarked evidence is the
    expensive error for this pipeline. It still appears in the evidence
    string so a reviewer can see it was found and discounted.
    """
    if not hmm_validated:
        hmm_fam = ""
    hmm_real = _is_real_call(hmm_fam)
    og_real = _is_real_call(og_fam)
    placement_real = placement_fam is not None and _is_real_call(placement_fam)

    if not (hmm_real or og_real or placement_real):
        return ("", "chemoreceptor-candidate", 0)

    # High: all three sources present and agree on the same family.
    if hmm_real and og_real and placement_real and hmm_fam == og_fam == placement_fam:
        return (hmm_fam, "non-chemoreceptor", 3)

    # Medium: HMM + OG (the annotation-based sources) agree.
    if hmm_real and og_real and hmm_fam == og_fam:
        return (hmm_fam, "likely-non-chemoreceptor", 2)

    return ("", "chemoreceptor-candidate", 1)


def _consensus_subfamily(consensus_family: str,
                         hmm_fam: str, hmm_sub: str,
                         og_fam: str, og_sub: str,
                         placement_fam: str | None, placement_sub: str | None
                         ) -> str:
    """Return the subfamily when ALL sources that agreed on
    consensus_family ALSO agree on the same non-empty subfamily.
    Otherwise empty string.
    """
    if not consensus_family:
        return ""
    sub_calls: list[str] = []
    if hmm_fam == consensus_family and hmm_sub:
        sub_calls.append(hmm_sub)
    if og_fam == consensus_family and og_sub:
        sub_calls.append(og_sub)
    if (placement_fam is not None
            and placement_fam == consensus_family
            and placement_sub):
        sub_calls.append(placement_sub)
    if not sub_calls:
        return ""
    if len(set(sub_calls)) == 1:
        return sub_calls[0]
    return ""


def _evidence_string(hmm_fam: str, og_fam: str,
                     placement_fam: str | None,
                     hmm_validated: bool = True) -> str:
    hmm_part = hmm_fam or "NA"
    if hmm_fam and not hmm_validated:
        hmm_part = f"{hmm_fam}(unvalidated-threshold)"
    parts = [f"hmm:{hmm_part}", f"og:{og_fam or 'NA'}"]
    if placement_fam is not None:
        parts.append(f"placement:{placement_fam or 'NA'}")
    return ";".join(parts)


def _read_tsv_safe(path: str | None) -> pd.DataFrame:
    """Read a TSV; return empty DataFrame if file missing or empty.
    Empty cells are kept as empty strings (not NaN) so the consensus
    logic can do string comparisons without NaN handling everywhere."""
    if not path or not Path(path).exists():
        return pd.DataFrame()
    if Path(path).stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", keep_default_na=False,
                           na_values=[""], dtype=str).fillna("")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def run_consensus(hmm_tsv: str, og_tsv: str,
                  placement_tsv: str | None, out_tsv: str
                  ) -> None:
    """Read inputs, compute per-candidate consensus, write output TSV."""
    hmm = _read_tsv_safe(hmm_tsv)
    og = _read_tsv_safe(og_tsv)
    # Distinguish "placement requested but absent" from "file empty/missing"
    placement = _read_tsv_safe(placement_tsv) if placement_tsv else None
    if placement is not None and placement.empty:
        # File was given but had no rows — treat as no placement evidence
        # for ALL candidates rather than falling back to "no placement
        # support".
        placement = None

    # Set of all candidate IDs across the inputs
    all_ids: set[str] = set()
    for df in (hmm, og, placement):
        if df is not None and not df.empty and "candidate_id" in df.columns:
            all_ids.update(df["candidate_id"].astype(str))

    # Build lookup dicts
    def _row(df: pd.DataFrame, cid: str) -> dict:
        if df is None or df.empty:
            return {}
        sel = df[df["candidate_id"].astype(str) == cid]
        if sel.empty:
            return {}
        return sel.iloc[0].to_dict()

    out_rows: list[dict] = []
    for cid in sorted(all_ids):
        h = _row(hmm, cid)
        o = _row(og, cid)
        p = _row(placement, cid) if placement is not None else {}

        hmm_fam = str(h.get("hmm_family", ""))
        hmm_sub = str(h.get("hmm_subfamily", ""))
        # Absent column => written before threshold provenance existed;
        # keep the historical meaning rather than silently dropping every
        # demotion on an old artifact.
        hmm_validated = (str(h.get("threshold_source", ""))
                         != UNVALIDATED_THRESHOLD_SOURCE)
        og_fam = str(o.get("og_vote_family", ""))
        og_sub = str(o.get("og_vote_subfamily", ""))
        if placement is not None:
            placement_fam = str(p.get("placement_family", "")) if p else ""
            placement_sub = str(p.get("placement_subfamily", "")) if p else ""
        else:
            placement_fam = None
            placement_sub = None

        consensus_family, classification, n_agree = _consensus(
            hmm_fam, og_fam, placement_fam, hmm_validated)
        subfamily = _consensus_subfamily(
            consensus_family, hmm_fam if hmm_validated else "", hmm_sub,
            og_fam, og_sub, placement_fam, placement_sub)
        confidence = ("high" if n_agree >= 3
                      else "medium" if n_agree >= 2
                      else "NA")
        evidence = _evidence_string(hmm_fam, og_fam, placement_fam,
                                    hmm_validated)

        out_rows.append({
            "candidate_id": cid,
            "classification": classification,
            "classification_confidence": confidence,
            "classification_family": consensus_family,
            "classification_subfamily": subfamily,
            "classification_evidence": evidence,
        })

    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(out_rows).to_csv(out_tsv, sep="\t", index=False)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--hmm-tsv", required=True,
                    help="HMM scan output TSV (Task 4.1)")
    ap.add_argument("--og-tsv", required=True,
                    help="OG vote output TSV (Task 4.2)")
    ap.add_argument("--placement-tsv", default="",
                    help="Phylogenetic placement output TSV (Task 4.3, "
                         "optional — without it, only HMM+OG 2-source "
                         "consensus is possible, capping at medium confidence)")
    ap.add_argument("--out", required=True,
                    help="Output TSV with consensus classification per candidate")
    args = ap.parse_args()
    run_consensus(args.hmm_tsv, args.og_tsv,
                  args.placement_tsv if args.placement_tsv else None,
                  args.out)

    df = pd.read_csv(args.out, sep="\t")
    n = len(df)
    n_high = int((df["classification_confidence"] == "high").sum())
    n_med = int((df["classification_confidence"] == "medium").sum())
    n_chem = int((df["classification"] == "chemoreceptor-candidate").sum())
    print(f"[consensus] {n} candidates total", file=sys.stderr)
    print(f"  non-chemoreceptor (high):     {n_high}", file=sys.stderr)
    print(f"  likely-non-chemo (medium):    {n_med}", file=sys.stderr)
    print(f"  chemoreceptor-candidate:      {n_chem}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
