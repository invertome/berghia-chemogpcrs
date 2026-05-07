#!/usr/bin/env python3
"""classify_via_og_vote.py — Orthogroup-vote classifier for the
non-chemoreceptor classification pipeline.

Phase 4 Task 4.2.

For each candidate, look up its orthogroup membership, vote on the
most common annotated family among OG members. Require:
  - ≥80% consensus (default; tunable)
  - ≥3 annotated members (default; tunable)

If both conditions met, the OG (and the candidate that belongs to it)
is classified into that family. Subfamily is included only if the
medium-granularity annotations also agree among the same members.
Otherwise the candidate is 'unclassified-og'.

Output: TSV with columns
    candidate_id, og_id, og_vote_family, og_vote_subfamily,
    n_members, n_annotated, consensus_fraction

Usage:
    python3 classify_via_og_vote.py \\
        --candidate-fasta path/to/candidates.fa \\
        --orthogroups-tsv path/to/Orthogroups.tsv \\
        --annotations-tsv references/non_chemo_gpcr/all_references.tsv \\
        --output-tsv results/classification/candidate_og_assignments.tsv
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter
from pathlib import Path


DEFAULT_CONSENSUS_THRESHOLD = 0.8
DEFAULT_MIN_ANNOTATED = 3


# ---- Inputs -------------------------------------------------------------

def load_reference_annotations(annotations_tsv: str) -> dict[str, tuple[str, str]]:
    """Read the curated reference TSV (Phase 1 output). Returns dict:
    accession -> (family, subfamily)."""
    out: dict[str, tuple[str, str]] = {}
    with open(annotations_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            acc = row.get("accession", "").strip()
            fam = row.get("family", "").strip()
            sub = row.get("subfamily", "").strip()
            if acc and fam:
                out[acc] = (fam, sub)
    return out


def load_og_members(orthogroups_tsv: str) -> dict[str, list[str]]:
    """Parse OrthoFinder Orthogroups.tsv -> {og_id: [member_ids]}."""
    og_members: dict[str, list[str]] = {}
    if not os.path.exists(orthogroups_tsv):
        return og_members
    with open(orthogroups_tsv) as f:
        next(f, None)  # header (species names)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            og_id = parts[0].strip()
            if not og_id:
                continue
            members: list[str] = []
            for cell in parts[1:]:
                cell = cell.strip()
                if not cell:
                    continue
                for g in cell.split(","):
                    g = g.strip()
                    if g:
                        members.append(g)
            og_members[og_id] = members
    return og_members


def build_member_to_og(og_members: dict[str, list[str]]
                       ) -> dict[str, str]:
    """Reverse mapping: gene_id -> OG id (each gene belongs to one OG)."""
    out: dict[str, str] = {}
    for og_id, members in og_members.items():
        for m in members:
            out[m] = og_id
    return out


# ---- Voting -------------------------------------------------------------

def vote_for_og(members: list[str],
                annotations: dict[str, tuple[str, str]],
                consensus_threshold: float = DEFAULT_CONSENSUS_THRESHOLD,
                min_annotated: int = DEFAULT_MIN_ANNOTATED
                ) -> dict:
    """Vote on the family + subfamily for an OG given its members and
    the reference annotations dict.

    Returns dict with keys: family, subfamily, n_annotated,
    consensus_fraction, top_count.
    """
    annotated = [(m, annotations[m]) for m in members if m in annotations]
    n_annotated = len(annotated)

    if n_annotated < min_annotated:
        return {
            "family": "unclassified-og",
            "subfamily": "",
            "n_annotated": n_annotated,
            "consensus_fraction": 0.0,
            "top_count": 0,
        }

    # Family vote
    families = [fam for _, (fam, _) in annotated]
    fam_counts = Counter(families)
    top_family, top_count = fam_counts.most_common(1)[0]
    consensus = top_count / n_annotated

    if consensus < consensus_threshold:
        return {
            "family": "unclassified-og",
            "subfamily": "",
            "n_annotated": n_annotated,
            "consensus_fraction": consensus,
            "top_count": top_count,
        }

    # Family is decided. Now check subfamily agreement among the
    # top-family members only.
    sub_for_top = [sub for _, (fam, sub) in annotated
                    if fam == top_family and sub]
    if sub_for_top:
        sub_counts = Counter(sub_for_top)
        top_sub, top_sub_count = sub_counts.most_common(1)[0]
        sub_consensus = top_sub_count / len(sub_for_top)
        # Require same threshold for subfamily call
        subfamily = top_sub if sub_consensus >= consensus_threshold else ""
    else:
        subfamily = ""

    return {
        "family": top_family,
        "subfamily": subfamily,
        "n_annotated": n_annotated,
        "consensus_fraction": consensus,
        "top_count": top_count,
    }


# ---- Orchestration ------------------------------------------------------

def classify_candidates_via_og(
    candidate_ids: list[str],
    orthogroups_tsv: str,
    annotations_tsv: str,
    output_tsv: str,
    consensus_threshold: float = DEFAULT_CONSENSUS_THRESHOLD,
    min_annotated: int = DEFAULT_MIN_ANNOTATED,
) -> None:
    """End-to-end: read inputs, run OG vote per candidate, write TSV."""
    annotations = load_reference_annotations(annotations_tsv)
    og_members = load_og_members(orthogroups_tsv)
    member_to_og = build_member_to_og(og_members)

    # Cache OG votes (same OG → same vote, save work)
    og_vote_cache: dict[str, dict] = {}

    Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(output_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "candidate_id", "og_id",
            "og_vote_family", "og_vote_subfamily",
            "n_members", "n_annotated", "consensus_fraction",
        ])
        for cid in candidate_ids:
            og_id = member_to_og.get(cid, "")
            if not og_id:
                w.writerow([cid, "", "unclassified-og", "", 0, 0, 0.0])
                continue
            if og_id not in og_vote_cache:
                members = og_members.get(og_id, [])
                # Exclude the candidate itself from the vote (we're trying to
                # classify it from its NEIGHBORS, not from itself).
                voting_members = [m for m in members if m != cid]
                og_vote_cache[og_id] = vote_for_og(
                    voting_members, annotations,
                    consensus_threshold=consensus_threshold,
                    min_annotated=min_annotated)
                og_vote_cache[og_id]["n_members"] = len(members)
            v = og_vote_cache[og_id]
            w.writerow([
                cid, og_id,
                v["family"], v["subfamily"],
                v["n_members"], v["n_annotated"],
                f"{v['consensus_fraction']:.3f}",
            ])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--candidate-fasta", required=True,
                    help="FASTA of candidate sequences (we read IDs only)")
    ap.add_argument("--orthogroups-tsv", required=True,
                    help="OrthoFinder Orthogroups.tsv")
    ap.add_argument("--annotations-tsv", required=True,
                    help="Reference annotations TSV "
                         "(references/non_chemo_gpcr/all_references.tsv)")
    ap.add_argument("--output-tsv", required=True)
    ap.add_argument("--consensus-threshold", type=float,
                    default=DEFAULT_CONSENSUS_THRESHOLD)
    ap.add_argument("--min-annotated", type=int,
                    default=DEFAULT_MIN_ANNOTATED)
    args = ap.parse_args()

    # Read candidate IDs
    candidate_ids: list[str] = []
    with open(args.candidate_fasta) as f:
        for line in f:
            if line.startswith(">"):
                candidate_ids.append(line[1:].split()[0])

    classify_candidates_via_og(
        candidate_ids=candidate_ids,
        orthogroups_tsv=args.orthogroups_tsv,
        annotations_tsv=args.annotations_tsv,
        output_tsv=args.output_tsv,
        consensus_threshold=args.consensus_threshold,
        min_annotated=args.min_annotated,
    )

    # Brief summary
    n_classified = 0
    n_total = 0
    with open(args.output_tsv) as f:
        next(f, None)  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                n_total += 1
                if parts[2] != "unclassified-og":
                    n_classified += 1
    print(f"[og-vote] DONE: {n_classified}/{n_total} candidates classified "
          f"-> {args.output_tsv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
