#!/usr/bin/env python3
"""identify_gpcrs.py — HMM-first GPCR identification for stage 02.

Bead -m1f restructure: stage 02 used to run TMbed on the full 86k-protein
Berghia transcriptome and post-filter. That was 99%+ wasted GPU compute
and forced an ad-hoc length filter to dodge ProtT5's quadratic-memory
slow-tail. The correct architecture is to HMM-filter FIRST (catches all
GPCR families: chemoreceptor LSEs + curated bioamine/peptide/opsin/lipid/
nucleotide/class-B/C/F + Pfam 7tm_1/2/3/F fallback), then run TMbed only
on the ~500-3000 GPCR-positive proteins.

This module owns the post-HMMER merge step:

  Inputs:
    - classification TSV (from scripts/classify_via_hmm.py)
    - hmmsearch --tblout against results/hmms/lse.hmm

  Outputs:
    - flat ID list of all GPCR-positive proteins (sorted, unique)
      -> downstream seqtk subseq + TMbed
    - census TSV with seq_id, family, subfamily, evalue, source
      -> follow-up Berghia GPCR/brain-expression paper

Precedence rule: classification (curated, specific family/subfamily)
beats lse-only (project-built chemoreceptor LSE OG HMMs). A protein
hit by both is annotated with the classification family but tagged
source='classification+lse' so the lse evidence isn't lost.
"""
from __future__ import annotations

import argparse
import csv
import sys
from typing import Optional


def _parse_classification(path: str) -> dict[str, dict]:
    """Return {seq_id: {family, subfamily, evalue, source}} for rows in
    the classification TSV that have a real family assignment.
    Drops 'unclassified-hmm' and empty-family rows."""
    out: dict[str, dict] = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            fam = (row.get("hmm_family") or "").strip()
            if not fam or fam == "unclassified-hmm":
                continue
            seq_id = row.get("candidate_id", "").strip()
            if not seq_id:
                continue
            ev_str = (row.get("evalue") or "").strip()
            try:
                ev = float(ev_str) if ev_str else None
            except ValueError:
                ev = None
            out[seq_id] = {
                "family": fam,
                "subfamily": (row.get("hmm_subfamily") or "").strip(),
                "evalue": ev,
                "source": "classification",
            }
    return out


def _parse_lse_tblout(path: str) -> dict[str, dict]:
    """Return {seq_id: {family='lse_chemoreceptor', subfamily=<OG>,
    evalue=<best>, source='lse'}} keeping the best (lowest-E) hit per
    query.

    HMMER --tblout columns: target query - - eval(full) score(full) ...
    """
    best: dict[str, tuple[str, float]] = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            target = parts[0]
            query = parts[2]
            try:
                ev = float(parts[4])
            except ValueError:
                continue
            cur = best.get(query)
            if cur is None or ev < cur[1]:
                best[query] = (target, ev)
    return {
        q: {
            "family": "lse_chemoreceptor",
            "subfamily": og,
            "evalue": ev,
            "source": "lse",
        }
        for q, (og, ev) in best.items()
    }


def merge_gpcr_evidence(
    classification_tsv: Optional[str],
    lse_tblout: Optional[str],
) -> dict[str, dict]:
    """Union classification + lse evidence per sequence.
    Classification family/subfamily wins on overlap; source tagged
    'classification+lse' when both sources hit the same seq_id."""
    cls = _parse_classification(classification_tsv) if classification_tsv else {}
    lse = _parse_lse_tblout(lse_tblout) if lse_tblout else {}
    merged: dict[str, dict] = {}
    for sid, rec in cls.items():
        merged[sid] = dict(rec)
    for sid, rec in lse.items():
        if sid in merged:
            merged[sid]["source"] = "classification+lse"
        else:
            merged[sid] = dict(rec)
    return merged


def write_ids(merged: dict[str, dict], output_path: str) -> None:
    """One sorted ID per line."""
    with open(output_path, "w") as out:
        for sid in sorted(merged):
            out.write(sid + "\n")


def write_census(merged: dict[str, dict], output_path: str) -> None:
    """TSV: seq_id, family, subfamily, evalue, source — sorted by seq_id."""
    with open(output_path, "w") as out:
        out.write("seq_id\tfamily\tsubfamily\tevalue\tsource\n")
        for sid in sorted(merged):
            r = merged[sid]
            ev = r.get("evalue")
            ev_str = f"{ev:.2e}" if isinstance(ev, float) else ""
            out.write(f"{sid}\t{r.get('family','')}\t{r.get('subfamily','')}"
                      f"\t{ev_str}\t{r.get('source','')}\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--classification-tsv",
                    help="Output from scripts/classify_via_hmm.py")
    ap.add_argument("--lse-tblout",
                    help="hmmsearch --tblout against results/hmms/lse.hmm")
    ap.add_argument("--ids-out", required=True,
                    help="Output: GPCR-positive sequence IDs, one per line")
    ap.add_argument("--census-out",
                    help="Output: census TSV (seq_id, family, subfamily, "
                         "evalue, source)")
    args = ap.parse_args()
    if not args.classification_tsv and not args.lse_tblout:
        print("identify_gpcrs.py: at least one of --classification-tsv or "
              "--lse-tblout must be provided", file=sys.stderr)
        return 2
    merged = merge_gpcr_evidence(args.classification_tsv, args.lse_tblout)
    write_ids(merged, args.ids_out)
    if args.census_out:
        write_census(merged, args.census_out)
    print(f"identify_gpcrs: {len(merged)} GPCR-positive proteins", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
