#!/usr/bin/env python3
"""identify_gpcrs.py — HMM-first GPCR identification for stage 02.

Bead -m1f restructure: stage 02 used to run TMbed on the full 86k-protein
Berghia transcriptome and post-filter. That was 99%+ wasted GPU compute
and forced an ad-hoc length filter to dodge ProtT5's quadratic-memory
slow-tail. The correct architecture is to HMM-filter FIRST, then run
TMbed only on the ~500-3000 GPCR-positive proteins.

This module owns the post-HMMER merge step. Three HMM evidence layers:

  - classification TSV (from scripts/classify_via_hmm.py): curated
    Swiss-Prot GPCR family HMMs (bioamine, peptide, opsin, lipid,
    nucleotide, class-B/C/F) + Pfam fallback (7tm_1/2/3, Frizzled).
    Per-family LOO-validated thresholds. NON-chemoreceptor families.
    -> family = specific family/subfamily name, source = 'classification'.

  - lse tblout: hmmsearch against results/hmms/lse.hmm — Nath et al.
    LINEAGE-SPECIFIC EXPANSION HMMs (one OG per cluster from
    references/nath_et_al/lse/). Expansion pattern alone does NOT
    confirm chemoreceptor function — these are uncharacterized GPCRs
    that happen to have expanded in specific lineages.
    -> family = 'lse', subfamily = OG name, source = 'lse'.

  - nath_ortholog tblout: hmmsearch against results/hmms/conserved.hmm —
    Nath et al. ONE-TO-ONE ORTHOLOGOUS GPCR HMMs (from
    references/nath_et_al/one_to_one_ortholog/). Broad metazoan GPCR
    orthologs; NOT chemoreceptor-specific.
    -> family = 'nath_ortholog', subfamily = OG name, source = 'nath_ortholog'.

(Future) curated_chemoreceptor: HMMs built from functionally
confirmed invertebrate chemoreceptor GPCRs (Cummins 2009 Aplysia,
deorphanized C. elegans SR families, etc.). Bead in progress.

Outputs:
  - flat ID list of all GPCR-positive proteins (sorted, unique)
    -> downstream seqtk subseq + TMbed
  - census TSV with seq_id, family, subfamily, evalue, source
    -> follow-up Berghia GPCR/brain-expression paper

Precedence rule: classification has named curated family assignments
so it wins on family/subfamily. lse and nath_ortholog only contribute
their tag when classification is absent. Source tag concatenates every
source that hit the protein (e.g., 'classification+lse+nath_ortholog'),
so no evidence is lost on the multi-hit case.
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


def _parse_hmm_tblout(path: str, family_tag: str) -> dict[str, dict]:
    """Generic HMMER --tblout parser keeping the best (lowest-E) hit per
    sequence. Used for both lse.hmm and conserved.hmm (Nath ortholog)
    scans, both run via hmmsearch.

    hmmsearch --tblout column order (DIFFERENT from hmmscan!):
      col 1 (target name)  = the SEQUENCE searched (FASTA query)
      col 3 (query name)   = the HMM profile name (DB record)
      col 5                = E-value (full sequence)

    family_tag is the value put in the 'family' field AND in 'source'
    (e.g., 'lse' or 'nath_ortholog'). Subfamily is the best-hit HMM name.
    """
    best: dict[str, tuple[str, float]] = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            seq_id = parts[0]    # hmmsearch: target = sequence
            hmm_name = parts[2]  # hmmsearch: query  = HMM profile
            try:
                ev = float(parts[4])
            except ValueError:
                continue
            cur = best.get(seq_id)
            if cur is None or ev < cur[1]:
                best[seq_id] = (hmm_name, ev)
    return {
        sid: {
            "family": family_tag,
            "subfamily": hmm_name,
            "evalue": ev,
            "source": family_tag,
        }
        for sid, (hmm_name, ev) in best.items()
    }


def _parse_curated_chemo_tblout(path: str) -> dict[str, dict]:
    """Curated invertebrate chemoreceptor HMM tblout. Each HMM in this
    library represents a single functionally validated chemoreceptor
    family (e.g., apgr from Aplysia Cummins 2009, schisto_chemo from
    deorphanized platyhelminth receptors, etc.). The hit's HMM name IS
    the family name. family='<hmm_name>', subfamily='<hmm_name>'
    (per-family granularity is what we have), source='curated_chemo'.

    This is the HIGHEST CONFIDENCE detection layer in identify_gpcrs.py
    because every HMM is built from functionally validated sequences,
    not just sequence-similarity-based annotation.

    hmmsearch --tblout column order: col 1 = SEQUENCE (target),
    col 3 = HMM profile (query). Keeps best (lowest-E) hit per sequence.
    """
    best: dict[str, tuple[str, float]] = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            seq_id = parts[0]
            hmm_name = parts[2]
            try:
                ev = float(parts[4])
            except ValueError:
                continue
            cur = best.get(seq_id)
            if cur is None or ev < cur[1]:
                best[seq_id] = (hmm_name, ev)
    return {
        sid: {
            "family": hmm_name,
            "subfamily": hmm_name,
            "evalue": ev,
            "source": "curated_chemo",
        }
        for sid, (hmm_name, ev) in best.items()
    }


def _parse_tiammat_tblout(path: str) -> dict[str, dict]:
    """TIAMMAT-revised GPCR HMM tblout (from hmmsearch). family='gpcr_7tm'
    (broad 7TM GPCR category), subfamily=<HMM name, e.g., 7tm_1_REV or
    7TM_GPCR_Sra_REV>, source='tiammat'.

    hmmsearch --tblout column order: col 1 = SEQUENCE (target),
    col 3 = HMM profile (query). Keeps best (lowest-E) hit per sequence.
    """
    best: dict[str, tuple[str, float]] = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            seq_id = parts[0]    # hmmsearch: target = sequence
            hmm_name = parts[2]  # hmmsearch: query  = HMM profile
            try:
                ev = float(parts[4])
            except ValueError:
                continue
            cur = best.get(seq_id)
            if cur is None or ev < cur[1]:
                best[seq_id] = (hmm_name, ev)
    return {
        sid: {
            "family": "gpcr_7tm",
            "subfamily": hmm_name,
            "evalue": ev,
            "source": "tiammat",
        }
        for sid, (hmm_name, ev) in best.items()
    }


def merge_gpcr_evidence(
    classification_tsv: Optional[str],
    lse_tblout: Optional[str],
    nath_ortholog_tblout: Optional[str] = None,
    tiammat_tblout: Optional[str] = None,
    curated_chemo_tblout: Optional[str] = None,
) -> dict[str, dict]:
    """Union curated_chemo + classification + tiammat + lse +
    nath_ortholog evidence per sequence.

    Precedence order on family/subfamily (most specific / highest
    confidence first):
      1. curated_chemo — functionally validated invertebrate
         chemoreceptors (Cummins 2009 Aplysia, deorphanized platyhelminth
         GPCRs, etc.). HIGHEST CONFIDENCE detection.
      2. classification — curated Swiss-Prot families (LOO-validated
         but non-chemoreceptor families like bioamine, peptide, opsin)
      3. tiammat — TIAMMAT-revised 7TM GPCR HMMs (mollusc-tuned Pfam +
         invertebrate Sr-style)
      4. lse — Nath et al. lineage-specific expansion OGs (generic tag)
      5. nath_ortholog — Nath et al. one_to_one_ortholog OGs (generic)

    Why curated_chemo beats classification: a protein hitting both
    apgr (functionally validated rhinophore chemoreceptor) AND
    aminergic (Swiss-Prot bioamine family signature) is clearly a
    chemoreceptor that happens to share some structural features with
    bioamine GPCRs — the chemoreceptor identity is the BIOLOGICALLY
    INTERESTING one for our paper.

    The source field concatenates every source that contributed evidence,
    e.g., 'curated_chemo+classification+tiammat+lse+nath_ortholog'.
    """
    cur = _parse_curated_chemo_tblout(curated_chemo_tblout) \
        if curated_chemo_tblout else {}
    cls = _parse_classification(classification_tsv) if classification_tsv else {}
    tia = _parse_tiammat_tblout(tiammat_tblout) if tiammat_tblout else {}
    lse = _parse_hmm_tblout(lse_tblout, "lse") if lse_tblout else {}
    nath = _parse_hmm_tblout(nath_ortholog_tblout, "nath_ortholog") \
        if nath_ortholog_tblout else {}

    merged: dict[str, dict] = {}
    sources_seen: dict[str, list[str]] = {}

    for source_name, source_dict in (("curated_chemo", cur),
                                      ("classification", cls),
                                      ("tiammat", tia),
                                      ("lse", lse),
                                      ("nath_ortholog", nath)):
        for sid, rec in source_dict.items():
            if sid not in merged:
                merged[sid] = dict(rec)
                sources_seen[sid] = [source_name]
            else:
                sources_seen[sid].append(source_name)
                # First-set wins family/subfamily — loop order enforces
                # the documented precedence.

    # Rewrite source field as the canonical concatenation
    for sid, srcs in sources_seen.items():
        merged[sid]["source"] = "+".join(srcs)
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
    ap.add_argument("--nath-ortholog-tblout",
                    help="hmmsearch --tblout against results/hmms/conserved.hmm "
                         "(Nath et al. one_to_one_ortholog GPCR-OG HMMs)")
    ap.add_argument("--tiammat-tblout",
                    help="hmmsearch --tblout against TIAMMAT-revised GPCR "
                         "HMMs (Pfam 7tm_* + invertebrate Sr-style "
                         "chemoreceptors, retrained on mollusc data)")
    ap.add_argument("--curated-chemo-tblout",
                    help="hmmsearch --tblout against curated invertebrate "
                         "chemoreceptor HMMs (functionally validated; "
                         "Cummins 2009 Aplysia + deorphanized lophotrochozoan "
                         "chemoreceptor families)")
    ap.add_argument("--ids-out", required=True,
                    help="Output: GPCR-positive sequence IDs, one per line")
    ap.add_argument("--census-out",
                    help="Output: census TSV (seq_id, family, subfamily, "
                         "evalue, source)")
    args = ap.parse_args()
    if not any((args.classification_tsv, args.lse_tblout,
                args.nath_ortholog_tblout, args.tiammat_tblout,
                args.curated_chemo_tblout)):
        print("identify_gpcrs.py: at least one of --classification-tsv, "
              "--lse-tblout, --nath-ortholog-tblout, --tiammat-tblout, "
              "or --curated-chemo-tblout must be provided", file=sys.stderr)
        return 2
    merged = merge_gpcr_evidence(args.classification_tsv, args.lse_tblout,
                                  nath_ortholog_tblout=args.nath_ortholog_tblout,
                                  tiammat_tblout=args.tiammat_tblout,
                                  curated_chemo_tblout=args.curated_chemo_tblout)
    write_ids(merged, args.ids_out)
    if args.census_out:
        write_census(merged, args.census_out)
    print(f"identify_gpcrs: {len(merged)} GPCR-positive proteins", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
