#!/usr/bin/env python3
# structural_evidence.py
# Purpose: Foldseek structural-evidence channel for the ML/PLM chemoreceptor
#   ranking reranker (Task 4, docs/plans/2026-07-01-ml-plm-chemoreceptor-
#   ranking.md). Turns each candidate's AlphaFold model into two honest
#   signals via Foldseek vs PDB+AFDB50+GPCRdb.
"""Parse and score Foldseek structural-similarity search results.

Council rule (bead berghia-chemogpcrs-875): cross-species structural
resemblance is used ONLY for RECALL (flagging structurally-novel
candidates -- no confident hit to any characterized family) and EXCLUSION
(corroborating that a candidate matches a known NON-chemoreceptor GPCR
family, e.g. an opsin or an aminergic receptor). It is explicitly NEVER
used as a positive "looks like a known chemoreceptor" score -- a confident
hit to a family that is not a recognized non-chemoreceptor family
("known_other") is not surfaced as evidence for chemoreceptor identity.

This module is a pure parser/scorer: it reads Foldseek `easy-search`
tabular output (produced on Unity by scripts/unity/run_foldseek_candidates.sh)
plus a target-id -> family-label map, and emits per-candidate flags. No
Foldseek binary is invoked here, and there are no import-time side effects.
"""

from __future__ import annotations

import os
from typing import Optional

# Family-label substrings recognized as NON-chemoreceptor GPCR families.
# A module constant (rather than hard-coded inline) so it is testable and
# overridable by callers without editing this module. Matches the family
# taxonomy the non-chemoreceptor classifier already uses (see
# scripts/_classification_labels.py) -- labels like "aminergic_5HT",
# "class-B-secretin", "class-F-frizzled".
NON_CHEMORECEPTOR_FAMILIES = {
    "peptide", "aminergic", "opsin", "lipid", "nucleotide",
    "class-B", "class-C", "class-F", "bioamine",
}

# foldseek easy-search --format-output "query,target,fident,alntmscore,evalue"
_EXPECTED_FIELDS = 5


def parse_foldseek(path: str) -> dict:
    """Parse Foldseek easy-search tabular output into a best-hit-per-query map.

    Expected columns (whitespace/tab-separated, no header):
        query, target, fident, alntmscore, evalue

    Keeps the hit with the highest alntmscore per query. Blank lines and
    malformed rows (wrong column count, non-numeric fident/alntmscore/evalue,
    empty query/target) are skipped rather than raising. A missing file
    returns {} rather than raising.

    Returns:
        {query_id: {"target": str, "fident": float, "alntmscore": float,
                    "evalue": float}}
    """
    if not path or not os.path.exists(path):
        return {}

    best: dict = {}
    with open(path) as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != _EXPECTED_FIELDS:
                continue
            query, target, fident_s, alntmscore_s, evalue_s = parts
            if not query or not target:
                continue
            try:
                fident = float(fident_s)
                alntmscore = float(alntmscore_s)
                evalue = float(evalue_s)
            except ValueError:
                continue
            current = best.get(query)
            if current is None or alntmscore > current["alntmscore"]:
                best[query] = {
                    "target": target,
                    "fident": fident,
                    "alntmscore": alntmscore,
                    "evalue": evalue,
                }
    return best


def classify_hit(best: Optional[dict], family_map: dict,
                  tm_threshold: float = 0.5) -> str:
    """Classify a candidate's best Foldseek hit into one of three honest states.

    Returns one of:
        "novel" -- best is None, or its alntmscore < tm_threshold (no
            confident structural match to anything characterized). This is
            the RECALL signal.
        "known_non_chemoreceptor" -- confident hit (alntmscore >=
            tm_threshold) whose target family (via family_map) contains one
            of the NON_CHEMORECEPTOR_FAMILIES tags. This is the EXCLUSION
            signal.
        "known_other" -- confident hit, but not a recognized
            non-chemoreceptor family. Deliberately NOT treated as "looks
            like a known chemoreceptor" -- it carries no positive weight.
    """
    if best is None or best["alntmscore"] < tm_threshold:
        return "novel"
    family = family_map.get(best["target"], "") or ""
    if any(tag in family for tag in NON_CHEMORECEPTOR_FAMILIES):
        return "known_non_chemoreceptor"
    return "known_other"


def structural_channel(foldseek_path: str, family_map: dict,
                        tm_threshold: float = 0.5) -> dict:
    """Build the per-candidate Foldseek structural-evidence channel.

    Only emits entries for queries present in the Foldseek output (those get
    has_struct_data=True). A candidate with no Foldseek row at all (e.g. no
    AlphaFold model was generated for it) is simply absent from the returned
    dict -- downstream integration (rank_candidates.py) is responsible for
    defaulting any candidate absent from these keys to has_struct_data=False,
    the same convention used by every other has_*_data signal column in
    ranked_candidates_sorted.csv. This function never fabricates a 0/False
    value for a candidate it has no evidence about.

    Returns:
        {candidate_id: {"struct_state": "novel"|"known_non_chemoreceptor"|
                        "known_other",
                        "struct_novelty": 1 if novel else 0,
                        "struct_nonchemo_corrob": 1 if known_non_chemoreceptor
                                                  else 0,
                        "has_struct_data": True}}
    """
    hits = parse_foldseek(foldseek_path)
    out = {}
    for query, best in hits.items():
        state = classify_hit(best, family_map, tm_threshold=tm_threshold)
        out[query] = {
            "struct_state": state,
            "struct_novelty": 1 if state == "novel" else 0,
            "struct_nonchemo_corrob": 1 if state == "known_non_chemoreceptor" else 0,
            "has_struct_data": True,
        }
    return out
