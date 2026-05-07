"""Shared helpers for the non-chemoreceptor classification pipeline.

Phase 4-5 of the feature has multiple scripts that all need to:
  1. Map an HMM file basename (no extension) to (family, subfamily) labels.
     Naming convention:
        coarse: '<family>'                 e.g. aminergic
        medium: '<family>_<subfamily>'    e.g. aminergic_5HT
     Family names contain hyphens (class-B-secretin, class-F-frizzled)
     but never underscores. Subfamily names CAN contain hyphens
     (NPY-NPF, vasopressin-oxytocin). Split on FIRST underscore only.

  2. Read candidate IDs from a FASTA (the first whitespace token after '>').

These helpers were duplicated across classify_via_hmm.py,
validate_classification_hmms.py, classify_via_og_vote.py, and
classify_via_placement.py. Consolidated here to fix H4 from the
2026-05-07 code review (preventing label-namespace drift between
modules — the kind of contract bug that bit us with the Pfam HMM
rename collision).
"""
from __future__ import annotations


def label_for_hmm(hmm_name: str) -> tuple[str, str]:
    """Parse an HMM name to (family, subfamily). Split on FIRST underscore."""
    if "_" not in hmm_name:
        return (hmm_name, "")
    family, _, subfamily = hmm_name.partition("_")
    return (family, subfamily)


def read_candidate_ids_from_fasta(fasta_path: str) -> list[str]:
    """Read every candidate ID (first whitespace token after '>') from a FASTA."""
    ids: list[str] = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                tok = line[1:].split()[0] if len(line) > 1 else ""
                if tok:
                    ids.append(tok)
    return ids
