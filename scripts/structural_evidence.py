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
import re
from typing import List, Optional

# Coarse NON-chemoreceptor GPCR family labels. A module constant (rather than
# hard-coded inline) so it is testable and overridable by callers without
# editing this module.
#
# These are the EXACT canonical coarse families the non-chemoreceptor
# classifier emits: an authoritative mirror of COARSE_FAMILIES in
# scripts/validate_classification_hmms.py (identical to the family list in
# scripts/curate_gpcr_references.py). Keep the three copies in sync. Using the
# classifier's own strings is what lets a Foldseek hit to e.g. the
# glycoprotein-hormone (TSHR/FSHR/LHCGR) or metabotropic-neurotransmitter
# family correctly corroborate EXCLUSION instead of silently leaking away as
# "known_other". family_is_non_chemoreceptor() matches these against a
# family_map label using that vocabulary's own grammar -- the whole label
# ("glycoprotein-hormone") or its "<coarse>_<subfamily>" coarse head
# ("aminergic_5HT" -> "aminergic"). The match is ANCHORED, never a substring
# test: `"opsin" in "Class A (Rhodopsin)"` is True, so a substring test would
# classify every class-A hit -- i.e. our chemoreceptor candidates -- as a
# non-chemoreceptor the moment a vocabulary using GPCRdb-style labels (e.g. via
# --gpcrdb-meta) reached this function. NB: "unclassified-gpcr" is the
# uncertain/candidate bucket, NOT a non-chemoreceptor family -- omitted on
# purpose.
NON_CHEMORECEPTOR_FAMILIES = {
    "aminergic", "peptide", "opsin", "lipid", "nucleotide",
    "metabotropic-neurotransmitter", "glycoprotein-hormone",
    "class-B-secretin", "class-C", "class-F-frizzled",
}

# foldseek easy-search --format-output "query,target,fident,alntmscore,evalue"
_EXPECTED_FIELDS = 5

# --------------------------------------------------------------------------- #
# Foldseek target-identifier parsing
#
# Foldseek targets are STRUCTURE identifiers, not sequence accessions. The
# family_map this module is handed is keyed on UniProt accessions (from
# references/anchors/anchor_set.tsv) and/or on literal target ids (from an
# optional --gpcrdb-meta table). Looking the raw target up directly therefore
# missed on every real hit, which silently disabled the exclusion signal
# entirely. target_keys() bridges the two by deriving the ordered set of keys a
# given target could legitimately be filed under, per DB family:
#
#   AFDB / AFDB50   AF-P31356-F1-model_v4   -> UniProt accession P31356
#                   (fragment number and model version both vary; a createdb'd
#                    directory may keep .cif/.pdb/.gz and append the chain
#                    foldseek assigned)
#   PDB             1f88_A                  -> 1f88_A, then 1f88
#                   (foldseek's prebuilt PDB db names entries <pdbid>_<chain>)
#   GPCRdb / local  <file basename>         -> the literal target, then the
#                    extension-stripped basename
#
# The raw target is always tried FIRST so an explicit meta mapping keyed on
# exactly what foldseek printed always wins over anything derived.
# --------------------------------------------------------------------------- #

# Compression / structure-file extensions, stripped repeatedly (".cif.gz").
_STRUCTURE_SUFFIXES = (".gz", ".bcif", ".mmcif", ".cif", ".pdb", ".ent")

# AF-<accession>-F<fragment>-model_v<version>, optionally + "_<chain>".
_AFDB_RE = re.compile(
    r"^AF-([A-Za-z0-9]+)-F\d+-model_v\d+(?:_[A-Za-z0-9]{1,4})?$")

# <4-char pdb id starting with a digit>_<chain>
_PDB_RE = re.compile(r"^(\d[A-Za-z0-9]{3})_([A-Za-z0-9]{1,4})$")


def _strip_structure_suffixes(name: str) -> str:
    """Repeatedly strip structure/compression extensions ("x.cif.gz" -> "x")."""
    while True:
        lowered = name.lower()
        for suffix in _STRUCTURE_SUFFIXES:
            if lowered.endswith(suffix) and len(name) > len(suffix):
                name = name[: -len(suffix)]
                break
        else:
            return name


def target_keys(target: Optional[str]) -> List[str]:
    """Ordered candidate family_map lookup keys for one Foldseek target id.

    The raw target comes first (so a meta table keyed on the literal target
    always wins), then the path basename, then the extension-stripped
    basename, then any DB-specific identifier parsed out of it: the UniProt
    accession for an AFDB/AFDB50 entry, the bare PDB id for a PDB entry.

    Deduplicated and order-stable. An empty/None target yields [].
    """
    raw = (target or "").strip()
    if not raw:
        return []

    keys: List[str] = []

    def add(key: str) -> None:
        if key and key not in keys:
            keys.append(key)

    add(raw)
    basename = raw.rsplit("/", 1)[-1]
    add(basename)
    stem = _strip_structure_suffixes(basename)
    add(stem)

    afdb = _AFDB_RE.match(stem)
    if afdb:
        # AFDB entries embed the UniProt accession the anchor set is keyed on.
        add(afdb.group(1))
        return keys

    pdb = _PDB_RE.match(stem)
    if pdb:
        pdb_id, chain = pdb.group(1), pdb.group(2)
        add(pdb_id)
        add(pdb_id.lower())
        add(pdb_id.upper())
        add(f"{pdb_id.lower()}_{chain}")
        add(f"{pdb_id.upper()}_{chain}")
    return keys


def family_for_target(target: Optional[str], family_map: dict) -> str:
    """Family label for a Foldseek target, or "" if none of its keys resolve.

    Tries target_keys() in order and returns the first non-empty label found.
    """
    for key in target_keys(target):
        label = (family_map.get(key) or "").strip()
        if label:
            return label
    return ""


def family_is_non_chemoreceptor(family: Optional[str],
                                 families: Optional[set] = None) -> bool:
    """True when a family label names a known NON-chemoreceptor GPCR family.

    ANCHORED, not a substring test (see NON_CHEMORECEPTOR_FAMILIES): the label
    must either be a family name outright ("opsin", "class-C") or have one as
    its "<coarse>_<subfamily>" coarse head ("aminergic_5HT" -> "aminergic"). A
    free-text label that merely CONTAINS a family name -- "Class A
    (Rhodopsin)" contains "opsin" -- is correctly rejected.
    """
    vocabulary = NON_CHEMORECEPTOR_FAMILIES if families is None else families
    label = (family or "").strip()
    if not label:
        return False
    if label in vocabulary:
        return True
    return label.split("_", 1)[0] in vocabulary


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
            tm_threshold) whose target resolves (via family_for_target, which
            parses the real Foldseek target-id formats) to one of the
            NON_CHEMORECEPTOR_FAMILIES. This is the EXCLUSION signal.
        "known_other" -- confident hit, but not a recognized
            non-chemoreceptor family. Deliberately NOT treated as "looks
            like a known chemoreceptor" -- it carries no positive weight.
    """
    if best is None or best["alntmscore"] < tm_threshold:
        return "novel"
    family = family_for_target(best["target"], family_map)
    if family_is_non_chemoreceptor(family):
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
