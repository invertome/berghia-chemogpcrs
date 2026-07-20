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
import sys
from typing import Dict, List, Optional, Set

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


# --------------------------------------------------------------------------- #
# Foldseek QUERY-identifier normalisation
#
# Symmetric with target_keys() above, and deliberately living in the same module
# so the two sides share _strip_structure_suffixes() and cannot drift.
#
# structural_channel() used the raw foldseek `query` field VERBATIM as the key
# of the channel it returns. A raw search-tool output field is not a join key.
#
# MEASURED, not assumed (foldseek 10.941cd33, this project's berghia-gpcr env on
# Unity, srun 61999969; flat query dir staged <candidate_id>.<ext> exactly as
# scripts/unity/run_foldseek_candidates.sh stages it):
#
#     BersteEVm000001t1.cif  (1 chain)   -> query "BersteEVm000001t1"
#     BersteEVm000002t1.pdb  (1 chain)   -> query "BersteEVm000002t1"
#     BersteEVm000003t1.cif  (2 chains)  -> queries "BersteEVm000003t1_A"
#                                                   "BersteEVm000003t1_B"
#
# So foldseek STRIPS the extension (structcreatedb.cpp: Util::remove_extension,
# applied twice for .gz/.zst) and APPENDS "_<chain>" once a file holds more than
# one chain, emitting one row PER CHAIN. The chain suffix is the real break: a
# multi-chain AF3 model (receptor + peptide/ligand/G-protein, or any complex)
# yields <cand_id>_A / <cand_id>_B, neither of which equals the bare candidate
# id, and one candidate silently becomes two orphan channel rows.
#
# target_keys() alone does not repair that: its _PDB_RE wants a 4-character id
# starting with a digit, so "BersteEVm000003t1_A" normalises to itself. And the
# chain strip cannot be applied blindly, because a candidate id may legitimately
# end in "_<alnum>" ("Berghia_scaffold_12"). It is therefore applied ONLY when
# the stripped form is corroborated against the real candidate id universe --
# resolution by identity, never by pattern-guess.
#
# NOTE: scripts/build_structural_channel.py carries an identical query_keys() /
# canonical_query_id() pair (it was fixed first, before these were hoisted here
# beside the target side). That module already imports _strip_structure_suffixes
# from here; its two local copies should be collapsed into an import of these,
# leaving one definition. Doing so is a change to a file this pass does not own.
# --------------------------------------------------------------------------- #

# Trailing "_<chain>" as foldseek appends it (mmCIF/PDB chain ids are short
# alphanumerics). Only ever consulted against a known candidate id set.
_CHAIN_SUFFIX_RE = re.compile(r"^(.+)_([A-Za-z0-9]{1,4})$")


def query_keys(query: Optional[str]) -> List[str]:
    """Ordered candidate join keys for one Foldseek query id.

    Mirrors target_keys() for the query side: raw -> path basename ->
    extension-stripped basename -> chain-stripped.

    Deduplicated and order-stable. An empty/None query yields [].
    """
    raw = (query or "").strip()
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

    chain = _CHAIN_SUFFIX_RE.match(stem)
    if chain:
        add(chain.group(1))
    return keys


def canonical_query_id(query: Optional[str],
                        candidate_ids: Optional[Set[str]] = None) -> Optional[str]:
    """Resolve a Foldseek query id into the candidate id namespace.

    With `candidate_ids`, returns the first of query_keys() that is a REAL
    candidate id, or None when none of them is -- an unresolvable query is
    reported, never rewritten into something that merely looks joinable.

    Without `candidate_ids` there is nothing to corroborate against, so only
    the unambiguous extension strip is applied and the chain suffix is left
    alone (stripping it on a guess could truncate a legitimate id).
    """
    keys = query_keys(query)
    if not keys:
        return None
    if candidate_ids is None:
        return _strip_structure_suffixes(keys[0].rsplit("/", 1)[-1]) or None
    for key in keys:
        if key in candidate_ids:
            return key
    return None


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
                        tm_threshold: float = 0.5,
                        candidate_ids: Optional[Set[str]] = None) -> dict:
    """Build the per-candidate Foldseek structural-evidence channel.

    Keys are CANDIDATE IDs, not raw Foldseek query fields. Each query is
    normalised through canonical_query_id() -- symmetric with the target
    side's family_for_target()/target_keys() -- so a multi-chain model's
    "<cand_id>_A" / "<cand_id>_B" rows collapse into ONE row for the
    candidate (strongest alntmscore wins) instead of becoming two orphan
    keys that join to nothing.

    `candidate_ids` is the real candidate id universe (e.g. the class-A
    candidate FASTA's headers). Passing it is strongly preferred: it is what
    lets a chain suffix be resolved rather than guessed, and it enables the
    zero-overlap assertion below. Without it the function stays
    backward-compatible (extension strip only, no chain strip, no assertion).

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

    Raises:
        ValueError: if `candidate_ids` is given but empty, or if there are
            hits and NOT ONE of them resolves to a candidate. Zero key
            overlap means the channel would left-join to nothing and the
            structural voter would go dormant while everything exits 0 --
            the exact silent failure this guard exists to make loud.
    """
    if candidate_ids is not None and not candidate_ids:
        raise ValueError(
            "structural_channel: candidate_ids is empty -- the candidate "
            "universe could not be read, so every Foldseek query would fail to "
            "resolve and the structural channel would silently join to nothing."
        )

    hits = parse_foldseek(foldseek_path)

    # Collapse each model's per-chain rows onto its candidate id first, so the
    # classification runs once per CANDIDATE rather than once per chain.
    resolved: Dict[str, dict] = {}
    unresolved: List[str] = []
    for query, best in hits.items():
        canonical = canonical_query_id(query, candidate_ids)
        if canonical is None:
            unresolved.append(query)
            continue
        current = resolved.get(canonical)
        if current is None or best["alntmscore"] > current["alntmscore"]:
            resolved[canonical] = best

    if hits and candidate_ids is not None and not resolved:
        raise ValueError(
            "structural_channel: ZERO of "
            f"{len(hits)} Foldseek queries resolved to a known candidate id. "
            "The structural channel would left-join to nothing and go silently "
            "dormant.\n"
            f"  saw (up to 5):      {sorted(unresolved)[:5]}\n"
            f"  expected ids like:  {sorted(candidate_ids)[:5]}\n"
            "  Likely cause: the Foldseek query dir was staged under the wrong "
            "id scheme (e.g. AF3 job names instead of candidate ids -- bead "
            "5ubd), or the candidate universe is a different candidate set."
        )

    if unresolved:
        print(
            f"[structural_channel] WARNING: {len(unresolved)} of "
            f"{len(hits)} Foldseek queries did not resolve to a candidate id "
            f"and were dropped: {sorted(unresolved)[:5]}"
            f"{' ...' if len(unresolved) > 5 else ''}",
            file=sys.stderr,
        )

    out = {}
    for candidate_id, best in resolved.items():
        state = classify_hit(best, family_map, tm_threshold=tm_threshold)
        out[candidate_id] = {
            "struct_state": state,
            "struct_novelty": 1 if state == "novel" else 0,
            "struct_nonchemo_corrob": 1 if state == "known_non_chemoreceptor" else 0,
            "has_struct_data": True,
        }
    return out
