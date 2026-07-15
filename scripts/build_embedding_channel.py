#!/usr/bin/env python3
"""build_embedding_channel.py — ESM-C embedding-channel PRODUCER (Glue task G2).

Wires the pure scorer (embedding_evidence.py: load_embeddings / centroid /
nearest_family / embedding_channel) to real inputs -- a precomputed candidate
embeddings .npz, a precomputed reference embeddings .npz, and a family/class-
labeled reference TSV (references/anchors/anchor_set.tsv, produced by
build_anchor_set.py / curate_gpcr_references.py) -- and writes the resulting
per-candidate channel out as a TSV that rank_aggregation.merge_evidence_channels
left-joins onto the ranking dataframe by id.

HARD GUARDRAIL (inherited from embedding_evidence.py's module contract): this
producer must NEVER build a chemoreceptor centroid or emit a positive
"similarity to a known chemoreceptor" key. The reference set consumed here is
the curated NON-chemoreceptor GPCR set by construction (every family in
anchor_set.tsv -- aminergic, class-A-other, class-B-secretin, class-C,
class-F-frizzled, lipid, opsin, peptide -- is a non-chemoreceptor family), so
every family present is eligible for `nonchemo_centroids`. As a defensive
belt-and-braces measure against future reference-set contamination,
`build_family_centroids` still drops any family whose label literally
contains "chemoreceptor" before computing anything. The emitted TSV's column
set is hardcoded to exactly the 5 keys `embedding_channel()` documents, so no
`*chemoreceptor_sim*` key can ever appear in this producer's output either
(see tests/unit/test_build_embedding_channel.py's guard test).

`classA_centroid` ("broad class-A-GPCR-ness", used only for RECALL) is
derived from each ref's family label via `family_to_class`, which reproduces
build_anchor_set.py's documented rule verbatim: "family -> class rule:
class-B*->B, class-C*->C, class-F*->F, else A". This keeps `ref_labels` a
plain {ref_id: family} mapping (no separate class column needed downstream of
`load_ref_labels`) while staying exactly consistent with how
build_anchor_set.py itself assigns GPCR classes.

Reference-embedding production: scripts/unity/run_esmc_reference_embeddings.sh
(Unity GPU; mirrors scripts/unity/run_esmc_embeddings.sh's corrected ESM-C
API). Candidate embeddings: scripts/unity/run_esmc_embeddings.sh. Both MUST
come from the same ESM-C checkpoint (e.g. esmc_300m) -- comparing embeddings
from different model sizes/checkpoints is meaningless.

Usage:
    python3 build_embedding_channel.py \\
        --candidate-npz results/ranking/embeddings/candidates_esmc300m.npz \\
        --ref-npz results/ranking/embeddings/reference_esmc300m.npz \\
        --ref-labels references/anchors/anchor_set.tsv \\
        --out results/ranking/embeddings/emb_channel.tsv
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, FrozenSet, List, Optional, Tuple

import numpy as np
import pandas as pd

from embedding_evidence import (
    centroid,
    embedding_channel,
    load_embeddings,
    mahalanobis_channel,
)

# Columns of the output TSV, in order. MUST stay exactly these 5 value keys
# (plus id) -- rank_aggregation.merge_evidence_channels reads
# emb_classA_sim / emb_nonchemo_sim / emb_nonchemo_family / emb_leakage_flag
# by these exact names, and embedding_evidence.embedding_channel() documents
# has_emb_data as the 5th per-candidate key.
TSV_COLUMNS: List[str] = [
    "id",
    "emb_nonchemo_sim",
    "emb_nonchemo_family",
    "emb_classA_sim",
    "has_emb_data",
    "emb_leakage_flag",
]

# Columns of the Mahalanobis (`maha`) channel TSV, in order. Emits emb_novelty
# (S_novel, the positive ranking axis) instead of the cosine channel's
# emb_nonchemo_sim/emb_classA_sim; rank_aggregation reads emb_novelty as a
# positive voter and emb_nonchemo_family as annotation.
MAHA_TSV_COLUMNS: List[str] = [
    "id",
    "emb_nonchemo_family",
    "emb_novelty",
    "has_emb_data",
    "emb_leakage_flag",
]

# family -> class rule, copied verbatim from build_anchor_set.py's docstring
# so classA membership is derived exactly the way the reference set itself
# was assembled, not re-guessed here.
_NON_A_CLASS_PREFIXES: Tuple[Tuple[str, str], ...] = (
    ("class-B", "B"),
    ("class-C", "C"),
    ("class-F", "F"),
)


def family_to_class(family: str) -> str:
    """GPCR class letter for a family label.

    Mirrors build_anchor_set.py's documented rule: "family -> class rule:
    class-B*->B, class-C*->C, class-F*->F, else A".
    """
    for prefix, gpcr_class in _NON_A_CLASS_PREFIXES:
        if family.startswith(prefix):
            return gpcr_class
    return "A"


# Class-A families whose members have NO known ligand/function. An orphan is not
# a *known non-chemoreceptor* -- it is an unknown, and some orphans may in fact be
# chemoreceptors (the reference set's own orphan bucket held a C. elegans str-33
# and an olfactory receptor before curation). They therefore must not define
# "known" space for novelty; they remain valid class-A anchors for the broad
# class-A recall centroid only.
UNCHARACTERIZED_FAMILIES: FrozenSet[str] = frozenset({"orphan"})


def novelty_reference_labels(ref_labels: Dict[str, str]) -> Dict[str, str]:
    """`ref_labels` restricted to CHARACTERIZED class-A families.

    S_novel is a MIN over family prototypes, so every extra prototype can only
    pull a candidate's novelty DOWN. Two kinds of reference must therefore be
    kept out of the novelty prototypes (and out of the tied covariance and
    background mean, which `mahalanobis_channel` derives from these same labels):

      * out-of-class families (class-B/C/F) -- they exist to support the upstream
        class-A gate, not to explain a class-A candidate away. Without this a
        divergent class-A candidate that happens to land nearest a class-C
        centroid is scored "known", which masks novelty exactly on the most
        divergent LSE candidates -- the ones we most want surfaced.
      * `UNCHARACTERIZED_FAMILIES` (orphan) -- see that constant.

    Novelty is thus "unlike any CHARACTERIZED non-chemoreceptor family".
    """
    return {
        ref_id: family
        for ref_id, family in ref_labels.items()
        if family_to_class(family) == "A" and family not in UNCHARACTERIZED_FAMILIES
    }


def load_ref_labels(path: str) -> Dict[str, str]:
    """{composite_id: family} from a labeled reference TSV (e.g. anchor_set.tsv).

    The reference embeddings .npz (scripts/unity/run_esmc_reference_embeddings.sh)
    is keyed by each sequence's FASTA header, which build_anchor_set.py writes
    as the COMPOSITE id ANCHOR_<class>_<tier>_<accession> (anchor_header(),
    build_anchor_set.py:128) -- e.g. 'ANCHOR_A_1_P31356'. anchor_set.tsv's
    `accession` column, by contrast, is the BARE accession ('P31356'). So this
    reconstructs the composite id per row -- mirroring build_anchor_set.py's
    anchor_header f-string VERBATIM -- so the keys match the .npz and
    build_family_centroids' join actually resolves. (Verified against source:
    bare accession vs the FASTA composite header overlaps 0/206; the
    reconstructed composite overlaps 206/206.)

    Uses the `accession`, `tier`, and `class` columns; any other columns
    (taxid, species, evidence, ...) are ignored, so this reads anchor_set.tsv's
    full schema directly. Missing file -> {} (no reference labels available),
    mirroring embedding_evidence.load_embeddings' missing-file contract --
    never a crash.
    """
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t")
    return {
        # mirrors build_anchor_set.anchor_header: f"ANCHOR_{klass}_{tier}_{accession}"
        f"ANCHOR_{cls}_{tier}_{acc}": family
        for acc, tier, cls, family in zip(
            df["accession"], df["tier"], df["class"], df["family"]
        )
    }


def build_family_centroids(
    ref_embeddings: Dict[str, np.ndarray], ref_labels: Dict[str, str]
) -> Tuple[Dict[str, np.ndarray], np.ndarray]:
    """(nonchemo_centroids, classA_centroid) from a labeled reference set.

    `ref_labels` is {ref_id: family} (e.g. from `load_ref_labels`). Every
    family is treated as non-chemoreceptor -- see the module HARD GUARDRAIL
    docstring above -- except any family whose label literally contains
    "chemoreceptor", which is dropped from BOTH outputs as a defensive
    contamination guard (never expected to fire against the real curated
    reference set).

    `nonchemo_centroids`: {family: centroid vector} for every distinct
    (non-dropped) family in `ref_labels`, via `embedding_evidence.centroid`.
    `classA_centroid`: centroid over every (non-dropped) ref id whose family
    maps to class "A" under `family_to_class`.

    ids in `ref_labels` absent from `ref_embeddings` are silently ignored
    (delegated to `centroid`'s existing contract); a family/class bucket
    with no resolvable embeddings gets a zero-vector centroid rather than
    raising.
    """
    clean_labels: Dict[str, str] = {}
    for ref_id, family in ref_labels.items():
        # Guard a malformed custom --ref-labels: a missing family cell reads as
        # NaN (a float), and NaN is truthy in Python so `family or ""` would NOT
        # catch it -- coerce any non-str to "" before .lower() to avoid an
        # AttributeError (never expected to fire against the real anchor_set.tsv).
        if not isinstance(family, str):
            family = ""
        if "chemoreceptor" in family.lower():
            print(
                f"[build_embedding_channel] WARNING: dropping ref {ref_id!r} "
                f"labeled family={family!r} -- contains 'chemoreceptor', which "
                "this producer must never centroid (see embedding_evidence.py's "
                "HARD RULE)",
                file=sys.stderr,
            )
            continue
        clean_labels[ref_id] = family

    families: Dict[str, List[str]] = {}
    for ref_id, family in clean_labels.items():
        families.setdefault(family, []).append(ref_id)
    nonchemo_centroids = {
        family: centroid(ref_embeddings, ids) for family, ids in families.items()
    }

    classA_ids = [
        ref_id for ref_id, family in clean_labels.items()
        if family_to_class(family) == "A"
    ]
    classA_centroid = centroid(ref_embeddings, classA_ids)

    return nonchemo_centroids, classA_centroid


def build_embedding_channel(
    candidate_npz: str, ref_npz: str, ref_labels_tsv: str
) -> Dict[str, Dict[str, object]]:
    """embedding_evidence.embedding_channel() for every candidate in
    `candidate_npz`, using family centroids built from `ref_npz` +
    `ref_labels_tsv`.

    Thin glue, delegating all scoring logic to embedding_evidence.py (this
    producer never computes a cosine similarity itself). Missing/empty
    inputs degrade gracefully to an empty channel (`{}`), never a crash --
    `load_embeddings` and `load_ref_labels` both already treat a missing
    path as "no data available".
    """
    candidate_embeddings = load_embeddings(candidate_npz)
    ref_embeddings = load_embeddings(ref_npz)
    ref_labels = load_ref_labels(ref_labels_tsv)
    nonchemo_centroids, classA_centroid = build_family_centroids(
        ref_embeddings, ref_labels
    )
    return embedding_channel(candidate_embeddings, nonchemo_centroids, classA_centroid)


def build_embedding_channel_maha(
    candidate_npz: str,
    ref_npz: str,
    ref_labels_tsv: str,
    k: int = 3,
    relative: bool = False,
) -> Dict[str, Dict[str, object]]:
    """Mahalanobis multi-prototype embedding channel (the Phase-0 scorer path).

    Thin glue over `embedding_evidence.mahalanobis_channel`: loads the candidate +
    reference embeddings + the family-labeled anchor TSV, then emits per candidate
    `emb_nonchemo_family` (annotation), `emb_novelty` (S_novel — the POSITIVE
    ranking axis), `has_emb_data`, `emb_leakage_flag`. Preferred with the ESM-C
    600M embeddings (the bake-off winner). Missing/empty inputs degrade to `{}`,
    never a crash (mirrors `build_embedding_channel`).
    """
    candidate_embeddings = load_embeddings(candidate_npz)
    ref_embeddings = load_embeddings(ref_npz)
    ref_labels = novelty_reference_labels(load_ref_labels(ref_labels_tsv))
    if not candidate_embeddings or not ref_embeddings or not ref_labels:
        return {}
    return mahalanobis_channel(
        candidate_embeddings, ref_embeddings, ref_labels, k=k, relative=relative
    )


def write_channel_tsv(
    channel: Dict[str, Dict[str, object]],
    path: str,
    columns: List[str] = TSV_COLUMNS,
) -> None:
    """Write a channel dict as a TSV keyed by id, with the given column order.

    `columns` defaults to the cosine `TSV_COLUMNS` (id + the 5 keys
    `embedding_evidence.embedding_channel()` documents); pass `MAHA_TSV_COLUMNS`
    for the Mahalanobis channel. Rows are sorted by id so re-runs diff cleanly.
    An empty `channel` still writes a header-only TSV.
    """
    rows = [{"id": candidate_id, **channel[candidate_id]} for candidate_id in sorted(channel)]
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)


def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(
        description="Build the ESM-C embedding-evidence channel TSV for candidate ranking."
    )
    parser.add_argument(
        "--candidate-npz", required=True,
        help="Candidate embeddings .npz (scripts/unity/run_esmc_embeddings.sh output)",
    )
    parser.add_argument(
        "--ref-npz", required=True,
        help="Reference embeddings .npz (scripts/unity/run_esmc_reference_embeddings.sh output)",
    )
    parser.add_argument(
        "--ref-labels", required=True,
        help="Family/class-labeled reference TSV (e.g. references/anchors/anchor_set.tsv)",
    )
    parser.add_argument("--out", required=True, help="Output channel TSV path")
    parser.add_argument(
        "--scorer", choices=["cosine", "maha"], default="cosine",
        help="cosine (default, unchanged) or maha (tied-covariance Mahalanobis + "
             "multi-prototype centroids; use with the ESM-C 600M embeddings)",
    )
    parser.add_argument(
        "--relative", action="store_true",
        help="maha only: relative Mahalanobis for emb_novelty (conservative "
             "variant); default is plain Mahalanobis",
    )
    parser.add_argument(
        "--k", type=int, default=3,
        help="maha only: k-means sub-prototypes per family (default 3)",
    )
    args = parser.parse_args(argv)

    if args.scorer == "maha":
        channel = build_embedding_channel_maha(
            args.candidate_npz, args.ref_npz, args.ref_labels,
            k=args.k, relative=args.relative,
        )
        write_channel_tsv(channel, args.out, columns=MAHA_TSV_COLUMNS)
    else:
        channel = build_embedding_channel(args.candidate_npz, args.ref_npz, args.ref_labels)
        write_channel_tsv(channel, args.out)
    print(
        f"[build_embedding_channel] scorer={args.scorer} wrote {len(channel)} "
        f"candidates -> {args.out}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
