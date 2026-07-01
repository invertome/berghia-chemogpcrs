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
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from embedding_evidence import centroid, embedding_channel, load_embeddings

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


def load_ref_labels(path: str) -> Dict[str, str]:
    """{accession: family} from a labeled reference TSV (e.g. anchor_set.tsv).

    Only the `accession` and `family` columns are used; any other columns
    (tier, taxid, species, class, evidence, ...) are ignored, so this reads
    anchor_set.tsv's full 7-column schema directly. Missing file -> {} (no
    reference labels available), mirroring embedding_evidence.load_embeddings'
    missing-file contract -- never a crash.
    """
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["accession"], df["family"]))


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


def write_channel_tsv(channel: Dict[str, Dict[str, object]], path: str) -> None:
    """Write `build_embedding_channel`'s output as a TSV keyed by id.

    Columns are exactly `TSV_COLUMNS` (id + the 5 keys
    `embedding_evidence.embedding_channel()` documents) -- this is the exact
    schema `rank_aggregation.merge_evidence_channels`'s `emb_tsv` argument
    expects (see its `value_cols` list). Rows are sorted by id so re-runs
    diff cleanly. An empty `channel` still writes a header-only TSV.
    """
    rows = [{"id": candidate_id, **channel[candidate_id]} for candidate_id in sorted(channel)]
    pd.DataFrame(rows, columns=TSV_COLUMNS).to_csv(path, sep="\t", index=False)


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
    args = parser.parse_args(argv)

    channel = build_embedding_channel(args.candidate_npz, args.ref_npz, args.ref_labels)
    write_channel_tsv(channel, args.out)
    print(
        f"[build_embedding_channel] wrote {len(channel)} candidates -> {args.out}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
