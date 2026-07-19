#!/usr/bin/env python3
"""add_embedding_columns.py — Augment ranked CSV with the embedding channel.

cw3.6 item 3. Left-joins the consensus (or single-model) embedding channel's
per-candidate columns into the ranked candidates CSV so ``emb_novelty`` is an
ALWAYS-PRESENT, sortable descriptive column in the ranked CSV (and therefore in
both ranked views), independent of RANK_METHOD.

Why this is needed: the embedding channel TSV (fusion_consensus.py or
build_embedding_channel.py) is otherwise only consumed by rank_candidates.py's
RANK_METHOD=rankagg voter (_merge_channels_if_present). On the default
RANK_METHOD=weighted path the channel is never merged, so novelty is invisible
to emit_ranked_views.py's discovery score and to reviewers. This augmenter is
the single, method-independent place that surfaces it — mirroring the
add_hcr_columns.py / add_og_coverage_columns.py / add_classification_columns.py
pattern (join in place, never re-score).

Columns added (from the channel's mahalanobis/consensus contract):
    emb_novelty            consensus embedding-novelty score (high = novel);
                           the sortable signal. Blank for candidates absent from
                           the channel, and an empty column when the channel is
                           dormant (no npz) — never a KeyError downstream.
    emb_nonchemo_family    nearest characterized family (descriptive context)
    has_emb_data           whether the candidate had embedding evidence

The channel is AUTHORITATIVE: any pre-existing emb_* columns (e.g. already
merged by the rankagg path) are dropped and re-joined from the channel TSV, so
the output has exactly one of each and no double-counting. emb_novelty is
guaranteed present even when the channel is absent, so downstream never breaks.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

# Channel contract columns we surface (emb_leakage_flag is a blanket True flag,
# not informative as a per-row column, so it is deliberately not carried over).
# emb_novelty_residual (A1, v4bs.2) is the phylogeny-residualized novelty: a
# dormant descriptive column, surfaced when the channel emits it (else absent).
EMB_COLUMNS = ["emb_novelty", "emb_novelty_residual", "emb_nonchemo_family",
               "has_emb_data"]
CHANNEL_ID = "id"


def _split_extra_spec(spec: str):
    """``[PREFIX:]PATH`` -> ``(prefix, path)``.

    A leading token is only treated as a prefix when it contains no ``/``, so
    absolute paths (and paths with colons) are never mis-split.
    """
    head = spec.split(":", 1)[0]
    if ":" in spec and head and "/" not in head:
        prefix, path = spec.split(":", 1)
        return prefix, path
    return "", spec


def add_embedding_columns(
    *,
    ranked_csv_path: str,
    channel_tsv_path: str,
    out_path: str,
    id_column: str = "id",
    extra_tsvs=None,
) -> int:
    """Join the embedding channel's emb_* columns into the ranked CSV in place.

    Returns the number of rows written. Guarantees an ``emb_novelty`` column
    exists in the output regardless of whether the channel TSV is present.

    ``extra_tsvs`` (A3/A4): additional id-keyed diagnostic TSVs to left-join, each
    given as ``[PREFIX:]PATH``. Producers whose column names are already
    self-describing (model_role_split -> emb_role_gap) need no prefix; those with
    generic names (embedding_nulls -> p_value, collapsed) pass a PREFIX so they are
    namespaced in the ranked CSV. A missing TSV is skipped (the producer is
    dormant), never an error — the same contract as the channel itself.
    """
    df = pd.read_csv(ranked_csv_path, dtype=str, keep_default_na=False)

    # Drop any pre-existing emb columns so the channel is authoritative (no dup).
    df = df.drop(columns=[c for c in EMB_COLUMNS if c in df.columns])

    if channel_tsv_path and os.path.exists(channel_tsv_path):
        ch = pd.read_csv(channel_tsv_path, sep="\t", dtype=str,
                         keep_default_na=False)
        if CHANNEL_ID in ch.columns:
            ch = ch.rename(columns={CHANNEL_ID: id_column})
            keep = [id_column] + [c for c in EMB_COLUMNS if c in ch.columns]
            df = df.merge(ch[keep], on=id_column, how="left")
            # Candidates absent from the channel come back as NaN -> blank.
            for c in EMB_COLUMNS:
                if c in df.columns:
                    df[c] = df[c].fillna("")
        else:
            print(f"WARN: embedding channel {channel_tsv_path!r} has no "
                  f"{CHANNEL_ID!r} column; emb columns not joined", file=sys.stderr)

    # A3/A4 diagnostic producers: left-join each extra TSV's columns.
    for spec in (extra_tsvs or []):
        prefix, path = _split_extra_spec(spec)
        if not path or not os.path.exists(path):
            continue                      # producer dormant -> nothing to join
        ex = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        if id_column not in ex.columns and CHANNEL_ID in ex.columns:
            ex = ex.rename(columns={CHANNEL_ID: id_column})
        if id_column not in ex.columns:
            print(f"WARN: extra TSV {path!r} has no {id_column!r} column; skipped",
                  file=sys.stderr)
            continue
        value_cols = [c for c in ex.columns if c != id_column]
        rename = {
            c: (f"{prefix}_{c}" if prefix and not c.startswith(f"{prefix}_") else c)
            for c in value_cols
        }
        ex = ex[[id_column] + value_cols].rename(columns=rename)
        # The producer is authoritative for its own columns (no duplicates).
        df = df.drop(columns=[c for c in rename.values() if c in df.columns])
        df = df.merge(ex, on=id_column, how="left")
        for c in rename.values():
            df[c] = df[c].fillna("")

    # Guarantee the always-present sortable column even when the channel is
    # dormant (no npz) or malformed.
    if "emb_novelty" not in df.columns:
        df["emb_novelty"] = ""

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    return len(df)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--ranked-csv", required=True, help="Input ranked CSV")
    ap.add_argument("--channel-tsv", required=True,
                    help="Embedding channel TSV (fusion_consensus.py / "
                         "build_embedding_channel.py output). Missing = the "
                         "emb_novelty column is added empty.")
    ap.add_argument("--out", required=True, help="Output CSV path")
    ap.add_argument("--id-column", default="id",
                    help="Candidate id column in the ranked CSV (default 'id')")
    ap.add_argument("--extra-tsv", action="append", default=[], metavar="[PREFIX:]PATH",
                    help="Additional id-keyed diagnostic TSV to left-join "
                         "(repeatable). PREFIX namespaces generic column names, "
                         "e.g. 'emb_null:...'. Missing file = dormant, not an error.")
    args = ap.parse_args()

    n_rows = add_embedding_columns(
        ranked_csv_path=args.ranked_csv,
        channel_tsv_path=args.channel_tsv,
        out_path=args.out,
        id_column=args.id_column,
        extra_tsvs=args.extra_tsv,
    )
    print(f"Wrote {args.out} (annotated {n_rows} rows)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
