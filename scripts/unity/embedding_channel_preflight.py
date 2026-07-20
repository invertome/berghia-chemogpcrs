#!/usr/bin/env python3
"""embedding_channel_preflight.py — hard-fail guards for the embedding channel rebuild.

Bead berghia-chemogpcrs-nsew. The embedding/novelty channel has three failure
modes that all exit 0 and emit a plausible-looking artifact. Each one silently
voids the run. This script asserts against all three BEFORE any scoring happens,
and exits non-zero (loudly) if any fires.

TRAP 1 — reference label overlap (the defect this rebuild fixes).
    load_ref_labels() rebuilds the composite key ANCHOR_<class>_<tier>_<acc>
    from the label TSV; the reference npz is keyed by the anchor FASTA header,
    which is that SAME composite format. So a wrong label file still produces
    format-valid keys and the join SUCCEEDS while resolving almost nothing.
    Measured with the old default (references/anchors/anchor_set.tsv, 206 rows):
    116/953 = 12%. With anchor_set_PROD.tsv: 953/953. There is no exception and
    no NaN — 837 reference vectors were just dropped.

TRAP 2 — deconfound key mismatch.
    Stage 07 runs fusion_consensus.py with --deconfound seq_len. seq_len is
    keyed on the candidate FASTA header's FIRST TOKEN (_read_lengths); novelty
    is keyed on the candidate .npz keys. confound_residuals keeps only ids
    present in EVERY confound (fusion_consensus.py:91), so zero overlap yields
    an empty id list -> empty residuals -> an RRA over nothing -> a header-only
    channel TSV -> has_emb_data=False for every candidate. Exit 0, one
    informational log line.

TRAP 3 — zero-vector prototypes.
    embedding_evidence.centroid() returns np.zeros(dim) for a family whose ids
    resolve to no embedding, and _cosine_similarity() converts a zero-norm
    vector into a plausible 0.0 rather than NaN. A family with no resolved
    members is therefore INVISIBLE downstream. This is not assertable as a
    single threshold (a genuinely small family is legitimate), so instead we
    REPORT per-family resolved member counts and hard-fail only on the
    degenerate case: a family present in the labels that resolves to zero
    embeddings.

Exit codes: 0 = all guards passed; 1 = a guard fired (details on stderr).
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from build_embedding_channel import (  # noqa: E402
    load_ref_labels,
    novelty_reference_labels,
)
from embedding_evidence import load_embeddings  # noqa: E402
from embedding_candidate_diagnostics import _read_lengths  # noqa: E402

REF_OVERLAP_MIN = 0.80
SEQLEN_OVERLAP_MIN = 0.95


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ref-labels", required=True)
    p.add_argument("--candidate-fasta", required=True)
    p.add_argument(
        "--models", nargs="+", required=True,
        help="tag:cand_npz:ref_npz specs (identity_tsv suffix tolerated/ignored)",
    )
    p.add_argument("--out-json", help="write the measured numbers here")
    a = p.parse_args()

    failures: list[str] = []
    report: dict = {"models": {}}

    ref_labels = load_ref_labels(a.ref_labels)
    novelty_labels = novelty_reference_labels(ref_labels)
    report["ref_labels_file"] = a.ref_labels
    report["ref_labels_rows"] = len(ref_labels)
    print(f"[preflight] ref labels: {a.ref_labels} -> {len(ref_labels)} composite keys "
          f"({len(novelty_labels)} after the characterized-class-A novelty filter)")

    if not ref_labels:
        failures.append(f"TRAP 1: --ref-labels {a.ref_labels} yielded ZERO labels "
                        "(missing file or unreadable schema)")

    seq_len = _read_lengths(a.candidate_fasta)
    report["seq_len_keys"] = len(seq_len)
    print(f"[preflight] seq_len: {a.candidate_fasta} -> {len(seq_len)} keys")

    for spec in a.models:
        parts = spec.split(":")
        tag, cand_npz, ref_npz = parts[0], parts[1], parts[2]
        cand = load_embeddings(cand_npz)
        ref = load_embeddings(ref_npz)
        m: dict = {"candidates": len(cand), "reference": len(ref)}
        print(f"\n[preflight] === model {tag} ===")
        print(f"[preflight] candidate npz: {len(cand)} keys | reference npz: {len(ref)} keys")

        if not cand or not ref:
            failures.append(f"TRAP: model {tag} has an EMPTY npz "
                            f"(cand={len(cand)}, ref={len(ref)})")
            report["models"][tag] = m
            continue

        # --- TRAP 1: reference label overlap ---------------------------------
        ov = set(ref_labels) & set(ref)
        frac = len(ov) / len(ref)
        m["ref_overlap"] = len(ov)
        m["ref_overlap_frac"] = round(frac, 4)
        print(f"[preflight] TRAP1 ref-label overlap: {len(ov)}/{len(ref)} = {frac:.4f} "
              f"(need > {REF_OVERLAP_MIN})")
        if frac <= REF_OVERLAP_MIN:
            failures.append(
                f"TRAP 1 FIRED ({tag}): reference-label overlap {len(ov)}/{len(ref)} "
                f"= {frac:.4f} <= {REF_OVERLAP_MIN}. The label TSV does not describe "
                f"the reference npz; family prototypes would be built from a "
                f"fraction of the references and the join would still 'succeed'."
            )

        # --- TRAP 2: deconfound (seq_len) key overlap ------------------------
        ov_len = set(seq_len) & set(cand)
        frac_len = len(ov_len) / len(cand)
        m["seqlen_overlap"] = len(ov_len)
        m["seqlen_overlap_frac"] = round(frac_len, 4)
        print(f"[preflight] TRAP2 seq_len/candidate overlap: {len(ov_len)}/{len(cand)} "
              f"= {frac_len:.4f} (need > {SEQLEN_OVERLAP_MIN})")
        if frac_len <= SEQLEN_OVERLAP_MIN:
            failures.append(
                f"TRAP 2 FIRED ({tag}): seq_len/candidate id overlap {len(ov_len)}/"
                f"{len(cand)} = {frac_len:.4f} <= {SEQLEN_OVERLAP_MIN}. "
                f"confound_residuals would keep only the intersection, so "
                f"--deconfound seq_len would residualize over "
                f"{len(ov_len)} candidates and the channel would come back "
                f"empty/partial with exit 0."
            )

        # --- TRAP 3: per-family resolved prototype members -------------------
        resolved = collections.Counter()
        labelled = collections.Counter()
        for rid, fam in novelty_labels.items():
            labelled[fam] += 1
            if rid in ref:
                resolved[fam] += 1
        m["family_resolved"] = dict(resolved)
        m["family_labelled"] = dict(labelled)
        print(f"[preflight] TRAP3 per-family resolved prototype members "
              f"(characterized class-A only):")
        print(f"[preflight]   {'family':<26} {'labelled':>9} {'resolved':>9}")
        for fam in sorted(labelled):
            flag = "  <-- DEGENERATE (zero-vector prototype)" if resolved[fam] == 0 else ""
            print(f"[preflight]   {fam:<26} {labelled[fam]:>9} {resolved[fam]:>9}{flag}")
            if resolved[fam] == 0:
                failures.append(
                    f"TRAP 3 FIRED ({tag}): family '{fam}' has {labelled[fam]} labelled "
                    f"references but ZERO resolve to an embedding. centroid() returns "
                    f"a zero vector and _cosine_similarity turns that into 0.0, not "
                    f"NaN — the family is silently invisible and any candidate that "
                    f"IS one scores maximally novel."
                )
        report["models"][tag] = m

    print()
    if failures:
        print("=" * 72, file=sys.stderr)
        print(f"PREFLIGHT FAILED — {len(failures)} guard(s) fired:", file=sys.stderr)
        for f in failures:
            print(f"  * {f}", file=sys.stderr)
        print("=" * 72, file=sys.stderr)
        rc = 1
    else:
        print("[preflight] ALL GUARDS PASSED — proceeding to score.")
        rc = 0

    report["failures"] = failures
    if a.out_json:
        with open(a.out_json, "w") as fh:
            json.dump(report, fh, indent=2, sort_keys=True)
        print(f"[preflight] wrote {a.out_json}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
