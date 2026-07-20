#!/usr/bin/env python3
"""embedding_family_assignment.py — is a candidate actually IN a known family?

Bead berghia-chemogpcrs-nsew, deliverable 4. `emb_nonchemo_family` names the
NEAREST non-chemoreceptor family prototype for every candidate unconditionally —
including candidates that are nowhere near any of them. Nearest-family alone
therefore cannot answer "are any Berghia candidates actually bioamine / opsin /
peptide receptors"; on its own it would label all 790 as something.

This script calibrates the assignment against the reference set's OWN spread.
For each family F it computes the distance from every labelled reference member
of F to F's prototypes, giving F's empirical within-family distance
distribution. A candidate's distance to F is then expressed as a percentile of
that distribution:

    pct <= 95   the candidate sits inside the envelope F's own members occupy
                -> a genuine "looks like an F" call
    pct  > 100  the candidate is further from F than EVERY reference F is
                -> nearest-family is an artifact of argmin, not membership

Specificity is reported alongside: the margin (in the same units) between the
best family and the runner-up, so a candidate that is ambiguously between two
families is visible as such.

LENGTH DECONFOUNDING OF THE ENVELOPE IS MANDATORY (fixed 2026-07-20).
---------------------------------------------------------------------------
This script previously scored RAW squared-Mahalanobis distance throughout.
Production does not use a raw novelty: 07_candidate_ranking.sh runs
`fusion_consensus.py --combiner rra --deconfound seq_len`, which residualizes
rank(novelty) on rank(seq_len). The artifact was demonstrably live — the
candidates inside proteinclip3b's envelope but not protrek's were significantly
LONGER (median 469.5 vs 376 aa, p=4.84e-16, rank-biserial +0.640).

But production deconfounds the novelty SCALAR, not the family call: in
`embedding_candidate_diagnostics.candidate_diagnostics`, `nearest_family` is a
raw argmin and `novelty` is the raw min-over-families distance; only that scalar
is then passed through `confound_residuals`. So production semantics are
"raw argmin, deconfounded magnitude", and that is what `--deconfound envelope`
(the default) reproduces:

    family call        -> RAW argmin          (matches `nearest_family`)
    envelope percentile-> DECONFOUNDED residual (matches the deconfounded novelty)

Residualizing every family column instead — the intuitive "deconfound
everything" reading — was tried and is DISQUALIFIED by measurement: reference
self-recovery collapses 98.4%→69.0% (proteinclip3b) and 99.8%→76.0% (protrek),
and cross-model family agreement falls 95.1%→40.6%. See the `--deconfound`
guard and the operator comment for the mechanism and the evidence trail.

Residualization detail. `fusion_consensus.confound_residuals` is imported and
reused verbatim — not reimplemented. Envelope membership compares a candidate's
distance against the REFERENCE members' distances, so both populations must land
on one common scale; the residualization is therefore fit over the POOLED set
(labelled references + candidates) over a single shared id ordering.

KNOWN, DELIBERATE DIVERGENCE FROM PRODUCTION (do not silently "align" it away).
The envelope test is a comparison production never performs, so exact
equivalence is not achievable — production has no deconfounded reference
distances at all. Two differences are therefore structural, not incidental:
  1. ID SET. Production residualizes over the 790 CANDIDATES only, because that
     is all it scores. The envelope test needs references on the same scale, so
     the residualization here is fit over the pooled 838 references + 790
     candidates (n=1628). Different rank basis, hence different residual values.
  2. COLUMN GROUPING. Production residualizes ONE column — the min-over-families
     distance. Here each family's distance column is residualized separately and
     the column of the raw-argmin family is read off. For any given candidate
     the value residualized IS its min-over-families distance, but it is ranked
     against other rows' distances-to-that-same-family rather than against other
     candidates' min distances.
What IS identical to production: the operator itself (`confound_residuals`,
rank-OLS of rank(score) on rank(seq_len)), the geometry (same prototypes, same
tied shrinkage precision, same min-over-prototypes distance), and the family
call (raw argmin).

`raw_argmin_family` is emitted alongside `best_family`; under `envelope` they
are identical by construction and the script hard-fails if they ever diverge.
It exists so downstream integrity checks can cross-compare against production's
raw `nearest_family` without firing on the deconfounding itself.

The same tied precision / multi-prototype construction as
embedding_evidence.mahalanobis_channel is used, so these distances are exactly
the ones the production novelty axis is built from — this is a readout of the
production geometry, not a parallel scoring scheme.

CAVEAT (must survive into any report): embedding proximity is NOT orthology. A
high-confidence call here is one line of evidence that warrants checking with
the dedicated three-source classifier (HMM scan + orthogroup vote +
phylogenetic placement), which has never been run. It is not a verdict.

INTERPRETIVE FRAME: the pipeline is SUBTRACTIVE. A confident EXCLUSION is the
load-bearing output; "unassigned" is a residual, NEVER a positive claim of
chemoreceptor identity.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from build_embedding_channel import (  # noqa: E402
    load_ref_labels,
    novelty_reference_labels,
)
from embedding_evidence import (  # noqa: E402
    family_prototypes,
    load_embeddings,
    mahalanobis_sq,
    shrinkage_precision,
)
from fusion_consensus import confound_residuals  # noqa: E402


def dist_to_family(vec: np.ndarray, protos: np.ndarray, precision: np.ndarray) -> float:
    return min(mahalanobis_sq(vec, p, precision) for p in np.atleast_2d(protos))


def read_fasta_lengths(path: str) -> dict:
    """{first whitespace token of header: residue count}."""
    lengths, cur, n = {}, None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    lengths[cur] = n
                cur, n = line[1:].strip().split()[0], 0
            else:
                n += len(line.strip())
    if cur is not None:
        lengths[cur] = n
    return lengths


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidate-npz", required=True)
    p.add_argument("--ref-npz", required=True)
    p.add_argument("--ref-labels", required=True)
    p.add_argument("--candidate-fasta", required=True,
                   help="class-A candidate FASTA; source of candidate seq_len "
                        "for the mandatory length deconfounding")
    p.add_argument("--anchor-fasta", required=True,
                   help="anchor FASTA; source of reference seq_len, so "
                        "references and candidates share one residual scale")
    p.add_argument("--deconfound", default="envelope",
                   help="how seq_len is residualized out. Only 'envelope' is "
                        "supported (the default): the family argmin stays RAW "
                        "and the envelope percentile is deconfounded. This is "
                        "production semantics. 'all-columns' (residualize "
                        "every family distance) and 'none' are both hard "
                        "failures -- see the guard for the evidence.")
    p.add_argument("--k", type=int, default=3)
    p.add_argument("--pct-threshold", type=float, default=95.0,
                   help="percentile of the family's own member-distance "
                        "distribution below which a candidate counts as inside "
                        "the family envelope (default 95)")
    p.add_argument("--out-tsv", required=True)
    a = p.parse_args()

    # ---- loud guard: reject BOTH failure modes, raw and over-corrected ------
    if a.deconfound == "all-columns":
        print(
            "FATAL: --deconfound 'all-columns' is DISQUALIFIED by measurement, "
            "not by preference.\n"
            "  Residualizing seq_len out of EVERY family-distance column "
            "destroys the family call:\n"
            "    reference self-recovery (does a labelled reference land "
            "nearest its OWN family?)\n"
            "      proteinclip3b  98.4% -> 69.0%   (-29.5 pp)\n"
            "      protrek        99.8% -> 76.0%   (-23.7 pp)\n"
            "    cross-model family agreement      95.1% -> 40.6%\n"
            "    peptide assignments               725/790 -> 19/790, though "
            "peptide is the LARGEST\n"
            "                                      reference family (386/838)\n"
            "  The labelled references stop classifying themselves, which is "
            "disqualifying on its own.\n"
            "  seq_len is strongly family-linked in the reference set "
            "(Kruskal-Wallis H=261.5,\n"
            "  p=1.4e-53), so a blanket residual removes real family signal "
            "along with the artifact.\n"
            "  Do NOT re-enable this to 'be more thorough' -- it is the "
            "over-correction, not the fix.\n"
            "  Evidence trail: results/ranking/diagnostics/"
            "deconfound_impact_audit_*_{raw,dec,hyb}.tsv",
            file=sys.stderr,
        )
        return 2
    if a.deconfound != "envelope":
        print(
            f"FATAL: --deconfound must be 'envelope' (got '{a.deconfound}').\n"
            "  A raw, wholly non-deconfounded assignment is NOT the quantity "
            "production uses:\n"
            "  07_candidate_ranking.sh runs `fusion_consensus.py --combiner "
            "rra --deconfound seq_len`.\n"
            "  The length artifact is real and measured -- envelope-XOR "
            "candidates were longer,\n"
            "  median 469.5 vs 376 aa, p=4.84e-16, rank-biserial +0.640.\n"
            "  To compare variants, use "
            "scripts/unity/embedding_deconfound_impact_audit.py.",
            file=sys.stderr,
        )
        return 2

    cand = load_embeddings(a.candidate_npz)
    ref = load_embeddings(a.ref_npz)
    labels = novelty_reference_labels(load_ref_labels(a.ref_labels))
    labels = {i: f for i, f in labels.items() if i in ref}
    print(f"[assign] candidates={len(cand)} references_resolved={len(labels)}")
    if not cand or not labels:
        print("[assign] FATAL: empty inputs", file=sys.stderr)
        return 1

    # ---- seq_len, asserted against REAL keys (a fixture cannot catch a wrong
    # ---- key; every candidate and every labelled reference must resolve) ----
    cand_len = read_fasta_lengths(a.candidate_fasta)
    ref_len = read_fasta_lengths(a.anchor_fasta)
    miss_c = [i for i in cand if i not in cand_len]
    miss_r = [i for i in labels if i not in ref_len]
    if miss_c or miss_r:
        print(f"[assign] FATAL: seq_len unresolved for {len(miss_c)} candidates "
              f"and {len(miss_r)} references (e.g. {(miss_c + miss_r)[:3]}); "
              "deconfounding cannot proceed on a partial key join",
              file=sys.stderr)
        return 1
    print(f"[assign] seq_len resolved: {len(cand)}/{len(cand)} candidates, "
          f"{len(labels)}/{len(labels)} references")

    protos = family_prototypes(ref, labels, a.k)
    by_family: dict = {}
    for i, f in labels.items():
        by_family.setdefault(f, []).append(np.asarray(ref[i], dtype=float))
    centered = np.vstack([np.array(v) - np.array(v).mean(axis=0)
                          for v in by_family.values()])
    precision = shrinkage_precision(centered)
    fams = sorted(protos)

    # ---- raw distances for references and candidates ----
    ref_ids, cand_ids = list(labels), list(cand)
    d_ref = pd.DataFrame(
        [{f: dist_to_family(np.asarray(ref[i], dtype=float), protos[f], precision)
          for f in fams} for i in ref_ids], index=ref_ids)
    d_cand = pd.DataFrame(
        [{f: dist_to_family(np.asarray(cand[i], dtype=float), protos[f], precision)
          for f in fams} for i in cand_ids], index=cand_ids)

    # ======================================================================
    # WHY THE ARGMIN STAYS RAW AND ONLY THE ENVELOPE IS DECONFOUNDED.
    # Do not "fix" this by residualizing every family column -- that was tried,
    # measured, and is disqualified (see the --deconfound guard above).
    #
    # These are two DIFFERENT comparisons and only one of them is length-biased:
    #
    #   * WHICH family is nearest is a WITHIN-sequence comparison. All seven
    #     family distances for a given candidate are computed from the SAME
    #     sequence, hence the same seq_len, so length enters all seven as a
    #     common factor and very largely CANCELS in the argmin. Residualizing
    #     each column separately destroys that cancellation: it re-ranks the
    #     families against each other using a length trend that had already
    #     cancelled, which is why it collapses reference self-recovery
    #     (98.4% -> 69.0%) and cross-model agreement (95.1% -> 40.6%).
    #
    #   * Whether a candidate falls INSIDE a family's envelope is a BETWEEN-
    #     sequence comparison -- this candidate against reference members of
    #     systematically different lengths. Nothing cancels here, so this is
    #     exactly where the length artifact lives, and this is the comparison
    #     that needs the residual.
    #
    # Hence: argmin from raw distance (matching production's `nearest_family`,
    # which fusion_consensus also leaves raw), envelope percentile from the
    # deconfounded residual (matching production's deconfounded novelty scalar).
    # ======================================================================
    raw_argmin = d_cand[fams].idxmin(axis=1)

    # ---- mandatory deconfounding over the pooled reference+candidate set ----
    pooled_len = {**{i: float(ref_len[i]) for i in ref_ids},
                  **{i: float(cand_len[i]) for i in cand_ids}}
    pooled_d = pd.concat([d_ref[fams], d_cand[fams]])
    if pooled_d.index.duplicated().any():
        print("[assign] FATAL: reference/candidate id collision in the pooled "
              "set; the residual scale would be corrupt", file=sys.stderr)
        return 1
    resid = confound_residuals({f: pooled_d[f].to_dict() for f in fams},
                               {"seq_len": pooled_len})
    if {len(resid[f]) for f in fams} != {len(pooled_d)}:
        print("[assign] FATAL: confound_residuals returned partial coverage",
              file=sys.stderr)
        return 1
    r_all = pd.DataFrame({f: pd.Series(resid[f]) for f in fams})
    r_ref, r_cand = r_all.loc[ref_ids], r_all.loc[cand_ids]
    print(f"[assign] envelope deconfounded on seq_len via "
          f"fusion_consensus.confound_residuals over pooled n={len(pooled_d)} "
          f"({len(ref_ids)} ref + {len(cand_ids)} cand) x {len(fams)} families; "
          "family argmin left RAW (production semantics)")

    # ---- envelopes, from the deconfounded reference residuals ----
    ref_fam = pd.Series(labels)
    envelope = {f: r_ref.loc[ref_fam[ref_fam == f].index, f].to_numpy()
                for f in fams}
    for f in fams:
        d = envelope[f]
        print(f"[assign] {f:<24} n={len(d):<5d} "
              f"median={np.median(d):9.2f} p95={np.percentile(d, 95):9.2f} "
              f"max={d.max():9.2f}   (deconfounded residual units)")

    rows = []
    for cid in cand_ids:
        # family choice: RAW distances (within-sequence comparison)
        draw = {f: float(d_cand.loc[cid, f]) for f in fams}
        order = sorted(fams, key=lambda f: draw[f])
        best, second = order[0], (order[1] if len(order) > 1 else None)
        # envelope test: DECONFOUNDED residual (between-sequence comparison)
        rv = float(r_cand.loc[cid, best])
        env = envelope[best]
        pct = float((env < rv).sum()) / len(env) * 100.0
        rows.append({
            "id": cid,
            "best_family": best,
            "best_dist_raw": draw[best],
            "best_resid_deconfounded": rv,
            "pct_within_family": pct,
            "inside_envelope": pct <= a.pct_threshold,
            "runner_up_family": second,
            "runner_up_dist_raw": draw[second] if second else np.nan,
            "margin_raw": (draw[second] - draw[best]) if second else np.nan,
            "raw_argmin_family": raw_argmin[cid],
            "seq_len": cand_len[cid],
            **{f"d_{f}": draw[f] for f in fams},
            **{f"resid_{f}": float(r_cand.loc[cid, f]) for f in fams},
        })

    df = pd.DataFrame(rows).sort_values(["inside_envelope", "pct_within_family"],
                                        ascending=[False, True])
    df.to_csv(a.out_tsv, sep="\t", index=False)
    print(f"\n[assign] wrote {a.out_tsv} ({len(df)} rows; family call = RAW "
          "argmin, envelope percentile deconfounded on seq_len)")

    n_shift = int((df["best_family"] != df["raw_argmin_family"]).sum())
    if n_shift:
        print(f"[assign] FATAL: best_family diverged from raw_argmin_family on "
              f"{n_shift} candidates; under 'envelope' deconfounding the family "
              "call MUST be the raw argmin", file=sys.stderr)
        return 1
    print(f"[assign] family call == raw argmin for all {len(df)} candidates "
          "(as required); only the envelope percentile is deconfounded")

    inside = df[df["inside_envelope"]]
    print(f"\n=== CANDIDATES INSIDE A KNOWN FAMILY ENVELOPE "
          f"(pct <= {a.pct_threshold}, deconfounded) : {len(inside)}/{len(df)} ===")
    if len(inside):
        print(inside["best_family"].value_counts().to_string())
        print(f"\n{'id':<34} {'family':<22} {'pct':>7} {'dist':>11} {'margin':>11}")
        for _, r in inside.sort_values("pct_within_family").head(40).iterrows():
            print(f"{r['id']:<34} {r['best_family']:<22} "
                  f"{r['pct_within_family']:7.2f} {r['best_dist_raw']:11.2f} "
                  f"{r['margin_raw']:11.2f}")
    else:
        print("(none — no candidate falls inside any family's own member spread)")

    print(f"\n=== nearest-family assignment for ALL {len(df)} candidates "
          "(RAW argmin, uncalibrated) ===")
    print(df["best_family"].value_counts().to_string())
    print("\npct_within_family distribution (100 = further than every real member):")
    s = df["pct_within_family"]
    for q in (0, 1, 5, 10, 25, 50, 75, 90, 100):
        print(f"  p{q}={np.percentile(s, q):.2f}")
    print(f"  == 100.00 (beyond every reference member): "
          f"{int((s >= 100.0).sum())}/{len(s)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
