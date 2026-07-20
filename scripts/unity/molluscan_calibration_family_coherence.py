#!/usr/bin/env python3
"""molluscan_calibration_family_coherence.py — what would happen to the
exclusion set if the phylogenetically incoherent reference families were fixed.

MEASURE ONLY. This script changes NOTHING in production. Which families exist
and where their boundaries fall is a design decision for the user; this
quantifies the consequences so that decision is made on numbers.

THE CENSUS FINDING BEING TESTED
-------------------------------
Three of the seven characterized class-A families the novelty channel uses are
not phylogenetically coherent for a molluscan query:

  chemokine     absent from all protostomes (Mollusca 0). Its envelope is
                vacuous -- there is no such thing as a molluscan chemokine
                receptor, so every candidate assigned to it is an argmin
                artifact.
  lipid         merges prostanoid (Mollusca 371) with cannabinoid / S1P / LPA /
                leukotriene (all 0). One coherent protostome group is averaged
                with four groups that do not exist outside deuterostomes.
  nucleotide    merges adenosine (Mollusca 41) with P2Y (0).

(glycoprotein-hormone merges the pre-bilaterian LGR group, Mollusca 157, with
the vertebrate TSHR/FSHR/LHCGR radiation; it is measured here as a fourth
variant, reported separately from the three the census named.)

THE MECHANISM, STATED UP FRONT SO THE NUMBERS ARE NOT OVER-READ
---------------------------------------------------------------
Novelty is a MIN over family prototypes, so the two edits push in opposite
directions and neither is a free improvement:

  * DROPPING a family removes prototypes, so every candidate's min-distance can
    only stay the same or INCREASE. But the candidates whose argmin WAS the
    dropped family are then re-assigned to their next-nearest family and tested
    against that family's envelope instead, which they may well fall inside. So
    even dropping is not monotone on the exclusion set.
  * SPLITTING a family adds prototypes (each subfamily gets its own k
    prototypes) and tightens each envelope onto a smaller, more homogeneous
    member set. Min-distance can only stay the same or DECREASE, but the
    envelope each candidate is tested against also narrows.

  In both directions the net effect is therefore NOT predictable from the
  mechanism and has to be measured -- which is the point of this script. The
  per-variant `gained`/`lost` counts are reported precisely because the net
  delta hides movement in both directions.

The tied shrinkage precision is re-estimated per variant, because it is derived
from the family-centered reference matrix and so is itself a function of the
family definition. Holding it fixed would measure a different, incoherent thing.

Everything else is the production geometry verbatim: k=3 prototypes,
min-over-prototypes squared Mahalanobis, RAW argmin family call, and
`fusion_consensus.confound_residuals` length deconfounding of the envelope
comparison over the pooled reference + candidate (+ optional null) set.
"""
from __future__ import annotations

import argparse
import json
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


def dist_to_family(vec, protos, precision) -> float:
    return min(mahalanobis_sq(vec, p, precision) for p in np.atleast_2d(protos))


def read_fasta_lengths(path: str) -> dict:
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


def build_variants(base_labels: dict, sub: dict) -> dict:
    """{variant name: {ref_id: family}}. `sub` maps ref_id -> subfamily."""
    def split(fam):
        out = {}
        for i, f in base_labels.items():
            if f == fam:
                s = sub.get(i)
                if s is None:
                    raise SystemExit(f"FATAL: no subfamily label for {i} ({fam})")
                out[i] = f"{fam}:{s}"
            else:
                out[i] = f
        return out

    def drop(labels, fam):
        return {i: f for i, f in labels.items() if f != fam}

    def lipid_binary(labels):
        # the census boundary: prostanoid (present in Mollusca) vs everything
        # else in `lipid` (absent from protostomes)
        out = {}
        for i, f in labels.items():
            if f.startswith("lipid:"):
                out[i] = ("lipid:prostanoid" if f == "lipid:prostanoid"
                          else "lipid:non-prostanoid")
            else:
                out[i] = f
        return out

    v = {}
    v["baseline"] = dict(base_labels)
    v["drop_chemokine"] = drop(base_labels, "chemokine")
    v["split_lipid"] = lipid_binary(split("lipid"))
    v["split_lipid_fine"] = split("lipid")
    v["split_nucleotide"] = split("nucleotide")
    v["split_glycoprotein_hormone"] = split("glycoprotein-hormone")

    comb = lipid_binary(split("lipid"))
    comb2 = {}
    for i, f in comb.items():
        if base_labels[i] == "nucleotide":
            comb2[i] = f"nucleotide:{sub[i]}"
        else:
            comb2[i] = f
    v["census_combined"] = drop(comb2, "chemokine")
    return v


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--model", required=True)
    p.add_argument("--ref-npz", required=True)
    p.add_argument("--candidate-npz", required=True)
    p.add_argument("--ref-labels", required=True)
    p.add_argument("--subfamily-tsv", required=True,
                   help="accession/family/subfamily, resolved from UniProt "
                        "protein names + gene names + curation provenance")
    p.add_argument("--candidate-fasta", required=True)
    p.add_argument("--anchor-fasta", required=True)
    p.add_argument("--null-npz", default="",
                   help="optional molluscan null; when given, each variant is "
                        "ALSO scored against the molluscan background")
    p.add_argument("--null-fasta", default="")
    p.add_argument("--k", type=int, default=3)
    p.add_argument("--pct-threshold", type=float, default=95.0)
    p.add_argument("--moll-threshold", type=float, default=5.0)
    p.add_argument("--out-prefix", required=True)
    a = p.parse_args()

    ref = load_embeddings(a.ref_npz)
    cand = load_embeddings(a.candidate_npz)
    base = novelty_reference_labels(load_ref_labels(a.ref_labels))
    base = {i: f for i, f in base.items() if i in ref}
    null = load_embeddings(a.null_npz) if a.null_npz else {}

    cand_len = read_fasta_lengths(a.candidate_fasta)
    ref_len = read_fasta_lengths(a.anchor_fasta)
    null_len = read_fasta_lengths(a.null_fasta) if a.null_fasta else {}

    # ---- subfamily labels, keyed the way the npz is keyed ----------------
    sdf = pd.read_csv(a.subfamily_tsv, sep="\t")
    sub = {f"ANCHOR_{c}_{t}_{acc}": s
           for acc, t, c, s in zip(sdf["accession"], sdf["tier"],
                                   sdf["class"], sdf["subfamily"])}
    # HARD key-overlap assertion on REAL data: a fixture cannot catch a wrong key
    split_fams = {"lipid", "nucleotide", "glycoprotein-hormone"}
    need = {i for i, f in base.items() if f in split_fams}
    resolved = need & set(sub)
    print(f"[coh] subfamily key overlap: {len(resolved)}/{len(need)} references "
          f"in the split-eligible families")
    if len(resolved) != len(need):
        print(f"[coh] FATAL: {len(need - resolved)} unresolved "
              f"(e.g. {sorted(need - resolved)[:3]})", file=sys.stderr)
        return 1

    variants = build_variants(base, sub)
    results, per_variant_calls = {}, {}

    for name, labels in variants.items():
        protos = family_prototypes(ref, labels, a.k)
        by_family: dict = {}
        for i, f in labels.items():
            by_family.setdefault(f, []).append(np.asarray(ref[i], dtype=float))
        # a k=3 prototype set on a 1-2 member family is degenerate; report it
        tiny = {f: len(v) for f, v in by_family.items() if len(v) < a.k}
        centered = np.vstack([np.array(v) - np.array(v).mean(axis=0)
                              for v in by_family.values()])
        precision = shrinkage_precision(centered)
        fams = sorted(protos)

        ref_ids, cand_ids = list(labels), list(cand)

        def dmat(store, ids):
            return pd.DataFrame(
                [{f: dist_to_family(np.asarray(store[i], dtype=float),
                                    protos[f], precision) for f in fams}
                 for i in ids], index=ids)

        d_ref, d_cand = dmat(ref, ref_ids), dmat(cand, cand_ids)
        frames = [d_ref[fams], d_cand[fams]]
        lens = {**{i: float(ref_len[i]) for i in ref_ids},
                **{i: float(cand_len[i]) for i in cand_ids}}
        null_ids = list(null)
        if null:
            d_null = dmat(null, null_ids)
            frames.append(d_null[fams])
            lens.update({i: float(null_len[i]) for i in null_ids})

        pooled = pd.concat(frames)
        r = confound_residuals({f: pooled[f].to_dict() for f in fams},
                               {"seq_len": lens})
        R = pd.DataFrame({f: pd.Series(r[f]) for f in fams})

        ref_fam = pd.Series(labels)
        env = {f: R.loc[ref_fam[ref_fam == f].index, f].to_numpy() for f in fams}
        argmin = d_cand[fams].idxmin(axis=1)

        calls = []
        for cid in cand_ids:
            best = argmin[cid]
            v = float(R.loc[cid, best])
            pv = float((env[best] < v).sum()) / len(env[best]) * 100.0
            rec = {"id": cid, "best_family": best, "pct_vert": pv,
                   "inside_vert": pv <= a.pct_threshold}
            if null:
                pop = R.loc[null_ids, best].to_numpy()
                pm = float((pop < v).sum()) / len(pop) * 100.0
                rec["pct_moll"] = pm
                rec["inside_moll"] = pm <= a.moll_threshold
            calls.append(rec)
        cdf = pd.DataFrame(calls).set_index("id")
        per_variant_calls[name] = cdf

        res = {"n_families": len(fams),
               "families": {f: int((ref_fam == f).sum()) for f in fams},
               "families_below_k": tiny,
               "inside_vert": int(cdf["inside_vert"].sum()),
               "inside_vert_by_family":
                   cdf.loc[cdf["inside_vert"], "best_family"].value_counts().to_dict(),
               "argmin_by_family": cdf["best_family"].value_counts().to_dict()}
        if null:
            res["inside_moll"] = int(cdf["inside_moll"].sum())
            res["inside_moll_by_family"] = (
                cdf.loc[cdf["inside_moll"], "best_family"].value_counts().to_dict())
        results[name] = res
        print(f"[coh] {name:<28} families={len(fams):<3} "
              f"inside_vert={res['inside_vert']:<4}"
              + (f" inside_moll={res['inside_moll']}" if null else "")
              + (f"  degenerate(<k members): {tiny}" if tiny else ""))

    # ---- deltas against baseline -----------------------------------------
    b = per_variant_calls["baseline"]
    print("\n" + "=" * 78)
    print(f"D4. EFFECT ON THE EXCLUSION SET  [{a.model}]")
    print("=" * 78)
    print(f"{'variant':<28} {'n_fam':>6} {'inside_vert':>12} {'delta':>7} "
          f"{'gained':>7} {'lost':>6}")
    for name, cdf in per_variant_calls.items():
        gained = int((cdf["inside_vert"] & ~b["inside_vert"]).sum())
        lost = int((~cdf["inside_vert"] & b["inside_vert"]).sum())
        results[name]["gained_vs_baseline"] = gained
        results[name]["lost_vs_baseline"] = lost
        print(f"{name:<28} {results[name]['n_families']:>6} "
              f"{results[name]['inside_vert']:>12} "
              f"{results[name]['inside_vert'] - results['baseline']['inside_vert']:>+7} "
              f"{gained:>7} {lost:>6}")
    if null:
        print(f"\n{'variant':<28} {'inside_moll':>12} {'delta':>7}")
        for name in per_variant_calls:
            print(f"{name:<28} {results[name]['inside_moll']:>12} "
                  f"{results[name]['inside_moll'] - results['baseline']['inside_moll']:>+7}")

    print("\n  candidates assigned to `chemokine` by raw argmin in the baseline "
          "(all of\n  which are argmin artifacts -- Mollusca has no chemokine "
          "receptors):")
    ch = b[b["best_family"] == "chemokine"]
    print(f"    n={len(ch)}, of which called inside the chemokine envelope: "
          f"{int(ch['inside_vert'].sum())}")

    for name, cdf in per_variant_calls.items():
        cdf.reset_index().to_csv(f"{a.out_prefix}_{name}.tsv", sep="\t", index=False)
    with open(f"{a.out_prefix}_summary.json", "w") as fh:
        json.dump({"model": a.model, "variants": results,
                   "pct_threshold": a.pct_threshold,
                   "moll_threshold": a.moll_threshold,
                   "note": "MEASURED ONLY -- no variant is applied to production"},
                  fh, indent=2)
    print(f"\n[coh] wrote {a.out_prefix}_*.tsv + _summary.json")
    print("[coh] MEASURED ONLY. No family definition was changed in production.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
