#!/usr/bin/env python3
"""molluscan_calibration_null.py — the molluscan background null, and what the
790 Berghia candidates look like when measured against it.

THE MEASURED PROBLEM THIS ADDRESSES
-----------------------------------
`embedding_family_assignment.py` decides whether a Berghia candidate is inside
non-chemoreceptor family F by comparing its distance to F against the
distribution of distances F's OWN reference members have to F. The reference
set is 1094 anchors, ~61 of them molluscan; per family (total/molluscan):
peptide 386/20, aminergic 212/22, opsin 105/4, lipid 70/1, chemokine 34/0,
nucleotide 23/0, glycoprotein-hormone 8/0. For four of the seven families the
within-family spread is therefore entirely vertebrate-to-vertebrate, while
every query distance is mollusc-to-vertebrate across ~550 My of divergence.
That is a systematic, directional inflation applied to every candidate
regardless of its true family membership.

Reference expansion cannot fix it: reviewed Swiss-Prot holds only 19 class-A
GPCRs for all of Lophotrochozoa and all 19 are already anchors. The ceiling is
exhausted. So the correction has to come from the other side -- calibrate
against a background of REAL molluscan class-A GPCRs and ask whether a Berghia
candidate is closer to family F than a typical mollusc is.

WHAT IS COMPUTED
----------------
The geometry is taken verbatim from the production readout
(`embedding_family_assignment.py` / `embedding_evidence.mahalanobis_channel`):
same k=3 family prototypes, same tied Ledoit-Wolf shrinkage precision, same
min-over-prototypes squared Mahalanobis distance, same RAW-argmin family call,
same `fusion_consensus.confound_residuals` length deconfounding of the
envelope comparison. The ONLY thing that changes is the population a
candidate's distance is compared against.

  vertebrate envelope test (status quo): pct of F's own reference members
      whose residual is below the candidate's. inside if pct <= 95.
  molluscan null test (this script): pct of the molluscan background whose
      residual to that SAME family is below the candidate's. "inside" now
      means unusually close to F for a mollusc -- a lower-tail call, reported
      over a threshold sweep rather than at one asserted cutoff.

RESIDUAL BASIS. `confound_residuals` is fit over the POOLED reference +
candidate + null set so all three populations land on one common scale. The
production readout pools reference + candidate only, so the basis differs and
the vertebrate-envelope numbers on this basis need not reproduce the published
91/16 exactly. Both are therefore computed and the basis-only shift is
reported, so no difference is silently attributed to the calibration.

INTERPRETIVE FRAME. The pipeline is SUBTRACTIVE: a confident EXCLUSION is the
load-bearing output and the unclassified residual is NEVER a positive claim of
chemoreceptor identity. Nothing here is a production ranking; these are
calibration/validation artifacts. There is no positive control for molluscan
odorant receptors and there cannot be one.
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

THRESHOLDS = (0.1, 1.0, 2.5, 5.0, 10.0)
CONVERGENCE_N = (250, 500, 1000, 2000, 4000, 6000, 8000, 10000)


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


def pct_below(pop: np.ndarray, v: float) -> float:
    """Percent of `pop` strictly below v -- the same statistic
    embedding_family_assignment.py uses for pct_within_family."""
    return float((pop < v).sum()) / len(pop) * 100.0


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--model", required=True, help="proteinclip3b | protrek")
    p.add_argument("--ref-npz", required=True)
    p.add_argument("--candidate-npz", required=True)
    p.add_argument("--null-npz", required=True)
    p.add_argument("--ref-labels", required=True)
    p.add_argument("--candidate-fasta", required=True)
    p.add_argument("--anchor-fasta", required=True)
    p.add_argument("--null-fasta", required=True)
    p.add_argument("--sample-tsv", required=True,
                   help="molluscan_calibration_sample.tsv: stratum + taxon scope")
    p.add_argument("--k", type=int, default=3)
    p.add_argument("--pct-threshold", type=float, default=95.0)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--seed", type=int, default=20260720)
    a = p.parse_args()

    rng = np.random.default_rng(a.seed)

    # ---- load ------------------------------------------------------------
    ref = load_embeddings(a.ref_npz)
    cand = load_embeddings(a.candidate_npz)
    null = load_embeddings(a.null_npz)
    labels = novelty_reference_labels(load_ref_labels(a.ref_labels))
    labels = {i: f for i, f in labels.items() if i in ref}
    if not (ref and cand and null and labels):
        print("[null] FATAL: empty input", file=sys.stderr)
        return 1

    meta = pd.read_csv(a.sample_tsv, sep="\t", dtype=str).set_index("seq_id")
    missing_meta = [i for i in null if i not in meta.index]
    if missing_meta:
        print(f"[null] FATAL: {len(missing_meta)} null ids have no sample-TSV "
              f"row (e.g. {missing_meta[:3]}); stratum/scope joins would be "
              "silently partial", file=sys.stderr)
        return 1

    cand_len = read_fasta_lengths(a.candidate_fasta)
    ref_len = read_fasta_lengths(a.anchor_fasta)
    null_len = read_fasta_lengths(a.null_fasta)
    for nm, keys, src in (("candidate", cand, cand_len),
                          ("reference", labels, ref_len),
                          ("null", null, null_len)):
        miss = [i for i in keys if i not in src]
        if miss:
            print(f"[null] FATAL: seq_len unresolved for {len(miss)} {nm} ids "
                  f"(e.g. {miss[:3]})", file=sys.stderr)
            return 1
    print(f"[null] model={a.model} refs={len(labels)} candidates={len(cand)} "
          f"null={len(null)}; all seq_len keys resolved")

    # ---- geometry (identical to the production readout) -------------------
    protos = family_prototypes(ref, labels, a.k)
    by_family: dict = {}
    for i, f in labels.items():
        by_family.setdefault(f, []).append(np.asarray(ref[i], dtype=float))
    centered = np.vstack([np.array(v) - np.array(v).mean(axis=0)
                          for v in by_family.values()])
    precision = shrinkage_precision(centered)
    fams = sorted(protos)
    print(f"[null] families ({len(fams)}): " +
          ", ".join(f"{f}={len(by_family[f])}" for f in fams))

    ref_ids, cand_ids, null_ids = list(labels), list(cand), list(null)

    def dmat(store, ids):
        return pd.DataFrame(
            [{f: dist_to_family(np.asarray(store[i], dtype=float), protos[f], precision)
              for f in fams} for i in ids], index=ids)

    d_ref, d_cand, d_null = dmat(ref, ref_ids), dmat(cand, cand_ids), dmat(null, null_ids)

    # ---- residual bases ---------------------------------------------------
    def residualize(frames, lens):
        pooled = pd.concat(frames)
        if pooled.index.duplicated().any():
            raise SystemExit("[null] FATAL: id collision in the pooled residual set")
        r = confound_residuals({f: pooled[f].to_dict() for f in fams},
                               {"seq_len": lens})
        if {len(r[f]) for f in fams} != {len(pooled)}:
            raise SystemExit("[null] FATAL: confound_residuals partial coverage")
        return pd.DataFrame({f: pd.Series(r[f]) for f in fams})

    L_ref = {i: float(ref_len[i]) for i in ref_ids}
    L_cand = {i: float(cand_len[i]) for i in cand_ids}
    L_null = {i: float(null_len[i]) for i in null_ids}

    # basis P = production readout (ref + cand) -> reproduces the published call
    rP = residualize([d_ref[fams], d_cand[fams]], {**L_ref, **L_cand})
    # basis C = calibration (ref + cand + null) -> one scale for all three
    rC = residualize([d_ref[fams], d_cand[fams], d_null[fams]],
                     {**L_ref, **L_cand, **L_null})
    print(f"[null] residual bases: production n={len(rP)}, calibration n={len(rC)}")

    ref_fam = pd.Series(labels)
    raw_argmin = d_cand[fams].idxmin(axis=1)

    def vert_envelope(rdf):
        return {f: rdf.loc[ref_fam[ref_fam == f].index, f].to_numpy() for f in fams}

    envP, envC = vert_envelope(rP), vert_envelope(rC)

    # ---- candidate table --------------------------------------------------
    rows = []
    for cid in cand_ids:
        best = raw_argmin[cid]
        vP = float(rP.loc[cid, best])
        vC = float(rC.loc[cid, best])
        rows.append({
            "id": cid,
            "best_family": best,
            "seq_len": cand_len[cid],
            "best_dist_raw": float(d_cand.loc[cid, best]),
            "pct_vert_envelope_prodbasis": pct_below(envP[best], vP),
            "pct_vert_envelope_calbasis": pct_below(envC[best], vC),
            "resid_calbasis": vC,
        })
    C = pd.DataFrame(rows).set_index("id")

    # molluscan-null percentile of each candidate, per stratum / scope
    null_meta = meta.loc[null_ids]
    scopes = {
        "mollusca_all": np.ones(len(null_ids), dtype=bool),
        "postgate_only": (null_meta["stratum"] == "postgate_classA").to_numpy(),
        "pregate_only": (null_meta["stratum"] == "pregate_only").to_numpy(),
        "gastropoda": (null_meta["is_gastropoda"] == "1").to_numpy(),
        "heterobranchia": (null_meta["is_heterobranchia"] == "1").to_numpy(),
        "bivalvia": (null_meta["tax_class"] == "Bivalvia").to_numpy(),
        "cephalopoda": (null_meta["tax_class"] == "Cephalopoda").to_numpy(),
        "gastropoda_postgate": ((null_meta["is_gastropoda"] == "1") &
                                (null_meta["stratum"] == "postgate_classA")).to_numpy(),
    }
    # detection-evidence strata bound the HMM step's chemoreceptor bias from
    # INSIDE the post-gate sample (see molluscan_calibration_augment_sample.py)
    if "detection_class" in null_meta.columns:
        for det in ("chemo_classA", "chemo_offclan", "generic_7tm1"):
            scopes[f"det_{det}"] = (null_meta["detection_class"] == det).to_numpy()
        # clan-clean sensitivity: drop the CL0176 Chemosens_recp sequences,
        # which classify_gpcr_by_class.py calls class A but which are not in
        # the rhodopsin clan this project gates on
        scopes["clan_clean_CL0192"] = (null_meta["detection_class"]
                                       != "chemo_offclan").to_numpy()
    rC_null = rC.loc[null_ids]
    for name, mask in scopes.items():
        n = int(mask.sum())
        if n < 30:
            print(f"[null] scope {name}: n={n} < 30, skipped")
            continue
        sub = rC_null[mask]
        C[f"pct_moll_{name}"] = [
            pct_below(sub[C.loc[cid, "best_family"]].to_numpy(),
                      float(C.loc[cid, "resid_calbasis"]))
            for cid in C.index]

    C["inside_vert_prodbasis"] = C["pct_vert_envelope_prodbasis"] <= a.pct_threshold
    C["inside_vert_calbasis"] = C["pct_vert_envelope_calbasis"] <= a.pct_threshold
    for q in THRESHOLDS:
        C[f"inside_moll_q{q}"] = C["pct_moll_mollusca_all"] <= q

    C.reset_index().to_csv(f"{a.out_prefix}_candidates.tsv", sep="\t", index=False)

    # ---- deliverable 1: the null itself ----------------------------------
    null_rows = []
    for j, nid in enumerate(null_ids):
        best = d_null.loc[nid, fams].idxmin()
        rec = {"seq_id": nid, "stratum": null_meta.loc[nid, "stratum"],
               "detection_class": (null_meta.loc[nid, "detection_class"]
                                   if "detection_class" in null_meta.columns else ""),
               "tax_class": null_meta.loc[nid, "tax_class"],
               "is_gastropoda": null_meta.loc[nid, "is_gastropoda"],
               "is_heterobranchia": null_meta.loc[nid, "is_heterobranchia"],
               "sample": null_meta.loc[nid, "sample"],
               "seq_len": null_len[nid], "best_family": best}
        for f in fams:
            rec[f"pct_vert_{f}"] = pct_below(envC[f], float(rC_null.loc[nid, f]))
        rec["pct_vert_best"] = rec[f"pct_vert_{best}"]
        null_rows.append(rec)
    N = pd.DataFrame(null_rows).set_index("seq_id")
    N.reset_index().to_csv(f"{a.out_prefix}_null.tsv", sep="\t", index=False)

    print("\n" + "=" * 78)
    print(f"D1. EXPECTED PERCENTILE OF A RANDOM MOLLUSCAN CLASS-A RECEPTOR "
          f"AGAINST EACH\n    VERTEBRATE FAMILY ENVELOPE  [{a.model}]")
    print("=" * 78)
    print(f"{'family':<22} {'n_ref':>6} {'median':>8} {'p25':>7} {'p75':>7} "
          f"{'>=100':>7} {'<=95':>7}")
    d1 = {}
    for f in fams:
        v = N[f"pct_vert_{f}"].to_numpy()
        d1[f] = {"n_ref_members": int((ref_fam == f).sum()),
                 "median_pct": float(np.median(v)),
                 "p25": float(np.percentile(v, 25)),
                 "p75": float(np.percentile(v, 75)),
                 "frac_at_100": float((v >= 100.0).mean()),
                 "frac_inside_95": float((v <= a.pct_threshold).mean())}
        s = d1[f]
        print(f"{f:<22} {s['n_ref_members']:>6} {s['median_pct']:>8.2f} "
              f"{s['p25']:>7.2f} {s['p75']:>7.2f} {100*s['frac_at_100']:>6.1f}% "
              f"{100*s['frac_inside_95']:>6.1f}%")
    print("\n  Read: 'median' is where a typical molluscan class-A receptor sits in the")
    print("  spread of family F's own vertebrate members. 100 = further from F than")
    print("  EVERY reference member of F. '<=95' is the fraction of ordinary molluscan")
    print("  receptors the status-quo rule would call INSIDE family F.")

    # ---- convergence -------------------------------------------------------
    print("\n" + "=" * 78)
    print("D1b. CONVERGENCE OF THE NULL (does the sample size support the estimate?)")
    print("=" * 78)
    conv = []
    for n in CONVERGENCE_N:
        if n > len(null_ids):
            continue
        reps = []
        for _ in range(50):
            idx = rng.choice(len(null_ids), size=n, replace=False)
            sub = rC_null.iloc[idx]
            reps.append([np.percentile(sub[f].to_numpy(), 5) for f in fams])
        arr = np.asarray(reps)
        conv.append({"n": n, **{f"sd_p5_{f}": float(arr[:, i].std(ddof=1))
                                for i, f in enumerate(fams)},
                     "max_sd_rel": float(np.max(arr.std(axis=0, ddof=1) /
                                                np.abs(arr.mean(axis=0))))})
        print(f"  n={n:>6}  max relative SD of the per-family 5th-percentile "
              f"estimate over 50 draws = {conv[-1]['max_sd_rel']:.4f}")
    pd.DataFrame(conv).to_csv(f"{a.out_prefix}_convergence.tsv", sep="\t", index=False)

    # ---- deliverable 2: rescore -------------------------------------------
    print("\n" + "=" * 78)
    print(f"D2. THE 790 BERGHIA CANDIDATES RESCORED AGAINST THE MOLLUSCAN NULL "
          f"[{a.model}]")
    print("=" * 78)
    n_vp = int(C["inside_vert_prodbasis"].sum())
    n_vc = int(C["inside_vert_calbasis"].sum())
    print(f"  inside a vertebrate family envelope, production residual basis : {n_vp}")
    print(f"  inside a vertebrate family envelope, calibration residual basis: {n_vc}")
    print(f"  (basis-only shift, attributable to pooling the null into the "
          f"residual fit: {n_vc - n_vp:+d})")
    print(f"\n  molluscan-null percentile of the candidates the STATUS QUO calls inside:")
    ins = C[C["inside_vert_prodbasis"]]
    if len(ins):
        v = ins["pct_moll_mollusca_all"]
        print(f"    n={len(ins)}  median={v.median():.2f}  p25={v.quantile(.25):.2f} "
              f" p75={v.quantile(.75):.2f}  frac<=5={float((v <= 5).mean()):.3f}")
    print(f"\n  {'threshold q':>12} {'inside_moll':>12} {'overlap w/ vert-inside':>24}")
    d2 = {"inside_vert_prodbasis": n_vp, "inside_vert_calbasis": n_vc}
    for q in THRESHOLDS:
        col = C[f"inside_moll_q{q}"]
        ov = int((col & C["inside_vert_prodbasis"]).sum())
        print(f"  {q:>12} {int(col.sum()):>12} {ov:>24}")
        d2[f"inside_moll_q{q}"] = int(col.sum())
        d2[f"overlap_q{q}"] = ov
    print("\n  full distribution of the candidates' molluscan-null percentile:")
    s = C["pct_moll_mollusca_all"]
    print("   " + "  ".join(f"p{q}={np.percentile(s, q):.1f}"
                            for q in (0, 1, 5, 10, 25, 50, 75, 90, 100)))

    # ---- deliverable 3: the enrichment caveat ------------------------------
    print("\n" + "=" * 78)
    print("D3. THE 6TM-CHEMO-GATE ENRICHMENT CAVEAT, QUANTIFIED")
    print("=" * 78)
    print("  The scan_record.tsv files DO retain the pre-gate set relative to the TM")
    print("  gate: every row is HMM-GPCR-positive and `passed_gate` records the >=6TM")
    print("  outcome. So the null can be computed both ways. What CANNOT be recovered")
    print("  is a pre-HMM-detection sample: the HMM stack itself (TIAMMAT + curated")
    print("  lophotrochozoan chemoreceptor HMMs) is chemoreceptor-biased, and no")
    print("  artifact of this run contains the sequences it rejected.")
    d3 = {}
    for f in fams:
        a_ = N.loc[N["stratum"] == "postgate_classA", f"pct_vert_{f}"].to_numpy()
        b_ = N.loc[N["stratum"] == "pregate_only", f"pct_vert_{f}"].to_numpy()
        d3[f] = {"postgate_median": float(np.median(a_)),
                 "pregate_median": float(np.median(b_)),
                 "postgate_p5": float(np.percentile(a_, 5)),
                 "pregate_p5": float(np.percentile(b_, 5)),
                 "n_postgate": int(len(a_)), "n_pregate": int(len(b_))}
    print(f"\n  {'family':<22} {'post med':>9} {'pre med':>9} {'post p5':>9} "
          f"{'pre p5':>9}")
    for f in fams:
        s = d3[f]
        print(f"  {f:<22} {s['postgate_median']:>9.2f} {s['pregate_median']:>9.2f} "
              f"{s['postgate_p5']:>9.2f} {s['pregate_p5']:>9.2f}")
    print(f"\n  n_postgate={d3[fams[0]]['n_postgate']}  "
          f"n_pregate_only={d3[fams[0]]['n_pregate']}")
    print("  seq_len (aa): post-gate median "
          f"{N.loc[N['stratum'] == 'postgate_classA', 'seq_len'].median():.0f}, "
          f"pre-gate-only median "
          f"{N.loc[N['stratum'] == 'pregate_only', 'seq_len'].median():.0f}")
    print("  NOTE: the pre-gate-only stratum is TM-FAILED, i.e. enriched for")
    print("  fragments and incomplete models. It is not a neutral sample either;")
    print("  it bounds the gate effect from the opposite side, it does not remove it.")

    # within-sample bound on the HMM detection step's chemoreceptor bias
    if "detection_class" in N.columns and (N["detection_class"] != "").any():
        print("\n  WITHIN-SAMPLE BOUND on the HMM detection step (which ran BEFORE")
        print("  the TM gate and used chemoreceptor-biased profiles). Post-gate")
        print("  sequences split by which Pfam actually hit (clan from InterPro):")
        print("    chemo_classA  = PF10324 Srw / PF05296 TAS2R   clan CL0192 GPCR_A")
        print("    generic_7tm1  = PF00001 rhodopsin             clan CL0192 GPCR_A")
        print("    chemo_offclan = PF08395 7tm_7 / PF02949 7tm_6 clan CL0176 "
              "(NOT class A;\n                    reported separately as "
              "contamination, not enrichment)")
        print(f"\n  {'family':<22} {'chemo med':>10} {'generic med':>12} "
              f"{'chemo p5':>10} {'generic p5':>11}")
        d3det = {}
        for f in fams:
            c_ = N.loc[N["detection_class"] == "chemo_classA", f"pct_vert_{f}"].to_numpy()
            g_ = N.loc[N["detection_class"] == "generic_7tm1", f"pct_vert_{f}"].to_numpy()
            if len(c_) < 30 or len(g_) < 30:
                continue
            d3det[f] = {"chemo_median": float(np.median(c_)),
                        "generic_median": float(np.median(g_)),
                        "chemo_p5": float(np.percentile(c_, 5)),
                        "generic_p5": float(np.percentile(g_, 5)),
                        "n_chemo": int(len(c_)), "n_generic": int(len(g_))}
            s = d3det[f]
            print(f"  {f:<22} {s['chemo_median']:>10.2f} {s['generic_median']:>12.2f} "
                  f"{s['chemo_p5']:>10.2f} {s['generic_p5']:>11.2f}")
        if d3det:
            print(f"\n  n_chemo_pfam={list(d3det.values())[0]['n_chemo']}  "
                  f"n_generic_7tm1={list(d3det.values())[0]['n_generic']}")
            print("  A SMALL gap here means the null is dominated by 'these are")
            print("  molluscs', not by the chemoreceptor enrichment. A LARGE gap")
            print("  means the enrichment is doing real work and the null is")
            print("  optimistic by roughly that much.")
        d3["_detection_contrast"] = d3det

    # candidate calls under each null stratum
    print("\n  candidate inside-counts under each null stratum (q<=5):")
    for name in ("mollusca_all", "postgate_only", "pregate_only"):
        col = f"pct_moll_{name}"
        if col in C:
            print(f"    {name:<16} {int((C[col] <= 5.0).sum()):>5}")

    # ---- deliverable 3b: phylogenetic scope -------------------------------
    print("\n" + "=" * 78)
    print("D3b. IS THE NULL DOMINATED BY PHYLOGENETIC DISTANCE WITHIN MOLLUSCA?")
    print("=" * 78)
    print(f"  {'scope':<22} {'n':>7} {'median pct_vert_best':>22}")
    d3b = {}
    for name, mask in scopes.items():
        if mask.sum() < 30:
            continue
        sub = N[mask]
        d3b[name] = {"n": int(mask.sum()),
                     "median_pct_vert_best": float(sub["pct_vert_best"].median())}
        print(f"  {name:<22} {int(mask.sum()):>7} "
              f"{sub['pct_vert_best'].median():>22.2f}")
    print(f"\n  candidate inside-counts (q<=5) under each scope's null:")
    for name in scopes:
        col = f"pct_moll_{name}"
        if col in C:
            n_in = int((C[col] <= 5.0).sum())
            d3b.setdefault(name, {})["candidates_inside_q5"] = n_in
            print(f"    {name:<22} {n_in:>5}")

    # per-family median by taxon class, so a bivalve-vs-gastropod gap is visible
    print("\n  median pct_vert by family x taxon class:")
    print(f"  {'family':<22} " + " ".join(f"{c:>12}" for c in
                                          ("Gastropoda", "Bivalvia", "Cephalopoda")))
    for f in fams:
        vals = []
        for c in ("Gastropoda", "Bivalvia", "Cephalopoda"):
            sub = N.loc[N["tax_class"] == c, f"pct_vert_{f}"]
            vals.append(f"{sub.median():>12.2f}" if len(sub) >= 30 else f"{'n/a':>12}")
        print(f"  {f:<22} " + " ".join(vals))

    summary = {"model": a.model, "n_ref": len(labels), "n_candidates": len(cand),
               "n_null": len(null), "families": fams,
               "d1_family_null": d1, "d2_rescore": d2,
               "d3_gate_enrichment": d3, "d3b_scope": d3b,
               "convergence": conv,
               "pct_threshold": a.pct_threshold, "seed": a.seed}
    with open(f"{a.out_prefix}_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"\n[null] wrote {a.out_prefix}_{{candidates,null,convergence}}.tsv "
          f"+ _summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
