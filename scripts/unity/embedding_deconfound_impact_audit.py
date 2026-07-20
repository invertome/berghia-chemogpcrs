#!/usr/bin/env python3
"""embedding_deconfound_impact_audit.py — quantify HOW WRONG the earlier RAW
family-assignment numbers were, against the deconfounded numbers that stand.

THIS IS AN AUDIT, NOT A PRODUCER. The canonical producer is
`embedding_family_assignment.py`, which deconfounds seq_len BY DEFAULT and
refuses to emit a raw result. This script exists only to measure the size of the
error the raw runs introduced, so that the two analyses already built on raw
distances can be corrected by a known amount. Deconfounded is not one option on
a menu: the deconfounded proteinclip3b+protrek consensus is a LOCKED decision,
and the raw column here is the WRONG answer shown for calibration.

WHY IT WAS WRONG. `embedding_family_assignment.py` used to compute RAW
squared-Mahalanobis distance. The production novelty channel
(07_candidate_ranking.sh) runs `fusion_consensus.py --combiner rra --deconfound
seq_len`, which residualizes rank(score) on rank(seq_len) before anything
downstream consumes it. Every envelope-membership and family call produced
before this fix was therefore computed on a quantity production does not use.
Direct evidence the length artifact was live: the candidates inside
proteinclip3b's envelope but not protrek's were significantly LONGER (median
469.5 vs 376 aa).

WHAT THIS DOES. The embedding geometry is reproduced VERBATIM from
embedding_family_assignment.py (same prototypes, same tied shrinkage precision,
same min-over-prototypes distance), then a deconfounded pass is added that
reuses `fusion_consensus.confound_residuals` — the exact production
residualizer, not a reimplementation.

THE ONE DESIGN DECISION, STATED OPENLY. Production residualizes over the 790
candidates alone, because that is all it scores. Envelope membership, however,
is a comparison BETWEEN a candidate's distance and the reference members' own
distances, so both populations must land on one common scale or the comparison
is meaningless. So the residualization here is fit over the POOLED set
(labelled references + candidates), once per family-distance column, using one
shared id ordering. Consequences, both real:
  * every family column is residualized over the SAME pooled set, so the
    residuals stay comparable across families and argmin remains well-defined;
  * the pooled rank basis mixes two populations, so if references of family F
    are systematically shorter/longer than candidates AND systematically closer
    to F, the length coefficient absorbs part of that. The reference-recovery
    check below is what detects it.

SIDE-EFFECT PROBES (flagged observations, NOT decision inputs). Deconfounding is
locked, so these do not gate anything; they are recorded so that any collateral
damage is documented with evidence rather than discovered later:
  1. Kruskal-Wallis of seq_len across the 7 reference families — is length
     correlated with family identity in the labelled data at all?
  2. Reference self-recovery: does a labelled reference still land nearest its
     OWN family, raw vs deconfounded? (Prototypes include the point being
     scored, so the absolute rate is optimistic; the RAW-vs-DECONFOUNDED DELTA
     is the interpretable quantity, since that bias is identical in both.)
A large negative delta would be worth flagging in its own right — it would mean
the residual is also removing family-linked length structure — but it is
reported as an observation, never as a reason to quote the raw numbers.

INTERPRETIVE FRAME (must survive into any report). The pipeline is SUBTRACTIVE.
It classifies every OTHER molluscan GPCR family confidently; what remains
unclassified is the enriched residual. There is no positive control for
molluscan odorant receptors and there cannot be one. A confident EXCLUSION is
the load-bearing output; "unassigned" is a residual, NEVER a positive claim of
chemoreceptor identity. A candidate that only ONE model excludes is a candidate
whose exclusion is NOT safe — it stays in the chemoreceptor pool.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

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

PCT_THRESHOLD = 95.0


def hdr(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def read_fasta_lengths(path: str) -> dict:
    """{first_whitespace_token_of_header: residue count}."""
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


def dist_to_family(vec: np.ndarray, protos: np.ndarray, precision: np.ndarray) -> float:
    return min(mahalanobis_sq(vec, p, precision) for p in np.atleast_2d(protos))


def mwu(name: str, a, b, la: str, lb: str) -> float:
    """Mann-Whitney U with medians and rank-biserial effect size. Returns p."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3:
        print(f"  {name:<30} NOT TESTABLE (n={len(a)}/{len(b)}, too few)")
        return float("nan")
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    rbc = 2.0 * u / (len(a) * len(b)) - 1.0
    print(f"  {name:<30} median {la}={np.median(a):9.2f}  {lb}={np.median(b):9.2f}"
          f"   n={len(a)}/{len(b)}  U={u:9.1f}  p={p:.3e}  rank-biserial={rbc:+.3f}")
    return float(p)


def build_geometry(cand_npz: str, ref_npz: str, ref_labels_tsv: str, k: int):
    """Prototypes + tied precision + raw distance frames. Verbatim geometry."""
    cand = load_embeddings(cand_npz)
    ref = load_embeddings(ref_npz)
    labels = novelty_reference_labels(load_ref_labels(ref_labels_tsv))
    labels = {i: f for i, f in labels.items() if i in ref}
    if not cand or not labels:
        raise SystemExit("FATAL: empty inputs")

    protos = family_prototypes(ref, labels, k)
    by_family: dict = {}
    for i, f in labels.items():
        by_family.setdefault(f, []).append(np.asarray(ref[i], dtype=float))
    centered = np.vstack([np.array(v) - np.array(v).mean(axis=0)
                          for v in by_family.values()])
    precision = shrinkage_precision(centered)
    fams = sorted(protos)

    def dist_frame(emb: dict, ids: list) -> pd.DataFrame:
        rows = []
        for i in ids:
            v = np.asarray(emb[i], dtype=float)
            rows.append({"id": i, **{f: dist_to_family(v, protos[f], precision)
                                     for f in fams}})
        return pd.DataFrame(rows).set_index("id")

    ref_ids = list(labels)
    cand_ids = list(cand)
    d_ref = dist_frame(ref, ref_ids)
    d_cand = dist_frame(cand, cand_ids)
    return fams, labels, d_ref, d_cand


def assign(d_cand: pd.DataFrame, d_ref: pd.DataFrame, labels: dict,
           fams: list) -> pd.DataFrame:
    """argmin family + envelope percentile against that family's own members."""
    ref_fam = pd.Series(labels)
    envelope = {f: d_ref.loc[ref_fam[ref_fam == f].index, f].to_numpy()
                for f in fams}
    arr = d_cand[fams].to_numpy()
    order = np.argsort(arr, axis=1)
    best_i, second_i = order[:, 0], order[:, 1]
    rows = []
    for r, (cid, bi, si) in enumerate(zip(d_cand.index, best_i, second_i)):
        best, second = fams[bi], fams[si]
        env = envelope[best]
        rows.append({
            "id": cid,
            "best_family": best,
            "best_dist": arr[r, bi],
            "pct_within_family": float((env < arr[r, bi]).sum()) / len(env) * 100.0,
            "runner_up_family": second,
            "margin": arr[r, si] - arr[r, bi],
        })
    out = pd.DataFrame(rows).set_index("id")
    out["inside_envelope"] = out["pct_within_family"] <= PCT_THRESHOLD
    return out


def assign_hybrid(d_cand_raw: pd.DataFrame, r_cand: pd.DataFrame,
                  r_ref: pd.DataFrame, labels: dict, fams: list) -> pd.DataFrame:
    """Family argmin from RAW distance; envelope percentile from the RESIDUAL.

    Rationale for the split, which the full-deconfound pass appears to violate:
      * choosing WHICH family a candidate is nearest is a WITHIN-sequence
        comparison — every family distance for that one candidate carries the
        same sequence length, so length is a common factor that largely cancels
        in the argmin. Residualizing each family column separately destroys
        that cancellation and re-ranks families against each other.
      * whether the candidate sits inside the family's envelope is a BETWEEN-
        sequence comparison (candidate vs reference members of different
        lengths) — this is where the length artifact actually operates, so this
        is the comparison that needs the residual.
    """
    ref_fam = pd.Series(labels)
    envelope = {f: r_ref.loc[ref_fam[ref_fam == f].index, f].to_numpy()
                for f in fams}
    raw = d_cand_raw[fams].to_numpy()
    order = np.argsort(raw, axis=1)
    rows = []
    for k, (cid, bi, si) in enumerate(zip(d_cand_raw.index,
                                          order[:, 0], order[:, 1])):
        best, second = fams[bi], fams[si]
        env = envelope[best]
        rv = float(r_cand.loc[cid, best])
        rows.append({
            "id": cid,
            "best_family": best,
            "best_dist": rv,
            "pct_within_family": float((env < rv).sum()) / len(env) * 100.0,
            "runner_up_family": second,
            "margin": raw[k, si] - raw[k, bi],
        })
    out = pd.DataFrame(rows).set_index("id")
    out["inside_envelope"] = out["pct_within_family"] <= PCT_THRESHOLD
    return out


def ref_self_recovery(d_ref: pd.DataFrame, labels: dict, fams: list) -> float:
    """Fraction of labelled references whose argmin family is their OWN."""
    pred = d_ref[fams].idxmin(axis=1)
    truth = pd.Series(labels).reindex(pred.index)
    return float((pred == truth).mean())


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--emb-dir", required=True)
    p.add_argument("--ref-labels", required=True)
    p.add_argument("--candidate-fasta", required=True)
    p.add_argument("--anchor-fasta", required=True)
    p.add_argument("--diag-dir", required=True)
    p.add_argument("--model-a", default="proteinclip3b")
    p.add_argument("--model-b", default="protrek")
    p.add_argument("--k", type=int, default=3)
    p.add_argument("--out-prefix", required=True)
    a = p.parse_args()

    # ---------------------------------------------------------------- lengths
    cand_len = read_fasta_lengths(a.candidate_fasta)
    ref_len = read_fasta_lengths(a.anchor_fasta)
    print(f"[len] candidate FASTA={len(cand_len)}  anchor FASTA={len(ref_len)}")

    results, frames = {}, {}
    for tag in (a.model_a, a.model_b):
        hdr(f"MODEL {tag}: geometry")
        fams, labels, d_ref, d_cand = build_geometry(
            os.path.join(a.emb_dir, f"candidates_{tag}_classA.npz"),
            os.path.join(a.emb_dir, f"reference_{tag}_PROD.npz"),
            a.ref_labels, a.k)
        print(f"  candidates={len(d_cand)}  references_resolved={len(labels)}  "
              f"families={len(fams)} {fams}")

        # -- key-overlap assertions on REAL data (never trust a fixture) -----
        miss_c = [i for i in d_cand.index if i not in cand_len]
        miss_r = [i for i in d_ref.index if i not in ref_len]
        if miss_c or miss_r:
            raise SystemExit(f"FATAL[{tag}]: seq_len unresolved for "
                             f"{len(miss_c)} candidates / {len(miss_r)} refs "
                             f"(e.g. {(miss_c + miss_r)[:3]})")
        print(f"  seq_len resolved: {len(d_cand)}/{len(d_cand)} candidates, "
              f"{len(d_ref)}/{len(d_ref)} references")

        # -- cross-check candidate seq_len against the production novelty TSV
        nv = pd.read_csv(os.path.join(a.diag_dir, f"novelty_{tag}_PROD.tsv"),
                         sep="\t").set_index("candidate_id")
        shared = [i for i in d_cand.index if i in nv.index]
        delta = int(sum(int(nv.loc[i, "seq_len"]) != cand_len[i] for i in shared))
        print(f"  seq_len cross-check vs novelty_{tag}_PROD.tsv: "
              f"{len(shared)} shared ids, {delta} mismatches")
        if delta:
            raise SystemExit(f"FATAL[{tag}]: FASTA seq_len disagrees with the "
                             f"production novelty TSV on {delta} candidates")

        # ------------------------------------------------------------- RAW --
        raw = assign(d_cand, d_ref, labels, fams)

        # Faithfulness: my geometry must reproduce the on-disk RAW artifact.
        disk = pd.read_csv(os.path.join(a.diag_dir,
                                        f"family_assignment_{tag}.tsv"),
                           sep="\t").set_index("id")
        common = raw.index.intersection(disk.index)
        dmax = float(np.abs(raw.loc[common, "best_dist"]
                            - disk.loc[common, "best_dist"]).max())
        fam_same = int((raw.loc[common, "best_family"]
                        == disk.loc[common, "best_family"]).sum())
        env_same = int((raw.loc[common, "inside_envelope"]
                        == disk.loc[common, "inside_envelope"].astype(bool)).sum())
        print(f"  [reproduce on-disk RAW] n={len(common)}  "
              f"max|Δbest_dist|={dmax:.3e}  family match={fam_same}/{len(common)}  "
              f"envelope match={env_same}/{len(common)}")
        if fam_same != len(common) or env_same != len(common):
            raise SystemExit(f"FATAL[{tag}]: could not reproduce the on-disk RAW "
                             "assignment; geometry differs, deconfounded numbers "
                             "would not be comparable")

        # --------------------------------------------------- DECONFOUNDED --
        # One pooled id set, one shared confound map, every family column
        # residualized over it -> residual scales comparable across families.
        pooled_len = {**{i: float(ref_len[i]) for i in d_ref.index},
                      **{i: float(cand_len[i]) for i in d_cand.index}}
        pooled_d = pd.concat([d_ref[fams], d_cand[fams]])
        if pooled_d.index.duplicated().any():
            raise SystemExit(f"FATAL[{tag}]: reference/candidate id collision")
        by_col = {f: pooled_d[f].to_dict() for f in fams}
        resid = confound_residuals(by_col, {"seq_len": pooled_len})
        cov = {f: len(resid[f]) for f in fams}
        if set(cov.values()) != {len(pooled_d)}:
            raise SystemExit(f"FATAL[{tag}]: residual coverage incomplete {cov}")
        print(f"  [deconfound] confound_residuals over pooled n={len(pooled_d)} "
              f"({len(d_ref)} ref + {len(d_cand)} cand) x {len(fams)} families")

        r_all = pd.DataFrame({f: pd.Series(resid[f]) for f in fams})
        r_ref = r_all.loc[d_ref.index]
        r_cand = r_all.loc[d_cand.index]
        dec = assign(r_cand, r_ref, labels, fams)
        hyb = assign_hybrid(d_cand, r_cand, r_ref, labels, fams)

        # ------------------------------------------- over-correction probes --
        rec_raw = ref_self_recovery(d_ref, labels, fams)
        rec_dec = ref_self_recovery(r_ref, labels, fams)
        lens_by_fam = {f: [ref_len[i] for i in d_ref.index if labels[i] == f]
                       for f in fams}
        kw_h, kw_p = stats.kruskal(*[lens_by_fam[f] for f in fams])
        print(f"  [over-correction] reference self-recovery "
              f"raw={rec_raw * 100:.1f}%  deconfounded={rec_dec * 100:.1f}%  "
              f"(Δ={100 * (rec_dec - rec_raw):+.1f} pp)")
        print(f"  [over-correction] seq_len across the {len(fams)} reference "
              f"families: Kruskal-Wallis H={kw_h:.2f} p={kw_p:.3e}")
        for f in fams:
            v = lens_by_fam[f]
            print(f"      {f:<24} n={len(v):<5d} median_len={np.median(v):7.1f}")

        results[tag] = {"raw": raw, "dec": dec, "hyb": hyb,
                        "fams": fams, "labels": labels,
                        "rec_raw": rec_raw, "rec_dec": rec_dec}
        frames[tag] = {"cand_len": cand_len}

    A, B = a.model_a, a.model_b
    fams = results[A]["fams"]
    ref_tbl = pd.read_csv(a.ref_labels, sep="\t")
    fam_n = pd.Series(results[A]["labels"]).value_counts().to_dict()

    # ======================================================== Q1 envelope ===
    MODE_LABEL = {"raw": "RAW (wrong: no deconfounding)",
                  "dec": "DECONFOUNDED (seq_len, every family column)",
                  "hyb": "HYBRID (raw argmin + deconfounded envelope)"}
    for mode in ("raw", "dec", "hyb"):
        label = MODE_LABEL[mode]
        hdr(f"Q1. ENVELOPE MEMBERSHIP — {label}")
        ia = results[A][mode]["inside_envelope"]
        ib = results[B][mode]["inside_envelope"].reindex(ia.index)
        both = ia & ib
        only_a = ia & ~ib
        only_b = ib & ~ia
        neither = ~ia & ~ib
        n = len(ia)
        print(f"  inside BOTH            {int(both.sum()):5d}")
        print(f"  inside {A} only  {int(only_a.sum()):5d}")
        print(f"  inside {B} only        {int(only_b.sum()):5d}")
        print(f"  inside NEITHER         {int(neither.sum()):5d}")
        print(f"  total                  {n:5d}")
        nested = int(only_b.sum()) == 0
        print(f"  {B} ⊂ {A} (strict nesting)? "
              f"{'YES' if nested else 'NO — nesting BROKEN'}")
        results[f"{mode}_sets"] = {"both": set(both[both].index),
                                   "only_a": set(only_a[only_a].index),
                                   "only_b": set(only_b[only_b].index),
                                   "xor": set(only_a[only_a].index)
                                   | set(only_b[only_b].index)}

    # ==================================================== Q2 both-inside ====
    hdr("Q2. THE BOTH-INSIDE SET — raw vs deconfounded")
    braw = results["raw_sets"]["both"]
    bdec = results["dec_sets"]["both"]
    kept, lost, new = sorted(braw & bdec), sorted(braw - bdec), sorted(bdec - braw)
    print(f"  RAW both-inside          n={len(braw)}")
    print(f"  DECONFOUNDED both-inside n={len(bdec)}")
    print(f"  survive={len(kept)}  lost={len(lost)}  new={len(new)}")

    def show(ids, title):
        print(f"\n  --- {title} ({len(ids)}) ---")
        if not ids:
            print("      (none)")
            return
        print(f"  {'id':<22} {'famA_raw':<12} {'famA_dec':<12} "
              f"{'pctA_raw':>9} {'pctA_dec':>9} {'pctB_raw':>9} {'pctB_dec':>9} "
              f"{'len':>5}")
        for i in ids:
            ra, da = results[A]["raw"].loc[i], results[A]["dec"].loc[i]
            rb, db = results[B]["raw"].loc[i], results[B]["dec"].loc[i]
            print(f"  {i:<22} {ra['best_family']:<12} {da['best_family']:<12} "
                  f"{ra['pct_within_family']:9.2f} {da['pct_within_family']:9.2f} "
                  f"{rb['pct_within_family']:9.2f} {db['pct_within_family']:9.2f} "
                  f"{frames[A]['cand_len'][i]:5d}")

    bhyb = set(results[A]["hyb"]["inside_envelope"][results[A]["hyb"]["inside_envelope"]].index) & \
           set(results[B]["hyb"]["inside_envelope"][results[B]["hyb"]["inside_envelope"]].index)
    print(f"  HYBRID both-inside       n={len(bhyb)}  "
          f"(vs raw: survive={len(braw & bhyb)} lost={len(braw - bhyb)} new={len(bhyb - braw)})")
    print("  HYBRID both-inside ids: " + ", ".join(sorted(bhyb)) if bhyb else "  (none)")
    show(kept, "SURVIVE deconfounding")
    show(lost, "LOST under deconfounding (exclusion no longer safe)")
    show(new, "NEW under deconfounding")

    # ================================================ Q3 family agreement ===
    for mode in ("raw", "dec", "hyb"):
        label = MODE_LABEL[mode]
        hdr(f"Q3. FAMILY ASSIGNMENT AGREEMENT across all candidates — {label}")
        fa = results[A][mode]["best_family"]
        fb = results[B][mode]["best_family"].reindex(fa.index)
        agree = fa == fb
        n = len(fa)
        print(f"  agreement {int(agree.sum())}/{n} = {100 * agree.mean():.1f}%")
        print(f"\n  {A} assignment counts:")
        print("   ", fa.value_counts().to_dict())
        print(f"  {B} assignment counts:")
        print("   ", fb.value_counts().to_dict())
        for f in fams:
            if int((fb == f).sum()) == 0:
                print(f"  NOTE: {B} assigns ZERO candidates to {f}")
            if int((fa == f).sum()) == 0:
                print(f"  NOTE: {A} assigns ZERO candidates to {f}")
        # disagreement vs reference family size
        size_dis = [fam_n[f] for f in fa[~agree]]
        size_agr = [fam_n[f] for f in fa[agree]]
        print(f"\n  Does disagreement track reference family size? "
              f"(size of the {A} call's family)")
        mwu("ref family size", size_dis, size_agr, "disagree", "agree")

    # ================================================ Q4 length bias in XOR =
    hdr("Q4. IS THE LENGTH BIAS IN THE ONE-SIDED (XOR) SET REMOVED?")
    idx = results[A]["raw"].index
    lens = pd.Series({i: frames[A]["cand_len"][i] for i in idx})
    q4 = {}
    for mode in ("raw", "dec", "hyb"):
        label = MODE_LABEL[mode]
        xor = results[f"{mode}_sets"]["xor"]
        oth = set(idx) - xor
        print(f"\n  {label}: XOR n={len(xor)}  other n={len(oth)}")
        q4[mode] = mwu("seq_len  XOR vs rest",
                       lens[sorted(xor)] if xor else [],
                       lens[sorted(oth)], "xor", "rest")
        oa = results[f"{mode}_sets"]["only_a"]
        ob = results[f"{mode}_sets"]["only_b"]
        if oa and ob:
            mwu("seq_len  A-only vs B-only",
                lens[sorted(oa)], lens[sorted(ob)], f"{A}", f"{B}")
        else:
            print(f"    (A-only n={len(oa)}, B-only n={len(ob)}: "
                  "one side empty, A-vs-B length test not defined)")

    # =========================================================== write out ==
    for tag in (A, B):
        for mode in ("raw", "dec", "hyb"):
            out = f"{a.out_prefix}_{tag}_{mode}.tsv"
            df = results[tag][mode].copy()
            df["seq_len"] = [frames[A]["cand_len"][i] for i in df.index]
            df.to_csv(out, sep="\t")
            print(f"[write] {out} ({len(df)} rows)")

    hdr("SUMMARY OF OVER-CORRECTION PROBES")
    for tag in (A, B):
        print(f"  {tag:<16} reference self-recovery "
              f"raw={results[tag]['rec_raw'] * 100:5.1f}%  "
              f"deconfounded={results[tag]['rec_dec'] * 100:5.1f}%  "
              f"Δ={100 * (results[tag]['rec_dec'] - results[tag]['rec_raw']):+.1f} pp")
    print(f"\n  reference rows in {os.path.basename(a.ref_labels)}: {len(ref_tbl)}; "
          f"resolved into novelty families: {sum(fam_n.values())}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
