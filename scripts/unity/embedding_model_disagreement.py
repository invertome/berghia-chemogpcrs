#!/usr/bin/env python3
"""embedding_model_disagreement.py — characterize where the two consensus
embedding models DISAGREE about non-chemoreceptor family membership.

The locked production scorer (EMB_SCORER=consensus) reports only the agreeing
subset. That subset is, by construction, the set on which the two models say the
same thing; the disagreements are unexamined.

INTERPRETIVE FRAME (must survive into any report). The pipeline is SUBTRACTIVE.
There is no positive control for molluscan odorant receptors and there cannot be
one. The method confidently classifies every OTHER molluscan GPCR family; what
remains unclassified is the enriched residual. A confident EXCLUSION is therefore
the load-bearing output, and "novel / unassigned" is a residual, NEVER a positive
claim of chemoreceptor identity.

Consequence for disagreement: a candidate that one model places inside a known
family envelope and the other does not is a candidate whose EXCLUSION IS NOT
SAFE. It stays in the chemoreceptor pool. Disagreement is not evidence of
chemoreceptor identity; it is absence of a safe exclusion.

Reads only the existing per-model diagnostics written by
embedding_family_assignment.py and the novelty TSVs. Computes nothing new about
the embedding geometry — this is a readout of the production artifacts.
"""
from __future__ import annotations

import argparse
import itertools
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

GENE_RE = re.compile(r"t\d+$")


def gene_of(tid: str) -> str:
    """BersteEVm007379t1 -> BersteEVm007379 (transcript suffix stripped)."""
    return GENE_RE.sub("", tid)


def load_pair(diag_dir: str, tag: str) -> pd.DataFrame:
    fa = pd.read_csv(os.path.join(diag_dir, f"family_assignment_{tag}.tsv"), sep="\t")
    # family_assignment now carries its own seq_len; drop it so the merge below
    # does not produce seq_len_x/seq_len_y. The novelty TSV stays the single
    # source, and the two are cross-checked for equality further down.
    fa_len = dict(zip(fa["id"], fa["seq_len"])) if "seq_len" in fa.columns else {}
    fa = fa.drop(columns=[c for c in ("seq_len",) if c in fa.columns])
    nv = pd.read_csv(os.path.join(diag_dir, f"novelty_{tag}_PROD.tsv"), sep="\t")
    nv = nv.rename(columns={"candidate_id": "id"})
    keep_nv = ["id", "novelty", "nearest_family", "seq_len", "identity_to_nearest"]
    df = fa.merge(nv[keep_nv], on="id", how="outer", indicator=True)
    if (df["_merge"] != "both").any():
        raise SystemExit(
            f"FATAL[{tag}]: family_assignment / novelty id sets differ "
            f"({(df['_merge'] != 'both').sum()} non-matching)"
        )
    if fa_len:
        bad = int(sum(int(fa_len[i]) != int(l)
                      for i, l in zip(df["id"], df["seq_len"])))
        if bad:
            raise SystemExit(
                f"FATAL[{tag}]: seq_len disagrees between family_assignment and "
                f"the novelty TSV on {bad} rows"
            )
    # Cross-check: the two producers must agree on the RAW argmin family.
    # Under the producer's `--deconfound envelope` default (production
    # semantics) the family call is the RAW argmin and only the envelope
    # percentile is deconfounded, so `best_family` == `raw_argmin_family` ==
    # the novelty TSV's `nearest_family`. `raw_argmin_family` is checked rather
    # than `best_family` so that this guard keeps testing GEOMETRY DRIFT
    # specifically, and would not start firing on the deconfounding itself if
    # the family call were ever changed.
    raw_col = "raw_argmin_family" if "raw_argmin_family" in df.columns else None
    if raw_col is None:
        raise SystemExit(
            f"FATAL[{tag}]: family_assignment_{tag}.tsv has no "
            "`raw_argmin_family` column — it was produced by the OLD raw "
            "embedding_family_assignment.py. Re-run run_family_assignment.sh; "
            "the raw numbers are not the ones production uses."
        )
    mismatch = int((df[raw_col] != df["nearest_family"]).sum())
    if mismatch:
        raise SystemExit(
            f"FATAL[{tag}]: raw_argmin_family != nearest_family for {mismatch} "
            "rows; the two diagnostics were built from different geometry"
        )
    df = df.drop(columns=["_merge", "nearest_family"])
    df["novelty_rank"] = df["novelty"].rank(ascending=False, method="average")
    return df


def hdr(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def mwu(name: str, a: np.ndarray, b: np.ndarray, la: str, lb: str) -> None:
    """Mann-Whitney U with medians and rank-biserial effect size."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3:
        print(f"  {name:<26} SKIPPED (n={len(a)}/{len(b)}, too few to test)")
        return
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    rbc = 2.0 * u / (len(a) * len(b)) - 1.0
    print(
        f"  {name:<26} median {la}={np.median(a):10.3f}  {lb}={np.median(b):10.3f}"
        f"   U={u:9.1f}  p={p:.3e}  rank-biserial={rbc:+.3f}"
    )


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--diag-dir", required=True)
    p.add_argument("--ref-labels", required=True)
    p.add_argument("--model-a", default="proteinclip3b")
    p.add_argument("--model-b", default="protrek")
    p.add_argument("--out-tsv", required=True)
    a = p.parse_args()

    A, B = a.model_a, a.model_b
    da = load_pair(a.diag_dir, A)
    db = load_pair(a.diag_dir, B)
    print(f"[disagree] {A}: {len(da)} candidates   {B}: {len(db)} candidates")

    m = da.merge(db, on="id", how="inner", suffixes=(f"_{A}", f"_{B}"))
    if len(m) != len(da) or len(m) != len(db):
        raise SystemExit(
            f"FATAL: candidate id sets differ: {A}={len(da)} {B}={len(db)} "
            f"intersection={len(m)}"
        )
    # seq_len / identity_to_nearest are model-independent; assert and collapse.
    for col in ("seq_len", "identity_to_nearest"):
        ca, cb = f"{col}_{A}", f"{col}_{B}"
        if not m[ca].equals(m[cb]):
            raise SystemExit(f"FATAL: {col} differs between models; not a shared covariate")
        m[col] = m[ca]
        m = m.drop(columns=[ca, cb])
    N = len(m)
    print(f"[disagree] merged on id: {N} candidates, 1:1")

    fam_a, fam_b = f"best_family_{A}", f"best_family_{B}"
    in_a, in_b = f"inside_envelope_{A}", f"inside_envelope_{B}"
    m[in_a] = m[in_a].astype(bool)
    m[in_b] = m[in_b].astype(bool)
    m["family_agree"] = m[fam_a] == m[fam_b]
    m["gene"] = m["id"].map(gene_of)

    # Reference family sizes (the labelled anchor set the envelopes are built from).
    ref = pd.read_csv(a.ref_labels, sep="\t")
    fam_n = ref["family"].value_counts().to_dict()
    print(f"[disagree] reference rows={len(ref)}  families={len(fam_n)}")

    # ------------------------------------------------------------------ 1
    hdr("1. AGREEMENT MATRIX over all %d candidates (nearest-family argmin)" % N)
    n_agree = int(m["family_agree"].sum())
    print(f"\nFamily agreement (argmin, envelope-independent): "
          f"{n_agree}/{N} = {100.0 * n_agree / N:.1f}%")
    print(f"Family disagreement:                              "
          f"{N - n_agree}/{N} = {100.0 * (N - n_agree) / N:.1f}%")

    ct = pd.crosstab(m[fam_a], m[fam_b], dropna=False)
    print(f"\nCross-tabulation  rows={A} (A)  cols={B} (B):")
    print(ct.to_string())

    print("\nPer-model marginal family assignment:")
    marg = pd.DataFrame({
        A: m[fam_a].value_counts(),
        B: m[fam_b].value_counts(),
    }).fillna(0).astype(int)
    marg["ref_n"] = pd.Series(fam_n)
    print(marg.to_string())

    print("\nMost-confused UNORDERED family pairs (A!=B), by count:")
    dis = m[~m["family_agree"]]
    pairs: dict = {}
    for _, r in dis.iterrows():
        key = tuple(sorted((r[fam_a], r[fam_b])))
        pairs[key] = pairs.get(key, 0) + 1
    for (f1, f2), c in sorted(pairs.items(), key=lambda kv: -kv[1]):
        print(f"  {f1:<22} <-> {f2:<22} {c:4d}  "
              f"({100.0 * c / max(N - n_agree, 1):5.1f}% of disagreements)")

    print("\nPer-family agreement rate (denominator = model A's assignment):")
    for f in sorted(m[fam_a].unique()):
        sub = m[m[fam_a] == f]
        print(f"  {f:<22} {int(sub['family_agree'].sum()):4d}/{len(sub):<4d} "
              f"= {100.0 * sub['family_agree'].mean():5.1f}%")

    # ------------------------------------------------------------------ 2
    hdr("2. ENVELOPE-MEMBERSHIP DISAGREEMENT (one model excludes, the other does not)")
    both_in = m[m[in_a] & m[in_b]]
    a_only = m[m[in_a] & ~m[in_b]]
    b_only = m[~m[in_a] & m[in_b]]
    neither = m[~m[in_a] & ~m[in_b]]
    print(f"\n  inside BOTH envelopes            : {len(both_in):4d}")
    print(f"  inside {A:<24} ONLY : {len(a_only):4d}")
    print(f"  inside {B:<24} ONLY : {len(b_only):4d}")
    print(f"  inside NEITHER (unassigned)      : {len(neither):4d}")
    print(f"  ---------------------------------------")
    print(f"  total                            : {N:4d}")
    print(f"\n  union inside >=1 envelope        : {len(both_in) + len(a_only) + len(b_only):4d}")
    print(f"  ENVELOPE DISAGREEMENT (XOR)      : {len(a_only) + len(b_only):4d}"
          f"  ({100.0 * (len(a_only) + len(b_only)) / N:.1f}% of all candidates)")

    tab = np.array([[len(both_in), len(a_only)], [len(b_only), len(neither)]])
    if tab.min() >= 0 and tab.sum() == N:
        _, pmc = stats.fisher_exact(tab) if tab.size == 4 else (None, np.nan)
        print(f"\n  2x2 envelope-membership table (Fisher p={pmc:.3e}) -- tests whether the "
              f"two\n  models' inside/outside calls are associated at all:")
        print(f"    {'':<14}{'B inside':>10}{'B outside':>11}")
        print(f"    {'A inside':<14}{len(both_in):>10}{len(a_only):>11}")
        print(f"    {'A outside':<14}{len(b_only):>10}{len(neither):>11}")

    for label, sub, inside_model, outside_model in (
        (f"INSIDE {A} ONLY  ({A} would exclude it; {B} would not)", a_only, A, B),
        (f"INSIDE {B} ONLY  ({B} would exclude it; {A} would not)", b_only, B, A),
    ):
        print(f"\n--- {label} : n={len(sub)} ---")
        if not len(sub):
            print("  (none)")
            continue
        print(f"  {'id':<24}{'inside-model call':<22}{'pct':>7}{'dist':>10}"
              f"   | {'other-model nearest':<22}{'pct':>7}{'dist':>10}")
        for _, r in sub.sort_values(f"pct_within_family_{inside_model}").iterrows():
            print(f"  {r['id']:<24}{r[f'best_family_{inside_model}']:<22}"
                  f"{r[f'pct_within_family_{inside_model}']:7.2f}"
                  f"{r[f'best_dist_raw_{inside_model}']:10.1f}"
                  f"   | {r[f'best_family_{outside_model}']:<22}"
                  f"{r[f'pct_within_family_{outside_model}']:7.2f}"
                  f"{r[f'best_dist_raw_{outside_model}']:10.1f}")
        print(f"  families ({inside_model} call): "
              f"{sub[f'best_family_{inside_model}'].value_counts().to_dict()}")
        print(f"  of these, the other model's ARGMIN family still matches: "
              f"{int(sub['family_agree'].sum())}/{len(sub)}")

    # ------------------------------------------------------------------ 3
    hdr("3. FAMILY DISAGREEMENT among candidates inside BOTH envelopes")
    both_dis = both_in[~both_in["family_agree"]]
    both_agr = both_in[both_in["family_agree"]]
    print(f"\n  inside both AND same family (the reported agreeing core): {len(both_agr)}"
          f"  [{both_agr['gene'].nunique()} genes]")
    print(f"  inside both BUT different family                        : {len(both_dis)}")
    if len(both_dis):
        print(f"\n  {'id':<24}{A + ' call':<24}{'pct':>7}{'dist':>10}"
              f"   | {B + ' call':<24}{'pct':>7}{'dist':>10}")
        for _, r in both_dis.iterrows():
            print(f"  {r['id']:<24}{r[fam_a]:<24}{r[f'pct_within_family_{A}']:7.2f}"
                  f"{r[f'best_dist_raw_{A}']:10.1f}"
                  f"   | {r[fam_b]:<24}{r[f'pct_within_family_{B}']:7.2f}"
                  f"{r[f'best_dist_raw_{B}']:10.1f}")
    else:
        print("  (none)")

    # The agreeing core, for provenance against previously reported numbers.
    union_in = m[m[in_a] | m[in_b]]
    union_agr = union_in[union_in["family_agree"]]
    print(f"\n  PROVENANCE of previously reported counts:")
    print(f"    inside >=1 envelope AND family-agreeing : {len(union_agr)} transcripts / "
          f"{union_agr['gene'].nunique()} genes")
    print(f"      breakdown: {union_agr[fam_a].value_counts().to_dict()}")
    print(f"    inside BOTH envelopes (Tier A)          : {len(both_in)} transcripts / "
          f"{both_in['gene'].nunique()} genes")

    # ------------------------------------------------------------------ 4
    hdr("4. NOVELTY-RANK DISAGREEMENT")
    na, nb = f"novelty_{A}", f"novelty_{B}"
    ra, rb = f"novelty_rank_{A}", f"novelty_rank_{B}"
    rho, prho = stats.spearmanr(m[na], m[nb])
    tau, ptau = stats.kendalltau(m[na], m[nb])
    pear, ppear = stats.pearsonr(m[na], m[nb])
    print(f"\n  Spearman rho = {rho:+.4f}  (p={prho:.3e})   n={N}")
    print(f"  Kendall  tau = {tau:+.4f}  (p={ptau:.3e})")
    print(f"  Pearson    r = {pear:+.4f}  (p={ppear:.3e})   [raw novelty is not "
          f"cross-model comparable in scale; rank is the meaningful axis]")

    m["rank_diff"] = m[ra] - m[rb]
    m["abs_rank_diff"] = m["rank_diff"].abs()
    print(f"\n  |rank difference| distribution over {N} candidates:")
    for q in (50, 75, 90, 95, 99, 100):
        print(f"    p{q:<3d} = {np.percentile(m['abs_rank_diff'], q):7.1f} ranks")
    print(f"    mean = {m['abs_rank_diff'].mean():.1f} ranks")

    # Top-k concordance: does the head of the ranking survive a model swap?
    for k in (10, 25, 50, 100):
        ta = set(m.nsmallest(k, ra)["id"])
        tb = set(m.nsmallest(k, rb)["id"])
        print(f"  top-{k:<4d} most-novel overlap: {len(ta & tb):3d}/{k}"
              f"  (Jaccard {len(ta & tb) / len(ta | tb):.3f})")

    print(f"\n--- 20 LARGEST rank divergences (highly novel to one model, "
          f"unremarkable to the other) ---")
    print(f"  {'id':<24}{'rank_' + A:>12}{'rank_' + B:>12}{'|diff|':>9}"
          f"  {'fam_' + A:<20}{'fam_' + B:<20}{'in_A':>6}{'in_B':>6}")
    for _, r in m.nlargest(20, "abs_rank_diff").iterrows():
        print(f"  {r['id']:<24}{r[ra]:12.0f}{r[rb]:12.0f}{r['abs_rank_diff']:9.0f}"
              f"  {r[fam_a]:<20}{r[fam_b]:<20}"
              f"{str(bool(r[in_a])):>6}{str(bool(r[in_b])):>6}")

    print(f"\n  Directional split among the 50 largest divergences:")
    top50 = m.nlargest(50, "abs_rank_diff")
    more_novel_a = int((top50["rank_diff"] < 0).sum())
    print(f"    ranked MORE novel by {A}: {more_novel_a}/50")
    print(f"    ranked MORE novel by {B}: {50 - more_novel_a}/50")

    # ----------------------------------------------------------------- 4b
    hdr("4b. DOES THE RANK DISAGREEMENT SURVIVE LENGTH-DECONFOUNDING?")
    print("\n  novelty_{tag}_PROD.tsv caches RAW per-model novelty; production\n"
          "  (fusion_consensus.py --deconfound seq_len) residualizes rank(novelty) on\n"
          "  rank(seq_len) BEFORE combining. Reproducing that residualization here tests\n"
          "  whether the disagreement above is a length artifact or survives it.")
    from fusion_consensus import confound_residuals  # noqa: E402

    nov_by_model = {
        A: dict(zip(m["id"], m[na])),
        B: dict(zip(m["id"], m[nb])),
    }
    confs = {"seq_len": dict(zip(m["id"], m["seq_len"].astype(float)))}
    resid = confound_residuals(nov_by_model, confs)
    m[f"resid_{A}"] = m["id"].map(resid[A])
    m[f"resid_{B}"] = m["id"].map(resid[B])
    cov = int(m[[f"resid_{A}", f"resid_{B}"]].notna().all(axis=1).sum())
    print(f"\n  candidates with a residual under both models: {cov}/{N}")
    m[f"resid_rank_{A}"] = m[f"resid_{A}"].rank(ascending=False, method="average")
    m[f"resid_rank_{B}"] = m[f"resid_{B}"].rank(ascending=False, method="average")
    rr, prr = stats.spearmanr(m[f"resid_{A}"], m[f"resid_{B}"])
    print(f"  Spearman rho (length-deconfounded) = {rr:+.4f} (p={prr:.3e})")
    print(f"  Spearman rho (raw, from section 4) = {rho:+.4f}")
    print(f"  change = {rr - rho:+.4f}")
    m["resid_abs_rank_diff"] = (m[f"resid_rank_{A}"] - m[f"resid_rank_{B}"]).abs()
    print(f"\n  |rank difference| after deconfounding: median="
          f"{np.median(m['resid_abs_rank_diff']):.1f}  mean="
          f"{m['resid_abs_rank_diff'].mean():.1f}"
          f"   (raw: median={np.median(m['abs_rank_diff']):.1f} "
          f"mean={m['abs_rank_diff'].mean():.1f})")
    for k in (10, 25, 50, 100):
        ta = set(m.nsmallest(k, f"resid_rank_{A}")["id"])
        tb = set(m.nsmallest(k, f"resid_rank_{B}")["id"])
        raw_a = set(m.nsmallest(k, ra)["id"])
        raw_b = set(m.nsmallest(k, rb)["id"])
        print(f"  top-{k:<4d} overlap deconfounded: {len(ta & tb):3d}/{k}"
              f"   (raw: {len(raw_a & raw_b):3d}/{k})")

    # ------------------------------------------------------------------ 5
    hdr("5. IS DISAGREEMENT STRUCTURED OR NOISE?")
    m["fam_n_a"] = m[fam_a].map(fam_n)
    m["fam_n_b"] = m[fam_b].map(fam_n)
    m["min_pct"] = m[[f"pct_within_family_{A}", f"pct_within_family_{B}"]].min(axis=1)
    m["env_xor"] = m[in_a] ^ m[in_b]

    print("\n(a) FAMILY disagreement (argmin differs) vs covariates "
          "[Mann-Whitney, two-sided]")
    d1 = m[~m["family_agree"]]
    d0 = m[m["family_agree"]]
    print(f"    n_disagree={len(d1)}  n_agree={len(d0)}")
    mwu("seq_len", d1["seq_len"], d0["seq_len"], "dis", "agr")
    mwu("identity_to_nearest", d1["identity_to_nearest"], d0["identity_to_nearest"], "dis", "agr")
    mwu(f"best_dist_raw ({A})", d1[f"best_dist_raw_{A}"], d0[f"best_dist_raw_{A}"], "dis", "agr")
    mwu(f"best_dist_raw ({B})", d1[f"best_dist_raw_{B}"], d0[f"best_dist_raw_{B}"], "dis", "agr")
    mwu(f"pct_within_family ({A})", d1[f"pct_within_family_{A}"],
        d0[f"pct_within_family_{A}"], "dis", "agr")
    mwu(f"margin_raw ({A})", d1[f"margin_raw_{A}"], d0[f"margin_raw_{A}"], "dis", "agr")
    mwu(f"margin_raw ({B})", d1[f"margin_raw_{B}"], d0[f"margin_raw_{B}"], "dis", "agr")
    mwu(f"ref family size ({A} call)", d1["fam_n_a"], d0["fam_n_a"], "dis", "agr")

    print("\n(b) ENVELOPE disagreement (XOR: inside for exactly one model) vs covariates")
    e1 = m[m["env_xor"]]
    e0 = m[~m["env_xor"]]
    print(f"    n_xor={len(e1)}  n_non_xor={len(e0)}")
    mwu("seq_len", e1["seq_len"], e0["seq_len"], "xor", "oth")
    mwu("identity_to_nearest", e1["identity_to_nearest"], e0["identity_to_nearest"], "xor", "oth")
    mwu("min pct_within_family", e1["min_pct"], e0["min_pct"], "xor", "oth")

    print("\n(c) CONTINUOUS rank divergence vs covariates [Spearman]")
    for col in ("seq_len", "identity_to_nearest", "min_pct",
                f"best_dist_raw_{A}", f"best_dist_raw_{B}",
                f"margin_raw_{A}", f"margin_raw_{B}", "fam_n_a"):
        v = m[col].astype(float)
        ok = np.isfinite(v)
        r_, p_ = stats.spearmanr(m.loc[ok, "abs_rank_diff"], v[ok])
        print(f"  |rank_diff| ~ {col:<24} rho={r_:+.4f}  p={p_:.3e}  n={int(ok.sum())}")

    print("\n(d) Distance-to-reference stratification: is disagreement just "
          "'far from everything'?")
    print("    Candidates binned by min_pct (percentile within the nearest family's "
          "own spread;\n    100 = further than EVERY reference member of that family).")
    bins = [(0, 50), (50, 90), (90, 95), (95, 99.999), (99.999, 100.001)]
    print(f"    {'min_pct bin':<18}{'n':>6}{'fam_disagree':>14}{'rate':>8}"
          f"{'env_xor':>10}{'rate':>8}{'median|rankdiff|':>18}")
    for lo, hi in bins:
        sub = m[(m["min_pct"] >= lo) & (m["min_pct"] < hi)]
        if not len(sub):
            continue
        lbl = f"[{lo:g}, {hi:g})" if hi < 100 else "== 100 (beyond all)"
        print(f"    {lbl:<18}{len(sub):>6}{int((~sub['family_agree']).sum()):>14}"
              f"{100.0 * (~sub['family_agree']).mean():7.1f}%"
              f"{int(sub['env_xor'].sum()):>10}{100.0 * sub['env_xor'].mean():7.1f}%"
              f"{np.median(sub['abs_rank_diff']):>18.1f}")

    print("\n(e) leakage flag")
    print("    emb_leakage_flag is hardcoded True for every candidate in "
          "embedding_evidence.py\n    (mahalanobis_channel / embedding_channel). It is a "
          "CONSTANT, not a variable, so it\n    CANNOT be correlated with disagreement. "
          "No test is possible from available artifacts.")

    # ------------------------------------------------------------------ out
    out_cols = ["id", "gene", "seq_len", "identity_to_nearest",
                f"resid_{A}", f"resid_{B}", f"resid_rank_{A}", f"resid_rank_{B}",
                "resid_abs_rank_diff",
                fam_a, fam_b, "family_agree",
                f"pct_within_family_{A}", f"pct_within_family_{B}",
                f"best_dist_raw_{A}", f"best_dist_raw_{B}",
                f"margin_raw_{A}", f"margin_raw_{B}",
                in_a, in_b, "env_xor",
                na, nb, ra, rb, "rank_diff", "abs_rank_diff", "min_pct"]
    m.sort_values("abs_rank_diff", ascending=False)[out_cols].to_csv(
        a.out_tsv, sep="\t", index=False)
    print(f"\n[disagree] wrote {a.out_tsv} ({len(m)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
