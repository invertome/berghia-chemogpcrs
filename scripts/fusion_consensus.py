#!/usr/bin/env python3
"""Fusion-consensus harness (bead cw3.16): combine per-model PLM novelty scores.

Score-space only. Pure functions (no torch/esm) so they unit-test without a GPU.
See docs/plans/2026-07-16-fusion-consensus-harness.md for the design + rationale.

The harness decides whether *combining* several PLM novelty scorers yields a
more robust, generalizable novelty channel than the single best model — with a
built-in confound-independence gate that can legitimately answer "no". Every
statistical primitive here is a pure function operating on plain
`{model: {candidate_id: value}}` mappings; the model set (and which models count
as "clean") is a parameter, nothing is hardcoded. IO glue (loading npz,
computing per-model novelty via `candidate_diagnostics`, writing the consensus
TSV) lives in `main()` and is imported lazily so the pure math stays lightweight.
"""
from __future__ import annotations

from collections import Counter
from itertools import combinations
from typing import Dict, Iterable, List, Mapping, Optional, Sequence

import numpy as np
from scipy.stats import beta, rankdata, spearmanr

NoveltyByModel = Mapping[str, Mapping[str, float]]


# ---------------------------------------------------------------------------
# Task 1 — consensus combiners (rank-average baseline + Robust Rank Aggregation)
# ---------------------------------------------------------------------------
def _descending_percentile_ranks(scores: Dict[str, float]) -> Dict[str, float]:
    """cid -> normalized rank in (0,1], smallest = most novel (highest score)."""
    ids = list(scores)
    r = rankdata([-scores[i] for i in ids], method="average")  # 1 = most novel
    n = len(ids)
    return {i: r[k] / n for k, i in enumerate(ids)}


def rank_average_consensus(novelty: NoveltyByModel) -> Dict[str, float]:
    """Mean cross-model *novelty* percentile (high = novel in more models)."""
    per_model = {m: _descending_percentile_ranks(dict(s)) for m, s in novelty.items()}
    out: Dict[str, float] = {}
    for cid in {c for s in novelty.values() for c in s}:
        pr = [1.0 - pm[cid] for pm in per_model.values() if cid in pm]  # novel = high
        out[cid] = float(np.mean(pr))
    return out


def robust_rank_aggregation(novelty: NoveltyByModel) -> Dict[str, float]:
    """RRA (Kolde et al. 2012, *Bioinformatics* 28:573): per-candidate p that
    its normalized ranks are consistently smaller (more novel) than the
    uniform-null expectation.

    For a candidate observed in ``m`` models, sort its normalized ranks
    ``r_(1) <= ... <= r_(m)``. Under the null (ranks ~ Uniform(0,1]) the k-th
    order statistic (1-based) is Beta(k, m-k+1)-distributed, so
    ``beta.cdf(r_(k), k, m-k+1)`` is the probability of seeing a rank at least
    that small by chance; ``rho = min_k`` of those tail probabilities is Kolde's
    RRA statistic. With the loop index ``k`` 0-based, the (k+1)-th order
    statistic uses Beta(k+1, m-k) — the SAME statistic
    :func:`rank_aggregation.rho_statistic` computes with a 1-based index
    (``betainc(j, m-j+1, r_(j))``), since ``betainc`` and ``beta.cdf`` are the
    same regularized incomplete beta and the 0-/1-based indices coincide.

    rho is then referred to its EXACT null distribution conditioned on ``m``
    (:func:`rank_aggregation.rho_null_pvalue`), NOT the Bonferroni bound
    ``min(rho * m, 1)`` this function used to return (bead sp4q / 8k8e). That
    bound CLIPS: every candidate with ``rho >= 1/m`` collapses onto a single
    p=1.0 tie block, flattening the whole low-novelty tail. The exact null is
    strictly increasing in rho, so it resolves that tail into distinct scores
    and never ties two candidates whose rank vectors differ — and it is the SAME
    estimator ``rank_aggregation`` applies to the downstream aggregator that
    rank-transforms this fusion output (``emb_novelty``), so the two are now
    consistent. LOWER = more consistently novel; p stays in (0,1] (every
    r_(k) >= 1/n > 0, so rho > 0).

    ``rho_null_pvalue`` is pure numpy/scipy but lives beside pandas-using code,
    so it is imported lazily (like the module's IO glue) to keep the pure-math
    top level lightweight.
    """
    from rank_aggregation import rho_null_pvalue  # single exact-null source (bead sp4q)
    per_model = {m: _descending_percentile_ranks(dict(s)) for m, s in novelty.items()}
    out: Dict[str, float] = {}
    for cid in {c for s in novelty.values() for c in s}:
        rs = sorted(pm[cid] for pm in per_model.values() if cid in pm)
        m = len(rs)
        # k-th smallest of m uniforms (1-based k=index+1): Beta(k, m-k+1)
        rho = min(beta.cdf(rs[k], k + 1, m - k) for k in range(m))
        out[cid] = float(rho_null_pvalue(rho, m))
    return out


# ---------------------------------------------------------------------------
# Task 2 — confound-independence gate
# ---------------------------------------------------------------------------
def confound_residuals(
    novelty: NoveltyByModel,
    confounds: Mapping[str, Mapping[str, float]],
) -> Dict[str, Dict[str, float]]:
    """Per model: residual of rank(novelty) after OLS on rank(each confound).

    Regressing on *ranks* makes the de-confounding monotone/robust (matches how
    the confound is measured with Spearman elsewhere). Only candidates with a
    value for every confound are used.
    """
    out: Dict[str, Dict[str, float]] = {}
    for model, scores in novelty.items():
        ids = [i for i in scores if all(i in confounds[c] for c in confounds)]
        y = rankdata([scores[i] for i in ids])
        X = np.column_stack(
            [np.ones(len(ids))]
            + [rankdata([confounds[c][i] for i in ids]) for c in confounds]
        )
        beta_hat, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta_hat
        out[model] = {i: float(resid[k]) for k, i in enumerate(ids)}
    return out


def residual_independence(
    residuals: Dict[str, Dict[str, float]],
    threshold: float = 0.7,
) -> Dict[str, object]:
    """Pairwise Spearman of the confound-removed residuals across models.

    Fusion only reduces error when model errors are partly *independent*. If
    every pairwise residual correlation is high (common-mode leakage/length
    bias), fusion is not justified — emit a STOP finding. ``independent`` is True
    when at least one pair falls below ``threshold`` (partly-independent
    confounds → GO).
    """
    models = list(residuals)
    pairs: Dict[tuple, float] = {}
    for a, b in combinations(models, 2):
        common = sorted(set(residuals[a]) & set(residuals[b]))
        rho, _ = spearmanr(
            [residuals[a][i] for i in common],
            [residuals[b][i] for i in common],
        )
        pairs[(a, b)] = float(rho)
    independent = any(abs(r) < threshold for r in pairs.values()) if pairs else True
    return {
        "pairwise": pairs,
        "independent": independent,
        "max_abs_rho": max((abs(r) for r in pairs.values()), default=0.0),
    }


# ---------------------------------------------------------------------------
# Task 3 — consensus family vote + reference LOO vote-famAcc
# ---------------------------------------------------------------------------
def _mode(labels: Iterable[str]) -> str:
    """Modal label with a deterministic alphabetical tie-break."""
    c = Counter(labels)
    top = max(c.values())
    return sorted(k for k, v in c.items() if v == top)[0]


def plurality_family(
    per_model_family: Mapping[str, Mapping[str, str]],
) -> Dict[str, str]:
    """Per candidate, the plurality-voted nearest family across models
    (alphabetical tie-break)."""
    out: Dict[str, str] = {}
    for cid in {c for f in per_model_family.values() for c in f}:
        votes = [f[cid] for f in per_model_family.values() if cid in f]
        out[cid] = _mode(votes)
    return out


def loo_vote_famacc(
    per_model_assignment: Mapping[str, Mapping[str, str]],
    true_family: Mapping[str, str],
) -> float:
    """Majority-vote family accuracy over held-out references.

    Applies the same per-model leave-one-out family assignment, majority-votes
    across models, and scores against the known reference family. This
    supervised anchor is what keeps the confound criterion from being gamed by
    injecting noise (which decorrelates the confound but destroys famAcc).
    """
    voted = plurality_family(per_model_assignment)
    ids = [i for i in voted if i in true_family]
    if not ids:
        return float("nan")
    return sum(voted[i] == true_family[i] for i in ids) / len(ids)


# ---------------------------------------------------------------------------
# Task 4 — dual selection criterion with bootstrap CIs
# ---------------------------------------------------------------------------
def bootstrap_ci_diff(
    values_a: Sequence[float],
    values_b: Sequence[float],
    stat,
    n: int = 2000,
    seed: int = 0,
    alpha: float = 0.05,
):
    """CI of ``stat(a) - stat(b)`` over paired bootstrap resamples.

    ``a`` and ``b`` share resampled indices so the difference respects the
    pairing (e.g. per-candidate confound ρ, or per-reference famAcc). Returns a
    ``(lo, hi)`` tuple at the ``(alpha/2, 1-alpha/2)`` quantiles.
    """
    rng = np.random.default_rng(seed)
    a, b = np.asarray(values_a), np.asarray(values_b)
    idx = np.arange(len(a))
    diffs = []
    for _ in range(n):
        s = rng.choice(idx, len(idx), replace=True)
        diffs.append(stat(a[s]) - stat(b[s]))
    lo, hi = np.quantile(diffs, [alpha / 2, 1 - alpha / 2])
    return float(lo), float(hi)


def _max_abs_confound(
    novelty_map: Mapping[str, float],
    confounds: Mapping[str, Mapping[str, float]],
) -> float:
    """max_confound |Spearman(novelty, confound)| over candidates scored on all
    confounds."""
    ids = [i for i in novelty_map if all(i in c for c in confounds.values())]
    y = [novelty_map[i] for i in ids]
    return max(
        abs(spearmanr(y, [c[i] for i in ids])[0]) for c in confounds.values()
    )


def _paired_confound_diff_ci(
    combiner_nov: Mapping[str, float],
    best_nov: Mapping[str, float],
    confounds: Mapping[str, Mapping[str, float]],
    n: int,
    seed: int,
    alpha: float,
):
    """CI of (combiner − best) max|confound ρ| over a paired candidate bootstrap.

    Each resample recomputes max_c|Spearman(novelty, confound_c)| for both the
    combiner and the best single model on the SAME resampled candidates, so the
    difference respects the pairing. Positive => combiner more confounded.
    """
    ids = [
        i for i in combiner_nov
        if i in best_nov and all(i in c for c in confounds.values())
    ]
    A = np.array([combiner_nov[i] for i in ids])
    B = np.array([best_nov[i] for i in ids])
    C = {name: np.array([vals[i] for i in ids]) for name, vals in confounds.items()}
    rng = np.random.default_rng(seed)
    idx = np.arange(len(ids))
    diffs = []
    for _ in range(n):
        s = rng.choice(idx, len(idx), replace=True)
        a = max(abs(spearmanr(A[s], C[name][s])[0]) for name in C)
        b = max(abs(spearmanr(B[s], C[name][s])[0]) for name in C)
        diffs.append(a - b)
    lo, hi = np.quantile(diffs, [alpha / 2, 1 - alpha / 2])
    return float(lo), float(hi)


def select_consensus(
    combiner_novelty: Mapping[str, Mapping[str, float]],
    combiner_ref_correct: Mapping[str, Mapping[str, int]],
    confounds: Mapping[str, Mapping[str, float]],
    best_single_novelty: Mapping[str, float],
    best_single_ref_correct: Mapping[str, int],
    n_boot: int = 2000,
    seed: int = 0,
    alpha: float = 0.05,
) -> Dict[str, object]:
    """Accept a combiner unless it is SIGNIFICANTLY worse than the best single
    model on either axis, judged by a paired bootstrap CI on the difference
    (so we don't chase noise-level ρ/famAcc deltas at large n):

      (a) known-class: reject if the CI of (combiner − best) reference
          vote-famAcc is entirely below 0 (robustly worse) → reason ``"famacc"``.
      (b) robustness:  reject if the CI of (combiner − best) max|confound ρ| is
          entirely above 0 (robustly more confounded) → reason ``"confound"``.

    famAcc is the supervised anchor and is checked first, so a noise combiner
    that decorrelates the confound but destroys family separability is rejected
    with reason ``"famacc"``. ``combiner_ref_correct`` / ``best_single_ref_correct``
    are per-reference 0/1 vote-correctness vectors (from the LOO family vote — the
    real vectors arrive with cw3.16.2; synthetic ones exercise it now). Among
    accepted combiners the highest point-famAcc is selected; if none pass, the
    finding is "no fusion improves on the single best model" (``selected`` None).
    Each verdict carries its point metrics and both CIs for auditability.
    """
    ref_ids = sorted(best_single_ref_correct)
    accepted: Dict[str, dict] = {}
    rejected: Dict[str, dict] = {}
    for name, nov in combiner_novelty.items():
        conf = _max_abs_confound(nov, confounds)
        common = [r for r in ref_ids if r in combiner_ref_correct[name]]
        comb_corr = np.array([combiner_ref_correct[name][r] for r in common])
        base_corr = np.array([best_single_ref_correct[r] for r in common])
        fam = float(np.mean(comb_corr)) if len(comb_corr) else float("nan")
        fam_lo, fam_hi = bootstrap_ci_diff(
            comb_corr, base_corr, np.mean, n=n_boot, seed=seed, alpha=alpha
        )
        conf_lo, conf_hi = _paired_confound_diff_ci(
            nov, best_single_novelty, confounds, n_boot, seed, alpha
        )
        # "confound"/"famacc" are the combiner's own POINT metrics; the "*_ci"
        # entries are CIs of the (combiner − best_single) DIFFERENCE, not of the
        # point metric — the decision below is on the difference CIs.
        info = {
            "confound": conf, "famacc": fam,
            "famacc_ci": (fam_lo, fam_hi), "confound_ci": (conf_lo, conf_hi),
        }
        if fam_hi < 0:                      # robustly worse famAcc
            rejected[name] = {"reason": "famacc", **info}
        elif conf_lo > 0:                   # robustly more confounded
            rejected[name] = {"reason": "confound", **info}
        else:
            accepted[name] = info
    selected = (
        max(accepted, key=lambda k: accepted[k]["famacc"]) if accepted else None
    )
    return {"selected": selected, "accepted": accepted, "rejected": rejected}


# ---------------------------------------------------------------------------
# Task 5 — contract-compatible consensus channel + CLI (IO glue)
# ---------------------------------------------------------------------------
def build_consensus_channel(
    per_model_novelty: NoveltyByModel,
    per_model_family: Mapping[str, Mapping[str, str]],
    combiner: str = "rra",
    deconfound: Optional[Mapping[str, Mapping[str, float]]] = None,
    residual_deconfound: Optional[Mapping[str, Mapping[str, float]]] = None,
) -> Dict[str, Dict[str, object]]:
    """Emit the `mahalanobis_channel` contract from a consensus of per-model
    scores — a drop-in replacement for a single-model embedding channel.

    ``emb_novelty`` = ``-log10(RRA p)`` (high = novel) for ``combiner="rra"``, or
    the mean novelty percentile for ``"rank_average"``. ``emb_nonchemo_family``
    = plurality vote of the per-model nearest families. ``has_emb_data`` and
    ``emb_leakage_flag`` are the blanket True flags the contract requires.

    ``deconfound`` (cw3.6): when given a ``{confound: {cid: value}}`` mapping,
    each model's novelty is residualized (rank-OLS, :func:`confound_residuals`)
    on those confounds BEFORE the combiner runs, so the emitted novelty is the
    length-deconfounded consensus. The RRA/rank combiners re-rank internally, so
    feeding residuals (higher residual = more novel than the confound predicts)
    is correct. Pass only ``seq_len`` here to strip the pervasive mean-pooling
    length artifact while KEEPING identity's added value. Candidates missing any
    deconfound value are dropped by the residualizer, so deconfound only on a
    fully-covered confound (seq_len covers every candidate).

    ``residual_deconfound`` (A1, v4bs.2): when given, ALSO emit a SECOND novelty
    ``emb_novelty_residual`` — the consensus of novelty additionally residualized
    on these confounds (e.g. ``tree_distance``), i.e. novelty *beyond*
    phylogenetic (and/or compositional) expectation. This is a DORMANT extra
    column/voter: ``emb_novelty`` (the production channel) is unchanged, so the
    production shortlist does not move; promoting the residual to a scored voter
    is a separate, gated decision. Candidates missing a residual confound value
    are dropped by the residualizer, so their ``emb_novelty_residual`` is None
    (partial coverage is expected — tree_distance only covers tree leaves).
    """
    fam = plurality_family(per_model_family)

    def _combine(confs: Optional[Mapping[str, Mapping[str, float]]]) -> Dict[str, float]:
        novelty = confound_residuals(per_model_novelty, confs) if confs else per_model_novelty
        if combiner == "rra":
            p = robust_rank_aggregation(novelty)
            return {c: float(-np.log10(max(p[c], 1e-300))) for c in p}
        return rank_average_consensus(novelty)

    score = _combine(deconfound)
    residual_score = _combine(residual_deconfound) if residual_deconfound else None

    out: Dict[str, Dict[str, object]] = {}
    for c in score:
        entry: Dict[str, object] = {
            "emb_novelty": score[c],
            "emb_nonchemo_family": fam.get(c),
            "has_emb_data": True,
            "emb_leakage_flag": True,
        }
        if residual_score is not None:
            entry["emb_novelty_residual"] = residual_score.get(c)  # None if dropped
        out[c] = entry
    return out


# --- CLI glue --------------------------------------------------------------
# Emitted per candidate; matches build_embedding_channel.MAHA_TSV_COLUMNS order
# so a consensus channel is a drop-in for a single-model one downstream.
CONSENSUS_TSV_COLUMNS: List[str] = [
    "id",
    "emb_nonchemo_family",
    "emb_novelty",
    "emb_novelty_residual",   # A1 (v4bs.2): dormant phylo/composition-residual novelty
    "has_emb_data",
    "emb_leakage_flag",
]


def _parse_model_spec(spec: str):
    """`tag:cand_npz:ref_npz[:identity_tsv]` -> (tag, cand_npz, ref_npz, identity_tsv)."""
    parts = spec.split(":")
    if len(parts) < 3:
        raise ValueError(
            f"--models entry {spec!r} must be tag:cand_npz:ref_npz[:identity_tsv]"
        )
    tag, cand, ref = parts[0], parts[1], parts[2]
    identity = parts[3] if len(parts) > 3 else ""
    return tag, cand, ref, identity


def _self_test() -> None:
    """Exercise gate -> combiners -> selection -> channel on synthetic in-memory
    inputs (no npz/GPU). Demonstrates the full decision procedure end-to-end."""
    rng = np.random.default_rng(0)
    ids = [f"c{i}" for i in range(60)]
    ident = {i: float(k) for k, i in enumerate(ids)}
    length = {i: float(rng.normal()) for i in ids}
    per_model_novelty = {
        m: {i: ident[i] * 0.05 + rng.normal(0, 5) for i in ids}
        for m in ("m1", "m2", "m3")
    }
    per_model_family = {
        m: {i: rng.choice(["pep", "amine", "opsin"]) for i in ids}
        for m in ("m1", "m2", "m3")
    }
    confounds = {"identity": ident, "length": length}

    resid = confound_residuals(per_model_novelty, confounds)
    gate = residual_independence(resid, threshold=0.7)
    rra = robust_rank_aggregation(per_model_novelty)  # exact-null RRA p (bead sp4q)
    ranks = rank_average_consensus(per_model_novelty)
    combiners = {
        "rra": {c: -float(np.log10(max(rra[c], 1e-300))) for c in rra},
        "rank_average": ranks,
    }
    refs = [f"r{i}" for i in range(50)]
    sel = select_consensus(
        combiner_novelty=combiners,
        combiner_ref_correct={
            "rra": {r: int(i < 43) for i, r in enumerate(refs)},
            "rank_average": {r: int(i < 41) for i, r in enumerate(refs)},
        },
        confounds=confounds,
        best_single_novelty=per_model_novelty["m1"],
        best_single_ref_correct={r: int(i < 40) for i, r in enumerate(refs)},
        n_boot=200,
        seed=0,
    )
    channel = build_consensus_channel(per_model_novelty, per_model_family)
    print(f"[self-test] gate independent={gate['independent']} "
          f"(max|resid ρ|={gate['max_abs_rho']:.3f})")
    print(f"[self-test] selected combiner={sel['selected']} "
          f"rejected={ {k: v['reason'] for k, v in sel['rejected'].items()} }")
    print(f"[self-test] channel emitted {len(channel)} candidates "
          f"(keys={sorted(next(iter(channel.values())))})")


def main(argv: Optional[Sequence[str]] = None) -> None:
    import argparse
    import os

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--models", nargs="+",
        help="one or more tag:cand_npz:ref_npz[:identity_tsv] specs",
    )
    parser.add_argument("--ref-labels", help="family/class-labeled reference TSV")
    parser.add_argument("--candidate-fasta", help="candidate FASTA (for seq_len)")
    parser.add_argument("--out", help="output consensus channel TSV")
    parser.add_argument(
        "--combiner", choices=["rra", "rank_average"], default="rra",
        help="channel combiner (default rra: -log10 RRA p)",
    )
    parser.add_argument(
        "--deconfound", default="",
        help="comma-separated confound names to residualize novelty on before "
             "combining (cw3.6 locked decision: 'seq_len' strips the mean-pooling "
             "length artifact while keeping identity's added value). Names: "
             "seq_len, identity_to_nearest. Empty = combine raw novelty.",
    )
    parser.add_argument(
        "--confound-tsv", action="append", default=[], metavar="NAME:PATH",
        help="extra confound source as NAME:PATH (a TSV whose first column is the "
             "candidate id and second column the value); repeatable. Registers "
             "NAME for use by --deconfound / --residual-confound. Producers: "
             "tree_distance_to_refs.py (tree_distance), composition_distance.py "
             "(composition).",
    )
    parser.add_argument(
        "--residual-confound", default="",
        help="comma-separated confound names to residualize a SECOND novelty on, "
             "emitted as the DORMANT emb_novelty_residual column (A1/v4bs.2: "
             "'tree_distance' = novelty beyond phylogenetic expectation; add "
             "'composition' for A2). The production emb_novelty is unchanged; the "
             "residual is not a scored voter until separately gated. Empty = off.",
    )
    parser.add_argument(
        "--diag-dir", default="results/ranking/diagnostics",
        help="where to cache per-model novelty TSVs",
    )
    parser.add_argument(
        "--indep-threshold", type=float, default=0.7,
        help="residual-correlation STOP/GO threshold (default 0.7)",
    )
    parser.add_argument(
        "--self-test", action="store_true",
        help="run the full decision procedure on synthetic in-memory inputs",
    )
    args = parser.parse_args(argv)

    if args.self_test:
        _self_test()
        return

    missing = [
        f for f, v in (
            ("--models", args.models), ("--ref-labels", args.ref_labels),
            ("--candidate-fasta", args.candidate_fasta), ("--out", args.out),
        ) if not v
    ]
    if missing:
        parser.error(f"the following are required (unless --self-test): {', '.join(missing)}")

    # Lazy IO imports so the pure math above stays torch/pandas-free for unit tests.
    import pandas as pd
    from embedding_evidence import load_embeddings
    from build_embedding_channel import load_ref_labels
    from embedding_candidate_diagnostics import (
        candidate_diagnostics, confound_report, _read_lengths, _read_identity,
    )

    ref_labels = load_ref_labels(args.ref_labels)
    seq_len = _read_lengths(args.candidate_fasta)
    os.makedirs(args.diag_dir, exist_ok=True)

    per_model_novelty: Dict[str, Dict[str, float]] = {}
    per_model_family: Dict[str, Dict[str, str]] = {}
    per_model_confound: Dict[str, Dict[str, float]] = {}
    identity_any: Dict[str, float] = {}

    for spec in args.models:
        tag, cand_npz, ref_npz, identity_tsv = _parse_model_spec(spec)
        cand = load_embeddings(cand_npz)
        ref = load_embeddings(ref_npz)
        identity = _read_identity(identity_tsv)
        identity_any.update(identity)
        df = candidate_diagnostics(cand, ref, ref_labels, seq_len, identity)
        cache = os.path.join(args.diag_dir, f"novelty_{tag}_PROD.tsv")
        out_df = df.copy()
        out_df.insert(0, "model", tag)
        out_df.to_csv(cache, sep="\t", index=False)
        per_model_novelty[tag] = dict(zip(df["candidate_id"], df["novelty"]))
        per_model_family[tag] = dict(zip(df["candidate_id"], df["nearest_family"]))
        rep = confound_report(df)
        per_model_confound[tag] = {c: float(r["spearman"]) for c, r in rep.items()}
        conf_str = ", ".join(
            f"{c}:{r['spearman']:+.3f}" for c, r in rep.items()
        )
        print(f"[fusion_consensus] {tag}: {len(df)} candidates cached -> {cache}"
              f"  (confounds {conf_str})")

    # Only gate on confounds with real coverage over the scored candidates. An
    # empty/low-coverage map (no identity TSV given, or mmseqs omits the no-hit
    # divergent candidates that matter most) would otherwise reduce the gate to a
    # meaningless nan verdict on an empty candidate set.
    cand_ids = {c for m in per_model_novelty.values() for c in m}
    raw_confounds = {
        "identity_to_nearest": identity_any,
        "seq_len": {i: float(v) for i, v in seq_len.items()},
    }
    # Register extra confound sources from --confound-tsv NAME:PATH (A1/A2:
    # tree_distance, composition). Each TSV's first column is the candidate id,
    # second the value. Registered here so they are usable by both the coverage
    # gate and --deconfound / --residual-confound.
    for spec in args.confound_tsv:
        if ":" not in spec:
            parser.error(f"--confound-tsv {spec!r} must be NAME:PATH")
        name, path = spec.split(":", 1)
        cdf = pd.read_csv(path, sep="\t")
        idcol, valcol = cdf.columns[0], cdf.columns[1]
        raw_confounds[name] = {str(k): float(v)
                               for k, v in zip(cdf[idcol], cdf[valcol])}
        print(f"[fusion_consensus] registered confound {name!r} "
              f"({len(raw_confounds[name])} candidates) from {path}")

    confounds: Dict[str, Dict[str, float]] = {}
    for name, cmap in raw_confounds.items():
        cov = len(cand_ids & set(cmap)) / len(cand_ids) if cand_ids else 0.0
        if cov >= 0.5:
            confounds[name] = cmap
        else:
            print(f"[fusion_consensus] WARNING: confound {name!r} coverage "
                  f"{cov:.0%} < 50% — excluded from the independence gate")

    if not confounds or len(per_model_novelty) < 2:
        print("[fusion_consensus] confound-independence gate SKIPPED "
              f"(usable confounds={list(confounds)}, models={len(per_model_novelty)})")
        gate = {"pairwise": {}, "independent": None, "max_abs_rho": float("nan")}
    else:
        resid = confound_residuals(per_model_novelty, confounds)
        gate = residual_independence(resid, threshold=args.indep_threshold)
        verdict = "GO (partly-independent confounds)" if gate["independent"] else \
            "STOP (common-mode confound — fusion unjustified)"
        print(f"[fusion_consensus] confound-independence gate: {verdict}")
        for (a, b), rho in gate["pairwise"].items():
            print(f"    residual ρ({a},{b}) = {rho:+.3f}")

    # cw3.6: residualize novelty on the requested confounds (seq_len) before the
    # combiner so the emitted emb_novelty is the length-deconfounded consensus.
    # Deconfound off the RAW confound maps (seq_len covers every candidate), not
    # the coverage-gated `confounds` used by the independence gate.
    deconfound_names = [s.strip() for s in args.deconfound.split(",") if s.strip()]
    deconfound_map = {n: raw_confounds[n] for n in deconfound_names if n in raw_confounds}
    unknown_dc = [n for n in deconfound_names if n not in raw_confounds]
    if unknown_dc:
        print(f"[fusion_consensus] WARNING: --deconfound name(s) {unknown_dc} not "
              f"in confounds {list(raw_confounds)} — ignored")
    if deconfound_map:
        print(f"[fusion_consensus] deconfounding novelty on {list(deconfound_map)} "
              f"before combining")

    # A1/A2 (v4bs.2/.3): a SECOND, dormant novelty residualized on the requested
    # confounds (e.g. tree_distance) -> emb_novelty_residual. Production emb_novelty
    # is unchanged, so the shortlist does not move; promotion to a scored voter is
    # a separate, gated decision.
    residual_names = [s.strip() for s in args.residual_confound.split(",") if s.strip()]
    residual_map = {n: raw_confounds[n] for n in residual_names if n in raw_confounds}
    unknown_res = [n for n in residual_names if n not in raw_confounds]
    if unknown_res:
        print(f"[fusion_consensus] WARNING: --residual-confound name(s) {unknown_res} "
              f"not in confounds {list(raw_confounds)} — ignored")
    if residual_map:
        print(f"[fusion_consensus] emitting emb_novelty_residual "
              f"(novelty residualized on {list(residual_map)})")

    channel = build_consensus_channel(
        per_model_novelty, per_model_family, combiner=args.combiner,
        deconfound=deconfound_map or None,
        residual_deconfound=residual_map or None,
    )
    rows = [{"id": c, **channel[c]} for c in sorted(channel)]
    pd.DataFrame(rows, columns=CONSENSUS_TSV_COLUMNS).to_csv(
        args.out, sep="\t", index=False
    )
    print(f"[fusion_consensus] combiner={args.combiner} wrote {len(channel)} "
          f"candidates -> {args.out}")

    # Persist the decision context (spec data-flow: "consensus TSV + summary").
    import json
    # cw3.6 audit (non-gating): Spearman of the EMITTED consensus emb_novelty vs
    # each confound. Verifies the deconfounding did its job — seq_len ρ should be
    # ~0 after --deconfound seq_len, while identity's added value may remain. This
    # is the confound arm of select_consensus recorded as an audit line only; the
    # famAcc arm (and any accept/reject gate) is deferred to cw3.16.2.
    emitted = {c: channel[c]["emb_novelty"] for c in channel}
    channel_confound: Dict[str, float] = {}
    for name, cmap in confounds.items():
        ids = [i for i in emitted if i in cmap]
        if len(ids) >= 3:
            rho, _ = spearmanr([emitted[i] for i in ids], [cmap[i] for i in ids])
            channel_confound[name] = float(rho)

    summary = {
        "models": list(per_model_novelty),
        "n_candidates": len(channel),
        "combiner": args.combiner,
        "deconfound": list(deconfound_map),
        "residual_confound": list(residual_map),
        "per_model_confound_spearman": per_model_confound,
        "consensus_channel_confound_spearman": channel_confound,
        "select_consensus_audit": {
            "status": "deferred",
            "note": "confound arm shown by consensus_channel_confound_spearman; "
                    "famAcc arm needs reference LOO vote-famAcc vectors (cw3.16.2). "
                    "Channel emitted per locked decision (scope A); not gated.",
        },
        "gate": {
            "independent": gate["independent"],
            "threshold": args.indep_threshold,
            "max_abs_resid_rho": gate["max_abs_rho"],
            "pairwise_resid_rho": {
                f"{a}|{b}": r for (a, b), r in gate["pairwise"].items()
            },
        },
    }
    summary_path = args.out.rsplit(".", 1)[0] + ".summary.json"
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[fusion_consensus] summary -> {summary_path}")


if __name__ == "__main__":
    main()
