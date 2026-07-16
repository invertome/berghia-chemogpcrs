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
    """RRA (Kolde et al. 2012, *Bioinformatics* 28:573): per-candidate corrected
    p that its normalized ranks are consistently smaller (more novel) than the
    uniform-null expectation.

    For a candidate observed in ``m`` models, sort its normalized ranks
    ``r_(1) <= ... <= r_(m)``. Under the null (ranks ~ Uniform(0,1]) the k-th
    order statistic (1-based) is Beta(k, m-k+1)-distributed, so
    ``beta.cdf(r_(k), k, m-k+1)`` is the probability of seeing a rank at least
    that small by chance. ``rho = min_k`` of those tail probabilities is the RRA
    statistic; a Bonferroni ``min(rho * m, 1)`` corrects for taking the minimum
    over m order statistics, yielding a p in (0,1] (every r_(k) >= 1/n > 0, so
    rho > 0). With the loop index ``k`` 0-based, the (k+1)-th order statistic
    uses Beta(k+1, m-k) — verified to keep both shape params strictly positive
    for every k in 0..m-1.
    """
    per_model = {m: _descending_percentile_ranks(dict(s)) for m, s in novelty.items()}
    out: Dict[str, float] = {}
    for cid in {c for s in novelty.values() for c in s}:
        rs = sorted(pm[cid] for pm in per_model.values() if cid in pm)
        m = len(rs)
        # k-th smallest of m uniforms (1-based k=index+1): Beta(k, m-k+1)
        rho = min(beta.cdf(rs[k], k + 1, m - k) for k in range(m))
        out[cid] = float(min(rho * m, 1.0))
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


def select_consensus(
    combiners: Mapping[str, Mapping[str, float]],
    combiner_famacc: Mapping[str, float],
    confounds: Mapping[str, Mapping[str, float]],
    best_single: Mapping[str, float],
    n_boot: int = 2000,
    seed: int = 0,
) -> Dict[str, object]:
    """Accept a combiner iff BOTH criteria beat the best single model:

      (a) robustness:  max|confound ρ| <= best_single["confound_rho"]
      (b) known-class: reference vote-famAcc >= best_single["famacc"]

    famAcc is checked first, so a noise combiner that decorrelates the confound
    (passes a) but destroys the supervised signal is rejected with reason
    ``"famacc"``. Among accepted combiners, the one with the best famAcc is
    selected; if none pass, the finding is "no fusion improves on the single
    best model" (``selected`` is None).
    """
    accepted: Dict[str, dict] = {}
    rejected: Dict[str, dict] = {}
    for name, nov in combiners.items():
        conf = _max_abs_confound(nov, confounds)
        fam = combiner_famacc[name]
        if fam < best_single["famacc"]:
            rejected[name] = {"reason": "famacc", "confound": conf, "famacc": fam}
        elif conf > best_single["confound_rho"]:
            rejected[name] = {"reason": "confound", "confound": conf, "famacc": fam}
        else:
            accepted[name] = {"confound": conf, "famacc": fam}
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
) -> Dict[str, Dict[str, object]]:
    """Emit the `mahalanobis_channel` contract from a consensus of per-model
    scores — a drop-in replacement for a single-model embedding channel.

    ``emb_novelty`` = ``-log10(RRA p)`` (high = novel) for ``combiner="rra"``, or
    the mean novelty percentile for ``"rank_average"``. ``emb_nonchemo_family``
    = plurality vote of the per-model nearest families. ``has_emb_data`` and
    ``emb_leakage_flag`` are the blanket True flags the contract requires.
    """
    fam = plurality_family(per_model_family)
    if combiner == "rra":
        p = robust_rank_aggregation(per_model_novelty)
        score = {c: float(-np.log10(max(p[c], 1e-300))) for c in p}
    else:
        score = rank_average_consensus(per_model_novelty)
    return {
        c: {
            "emb_novelty": score[c],
            "emb_nonchemo_family": fam.get(c),
            "has_emb_data": True,
            "emb_leakage_flag": True,
        }
        for c in score
    }


# --- CLI glue --------------------------------------------------------------
# Emitted per candidate; matches build_embedding_channel.MAHA_TSV_COLUMNS order
# so a consensus channel is a drop-in for a single-model one downstream.
CONSENSUS_TSV_COLUMNS: List[str] = [
    "id",
    "emb_nonchemo_family",
    "emb_novelty",
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
    rra = robust_rank_aggregation(per_model_novelty)
    ranks = rank_average_consensus(per_model_novelty)
    combiners = {
        "rra": {c: -float(np.log10(max(rra[c], 1e-300))) for c in rra},
        "rank_average": ranks,
    }
    sel = select_consensus(
        combiners=combiners,
        combiner_famacc={"rra": 0.85, "rank_average": 0.80},
        confounds=confounds,
        best_single={"confound_rho": 0.5, "famacc": 0.80},
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
        conf_str = ", ".join(
            f"{c}:{r['spearman']:+.3f}" for c, r in rep.items()
        )
        print(f"[fusion_consensus] {tag}: {len(df)} candidates cached -> {cache}"
              f"  (confounds {conf_str})")

    confounds = {
        "identity_to_nearest": identity_any,
        "seq_len": {i: float(v) for i, v in seq_len.items()},
    }
    resid = confound_residuals(per_model_novelty, confounds)
    gate = residual_independence(resid, threshold=args.indep_threshold)
    verdict = "GO (partly-independent confounds)" if gate["independent"] else \
        "STOP (common-mode confound — fusion unjustified)"
    print(f"[fusion_consensus] confound-independence gate: {verdict}")
    for (a, b), rho in gate["pairwise"].items():
        print(f"    residual ρ({a},{b}) = {rho:+.3f}")

    channel = build_consensus_channel(
        per_model_novelty, per_model_family, combiner=args.combiner
    )
    rows = [{"id": c, **channel[c]} for c in sorted(channel)]
    pd.DataFrame(rows, columns=CONSENSUS_TSV_COLUMNS).to_csv(
        args.out, sep="\t", index=False
    )
    print(f"[fusion_consensus] combiner={args.combiner} wrote {len(channel)} "
          f"candidates -> {args.out}")


if __name__ == "__main__":
    main()
