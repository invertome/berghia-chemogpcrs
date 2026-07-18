#!/usr/bin/env python3
"""embedding_nulls.py — A4 residue-shuffle embedding null (epic v4bs, task v4bs.5).

Turns a candidate's raw embedding-novelty into a per-candidate p-value against a
NULL in which the sequence's residues are shuffled. If novelty reflects genuine
evolutionary/functional signal (not mere amino-acid composition), re-embedding a
residue-shuffled copy should COLLAPSE novelty toward the reference centroid, so
the observed novelty sits far in the upper tail of the shuffled-null distribution.

Scope: this module is the LOCAL, pure, testable harness only. It (1) GENERATES the
shuffled sequences to be re-embedded elsewhere (GPU/HPC) and (2) CONSUMES the
already-computed null novelties to produce a right-tail p-value and a collapse
diagnostic. It never embeds anything and imports no ranking/embedding code.

The empirical p uses the never-zero (b+1)/(m+1) rule (Phipson & Smyth 2010),
matching scripts/permutation_null.py's convention.
"""
from __future__ import annotations

from typing import Dict, Sequence

import numpy as np


def shuffle_residues(seq: str, frac: float, rng: "np.random.Generator") -> str:
    """Permute the residues among a random ``frac`` of positions; leave the rest.

    A fraction ``frac`` of positions is chosen without replacement, and their
    residues are permuted among ONLY those positions (a bijection over the
    selected set). Length and the residue multiset are therefore preserved for
    any ``frac``. ``frac=0.0`` returns ``seq`` unchanged (identity); ``frac=1.0``
    permutes every position. Randomness flows exclusively through the passed-in
    numpy ``Generator`` ``rng`` (no global state), so output is seed-deterministic.

    Note: because a permutation can map elements to their own positions, the
    output MAY equal the input with low probability (e.g. very short or
    low-complexity sequences); the length/multiset invariants always hold.
    """
    if frac == 0.0:
        return seq
    L = len(seq)
    n_sel = int(round(frac * L))
    if n_sel < 2:
        return seq  # 0 or 1 position selected -> permutation is the identity
    sel = rng.choice(L, size=n_sel, replace=False)
    perm = rng.permutation(n_sel)
    chars = list(seq)
    for dst, src in zip(sel, sel[perm]):
        chars[dst] = seq[src]
    return "".join(chars)


def empirical_pvalue(observed: float, null_values: Sequence[float]) -> float:
    """Never-zero right-tail empirical p: ``(b + 1) / (m + 1)`` with
    ``b = #{v in null_values : v >= observed}`` and ``m = len(null_values)``.

    Higher novelty is more extreme, so this is a right-tail test. The +1 in
    numerator and denominator is the Phipson & Smyth (2010) correction that
    forbids an exact-zero p-value, matching scripts/permutation_null.py's
    convention (re-derived here to avoid importing that module's pipeline deps).
    """
    m = len(null_values)
    b = int(sum(1 for v in null_values if v >= observed))
    return (b + 1) / (m + 1)


def novelty_collapse(observed: float, null_values: Sequence[float]) -> Dict[str, float]:
    """Diagnostic for whether observed novelty stands above shuffle noise.

    Returns ``null_mean``, ``null_std`` (population std), ``excess`` =
    ``(observed - null_mean) / null_std`` (a z-like score of how many null std
    the observed novelty exceeds the shuffled-null mean), and ``collapsed`` =
    ``observed <= null_mean`` (the shuffled null reproduces or exceeds the
    observed novelty, so the signal is NOT above shuffle noise).

    Zero-variance guard: when ``null_std == 0`` the z-score is undefined, so
    ``excess`` is reported as ``0.0`` when ``observed <= null_mean`` (no excess)
    and ``+inf`` when ``observed > null_mean`` (an unbounded departure from a
    point-mass null); ``collapsed`` is unaffected.
    """
    arr = np.asarray(list(null_values), dtype=float)
    null_mean = float(arr.mean())
    null_std = float(arr.std())
    collapsed = observed <= null_mean
    if null_std == 0.0:
        excess = 0.0 if collapsed else float("inf")
    else:
        excess = (observed - null_mean) / null_std
    return {
        "null_mean": null_mean,
        "null_std": null_std,
        "excess": float(excess),
        "collapsed": bool(collapsed),
    }


def main(argv=None):
    """CLI: p-value + collapse diagnostic from precomputed observed/null novelties.

    ``--observed-tsv`` has columns ``candidate_id``, ``novelty``. Each
    ``--null-tsv`` file is one shuffle-replicate's per-candidate novelty (same
    schema); pass one or more. For every candidate in the observed TSV, the null
    novelties are gathered across the replicate files and reduced to an empirical
    p-value and a collapse diagnostic. Writes a TSV with columns
    ``id, novelty, null_mean, p_value, collapsed``.
    """
    import argparse

    import pandas as pd

    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--observed-tsv", required=True)
    ap.add_argument("--null-tsv", nargs="+", required=True,
                    help="one or more per-replicate per-candidate novelty TSVs")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    observed = pd.read_csv(a.observed_tsv, sep="\t")
    null_by_id: Dict[str, list] = {}
    for path in a.null_tsv:
        rep = pd.read_csv(path, sep="\t")
        for cid, nov in zip(rep["candidate_id"], rep["novelty"]):
            null_by_id.setdefault(cid, []).append(float(nov))

    rows = []
    for cid, nov in zip(observed["candidate_id"], observed["novelty"]):
        nulls = null_by_id.get(cid, [])
        diag = novelty_collapse(float(nov), nulls) if nulls else {
            "null_mean": float("nan"), "collapsed": False}
        rows.append({
            "id": cid,
            "novelty": float(nov),
            "null_mean": diag["null_mean"],
            "p_value": empirical_pvalue(float(nov), nulls) if nulls else float("nan"),
            "collapsed": diag["collapsed"],
        })

    out_df = pd.DataFrame(rows, columns=["id", "novelty", "null_mean", "p_value", "collapsed"])
    out_df.to_csv(a.out, sep="\t", index=False)
    print(f"[embedding_nulls] wrote {len(out_df)} candidates -> {a.out}")


if __name__ == "__main__":
    main()
