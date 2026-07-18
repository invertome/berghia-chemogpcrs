#!/usr/bin/env python3
"""Model-role split producer (epic v4bs, task v4bs.4 / A3).

Embedding distance-to-reference-centroid conflates PHYLOGENY with FUNCTION.
Different PLMs specialize: ESM-2-style embeddings track phylogenetic distance,
while text/function-aligned models (ProTrek, ProteinCLIP) track functional
novelty. This module turns that specialization into a signal: it converts each
model's per-candidate novelty score into a within-model novelty percentile, then
emits the GAP between a "function" model and a "phylo" model. A positive gap
marks a candidate that is functionally novel BEYOND its phylogenetic expectation
(a functional shift without much sequence divergence — genuinely interesting); a
negative gap marks a phylo-distant but functionally typical distant homolog.

A dormant descriptive producer (like the cw3.6 embedding channel): pure,
GPU-free functions plus a thin CLI. It is NOT wired into any production ranking
voter set here — the controller decides that later. The math operates on plain
``{candidate_id: score}`` mappings; pandas is used only inside ``main`` (lazy
import) so the pure functions stay importable without pandas/torch.
"""
from __future__ import annotations

from typing import Dict, List, Mapping, Optional, Sequence

import numpy as np
from scipy.stats import rankdata


def novelty_percentiles(scores: Mapping[str, float]) -> Dict[str, float]:
    """cid -> novelty percentile in (0, 1], HIGH = more novel.

    Higher input score = more novel, so scores are ranked ASCENDING (rank 1 = the
    least-novel, lowest score) and divided by ``n``. The maximum-score candidate
    therefore maps to ``n / n = 1.0`` and every value lands in ``(0, 1]``. Ties
    receive the average of their ranks (``scipy.rankdata(method="average")``), so
    equal scores always get identical, deterministic percentiles.
    """
    ids = list(scores)
    r = rankdata([scores[i] for i in ids], method="average")  # 1 = least novel
    n = len(ids)
    return {i: float(r[k] / n) for k, i in enumerate(ids)}


def role_gap(
    phylo_novelty: Mapping[str, float],
    function_novelty: Mapping[str, float],
) -> Dict[str, float]:
    """cid -> ``function_novelty_percentile - phylo_novelty_percentile``.

    Each model's novelty is converted to a percentile WITHIN that model's own
    candidate set (:func:`novelty_percentiles`, rank-based so it is scale-free),
    then the gap is emitted only for candidates present in BOTH maps. The sign
    is the function-vs-phylogeny disentanglement signal:

      * ``> 0`` — novel in function beyond phylogenetic expectation (a functional
        shift without much sequence divergence: genuinely interesting);
      * ``~ 0`` — the two roles agree;
      * ``< 0`` — phylo-distant but functionally typical (just a distant homolog).
    """
    phylo_pct = novelty_percentiles(phylo_novelty)
    function_pct = novelty_percentiles(function_novelty)
    common = set(phylo_pct) & set(function_pct)
    return {cid: function_pct[cid] - phylo_pct[cid] for cid in common}


# --- CLI glue --------------------------------------------------------------
# One row per intersection candidate, sorted by emb_role_gap descending.
ROLE_SPLIT_TSV_COLUMNS: List[str] = [
    "id",
    "emb_novelty_phylo",
    "emb_novelty_function",
    "emb_role_gap",
]


def main(argv: Optional[Sequence[str]] = None) -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--phylo-tsv", required=True,
        help="phylo-model novelty TSV (columns candidate_id, novelty — the "
             "schema of results/ranking/diagnostics/novelty_<tag>_PROD.tsv)",
    )
    parser.add_argument(
        "--function-tsv", required=True,
        help="function-model novelty TSV (same candidate_id, novelty schema)",
    )
    parser.add_argument("--out", required=True, help="output role-split TSV")
    args = parser.parse_args(argv)

    # Lazy IO import so the pure functions above stay pandas-free for unit tests.
    import pandas as pd

    def _load(path: str) -> Dict[str, float]:
        df = pd.read_csv(path, sep="\t")
        return dict(zip(df["candidate_id"], df["novelty"].astype(float)))

    phylo = _load(args.phylo_tsv)
    function = _load(args.function_tsv)
    phylo_pct = novelty_percentiles(phylo)
    function_pct = novelty_percentiles(function)
    gap = role_gap(phylo, function)

    rows = [
        {
            "id": cid,
            "emb_novelty_phylo": phylo_pct[cid],
            "emb_novelty_function": function_pct[cid],
            "emb_role_gap": gap[cid],
        }
        for cid in gap
    ]
    out_df = pd.DataFrame(rows, columns=ROLE_SPLIT_TSV_COLUMNS)
    out_df = out_df.sort_values("emb_role_gap", ascending=False, kind="mergesort")
    out_df.to_csv(args.out, sep="\t", index=False)
    print(f"[model_role_split] wrote {len(out_df)} candidates "
          f"(phylo∩function) -> {args.out}")


if __name__ == "__main__":
    main()
