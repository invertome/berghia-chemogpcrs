"""Pure helpers for Aplysia retrospective validation (bead -bdu).

Used by validate_against_aplysia.sh / .py to compute recall@N of known
chemoreceptor genes from a labeled CSV against a pipeline-produced
ranked CSV.
"""
from __future__ import annotations

from typing import Iterable, Optional, Sequence

import pandas as pd


def recall_at_n(ranked_ids: Sequence[str],
                truth_ids: Iterable[str],
                n: int) -> float:
    """Fraction of truth IDs found in the top-N of ranked_ids.

    truth_ids is treated as a set; duplicates are ignored.
    """
    truth = set(truth_ids)
    if not truth:
        return 0.0
    n = min(max(int(n), 0), len(ranked_ids))
    top_set = set(ranked_ids[:n])
    return len(top_set & truth) / len(truth)


def precision_at_n(ranked_ids: Sequence[str],
                   truth_ids: Iterable[str],
                   n: int) -> float:
    """Fraction of top-N candidates that are in the truth set."""
    truth = set(truth_ids)
    n = min(max(int(n), 0), len(ranked_ids))
    if n == 0:
        return 0.0
    top_set = set(ranked_ids[:n])
    return len(top_set & truth) / n


def average_rank(ranked_ids: Sequence[str],
                 truth_ids: Iterable[str]) -> Optional[float]:
    """Average 1-based rank of truth IDs in ranked_ids; None if none present."""
    rank_of: dict[str, int] = {gid: i + 1 for i, gid in enumerate(ranked_ids)}
    found = [rank_of[t] for t in set(truth_ids) if t in rank_of]
    if not found:
        return None
    return sum(found) / len(found)


def recall_curve(ranked_ids: Sequence[str],
                 truth_ids: Iterable[str],
                 ns: Iterable[int] = (10, 25, 50, 100, 250, 500)) -> pd.DataFrame:
    """Return recall + precision at each N as a DataFrame."""
    truth = set(truth_ids)
    n_total = len(truth)
    rows = []
    for n in ns:
        rows.append({
            "N": int(n),
            "recall": recall_at_n(ranked_ids, truth, n),
            "precision": precision_at_n(ranked_ids, truth, n),
            "n_truth_total": n_total,
            "n_truth_found_in_topN": int(round(recall_at_n(ranked_ids, truth, n) * n_total)),
        })
    return pd.DataFrame(rows)
