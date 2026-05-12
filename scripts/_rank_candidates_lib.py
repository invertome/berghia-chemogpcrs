"""Pure, importable scoring/statistics functions used by rank_candidates.py.

This module exists separately from rank_candidates.py so that unit tests can
exercise the functions without triggering the heavy module-level data loading
and main-loop execution that rank_candidates.py performs at import time.

All functions here MUST be free of module-level side effects: no file I/O,
no env var reads at import time, no argparse, no global mutable state.

Authored as part of the 2026-05-01 review-fix sweep (see
docs/plans/2026-05-01-pipeline-review-fixes.md). Bead refs in docstrings.
"""
from __future__ import annotations

import math
import re
from typing import Iterable, Mapping, Optional, Sequence

import numpy as np


# -----------------------------------------------------------------------------
# Reference categorization (bead -wux part 2)
# -----------------------------------------------------------------------------

# Word-boundary regex avoids the 'or' substring bug that previously matched
# Orexin, Orphanin, Origin, Adrenergic, etc. Keeps olfactory-receptor naming
# conventions: OR1A1, OLFR123, V1R, V2R, T1R, T2R, TAAR, etc.
# Custom word-boundary that treats underscore as a separator (Python's \b
# treats _ as a word char, which would prevent matching e.g. "OR1A1_HUMAN").
_CHEMORECEPTOR_PATTERNS = re.compile(
    r"(?<![A-Za-z0-9])("
    r"olfr\d*"                           # mammalian olfactory receptor prefix (Olfr, OLFR1234)
    r"|or\d+[a-z]?\d*"                   # OR1A1, Or22a, OR7D4, etc.
    r"|olfactory\s*receptor"             # "Olfactory receptor 7A1"
    r"|odorant\s*receptor"
    r"|taste\s*receptor|tasr\d*|t[12]r\d*"
    r"|vomero(?:nasal)?(?:\s*receptor)?|v[12]r\d*"
    r"|formyl[\s_-]?peptide(?:\s*receptor)?|fpr\d*"
    r"|trace[\s_-]?amine(?:[\s_-]?associated)?(?:\s*receptor)?|taar\d*"
    r"|gustatory(?:\s*receptor)?"
    r"|chemosensory(?:\s*receptor)?"
    r"|chemoreceptor"
    r")(?![A-Za-z0-9])",
    re.IGNORECASE,
)


def categorize_reference(
    ref_name: str,
    *,
    chemoreceptor_weight: float = 2.0,
    other_gpcr_weight: float = 1.0,
    explicit_category_map: Optional[Mapping[str, float]] = None,
    explicit_chemoreceptor_set: Optional[Iterable[str]] = None,
) -> float:
    """Return reference weight based on category.

    Priority:
      1. Explicit category map (e.g. from ref_categories_final.csv)
      2. Explicit chemoreceptor set (e.g. legacy JSON-loaded list)
      3. Word-boundary keyword regex (no 'or' substring false-positives)

    Args:
        ref_name: header/name of the reference sequence
        chemoreceptor_weight: weight to return for chemoreceptor refs
        other_gpcr_weight: weight to return for non-chemoreceptor GPCR refs
        explicit_category_map: optional name -> weight mapping (highest priority)
        explicit_chemoreceptor_set: optional collection of names known to be
            chemoreceptors (second priority)

    Returns:
        Weight (float) to apply when this reference is used in scoring.
    """
    if not ref_name:
        return other_gpcr_weight
    if explicit_category_map and ref_name in explicit_category_map:
        return float(explicit_category_map[ref_name])
    if explicit_chemoreceptor_set and ref_name in explicit_chemoreceptor_set:
        return float(chemoreceptor_weight)
    if _CHEMORECEPTOR_PATTERNS.search(ref_name):
        return float(chemoreceptor_weight)
    return float(other_gpcr_weight)


# -----------------------------------------------------------------------------
# Benjamini-Hochberg FDR (bead -wux part 1)
# -----------------------------------------------------------------------------

def benjamini_hochberg(pvalues: Sequence[float]) -> list[float]:
    """Benjamini-Hochberg FDR-corrected q-values, preserving input order.

    Uses statsmodels.stats.multitest.multipletests as the reference
    implementation; replaces a previously hand-rolled version that had a
    rank-indexing bug producing non-monotonic q-values.

    NaNs in input are propagated to NaN q-values (excluded from correction).

    Args:
        pvalues: list/array of p-values in [0, 1] (NaN allowed).

    Returns:
        List of q-values aligned to the input order.
    """
    n = len(pvalues)
    if n == 0:
        return []
    from statsmodels.stats.multitest import multipletests  # local import: only when used

    arr = np.asarray(pvalues, dtype=float)
    nan_mask = np.isnan(arr)
    if nan_mask.all():
        return [float("nan")] * n
    out = np.full(n, np.nan, dtype=float)
    ok = ~nan_mask
    _, q_ok, _, _ = multipletests(arr[ok], method="fdr_bh")
    out[ok] = q_ok
    return out.tolist()


# -----------------------------------------------------------------------------
# aBSREL omega rate-class extraction (bead -ea9 part 2)
# -----------------------------------------------------------------------------

def extract_branch_omega(branch_data: dict) -> dict:
    """Extract omega rate-class statistics from a single aBSREL branch entry.

    aBSREL fits a mixture of rate classes per branch; reporting the weighted
    mean ω destroys the episodic-positive-selection signal aBSREL is designed
    to detect. This function reports both omega_max (max across rate classes
    -- the chemoreceptor-relevant signal) and omega_mean (legacy compatibility).

    Args:
        branch_data: a single branch's dict from aBSREL JSON output, expected to
            contain 'Rate Distributions' as a list of [omega, weight] pairs.

    Returns:
        Dict with keys: omega_max, omega_mean, weight_at_max, n_rate_classes.
        Returns NaNs and zero rate-class count if input is missing/empty.
    """
    rate_classes = branch_data.get("Rate Distributions", []) or []
    if not rate_classes:
        return {
            "omega_max": float("nan"),
            "omega_mean": float("nan"),
            "weight_at_max": float("nan"),
            "n_rate_classes": 0,
        }
    pairs = [(float(rc[0]), float(rc[1])) for rc in rate_classes]
    omegas = [p[0] for p in pairs]
    weights = [p[1] for p in pairs]
    omega_max = max(omegas)
    weight_at_max = weights[omegas.index(omega_max)]
    omega_mean = sum(o * w for o, w in pairs)
    return {
        "omega_max": omega_max,
        "omega_mean": omega_mean,
        "weight_at_max": weight_at_max,
        "n_rate_classes": len(pairs),
    }


# -----------------------------------------------------------------------------
# Selection scoring (bead -ea9 part 1)
# -----------------------------------------------------------------------------

def get_selection_scores(
    *,
    omega: float,
    p_corrected: float,
    purifying_weight: float = 0.0,
    positive_weight: float = 1.0,
    significance_threshold: float = 0.05,
    significance_boost: float = 1.5,
) -> dict:
    """Compute purifying and positive selection scores from a branch's omega.

    Critical fix (bead -ea9): default ``purifying_weight`` is 0 because
    chemoreceptor identification rewards diversifying selection on
    extracellular loops, not whole-gene purifying selection. Set
    ``purifying_weight`` > 0 only when looking for conserved-function GPCRs.

    Args:
        omega: the dN/dS value to score (use omega_max from
            ``extract_branch_omega`` for the chemoreceptor-relevant signal).
        p_corrected: BH-FDR-corrected p-value from aBSREL/BUSTED.
        purifying_weight: scaling for the purifying-selection axis (default 0).
        positive_weight: scaling for the positive-selection axis (default 1).
        significance_threshold: q-value threshold for the significance boost.
        significance_boost: multiplicative boost when significant.

    Returns:
        Dict with keys: purifying_score, positive_score, is_significant.
        NaN omega -> zero scores.
    """
    if omega is None or (isinstance(omega, float) and math.isnan(omega)):
        return {"purifying_score": 0.0, "positive_score": 0.0, "is_significant": False}

    is_significant = (
        p_corrected is not None
        and not (isinstance(p_corrected, float) and math.isnan(p_corrected))
        and float(p_corrected) < significance_threshold
    )

    # Use log10 magnitude so "ω = 0.05" and "ω = 20" both produce log10 magnitude ~1.3
    if omega < 1.0:
        purifying_raw = abs(math.log10(max(omega, 1e-9)))
        positive_raw = 0.0
    elif omega > 1.0:
        purifying_raw = 0.0
        positive_raw = math.log10(omega)
    else:
        purifying_raw = 0.0
        positive_raw = 0.0

    purifying_score = purifying_weight * purifying_raw
    positive_score = positive_weight * positive_raw
    if is_significant:
        purifying_score *= significance_boost
        positive_score *= significance_boost
    return {
        "purifying_score": float(purifying_score),
        "positive_score": float(positive_score),
        "is_significant": bool(is_significant),
    }


# -----------------------------------------------------------------------------
# Composite score with evidence-completeness multiplier (bead -ce4)
# -----------------------------------------------------------------------------

def calculate_fair_rank_score(
    scores: Mapping[str, Optional[float]],
    weights: Mapping[str, float],
    *,
    completeness_floor: float = 0.4,
    return_diagnostics: bool = False,
):
    """Combine per-axis scores into a single composite, penalizing missing data.

    The previous implementation normalized by *available* weight only, so a
    candidate with one strong signal could rank equal to a candidate with all
    signals at the same per-axis level. This version multiplies the
    available-weight-normalized score by ``evidence_completeness`` (the
    fraction of weight that was actually contributed), floored at
    ``completeness_floor`` so a single very strong signal is not zeroed.

    Args:
        scores: name -> score in [0, 1] or None for missing.
        weights: name -> weight (>= 0).
        completeness_floor: minimum multiplier for evidence_completeness.
        return_diagnostics: if True, return dict with score + diagnostics.

    Returns:
        Float in [0, 1] (or dict if return_diagnostics).
    """
    if not weights:
        out = {"score": 0.0, "evidence_completeness": 0.0, "available_weight": 0.0}
        return out if return_diagnostics else 0.0

    total_weight = float(sum(max(w, 0.0) for w in weights.values()))
    if total_weight <= 0.0:
        out = {"score": 0.0, "evidence_completeness": 0.0, "available_weight": 0.0}
        return out if return_diagnostics else 0.0

    avail_weight = 0.0
    weighted_sum = 0.0
    for k, w in weights.items():
        v = scores.get(k)
        if v is None:
            continue
        if isinstance(v, float) and math.isnan(v):
            continue
        avail_weight += float(w)
        weighted_sum += float(w) * float(v)

    if avail_weight <= 0.0:
        out = {"score": 0.0, "evidence_completeness": 0.0, "available_weight": 0.0}
        return out if return_diagnostics else 0.0

    completeness_raw = avail_weight / total_weight
    completeness = max(completeness_floor, completeness_raw)
    score = (weighted_sum / avail_weight) * completeness
    out = {
        "score": float(score),
        "evidence_completeness": float(completeness),
        "evidence_completeness_raw": float(completeness_raw),
        "available_weight": float(avail_weight),
    }
    return out if return_diagnostics else float(score)


# -----------------------------------------------------------------------------
# Synteny score normalization (bead -ce4 part 2 / -mqt)
# -----------------------------------------------------------------------------

def normalize_synteny_counts(
    counts: Mapping[str, int],
    *,
    min_max_anchors: int = 5,
) -> dict[str, Optional[float]]:
    """Normalize per-candidate synteny anchor counts to a [0, 1] score.

    The previous implementation divided by ``max(counts)``, which produced
    a score of 1.0 for any candidate with one anchor when the dataset-wide
    max was also one — degenerate. This version uses log1p scaling and
    returns ``None`` (missing) for every candidate when the dataset-wide
    max is below ``min_max_anchors``, so the composite-score function's
    evidence-completeness multiplier handles it correctly.

    Args:
        counts: name -> integer anchor count.
        min_max_anchors: dataset-wide max below which all scores are None.

    Returns:
        Dict mapping the same names to ``Optional[float]`` scores in (0, 1].
    """
    if not counts:
        return {}
    max_count = max(counts.values()) if counts else 0
    if max_count < min_max_anchors:
        return {k: None for k in counts}
    denom = math.log1p(float(max_count))
    if denom <= 0.0:
        return {k: None for k in counts}
    return {k: (math.log1p(float(v)) / denom if v > 0 else 0.0) for k, v in counts.items()}


# -----------------------------------------------------------------------------
# MEME concordance loader (bead -7cy step 2)
# -----------------------------------------------------------------------------

def load_meme_concordance(path: str) -> dict:
    """Read parse_meme.py's dual-mode per-OG concordance CSV.

    The CSV is the output of ``parse_meme.py --lenient-json ...`` (one row
    per OG). Each row carries strict / lenient / high_confidence /
    lenient_only / strict_only counts plus the alignment_robustness_index.

    Returns dict og_name -> {high_confidence_sites_n, lenient_only_sites_n,
    strict_only_sites_n, n_strict_positive_sites, n_lenient_positive_sites,
    alignment_robustness_index}. Missing file returns {}; blank og_name
    rows are skipped.
    """
    import csv as _csv
    import os as _os
    out: dict = {}
    if not _os.path.exists(path):
        return out
    with open(path, newline="") as f:
        reader = _csv.DictReader(f)
        for row in reader:
            og = (row.get("og_name") or "").strip()
            if not og:
                continue
            def _i(k):
                v = row.get(k, "")
                try:
                    return int(v) if v not in ("", None) else 0
                except (TypeError, ValueError):
                    return 0
            def _f(k):
                v = row.get(k, "")
                try:
                    return float(v) if v not in ("", None) else 0.0
                except (TypeError, ValueError):
                    return 0.0
            out[og] = {
                "n_strict_positive_sites": _i("n_strict_positive_sites"),
                "n_lenient_positive_sites": _i("n_lenient_positive_sites"),
                "high_confidence_sites_n": _i("high_confidence_sites_n"),
                "lenient_only_sites_n": _i("lenient_only_sites_n"),
                "strict_only_sites_n": _i("strict_only_sites_n"),
                "alignment_robustness_index": _f("alignment_robustness_index"),
            }
    return out
