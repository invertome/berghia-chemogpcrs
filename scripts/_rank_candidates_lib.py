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
import os
import re
import sys
from typing import Iterable, Mapping, Optional, Sequence

import numpy as np


# -----------------------------------------------------------------------------
# `lse_depth` -> `lse_divergence` rename (2026-07-20) — persisted-schema compat
# -----------------------------------------------------------------------------
#
# The axis was named for nesting depth but measures CUMULATIVE BRANCH LENGTH
# (divergence = rate x time). Its companion `lse_nesting_depth` measures the
# topological quantity. The old name therefore denoted the wrong thing and was
# renamed; `lse_nesting_depth` is deliberately NOT touched.
#
# The rename crosses a persisted schema (written CSV columns) and a public
# configuration surface (env vars), so the compatibility policy is explicit:
#
#   READ  legacy names are still accepted, via the two resolvers below, and
#         every acceptance prints a one-time DEPRECATION line naming the file
#         or variable. A legacy read is therefore never SILENT -- which is the
#         whole hazard the rename would otherwise introduce, since a ranked CSV
#         written before the rename would otherwise present no `lse_divergence`
#         column at all and the axis would quietly vanish from the ranking.
#   WRITE only the canonical names are ever emitted. Nothing writes `lse_depth_*`
#         again, so the legacy surface is strictly shrinking.
#
# Stale files are NOT rejected: the legacy column holds exactly the same
# measurement under a wrong label, so refusing it would discard good data and
# break the signal-independence audit's read of the PRIOR run's ranked CSV
# (07_candidate_ranking.sh) on the first post-rename run. Absence of BOTH names
# remains "unavailable" -- never zero.
LEGACY_COLUMN_ALIASES = {
    "lse_divergence_score": "lse_depth_score",
    "lse_divergence_score_norm": "lse_depth_score_norm",
    "has_lse_divergence_data": "has_lse_depth_data",
}

LEGACY_ENV_ALIASES = {
    "LSE_DIVERGENCE_WEIGHT": "LSE_DEPTH_WEIGHT",
    "LSE_DIVERGENCE_PERCENTILE": "LSE_DEPTH_PERCENTILE",
}

_LEGACY_WARNED: set = set()


def _warn_once(key, message) -> None:
    """Emit ``message`` to stderr the first time ``key`` is seen this process."""
    if key in _LEGACY_WARNED:
        return
    _LEGACY_WARNED.add(key)
    print(message, file=sys.stderr)


def apply_legacy_column_aliases(df, source="ranked CSV"):
    """Rename any pre-rename `lse_depth_*` columns in ``df`` to their canonical
    `lse_divergence_*` names, announcing each one.

    Call this ONCE immediately after reading a ranked CSV, rather than
    threading a fallback through every per-signal lookup: a single
    normalization point cannot miss a consumer, and it makes the canonical
    name the only thing any downstream code has to know about.

    A canonical column already present WINS -- the legacy column is then left
    untouched rather than overwriting it, so a hand-merged file carrying both
    cannot silently substitute the stale copy.

    Returns ``df`` (renamed in place via ``DataFrame.rename``; the caller must
    use the return value).
    """
    if df is None or not hasattr(df, "columns"):
        return df
    present = set(df.columns)
    mapping = {
        legacy: canonical
        for canonical, legacy in LEGACY_COLUMN_ALIASES.items()
        if legacy in present and canonical not in present
    }
    if not mapping:
        return df
    for legacy, canonical in sorted(mapping.items()):
        _warn_once(
            ("column", source, legacy),
            f"DEPRECATION: {source} carries the pre-rename column '{legacy}'; "
            f"reading it as '{canonical}'. The lse_depth axis was renamed to "
            f"lse_divergence (it measures cumulative branch length, not nesting "
            f"depth). Re-run stage 07 to emit the canonical schema.",
        )
    return df.rename(columns=mapping)


def getenv_renamed(name, default, *, cast=float):
    """Read env var ``name``, falling back to its pre-rename alias.

    Precedence is canonical-first: if both are set the canonical name wins and
    the legacy one is reported as ignored, so a half-migrated environment can
    never silently apply the value the operator did NOT intend. A legacy-only
    read is honoured and announced. Neither set -> ``default``.
    """
    legacy_name = LEGACY_ENV_ALIASES.get(name)
    canonical = os.getenv(name)
    legacy = os.getenv(legacy_name) if legacy_name else None

    if canonical is not None:
        if legacy is not None and legacy != canonical:
            _warn_once(
                ("env-both", name),
                f"DEPRECATION: both {name}={canonical} and the pre-rename "
                f"{legacy_name}={legacy} are set; using {name}. Unset "
                f"{legacy_name}.",
            )
        return cast(canonical)
    if legacy is not None:
        _warn_once(
            ("env", legacy_name),
            f"DEPRECATION: {legacy_name} is the pre-rename name of {name} and "
            f"is still honoured, but will not be read forever. Rename it to "
            f"{name}.",
        )
        return cast(legacy)
    return cast(default)


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
# Branch-support parsing (bead -i3w9)
# -----------------------------------------------------------------------------
#
# rank_candidates.py loads trees with ``Tree(path, format=1)``. ete3's format=1
# parser puts the Newick internal LABEL into ``node.name`` and leaves
# ``node.support`` at the constructor default of 1.0 -- it does NOT parse the
# label into a number. Every consumer used to read ``node.support`` and compare
# it to BOOTSTRAP_THRESHOLD (70), so on the repo's own class-A tree all 437
# internal nodes reported 1.0 and none passed. The real support is in the label.
#
# Label forms measured in this repository on 2026-07-20 (not assumed):
#
#   producer                                        example label   meaning
#   ----------------------------------------------  --------------  ---------------
#   IQ-TREE `-B N -alrt N` (production, 04 line 63)  '0.777/97'     SH-aLRT/UFBoot
#   IQ-TREE `-alrt N -abayes -B N`                   '88.3/.99/95'  UFBoot LAST
#   IQ-TREE .contree                                 '100'          UFBoot percent
#   FastTree                                         '0.768'        proportion
#   OrthoFinder Gene_Trees                           '1', '0.999'   proportion
#   OrthoFinder Resolved_Gene_Trees                  'n5'           NODE NAME
#   unlabelled OG tree (e.g. OG0000339.treefile)     ''             nothing
#
# Which component is the support? For composite labels, the LAST one: that is
# UFBoot, the value BOOTSTRAP_THRESHOLD=70 conventionally refers to, and
# IQ-TREE always writes it as an integer percent. SH-aLRT is deliberately NOT
# used -- it is emitted on a 0-1 scale by some invocations and 0-100 by others
# (this repo's own treefile is 0-1: 0.229, 0.777, 1.000), so thresholding it at
# 70 would be scale-dependent. UFBoot has no such ambiguity.

_SUPPORT_SEPARATOR = "/"


def parse_support_label(label: Optional[str]):
    """Parse one Newick internal-node label into a raw branch-support number.

    Args:
        label: the internal node label (ete3 format=1 puts this in ``node.name``).

    Returns:
        ``None`` if the label carries no support value at all -- an empty label,
        or a label that is not numeric (OrthoFinder's ``n5`` node names). The
        caller must treat this as "not measured", never as a default value.

        ``(value, True)`` when the value is known to be on the 0-100 percent
        scale (composite IQ-TREE labels, whose trailing UFBoot field always is).

        ``(value, False)`` for a bare number, whose scale cannot be decided from
        one label -- a bare ``1`` is 100% in an OrthoFinder Gene_Tree (siblings
        are 0.999) but 1% in an IQ-TREE .contree. Resolve those together with
        :func:`resolve_support_scale`, which looks at the whole tree.
    """
    if label is None:
        return None
    text = str(label).strip()
    if not text:
        return None

    if _SUPPORT_SEPARATOR in text:
        parts = [p.strip() for p in text.split(_SUPPORT_SEPARATOR)]
        if any(not p for p in parts):
            return None                      # malformed, e.g. '0.9/' or '//'
        try:
            values = [float(p) for p in parts]
        except ValueError:
            return None
        return (values[-1], True)            # UFBoot is last, always percent

    try:
        return (float(text), False)
    except ValueError:
        return None                          # 'n5', 'abc', ...


def resolve_support_scale(raw_values: Sequence[float]) -> float:
    """Multiplier converting bare support values to the 0-100 percent scale.

    Decided once per tree rather than per label, because a single bare value is
    genuinely ambiguous. If every bare value in the tree is <= 1.0 the tree is
    proportion-scaled (FastTree, OrthoFinder Gene_Trees) and must be multiplied
    by 100; if any exceeds 1.0 the tree is already percent-scaled (IQ-TREE
    .contree). Without this, a perfectly supported FastTree node (1.000) would
    read as 1% and fail a threshold of 70 -- the same failure mode as the ete3
    default this bead fixes.
    """
    values = [v for v in raw_values if v is not None]
    if not values:
        return 1.0
    return 100.0 if max(values) <= 1.0 else 1.0


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
    dnds_reliability: float = 1.0,
    dnds_axes: Iterable[str] = ("positive", "purifying"),
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

    # dN/dS reliability (bead -8st): scale the positive/purifying axes' weight
    # by dnds_reliability in [0, 1] everywhere — the score AND the completeness
    # denominator (total_weight) — so an under-supported omega estimate cleanly
    # falls out and the candidate is judged on its other axes (clean-removal).
    _dnds = set(dnds_axes)
    _rw = max(0.0, min(1.0, float(dnds_reliability)))

    def _eff_weight(name, w):
        w = max(float(w), 0.0)
        return w * _rw if name in _dnds else w

    total_weight = float(sum(_eff_weight(k, w) for k, w in weights.items()))
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
        ew = _eff_weight(k, w)
        avail_weight += ew
        weighted_sum += ew * float(v)

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
