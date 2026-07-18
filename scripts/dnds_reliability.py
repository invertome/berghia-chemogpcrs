#!/usr/bin/env python3
"""dnds_reliability.py — reliability-aware shrinkage of a selection signal.

A per-candidate dN/dS positive/purifying-selection estimate is UNDERPOWERED for
reference-poor divergent candidates (aBSREL/BUSTED power scales with divergence
and #informative branches; Anisimova, Bielawski & Yang 2001). "No positive
selection" there is a power artifact, not conservation. This shrinks each
candidate's signal toward the cohort NEUTRAL (median) in proportion to
(1 - reliability), so an underpowered candidate contributes a neutral value
(never penalized) and a fully-powered one is unchanged. Reliability in [0,1] is
the existing per-OG proxy min(1, n_ref_cds / DNDS_RELIABILITY_FULL).
"""
from __future__ import annotations

from statistics import median
from typing import Dict, Mapping


def reliability_shrink(
    values: Mapping[str, float],
    reliability: Mapping[str, float],
    neutral: float | None = None,
) -> Dict[str, float]:
    """Return values shrunk toward ``neutral`` (default: their median) by each
    candidate's reliability. adj = neutral + w * (value - neutral); w defaults
    to 1.0 for any candidate absent from ``reliability``; w is clamped to [0,1]."""
    if not values:
        return {}
    neutral = median(values.values()) if neutral is None else float(neutral)
    out: Dict[str, float] = {}
    for cid, v in values.items():
        w = reliability.get(cid, 1.0)
        try:
            w = float(w)
        except (TypeError, ValueError):
            w = 1.0
        w = min(1.0, max(0.0, w))
        # Lerp form of ``neutral + w*(v-neutral)``; exact at the endpoints
        # (w=1 -> v, w=0 -> neutral) so full/zero reliability round-trips
        # without floating-point drift.
        out[cid] = float(v) * w + neutral * (1.0 - w)
    return out
