"""Pure helpers for HCR probe-design diagnostic columns (bead -xqz).

The wet-lab pipeline validates each ranked candidate one at a time via HCR
(in situ hybridization chain reaction). HCR probes are typically
~30–50 nt, so probe specificity depends on local sequence divergence
between a candidate and its closest paralogs. This module computes
diagnostic columns that help wet-lab staff pick which candidates to
design probes against:

- ``cds_length_bp`` and ``protein_length_aa``
- ``paralog_minimum_identity``: highest pairwise % identity to any other
  GPCR over a 300 aa local window (proxy for HCR probe specificity).
  Higher = harder to design specific probes.
- ``hcr_probe_friendly``: boolean, True iff the candidate is reasonably
  long AND its closest paralog is divergent enough that probe design
  is plausible.

The thresholds below are PROXIES, not validated cutoffs from a
benchmarked HCR design study. Document them as such in the manuscript.

Library design rule: keep this module pure and importable; rank_candidates.py
orchestrates the I/O and calls these functions.
"""
from __future__ import annotations

from typing import Iterable, Mapping, Optional, Sequence


# Probe-friendliness thresholds (PROXIES — not benchmark-validated).
# Documented for transparency: a candidate is "probe-friendly" if its
# closest paralog has <= 80% identity over a sliding 300-aa window AND
# the CDS is long enough to permit several non-overlapping probe pairs.
HCR_DEFAULT_MAX_PARALOG_IDENTITY = 0.80
HCR_DEFAULT_MIN_CDS_LENGTH_BP = 600


def compute_pairwise_identity(seq_a: str, seq_b: str) -> float:
    """Fraction of matching residues over the aligned overlap length.

    Both sequences should already be aligned columnar (same length); this
    function does NOT realign. Gaps on either side are excluded from
    the denominator. Returns 0.0 if there is no overlap.
    """
    if len(seq_a) != len(seq_b):
        raise ValueError("compute_pairwise_identity expects aligned sequences "
                         f"of equal length (got {len(seq_a)} and {len(seq_b)})")
    matches = 0
    counted = 0
    for a, b in zip(seq_a, seq_b):
        if a in ("-", ".") or b in ("-", "."):
            continue
        counted += 1
        if a == b:
            matches += 1
    return matches / counted if counted else 0.0


def windowed_max_identity(seq_a: str, seq_b: str,
                          window_aa: int = 300,
                          step_aa: int = 50) -> float:
    """Maximum pairwise identity over any sliding window of length ``window_aa``.

    Returns 0.0 if both sequences are shorter than the window or contain
    no overlap.
    """
    n = min(len(seq_a), len(seq_b))
    if n == 0:
        return 0.0
    if n < window_aa:
        return compute_pairwise_identity(seq_a[:n], seq_b[:n])
    best = 0.0
    for start in range(0, n - window_aa + 1, max(1, step_aa)):
        end = start + window_aa
        ident = compute_pairwise_identity(seq_a[start:end], seq_b[start:end])
        if ident > best:
            best = ident
    return best


def closest_paralog_identity(target_id: str,
                             aligned_seqs: Mapping[str, str],
                             window_aa: int = 300,
                             step_aa: int = 50) -> tuple[Optional[str], float]:
    """For ``target_id``, find the paralog with the highest windowed identity.

    Returns (closest_paralog_id, max_identity). Returns (None, 0.0) if the
    target is missing or there are no other sequences to compare against.
    """
    if target_id not in aligned_seqs:
        return None, 0.0
    target = aligned_seqs[target_id]
    best_id: Optional[str] = None
    best_ident = 0.0
    for other_id, other in aligned_seqs.items():
        if other_id == target_id:
            continue
        ident = windowed_max_identity(target, other,
                                      window_aa=window_aa, step_aa=step_aa)
        if ident > best_ident:
            best_ident = ident
            best_id = other_id
    return best_id, best_ident


def is_hcr_probe_friendly(*,
                          cds_length_bp: Optional[int],
                          paralog_min_identity: Optional[float],
                          tandem_cluster_size: Optional[int] = None,
                          max_paralog_identity: float = HCR_DEFAULT_MAX_PARALOG_IDENTITY,
                          min_cds_length_bp: int = HCR_DEFAULT_MIN_CDS_LENGTH_BP,
                          tandem_cluster_warning: int = 5) -> bool:
    """Return True iff the candidate is plausibly suitable for HCR probe design.

    Rules:
      * cds_length_bp >= min_cds_length_bp
      * paralog_min_identity <= max_paralog_identity
        (or unknown — treat unknown as friendly so we don't exclude
        candidates that simply have no aligned paralog)
      * If tandem_cluster_size is provided and >= tandem_cluster_warning,
        the candidate is flagged unfriendly because probe-binding regions
        are likely to be shared across cluster members.
    """
    if cds_length_bp is None or cds_length_bp < min_cds_length_bp:
        return False
    if (paralog_min_identity is not None
            and paralog_min_identity > max_paralog_identity):
        return False
    if (tandem_cluster_size is not None
            and tandem_cluster_size >= tandem_cluster_warning):
        return False
    return True
