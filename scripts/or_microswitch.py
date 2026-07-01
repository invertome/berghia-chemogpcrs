#!/usr/bin/env python3
# or_microswitch.py
# Purpose: cheap OR-microswitch structural flag for the ML/PLM chemoreceptor
#   ranking reranker (Task 6, docs/plans/2026-07-01-ml-plm-chemoreceptor-
#   ranking.md), from Ballesteros-Weinstein (BW) numbered residues.
"""Score the "altered TM6 toggle" chemoreceptor structural fingerprint.

Ballesteros-Weinstein (BW) position 6.48 is the class-A GPCR "toggle switch"
residue, canonically Trp (W), anchoring the CWxP motif that spans TM6
positions 6.47-6.50 (Cys6.47-Trp6.48-x6.49-Pro6.50). A documented structural
signature in a subset of olfactory/chemosensory-receptor-like GPCRs (hence
"OR-microswitch") is a SUBSTITUTED toggle residue -- Tyr (Y) or Phe (F) at
6.48 -- while the flanking Pro6.50 that anchors the motif stays conserved.
This module scores that altered-motif signature from a per-candidate,
already-BW-numbered residue map.

The upstream step that assigns Ballesteros-Weinstein numbers to an
AlphaFold model's residues (structural alignment to a numbered class-A GPCR
template) is OUT OF SCOPE here -- it is a Unity/upstream step. This module is
a pure, cheap classifier over an already-produced ``{bw_position: residue}``
map for one candidate; it never touches a structure file directly.
"""
from __future__ import annotations

from typing import Dict, Optional

# BW position of the canonical class-A "toggle switch" residue (Trp6.48 in
# most class-A GPCRs), and the invariant downstream proline that anchors the
# CWxP motif (Cys6.47-Trp6.48-x6.49-Pro6.50). Both must be present in the
# input map, or there isn't enough BW-numbered structure to call the flag
# either way.
TOGGLE_POSITION = "6.48"
MOTIF_ANCHOR_POSITION = "6.50"
ALTERED_TOGGLE_RESIDUES = {"Y", "F"}
CONSERVED_ANCHOR_RESIDUE = "P"


def or_microswitch_flag(bw_residues: Optional[Dict[str, str]]) -> Optional[int]:
    """1/0/None flag for the altered-TM6-toggle chemoreceptor microswitch.

    Args:
        bw_residues: a per-candidate ``{bw_position: residue}`` map, e.g.
            ``{"6.48": "Y", "6.50": "P", ...}``.

    Returns:
        1 -- position 6.48 is Tyr or Phe (the altered toggle) AND 6.50 is
            the conserved Pro anchoring the CWxP motif: the documented
            "olfactory-GPCR" structural fingerprint.
        0 -- both required positions are present but the pattern above is
            not met (e.g. the canonical Trp6.48, or an altered toggle
            without the conserved anchor).
        None -- `bw_residues` is falsy, or either required BW position
            (6.48 and/or 6.50) is absent -- insufficient BW-numbered
            residues to call the flag either way.
    """
    if not bw_residues:
        return None
    if TOGGLE_POSITION not in bw_residues or MOTIF_ANCHOR_POSITION not in bw_residues:
        return None
    toggle = str(bw_residues[TOGGLE_POSITION]).strip().upper()
    anchor = str(bw_residues[MOTIF_ANCHOR_POSITION]).strip().upper()
    altered = toggle in ALTERED_TOGGLE_RESIDUES and anchor == CONSERVED_ANCHOR_RESIDUE
    return 1 if altered else 0
