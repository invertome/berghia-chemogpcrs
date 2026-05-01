#!/usr/bin/env python3
"""DEPRECATED: binding_site_prediction.py has been quarantined.

Bead -hg1: the previous implementation had two fatal flaws:

1. The Ballesteros-Weinstein numbering dictionary at module top was DEAD CODE
   — the script never actually aligned candidate sequences to the canonical
   Class A topology, so positions like 3.32 were never mapped to any residue.

2. ``is_pocket_facing = dist_to_center < median(distances)`` made
   "pocket-facing" mathematically equivalent to "the closer half of TM
   residues" by definition. Cavity volume estimates were thus inflated.

The correct rebuild requires either:
- GPCRdb structure-based BW alignment (e.g. via gpcrdb-data-cli), OR
- fpocket cavity detection on AlphaFold structures, OR
- TM-aligner (Hildebrand 2009) + reference Class A inactive-state structure

Until then, downstream code should treat any binding_site_* columns in CSVs
as missing data (not zero), so the composite-score evidence-completeness
multiplier handles them correctly.

Original implementation kept under scripts/legacy/ for reference.
"""
import sys

print(
    "ERROR: binding_site_prediction.py is quarantined (bead -hg1).\n"
    "       The previous implementation produced misleading binding-pocket\n"
    "       calls (Ballesteros-Weinstein numbering was dead code; pocket-\n"
    "       facing was a tautology). Rebuild required.\n"
    "\n"
    "       Original implementation at scripts/legacy/binding_site_prediction.py.\n"
    "       Track the rebuild via beads issue berghia-chemogpcrs-hg1.\n",
    file=sys.stderr,
)
sys.exit(2)
