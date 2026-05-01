#!/usr/bin/env python3
"""DEPRECATED: convergent_evolution.py has been quarantined.

Bead -0ku: the previous implementation had three fatal flaws:

1. Ancestral-state reconstruction was a parsimony-style postorder majority
   vote across direct children, with NO branch lengths and NO substitution
   model. The pipeline already runs FastML / IQ-TREE --ancestral in step 05;
   this duplicate, inferior implementation should never have been called.

2. When no focal_species was provided, the script split FASTA leaves into
   first-half vs. second-half by alphabetical order. Convergent-evolution
   claims require *phylogenetically independent* lineages, not arbitrary
   text-order splits. This actively manufactured false convergence signal.

3. Root tree.name often defaults to empty string (falsy in Python), which
   bypassed the derived-state filter, flagging any convergent AA as derived.

The correct rebuild should:
- consume pre-computed FastML / IQ-TREE --ancestral sequences (do not
  re-implement ASR);
- require explicit focal-clade specification, OR derive monophyletic clades
  from tree topology;
- properly handle root naming and derived-state checks.

Original implementation kept under scripts/legacy/ for reference.
"""
import sys

print(
    "ERROR: convergent_evolution.py is quarantined (bead -0ku).\n"
    "       The previous implementation produced fabricated convergence\n"
    "       signal (parsimony ASR without branch lengths; arbitrary\n"
    "       half/half clade split). Rebuild required.\n"
    "\n"
    "       Original implementation at scripts/legacy/convergent_evolution.py.\n"
    "       Track the rebuild via beads issue berghia-chemogpcrs-0ku.\n",
    file=sys.stderr,
)
sys.exit(2)
