"""Load clade-selection policies from references/species_tree/clade_policies.tsv.

Replaces the hardcoded POLICIES list in the phase1d extension builder. Reuses
the ClaadePolicy dataclass so existing selection logic is unchanged.
"""
from __future__ import annotations

import csv
from pathlib import Path

from build_species_tree_phase1d_extension_inventory import ClaadePolicy

# Anchor to the repo (scripts/ → repo root) so the default resolves regardless
# of the caller's working directory (this constant is imported by other modules).
DEFAULT_POLICY_CONFIG = (
    Path(__file__).resolve().parent.parent / "references/species_tree/clade_policies.tsv"
)


def load_clade_policies(path: str | Path = DEFAULT_POLICY_CONFIG) -> list[ClaadePolicy]:
    """Parse the clade-policy TSV into ClaadePolicy objects, preserving order.

    Columns: clade_name, policy_class, min_assembly_level, max_count.
    Empty max_count cell → None (no cap). require_annotation is always False
    (Phase 1f BRAKER4-annotates, so unannotated GenBank assemblies are kept).
    """
    policies: list[ClaadePolicy] = []
    with Path(path).open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            cap = (row.get("max_count") or "").strip()
            policies.append(ClaadePolicy(
                clade_name=row["clade_name"].strip(),
                policy_class=row["policy_class"].strip(),
                min_assembly_level=(row.get("min_assembly_level") or "").strip(),
                require_annotation=False,
                max_count=int(cap) if cap else None,
            ))
    return policies
