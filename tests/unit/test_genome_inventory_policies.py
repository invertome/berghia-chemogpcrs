"""Tests for scripts/genome_inventory_policies.py — clade-policy config loader."""
from __future__ import annotations

from pathlib import Path

import genome_inventory_policies as gip
from build_species_tree_phase1d_extension_inventory import POLICIES


def test_loads_fields(tmp_path: Path) -> None:
    cfg = tmp_path / "p.tsv"
    cfg.write_text(
        "clade_name\tpolicy_class\tmin_assembly_level\tmax_count\n"
        "Heterobranchia\theterobranchia\t\t\n"
        "Annelida\toutgroup_annelida\tChromosome\t2\n"
    )
    pols = gip.load_clade_policies(cfg)
    assert len(pols) == 2
    assert pols[0].clade_name == "Heterobranchia"
    assert pols[0].min_assembly_level == ""
    assert pols[0].max_count is None
    assert pols[0].require_annotation is False
    assert pols[1].max_count == 2
    assert pols[1].min_assembly_level == "Chromosome"


def test_config_matches_legacy_hardcoded_policies() -> None:
    """The shipped clade_policies.tsv must reproduce the legacy POLICIES list."""
    loaded = gip.load_clade_policies(gip.DEFAULT_POLICY_CONFIG)
    assert loaded == list(POLICIES)
