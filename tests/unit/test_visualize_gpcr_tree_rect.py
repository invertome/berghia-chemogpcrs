"""Tests for the pure (non-plotting) logic in visualize_gpcr_tree_rect.py.

The rendering functions are side-effecting matplotlib code exercised by
actually generating the figures; here we cover the deterministic helpers
that decide *what* gets drawn:

    1. is_berghia        — leaf-name → Berghia membership.
    2. load_ranking      — CSV → (1-based rank, score) maps.
    3. select_lse_clade  — picks the Berghia-dominated clade holding the
                           most Tier-1 candidates (ties → smallest clade).
"""
from __future__ import annotations

import io

from Bio import Phylo

import visualize_gpcr_tree_rect as vr


def test_is_berghia():
    assert vr.is_berghia("BersteEVm011558t16")
    assert vr.is_berghia("TRINITY_DN1_c0_g1_i1")
    assert not vr.is_berghia("ref_alvmar_JAB_55")
    assert not vr.is_berghia("")
    assert not vr.is_berghia(None)


def test_load_ranking(tmp_path):
    csv = tmp_path / "ranked.csv"
    csv.write_text(
        "id,rank_score\n"
        "BersteEVm011558t16,6.93\n"
        "BersteEVm010564t1,6.33\n"
        "BersteEVm009991t7,6.29\n"
    )
    rank_of, score_of = vr.load_ranking(str(csv))
    assert rank_of["BersteEVm011558t16"] == 1          # 1-based, row order
    assert rank_of["BersteEVm009991t7"] == 3
    assert score_of["BersteEVm010564t1"] == 6.33


def _tree(newick):
    return Phylo.read(io.StringIO(newick), "newick")


def test_select_lse_clade_picks_berghia_clade():
    # Clade A = 4 Berghia tips; clade B = 3 gastropod refs.
    t = _tree("((Berste1,Berste2,Berste3,Berste4),"
              "(ref_alvmar_1,ref_alvmar_2,ref_alvmar_3));")
    tier1 = {"Berste1", "Berste2"}
    clade = vr.select_lse_clade(t, tier1, min_berghia_frac=0.85, min_size=3)
    names = sorted(l.name for l in clade.get_terminals())
    assert names == ["Berste1", "Berste2", "Berste3", "Berste4"]


def test_select_lse_clade_prefers_more_tier1_then_smaller():
    # Two pure-Berghia clades; the one with more Tier-1 wins even if smaller.
    t = _tree("((Berste1,Berste2,Berste3,Berste4,Berste5),"
              "(Berste6,Berste7,Berste8));")
    tier1 = {"Berste6", "Berste7"}            # both Tier-1 live in the 3-tip clade
    clade = vr.select_lse_clade(t, tier1, min_berghia_frac=0.85, min_size=3)
    names = sorted(l.name for l in clade.get_terminals())
    assert names == ["Berste6", "Berste7", "Berste8"]


def test_select_lse_clade_respects_min_berghia_frac():
    # A mostly-non-Berghia clade must be rejected despite holding a Tier-1 tip.
    t = _tree("((Berste1,ref_a_1,ref_a_2,ref_a_3),"
              "(Berste2,Berste3,Berste4,Berste5));")
    tier1 = {"Berste1", "Berste2"}
    clade = vr.select_lse_clade(t, tier1, min_berghia_frac=0.85, min_size=3)
    names = sorted(l.name for l in clade.get_terminals())
    assert names == ["Berste2", "Berste3", "Berste4", "Berste5"]
