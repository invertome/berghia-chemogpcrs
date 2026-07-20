"""Regression guard for the PURIFYING_WEIGHT code default in rank_candidates.py.

Adopted decision (May 2026): PURIFYING_WEIGHT defaults to 0 — chemoreceptor
identification rewards positive/diversifying selection, not whole-gene purifying
selection, so conserved housekeeping GPCRs must not rank highly. config.sh
exports 0, but the CLI's own getenv fallback must also be 0 so a standalone
``python3 rank_candidates.py ...`` (config.sh not sourced) does not silently
reintroduce purifying credit. rank_candidates.py is not import-safe (top-level
``sys.argv``), so this guards the default at the source level.
"""
from __future__ import annotations

import re
from pathlib import Path

RANK = Path(__file__).resolve().parent.parent.parent / "scripts" / "rank_candidates.py"


def test_purifying_weight_default_is_zero() -> None:
    src = RANK.read_text()
    m = re.search(
        r"PURIFYING_WEIGHT\s*=\s*float\(\s*os\.getenv\(\s*['\"]PURIFYING_WEIGHT['\"]\s*,\s*([0-9.]+)\s*\)",
        src,
    )
    assert m is not None, "PURIFYING_WEIGHT getenv default line not found in rank_candidates.py"
    assert float(m.group(1)) == 0.0, (
        f"PURIFYING_WEIGHT code default must be 0 (adopted decision), got {m.group(1)}"
    )


def test_production_wrapper_forwards_dnds_reliability() -> None:
    """The production scorer wrapper calculate_fair_rank_score(row) must forward
    the row's per-OG dN/dS reliability weight into the fair-rank scorer (bead
    -8st). MEDIUM-3 was exactly this wiring being absent. rank_candidates.py is
    not import-safe, so this guards the wiring at the source level."""
    src = RANK.read_text()
    start = src.find("def calculate_fair_rank_score(row):")
    assert start != -1, "production wrapper calculate_fair_rank_score(row) not found"
    end = src.find("\n\ndef ", start + 1)
    body = src[start:end] if end != -1 else src[start:]
    assert "dnds_reliability_weight" in body, "wrapper does not read the reliability weight"
    assert "dnds_reliability=" in body, "wrapper does not pass dnds_reliability to the scorer"


def test_phylo_and_lse_depth_gated_by_tree_membership() -> None:
    """o98: phylo/lse_depth are class-A chemoreceptor-tree signals, so they must
    be dropped (None) for candidates not in that tree (class B/C/F, unclassified,
    or too divergent to be pooled) — those candidates then fall out of
    fair-scoring instead of carrying a full-weight present-zero that drags their
    composite down. rank_candidates.py is not import-safe; guard the wiring here.
    (Berghia GPCR chemoreceptors are class-A rhodopsin-like; class-C
    chemoreceptors are a vertebrate-only innovation, so class A is the only
    scored chemoreceptor tree.)"""
    src = RANK.read_text()
    assert "'has_phylo_data'" in src, \
        "row dict must record class-A-tree membership as has_phylo_data"
    start = src.find("def calculate_fair_rank_score(row):")
    assert start != -1
    end = src.find("\n\ndef ", start + 1)
    body = src[start:end] if end != -1 else src[start:]
    assert "'phylo': row.get('phylo_score_norm') if row.get('has_phylo_data') else None" in body
    assert "'lse_depth': row.get('lse_depth_score_norm') if row.get('has_phylo_data') else None" in body
