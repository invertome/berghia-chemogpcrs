"""Tests for scripts/classify_consensus.py.

Phase 4 Task 4.4 — combines HMM scan (Task 4.1), OG vote (Task 4.2),
and phylogenetic placement (Task 4.3) into a single per-candidate
classification with confidence stratification.

Consensus logic:
  - 3-of-3 family agreement (HMM + OG + placement)
      → 'non-chemoreceptor' (high confidence)
  - 2-of-3 family agreement (HMM + OG agree; placement disagrees or absent)
      → 'likely-non-chemoreceptor' (medium confidence)
  - <2 agreement
      → 'chemoreceptor-candidate' (default, no override)

Subfamily emitted only when ALL agreeing sources also agree on subfamily.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

import classify_consensus as cc


def _make_inputs(tmp_path: Path,
                 hmm_rows: list[dict],
                 og_rows: list[dict],
                 placement_rows: list[dict] | None = None
                 ) -> tuple[str, str, str | None]:
    """Build the three input TSVs and return paths."""
    hmm_path = tmp_path / "hmm.tsv"
    pd.DataFrame(hmm_rows).to_csv(hmm_path, sep="\t", index=False)
    og_path = tmp_path / "og.tsv"
    pd.DataFrame(og_rows).to_csv(og_path, sep="\t", index=False)
    if placement_rows is None:
        return str(hmm_path), str(og_path), None
    placement_path = tmp_path / "placement.tsv"
    pd.DataFrame(placement_rows).to_csv(placement_path, sep="\t", index=False)
    return str(hmm_path), str(og_path), str(placement_path)


# ---- 3-of-3 agreement ---------------------------------------------------

def test_three_way_agreement_high_confidence(tmp_path: Path) -> None:
    """All three sources agree on aminergic/5HT → high-confidence non-chemo."""
    hmm = [{"candidate_id": "g1", "hmm_family": "aminergic",
            "hmm_subfamily": "5HT", "evalue": "1e-150",
            "score": "500", "hmm_target": "aminergic_5HT"}]
    og = [{"candidate_id": "g1", "og_id": "OG0001",
           "og_vote_family": "aminergic", "og_vote_subfamily": "5HT",
           "n_members": 5, "n_annotated": 5, "consensus_fraction": "1.000"}]
    placement = [{"candidate_id": "g1", "placement_family": "aminergic",
                  "placement_subfamily": "5HT", "placement_lwr": "0.95"}]
    h, o, p = _make_inputs(tmp_path, hmm, og, placement)
    out = tmp_path / "consensus.tsv"
    cc.run_consensus(h, o, p, str(out))
    df = pd.read_csv(out, sep="\t", keep_default_na=False, dtype=str)
    row = df[df["candidate_id"] == "g1"].iloc[0]
    assert row["classification"] == "non-chemoreceptor"
    assert row["classification_confidence"] == "high"
    assert row["classification_family"] == "aminergic"
    assert row["classification_subfamily"] == "5HT"


# ---- 2-of-3 agreement (HMM + OG, no placement support) -----------------

def test_two_way_agreement_no_placement_medium_confidence(tmp_path: Path) -> None:
    """HMM + OG agree; placement disagrees → likely-non-chemoreceptor."""
    hmm = [{"candidate_id": "g1", "hmm_family": "aminergic",
            "hmm_subfamily": "5HT", "evalue": "1e-150", "score": "500",
            "hmm_target": "aminergic_5HT"}]
    og = [{"candidate_id": "g1", "og_id": "OG0001",
           "og_vote_family": "aminergic", "og_vote_subfamily": "5HT",
           "n_members": 5, "n_annotated": 5, "consensus_fraction": "1.000"}]
    placement = [{"candidate_id": "g1", "placement_family": "peptide",
                  "placement_subfamily": "tachykinin", "placement_lwr": "0.85"}]
    h, o, p = _make_inputs(tmp_path, hmm, og, placement)
    out = tmp_path / "consensus.tsv"
    cc.run_consensus(h, o, p, str(out))
    df = pd.read_csv(out, sep="\t", keep_default_na=False, dtype=str)
    row = df[df["candidate_id"] == "g1"].iloc[0]
    assert row["classification"] == "likely-non-chemoreceptor"
    assert row["classification_confidence"] == "medium"
    assert row["classification_family"] == "aminergic"


def test_two_way_agreement_placement_absent_medium_confidence(tmp_path: Path) -> None:
    """HMM + OG agree; placement file not provided → likely-non-chemo."""
    hmm = [{"candidate_id": "g1", "hmm_family": "peptide",
            "hmm_subfamily": "tachykinin", "evalue": "1e-100", "score": "350",
            "hmm_target": "peptide_tachykinin"}]
    og = [{"candidate_id": "g1", "og_id": "OG0001",
           "og_vote_family": "peptide", "og_vote_subfamily": "tachykinin",
           "n_members": 5, "n_annotated": 5, "consensus_fraction": "1.000"}]
    h, o, _ = _make_inputs(tmp_path, hmm, og, None)
    out = tmp_path / "consensus.tsv"
    cc.run_consensus(h, o, None, str(out))
    df = pd.read_csv(out, sep="\t", keep_default_na=False, dtype=str)
    row = df[df["candidate_id"] == "g1"].iloc[0]
    assert row["classification"] == "likely-non-chemoreceptor"
    assert row["classification_confidence"] == "medium"


# ---- 1-of-3 agreement (default to chemoreceptor-candidate) -------------

def test_only_hmm_agreement_chemoreceptor_default(tmp_path: Path) -> None:
    """Only HMM classifies; OG can't (insufficient annotated members) →
    candidate stays as chemoreceptor-candidate."""
    hmm = [{"candidate_id": "g1", "hmm_family": "aminergic",
            "hmm_subfamily": "5HT", "evalue": "1e-150", "score": "500",
            "hmm_target": "aminergic_5HT"}]
    og = [{"candidate_id": "g1", "og_id": "OG0001",
           "og_vote_family": "unclassified-og", "og_vote_subfamily": "",
           "n_members": 2, "n_annotated": 1, "consensus_fraction": "0.0"}]
    placement = [{"candidate_id": "g1", "placement_family": "unclassified-placement",
                  "placement_subfamily": "", "placement_lwr": "0.0"}]
    h, o, p = _make_inputs(tmp_path, hmm, og, placement)
    out = tmp_path / "consensus.tsv"
    cc.run_consensus(h, o, p, str(out))
    df = pd.read_csv(out, sep="\t", keep_default_na=False, dtype=str)
    row = df[df["candidate_id"] == "g1"].iloc[0]
    assert row["classification"] == "chemoreceptor-candidate"
    assert row["classification_confidence"] == "NA"


def test_no_classifier_agreement_chemoreceptor_default(tmp_path: Path) -> None:
    """All three sources disagree → chemoreceptor-candidate."""
    hmm = [{"candidate_id": "g1", "hmm_family": "aminergic",
            "hmm_subfamily": "5HT", "evalue": "1e-150", "score": "500",
            "hmm_target": "aminergic_5HT"}]
    og = [{"candidate_id": "g1", "og_id": "OG0001",
           "og_vote_family": "peptide", "og_vote_subfamily": "tachykinin",
           "n_members": 5, "n_annotated": 5, "consensus_fraction": "1.000"}]
    placement = [{"candidate_id": "g1", "placement_family": "opsin",
                  "placement_subfamily": "", "placement_lwr": "0.85"}]
    h, o, p = _make_inputs(tmp_path, hmm, og, placement)
    out = tmp_path / "consensus.tsv"
    cc.run_consensus(h, o, p, str(out))
    df = pd.read_csv(out, sep="\t", keep_default_na=False, dtype=str)
    row = df[df["candidate_id"] == "g1"].iloc[0]
    assert row["classification"] == "chemoreceptor-candidate"


# ---- subfamily handling -------------------------------------------------

def test_subfamily_consensus_when_all_agree(tmp_path: Path) -> None:
    """All sources agree on family AND subfamily → emit subfamily."""
    hmm = [{"candidate_id": "g1", "hmm_family": "aminergic",
            "hmm_subfamily": "5HT", "evalue": "1e-150", "score": "500",
            "hmm_target": "aminergic_5HT"}]
    og = [{"candidate_id": "g1", "og_id": "OG0001",
           "og_vote_family": "aminergic", "og_vote_subfamily": "5HT",
           "n_members": 5, "n_annotated": 5, "consensus_fraction": "1.000"}]
    placement = [{"candidate_id": "g1", "placement_family": "aminergic",
                  "placement_subfamily": "5HT", "placement_lwr": "0.95"}]
    h, o, p = _make_inputs(tmp_path, hmm, og, placement)
    out = tmp_path / "consensus.tsv"
    cc.run_consensus(h, o, p, str(out))
    df = pd.read_csv(out, sep="\t", keep_default_na=False, dtype=str)
    assert df.iloc[0]["classification_subfamily"] == "5HT"


def test_subfamily_disagreement_falls_back_to_empty(tmp_path: Path) -> None:
    """Family agrees (aminergic) but subfamilies disagree (5HT vs dopamine)
    → emit family only, subfamily empty."""
    hmm = [{"candidate_id": "g1", "hmm_family": "aminergic",
            "hmm_subfamily": "5HT", "evalue": "1e-150", "score": "500",
            "hmm_target": "aminergic_5HT"}]
    og = [{"candidate_id": "g1", "og_id": "OG0001",
           "og_vote_family": "aminergic", "og_vote_subfamily": "dopamine",
           "n_members": 5, "n_annotated": 5, "consensus_fraction": "1.000"}]
    placement = [{"candidate_id": "g1", "placement_family": "aminergic",
                  "placement_subfamily": "", "placement_lwr": "0.85"}]
    h, o, p = _make_inputs(tmp_path, hmm, og, placement)
    out = tmp_path / "consensus.tsv"
    cc.run_consensus(h, o, p, str(out))
    df = pd.read_csv(out, sep="\t", keep_default_na=False, dtype=str)
    row = df.iloc[0]
    assert row["classification_family"] == "aminergic"
    # subfamilies disagreed (5HT vs dopamine vs missing) → empty
    assert row["classification_subfamily"] == ""


# ---- evidence string ----------------------------------------------------

def test_evidence_string_format(tmp_path: Path) -> None:
    """Evidence column shows source:label per source, semicolon-separated."""
    hmm = [{"candidate_id": "g1", "hmm_family": "aminergic",
            "hmm_subfamily": "5HT", "evalue": "1e-150", "score": "500",
            "hmm_target": "aminergic_5HT"}]
    og = [{"candidate_id": "g1", "og_id": "OG0001",
           "og_vote_family": "aminergic", "og_vote_subfamily": "5HT",
           "n_members": 5, "n_annotated": 5, "consensus_fraction": "1.000"}]
    placement = [{"candidate_id": "g1", "placement_family": "aminergic",
                  "placement_subfamily": "5HT", "placement_lwr": "0.95"}]
    h, o, p = _make_inputs(tmp_path, hmm, og, placement)
    out = tmp_path / "consensus.tsv"
    cc.run_consensus(h, o, p, str(out))
    df = pd.read_csv(out, sep="\t", keep_default_na=False, dtype=str)
    ev = df.iloc[0]["classification_evidence"]
    assert "hmm:aminergic" in ev
    assert "og:aminergic" in ev
    assert "placement:aminergic" in ev


# ---- candidate present in HMM only ------------------------------------

def test_candidate_in_hmm_but_not_og(tmp_path: Path) -> None:
    """Some candidates are in HMM TSV but not OG TSV (e.g. orphan with no OG).
    Should still produce a row, default to chemoreceptor-candidate."""
    hmm = [{"candidate_id": "orphan", "hmm_family": "aminergic",
            "hmm_subfamily": "", "evalue": "1e-150", "score": "500",
            "hmm_target": "aminergic"}]
    og: list[dict] = []
    h, o, _ = _make_inputs(tmp_path, hmm, og, None)
    out = tmp_path / "consensus.tsv"
    cc.run_consensus(h, o, None, str(out))
    df = pd.read_csv(out, sep="\t", keep_default_na=False, dtype=str)
    row = df[df["candidate_id"] == "orphan"].iloc[0]
    # Only HMM signal → not enough for non-chemo call
    assert row["classification"] == "chemoreceptor-candidate"
