"""Tests for the seed-calibrated Apgr family recovery.

The recovery turns three probes into 159. Two things can go silently wrong and
both are checked here on the REAL artifacts, not fixtures:

  * a threshold that drifted away from the seeds it was supposed to be
    calibrated on, and
  * a recovered "new member" that is really a seed under a RefSeq accession,
    which inflates the count and fakes independence while every accession-based
    check passes.
"""
import csv
from pathlib import Path

import pytest

from refexp6_recover_apgr_family import (
    APLYSIA_CALIFORNICA_TAXID,
    ASSEMBLY,
    CALIBRATION_BAR,
    PFAM_GA,
    SEED_VS_SEED_BITSCORES,
    read_fasta,
)

REPO = Path(__file__).resolve().parents[2]
REC_TSV = REPO / "references/anchors/apgr_recovered_members.tsv"
REC_FASTA = REPO / "references/anchors/apgr_recovered_members.fasta"
PROBE_TSV = REPO / "references/anchors/chemoreceptor_probe_set.tsv"
PROBE_FASTA = REPO / "references/anchors/chemoreceptor_probe_set.fasta"

needs_artifacts = pytest.mark.skipif(
    not (REC_TSV.exists() and PROBE_TSV.exists()),
    reason="recovered/probe artifacts absent (references/ is gitignored)",
)


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# --- the calibration is derived from the seeds, not invented ----------------

def test_calibration_bar_is_the_minimum_seed_to_seed_score():
    """The bar must BE a measured seed-seed distance, not a round number.

    If someone replaces it with a hand-picked cutoff this fails, which is the
    point: a threshold chosen independently of the family is not calibrated to
    the family.
    """
    assert CALIBRATION_BAR == min(SEED_VS_SEED_BITSCORES)
    assert CALIBRATION_BAR in SEED_VS_SEED_BITSCORES


def test_seed_scores_are_a_full_all_vs_all():
    """Three seeds give six ordered off-diagonal pairs."""
    assert len(SEED_VS_SEED_BITSCORES) == 6


def test_pfam_gathering_thresholds_are_the_real_ones():
    """Pinned from the Pfam-A.hmm actually searched (HMMER 3.4)."""
    assert PFAM_GA["PF00001"] == 30.5
    assert PFAM_GA["PF10324"] == 22.2
    assert PFAM_GA["PF10328"] == 23.8


def test_calibration_bar_is_far_above_every_pfam_threshold():
    """The bar does real work beyond bare family membership.

    Srw GA alone returns 266 proteins, far more than Cummins's 90, because Srw
    holds peptide receptors too. The bar is what removes them.
    """
    assert CALIBRATION_BAR > max(PFAM_GA.values())


# --- the real recovered set -------------------------------------------------

@needs_artifacts
def test_every_recovered_member_clears_the_calibration_bar():
    rows = read_tsv(REC_TSV)
    assert rows
    for r in rows:
        assert float(r["best_seed_bitscore"]) >= CALIBRATION_BAR, r


@needs_artifacts
def test_every_recovered_member_beats_its_best_peptide_score():
    """The discrimination that separates the target family from the reference
    family it shares a Pfam entry with."""
    for r in read_tsv(REC_TSV):
        assert float(r["best_seed_bitscore"]) > float(r["best_peptide_bitscore"]), r


@needs_artifacts
def test_every_recovered_member_clears_both_family_profiles():
    for r in read_tsv(REC_TSV):
        assert float(r["srw_bitscore"]) >= PFAM_GA["PF10324"], r
        assert float(r["tm7_bitscore"]) >= PFAM_GA["PF00001"], r


@needs_artifacts
def test_every_recovered_member_is_aplysia_californica():
    """Verified against NCBI at build time; re-asserted on the artifact."""
    for r in read_tsv(REC_TSV):
        assert r["taxid"] == APLYSIA_CALIFORNICA_TAXID, r
        assert r["organism"] == "Aplysia californica", r
        assert r["assembly"] == ASSEMBLY, r


@needs_artifacts
def test_recovered_accessions_and_sequences_are_both_unique():
    """Accession uniqueness is not enough -- the same protein can appear twice
    under two accessions and pass every key-based check."""
    rows = read_tsv(REC_TSV)
    accs = [r["refseq_accession"] for r in rows]
    assert len(set(accs)) == len(accs)
    seqs = read_fasta(REC_FASTA)
    assert len(seqs) == len(rows)
    assert len(set(seqs.values())) == len(seqs), "duplicate sequences"


@needs_artifacts
def test_recovered_lengths_match_the_manifest():
    seqs = read_fasta(REC_FASTA)
    for r in read_tsv(REC_TSV):
        assert len(seqs[r["refseq_accession"]]) == int(r["sequence_length"])


@needs_artifacts
def test_no_recovered_member_duplicates_a_published_probe_sequence():
    """The load-bearing check.

    Two AplCal3.0 entries ARE seeds under RefSeq accessions
    (NP_001191494.1 == C5H675, NP_001191495.1 == C5H674). Recovering them as
    "new members" would inflate the count and make correlated evidence look
    independent. Only a sequence comparison catches it.
    """
    rec = read_fasta(REC_FASTA)
    probe = read_fasta(PROBE_FASTA)
    published = {s for pid, s in probe.items()
                 if pid in {"PROBE_A_C5H877", "PROBE_A_C5H675",
                            "PROBE_A_C5H674", "PROBE_A_C6FGJ8"}}
    assert published, "no published probe sequences resolved"
    overlap = [a for a, s in rec.items() if s in published]
    assert not overlap, f"recovered members duplicate a published probe: {overlap}"


@needs_artifacts
def test_the_two_known_refseq_twins_are_absent_from_the_recovered_set():
    accs = {r["refseq_accession"] for r in read_tsv(REC_TSV)}
    assert "NP_001191494.1" not in accs
    assert "NP_001191495.1" not in accs
