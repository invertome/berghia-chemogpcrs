"""Pins the corrected exclusion standard in the family ratification.

Q9NJC9 was excluded on its protein name ("RHO G-protein coupled receptor"
commits to no family) while every ADMISSION was required to show a curated
InterPro family entry. Two standards, and the looser one was wrong: the record
carries better opsin evidence than any of the nine opsins that were admitted.

These tests exist so the ruling cannot silently revert and so the exclusions
are pinned to the evidence that actually justifies them rather than to a name.
"""
import csv
from pathlib import Path

import pytest

from refexp4_ratify_held_families import (
    ADMISSION_NOTES,
    ADMIT,
    EXCLUDE,
    OPSIN_SIGNATURES,
    PEPTIDE_SIGNATURES,
    PROBE,
)

REPO = Path(__file__).resolve().parents[2]
ANCHORS = REPO / "references/anchors/anchor_set_PROD.tsv"
DISPOSITION = REPO / "references/anchors/anchor_set_PROD_held_disposition.tsv"

needs_artifacts = pytest.mark.skipif(
    not ANCHORS.exists(),
    reason="anchor table absent (references/ is gitignored)")


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# --- the reversal ----------------------------------------------------------

def test_q9njc9_is_admitted_as_opsin():
    assert "Q9NJC9" in ADMIT
    family, organism, name, required = ADMIT["Q9NJC9"]
    assert family == "opsin"
    assert organism == "Schistosoma mansoni"
    assert required is OPSIN_SIGNATURES


def test_q9njc9_is_no_longer_excluded():
    assert "Q9NJC9" not in EXCLUDE


def test_the_reversal_is_recorded_rather_than_quietly_applied():
    """A changed ruling that leaves no trace is indistinguishable from one that
    was always the call."""
    note = ADMISSION_NOTES.get("Q9NJC9", "")
    assert "REVERSING" in note
    assert "IPR050125" in note and "IPR027430" in note


def test_the_four_remaining_exclusions_are_unchanged():
    assert set(EXCLUDE) == {"Q962I3", "C9K4W2", "F2VWU2", "C0M0N9"}


def test_every_held_entry_still_has_exactly_one_ruling():
    """22 held entries, each accounted for once. Overlapping rulings would let
    an entry be admitted and excluded at the same time."""
    assert len(ADMIT) + len(EXCLUDE) + len(PROBE) == 22
    assert not set(ADMIT) & set(EXCLUDE)
    assert not set(ADMIT) & set(PROBE)
    assert not set(EXCLUDE) & set(PROBE)


def test_generic_rhodopsin_entries_are_accepted_nowhere():
    """IPR000276/IPR017452 are on all 22 held records, so accepting either
    would turn the curated check into a rubber stamp."""
    for sigs in (OPSIN_SIGNATURES, PEPTIDE_SIGNATURES):
        assert "IPR000276" not in sigs
        assert "IPR017452" not in sigs


# --- the audited exclusions, pinned to evidence not names -------------------
#
# Re-audited against live UniProt + the InterPro API on 2026-07-21. Entry TYPES
# were resolved from InterPro, not assumed. None of the four carries an
# InterPro FAMILY-type entry beyond the generic pair, none has a SIMILARITY
# statement, and C9K4W2 has no Pfam hit at all.
AUDITED_EXCLUSIONS = {
    "Q962I3":  {"pfam": "PF00001", "family_signatures": 0, "similarity": False},
    "C9K4W2":  {"pfam": "",        "family_signatures": 0, "similarity": False},
    "F2VWU2":  {"pfam": "PF00001", "family_signatures": 0, "similarity": False},
    "C0M0N9":  {"pfam": "PF00001", "family_signatures": 0, "similarity": False},
}


def test_audited_exclusions_cover_exactly_the_remaining_four():
    assert set(AUDITED_EXCLUSIONS) == set(EXCLUDE)


def test_no_audited_exclusion_had_family_level_evidence():
    """This is what justifies the exclusions now -- not the protein name."""
    for acc, found in AUDITED_EXCLUSIONS.items():
        assert found["family_signatures"] == 0, acc


@needs_artifacts
def test_no_excluded_accession_reached_the_reference_set():
    accs = {r["accession"] for r in read_tsv(ANCHORS)}
    assert not set(EXCLUDE) & accs


@needs_artifacts
def test_q9njc9_is_in_the_reference_set_as_opsin():
    rows = [r for r in read_tsv(ANCHORS) if r["accession"] == "Q9NJC9"]
    assert len(rows) == 1
    assert rows[0]["family"] == "opsin"
    assert rows[0]["class"] == "A"
    assert rows[0]["taxid"] == "6183"


@needs_artifacts
def test_no_probe_accession_reached_the_reference_set():
    """The probe entries are target class and must never be admitted."""
    accs = {r["accession"] for r in read_tsv(ANCHORS)}
    assert not set(PROBE) & accs


@needs_artifacts
def test_disposition_records_q9njc9_as_admitted():
    if not DISPOSITION.exists():
        pytest.skip("disposition table absent")
    rows = {r["accession"]: r for r in read_tsv(DISPOSITION)}
    r = rows["Q9NJC9"]
    assert r["disposition"] == "admitted"
    assert r["assigned_family"] == "opsin"
    assert "IPR050125" in r["family_basis"]
    assert "REVERSING" in r["reason"]
