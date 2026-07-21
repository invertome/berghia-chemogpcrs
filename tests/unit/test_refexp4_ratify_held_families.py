"""Tests for the ratification of the held ``unclassified`` reference entries.

A family label is a novelty prototype, so assigning one is a scientific call
with a measurable effect on every candidate's score. The rule this module
enforces is that the call rests on a curated InterPro FAMILY entry, never on the
protein name -- a name is an author's label, and the reference set has already
been corrupted once by a name-driven assignment (the 'gonadotropin' substring
collision that filed two vasopressin/oxytocin receptors as glycoprotein-hormone).
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
    verify_admissions,
)

REPO = Path(__file__).resolve().parents[2]
HELD_TSV = REPO / "references/anchors/anchor_set_PROD_held_pending_family.tsv"

# The generic rhodopsin-like 7TM entries. EVERY one of the 22 held entries
# carries both, so neither discriminates anything.
GENERIC_7TM = {"IPR000276", "IPR017452"}

needs_held = pytest.mark.skipif(
    not HELD_TSV.exists(),
    reason="held table absent (references/ is gitignored)",
)


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# --- the ruling is a partition ----------------------------------------------

def test_the_three_buckets_are_disjoint():
    assert not set(ADMIT) & set(EXCLUDE)
    assert not set(ADMIT) & set(PROBE)
    assert not set(EXCLUDE) & set(PROBE)


def test_bucket_sizes_are_pinned():
    """15/4/3 since 2026-07-21: Q9NJC9 moved from EXCLUDE to ADMIT when the
    exclusions were re-audited under the curated standard the admissions face.
    See tests/unit/test_refexp4_exclusion_audit.py."""
    assert (len(ADMIT), len(EXCLUDE), len(PROBE)) == (15, 4, 3)


@needs_held
def test_every_held_entry_has_exactly_one_ruling():
    """A forgotten held entry would stay out with no record of why."""
    held = {r["accession"] for r in read_tsv(HELD_TSV)}
    ruled = set(ADMIT) | set(EXCLUDE) | set(PROBE)
    assert held == ruled


def test_the_aplysia_chemosensory_receptors_are_probes_not_admissions():
    """The load-bearing ruling: target class never becomes a reference."""
    for acc in ("C5H877", "C5H675", "C5H674"):
        assert acc in PROBE
        assert acc not in ADMIT


# --- the family assignment rests on a curated signature ----------------------

def test_no_admission_accepts_a_generic_7tm_signature():
    """Accepting IPR000276 would make the check a rubber stamp."""
    for acc, (_, _, _, required) in ADMIT.items():
        assert not (required & GENERIC_7TM), acc
        assert required


def test_admitted_families_are_only_the_two_ratified_ones():
    assert {f for f, _, _, _ in ADMIT.values()} == {"opsin", "peptide"}


def test_each_family_requires_its_own_signature_set():
    for acc, (family, _, _, required) in ADMIT.items():
        expected = OPSIN_SIGNATURES if family == "opsin" else PEPTIDE_SIGNATURES
        assert required == expected, acc


# --- verify_admissions refuses every way the assignment can be wrong ---------

def _record(**over):
    rec = {
        "accession": "A0A0H5ANL2", "protein_name": "Retinochrome1",
        "organism": "Idiosepius paradoxus", "taxid": "294707", "length": 301,
        "interpro": frozenset({"IPR050125", "IPR000276"}),
    }
    rec.update(over)
    return rec


def _held(**over):
    row = {"accession": "A0A0H5ANL2", "sequence_length": "301",
           "evidence_tier": "published-not-deorphanized", "phylum": "Mollusca",
           "organism": "Idiosepius paradoxus"}
    row.update(over)
    return row


@pytest.fixture
def one_admission(monkeypatch):
    monkeypatch.setattr(
        "refexp4_ratify_held_families.ADMIT",
        {"A0A0H5ANL2": ("opsin", "Idiosepius paradoxus", "Retinochrome1",
                        OPSIN_SIGNATURES)},
    )


def test_a_correct_admission_verifies(one_admission):
    verified, errors = verify_admissions(
        {"A0A0H5ANL2": _record()}, {"A0A0H5ANL2": _held()},
        {"A0A0H5ANL2": "M" * 301})
    assert not errors
    assert len(verified) == 1
    assert verified[0]["interpro_hits"] == ["IPR050125"]


def test_a_wrong_organism_is_refused(one_admission):
    """A resolving accession is not proof it is the right record."""
    verified, errors = verify_admissions(
        {"A0A0H5ANL2": _record(organism="Sepia officinalis")},
        {"A0A0H5ANL2": _held()}, {"A0A0H5ANL2": "M" * 301})
    assert not verified
    assert "organism" in errors[0]


def test_a_wrong_protein_name_is_refused(one_admission):
    verified, errors = verify_admissions(
        {"A0A0H5ANL2": _record(protein_name="Rhodopsin")},
        {"A0A0H5ANL2": _held()}, {"A0A0H5ANL2": "M" * 301})
    assert not verified
    assert "protein name" in errors[0]


def test_a_name_only_assignment_is_refused(one_admission):
    """The core rule: the record carries only the generic 7TM entries.

    The protein is still called 'Retinochrome1' -- which is why assigning from
    the name would have let this through.
    """
    verified, errors = verify_admissions(
        {"A0A0H5ANL2": _record(interpro=frozenset(GENERIC_7TM))},
        {"A0A0H5ANL2": _held()}, {"A0A0H5ANL2": "M" * 301})
    assert not verified
    assert "NOT corroborated" in errors[0]


def test_a_contradicting_family_signature_is_refused(one_admission):
    """An opsin admission whose record carries only a peptide-family entry."""
    verified, errors = verify_admissions(
        {"A0A0H5ANL2": _record(interpro=frozenset({"IPR019427", "IPR000276"}))},
        {"A0A0H5ANL2": _held()}, {"A0A0H5ANL2": "M" * 301})
    assert not verified
    assert "NOT corroborated" in errors[0]


def test_a_length_drift_against_the_held_table_is_refused(one_admission):
    verified, errors = verify_admissions(
        {"A0A0H5ANL2": _record(length=302)}, {"A0A0H5ANL2": _held()},
        {"A0A0H5ANL2": "M" * 302})
    assert not verified
    assert "held length" in errors[0]


def test_a_truncated_or_doubled_sequence_is_refused(one_admission):
    """An accumulate-by-accession bug has silently doubled sequences before."""
    verified, errors = verify_admissions(
        {"A0A0H5ANL2": _record()}, {"A0A0H5ANL2": _held()},
        {"A0A0H5ANL2": "M" * 602})
    assert not verified
    assert "FASTA sequence length" in errors[0]


def test_a_missing_uniprot_record_is_refused(one_admission):
    verified, errors = verify_admissions(
        {}, {"A0A0H5ANL2": _held()}, {"A0A0H5ANL2": "M" * 301})
    assert not verified
    assert "no record" in errors[0]


# --- the exclusions ----------------------------------------------------------

def test_the_opsin_evidence_on_Q9NJC9_is_still_recorded_after_its_admission():
    """The live question this test was written to keep open has been ANSWERED.

    Q9NJC9 used to be excluded, and this test asserted its opsin evidence was
    recorded in the exclusion reason so the call stayed reviewable. Reviewing
    it is exactly what happened: re-audited under the curated standard the
    admissions face, the record carries IPR050125 plus the IPR027430 retinal
    binding site that none of the nine admitted opsins has, and it was admitted
    on 2026-07-21.

    The requirement is unchanged in substance -- the evidence must still be
    written down where a reader will find it -- so the assertion now points at
    the admission note instead of the exclusion reason.
    """
    assert "Q9NJC9" not in EXCLUDE
    note = ADMISSION_NOTES["Q9NJC9"]
    assert "IPR050125" in note and "IPR027430" in note
    assert "REVERSING" in note


def test_every_exclusion_states_a_reason():
    for acc, (reason, name) in EXCLUDE.items():
        assert len(reason) > 30, acc
        assert name
