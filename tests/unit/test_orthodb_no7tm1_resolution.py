"""The 7tm_1-miss resolver, and the column-label trap it walked into.

The resolver decides whether a class-A anchor that does not hit PF00001 is a
genuinely divergent class-A receptor (ACKR1/DARC) or class-B/C contamination.
It reads UniProt's curated family statement to decide.

UniProt's TSV header labels are NOT the field names requested: `cc_similarity`
comes back as "Sequence similarities". Read under any other label it yields an
empty string for EVERY record, which silently demotes the resolver to a
Pfam-only test -- exactly the domain-only criterion it exists to replace, and
one that resolves every 7tm_1 miss as "not class A". Observed live: all five
surviving anchors, including both ACKR1 entries, resolved to UNKNOWN and would
have been evicted.

That failure is invisible in a fixture keyed on the assumed label, so the tests
below pin the label itself and the guard that asserts it.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

import orthodb_resolve_no7tm1 as resolver  # noqa: E402
from curate_gpcr_references import gpcr_class_from_evidence  # noqa: E402


def test_curated_similarity_is_read_under_uniprots_actual_header_label():
    """The label is the defect. Pin it explicitly."""
    assert resolver.COLUMN_LABELS["curated_similarity"] == "Sequence similarities"


def test_requested_field_name_is_not_the_header_label():
    """Guard against 'fixing' this by assuming request name == header label."""
    assert "cc_similarity" in resolver.FIELDS
    assert resolver.COLUMN_LABELS["curated_similarity"] != "cc_similarity"


def test_every_declared_label_is_requested_as_a_field():
    """A label with no corresponding requested field can never be populated."""
    request_to_label = {
        "accession": "Entry",
        "cc_similarity": "Sequence similarities",
        "xref_pfam": "Pfam",
        "reviewed": "Reviewed",
        "protein_name": "Protein names",
        "organism_name": "Organism",
        "organism_id": "Organism (ID)",
        "length": "Length",
    }
    for field, label in request_to_label.items():
        assert field in resolver.FIELDS, f"{field} is not requested"
        assert label in resolver.COLUMN_LABELS.values(), f"{label} is not mapped"


def test_missing_column_raises_rather_than_resolving_on_empty_evidence(monkeypatch):
    """A header without the curated column must abort, not resolve Pfam-only.

    This is the whole point: the broken run produced a complete, plausible
    verdict table. Silence is the failure mode, so absence must be loud.
    """
    header = "\t".join(["Entry", "Pfam", "Reviewed"])
    body = "\t".join(["Q8IWP5", "PF00002;", "reviewed"])
    payload = f"{header}\n{body}\n"

    class FakeResponse:
        def read(self):
            return payload.encode()

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    monkeypatch.setattr(resolver, "urlopen", lambda *a, **k: FakeResponse())

    with pytest.raises(SystemExit) as excinfo:
        resolver.fetch(["Q8IWP5"])
    assert "Sequence similarities" in str(excinfo.value)


def test_ackr1_resolves_to_class_a_from_curation_despite_no_pfam_signature():
    """The known divergent case, with its real curated statement.

    ACKR1/DARC has no PF00001 hit. Curation still places it in family 1, and
    curation is authoritative, so it must be RETAINED.
    """
    curated = ("SIMILARITY: Belongs to the G-protein coupled receptor 1 family. "
               "Atypical chemokine receptor subfamily. "
               "{ECO:0000256|ARBA:ARBA00008790}.")
    assert gpcr_class_from_evidence(curated, "") == "A"


def test_a_7tm1_miss_with_no_curation_does_not_default_to_class_a():
    """No evidence must resolve to UNKNOWN, never to a silent 'A'."""
    assert gpcr_class_from_evidence("", "") == "UNKNOWN"


def test_curation_outranks_signature_in_both_directions():
    """The old test erred both ways; the resolver must not."""
    b_curated = "SIMILARITY: Belongs to the G protein-coupled receptor 2 family."
    assert gpcr_class_from_evidence(b_curated, "PF00001;") == "B"
    a_curated = "SIMILARITY: Belongs to the G-protein coupled receptor 1 family."
    assert gpcr_class_from_evidence(a_curated, "PF00002;") == "A"
