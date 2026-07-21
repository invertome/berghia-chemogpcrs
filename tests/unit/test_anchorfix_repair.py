"""Tests for the anchor-set class/taxid repair (``anchorfix_class_reconcile``).

The repair edits ``anchor_set_PROD.tsv`` in place. Three things make that
dangerous, and each has a test here:

* **Identity.** The accession is the identity, and the derived FASTA/npz key is
  ``ANCHOR_<class>_<tier>_<accession>``. Rewriting ``class`` would re-mint that
  key, so the repair EVICTS non-conforming rows instead of relabelling them, and
  every retained row must keep its composite id byte-identical.
* **Silent doubling.** Accumulating by accession has previously duplicated
  sequences while passing every count check, so integrity is asserted on the
  sequence, not on the row count.
* **Drift.** A correction keyed only on accession would fire against a file that
  has since changed underneath it, so each correction also asserts the current
  value it expects to replace.
"""
import csv
from pathlib import Path

import pytest

from anchorfix_class_reconcile import (
    apply_family_corrections,
    composite_id,
    fill_taxids,
    partition_by_resolved_class,
    verify_written_table,
    write_table_atomic,
)

_COLS = ["accession", "tier", "taxid", "species", "family", "class", "evidence"]


def _row(acc, tier="1", taxid="", species="", family="peptide", klass="A",
         evidence="reviewed"):
    return dict(zip(_COLS, [acc, tier, taxid, species, family, klass, evidence]))


# --- eviction preserves the identity of everything it retains ----------------

def test_partition_evicts_only_rows_whose_resolved_class_differs():
    rows = [_row("P1"), _row("P2"), _row("P3")]
    resolved = {"P1": "A", "P2": "B", "P3": "A"}
    keep, evict = partition_by_resolved_class(rows, resolved, envelope="A")
    assert [r["accession"] for r in keep] == ["P1", "P3"]
    assert [r["accession"] for r in evict] == ["P2"]


def test_unknown_class_is_evicted_not_kept():
    """UNKNOWN must never be treated as 'probably A'."""
    rows = [_row("P1"), _row("P2")]
    keep, evict = partition_by_resolved_class(
        rows, {"P1": "A", "P2": "UNKNOWN"}, envelope="A")
    assert [r["accession"] for r in evict] == ["P2"]


def test_rows_outside_the_envelope_are_untouched():
    """Repairing the class-A envelope must not disturb the B/C/F rows."""
    rows = [_row("P1", klass="A"), _row("P2", klass="B"), _row("P3", klass="C")]
    keep, evict = partition_by_resolved_class(
        rows, {"P1": "A", "P2": "C", "P3": "C"}, envelope="A")
    assert [r["accession"] for r in keep] == ["P1", "P2", "P3"]
    assert evict == []


def test_retained_rows_keep_composite_id_and_order():
    rows = [_row("P1", tier="3"), _row("P2"), _row("P3", tier="8")]
    before = [composite_id(r) for r in rows]
    keep, _ = partition_by_resolved_class(
        rows, {"P1": "A", "P2": "B", "P3": "A"}, envelope="A")
    assert [composite_id(r) for r in keep] == [before[0], before[2]]
    assert composite_id(keep[0]) == "ANCHOR_A_3_P1"


def test_partition_is_idempotent():
    """A second pass over an already-partitioned table evicts nothing more."""
    rows = [_row("P1"), _row("P2"), _row("P3")]
    resolved = {"P1": "A", "P2": "B", "P3": "A"}
    keep, evict = partition_by_resolved_class(rows, resolved, envelope="A")
    keep2, evict2 = partition_by_resolved_class(keep, resolved, envelope="A")
    assert evict2 == []
    assert [r["accession"] for r in keep2] == [r["accession"] for r in keep]


def test_fill_taxids_is_idempotent():
    rows = [_row("P1", taxid=""), _row("P2", taxid="")]
    src = {"P1": "6523", "P2": "9606"}
    assert fill_taxids(rows, src) == 2
    assert fill_taxids(rows, src) == 0
    assert [r["taxid"] for r in rows] == ["6523", "9606"]


def test_partition_raises_when_a_row_has_no_resolved_class():
    """A missing lookup is a broken join, not an implicit keep."""
    rows = [_row("P1"), _row("P2")]
    with pytest.raises(ValueError, match="P2"):
        partition_by_resolved_class(rows, {"P1": "A"}, envelope="A")


# --- taxid backfill ----------------------------------------------------------

def test_fill_taxids_populates_blanks_only():
    rows = [_row("P1", taxid=""), _row("P2", taxid="9606")]
    n = fill_taxids(rows, {"P1": "6523", "P2": "9606"})
    assert n == 1
    assert [r["taxid"] for r in rows] == ["6523", "9606"]


def test_fill_taxids_raises_on_disagreement_with_source():
    """An existing taxid that contradicts UniProt is a data defect; the repair
    must stop rather than silently overwrite or silently keep it."""
    rows = [_row("P1", taxid="9606")]
    with pytest.raises(ValueError, match="P1"):
        fill_taxids(rows, {"P1": "10090"})


def test_fill_taxids_raises_when_source_lacks_the_accession():
    rows = [_row("P1", taxid="")]
    with pytest.raises(ValueError, match="P1"):
        fill_taxids(rows, {})


def test_fill_taxids_never_leaves_a_blank():
    rows = [_row(a, taxid="") for a in ("P1", "P2", "P3")]
    fill_taxids(rows, {"P1": "1", "P2": "2", "P3": "3"})
    assert all(r["taxid"] for r in rows)


# --- family corrections ------------------------------------------------------

def test_family_correction_applies_and_leaves_class_alone():
    rows = [_row("Q75W84", family="orphan", klass="A")]
    applied = apply_family_corrections(
        rows, {"Q75W84": ("orphan", "peptide")})
    assert applied == 1
    assert rows[0]["family"] == "peptide"
    assert rows[0]["class"] == "A"
    assert composite_id(rows[0]) == "ANCHOR_A_1_Q75W84"


def test_family_correction_is_idempotent():
    """Re-running the repair against an already-repaired file must converge.
    The file is unversioned (references/** is gitignored), so a second run is a
    realistic recovery action after an interrupted first run."""
    rows = [_row("Q75W84", family="orphan")]
    corr = {"Q75W84": ("orphan", "peptide")}
    assert apply_family_corrections(rows, corr) == 1
    assert apply_family_corrections(rows, corr) == 0
    assert rows[0]["family"] == "peptide"


def test_family_correction_refuses_when_current_value_drifted():
    rows = [_row("Q75W84", family="aminergic")]
    with pytest.raises(ValueError, match="Q75W84"):
        apply_family_corrections(rows, {"Q75W84": ("orphan", "peptide")})


def test_family_correction_refuses_unknown_accession():
    rows = [_row("P1", family="orphan")]
    with pytest.raises(ValueError, match="ZZZZZZ"):
        apply_family_corrections(rows, {"ZZZZZZ": ("orphan", "peptide")})


def test_corrections_never_change_class_so_no_key_moves():
    """Every shipped correction is family-only. If one ever changed `class` it
    would re-mint the composite id of a retained row."""
    from anchorfix_class_reconcile import VERIFIED_FAMILY_CORRECTIONS
    rows = [_row(acc, family=expected, klass="A")
            for acc, (expected, _) in VERIFIED_FAMILY_CORRECTIONS.items()]
    before = [composite_id(r) for r in rows]
    apply_family_corrections(rows, VERIFIED_FAMILY_CORRECTIONS)
    assert [composite_id(r) for r in rows] == before
    assert all(r["class"] == "A" for r in rows)


def test_p46023_is_not_a_correction_target():
    """GRL101 has only the bare class-A family and no curated subfamily, so
    promoting it out of `orphan` would be inference, not curation."""
    from anchorfix_class_reconcile import VERIFIED_FAMILY_CORRECTIONS
    assert "P46023" not in VERIFIED_FAMILY_CORRECTIONS


def test_family_correction_does_not_touch_other_rows():
    rows = [_row("P1", family="orphan"), _row("Q75W84", family="orphan")]
    apply_family_corrections(rows, {"Q75W84": ("orphan", "peptide")})
    assert rows[0]["family"] == "orphan"


# --- atomic write + read-back verification -----------------------------------

def test_write_is_atomic_and_round_trips(tmp_path: Path):
    out = tmp_path / "prod.tsv"
    rows = [_row("P1", taxid="9606"), _row("P2", taxid="7227")]
    write_table_atomic(str(out), rows, _COLS)
    assert not list(tmp_path.glob("*.tmp"))
    back = list(csv.DictReader(out.open(), delimiter="\t"))
    assert [r["accession"] for r in back] == ["P1", "P2"]
    assert list(back[0].keys()) == _COLS


def test_write_leaves_original_intact_when_serialisation_fails(tmp_path: Path):
    """An unexpected column must abort the write, not append itself silently
    or truncate the file it was replacing."""
    out = tmp_path / "prod.tsv"
    out.write_text("sentinel\n")
    rogue = dict(_row("P1"), unexpected_column="x")
    with pytest.raises(ValueError):
        write_table_atomic(str(out), [rogue], _COLS)
    assert out.read_text() == "sentinel\n"
    assert not list(tmp_path.glob("*.tmp"))


def test_verify_written_table_detects_duplicate_keys(tmp_path: Path):
    out = tmp_path / "prod.tsv"
    out.write_text("\t".join(_COLS) + "\nP1\t1\t9606\t\tpeptide\tA\tx\n"
                                      "P1\t1\t9606\t\tpeptide\tA\tx\n")
    with pytest.raises(ValueError, match="duplicate"):
        verify_written_table(str(out), expected_accessions=["P1", "P1"])


def test_verify_written_table_detects_row_count_drift(tmp_path: Path):
    out = tmp_path / "prod.tsv"
    out.write_text("\t".join(_COLS) + "\nP1\t1\t9606\t\tpeptide\tA\tx\n")
    with pytest.raises(ValueError, match="row count"):
        verify_written_table(str(out), expected_accessions=["P1", "P2"])


def test_verify_written_table_detects_reordering(tmp_path: Path):
    """Order is part of the contract: a reordered file would silently change
    every positional join downstream."""
    out = tmp_path / "prod.tsv"
    out.write_text("\t".join(_COLS) + "\nP2\t1\t1\t\tpeptide\tA\tx\n"
                                      "P1\t1\t1\t\tpeptide\tA\tx\n")
    with pytest.raises(ValueError, match="order"):
        verify_written_table(str(out), expected_accessions=["P1", "P2"])


def test_verify_written_table_honours_an_alternate_key_column(tmp_path: Path):
    """The UniProt cache keys on 'queried_accession', not 'accession'."""
    out = tmp_path / "cache.tsv"
    out.write_text("queried_accession\tEntry\nP1\tP1\nP2\tP2\n")
    verify_written_table(str(out), ["P1", "P2"],
                         key_column="queried_accession")
    with pytest.raises(ValueError, match="key column"):
        verify_written_table(str(out), ["P1", "P2"])


def test_verify_written_table_accepts_a_correct_write(tmp_path: Path):
    out = tmp_path / "prod.tsv"
    rows = [_row("P1", taxid="9606"), _row("P2", taxid="7227")]
    write_table_atomic(str(out), rows, _COLS)
    verify_written_table(str(out), expected_accessions=["P1", "P2"])
