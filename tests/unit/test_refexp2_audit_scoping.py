"""The anchor audit can only reach UniProt-resolvable accessions.

436 of the 1410 class-A anchors are tier-9 OrthoDB gene ids. Querying
``accession:10228_0_00188c`` against UniProt returns HTTP 400 and kills the
whole batch, so the audit must be scoped to the resolvable subset and must
ITEMISE what it skipped rather than reporting a bare count.

The partition guard is tested in BOTH directions: a guard tuned only against
false positives ends up too strict, one tuned only against false negatives too
loose, and both failure modes have happened on this project.
"""
import csv
from pathlib import Path

import pytest

from refexp2_evidence_gate import is_uniprot_accession, partition_by_resolvability

REPO = Path(__file__).resolve().parents[2]
ANCHORS = REPO / "references/anchors/anchor_set_PROD.tsv"


@pytest.mark.parametrize("acc", [
    "P31356",      # [OPQ] form
    "O15973",
    "Q9NYQ7",
    "A0A0D3MML9",  # long [A-NR-Z] form, two trailing groups
    "E7CCA9",
    "V5UZ27",
])
def test_must_accept_real_uniprot_accessions(acc):
    assert is_uniprot_accession(acc), f"{acc} is a real UniProt accession"


@pytest.mark.parametrize("acc", [
    "10228_0_00188c",    # OrthoDB gene id, the case that 400s
    "400682_1_001cae",
    "2023355_0_0028d8",
    "",
    "   ",
])
def test_must_reject_non_uniprot_identifiers(acc):
    assert not is_uniprot_accession(acc)


def test_partition_loses_nothing():
    """Every input lands in exactly one side. A silent drop here would shrink
    the audit denominator without any signal."""
    rows = [
        {"accession": "P31356"},
        {"accession": "10228_0_00188c"},
        {"accession": "A0A0D3MML9"},
        {"accession": "400682_1_001cae"},
    ]
    resolvable, skipped = partition_by_resolvability(rows)
    assert [r["accession"] for r in resolvable] == ["P31356", "A0A0D3MML9"]
    assert [r["accession"] for r in skipped] == ["10228_0_00188c", "400682_1_001cae"]
    assert len(resolvable) + len(skipped) == len(rows)


def test_partition_preserves_the_whole_row():
    """The skipped list is itemised to a file, so it must carry enough to
    identify each anchor, not just its accession."""
    rows = [{"accession": "10228_0_00188c", "tier": "9", "family": "opsin"}]
    _, skipped = partition_by_resolvability(rows)
    assert skipped[0]["tier"] == "9" and skipped[0]["family"] == "opsin"


@pytest.mark.skipif(not ANCHORS.exists(), reason="production anchor set absent")
def test_real_anchor_set_partitions_into_the_measured_split():
    """Real-data guard. The counts are the ones that made the audit 400."""
    with open(ANCHORS, newline="") as fh:
        class_a = [r for r in csv.DictReader(fh, delimiter="\t") if r["class"] == "A"]
    resolvable, skipped = partition_by_resolvability(class_a)
    assert len(class_a) == 1410
    assert len(resolvable) == 974
    assert len(skipped) == 436
    # every skipped anchor is tier 9; if a non-harvest tier ever lands here it
    # means a new unresolvable namespace entered the set unnoticed.
    assert {r["tier"] for r in skipped} == {"9"}
