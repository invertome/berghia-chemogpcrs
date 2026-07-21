"""The class-A anchors are validated by TWO paths, and neither may be defaulted.

UniProt-accession anchors are gated by the refexp2 audit. The tier-9 OrthoDB
harvest anchors have no UniProt record at all, so they are validated instead by
``harvest_final.tsv`` (organism check, 7tm_1 Pfam evidence, TM coverage, gate
verdict). An anchor covered by neither is an EXPLICIT unknown, never an
optimistic pass.

The join between the two namespaces differs by a single separator, and the
unnormalised join silently returns zero rows, so the normalisation is asserted
here on real data rather than assumed.
"""
import csv
import re
from pathlib import Path

import pytest

from orthodb_snapshot_anchors import (
    build_snapshot,
    gene_id_to_accession,
    load_harvest_verdicts,
)

REPO = Path(__file__).resolve().parents[2]
ANCHORS = REPO / "references/anchors/anchor_set_PROD.tsv"
HARVEST = REPO / "results/ranking/diagnostics/orthodb/harvest_final.tsv"
UNIPROT_ACC = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")

_HARVEST_HEADER = ["gene_id", "verdict", "quality_verdict", "organism_verified"]


def _harvest(tmp_path: Path, rows) -> Path:
    p = tmp_path / "harvest_final.tsv"
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=_HARVEST_HEADER, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return p


def _ok(gene_id):
    return {"gene_id": gene_id, "verdict": "verified_class_A",
            "quality_verdict": "PASS", "organism_verified": "True"}


def test_gene_id_separator_is_normalised_to_the_anchor_form():
    """The last ':' becomes '_'; nothing else about the id changes."""
    assert gene_id_to_accession("10228_0:00188c") == "10228_0_00188c"
    assert gene_id_to_accession("400682_1:001cae") == "400682_1_001cae"


def test_uniprot_accession_is_left_untouched():
    assert gene_id_to_accession("P31356") == "P31356"


def test_harvest_verdict_requires_all_three_checks(tmp_path: Path):
    """A single failing check disqualifies; they are ANDed, not sampled."""
    p = _harvest(tmp_path, [
        _ok("1_0:aaa"),
        {**_ok("1_0:bbb"), "verdict": "dropped"},
        {**_ok("1_0:ccc"), "quality_verdict": "FAIL"},
        {**_ok("1_0:ddd"), "organism_verified": "False"},
    ])
    v = load_harvest_verdicts(p)
    # Each assertion names its own guard. These are SIBLING guards over the same
    # row, so a bare `assert x is False` makes every one of them fail with an
    # identical diff and you cannot tell which check was removed.
    assert v["1_0_aaa"] is True, "a row passing all three checks was rejected"
    assert v["1_0_bbb"] is False, "guard 'verdict == verified_class_A' did not fire"
    assert v["1_0_ccc"] is False, "guard 'quality_verdict == PASS' did not fire"
    assert v["1_0_ddd"] is False, "guard 'organism_verified == True' did not fire"


def test_anchor_in_neither_source_is_explicitly_unvalidated(tmp_path: Path):
    """The load-bearing case: no source must NOT mean 'passes'.

    The previous behaviour was audit.get(acc, True), which collapsed 'never
    checked' into 'verified good'.
    """
    anchors = [{"class": "A", "accession": "Q00001", "family": "opsin",
                "taxid": "9606", "species": "Homo sapiens", "tier": "10",
                "evidence": "literature-gated"}]
    snap = build_snapshot(anchors, audit={}, harvest={})
    assert snap[0]["validation_source"] == "none"
    assert snap[0]["pass_class_a"] == "False"
    assert snap[0]["use_primary"] == "False"


def test_each_source_validates_its_own_anchors(tmp_path: Path):
    anchors = [
        {"class": "A", "accession": "P31356", "family": "aminergic", "taxid": "9606",
         "species": "Homo sapiens", "tier": "1", "evidence": "reviewed"},
        {"class": "A", "accession": "10228_0_00188c", "family": "opsin", "taxid": "10228",
         "species": "x", "tier": "9", "evidence": "orthodb-harvest"},
    ]
    snap = build_snapshot(anchors, audit={"P31356": True},
                          harvest={"10228_0_00188c": True})
    by = {r["accession"]: r for r in snap}
    assert by["P31356"]["validation_source"] == "uniprot_audit"
    assert by["10228_0_00188c"]["validation_source"] == "orthodb_harvest"
    assert all(r["use_primary"] == "True" for r in snap)


def test_a_failing_audit_still_beats_a_passing_harvest(tmp_path: Path):
    """Precedence is audit-first: a curated FALSE is not overridden."""
    anchors = [{"class": "A", "accession": "P31356", "family": "opsin", "taxid": "9606",
                "species": "x", "tier": "1", "evidence": "reviewed"}]
    snap = build_snapshot(anchors, audit={"P31356": False}, harvest={"P31356": True})
    assert snap[0]["validation_source"] == "uniprot_audit"
    assert snap[0]["pass_class_a"] == "False"


@pytest.mark.skipif(not (ANCHORS.exists() and HARVEST.exists()),
                    reason="production anchor/harvest artifacts absent")
def test_every_orthodb_anchor_resolves_against_the_real_harvest():
    """Real-data key-overlap guard.

    A fixture cannot catch a wrong key: it encodes the key it assumes. The
    unnormalised join returns ZERO rows here, which would read as 'no
    validation exists' rather than as a broken join.
    """
    with open(ANCHORS, newline="") as fh:
        anchors = [r for r in csv.DictReader(fh, delimiter="\t") if r["class"] == "A"]
    odb = [r["accession"] for r in anchors if not UNIPROT_ACC.match(r["accession"])]
    assert odb, "expected tier-9 OrthoDB anchors in the production set"

    verdicts = load_harvest_verdicts(HARVEST)
    resolved = [a for a in odb if a in verdicts]
    assert len(resolved) == len(odb), (
        f"only {len(resolved)}/{len(odb)} OrthoDB anchors resolved against the "
        "harvest validation table")
    assert all(verdicts[a] for a in resolved), "an admitted anchor failed its own gate"
