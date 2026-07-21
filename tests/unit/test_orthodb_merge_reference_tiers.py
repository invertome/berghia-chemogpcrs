"""Invariants of the two-tier reference merge.

The merge appends to a file that is gitignored, keyed on identifiers that other
artifacts already reference. Every failure mode here is silent: a re-minted id
still looks like an id, a doubled sequence still passes a row count, and a
colon in a label still writes a valid-looking FASTA header while corrupting
every Newick tree built from it.
"""

import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

import orthodb_merge_reference_tiers as merge  # noqa: E402


ANCHOR_COLUMNS = ["accession", "tier", "taxid", "species", "family", "class",
                  "evidence"]


def _anchor(acc, tier="1", klass="A", family="opsin", taxid="9606"):
    return {"accession": acc, "tier": tier, "taxid": taxid,
            "species": "Homo sapiens", "family": family, "class": klass,
            "evidence": "reviewed"}


def _harvest(gene_id, family="opsin", taxid="6500", verdict="PASS"):
    return {
        "gene_id": gene_id, "family": family, "taxid": taxid,
        "species": "Aplysia californica", "og_id": "OG1", "og_name": "opsins",
        "og_seeds": "3", "og_members": "900", "pf00001_evalue": "1e-40",
        "tm7_coverage": "0.99", "stratum": "invertebrate",
        "quality_verdict": verdict, "verdict": "verified_class_A",
        "sequence": "MEEPLA" * 60,
    }


def _verdict(acc, family="peptide", failed_at=""):
    return {
        "accession": acc, "protein_name": "Myomodulin receptor", "family": family,
        "organism": "Platynereis dumerilii", "taxid": "6359", "phylum": "Annelida",
        "reviewed": "unreviewed", "sequence_length": "360",
        "class_a_evidence": "curated:GPCR-1-family", "evidence_tier":
        "functionally-characterized", "failed_at": failed_at,
        "primary_pmids": "12345678",
    }


# --- identifier construction ------------------------------------------------

def test_composite_id_embeds_class_and_tier():
    """The id is the join key into the sequence store and the npz."""
    assert merge.composite_id(_anchor("P31356", tier="4", klass="A")) == \
        "ANCHOR_A_4_P31356"


def test_orthodb_gene_ids_have_their_colon_replaced():
    """A colon terminates a taxon label in Newick and would corrupt trees."""
    rows, prov, seqs = merge.build_orthodb_rows([_harvest("6500_0:00122c")], "now")
    assert rows[0]["accession"] == "6500_0_00122c"
    assert ":" not in merge.composite_id(rows[0])
    # the untransformed id must survive in provenance, or the row is unlinkable
    assert prov[0]["source_id"] == "6500_0:00122c"


def test_colon_replacement_cannot_collide():
    """Two distinct gene ids must not sanitize onto one accession."""
    rows, _, _ = merge.build_orthodb_rows(
        [_harvest("6500_0:00122c"), _harvest("6500_0:00122d")], "now")
    accs = [r["accession"] for r in rows]
    assert len(set(accs)) == len(accs)


# --- what is admitted -------------------------------------------------------

def test_only_quality_pass_entries_are_admitted():
    rows, _, _ = merge.build_orthodb_rows(
        [_harvest("1_0:a", verdict="PASS"), _harvest("1_0:b", verdict="FAIL")],
        "now")
    assert [r["accession"] for r in rows] == ["1_0_a"]


def test_a_pass_row_whose_class_was_not_verified_is_refused():
    """Class must be verified per entry, never inherited from the orthogroup."""
    bad = _harvest("1_0:a")
    bad["verdict"] = "dropped"
    with pytest.raises(SystemExit, match="not verified per entry"):
        merge.build_orthodb_rows([bad], "now")


def test_failed_literature_entries_are_not_admitted():
    rows, _, _, _ = merge.build_literature_rows(
        [_verdict("A1"), _verdict("A2", failed_at="literature")],
        {"A1": "MKV" * 120, "A2": "MKV" * 120}, "now")
    assert [r["accession"] for r in rows] == ["A1"]


def test_a_passing_entry_with_no_sequence_is_refused():
    with pytest.raises(SystemExit, match="no sequence"):
        merge.build_literature_rows([_verdict("A1")], {}, "now")


def test_held_families_are_withheld_and_reported():
    """A held family must not silently vanish -- it comes back for ratification."""
    rows, _, _, held = merge.build_literature_rows(
        [_verdict("A1", family="peptide"), _verdict("A2", family="unclassified")],
        {"A1": "MKV" * 120, "A2": "MKV" * 120}, "now",
        hold_families=frozenset({"unclassified"}))
    assert [r["accession"] for r in rows] == ["A1"]
    assert [r["accession"] for r in held] == ["A2"]


def test_tiers_are_distinct_so_either_can_be_reverted_alone():
    assert merge.TIER_ORTHODB != merge.TIER_LITERATURE
    assert set(merge.TIER_LABELS) == {merge.TIER_ORTHODB, merge.TIER_LITERATURE}


# --- end-to-end invariants --------------------------------------------------

def _write(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _run(tmp_path, existing, harvest, verdicts, fasta, extra=()):
    anchors = tmp_path / "anchor_set_PROD.tsv"
    _write(anchors, existing, ANCHOR_COLUMNS)
    hpath = tmp_path / "harvest_final.tsv"
    _write(hpath, harvest, list(harvest[0].keys()))
    vpath = tmp_path / "verdicts.tsv"
    _write(vpath, verdicts, list(verdicts[0].keys()))
    fpath = tmp_path / "passing.fasta"
    fpath.write_text("".join(f">{k}\n{v}\n" for k, v in fasta.items()))
    rc = merge.main([
        "--anchors", str(anchors), "--harvest", str(hpath),
        "--literature", str(vpath), "--literature-fasta", str(fpath),
        "--outdir", str(tmp_path),
        "--existing-sequences", str(tmp_path / "absent.tsv"),
        *extra,
    ])
    return rc, anchors


def test_existing_rows_are_preserved_byte_identical_and_in_order(tmp_path):
    existing = [_anchor("P1"), _anchor("P2", tier="3"), _anchor("P3", klass="B",
                                                                family="class-B-secretin")]
    rc, anchors = _run(tmp_path, existing, [_harvest("1_0:a")], [_verdict("A1")],
                       {"A1": "MKV" * 120})
    assert rc == 0
    out = list(csv.DictReader(open(anchors), delimiter="\t"))
    assert out[:3] == existing, "pre-existing rows were altered or reordered"


def test_new_rows_are_appended_never_interleaved(tmp_path):
    existing = [_anchor("P1"), _anchor("P2")]
    rc, anchors = _run(tmp_path, existing, [_harvest("1_0:a")], [_verdict("A1")],
                       {"A1": "MKV" * 120})
    out = list(csv.DictReader(open(anchors), delimiter="\t"))
    assert len(out) == 4
    assert [r["tier"] for r in out[2:]] == [merge.TIER_ORTHODB, merge.TIER_LITERATURE]


def test_an_accession_already_present_is_refused(tmp_path):
    """Re-adding an accession would mint a SECOND id for one protein."""
    existing = [_anchor("A1")]
    with pytest.raises(SystemExit, match="refusing to merge"):
        _run(tmp_path, existing, [_harvest("1_0:a")], [_verdict("A1")],
             {"A1": "MKV" * 120})


def test_duplicate_ids_in_the_existing_table_block_the_merge(tmp_path):
    """Never merge onto a table whose keys are already ambiguous."""
    existing = [_anchor("P1"), _anchor("P1")]
    with pytest.raises(SystemExit, match="duplicate"):
        _run(tmp_path, existing, [_harvest("1_0:a")], [_verdict("A1")],
             {"A1": "MKV" * 120})


def test_provenance_covers_every_new_row_and_dangles_on_none(tmp_path):
    rc, anchors = _run(tmp_path, [_anchor("P1")],
                       [_harvest("1_0:a"), _harvest("1_0:b")],
                       [_verdict("A1")], {"A1": "MKV" * 120})
    assert rc == 0
    out = list(csv.DictReader(open(anchors), delimiter="\t"))
    prov = list(csv.DictReader(
        open(tmp_path / "anchor_set_PROD_tier_provenance.tsv"), delimiter="\t"))
    ids = {merge.composite_id(r) for r in out}
    prov_ids = {r["composite_id"] for r in prov}
    assert prov_ids <= ids, "provenance names an anchor that does not exist"
    assert len(prov) == 3


def test_recorded_length_disagreeing_with_the_written_sequence_fails(tmp_path):
    """The doubling bug: counts pass, only the length reveals it."""
    v = _verdict("A1")
    v["sequence_length"] = "999"          # source says 999
    rc, _ = _run(tmp_path, [_anchor("P1")], [_harvest("1_0:a")], [v],
                 {"A1": "MKV" * 120})     # but 360 residues are written
    assert rc == 1, "a length mismatch must fail the merge, not warn"


def test_snapshot_is_taken_before_writing(tmp_path):
    _run(tmp_path, [_anchor("P1")], [_harvest("1_0:a")], [_verdict("A1")],
         {"A1": "MKV" * 120})
    snaps = list(tmp_path.glob(".merge_snapshot_*"))
    assert snaps, "no snapshot -- references/ is gitignored and unrecoverable"
    assert (snaps[0] / "anchor_set_PROD.tsv").exists()


def test_write_atomic_refuses_an_empty_table(tmp_path):
    with pytest.raises(SystemExit, match="empty table"):
        merge.write_atomic(tmp_path / "x.tsv", [], ANCHOR_COLUMNS)


def test_merged_table_keeps_exactly_the_seven_column_schema(tmp_path):
    """Positional readers of this file must not break."""
    rc, anchors = _run(tmp_path, [_anchor("P1")], [_harvest("1_0:a")],
                       [_verdict("A1")], {"A1": "MKV" * 120})
    header = open(anchors).readline().rstrip("\n").split("\t")
    assert header == ANCHOR_COLUMNS
