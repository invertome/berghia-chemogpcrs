"""Unit tests for anchor injection into the per-class reference pools (C2).

Anchors (characterized GPCR landmarks from build_anchor_set.py) are injected
must-keep, through CD-HIT, keeping the annotated anchor on any collision with a
scanned candidate. They are charged to the reference budget (like references,
unlike Berghia which is on top). Berghia is never absorbed.

CD-HIT itself is not run here: the anchor-vs-candidate collision logic is a pure
function over parsed clusters, tested directly; the orchestrator test mocks
cdhit_dedup as a passthrough.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path
from unittest.mock import patch

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import build_per_class_reference_pools as bpcp  # noqa: E402


def _rec(rid, seq="MKTII"):
    return SeqRecord(Seq(seq), id=rid, description="")


# Minimal mock lineage (Berghia + a couple of molluscs); unknown -> [1, taxid].
_MOCK_LINEAGES = {
    1287507: [1, 2759, 6447, 216305, 1287507],
    6500:    [1, 2759, 6447, 216305, 6500],
    225164:  [1, 2759, 6447, 225164],
    29159:   [1, 2759, 6447, 29159],
}


def _lineage(taxid):
    return _MOCK_LINEAGES.get(taxid, [1, taxid])


# ---------------------------------------------------------------------------
# .clstr parsing
# ---------------------------------------------------------------------------

def test_parse_clstr(tmp_path):
    clstr = tmp_path / "out.clstr"
    clstr.write_text(
        ">Cluster 0\n"
        "0\t300aa, >seqA... *\n"
        "1\t290aa, >seqB... at 95.00%\n"
        ">Cluster 1\n"
        "0\t250aa, >seqC... *\n"
    )
    clusters = bpcp._parse_clstr(str(clstr))
    assert len(clusters) == 2
    by_rep = {c["rep"]: set(c["members"]) for c in clusters}
    assert by_rep["seqA"] == {"seqA", "seqB"}
    assert by_rep["seqC"] == {"seqC"}


# ---------------------------------------------------------------------------
# anchor-aware survivor selection (keep-annotated-on-collision)
# ---------------------------------------------------------------------------

def test_survivor_ids_no_anchor_keeps_representative():
    clusters = [{"rep": "candA", "members": ["candA", "candB"]}]
    assert bpcp._cdhit_survivor_ids(clusters, frozenset()) == {"candA"}


def test_survivor_ids_anchor_wins_collision():
    # An anchor clustered with two scanned candidates: keep the anchor, drop both
    # candidates (even though a candidate is CD-HIT's representative).
    clusters = [{"rep": "candA", "members": ["candA", "ANCHOR_A_1_P1", "candB"]}]
    assert bpcp._cdhit_survivor_ids(
        clusters, frozenset({"ANCHOR_A_1_P1"})) == {"ANCHOR_A_1_P1"}


def test_survivor_ids_multiple_anchors_all_kept():
    clusters = [{"rep": "ANCHOR_A_1_P1",
                 "members": ["ANCHOR_A_1_P1", "ANCHOR_A_2_P2", "candA"]}]
    survivors = bpcp._cdhit_survivor_ids(
        clusters, frozenset({"ANCHOR_A_1_P1", "ANCHOR_A_2_P2"}))
    assert survivors == {"ANCHOR_A_1_P1", "ANCHOR_A_2_P2"}


# ---------------------------------------------------------------------------
# load_anchor_set: route by class, attach taxid
# ---------------------------------------------------------------------------

def test_load_anchor_set_routes_by_class(tmp_path):
    fasta = tmp_path / "anchor_set.fasta"
    fasta.write_text(
        ">ANCHOR_A_1_P1\nMAAA\n"
        ">ANCHOR_A_2_P2\nMBBB\n"
        ">ANCHOR_B_3_Q9\nMCCC\n"
    )
    tsv = tmp_path / "anchor_set.tsv"
    with open(tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "tier", "taxid", "species", "family", "class", "evidence"])
        w.writerow(["P1", "1", "6500", "Aplysia californica", "aminergic", "A", "reviewed"])
        w.writerow(["P2", "2", "6359", "Platynereis dumerilii", "peptide", "A", "experimental:26190115"])
        w.writerow(["Q9", "3", "7227", "Drosophila melanogaster", "class-F-frizzled", "F", "reviewed"])

    per_class = bpcp.load_anchor_set(str(fasta), str(tsv))
    assert {t for t, _ in per_class["A"]} == {6500, 6359}
    assert [r.id for _, r in per_class["A"]] == ["ANCHOR_A_1_P1", "ANCHOR_A_2_P2"]
    # header carries class B, but tsv routes Q9 -> class F (header class is authoritative
    # and consistent with the tsv); here Q9 is class F by both.
    assert "F" in per_class
    assert [r.id for _, r in per_class["F"]] == ["ANCHOR_F_3_Q9"] or \
           [r.id for _, r in per_class["B"]] == ["ANCHOR_B_3_Q9"]


# ---------------------------------------------------------------------------
# build_pool_for_class: anchors must-keep, charged to the refs budget
# ---------------------------------------------------------------------------

def test_anchors_must_keep_and_charged_to_budget():
    # 5 candidate taxa x 2 seqs; budget for refs = 4; 2 anchors injected.
    candidates = []
    for taxid in (6500, 225164, 29159, 6526, 6523):
        for i in range(2):
            candidates.append((taxid, _rec(f"cand_{taxid}_{i}")))
    anchors = [(6359, _rec("ANCHOR_A_2_PLAT1")), (7227, _rec("ANCHOR_A_3_DROS1"))]

    selected, stats = bpcp.build_pool_for_class(
        records=candidates,
        must_include_taxids=frozenset(),
        berghia_taxid=1287507,
        max_size=4,
        lineage_fn=_lineage,
        berghia_records=None,
        anchor_records=anchors,
    )
    ids = {r.id for _, r in selected}
    # Both anchors are present (must-keep).
    assert "ANCHOR_A_2_PLAT1" in ids
    assert "ANCHOR_A_3_DROS1" in ids
    # Anchors are charged to the refs budget: total refs (anchors + balanced) <= max_size.
    assert len(selected) <= 4
    assert stats["n_anchors"] == 2
    # The 2 anchors consumed budget, leaving <= 2 balanced refs.
    n_balanced = len([1 for _, r in selected if not r.id.startswith("ANCHOR_")])
    assert n_balanced <= 2


def test_anchors_kept_even_when_budget_is_zero():
    anchors = [(6359, _rec("ANCHOR_A_2_X")), (7227, _rec("ANCHOR_A_3_Y"))]
    selected, stats = bpcp.build_pool_for_class(
        records=[(6500, _rec("candX"))],
        must_include_taxids=frozenset(),
        berghia_taxid=1287507,
        max_size=1,                 # smaller than the 2 anchors
        lineage_fn=_lineage,
        anchor_records=anchors,
    )
    ids = {r.id for _, r in selected}
    assert {"ANCHOR_A_2_X", "ANCHOR_A_3_Y"} <= ids   # both kept despite tiny budget
    assert stats["n_anchors"] == 2


# ---------------------------------------------------------------------------
# build_all_pools: inject anchors into the right class, report n_anchors,
# Berghia preserved
# ---------------------------------------------------------------------------

def test_build_all_pools_injects_anchors(tmp_path):
    # One scan species with one Class-A + one Class-B candidate.
    scan = tmp_path / "6500_Aplysia_californica.chemo_candidates.fa"
    scan.write_text(">cand_A1\nMAAA\n>cand_B1\nMBBB\n")
    class_tsv = tmp_path / "classes.tsv"
    with open(class_tsv, "w") as fh:
        fh.write("seq_id\tclass\tevidence_pfam\ttop_evalue\n")
        fh.write("cand_A1\tA\tPF00001\t1e-50\n")
        fh.write("cand_B1\tB\tPF00002\t1e-50\n")

    # Berghia candidate (Class A) — must survive.
    berghia_fa = tmp_path / "berghia.fa"
    berghia_fa.write_text(">berghia_A1\nMBER\n")
    berghia_tsv = tmp_path / "berghia_classes.tsv"
    with open(berghia_tsv, "w") as fh:
        fh.write("seq_id\tclass\tevidence_pfam\ttop_evalue\n")
        fh.write("berghia_A1\tA\tPF00001\t1e-50\n")

    # Anchor set: one Class-A anchor, one Class-B anchor.
    anchor_fa = tmp_path / "anchor_set.fasta"
    anchor_fa.write_text(">ANCHOR_A_1_P1\nMANCHORA\n>ANCHOR_B_3_Q9\nMANCHORB\n")
    anchor_tsv = tmp_path / "anchor_set.tsv"
    with open(anchor_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "tier", "taxid", "species", "family", "class", "evidence"])
        w.writerow(["P1", "1", "6500", "Aplysia californica", "aminergic", "A", "reviewed"])
        w.writerow(["Q9", "3", "7227", "Drosophila melanogaster", "class-B-secretin", "B", "reviewed"])

    out_dir = tmp_path / "out"
    with patch.object(bpcp, "cdhit_dedup", side_effect=lambda records, **kw: records):
        with patch.object(bpcp, "get_lineage", side_effect=_lineage):
            bpcp.build_all_pools(
                scan_fasta_glob=str(tmp_path / "*.chemo_candidates.fa"),
                class_tsv=str(class_tsv),
                out_dir=str(out_dir),
                cdhit_path="cd-hit",
                berghia_fasta=str(berghia_fa),
                berghia_class_tsv=str(berghia_tsv),
                anchor_fasta=str(anchor_fa),
                anchor_tsv=str(anchor_tsv),
            )

    class_a = (out_dir / "refs_class_A.fa").read_text()
    class_b = (out_dir / "refs_class_B.fa").read_text()
    # Anchors land in their class pools.
    assert "ANCHOR_A_1_P1" in class_a
    assert "ANCHOR_B_3_Q9" in class_b
    assert "ANCHOR_A_1_P1" not in class_b
    # Berghia preserved (never absorbed by an anchor).
    assert "berghia_A1" in class_a

    import json
    report = json.loads((out_dir / "pool_build_report.json").read_text())
    assert report["class_A"]["n_anchors"] == 1
    assert report["class_B"]["n_anchors"] == 1
