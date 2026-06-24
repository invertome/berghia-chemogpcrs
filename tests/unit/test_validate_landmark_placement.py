"""Unit tests for validate_landmark_placement.py — Task 0: landmark loading
and family-vocab normalizer."""
import json
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))
import validate_landmark_placement as vlp  # noqa: E402

ANCHOR_TSV = "accession\ttier\ttaxid\tspecies\tfamily\tclass\tevidence\n" \
    "P1\t1\t6637\tsp\topsin\tA\treviewed\n" \
    "Q2\t2\t6340\tPlatynereis\taminergic\tA\texperimental\n" \
    "Q3\t3\t7227\tDrosophila\tpeptide\tA\treviewed\n" \
    "Z9\t2\t6340\tPlatynereis\tclass-F-frizzled\tF\texperimental\n"
ANCHOR_FA = ">ANCHOR_A_2_Q2\nMKVL\n>ANCHOR_A_3_Q3\nMMMM\n>ANCHOR_F_2_Z9\nWWWW\n>ANCHOR_A_1_P1\nAAAA\n"

def test_load_landmarks_filters_tier23_and_class(tmp_path):
    t = tmp_path / "a.tsv"; t.write_text(ANCHOR_TSV)
    f = tmp_path / "a.fa"; f.write_text(ANCHOR_FA)
    got = vlp.load_landmarks(str(t), str(f), "A")
    ids = sorted(r["id"] for r in got)
    assert ids == ["ANCHOR_A_2_Q2", "ANCHOR_A_3_Q3"]      # tier1 P1 and class-F Z9 excluded
    fam = {r["id"]: r["family"] for r in got}
    assert fam["ANCHOR_A_2_Q2"] == "aminergic"
    seqs = {r["id"]: r["seq"] for r in got}
    assert seqs["ANCHOR_A_2_Q2"] == "MKVL"
    assert seqs["ANCHOR_A_3_Q3"] == "MMMM"

def test_load_landmarks_warns_on_missing_seq(tmp_path):
    t = tmp_path / "a.tsv"; t.write_text(ANCHOR_TSV)
    f = tmp_path / "a.fa"; f.write_text(">ANCHOR_A_2_Q2\nMKVL\n")  # ANCHOR_A_3_Q3 absent
    with pytest.warns(UserWarning):
        vlp.load_landmarks(str(t), str(f), "A")

def test_normalize_family_identity_and_unknown():
    assert vlp.normalize_family("class-F-frizzled") == "class-f-frizzled"
    assert vlp.normalize_family(" Opsin ") == "opsin"
    assert vlp.normalize_family("") == "?unknown?"
    assert vlp.normalize_family(None) == "?unknown?"


# ---- Task 1: placement runner -------------------------------------------

JPLACE = json.dumps({
    "version": 3,
    "tree": "((Berg1:0.1{0},Berg2:0.1{1}):0.1{2},(RefM1:0.2{3},LandX:0.2{4}):0.1{5}):0{6};",
    "fields": ["edge_num", "like_weight_ratio"],
    "placements": [{"n": ["ANCHOR_A_2_Q2"], "p": [[3, 0.95], [4, 0.05]]}],
})


def test_place_landmarks_parses_edge_and_leaves(tmp_path, monkeypatch):
    jp = tmp_path / "p.jplace"
    jp.write_text(JPLACE)
    monkeypatch.setattr(vlp, "run_epang_placement", lambda *a, **k: str(jp))
    placements, edge_to_leaves = vlp.place_landmarks(
        ref_msa="ignored", ref_tree="ignored",
        landmarks=[{"id": "ANCHOR_A_2_Q2", "family": "aminergic", "seq": "MK"}],
        outdir=str(tmp_path))
    assert placements["ANCHOR_A_2_Q2"]["edge"] == 3          # best LWR (0.95 > 0.05)
    assert abs(placements["ANCHOR_A_2_Q2"]["lwr"] - 0.95) < 1e-9
    assert edge_to_leaves[3] == ["RefM1"]                    # edge 3 -> its leaf
