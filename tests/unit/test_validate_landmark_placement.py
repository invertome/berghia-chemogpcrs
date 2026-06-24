"""Unit tests for validate_landmark_placement.py — Task 0: landmark loading
and family-vocab normalizer."""
import json
import os
import sys

import pytest
from ete3 import Tree

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


# ---- Task 2: per-candidate + clade-consensus functional reads ---------------


def test_nearest_landmark_assigns_family_and_gates_lwr():
    # clean reference tree: NO landmark leaf — the landmark is *placed* on an edge.
    t = Tree("((Berg1:1,Berg2:1)100:1,(RefM1:1,RefM2:1)90:1);", format=0)
    # landmark LM1 placed on edge 7 == the branch above the (Berg1,Berg2) clade,
    # so its placement node is MRCA(Berg1,Berg2). edge_to_leaves comes from place_landmarks().
    edge_to_leaves = {7: ["Berg1", "Berg2"]}
    fam = {"LM1": "aminergic"}
    low = vlp.nearest_landmark_per_candidate(
        t, {"LM1": {"edge": 7, "lwr": 0.5}}, fam, ["Berg1", "Berg2"],
        lwr_min=0.80, edge_to_leaves=edge_to_leaves)
    assert low == []                                   # gated out by lwr
    rows = vlp.nearest_landmark_per_candidate(
        t, {"LM1": {"edge": 7, "lwr": 0.95}}, fam, ["Berg1", "Berg2"],
        lwr_min=0.80, edge_to_leaves=edge_to_leaves)
    assert {r["candidate"] for r in rows} == {"Berg1", "Berg2"}
    assert all(r["family"] == "aminergic" for r in rows)


def test_clade_consensus_majority_over_supported_berghia_clade():
    # (Berg1,Berg2,Berg3) form a supported all-Berghia clade; RefM1 outside.
    t = Tree("(((Berg1:1,Berg2:1):1,Berg3:1)95:1,RefM1:1);", format=0)
    rows = [
        {"candidate": "Berg1", "family": "aminergic", "landmark": "L1", "distance": 1, "lwr": 0.9},
        {"candidate": "Berg2", "family": "aminergic", "landmark": "L1", "distance": 1, "lwr": 0.9},
        {"candidate": "Berg3", "family": "peptide",   "landmark": "L2", "distance": 2, "lwr": 0.9},
    ]
    out = vlp.clade_consensus(t, rows, ["Berg1", "Berg2", "Berg3"], supp_min=80)
    # the supported all-Berghia clade -> majority aminergic (2 vs 1)
    top = [c for c in out if set(c["clade_leaves"]) == {"Berg1", "Berg2", "Berg3"}]
    assert top and top[0]["consensus_family"] == "aminergic" and top[0]["n_called"] == 3


def test_nearest_landmark_tiebreak_is_alphabetical():
    t = Tree("((Berg1:1,Berg2:1)100:1,(RefM1:1,RefM2:1)90:1);", format=0)
    edge_to_leaves = {7: ["Berg1", "Berg2"]}          # both landmarks on the same edge
    placements = {"ZZZ": {"edge": 7, "lwr": 0.95}, "AAA": {"edge": 7, "lwr": 0.95}}
    fam = {"ZZZ": "peptide", "AAA": "aminergic"}
    rows = vlp.nearest_landmark_per_candidate(t, placements, fam, ["Berg1"], edge_to_leaves=edge_to_leaves)
    assert len(rows) == 1
    assert rows[0]["landmark"] == "AAA"               # alphabetically-first wins the tie
    assert rows[0]["family"] == "aminergic"
