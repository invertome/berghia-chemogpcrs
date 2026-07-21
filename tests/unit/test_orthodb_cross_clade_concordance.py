"""Unit tests for the cross-clade family concordance classifier.

This test is the centrepiece of the OrthoDB validation: it decides whether an
orthogroup correctly connects a characterized vertebrate receptor to a
characterized invertebrate receptor of the same family. The three buckets carry
very different meanings (right call / confident wrong call / no call), so the
classifier must never collapse one into another.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "orthodb_cross_clade_concordance",
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "orthodb_cross_clade_concordance.py",
)
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


# ------------------------------------------------------------- modal family --

def test_modal_family_picks_the_majority():
    fam, n, unique = mod.modal_family(["opsin", "opsin", "peptide"])
    assert (fam, n, unique) == ("opsin", 2, True)


def test_modal_family_reports_ties_as_non_unique():
    """A tie means the vertebrate anchors disagree among themselves, which is
    a different finding from a clean modal family and must not be hidden."""
    _, _, unique = mod.modal_family(["opsin", "peptide"])
    assert unique is False


def test_modal_family_rejects_empty():
    with pytest.raises(ValueError):
        mod.modal_family([])


# ------------------------------------------------------------ the 3 buckets --

def test_agree_when_modal_vertebrate_family_matches():
    bucket, mode = mod.classify_anchor("opsin", ["opsin", "opsin", "peptide"])
    assert bucket == "AGREE"
    assert mode == "opsin"


def test_disagree_when_modal_vertebrate_family_differs():
    """The dangerous outcome: the route would make a confident wrong call."""
    bucket, mode = mod.classify_anchor("opsin", ["peptide", "peptide"])
    assert bucket == "DISAGREE"
    assert mode == "peptide"


def test_not_spanned_when_no_vertebrate_anchor_shares_the_orthogroup():
    """No answer is not the same as a wrong answer and gets its own bucket."""
    bucket, mode = mod.classify_anchor("opsin", [])
    assert bucket == "NOT_SPANNED"
    assert mode == ""


def test_tied_vertebrate_families_are_not_silently_called_agree():
    bucket, _ = mod.classify_anchor("opsin", ["opsin", "peptide"])
    assert bucket == "AMBIGUOUS_VERTEBRATE_MODE"


def test_single_vertebrate_anchor_is_enough_to_span():
    bucket, mode = mod.classify_anchor("aminergic", ["aminergic"])
    assert bucket == "AGREE"
    assert mode == "aminergic"


def test_minority_match_still_counts_as_disagree():
    """Matching a non-modal family is a disagreement: the harvest would take
    the modal call, not the one that happens to be right."""
    bucket, mode = mod.classify_anchor("opsin", ["peptide", "peptide", "opsin"])
    assert bucket == "DISAGREE"
    assert mode == "peptide"


# ------------------------------------------------------- end-to-end on a level --

def _fixture(tmp_path):
    """Two orthogroups at Metazoa level:
    OG_A: molluscan opsin + 2 vertebrate opsins           -> AGREE
    OG_B: molluscan peptide + 2 vertebrate aminergics     -> DISAGREE
    OG_C: molluscan lipid, alone                          -> NOT_SPANNED
    """
    rows = [
        ("M1", "opsin", "Mollusca", "True"),
        ("M2", "peptide", "Mollusca", "True"),
        ("M3", "lipid", "Mollusca", "True"),
        ("V1", "opsin", "Vertebrata", "True"),
        ("V2", "opsin", "Vertebrata", "True"),
        ("V3", "aminergic", "Vertebrata", "True"),
        ("V4", "aminergic", "Vertebrata", "True"),
    ]
    audit = tmp_path / "mapping_audit.tsv"
    audit.write_text(
        "accession\tfamily\tclade\tuse_primary\n"
        + "".join(f"{a}\t{f}\t{c}\t{p}\n" for a, f, c, p in rows)
    )
    (tmp_path / "anchor_gene_map.tsv").write_text(
        "".join(f"{a}\tg_{a}\t{a}\n" for a, _, _, _ in rows)
    )
    og = {"M1": "OG_A", "V1": "OG_A", "V2": "OG_A",
          "M2": "OG_B", "V3": "OG_B", "V4": "OG_B",
          "M3": "OG_C"}
    (tmp_path / "anchor_og2genes.tsv").write_text(
        "".join(f"{o}\tg_{a}\n" for a, o in og.items())
    )
    (tmp_path / "anchor_og_meta.tsv").write_text(
        "OG_A\t33208\tgroupA\nOG_B\t33208\tgroupB\nOG_C\t33208\tgroupC\n"
    )
    (tmp_path / "odb_levels.tsv").write_text("33208\tMetazoa\t1\t1\t1\n")
    return tmp_path


def test_end_to_end_buckets_on_a_synthetic_level(tmp_path):
    data = mod.load(_fixture(tmp_path))
    rows, summary, confusion = mod.analyse_level("33208", data)
    buckets = {r["accession"]: r["bucket"] for r in rows}
    assert buckets == {"M1": "AGREE", "M2": "DISAGREE", "M3": "NOT_SPANNED"}
    assert summary["invertebrate_anchors_tested"] == 3
    assert confusion[("peptide", "aminergic")] == 1


def test_vertebrate_anchors_are_not_themselves_tested(tmp_path):
    """The test is about invertebrate anchors; vertebrates are the reference."""
    data = mod.load(_fixture(tmp_path))
    rows, _, _ = mod.analyse_level("33208", data)
    assert all(r["clade"] != "Vertebrata" for r in rows)
    assert len(rows) == 3


def test_per_family_breakdown_separates_the_outcomes(tmp_path):
    data = mod.load(_fixture(tmp_path))
    _, summary, _ = mod.analyse_level("33208", data)
    assert summary["by_family"]["opsin"] == {"AGREE": 1}
    assert summary["by_family"]["peptide"] == {"DISAGREE": 1}
    assert summary["by_family"]["lipid"] == {"NOT_SPANNED": 1}


def test_level_with_no_invertebrate_anchors_yields_nothing(tmp_path):
    """At Vertebrata level the test is vacuous by construction and must report
    zero tested rather than a misleading 100% agreement."""
    d = _fixture(tmp_path)
    (d / "anchor_og_meta.tsv").write_text(
        "OG_A\t7742\tgroupA\nOG_B\t7742\tgroupB\nOG_C\t7742\tgroupC\n"
    )
    (d / "mapping_audit.tsv").write_text(
        "accession\tfamily\tclade\tuse_primary\n"
        "V1\topsin\tVertebrata\tTrue\nV2\topsin\tVertebrata\tTrue\n"
    )
    (d / "anchor_gene_map.tsv").write_text("V1\tg_V1\tV1\nV2\tg_V2\tV2\n")
    (d / "anchor_og2genes.tsv").write_text("OG_A\tg_V1\nOG_A\tg_V2\n")
    data = mod.load(d)
    rows, summary, _ = mod.analyse_level("7742", data)
    assert rows == []
    assert summary["invertebrate_anchors_tested"] == 0


def test_non_primary_anchors_are_excluded(tmp_path):
    d = _fixture(tmp_path)
    txt = (d / "mapping_audit.tsv").read_text().replace(
        "M1\topsin\tMollusca\tTrue", "M1\topsin\tMollusca\tFalse"
    )
    (d / "mapping_audit.tsv").write_text(txt)
    data = mod.load(d)
    rows, _, _ = mod.analyse_level("33208", data)
    assert "M1" not in {r["accession"] for r in rows}


def test_load_aborts_when_inputs_are_missing(tmp_path):
    with pytest.raises(SystemExit):
        mod.load(tmp_path)
