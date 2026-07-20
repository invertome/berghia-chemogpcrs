"""Guard: the 'unclassified' sentinel is not a class and must not get a vote.

``classify_gpcr_by_class.py`` emits ``class in {A, B, C, F, unclassified}``.
``build_og_class_majority`` read that column and counted every non-empty value
as a class vote, so the sentinel -- which means "we could not tell" -- was
tallied alongside real evidence. Two consequences:

  1. It INFLATED THE DENOMINATOR. With ``--min-fraction``, an orthogroup of
     6 class-A members and 4 unclassified scored 6/10 = 0.60 rather than
     6/6 = 1.00, so unanimous real evidence could be discarded as a weak
     majority.

  2. It COULD WIN. 3 class-A members against 7 unclassified made
     ``unclassified`` the plurality winner, and the producer then wrote a
     confident per-orthogroup class assignment whose content is "we do not
     know" -- while the only actual evidence said class A.

The fix excludes the sentinel from the vote, and makes an orthogroup with no
real votes say ``unclassified`` EXPLICITLY. Previously such an orthogroup was
simply omitted from the table, so "no class evidence" was indistinguishable
from "this orthogroup was never considered" or "the table is truncated" --
stage 04 and stage 05 both inferred it from a failed lookup.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

import build_og_class_majority as bom  # noqa: E402


# ---------------------------------------------------------------------------
# 1. The sentinel never wins a vote
# ---------------------------------------------------------------------------

def test_sentinel_does_not_win_a_plurality():
    """7 'unclassified' vs 3 real class-A votes: A is the only evidence."""
    assert bom.majority_class(["unclassified"] * 7 + ["A"] * 3) == "A"


def test_sentinel_only_is_not_a_class():
    assert bom.majority_class(["unclassified", "unclassified"]) is None


def test_sentinel_is_case_insensitive():
    assert bom.majority_class(["Unclassified", "UNCLASSIFIED", "C"]) == "C"


def test_real_votes_still_decide_normally():
    assert bom.majority_class(["A", "A", "C"]) == "A"
    assert bom.majority_class(["C", "A"]) == "A"      # tie -> CLASS_ORDER
    assert bom.majority_class([]) is None


# ---------------------------------------------------------------------------
# 2. The sentinel is out of the min-fraction denominator
# ---------------------------------------------------------------------------

def test_sentinel_does_not_dilute_min_fraction():
    """6 real A + 4 sentinel is UNANIMOUS real evidence (6/6), not 0.60."""
    ogs = {"OG0": [f"s{i}" for i in range(10)]}
    classes = {f"s{i}": "A" for i in range(6)}
    classes.update({f"s{i}": "unclassified" for i in range(6, 10)})
    rows = dict(bom.build(ogs, classes, min_fraction=0.9))
    assert rows["OG0"] == "A"


def test_genuinely_weak_majority_is_still_refused():
    """Guard against over-correction: a real 1/4 split must still fail 0.5."""
    ogs = {"OG0": ["s1", "s2", "s3", "s4"]}
    classes = {"s1": "A", "s2": "C", "s3": "F", "s4": "B"}
    rows = dict(bom.build(ogs, classes, min_fraction=0.5))
    assert rows["OG0"] == "unclassified"


# ---------------------------------------------------------------------------
# 3. "No class evidence" is stated, not inferred from absence
# ---------------------------------------------------------------------------

def test_og_with_no_real_votes_says_unclassified_explicitly():
    ogs = {"OG0": ["s1", "s2"], "OG1": ["x1", "x2"]}
    classes = {"s1": "A", "s2": "A", "x1": "unclassified", "x2": "unclassified"}
    rows = dict(bom.build(ogs, classes, min_fraction=0.0))
    assert rows["OG0"] == "A"
    assert rows["OG1"] == "unclassified", "must be stated, not inferred"


def test_og_with_no_classified_members_at_all_is_listed():
    ogs = {"OG0": ["s1"], "OG1": ["absent_from_class_map"]}
    rows = dict(bom.build(ogs, {"s1": "C"}, min_fraction=0.0))
    assert rows["OG1"] == "unclassified"


def test_every_orthogroup_appears_exactly_once():
    """The table is now total over the orthogroups it was given, so a failed
    lookup downstream means a genuinely missing/!truncated table."""
    ogs = {"OG0": ["s1"], "OG1": ["s2"], "OG2": []}
    rows = bom.build(ogs, {"s1": "A"}, min_fraction=0.0)
    assert [og for og, _ in rows] == ["OG0", "OG1", "OG2"]
    assert len(rows) == len(set(og for og, _ in rows))


def test_empty_orthogroup_is_unclassified_not_dropped():
    rows = dict(bom.build({"OG0": []}, {"s1": "A"}, min_fraction=0.0))
    assert rows["OG0"] == "unclassified"


# ---------------------------------------------------------------------------
# 4. The value written is the one stage 04/05 route on
# ---------------------------------------------------------------------------

def test_written_sentinel_matches_the_stage_routing_token(tmp_path: Path):
    """Stage 04/05 set OG_CLASS='unclassified' and read
    results/phylogenies/protein/class_${OG_CLASS}/. The token must match
    exactly, or the explicit row routes somewhere that does not exist."""
    out = tmp_path / "og_class_majority.tsv"
    bom.write_tsv(bom.build({"OG0": ["x"]}, {}, min_fraction=0.0), out)
    lines = out.read_text().splitlines()
    assert lines[0] == "orthogroup\tclass"
    assert lines[1] == "OG0\tunclassified"
    assert b"\r" not in out.read_bytes(), "awk col2 match breaks on trailing CR"


def test_cli_emits_explicit_unclassified_rows(tmp_path: Path):
    ogf = tmp_path / "Orthogroups.tsv"
    ogf.write_text("Orthogroup\tsp1\tsp2\n"
                   "OG0\ts1, s2\ts3\n"
                   "OG1\tu1\t\n")
    cf = tmp_path / "candidate_classes.tsv"
    cf.write_text("seq_id\tclass\ns1\tA\ns2\tA\ns3\tF\nu1\tunclassified\n")
    out = tmp_path / "og_class_majority.tsv"
    assert bom.main(["--orthogroups", str(ogf), "--classes", str(cf),
                     "--out", str(out)]) == 0
    rows = dict(l.split("\t") for l in out.read_text().splitlines()[1:])
    assert rows == {"OG0": "A", "OG1": "unclassified"}


def test_missing_inputs_still_write_nothing(tmp_path: Path):
    """Unchanged contract: a table that cannot be built must NOT be written as
    an all-unclassified table -- stage 04 keeps its own fallback instead."""
    out = tmp_path / "og_class_majority.tsv"
    assert bom.main(["--orthogroups", str(tmp_path / "nope.tsv"),
                     "--classes", str(tmp_path / "alsono.tsv"),
                     "--out", str(out)]) == 0
    assert not out.exists()
