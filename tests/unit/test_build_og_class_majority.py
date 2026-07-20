"""Tests for scripts/build_og_class_majority.py.

Produces results/classification/og_class_majority.tsv (orthogroup<TAB>class),
consumed by stage 04's per-OG routing (awk col1=og, col2=class, header skipped).

Inputs:
  - OrthoFinder Orthogroups.tsv (OG_id<TAB>sp1 members<TAB>sp2 members ...,
    members comma-space separated within a cell).
  - One or more classifier TSVs (seq_id<TAB>class<TAB>...; class in A/B/C/F).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

import build_og_class_majority as bom


# ---------------------------------------------------------------------------
# parse_class_tsv
# ---------------------------------------------------------------------------

class TestParseClassTsv:
    def test_reads_seq_to_class(self, tmp_path: Path) -> None:
        f = tmp_path / "classes.tsv"
        f.write_text("seq_id\tclass\tevidence_pfam\n"
                     "s1\tA\tPF00001\n"
                     "s2\tC\tPF00003\n")
        m = bom.parse_class_tsv(f)
        assert m == {"s1": "A", "s2": "C"}

    def test_skips_empty_class(self, tmp_path: Path) -> None:
        f = tmp_path / "classes.tsv"
        f.write_text("seq_id\tclass\n"
                     "s1\tA\n"
                     "s2\t\n")          # unclassified → omitted
        m = bom.parse_class_tsv(f)
        assert m == {"s1": "A"}

    def test_multiple_files_merge(self, tmp_path: Path) -> None:
        a = tmp_path / "a.tsv"; a.write_text("seq_id\tclass\ns1\tA\n")
        b = tmp_path / "b.tsv"; b.write_text("seq_id\tclass\nberghia_7\tB\n")
        m = bom.parse_class_tsv(a)
        m.update(bom.parse_class_tsv(b))
        assert m == {"s1": "A", "berghia_7": "B"}


# ---------------------------------------------------------------------------
# parse_orthogroups_tsv
# ---------------------------------------------------------------------------

class TestParseOrthogroupsTsv:
    def test_members_across_species_columns(self, tmp_path: Path) -> None:
        f = tmp_path / "Orthogroups.tsv"
        f.write_text(
            "Orthogroup\tsp1\tsp2\tsp3\n"
            "OG0000000\ts1, s2\ts3\t\n"
            "OG0000001\ts4\t\ts5, s6\n"
        )
        og = bom.parse_orthogroups_tsv(f)
        assert og["OG0000000"] == ["s1", "s2", "s3"]
        assert og["OG0000001"] == ["s4", "s5", "s6"]

    def test_empty_cells_ignored(self, tmp_path: Path) -> None:
        f = tmp_path / "Orthogroups.tsv"
        f.write_text("Orthogroup\tsp1\tsp2\nOG0\t\t\n")
        og = bom.parse_orthogroups_tsv(f)
        assert og["OG0"] == []


# ---------------------------------------------------------------------------
# majority_class
# ---------------------------------------------------------------------------

class TestMajorityClass:
    def test_plurality_wins(self) -> None:
        assert bom.majority_class(["A", "A", "C"]) == "A"

    def test_tie_broken_by_class_order(self) -> None:
        # A,B,C,F order → A beats C on a tie
        assert bom.majority_class(["C", "A"]) == "A"
        assert bom.majority_class(["F", "B"]) == "B"

    def test_empty_returns_none(self) -> None:
        assert bom.majority_class([]) is None


# ---------------------------------------------------------------------------
# build + write
# ---------------------------------------------------------------------------

class TestBuild:
    def test_majority_per_og_and_states_unclassified(self) -> None:
        """An OG with no class evidence is now STATED 'unclassified' rather
        than omitted, so a failed lookup downstream means a missing/truncated
        table instead of doubling as 'no evidence'."""
        ogs = {"OG0": ["s1", "s2", "s3"], "OG1": ["x1", "x2"]}
        classes = {"s1": "A", "s2": "A", "s3": "C"}   # OG1 members unclassified
        rows = bom.build(ogs, classes, min_fraction=0.0)
        assert rows == [("OG0", "A"), ("OG1", "unclassified")]

    def test_min_fraction_reports_weak_majority_as_unclassified(self) -> None:
        ogs = {"OG0": ["s1", "s2", "s3", "s4"]}
        classes = {"s1": "A", "s2": "C", "s3": "F", "s4": "B"}  # 1/4 each, max frac 0.25
        rows = bom.build(ogs, classes, min_fraction=0.5)
        assert rows == [("OG0", "unclassified")]       # no class reaches 50%


class TestWriteTsv:
    def test_output_format_header_and_lf(self, tmp_path: Path) -> None:
        out = tmp_path / "og_class_majority.tsv"
        bom.write_tsv([("OG0", "A"), ("OG1", "C")], out)
        raw = out.read_bytes()
        assert b"\r" not in raw, "must use LF — bash awk col2 match breaks on trailing CR"
        lines = out.read_text().splitlines()
        assert lines[0] == "orthogroup\tclass"
        assert lines[1] == "OG0\tA"
        assert lines[2] == "OG1\tC"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

class TestCLI:
    def test_end_to_end(self, tmp_path: Path) -> None:
        ogf = tmp_path / "Orthogroups.tsv"
        ogf.write_text("Orthogroup\tsp1\tsp2\n"
                       "OG0\ts1, s2\ts3\n"
                       "OG1\tu1\t\n")
        cf = tmp_path / "candidate_classes.tsv"
        cf.write_text("seq_id\tclass\ns1\tA\ns2\tA\ns3\tF\n")  # u1 unclassified
        out = tmp_path / "og_class_majority.tsv"
        rc = bom.main(["--orthogroups", str(ogf), "--classes", str(cf), "--out", str(out)])
        assert rc == 0
        text = out.read_text()
        assert "OG0\tA" in text
        # OG1 has no classified member -> stated explicitly, not omitted.
        assert "OG1\tunclassified" in text

    def test_missing_inputs_writes_nothing(self, tmp_path: Path) -> None:
        """Missing Orthogroups → exit 0 WITHOUT creating the file, so stage 04
        falls back to its class_A default rather than routing all OGs to
        'unclassified' on an empty file."""
        out = tmp_path / "og_class_majority.tsv"
        rc = bom.main([
            "--orthogroups", str(tmp_path / "nope.tsv"),
            "--classes", str(tmp_path / "alsono.tsv"),
            "--out", str(out),
        ])
        assert rc == 0
        assert not out.exists()
