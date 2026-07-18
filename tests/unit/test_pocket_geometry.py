"""Unit tests for scripts/pocket_geometry.py.

Receptor-intrinsic small-molecule-vs-peptide discriminant (Foster et al.
2019, Cell 177:1933): peptide-binding class-A GPCRs have a LONG ECL2
(> ~20 residues) and a large orthosteric cavity; small-molecule / odorant
receptors have a short ECL2 and a small enclosed pocket.

These tests drive the two pure, cleanly-TDD-able pieces hard:
  * ecl2_length -- parse the 7 TM helices out of a TMbed per-residue
    topology string (contiguous runs of helix labels) and return the length
    of the loop between the 4th and 5th helix (ECL2).
  * small_molecule_vs_peptide_score -- provisional score (higher = more
    small-molecule-like) from ECL2 length + pocket volume.
  * parse_fpocket_info -- pull the top-ranked pocket's volume out of an
    fpocket <name>_info.txt (pure text parse, no binary needed).

Structure parsing (pocket_volume) is exercised on a tiny synthetic PDB /
mmCIF written in the test: we assert it parses and returns a number, NOT an
exact volume (the proxy is provisional and flagged for calibration).
"""
from __future__ import annotations

import math
from pathlib import Path

import pytest

# conftest.py adds scripts/ to sys.path
import pocket_geometry as pg


# --------------------------------------------------------------------------
# helpers to build synthetic TMbed topology strings
# --------------------------------------------------------------------------
def build_gpcr_topology(
    nterm=3, tm=(20, 20, 20, 20, 20, 20, 20),
    loops=(4, 8, 6, 15, 5, 9), cterm=10,
):
    """Build a canonical 7-TM GPCR TMbed topology string.

    Layout (N-term extracellular, alternating helix orientation h/H):
        o*nterm  h*tm1  i*ICL1  H*tm2  o*ECL1  h*tm3  i*ICL2  H*tm4
        o*ECL2   h*tm5  i*ICL3  H*tm6  o*ECL3  h*tm7  i*cterm
    loops = (ICL1, ECL1, ICL2, ECL2, ICL3, ECL3).
    """
    icl1, ecl1, icl2, ecl2, icl3, ecl3 = loops
    helix = "hHhHhHh"  # alternating orientation, both count as helix
    parts = ["o" * nterm]
    loop_seq = ["i" * icl1, "o" * ecl1, "i" * icl2, "o" * ecl2,
                "i" * icl3, "o" * ecl3, "i" * cterm]
    for k in range(7):
        parts.append(helix[k] * tm[k])
        parts.append(loop_seq[k])
    return "".join(parts)


# --------------------------------------------------------------------------
# ecl2_length
# --------------------------------------------------------------------------
def test_ecl2_length_canonical_seven_tm():
    topo = build_gpcr_topology(loops=(4, 8, 6, 15, 5, 9))
    assert pg.ecl2_length(topo) == 15


def test_ecl2_length_picks_the_right_loop_not_ecl1_or_ecl3():
    # ECL1=8, ECL2=13, ECL3=9 all distinct -> must return the middle one.
    topo = build_gpcr_topology(loops=(4, 8, 6, 13, 5, 9))
    assert pg.ecl2_length(topo) == 13


def test_ecl2_length_counts_both_upper_and_lower_case_helix_as_helix():
    # All helices upper-case H -> still 7 segments, ECL2 still found.
    topo = build_gpcr_topology(loops=(4, 8, 6, 17, 5, 9)).replace("h", "H")
    assert pg.ecl2_length(topo) == 17


def test_ecl2_length_short_peptide_vs_small_molecule_regimes():
    short = build_gpcr_topology(loops=(4, 8, 6, 7, 5, 9))
    long = build_gpcr_topology(loops=(4, 8, 6, 28, 5, 9))
    assert pg.ecl2_length(short) == 7
    assert pg.ecl2_length(long) == 28


def test_ecl2_length_minimum_single_residue_loop():
    # Consecutive helix runs are always separated by >=1 non-helix residue.
    topo = build_gpcr_topology(loops=(4, 8, 6, 1, 5, 9))
    assert pg.ecl2_length(topo) == 1


def test_ecl2_length_ignores_extra_tm_helices_beyond_seven():
    # 8 helices: ECL2 is still the loop between the 4th and 5th, unchanged.
    topo = build_gpcr_topology(loops=(4, 8, 6, 11, 5, 9)) + "o" * 3 + "H" * 20
    assert pg.ecl2_length(topo) == 11


def test_ecl2_length_tolerates_trailing_whitespace_and_newline():
    topo = build_gpcr_topology(loops=(4, 8, 6, 14, 5, 9)) + "\n"
    assert pg.ecl2_length(topo) == 14


def test_ecl2_length_returns_none_when_fewer_than_five_tm():
    # Only 4 helix segments -> ECL2 (between TM4 and TM5) is undefined.
    topo = "ooo" + "h" * 20 + "iiii" + "H" * 20 + "oooo" + "h" * 20 + "iii" + "H" * 20 + "oo"
    assert pg.ecl2_length(topo) is None


def test_ecl2_length_returns_none_on_empty_string():
    assert pg.ecl2_length("") is None


def test_ecl2_length_returns_none_on_none():
    assert pg.ecl2_length(None) is None


def test_ecl2_length_returns_none_when_no_helix_labels_at_all():
    assert pg.ecl2_length("oooiiiooo") is None


# --------------------------------------------------------------------------
# parse_tm_segments (the parser under ecl2_length)
# --------------------------------------------------------------------------
def test_parse_tm_segments_finds_seven_contiguous_runs():
    topo = build_gpcr_topology()
    segs = pg.parse_tm_segments(topo)
    assert len(segs) == 7
    # each segment is a (start, end) inclusive pair, monotonically ordered
    for start, end in segs:
        assert start <= end
    assert all(segs[i][1] < segs[i + 1][0] for i in range(len(segs) - 1))


def test_parse_tm_segments_empty_for_no_helices():
    assert pg.parse_tm_segments("oooiii") == []


# --------------------------------------------------------------------------
# small_molecule_vs_peptide_score
# --------------------------------------------------------------------------
def test_score_short_ecl2_small_pocket_is_small_molecule_like():
    s = pg.small_molecule_vs_peptide_score(5, 400.0)
    assert s > 0.5
    assert 0.0 <= s <= 1.0


def test_score_long_ecl2_large_pocket_is_peptide_like():
    s = pg.small_molecule_vs_peptide_score(30, 1500.0)
    assert s < 0.5
    assert 0.0 <= s <= 1.0


def test_score_is_monotonic_decreasing_in_ecl2_length():
    hi = pg.small_molecule_vs_peptide_score(5, 800.0)
    lo = pg.small_molecule_vs_peptide_score(25, 800.0)
    assert hi >= lo


def test_score_is_monotonic_decreasing_in_pocket_volume():
    hi = pg.small_molecule_vs_peptide_score(12, 400.0)
    lo = pg.small_molecule_vs_peptide_score(12, 1400.0)
    assert hi >= lo


def test_score_uses_ecl2_alone_when_pocket_missing():
    s = pg.small_molecule_vs_peptide_score(5, None)
    assert s is not None and 0.0 <= s <= 1.0
    assert s > 0.5


def test_score_uses_pocket_alone_when_ecl2_missing():
    s = pg.small_molecule_vs_peptide_score(None, 400.0)
    assert s is not None and 0.0 <= s <= 1.0
    assert s > 0.5


def test_score_none_when_both_missing():
    assert pg.small_molecule_vs_peptide_score(None, None) is None


# --------------------------------------------------------------------------
# parse_fpocket_info (pure text parser, no fpocket binary required)
# --------------------------------------------------------------------------
FPOCKET_INFO = """Pocket 1 :
\tScore : \t0.834
\tDruggability Score : \t0.912
\tNumber of Alpha Spheres : \t142
\tTotal SASA : \t321.5
\tVolume : \t1053.42

Pocket 2 :
\tScore : \t0.201
\tDruggability Score : \t0.050
\tNumber of Alpha Spheres : \t44
\tVolume : \t402.10
"""


def test_parse_fpocket_info_returns_top_pocket_volume():
    vol = pg.parse_fpocket_info(FPOCKET_INFO)
    assert vol == pytest.approx(1053.42)


def test_parse_fpocket_info_none_when_no_pocket():
    assert pg.parse_fpocket_info("no pockets found\n") is None


# --------------------------------------------------------------------------
# pocket_volume (structure parsing) on a tiny synthetic structure
# --------------------------------------------------------------------------
def _write_tiny_pdb(path: Path, n=14):
    """A small non-coplanar CA-only 'structure' (a helical point cloud)."""
    lines = []
    for i in range(n):
        x = 6.0 * math.cos(i * 0.9)
        y = 6.0 * math.sin(i * 0.9)
        z = 1.5 * i
        lines.append(
            "ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C"
            % (i + 1, i + 1, x, y, z)
        )
    lines.append("TER")
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def test_pocket_volume_proxy_on_pdb_returns_positive_number(tmp_path):
    pdb = tmp_path / "cand.pdb"
    _write_tiny_pdb(pdb)
    vol, method = pg.pocket_volume(str(pdb))
    assert isinstance(vol, float)
    assert vol > 0.0
    # fpocket is not on PATH in the test env -> proxy method must be used
    assert method.startswith("proxy")


def test_pocket_volume_handles_mmcif(tmp_path):
    from Bio.PDB import PDBParser, MMCIFIO

    pdb = tmp_path / "cand.pdb"
    _write_tiny_pdb(pdb)
    struct = PDBParser(QUIET=True).get_structure("cand", str(pdb))
    cif = tmp_path / "cand.cif"
    io = MMCIFIO()
    io.set_structure(struct)
    io.save(str(cif))

    vol, method = pg.pocket_volume(str(cif))
    assert isinstance(vol, float)
    assert vol > 0.0


def test_pocket_volume_missing_file_returns_none(tmp_path):
    vol, method = pg.pocket_volume(str(tmp_path / "does_not_exist.pdb"))
    assert vol is None
    assert method == "none"


# --------------------------------------------------------------------------
# parse_topology_file
# --------------------------------------------------------------------------
def test_parse_topology_file_reads_three_line_records(tmp_path):
    topo1 = build_gpcr_topology(loops=(4, 8, 6, 12, 5, 9))
    seq1 = "A" * len(topo1)
    topo2 = build_gpcr_topology(loops=(4, 8, 6, 22, 5, 9))
    seq2 = "A" * len(topo2)
    f = tmp_path / "topologies.3line"
    f.write_text(
        ">cand_A\n" + seq1 + "\n" + topo1 + "\n"
        + ">cand_B some description\n" + seq2 + "\n" + topo2 + "\n"
    )
    parsed = pg.parse_topology_file(str(f))
    assert set(parsed) == {"cand_A", "cand_B"}
    assert pg.ecl2_length(parsed["cand_A"]) == 12
    assert pg.ecl2_length(parsed["cand_B"]) == 22


# --------------------------------------------------------------------------
# main() CLI end-to-end (proxy path, fpocket absent)
# --------------------------------------------------------------------------
def test_main_emits_tsv_with_expected_columns(tmp_path):
    import csv

    topo = build_gpcr_topology(loops=(4, 8, 6, 15, 5, 9))
    seq = "A" * len(topo)
    topo_file = tmp_path / "cand.3line"
    topo_file.write_text(">cand_1\n" + seq + "\n" + topo + "\n")

    pdb = tmp_path / "cand_1.pdb"
    _write_tiny_pdb(pdb)

    out = tmp_path / "pocket_geometry.tsv"
    rc = pg.main([
        "--topology", str(topo_file),
        "--structures", str(pdb),
        "--output", str(out),
    ])
    assert rc == 0
    assert out.exists()

    with open(out) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert set(row) == {
        "candidate_id", "ecl2_length", "pocket_volume",
        "sm_vs_peptide_score", "method_used",
    }
    assert row["candidate_id"] == "cand_1"
    assert row["ecl2_length"] == "15"
    assert float(row["pocket_volume"]) > 0.0
    assert row["method_used"].startswith("proxy")
    # score is present and in [0,1]
    assert 0.0 <= float(row["sm_vs_peptide_score"]) <= 1.0


def test_main_topology_only_when_no_structure_found(tmp_path):
    import csv

    topo = build_gpcr_topology(loops=(4, 8, 6, 9, 5, 9))
    seq = "A" * len(topo)
    topo_file = tmp_path / "cand.3line"
    topo_file.write_text(">solo\n" + seq + "\n" + topo + "\n")

    empty_struct_dir = tmp_path / "structs"
    empty_struct_dir.mkdir()

    out = tmp_path / "out.tsv"
    rc = pg.main([
        "--topology", str(topo_file),
        "--structures", str(empty_struct_dir),
        "--output", str(out),
    ])
    assert rc == 0
    with open(out) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert row["candidate_id"] == "solo"
    assert row["ecl2_length"] == "9"
    assert row["pocket_volume"] == ""      # no structure
    assert row["method_used"] == "none"
    # score still computable from ECL2 alone
    assert row["sm_vs_peptide_score"] != ""
