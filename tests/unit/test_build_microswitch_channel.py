"""Tests for scripts/build_microswitch_channel.py.

Glue task G5 -- the OR-microswitch BW-numbering channel PRODUCER.
scripts/or_microswitch.py is the pure SCORER (already tested in
test_or_microswitch.py): given an already-BW-numbered {bw_position:
residue} map for one candidate, it returns 1/0/None. This module is the
PRODUCER one level up -- it builds that map via STRUCTURAL NUMBER-TRANSFER,
superposing each candidate's AlphaFold model onto a reference class-A GPCR
with known Ballesteros-Weinstein numbering (bovine rhodopsin: BW 6.48 =
Trp265, BW 6.50 = Pro267) and reading off the candidate residues that
structurally correspond to those two reference positions.

Honesty gate under test throughout: a candidate whose structural alignment
to the reference is too poor to trust (TM-score below threshold), or whose
correspondence is gapped at a required BW position, must come out
has_or_microswitch_data=False -- never a fabricated 0/1 call.

No real structures or aligner binaries are used here -- every test parses a
small hand-written fixture of the aligner-output report format documented
at the top of build_microswitch_channel.py (produced in production by
scripts/unity/run_microswitch_bw_transfer.sh from TM-align).

Coverage:
    - parse_alignment: (tmscore, correspondence) from a fixture report;
      gapped positions excluded; missing file -> (None, {})
    - transfer_bw_positions: correct residue at each BW position;
      below/at the min_tmscore boundary; gapped required position -> None;
      reference_bw is genuinely overridable
    - build_microswitch_channel: end-to-end parse -> transfer -> score,
      for both the flag=1 and flag=0 cases, and the has-no-data case
    - write_channel_tsv: output TSV columns EXACTLY match
      rank_aggregation.merge_evidence_channels' expected microswitch
      columns, unreliable rows are never written, and it genuinely
      interops with merge_evidence_channels (not just column-name matching)
    - resolve_alignment_paths / main(): the --alignments dir-or-files CLI
      plumbing
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

import rank_aggregation as ra
from build_microswitch_channel import (
    REFERENCE_BW,
    build_microswitch_channel,
    main,
    parse_alignment,
    resolve_alignment_paths,
    transfer_bw_positions,
    write_channel_tsv,
)


def _write_alignment(tmp_path: Path, name: str, tmscore, rows: list) -> Path:
    """Write a small BW-transfer alignment-report fixture.

    rows: list of (ref_resnum, cand_resnum, cand_residue) triples. Use
    ("-", "-") for cand_resnum/cand_residue to mark a gapped reference
    position (no structurally-corresponding candidate residue).
    """
    lines = [f"TM-score\t{tmscore}", "ref_resnum\tcand_resnum\tcand_residue"]
    for ref_resnum, cand_resnum, cand_residue in rows:
        lines.append(f"{ref_resnum}\t{cand_resnum}\t{cand_residue}")
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return path


# A small window around the two REFERENCE_BW positions (264-267: the CWxP
# toggle motif). Candidate has the documented altered-toggle signature (Y
# aligned to reference 265, P aligned to reference 267).
GOOD_ALIGNMENT_ROWS = [
    (263, 261, "C"),
    (264, 262, "W"),
    (265, 263, "Y"),
    (266, 264, "L"),
    (267, 265, "P"),
]

# Canonical (unaltered) toggle: Trp aligned to reference 265.
CANONICAL_ALIGNMENT_ROWS = [
    (265, 263, "W"),
    (266, 264, "L"),
    (267, 265, "P"),
]


# --------------------------------------------------------------------------- #
# parse_alignment
# --------------------------------------------------------------------------- #

def test_parse_alignment_reads_tmscore_and_correspondence(tmp_path):
    path = _write_alignment(tmp_path, "good.txt", 0.85, GOOD_ALIGNMENT_ROWS)
    tmscore, correspondence = parse_alignment(str(path))
    assert tmscore == pytest.approx(0.85)
    assert correspondence[264] == (262, "W")
    assert correspondence[265] == (263, "Y")
    assert correspondence[267] == (265, "P")


def test_parse_alignment_excludes_gapped_positions_from_correspondence(tmp_path):
    rows = [
        (264, 262, "W"),
        (265, "-", "-"),  # gapped -- no structurally-corresponding residue
        (266, 264, "L"),
        (267, 265, "P"),
    ]
    path = _write_alignment(tmp_path, "gapped.txt", 0.85, rows)
    tmscore, correspondence = parse_alignment(str(path))
    assert 265 not in correspondence
    assert correspondence[267] == (265, "P")


def test_parse_alignment_missing_file_returns_none_and_empty(tmp_path):
    missing = tmp_path / "does_not_exist.txt"
    tmscore, correspondence = parse_alignment(str(missing))
    assert tmscore is None
    assert correspondence == {}


def test_parse_alignment_skips_blank_and_comment_lines(tmp_path):
    lines = [
        "# BW-transfer alignment report for cand1",
        "TM-score\t0.85",
        "",
        "ref_resnum\tcand_resnum\tcand_residue",
        "265\t263\tY",
        "267\t265\tP",
    ]
    path = tmp_path / "commented.txt"
    path.write_text("\n".join(lines) + "\n")
    tmscore, correspondence = parse_alignment(str(path))
    assert tmscore == pytest.approx(0.85)
    assert correspondence[265] == (263, "Y")


# --------------------------------------------------------------------------- #
# transfer_bw_positions
# --------------------------------------------------------------------------- #

def test_transfer_bw_positions_extracts_correct_residues(tmp_path):
    path = _write_alignment(tmp_path, "good.txt", 0.85, GOOD_ALIGNMENT_ROWS)
    tmscore, correspondence = parse_alignment(str(path))
    assert transfer_bw_positions(tmscore, correspondence) == {"6.48": "Y", "6.50": "P"}


def test_transfer_bw_positions_below_min_tmscore_is_none(tmp_path):
    path = _write_alignment(tmp_path, "low_tm.txt", 0.2, GOOD_ALIGNMENT_ROWS)
    tmscore, correspondence = parse_alignment(str(path))
    assert transfer_bw_positions(tmscore, correspondence, min_tmscore=0.4) is None


def test_transfer_bw_positions_at_min_tmscore_boundary_passes(tmp_path):
    """tmscore == min_tmscore is NOT below threshold -- the boundary passes
    (mirrors the "< threshold" wording, strict less-than)."""
    path = _write_alignment(tmp_path, "boundary.txt", 0.4, GOOD_ALIGNMENT_ROWS)
    tmscore, correspondence = parse_alignment(str(path))
    assert transfer_bw_positions(tmscore, correspondence, min_tmscore=0.4) is not None


def test_transfer_bw_positions_none_tmscore_is_none(tmp_path):
    missing = tmp_path / "does_not_exist.txt"
    tmscore, correspondence = parse_alignment(str(missing))
    assert transfer_bw_positions(tmscore, correspondence) is None


def test_transfer_bw_positions_gapped_required_position_is_none(tmp_path):
    rows = [
        (264, 262, "W"),
        (265, "-", "-"),  # 6.48 gapped
        (266, 264, "L"),
        (267, 265, "P"),
    ]
    path = _write_alignment(tmp_path, "gapped_648.txt", 0.85, rows)
    tmscore, correspondence = parse_alignment(str(path))
    assert transfer_bw_positions(tmscore, correspondence) is None


def test_transfer_bw_positions_honors_custom_reference_bw(tmp_path):
    """REFERENCE_BW is documented as overridable -- a caller can target an
    entirely different reference residue-number -> BW mapping."""
    rows = [(100, 50, "Y"), (102, 52, "P")]
    path = _write_alignment(tmp_path, "custom.txt", 0.9, rows)
    tmscore, correspondence = parse_alignment(str(path))
    custom_bw = {100: "6.48", 102: "6.50"}
    bw_map = transfer_bw_positions(tmscore, correspondence, reference_bw=custom_bw)
    assert bw_map == {"6.48": "Y", "6.50": "P"}


def test_reference_bw_default_is_rhodopsin_265_267():
    assert REFERENCE_BW == {265: "6.48", 267: "6.50"}


# --------------------------------------------------------------------------- #
# build_microswitch_channel
# --------------------------------------------------------------------------- #

def test_build_channel_good_alignment_flags_altered_toggle(tmp_path):
    path = _write_alignment(tmp_path, "cand1.txt", 0.85, GOOD_ALIGNMENT_ROWS)
    channel = build_microswitch_channel({"cand1": str(path)})
    assert channel["cand1"] == {"or_microswitch": 1, "has_or_microswitch_data": True}


def test_build_channel_canonical_toggle_flags_zero_not_none(tmp_path):
    path = _write_alignment(tmp_path, "cand2.txt", 0.85, CANONICAL_ALIGNMENT_ROWS)
    channel = build_microswitch_channel({"cand2": str(path)})
    assert channel["cand2"] == {"or_microswitch": 0, "has_or_microswitch_data": True}


def test_build_channel_low_tmscore_has_no_data(tmp_path):
    path = _write_alignment(tmp_path, "cand3.txt", 0.2, GOOD_ALIGNMENT_ROWS)
    channel = build_microswitch_channel({"cand3": str(path)}, min_tmscore=0.4)
    assert channel["cand3"] == {"or_microswitch": None, "has_or_microswitch_data": False}


def test_build_channel_missing_alignment_file_has_no_data(tmp_path):
    missing = tmp_path / "does_not_exist.txt"
    channel = build_microswitch_channel({"cand4": str(missing)})
    assert channel["cand4"] == {"or_microswitch": None, "has_or_microswitch_data": False}


def test_build_channel_processes_multiple_candidates_independently(tmp_path):
    good = _write_alignment(tmp_path, "good.txt", 0.85, GOOD_ALIGNMENT_ROWS)
    bad = _write_alignment(tmp_path, "bad.txt", 0.1, GOOD_ALIGNMENT_ROWS)
    channel = build_microswitch_channel({"good": str(good), "bad": str(bad)})
    assert channel["good"]["has_or_microswitch_data"] is True
    assert channel["good"]["or_microswitch"] == 1
    assert channel["bad"]["has_or_microswitch_data"] is False
    assert channel["bad"]["or_microswitch"] is None


# --------------------------------------------------------------------------- #
# write_channel_tsv -- must match merge_evidence_channels' microswitch schema
# --------------------------------------------------------------------------- #

_EXPECTED_TSV_COLUMNS = ["id", "or_microswitch", "has_or_microswitch_data"]


def test_write_channel_tsv_columns_exactly_match_merge_evidence_channels(tmp_path):
    channel = {"cand1": {"or_microswitch": 1, "has_or_microswitch_data": True}}
    out = tmp_path / "microswitch.tsv"
    write_channel_tsv(channel, str(out))
    df = pd.read_csv(out, sep="\t")
    assert list(df.columns) == _EXPECTED_TSV_COLUMNS


def test_write_channel_tsv_empty_channel_writes_header_only(tmp_path):
    out = tmp_path / "microswitch.tsv"
    write_channel_tsv({}, str(out))
    df = pd.read_csv(out, sep="\t")
    assert list(df.columns) == _EXPECTED_TSV_COLUMNS
    assert len(df) == 0


def test_write_channel_tsv_omits_rows_without_data(tmp_path):
    """A candidate with has_or_microswitch_data=False must be OMITTED from
    the TSV, not written with a False value: merge_evidence_channels
    derives has_or_microswitch_data from row PRESENCE
    (id.isin(chan.index)), not from a column value (its microswitch
    value_cols is just ["or_microswitch"]) -- an included False-flagged row
    would be silently promoted to True downstream."""
    channel = {
        "reliable": {"or_microswitch": 0, "has_or_microswitch_data": True},
        "unreliable": {"or_microswitch": None, "has_or_microswitch_data": False},
    }
    out = tmp_path / "microswitch.tsv"
    write_channel_tsv(channel, str(out))
    df = pd.read_csv(out, sep="\t")
    assert list(df["id"]) == ["reliable"]


def test_write_channel_tsv_rows_sorted_by_id(tmp_path):
    channel = {
        "z_cand": {"or_microswitch": 0, "has_or_microswitch_data": True},
        "a_cand": {"or_microswitch": 1, "has_or_microswitch_data": True},
    }
    out = tmp_path / "microswitch.tsv"
    write_channel_tsv(channel, str(out))
    df = pd.read_csv(out, sep="\t")
    assert list(df["id"]) == ["a_cand", "z_cand"]


def test_write_channel_tsv_interops_with_merge_evidence_channels(tmp_path):
    """Real proof of interop, not just column-name matching: feed the
    written TSV straight into the real rank_aggregation.
    merge_evidence_channels and confirm the left-join + has_data flag
    behave as documented there -- True only for the reliable candidate,
    False for both the unreliable (present-but-False) and unseen ids."""
    channel = {
        "reliable": {"or_microswitch": 1, "has_or_microswitch_data": True},
        "unreliable": {"or_microswitch": None, "has_or_microswitch_data": False},
    }
    out = tmp_path / "microswitch.tsv"
    write_channel_tsv(channel, str(out))

    df = pd.DataFrame({"id": ["reliable", "unreliable", "unseen"]})
    merged = ra.merge_evidence_channels(df, microswitch_tsv=str(out))
    result = merged.set_index("id")

    assert list(merged["has_or_microswitch_data"]) == [True, False, False]
    assert result.loc["reliable", "or_microswitch"] == 1
    assert pd.isna(result.loc["unreliable", "or_microswitch"])


# --------------------------------------------------------------------------- #
# resolve_alignment_paths
# --------------------------------------------------------------------------- #

def test_resolve_alignment_paths_from_directory(tmp_path):
    align_dir = tmp_path / "alignments"
    align_dir.mkdir()
    _write_alignment(align_dir, "cand1.txt", 0.85, GOOD_ALIGNMENT_ROWS)
    _write_alignment(align_dir, "cand2.txt", 0.85, GOOD_ALIGNMENT_ROWS)
    resolved = resolve_alignment_paths([str(align_dir)])
    assert set(resolved) == {"cand1", "cand2"}
    assert resolved["cand1"] == str(align_dir / "cand1.txt")


def test_resolve_alignment_paths_from_explicit_files(tmp_path):
    p1 = _write_alignment(tmp_path, "cand1.txt", 0.85, GOOD_ALIGNMENT_ROWS)
    p2 = _write_alignment(tmp_path, "cand2.txt", 0.85, GOOD_ALIGNMENT_ROWS)
    resolved = resolve_alignment_paths([str(p1), str(p2)])
    assert resolved == {"cand1": str(p1), "cand2": str(p2)}


def test_resolve_alignment_paths_skips_nonexistent_items(tmp_path):
    p1 = _write_alignment(tmp_path, "cand1.txt", 0.85, GOOD_ALIGNMENT_ROWS)
    resolved = resolve_alignment_paths([str(p1), str(tmp_path / "nope.txt")])
    assert set(resolved) == {"cand1"}


def test_resolve_alignment_paths_strips_bw_report_compound_suffix_from_directory(tmp_path):
    """Real production convention: run_microswitch_bw_transfer.sh writes each
    report as '<candidate_id>.bw_report.txt'. The recovered id MUST be the
    bare candidate id ('BersteEVm001t1'), NOT 'BersteEVm001t1.bw_report' --
    os.path.splitext strips only the LAST extension, leaving the stray
    '.bw_report' fragment, which then never matches a real candidate id in
    the ranking df, so merge_evidence_channels joins ZERO rows and the whole
    channel silently dies in production (the exact bug this guards)."""
    align_dir = tmp_path / "alignments"
    align_dir.mkdir()
    _write_alignment(
        align_dir, "BersteEVm001t1.bw_report.txt", 0.85, GOOD_ALIGNMENT_ROWS
    )
    resolved = resolve_alignment_paths([str(align_dir)])
    assert set(resolved) == {"BersteEVm001t1"}
    assert resolved["BersteEVm001t1"] == str(align_dir / "BersteEVm001t1.bw_report.txt")


def test_resolve_alignment_paths_strips_bw_report_compound_suffix_from_explicit_file(tmp_path):
    p = _write_alignment(
        tmp_path, "BersteEVm001t1.bw_report.txt", 0.85, GOOD_ALIGNMENT_ROWS
    )
    resolved = resolve_alignment_paths([str(p)])
    assert set(resolved) == {"BersteEVm001t1"}
    assert resolved["BersteEVm001t1"] == str(p)


# --------------------------------------------------------------------------- #
# main() CLI
# --------------------------------------------------------------------------- #

def test_main_cli_writes_expected_tsv(tmp_path):
    align_dir = tmp_path / "alignments"
    align_dir.mkdir()
    _write_alignment(align_dir, "cand1.txt", 0.85, GOOD_ALIGNMENT_ROWS)
    _write_alignment(align_dir, "cand2.txt", 0.1, GOOD_ALIGNMENT_ROWS)
    out = tmp_path / "out.tsv"

    rc = main(["--alignments", str(align_dir), "--out", str(out)])

    assert rc == 0
    df = pd.read_csv(out, sep="\t")
    assert list(df.columns) == _EXPECTED_TSV_COLUMNS
    assert list(df["id"]) == ["cand1"]


def test_main_cli_respects_min_tmscore_override(tmp_path):
    align_dir = tmp_path / "alignments"
    align_dir.mkdir()
    _write_alignment(align_dir, "cand1.txt", 0.5, GOOD_ALIGNMENT_ROWS)
    out = tmp_path / "out.tsv"

    rc = main(
        ["--alignments", str(align_dir), "--out", str(out), "--min-tmscore", "0.9"]
    )

    assert rc == 0
    df = pd.read_csv(out, sep="\t")
    assert df.empty
