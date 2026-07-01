#!/usr/bin/env python3
# build_microswitch_channel.py
# Purpose: PRODUCER for the OR-microswitch BW-numbering evidence channel
#   (Glue task G5, docs/plans/2026-07-01-ml-plm-chemoreceptor-ranking.md).
"""Structural Ballesteros-Weinstein (BW) number-transfer, then score.

scripts/or_microswitch.py is the SCORER: a pure classifier over an
already-BW-numbered ``{bw_position: residue}`` map for one candidate (e.g.
``{"6.48": "Y", "6.50": "P"}`` -> 1). This module is the PRODUCER one level
up the pipeline -- it builds that map in the first place, via STRUCTURAL
NUMBER-TRANSFER: superpose each candidate's AlphaFold model onto a
reference class-A GPCR with known BW numbering (bovine rhodopsin: BW 6.48 =
Trp265, BW 6.50 = Pro267 -- the CWxP toggle motif), then read off the
candidate residues that structurally correspond to those reference
positions.

Pipeline:
    1. parse_alignment           -- one candidate's structural-aligner
                                     output -> (tmscore, correspondence)
    2. transfer_bw_positions     -- project REFERENCE_BW onto the candidate
                                     via `correspondence`, gated on tmscore
    3. or_microswitch.or_microswitch_flag -- score the transferred
                                     {bw: residue} map
    4. build_microswitch_channel -- steps 1-3 for every candidate
    5. write_channel_tsv         -- emit the TSV
                                     rank_aggregation.merge_evidence_channels
                                     left-joins onto the ranking df by id

Aligner-output format targeted (one file per candidate, produced by
scripts/unity/run_microswitch_bw_transfer.sh -- see that script's header
comment for exactly how it derives this from TM-align's superposition):

    TM-score\t<float>
    ref_resnum\tcand_resnum\tcand_residue
    <ref_resnum>\t<cand_resnum>\t<cand_residue>
    ...

Line 1 is the global structural-alignment quality score (TM-align's score
normalized by the REFERENCE structure's length -- TM-align's own stdout
recommends this normalization when comparing structures of unequal
length). Line 2 is a fixed header. Every remaining line is one reference
residue's correspondence: `ref_resnum` is the reference (rhodopsin)
structure's residue number; `cand_resnum`/`cand_residue` are the
structurally-nearest candidate residue's number and one-letter amino acid,
or literally "-"/"-" if that reference position has no aligned candidate
residue (a gap -- e.g. a disordered loop or a genuinely divergent region).
Blank lines and lines starting with "#" are ignored.

Honesty gate (the point of this whole module): a candidate whose structure
is too divergent from the rhodopsin reference to be numbered reliably
(tmscore < min_tmscore), or whose correspondence is gapped at either
required BW position, gets has_or_microswitch_data=False -- never a
fabricated 0/1 call. No aligner binary is invoked here and there are no
import-time side effects; scripts/unity/run_microswitch_bw_transfer.sh is
the only place TM-align actually runs.
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Optional, Tuple

import pandas as pd

from or_microswitch import or_microswitch_flag

# Reference residue-number (bovine rhodopsin) -> Ballesteros-Weinstein
# position. Trp6.48 is the class-A "toggle switch" residue; Pro6.50 anchors
# the CWxP motif (Cys6.47-Trp6.48-x6.49-Pro6.50; see or_microswitch.py).
# Documented module constant, overridable per-call via
# transfer_bw_positions'/build_microswitch_channel's `reference_bw`
# parameter -- e.g. to target a differently-numbered reference structure.
REFERENCE_BW: Dict[int, str] = {265: "6.48", 267: "6.50"}

# The exact columns rank_aggregation.merge_evidence_channels() reads from a
# microswitch TSV: the id_col join key ("id") plus its microswitch
# value_cols (["or_microswitch"]). has_or_microswitch_data is written too,
# for human-readable provenance / parity with build_microswitch_channel()'s
# own per-candidate schema -- mirroring build_structural_channel.py's
# CHANNEL_TSV_COLUMNS -- even though merge_evidence_channels recomputes
# that flag itself from row PRESENCE (id.isin(chan.index)), never from this
# column's value. Because of that row-presence semantics, write_channel_tsv
# below only ever writes rows with has_or_microswitch_data=True: see its
# docstring for why that filter matters here specifically.
CHANNEL_TSV_COLUMNS: List[str] = ["id", "or_microswitch", "has_or_microswitch_data"]

Correspondence = Dict[int, Tuple[int, str]]


def parse_alignment(path: str) -> Tuple[Optional[float], Correspondence]:
    """Parse one candidate's structural-alignment report into (tmscore, correspondence).

    See the module docstring for the exact target format. `correspondence`
    maps reference (rhodopsin) residue-number -> (candidate residue-number,
    candidate one-letter residue), for every reference position that IS
    aligned to a candidate residue. A reference position marked "-"/"-"
    (gapped) is simply absent from the returned dict, same as a position
    never mentioned at all -- callers never need to distinguish the two.

    A missing/nonexistent `path` returns (None, {}) -- "no alignment data
    available" -- rather than raising, mirroring parse_foldseek() /
    load_embeddings()'s existing missing-file contract elsewhere in this
    codebase.
    """
    if not path or not os.path.exists(path):
        return None, {}

    tmscore: Optional[float] = None
    correspondence: Correspondence = {}
    with open(path) as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if line.lower().startswith("tm-score"):
                parts = line.split("\t") if "\t" in line else line.split()
                if len(parts) == 2:
                    try:
                        tmscore = float(parts[1])
                    except ValueError:
                        pass
                continue
            if line.lower().startswith("ref_resnum"):
                continue  # header line
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) != 3:
                continue  # malformed row -- skip rather than raise
            ref_s, cand_s, cand_res = parts
            if cand_s == "-" or cand_res == "-":
                continue  # gapped reference position
            try:
                ref_num = int(ref_s)
                cand_num = int(cand_s)
            except ValueError:
                continue
            if len(cand_res) != 1:
                continue
            correspondence[ref_num] = (cand_num, cand_res.upper())
    return tmscore, correspondence


def transfer_bw_positions(
    tmscore: Optional[float],
    correspondence: Correspondence,
    reference_bw: Dict[int, str] = REFERENCE_BW,
    min_tmscore: float = 0.4,
) -> Optional[Dict[str, str]]:
    """Project `reference_bw` onto a candidate via its `correspondence`.

    Returns ``{bw_position: candidate_residue}`` (e.g. ``{"6.48": "Y",
    "6.50": "P"}``) -- exactly the input shape
    or_microswitch.or_microswitch_flag() expects -- when the structural
    alignment is trustworthy enough to number reliably: `tmscore` is not
    below `min_tmscore`, AND every reference residue number in
    `reference_bw` is present (aligned, not gapped) in `correspondence`.

    Returns None -- "can't number this candidate reliably" -- if `tmscore`
    is missing or below `min_tmscore`, OR if any required reference BW
    position is gapped/absent from `correspondence`. This is all-or-nothing
    by design: or_microswitch_flag() itself needs every required BW
    position to make a definitive call, so a partial map is no more useful
    than no map -- returning None here (rather than a partial dict) is what
    makes the downstream has_or_microswitch_data=False gate mean something.
    """
    if tmscore is None or tmscore < min_tmscore:
        return None
    bw_map: Dict[str, str] = {}
    for ref_resnum, bw_position in reference_bw.items():
        if ref_resnum not in correspondence:
            return None
        _cand_resnum, cand_residue = correspondence[ref_resnum]
        bw_map[bw_position] = cand_residue
    return bw_map


def build_microswitch_channel(
    candidate_alignment_paths: Dict[str, str],
    reference_bw: Dict[int, str] = REFERENCE_BW,
    min_tmscore: float = 0.4,
) -> Dict[str, Dict[str, object]]:
    """parse_alignment -> transfer_bw_positions -> or_microswitch_flag, per candidate.

    Args:
        candidate_alignment_paths: {candidate_id: alignment_report_path}.
        reference_bw: see REFERENCE_BW.
        min_tmscore: see transfer_bw_positions.

    Unlike structural_evidence.structural_channel() / embedding_evidence.
    embedding_channel() (which only emit entries for candidates they have
    usable data for, leaving "no data" to be inferred from omission), this
    function emits an entry for EVERY candidate in
    `candidate_alignment_paths`: has_or_microswitch_data is an explicit
    True/False verdict here, because "we attempted structural
    number-transfer for this candidate and its alignment was too divergent
    to trust" is itself a meaningful result worth recording -- as opposed
    to struct/emb's "this candidate simply has no Foldseek/embedding row",
    which can't distinguish attempted-and-failed from never-attempted.

    Returns:
        {candidate_id: {"or_microswitch": 1|0|None,
                        "has_or_microswitch_data": bool}}
        has_or_microswitch_data is False whenever transfer_bw_positions()
        returns None (untrustworthy alignment or a gapped required BW
        position), and also False in the defensive edge case where
        or_microswitch_flag() itself still returns None despite a non-None
        bw_map.
    """
    channel: Dict[str, Dict[str, object]] = {}
    for candidate_id, path in candidate_alignment_paths.items():
        tmscore, correspondence = parse_alignment(path)
        bw_map = transfer_bw_positions(
            tmscore, correspondence, reference_bw=reference_bw, min_tmscore=min_tmscore
        )
        flag = or_microswitch_flag(bw_map)
        channel[candidate_id] = {
            "or_microswitch": flag,
            "has_or_microswitch_data": bw_map is not None and flag is not None,
        }
    return channel


def write_channel_tsv(channel: Dict[str, Dict[str, object]], path: str) -> None:
    """Write the microswitch channel TSV that merge_evidence_channels() reads.

    Columns are exactly CHANNEL_TSV_COLUMNS. Only rows with
    has_or_microswitch_data=True are written: rank_aggregation.
    merge_evidence_channels derives has_or_microswitch_data from row
    PRESENCE in this file (`id.isin(chan.index)`), not from a column value
    -- its microswitch value_cols is just ``["or_microswitch"]`` -- so an
    included False-flagged row would be silently promoted to True
    downstream. build_microswitch_channel() (unlike structural_channel()/
    embedding_channel()) DOES surface explicit False entries in its
    returned dict, so this filter is what keeps that honest-dormancy
    distinction intact all the way through to the join. Rows are sorted by
    candidate id for deterministic output; an empty/all-unreliable channel
    still writes a valid header-only TSV.
    """
    rows = [
        {"id": candidate_id, **channel[candidate_id]}
        for candidate_id in sorted(channel)
        if channel[candidate_id].get("has_or_microswitch_data")
    ]
    pd.DataFrame(rows, columns=CHANNEL_TSV_COLUMNS).to_csv(path, sep="\t", index=False)


# The COMPOUND suffix run_microswitch_bw_transfer.sh writes each
# per-candidate report under: "<candidate_id>.bw_report.txt". Recovering
# the candidate id needs this whole suffix stripped, not just the trailing
# ".txt" -- os.path.splitext removes only the last extension, leaving a
# stray ".bw_report" on the id, which then never matches a real candidate
# id in the ranking df and silently deadens the entire channel through
# merge_evidence_channels' id-presence join (has_or_microswitch_data=False
# for every candidate).
_REPORT_SUFFIX = ".bw_report.txt"


def _candidate_id_from_name(name: str) -> str:
    """Candidate id from an alignment-report filename.

    Strips run_microswitch_bw_transfer.sh's full "<id>.bw_report.txt"
    compound suffix when present (so 'BersteEVm001t1.bw_report.txt' ->
    'BersteEVm001t1', NOT 'BersteEVm001t1.bw_report'); falls back to plain
    single-extension stripping for any other input filename.
    """
    if name.endswith(_REPORT_SUFFIX):
        return name[: -len(_REPORT_SUFFIX)]
    return os.path.splitext(name)[0]


def resolve_alignment_paths(inputs: List[str]) -> Dict[str, str]:
    """{candidate_id: path} from `--alignments` CLI inputs.

    Each item in `inputs` can be a directory (every regular file directly
    inside it is taken as one candidate's alignment report, non-recursive)
    or a single file path. candidate_id is recovered via
    `_candidate_id_from_name`: the part BEFORE the ".bw_report.txt" suffix
    run_microswitch_bw_transfer.sh writes (e.g.
    'BersteEVm001t1.bw_report.txt' -> 'BersteEVm001t1'), so it matches the
    real candidate ids in the ranking df; any other filename falls back to
    its plain single-extension stem. Items that are neither an existing
    file nor an existing directory are silently skipped.
    """
    paths: Dict[str, str] = {}
    for item in inputs:
        if os.path.isdir(item):
            for name in sorted(os.listdir(item)):
                full = os.path.join(item, name)
                if os.path.isfile(full):
                    paths[_candidate_id_from_name(name)] = full
        elif os.path.isfile(item):
            paths[_candidate_id_from_name(os.path.basename(item))] = item
    return paths


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build the OR-microswitch BW-numbering evidence channel "
                     "TSV for candidate ranking (structural number-transfer "
                     "from a bovine-rhodopsin BW-numbered reference)."
    )
    parser.add_argument(
        "--alignments", nargs="+", required=True,
        help="One or more per-candidate structural-alignment report files "
             "and/or directories of them (scripts/unity/"
             "run_microswitch_bw_transfer.sh output)",
    )
    parser.add_argument("--out", required=True, help="Output channel TSV path")
    parser.add_argument(
        "--min-tmscore", type=float, default=0.4,
        help="Minimum TM-score to trust the structural number-transfer "
             "(default: 0.4)",
    )
    args = parser.parse_args(argv)

    candidate_paths = resolve_alignment_paths(args.alignments)
    channel = build_microswitch_channel(candidate_paths, min_tmscore=args.min_tmscore)
    write_channel_tsv(channel, args.out)
    n_data = sum(1 for r in channel.values() if r["has_or_microswitch_data"])
    n_flagged = sum(1 for r in channel.values() if r["or_microswitch"] == 1)
    print(
        f"[build_microswitch_channel] {len(channel)} candidates -> {args.out} "
        f"({n_data} reliably numbered, {n_flagged} flagged)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
