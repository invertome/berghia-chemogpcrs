#!/usr/bin/env python3
"""_scan_proteome_filter_helper.py

Python filter/emit core for scan_proteome_for_chemoreceptors.sh.

Called by the bash wrapper after HMM scan + TMbed have run:
    python3 _scan_proteome_filter_helper.py \
        --tmbed-prediction <prediction> \
        --hmm-positive-ids <ids_file> \
        --input-fa <gpcr_positive.aa> \
        --out-fa <candidates.fa> \
        --out-tsv <scan_record.tsv> \
        [--min-tm 6] [--min-confidence 0.5] [--force]

Exit codes:
    0  normal exit (outputs written, or skipped)
    1  error
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from Bio import SeqIO

# ---------------------------------------------------------------------------
# Shared constants (keep in sync with build_braker4_samples_csv.py)
# ---------------------------------------------------------------------------

_NON_ID_CHARS = re.compile(r"[^A-Za-z0-9_]+")
_MULTI_UNDERSCORE = re.compile(r"_+")


# ---------------------------------------------------------------------------
# Public types
# ---------------------------------------------------------------------------

@dataclass
class PredictionRow:
    """One row from the TMbed prediction file (columns 3 + 5 carry the signal)."""
    tm_count: int
    tm_confidence: float


# ---------------------------------------------------------------------------
# 1. sanitize_sample_name
# ---------------------------------------------------------------------------

def sanitize_sample_name(name: str) -> str:
    """Normalize to [A-Za-z0-9_] — collapse non-id char runs to '_',
    dedupe underscore runs, strip leading/trailing underscores.

    Matches the identical function in build_braker4_samples_csv.py and
    consolidate_proteomes_for_genome_wide_og.py — do NOT diverge.
    """
    s = _NON_ID_CHARS.sub("_", name.strip())
    s = _MULTI_UNDERSCORE.sub("_", s)
    return s.strip("_")


# ---------------------------------------------------------------------------
# 2. parse_tmbed_prediction
# ---------------------------------------------------------------------------

def parse_tmbed_prediction(prediction_path: Path) -> dict[str, PredictionRow]:
    """Parse TMbed prediction file; return dict[seq_id → PredictionRow].

    File format (tab-separated, no header):
        col 1 = seq_id
        col 3 = confidence  (0-1 float)
        col 5 = TM-region count  (int)

    This exactly mirrors the stage-02 awk filter:
        awk -F"\\t" 'NF >= 5 && $5+0 >= min_tm && $3+0 >= min_conf {print $1}'

    Rows with fewer than 5 fields are silently skipped (same as the awk NF≥5 guard).
    """
    result: dict[str, PredictionRow] = {}
    prediction_path = Path(prediction_path)
    if not prediction_path.exists() or prediction_path.stat().st_size == 0:
        return result

    with prediction_path.open() as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            seq_id = parts[0]
            try:
                confidence = float(parts[2])
                tm_count = int(float(parts[4]))  # int(float()) tolerates "6.0"
            except (ValueError, IndexError):
                continue
            result[seq_id] = PredictionRow(tm_count=tm_count, tm_confidence=confidence)
    return result


# ---------------------------------------------------------------------------
# 3. apply_filter
# ---------------------------------------------------------------------------

def apply_filter(
    predictions: dict[str, PredictionRow],
    hmm_positive_set: set[str],
    min_tm: int,
    min_confidence: float,
) -> set[str]:
    """Return IDs that pass all three gates (mirrors stage-02 awk + HMM filter).

    Gates applied in order (AND):
      1. seq_id is in hmm_positive_set  (HMM-first GPCR filter)
      2. tm_count >= min_tm             ($5+0 >= min_tm)
      3. tm_confidence >= min_confidence ($3+0 >= min_conf)
    """
    passed: set[str] = set()
    for seq_id, row in predictions.items():
        if (
            seq_id in hmm_positive_set
            and row.tm_count >= min_tm
            and row.tm_confidence >= min_confidence
        ):
            passed.add(seq_id)
    return passed


# ---------------------------------------------------------------------------
# 4. emit_scan_record
# ---------------------------------------------------------------------------

_TSV_COLUMNS = ("seq_id", "gpcr_positive", "tm_count", "tm_confidence",
                "passed_gate", "notes")


def emit_scan_record(
    predictions: dict[str, PredictionRow],
    hmm_positive_set: set[str],
    passed_ids: set[str],
    out_tsv: Path,
    min_tm: int = 6,
    min_confidence: float = 0.5,
) -> None:
    """Write the 6-column scan record TSV — one row per sequence in predictions.

    Columns (fixed order per AC spec):
        seq_id, gpcr_positive, tm_count, tm_confidence, passed_gate, notes
    """
    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(_TSV_COLUMNS)
        for seq_id, row in sorted(predictions.items()):
            gpcr_positive = 1 if seq_id in hmm_positive_set else 0
            passed_gate = 1 if seq_id in passed_ids else 0
            notes = _build_notes(
                seq_id, row, hmm_positive_set, passed_ids, min_tm, min_confidence
            )
            writer.writerow([
                seq_id,
                gpcr_positive,
                row.tm_count,
                f"{row.tm_confidence:.4f}",
                passed_gate,
                notes,
            ])


def _build_notes(
    seq_id: str,
    row: PredictionRow,
    hmm_positive_set: set[str],
    passed_ids: set[str],
    min_tm: int = 6,
    min_confidence: float = 0.5,
) -> str:
    """Human-readable reason why a sequence was rejected (empty if passed).

    Checks filters in order (AND logic):
      1. HMM positivity
      2. TM count against min_tm threshold
      3. TM confidence against min_confidence threshold
    """
    if seq_id in passed_ids:
        return ""
    if seq_id not in hmm_positive_set:
        return "hmm_negative"
    if row.tm_count < min_tm:
        return f"tm_count={row.tm_count}_below_min={min_tm}"
    if row.tm_confidence < min_confidence:
        return f"tm_confidence={row.tm_confidence:.2f}_below_min={min_confidence:.2f}"
    return ""


# ---------------------------------------------------------------------------
# 5. emit_candidate_fasta
# ---------------------------------------------------------------------------

def emit_candidate_fasta(
    input_fa: Path,
    passed_ids: set[str],
    out_fa: Path,
) -> int:
    """Write FASTA subset containing only sequences whose ID is in passed_ids.

    Matches on the first whitespace-delimited token of each header line
    (standard FASTA ID convention, same as seqtk subseq).

    Returns the count of sequences written.
    """
    out_fa = Path(out_fa)
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with out_fa.open("w") as fh:
        for rec in SeqIO.parse(str(input_fa), "fasta"):
            if rec.id in passed_ids:
                SeqIO.write(rec, fh, "fasta")
                written += 1
    return written


# ---------------------------------------------------------------------------
# 6. check_outputs_exist (idempotency)
# ---------------------------------------------------------------------------

def check_outputs_exist(
    out_fa: Path,
    out_tsv: Path,
    force: bool,
) -> Optional[str]:
    """Return 'skipped' if outputs appear complete and force=False; else None.

    "Complete" means:
      - both out_fa and out_tsv exist, AND
      - either out_fa contains at least one sequence (>0 bytes with a '>') OR
        a sibling zero_candidates flag file exists next to out_fa.

    force=True always returns None (caller must re-run).
    """
    if force:
        return None

    out_fa = Path(out_fa)
    out_tsv = Path(out_tsv)

    if not out_fa.exists() or not out_tsv.exists():
        return None

    # Check for zero_candidates flag (stem = out_fa without .chemo_candidates.fa)
    stem = out_fa.name
    for suffix in (".chemo_candidates.fa", ".fa"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    zero_flag = out_fa.parent / f"{stem}.zero_candidates"

    has_candidates = out_fa.stat().st_size > 0 and _fasta_has_sequences(out_fa)
    has_zero_flag = zero_flag.exists()

    if has_candidates or has_zero_flag:
        return "skipped"

    return None


def _fasta_has_sequences(fa: Path) -> bool:
    """Return True if file contains at least one FASTA header."""
    try:
        with fa.open() as fh:
            for line in fh:
                if line.startswith(">"):
                    return True
    except OSError:
        pass
    return False


# ---------------------------------------------------------------------------
# 7. CLI
# ---------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    """Build and return the argument parser for the filter helper."""
    p = argparse.ArgumentParser(
        description="Apply HMM + TMbed filters and emit outputs for one proteome scan.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--tmbed-prediction", required=True, type=Path, metavar="FILE",
        help="TMbed prediction file (tab-separated, ≥5 fields per row)",
    )
    p.add_argument(
        "--hmm-positive-ids", required=True, type=Path, metavar="FILE",
        help="Text file (one ID per line) of HMM-positive sequences",
    )
    p.add_argument(
        "--input-fa", required=True, type=Path, metavar="FASTA",
        help="Length-filtered GPCR-positive input FASTA (TMbed was run on this)",
    )
    p.add_argument(
        "--out-fa", required=True, type=Path, metavar="FASTA",
        help="Output candidate FASTA path",
    )
    p.add_argument(
        "--out-tsv", required=True, type=Path, metavar="TSV",
        help="Output scan_record TSV path",
    )
    p.add_argument(
        "--min-tm", type=int, default=6, metavar="N",
        help="Minimum TM regions required (≥N passes)",
    )
    p.add_argument(
        "--min-confidence", type=float, default=0.5, metavar="FLOAT",
        help="Minimum TMbed confidence score required (≥threshold passes)",
    )
    p.add_argument(
        "--force", action="store_true",
        help="Re-run even if output files already exist",
    )
    return p


def _load_hmm_positive_set(ids_file: Path) -> set[str]:
    """Read a one-ID-per-line text file; return a set of stripped IDs."""
    ids: set[str] = set()
    if not ids_file.exists():
        return ids
    with ids_file.open() as fh:
        for line in fh:
            stripped = line.strip()
            if stripped:
                ids.add(stripped)
    return ids


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    # Idempotency check
    skip = check_outputs_exist(args.out_fa, args.out_tsv, force=args.force)
    if skip == "skipped":
        print(
            f"[scan_filter] skipped: outputs already exist at {args.out_fa}",
            file=sys.stderr,
        )
        return 0

    # Load TMbed predictions
    predictions = parse_tmbed_prediction(args.tmbed_prediction)

    # Load HMM-positive IDs
    hmm_positive_set = _load_hmm_positive_set(args.hmm_positive_ids)

    # Apply filter
    passed_ids = apply_filter(
        predictions, hmm_positive_set,
        min_tm=args.min_tm,
        min_confidence=args.min_confidence,
    )

    # Emit scan record TSV (all sequences that went through TMbed)
    emit_scan_record(
        predictions=predictions,
        hmm_positive_set=hmm_positive_set,
        passed_ids=passed_ids,
        out_tsv=args.out_tsv,
        min_tm=args.min_tm,
        min_confidence=args.min_confidence,
    )

    # Emit candidate FASTA
    n_written = emit_candidate_fasta(args.input_fa, passed_ids, args.out_fa)

    # If zero candidates, write flag file to allow idempotent skipping next run
    if n_written == 0:
        stem = args.out_fa.name
        for suffix in (".chemo_candidates.fa", ".fa"):
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
                break
        zero_flag = args.out_fa.parent / f"{stem}.zero_candidates"
        zero_flag.touch()
        print(
            f"[scan_filter] WARNING: 0 candidates passed all gates; "
            f"wrote zero_candidates flag at {zero_flag}",
            file=sys.stderr,
        )

    print(
        f"[scan_filter] {n_written} candidates written to {args.out_fa}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
