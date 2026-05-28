"""Tests for scripts/_scan_proteome_filter_helper.py.

Covers 7+ acceptance-criteria test groups:
  1. sanitize_sample_name matches existing convention
  2. parse_tmbed_prediction — correct parsing of TMbed prediction TSV format
  3. apply_filter — gates (≥min_tm AND ≥min_confidence AND in hmm_positive_set)
  4. emit_scan_record — 6-column TSV, one row per input sequence
  5. emit_candidate_fasta — FASTA subset for passed IDs
  6. Idempotency — existing outputs with non-zero candidates returns "skipped"
  7. CLI argparse for the Python helper

Import style matches existing unit tests: conftest.py adds scripts/ to sys.path.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pytest

import _scan_proteome_filter_helper as helper


# ---------------------------------------------------------------------------
# 1. sanitize_sample_name — must match build_braker4_samples_csv convention
# ---------------------------------------------------------------------------

class TestSanitizeSampleName:
    """[^A-Za-z0-9_]+ → '_', dedupe runs, strip leading/trailing '_'."""

    def test_space_becomes_underscore(self) -> None:
        assert helper.sanitize_sample_name("Berghia stephanieae") == "Berghia_stephanieae"

    def test_already_clean_unchanged(self) -> None:
        assert helper.sanitize_sample_name("1287507_berghia_stephanieae") == "1287507_berghia_stephanieae"

    def test_hyphen_becomes_underscore(self) -> None:
        assert helper.sanitize_sample_name("Aeolidia-papillosa") == "Aeolidia_papillosa"

    def test_multiple_spaces_collapsed(self) -> None:
        assert helper.sanitize_sample_name("Genus   species") == "Genus_species"

    def test_leading_trailing_whitespace_stripped(self) -> None:
        assert helper.sanitize_sample_name("  trimme  ") == "trimme"

    def test_dot_and_comma_replaced(self) -> None:
        result = helper.sanitize_sample_name("Ctena cf. galapagana, sp. nov")
        assert result == "Ctena_cf_galapagana_sp_nov"


# ---------------------------------------------------------------------------
# 2. parse_tmbed_prediction — TMbed TSV format
#    col 1 = id, col 3 = confidence, col 5 = TM count (1-based column indices)
#    as used in stage 02 line 67 awk: $1=id, $3=confidence, $5=tm_count
# ---------------------------------------------------------------------------

class TestParseTmbedPrediction:
    """parse_tmbed_prediction(path) → dict[seq_id, PredictionRow]."""

    def _write_prediction(self, path: Path, rows: list[str]) -> None:
        path.write_text("\n".join(rows) + "\n")

    def test_basic_parse_returns_all_rows(self, tmp_path: Path) -> None:
        pred = tmp_path / "prediction"
        self._write_prediction(pred, [
            "seq1\tTM\t0.9\tsome\t7\textra",
            "seq2\tTM\t0.4\tsome\t3\textra",
            "seq3\tSP\t0.8\tsome\t0\textra",
        ])
        result = helper.parse_tmbed_prediction(pred)
        assert len(result) == 3
        assert result["seq1"].tm_count == 7
        assert result["seq1"].tm_confidence == pytest.approx(0.9)
        assert result["seq2"].tm_count == 3
        assert result["seq3"].tm_count == 0

    def test_rows_with_fewer_than_5_fields_skipped(self, tmp_path: Path) -> None:
        pred = tmp_path / "prediction"
        self._write_prediction(pred, [
            "seq_ok\tTM\t0.9\tsome\t6",          # exactly 5 — should be kept
            "seq_short\tTM\t0.9\t3",              # only 4 fields — should be skipped
            "seq_also_ok\tTM\t0.7\tsome\t8\tex",  # 6 fields — should be kept
        ])
        result = helper.parse_tmbed_prediction(pred)
        assert "seq_ok" in result
        assert "seq_also_ok" in result
        assert "seq_short" not in result

    def test_empty_prediction_file(self, tmp_path: Path) -> None:
        pred = tmp_path / "prediction"
        pred.write_text("")
        result = helper.parse_tmbed_prediction(pred)
        assert result == {}

    def test_fields_are_tab_separated(self, tmp_path: Path) -> None:
        pred = tmp_path / "prediction"
        pred.write_text("myseq\tTM\t0.75\tother\t5\n")
        result = helper.parse_tmbed_prediction(pred)
        assert "myseq" in result
        assert result["myseq"].tm_count == 5
        assert result["myseq"].tm_confidence == pytest.approx(0.75)


# ---------------------------------------------------------------------------
# 3. apply_filter — gates: ≥min_tm AND ≥min_confidence AND in hmm_positive_set
# ---------------------------------------------------------------------------

class TestApplyFilter:
    """apply_filter(predictions, hmm_positive_set, min_tm, min_confidence)
    → set[str] of IDs that pass all gates."""

    def _make_predictions(self) -> dict:
        return {
            "high_tm_high_conf":  helper.PredictionRow(tm_count=7, tm_confidence=0.9),
            "low_tm_high_conf":   helper.PredictionRow(tm_count=3, tm_confidence=0.9),
            "high_tm_low_conf":   helper.PredictionRow(tm_count=7, tm_confidence=0.3),
            "exact_min_tm":       helper.PredictionRow(tm_count=6, tm_confidence=0.5),
            "not_in_hmm_set":     helper.PredictionRow(tm_count=8, tm_confidence=0.95),
        }

    def test_passes_when_all_gates_met(self) -> None:
        preds = self._make_predictions()
        hmm_pos = {"high_tm_high_conf", "low_tm_high_conf", "high_tm_low_conf",
                   "exact_min_tm", "extra"}
        passed = helper.apply_filter(preds, hmm_pos, min_tm=6, min_confidence=0.5)
        assert "high_tm_high_conf" in passed
        assert "exact_min_tm" in passed

    def test_fails_when_tm_below_min(self) -> None:
        preds = self._make_predictions()
        hmm_pos = {"low_tm_high_conf"}
        passed = helper.apply_filter(preds, hmm_pos, min_tm=6, min_confidence=0.5)
        assert "low_tm_high_conf" not in passed

    def test_fails_when_confidence_below_min(self) -> None:
        preds = self._make_predictions()
        hmm_pos = {"high_tm_low_conf"}
        passed = helper.apply_filter(preds, hmm_pos, min_tm=6, min_confidence=0.5)
        assert "high_tm_low_conf" not in passed

    def test_fails_when_not_in_hmm_positive_set(self) -> None:
        preds = self._make_predictions()
        # "not_in_hmm_set" has good TM + confidence but NOT in hmm_pos
        hmm_pos = {"high_tm_high_conf"}  # only this one in set
        passed = helper.apply_filter(preds, hmm_pos, min_tm=6, min_confidence=0.5)
        assert "not_in_hmm_set" not in passed

    def test_empty_predictions_returns_empty(self) -> None:
        passed = helper.apply_filter({}, {"anything"}, min_tm=6, min_confidence=0.5)
        assert passed == set()

    def test_empty_hmm_set_returns_empty(self) -> None:
        preds = self._make_predictions()
        passed = helper.apply_filter(preds, set(), min_tm=6, min_confidence=0.5)
        assert passed == set()


# ---------------------------------------------------------------------------
# 4. emit_scan_record — 6-column TSV, one row per input sequence
# ---------------------------------------------------------------------------

class TestEmitScanRecord:
    """emit_scan_record writes exactly 6 columns per the AC spec:
    seq_id, gpcr_positive, tm_count, tm_confidence, passed_gate, notes.
    One row per sequence in predictions; HMM-negative IDs get defaults.
    """

    def _make_predictions(self) -> dict:
        return {
            "s1": helper.PredictionRow(tm_count=7, tm_confidence=0.9),
            "s2": helper.PredictionRow(tm_count=2, tm_confidence=0.3),
        }

    def test_correct_columns_present(self, tmp_path: Path) -> None:
        out = tmp_path / "scan_record.tsv"
        helper.emit_scan_record(
            predictions=self._make_predictions(),
            hmm_positive_set={"s1", "s2"},
            passed_ids={"s1"},
            out_tsv=out,
        )
        with out.open() as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
        assert set(reader.fieldnames) == {
            "seq_id", "gpcr_positive", "tm_count", "tm_confidence",
            "passed_gate", "notes"
        }

    def test_passed_id_has_passed_gate_1(self, tmp_path: Path) -> None:
        out = tmp_path / "scan_record.tsv"
        helper.emit_scan_record(
            predictions=self._make_predictions(),
            hmm_positive_set={"s1", "s2"},
            passed_ids={"s1"},
            out_tsv=out,
        )
        rows = _read_tsv(out)
        row = {r["seq_id"]: r for r in rows}
        assert row["s1"]["passed_gate"] == "1"
        assert row["s2"]["passed_gate"] == "0"

    def test_gpcr_positive_reflects_hmm_set(self, tmp_path: Path) -> None:
        out = tmp_path / "scan_record.tsv"
        helper.emit_scan_record(
            predictions=self._make_predictions(),
            hmm_positive_set={"s1"},  # s2 NOT in hmm set
            passed_ids={"s1"},
            out_tsv=out,
        )
        rows = _read_tsv(out)
        row = {r["seq_id"]: r for r in rows}
        assert row["s1"]["gpcr_positive"] == "1"
        assert row["s2"]["gpcr_positive"] == "0"

    def test_one_row_per_sequence(self, tmp_path: Path) -> None:
        out = tmp_path / "scan_record.tsv"
        helper.emit_scan_record(
            predictions=self._make_predictions(),
            hmm_positive_set={"s1"},
            passed_ids={"s1"},
            out_tsv=out,
        )
        rows = _read_tsv(out)
        assert len(rows) == 2  # s1 + s2


# ---------------------------------------------------------------------------
# 4b. _build_notes — diagnostic notes function
# ---------------------------------------------------------------------------

class TestBuildNotes:
    """_build_notes generates human-readable diagnostics about rejection reasons."""

    def test_build_notes_returns_empty_for_passed_sequence(self) -> None:
        """A sequence in passed_ids gets empty notes."""
        row = helper.PredictionRow(tm_count=3, tm_confidence=0.1)
        notes = helper._build_notes(
            "seq1", row,
            hmm_positive_set={"seq1"},
            passed_ids={"seq1"},
            min_tm=6, min_confidence=0.5,
        )
        assert notes == ""

    def test_build_notes_reports_hmm_negative(self) -> None:
        """A sequence not in hmm_positive_set gets 'hmm_negative'."""
        row = helper.PredictionRow(tm_count=7, tm_confidence=0.9)
        notes = helper._build_notes(
            "seq1", row,
            hmm_positive_set=set(),  # seq1 NOT in set
            passed_ids=set(),
            min_tm=6, min_confidence=0.5,
        )
        assert notes == "hmm_negative"

    def test_build_notes_uses_actual_min_tm_threshold(self) -> None:
        """When min_tm=7 and tm_count=6, should report 'tm_count=6_below_min=7'."""
        row = helper.PredictionRow(tm_count=6, tm_confidence=0.9)
        notes = helper._build_notes(
            "seq1", row,
            hmm_positive_set={"seq1"},  # HMM positive
            passed_ids=set(),  # not passed
            min_tm=7,  # threshold is 7
            min_confidence=0.5,
        )
        assert notes == "tm_count=6_below_min=7"

    def test_build_notes_reports_confidence_failure(self) -> None:
        """When min_confidence=0.5 and tm_confidence=0.3, should report confidence failure."""
        row = helper.PredictionRow(tm_count=8, tm_confidence=0.3)
        notes = helper._build_notes(
            "seq1", row,
            hmm_positive_set={"seq1"},
            passed_ids=set(),
            min_tm=6,
            min_confidence=0.5,
        )
        assert notes == "tm_confidence=0.30_below_min=0.50"


# ---------------------------------------------------------------------------
# 5. emit_candidate_fasta — FASTA subset for passed IDs
# ---------------------------------------------------------------------------

class TestEmitCandidateFasta:
    """emit_candidate_fasta(input_fa, passed_ids, out_fa) writes FASTA subset."""

    def test_subset_contains_only_passed_ids(self, tmp_path: Path) -> None:
        src = tmp_path / "input.fa"
        src.write_text(">seqA\nMSEQA\n>seqB\nMSEQB\n>seqC\nMSEQC\n")
        out = tmp_path / "out.fa"
        helper.emit_candidate_fasta(src, {"seqA", "seqC"}, out)
        content = out.read_text()
        assert ">seqA" in content
        assert ">seqC" in content
        assert ">seqB" not in content

    def test_sequences_preserved_correctly(self, tmp_path: Path) -> None:
        src = tmp_path / "input.fa"
        src.write_text(">prot1\nMSEQ\n>prot2\nACGT\n")
        out = tmp_path / "out.fa"
        helper.emit_candidate_fasta(src, {"prot1"}, out)
        content = out.read_text()
        assert "MSEQ" in content
        assert "ACGT" not in content

    def test_empty_passed_ids_writes_empty_fasta(self, tmp_path: Path) -> None:
        src = tmp_path / "input.fa"
        src.write_text(">s1\nM\n>s2\nA\n")
        out = tmp_path / "out.fa"
        helper.emit_candidate_fasta(src, set(), out)
        # File exists but no sequences
        content = out.read_text()
        assert ">" not in content

    def test_creates_output_parent_dirs(self, tmp_path: Path) -> None:
        src = tmp_path / "input.fa"
        src.write_text(">x\nM\n")
        out = tmp_path / "deep" / "nested" / "out.fa"
        helper.emit_candidate_fasta(src, {"x"}, out)
        assert out.exists()

    def test_fasta_id_matched_by_first_token(self, tmp_path: Path) -> None:
        """Header '>seq1 some description' should match on 'seq1'."""
        src = tmp_path / "input.fa"
        src.write_text(">seq1 description here\nMSTAR\n>seq2 other\nMFOO\n")
        out = tmp_path / "out.fa"
        helper.emit_candidate_fasta(src, {"seq1"}, out)
        content = out.read_text()
        assert ">seq1" in content
        assert "MSTAR" in content
        assert ">seq2" not in content


# ---------------------------------------------------------------------------
# 6. Idempotency — existing outputs with non-zero candidates returns "skipped"
# ---------------------------------------------------------------------------

class TestIdempotency:
    """check_outputs_exist returns 'skipped' when both target files exist and
    candidate FASTA has at least one sequence (or zero_candidates flag exists)."""

    def test_returns_skipped_when_both_outputs_exist_with_candidates(
        self, tmp_path: Path
    ) -> None:
        out_fa = tmp_path / "1287507_berghia.chemo_candidates.fa"
        out_tsv = tmp_path / "1287507_berghia.scan_record.tsv"
        out_fa.write_text(">cand1\nMSEQ\n")
        out_tsv.write_text("seq_id\tpassed_gate\n")
        result = helper.check_outputs_exist(out_fa, out_tsv, force=False)
        assert result == "skipped"

    def test_returns_skipped_when_zero_candidates_flag_exists(
        self, tmp_path: Path
    ) -> None:
        out_fa = tmp_path / "42_empty.chemo_candidates.fa"
        out_tsv = tmp_path / "42_empty.scan_record.tsv"
        zero_flag = tmp_path / "42_empty.zero_candidates"
        out_fa.write_text("")  # empty FASTA
        out_tsv.write_text("seq_id\tpassed_gate\n")
        zero_flag.touch()
        result = helper.check_outputs_exist(out_fa, out_tsv, force=False)
        assert result == "skipped"

    def test_returns_none_when_outputs_missing(self, tmp_path: Path) -> None:
        out_fa = tmp_path / "missing.chemo_candidates.fa"
        out_tsv = tmp_path / "missing.scan_record.tsv"
        result = helper.check_outputs_exist(out_fa, out_tsv, force=False)
        assert result is None

    def test_force_overrides_skip(self, tmp_path: Path) -> None:
        out_fa = tmp_path / "1_A.chemo_candidates.fa"
        out_tsv = tmp_path / "1_A.scan_record.tsv"
        out_fa.write_text(">cand1\nM\n")
        out_tsv.write_text("seq_id\tpassed_gate\n")
        result = helper.check_outputs_exist(out_fa, out_tsv, force=True)
        assert result is None  # force=True → do NOT skip

    def test_returns_none_when_fasta_empty_and_no_flag(
        self, tmp_path: Path
    ) -> None:
        """Empty FASTA + no zero_candidates flag means run didn't finish."""
        out_fa = tmp_path / "1_A.chemo_candidates.fa"
        out_tsv = tmp_path / "1_A.scan_record.tsv"
        out_fa.write_text("")
        out_tsv.write_text("seq_id\tpassed_gate\n")
        result = helper.check_outputs_exist(out_fa, out_tsv, force=False)
        assert result is None


# ---------------------------------------------------------------------------
# 7. CLI argparse
# ---------------------------------------------------------------------------

class TestCLI:
    """The helper accepts --tmbed-prediction, --hmm-positive-ids, --input-fa,
    --out-fa, --out-tsv, --min-tm, --min-confidence, --force."""

    def test_required_args_parsed(self, tmp_path: Path) -> None:
        pred = tmp_path / "prediction"
        ids = tmp_path / "ids.txt"
        fa = tmp_path / "input.fa"
        out_fa = tmp_path / "out.fa"
        out_tsv = tmp_path / "out.tsv"
        # Create dummy files so argparse type=Path exists-check passes
        pred.write_text("")
        ids.write_text("")
        fa.write_text("")
        args = helper.build_arg_parser().parse_args([
            "--tmbed-prediction", str(pred),
            "--hmm-positive-ids", str(ids),
            "--input-fa", str(fa),
            "--out-fa", str(out_fa),
            "--out-tsv", str(out_tsv),
        ])
        assert args.tmbed_prediction == pred
        assert args.hmm_positive_ids == ids
        assert args.input_fa == fa
        assert args.out_fa == out_fa
        assert args.out_tsv == out_tsv
        assert args.min_tm == 6       # default
        assert args.min_confidence == pytest.approx(0.5)  # default
        assert args.force is False    # default

    def test_optional_args_override_defaults(self, tmp_path: Path) -> None:
        pred = tmp_path / "prediction"
        ids = tmp_path / "ids.txt"
        fa = tmp_path / "input.fa"
        pred.write_text(""); ids.write_text(""); fa.write_text("")
        args = helper.build_arg_parser().parse_args([
            "--tmbed-prediction", str(pred),
            "--hmm-positive-ids", str(ids),
            "--input-fa", str(fa),
            "--out-fa", str(tmp_path / "out.fa"),
            "--out-tsv", str(tmp_path / "out.tsv"),
            "--min-tm", "7",
            "--min-confidence", "0.7",
            "--force",
        ])
        assert args.min_tm == 7
        assert args.min_confidence == pytest.approx(0.7)
        assert args.force is True


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_tsv(path: Path) -> list[dict]:
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))
