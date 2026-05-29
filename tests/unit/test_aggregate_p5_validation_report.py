"""Tests for scripts/aggregate_p5_validation_report.py.

P5 Phase 1a validation infrastructure — aggregation + report generator.

Reads:
  - scan/<taxid>_<binomial>.scan_record.tsv (per-species scan stats)
  - classify/class_phase1a.tsv  (P1 classifier output for Phase 1a candidates)
  - classify/class_berghia.tsv  (P1 classifier output for Berghia 888)
  - pools/pool_build_report.json (P2 pool builder output)

Emits:
  - validation_report.json
  - validation_report.md

Idempotent: skips if outputs already exist unless --force.
Partial: works even if classify or pool outputs are missing.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

import aggregate_p5_validation_report as agg


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_scan_record(path: Path, rows: list[dict]) -> None:
    """Write a scan_record.tsv with header: seq_id, gpcr_positive, tm_count,
    tm_confidence, passed_gate, notes.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["seq_id", "gpcr_positive", "tm_count", "tm_confidence",
            "passed_gate", "notes"]
    with path.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def _write_class_tsv(path: Path, rows: list[dict]) -> None:
    """Write a classify_gpcr_by_class.py output TSV:
    seq_id, class, evidence_pfam, evidence_family_hmm, top_evalue.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["seq_id", "class", "evidence_pfam", "evidence_family_hmm", "top_evalue"]
    with path.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def _write_pool_report(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(data, f, indent=2)


MINIMAL_POOL_REPORT = {
    "class_A": {
        "n_input_refs": 20,
        "n_after_cdhit": 18,
        "n_berghia_included": 5,
        "n_must_include": 3,
        "n_subsampled": 10,
        "n_output": 15,
        "n_refs_target": 15,
        "total_budget_for_class": 20,
    },
    "class_B": {
        "n_input_refs": 5,
        "n_after_cdhit": 5,
        "n_berghia_included": 1,
        "n_must_include": 1,
        "n_subsampled": 2,
        "n_output": 3,
        "n_refs_target": 5,
        "total_budget_for_class": 10,
    },
    "class_C": {
        "n_input_refs": 3,
        "n_after_cdhit": 3,
        "n_berghia_included": 0,
        "n_must_include": 0,
        "n_subsampled": 2,
        "n_output": 2,
        "n_refs_target": 5,
        "total_budget_for_class": 10,
    },
    "class_F": {
        "n_input_refs": 2,
        "n_after_cdhit": 2,
        "n_berghia_included": 0,
        "n_must_include": 0,
        "n_subsampled": 1,
        "n_output": 1,
        "n_refs_target": 5,
        "total_budget_for_class": 10,
    },
    "unclassified": {"n": 7},
    "berghia_included": {"taxid": 1287507, "n_per_class": {"A": 5, "B": 1, "C": 0, "F": 0}, "n_total": 6},
}


# ---------------------------------------------------------------------------
# parse_scan_records
# ---------------------------------------------------------------------------

class TestParseScanRecords:
    def test_counts_correct_for_single_file(self, tmp_path: Path) -> None:
        scan_dir = tmp_path / "scan"
        _write_scan_record(
            scan_dir / "100_Aaa_bbb.scan_record.tsv",
            [
                {"seq_id": "s1", "gpcr_positive": "1", "passed_gate": "1"},
                {"seq_id": "s2", "gpcr_positive": "1", "passed_gate": "0"},
                {"seq_id": "s3", "gpcr_positive": "0", "passed_gate": "0"},
            ],
        )
        records = agg.parse_scan_records(scan_dir)
        assert len(records) == 1
        stats = records[0]
        assert stats["n_input"] == 3
        assert stats["n_hmm_positive"] == 2
        assert stats["n_passed_gate"] == 1

    def test_taxid_and_binomial_extracted_from_filename(
        self, tmp_path: Path
    ) -> None:
        scan_dir = tmp_path / "scan"
        _write_scan_record(
            scan_dir / "6183_Schistosoma_mansoni.scan_record.tsv",
            [{"seq_id": "x", "gpcr_positive": "1", "passed_gate": "1"}],
        )
        records = agg.parse_scan_records(scan_dir)
        assert records[0]["taxid"] == "6183"
        assert "Schistosoma" in records[0]["binomial"]

    def test_empty_scan_dir_returns_empty_list(self, tmp_path: Path) -> None:
        records = agg.parse_scan_records(tmp_path / "scan_missing")
        assert records == []

    def test_multiple_files_aggregated_separately(self, tmp_path: Path) -> None:
        scan_dir = tmp_path / "scan"
        for i in range(3):
            _write_scan_record(
                scan_dir / f"10{i}_Species_sp{i}.scan_record.tsv",
                [{"seq_id": f"s{j}", "gpcr_positive": "1", "passed_gate": "1"}
                 for j in range(i + 1)],
            )
        records = agg.parse_scan_records(scan_dir)
        assert len(records) == 3
        total_input = sum(r["n_input"] for r in records)
        assert total_input == 1 + 2 + 3


# ---------------------------------------------------------------------------
# parse_class_distribution
# ---------------------------------------------------------------------------

class TestParseClassDistribution:
    def test_counts_class_labels(self, tmp_path: Path) -> None:
        tsv = tmp_path / "class.tsv"
        _write_class_tsv(tsv, [
            {"seq_id": "a", "class": "A"},
            {"seq_id": "b", "class": "A"},
            {"seq_id": "c", "class": "B"},
            {"seq_id": "d", "class": "unclassified"},
        ])
        dist = agg.parse_class_distribution(tsv)
        assert dist["A"] == 2
        assert dist["B"] == 1
        assert dist["unclassified"] == 1
        assert dist.get("C", 0) == 0

    def test_missing_file_returns_none(self, tmp_path: Path) -> None:
        result = agg.parse_class_distribution(tmp_path / "missing.tsv")
        assert result is None

    def test_empty_file_body_returns_zero_counts(self, tmp_path: Path) -> None:
        tsv = tmp_path / "class.tsv"
        _write_class_tsv(tsv, [])
        dist = agg.parse_class_distribution(tsv)
        assert dist is not None
        assert sum(dist.values()) == 0


# ---------------------------------------------------------------------------
# load_pool_report
# ---------------------------------------------------------------------------

class TestLoadPoolReport:
    def test_returns_dict_for_valid_json(self, tmp_path: Path) -> None:
        rp = tmp_path / "pools" / "pool_build_report.json"
        _write_pool_report(rp, MINIMAL_POOL_REPORT)
        data = agg.load_pool_report(rp)
        assert data is not None
        assert "class_A" in data

    def test_returns_none_for_missing_file(self, tmp_path: Path) -> None:
        result = agg.load_pool_report(tmp_path / "missing.json")
        assert result is None


# ---------------------------------------------------------------------------
# build_validation_report
# ---------------------------------------------------------------------------

class TestBuildValidationReport:
    def _setup(self, tmp_path: Path) -> dict:
        """Build a minimal valid dataset and call build_validation_report."""
        scan_dir = tmp_path / "scan"
        classify_dir = tmp_path / "classify"
        pool_dir = tmp_path / "pools"

        # Two species in scan
        for taxid, sp in [("100", "Alpha_beta"), ("101", "Gamma_delta")]:
            _write_scan_record(
                scan_dir / f"{taxid}_{sp}.scan_record.tsv",
                [
                    {"seq_id": "s1", "gpcr_positive": "1", "passed_gate": "1"},
                    {"seq_id": "s2", "gpcr_positive": "0", "passed_gate": "0"},
                ],
            )

        _write_class_tsv(
            classify_dir / "class_phase1a.tsv",
            [{"seq_id": "s1", "class": "A"}, {"seq_id": "s2", "class": "B"}],
        )
        _write_class_tsv(
            classify_dir / "class_berghia.tsv",
            [{"seq_id": "b1", "class": "A"}, {"seq_id": "b2", "class": "unclassified"}],
        )
        _write_pool_report(pool_dir / "pool_build_report.json", MINIMAL_POOL_REPORT)

        return agg.build_validation_report(
            scan_dir=scan_dir,
            classify_dir=classify_dir,
            pool_report_path=pool_dir / "pool_build_report.json",
        )

    def test_report_has_per_species_section(self, tmp_path: Path) -> None:
        report = self._setup(tmp_path)
        assert "per_species" in report
        assert len(report["per_species"]) == 2

    def test_report_has_class_distribution(self, tmp_path: Path) -> None:
        report = self._setup(tmp_path)
        assert "class_distribution" in report
        assert "phase1a" in report["class_distribution"]
        assert "berghia" in report["class_distribution"]

    def test_report_has_pool_composition(self, tmp_path: Path) -> None:
        report = self._setup(tmp_path)
        assert "pool_composition" in report

    def test_report_has_summary(self, tmp_path: Path) -> None:
        report = self._setup(tmp_path)
        assert "summary" in report
        assert "n_species" in report["summary"]
        assert report["summary"]["n_species"] == 2

    def test_partial_run_missing_classify(self, tmp_path: Path) -> None:
        """classify dir absent → report still emits with missing=True flags."""
        scan_dir = tmp_path / "scan"
        _write_scan_record(
            scan_dir / "100_Alpha_beta.scan_record.tsv",
            [{"seq_id": "s1", "gpcr_positive": "1", "passed_gate": "1"}],
        )
        report = agg.build_validation_report(
            scan_dir=scan_dir,
            classify_dir=tmp_path / "missing_classify",
            pool_report_path=tmp_path / "missing_pools" / "pool_build_report.json",
        )
        assert report["class_distribution"]["phase1a"] is None
        assert report["class_distribution"]["berghia"] is None
        assert report["pool_composition"] is None

    def test_per_species_stats_are_correct(self, tmp_path: Path) -> None:
        scan_dir = tmp_path / "scan"
        _write_scan_record(
            scan_dir / "200_Yyy_zzz.scan_record.tsv",
            [
                {"seq_id": "s1", "gpcr_positive": "1", "passed_gate": "1"},
                {"seq_id": "s2", "gpcr_positive": "1", "passed_gate": "0"},
                {"seq_id": "s3", "gpcr_positive": "0", "passed_gate": "0"},
            ],
        )
        report = agg.build_validation_report(
            scan_dir=scan_dir,
            classify_dir=tmp_path / "classify_missing",
            pool_report_path=tmp_path / "pool_missing.json",
        )
        sp = report["per_species"][0]
        assert sp["n_input"] == 3
        assert sp["n_hmm_positive"] == 2
        assert sp["n_passed_gate"] == 1


# ---------------------------------------------------------------------------
# write_json_report + write_md_report
# ---------------------------------------------------------------------------

class TestWriteReports:
    def _minimal_report(self) -> dict:
        return {
            "per_species": [
                {"taxid": "100", "binomial": "Alpha_beta",
                 "n_input": 5, "n_hmm_positive": 3, "n_passed_gate": 2},
            ],
            "class_distribution": {
                "phase1a": {"A": 2, "B": 1, "C": 0, "F": 0, "unclassified": 0},
                "berghia": {"A": 5, "B": 2, "C": 0, "F": 0, "unclassified": 1},
                "combined": {"A": 7, "B": 3, "C": 0, "F": 0, "unclassified": 1},
            },
            "pool_composition": {"class_A": {"n_output": 15}},
            "summary": {"n_species": 1, "total_candidates": 2},
        }

    def test_json_output_is_valid_json(self, tmp_path: Path) -> None:
        out = tmp_path / "report.json"
        agg.write_json_report(self._minimal_report(), out)
        data = json.loads(out.read_text())
        assert "per_species" in data

    def test_md_output_contains_section_headers(self, tmp_path: Path) -> None:
        out = tmp_path / "report.md"
        agg.write_md_report(self._minimal_report(), out)
        text = out.read_text()
        assert "Per-species" in text or "per-species" in text.lower()
        assert "Class distribution" in text or "class distribution" in text.lower()
        assert "Pool composition" in text or "pool composition" in text.lower()
        assert "summary" in text.lower() or "Summary" in text

    def test_md_output_contains_species_data(self, tmp_path: Path) -> None:
        out = tmp_path / "report.md"
        agg.write_md_report(self._minimal_report(), out)
        text = out.read_text()
        assert "100" in text  # taxid present


# ---------------------------------------------------------------------------
# main() idempotency
# ---------------------------------------------------------------------------

class TestMainIdempotency:
    def _make_valid_inputs(self, tmp_path: Path) -> None:
        scan_dir = tmp_path / "scan"
        classify_dir = tmp_path / "classify"
        pool_dir = tmp_path / "pools"
        _write_scan_record(
            scan_dir / "100_Aaa_bbb.scan_record.tsv",
            [{"seq_id": "s1", "gpcr_positive": "1", "passed_gate": "1"}],
        )
        _write_class_tsv(
            classify_dir / "class_phase1a.tsv",
            [{"seq_id": "s1", "class": "A"}],
        )
        _write_class_tsv(
            classify_dir / "class_berghia.tsv",
            [{"seq_id": "b1", "class": "A"}],
        )
        _write_pool_report(pool_dir / "pool_build_report.json", MINIMAL_POOL_REPORT)

    def test_main_returns_0(self, tmp_path: Path) -> None:
        self._make_valid_inputs(tmp_path)
        rc = agg.main([
            "--scan-dir", str(tmp_path / "scan"),
            "--classify-dir", str(tmp_path / "classify"),
            "--pool-report", str(tmp_path / "pools" / "pool_build_report.json"),
            "--out-json", str(tmp_path / "report.json"),
            "--out-md", str(tmp_path / "report.md"),
        ])
        assert rc == 0
        assert (tmp_path / "report.json").exists()
        assert (tmp_path / "report.md").exists()

    def test_main_skips_if_outputs_exist(self, tmp_path: Path) -> None:
        """Second call without --force should not overwrite."""
        self._make_valid_inputs(tmp_path)
        out_json = tmp_path / "report.json"
        out_md = tmp_path / "report.md"
        out_json.write_text('{"existing": true}')
        out_md.write_text("# existing\n")

        rc = agg.main([
            "--scan-dir", str(tmp_path / "scan"),
            "--classify-dir", str(tmp_path / "classify"),
            "--pool-report", str(tmp_path / "pools" / "pool_build_report.json"),
            "--out-json", str(out_json),
            "--out-md", str(out_md),
        ])
        assert rc == 0
        assert json.loads(out_json.read_text()) == {"existing": True}

    def test_force_overwrites_existing_outputs(self, tmp_path: Path) -> None:
        self._make_valid_inputs(tmp_path)
        out_json = tmp_path / "report.json"
        out_md = tmp_path / "report.md"
        out_json.write_text('{"stale": true}')
        out_md.write_text("# stale\n")

        rc = agg.main([
            "--scan-dir", str(tmp_path / "scan"),
            "--classify-dir", str(tmp_path / "classify"),
            "--pool-report", str(tmp_path / "pools" / "pool_build_report.json"),
            "--out-json", str(out_json),
            "--out-md", str(out_md),
            "--force",
        ])
        assert rc == 0
        data = json.loads(out_json.read_text())
        assert "per_species" in data
