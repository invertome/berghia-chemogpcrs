"""Tests for scripts/classify_via_hmm.py.

Phase 4 Task 4.1 — HMM scan classifier. Runs hmmscan on candidate
sequences against the custom HMM library (Phase 2) plus Pfam fallback
HMMs, applies per-HMM E-value thresholds from the LOO benchmark, and
outputs per-candidate family/subfamily classifications.

Tests cover the local (no-subprocess) parts:
  - hmmscan tblout parsing (per-query hits sorted by E-value)
  - per-HMM threshold loading from loo_metrics.tsv
  - per-candidate best-hit assignment given thresholds + label map
  - HMM name -> (family, subfamily) parsing (reuse from Task 2.2)
  - unclassified path when best hit doesn't pass threshold
"""
from __future__ import annotations

from pathlib import Path

import classify_via_hmm as cvh


# ---- parse_hmmscan_tblout -----------------------------------------------

EXAMPLE_TBL = """\
#                                                               --- full sequence ---- --- best 1 domain ----
# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias
# ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----
aminergic            -          bste_gene_001        -           1.2e-150 510.5   0.0   2.1e-150 509.8   0.0
peptide              -          bste_gene_001        -           4.5e-50  175.3   0.0   8.2e-50  174.6   0.0
opsin                -          bste_gene_002        -           3.8e-200 670.9   0.0   3.8e-200 670.9   0.0
# end
"""


def test_parse_hmmscan_tblout_groups_by_query(tmp_path: Path) -> None:
    """Output: dict query_id -> list of (hmm_name, evalue, score) sorted
    ascending by E-value."""
    tbl = tmp_path / "scan.tbl"
    tbl.write_text(EXAMPLE_TBL)
    hits = cvh.parse_hmmscan_tblout(str(tbl))
    assert "bste_gene_001" in hits
    assert "bste_gene_002" in hits
    # bste_gene_001 has 2 hits; sorted by E-value (best first)
    assert hits["bste_gene_001"][0][0] == "aminergic"  # 1.2e-150 first
    assert hits["bste_gene_001"][0][1] == 1.2e-150
    assert hits["bste_gene_001"][1][0] == "peptide"
    # bste_gene_002 has 1 hit
    assert len(hits["bste_gene_002"]) == 1
    assert hits["bste_gene_002"][0][0] == "opsin"


def test_parse_hmmscan_tblout_skips_comments(tmp_path: Path) -> None:
    tbl = tmp_path / "empty.tbl"
    tbl.write_text("# only comments\n# nothing else\n")
    assert cvh.parse_hmmscan_tblout(str(tbl)) == {}


def test_parse_hmmscan_tblout_handles_missing_file(tmp_path: Path) -> None:
    """Gracefully handle missing tblout (e.g., hmmscan failed silently)."""
    assert cvh.parse_hmmscan_tblout(str(tmp_path / "missing.tbl")) == {}


# ---- load_thresholds ----------------------------------------------------

def test_load_thresholds_from_loo_metrics(tmp_path: Path) -> None:
    """Read per-family E-value thresholds from the LOO output."""
    f = tmp_path / "loo_metrics.tsv"
    f.write_text(
        "family\tn_total\ttp\tfn\tfp\trecall\tprecision\tevalue_threshold\n"
        "aminergic\t165\t165\t0\t0\t1.000\t1.000\t2.4e-105\n"
        "peptide\t189\t187\t2\t0\t0.989\t1.000\t2.7e-58\n"
    )
    thresholds = cvh.load_thresholds(str(f))
    assert thresholds["aminergic"] == 2.4e-105
    assert thresholds["peptide"] == 2.7e-58


def test_load_thresholds_handles_missing_file(tmp_path: Path) -> None:
    """If LOO metrics aren't available, return empty dict (caller decides
    fallback — typically a conservative default like 1e-10)."""
    assert cvh.load_thresholds(str(tmp_path / "missing.tsv")) == {}


# ---- assign_classification ----------------------------------------------

def test_assign_classification_picks_best_passing_hit() -> None:
    """Best hit wins, but it must pass that HMM's threshold."""
    # bste_gene_001 hits:
    #   aminergic_5HT 1e-200 (best, but threshold for aminergic is 1e-105)
    #   peptide       1e-50  (worse, threshold for peptide is 1e-58 — fails)
    hits = [
        ("aminergic_5HT", 1e-200, 700.0),
        ("peptide", 1e-50, 175.0),
    ]
    thresholds = {"aminergic": 1e-105, "peptide": 1e-58}
    result = cvh.assign_classification(hits, thresholds)
    assert result["family"] == "aminergic"
    assert result["subfamily"] == "5HT"
    assert result["evalue"] == 1e-200


def test_assign_classification_falls_through_when_best_fails() -> None:
    """If the best hit's E-value FAILS its threshold, look at the next-best
    that PASSES."""
    hits = [
        ("aminergic", 1e-50, 200.0),  # but threshold is 1e-105 → fails
        ("peptide", 1e-100, 350.0),    # threshold is 1e-58 → passes
    ]
    thresholds = {"aminergic": 1e-105, "peptide": 1e-58}
    result = cvh.assign_classification(hits, thresholds)
    assert result["family"] == "peptide"
    assert result["evalue"] == 1e-100


def test_assign_classification_returns_unclassified_when_all_fail() -> None:
    hits = [
        ("aminergic", 1e-50, 200.0),  # threshold 1e-105
        ("peptide", 1e-30, 100.0),     # threshold 1e-58
    ]
    thresholds = {"aminergic": 1e-105, "peptide": 1e-58}
    result = cvh.assign_classification(hits, thresholds)
    assert result["family"] == "unclassified-hmm"
    assert result["subfamily"] == ""


def test_assign_classification_returns_unclassified_when_no_hits() -> None:
    result = cvh.assign_classification([], {})
    assert result["family"] == "unclassified-hmm"


def test_assign_classification_default_threshold_when_family_missing() -> None:
    """If a family isn't in the thresholds dict (e.g. Pfam fallback HMM),
    use a conservative default (1e-10) so it can still pass."""
    hits = [
        ("class-A-7tm", 1e-30, 100.0),  # Pfam fallback, no LOO threshold
    ]
    thresholds = {"aminergic": 1e-105}  # no class-A-7tm key
    result = cvh.assign_classification(hits, thresholds)
    # Should pass with default 1e-10
    assert result["family"] == "class-A-7tm"


# ---- HMM name parsing (reused) ------------------------------------------

def test_label_for_hmm_uses_same_logic_as_validation() -> None:
    """Coarse: 'aminergic' -> ('aminergic', '')
       Medium: 'aminergic_5HT' -> ('aminergic', '5HT')
       Hyphen-rich: 'class-B-secretin' -> ('class-B-secretin', '')
       Compound subfamily: 'peptide_NPY-NPF' -> ('peptide', 'NPY-NPF')"""
    assert cvh.label_for_hmm("aminergic") == ("aminergic", "")
    assert cvh.label_for_hmm("aminergic_5HT") == ("aminergic", "5HT")
    assert cvh.label_for_hmm("class-B-secretin") == ("class-B-secretin", "")
    assert cvh.label_for_hmm("peptide_NPY-NPF") == ("peptide", "NPY-NPF")
