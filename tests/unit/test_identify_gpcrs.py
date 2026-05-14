"""Tests for scripts/identify_gpcrs.py.

Bead -m1f restructure follow-up: stage 02 now does HMM-first GPCR
identification (against lse.hmm + classification HMMs + Pfam fallback)
BEFORE running TMbed, so TMbed only sees the ~500-3000 GPCR-positive
proteins instead of the full 86k-protein transcriptome.

identify_gpcrs.py owns the merge: given the classification TSV
(produced by classify_via_hmm.py) and an hmmsearch --tblout against
lse.hmm, it unions the GPCR-positive IDs and writes:
  - a flat ID list for downstream extract+TMbed
  - a census TSV with family/subfamily/source for the follow-up paper
"""
from __future__ import annotations

from pathlib import Path

import identify_gpcrs as ig


# ---- merge_gpcr_evidence -------------------------------------------------

CLASS_TSV = """\
candidate_id\thmm_family\thmm_subfamily\tevalue\tscore\tevidence
bste_gene_001\taminergic\t5HT\t1.2e-150\t510.5\thmm
bste_gene_002\topsin\t\t3.8e-200\t670.9\thmm
bste_gene_010\tunclassified-hmm\t\t\t\thmm
bste_gene_020\t\t\t\t\thmm
"""

LSE_TBLOUT = """\
#                                                               --- full sequence ---- --- best 1 domain ----
# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias
# ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----
OG_lse_00042       -          bste_gene_100        -           1.5e-80   272.4   0.0   1.5e-80   272.4   0.0
OG_lse_00042       -          bste_gene_001        -           2.3e-12    45.7   0.0   2.3e-12    45.7   0.0
OG_lse_00099       -          bste_gene_200        -           4.7e-15    52.1   0.0   4.7e-15    52.1   0.0
# end
"""

# Nath et al. one_to_one_ortholog GPCR-OG HMMs. NOT functionally confirmed
# chemoreceptors — these are broad metazoan GPCR orthologs. They give
# additional GPCR-detection coverage on top of classification + lse, but
# carry no chemoreceptor-specific labeling.
NATH_ORTHOLOG_TBLOUT = """\
#                                                               --- full sequence ---- --- best 1 domain ----
# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias
# ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----
OG_nath_00007      -          bste_gene_300        -           5.0e-90   305.0   0.0   5.0e-90   305.0   0.0
OG_nath_00007      -          bste_gene_100        -           1.0e-30   100.2   0.0   1.0e-30   100.2   0.0
OG_nath_00050      -          bste_gene_400        -           2.0e-50   180.5   0.0   2.0e-50   180.5   0.0
# end
"""


def test_merge_classification_only(tmp_path: Path) -> None:
    """When only classification TSV present, returns those classified IDs.
    'unclassified-hmm' and empty-family entries are dropped."""
    cls = tmp_path / "cls.tsv"
    cls.write_text(CLASS_TSV)
    merged = ig.merge_gpcr_evidence(str(cls), None)
    # bste_gene_001 and 002 are classified; 010 and 020 are not GPCR-positive
    assert set(merged.keys()) == {"bste_gene_001", "bste_gene_002"}
    assert merged["bste_gene_001"]["family"] == "aminergic"
    assert merged["bste_gene_001"]["subfamily"] == "5HT"
    assert merged["bste_gene_001"]["source"] == "classification"


def test_merge_lse_only(tmp_path: Path) -> None:
    """When only LSE tblout present, returns those IDs with family='lse'
    (lineage-specific expansion — NOT functionally confirmed chemoreceptor)
    and the best-hit OG as subfamily."""
    lse = tmp_path / "lse.tbl"
    lse.write_text(LSE_TBLOUT)
    merged = ig.merge_gpcr_evidence(None, str(lse))
    # 3 unique queries: 100, 001, 200
    assert set(merged.keys()) == {"bste_gene_100", "bste_gene_001", "bste_gene_200"}
    assert merged["bste_gene_100"]["family"] == "lse"
    assert merged["bste_gene_100"]["subfamily"] == "OG_lse_00042"
    assert merged["bste_gene_100"]["source"] == "lse"


def test_merge_nath_ortholog_only(tmp_path: Path) -> None:
    """When only nath-ortholog tblout present, returns those IDs with
    family='nath_ortholog' (broad metazoan GPCR — NOT chemoreceptor-confirmed)
    and the best-hit OG as subfamily."""
    nath = tmp_path / "nath.tbl"
    nath.write_text(NATH_ORTHOLOG_TBLOUT)
    merged = ig.merge_gpcr_evidence(None, None, nath_ortholog_tblout=str(nath))
    # 3 unique queries: 300, 100, 400
    assert set(merged.keys()) == {"bste_gene_300", "bste_gene_100", "bste_gene_400"}
    assert merged["bste_gene_300"]["family"] == "nath_ortholog"
    assert merged["bste_gene_300"]["subfamily"] == "OG_nath_00007"
    assert merged["bste_gene_300"]["source"] == "nath_ortholog"


def test_merge_classification_wins_when_both(tmp_path: Path) -> None:
    """When a protein is in BOTH classification and lse tblout,
    classification takes precedence (it's the more specific assignment).
    Source reflects the dual evidence."""
    cls = tmp_path / "cls.tsv"
    cls.write_text(CLASS_TSV)
    lse = tmp_path / "lse.tbl"
    lse.write_text(LSE_TBLOUT)
    merged = ig.merge_gpcr_evidence(str(cls), str(lse))
    # bste_gene_001 is in BOTH — classification wins
    assert merged["bste_gene_001"]["family"] == "aminergic"
    assert merged["bste_gene_001"]["subfamily"] == "5HT"
    assert merged["bste_gene_001"]["source"] == "classification+lse"
    # bste_gene_100 / 200 are only in lse
    assert merged["bste_gene_100"]["source"] == "lse"
    # bste_gene_002 is only in classification
    assert merged["bste_gene_002"]["source"] == "classification"
    # Total unique GPCR-positive: 001, 002, 100, 200 (NOT 010, 020 — unclassified)
    assert set(merged.keys()) == {"bste_gene_001", "bste_gene_002",
                                  "bste_gene_100", "bste_gene_200"}


def test_merge_three_sources_precedence(tmp_path: Path) -> None:
    """Triple-evidence merge: classification > {lse, nath_ortholog}.
    Source tag should reflect every source that hit the protein. The
    Nath ortholog set is broad metazoan GPCRs, not chemoreceptors —
    its family label is 'nath_ortholog', not 'conserved_chemoreceptor'."""
    cls = tmp_path / "cls.tsv"
    cls.write_text(CLASS_TSV)
    lse = tmp_path / "lse.tbl"
    lse.write_text(LSE_TBLOUT)
    nath = tmp_path / "nath.tbl"
    nath.write_text(NATH_ORTHOLOG_TBLOUT)
    merged = ig.merge_gpcr_evidence(str(cls), str(lse), nath_ortholog_tblout=str(nath))
    # bste_gene_100 hits lse + nath_ortholog (not classification) -> lse family,
    # source 'lse+nath_ortholog'
    assert merged["bste_gene_100"]["family"] == "lse"
    assert merged["bste_gene_100"]["source"] == "lse+nath_ortholog"
    # bste_gene_300 only in nath_ortholog
    assert merged["bste_gene_300"]["family"] == "nath_ortholog"
    assert merged["bste_gene_300"]["source"] == "nath_ortholog"
    # bste_gene_001 in classification + lse (not nath_ortholog) -> classification wins
    assert merged["bste_gene_001"]["family"] == "aminergic"
    assert merged["bste_gene_001"]["source"] == "classification+lse"
    # bste_gene_400 only in nath_ortholog
    assert merged["bste_gene_400"]["family"] == "nath_ortholog"


def test_merge_neither_returns_empty(tmp_path: Path) -> None:
    """All inputs missing/None -> empty dict."""
    assert ig.merge_gpcr_evidence(None, None) == {}
    assert ig.merge_gpcr_evidence(None, None, nath_ortholog_tblout=None) == {}


# ---- write_ids -----------------------------------------------------------

def test_write_ids_sorted_unique(tmp_path: Path) -> None:
    merged = {
        "bste_gene_002": {"family": "opsin", "subfamily": "", "source": "classification"},
        "bste_gene_001": {"family": "aminergic", "subfamily": "5HT", "source": "classification"},
        "bste_gene_010": {"family": "lse_chemoreceptor", "subfamily": "OG_lse_007", "source": "lse"},
    }
    out = tmp_path / "ids.txt"
    ig.write_ids(merged, str(out))
    content = out.read_text().splitlines()
    assert content == ["bste_gene_001", "bste_gene_002", "bste_gene_010"]


# ---- write_census --------------------------------------------------------

def test_write_census_has_header_and_rows(tmp_path: Path) -> None:
    merged = {
        "bste_gene_001": {
            "family": "aminergic", "subfamily": "5HT",
            "evalue": 1.2e-150, "source": "classification",
        },
        "bste_gene_100": {
            "family": "lse_chemoreceptor", "subfamily": "OG_lse_00042",
            "evalue": 1.5e-80, "source": "lse",
        },
    }
    out = tmp_path / "census.tsv"
    ig.write_census(merged, str(out))
    lines = out.read_text().splitlines()
    assert lines[0].startswith("seq_id\t")
    assert "family" in lines[0]
    assert "subfamily" in lines[0]
    assert "evalue" in lines[0]
    assert "source" in lines[0]
    # Two data rows, sorted by seq_id
    assert lines[1].startswith("bste_gene_001\t")
    assert "aminergic" in lines[1]
    assert "5HT" in lines[1]
    assert lines[2].startswith("bste_gene_100\t")
    assert "lse_chemoreceptor" in lines[2]
