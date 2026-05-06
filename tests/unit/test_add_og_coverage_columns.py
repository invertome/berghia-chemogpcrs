"""Tests for scripts/add_og_coverage_columns.py.

The transparency-column script adds per-orthogroup reference-CDS coverage
information to the ranked candidates CSV so reviewers can see which
candidates' dN/dS estimates were computed from sparse reference data.

Adds three columns:
    og_n_ref_cds          — number of OG members that have a CDS in
                            all_references_cds.fna
    og_n_total            — total number of OG members
    og_dnds_reliability   — 'high' (>=10 ref CDS) / 'medium' (>=5) / 'low'

This does NOT change ranking. It only annotates each row with the
data-quality context for its dN/dS contribution.

Bug context: the May 2026 review uncovered that bad-CDS rates are non-
uniform across species (cate 81%, aplcal 68%, arvu 58% bad). The OGs
where dN/dS would be most informative (chemoreceptor LSE expansions)
have the sparsest reference CDS coverage. Without these transparency
columns, reviewers can't tell if a top-ranked candidate's dN/dS score
came from a well-sampled OG or a degenerate one.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

# conftest.py adds scripts/ to sys.path
import add_og_coverage_columns as aoc


# ---- reliability_flag ---------------------------------------------------

@pytest.mark.parametrize("n,expected", [
    (0,  "low"),
    (1,  "low"),
    (4,  "low"),
    (5,  "medium"),
    (9,  "medium"),
    (10, "high"),
    (50, "high"),
])
def test_reliability_flag_thresholds(n: int, expected: str) -> None:
    """Reliability bins: <5 = low, 5-9 = medium, >=10 = high."""
    assert aoc.reliability_flag(n) == expected


# ---- parse_cds_ids ------------------------------------------------------

def test_parse_cds_ids_extracts_first_token(tmp_path: Path) -> None:
    """CDS FASTA headers: ID is the first whitespace token after '>'."""
    fa = tmp_path / "cds.fna"
    fa.write_text(
        ">cate_KB300474.1_27 some metadata\nATGAAA\n"
        ">aplcal_NP_001191513.1\nATGAAACTG\n"
        ">phaust_GCAxx_42 [gene=foo]\nATGGGG\n"
    )
    ids = aoc.parse_cds_ids(str(fa))
    assert ids == {
        "cate_KB300474.1_27",
        "aplcal_NP_001191513.1",
        "phaust_GCAxx_42",
    }


def test_parse_cds_ids_handles_missing_file(tmp_path: Path) -> None:
    """If the CDS FASTA doesn't exist, return an empty set (graceful)."""
    nonexistent = tmp_path / "no_cds.fna"
    assert aoc.parse_cds_ids(str(nonexistent)) == set()


# ---- load_og_members ----------------------------------------------------

def _write_orthogroups_tsv(path: Path, content: str) -> None:
    path.write_text(content)


def test_load_og_members_orthofinder_format(tmp_path: Path) -> None:
    """OrthoFinder Orthogroups.tsv: tab-separated; column 0 is OG id; remaining
    columns are comma-separated gene lists per species (one column = one species)."""
    tsv = tmp_path / "Orthogroups.tsv"
    _write_orthogroups_tsv(tsv,
        "Orthogroup\tspeciesA\tspeciesB\tspeciesC\n"
        "OG0000001\tgene1, gene2\tgene3\tgene4, gene5\n"
        "OG0000002\tgene6\t\tgene7\n"
    )
    og = aoc.load_og_members(str(tsv))
    assert og["OG0000001"] == ["gene1", "gene2", "gene3", "gene4", "gene5"]
    assert og["OG0000002"] == ["gene6", "gene7"]


def test_load_og_members_handles_empty_cells(tmp_path: Path) -> None:
    """Empty cells (species with no member in OG) should not contribute
    spurious empty-string members."""
    tsv = tmp_path / "Orthogroups.tsv"
    _write_orthogroups_tsv(tsv,
        "Orthogroup\tspA\tspB\n"
        "OG0000001\t\tgene1\n"
    )
    og = aoc.load_og_members(str(tsv))
    assert og["OG0000001"] == ["gene1"]


# ---- end-to-end on a small DataFrame -----------------------------------

def test_add_columns_end_to_end(tmp_path: Path) -> None:
    """Smoke test of the full pipeline: ranked CSV in -> ranked CSV out
    with og_n_ref_cds + og_n_total + og_dnds_reliability columns added."""
    cds_fa = tmp_path / "all_references_cds.fna"
    cds_fa.write_text(
        ">aplcal_NP_001191513.1\nATGAAA\n"
        ">aplcal_NP_001191494.1\nATGAAA\n"
        ">aplcal_NP_001191495.1\nATGAAA\n"
        ">cate_KB300474.1_27\nATGAAA\n"
        ">phaust_GCAxx_42\nATGAAA\n"
    )
    og_tsv = tmp_path / "Orthogroups.tsv"
    og_tsv.write_text(
        "Orthogroup\tA\tB\tC\n"
        # OG with all 5 reference CDS present + 2 berghia (high coverage)
        "OG0000001\taplcal_NP_001191513.1, aplcal_NP_001191494.1, aplcal_NP_001191495.1, "
        "cate_KB300474.1_27, phaust_GCAxx_42\tbste_gene_a, bste_gene_b\t\n"
        # OG with 0 reference CDS hits (low coverage; only Berghia genes)
        "OG0000002\tbste_gene_c\t\t\n"
    )
    ranked_csv = tmp_path / "ranked.csv"
    ranked_csv.write_text(
        "id,orthogroup,rank_score\n"
        "bste_gene_a,OG0000001,0.95\n"
        "bste_gene_b,OG0000001,0.92\n"
        "bste_gene_c,OG0000002,0.80\n"
    )
    out_csv = tmp_path / "ranked_with_coverage.csv"

    aoc.add_coverage_columns(
        ranked_csv_path=str(ranked_csv),
        cds_fasta_path=str(cds_fa),
        orthogroups_tsv_path=str(og_tsv),
        out_path=str(out_csv),
    )

    df = pd.read_csv(out_csv)
    assert list(df.columns)[-3:] == ["og_n_ref_cds", "og_n_total", "og_dnds_reliability"]
    # OG0000001: 5 ref CDS hits + 2 berghia = 7 total
    row1 = df[df["id"] == "bste_gene_a"].iloc[0]
    assert row1["og_n_ref_cds"] == 5
    assert row1["og_n_total"] == 7
    assert row1["og_dnds_reliability"] == "medium"  # 5-9 = medium
    # OG0000002: 0 ref CDS, 1 berghia
    row3 = df[df["id"] == "bste_gene_c"].iloc[0]
    assert row3["og_n_ref_cds"] == 0
    assert row3["og_n_total"] == 1
    assert row3["og_dnds_reliability"] == "low"


def test_missing_orthogroup_for_candidate(tmp_path: Path) -> None:
    """If a candidate's OG isn't in the Orthogroups.tsv (orphan gene),
    coverage = 0, total = 0, reliability = low. Don't crash."""
    cds_fa = tmp_path / "cds.fna"
    cds_fa.write_text(">x\nATG\n")
    og_tsv = tmp_path / "og.tsv"
    og_tsv.write_text("Orthogroup\tA\nOG_other\tx\n")
    ranked_csv = tmp_path / "ranked.csv"
    ranked_csv.write_text("id,orthogroup,rank_score\norphan,OG_missing,0.5\n")
    out_csv = tmp_path / "out.csv"

    aoc.add_coverage_columns(
        ranked_csv_path=str(ranked_csv),
        cds_fasta_path=str(cds_fa),
        orthogroups_tsv_path=str(og_tsv),
        out_path=str(out_csv),
    )

    df = pd.read_csv(out_csv)
    assert df.iloc[0]["og_n_ref_cds"] == 0
    assert df.iloc[0]["og_n_total"] == 0
    assert df.iloc[0]["og_dnds_reliability"] == "low"


def test_ranked_csv_without_orthogroup_column(tmp_path: Path) -> None:
    """If the ranked CSV doesn't have an 'orthogroup' column, the script
    must still produce output but with NaN coverage and 'low' reliability
    (and a stderr warning)."""
    cds_fa = tmp_path / "cds.fna"; cds_fa.write_text(">x\nATG\n")
    og_tsv = tmp_path / "og.tsv"; og_tsv.write_text("Orthogroup\tA\nOG1\tx\n")
    ranked_csv = tmp_path / "ranked.csv"
    ranked_csv.write_text("id,rank_score\ngene_q,0.5\n")  # no orthogroup col
    out_csv = tmp_path / "out.csv"

    aoc.add_coverage_columns(
        ranked_csv_path=str(ranked_csv),
        cds_fasta_path=str(cds_fa),
        orthogroups_tsv_path=str(og_tsv),
        out_path=str(out_csv),
    )

    df = pd.read_csv(out_csv)
    assert df.iloc[0]["og_n_ref_cds"] == 0
    assert df.iloc[0]["og_dnds_reliability"] == "low"
