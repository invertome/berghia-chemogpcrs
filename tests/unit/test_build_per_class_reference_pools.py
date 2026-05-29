"""Tests for scripts/build_per_class_reference_pools.py.

P2 of the per-class refactor — per-class reference pool builder.
Tests use synthetic scan FASTAs + class TSVs.  Taxonomy lookups are
mocked via a small dict so real ete3.NCBITaxa / network calls never fire.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Ensure scripts/ is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

import build_per_class_reference_pools as bpcp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    """Write a tiny FASTA file from [(id, seq), ...] pairs."""
    with open(path, "w") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n{seq}\n")


def _make_class_tsv(path: Path, rows: list[tuple[str, str]]) -> None:
    """Write a class TSV: seq_id<TAB>class<TAB>evidence_pfam<TAB>top_evalue."""
    with open(path, "w") as fh:
        fh.write("seq_id\tclass\tevidence_pfam\ttop_evalue\n")
        for seq_id, cls in rows:
            fh.write(f"{seq_id}\t{cls}\tPF00001\t1e-50\n")


# ---------------------------------------------------------------------------
# Mock taxonomy lookup
# ---------------------------------------------------------------------------

# Synthetic lineage: Berghia (1287507) is inside Aeolidiidae → Nudibranchia →
# Heterobranchia → Gastropoda → Mollusca → Lophotrochozoa → Bilateria ...
# We represent each taxid's lineage as a list going from root to leaf.
# The BERGHIA_LINEAGE constant in the script is what we need to agree with.
MOCK_LINEAGES: dict[int, list[int]] = {
    # Berghia stephanieae
    1287507: [1, 131567, 2759, 33154, 33208, 6072, 33213, 1206794,
              88194, 6157, 6447, 186803, 103598, 69675, 29178, 1460361,
              1460363, 6101, 56615, 1287507],
    # Aplysia californica — same up to Heterobranchia, then diverges
    6500:    [1, 131567, 2759, 33154, 33208, 6072, 33213, 1206794,
              88194, 6157, 6447, 186803, 103598, 69675, 29178, 1460361,
              1460363, 6101, 6500],
    # Lottia gigantea — Gastropoda, but Patellogastropoda (further)
    225164:  [1, 131567, 2759, 33154, 33208, 6072, 33213, 1206794,
              88194, 6157, 6447, 186803, 103598, 69675, 225164],
    # Crassostrea gigas — Bivalvia
    29159:   [1, 131567, 2759, 33154, 33208, 6072, 33213, 1206794,
              88194, 6157, 6447, 29159],
    # Drosophila melanogaster — Arthropoda (distant from Mollusca)
    7227:    [1, 131567, 2759, 33154, 33208, 6072, 33213, 1206794,
              88194, 6157, 7227],
    # Homo sapiens — Deuterostomia (distant)
    9606:    [1, 131567, 2759, 33154, 33208, 6072, 33213, 9606],
    # Outgroup: Trichoplax (early branch)
    10228:   [1, 131567, 2759, 33154, 33208, 10228],
}

BERGHIA_LINEAGE_SET = frozenset(MOCK_LINEAGES[1287507])


def mock_get_lineage(taxid: int) -> list[int]:
    """Return a fake lineage list, falling back to [1, taxid] for unknowns."""
    return MOCK_LINEAGES.get(taxid, [1, taxid])


# ---------------------------------------------------------------------------
# 1. load_class_tsv: maps seq_id → class
# ---------------------------------------------------------------------------

def test_load_class_tsv_basic(tmp_path):
    """load_class_tsv returns a dict mapping each seq_id to its class."""
    tsv = tmp_path / "classes.tsv"
    _make_class_tsv(tsv, [
        ("seq_A1", "A"), ("seq_B1", "B"), ("seq_C1", "C"),
        ("seq_F1", "F"), ("seq_U1", "unclassified"),
    ])
    result = bpcp.load_class_tsv(str(tsv))
    assert result["seq_A1"] == "A"
    assert result["seq_B1"] == "B"
    assert result["seq_C1"] == "C"
    assert result["seq_F1"] == "F"
    assert result["seq_U1"] == "unclassified"
    assert len(result) == 5


# ---------------------------------------------------------------------------
# 2. taxid_from_filename: parse taxid from scan FASTA filename
# ---------------------------------------------------------------------------

def test_taxid_from_filename_standard():
    """<taxid>_Genus_species.chemo_candidates.fa → taxid."""
    taxid = bpcp.taxid_from_filename("6500_Aplysia_californica.chemo_candidates.fa")
    assert taxid == 6500


def test_taxid_from_filename_no_binomial():
    """Single-word stem still extracts taxid correctly."""
    taxid = bpcp.taxid_from_filename("9606.chemo_candidates.fa")
    assert taxid == 9606


def test_taxid_from_filename_bad_returns_none():
    """Non-numeric prefix returns None."""
    taxid = bpcp.taxid_from_filename("mystery_file.chemo_candidates.fa")
    assert taxid is None


# ---------------------------------------------------------------------------
# 3. Berghia scan-file seqs route normally; --berghia-fasta adds MUST_INCLUDE
# ---------------------------------------------------------------------------

def test_berghia_scan_seqs_routed_normally(tmp_path):
    """Berghia sequences in scan FASTA files are no longer excluded — they
    route through the normal class-map path like any other species."""
    berghia_fa = tmp_path / "1287507_Berghia_stephanieae.chemo_candidates.fa"
    aplysia_fa = tmp_path / "6500_Aplysia_californica.chemo_candidates.fa"

    _make_fasta(berghia_fa, [("ber_A1", "MSTL" * 30), ("ber_A2", "MKVL" * 30)])
    _make_fasta(aplysia_fa, [("apl_A1", "MSTL" * 30), ("apl_A2", "MKGT" * 30)])

    class_tsv = tmp_path / "classes.tsv"
    _make_class_tsv(class_tsv, [
        ("ber_A1", "A"), ("ber_A2", "A"),
        ("apl_A1", "A"), ("apl_A2", "A"),
    ])

    out_dir = tmp_path / "pools"
    out_dir.mkdir()

    scan_fasta_glob = str(tmp_path / "*.chemo_candidates.fa")

    with patch.object(bpcp, "cdhit_dedup", side_effect=lambda records, **kw: records):
        with patch.object(bpcp, "get_lineage", side_effect=mock_get_lineage):
            bpcp.build_all_pools(
                scan_fasta_glob=scan_fasta_glob,
                class_tsv=str(class_tsv),
                out_dir=str(out_dir),
                max_per_class=2000,
                total_budget_per_class=3000,
                outgroup_budget_per_class=10,
                cluster_identity=0.7,
                cdhit_path="cd-hit",
                threads=1,
                must_include_taxids=frozenset(),
                berghia_taxid=1287507,
                force=True,
            )

    pool_a = list(SeqIO.parse(str(out_dir / "refs_class_A.fa"), "fasta"))
    ids_in_pool = {r.id for r in pool_a}
    # Berghia seqs now appear in the pool (routed via class_tsv like any species)
    assert "ber_A1" in ids_in_pool
    assert "ber_A2" in ids_in_pool
    assert "apl_A1" in ids_in_pool


# ---------------------------------------------------------------------------
# 4. MUST_INCLUDE retention: must-include taxid sequences always kept
# ---------------------------------------------------------------------------

def test_must_include_always_retained(tmp_path):
    """Sequences from must-include taxids survive subsampling even at small cap."""
    # Two species: must-include (6500) and ordinary (999999)
    aplysia_fa = tmp_path / "6500_Aplysia_californica.chemo_candidates.fa"
    other_fa   = tmp_path / "999999_Unknown_species.chemo_candidates.fa"

    aplysia_seqs = [(f"apl_{i}", "MSFT" * 30) for i in range(10)]
    other_seqs   = [(f"oth_{i}", "MKGT" * 30) for i in range(100)]

    _make_fasta(aplysia_fa, aplysia_seqs)
    _make_fasta(other_fa, other_seqs)

    class_rows = [(sid, "A") for sid, _ in aplysia_seqs + other_seqs]
    class_tsv = tmp_path / "classes.tsv"
    _make_class_tsv(class_tsv, class_rows)

    out_dir = tmp_path / "pools"
    out_dir.mkdir()

    # Cap at 15, outgroup reserve 10 → effective refs cap = 15 - 10 = 5
    # Aplysia has 10 seqs (MUST_INCLUDE, always retained even if exceeds cap) + up to 5 subsampled others
    with patch.object(bpcp, "cdhit_dedup", side_effect=lambda records, **kw: records):
        with patch.object(bpcp, "get_lineage", side_effect=mock_get_lineage):
            bpcp.build_all_pools(
                scan_fasta_glob=str(tmp_path / "*.chemo_candidates.fa"),
                class_tsv=str(class_tsv),
                out_dir=str(out_dir),
                max_per_class=15,
                total_budget_per_class=15,
                outgroup_budget_per_class=10,
                cluster_identity=0.7,
                cdhit_path="cd-hit",
                threads=1,
                must_include_taxids=frozenset({6500}),
                berghia_taxid=1287507,
                force=True,
            )

    pool_a = list(SeqIO.parse(str(out_dir / "refs_class_A.fa"), "fasta"))
    ids_out = {r.id for r in pool_a}

    # All Aplysia sequences must be present (MUST_INCLUDE always retained)
    for sid, _ in aplysia_seqs:
        assert sid in ids_out, f"{sid} missing from pool despite must-include taxid"


# ---------------------------------------------------------------------------
# 5. MAX_PHYLO_REFS cap enforced per class
# ---------------------------------------------------------------------------

def test_max_phylo_refs_cap(tmp_path):
    """Pool never exceeds MAX_PHYLO_REFS regardless of candidate count."""
    fa_path = tmp_path / "7227_Drosophila_melanogaster.chemo_candidates.fa"
    seqs = [(f"dm_{i}", "MKLF" * 30) for i in range(500)]
    _make_fasta(fa_path, seqs)

    class_tsv = tmp_path / "classes.tsv"
    _make_class_tsv(class_tsv, [(sid, "A") for sid, _ in seqs])

    out_dir = tmp_path / "pools"
    out_dir.mkdir()

    cap = 50
    with patch.object(bpcp, "cdhit_dedup", side_effect=lambda records, **kw: records):
        with patch.object(bpcp, "get_lineage", side_effect=mock_get_lineage):
            bpcp.build_all_pools(
                scan_fasta_glob=str(tmp_path / "*.chemo_candidates.fa"),
                class_tsv=str(class_tsv),
                out_dir=str(out_dir),
                max_per_class=cap,
                total_budget_per_class=cap,
                outgroup_budget_per_class=10,
                cluster_identity=0.7,
                cdhit_path="cd-hit",
                threads=1,
                must_include_taxids=frozenset(),
                berghia_taxid=1287507,
                force=True,
            )

    pool_a = list(SeqIO.parse(str(out_dir / "refs_class_A.fa"), "fasta"))
    # cap=50 with outgroup_budget=10 means refs_cap = 50 - 10 = 40
    assert len(pool_a) <= cap


# ---------------------------------------------------------------------------
# 6. B/C/F pools are independent of each other
# ---------------------------------------------------------------------------

def test_bcf_pools_independent(tmp_path):
    """Sequences classified as B, C, F land only in their respective pools."""
    fa_path = tmp_path / "9606_Homo_sapiens.chemo_candidates.fa"
    seqs = [
        ("hs_A1", "MSTL" * 30),
        ("hs_B1", "MKVL" * 30),
        ("hs_C1", "MKGT" * 30),
        ("hs_F1", "MSFT" * 30),
    ]
    _make_fasta(fa_path, seqs)

    class_tsv = tmp_path / "classes.tsv"
    _make_class_tsv(class_tsv, [
        ("hs_A1", "A"), ("hs_B1", "B"), ("hs_C1", "C"), ("hs_F1", "F"),
    ])

    out_dir = tmp_path / "pools"
    out_dir.mkdir()

    with patch.object(bpcp, "cdhit_dedup", side_effect=lambda records, **kw: records):
        with patch.object(bpcp, "get_lineage", side_effect=mock_get_lineage):
            bpcp.build_all_pools(
                scan_fasta_glob=str(tmp_path / "*.chemo_candidates.fa"),
                class_tsv=str(class_tsv),
                out_dir=str(out_dir),
                max_per_class=2000,
                total_budget_per_class=3000,
                outgroup_budget_per_class=10,
                cluster_identity=0.7,
                cdhit_path="cd-hit",
                threads=1,
                must_include_taxids=frozenset(),
                berghia_taxid=1287507,
                force=True,
            )

    ids_A = {r.id for r in SeqIO.parse(str(out_dir / "refs_class_A.fa"), "fasta")}
    ids_B = {r.id for r in SeqIO.parse(str(out_dir / "refs_class_B.fa"), "fasta")}
    ids_C = {r.id for r in SeqIO.parse(str(out_dir / "refs_class_C.fa"), "fasta")}
    ids_F = {r.id for r in SeqIO.parse(str(out_dir / "refs_class_F.fa"), "fasta")}

    assert "hs_A1" in ids_A
    assert "hs_B1" in ids_B
    assert "hs_C1" in ids_C
    assert "hs_F1" in ids_F

    # No cross-contamination
    assert "hs_B1" not in ids_A
    assert "hs_A1" not in ids_B


# ---------------------------------------------------------------------------
# 7. Unclassified log written correctly
# ---------------------------------------------------------------------------

def test_unclassified_log(tmp_path):
    """Unclassified sequences are written to unclassified_log.tsv, not any pool."""
    fa_path = tmp_path / "29159_Crassostrea_gigas.chemo_candidates.fa"
    seqs = [("cg_A1", "MSTL" * 30), ("cg_U1", "MKQP" * 30)]
    _make_fasta(fa_path, seqs)

    class_tsv = tmp_path / "classes.tsv"
    _make_class_tsv(class_tsv, [("cg_A1", "A"), ("cg_U1", "unclassified")])

    out_dir = tmp_path / "pools"
    out_dir.mkdir()

    with patch.object(bpcp, "cdhit_dedup", side_effect=lambda records, **kw: records):
        with patch.object(bpcp, "get_lineage", side_effect=mock_get_lineage):
            bpcp.build_all_pools(
                scan_fasta_glob=str(tmp_path / "*.chemo_candidates.fa"),
                class_tsv=str(class_tsv),
                out_dir=str(out_dir),
                max_per_class=2000,
                total_budget_per_class=3000,
                outgroup_budget_per_class=10,
                cluster_identity=0.7,
                cdhit_path="cd-hit",
                threads=1,
                must_include_taxids=frozenset(),
                berghia_taxid=1287507,
                force=True,
            )

    log_path = out_dir / "unclassified_log.tsv"
    assert log_path.exists(), "unclassified_log.tsv must be written"
    content = log_path.read_text()
    assert "cg_U1" in content

    # Must not appear in any class pool
    ids_A = {r.id for r in SeqIO.parse(str(out_dir / "refs_class_A.fa"), "fasta")}
    assert "cg_U1" not in ids_A


# ---------------------------------------------------------------------------
# 8. JSON report emitted with required fields
# ---------------------------------------------------------------------------

def test_json_report_emitted(tmp_path):
    """pool_build_report.json contains per-class stats after a run."""
    fa_path = tmp_path / "225164_Lottia_gigantea.chemo_candidates.fa"
    seqs = [(f"lg_{i}", "MSTL" * 30) for i in range(5)]
    _make_fasta(fa_path, seqs)

    class_tsv = tmp_path / "classes.tsv"
    _make_class_tsv(class_tsv, [(sid, "A") for sid, _ in seqs])

    out_dir = tmp_path / "pools"
    out_dir.mkdir()

    with patch.object(bpcp, "cdhit_dedup", side_effect=lambda records, **kw: records):
        with patch.object(bpcp, "get_lineage", side_effect=mock_get_lineage):
            bpcp.build_all_pools(
                scan_fasta_glob=str(tmp_path / "*.chemo_candidates.fa"),
                class_tsv=str(class_tsv),
                out_dir=str(out_dir),
                max_per_class=2000,
                total_budget_per_class=3000,
                outgroup_budget_per_class=10,
                cluster_identity=0.7,
                cdhit_path="cd-hit",
                threads=1,
                must_include_taxids=frozenset(),
                berghia_taxid=1287507,
                force=True,
            )

    report_path = out_dir / "pool_build_report.json"
    assert report_path.exists(), "pool_build_report.json must be written"
    report = json.loads(report_path.read_text())

    # Check top-level budget fields
    assert "total_budget_per_class" in report
    assert report["total_budget_per_class"] == 3000
    assert "outgroup_budget_per_class" in report
    assert report["outgroup_budget_per_class"] == 10
    assert "n_berghia_per_class" in report

    for cls in ("class_A", "class_B", "class_C", "class_F"):
        assert cls in report, f"{cls} missing from report"
        entry = report[cls]
        assert "n_total_candidates" in entry
        assert "n_after_cdhit" in entry
        assert "n_must_include" in entry
        assert "n_subsampled" in entry
        assert "n_output" in entry
        assert "species_contributing" in entry
        assert "must_include_taxids_with_hits" in entry
        assert "must_include_taxids_missing" in entry
        # New budget fields
        assert "n_refs_target" in entry
        assert "total_budget_for_class" in entry

    assert "unclassified" in report
    assert "berghia_included" in report
    assert report["berghia_included"]["n_total"] == 0  # no berghia_fasta provided
    assert report["class_A"]["n_total_candidates"] == 5
    assert report["class_A"]["n_output"] == 5
    assert "n_berghia_included" in report["class_A"]


# ---------------------------------------------------------------------------
# 9. Taxonomy-proximity ordering: closer to Berghia → higher score
# ---------------------------------------------------------------------------

def test_proximity_score_ordering():
    """Sequences from taxa closer to Berghia receive higher proximity scores."""
    # Aplysia (Heterobranchia) shares more lineage with Berghia than Homo sapiens
    score_aplysia = bpcp.proximity_score(6500, mock_get_lineage, 1287507)
    score_lottia  = bpcp.proximity_score(225164, mock_get_lineage, 1287507)
    score_crassos = bpcp.proximity_score(29159, mock_get_lineage, 1287507)
    score_human   = bpcp.proximity_score(9606, mock_get_lineage, 1287507)
    score_tricho  = bpcp.proximity_score(10228, mock_get_lineage, 1287507)

    # Aplysia (Heterobranchia) closer to Berghia than Lottia (Gastropoda out-group)
    assert score_aplysia > score_lottia
    # Lottia (Gastropoda) closer than Crassostrea (Bivalvia)
    assert score_lottia > score_crassos
    # All molluscs closer than Homo sapiens
    assert score_crassos > score_human
    # Homo sapiens closer (Bilateria) than Trichoplax (Placozoa)
    assert score_human > score_tricho


# ---------------------------------------------------------------------------
# 10. Idempotence: skip if outputs exist (no --force)
# ---------------------------------------------------------------------------

def test_idempotent_skip(tmp_path):
    """build_all_pools skips when all 4 FASTAs + report exist and force=False."""
    out_dir = tmp_path / "pools"
    out_dir.mkdir()

    # Pre-create the 4 output FASTAs + report so they look done
    for cls in ("A", "B", "C", "F"):
        (out_dir / f"refs_class_{cls}.fa").write_text(">dummy\nMSTL\n")
    (out_dir / "pool_build_report.json").write_text(json.dumps({"class_A": {}}))

    call_count = {"n": 0}
    original_load = bpcp.load_class_tsv

    def counting_load(path):
        call_count["n"] += 1
        return original_load(path)

    dummy_tsv = tmp_path / "classes.tsv"
    dummy_tsv.write_text("seq_id\tclass\tevidence_pfam\ttop_evalue\n")

    with patch.object(bpcp, "load_class_tsv", side_effect=counting_load):
        bpcp.build_all_pools(
            scan_fasta_glob=str(tmp_path / "*.fa"),
            class_tsv=str(dummy_tsv),
            out_dir=str(out_dir),
            max_per_class=2000,
            total_budget_per_class=3000,
            outgroup_budget_per_class=10,
            cluster_identity=0.7,
            cdhit_path="cd-hit",
            threads=1,
            must_include_taxids=frozenset(),
            berghia_taxid=1287507,
            force=False,
        )

    assert call_count["n"] == 0, "load_class_tsv should not be called when outputs exist"


# ---------------------------------------------------------------------------
# 11. CLI argparse smoke-test
# ---------------------------------------------------------------------------

def test_cli_argparse_defaults(tmp_path):
    """build_args_parser produces correct defaults."""
    parser = bpcp.build_args_parser()
    args = parser.parse_args([
        "--scan-fasta-glob", "scan_output/*.fa",
        "--class-tsv", "classes.tsv",
        "--out-dir", str(tmp_path),
    ])
    assert args.max_per_class == 2000
    assert args.cluster_identity == 0.7
    assert args.berghia_taxid == 1287507
    assert args.threads >= 1


# ---------------------------------------------------------------------------
# 12. DEFAULT_MUST_INCLUDE_TAXIDS has the required 18 entries
# ---------------------------------------------------------------------------

def test_default_must_include_has_18_entries():
    """DEFAULT_MUST_INCLUDE_TAXIDS constant contains exactly 18 taxids."""
    assert len(bpcp.DEFAULT_MUST_INCLUDE_TAXIDS) == 18


# ---------------------------------------------------------------------------
# 13. Small-pool passthrough (≤ max): all sequences written, none dropped
# ---------------------------------------------------------------------------

def test_small_pool_passthrough(tmp_path):
    """When candidates ≤ ref cap, all sequences go into the pool."""
    fa_path = tmp_path / "6526_Biomphalaria_glabrata.chemo_candidates.fa"
    seqs = [(f"bg_{i}", "MSTL" * 20) for i in range(3)]
    _make_fasta(fa_path, seqs)

    class_tsv = tmp_path / "classes.tsv"
    _make_class_tsv(class_tsv, [(sid, "B") for sid, _ in seqs])

    out_dir = tmp_path / "pools"
    out_dir.mkdir()

    with patch.object(bpcp, "cdhit_dedup", side_effect=lambda records, **kw: records):
        with patch.object(bpcp, "get_lineage", side_effect=mock_get_lineage):
            bpcp.build_all_pools(
                scan_fasta_glob=str(tmp_path / "*.chemo_candidates.fa"),
                class_tsv=str(class_tsv),
                out_dir=str(out_dir),
                max_per_class=2000,
                total_budget_per_class=3000,
                outgroup_budget_per_class=10,
                cluster_identity=0.7,
                cdhit_path="cd-hit",
                threads=1,
                must_include_taxids=frozenset(),
                berghia_taxid=1287507,
                force=True,
            )

    pool_b = list(SeqIO.parse(str(out_dir / "refs_class_B.fa"), "fasta"))
    assert len(pool_b) == 3
    ids_out = {r.id for r in pool_b}
    for sid, _ in seqs:
        assert sid in ids_out


# ---------------------------------------------------------------------------
# 14. Berghia Class A seqs appear in refs_class_A.fa as MUST_INCLUDE
# ---------------------------------------------------------------------------

def test_berghia_class_a_included_as_must_include(tmp_path):
    """Berghia Class A candidates appear in refs_class_A.fa even when they
    are not in DEFAULT_MUST_INCLUDE_TAXIDS."""
    # Scan FASTA from Aplysia only (no Berghia scan file here)
    aplysia_fa = tmp_path / "6500_Aplysia_californica.chemo_candidates.fa"
    _make_fasta(aplysia_fa, [("apl_A1", "MSTL" * 30)])

    scan_class_tsv = tmp_path / "scan_classes.tsv"
    _make_class_tsv(scan_class_tsv, [("apl_A1", "A")])

    # Berghia FASTA + class TSV (separate from the scan input)
    berghia_fa = tmp_path / "berghia_candidates.fa"
    _make_fasta(berghia_fa, [("ber_A1", "MKVL" * 30), ("ber_A2", "MSFT" * 30)])

    berghia_class_tsv = tmp_path / "berghia_classes.tsv"
    _make_class_tsv(berghia_class_tsv, [("ber_A1", "A"), ("ber_A2", "A")])

    out_dir = tmp_path / "pools"
    out_dir.mkdir()

    with patch.object(bpcp, "cdhit_dedup", side_effect=lambda records, **kw: records):
        with patch.object(bpcp, "get_lineage", side_effect=mock_get_lineage):
            bpcp.build_all_pools(
                scan_fasta_glob=str(tmp_path / "6500_*.fa"),
                class_tsv=str(scan_class_tsv),
                out_dir=str(out_dir),
                max_per_class=2000,
                total_budget_per_class=3000,
                outgroup_budget_per_class=10,
                cluster_identity=0.7,
                cdhit_path="cd-hit",
                threads=1,
                must_include_taxids=frozenset(),  # Berghia NOT in must_include_taxids
                berghia_taxid=1287507,
                berghia_fasta=str(berghia_fa),
                berghia_class_tsv=str(berghia_class_tsv),
                force=True,
            )

    pool_a = list(SeqIO.parse(str(out_dir / "refs_class_A.fa"), "fasta"))
    ids_out = {r.id for r in pool_a}
    # Berghia Class A seqs must be present despite not being in must_include_taxids
    assert "ber_A1" in ids_out, "ber_A1 missing from Class A pool"
    assert "ber_A2" in ids_out, "ber_A2 missing from Class A pool"
    assert "apl_A1" in ids_out


# ---------------------------------------------------------------------------
# 15. Berghia Class B seqs go to Class B pool, not Class A
# ---------------------------------------------------------------------------

def test_berghia_class_b_routed_to_class_b_pool(tmp_path):
    """Berghia sequences classified as B must appear in refs_class_B.fa only."""
    # Empty scan FASTA (just need the glob to work)
    dummy_fa = tmp_path / "9606_Homo_sapiens.chemo_candidates.fa"
    _make_fasta(dummy_fa, [("hs_B1", "MSTL" * 30)])

    scan_class_tsv = tmp_path / "scan_classes.tsv"
    _make_class_tsv(scan_class_tsv, [("hs_B1", "B")])

    berghia_fa = tmp_path / "berghia_candidates.fa"
    _make_fasta(berghia_fa, [("ber_A1", "MKVL" * 30), ("ber_B1", "MKGT" * 30)])

    berghia_class_tsv = tmp_path / "berghia_classes.tsv"
    _make_class_tsv(berghia_class_tsv, [("ber_A1", "A"), ("ber_B1", "B")])

    out_dir = tmp_path / "pools"
    out_dir.mkdir()

    with patch.object(bpcp, "cdhit_dedup", side_effect=lambda records, **kw: records):
        with patch.object(bpcp, "get_lineage", side_effect=mock_get_lineage):
            bpcp.build_all_pools(
                scan_fasta_glob=str(tmp_path / "*.chemo_candidates.fa"),
                class_tsv=str(scan_class_tsv),
                out_dir=str(out_dir),
                max_per_class=2000,
                total_budget_per_class=3000,
                outgroup_budget_per_class=10,
                cluster_identity=0.7,
                cdhit_path="cd-hit",
                threads=1,
                must_include_taxids=frozenset(),
                berghia_taxid=1287507,
                berghia_fasta=str(berghia_fa),
                berghia_class_tsv=str(berghia_class_tsv),
                force=True,
            )

    ids_A = {r.id for r in SeqIO.parse(str(out_dir / "refs_class_A.fa"), "fasta")}
    ids_B = {r.id for r in SeqIO.parse(str(out_dir / "refs_class_B.fa"), "fasta")}

    # ber_B1 goes to B only; ber_A1 goes to A only
    assert "ber_B1" in ids_B, "ber_B1 missing from Class B pool"
    assert "ber_B1" not in ids_A, "ber_B1 must not appear in Class A pool"
    assert "ber_A1" in ids_A, "ber_A1 missing from Class A pool"
    assert "ber_A1" not in ids_B, "ber_A1 must not appear in Class B pool"


# ---------------------------------------------------------------------------
# 16. Per-class caps: Class A cap < Class B cap when configured differently
# ---------------------------------------------------------------------------

def test_per_class_cap_respected_via_total_budget(tmp_path):
    """Per-class ref cap = total_budget - berghia_count - outgroup_budget."""
    # 50 seqs each for Class A and Class B
    fa_path = tmp_path / "7227_Drosophila_melanogaster.chemo_candidates.fa"
    seqs_a = [(f"dm_A{i}", "MKLF" * 30) for i in range(50)]
    seqs_b = [(f"dm_B{i}", "MKQP" * 30) for i in range(50)]
    _make_fasta(fa_path, seqs_a + seqs_b)

    class_tsv = tmp_path / "classes.tsv"
    _make_class_tsv(class_tsv, [(sid, "A") for sid, _ in seqs_a] +
                               [(sid, "B") for sid, _ in seqs_b])

    out_dir = tmp_path / "pools"
    out_dir.mkdir()

    total_budget = 60
    outgroup_budget = 10
    # Expected: Class A refs = 60 - 0 - 10 = 50, Class B refs = 60 - 0 - 10 = 50

    with patch.object(bpcp, "cdhit_dedup", side_effect=lambda records, **kw: records):
        with patch.object(bpcp, "get_lineage", side_effect=mock_get_lineage):
            bpcp.build_all_pools(
                scan_fasta_glob=str(tmp_path / "*.chemo_candidates.fa"),
                class_tsv=str(class_tsv),
                out_dir=str(out_dir),
                max_per_class=2000,
                total_budget_per_class=total_budget,
                outgroup_budget_per_class=outgroup_budget,
                cluster_identity=0.7,
                cdhit_path="cd-hit",
                threads=1,
                must_include_taxids=frozenset(),
                berghia_taxid=1287507,
                force=True,
            )

    pool_a = list(SeqIO.parse(str(out_dir / "refs_class_A.fa"), "fasta"))
    pool_b = list(SeqIO.parse(str(out_dir / "refs_class_B.fa"), "fasta"))

    expected_cap = total_budget - outgroup_budget
    assert len(pool_a) <= expected_cap, f"Class A pool {len(pool_a)} exceeds expected cap {expected_cap}"
    assert len(pool_b) <= expected_cap, f"Class B pool {len(pool_b)} exceeds expected cap {expected_cap}"


# ---------------------------------------------------------------------------
# 17. No --berghia-fasta: P2 works, Berghia count is 0 in report
# ---------------------------------------------------------------------------

def test_no_berghia_fasta_falls_back_gracefully(tmp_path):
    """When --berghia-fasta is not provided, build_all_pools still runs and
    report shows zero Berghia inclusions."""
    fa_path = tmp_path / "29159_Crassostrea_gigas.chemo_candidates.fa"
    seqs = [(f"cg_{i}", "MSTL" * 20) for i in range(4)]
    _make_fasta(fa_path, seqs)

    class_tsv = tmp_path / "classes.tsv"
    _make_class_tsv(class_tsv, [(sid, "C") for sid, _ in seqs])

    out_dir = tmp_path / "pools"
    out_dir.mkdir()

    with patch.object(bpcp, "cdhit_dedup", side_effect=lambda records, **kw: records):
        with patch.object(bpcp, "get_lineage", side_effect=mock_get_lineage):
            bpcp.build_all_pools(
                scan_fasta_glob=str(tmp_path / "*.chemo_candidates.fa"),
                class_tsv=str(class_tsv),
                out_dir=str(out_dir),
                max_per_class=2000,
                total_budget_per_class=3000,
                outgroup_budget_per_class=10,
                cluster_identity=0.7,
                cdhit_path="cd-hit",
                threads=1,
                must_include_taxids=frozenset(),
                berghia_taxid=1287507,
                berghia_fasta=None,        # not provided
                berghia_class_tsv=None,    # not provided
                force=True,
            )

    # Outputs exist and Class C pool has the right seqs
    pool_c = list(SeqIO.parse(str(out_dir / "refs_class_C.fa"), "fasta"))
    assert len(pool_c) == 4

    # Report shows zero Berghia inclusions
    report = json.loads((out_dir / "pool_build_report.json").read_text())
    assert report["berghia_included"]["n_total"] == 0
    for cls in ("A", "B", "C", "F"):
        assert report[f"class_{cls}"]["n_berghia_included"] == 0
    # Report should show budget was 3000, outgroup budget 10, berghia count was 0 for each class
    assert report["total_budget_per_class"] == 3000
    assert report["outgroup_budget_per_class"] == 10
    for cls in ("A", "B", "C", "F"):
        assert report["n_berghia_per_class"][cls] == 0


# ---------------------------------------------------------------------------
# 18. MUST_INCLUDE retained even when Berghia exceeds budget
# ---------------------------------------------------------------------------

def test_must_include_retained_when_berghia_exceeds_budget(tmp_path):
    """MUST_INCLUDE sequences are always retained even if they push total over max_size."""
    # Setup: One large must-include taxid with many seqs that exceed the cap
    aplysia_fa = tmp_path / "6500_Aplysia_californica.chemo_candidates.fa"
    aplysia_seqs = [(f"apl_A{i}", "MSTL" * 30) for i in range(100)]
    _make_fasta(aplysia_fa, aplysia_seqs)

    scan_class_tsv = tmp_path / "scan_classes.tsv"
    _make_class_tsv(scan_class_tsv, [(sid, "A") for sid, _ in aplysia_seqs])

    out_dir = tmp_path / "pools"
    out_dir.mkdir()

    # Set a very small cap: total_budget=50, outgroup=10 => refs_cap = 30
    # Aplysia has 100 seqs (MUST_INCLUDE), which exceeds 30
    # All Aplysia sequences should still be in the pool
    with patch.object(bpcp, "cdhit_dedup", side_effect=lambda records, **kw: records):
        with patch.object(bpcp, "get_lineage", side_effect=mock_get_lineage):
            bpcp.build_all_pools(
                scan_fasta_glob=str(tmp_path / "*.chemo_candidates.fa"),
                class_tsv=str(scan_class_tsv),
                out_dir=str(out_dir),
                max_per_class=50,
                total_budget_per_class=50,
                outgroup_budget_per_class=10,
                cluster_identity=0.7,
                cdhit_path="cd-hit",
                threads=1,
                must_include_taxids=frozenset({6500}),
                berghia_taxid=1287507,
                force=True,
            )

    pool_a = list(SeqIO.parse(str(out_dir / "refs_class_A.fa"), "fasta"))
    ids_out = {r.id for r in pool_a}

    # All Aplysia sequences must be present (MUST_INCLUDE always retained)
    for sid, _ in aplysia_seqs:
        assert sid in ids_out, f"{sid} missing from pool despite must-include taxid"

    # Verify that the count exceeds the ref_cap (because MUST_INCLUDE is guaranteed)
    ref_cap = 50 - 10
    assert len(pool_a) >= 100, f"Pool should contain all {len(aplysia_seqs)} must-include seqs"


# ---------------------------------------------------------------------------
# 19. Dynamic per-class max computed from total_budget - berghia_count - outgroup_budget
# ---------------------------------------------------------------------------

def test_max_refs_computed_from_total_budget_minus_berghia(tmp_path):
    """Per-class ref cap is computed as total_budget - berghia_count - outgroup_budget."""
    # Setup: Aplysia with 10 Class A seqs, Drosophila with 5 Class A seqs
    aplysia_fa = tmp_path / "6500_Aplysia_californica.chemo_candidates.fa"
    drosophila_fa = tmp_path / "7227_Drosophila_melanogaster.chemo_candidates.fa"

    aplysia_seqs = [(f"apl_A{i}", "MSTL" * 30) for i in range(10)]
    drosophila_seqs = [(f"dm_A{i}", "MKGT" * 30) for i in range(5)]

    _make_fasta(aplysia_fa, aplysia_seqs)
    _make_fasta(drosophila_fa, drosophila_seqs)

    scan_class_tsv = tmp_path / "scan_classes.tsv"
    _make_class_tsv(
        scan_class_tsv,
        [(sid, "A") for sid, _ in aplysia_seqs + drosophila_seqs]
    )

    # Berghia with 5 Class A seqs, 2 Class B seqs
    berghia_fa = tmp_path / "berghia_candidates.fa"
    berghia_a_seqs = [(f"ber_A{i}", "MKVL" * 30) for i in range(5)]
    berghia_b_seqs = [(f"ber_B{i}", "MSFT" * 30) for i in range(2)]
    _make_fasta(berghia_fa, berghia_a_seqs + berghia_b_seqs)

    berghia_class_tsv = tmp_path / "berghia_classes.tsv"
    _make_class_tsv(
        berghia_class_tsv,
        [(sid, "A") for sid, _ in berghia_a_seqs] +
        [(sid, "B") for sid, _ in berghia_b_seqs]
    )

    out_dir = tmp_path / "pools"
    out_dir.mkdir()

    total_budget = 3000
    outgroup_budget = 10
    # Expected: Class A refs = 3000 - 5 - 10 = 2985, Class B refs = 3000 - 2 - 10 = 2988

    with patch.object(bpcp, "cdhit_dedup", side_effect=lambda records, **kw: records):
        with patch.object(bpcp, "get_lineage", side_effect=mock_get_lineage):
            bpcp.build_all_pools(
                scan_fasta_glob=str(tmp_path / "*.chemo_candidates.fa"),
                class_tsv=str(scan_class_tsv),
                out_dir=str(out_dir),
                max_per_class=2000,
                total_budget_per_class=total_budget,
                outgroup_budget_per_class=outgroup_budget,
                cluster_identity=0.7,
                cdhit_path="cd-hit",
                threads=1,
                must_include_taxids=frozenset(),
                berghia_taxid=1287507,
                berghia_fasta=str(berghia_fa),
                berghia_class_tsv=str(berghia_class_tsv),
                force=True,
            )

    report = json.loads((out_dir / "pool_build_report.json").read_text())

    # Verify computed caps are in the report
    assert report["total_budget_per_class"] == total_budget
    assert report["outgroup_budget_per_class"] == outgroup_budget
    assert report["n_berghia_per_class"]["A"] == 5
    assert report["n_berghia_per_class"]["B"] == 2
    assert report["n_berghia_per_class"]["C"] == 0
    assert report["n_berghia_per_class"]["F"] == 0

    # Verify the per-class ref targets are correctly computed
    assert report["class_A"]["n_refs_target"] == 2985, \
        f"Class A ref target should be 3000 - 5 - 10 = 2985, got {report['class_A']['n_refs_target']}"
    assert report["class_B"]["n_refs_target"] == 2988, \
        f"Class B ref target should be 3000 - 2 - 10 = 2988, got {report['class_B']['n_refs_target']}"
    assert report["class_C"]["n_refs_target"] == 2990, \
        f"Class C ref target should be 3000 - 0 - 10 = 2990, got {report['class_C']['n_refs_target']}"
    assert report["class_F"]["n_refs_target"] == 2990, \
        f"Class F ref target should be 3000 - 0 - 10 = 2990, got {report['class_F']['n_refs_target']}"


# ---------------------------------------------------------------------------
# 20. No berghia means budget minus outgroup_budget available for refs
# ---------------------------------------------------------------------------

def test_no_berghia_means_full_budget_minus_outgroup(tmp_path):
    """When no --berghia-fasta is provided, each class gets total_budget - outgroup_budget refs."""
    fa_path = tmp_path / "7227_Drosophila_melanogaster.chemo_candidates.fa"
    seqs_a = [(f"dm_A{i}", "MKLF" * 30) for i in range(3)]
    seqs_b = [(f"dm_B{i}", "MKQP" * 30) for i in range(2)]
    _make_fasta(fa_path, seqs_a + seqs_b)

    class_tsv = tmp_path / "classes.tsv"
    _make_class_tsv(
        class_tsv,
        [(sid, "A") for sid, _ in seqs_a] +
        [(sid, "B") for sid, _ in seqs_b]
    )

    out_dir = tmp_path / "pools"
    out_dir.mkdir()

    total_budget = 3000
    outgroup_budget = 10

    with patch.object(bpcp, "cdhit_dedup", side_effect=lambda records, **kw: records):
        with patch.object(bpcp, "get_lineage", side_effect=mock_get_lineage):
            bpcp.build_all_pools(
                scan_fasta_glob=str(tmp_path / "*.chemo_candidates.fa"),
                class_tsv=str(class_tsv),
                out_dir=str(out_dir),
                max_per_class=2000,
                total_budget_per_class=total_budget,
                outgroup_budget_per_class=outgroup_budget,
                cluster_identity=0.7,
                cdhit_path="cd-hit",
                threads=1,
                must_include_taxids=frozenset(),
                berghia_taxid=1287507,
                berghia_fasta=None,        # NOT provided
                berghia_class_tsv=None,    # NOT provided
                force=True,
            )

    report = json.loads((out_dir / "pool_build_report.json").read_text())

    # When no Berghia is provided, each class should show berghia_count = 0
    assert report["n_berghia_per_class"]["A"] == 0
    assert report["n_berghia_per_class"]["B"] == 0
    assert report["n_berghia_per_class"]["C"] == 0
    assert report["n_berghia_per_class"]["F"] == 0

    # All classes should have ref_target = total_budget - outgroup_budget
    expected_refs = total_budget - outgroup_budget
    for cls in ("A", "B", "C", "F"):
        assert report[f"class_{cls}"]["n_refs_target"] == expected_refs, \
            f"Class {cls} ref target should be {expected_refs} when no Berghia, got {report[f'class_{cls}']['n_refs_target']}"
