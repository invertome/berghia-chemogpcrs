"""Tests for scripts/calibrate_reconcile_thresholds.py — Task 8.

The calibrator measures the transcript->genome %id/%cov distribution on the
Berghia BUSCO single-copy gene set (the "true same-gene" distribution),
estimates the paralog-family %id ceiling, and recommends a (min_id, min_cov)
accept bar that sits BELOW the true-same-gene floor but ABOVE the paralog
ceiling (design §9). These tests drive the pure computation + the --dry-run
path; the REAL run against Unity alignment files is done separately.
"""
from __future__ import annotations

import calibrate_reconcile_thresholds as cal
import reconcile_candidates as rc


def _pl(query, chrom, start, end, pid, cov=95.0, strand="+"):
    return rc.Placement(query, chrom, start, end, strand, pid, cov, "minimap2", 60.0)


# ---- _percentile (linear interpolation) -----------------------------------

def test_percentile_endpoints_and_median():
    vals = [10.0, 20.0, 30.0, 40.0, 50.0]
    assert cal._percentile(vals, 0) == 10.0
    assert cal._percentile(vals, 100) == 50.0
    assert cal._percentile(vals, 50) == 30.0


def test_percentile_linear_interpolates_between_ranks():
    # 5th percentile of two points sits 5% of the way from the low value.
    assert abs(cal._percentile([0.0, 100.0], 5) - 5.0) < 1e-9


def test_percentile_single_value():
    assert cal._percentile([42.0], 5) == 42.0


# ---- _distinct_locus_reps: self vs paralog cross-maps ----------------------

def test_distinct_locus_reps_orders_self_then_paralogs():
    # One transcript: near-perfect self hit + two lower paralog cross-maps at
    # distinct loci -> one representative per locus, highest identity first.
    placements = [
        _pl("t", "chr1", 1000, 2000, 99.0),   # self
        _pl("t", "chr1", 5000, 6000, 90.0),   # paralog
        _pl("t", "chr1", 9000, 10000, 88.0),  # paralog
    ]
    reps = cal._distinct_locus_reps(placements)
    assert [r.pct_identity for r in reps] == [99.0, 90.0, 88.0]


def test_distinct_locus_reps_folds_overlapping_hits_into_one_locus():
    # Two overlapping hits at ONE locus (99, 97) must collapse to a single
    # representative (the 99); only the genuinely distinct locus (90) is a
    # second rep. A within-locus second alignment is NOT a paralog.
    placements = [
        _pl("t", "chr1", 1000, 2000, 99.0),
        _pl("t", "chr1", 1500, 2500, 97.0),   # overlaps the 99 locus
        _pl("t", "chr1", 5000, 6000, 90.0),   # distinct locus
    ]
    reps = cal._distinct_locus_reps(placements)
    assert [r.pct_identity for r in reps] == [99.0, 90.0]


# ---- same_gene_reps / paralog_identities -----------------------------------

def test_same_gene_reps_takes_self_placement_per_id():
    # For a single-copy (BUSCO) transcript, the self placement is its best
    # (highest-identity) locus. Absent ids contribute nothing.
    placements_by_query = {
        "busco1": [_pl("busco1", "chr1", 1000, 2000, 99.5, cov=98.0),
                   _pl("busco1", "chr1", 8000, 9000, 92.0, cov=97.0)],
    }
    reps = cal.same_gene_reps(placements_by_query, ["busco1", "absent"])
    assert len(reps) == 1
    assert reps[0].pct_identity == 99.5 and reps[0].pct_coverage == 98.0


def test_paralog_identities_excludes_the_self_locus():
    # For a paralog-family transcript, the self (best) locus is dropped; only
    # the cross-map identities form the paralog distribution.
    placements_by_query = {
        "para1": [_pl("para1", "chr1", 1000, 2000, 98.0),   # self
                  _pl("para1", "chr1", 5000, 6000, 90.0),   # paralog
                  _pl("para1", "chr1", 9000, 10000, 88.0)],  # paralog
    }
    ids = cal.paralog_identities(placements_by_query, ["para1"])
    assert sorted(ids) == [88.0, 90.0]


# ---- compute_separation ----------------------------------------------------

def test_compute_separation_recommends_between_ceiling_and_floor():
    same_reps = [_pl("b", "c", 1, 2, pid, cov=cov) for pid, cov in
                 [(98.0, 95.0), (99.0, 96.0), (100.0, 97.0),
                  (97.9, 98.0), (98.5, 99.0)]]
    paralog_ids = [88.0, 90.0, 91.0, 92.0]
    res = cal.compute_separation(same_reps, paralog_ids)
    assert res.clean_separation is True
    # The accept bar sits strictly between the paralog ceiling and the
    # true-same-gene identity floor (p05) — the whole point of the calibration.
    assert res.paralog_ceiling < res.rec_min_id < res.id_dist.p05
    assert res.separation_margin > 0
    # Recommended coverage bar is at/below the true-same-gene coverage floor.
    assert res.rec_min_cov <= res.cov_dist.p05
    assert res.n_same_gene == 5 and res.n_paralog == 4


def test_compute_separation_flags_no_clean_separation():
    # Paralog identities that reach INTO the true-same-gene range -> no clean
    # bar exists; the calibrator refuses to invent one and falls back to the
    # working hypothesis, flagged so the operator sees it.
    same_reps = [_pl("b", "c", 1, 2, pid, cov=cov) for pid, cov in
                 [(96.0, 95.0), (97.0, 96.0), (98.0, 97.0)]]
    paralog_ids = [97.0, 98.0, 99.0]
    res = cal.compute_separation(same_reps, paralog_ids)
    assert res.clean_separation is False
    assert res.rec_min_id == cal.WORKING_MIN_ID
    assert res.rec_min_cov == cal.WORKING_MIN_COV


def test_compute_separation_requires_same_gene_data():
    # No BUSCO self-hits at all -> cannot calibrate; fail loud rather than
    # emit a meaningless bar.
    import pytest
    with pytest.raises(ValueError):
        cal.compute_separation([], [90.0])


# ---- format_table ----------------------------------------------------------

def test_format_table_reports_floor_ceiling_and_recommendation():
    same_reps = [_pl("b", "c", 1, 2, pid, cov=97.0) for pid in
                 [98.0, 99.0, 100.0, 97.9, 98.5]]
    res = cal.compute_separation(same_reps, [88.0, 90.0, 91.0, 92.0])
    table = cal.format_table(res)
    # Numbers are cited, and the two decision quantities are named.
    assert "paralog" in table.lower()
    assert "min_id" in table and "min_cov" in table
    assert f"{res.rec_min_id:g}" in table


# ---- --dry-run CLI ---------------------------------------------------------

def test_dry_run_prints_table_and_returns_zero(capsys):
    rc_code = cal.main(["--dry-run"])
    assert rc_code == 0
    out = capsys.readouterr().out
    assert "min_id" in out and "min_cov" in out
    assert "paralog" in out.lower()


def test_dry_run_example_shows_clean_separation():
    # The synthetic dry-run distribution is a WORKED EXAMPLE of the design's
    # separation logic: it must exhibit a clean gap (the recommended min_id
    # lands between the paralog ceiling and the true-same-gene floor).
    placements, busco_ids, paralog_ids = cal.synthetic_inputs()
    res = cal.compute_separation(cal.same_gene_reps(placements, busco_ids),
                                 cal.paralog_identities(placements, paralog_ids))
    assert res.clean_separation is True
    assert res.paralog_ceiling < res.rec_min_id < res.id_dist.p05


# ---- Minor 1: distinct-locus logic must MATCH the module it calibrates ------

def test_distinct_locus_reps_transitive_chain_is_one_locus():
    # A–B–C overlap CHAIN: A overlaps B, B overlaps C, but A is DISJOINT from C.
    # The module clusters by connected components (transitive), so the chain is
    # ONE locus -> one representative (A@99). A greedy pairwise scan would wrongly
    # keep C@98 as a second "distinct locus", inflating the paralog ceiling with
    # a within-gene alignment. Pins that the calibrator uses the SAME distinct-
    # locus definition as reconcile_candidates.py.
    placements = [
        _pl("t", "chr1", 1000, 2000, 99.0),   # A
        _pl("t", "chr1", 1990, 3000, 95.0),   # B bridges A and C
        _pl("t", "chr1", 2990, 4000, 98.0),   # C (disjoint from A)
    ]
    reps = cal._distinct_locus_reps(placements)
    assert [r.pct_identity for r in reps] == [99.0]


def test_distinct_locus_reps_matches_module_clustering():
    # Independent cross-check: the reps' loci equal the module's own clustering
    # of the same placements (one rep per module cluster).
    placements = [
        _pl("t", "chr1", 1000, 2000, 99.0),
        _pl("t", "chr1", 1990, 3000, 95.0),
        _pl("t", "chr1", 8000, 9000, 90.0),   # a genuinely separate locus
    ]
    reps = cal._distinct_locus_reps(placements)
    assert len(reps) == len(rc._cluster_placements(placements))
    assert [r.pct_identity for r in reps] == [99.0, 90.0]


# ---- Minor 2: strict-between invariant / degenerate near-touch gap ----------

def test_compute_separation_subthreshold_gap_is_not_clean():
    # ceiling 95.0, floor 95.05 -> gap 0.05 < the minimum separation gap: the
    # distributions essentially touch. Without the min-gap guard, midpoint
    # rounding (round(95.025, 1) == 95.0) would emit rec_min_id == ceiling yet
    # declare clean_separation True — violating paralog_ceiling < min_id < floor.
    same_reps = [_pl("b", "c", 1, 2, 95.05, cov=95.0) for _ in range(4)]
    paralog_ids = [95.0, 95.0, 95.0, 95.0]
    res = cal.compute_separation(same_reps, paralog_ids)
    assert res.clean_separation is False
    assert res.rec_min_id == cal.WORKING_MIN_ID


def test_compute_separation_clean_recommendation_is_strictly_interior():
    # Whenever a clean bar is declared, the recommended min_id must lie STRICTLY
    # between the paralog ceiling and the true-same-gene floor (never on either).
    same_reps = [_pl("b", "c", 1, 2, pid, cov=97.0) for pid in
                 [98.0, 99.0, 100.0, 97.9, 98.5]]
    res = cal.compute_separation(same_reps, [88.0, 90.0, 91.0, 92.0])
    assert res.clean_separation is True
    assert res.paralog_ceiling < res.rec_min_id < res.id_dist.p05


# ---- Minor 5: empty-paralog branch ------------------------------------------

def test_compute_separation_no_paralog_data_is_not_clean():
    same_reps = [_pl("b", "c", 1, 2, 99.0, cov=98.0)]
    res = cal.compute_separation(same_reps, [])
    assert res.paralog_ceiling is None
    assert res.separation_margin is None
    assert res.clean_separation is False
    assert res.rec_min_id == cal.WORKING_MIN_ID
    assert res.rec_min_cov == cal.WORKING_MIN_COV


# ---- Minor 5: _read_ids parsing ---------------------------------------------

def test_read_ids_skips_comments_blanks_and_takes_first_token(tmp_path):
    f = tmp_path / "ids.txt"
    f.write_text("# comment\n\nq1\nq2 extra tokens\n  q3  \n")
    assert cal._read_ids(str(f)) == ["q1", "q2", "q3"]


# ---- Minor 3: partial id-match warning --------------------------------------

def test_matched_fraction_warning_levels():
    assert cal.matched_fraction_warning("same-gene", 5, 5) is None
    assert cal.matched_fraction_warning("same-gene", 0, 0) is None
    assert "matched 3 of 4 requested same-gene ids" in \
        cal.matched_fraction_warning("same-gene", 3, 4)
    assert "none" in cal.matched_fraction_warning("paralog", 0, 4).lower()
    assert "low" in cal.matched_fraction_warning("paralog", 1, 4).lower()


def test_real_mode_warns_on_partial_id_match(tmp_path, capsys):
    # q3 is requested but absent from the alignment -> a silent tiny-distribution
    # bar would otherwise result. Also exercises _read_ids (comment + blank).
    paf = tmp_path / "aln.paf"
    # Identity comes from the de:f divergence tag (single-exon synthetic rows,
    # so de:f matches the old nmatch/alen value); only the id-match count is
    # asserted below.
    paf.write_text(
        "q1\t1000\t0\t1000\t+\tchr1\t100000\t5000\t6000\t990\t1000\t60\tde:f:0.01\n"
        "q2\t1000\t0\t900\t+\tchr1\t100000\t20000\t20900\t880\t900\t60\tde:f:0.0222\n")
    ids = tmp_path / "busco.txt"
    ids.write_text("# single-copy ids\n\nq1\nq2\nq3\n")
    code = cal.main(["--minimap2-paf", str(paf),
                     "--busco-single-copy-ids", str(ids)])
    assert code == 0
    err = capsys.readouterr().err
    assert "matched 2 of 3 requested same-gene ids" in err


# ---- Minor 4: dry-run --out must be marked synthetic -------------------------

def test_dry_run_out_file_marked_synthetic_and_sets_vars(tmp_path):
    out = tmp_path / "rec.sh"
    assert cal.main(["--dry-run", "--out", str(out)]) == 0
    text = out.read_text()
    # Clearly marked NOT for production...
    assert "DRY RUN synthetic" in text and "NOT for production" in text
    # ...yet still sets both vars and cites the decision numbers (closes the
    # "--out content never asserted" gap).
    assert "GENOME_TRACK_MIN_ID=" in text and "GENOME_TRACK_MIN_COV=" in text
    assert "floor" in text and "ceiling" in text and "clean_separation" in text
