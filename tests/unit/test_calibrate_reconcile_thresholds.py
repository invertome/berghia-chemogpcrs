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
