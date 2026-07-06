import pytest
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import calibrate_genome_track_margin as cm
import reconcile_candidates as rc

def test_locus_join():
    h = ">lcl|NC_088360.1_cds_XP_078753709.1_1 [gene=L] [location=join(4611..4721,5584..5704,6641..6776)] [gbkey=CDS]"
    assert cm.parse_cds_true_locus(h) == ("NC_088360.1", 4611, 6776)

def test_locus_complement_partial():
    h = ">lcl|NC_088360.1_cds_XP_1.1_2 [location=complement(join(40140..40160,48651..>48764))]"
    assert cm.parse_cds_true_locus(h) == ("NC_088360.1", 40140, 48764)

def test_locus_single_exon():
    h = ">lcl|NW_999.1_cds_XP_9.1_3 [location=100..250]"
    assert cm.parse_cds_true_locus(h) == ("NW_999.1", 100, 250)

def test_locus_missing_fields():
    with pytest.raises(ValueError):
        cm.parse_cds_true_locus(">lcl|NC_1.1_cds_XP_9.1_1 [gene=x]")   # no [location=]
    with pytest.raises(ValueError):
        cm.parse_cds_true_locus(">bad_header_no_cds [location=1..2]")   # no _cds_


def _P(chrom, start, end, pid):
    return rc.Placement("q", chrom, start, end, "+", pid, 100.0, "minimap2", 60.0)

def test_margins_single_locus_inf():
    reps = cm.distinct_locus_reps([_P("chrA", 100, 200, 99.9)])
    obs = cm.query_margins(reps, ("chrA", 100, 200))
    assert obs.tp == cm.INF and obs.fp is None

def test_margins_tp_plus_synthetic_fp():
    pls = [_P("chrA", 100, 200, 99.9), _P("chrA", 5000, 5100, 97.0), _P("chrA", 9000, 9100, 95.0)]
    reps = cm.distinct_locus_reps(pls)
    obs = cm.query_margins(reps, ("chrA", 100, 200))            # true = top locus
    assert obs.tp == pytest.approx(2.9)                          # 99.9 - 97.0
    assert obs.fp == pytest.approx(2.0) and obs.fp_kind == "synthetic"  # 97.0 - 95.0

def test_margins_natural_fp_when_paralog_outranks_true():
    pls = [_P("chrA", 5000, 5100, 99.0), _P("chrA", 100, 200, 98.0)]  # true locus is the 98.0 one
    reps = cm.distinct_locus_reps(pls)
    obs = cm.query_margins(reps, ("chrA", 100, 200))
    assert obs.tp is None
    assert obs.fp == pytest.approx(1.0) and obs.fp_kind == "natural"

def test_margins_busco_mode_rep0_is_truth():
    pls = [_P("chrA", 100, 200, 99.5), _P("chrA", 5000, 5100, 96.0)]
    reps = cm.distinct_locus_reps(pls)
    obs = cm.query_margins(reps, None)
    assert obs.tp == pytest.approx(3.5) and obs.fp is None       # only 1 paralog -> no synthetic FP


def test_roc_inf_tp_always_retained():
    pts = cm.roc_points([cm.INF, cm.INF, 5.0], [1.0], [0.0, 3.0, 6.0])
    ret_at_6 = [p for p in pts if p.m == 6.0][0].retention
    assert ret_at_6 == pytest.approx(2/3)      # the two inf margins still retained at m=6

def test_knee_clean_separation():
    tp = [8.0, 9.0, 10.0, 11.0]                # TP floor 8
    fp = [1.0, 2.0, 3.0]                        # FP ceiling 3
    rec = cm.recommend(tp, fp)
    assert 3.0 < rec.min_margin < 8.0 and not rec.split_advisory

def test_knee_overlap_flags_split():
    tp = [1.0, 2.0, 3.0, 4.0]                  # TP and FP overlap heavily
    fp = [1.0, 2.0, 3.0, 4.0]
    rec = cm.recommend(tp, fp)
    assert rec.split_advisory                    # no clean knee


def test_cli_dry_run_reproduces_clean_knee(capsys):
    rc_code = cm.main(["--dry-run"])
    out = capsys.readouterr().out
    assert rc_code == 0 and "GENOME_TRACK_MIN_MARGIN=" in out
    val = float([l for l in out.splitlines() if "GENOME_TRACK_MIN_MARGIN=" in l][0].split("=")[1])
    assert 3.0 < val < 8.0

def test_cli_writes_sourceable_recommendation(tmp_path):
    out = tmp_path / "rec.sh"
    cm.main(["--dry-run", "--out-recommendation", str(out)])
    text = out.read_text()
    assert text.startswith("#") and "GENOME_TRACK_MIN_MARGIN=" in text
