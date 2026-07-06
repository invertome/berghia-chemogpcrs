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


# ---- review fix 1: grid pads one full step beyond the max data value --------

def test_grid_last_point_exceeds_max_margin():
    # 8.1/0.1 == 80.9999... in float; int() would truncate and clip the grid at
    # 8.0 (== max data), losing the tail of the curve. round() keeps the pad.
    grid = cm._grid([8.0], [3.0])
    assert grid[-1] > 8.0


# ---- review fix 2: true locus overlaps NONE of the placements ---------------

def test_margins_true_locus_matches_no_placement():
    reps = cm.distinct_locus_reps([_P("chrA", 5000, 5100, 99.0), _P("chrA", 9000, 9100, 95.0)])
    obs = cm.query_margins(reps, ("chrX", 1, 2))          # true gene mapped nowhere here
    assert obs.tp is None
    assert obs.fp == pytest.approx(4.0) and obs.fp_kind == "synthetic"


# ---- review fix 3: real (non-dry-run) file-mode CLI, both modes -------------

def _paf_rows():
    # Gene A: true locus (99) + two distinct paralog loci (96, 94) -> TP margin
    #   3.0 and synthetic FP margin 2.0. Gene B: one confident locus -> TP inf.
    # Query names are the full lcl| CDS ids so they join the CDS FASTA headers.
    rows = [
        ("lcl|NC_1.1_cds_XPA_1", 100, 0, 100, "+", "NC_1.1", 1000000, 100, 200, 99, 100, 60),
        ("lcl|NC_1.1_cds_XPA_1", 100, 0, 100, "+", "NC_1.1", 1000000, 5000, 5100, 96, 100, 60),
        ("lcl|NC_1.1_cds_XPA_1", 100, 0, 100, "+", "NC_1.1", 1000000, 9000, 9100, 94, 100, 60),
        ("lcl|NC_1.1_cds_XPB_1", 100, 0, 100, "+", "NC_1.1", 1000000, 20000, 20100, 98, 100, 60),
    ]
    return "\n".join("\t".join(str(x) for x in r) for r in rows) + "\n"

def _cds_fasta():
    return (">lcl|NC_1.1_cds_XPA_1 [location=join(100..150,160..200)]\nMKAA\n"
            ">lcl|NC_1.1_cds_XPB_1 [location=20000..20100]\nMKBB\n")

def test_cli_refseq_mode_file_inputs(tmp_path, capsys):
    paf = tmp_path / "aln.paf"; paf.write_text(_paf_rows())
    cds = tmp_path / "gpcr_cds.fna"; cds.write_text(_cds_fasta())
    code = cm.main(["--mode", "refseq", "--minimap2-paf", str(paf), "--refseq-cds", str(cds)])
    out = capsys.readouterr().out
    assert code == 0
    assert "Margin-gate calibration" in out and "GENOME_TRACK_MIN_MARGIN=" in out

def test_cli_busco_mode_file_inputs(tmp_path, capsys):
    paf = tmp_path / "aln.paf"; paf.write_text(_paf_rows())
    ids = tmp_path / "busco.ids"
    ids.write_text("# busco single-copy\n\nlcl|NC_1.1_cds_XPA_1\nlcl|NC_1.1_cds_XPB_1\n")
    code = cm.main(["--mode", "busco", "--minimap2-paf", str(paf), "--busco-ids", str(ids)])
    out = capsys.readouterr().out
    assert code == 0 and "GENOME_TRACK_MIN_MARGIN=" in out


# ---- review fix 4: collect_margins match-fraction safeguard -----------------

def test_collect_margins_warns_on_partial_id_match(capsys):
    # 1 of 4 requested CDS ids joins the alignment -> a silent tiny-sample knee
    # would otherwise result. Surfaces the requested denominator (4), not n_tp.
    by_q = {"lcl|NC_1.1_cds_XPA_1": [_P("chrA", 100, 200, 99.0)]}
    true_loci = {"lcl|NC_1.1_cds_XPA_1": ("chrA", 100, 200),
                 "lcl|NC_1.1_cds_XPB_1": ("chrB", 1, 2),
                 "lcl|NC_1.1_cds_XPC_1": ("chrC", 1, 2),
                 "lcl|NC_1.1_cds_XPD_1": ("chrD", 1, 2)}
    cm.collect_margins(by_q, true_loci)
    err = capsys.readouterr().err
    assert "matched 1 of 4 requested refseq ids" in err


# ---- review fix 5: overlap -> split_advisory + GENOME_TRACK_LOW_MARGIN -------

def test_recommend_overlap_splits_and_emits_low_margin(tmp_path):
    rec = cm.recommend([1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0])
    assert rec.split_advisory is True and rec.low_margin_threshold is not None
    assert "GENOME_TRACK_LOW_MARGIN=" in cm.format_table(rec)
    out = tmp_path / "rec.sh"
    cm.write_recommendation(str(out), rec)
    assert "GENOME_TRACK_LOW_MARGIN=" in out.read_text()


# ---- review fix 6: write_figure smoke test (matplotlib installed) -----------

def test_write_figure_creates_png(tmp_path):
    tp, fp = cm._synthetic()
    rec = cm.recommend(tp, fp)
    out = tmp_path / "roc.png"
    cm.write_figure(str(out), tp, fp, rec)
    assert out.exists() and out.stat().st_size > 0


# ---- BUSCO-run fix: knee is built on the multi-locus (finite) population -----

def test_recommend_knee_ignores_infinite_single_locus_tp():
    # 10 single-locus genes (margin inf) must NOT pull the knee up; the knee is
    # built from the multi-locus TP {8,9} vs FP {1,2} (design §2). Without this
    # the inf-dominated retention pushes the knee to an artifactual ~99.5.
    rec = cm.recommend([cm.INF] * 10 + [8.0, 9.0], [1.0, 2.0])
    assert 2.0 < rec.min_margin < 8.0
    assert rec.n_tp == 12 and rec.n_tp_multi == 2
    expected = (10 + sum(1 for t in [8.0, 9.0] if t >= rec.min_margin)) / 12
    assert rec.full_pop_retention == pytest.approx(expected)

def test_recommend_raises_without_multi_locus_tp():
    with pytest.raises(ValueError):
        cm.recommend([cm.INF, cm.INF], [1.0])          # no finite-margin TP

def test_recommend_raises_without_fp():
    with pytest.raises(ValueError):
        cm.recommend([8.0], [])                         # no FP samples
