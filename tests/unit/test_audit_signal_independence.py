import json
import numpy as np
import pandas as pd
import audit_signal_ranking_independence as au

def _csv(tmp_path):
    df = pd.DataFrame({
        "id": [f"g{i}" for i in range(6)],
        "phylo_score":      [1, 2, 3, 4, 5, 6],
        "has_phylo_data":   [True]*6,
        "og_confidence_score":   [10, 20, 30, 40, 50, 60],   # monotone in phylo -> rho=1
        "has_og_confidence_data":[True]*6,
        "synteny_score":    [6, 1, 5, 2, 4, 3],              # scrambled -> independent
        "has_synteny_data": [True]*6,
    })
    p = tmp_path / "ranked.csv"; df.to_csv(p, index=False); return str(p)

def test_matrix_shape_and_alignment(tmp_path):
    m = au.load_signal_matrix(_csv(tmp_path))
    assert "phylo_score" in m.columns and "og_confidence_score" in m.columns
    assert len(m) == 6

def test_rank_correlation_detects_redundancy(tmp_path):
    m = au.load_signal_matrix(_csv(tmp_path))
    corr = au.rank_correlation(m)
    assert corr.loc["phylo_score", "og_confidence_score"] > 0.99
    assert abs(corr.loc["phylo_score", "synteny_score"]) < 0.9

def test_grouping_clusters_redundant_signals(tmp_path):
    m = au.load_signal_matrix(_csv(tmp_path))
    groups = au.group_correlated_signals(au.rank_correlation(m), threshold=0.7)
    assert any({"phylo_score", "og_confidence_score"} <= set(g) for g in groups)
    assert any(g == ["synteny_score"] for g in groups)

def test_write_report_emits_groups_json(tmp_path):
    m = au.load_signal_matrix(_csv(tmp_path))
    corr = au.rank_correlation(m); groups = au.group_correlated_signals(corr, 0.7)
    au.write_report(corr, groups, str(tmp_path / "sig"))
    data = json.loads((tmp_path / "sig_groups.json").read_text())
    assert data["threshold"] == 0.7 and isinstance(data["groups"], list)

def test_missing_signal_column_is_skipped(tmp_path):
    # a CSV with only two signal columns must not crash
    df = pd.DataFrame({"id": ["a","b","c"], "phylo_score": [1,2,3],
                       "has_phylo_data": [True]*3})
    p = tmp_path / "sparse.csv"; df.to_csv(p, index=False)
    m = au.load_signal_matrix(str(p))
    assert "phylo_score" in m.columns

def test_production_flag_names_mask_gprotein_and_ecl(tmp_path):
    # regression: gprotein_coexpr_score / ecl_divergence_score pair with the
    # REAL pipeline flags has_gprotein_data / has_ecl_data (NOT the naive
    # has_<name>_data suffix-swap). Rows flagged False must become NaN, else
    # their 0.0 no-data sentinels leak in as real scores.
    df = pd.DataFrame({
        "id": ["a", "b", "c", "d"],
        "gprotein_coexpr_score": [0.0, 0.5, 0.0, 0.9],
        "has_gprotein_data":     [False, True, False, True],
        "ecl_divergence_score":  [0.0, 0.0, 0.7, 0.3],
        "has_ecl_data":          [False, False, True, True],
    })
    p = tmp_path / "prod.csv"; df.to_csv(p, index=False)
    m = au.load_signal_matrix(str(p))
    g = m["gprotein_coexpr_score"]
    e = m["ecl_divergence_score"]
    assert np.isnan(g.loc["a"]) and np.isnan(g.loc["c"])   # flagged False -> NaN
    assert g.loc["b"] == 0.5 and g.loc["d"] == 0.9          # flagged True -> kept
    assert np.isnan(e.loc["a"]) and np.isnan(e.loc["b"])   # flagged False -> NaN
    assert e.loc["c"] == 0.7 and e.loc["d"] == 0.3          # flagged True -> kept

def test_tandem_cluster_signal_loaded_and_masked(tmp_path):
    # the 12th live scorer signal (bead -ar8, weight 2.5) must be audited too.
    # flag has_tandem_cluster_data derives via the plain suffix-swap.
    df = pd.DataFrame({
        "id": ["a", "b", "c", "d"],
        "tandem_cluster_score":    [0.0, 0.8, 0.0, 0.4],
        "has_tandem_cluster_data": [False, True, False, True],
    })
    p = tmp_path / "tandem.csv"; df.to_csv(p, index=False)
    m = au.load_signal_matrix(str(p))
    assert "tandem_cluster_score" in m.columns
    t = m["tandem_cluster_score"]
    assert np.isnan(t.loc["a"]) and np.isnan(t.loc["c"])   # flagged False -> NaN
    assert t.loc["b"] == 0.8 and t.loc["d"] == 0.4          # flagged True -> kept
