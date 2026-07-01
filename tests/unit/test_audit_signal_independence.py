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
