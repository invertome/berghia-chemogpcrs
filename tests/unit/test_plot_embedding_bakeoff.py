"""Light tests for the bake-off ranking figure builder.

Per the diagnostics-suite design, plotting tests assert structure (right models,
file written), not pixels. The scientific claims are tested in the metric
modules; here we only guard that the figure faithfully renders the tidy table.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import pandas as pd
import pytest

from plot_embedding_bakeoff import PROD_CONFIG, ranking_frame, ranking_figure


def _metrics():
    # two models x two refsets, production config + a decoy cosine row that must
    # be ignored so the figure reflects the maha/raw/multi production scorer only.
    rows = []
    for model, full, ver in [("protrek", 0.794, 0.791), ("esmc600m", 0.632, 0.645)]:
        for refset, fa in [("FULL", full), ("VERIFIED", ver)]:
            rows.append(dict(model=model, refset=refset, n_refs=1094, n_families=11,
                             dim=1024, score="maha", calib="raw", cent="multi",
                             lofo_auroc=0.99, invariance=1.0, loo_famacc=fa))
            rows.append(dict(model=model, refset=refset, n_refs=1094, n_families=11,
                             dim=1024, score="cos", calib="raw", cent="single",
                             lofo_auroc=0.9, invariance=0.8, loo_famacc=0.5))
    return pd.DataFrame(rows)


def test_ranking_frame_selects_production_config_only():
    rf = ranking_frame(_metrics())
    assert set(rf["score"]) == {PROD_CONFIG[0]}       # only maha rows survive
    assert len(rf) == 4                               # 2 models x 2 refsets, no cosine


def test_ranking_frame_orders_models_by_famacc_descending():
    rf = ranking_frame(_metrics())
    # order derived from the ranking metric (VERIFIED famAcc, fallback FULL)
    assert list(rf.drop_duplicates("model")["model"]) == ["protrek", "esmc600m"]


def test_ranking_figure_has_one_row_per_model():
    fig = ranking_figure(_metrics())
    ax = fig.axes[0]
    labels = [t.get_text() for t in ax.get_yticklabels() if t.get_text()]
    assert set(labels) == {"protrek", "esmc600m"}


def test_main_writes_outputs(tmp_path):
    from plot_embedding_bakeoff import main
    tsv = tmp_path / "m.tsv"
    _metrics().to_csv(tsv, sep="\t", index=False)
    out = tmp_path / "fig_rank"
    main([str(tsv), "--out", str(out)])
    assert (tmp_path / "fig_rank.pdf").exists()
    assert (tmp_path / "fig_rank.png").exists()


def test_empty_metrics_raises_clearly(tmp_path):
    tsv = tmp_path / "empty.tsv"
    pd.DataFrame(columns=_metrics().columns).to_csv(tsv, sep="\t", index=False)
    with pytest.raises(SystemExit):
        from plot_embedding_bakeoff import main
        main([str(tsv), "--out", str(tmp_path / "x")])
