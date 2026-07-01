"""Audit statistical redundancy among the per-signal ranking scores.

The candidate ranking combines ~11 signals, several of which the code derives
from the SAME OrthoFinder orthogroup / gene tree (phylo distance, OG confidence,
and the dN/dS OG-reliability weight). Rank aggregation must not let one shared
confound count as several independent votes. This emits the correlation matrix
plus a grouping consumed downstream.
"""
from __future__ import annotations
import json, sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SIGNAL_COLUMNS = [
    "phylo_score", "purifying_score", "positive_score", "synteny_score",
    "expression_score", "lse_depth_score", "chemosensory_expr_score",
    "gprotein_coexpr_score", "ecl_divergence_score", "expansion_score",
    "og_confidence_score",
    "tandem_cluster_score",  # flag suffix-swaps correctly (has_tandem_cluster_data)
    # Task 6: structural/embedding/microswitch evidence channels
    # (rank_aggregation.py SIGNAL_SPEC). None of these follow the
    # "<name>_score" naming convention, so each needs a FLAG_OVERRIDES entry.
    "struct_novelty", "struct_nonchemo_corrob",
    "emb_classA_sim", "emb_nonchemo_sim",
    "or_microswitch",
]

# Signals whose has_*_data flag does not follow the "<name>_score" ->
# "has_<name>_data" suffix-swap. The production ranked CSV names these two
# flags after the shorter data-source key, not the full score column.
FLAG_OVERRIDES = {
    "gprotein_coexpr_score": "has_gprotein_data",
    "ecl_divergence_score": "has_ecl_data",
    "struct_novelty": "has_struct_data",
    "struct_nonchemo_corrob": "has_struct_data",
    "emb_classA_sim": "has_emb_data",
    "emb_nonchemo_sim": "has_emb_data",
    "or_microswitch": "has_or_microswitch_data",
}

def load_signal_matrix(csv_path):
    df = pd.read_csv(csv_path)
    cols = [c for c in SIGNAL_COLUMNS if c in df.columns]
    for c in SIGNAL_COLUMNS:
        if c not in df.columns:
            print(f"[audit] warning: signal column {c} absent - skipped", file=sys.stderr)
    m = df[cols].astype(float).copy()
    for c in cols:
        flag = FLAG_OVERRIDES.get(c, "has_" + c.replace("_score", "_data"))
        if flag in df.columns:
            m.loc[~df[flag].astype(bool), c] = np.nan
    if "id" in df.columns:
        m.index = df["id"].values
    return m[sorted(m.columns)]

def rank_correlation(matrix):
    cols = list(matrix.columns)
    out = pd.DataFrame(np.eye(len(cols)), index=cols, columns=cols)
    for i, a in enumerate(cols):
        for b in cols[i + 1:]:
            pair = matrix[[a, b]].dropna()
            rho = spearmanr(pair[a], pair[b]).correlation if len(pair) >= 3 else np.nan
            rho = 0.0 if rho is None or (isinstance(rho, float) and np.isnan(rho)) else float(rho)
            out.loc[a, b] = out.loc[b, a] = rho
    return out

def group_correlated_signals(corr, threshold=0.7):
    cols = list(corr.columns)
    parent = {c: c for c in cols}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for i, a in enumerate(cols):
        for b in cols[i + 1:]:
            if abs(corr.loc[a, b]) >= threshold:
                parent[find(a)] = find(b)
    groups = {}
    for c in cols:
        groups.setdefault(find(c), []).append(c)
    return [sorted(g) for g in sorted(groups.values(), key=lambda g: sorted(g)[0])]

def write_report(corr, groups, out_prefix, threshold=0.7):
    corr.to_csv(out_prefix + ".tsv", sep="\t")
    with open(out_prefix + "_groups.json", "w") as fh:
        json.dump({"groups": groups, "threshold": threshold}, fh, indent=2)
    lines = ["# Signal-independence audit", "",
             f"Grouping threshold |rho| >= {threshold}", "", "## Groups"]
    for g in groups:
        note = " (REDUNDANT - counts once in aggregation)" if len(g) > 1 else ""
        lines.append(f"- {', '.join(g)}{note}")
    with open(out_prefix + ".md", "w") as fh:
        fh.write("\n".join(lines) + "\n")

def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--ranked-csv", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--threshold", type=float, default=0.7)
    a = ap.parse_args(argv)
    m = load_signal_matrix(a.ranked_csv)
    corr = rank_correlation(m)
    groups = group_correlated_signals(corr, a.threshold)
    write_report(corr, groups, a.out_prefix, a.threshold)
    print(f"[audit] {len(m.columns)} signals -> {len(groups)} groups")

if __name__ == "__main__":
    main()
