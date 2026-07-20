#!/usr/bin/env python3
# plot_ranking.py
# Purpose: Generate publication-quality visualizations of ranked GPCR candidates.
# Inputs: Ranked candidates CSV ($1), output plot file ($2)
# Outputs: Multi-panel figure with score breakdowns (${output_file}.png, ${output_file}.svg, ${output_file}.pdf)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
import sys
import os

# Set publication-quality defaults
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2

# Color palette for score components
SCORE_COLORS = {
    'phylo_score': '#2ca02c',       # Green - phylogenetic proximity
    'purifying_score': '#1f77b4',   # Blue - purifying selection
    'positive_score': '#d62728',    # Red - positive selection
    'synteny_score': '#9467bd',     # Purple - synteny
    'expression_score': '#ff7f0e',  # Orange - expression
    'lse_divergence_score': '#8c564b'    # Brown - LSE divergence
}

# Rank-confidence tiers from scripts/rank_confidence.py: the signal bootstrap's
# P(candidate is in the top k), binned high (>=0.8) / plausible (>=0.2) / tail.
# NOT the legacy evidence-completeness grade, which is no longer displayed.
#
# `tail` is grey, not red: it is a true statistical statement (this candidate
# seldom reaches the top k), and since only k candidates fit in the top k, most
# of the list is necessarily tail. Red would paint the expected majority
# outcome as an error.
RANK_TIER_ORDER = ['high', 'plausible', 'tail']
TIER_COLORS = {
    'high': '#2ca02c',       # green  - stably in the top k
    'plausible': '#ff7f0e',  # orange - membership depends on which signals are drawn
    'tail': '#7f7f7f'        # grey   - seldom in the top k; not a failure
}
# Anything else (blank cell = no signal covers the candidate; absent column =
# rank_confidence.py did not run) is "not scored" -- never binned as tail.
TIER_NA_COLOR = '#c7c7c7'
TIER_NA_LABEL = 'not scored'


def rank_tier_series(frame):
    """Normalized rank_tier values; unrecognized/blank -> TIER_NA_LABEL."""
    raw = frame['rank_tier'].fillna('').astype(str).str.strip()
    return raw.where(raw.isin(RANK_TIER_ORDER), TIER_NA_LABEL)

# Command-line arguments
ranking_file = sys.argv[1]
output_file = sys.argv[2]

# Load ranked data. Accept the pre-rename `lse_depth_*` schema under its
# canonical `lse_divergence_*` names, announced on stderr -- otherwise plotting
# a ranked CSV written before the rename would silently drop the divergence
# component out of the stacked bar and the violin panel, which reads as "this
# axis contributed nothing" rather than "this file predates the rename".
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _rank_candidates_lib import (apply_legacy_column_aliases,  # noqa: E402
                                  getenv_renamed)

df = pd.read_csv(ranking_file)
df = apply_legacy_column_aliases(df, source=f"ranked CSV {ranking_file}")

# Check if we have the enhanced columns
has_components = 'phylo_score' in df.columns

# Create figure with multiple panels
if has_components:
    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(3, 2, height_ratios=[1.2, 1, 1], hspace=0.35, wspace=0.25)
else:
    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(1, 1)

# --- Panel A: Top Candidates with Score Breakdown (Stacked Bar) ---
ax1 = fig.add_subplot(gs[0, :] if has_components else gs[0, 0])

top_n = min(15, len(df))
top_df = df.head(top_n).copy()

if has_components:
    # Get weights from environment or use defaults
    weights = {
        'phylo_score': float(os.getenv('PHYLO_WEIGHT', 2)),
        'purifying_score': float(os.getenv('PURIFYING_WEIGHT', 1)),
        'positive_score': float(os.getenv('POSITIVE_WEIGHT', 1)),
        'synteny_score': float(os.getenv('SYNTENY_WEIGHT', 3)),
        'expression_score': float(os.getenv('EXPR_WEIGHT', 1)),
        # Honours the pre-rename LSE_DEPTH_WEIGHT (announced) so a legacy
        # environment still labels the panel with the weight actually in force.
        'lse_divergence_score': getenv_renamed('LSE_DIVERGENCE_WEIGHT', 1)
    }

    # Normalize scores for stacking (use normalized versions if available)
    score_cols = []
    for col in ['phylo_score', 'purifying_score', 'positive_score',
                'synteny_score', 'expression_score', 'lse_divergence_score']:
        norm_col = f'{col}_norm'
        if norm_col in df.columns:
            score_cols.append(norm_col)
        elif col in df.columns:
            score_cols.append(col)

    # Create stacked bar data
    x = np.arange(top_n)
    bottom = np.zeros(top_n)

    component_labels = {
        'phylo_score_norm': f'Phylogenetic (w={weights["phylo_score"]})',
        'purifying_score_norm': f'Purifying sel. (w={weights["purifying_score"]})',
        'positive_score_norm': f'Positive sel. (w={weights["positive_score"]})',
        'synteny_score_norm': f'Synteny (w={weights["synteny_score"]})',
        'expression_score_norm': f'Expression (w={weights["expression_score"]})',
        'lse_divergence_score_norm': f'LSE divergence (w={weights["lse_divergence_score"]})',
        'phylo_score': 'Phylogenetic',
        'purifying_score': 'Purifying sel.',
        'positive_score': 'Positive sel.',
        'synteny_score': 'Synteny',
        'expression_score': 'Expression',
        'lse_divergence_score': 'LSE divergence'
    }

    for col in score_cols:
        base_col = col.replace('_norm', '')
        weight = weights.get(base_col, 1)
        values = top_df[col].values * weight if col in top_df.columns else np.zeros(top_n)
        color = SCORE_COLORS.get(base_col, '#7f7f7f')
        label = component_labels.get(col, col)
        ax1.bar(x, values, bottom=bottom, label=label, color=color, edgecolor='white', linewidth=0.5)
        bottom += values

    # Add rank-confidence tier indicators
    if 'rank_tier' in top_df.columns:
        tiers = rank_tier_series(top_df)
        for i, tier in enumerate(tiers.values):
            color = TIER_COLORS.get(tier, TIER_NA_COLOR)
            ax1.plot(i, bottom[i] + 0.1, marker='o', markersize=8, color=color,
                    markeredgecolor='black', markeredgewidth=0.5)

    ax1.legend(loc='upper right', framealpha=0.9, fontsize=9)
else:
    # Simple bar plot if no component data
    sns.barplot(data=top_df, x='id', y='rank_score', ax=ax1, palette='viridis')

ax1.set_xticks(range(top_n))
ax1.set_xticklabels(top_df['id'].values, rotation=45, ha='right', fontsize=9)
ax1.set_xlabel('Candidate ID', fontsize=11)
ax1.set_ylabel('Weighted Score', fontsize=11)
ax1.set_title('A. Top Ranked GPCR Candidates (Score Breakdown)', fontsize=12, fontweight='bold', loc='left')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

if has_components:
    # --- Panel B: Score Distribution by Component ---
    ax2 = fig.add_subplot(gs[1, 0])

    # Melt data for violin plot
    score_data = []
    for col in ['phylo_score', 'purifying_score', 'positive_score',
                'synteny_score', 'expression_score', 'lse_divergence_score']:
        if col in df.columns:
            for val in df[col].values:
                score_data.append({'Component': col.replace('_score', '').replace('_', ' ').title(),
                                   'Score': val})

    if score_data:
        score_df = pd.DataFrame(score_data)
        palette = [SCORE_COLORS.get(f"{c.lower().replace(' ', '_')}_score", '#7f7f7f')
                   for c in score_df['Component'].unique()]
        sns.violinplot(data=score_df, x='Component', y='Score', ax=ax2, palette=palette, inner='box')
        ax2.set_xticklabels(ax2.get_xticklabels(), rotation=30, ha='right')

    ax2.set_xlabel('')
    ax2.set_ylabel('Raw Score', fontsize=11)
    ax2.set_title('B. Score Distribution by Component', fontsize=12, fontweight='bold', loc='left')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # --- Panel C: Rank-Confidence Tier Distribution ---
    ax3 = fig.add_subplot(gs[1, 1])

    if 'rank_tier' in df.columns:
        tiers = rank_tier_series(df)
        tier_order = RANK_TIER_ORDER + [TIER_NA_LABEL]
        tier_counts = tiers.value_counts().reindex(tier_order).fillna(0)
        colors = [TIER_COLORS.get(t, TIER_NA_COLOR) for t in tier_order]

        bars = ax3.bar(tier_order, tier_counts.values, color=colors, edgecolor='black', linewidth=0.8)

        # Add count labels
        for bar, count in zip(bars, tier_counts.values):
            ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    f'{int(count)}', ha='center', va='bottom', fontsize=10, fontweight='bold')

        topk = os.getenv('RANK_TOPK', '20')
        ax3.set_xlabel(f'Rank tier: bootstrap P(in top {topk})', fontsize=11)
        ax3.set_ylabel('Number of Candidates', fontsize=11)
        ax3.set_title('C. Rank-Confidence Tier Distribution', fontsize=12, fontweight='bold', loc='left')
    else:
        # rank_confidence.py is invoked non-fatally by stage 07, so the column
        # can legitimately be missing. Say so, rather than drawing an empty or
        # all-tail bar chart that would read as a real result.
        ax3.text(0.5, 0.5,
                 'Rank confidence not computed\n(rank_confidence.py did not run)',
                 ha='center', va='center', transform=ax3.transAxes, fontsize=10)
        ax3.set_title('C. Rank-Confidence Tier Distribution', fontsize=12, fontweight='bold', loc='left')
        ax3.set_xticks([])
        ax3.set_yticks([])

    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    # --- Panel D: Evidence Completeness vs Rank Score ---
    ax4 = fig.add_subplot(gs[2, 0])

    if 'evidence_completeness' in df.columns:
        scatter = ax4.scatter(df['evidence_completeness'], df['rank_score'],
                             c=df['rank_score'], cmap='viridis',
                             s=50, alpha=0.7, edgecolors='black', linewidths=0.5)

        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax4)
        cbar.set_label('Rank Score', fontsize=10)

        # Highlight top candidates
        top_5 = df.head(5)
        ax4.scatter(top_5['evidence_completeness'], top_5['rank_score'],
                   s=100, facecolors='none', edgecolors='red', linewidths=2, label='Top 5')

        # Add correlation
        corr = df['evidence_completeness'].corr(df['rank_score'])
        ax4.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax4.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        ax4.set_xlabel('Evidence Completeness', fontsize=11)
        ax4.set_ylabel('Rank Score', fontsize=11)
        ax4.legend(loc='lower right', fontsize=9)
    else:
        ax4.text(0.5, 0.5, 'No evidence completeness data', ha='center', va='center', transform=ax4.transAxes)

    ax4.set_title('D. Evidence Completeness vs. Rank Score', fontsize=12, fontweight='bold', loc='left')
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)

    # --- Panel E: Selection Pressure Summary ---
    ax5 = fig.add_subplot(gs[2, 1])

    if 'purifying_score' in df.columns and 'positive_score' in df.columns:
        # Categorize by selection type
        df['selection_type'] = 'Neutral'
        df.loc[df['purifying_score'] > df['purifying_score'].median(), 'selection_type'] = 'Purifying'
        df.loc[df['positive_score'] > df['positive_score'].median(), 'selection_type'] = 'Positive'

        selection_counts = df['selection_type'].value_counts()
        colors_sel = {'Purifying': '#1f77b4', 'Positive': '#d62728', 'Neutral': '#7f7f7f'}

        wedges, texts, autotexts = ax5.pie(
            selection_counts.values,
            labels=selection_counts.index,
            colors=[colors_sel.get(s, '#7f7f7f') for s in selection_counts.index],
            autopct='%1.1f%%',
            startangle=90,
            explode=[0.02] * len(selection_counts),
            wedgeprops=dict(edgecolor='white', linewidth=1.5)
        )

        # Add significant selection count if available
        if 'selection_significant' in df.columns:
            sig_count = df['selection_significant'].sum()
            ax5.text(0, -1.3, f'{int(sig_count)} candidates with significant selection (FDR < 0.05)',
                    ha='center', fontsize=10, style='italic')

        ax5.set_title('E. Selection Pressure Distribution', fontsize=12, fontweight='bold', loc='left')
    else:
        ax5.text(0.5, 0.5, 'No selection data', ha='center', va='center', transform=ax5.transAxes)
        ax5.set_title('E. Selection Pressure Distribution', fontsize=12, fontweight='bold', loc='left')

# Add overall title
fig.suptitle('GPCR Candidate Ranking Summary', fontsize=14, fontweight='bold', y=0.98)

# Save in multiple formats
plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(f"{output_file}.svg", format='svg', bbox_inches='tight', facecolor='white')
plt.savefig(f"{output_file}.pdf", format='pdf', bbox_inches='tight', facecolor='white')
plt.close()

# --- Generate Summary Statistics Table ---
summary_file = f"{output_file}_summary.csv"
summary_stats = {
    'Metric': [],
    'Value': []
}

summary_stats['Metric'].append('Total candidates')
summary_stats['Value'].append(len(df))

if 'rank_tier' in df.columns:
    tier_values = rank_tier_series(df)
    for tier in RANK_TIER_ORDER + [TIER_NA_LABEL]:
        summary_stats['Metric'].append(f'Rank tier: {tier}')
        summary_stats['Value'].append(int((tier_values == tier).sum()))
else:
    summary_stats['Metric'].append('Rank tier')
    summary_stats['Value'].append('not computed')

if 'selection_significant' in df.columns:
    summary_stats['Metric'].append('Significant selection')
    summary_stats['Value'].append(int(df['selection_significant'].sum()))

if 'synteny_score' in df.columns:
    summary_stats['Metric'].append('With synteny support')
    summary_stats['Value'].append(int((df['synteny_score'] > 0).sum()))

if 'has_expression_data' in df.columns:
    summary_stats['Metric'].append('With expression data')
    summary_stats['Value'].append(int(df['has_expression_data'].sum()))

summary_stats['Metric'].append('Mean rank score')
summary_stats['Value'].append(f"{df['rank_score'].mean():.3f}")

summary_stats['Metric'].append('Max rank score')
summary_stats['Value'].append(f"{df['rank_score'].max():.3f}")

pd.DataFrame(summary_stats).to_csv(summary_file, index=False)

print(f"Generated ranking visualizations:")
print(f"  PNG: {output_file}.png")
print(f"  SVG: {output_file}.svg")
print(f"  PDF: {output_file}.pdf")
print(f"  Summary: {summary_file}")
