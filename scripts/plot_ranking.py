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
    'lse_depth_score': '#8c564b'    # Brown - LSE depth
}

TIER_COLORS = {
    'High': '#2ca02c',
    'Medium': '#ff7f0e',
    'Low': '#d62728'
}

# Command-line arguments
ranking_file = sys.argv[1]
output_file = sys.argv[2]

# Load ranked data
df = pd.read_csv(ranking_file)

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
        'lse_depth_score': float(os.getenv('LSE_DEPTH_WEIGHT', 1))
    }

    # Normalize scores for stacking (use normalized versions if available)
    score_cols = []
    for col in ['phylo_score', 'purifying_score', 'positive_score',
                'synteny_score', 'expression_score', 'lse_depth_score']:
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
        'lse_depth_score_norm': f'LSE depth (w={weights["lse_depth_score"]})',
        'phylo_score': 'Phylogenetic',
        'purifying_score': 'Purifying sel.',
        'positive_score': 'Positive sel.',
        'synteny_score': 'Synteny',
        'expression_score': 'Expression',
        'lse_depth_score': 'LSE depth'
    }

    for col in score_cols:
        base_col = col.replace('_norm', '')
        weight = weights.get(base_col, 1)
        values = top_df[col].values * weight if col in top_df.columns else np.zeros(top_n)
        color = SCORE_COLORS.get(base_col, '#7f7f7f')
        label = component_labels.get(col, col)
        ax1.bar(x, values, bottom=bottom, label=label, color=color, edgecolor='white', linewidth=0.5)
        bottom += values

    # Add confidence tier indicators
    if 'confidence_tier' in top_df.columns:
        for i, (_, row) in enumerate(top_df.iterrows()):
            tier = row.get('confidence_tier', 'Low')
            color = TIER_COLORS.get(tier, '#7f7f7f')
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
                'synteny_score', 'expression_score', 'lse_depth_score']:
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

    # --- Panel C: Confidence Tier Distribution ---
    ax3 = fig.add_subplot(gs[1, 1])

    if 'confidence_tier' in df.columns:
        tier_counts = df['confidence_tier'].value_counts()
        tier_order = ['High', 'Medium', 'Low']
        tier_counts = tier_counts.reindex(tier_order).fillna(0)
        colors = [TIER_COLORS.get(t, '#7f7f7f') for t in tier_order]

        bars = ax3.bar(tier_order, tier_counts.values, color=colors, edgecolor='black', linewidth=0.8)

        # Add count labels
        for bar, count in zip(bars, tier_counts.values):
            ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    f'{int(count)}', ha='center', va='bottom', fontsize=10, fontweight='bold')

        ax3.set_xlabel('Confidence Tier', fontsize=11)
        ax3.set_ylabel('Number of Candidates', fontsize=11)
        ax3.set_title('C. Candidate Confidence Distribution', fontsize=12, fontweight='bold', loc='left')
    else:
        ax3.text(0.5, 0.5, 'No confidence tier data', ha='center', va='center', transform=ax3.transAxes)

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

if 'confidence_tier' in df.columns:
    for tier in ['High', 'Medium', 'Low']:
        count = len(df[df['confidence_tier'] == tier])
        summary_stats['Metric'].append(f'{tier} confidence')
        summary_stats['Value'].append(count)

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
