#!/usr/bin/env python3
# plot_selective_pressure.py
# Purpose: Generate publication-quality selective pressure visualizations.
# Inputs: aBSREL results CSV ($1), output plot file ($2), [ranked candidates CSV ($3)]
# Outputs: Multi-panel figure (${output_file}.png, ${output_file}.svg, ${output_file}.pdf)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import sys
import os

# Set publication-quality defaults
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2

# Get thresholds from environment
FDR_THRESHOLD = float(os.getenv('ABSREL_FDR_THRESHOLD', 0.05))

# Command-line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]
ranked_file = sys.argv[3] if len(sys.argv) > 3 else None

# Try to load aBSREL results (new format) or LRT results (old format)
df = None
is_absrel = False

# Try reading as aBSREL CSV first
try:
    df = pd.read_csv(input_file)
    if 'omega' in df.columns or 'dnds' in df.columns or 'dN/dS' in df.columns:
        is_absrel = True
        # Standardize column names
        if 'dnds' in df.columns:
            df['omega'] = df['dnds']
        elif 'dN/dS' in df.columns:
            df['omega'] = df['dN/dS']
        if 'pvalue' in df.columns:
            df['p_value'] = df['pvalue']
        elif 'p-value' in df.columns:
            df['p_value'] = df['p-value']
except Exception:
    pass

# Fall back to old LRT format
if df is None or df.empty:
    try:
        df = pd.read_csv(input_file, names=['base', 'chi2_stat', 'p_value'])
        is_absrel = False
    except Exception as e:
        print(f"Error: Could not parse {input_file}: {e}", file=sys.stderr)
        # Create empty plot
        plt.figure(figsize=(10, 6))
        plt.text(0.5, 0.5, "No Selective Pressure Data", ha='center', va='center', fontsize=14)
        plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{output_file}.svg", format='svg', bbox_inches='tight')
        plt.close()
        sys.exit(0)

if df.empty:
    print(f"Warning: No data in {input_file}", file=sys.stderr)
    plt.figure(figsize=(10, 6))
    plt.text(0.5, 0.5, "No Selective Pressure Data", ha='center', va='center', fontsize=14)
    plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_file}.svg", format='svg', bbox_inches='tight')
    plt.close()
    sys.exit(0)

# Load ranked candidates if available (for annotation)
ranked_df = None
if ranked_file and os.path.exists(ranked_file):
    try:
        ranked_df = pd.read_csv(ranked_file)
    except Exception:
        pass

# Create figure
if is_absrel:
    fig = plt.figure(figsize=(14, 10))
    gs = gridspec.GridSpec(2, 2, hspace=0.3, wspace=0.25)
else:
    fig = plt.figure(figsize=(10, 6))
    gs = gridspec.GridSpec(1, 1)

if is_absrel:
    # Prepare data
    df['log_omega'] = np.log10(df['omega'].clip(lower=0.001, upper=100))
    df['neg_log_p'] = -np.log10(df['p_value'].clip(lower=1e-300))

    # Determine significance
    df['significant'] = df['p_value'] < FDR_THRESHOLD
    df['selection_type'] = 'Neutral'
    df.loc[(df['omega'] < 1) & df['significant'], 'selection_type'] = 'Purifying'
    df.loc[(df['omega'] > 1) & df['significant'], 'selection_type'] = 'Positive'

    # Color mapping
    colors = {'Positive': '#d62728', 'Purifying': '#1f77b4', 'Neutral': '#7f7f7f'}

    # --- Panel A: Volcano Plot (log(omega) vs -log10(p-value)) ---
    ax1 = fig.add_subplot(gs[0, 0])

    for sel_type in ['Neutral', 'Purifying', 'Positive']:
        subset = df[df['selection_type'] == sel_type]
        ax1.scatter(subset['log_omega'], subset['neg_log_p'],
                   c=colors[sel_type], label=f'{sel_type} (n={len(subset)})',
                   alpha=0.6, s=40, edgecolors='white', linewidths=0.3)

    # Add significance threshold line
    ax1.axhline(y=-np.log10(FDR_THRESHOLD), color='gray', linestyle='--', linewidth=1,
               label=f'FDR = {FDR_THRESHOLD}')

    # Add omega=1 line (neutral)
    ax1.axvline(x=0, color='gray', linestyle=':', linewidth=1, alpha=0.5)

    # Label top significant candidates
    if ranked_df is not None and 'id' in df.columns:
        top_ids = ranked_df.head(5)['id'].values
        for _, row in df[df['id'].isin(top_ids)].iterrows():
            ax1.annotate(row['id'], (row['log_omega'], row['neg_log_p']),
                        fontsize=8, alpha=0.8,
                        xytext=(5, 5), textcoords='offset points')

    ax1.set_xlabel('log₁₀(ω)', fontsize=11)
    ax1.set_ylabel('-log₁₀(p-value)', fontsize=11)
    ax1.set_title('A. Selective Pressure Volcano Plot', fontsize=12, fontweight='bold', loc='left')
    ax1.legend(loc='upper right', fontsize=9, framealpha=0.9)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # --- Panel B: Omega Distribution ---
    ax2 = fig.add_subplot(gs[0, 1])

    # Histogram of omega values
    omega_vals = df['omega'].clip(upper=10)  # Clip for visualization
    bins = np.linspace(0, 3, 31)

    ax2.hist(omega_vals[omega_vals <= 1], bins=bins, color=colors['Purifying'],
            alpha=0.7, label='ω ≤ 1 (Purifying)', edgecolor='white')
    ax2.hist(omega_vals[omega_vals > 1], bins=bins, color=colors['Positive'],
            alpha=0.7, label='ω > 1 (Positive)', edgecolor='white')

    ax2.axvline(x=1, color='black', linestyle='--', linewidth=1.5, label='ω = 1 (Neutral)')

    # Add summary stats
    median_omega = df['omega'].median()
    ax2.axvline(x=median_omega, color='green', linestyle=':', linewidth=1.5,
               label=f'Median ω = {median_omega:.2f}')

    ax2.set_xlabel('ω (dN/dS)', fontsize=11)
    ax2.set_ylabel('Count', fontsize=11)
    ax2.set_title('B. Distribution of Selection Coefficients', fontsize=12, fontweight='bold', loc='left')
    ax2.legend(loc='upper right', fontsize=9, framealpha=0.9)
    ax2.set_xlim(0, 3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # --- Panel C: Selection Type Breakdown ---
    ax3 = fig.add_subplot(gs[1, 0])

    selection_counts = df['selection_type'].value_counts()
    sel_order = ['Purifying', 'Neutral', 'Positive']
    selection_counts = selection_counts.reindex(sel_order).fillna(0)

    bars = ax3.bar(sel_order, selection_counts.values,
                  color=[colors[s] for s in sel_order],
                  edgecolor='black', linewidth=0.8)

    # Add percentage labels
    total = selection_counts.sum()
    for bar, count in zip(bars, selection_counts.values):
        pct = count / total * 100 if total > 0 else 0
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{int(count)}\n({pct:.1f}%)', ha='center', va='bottom', fontsize=10)

    ax3.set_xlabel('Selection Type', fontsize=11)
    ax3.set_ylabel('Number of Branches', fontsize=11)
    ax3.set_title('C. Selection Classification Summary', fontsize=12, fontweight='bold', loc='left')
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    # --- Panel D: P-value Distribution ---
    ax4 = fig.add_subplot(gs[1, 1])

    # Q-Q style: expected vs observed p-values
    sorted_p = np.sort(df['p_value'].values)
    expected_p = np.arange(1, len(sorted_p) + 1) / len(sorted_p)

    ax4.scatter(expected_p, sorted_p, alpha=0.5, s=20, c='steelblue', edgecolors='none')
    ax4.plot([0, 1], [0, 1], 'k--', linewidth=1, label='Expected (null)')

    # Highlight deviation from null
    ax4.fill_between([0, 1], [0, 1], [1, 1], alpha=0.1, color='red',
                     label='Enrichment region')

    ax4.set_xlabel('Expected p-value', fontsize=11)
    ax4.set_ylabel('Observed p-value', fontsize=11)
    ax4.set_title('D. P-value Calibration (Q-Q Plot)', fontsize=12, fontweight='bold', loc='left')
    ax4.legend(loc='lower right', fontsize=9)
    ax4.set_xlim(0, 1)
    ax4.set_ylim(0, 1)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)

else:
    # Old LRT format - simple scatter plot
    ax = fig.add_subplot(gs[0, 0])

    # Calculate -log10(p-value) safely
    neg_log_p = df['p_value'].apply(lambda x: -np.log10(x) if x > 0 else 0)

    scatter = ax.scatter(df['chi2_stat'], neg_log_p, c=neg_log_p, cmap='viridis',
                        alpha=0.7, s=50, edgecolors='white', linewidths=0.5)

    plt.colorbar(scatter, ax=ax, label='-log₁₀(p-value)')

    # Add significance threshold
    ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', linewidth=1,
              label='p = 0.05')

    ax.set_xlabel('Chi-squared Statistic', fontsize=11)
    ax.set_ylabel('-log₁₀(p-value)', fontsize=11)
    ax.set_title('Selective Pressure Analysis (LRT)', fontsize=12, fontweight='bold')
    ax.legend(loc='upper right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# Add overall title
fig.suptitle('Selective Pressure Analysis', fontsize=14, fontweight='bold', y=0.98)

# Save in multiple formats
plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(f"{output_file}.svg", format='svg', bbox_inches='tight', facecolor='white')
plt.savefig(f"{output_file}.pdf", format='pdf', bbox_inches='tight', facecolor='white')
plt.close()

# --- Generate Summary Statistics ---
summary_file = f"{output_file}_summary.txt"
with open(summary_file, 'w') as f:
    f.write("Selective Pressure Analysis Summary\n")
    f.write("=" * 40 + "\n\n")

    f.write(f"Total branches analyzed: {len(df)}\n\n")

    if is_absrel:
        f.write("Selection Classification:\n")
        for sel_type in ['Purifying', 'Neutral', 'Positive']:
            count = len(df[df['selection_type'] == sel_type])
            pct = count / len(df) * 100 if len(df) > 0 else 0
            f.write(f"  {sel_type}: {count} ({pct:.1f}%)\n")

        f.write(f"\nOmega (dN/dS) Statistics:\n")
        f.write(f"  Mean: {df['omega'].mean():.4f}\n")
        f.write(f"  Median: {df['omega'].median():.4f}\n")
        f.write(f"  Min: {df['omega'].min():.4f}\n")
        f.write(f"  Max: {df['omega'].max():.4f}\n")

        sig_count = df['significant'].sum()
        f.write(f"\nSignificant branches (p < {FDR_THRESHOLD}): {sig_count} ({sig_count/len(df)*100:.1f}%)\n")

print(f"Generated selective pressure visualizations:")
print(f"  PNG: {output_file}.png")
print(f"  SVG: {output_file}.svg")
print(f"  PDF: {output_file}.pdf")
print(f"  Summary: {summary_file}")
