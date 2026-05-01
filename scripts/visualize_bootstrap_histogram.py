#!/usr/bin/env python3
"""Bootstrap support histogram from IQ-TREE consensus tree."""

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--contree', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--theme', choices=['light', 'dark'], default='light')
    parser.add_argument('--dpi', type=int, default=300)
    args = parser.parse_args()

    dark = args.theme == 'dark'
    bg = '#000000' if dark else 'white'
    text = '#e0e0e0' if dark else 'black'
    grid_color = '#444466' if dark else '#cccccc'

    tree = Phylo.read(args.contree, 'newick')
    supports = []
    for clade in tree.find_clades():
        if not clade.is_terminal() and clade.confidence is not None:
            supports.append(clade.confidence)
    supports = np.array(supports)

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor(bg)
    ax.set_facecolor(bg)

    bins = np.arange(0, 105, 5)
    counts, edges, patches = ax.hist(supports, bins=bins, edgecolor=bg, linewidth=0.5)

    for patch, left_edge in zip(patches, edges[:-1]):
        if left_edge >= 95:
            patch.set_facecolor('#4daf4a')
        elif left_edge >= 70:
            patch.set_facecolor('#ffbf00' if dark else '#ff7f00')
        else:
            patch.set_facecolor('#e41a1c')

    ax.axvline(70, color='#aaaaaa', linestyle='--', linewidth=1.2, alpha=0.8)
    ax.axvline(95, color='#aaaaaa', linestyle='--', linewidth=1.2, alpha=0.8)

    strong = (supports >= 95).sum()
    moderate = ((supports >= 70) & (supports < 95)).sum()
    weak = (supports < 70).sum()

    ax.legend([
        plt.Rectangle((0, 0), 1, 1, fc='#4daf4a'),
        plt.Rectangle((0, 0), 1, 1, fc='#ffbf00' if dark else '#ff7f00'),
        plt.Rectangle((0, 0), 1, 1, fc='#e41a1c'),
    ], [
        f'Strong (>=95): {strong} ({strong/len(supports)*100:.1f}%)',
        f'Moderate (70-94): {moderate} ({moderate/len(supports)*100:.1f}%)',
        f'Weak (<70): {weak} ({weak/len(supports)*100:.1f}%)',
    ], fontsize=11, facecolor=bg, edgecolor=grid_color, labelcolor=text)

    ax.set_xlabel('UFBoot Support Value', fontsize=13, color=text)
    ax.set_ylabel('Number of Internal Nodes', fontsize=13, color=text)
    ax.set_title(f'Distribution of Bootstrap Support Values\n'
                 f'Consensus tree | {len(supports)} internal nodes | '
                 f'Median={np.median(supports):.0f} | Mean={supports.mean():.1f}',
                 fontsize=14, fontweight='bold', color=text)

    ax.tick_params(colors=text)
    ax.spines['bottom'].set_color(text)
    ax.spines['left'].set_color(text)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3, color=grid_color)

    plt.tight_layout()
    plt.savefig(args.output, dpi=args.dpi, bbox_inches='tight', facecolor=bg)
    plt.close()
    print(f'Saved: {args.output}')


if __name__ == '__main__':
    main()
