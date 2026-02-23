#!/usr/bin/env python3
# plot_synteny.py
# Purpose: Generate multi-level synteny plots using MCScanX output.
# Inputs: Synteny dir ($1), output prefix ($2)
# Outputs: Synteny plots in ${output_prefix}_${level}.png for each taxonomic level
# Logic: Plots synteny blocks at Aeolids, Nudibranchs, Gastropods levels based on available genomes.
# Author: Jorge L. Perez-Moreno, Ph.D.

import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import re
from collections import defaultdict

synteny_dir = sys.argv[1]
output_prefix = sys.argv[2]

# Taxonomic levels from config.sh (hardcoded for simplicity, could parse config.sh)
levels = ["Aeolids", "Nudibranchs", "Gastropods"]


def parse_collinearity_file(filepath):
    """
    Parse MCScanX .collinearity file format.

    Format example:
    ## Alignment 0: score=1234.0 e_value=1e-100 N=10 chr1&chr2 plus
    0-  0:	gene1	gene2
    0-  1:	gene3	gene4
    ...
    ## Alignment 1: score=567.0 e_value=1e-50 N=5 chr3&chr4 minus
    1-  0:	gene5	gene6
    ...

    Returns list of synteny blocks, each with metadata and gene pairs.
    """
    blocks = []
    current_block = None

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()

            # Skip empty lines
            if not line:
                continue

            # Parse alignment header
            if line.startswith('## Alignment'):
                # Save previous block if exists
                if current_block and current_block['pairs']:
                    blocks.append(current_block)

                # Parse header: ## Alignment N: score=X e_value=Y N=Z chr1&chr2 orientation
                # Note: Using [\d.eE+\-]+ for scientific notation (e.g., 1.5e-10)
                header_match = re.match(
                    r'## Alignment (\d+): score=([\d.eE+\-]+) e_value=([\d.eE+\-]+) N=(\d+)\s+(\S+)&(\S+)\s+(\S+)',
                    line
                )
                if header_match:
                    current_block = {
                        'id': int(header_match.group(1)),
                        'score': float(header_match.group(2)),
                        'e_value': float(header_match.group(3)),
                        'n_genes': int(header_match.group(4)),
                        'chr1': header_match.group(5),
                        'chr2': header_match.group(6),
                        'orientation': header_match.group(7),
                        'pairs': []
                    }
                else:
                    # Fallback for simpler header format
                    current_block = {
                        'id': len(blocks),
                        'score': 0,
                        'e_value': 1,
                        'n_genes': 0,
                        'chr1': 'unknown',
                        'chr2': 'unknown',
                        'orientation': 'plus',
                        'pairs': []
                    }

            # Parse gene pair lines
            elif current_block is not None and not line.startswith('#'):
                # Format: N-  M:	gene1	gene2
                # or simpler: gene1	gene2
                parts = line.split()

                # Try to extract gene pair
                if len(parts) >= 2:
                    # Check if first part contains alignment index (e.g., "0-  0:")
                    if ':' in parts[0]:
                        # Format with index: "0-  0:\tgene1\tgene2" -> parts might be split
                        # Find the gene names (usually last two non-empty fields)
                        gene_parts = [p for p in parts if p and not p.endswith(':')]
                        if len(gene_parts) >= 2:
                            gene1, gene2 = gene_parts[-2], gene_parts[-1]
                            current_block['pairs'].append((gene1, gene2))
                    else:
                        # Simple format: gene1 gene2
                        gene1, gene2 = parts[0], parts[1]
                        current_block['pairs'].append((gene1, gene2))

    # Don't forget the last block
    if current_block and current_block['pairs']:
        blocks.append(current_block)

    return blocks


def get_gene_positions(blocks):
    """
    Assign positions to genes for plotting.
    Returns dict mapping gene names to their chromosome and position.
    """
    chr_genes = defaultdict(list)

    for block in blocks:
        for gene1, gene2 in block['pairs']:
            if gene1 not in chr_genes[block['chr1']]:
                chr_genes[block['chr1']].append(gene1)
            if gene2 not in chr_genes[block['chr2']]:
                chr_genes[block['chr2']].append(gene2)

    # Assign positions based on order
    gene_positions = {}
    for chrom, genes in chr_genes.items():
        for i, gene in enumerate(genes):
            gene_positions[gene] = {'chr': chrom, 'pos': i}

    return gene_positions, chr_genes


def plot_synteny_blocks(blocks, output_file, title="Synteny Plot"):
    """
    Create a synteny dot plot from parsed blocks.
    """
    if not blocks:
        # Create empty plot with message
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, 'No synteny blocks found', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title)
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        return

    gene_positions, chr_genes = get_gene_positions(blocks)

    # Get unique chromosomes
    all_chroms = list(chr_genes.keys())
    if len(all_chroms) < 2:
        # Not enough chromosomes for comparison
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, 'Insufficient data for synteny comparison',
                ha='center', va='center', transform=ax.transAxes, fontsize=14)
        ax.set_title(title)
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        return

    fig, ax = plt.subplots(figsize=(12, 8))

    # Color palette for different blocks
    colors = plt.cm.tab20(range(min(len(blocks), 20)))

    # Plot each synteny block
    for i, block in enumerate(blocks):
        color = colors[i % len(colors)]

        for gene1, gene2 in block['pairs']:
            if gene1 in gene_positions and gene2 in gene_positions:
                pos1 = gene_positions[gene1]['pos']
                pos2 = gene_positions[gene2]['pos']

                # Plot connection line
                ax.plot([pos1, pos2], [0, 1], color=color, alpha=0.6, linewidth=0.8)

                # Plot gene points
                ax.scatter(pos1, 0, color=color, s=20, zorder=5)
                ax.scatter(pos2, 1, color=color, s=20, zorder=5)

    # Customize plot
    ax.set_ylim(-0.1, 1.1)
    ax.set_yticks([0, 1])

    # Label chromosomes
    if len(all_chroms) >= 2:
        ax.set_yticklabels([all_chroms[0], all_chroms[1] if len(all_chroms) > 1 else 'Genome 2'])

    ax.set_xlabel("Gene Position")
    ax.set_ylabel("Genome")
    ax.set_title(title)

    # Add legend for block count
    legend_text = f"{len(blocks)} synteny blocks, {sum(len(b['pairs']) for b in blocks)} gene pairs"
    ax.text(0.02, 0.98, legend_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()


# Process each taxonomic level
for level in levels:
    all_blocks = []

    # Find collinearity files for this level
    for file in os.listdir(synteny_dir):
        if file.endswith('.collinearity'):
            filepath = os.path.join(synteny_dir, file)
            try:
                blocks = parse_collinearity_file(filepath)
                all_blocks.extend(blocks)
            except Exception as e:
                print(f"Warning: Could not parse {file}: {e}", file=sys.stderr)

    # Create plot for this level
    output_file = f"{output_prefix}_{level}.png"
    plot_synteny_blocks(all_blocks, output_file, title=f"Synteny at {level} Level")
    print(f"Created {output_file} with {len(all_blocks)} blocks")
