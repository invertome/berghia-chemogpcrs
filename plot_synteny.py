#!/usr/bin/env python3
# plot_synteny.py
# Purpose: Generate multi-level synteny plots using MCScanX output.
# Inputs: Synteny dir ($1), output prefix ($2)
# Outputs: Synteny plots in ${output_prefix}_${level}.png for each taxonomic level
# Logic: Plots synteny blocks at Aeolids, Nudibranchs, Gastropods levels based on available genomes.
# Author: Jorge L. Perez-Moreno, Ph.D.

import sys
import matplotlib.pyplot as plt
import os

synteny_dir = sys.argv[1]
output_prefix = sys.argv[2]

# Taxonomic levels from config.sh (hardcoded for simplicity, could parse config.sh)
levels = ["Aeolids", "Nudibranchs", "Gastropods"]

for level in levels:
    plt.figure(figsize=(10, 6))
    plt.title(f"Synteny at {level} Level")
    # Placeholder for actual synteny plotting (requires parsing MCScanX output)
    for file in os.listdir(synteny_dir):
        if file.endswith('.collinearity'):
            with open(os.path.join(synteny_dir, file), 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        parts = line.split()
                        # Example plotting (simplified)
                        plt.plot([int(parts[2]), int(parts[5])], [int(parts[3]), int(parts[6])], 'b-')
    plt.xlabel("Genomic Position (Genome 1)")
    plt.ylabel("Genomic Position (Genome 2)")
    plt.savefig(f"{output_prefix}_{level}.png", dpi=300)
    plt.close()
