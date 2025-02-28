#!/usr/bin/env python3
# rank_candidates.py
# Purpose: Rank GPCR candidates using integrated data (phylogeny, selection, synteny, expression) with weighted scoring.
# Inputs: Candidate IDs ($1), expression data ($2), phylogeny directory ($3), selective pressure directory ($4), synteny IDs ($5), output CSV ($6)
# Outputs: Ranked candidates CSV with scores
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import pandas as pd
import os
import sys
from ete3 import Tree

# Command-line arguments
candidates_file = sys.argv[1]  # File with candidate GPCR IDs
expression_file = sys.argv[2]  # Expression data file
phylo_dir = sys.argv[3]        # Directory with phylogenetic trees
selective_dir = sys.argv[4]    # Directory with selective pressure results
synteny_file = sys.argv[5]     # File with synteny IDs
output_file = sys.argv[6]      # Output ranked candidates CSV

# Load candidate IDs
candidates = pd.read_csv(candidates_file, names=['id'])
# Load expression data if available, otherwise create empty DataFrame
expr_data = pd.read_csv(expression_file, names=['id', 'weight'], sep=',') if os.path.exists(expression_file) else pd.DataFrame(columns=['id', 'weight'])

# Initialize score lists
phylo_scores, dnds_scores, synteny_scores, expr_scores, lse_depths = [], [], [], [], []

# Score each candidate
for id in candidates['id']:
    # Phylogeny score: count of trees containing the ID
    phylo_score = sum(1 for f in os.listdir(phylo_dir) if os.path.isfile(f) and 'treefile' in f and id in open(os.path.join(phylo_dir, f)).read())
    # dN/dS score: 10 if dN/dS < 1 (purifying selection), else 0
    dnds = next((float(line.split()[-1]) for f in os.listdir(selective_dir) if 'codeml' in f and id in f for line in open(os.path.join(selective_dir, f)) if 'omega' in line), 0)
    dnds_score = 10 if dnds < 1 else 0
    # Synteny score: count of synteny blocks containing the ID
    synteny_score = sum(1 for line in open(synteny_file) if id in line) if os.path.exists(synteny_file) else 0
    # Expression score: sum of weighted TPM values
    expr_score = expr_data[expr_data['id'] == id]['weight'].sum() if id in expr_data['id'].values else 0
    # LSE depth score: 10 if depth > 0.5, else 0
    tree_file = next((f for f in os.listdir(phylo_dir) if id in f and f.endswith('.treefile')), None)
    lse_depth = Tree(os.path.join(phylo_dir, tree_file)).search_nodes(name=id)[0].get_distance(t) if tree_file else 0
    lse_depth_score = 10 if lse_depth > 0.5 else 0
    
    # Append scores to lists
    phylo_scores.append(phylo_score)
    dnds_scores.append(dnds_score)
    synteny_scores.append(synteny_score)
    expr_scores.append(expr_score)
    lse_depths.append(lse_depth_score)

# Create DataFrame with all scores
df = pd.DataFrame({
    'id': candidates['id'],
    'phylo_score': phylo_scores,
    'dnds_score': dnds_scores,
    'synteny_score': synteny_scores,
    'expression_score': expr_scores,
    'lse_depth_score': lse_depths
})

# Calculate overall rank score with weights
df['rank_score'] = (df['phylo_score'] * 2) + df['dnds_score'] + (df['synteny_score'] * 3) + df['expression_score'] + df['lse_depth_score']

# Sort by rank score in descending order and save to CSV
df.sort_values('rank_score', ascending=False).to_csv(output_file, index=False)
