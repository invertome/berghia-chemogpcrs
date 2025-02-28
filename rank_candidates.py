#!/usr/bin/env python3
# rank_candidates.py
# Purpose: Rank GPCR candidates using integrated data and expression weights.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import pandas as pd
import os
import sys
from ete3 import Tree

candidates_file = sys.argv[1]
expression_file = sys.argv[2]
phylo_dir = sys.argv[3]
selective_dir = sys.argv[4]
synteny_file = sys.argv[5]
output_file = sys.argv[6]

candidates = pd.read_csv(candidates_file, names=['id'])
expr_data = pd.read_csv(expression_file, names=['id', 'weight'], sep=',') if os.path.exists(expression_file) else pd.DataFrame(columns=['id', 'weight'])

phylo_scores = []
dnds_scores = []
synteny_scores = []
expr_scores = []
lse_depths = []

for id in candidates['id']:
    phylo_score = sum(1 for f in os.listdir(phylo_dir) if os.path.isfile(f) and 'treefile' in f and id in open(os.path.join(phylo_dir, f)).read())
    dnds = next((float(line.split()[-1]) for f in os.listdir(selective_dir) if 'codeml' in f and id in f for line in open(os.path.join(selective_dir, f)) if 'omega' in line), 0)
    synteny_score = sum(1 for line in open(synteny_file) if id in line) if os.path.exists(synteny_file) else 0
    expr_score = expr_data[expr_data['id'] == id]['weight'].sum() if id in expr_data['id'].values else 0
    tree_file = next((f for f in os.listdir(phylo_dir) if id in f and f.endswith('.treefile')), None)
    lse_depth = Tree(os.path.join(phylo_dir, tree_file)).search_nodes(name=id)[0].get_distance(t) if tree_file else 0
    phylo_scores.append(phylo_score)
    dnds_scores.append(dnds)
    synteny_scores.append(synteny_score)
    expr_scores.append(expr_score)
    lse_depths.append(lse_depth)

df = pd.DataFrame({
    'id': candidates['id'],
    'phylo_score': phylo_scores,
    'dnds': dnds_scores,
    'synteny_score': synteny_scores,
    'expression_score': expr_scores,
    'lse_depth': lse_depths
})
df['rank_score'] = df['phylo_score'] * 2 + df['dnds'].apply(lambda x: 10 if x < 1 else 0) + df['synteny_score'] * 3 + df['expression_score'] + df['lse_depth'].apply(lambda x: 10 if x > 0.5 else 0)
df.sort_values('rank_score', ascending=False).to_csv(output_file, index=False)
