#!/usr/bin/env python3
# rank_candidates.py
# Purpose: Rank GPCR candidates using phylogenetic proximity, dN/dS, expression, and synteny.
# Inputs: Candidate IDs ($1), expression data ($2), phylogeny dir ($3), selective pressure dir ($4), synteny IDs ($5), output CSV ($6)
# Outputs: Ranked candidates CSV
# Author: Jorge L. Perez-Moreno, Ph.D.

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
expr_data = pd.read_csv(expression_file, names=['id', 'weight']) if os.path.exists(expression_file) else pd.DataFrame(columns=['id', 'weight'])
synteny_ids = pd.read_csv(synteny_file, names=['id'])['id'].tolist() if os.path.exists(synteny_file) else []

t = Tree(f"{phylo_dir}/all_berghia_refs.treefile")
ref_ids = [leaf.name for leaf in t if leaf.name.startswith('ref_')]

# Load aBSREL dN/dS results
dnds_data = {}
absrel_files = ['absrel_results.csv', 'absrel_results_lse.csv']
for absrel_file in absrel_files:
    absrel_path = os.path.join(selective_dir, absrel_file)
    if os.path.exists(absrel_path):
        try:
            absrel_df = pd.read_csv(absrel_path)
            # Expected columns: branch_id, omega, p_value, corrected_p_value
            for _, row in absrel_df.iterrows():
                branch_id = row.get('branch_id', row.get('id', ''))
                omega = row.get('omega', row.get('dnds', row.get('dN/dS', None)))
                if pd.notna(omega) and branch_id:
                    # Store omega value; for branches under positive selection (p < 0.05),
                    # omega > 1 indicates diversifying selection
                    dnds_data[branch_id] = float(omega)
        except Exception as e:
            print(f"Warning: Could not parse {absrel_path}: {e}", file=sys.stderr)

def min_distance_to_refs(node_name):
    if node_name in ref_ids:
        return 0
    distances = [t.get_distance(node_name, ref) for ref in ref_ids]
    return min(distances) if distances else float('inf')

def get_dnds_score(candidate_id, dnds_data):
    """
    Calculate dN/dS score for ranking.

    Biological interpretation:
    - dN/dS < 1: purifying selection (conserved function)
    - dN/dS = 1: neutral evolution
    - dN/dS > 1: positive/diversifying selection (potential functional divergence)

    For chemoreceptor discovery, we want candidates showing evidence of:
    1. Strong purifying selection (conserved chemoreceptor function) OR
    2. Positive selection (potential novel/specialized function)

    Score: distance from neutrality (|omega - 1|) weighted by direction
    """
    if candidate_id in dnds_data:
        omega = dnds_data[candidate_id]
        # Score candidates far from neutral evolution higher
        # Both strong purifying (omega << 1) and positive selection (omega > 1) are interesting
        if omega > 1:
            # Positive selection - potentially diversifying/novel function
            return omega  # Higher omega = higher score
        else:
            # Purifying selection - conserved function
            # Transform so strong purifying selection scores well: 1/(omega + 0.01)
            return 1 / (omega + 0.01)
    # Default: neutral assumption (no data available)
    return 1.0

phylo_scores, dnds_scores, synteny_scores, expr_scores, lse_depths = [], [], [], [], []
for id in candidates['id']:
    phylo_score = 1 / (min_distance_to_refs(id) + 1e-5)
    dnds_score = get_dnds_score(id, dnds_data)
    synteny_score = 1 if id in synteny_ids else 0
    expr_score = expr_data[expr_data['id'] == id]['weight'].sum() if id in expr_data['id'].values else 0
    lse_depth = t.search_nodes(name=id)[0].get_distance(t) if id in t else 0
    lse_depth_score = lse_depth if lse_depth > 0.5 else 0

    phylo_scores.append(phylo_score)
    dnds_scores.append(dnds_score)
    synteny_scores.append(synteny_score)
    expr_scores.append(expr_score)
    lse_depths.append(lse_depth_score)

df = pd.DataFrame({
    'id': candidates['id'], 'phylo_score': phylo_scores, 'dnds_score': dnds_scores,
    'synteny_score': synteny_scores, 'expression_score': expr_scores, 'lse_depth_score': lse_depths
})

for col in ['phylo_score', 'dnds_score', 'expression_score', 'lse_depth_score']:
    df[col] = (df[col] - df[col].min()) / (df[col].max() - df[col].min()) if df[col].max() > df[col].min() else 0

phylo_weight, dnds_weight, synteny_weight, expr_weight, lse_depth_weight = map(float, [
    os.getenv('PHYLO_WEIGHT', 2), os.getenv('DNDS_WEIGHT', 1), os.getenv('SYNTENY_WEIGHT', 3),
    os.getenv('EXPR_WEIGHT', 1), os.getenv('LSE_DEPTH_WEIGHT', 1)
])

df['rank_score'] = (df['phylo_score'] * phylo_weight) + (df['dnds_score'] * dnds_weight) + \
                   (df['synteny_score'] * synteny_weight) + (df['expression_score'] * expr_weight) + \
                   (df['lse_depth_score'] * lse_depth_weight)

df.sort_values('rank_score', ascending=False).to_csv(output_file, index=False)
