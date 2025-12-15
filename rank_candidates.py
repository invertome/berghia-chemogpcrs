#!/usr/bin/env python3
# rank_candidates.py
# Purpose: Rank GPCR candidates using phylogenetic proximity, dN/dS, expression, and synteny.
# Inputs: Candidate IDs ($1), expression data ($2), phylogeny dir ($3), selective pressure dir ($4),
#         synteny dir ($5), output CSV ($6), [reference categories JSON ($7)]
# Outputs: Ranked candidates CSV with evidence completeness scores
# Author: Jorge L. Perez-Moreno, Ph.D.

"""
Improved Ranking Algorithm

Key improvements over original:
1. Separate purifying and positive selection scores (not combined |log(omega)|)
2. Weighted reference distances (chemoreceptors weighted higher than other GPCRs)
3. Quantitative synteny scoring (uses anchor counts, not just binary)
4. Proper missing data handling (tracks evidence completeness)
5. Configurable via environment variables
"""

import pandas as pd
import numpy as np
import os
import sys
import json
from pathlib import Path
from ete3 import Tree

# --- Input Arguments ---
candidates_file = sys.argv[1]
expression_file = sys.argv[2]
phylo_dir = sys.argv[3]
selective_dir = sys.argv[4]
synteny_dir = sys.argv[5]  # Changed: now takes directory instead of single file
output_file = sys.argv[6]
ref_categories_file = sys.argv[7] if len(sys.argv) > 7 else None

# --- Configuration from Environment ---
# Ranking weights
PHYLO_WEIGHT = float(os.getenv('PHYLO_WEIGHT', 2))
PURIFYING_WEIGHT = float(os.getenv('PURIFYING_WEIGHT', 1))  # New: separate from positive
POSITIVE_WEIGHT = float(os.getenv('POSITIVE_WEIGHT', 1))    # New: separate from purifying
SYNTENY_WEIGHT = float(os.getenv('SYNTENY_WEIGHT', 3))
EXPR_WEIGHT = float(os.getenv('EXPR_WEIGHT', 1))
LSE_DEPTH_WEIGHT = float(os.getenv('LSE_DEPTH_WEIGHT', 1))

# Reference weighting
CHEMORECEPTOR_REF_WEIGHT = float(os.getenv('CHEMORECEPTOR_REF_WEIGHT', 2.0))  # New
OTHER_GPCR_REF_WEIGHT = float(os.getenv('OTHER_GPCR_REF_WEIGHT', 1.0))        # New

# Statistical thresholds
ABSREL_FDR_THRESHOLD = float(os.getenv('ABSREL_FDR_THRESHOLD', 0.05))  # New: for significant selection
BOOTSTRAP_THRESHOLD = float(os.getenv('BOOTSTRAP_THRESHOLD', 70))       # New: minimum support

# LSE threshold percentile
LSE_DEPTH_PERCENTILE = float(os.getenv('LSE_DEPTH_PERCENTILE', 75))

# Local database directory (for reference categories)
LOCAL_DB_DIR = os.getenv('LOCAL_DB_DIR', '')

# --- Load Input Data ---
candidates = pd.read_csv(candidates_file, names=['id'])

# Expression data (handle missing gracefully)
if os.path.exists(expression_file):
    expr_data = pd.read_csv(expression_file, names=['id', 'weight'])
else:
    expr_data = pd.DataFrame(columns=['id', 'weight'])
    print(f"Warning: Expression file not found: {expression_file}", file=sys.stderr)

# --- Load Reference Categories ---
chemoreceptor_refs = set()
other_gpcr_refs = set()

# Try to load from provided file or local database
ref_cat_path = None
if ref_categories_file and os.path.exists(ref_categories_file):
    ref_cat_path = ref_categories_file
elif LOCAL_DB_DIR:
    potential_path = Path(LOCAL_DB_DIR) / "gpcrdb" / "reference_categories.json"
    if potential_path.exists():
        ref_cat_path = str(potential_path)

if ref_cat_path:
    try:
        with open(ref_cat_path) as f:
            ref_cats = json.load(f)
            chemoreceptor_refs = {r.get('entry_name', '') for r in ref_cats.get('chemoreceptors', [])}
            other_gpcr_refs = {r.get('entry_name', '') for r in ref_cats.get('other_gpcrs', [])}
            print(f"Loaded reference categories: {len(chemoreceptor_refs)} chemoreceptors, "
                  f"{len(other_gpcr_refs)} other GPCRs", file=sys.stderr)
    except Exception as e:
        print(f"Warning: Could not load reference categories: {e}", file=sys.stderr)

# --- Load Phylogenetic Tree ---
tree_file = f"{phylo_dir}/all_berghia_refs.treefile"
if not os.path.exists(tree_file):
    print(f"Error: Tree file not found: {tree_file}", file=sys.stderr)
    sys.exit(1)

t = Tree(tree_file, format=1)  # format=1 for IQ-TREE output with support values

# Identify reference sequences and categorize them
ref_ids = [leaf.name for leaf in t if leaf.name.startswith('ref_')]


def categorize_reference(ref_name):
    """
    Categorize a reference as chemoreceptor or other GPCR.

    Returns weight for this reference type.
    """
    # Check against loaded categories
    if ref_name in chemoreceptor_refs:
        return CHEMORECEPTOR_REF_WEIGHT

    # Heuristic: check for chemoreceptor keywords in name
    chemoreceptor_keywords = ['olfr', 'or', 'taste', 'gust', 'odorant', 'chem', 'vomero']
    name_lower = ref_name.lower()
    if any(kw in name_lower for kw in chemoreceptor_keywords):
        return CHEMORECEPTOR_REF_WEIGHT

    return OTHER_GPCR_REF_WEIGHT


# --- Load aBSREL dN/dS Results with FDR Correction ---
def load_absrel_with_fdr(selective_dir):
    """
    Load aBSREL results and apply Benjamini-Hochberg FDR correction.

    Returns dict: branch_id -> {omega, p_value, corrected_p, significant}
    """
    dnds_data = {}
    absrel_files = ['absrel_results.csv', 'absrel_results_lse.csv']

    all_pvalues = []
    all_entries = []

    for absrel_file in absrel_files:
        absrel_path = os.path.join(selective_dir, absrel_file)
        if os.path.exists(absrel_path):
            try:
                absrel_df = pd.read_csv(absrel_path)
                for _, row in absrel_df.iterrows():
                    branch_id = row.get('branch_id', row.get('id', ''))
                    omega = row.get('omega', row.get('dnds', row.get('dN/dS', None)))
                    p_value = row.get('p_value', row.get('pvalue', row.get('p-value', 1.0)))

                    if pd.notna(omega) and branch_id:
                        all_pvalues.append(float(p_value) if pd.notna(p_value) else 1.0)
                        all_entries.append({
                            'branch_id': str(branch_id),
                            'omega': float(omega),
                            'p_value': float(p_value) if pd.notna(p_value) else 1.0
                        })
            except Exception as e:
                print(f"Warning: Could not parse {absrel_path}: {e}", file=sys.stderr)

    # Apply Benjamini-Hochberg FDR correction
    if all_pvalues:
        corrected_pvalues = benjamini_hochberg(all_pvalues)

        for entry, corrected_p in zip(all_entries, corrected_pvalues):
            entry['corrected_p'] = corrected_p
            entry['significant'] = corrected_p < ABSREL_FDR_THRESHOLD
            dnds_data[entry['branch_id']] = entry

    return dnds_data


def benjamini_hochberg(pvalues):
    """
    Apply Benjamini-Hochberg FDR correction to a list of p-values.

    Returns list of corrected p-values in the same order.
    """
    n = len(pvalues)
    if n == 0:
        return []

    # Sort p-values and track original indices
    indexed_pvals = [(p, i) for i, p in enumerate(pvalues)]
    indexed_pvals.sort(key=lambda x: x[0])

    # Calculate corrected p-values
    corrected = [0.0] * n
    prev_corrected = 1.0

    for rank, (pval, orig_idx) in enumerate(reversed(indexed_pvals), 1):
        # BH formula: p_corrected = p * n / rank
        # But we need to ensure monotonicity
        corrected_p = min(pval * n / (n - rank + 1), prev_corrected)
        corrected[orig_idx] = min(corrected_p, 1.0)
        prev_corrected = corrected_p

    return corrected


dnds_data = load_absrel_with_fdr(selective_dir)

# --- Load Synteny Data (Quantitative) ---
def load_synteny_scores(synteny_dir):
    """
    Load synteny data and compute quantitative scores.

    Returns dict: gene_id -> synteny_score (0.0 to 1.0)
    """
    synteny_scores = {}
    anchor_counts = {}
    max_anchors = 1

    # Load synteny IDs file
    synteny_ids_file = os.path.join(synteny_dir, 'synteny_ids.txt')
    if os.path.exists(synteny_ids_file):
        with open(synteny_ids_file) as f:
            for line in f:
                gene_id = line.strip()
                if gene_id:
                    anchor_counts[gene_id] = anchor_counts.get(gene_id, 0) + 1
                    max_anchors = max(max_anchors, anchor_counts[gene_id])

    # Load collinearity files for more detailed scoring
    for collin_file in Path(synteny_dir).glob('*.collinearity'):
        try:
            with open(collin_file) as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        gene1, gene2 = parts[1], parts[2]
                        anchor_counts[gene1] = anchor_counts.get(gene1, 0) + 1
                        anchor_counts[gene2] = anchor_counts.get(gene2, 0) + 1
                        max_anchors = max(max_anchors, anchor_counts[gene1], anchor_counts[gene2])
        except Exception as e:
            print(f"Warning: Could not parse {collin_file}: {e}", file=sys.stderr)

    # Normalize to 0-1 range
    for gene_id, count in anchor_counts.items():
        synteny_scores[gene_id] = count / max_anchors

    return synteny_scores


synteny_scores = load_synteny_scores(synteny_dir)

# --- Scoring Functions ---

def weighted_distance_to_refs(node_name, tree, ref_ids):
    """
    Calculate weighted phylogenetic distance to references.

    Chemoreceptor references are weighted higher than other GPCRs.
    Returns inverse weighted distance (higher = closer to important refs).
    """
    if node_name in ref_ids:
        return float('inf')  # Reference itself

    try:
        total_weighted_inverse = 0.0
        for ref in ref_ids:
            try:
                distance = tree.get_distance(node_name, ref)
                weight = categorize_reference(ref)
                # Weighted inverse distance
                total_weighted_inverse += weight / (distance + 1e-6)
            except Exception:
                continue

        return total_weighted_inverse if total_weighted_inverse > 0 else 0.0
    except Exception:
        return 0.0


def get_selection_scores(candidate_id, dnds_data):
    """
    Calculate separate purifying and positive selection scores.

    Returns: (purifying_score, positive_score, is_significant)

    Biological interpretation:
    - Purifying selection (omega < 1): Conserved function, constrained evolution
    - Positive selection (omega > 1): Adaptive evolution, potential novel function
    - Neutral (omega ~ 1): Relaxed constraint or pseudogenization

    For chemoreceptor discovery:
    - Strong purifying: likely functional receptor with conserved role
    - Strong positive: potential lineage-specific adaptation/specialization
    """
    if candidate_id not in dnds_data:
        return 0.0, 0.0, False

    entry = dnds_data[candidate_id]
    omega = entry['omega']
    is_significant = entry.get('significant', False)

    # Clamp omega to avoid extreme values
    omega = max(omega, 0.001)
    omega = min(omega, 100.0)

    if omega < 1.0:
        # Purifying selection: score increases as omega decreases
        purifying_score = abs(np.log(omega))
        positive_score = 0.0
    else:
        # Positive selection: score increases as omega increases
        purifying_score = 0.0
        positive_score = np.log(omega)

    # Boost scores if statistically significant
    if is_significant:
        purifying_score *= 1.5
        positive_score *= 1.5

    return purifying_score, positive_score, is_significant


def get_synteny_score(candidate_id, synteny_scores):
    """
    Get quantitative synteny score for a candidate.

    Returns score between 0.0 and 1.0.
    """
    return synteny_scores.get(candidate_id, 0.0)


def get_expression_score(candidate_id, expr_data):
    """
    Get expression score, handling missing data appropriately.

    Returns: (score, has_data)
    """
    if candidate_id in expr_data['id'].values:
        score = expr_data[expr_data['id'] == candidate_id]['weight'].sum()
        return score, True
    return 0.0, False


def get_lse_depth_score(candidate_id, tree, threshold):
    """
    Calculate LSE depth score based on tree position.

    Returns: (depth_score, raw_depth)
    """
    try:
        nodes = tree.search_nodes(name=candidate_id)
        if nodes:
            depth = nodes[0].get_distance(tree)
            # Only score if above threshold (deep in tree)
            if depth > threshold:
                return depth, depth
    except Exception:
        pass
    return 0.0, 0.0


# --- Calculate Scores for All Candidates ---

# First pass: collect all depths for percentile calculation
all_depths = []
for cand_id in candidates['id']:
    try:
        nodes = t.search_nodes(name=cand_id)
        if nodes:
            all_depths.append(nodes[0].get_distance(t))
    except Exception:
        pass

lse_threshold = np.percentile(all_depths, LSE_DEPTH_PERCENTILE) if all_depths else 0.5

# Second pass: calculate all scores
results = []

for cand_id in candidates['id']:
    # Phylogenetic score (weighted distance to refs)
    phylo_score = weighted_distance_to_refs(cand_id, t, ref_ids)

    # Selection scores (separate purifying and positive)
    purifying_score, positive_score, selection_significant = get_selection_scores(cand_id, dnds_data)

    # Synteny score (quantitative)
    synteny_score = get_synteny_score(cand_id, synteny_scores)

    # Expression score
    expr_score, has_expression = get_expression_score(cand_id, expr_data)

    # LSE depth score
    lse_score, raw_depth = get_lse_depth_score(cand_id, t, lse_threshold)

    # Track evidence completeness
    evidence_count = sum([
        phylo_score > 0,
        purifying_score > 0 or positive_score > 0,
        synteny_score > 0,
        has_expression,
        lse_score > 0
    ])
    evidence_completeness = evidence_count / 5.0

    results.append({
        'id': cand_id,
        'phylo_score': phylo_score,
        'purifying_score': purifying_score,
        'positive_score': positive_score,
        'selection_significant': selection_significant,
        'synteny_score': synteny_score,
        'expression_score': expr_score,
        'has_expression_data': has_expression,
        'lse_depth_score': lse_score,
        'raw_tree_depth': raw_depth,
        'evidence_completeness': evidence_completeness
    })

df = pd.DataFrame(results)

# --- Normalize Scores ---
# Normalize continuous scores to [0,1] range for fair comparison
normalize_cols = ['phylo_score', 'purifying_score', 'positive_score', 'expression_score', 'lse_depth_score']

for col in normalize_cols:
    col_min = df[col].min()
    col_max = df[col].max()
    col_range = col_max - col_min
    if col_range > 0:
        df[f'{col}_norm'] = (df[col] - col_min) / col_range
    else:
        df[f'{col}_norm'] = 0.0

# Synteny is already 0-1
df['synteny_score_norm'] = df['synteny_score']

# --- Calculate Final Rank Score ---
df['rank_score'] = (
    df['phylo_score_norm'] * PHYLO_WEIGHT +
    df['purifying_score_norm'] * PURIFYING_WEIGHT +
    df['positive_score_norm'] * POSITIVE_WEIGHT +
    df['synteny_score_norm'] * SYNTENY_WEIGHT +
    df['expression_score_norm'] * EXPR_WEIGHT +
    df['lse_depth_score_norm'] * LSE_DEPTH_WEIGHT
)

# --- Add Confidence Tier ---
def assign_confidence_tier(row):
    """
    Assign confidence tier based on evidence completeness and score.

    Tiers:
    - High: Strong evidence across multiple sources
    - Medium: Good evidence but some missing data
    - Low: Limited evidence or low scores
    """
    completeness = row['evidence_completeness']
    score = row['rank_score']

    # Normalize rank score for tier assignment
    max_possible = PHYLO_WEIGHT + PURIFYING_WEIGHT + POSITIVE_WEIGHT + SYNTENY_WEIGHT + EXPR_WEIGHT + LSE_DEPTH_WEIGHT
    score_pct = score / max_possible if max_possible > 0 else 0

    if completeness >= 0.6 and score_pct >= 0.5:
        return 'High'
    elif completeness >= 0.4 or score_pct >= 0.3:
        return 'Medium'
    else:
        return 'Low'


df['confidence_tier'] = df.apply(assign_confidence_tier, axis=1)

# --- Sort and Save ---
df_sorted = df.sort_values('rank_score', ascending=False)

# Select columns for output
output_cols = [
    'id', 'rank_score', 'confidence_tier', 'evidence_completeness',
    'phylo_score', 'purifying_score', 'positive_score', 'selection_significant',
    'synteny_score', 'expression_score', 'has_expression_data',
    'lse_depth_score', 'raw_tree_depth'
]

df_sorted[output_cols].to_csv(output_file, index=False)

# --- Summary Statistics ---
print(f"\nRanking Summary:", file=sys.stderr)
print(f"  Total candidates: {len(df)}", file=sys.stderr)
print(f"  High confidence: {len(df[df['confidence_tier'] == 'High'])}", file=sys.stderr)
print(f"  Medium confidence: {len(df[df['confidence_tier'] == 'Medium'])}", file=sys.stderr)
print(f"  Low confidence: {len(df[df['confidence_tier'] == 'Low'])}", file=sys.stderr)
print(f"  With significant selection: {df['selection_significant'].sum()}", file=sys.stderr)
print(f"  With synteny support: {(df['synteny_score'] > 0).sum()}", file=sys.stderr)
print(f"  With expression data: {df['has_expression_data'].sum()}", file=sys.stderr)
print(f"\nOutput written to: {output_file}", file=sys.stderr)
