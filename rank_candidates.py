#!/usr/bin/env python3
# rank_candidates.py
# Purpose: Rank GPCR candidates using phylogenetic proximity, dN/dS, expression, and synteny.
# Inputs: Candidate IDs ($1), expression data ($2), phylogeny dir ($3), selective pressure dir ($4),
#         synteny dir ($5), output CSV ($6), [reference categories JSON ($7)]
# Outputs: Ranked candidates CSV with evidence completeness scores and sensitivity analysis
# Author: Jorge L. Perez-Moreno, Ph.D.

"""
Improved Ranking Algorithm with Sensitivity Analysis

Key improvements over original:
1. Separate purifying and positive selection scores (not combined |log(omega)|)
2. Weighted reference distances (chemoreceptors weighted higher than other GPCRs)
3. Quantitative synteny scoring (uses anchor counts, not just binary)
4. Proper missing data handling (tracks evidence completeness)
5. Configurable via environment variables
6. Sensitivity analysis: tests weight parameter variations, reports ranking stability
7. Cross-validation support for reference classification validation
"""

import pandas as pd
import numpy as np
import os
import sys
import json
import itertools
from pathlib import Path
from scipy import stats
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

# NEW: Phase 1 - Chemosensory expression weight
CHEMOSENSORY_EXPR_WEIGHT = float(os.getenv('CHEMOSENSORY_EXPR_WEIGHT', 3))
CHEMOSENSORY_HARD_FILTER = os.getenv('CHEMOSENSORY_HARD_FILTER', 'false').lower() == 'true'

# NEW: Phase 2 - G-protein co-expression weight
GPROTEIN_COEXPR_WEIGHT = float(os.getenv('GPROTEIN_COEXPR_WEIGHT', 2))

# NEW: Phase 3 - ECL divergence weight
ECL_DIVERGENCE_WEIGHT = float(os.getenv('ECL_DIVERGENCE_WEIGHT', 1.5))

# NEW: Phase 5 - CAFE expansion weight
EXPANSION_WEIGHT = float(os.getenv('EXPANSION_WEIGHT', 1.5))
PREFERRED_EXPANSION_LEVELS = os.getenv('PREFERRED_EXPANSION_LEVELS', 'Aeolid-specific').split(',')

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

# Sensitivity analysis settings
RUN_SENSITIVITY = os.getenv('RUN_SENSITIVITY', 'true').lower() == 'true'
SENSITIVITY_PERTURBATION = float(os.getenv('SENSITIVITY_PERTURBATION', 0.5))  # +/- 50% weight variation
SENSITIVITY_ITERATIONS = int(os.getenv('SENSITIVITY_ITERATIONS', 100))  # Monte Carlo iterations

# Cross-validation settings
RUN_CROSSVAL = os.getenv('RUN_CROSSVAL', 'false').lower() == 'true'
CROSSVAL_FOLDS = int(os.getenv('CROSSVAL_FOLDS', 5))

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
    Only applies correction to p-values that aren't already corrected by HyPhy.

    Returns dict: branch_id -> {omega, p_value, corrected_p, significant}
    """
    dnds_data = {}
    absrel_files = ['absrel_results.csv', 'absrel_results_lse.csv']

    # Separate entries by whether they need FDR correction
    needs_correction = []  # (p_value, entry) pairs for uncorrected p-values
    already_corrected = []  # entries with HyPhy-corrected p-values

    for absrel_file in absrel_files:
        absrel_path = os.path.join(selective_dir, absrel_file)
        if os.path.exists(absrel_path):
            try:
                absrel_df = pd.read_csv(absrel_path)
                for _, row in absrel_df.iterrows():
                    branch_id = row.get('branch_id', row.get('id', ''))
                    omega = row.get('omega', row.get('dnds', row.get('dN/dS', None)))
                    p_value = row.get('p_value', row.get('pvalue', row.get('p-value', 1.0)))
                    is_already_corrected = row.get('is_already_corrected', 0)
                    corrected_p_from_hyphy = row.get('corrected_p_value')

                    if pd.notna(omega) and branch_id:
                        entry = {
                            'branch_id': str(branch_id),
                            'omega': float(omega),
                            'p_value': float(p_value) if pd.notna(p_value) else 1.0
                        }

                        # Check if HyPhy already provided corrected p-value
                        if is_already_corrected == 1 and pd.notna(corrected_p_from_hyphy):
                            entry['corrected_p'] = float(corrected_p_from_hyphy)
                            entry['significant'] = entry['corrected_p'] < ABSREL_FDR_THRESHOLD
                            already_corrected.append(entry)
                        else:
                            needs_correction.append((entry['p_value'], entry))
            except Exception as e:
                print(f"Warning: Could not parse {absrel_path}: {e}", file=sys.stderr)

    # Apply Benjamini-Hochberg FDR correction only to uncorrected p-values
    if needs_correction:
        pvalues = [p for p, _ in needs_correction]
        corrected_pvalues = benjamini_hochberg(pvalues)

        for (_, entry), corrected_p in zip(needs_correction, corrected_pvalues):
            entry['corrected_p'] = corrected_p
            entry['significant'] = corrected_p < ABSREL_FDR_THRESHOLD
            dnds_data[entry['branch_id']] = entry

    # Add already-corrected entries
    for entry in already_corrected:
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

# --- Derive results directory from phylo_dir ---
# phylo_dir is typically {RESULTS_DIR}/phylogenies/protein
results_dir = str(Path(phylo_dir).parent.parent)

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

    Returns: (score, has_data)
        - score: float between 0.0 and 1.0
        - has_data: bool indicating if synteny data exists for this candidate
    """
    if candidate_id in synteny_scores:
        return synteny_scores[candidate_id], True
    return 0.0, False


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


# --- NEW: Phase 1 - Chemosensory Expression Scoring ---

def load_chemosensory_expression(results_dir):
    """
    Load chemosensory expression summary from process_expression.py output.

    Returns dict: gene_id -> {tau_index, enrichment, is_specific, has_data}
    """
    expr_file = os.path.join(results_dir, 'expression', 'expression_summary.csv')
    expr_data = {}

    if not os.path.exists(expr_file):
        return expr_data

    try:
        df = pd.read_csv(expr_file)
        for _, row in df.iterrows():
            gene_id = str(row['gene_id'])
            expr_data[gene_id] = {
                'tau_index': row.get('tau_index', 0.0),
                'chemosensory_enrichment': row.get('chemosensory_enrichment', 0.0),
                'mean_chemosensory_tpm': row.get('mean_chemosensory_tpm', 0.0),
                'expressed_in_chemosensory': row.get('expressed_in_chemosensory', False),
                'is_chemosensory_specific': row.get('is_chemosensory_specific', False),
                'has_data': True
            }
        print(f"Loaded chemosensory expression data for {len(expr_data)} genes", file=sys.stderr)
    except Exception as e:
        print(f"Warning: Could not load chemosensory expression: {e}", file=sys.stderr)

    return expr_data


def get_chemosensory_expression_score(candidate_id, chemo_expr_data):
    """
    Calculate chemosensory expression score.

    Components:
    1. Is expressed in any chemosensory tissue (binary bonus)
    2. Tau index (tissue specificity)
    3. Chemosensory enrichment (fold-change vs other tissues)

    Returns: (score, has_data)
    """
    if candidate_id not in chemo_expr_data:
        return 0.0, False

    data = chemo_expr_data[candidate_id]

    if not data.get('has_data', False):
        return 0.0, False

    score = 0.0

    # Component 1: Expressed in chemosensory tissue
    if data.get('expressed_in_chemosensory', False):
        score += 0.3

    # Component 2: Tissue specificity (tau index)
    tau = data.get('tau_index', 0.0)
    if pd.notna(tau):
        score += tau * 0.3  # Max 0.3 contribution

    # Component 3: Chemosensory enrichment (log fold-change, capped)
    enrichment = data.get('chemosensory_enrichment', 0.0)
    if enrichment > 1:
        # Log2 fold-change, capped at 0.4 contribution
        log_enrichment = min(np.log2(enrichment) / 5, 1.0)  # Cap at 5-fold (log2=2.32)
        score += log_enrichment * 0.4

    # Bonus for chemosensory-specific expression
    if data.get('is_chemosensory_specific', False):
        score *= 1.5  # 50% bonus

    return score, True


# --- NEW: Phase 2 - G-protein Co-expression Scoring ---

def load_gprotein_coexpression(results_dir):
    """
    Load G-protein co-expression data.

    Returns dict: gpcr_id -> {gprotein_class, correlation, coexpr_score, has_data}
    """
    coexpr_file = os.path.join(results_dir, 'gproteins', 'gprotein_coexpression.csv')
    coexpr_data = {}

    if not os.path.exists(coexpr_file):
        return coexpr_data

    try:
        df = pd.read_csv(coexpr_file)
        for _, row in df.iterrows():
            gpcr_id = str(row['gpcr_id'])
            coexpr_data[gpcr_id] = {
                'best_gprotein': row.get('best_coexpr_gprotein', ''),
                'gprotein_class': row.get('gprotein_class', ''),
                'correlation': row.get('correlation', 0.0),
                'coexpr_tissue': row.get('coexpr_tissue', ''),
                'coexpr_score': row.get('coexpr_score', 0.0),
                'has_data': True
            }
        print(f"Loaded G-protein co-expression data for {len(coexpr_data)} GPCRs", file=sys.stderr)
    except Exception as e:
        print(f"Warning: Could not load G-protein co-expression: {e}", file=sys.stderr)

    return coexpr_data


def get_gprotein_coexpr_score(candidate_id, coexpr_data, olfactory_classes=None):
    """
    Calculate G-protein co-expression score.

    Prioritizes co-expression with olfactory G-proteins (Golf, Gi, Go).

    Returns: (score, has_data)
    """
    if olfactory_classes is None:
        olfactory_classes = ['Golf', 'Gi', 'Go']

    if candidate_id not in coexpr_data:
        return 0.0, False

    data = coexpr_data[candidate_id]

    if not data.get('has_data', False):
        return 0.0, False

    score = data.get('coexpr_score', 0.0)

    # Bonus for olfactory G-protein class
    gprotein_class = data.get('gprotein_class', '')
    if gprotein_class in olfactory_classes:
        score *= 1.5  # 50% bonus for olfactory G-proteins

    return score, True


# --- NEW: Phase 3 - ECL Divergence Scoring ---

def load_ecl_divergence(results_dir):
    """
    Load ECL divergence analysis results.

    Returns dict: gene_id -> {ecl_tm_ratio, ecl_divergence, tm_divergence, has_data}
    """
    ecl_file = os.path.join(results_dir, 'ecl_analysis', 'ecl_divergence.csv')
    ecl_data = {}

    if not os.path.exists(ecl_file):
        return ecl_data

    try:
        df = pd.read_csv(ecl_file)
        for _, row in df.iterrows():
            gene_id = str(row['gene_id'])
            ecl_data[gene_id] = {
                'ecl_tm_ratio': row.get('ecl_tm_ratio', 0.0),
                'ecl_divergence': row.get('ecl_divergence', 0.0),
                'tm_divergence': row.get('tm_divergence', 0.0),
                'ecl_divergence_score': row.get('ecl_divergence_score', 0.0),
                'has_data': True
            }
        print(f"Loaded ECL divergence data for {len(ecl_data)} genes", file=sys.stderr)
    except Exception as e:
        print(f"Warning: Could not load ECL divergence: {e}", file=sys.stderr)

    return ecl_data


def get_ecl_divergence_score(candidate_id, ecl_data, ratio_threshold=2.0):
    """
    Calculate ECL divergence score.

    High ratio (ECL divergent, TM conserved) suggests ligand-binding diversification.

    Returns: (score, has_data)
    """
    if candidate_id not in ecl_data:
        return 0.0, False

    data = ecl_data[candidate_id]

    if not data.get('has_data', False):
        return 0.0, False

    ratio = data.get('ecl_tm_ratio', 0.0)

    # Score based on ratio exceeding threshold
    if ratio >= ratio_threshold:
        # Normalize: ratio of 2 = base score, higher = better
        score = min(ratio / ratio_threshold, 3.0)  # Cap at 3x threshold
    else:
        # Below threshold, reduced score
        score = ratio / ratio_threshold * 0.5

    return score, True


# --- NEW: Phase 5 - CAFE Expansion Scoring ---

def load_cafe_expansion(results_dir):
    """
    Load CAFE expansion interpretation results.

    Returns dict: orthogroup -> {taxonomic_level, pvalue, expansion_fold, has_data}
    """
    cafe_file = os.path.join(results_dir, 'cafe', 'expansion_interpretation.csv')
    cafe_data = {}

    if not os.path.exists(cafe_file):
        return cafe_data

    try:
        df = pd.read_csv(cafe_file)
        for _, row in df.iterrows():
            og = str(row['orthogroup'])
            cafe_data[og] = {
                'taxonomic_level': row.get('taxonomic_level', ''),
                'expansion_pvalue': row.get('expansion_pvalue', 1.0),
                'expansion_fold': row.get('expansion_fold', 1.0),
                'berghia_copies': row.get('berghia_copies', 0),
                'ancestral_copies': row.get('ancestral_copies', 0),
                'has_data': True
            }
        print(f"Loaded CAFE expansion data for {len(cafe_data)} orthogroups", file=sys.stderr)
    except Exception as e:
        print(f"Warning: Could not load CAFE expansion: {e}", file=sys.stderr)

    return cafe_data


def load_gene_to_orthogroup(results_dir):
    """
    Load gene to orthogroup mapping.

    Returns dict: gene_id -> orthogroup
    """
    # Try multiple possible file locations
    possible_files = [
        os.path.join(results_dir, 'orthofinder', 'Orthogroups', 'Orthogroups.tsv'),
        os.path.join(results_dir, 'orthofinder', 'Orthogroups.tsv'),
        os.path.join(results_dir, 'orthology', 'gene_orthogroups.csv')
    ]

    gene_to_og = {}

    for og_file in possible_files:
        if os.path.exists(og_file):
            try:
                if og_file.endswith('.tsv'):
                    df = pd.read_csv(og_file, sep='\t')
                    # OrthoFinder format: Orthogroup, Species1, Species2, ...
                    for _, row in df.iterrows():
                        og = row.iloc[0]
                        for col in df.columns[1:]:
                            genes = str(row[col]).split(', ')
                            for gene in genes:
                                if gene and gene != 'nan':
                                    gene_to_og[gene.strip()] = og
                else:
                    df = pd.read_csv(og_file)
                    for _, row in df.iterrows():
                        gene_to_og[str(row['gene_id'])] = str(row['orthogroup'])
                print(f"Loaded gene-orthogroup mapping: {len(gene_to_og)} genes", file=sys.stderr)
                break
            except Exception as e:
                print(f"Warning: Could not parse {og_file}: {e}", file=sys.stderr)

    return gene_to_og


def get_expansion_score(candidate_id, cafe_data, gene_to_og, preferred_levels=None):
    """
    Calculate CAFE expansion score.

    Prioritizes expansions at preferred taxonomic levels (e.g., Aeolid-specific).

    Returns: (score, has_data)
    """
    if preferred_levels is None:
        preferred_levels = ['Aeolid-specific']

    # Get orthogroup for this candidate
    og = gene_to_og.get(candidate_id)
    if not og or og not in cafe_data:
        return 0.0, False

    data = cafe_data[og]

    if not data.get('has_data', False):
        return 0.0, False

    pvalue = data.get('expansion_pvalue', 1.0)

    # Base score from significance (-log10 p-value)
    if pvalue > 0 and pvalue < 1:
        base_score = -np.log10(pvalue)
    else:
        base_score = 0.0

    # Bonus for preferred taxonomic level
    tax_level = data.get('taxonomic_level', '')
    if tax_level in preferred_levels:
        base_score *= 1.5  # 50% bonus

    # Additional bonus for high expansion fold
    fold = data.get('expansion_fold', 1.0)
    if fold > 2:
        base_score *= min(fold / 2, 2.0)  # Up to 2x bonus

    return base_score, True


# --- Sensitivity Analysis Functions ---

def calculate_rank_score(df, weights):
    """
    Calculate rank scores with given weight parameters.

    Uses fair scoring that handles missing data (synteny, expression, etc.) by
    normalizing per-candidate based on available evidence.

    Args:
        df: DataFrame with normalized score columns and has_*_data flags
        weights: dict with keys for all weight parameters

    Returns:
        Series of rank scores
    """
    max_possible_weight = sum(weights.values())

    def calc_row(row):
        # Base scores always available
        score = (
            row['phylo_score_norm'] * weights.get('phylo', 0) +
            row['purifying_score_norm'] * weights.get('purifying', 0) +
            row['positive_score_norm'] * weights.get('positive', 0) +
            row['lse_depth_score_norm'] * weights.get('lse_depth', 0)
        )
        total_weight = (
            weights.get('phylo', 0) + weights.get('purifying', 0) +
            weights.get('positive', 0) + weights.get('lse_depth', 0)
        )

        # Conditionally add synteny and expression
        if row.get('has_synteny_data', False):
            score += row['synteny_score_norm'] * weights.get('synteny', 0)
            total_weight += weights.get('synteny', 0)
        if row.get('has_expression_data', False):
            score += row['expression_score_norm'] * weights.get('expr', 0)
            total_weight += weights.get('expr', 0)

        # NEW: Conditionally add Phase 1-5 scores
        if row.get('has_chemosensory_expr_data', False):
            score += row.get('chemosensory_expr_score_norm', 0) * weights.get('chemosensory_expr', 0)
            total_weight += weights.get('chemosensory_expr', 0)
        if row.get('has_gprotein_data', False):
            score += row.get('gprotein_coexpr_score_norm', 0) * weights.get('gprotein_coexpr', 0)
            total_weight += weights.get('gprotein_coexpr', 0)
        if row.get('has_ecl_data', False):
            score += row.get('ecl_divergence_score_norm', 0) * weights.get('ecl_divergence', 0)
            total_weight += weights.get('ecl_divergence', 0)
        if row.get('has_expansion_data', False):
            score += row.get('expansion_score_norm', 0) * weights.get('expansion', 0)
            total_weight += weights.get('expansion', 0)

        # Normalize to max possible
        if total_weight > 0:
            return (score / total_weight) * max_possible_weight
        return 0.0

    return df.apply(calc_row, axis=1)


def run_sensitivity_analysis(df, base_weights, perturbation=0.5, n_iterations=100):
    """
    Run Monte Carlo sensitivity analysis on weight parameters.

    Tests how robust rankings are to weight variations.

    Args:
        df: DataFrame with normalized score columns
        base_weights: dict of base weight values
        perturbation: fraction of weight to vary (+/- this amount)
        n_iterations: number of Monte Carlo iterations

    Returns:
        dict with sensitivity analysis results
    """
    np.random.seed(42)  # For reproducibility

    # Calculate baseline rankings
    baseline_scores = calculate_rank_score(df, base_weights)
    baseline_ranks = baseline_scores.rank(ascending=False)

    # Store rank variations for each candidate
    rank_variations = {cid: [] for cid in df['id']}
    weight_samples = []

    for _ in range(n_iterations):
        # Perturb weights randomly within range
        perturbed_weights = {}
        for key, base_val in base_weights.items():
            min_val = base_val * (1 - perturbation)
            max_val = base_val * (1 + perturbation)
            perturbed_weights[key] = np.random.uniform(min_val, max_val)

        weight_samples.append(perturbed_weights.copy())

        # Calculate perturbed scores and ranks
        perturbed_scores = calculate_rank_score(df, perturbed_weights)
        perturbed_ranks = perturbed_scores.rank(ascending=False)

        # Record rank for each candidate
        for idx, cid in enumerate(df['id']):
            rank_variations[cid].append(perturbed_ranks.iloc[idx])

    # Calculate stability metrics for each candidate
    stability_results = []
    for cid in df['id']:
        ranks = rank_variations[cid]
        baseline_rank = baseline_ranks[df['id'] == cid].values[0]

        stability_results.append({
            'id': cid,
            'baseline_rank': int(baseline_rank),
            'mean_rank': np.mean(ranks),
            'std_rank': np.std(ranks),
            'min_rank': int(np.min(ranks)),
            'max_rank': int(np.max(ranks)),
            'rank_range': int(np.max(ranks) - np.min(ranks)),
            'rank_stability': 1.0 - (np.std(ranks) / len(df))  # Higher = more stable
        })

    stability_df = pd.DataFrame(stability_results)

    # Calculate global metrics
    overall_stability = stability_df['rank_stability'].mean()
    top10_stable = stability_df.nsmallest(10, 'baseline_rank')['rank_stability'].mean()

    # Calculate weight importance via correlation
    weight_importance = {}
    for key in base_weights.keys():
        weight_values = [ws[key] for ws in weight_samples]
        # Correlate weight values with top candidate's rank
        top_cid = df.loc[baseline_scores.idxmax(), 'id']
        top_ranks = rank_variations[top_cid]
        corr, _ = stats.spearmanr(weight_values, top_ranks)
        weight_importance[key] = abs(corr) if not np.isnan(corr) else 0.0

    return {
        'stability_df': stability_df,
        'overall_stability': overall_stability,
        'top10_stability': top10_stable,
        'weight_importance': weight_importance,
        'n_iterations': n_iterations,
        'perturbation': perturbation
    }


def run_leave_one_out_crossval(df, tree, ref_ids, base_weights, n_folds=5):
    """
    Run k-fold cross-validation on reference classification.

    Tests how well the scoring function identifies known chemoreceptors
    when they're held out from the reference set.

    Args:
        df: DataFrame with candidate scores
        tree: phylogenetic tree
        ref_ids: list of reference sequence IDs
        base_weights: weight parameters
        n_folds: number of CV folds

    Returns:
        dict with cross-validation metrics
    """
    if len(ref_ids) < n_folds:
        print(f"Warning: Not enough references ({len(ref_ids)}) for {n_folds}-fold CV",
              file=sys.stderr)
        return None

    # Shuffle and split references into folds
    np.random.seed(42)
    shuffled_refs = np.random.permutation(ref_ids)
    fold_size = len(shuffled_refs) // n_folds

    cv_results = []

    for fold in range(n_folds):
        # Hold out this fold's references
        start_idx = fold * fold_size
        end_idx = start_idx + fold_size if fold < n_folds - 1 else len(shuffled_refs)
        held_out = set(shuffled_refs[start_idx:end_idx])
        training_refs = [r for r in ref_ids if r not in held_out]

        # Recalculate phylo scores with reduced reference set
        phylo_scores_cv = []
        for cid in df['id']:
            try:
                total_weighted_inverse = 0.0
                for ref in training_refs:
                    try:
                        distance = tree.get_distance(cid, ref)
                        weight = categorize_reference(ref)
                        total_weighted_inverse += weight / (distance + 1e-6)
                    except Exception:
                        continue
                phylo_scores_cv.append(total_weighted_inverse)
            except Exception:
                phylo_scores_cv.append(0.0)

        # Normalize and calculate rank scores
        phylo_array = np.array(phylo_scores_cv)
        if phylo_array.max() > phylo_array.min():
            phylo_norm = (phylo_array - phylo_array.min()) / (phylo_array.max() - phylo_array.min())
        else:
            phylo_norm = np.zeros_like(phylo_array)

        # Calculate scores with CV phylo scores
        cv_scores = (
            phylo_norm * base_weights['phylo'] +
            df['purifying_score_norm'].values * base_weights['purifying'] +
            df['positive_score_norm'].values * base_weights['positive'] +
            df['synteny_score_norm'].values * base_weights['synteny'] +
            df['expression_score_norm'].values * base_weights['expr'] +
            df['lse_depth_score_norm'].values * base_weights['lse_depth']
        )

        cv_ranks = pd.Series(cv_scores).rank(ascending=False)

        # Check how well held-out references rank
        for held_ref in held_out:
            if held_ref in df['id'].values:
                ref_idx = df[df['id'] == held_ref].index[0]
                ref_rank = cv_ranks.iloc[ref_idx]
                ref_percentile = 100 * (1 - ref_rank / len(df))
                cv_results.append({
                    'fold': fold,
                    'held_out_ref': held_ref,
                    'rank': int(ref_rank),
                    'percentile': ref_percentile,
                    'in_top_10pct': ref_percentile >= 90
                })

    if not cv_results:
        return None

    cv_df = pd.DataFrame(cv_results)

    return {
        'cv_results': cv_df,
        'mean_percentile': cv_df['percentile'].mean(),
        'std_percentile': cv_df['percentile'].std(),
        'top_10pct_rate': cv_df['in_top_10pct'].mean(),
        'n_folds': n_folds,
        'n_refs_tested': len(cv_df)
    }


# --- Calculate Scores for All Candidates ---

# First pass: collect all depths for percentile calculation
# Only include nodes on paths with sufficient bootstrap support (H12 fix)
all_depths = []
for cand_id in candidates['id']:
    try:
        nodes = t.search_nodes(name=cand_id)
        if nodes:
            node = nodes[0]
            # Check if the path to root has sufficient bootstrap support
            # Walk up the tree and verify bootstrap values
            has_good_support = True
            current = node
            while current.up:
                # IQ-TREE stores bootstrap as support attribute
                support = getattr(current.up, 'support', None)
                if support is not None and support < BOOTSTRAP_THRESHOLD:
                    has_good_support = False
                    break
                current = current.up

            if has_good_support:
                all_depths.append(node.get_distance(t))
    except Exception:
        pass

lse_threshold = np.percentile(all_depths, LSE_DEPTH_PERCENTILE) if all_depths else 0.5

# --- Load NEW Data Sources (Phases 1-5) ---
# Phase 1: Chemosensory expression data
chemo_expr_data = load_chemosensory_expression(results_dir)

# Phase 2: G-protein co-expression data
gprotein_coexpr_data = load_gprotein_coexpression(results_dir)

# Phase 3: ECL divergence data
ecl_divergence_data = load_ecl_divergence(results_dir)

# Phase 5: CAFE expansion data
cafe_expansion_data = load_cafe_expansion(results_dir)
gene_to_orthogroup = load_gene_to_orthogroup(results_dir)

# ECL divergence threshold from environment
ECL_RATIO_THRESHOLD = float(os.getenv('ECL_CONSERVATION_RATIO_THRESHOLD', 2.0))

# Second pass: calculate all scores
results = []

for cand_id in candidates['id']:
    # Phylogenetic score (weighted distance to refs)
    phylo_score = weighted_distance_to_refs(cand_id, t, ref_ids)

    # Selection scores (separate purifying and positive)
    purifying_score, positive_score, selection_significant = get_selection_scores(cand_id, dnds_data)

    # Synteny score (quantitative)
    synteny_score, has_synteny = get_synteny_score(cand_id, synteny_scores)

    # Expression score
    expr_score, has_expression = get_expression_score(cand_id, expr_data)

    # LSE depth score
    lse_score, raw_depth = get_lse_depth_score(cand_id, t, lse_threshold)

    # NEW: Phase 1 - Chemosensory expression score
    chemo_expr_score, has_chemo_expr = get_chemosensory_expression_score(cand_id, chemo_expr_data)

    # NEW: Phase 2 - G-protein co-expression score
    gprotein_score, has_gprotein = get_gprotein_coexpr_score(cand_id, gprotein_coexpr_data)

    # NEW: Phase 3 - ECL divergence score
    ecl_score, has_ecl = get_ecl_divergence_score(cand_id, ecl_divergence_data, ECL_RATIO_THRESHOLD)

    # NEW: Phase 5 - CAFE expansion score
    expansion_score, has_expansion = get_expansion_score(
        cand_id, cafe_expansion_data, gene_to_orthogroup, PREFERRED_EXPANSION_LEVELS
    )

    # Track evidence completeness (based on data availability, not just score > 0)
    # Updated to include new score sources
    evidence_sources = [
        phylo_score > 0,
        purifying_score > 0 or positive_score > 0,
        has_synteny,
        has_expression,
        lse_score > 0,
        has_chemo_expr,
        has_gprotein,
        has_ecl,
        has_expansion
    ]
    evidence_count = sum(evidence_sources)
    evidence_completeness = evidence_count / len(evidence_sources)

    results.append({
        'id': cand_id,
        'phylo_score': phylo_score,
        'purifying_score': purifying_score,
        'positive_score': positive_score,
        'selection_significant': selection_significant,
        'synteny_score': synteny_score,
        'has_synteny_data': has_synteny,
        'expression_score': expr_score,
        'has_expression_data': has_expression,
        'lse_depth_score': lse_score,
        'raw_tree_depth': raw_depth,
        # NEW scores
        'chemosensory_expr_score': chemo_expr_score,
        'has_chemosensory_expr_data': has_chemo_expr,
        'gprotein_coexpr_score': gprotein_score,
        'has_gprotein_data': has_gprotein,
        'ecl_divergence_score': ecl_score,
        'has_ecl_data': has_ecl,
        'expansion_score': expansion_score,
        'has_expansion_data': has_expansion,
        'evidence_completeness': evidence_completeness
    })

df = pd.DataFrame(results)

# --- Normalize Scores ---
# Normalize continuous scores to [0,1] range for fair comparison
normalize_cols = [
    'phylo_score', 'purifying_score', 'positive_score', 'expression_score', 'lse_depth_score',
    # NEW score columns
    'chemosensory_expr_score', 'gprotein_coexpr_score', 'ecl_divergence_score', 'expansion_score'
]

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

# --- Calculate Final Rank Score with Fair Missing Data Handling ---
# For candidates missing synteny or expression data, we normalize by available weights
# This prevents penalizing candidates from species without genome assemblies


def calculate_fair_rank_score(row):
    """
    Calculate rank score with per-candidate weight normalization.

    Candidates missing data (synteny, expression, chemosensory, etc.) are scored
    only on available evidence, then normalized to be comparable with candidates
    having full data.

    Updated to include all new Phase 1-5 scoring components.
    """
    # Base scores that are always available
    score = (
        row['phylo_score_norm'] * PHYLO_WEIGHT +
        row['purifying_score_norm'] * PURIFYING_WEIGHT +
        row['positive_score_norm'] * POSITIVE_WEIGHT +
        row['lse_depth_score_norm'] * LSE_DEPTH_WEIGHT
    )
    total_weight = PHYLO_WEIGHT + PURIFYING_WEIGHT + POSITIVE_WEIGHT + LSE_DEPTH_WEIGHT

    # Add synteny if data is available
    if row['has_synteny_data']:
        score += row['synteny_score_norm'] * SYNTENY_WEIGHT
        total_weight += SYNTENY_WEIGHT

    # Add expression if data is available
    if row['has_expression_data']:
        score += row['expression_score_norm'] * EXPR_WEIGHT
        total_weight += EXPR_WEIGHT

    # NEW: Add chemosensory expression if data is available (Phase 1)
    if row.get('has_chemosensory_expr_data', False):
        score += row['chemosensory_expr_score_norm'] * CHEMOSENSORY_EXPR_WEIGHT
        total_weight += CHEMOSENSORY_EXPR_WEIGHT

    # NEW: Add G-protein co-expression if data is available (Phase 2)
    if row.get('has_gprotein_data', False):
        score += row['gprotein_coexpr_score_norm'] * GPROTEIN_COEXPR_WEIGHT
        total_weight += GPROTEIN_COEXPR_WEIGHT

    # NEW: Add ECL divergence if data is available (Phase 3)
    if row.get('has_ecl_data', False):
        score += row['ecl_divergence_score_norm'] * ECL_DIVERGENCE_WEIGHT
        total_weight += ECL_DIVERGENCE_WEIGHT

    # NEW: Add CAFE expansion if data is available (Phase 5)
    if row.get('has_expansion_data', False):
        score += row['expansion_score_norm'] * EXPANSION_WEIGHT
        total_weight += EXPANSION_WEIGHT

    # Normalize to max possible score (as if all weights were available)
    max_possible_weight = (
        PHYLO_WEIGHT + PURIFYING_WEIGHT + POSITIVE_WEIGHT +
        SYNTENY_WEIGHT + EXPR_WEIGHT + LSE_DEPTH_WEIGHT +
        CHEMOSENSORY_EXPR_WEIGHT + GPROTEIN_COEXPR_WEIGHT +
        ECL_DIVERGENCE_WEIGHT + EXPANSION_WEIGHT
    )

    # Scale score to be comparable across candidates with different data availability
    if total_weight > 0:
        normalized_score = (score / total_weight) * max_possible_weight
    else:
        normalized_score = 0.0

    return normalized_score


df['rank_score'] = df.apply(calculate_fair_rank_score, axis=1)

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

    # Normalize rank score for tier assignment (updated with new weights)
    max_possible = (
        PHYLO_WEIGHT + PURIFYING_WEIGHT + POSITIVE_WEIGHT + SYNTENY_WEIGHT +
        EXPR_WEIGHT + LSE_DEPTH_WEIGHT + CHEMOSENSORY_EXPR_WEIGHT +
        GPROTEIN_COEXPR_WEIGHT + ECL_DIVERGENCE_WEIGHT + EXPANSION_WEIGHT
    )
    score_pct = score / max_possible if max_possible > 0 else 0

    if completeness >= 0.6 and score_pct >= 0.5:
        return 'High'
    elif completeness >= 0.4 or score_pct >= 0.3:
        return 'Medium'
    else:
        return 'Low'


df['confidence_tier'] = df.apply(assign_confidence_tier, axis=1)

# --- Apply Hard Filter if Enabled (Phase 1) ---
if CHEMOSENSORY_HARD_FILTER and 'has_chemosensory_expr_data' in df.columns:
    original_count = len(df)
    # Only keep candidates with chemosensory expression
    df = df[df['has_chemosensory_expr_data'] & (df['chemosensory_expr_score'] > 0)]
    filtered_count = len(df)
    if original_count > filtered_count:
        print(f"CHEMOSENSORY_HARD_FILTER: Removed {original_count - filtered_count} candidates "
              f"without chemosensory expression", file=sys.stderr)

# --- Sort and Save ---
df_sorted = df.sort_values('rank_score', ascending=False)

# Select columns for output (updated with new score columns)
output_cols = [
    'id', 'rank_score', 'confidence_tier', 'evidence_completeness',
    'phylo_score', 'purifying_score', 'positive_score', 'selection_significant',
    'synteny_score', 'has_synteny_data', 'expression_score', 'has_expression_data',
    'lse_depth_score', 'raw_tree_depth',
    # NEW: Phase 1-5 scores
    'chemosensory_expr_score', 'has_chemosensory_expr_data',
    'gprotein_coexpr_score', 'has_gprotein_data',
    'ecl_divergence_score', 'has_ecl_data',
    'expansion_score', 'has_expansion_data'
]

df_sorted[output_cols].to_csv(output_file, index=False)

# --- Run Sensitivity Analysis ---
sensitivity_results = None
if RUN_SENSITIVITY:
    print(f"\nRunning sensitivity analysis ({SENSITIVITY_ITERATIONS} iterations)...", file=sys.stderr)

    base_weights = {
        'phylo': PHYLO_WEIGHT,
        'purifying': PURIFYING_WEIGHT,
        'positive': POSITIVE_WEIGHT,
        'synteny': SYNTENY_WEIGHT,
        'expr': EXPR_WEIGHT,
        'lse_depth': LSE_DEPTH_WEIGHT,
        # NEW: Phase 1-5 weights
        'chemosensory_expr': CHEMOSENSORY_EXPR_WEIGHT,
        'gprotein_coexpr': GPROTEIN_COEXPR_WEIGHT,
        'ecl_divergence': ECL_DIVERGENCE_WEIGHT,
        'expansion': EXPANSION_WEIGHT
    }

    sensitivity_results = run_sensitivity_analysis(
        df, base_weights,
        perturbation=SENSITIVITY_PERTURBATION,
        n_iterations=SENSITIVITY_ITERATIONS
    )

    # Add stability metrics to main output
    stability_df = sensitivity_results['stability_df']
    df_sorted = df_sorted.merge(
        stability_df[['id', 'rank_stability', 'std_rank', 'rank_range']],
        on='id',
        how='left'
    )

    # Write sensitivity analysis results
    output_dir = Path(output_file).parent
    sensitivity_file = output_dir / 'sensitivity_analysis.csv'
    stability_df.to_csv(sensitivity_file, index=False)

    # Write weight importance
    importance_file = output_dir / 'weight_importance.json'
    with open(importance_file, 'w') as f:
        json.dump({
            'weight_importance': sensitivity_results['weight_importance'],
            'overall_stability': sensitivity_results['overall_stability'],
            'top10_stability': sensitivity_results['top10_stability'],
            'perturbation': SENSITIVITY_PERTURBATION,
            'n_iterations': SENSITIVITY_ITERATIONS
        }, f, indent=2)

    print(f"  Overall ranking stability: {sensitivity_results['overall_stability']:.3f}", file=sys.stderr)
    print(f"  Top-10 candidates stability: {sensitivity_results['top10_stability']:.3f}", file=sys.stderr)
    print(f"  Weight importance (Spearman |rho|):", file=sys.stderr)
    for weight, importance in sorted(sensitivity_results['weight_importance'].items(),
                                     key=lambda x: x[1], reverse=True):
        print(f"    {weight}: {importance:.3f}", file=sys.stderr)

# --- Run Cross-Validation ---
cv_results = None
if RUN_CROSSVAL and len(ref_ids) >= CROSSVAL_FOLDS:
    print(f"\nRunning {CROSSVAL_FOLDS}-fold cross-validation on reference classification...",
          file=sys.stderr)

    base_weights = {
        'phylo': PHYLO_WEIGHT,
        'purifying': PURIFYING_WEIGHT,
        'positive': POSITIVE_WEIGHT,
        'synteny': SYNTENY_WEIGHT,
        'expr': EXPR_WEIGHT,
        'lse_depth': LSE_DEPTH_WEIGHT
    }

    cv_results = run_leave_one_out_crossval(
        df, t, ref_ids, base_weights,
        n_folds=CROSSVAL_FOLDS
    )

    if cv_results:
        # Write CV results
        output_dir = Path(output_file).parent
        cv_file = output_dir / 'crossval_results.csv'
        cv_results['cv_results'].to_csv(cv_file, index=False)

        cv_summary_file = output_dir / 'crossval_summary.json'
        with open(cv_summary_file, 'w') as f:
            json.dump({
                'mean_percentile': cv_results['mean_percentile'],
                'std_percentile': cv_results['std_percentile'],
                'top_10pct_recovery_rate': cv_results['top_10pct_rate'],
                'n_folds': cv_results['n_folds'],
                'n_refs_tested': cv_results['n_refs_tested']
            }, f, indent=2)

        print(f"  Mean held-out reference percentile: {cv_results['mean_percentile']:.1f}%",
              file=sys.stderr)
        print(f"  Top-10% recovery rate: {cv_results['top_10pct_rate']*100:.1f}%",
              file=sys.stderr)

# --- Update output columns if sensitivity was run ---
output_cols = [
    'id', 'rank_score', 'confidence_tier', 'evidence_completeness',
    'phylo_score', 'purifying_score', 'positive_score', 'selection_significant',
    'synteny_score', 'has_synteny_data', 'expression_score', 'has_expression_data',
    'lse_depth_score', 'raw_tree_depth'
]

if sensitivity_results and 'rank_stability' in df_sorted.columns:
    output_cols.extend(['rank_stability', 'std_rank', 'rank_range'])

# Re-save with additional columns
df_sorted[output_cols].to_csv(output_file, index=False)

# --- Summary Statistics ---
print(f"\nRanking Summary:", file=sys.stderr)
print(f"  Total candidates: {len(df)}", file=sys.stderr)
print(f"  High confidence: {len(df[df['confidence_tier'] == 'High'])}", file=sys.stderr)
print(f"  Medium confidence: {len(df[df['confidence_tier'] == 'Medium'])}", file=sys.stderr)
print(f"  Low confidence: {len(df[df['confidence_tier'] == 'Low'])}", file=sys.stderr)
print(f"  With significant selection: {df['selection_significant'].sum()}", file=sys.stderr)
print(f"  With synteny data: {df['has_synteny_data'].sum()} (score > 0: {(df['synteny_score'] > 0).sum()})", file=sys.stderr)
print(f"  With expression data: {df['has_expression_data'].sum()}", file=sys.stderr)

# Warn if synteny data is only available for a subset
synteny_coverage = df['has_synteny_data'].mean()
if 0 < synteny_coverage < 1.0:
    print(f"\nNote: Synteny data available for {synteny_coverage*100:.1f}% of candidates.", file=sys.stderr)
    print(f"  Scores are normalized per-candidate to avoid penalizing species without genomes.", file=sys.stderr)

if sensitivity_results:
    print(f"\nSensitivity analysis written to: {output_dir}/sensitivity_analysis.csv", file=sys.stderr)
    print(f"Weight importance written to: {output_dir}/weight_importance.json", file=sys.stderr)

if cv_results:
    print(f"Cross-validation results written to: {output_dir}/crossval_results.csv", file=sys.stderr)

# --- Score Component Correlation Analysis (M27) ---
# Analyze correlations between score components to detect redundancy or unexpected relationships
score_cols_for_corr = ['phylo_score', 'purifying_score', 'positive_score',
                       'synteny_score', 'expression_score', 'lse_depth_score']
available_cols = [col for col in score_cols_for_corr if col in df.columns]

if len(available_cols) >= 2:
    from scipy.stats import spearmanr
    import numpy as np

    output_dir = Path(output_file).parent
    corr_file = output_dir / 'score_correlations.csv'

    # Compute Spearman correlations between all score pairs
    corr_data = []
    for i, col1 in enumerate(available_cols):
        for col2 in available_cols[i+1:]:
            # Get non-null pairs
            mask = df[col1].notna() & df[col2].notna()
            if mask.sum() >= 5:  # Need at least 5 pairs for meaningful correlation
                rho, p_value = spearmanr(df.loc[mask, col1], df.loc[mask, col2])
                corr_data.append({
                    'component_1': col1.replace('_score', ''),
                    'component_2': col2.replace('_score', ''),
                    'spearman_rho': round(rho, 4),
                    'p_value': round(p_value, 6),
                    'n_pairs': int(mask.sum()),
                    'significant': p_value < 0.05
                })

    if corr_data:
        corr_df = pd.DataFrame(corr_data)
        corr_df.to_csv(corr_file, index=False)

        # Flag highly correlated pairs (potential redundancy)
        high_corr = corr_df[(corr_df['spearman_rho'].abs() > 0.7) & corr_df['significant']]
        if len(high_corr) > 0:
            print(f"\nWarning: Highly correlated score components detected:", file=sys.stderr)
            for _, row in high_corr.iterrows():
                print(f"  {row['component_1']} <-> {row['component_2']}: rho={row['spearman_rho']:.2f}",
                      file=sys.stderr)
            print("  Consider adjusting weights if this represents redundant information.", file=sys.stderr)

        print(f"Score correlations written to: {corr_file}", file=sys.stderr)

print(f"\nOutput written to: {output_file}", file=sys.stderr)
