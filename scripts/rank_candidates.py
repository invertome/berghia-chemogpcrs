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

# Import unit-tested pure functions from the library module
# (see tests/unit/test_ranking_lib.py for coverage)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _rank_candidates_lib import (
    benjamini_hochberg as _bh_corrected,
    categorize_reference as _categorize_reference_corrected,
    extract_branch_omega,
    get_selection_scores as _get_selection_scores_corrected,
    calculate_fair_rank_score as _calculate_fair_rank_score_corrected,
    normalize_synteny_counts,
)

# Conditional ete3 import - may fail on Python 3.13+ due to removed cgi module
try:
    from ete3 import Tree
    ETE3_AVAILABLE = True
except ImportError as e:
    print(f"Warning: ete3 not available ({e}). Phylogenetic features will be limited.",
          file=sys.stderr)
    Tree = None
    ETE3_AVAILABLE = False

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

# Bead -8st: per-OG dN/dS reliability gating. The dN/dS axis is built from
# aBSREL run on a per-OG codon alignment that depends on having recovered
# CDS for the orthogroup's reference members. When few references have CDS
# (chemoreceptor LSE expansions where miniprot fails on most paralogs), the
# omega estimate is statistically unreliable AND there is no way to tell
# from the omega alone whether it's a real signal. We multiply the dN/dS
# axis contribution by min(1.0, n_ref_cds_in_og / DNDS_RELIABILITY_FULL),
# applied to BOTH the score AND the total_weight denominator so a candidate
# in an under-supported OG falls back to the other axes via fair-scoring,
# rather than getting a score artificially shaped by noise.
#
# 5 references is the minimum for a defensible omega estimate (matches the
# "medium" threshold in scripts/add_og_coverage_columns.py); above 10 the
# estimate is robust ("high"). The default fully unlocks the dN/dS axis at
# n_ref_cds == DNDS_RELIABILITY_FULL.
DNDS_RELIABILITY_FULL = float(os.getenv('DNDS_RELIABILITY_FULL', 5))
SYNTENY_WEIGHT = float(os.getenv('SYNTENY_WEIGHT', 3))
EXPR_WEIGHT = float(os.getenv('EXPR_WEIGHT', 1))
LSE_DEPTH_WEIGHT = float(os.getenv('LSE_DEPTH_WEIGHT', 1))

# NEW: Phase 1 - Chemosensory expression weight
CHEMOSENSORY_EXPR_WEIGHT = float(os.getenv('CHEMOSENSORY_EXPR_WEIGHT', 3))
CHEMOSENSORY_HARD_FILTER = os.getenv('CHEMOSENSORY_HARD_FILTER', 'false').lower() == 'true'

# Tissue-specific weights for chemosensory scoring (rhinophore prioritized)
# Format: comma-separated tissue:weight pairs
CHEMOSENSORY_TISSUE_WEIGHTS_STR = os.getenv('CHEMOSENSORY_TISSUE_WEIGHTS', 'rhinophore:2.0,oral-tentacle:1.0')
CHEMOSENSORY_TISSUE_WEIGHTS = {}
for pair in CHEMOSENSORY_TISSUE_WEIGHTS_STR.split(','):
    tissue, weight = pair.strip().split(':')
    CHEMOSENSORY_TISSUE_WEIGHTS[tissue.strip()] = float(weight.strip())

# NEW: Phase 2 - G-protein co-expression weight
GPROTEIN_COEXPR_WEIGHT = float(os.getenv('GPROTEIN_COEXPR_WEIGHT', 2))

# NEW: Phase 3 - ECL divergence weight
ECL_DIVERGENCE_WEIGHT = float(os.getenv('ECL_DIVERGENCE_WEIGHT', 1.5))

# NEW: Phase 5 - CAFE expansion weight
EXPANSION_WEIGHT = float(os.getenv('EXPANSION_WEIGHT', 1.5))
PREFERRED_EXPANSION_LEVELS = os.getenv('PREFERRED_EXPANSION_LEVELS', 'Aeolid-specific').split(',')

# Bead -ar8: Intra-genome tandem-cluster weight. The field's signature
# evidence type for chemoreceptor claims (Cummins 2009, Vogt 2023, Robertson
# 2015, Nath 2025). Especially important since rhinophore RNA-seq has known
# confounds (starvation + low depth, bead -qu9). Default weight set high.
TANDEM_CLUSTER_WEIGHT = float(os.getenv('TANDEM_CLUSTER_WEIGHT', 2.5))
TANDEM_CLUSTERS_FILE = os.getenv('TANDEM_CLUSTERS_FILE', '')

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
# Try legacy 2-column expression_file first; fall back to expression_summary.csv
if os.path.exists(expression_file):
    expr_data = pd.read_csv(expression_file, names=['id', 'weight'])
else:
    expr_data = pd.DataFrame(columns=['id', 'weight'])
    print(f"Note: Legacy expression file not found: {expression_file}. "
          "Will use expression_summary.csv if available.", file=sys.stderr)

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
tree_filename = os.environ.get("PHYLO_TREE_FILENAME", "all_berghia_refs.treefile")
tree_file = f"{phylo_dir}/{tree_filename}"
t = None
ref_ids = []

if not os.path.exists(tree_file):
    print(f"Warning: Tree file not found: {tree_file}. Phylogenetic scoring disabled.",
          file=sys.stderr)
elif not ETE3_AVAILABLE:
    print(f"Warning: ete3 not available. Phylogenetic scoring disabled.", file=sys.stderr)
else:
    try:
        t = Tree(tree_file, format=1)  # format=1 for IQ-TREE output with support values
        # Identify reference sequences — exclude outgroup (used only for rooting)
        ref_ids = [leaf.name for leaf in t
                   if leaf.name.startswith('ref_') and not leaf.name.startswith('ref_outgroup_')]
        outgroup_ids = [leaf.name for leaf in t if leaf.name.startswith('ref_outgroup_')]
        if outgroup_ids:
            print(f"Excluded {len(outgroup_ids)} outgroup sequence(s) from distance calculations",
                  file=sys.stderr)
    except Exception as e:
        print(f"Warning: Could not load tree: {e}. Phylogenetic scoring disabled.",
              file=sys.stderr)
        t = None
        ref_ids = []


# --- Load Reference Categories from CSV ---
# Maps ref_id -> weight based on LSE/one_to_one category
import csv as _csv
ref_category_weights = {}  # ref_id -> weight (from ref_categories_final.csv)

REF_CATEGORIES_CSV = os.getenv('REF_CATEGORIES_CSV', '')
if not REF_CATEGORIES_CSV:
    # Bead -o7k: results_dir is computed later from phylo_dir; here we must
    # resolve a temporary candidate from $RESULTS_DIR (env) or by walking up
    # from phylo_dir.
    _provisional_results = os.environ.get('RESULTS_DIR', '')
    if not _provisional_results:
        _p = Path(phylo_dir)
        while _p != _p.parent:
            if (_p / 'reference_sequences').exists() or (_p / 'ranking').exists():
                _provisional_results = str(_p)
                break
            _p = _p.parent
        if not _provisional_results:
            _provisional_results = str(Path(phylo_dir).parent.parent)
    potential = os.path.join(_provisional_results, 'reference_sequences', 'ref_categories_final.csv')
    if os.path.exists(potential):
        REF_CATEGORIES_CSV = potential

if REF_CATEGORIES_CSV and os.path.exists(REF_CATEGORIES_CSV):
    with open(REF_CATEGORIES_CSV) as _f:
        for row in _csv.DictReader(_f):
            ref_id = row['ref_id']
            category = row['category']
            if category == 'lse':
                ref_category_weights[ref_id] = CHEMORECEPTOR_REF_WEIGHT
            elif category == 'outgroup':
                pass  # Outgroup excluded from distance calculations
            else:  # one_to_one_ortholog
                ref_category_weights[ref_id] = OTHER_GPCR_REF_WEIGHT
    lse_count = sum(1 for v in ref_category_weights.values() if v == CHEMORECEPTOR_REF_WEIGHT)
    oto_count = sum(1 for v in ref_category_weights.values() if v == OTHER_GPCR_REF_WEIGHT)
    print(f"Loaded {len(ref_category_weights)} reference categories "
          f"({lse_count} LSE, {oto_count} one-to-one)",
          file=sys.stderr)


def categorize_reference(ref_name):
    """
    Categorize a reference as chemoreceptor or other GPCR.

    Returns weight for this reference type.
    Priority: (1) ref_categories_final.csv, (2) loaded JSON categories, (3) keyword heuristic.

    Delegates to the unit-tested _rank_candidates_lib.categorize_reference,
    which uses a word-boundary regex (no 'or' substring false-positives).
    See bead -wux for the previous bug.
    """
    return _categorize_reference_corrected(
        ref_name,
        chemoreceptor_weight=CHEMORECEPTOR_REF_WEIGHT,
        other_gpcr_weight=OTHER_GPCR_REF_WEIGHT,
        explicit_category_map=ref_category_weights,
        explicit_chemoreceptor_set=chemoreceptor_refs,
    )


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

                    # Read new omega_max / omega_mean / weight_at_max columns
                    # (added 2026-05 to preserve aBSREL episodic-selection
                    # signal; bead -ea9). Falls back to legacy 'omega' column.
                    omega_max = row.get('omega_max', omega)
                    omega_mean = row.get('omega_mean', omega)
                    weight_at_max = row.get('weight_at_max', float('nan'))

                    if pd.notna(omega) and branch_id:
                        entry = {
                            'branch_id': str(branch_id),
                            'omega': float(omega),                       # legacy
                            'omega_max': float(omega_max) if pd.notna(omega_max) else float(omega),
                            'omega_mean': float(omega_mean) if pd.notna(omega_mean) else float(omega),
                            'weight_at_max': float(weight_at_max) if pd.notna(weight_at_max) else float('nan'),
                            'p_value': float(p_value) if pd.notna(p_value) else 1.0,
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

    Returns list of corrected q-values in the same order. Delegates to
    statsmodels.stats.multitest.multipletests via _rank_candidates_lib.
    The previous hand-rolled implementation had a rank-indexing bug
    producing non-monotonic q-values (bead -wux).
    """
    return _bh_corrected(pvalues)


dnds_data = load_absrel_with_fdr(selective_dir)


# --- Load BUSTED-S / BUSTED-MH / MEME OG-level signals (bead -urk) ---
def load_busted_signals(selective_dir):
    """Load per-OG BUSTED-S, BUSTED-MH, and MEME aggregate CSVs.

    Returns dict: og_name -> {
        busted_s_p, busted_s_significant,
        busted_mh_p, busted_mh_significant,
        meme_n_episodic, meme_fraction_episodic, meme_n_total,
    }
    Missing files / OGs default to NaN / False.
    """
    out = {}
    for fname, key_p, key_sig in [
        ('busted_s_results.csv', 'busted_s_p', 'busted_s_significant'),
        ('busted_mh_results.csv', 'busted_mh_p', 'busted_mh_significant'),
    ]:
        path = os.path.join(selective_dir, fname)
        if not os.path.exists(path):
            continue
        try:
            df_b = pd.read_csv(path)
        except Exception as e:
            print(f"Warning: could not parse {path}: {e}", file=sys.stderr)
            continue
        for _, r in df_b.iterrows():
            og = r.get('og_name', '')
            if not og or pd.isna(og):
                continue
            entry = out.setdefault(og, {})
            entry[key_p] = float(r.get('p_value')) if pd.notna(r.get('p_value')) else float('nan')
            entry[key_sig] = bool(int(r.get('is_significant', 0)))

    meme_path = os.path.join(selective_dir, 'meme_results.csv')
    if os.path.exists(meme_path):
        try:
            df_m = pd.read_csv(meme_path)
            for _, r in df_m.iterrows():
                og = r.get('og_name', '')
                if not og or pd.isna(og):
                    continue
                entry = out.setdefault(og, {})
                entry['meme_n_total'] = int(r.get('n_sites_total', 0) or 0)
                entry['meme_n_episodic'] = int(r.get('n_sites_episodic', 0) or 0)
                entry['meme_fraction_episodic'] = float(r.get('fraction_episodic', 0.0) or 0.0)
        except Exception as e:
            print(f"Warning: could not parse {meme_path}: {e}", file=sys.stderr)

    if out:
        n_busted_s_sig = sum(1 for v in out.values() if v.get('busted_s_significant'))
        n_busted_mh_sig = sum(1 for v in out.values() if v.get('busted_mh_significant'))
        n_meme_pos = sum(1 for v in out.values() if v.get('meme_n_episodic', 0) > 0)
        print(f"Loaded BUSTED/MEME signals for {len(out)} OGs "
              f"(BUSTED-S sig: {n_busted_s_sig}, BUSTED-MH sig: {n_busted_mh_sig}, "
              f"MEME with >=1 episodic site: {n_meme_pos})", file=sys.stderr)
    return out


busted_meme_signals = load_busted_signals(selective_dir)


# --- Load CDS provenance manifest (bead -325) ---
def load_cds_provenance_map(path):
    """Return dict seq_id -> source ('native' | 'miniprot' | 'unknown')."""
    if not path or not os.path.exists(path):
        return {}
    try:
        df_prov = pd.read_csv(path)
    except Exception as e:
        print(f"Warning: could not parse CDS provenance {path}: {e}",
              file=sys.stderr)
        return {}
    out = {}
    for _, r in df_prov.iterrows():
        sid = r.get('seq_id', '')
        src = r.get('source', '')
        if sid and not pd.isna(sid):
            out[str(sid)] = str(src) if src and not pd.isna(src) else 'unknown'
    if out:
        n_miniprot = sum(1 for v in out.values() if v == 'miniprot')
        n_native = sum(1 for v in out.values() if v == 'native')
        print(f"Loaded CDS provenance for {len(out)} sequences "
              f"(native={n_native}, miniprot={n_miniprot})", file=sys.stderr)
    return out


CDS_PROVENANCE_CSV = os.environ.get(
    'CDS_PROVENANCE_CSV',
    os.path.join(os.environ.get('RESULTS_DIR', ''),
                 'reference_sequences', 'cds_provenance.csv') if os.environ.get('RESULTS_DIR') else ''
)
cds_provenance_map = load_cds_provenance_map(CDS_PROVENANCE_CSV)

# --- Load Synteny Data (Quantitative) ---
def load_synteny_scores(synteny_dir):
    """Load synteny anchor counts and produce per-gene scores in (0, 1].

    Bead -e59 fix: prefer the JCVI MCscan per-gene CSV (canonical) when
    available; fall back to MCScanX `.collinearity` files for legacy runs.
    Drops the previous synteny_ids.txt line-count input (counted any
    minimap2 hit, not actual collinearity blocks — bead -ce4/-mqt).

    Bead -ce4 fix: log-scale via _rank_candidates_lib.normalize_synteny_counts;
    returns ``None`` per gene when max anchor count is degenerate (<5),
    so the composite-score evidence-completeness multiplier handles it
    correctly (no spurious 1.0 scores for trivial mapping hits).
    """
    counts: dict = {}

    # Preferred: JCVI per-gene anchors CSV (produced by 06_synteny_and_mapping.sh)
    jcvi_csv_env = os.getenv("SYNTENY_ANCHORS_CSV", "")
    candidate_csvs = []
    if jcvi_csv_env:
        candidate_csvs.append(jcvi_csv_env)
    if synteny_dir:
        candidate_csvs.append(os.path.join(synteny_dir, "jcvi_anchors_per_gene.csv"))

    for csv_path in candidate_csvs:
        if not csv_path or not os.path.exists(csv_path):
            continue
        try:
            df_anch = pd.read_csv(csv_path)
        except Exception as e:
            print(f"Warning: Could not parse JCVI anchors CSV {csv_path}: {e}",
                  file=sys.stderr)
            continue
        if "candidate_id" not in df_anch.columns:
            continue
        # Aggregate over all target species (per-Berghia-gene n_anchor_blocks
        # summed across targets if multiple). Use total_anchor_genes as the
        # raw count proxy.
        agg_col = "total_anchor_genes" if "total_anchor_genes" in df_anch.columns \
            else "n_anchor_blocks"
        for _, row in df_anch.iterrows():
            gid = row.get("candidate_id")
            if not gid or pd.isna(gid):
                continue
            v = row.get(agg_col, 0)
            if pd.isna(v):
                continue
            counts[gid] = counts.get(gid, 0) + int(v)
        if counts:
            print(f"Loaded JCVI synteny anchors for {len(counts)} candidates "
                  f"from {csv_path}", file=sys.stderr)
            break  # first non-empty CSV wins

    # Legacy fallback: MCScanX .collinearity files
    if not counts and synteny_dir:
        for collin_file in Path(synteny_dir).glob('*.collinearity'):
            try:
                with open(collin_file) as f:
                    for line in f:
                        if line.startswith('#') or not line.strip():
                            continue
                        parts = line.strip().split()
                        if len(parts) >= 3:
                            gene1, gene2 = parts[1], parts[2]
                            counts[gene1] = counts.get(gene1, 0) + 1
                            counts[gene2] = counts.get(gene2, 0) + 1
            except Exception as e:
                print(f"Warning: Could not parse {collin_file}: {e}", file=sys.stderr)

    if not counts:
        return {}

    # Use the unit-tested log-scale normalizer; returns None per gene when
    # the dataset is degenerate (max < 5 anchors).
    normalized = normalize_synteny_counts(counts, min_max_anchors=5)
    # Filter out None values; the composite-score function treats absence
    # as "no synteny data" via has_synteny_data flag, so we only return
    # genes with concrete scores here.
    return {k: v for k, v in normalized.items() if v is not None}


synteny_scores = load_synteny_scores(synteny_dir)


# --- Tandem-cluster scoring (bead -ar8) ---
def load_tandem_cluster_data(path):
    """Load tandem-cluster CSV produced by compute_tandem_clusters.py.

    Returns dict: candidate_id -> (cluster_size, cluster_id) (None if missing).
    Also computes a per-candidate score in [0, 1] = log1p(size)/log1p(20).
    """
    if not path or not os.path.exists(path):
        return {}, {}
    try:
        df_tc = pd.read_csv(path)
    except Exception as e:
        print(f"Warning: could not read tandem-cluster CSV at {path}: {e}",
              file=sys.stderr)
        return {}, {}
    info = {}
    scores = {}
    for _, r in df_tc.iterrows():
        cid = r.get('candidate_id', '')
        if not cid or pd.isna(cid):
            continue
        size = r.get('tandem_cluster_size')
        cluster_id = r.get('tandem_cluster_id', '') or None
        if pd.isna(size):
            info[cid] = (None, cluster_id)
            continue
        info[cid] = (int(size), cluster_id if cluster_id else None)
        # Score: log1p-scale capped at cluster size 20
        scores[cid] = float(np.log1p(int(size)) / np.log1p(20.0))
    return info, scores


tandem_cluster_info, tandem_cluster_scores = load_tandem_cluster_data(TANDEM_CLUSTERS_FILE)
if tandem_cluster_info:
    print(f"Loaded tandem-cluster info for {len(tandem_cluster_info)} candidates "
          f"({sum(1 for v in tandem_cluster_info.values() if v[0] and v[0] >= 3)} "
          f"in clusters of size >=3)", file=sys.stderr)


# --- Derive results directory from phylo_dir ---
# Derive results_dir: walk up from phylo_dir until we find a known results structure
# Handles both {RESULTS_DIR}/phylogenies/protein and .../protein/v2 etc.
results_dir = os.environ.get("RESULTS_DIR", "")
if not results_dir:
    # Fallback: walk up from phylo_dir looking for ranking/ or chemogpcrs/ siblings
    _p = Path(phylo_dir)
    while _p != _p.parent:
        if (_p / "ranking").exists() or (_p / "chemogpcrs").exists():
            results_dir = str(_p)
            break
        _p = _p.parent
    if not results_dir:
        results_dir = str(Path(phylo_dir).parent.parent)

# --- Scoring Functions ---

def path_bootstrap_confidence(tree, name_a, name_b, threshold):
    """
    Calculate a confidence multiplier for the path between two nodes
    based on bootstrap support values of internal nodes along the path.

    Returns a value in [0, 1]:
    - 1.0 if all branches on the path have support >= threshold (or no support data)
    - Reduced proportionally for each low-support branch
    - 0.0 if the path cannot be determined
    """
    try:
        ancestor = tree.get_common_ancestor(name_a, name_b)
    except Exception:
        return 0.0

    # Walk from each leaf up to the common ancestor, collecting support values
    supports = []
    for name in (name_a, name_b):
        try:
            nodes = tree.search_nodes(name=name)
            if not nodes:
                continue
            current = nodes[0]
            while current != ancestor and current.up:
                support = getattr(current.up, 'support', None)
                if support is not None:
                    supports.append(support)
                current = current.up
        except Exception:
            continue

    if not supports:
        return 1.0  # No bootstrap data available, don't penalize

    # Fraction of branches meeting threshold
    good = sum(1 for s in supports if s >= threshold)
    return good / len(supports)


# --- Efficient Phylo Scoring ---
# Precompute k-nearest distances and cache bootstrap values to avoid
# O(candidates * refs) tree traversals per candidate.

K_NEAREST_REFS = int(os.getenv('K_NEAREST_REFS', 50))


def precompute_phylo_data(tree, candidate_names, ref_ids, k_nearest=50):
    """Precompute distance matrix and bootstrap data for efficient scoring.

    Returns:
        distances: dict of candidate_name -> list of (distance, ref_name) sorted ascending, truncated to k
        node_supports: dict of node_id -> support value for all internal nodes
        leaf_lookup: dict of leaf_name -> ete3 node
    """
    # Cache all internal node support values by node id
    node_supports = {}
    for node in tree.traverse():
        if not node.is_leaf():
            node_supports[id(node)] = getattr(node, 'support', None)

    # Build leaf name -> node lookup once
    leaf_lookup = {leaf.name: leaf for leaf in tree}

    # Precompute k-nearest ref distances for each candidate
    distances = {}
    for cand_name in candidate_names:
        if cand_name not in leaf_lookup:
            distances[cand_name] = []
            continue

        ref_dists = []
        for ref_name in ref_ids:
            if ref_name not in leaf_lookup:
                continue
            try:
                dist = tree.get_distance(cand_name, ref_name)
                ref_dists.append((dist, ref_name))
            except Exception:
                continue

        # Sort and keep only k nearest
        ref_dists.sort(key=lambda x: x[0])
        distances[cand_name] = ref_dists[:k_nearest]

    return distances, node_supports, leaf_lookup


def fast_path_bootstrap_confidence(leaf_a, leaf_b, leaf_lookup, node_supports, threshold):
    """Fast bootstrap confidence using cached node supports."""
    if leaf_a not in leaf_lookup or leaf_b not in leaf_lookup:
        return 0.0

    try:
        node_a = leaf_lookup[leaf_a]
        node_b = leaf_lookup[leaf_b]
        ancestor = node_a.get_common_ancestor(node_b)
    except Exception:
        return 0.0

    supports = []
    for node in (node_a, node_b):
        current = node
        while current != ancestor and current.up:
            sup = node_supports.get(id(current.up))
            if sup is not None:
                supports.append(sup)
            current = current.up

    if not supports:
        return 1.0

    good = sum(1 for s in supports if s >= threshold)
    return good / len(supports)


def weighted_distance_to_refs(node_name, tree, ref_ids):
    """
    Calculate weighted phylogenetic distance to references.

    Chemoreceptor references are weighted higher than other GPCRs.
    Paths with low bootstrap support are penalized proportionally.
    Returns inverse weighted distance (higher = closer to important refs).

    NOTE: This is the legacy function used when precomputed data is not available.
    The fast path uses weighted_distance_to_refs_fast() instead.
    """
    # Guard for when tree is not available
    if tree is None or not ref_ids:
        return 0.0

    if node_name in ref_ids:
        return float('inf')  # Reference itself

    try:
        total_weighted_inverse = 0.0
        for ref in ref_ids:
            try:
                distance = tree.get_distance(node_name, ref)
                weight = categorize_reference(ref)
                # Penalize paths with low bootstrap support
                confidence = path_bootstrap_confidence(
                    tree, node_name, ref, BOOTSTRAP_THRESHOLD
                )
                # Weighted inverse distance scaled by path confidence
                total_weighted_inverse += confidence * weight / (distance + 1e-6)
            except Exception:
                continue

        return total_weighted_inverse if total_weighted_inverse > 0 else 0.0
    except Exception:
        return 0.0


def weighted_distance_to_refs_fast(node_name, precomputed_distances, leaf_lookup,
                                    node_supports):
    """Efficient phylo scoring using precomputed k-nearest distances.

    Uses cached distances and bootstrap values instead of per-pair tree traversal.
    """
    if node_name not in precomputed_distances:
        return 0.0

    nearest = precomputed_distances[node_name]
    if not nearest:
        return 0.0

    total = 0.0
    for dist, ref_name in nearest:
        weight = categorize_reference(ref_name)
        confidence = fast_path_bootstrap_confidence(
            node_name, ref_name, leaf_lookup, node_supports, BOOTSTRAP_THRESHOLD
        )
        total += confidence * weight / (dist + 1e-6)

    return total


def get_selection_scores(candidate_id, dnds_data):
    """
    Calculate separate purifying and positive selection scores.

    Returns: (purifying_score, positive_score, is_significant)

    Biological interpretation for chemoreceptor identification (bead -ea9):
    - Default PURIFYING_WEIGHT=0 (set in config.sh) because chemoreceptor
      discovery rewards diversifying selection on extracellular loops, NOT
      whole-gene purifying selection. A conserved housekeeping GPCR with
      ω≪1 should not rank highly for the chemoreceptor question.
    - Strong positive selection (ω≫1) is the chemoreceptor-relevant signal.
    - Set PURIFYING_WEIGHT>0 only when looking for conserved-function GPCRs.

    Uses omega_max from the aBSREL rate-class mixture (preserves the episodic
    positive-selection signal that aBSREL detects), falling back to legacy
    'omega' (weighted-mean) if omega_max is not available.
    """
    if candidate_id not in dnds_data:
        return 0.0, 0.0, False

    entry = dnds_data[candidate_id]
    # Prefer omega_max (max across rate classes) to preserve episodic signal.
    # Falls back to 'omega' for legacy aBSREL CSVs without the new column.
    omega = entry.get('omega_max', entry.get('omega', float('nan')))
    is_significant = entry.get('significant', False)
    corrected_p = entry.get('corrected_p', 1.0 if not is_significant else 0.0)

    # Clamp to avoid extreme log values from edge-case mixture parameters.
    if not (isinstance(omega, float) and np.isnan(omega)):
        omega = max(min(float(omega), 100.0), 0.001)

    # Pass weights=1 to get the unweighted raw signal (with the significance
    # boost still applied). The caller in calculate_rank_score multiplies by
    # PURIFYING_WEIGHT / POSITIVE_WEIGHT exactly once. Note: lib uses log10
    # (more standard for ω); previous code used natural log. If you have
    # tuned weights from before this change, multiply them by ~2.3 (=ln 10)
    # to preserve absolute scale, OR retune. With PURIFYING_WEIGHT=0 (new
    # default), purifying contributes nothing regardless.
    result = _get_selection_scores_corrected(
        omega=omega,
        p_corrected=corrected_p,
        purifying_weight=1.0,
        positive_weight=1.0,
        significance_threshold=ABSREL_FDR_THRESHOLD,
        significance_boost=1.5,
    )
    return result["purifying_score"], result["positive_score"], result["is_significant"]


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


# --- Orthogroup Bootstrap Confidence Score ---
OG_CONFIDENCE_WEIGHT = float(os.getenv('OG_CONFIDENCE_WEIGHT', 1))


def load_og_trees(results_dir):
    """Load OrthoFinder resolved gene trees.

    Returns dict: og_name -> ete3.Tree
    """
    if not ETE3_AVAILABLE:
        return {}

    og_trees = {}
    og_tree_dirs = list(Path(results_dir).glob(
        'orthogroups/input/OrthoFinder/Results_*/Resolved_Gene_Trees'
    ))

    for tree_dir in og_tree_dirs:
        for tree_file in tree_dir.glob('OG*_tree.txt'):
            og_name = tree_file.stem.replace('_tree', '')
            try:
                og_trees[og_name] = Tree(str(tree_file), format=1)
            except Exception:
                try:
                    og_trees[og_name] = Tree(str(tree_file), format=0)
                except Exception:
                    continue

    return og_trees


def get_og_confidence_score(candidate_id, gene_to_og, og_trees, bootstrap_threshold):
    """Compute bootstrap confidence from orthogroup gene tree.

    Finds nearest reference in OG tree, returns fraction of path nodes
    with bootstrap >= threshold.

    Returns: (score, has_data)
    """
    og = gene_to_og.get(candidate_id)
    if not og or og not in og_trees:
        return 0.0, False

    tree = og_trees[og]
    leaves = [l.name for l in tree]

    # Find this candidate's leaf (Berghia IDs in OG trees use full OrthoFinder
    # format). Bead -mqt: previous version used substring match (`candidate_id in l`)
    # which incorrectly matched any prefix-collision (e.g. 'TRINITY_DN1' would
    # match 'TRINITY_DN10'). Use exact match plus underscore-separated prefix.
    cand_leaves = [l for l in leaves
                   if l == candidate_id or l.startswith(candidate_id + '_')]
    if not cand_leaves:
        return 0.0, False
    cand_leaf = cand_leaves[0]

    # Find reference leaves (exclude outgroup)
    ref_leaves = [l for l in leaves
                  if (l.startswith('ref_lse_') or l.startswith('ref_one_to_one_ortholog_'))]
    if not ref_leaves:
        return 0.0, False

    # Find nearest reference by distance
    best_dist = float('inf')
    best_ref = None
    for ref_leaf in ref_leaves:
        try:
            dist = tree.get_distance(cand_leaf, ref_leaf)
            if dist < best_dist:
                best_dist = dist
                best_ref = ref_leaf
        except Exception:
            continue

    if best_ref is None:
        return 0.0, False

    # Compute bootstrap confidence along path to nearest ref
    try:
        ancestor = tree.get_common_ancestor(cand_leaf, best_ref)
    except Exception:
        return 0.0, False

    supports = []
    for name in (cand_leaf, best_ref):
        nodes = tree.search_nodes(name=name)
        if not nodes:
            continue
        current = nodes[0]
        while current != ancestor and current.up:
            sup = getattr(current.up, 'support', None)
            if sup is not None:
                supports.append(sup)
            current = current.up

    if not supports:
        return 1.0, True  # No bootstrap data — don't penalize

    good = sum(1 for s in supports if s >= bootstrap_threshold)
    return good / len(supports), True


# --- NEW: Phase 1 - Chemosensory Expression Scoring ---

def load_chemosensory_expression(results_dir):
    """
    Load chemosensory expression summary from process_expression.py output.

    Returns dict: gene_id -> {tau_index, enrichment, is_specific, tissue_tpms, has_data}
    """
    expr_file = os.path.join(results_dir, 'expression', 'expression_summary.csv')
    expr_data = {}

    if not os.path.exists(expr_file):
        return expr_data

    try:
        df = pd.read_csv(expr_file)
        # Detect per-tissue TPM columns (end with _tpm)
        tpm_cols = [c for c in df.columns if c.endswith('_tpm')]
        for _, row in df.iterrows():
            gene_id = str(row['gene_id'])
            tissue_tpms = {}
            for col in tpm_cols:
                tissue_name = col.replace('_tpm', '')
                val = row.get(col, 0.0)
                tissue_tpms[tissue_name] = val if pd.notna(val) else 0.0
            expr_data[gene_id] = {
                'tau_index': row.get('tau_index', 0.0),
                'chemosensory_enrichment': row.get('chemosensory_enrichment', 0.0),
                'mean_chemosensory_tpm': row.get('mean_chemosensory_tpm', 0.0),
                'expressed_in_chemosensory': row.get('expressed_in_chemosensory', False),
                'is_chemosensory_specific': row.get('is_chemosensory_specific', False),
                'tissue_tpms': tissue_tpms,
                'has_data': True
            }
        print(f"Loaded chemosensory expression data for {len(expr_data)} genes", file=sys.stderr)
    except Exception as e:
        print(f"Warning: Could not load chemosensory expression: {e}", file=sys.stderr)

    return expr_data


def get_chemosensory_expression_score(candidate_id, chemo_expr_data,
                                      tissue_weights=None):
    """
    Calculate chemosensory expression score with tissue-specific weighting.

    Components:
    1. Tissue-weighted chemosensory expression (rhinophore > oral-tentacle)
    2. Tau index (tissue specificity)
    3. Chemosensory enrichment (fold-change vs other tissues)

    Args:
        tissue_weights: dict of tissue_name -> weight (e.g. rhinophore:2.0, oral-tentacle:1.0)

    Returns: (score, has_data)
    """
    if candidate_id not in chemo_expr_data:
        return 0.0, False

    data = chemo_expr_data[candidate_id]

    if not data.get('has_data', False):
        return 0.0, False

    if tissue_weights is None:
        tissue_weights = {'rhinophore': 2.0, 'oral-tentacle': 1.0}

    score = 0.0
    tissue_tpms = data.get('tissue_tpms', {})

    # Component 1: Tissue-weighted chemosensory expression (0-0.4)
    # Weight each chemosensory tissue's contribution by its importance
    weighted_tpm = 0.0
    total_weight = 0.0
    for tissue, weight in tissue_weights.items():
        tpm = tissue_tpms.get(tissue, 0.0)
        weighted_tpm += tpm * weight
        total_weight += weight

    if total_weight > 0 and weighted_tpm > 0:
        # Log-scale weighted TPM, normalize to 0-0.4 range
        # log2(1 + weighted_tpm/total_weight) capped at ~7 (128 TPM)
        norm_tpm = weighted_tpm / total_weight
        log_tpm = min(np.log2(1 + norm_tpm) / 7, 1.0)
        score += log_tpm * 0.4

    # Component 2: Tissue specificity (tau index, 0-0.3)
    tau = data.get('tau_index', 0.0)
    if pd.notna(tau):
        score += tau * 0.3

    # Component 3: Chemosensory enrichment (log fold-change, 0-0.3)
    enrichment = data.get('chemosensory_enrichment', 0.0)
    if enrichment > 1:
        log_enrichment = min(np.log2(enrichment) / 5, 1.0)
        score += log_enrichment * 0.3

    # Bonus for chemosensory-specific expression (tau >= threshold AND max in chemo tissue)
    if data.get('is_chemosensory_specific', False):
        # Extra bonus if rhinophore is the dominant tissue
        rhino_tpm = tissue_tpms.get('rhinophore', 0.0)
        other_chemo_tpms = [tissue_tpms.get(t, 0.0) for t in tissue_weights if t != 'rhinophore']
        max_other_chemo = max(other_chemo_tpms) if other_chemo_tpms else 0.0

        if rhino_tpm > 0 and rhino_tpm >= max_other_chemo:
            score *= 1.75  # 75% bonus for rhinophore-dominant
        else:
            score *= 1.5   # 50% bonus for other chemosensory-specific

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


def load_og_dnds_reliability(results_dir):
    """Compute per-OG dN/dS-axis reliability weights in [0, 1].

    Bead -8st (2026-05-08). For each orthogroup, count how many of its
    members have a CDS in the merged reference CDS file. Weight is
    ``min(1.0, n_ref_cds / DNDS_RELIABILITY_FULL)``. Below the threshold
    the dN/dS axis still contributes proportionally; above it the axis
    is fully unlocked. The default threshold (5) matches the "medium"
    reliability tier in scripts/add_og_coverage_columns.py.

    Returns dict: og_id -> float in [0, 1]. Empty dict if either input
    is missing — callers must default to 0.0 in that case so the dN/dS
    axis falls out of fair-scoring rather than getting full weight on
    unverified data.
    """
    cds_path = os.path.join(results_dir, 'reference_sequences', 'cds',
                            'all_references_cds.fna')
    og_tsv = None
    for candidate in (
            os.path.join(results_dir, 'orthogroups', 'OrthoFinder'),
            os.path.join(results_dir, 'orthofinder')):
        if not os.path.isdir(candidate):
            continue
        for root, _, files in os.walk(candidate):
            if 'Orthogroups.tsv' in files and root.endswith('Orthogroups'):
                og_tsv = os.path.join(root, 'Orthogroups.tsv')
                break
        if og_tsv:
            break

    if not (os.path.exists(cds_path) and og_tsv):
        return {}

    # Reuse the helpers from add_og_coverage_columns.py — same data sources,
    # same parsing, single source of truth for "what counts as a recovered
    # reference CDS for an OG".
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from add_og_coverage_columns import parse_cds_ids, load_og_members
    except ImportError:
        return {}

    cds_ids = parse_cds_ids(cds_path)
    og_members = load_og_members(og_tsv)
    weights: dict = {}
    full = max(DNDS_RELIABILITY_FULL, 1.0)
    for og_id, members in og_members.items():
        n_ref_cds = sum(1 for m in members if m in cds_ids)
        weights[og_id] = min(1.0, n_ref_cds / full)
    print(f"Computed dN/dS reliability weights for {len(weights)} orthogroups "
          f"(full-credit threshold = {DNDS_RELIABILITY_FULL} ref CDS)",
          file=sys.stderr)
    return weights


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
        # Bead -8st: per-OG dN/dS reliability multiplier in [0, 1].
        # Falls back to 1.0 if the column wasn't populated (e.g. legacy
        # CSVs run through downstream consumers) so we don't silently
        # zero out the dN/dS axis on inputs that predate the multiplier.
        dnds_rw = row.get('dnds_reliability_weight', 1.0)
        try:
            dnds_rw = float(dnds_rw)
        except (TypeError, ValueError):
            dnds_rw = 1.0

        # Base scores always available
        score = (
            row['phylo_score_norm'] * weights.get('phylo', 0) +
            row['purifying_score_norm'] * weights.get('purifying', 0) * dnds_rw +
            row['positive_score_norm'] * weights.get('positive', 0) * dnds_rw +
            row['lse_depth_score_norm'] * weights.get('lse_depth', 0)
        )
        total_weight = (
            weights.get('phylo', 0) +
            weights.get('purifying', 0) * dnds_rw +
            weights.get('positive', 0) * dnds_rw +
            weights.get('lse_depth', 0)
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
        if row.get('has_og_confidence_data', False):
            score += row.get('og_confidence_score_norm', 0) * weights.get('og_confidence', 0)
            total_weight += weights.get('og_confidence', 0)

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
    # Guard for when tree is not available
    if tree is None or not ref_ids:
        print("Warning: Cross-validation requires phylogenetic tree. Skipping.",
              file=sys.stderr)
        return None

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

# Precompute phylo data for efficient scoring
precomputed_distances = {}
node_supports = {}
leaf_lookup = {}
if t is not None and ref_ids:
    print(f"Precomputing phylo data: {len(candidates)} candidates x {len(ref_ids)} refs "
          f"(k={K_NEAREST_REFS})...", file=sys.stderr)
    precomputed_distances, node_supports, leaf_lookup = precompute_phylo_data(
        t, list(candidates['id']), ref_ids, k_nearest=K_NEAREST_REFS
    )
    print(f"Precomputation complete: {len(precomputed_distances)} candidates indexed",
          file=sys.stderr)

# First pass: collect all depths for percentile calculation
# Only include nodes on paths with sufficient bootstrap support (H12 fix)
all_depths = []
if t is not None:
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
else:
    print("Warning: Skipping depth calculation - tree not available.", file=sys.stderr)

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

# Bead -8st: per-OG dN/dS-axis reliability weights. Used downstream to
# scale the purifying/positive contributions in the composite so that
# under-supported OGs (chemoreceptor LSE expansions where miniprot
# typically fails on most paralogs) don't get full dN/dS weight on
# what is statistically just noise.
og_dnds_reliability = load_og_dnds_reliability(results_dir)

# Orthogroup gene trees for bootstrap confidence scoring
og_trees = {}
if OG_CONFIDENCE_WEIGHT > 0 and gene_to_orthogroup:
    og_trees = load_og_trees(results_dir)
    if og_trees:
        print(f"Loaded {len(og_trees)} orthogroup gene trees for confidence scoring",
              file=sys.stderr)
    else:
        print("No orthogroup gene trees found for confidence scoring", file=sys.stderr)

# ECL divergence threshold from environment
ECL_RATIO_THRESHOLD = float(os.getenv('ECL_CONSERVATION_RATIO_THRESHOLD', 2.0))

# Second pass: calculate all scores
results = []

for cand_id in candidates['id']:
    # Phylogenetic score — use fast precomputed path when available
    if precomputed_distances:
        phylo_score = weighted_distance_to_refs_fast(
            cand_id, precomputed_distances, leaf_lookup, node_supports
        )
    else:
        phylo_score = weighted_distance_to_refs(cand_id, t, ref_ids)

    # Selection scores (separate purifying and positive)
    purifying_score, positive_score, selection_significant = get_selection_scores(cand_id, dnds_data)

    # Bead -urk: BUSTED-S/MH gene-wide signal + MEME site-level signal for
    # this candidate's OG. These are SUPPLEMENTARY to the per-branch aBSREL
    # signal — they don't get their own axis weight (which would compound
    # with positive_score), but a small multiplicative boost (up to 1.3x) is
    # applied to positive_score when the gene-wide BUSTED tests independently
    # confirm aBSREL's branch-level signal. This is the "two-test agreement"
    # heuristic — credible only when both tests fire.
    cand_og = gene_to_orthogroup.get(cand_id, None)
    # Bead -325: CDS provenance — Berghia gets 'native' (RefSeq genome
    # annotation); reference candidates get whatever the manifest says.
    if cds_provenance_map:
        cds_source = cds_provenance_map.get(cand_id, 'native')
    else:
        cds_source = 'unknown'
    busted_s_p = float('nan')
    busted_mh_p = float('nan')
    busted_s_sig = False
    busted_mh_sig = False
    meme_n_episodic = 0
    meme_fraction_episodic = 0.0
    if cand_og and cand_og in busted_meme_signals:
        sig = busted_meme_signals[cand_og]
        busted_s_p = float(sig.get('busted_s_p', float('nan')))
        busted_mh_p = float(sig.get('busted_mh_p', float('nan')))
        busted_s_sig = bool(sig.get('busted_s_significant', False))
        busted_mh_sig = bool(sig.get('busted_mh_significant', False))
        meme_n_episodic = int(sig.get('meme_n_episodic', 0))
        meme_fraction_episodic = float(sig.get('meme_fraction_episodic', 0.0))
        # Boost: BUSTED-MH is the strictest test; positive selection that
        # survives MH correction is the most credible. aBSREL alone can be
        # tricked by short branches; BUSTED + MEME confirmation is the
        # paper-publishable signal.
        if selection_significant and busted_mh_sig:
            positive_score *= 1.3
        elif selection_significant and busted_s_sig:
            positive_score *= 1.15

    # Synteny score (quantitative)
    synteny_score, has_synteny = get_synteny_score(cand_id, synteny_scores)

    # Expression score (from legacy file or expression_summary.csv)
    expr_score, has_expression = get_expression_score(cand_id, expr_data)
    # If legacy expression file was missing, derive from expression_summary.csv
    if not has_expression and cand_id in chemo_expr_data:
        tissue_tpms = chemo_expr_data[cand_id].get('tissue_tpms', {})
        mean_tpm = np.mean(list(tissue_tpms.values())) if tissue_tpms else 0.0
        if mean_tpm > 0:
            expr_score = np.log2(1 + mean_tpm)
            has_expression = True

    # LSE depth score
    lse_score, raw_depth = get_lse_depth_score(cand_id, t, lse_threshold)

    # Phase 1 - Chemosensory expression score (tissue-weighted: rhinophore > oral-tentacle)
    chemo_expr_score, has_chemo_expr = get_chemosensory_expression_score(
        cand_id, chemo_expr_data, tissue_weights=CHEMOSENSORY_TISSUE_WEIGHTS)

    # NEW: Phase 2 - G-protein co-expression score
    gprotein_score, has_gprotein = get_gprotein_coexpr_score(cand_id, gprotein_coexpr_data)

    # NEW: Phase 3 - ECL divergence score
    ecl_score, has_ecl = get_ecl_divergence_score(cand_id, ecl_divergence_data, ECL_RATIO_THRESHOLD)

    # NEW: Phase 5 - CAFE expansion score
    expansion_score, has_expansion = get_expansion_score(
        cand_id, cafe_expansion_data, gene_to_orthogroup, PREFERRED_EXPANSION_LEVELS
    )

    # Orthogroup bootstrap confidence score
    og_conf, has_og_data = get_og_confidence_score(
        cand_id, gene_to_orthogroup, og_trees, BOOTSTRAP_THRESHOLD
    )

    # Bead -ar8: tandem-cluster score
    tandem_score = tandem_cluster_scores.get(cand_id, 0.0)
    has_tandem = cand_id in tandem_cluster_info
    tandem_size, tandem_cluster_label = tandem_cluster_info.get(cand_id, (None, None))

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
        has_expansion,
        has_og_data,
        has_tandem,
    ]
    evidence_count = sum(evidence_sources)
    evidence_completeness = evidence_count / len(evidence_sources)

    results.append({
        'id': cand_id,
        'phylo_score': phylo_score,
        'purifying_score': purifying_score,
        'positive_score': positive_score,
        'selection_significant': selection_significant,
        # Bead -urk: BUSTED/MEME diagnostic columns (read from selection-stack
        # outputs; used to boost positive_score on two-test agreement).
        'busted_s_p': busted_s_p,
        'busted_s_significant': busted_s_sig,
        'busted_mh_p': busted_mh_p,
        'busted_mh_significant': busted_mh_sig,
        'meme_n_episodic_sites': meme_n_episodic,
        'meme_fraction_episodic_sites': meme_fraction_episodic,
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
        'og_confidence_score': og_conf,
        'has_og_confidence_data': has_og_data,
        # Bead -ar8: tandem cluster
        'tandem_cluster_score': tandem_score,
        'tandem_cluster_size': tandem_size if tandem_size is not None else 0,
        'tandem_cluster_id': tandem_cluster_label or '',
        'has_tandem_cluster_data': has_tandem,
        # Bead -325: CDS provenance
        'cds_source': cds_source,
        # Bead -8st: dN/dS-axis reliability weight derived from how many
        # of this OG's members have a CDS in the reference set (linear
        # ramp to full credit at DNDS_RELIABILITY_FULL refs). Applied as
        # a multiplier on the dN/dS contribution in calculate_rank_score.
        # 0.0 when no OG can be resolved or no reliability map was loaded
        # — that takes dN/dS out of fair-scoring entirely for that row.
        'dnds_reliability_weight': (
            og_dnds_reliability.get(cand_og, 0.0) if cand_og else 0.0),
        'evidence_completeness': evidence_completeness
    })

df = pd.DataFrame(results)

# --- Normalize Scores ---
# Normalize continuous scores to [0,1] range for fair comparison
normalize_cols = [
    'phylo_score', 'purifying_score', 'positive_score', 'expression_score', 'lse_depth_score',
    # NEW score columns
    'chemosensory_expr_score', 'gprotein_coexpr_score', 'ecl_divergence_score', 'expansion_score',
    'og_confidence_score'
]

for col in normalize_cols:
    if col not in df.columns:
        df[col] = 0.0
    col_min = df[col].min()
    col_max = df[col].max()
    col_range = col_max - col_min
    if col_range > 0:
        df[f'{col}_norm'] = (df[col] - col_min) / col_range
    else:
        df[f'{col}_norm'] = 0.0

# Synteny is already 0-1
if 'synteny_score' not in df.columns:
    df['synteny_score'] = 0.0
df['synteny_score_norm'] = df['synteny_score']

# Tandem-cluster score is already 0-1 (log1p-scaled in load_tandem_cluster_data)
if 'tandem_cluster_score' not in df.columns:
    df['tandem_cluster_score'] = 0.0
df['tandem_cluster_score_norm'] = df['tandem_cluster_score']

# --- Calculate Final Rank Score with Fair Missing Data Handling ---
# For candidates missing synteny or expression data, we normalize by available weights
# This prevents penalizing candidates from species without genome assemblies


def calculate_fair_rank_score(row):
    """
    Calculate rank score with per-candidate weight normalization.

    Bead -ce4 fix: previously this function multiplied
    (sum_weighted / available_weight) * max_possible_weight, which let a
    candidate with one strong signal rank equal to a candidate with all
    signals at the same per-axis score. The corrected behavior multiplies
    by evidence_completeness (the fraction of weight that was available),
    floored at 0.4 so a single very strong signal isn't zeroed.

    Includes all Phase 1-5 scoring components. Missing data is signaled
    by has_*_data flags being False; the corresponding axis is then
    excluded from BOTH the numerator and the available-weight denominator.
    """
    # Build scores/weights dicts; missing axes set to None so the lib's
    # evidence-completeness multiplier applies correctly.
    scores = {
        'phylo': row.get('phylo_score_norm'),
        'purifying': row.get('purifying_score_norm'),
        'positive': row.get('positive_score_norm'),
        'lse_depth': row.get('lse_depth_score_norm'),
        'synteny': row.get('synteny_score_norm') if row.get('has_synteny_data') else None,
        'expression': row.get('expression_score_norm') if row.get('has_expression_data') else None,
        'chemosensory_expr': row.get('chemosensory_expr_score_norm') if row.get('has_chemosensory_expr_data') else None,
        'gprotein_coexpr': row.get('gprotein_coexpr_score_norm') if row.get('has_gprotein_data') else None,
        'ecl_divergence': row.get('ecl_divergence_score_norm') if row.get('has_ecl_data') else None,
        'expansion': row.get('expansion_score_norm') if row.get('has_expansion_data') else None,
        'og_confidence': row.get('og_confidence_score_norm') if row.get('has_og_confidence_data') else None,
        # Bead -ar8: tandem-cluster score (intra-genome paralog clustering)
        'tandem_cluster': row.get('tandem_cluster_score_norm') if row.get('has_tandem_cluster_data') else None,
    }
    weights = {
        'phylo': PHYLO_WEIGHT,
        'purifying': PURIFYING_WEIGHT,
        'positive': POSITIVE_WEIGHT,
        'lse_depth': LSE_DEPTH_WEIGHT,
        'synteny': SYNTENY_WEIGHT,
        'expression': EXPR_WEIGHT,
        'chemosensory_expr': CHEMOSENSORY_EXPR_WEIGHT,
        'gprotein_coexpr': GPROTEIN_COEXPR_WEIGHT,
        'ecl_divergence': ECL_DIVERGENCE_WEIGHT,
        'expansion': EXPANSION_WEIGHT,
        'og_confidence': OG_CONFIDENCE_WEIGHT,
        'tandem_cluster': TANDEM_CLUSTER_WEIGHT,
    }
    return _calculate_fair_rank_score_corrected(scores, weights, completeness_floor=0.4)


df['rank_score'] = df.apply(calculate_fair_rank_score, axis=1, result_type='reduce')

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
        GPROTEIN_COEXPR_WEIGHT + ECL_DIVERGENCE_WEIGHT + EXPANSION_WEIGHT +
        OG_CONFIDENCE_WEIGHT
    )
    score_pct = score / max_possible if max_possible > 0 else 0

    if completeness >= 0.6 and score_pct >= 0.5:
        return 'High'
    elif completeness >= 0.4 or score_pct >= 0.3:
        return 'Medium'
    else:
        return 'Low'


df['confidence_tier'] = df.apply(assign_confidence_tier, axis=1, result_type='reduce')

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
    # Bead -urk: BUSTED/MEME diagnostic columns
    'busted_s_p', 'busted_s_significant', 'busted_mh_p', 'busted_mh_significant',
    'meme_n_episodic_sites', 'meme_fraction_episodic_sites',
    'synteny_score', 'has_synteny_data', 'expression_score', 'has_expression_data',
    'lse_depth_score', 'raw_tree_depth',
    # NEW: Phase 1-5 scores
    'chemosensory_expr_score', 'has_chemosensory_expr_data',
    'gprotein_coexpr_score', 'has_gprotein_data',
    'ecl_divergence_score', 'has_ecl_data',
    'expansion_score', 'has_expansion_data',
    'og_confidence_score', 'has_og_confidence_data',
    # Bead -ar8: tandem cluster columns
    'tandem_cluster_score', 'tandem_cluster_size', 'tandem_cluster_id',
    'has_tandem_cluster_data',
    # Bead -325: CDS provenance (native vs miniprot-recovered)
    'cds_source',
]

# Ensure all output columns exist (fill missing with defaults)
for col in output_cols:
    if col not in df_sorted.columns:
        if col.startswith('has_'):
            df_sorted[col] = False
        else:
            df_sorted[col] = 0.0

df_sorted[output_cols].to_csv(output_file, index=False)

# --- Run Sensitivity Analysis ---
sensitivity_results = None
if RUN_SENSITIVITY:
    print(f"\nRunning sensitivity analysis ({SENSITIVITY_ITERATIONS} iterations)...", file=sys.stderr)

    # Bead -j6f: complete weight dict (previously omitted half the
    # contributing axes from sensitivity analysis).
    base_weights = {
        'phylo': PHYLO_WEIGHT,
        'purifying': PURIFYING_WEIGHT,
        'positive': POSITIVE_WEIGHT,
        'synteny': SYNTENY_WEIGHT,
        'expression': EXPR_WEIGHT,
        'lse_depth': LSE_DEPTH_WEIGHT,
        'chemosensory_expr': CHEMOSENSORY_EXPR_WEIGHT,
        'gprotein_coexpr': GPROTEIN_COEXPR_WEIGHT,
        'ecl_divergence': ECL_DIVERGENCE_WEIGHT,
        'expansion': EXPANSION_WEIGHT,
        'og_confidence': OG_CONFIDENCE_WEIGHT,
        'tandem_cluster': TANDEM_CLUSTER_WEIGHT,
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

    # Bead -j6f: Latin Hypercube sweep over weight ORDERINGS (not just
    # ±50% magnitude perturbation). Generates n_lhs random weight schemes
    # uniformly from [0.5, 3.0] for each component, recomputes the rank,
    # and reports top-50 Jaccard-similarity stability across schemes.
    try:
        from scipy.stats import qmc as _qmc
        n_lhs = int(os.getenv('SENSITIVITY_LHS_N', 100))
        if n_lhs > 0 and len(df) > 0:
            sampler = _qmc.LatinHypercube(d=len(base_weights), seed=42)
            sample = sampler.random(n=n_lhs)
            lo, hi = 0.5, 3.0
            lhs_weights_matrix = lo + sample * (hi - lo)
            keys = list(base_weights.keys())
            top_n = min(50, len(df))
            top_sets = []
            for row in lhs_weights_matrix:
                w = dict(zip(keys, row.tolist()))
                # Recompute rank under this weight scheme
                tmp_scores = df.apply(
                    lambda r: _calculate_fair_rank_score_corrected(
                        scores={
                            'phylo': r.get('phylo_score_norm'),
                            'purifying': r.get('purifying_score_norm'),
                            'positive': r.get('positive_score_norm'),
                            'lse_depth': r.get('lse_depth_score_norm'),
                            'synteny': r.get('synteny_score_norm') if r.get('has_synteny_data') else None,
                            'expression': r.get('expression_score_norm') if r.get('has_expression_data') else None,
                            'chemosensory_expr': r.get('chemosensory_expr_score_norm') if r.get('has_chemosensory_expr_data') else None,
                            'gprotein_coexpr': r.get('gprotein_coexpr_score_norm') if r.get('has_gprotein_data') else None,
                            'ecl_divergence': r.get('ecl_divergence_score_norm') if r.get('has_ecl_data') else None,
                            'expansion': r.get('expansion_score_norm') if r.get('has_expansion_data') else None,
                            'og_confidence': r.get('og_confidence_score_norm') if r.get('has_og_confidence_data') else None,
                            'tandem_cluster': r.get('tandem_cluster_score_norm') if r.get('has_tandem_cluster_data') else None,
                        },
                        weights=w,
                        completeness_floor=0.4,
                    ),
                    axis=1,
                )
                tmp_top = set(df.assign(_s=tmp_scores)
                              .nlargest(top_n, '_s')['id'].tolist())
                top_sets.append(tmp_top)

            # Pairwise Jaccard
            jaccards = []
            for i in range(len(top_sets)):
                for j in range(i + 1, len(top_sets)):
                    a, b = top_sets[i], top_sets[j]
                    if a or b:
                        jaccards.append(len(a & b) / len(a | b))
            median_jaccard = float(np.median(jaccards)) if jaccards else float('nan')
            min_jaccard = float(np.min(jaccards)) if jaccards else float('nan')

            # Frequency table: how often each candidate appears in any top-N
            freq = {}
            for s in top_sets:
                for cid in s:
                    freq[cid] = freq.get(cid, 0) + 1
            freq_df = pd.DataFrame([(k, v / len(top_sets)) for k, v in freq.items()],
                                   columns=['id', 'lhs_top50_frequency'])

            output_dir = Path(output_file).parent
            with open(output_dir / 'weight_sensitivity.json', 'w') as f:
                json.dump({
                    'method': 'latin_hypercube',
                    'n_schemes': n_lhs,
                    'weight_range_low': lo,
                    'weight_range_high': hi,
                    'top_n': top_n,
                    'median_jaccard_topN': median_jaccard,
                    'min_jaccard_topN': min_jaccard,
                    'n_pairs_compared': len(jaccards),
                }, f, indent=2)
            freq_df.to_csv(output_dir / 'weight_sensitivity_topN_freq.csv', index=False)
            print(f"  LHS sensitivity: median top-{top_n} Jaccard={median_jaccard:.3f} "
                  f"(min={min_jaccard:.3f}, n_pairs={len(jaccards)})",
                  file=sys.stderr)
    except Exception as e:
        print(f"  LHS sensitivity skipped: {e}", file=sys.stderr)

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
    # Bead -urk: BUSTED/MEME diagnostic columns
    'busted_s_p', 'busted_s_significant', 'busted_mh_p', 'busted_mh_significant',
    'meme_n_episodic_sites', 'meme_fraction_episodic_sites',
    'synteny_score', 'has_synteny_data', 'expression_score', 'has_expression_data',
    'lse_depth_score', 'raw_tree_depth',
    # NEW: Phase 1-5 scores
    'chemosensory_expr_score', 'has_chemosensory_expr_data',
    'gprotein_coexpr_score', 'has_gprotein_data',
    'ecl_divergence_score', 'has_ecl_data',
    'expansion_score', 'has_expansion_data',
    'og_confidence_score', 'has_og_confidence_data',
    # Bead -ar8: tandem cluster columns
    'tandem_cluster_score', 'tandem_cluster_size', 'tandem_cluster_id',
    'has_tandem_cluster_data',
    # Bead -325: CDS provenance (native vs miniprot-recovered)
    'cds_source',
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
