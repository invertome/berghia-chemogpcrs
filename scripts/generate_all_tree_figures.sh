#!/bin/bash
# Generate all tree figure variants: PNG+PDF, light+dark themes
# Run from project root.
#
# Defaults to the March 2026 preliminary v2 phylogeny (historical paper
# artifact). Override via env vars for the May-2026+ production run:
#   TREE=results/phylogenies/protein/all_berghia_refs.treefile \
#   CONTREE=results/phylogenies/protein/all_berghia_refs.contree \
#   RANKING=results/ranking/ranked_candidates_sorted.csv \
#   CLASSIFICATION=results/ranking/berghia_gpcr_classification.csv \
#   OUTDIR=results/figures \
#   bash scripts/generate_all_tree_figures.sh

set -uo pipefail

TREE="${TREE:-preliminary/results/phylogenies/protein/v2/gpcrs.treefile}"
CONTREE="${CONTREE:-preliminary/results/phylogenies/protein/v2/gpcrs.contree}"
RANKING="${RANKING:-preliminary/results/ranking/ranked_candidates_sorted.csv}"
CLASSIFICATION="${CLASSIFICATION:-preliminary/results/ranking/berghia_gpcr_classification.csv}"
OUTDIR="${OUTDIR:-preliminary/results/phylogenies/protein/v2/figures}"
mkdir -p "$OUTDIR"

echo "Generating all tree figure variants..."

# Taxonomy tree
for theme in light dark; do
    for fmt in png pdf; do
        echo "  taxonomy_${theme}.${fmt}"
        python3 scripts/visualize_gpcr_tree.py \
            --tree "$TREE" \
            --ranking "$RANKING" \
            --output "$OUTDIR/gpcrs_taxonomy_${theme}.${fmt}" \
            --theme "$theme" \
            --dpi 300
    done
done

# Functional classification tree
for theme in light dark; do
    for fmt in png pdf; do
        echo "  classification_${theme}.${fmt}"
        python3 scripts/visualize_gpcr_classification.py \
            --tree "$TREE" \
            --ranking "$RANKING" \
            --classification "$CLASSIFICATION" \
            --output "$OUTDIR/gpcrs_classification_${theme}.${fmt}" \
            --theme "$theme" \
            --dpi 300
    done
done

# Bootstrap histogram
for theme in light dark; do
    for fmt in png pdf; do
        echo "  bootstrap_${theme}.${fmt}"
        python3 scripts/visualize_bootstrap_histogram.py \
            --contree "$CONTREE" \
            --output "$OUTDIR/gpcrs_bootstrap_${theme}.${fmt}" \
            --theme "$theme" \
            --dpi 300
    done
done

echo "Done. Figures in: $OUTDIR/"
ls -lh "$OUTDIR/"
