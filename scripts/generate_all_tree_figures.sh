#!/bin/bash
# Generate all tree figure variants: PNG+PDF, light+dark themes
# Run from project root

set -uo pipefail

TREE="preliminary/results/phylogenies/protein/v2/gpcrs.treefile"
RANKING="preliminary/results/ranking/ranked_candidates_sorted.csv"
CLASSIFICATION="preliminary/results/ranking/berghia_gpcr_classification.csv"
OUTDIR="preliminary/results/phylogenies/protein/v2/figures"
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
            --contree "preliminary/results/phylogenies/protein/v2/gpcrs.contree" \
            --output "$OUTDIR/gpcrs_bootstrap_${theme}.${fmt}" \
            --theme "$theme" \
            --dpi 300
    done
done

echo "Done. Figures in: $OUTDIR/"
ls -lh "$OUTDIR/"
