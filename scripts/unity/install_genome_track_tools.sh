#!/bin/bash
# Install the genome-track reconciliation tools (bead: genome-track stage 02c).
set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh" 2>/dev/null || \
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
ENV_NAME="${1:-berghia-gpcr}"
mamba install -y -n "$ENV_NAME" -c bioconda -c conda-forge gmap gffcompare bedtools miniprot blast
conda activate "$ENV_NAME"
for t in gmap gmap_build gffcompare bedtools miniprot blastp makeblastdb minimap2 diamond; do
    printf "%-12s %s\n" "$t" "$(command -v "$t" || echo MISSING)"
done
