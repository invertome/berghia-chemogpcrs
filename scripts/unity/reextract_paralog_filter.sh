#!/bin/bash
#SBATCH --job-name=reextract_paralog
#SBATCH --partition=cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=results/logs/reextract_paralog_%j.out
#SBATCH --error=results/logs/reextract_paralog_%j.err
#
# Bead -8st follow-up: apply the tightened paralog-discrimination thresholds
# (commit 9e2ab2d) to the existing miniprot recovery WITHOUT re-running
# miniprot. Re-uses cached GFFs in results/reference_sequences/cds_recovered/
# and cached genomes in genomes/ref_assemblies/. Runs extract_cds_from_gff
# under:
#   CDS_PROTEIN_IDENTITY_MIN=0.95   (was 0.50; paralog-discrimination)
#   CDS_AMBIGUITY_MARGIN=0.05       (best/2nd-best identity gate)
#   CDS_OVERLAP_FRACTION=0.5        (genomic-locus dedup)
#
# After this finishes, stage 01's merge step picks up the cleaned
# cds_recovered/ output and concatenates into all_references_cds.fna.

set -eo pipefail
cd /scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05

source ~/.miniconda3/etc/profile.d/conda.sh
conda activate berghia-gpcr

# --cds-file /dev/null forces get_missing_proteins to treat every protein
# as missing (cds_ids set is empty), so we re-process all of them. The
# script's --reextract-only flag then short-circuits the miniprot run and
# just re-runs extract_cds_from_gff against the cached GFFs.
mkdir -p results/logs

python3 scripts/recover_cds_from_assemblies.py \
    --og-dir references/nath_et_al \
    --cds-file /dev/null \
    --output-dir results/reference_sequences/cds_recovered \
    --genome-dir genomes/ref_assemblies \
    --threads ${SLURM_CPUS_PER_TASK:-4} \
    --keep-genomes \
    --reextract-only \
    --manifest results/reference_sequences/cds_provenance.csv

echo "[reextract] DONE — cds_recovered/ now contains paralog-filtered CDS"
echo "[reextract] sequence count by species:"
for f in results/reference_sequences/cds_recovered/*_cds.fna; do
    [ -f "$f" ] && echo "  $(basename "$f"): $(grep -c '^>' "$f")"
done
