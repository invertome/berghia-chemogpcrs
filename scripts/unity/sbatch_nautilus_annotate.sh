#!/bin/bash
#SBATCH --job-name=nautilus_annot
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --time=5-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --output=logs/nautilus_annotate-%j.out
#SBATCH --error=logs/nautilus_annotate-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jperezmoreno@umass.edu
#SBATCH --account=pi_pkatz_umass_edu
#
# One-off recovery of Nautilus macromphalus (34576), which failed default
# BRAKER-EP (GeneMark-EP div-by-zero: ProtHint gave only ~85 hits from OrthoDB
# on this deep-branching lineage). MODE env var picks the approach:
#   MODE=ep_enriched -> protein_fasta = OrthoDB metazoa + N. pompilius (same
#                       genus, 30,895 proteins) -> abundant ProtHint hits -> EP
#   MODE=ep_mollusca -> protein_fasta = OrthoDB metazoa + MolluscaGenes mollusca_aa
#                       (200+ mollusc spp, ~57M seqs): the STANDARDIZED evidence
#   MODE=es          -> no protein evidence -> GeneMark-ES ab initio
# Runs in an ISOLATED RUN_DIR + real (non-shared) output so it cannot collide
# with the main array or the other mode. Submit: MODE=es sbatch this; MODE=ep_enriched sbatch this.

set -eo pipefail
mkdir -p logs
source ~/.miniconda3/etc/profile.d/conda.sh
set +u; conda activate braker4_runner; set -u 2>/dev/null || true

REPO_ROOT="/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05"
cd "$REPO_ROOT"
MODE="${MODE:?set MODE=ep_mollusca, ep_enriched, or es}"
SAMPLE="34576_Nautilus_macromphalus"

GENOME="${REPO_ROOT}/species_tree_data/braker4_genomes/${SAMPLE}.fasta"
ODB="${REPO_ROOT}/species_tree_data/orthodb/odb12_metazoa.fa"
POMPILIUS="${REPO_ROOT}/references/species_tree/cache/proteomes_braker4/34573_Nautilus_pompilius.aa.fna"
SIF_CACHE="${REPO_ROOT}/species_tree_data/braker4_run/.singularity_cache"
CONFIG_INI="${REPO_ROOT}/scripts/unity/braker4_config_ep_metazoa.ini"
SNAKEFILE="${REPO_ROOT}/external/braker4/Snakefile"
POSTPROCESS="${REPO_ROOT}/scripts/postprocess_braker4_outputs.py"

RUN_DIR="${REPO_ROOT}/species_tree_data/braker4_run/array_runs/nautilus_${MODE}"
rm -rf "$RUN_DIR"; mkdir -p "$RUN_DIR/output" "$RUN_DIR/augustus_config" "$RUN_DIR/cache"
export AUGUSTUS_CONFIG_PATH="${RUN_DIR}/augustus_config"
export BRAKER4_CONFIG="$CONFIG_INI"
export APPTAINER_BIND="${REPO_ROOT}:${REPO_ROOT}"

PROT=""
if [ "$MODE" = "ep_enriched" ]; then
    PROT="${RUN_DIR}/nautilus_enriched_proteins.fa"
    cat "$ODB" "$POMPILIUS" > "$PROT"
    echo "[nautilus] enriched protein set: $(grep -c '^>' "$PROT") seqs (OrthoDB + pompilius)"
elif [ "$MODE" = "ep_mollusca" ]; then
    # Standardized evidence: OrthoDB metazoa + MolluscaGenes mollusca_aa (200+ spp).
    # Use the prebuilt combined DB directly (no per-run copy of ~25 GB).
    PROT="${REPO_ROOT}/species_tree_data/orthodb/odb12_metazoa_plus_mollusca_aa.fa"
    [ -s "$PROT" ] || { echo "[nautilus] ERROR: combined DB missing: $PROT (build via scripts/unity/build_braker4_prothint_db.sh)"; exit 1; }
    echo "[nautilus] mollusca-enriched protein set (standardized): $PROT"
fi

# samples.csv (14 cols), built by joining an explicit field array (no miscount)
HDR="sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf"
echo "$HDR" > "${RUN_DIR}/samples.csv"
fields=("$SAMPLE" "$GENOME" "" "$PROT" "" "" "" "" "" "" "" "" "metazoa_odb12" "")
(IFS=,; echo "${fields[*]}") >> "${RUN_DIR}/samples.csv"
echo "[nautilus] samples.csv:"; cat "${RUN_DIR}/samples.csv"

cd "$RUN_DIR"
set +e
snakemake \
    --use-singularity \
    --singularity-prefix "$SIF_CACHE" \
    --singularity-args "-B ${REPO_ROOT}:${REPO_ROOT}" \
    --executor local \
    --cores "${SLURM_CPUS_PER_TASK:-24}" \
    --keep-going --rerun-incomplete \
    --snakefile "$SNAKEFILE"
RC=$?
set -e
echo "[nautilus] snakemake rc=$RC (mode=$MODE)"
cd "$REPO_ROOT"

python "$POSTPROCESS" \
    --braker4-output "${RUN_DIR}/output" \
    --samples-csv "${RUN_DIR}/samples.csv" \
    --cache-dir "${RUN_DIR}/cache" \
    --report "${RUN_DIR}/postprocess_report.tsv" || echo "[nautilus] postprocess non-zero (mode=$MODE)"

echo "=== ${SAMPLE} mode=${MODE} result ==="
ls -la "${RUN_DIR}/cache/"*.aa.fna 2>/dev/null || echo "  NO proteome produced"
grep -c '^>' "${RUN_DIR}/cache/${SAMPLE}.aa.fna" 2>/dev/null | sed 's/^/  protein count: /' || true
