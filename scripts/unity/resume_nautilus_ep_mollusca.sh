#!/bin/bash
# resume_nautilus_ep_mollusca.sh — validate the GeneMark-EP -> ES fallback patch
# by RESUMING the existing nautilus_ep_mollusca BRAKER4 run.
#
# The pilot (job 60716113) completed masking/compleasm/GeneMark-ES/ProtHint, then
# GeneMark-EP crashed (parse_ET div-by-zero) -> no annotation. With the EP rule
# now patched (scripts/unity/patch_braker4_run_genemark_ep.py) to fall back to the
# GeneMark-ES genes, re-running snakemake in the SAME workdir re-runs only the
# (now non-fatal) EP rule -> ES fallback -> AUGUSTUS/TSEBRA -> braker.aa. The slow
# upstream steps are reused (~hours, not the full ~1.5 days).
#
# Submit:  sbatch scripts/unity/resume_nautilus_ep_mollusca.sh
#
#SBATCH --job-name=naut_resume
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --output=logs/naut_resume-%j.out
#SBATCH --error=logs/naut_resume-%j.err
#SBATCH --account=pi_pkatz_umass_edu

set -eo pipefail
mkdir -p logs
source ~/.miniconda3/etc/profile.d/conda.sh
set +u; conda activate braker4_runner; set -u 2>/dev/null || true

REPO_ROOT="/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05"
cd "$REPO_ROOT"
SAMPLE="34576_Nautilus_macromphalus"
RUN_DIR="${REPO_ROOT}/species_tree_data/braker4_run/array_runs/nautilus_ep_mollusca"
[ -d "$RUN_DIR" ] || { echo "ERROR: no existing run dir $RUN_DIR (nothing to resume)"; exit 1; }

# Re-assert the patch is live before resuming (idempotent).
python3 scripts/unity/patch_braker4_run_genemark_ep.py >/dev/null 2>&1 || true

SIF_CACHE="${REPO_ROOT}/species_tree_data/braker4_run/.singularity_cache"
CONFIG_INI="${REPO_ROOT}/scripts/unity/braker4_config_ep_metazoa.ini"
SNAKEFILE="${REPO_ROOT}/external/braker4/Snakefile"
POSTPROCESS="${REPO_ROOT}/scripts/postprocess_braker4_outputs.py"
export AUGUSTUS_CONFIG_PATH="${RUN_DIR}/augustus_config"
export BRAKER4_CONFIG="$CONFIG_INI"
export APPTAINER_BIND="${REPO_ROOT}:${REPO_ROOT}"

# Clear the failed GeneMark-EP output so snakemake re-runs (now-patched) EP only.
rm -f "${RUN_DIR}/output/${SAMPLE}/GeneMark-EP/genemark.gtf"

cd "$RUN_DIR"
set +e
snakemake \
    --use-singularity \
    --singularity-prefix "$SIF_CACHE" \
    --singularity-args "-B ${REPO_ROOT}:${REPO_ROOT}" \
    --executor local \
    --cores "${SLURM_CPUS_PER_TASK:-24}" \
    --keep-going --rerun-incomplete --rerun-triggers mtime \
    --snakefile "$SNAKEFILE"
RC=$?
set -e
echo "[resume] snakemake rc=$RC"
cd "$REPO_ROOT"

python "$POSTPROCESS" \
    --braker4-output "${RUN_DIR}/output" \
    --samples-csv "${RUN_DIR}/samples.csv" \
    --cache-dir "${RUN_DIR}/cache" \
    --report "${RUN_DIR}/postprocess_report.tsv" || echo "[resume] postprocess non-zero"

echo "=== EP fallback message (should mention GeneMark-ES) ==="
grep -i "falling back to GeneMark-ES" "${RUN_DIR}/logs/${SAMPLE}/genemark_ep/genemark_ep.log" 2>/dev/null | head -1 || echo "(no fallback message — check EP log)"
echo "=== result ==="
cat "${RUN_DIR}/postprocess_report.tsv" 2>/dev/null
ls -la "${RUN_DIR}/cache/"*.aa.fna 2>/dev/null && echo "proteins: $(grep -c '^>' "${RUN_DIR}/cache/${SAMPLE}.aa.fna" 2>/dev/null)" || echo "NO proteome produced"
