#!/bin/bash
# sweep_phase1f_incomplete.sh — Phase-1f completeness reconciliation.
#
# The main BRAKER array processes samples_full[task_id] for a fixed task range.
# That leaves gaps when: (a) the buildable set changes mid-run (genomes added
# after submit shift samples_full, pushing species past the array window);
# (b) tasks fail (e.g. a stale snakemake lock in a reused RUN_DIR); (c) the
# array range simply undershoots the current row count. This sweep closes all
# of those the same way: rebuild samples_full, then submit ONE BRAKER job per
# species whose proteome is not yet cached.
#
# Each job runs in MANUAL mode (SAMPLES_LIMIT_REGEX='^<sample>$'), so its RUN_DIR
# is keyed by sample name (manual_<sample>), NOT task id — it never collides with
# the array's task_<n> RUN_DIRs or their stale locks.
#
# Idempotent: a species with a cached proteome is skipped, so re-running the
# sweep only ever submits the still-missing remainder.
#
# Usage (run when the main array has finished, or as a dependency job):
#   sbatch --dependency=afterany:<ARRAY_JOBID> scripts/unity/sweep_phase1f_incomplete.sh
#
#SBATCH --job-name=phase1f_sweep
#SBATCH --partition=cpu
#SBATCH --time=30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --output=logs/phase1f_sweep-%j.out
#SBATCH --error=logs/phase1f_sweep-%j.err

set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr

REPO_ROOT="${REPO_ROOT:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
cd "${REPO_ROOT}"
mkdir -p logs

MANIFEST="references/species_tree/genome_inventory.tsv"
GENOME_CACHE="species_tree_data/braker4_genomes"
PROTEIN_DB="species_tree_data/orthodb/odb12_metazoa_plus_mollusca_aa.linclust0.9_rep_seq.fasta"
PROT_CACHE="references/species_tree/cache/proteomes_braker4"
WRAPPER="sbatch_run_species_tree_phase1f_local.sh"
THROTTLE="${SWEEP_THROTTLE:-40}"        # cap concurrent sweep jobs

MAX_SUBMIT="${SWEEP_MAX_SUBMIT:-80}"    # flood guard: abort if more than this look missing

SF="$(mktemp)"
python3 scripts/build_braker4_samples_csv.py \
    --manifest "$MANIFEST" --genome-cache "$GENOME_CACHE" \
    --protein-db "$PROTEIN_DB" --out "$SF"

# First pass: collect the still-missing species (proteome not cached).
total=0; done=0
missing=()
while IFS=, read -r sample _; do
    [ "$sample" = "sample_name" ] && continue
    [ -z "$sample" ] && continue
    total=$((total + 1))
    if [ -s "${PROT_CACHE}/${sample}.aa.fna" ]; then
        done=$((done + 1))
    else
        missing+=("$sample")
    fi
done < "$SF"
rm -f "$SF"

echo "phase1f sweep: ${total} buildable, ${done} already cached, ${#missing[@]} missing"

# Flood guard: a huge 'missing' count means a path/config bug (e.g. wrong
# proteome cache), not a real reconciliation — do not carpet-bomb the queue.
if [ "${#missing[@]}" -gt "$MAX_SUBMIT" ]; then
    echo "ERROR: ${#missing[@]} missing exceeds SWEEP_MAX_SUBMIT=${MAX_SUBMIT}; aborting." >&2
    echo "       Check PROT_CACHE=${PROT_CACHE} and the manifest before re-running." >&2
    exit 1
fi

# Each missing species runs in MANUAL mode (no --array) -> RUN_DIR manual_<sample>,
# isolated from the array's task_<n> dirs.
submitted=0
for sample in "${missing[@]}"; do
    sbatch --job-name="p1f_sweep_${sample}" \
        --export=ALL,SAMPLES_LIMIT_REGEX="^${sample}$" \
        "$WRAPPER" >/dev/null
    submitted=$((submitted + 1))
done

echo "=== phase1f sweep: submitted ${submitted} reconciliation job(s) ==="
