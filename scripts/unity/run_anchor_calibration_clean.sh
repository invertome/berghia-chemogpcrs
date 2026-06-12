#!/bin/bash
# run_anchor_calibration_clean.sh — DASH-free C3 calibration re-run, all 4 classes.
#
# The first calibration was invalidated: MAFFT-DASH injected 'DASH|<pdb>_<chain>||...'
# structural-homolog rows into every per-class alignment/tree (fixed by adding
# --originalseqonly to the --dash aligners, commit 1f3b0af), and the eval crashed
# on 3-part aBayes/SH-aLRT/UFBoot support (fixed, commit 7685660). The POOLS are
# clean (DASH is injected at the alignment step, not pool build), so this REUSES
# the existing refs_class_*.fa and rebuilds DASH-free trees + verdicts. Class A no
# longer carries ~239 DASH rows, so 128G clears the earlier 64G OOM.
#
# Submit:  sbatch --array=1-4 scripts/unity/run_anchor_calibration_clean.sh
#          (task 1=A, 2=B, 3=C, 4=F)
#
#SBATCH --job-name=anchor_cal_clean
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --output=logs/anchor_cal_clean-%A_%a.out
#SBATCH --error=logs/anchor_cal_clean-%A_%a.err
#SBATCH --account=pi_pkatz_umass_edu

set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
set -u

REPO_ROOT="${REPO_ROOT:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
cd "${REPO_ROOT}"
# shellcheck source=/dev/null
source "${REPO_ROOT}/config.sh"
# shellcheck source=/dev/null
source "${REPO_ROOT}/functions.sh"

# --- Array task -> GPCR class (1=A, 2=B, 3=C, 4=F) --------------------------
CLASSES=(_ A B C F)
CLS="${CLASSES[${SLURM_ARRAY_TASK_ID:?submit with --array=1-4}]}"

IQTREE_BOOT_FLAGS="-B ${IQTREE_BOOTSTRAP} -alrt 1000"
if [ "${IQTREE_TBE:-1}" = "1" ]; then
    IQTREE_BOOT_FLAGS+=" --tbe"
fi
CPUS="${CPUS:-${SLURM_CPUS_PER_TASK:-16}}"

CAL="${RESULTS_DIR}/p5_phase1a_validation/anchor_calibration"
TREES_DIR="${CAL}/trees_clean"
VERD_DIR="${CAL}/verdicts_clean"
export LOGS_DIR="${CAL}/logs_clean"     # isolate filter-stack tool stderr
mkdir -p "${TREES_DIR}" "${VERD_DIR}" "${LOGS_DIR}"

with_fa="${CAL}/pools/refs_class_${CLS}.fa"
[ -s "$with_fa" ] || { log --level=ERROR "anchor_cal_clean: class ${CLS} pool missing: $with_fa"; exit 1; }

build_tree() {                       # $1=input_fa  $2=prefix(no ext)
    local in_fa="$1" prefix="$2" tag
    tag="$(basename "$prefix")"
    local wd; wd="$(dirname "$prefix")"
    run_alignment_filter_stack "$in_fa" "${prefix}_aligned.fa" \
        "${wd}/_fs_${tag}" "$tag" "$CPUS"
    # Fail-fast: the --originalseqonly fix must keep DASH homologs out of output.
    if grep -q '^>DASH|' "${prefix}_aligned.fa" 2>/dev/null; then
        log --level=ERROR "anchor_cal_clean: DASH| rows still in ${prefix}_aligned.fa — --originalseqonly not applied; aborting"
        exit 3
    fi
    # shellcheck disable=SC2086
    ${CLIPKIT} "${prefix}_aligned.fa" -m kpic-smart-gap -o "${prefix}_trimmed.fa"
    python3 "${SCRIPTS_DIR}/drop_near_all_gap_rows.py" \
        --input "${prefix}_trimmed.fa" --output "${prefix}_trimmed.fa" || true
    # shellcheck disable=SC2086
    ${FASTTREE} -seed "${FASTTREE_SEED}" -lg -gamma "${prefix}_trimmed.fa" \
        > "${prefix}_fasttree.tre" 2>/dev/null || true
    local seed_arg=""
    [ -s "${prefix}_fasttree.tre" ] && seed_arg="-t ${prefix}_fasttree.tre"
    # shellcheck disable=SC2086
    ${IQTREE} -s "${prefix}_trimmed.fa" -st AA \
        -m "${IQTREE_MODEL_FIND}" -mset "${IQTREE_MODEL_SET}" \
        ${IQTREE_BOOT_FLAGS} -seed "${IQTREE_SEED}" -T "${CPUS}" \
        ${seed_arg} --prefix "${prefix}"
}

# WITHOUT out-group anchors = strip tier-2/3 anchor records (identical in-group).
without_fa="${TREES_DIR}/refs_class_${CLS}_noOG.fa"
awk '/^>/{og=($0 ~ /^>ANCHOR_[A-F]_[23]_/); if(!og)print; next}{if(!og)print}' \
    "$with_fa" > "$without_fa"

log "anchor_cal_clean: class ${CLS} — building DASH-free WITH / WITHOUT trees"
build_tree "$with_fa"    "${TREES_DIR}/with_class_${CLS}"
build_tree "$without_fa" "${TREES_DIR}/without_class_${CLS}"

with_tree="${TREES_DIR}/with_class_${CLS}.treefile"
without_tree="${TREES_DIR}/without_class_${CLS}.treefile"
if [ -s "$with_tree" ] && [ -s "$without_tree" ]; then
    python3 scripts/evaluate_anchor_divergence.py \
        --class "${CLS}" \
        --tree-without "$without_tree" \
        --tree-with "$with_tree" \
        --pool-members "${CAL}/pools/pool_members_class_${CLS}.tsv" \
        --out "${VERD_DIR}/verdict_class_${CLS}.json"
    log "anchor_cal_clean: class ${CLS} verdict -> ${VERD_DIR}/verdict_class_${CLS}.json"
else
    log --level=WARN "anchor_cal_clean: class ${CLS} tree(s) missing, no verdict"
fi
