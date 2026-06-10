#!/bin/bash
# run_anchor_calibration_bcf.sh — C3 calibration for classes B/C/F (parallel side-array).
#
# Companion to run_anchor_calibration.sh (bead berghia-chemogpcrs-521.4). The
# monolithic driver (job 60689713) builds all four classes SEQUENTIALLY and will
# only reach class A inside its 48h wall (class-A tree pair ~39h alone). This
# array builds the WITH/WITHOUT trees + verdicts for the three remaining classes
# IN PARALLEL, REUSING the per-class pools that the monolith already built
# (results/.../anchor_calibration/pools/, written 01:22).
#
# Collision-free with the running monolith: all outputs go to DISJOINT dirs
# (trees_bcf/, verdicts_bcf/, logs_bcf/). The pools/ dir is read-only here.
# Everything else (CD-HIT identity, filter stack, IQ-TREE flags, evaluate logic)
# is identical to the monolithic driver — this only parallelizes the remaining
# classes, it changes no calibration parameter.
#
# Submit:  sbatch --array=1-3 scripts/unity/run_anchor_calibration_bcf.sh
#          (task 1 = class B, 2 = class C, 3 = class F)
#
#SBATCH --job-name=anchor_cal_bcf
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=logs/anchor_cal_bcf-%A_%a.out
#SBATCH --error=logs/anchor_cal_bcf-%A_%a.err
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

# --- Array task -> GPCR class (1=B, 2=C, 3=F) -------------------------------
CLASSES=(_ B C F)
CLS="${CLASSES[${SLURM_ARRAY_TASK_ID:?submit with --array=1-3}]}"

# Canonical UFBoot + SH-aLRT (+TBE) — identical to run_anchor_calibration.sh.
IQTREE_BOOT_FLAGS="-B ${IQTREE_BOOTSTRAP} -alrt 1000"
if [ "${IQTREE_TBE:-1}" = "1" ]; then
    IQTREE_BOOT_FLAGS+=" --tbe"
fi
CPUS="${CPUS:-${SLURM_CPUS_PER_TASK:-16}}"

CAL="${RESULTS_DIR}/p5_phase1a_validation/anchor_calibration"
TREES_DIR="${CAL}/trees_bcf"
VERD_DIR="${CAL}/verdicts_bcf"
# Isolate filter-stack tool stderr: functions.sh keys *_${tag}.err on $LOGS_DIR,
# so override it to a side-job dir to stay fully disjoint from the monolith.
export LOGS_DIR="${CAL}/logs_bcf"
mkdir -p "${TREES_DIR}" "${VERD_DIR}" "${LOGS_DIR}"

with_fa="${CAL}/pools/refs_class_${CLS}.fa"
[ -s "$with_fa" ] || { log --level=ERROR "anchor_cal_bcf: class ${CLS} pool missing: $with_fa"; exit 1; }

# build_tree — copied verbatim from run_anchor_calibration.sh (md5-matched driver)
build_tree() {                       # $1=input_fa  $2=prefix(no ext)
    local in_fa="$1" prefix="$2" tag
    tag="$(basename "$prefix")"
    local wd; wd="$(dirname "$prefix")"
    run_alignment_filter_stack "$in_fa" "${prefix}_aligned.fa" \
        "${wd}/_fs_${tag}" "$tag" "$CPUS"
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

log "anchor_cal_bcf: class ${CLS} — building WITH / WITHOUT trees (parallel side-array)"
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
    log "anchor_cal_bcf: class ${CLS} verdict -> ${VERD_DIR}/verdict_class_${CLS}.json"
else
    log --level=WARN "anchor_cal_bcf: class ${CLS} tree(s) missing, no verdict"
fi
