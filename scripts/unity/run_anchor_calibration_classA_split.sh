#!/bin/bash
# run_anchor_calibration_classA_split.sh — Class A C3 calibration, with/without
# split into two independent jobs so each tree gets its own 14-day wall.
#
# WHY a separate script for class A:
#   The 4-class clean array (run_anchor_calibration_clean.sh) builds BOTH the
#   with-anchor and without-anchor trees inside ONE 48h task. For class A (949
#   taxa) that is infeasible: ModelFinder alone (-m MFP, 308 protein models on
#   949 taxa) ran >48h and never finished — the 2026-06-14 task 60770176_1 hit
#   the 48h wall still grinding through the +R series, ~100/308 models in, and
#   never reached the ML search or bootstrap. B/C/F are 59-277 taxa and finish
#   easily, so they stay on the 48h array. Class A is split here: task 1 = with,
#   task 2 = without, each --time=14d --mem=384G.
#
#   METHOD IS IDENTICAL to B/C/F: same -m MFP -mset (full 308-model search),
#   same UFBoot 1000 + SH-aLRT 1000 + TBE, same FastTree seed, same filter
#   stack. The only change is job partitioning + resources — no model, sampling,
#   or threshold change.
#
#   RESUME: the with-tree's DASH-free alignment (with_class_A_trimmed.fa) and
#   ModelFinder checkpoint (with_class_A.model.gz) already exist from the failed
#   run. build_tree reuses the alignment as-is and does NOT pass -redo, so
#   IQ-TREE resumes ModelFinder from where it stopped (~model 100/308). The
#   without-tree builds the full chain fresh and writes its own checkpoint.
#
# Submit:
#   AID=$(sbatch --parsable --array=1-2 scripts/unity/run_anchor_calibration_classA_split.sh)
#   sbatch --dependency=afterok:$AID scripts/unity/run_anchor_calibration_classA_agg.sh
#   (task 1 = with-anchor tree, task 2 = without-anchor tree)
#
#SBATCH --job-name=anchor_cal_A_split
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --time=14-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=384G
#SBATCH --output=logs/anchor_cal_A_split-%A_%a.out
#SBATCH --error=logs/anchor_cal_A_split-%A_%a.err
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

CLS="A"   # this script is class-A only

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

# --- Array task -> which tree (1=with, 2=without) --------------------------
case "${SLURM_ARRAY_TASK_ID:?submit with --array=1-2}" in
    1) WHICH=with ;;
    2) WHICH=without ;;
    *) log --level=ERROR "anchor_cal_A_split: --array id must be 1 (with) or 2 (without)"; exit 1 ;;
esac

build_tree() {                       # $1=input_fa  $2=prefix(no ext)
    local in_fa="$1" prefix="$2" tag
    tag="$(basename "$prefix")"
    local wd; wd="$(dirname "$prefix")"
    if [ -s "${prefix}_trimmed.fa" ]; then
        log "anchor_cal_A_split: ${tag} — reusing existing DASH-free trimmed alignment (resume); skipping filter stack"
        # Robustness check still applies on the reused alignment.
        if grep -q '^>DASH|' "${prefix}_trimmed.fa" 2>/dev/null; then
            log --level=ERROR "anchor_cal_A_split: DASH| rows in reused ${prefix}_trimmed.fa — aborting"
            exit 3
        fi
    else
        run_alignment_filter_stack "$in_fa" "${prefix}_aligned.fa" \
            "${wd}/_fs_${tag}" "$tag" "$CPUS"
        # Fail-fast: the --originalseqonly fix must keep DASH homologs out of output.
        if grep -q '^>DASH|' "${prefix}_aligned.fa" 2>/dev/null; then
            log --level=ERROR "anchor_cal_A_split: DASH| rows still in ${prefix}_aligned.fa — --originalseqonly not applied; aborting"
            exit 3
        fi
        # shellcheck disable=SC2086
        ${CLIPKIT} "${prefix}_aligned.fa" -m kpic-smart-gap -o "${prefix}_trimmed.fa"
        python3 "${SCRIPTS_DIR}/drop_near_all_gap_rows.py" \
            --input "${prefix}_trimmed.fa" --output "${prefix}_trimmed.fa" || true
    fi
    if [ ! -s "${prefix}_fasttree.tre" ]; then
        # shellcheck disable=SC2086
        ${FASTTREE} -seed "${FASTTREE_SEED}" -lg -gamma "${prefix}_trimmed.fa" \
            > "${prefix}_fasttree.tre" 2>/dev/null || true
    fi
    local seed_arg=""
    [ -s "${prefix}_fasttree.tre" ] && seed_arg="-t ${prefix}_fasttree.tre"
    # NO -redo: IQ-TREE resumes from ${prefix}.model.gz / .ckp.gz if present.
    # shellcheck disable=SC2086
    ${IQTREE} -s "${prefix}_trimmed.fa" -st AA \
        -m "${IQTREE_MODEL_FIND}" -mset "${IQTREE_MODEL_SET}" \
        ${IQTREE_BOOT_FLAGS} -seed "${IQTREE_SEED}" -T "${CPUS}" \
        ${seed_arg} --prefix "${prefix}"
}

if [ "${WHICH}" = "with" ]; then
    with_fa="${CAL}/pools/refs_class_${CLS}.fa"
    [ -s "$with_fa" ] || { log --level=ERROR "anchor_cal_A_split: class ${CLS} pool missing: $with_fa"; exit 1; }
    log "anchor_cal_A_split: class ${CLS} — building WITH-anchor tree (resume-aware)"
    build_tree "$with_fa" "${TREES_DIR}/with_class_${CLS}"
else
    # WITHOUT out-group anchors = strip tier-2/3 anchor records (identical in-group).
    with_fa="${CAL}/pools/refs_class_${CLS}.fa"
    without_fa="${TREES_DIR}/refs_class_${CLS}_noOG.fa"
    if [ ! -s "$without_fa" ]; then
        [ -s "$with_fa" ] || { log --level=ERROR "anchor_cal_A_split: class ${CLS} pool missing: $with_fa"; exit 1; }
        awk '/^>/{og=($0 ~ /^>ANCHOR_[A-F]_[23]_/); if(!og)print; next}{if(!og)print}' \
            "$with_fa" > "$without_fa"
    fi
    log "anchor_cal_A_split: class ${CLS} — building WITHOUT-anchor tree (resume-aware)"
    build_tree "$without_fa" "${TREES_DIR}/without_class_${CLS}"
fi

log "anchor_cal_A_split: class ${CLS} ${WHICH}-tree task done; verdict computed by classA_agg once both trees exist"
