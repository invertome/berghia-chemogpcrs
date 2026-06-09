#!/bin/bash
# run_anchor_calibration.sh — C3 anchor-divergence calibration pilot (Unity).
#
# Bead berghia-chemogpcrs-521.4. For each GPCR class, build a per-class pool WITH
# all anchors, derive the WITHOUT-out-group-anchors counterpart by removing the
# tier-2/3 anchor tips (so the in-group is IDENTICAL between the two — a clean RF
# comparison), build a tree for each with the canonical stage-04 stack
# (PREQUAL -> MAFFT ensemble -> CLOAK -> TAPER -> ClipKit -> FastTree seed ->
# IQ-TREE UFBoot+SH-aLRT+TBE), then run evaluate_anchor_divergence.py on the pair
# and aggregate the per-class verdicts into a REPORT (report-only; nothing is
# applied to config — review the trees, then set ANCHOR_OUTGROUP_CLASSES).
#
# Reuses the existing Phase-1a scan + classify outputs (results/p5_phase1a_validation).
#
# Pilot scaling (spec §8 — a fast LBA read, not the final P7 run):
#   PILOT_BERGHIA_MAX  cap on focal-species tips (seeded subsample; 0 = keep all)
#   PILOT_FLOOR/PILOT_CAP  reference per-species floor / per-taxon cap (thin the
#                          reference backdrop). Anchors are always kept in full.
# Keeping all focal-species paralogs is the safer (signal-preserving) default at
# the cost of a heavier Class-A run; cap PILOT_BERGHIA_MAX for a faster pilot.
#
# CD-HIT identity is forced to 0.9 here (the C4 value) regardless of the repo's
# config.sh, so the calibration matches the intended production setting.
#
# Submit:  sbatch scripts/unity/run_anchor_calibration.sh
#
#SBATCH --job-name=anchor_calibration
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/anchor_calibration-%j.out
#SBATCH --error=logs/anchor_calibration-%j.err
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

# IQ-TREE bootstrap flags (defined in 04_phylogenetic_analysis.sh, not config.sh —
# replicate the canonical UFBoot + SH-aLRT (+TBE) combination here).
IQTREE_BOOT_FLAGS="-B ${IQTREE_BOOTSTRAP} -alrt 1000"
if [ "${IQTREE_TBE:-1}" = "1" ]; then
    IQTREE_BOOT_FLAGS+=" --tbe"
fi

# --- Calibration knobs (override the repo defaults for this run) ------------
export REF_CLUSTER_IDENTITY=0.9                     # C4 value, independent of config.sh
# Moderate pilot defaults → Class-A tree ~1000 tips (anchors ~177 + focal ~400 +
# refs ~450): a sensitive-enough LBA read per spec §8 without the full-run cost.
# The final P7 build uses the full pool (all focal-species paralogs, full ref
# budget); the full pre-P7 confirmation re-checks any clean verdict on all tips.
export REF_FLOOR_PER_SPECIES="${PILOT_FLOOR:-1}"
export REF_CAP_PER_TAXON="${PILOT_CAP:-5}"
PILOT_BERGHIA_MAX="${PILOT_BERGHIA_MAX:-400}"       # 0 = keep all focal-species tips
CPUS="${CPUS:-${SLURM_CPUS_PER_TASK:-16}}"

# --- Reuse Phase-1a P5 scan + classify outputs ------------------------------
WORK_DIR="${RESULTS_DIR}/p5_phase1a_validation"
SCAN_DIR="${WORK_DIR}/scan"
CLASSIFY_DIR="${WORK_DIR}/classify"
CLASS_PHASE1A="${CLASSIFY_DIR}/class_phase1a.tsv"
CLASS_BERGHIA="${CLASSIFY_DIR}/class_berghia.tsv"
BERGHIA_FA="${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa"

CAL="${WORK_DIR}/anchor_calibration"
mkdir -p "${CAL}/pools" "${CAL}/trees" "${CAL}/verdicts" logs references/anchors

for f in "${CLASS_PHASE1A}" "${CLASS_BERGHIA}" "${BERGHIA_FA}"; do
    [ -s "$f" ] || { log --level=ERROR "anchor_calibration: missing input $f"; exit 1; }
done

# --- 0. Anchor set (one-time prep; idempotent) ------------------------------
if [ ! -s references/anchors/anchor_set.fasta ]; then
    log "anchor_calibration: building anchor set"
    python3 scripts/build_anchor_set.py --out-dir references/anchors
fi

# --- 1. Optional focal-species subsample for a faster pilot ------------------
berghia_in="${BERGHIA_FA}"
if [ "${PILOT_BERGHIA_MAX}" -gt 0 ]; then
    berghia_in="${CAL}/berghia_pilot.fa"
    python3 - "$BERGHIA_FA" "$berghia_in" "$PILOT_BERGHIA_MAX" <<'PY'
import sys, random
from Bio import SeqIO
src, dst, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
recs = list(SeqIO.parse(src, "fasta"))
random.Random(12345).shuffle(recs)
SeqIO.write(recs[:n], dst, "fasta")
print(f"[pilot] focal-species subsample: {min(n,len(recs))}/{len(recs)}", file=sys.stderr)
PY
fi

# --- 2. Build the WITH-all-anchors per-class pools (all 4 classes) -----------
log "anchor_calibration: building per-class pools (CD-HIT ${REF_CLUSTER_IDENTITY}, floor ${REF_FLOOR_PER_SPECIES}, cap ${REF_CAP_PER_TAXON})"
python3 scripts/build_per_class_reference_pools.py \
    --scan-fasta-glob "${SCAN_DIR}/*.chemo_candidates.fa" \
    --class-tsv "${CLASS_PHASE1A}" \
    --berghia-fasta "${berghia_in}" \
    --berghia-class-tsv "${CLASS_BERGHIA}" \
    --anchor-fasta references/anchors/anchor_set.fasta \
    --anchor-tsv references/anchors/anchor_set.tsv \
    --out-dir "${CAL}/pools" \
    --threads "${CPUS}" \
    --force

# --- 3. Per class: derive WITHOUT pool, build both trees, evaluate -----------
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

for CLS in A B C F; do
    with_fa="${CAL}/pools/refs_class_${CLS}.fa"
    [ -s "$with_fa" ] || { log "anchor_calibration: class ${CLS} pool empty, skipping"; continue; }

    # WITHOUT out-group anchors = remove tier-2/3 anchor records (identical in-group)
    without_fa="${CAL}/pools/refs_class_${CLS}_noOG.fa"
    awk '/^>/{og=($0 ~ /^>ANCHOR_[A-F]_[23]_/); if(!og)print; next}{if(!og)print}' \
        "$with_fa" > "$without_fa"

    log "anchor_calibration: class ${CLS} — building WITH / WITHOUT trees"
    build_tree "$with_fa"    "${CAL}/trees/with_class_${CLS}"
    build_tree "$without_fa" "${CAL}/trees/without_class_${CLS}"

    with_tree="${CAL}/trees/with_class_${CLS}.treefile"
    without_tree="${CAL}/trees/without_class_${CLS}.treefile"
    if [ -s "$with_tree" ] && [ -s "$without_tree" ]; then
        python3 scripts/evaluate_anchor_divergence.py \
            --class "${CLS}" \
            --tree-without "$without_tree" \
            --tree-with "$with_tree" \
            --pool-members "${CAL}/pools/pool_members_class_${CLS}.tsv" \
            --out "${CAL}/verdicts/verdict_class_${CLS}.json"
    else
        log --level=WARN "anchor_calibration: class ${CLS} tree(s) missing, no verdict"
    fi
done

# --- 4. Aggregate verdicts into a report (report-only) ----------------------
python3 scripts/aggregate_anchor_verdicts.py \
    --in-dir "${CAL}/verdicts" \
    --out-json "${CAL}/anchor_calibration_report.json" \
    --out-md "${CAL}/anchor_calibration_report.md"

log "anchor_calibration: DONE — review ${CAL}/anchor_calibration_report.md"
