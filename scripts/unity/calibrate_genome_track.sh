#!/bin/bash
# calibrate_genome_track.sh — empirical GENOME_TRACK_MIN_ID / _MIN_COV
# calibration (bead berghia-chemogpcrs-elh; Task 8 of the genome-track
# reconciliation plan).
#
# Tracked formalization of the ad-hoc wrapper that was staged Unity-side at
# /scratch3/.../jorge/calibrate_genome_track.sh and run as job 61394938 (the
# id/cov calibration). Its outputs seed the margin calibration too — the
# tracked calibrate_genome_track_margin.sh REUSES the gmap index + alignments
# this script writes under results/reconciliation/calibration/.
#
# Method (design 2026-06-30 S4/S9): map the Berghia BUSCO single-copy
# transcripts to the Berghia genome = the "true same-gene" %id/%cov floor; map
# the stage-02 chemoreceptor candidate transcripts = the paralog cross-map %id
# ceiling; scripts/calibrate_reconcile_thresholds.py recommends a bar between
# them. Alignments are produced EXACTLY as stage 02c does (minimap2 -x splice +
# gmap -f gff3_gene) over the UNION of both id sets, then the calibrator selects
# each distribution by id from the combined placement pool.
#
# Prereq: gmap/gmap_build/seqtk/minimap2 in the berghia-gpcr env
#   -> if a tool is missing: sbatch scripts/unity/install_genome_track_tools.sh
#
# No config.sh edit here — the USER-GATE reads the recommendation files this
# script writes (results/reconciliation/calibrated_cutoffs.sh +
# separation_table.txt) and adopts GENOME_TRACK_MIN_ID/_MIN_COV by hand
# (kept 95/90 when there is no clean separation, as job 61394938 found).
#
# Usage (submit from the repo root on Unity):
#   sbatch scripts/unity/calibrate_genome_track.sh
#
#SBATCH --job-name=gtrack_calib
#SBATCH --partition=cpu
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --output=logs/gtrack_calib-%j.out
#SBATCH --error=logs/gtrack_calib-%j.err

# set -eo pipefail ONLY — NOT set -u. config.sh/functions.sh reference several
# optionally-set vars (e.g. PROVENANCE_FILE, guarded with `[ -z ... ] && return`
# — safe without -u, fatal with it), and set -u crashed this exact calibration
# once already (job 61394938's first attempt). Activate conda BEFORE any -u for
# the same reason (the CONDA_BACKUP_CXX deactivate-hook trap).
set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr

REPO_ROOT="${REPO_ROOT:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
cd "${REPO_ROOT}"
source "${REPO_ROOT}/config.sh"
source "${REPO_ROOT}/functions.sh"

# Honor the SLURM allocation for tool threads (falls back to CPUS for local runs).
THREADS="${SLURM_CPUS_PER_TASK:-${CPUS:-8}}"

# --- Resolve inputs (all verified present at merge 2e9b8a4) -------------------
SC_DIR="${RESULTS_DIR}/busco/single_copy/${BERGHIA_FILE_PREFIX}.aa"   # 3538 single-copy .faa
CAND_FA="${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa"             # 888 stage-02 candidates
MRNA="${BERGHIA_TRANSCRIPTOME_MRNA}"                                  # EvidentialGene .mrna (86026)
GENOME_FA="${GENOME}"                                                 # RefSeq genome FASTA

WORK="${RESULTS_DIR}/reconciliation/calibration"
OUT_CUTOFFS="${RESULTS_DIR}/reconciliation/calibrated_cutoffs.sh"
mkdir -p "${WORK}" "${RESULTS_DIR}/reconciliation" "${LOGS_DIR}"

# --- Pre-flight: tools + inputs ---------------------------------------------
for t in "${GMAP_BUILD}" "${GMAP}" "${MINIMAP2}" "${SEQTK}"; do
    command -v "$t" >/dev/null 2>&1 || {
        echo "ERROR: required tool '$t' not on PATH in env berghia-gpcr." >&2
        echo "       Install the genome-track tools first:" >&2
        echo "       sbatch scripts/unity/install_genome_track_tools.sh" >&2
        exit 1
    }
done
check_file "${GENOME_FA}" "${MRNA}" "${CAND_FA}"
check_dir "${SC_DIR}"

# --- (a) id lists ------------------------------------------------------------
# BUSCO single-copy transcript ids: first token of each .faa header.
BUSCO_IDS="${WORK}/busco_single_copy.ids"
grep -h '^>' "${SC_DIR}"/*.faa | sed 's/^>//' | awk '{print $1}' | sort -u > "${BUSCO_IDS}"
# Paralog-family (chemoreceptor candidate) transcript ids.
CAND_IDS="${WORK}/candidates.ids"
grep '^>' "${CAND_FA}" | sed 's/^>//' | awk '{print $1}' | sort -u > "${CAND_IDS}"
log "BUSCO single-copy ids: $(wc -l < "${BUSCO_IDS}")"
log "candidate ids:        $(wc -l < "${CAND_IDS}")"

# --- (b) subset the .mrna to the UNION (align both id sets once) --------------
SUBSET_IDS="${WORK}/subset.ids"
SUBSET_MRNA="${WORK}/subset.mrna"
cat "${BUSCO_IDS}" "${CAND_IDS}" | sort -u > "${SUBSET_IDS}"
run_command "seqtk_subseq" --stdout="${SUBSET_MRNA}" \
    ${SEQTK} subseq "${MRNA}" "${SUBSET_IDS}"
log "subset transcripts:   $(grep -c '^>' "${SUBSET_MRNA}")"

# --- (c) build the GMAP index once (same as stage 02c) -----------------------
GMAP_DB_DIR="${WORK}/gmap_index"
GMAP_DB_NAME="berghia_genome"                            # matches 02c_genome_reconcile.sh
mkdir -p "${GMAP_DB_DIR}"
run_command "gmap_build" ${GMAP_BUILD} -D "${GMAP_DB_DIR}" -d "${GMAP_DB_NAME}" -t "${THREADS}" "${GENOME_FA}"

# --- (d) transcript->genome alignments (exactly as stage 02c) ----------------
MINIMAP_PAF="${WORK}/minimap2.paf"
GMAP_GFF="${WORK}/gmap.gff3"
run_command "minimap2_splice" --stdout="${MINIMAP_PAF}" \
    ${MINIMAP2} -x splice -t "${THREADS}" "${GENOME_FA}" "${SUBSET_MRNA}"
run_command "gmap_align" --stdout="${GMAP_GFF}" \
    ${GMAP} -D "${GMAP_DB_DIR}" -d "${GMAP_DB_NAME}" -f gff3_gene -t "${THREADS}" "${SUBSET_MRNA}"

# --- (e) calibrate -----------------------------------------------------------
# calibrate_reconcile_thresholds.py `import reconcile_candidates` -> scripts/ on PYTHONPATH.
SEP_TABLE="${WORK}/separation_table.txt"
PYTHONPATH="${SCRIPTS_DIR}${PYTHONPATH:+:${PYTHONPATH}}" \
python3 "${SCRIPTS_DIR}/calibrate_reconcile_thresholds.py" \
    --busco-single-copy-ids "${BUSCO_IDS}" \
    --paralog-ids "${CAND_IDS}" \
    --minimap2-paf "${MINIMAP_PAF}" \
    --gmap-gff "${GMAP_GFF}" \
    --out "${OUT_CUTOFFS}" \
    | tee "${SEP_TABLE}"

log "id/cov calibration complete."
log "separation table  : ${SEP_TABLE}  (also echoed above / in the .out log)"
log "recommended cutoffs: ${OUT_CUTOFFS}"
