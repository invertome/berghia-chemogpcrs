#!/bin/bash
# calibrate_genome_track_margin.sh — empirical GENOME_TRACK_MIN_MARGIN
# calibration (bead berghia-chemogpcrs-elh; Task 5 of the margin-calibration
# plan).
#
# Formalizes the ad-hoc calibrate_genome_track.sh (staged Unity-side at
# /scratch3/.../jorge/calibrate_genome_track.sh for the earlier id/cov
# calibration, job 61394938) as a tracked wrapper for the margin gate:
#   (a) selects the RefSeq GPCRs with the pipeline's own HMM detector
#       (identify_gpcr_candidates — HMM layer only, no TMbed; functions.sh),
#   (b) subsets the RefSeq CDS to those genes by protein_id and maps ONLY that
#       subset to the genome with minimap2 -x splice + gmap -f gff3_gene,
#       REUSING the id/cov calibration's prebuilt gmap index, and
#   (c) runs scripts/calibrate_genome_track_margin.py in --mode refseq
#       (on-target) and --mode busco (cross-check, reusing the id/cov
#       calibration's on-disk alignments), writing tables + figures +
#       sourceable recommendations under
#       results/reconciliation/margin_calibration/.
#
# Design: docs/plans/2026-07-05-genome-track-margin-calibration-design.md (S3a, S7)
# Plan:   docs/plans/2026-07-05-genome-track-margin-calibration-plan.md (Task 5)
#
# No config.sh edit here — Task 6 (USER-GATE) reads the recommendation files
# this script writes and picks GENOME_TRACK_MIN_MARGIN (± a tier split) by hand.
#
# Usage (submit from the repo root on Unity):
#   sbatch scripts/unity/calibrate_genome_track_margin.sh
#
#SBATCH --job-name=gtrack_margin
#SBATCH --partition=cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --output=logs/gtrack_margin-%j.out
#SBATCH --error=logs/gtrack_margin-%j.err

# set -eo pipefail ONLY — NOT set -u. config.sh/functions.sh reference several
# optionally-set vars (e.g. PROVENANCE_FILE, guarded with `[ -z ... ] && return`
# — safe without -u, fatal with it), and set -u crashed this exact calibration
# lineage once already (job 61394938's first attempt). See
# 02c_genome_reconcile.sh's own comment on the same point.
set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr

REPO_ROOT="${REPO_ROOT:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
cd "${REPO_ROOT}"
source "${REPO_ROOT}/config.sh"
source "${REPO_ROOT}/functions.sh"

# Honor the SLURM allocation for tool threads (falls back to CPUS for local runs).
THREADS="${SLURM_CPUS_PER_TASK:-${CPUS:-8}}"

CAL="${RESULTS_DIR}/reconciliation/calibration"          # id/cov calibration artifacts (reused, not rebuilt)
OUT="${RESULTS_DIR}/reconciliation/margin_calibration"    # this script's own outputs
GMAP_DB_NAME="berghia_genome"                             # matches 02c_genome_reconcile.sh + the id/cov calibration build
GPCR_FA="$OUT/refseq_gpcrs.fa"
GPCR_CDS="$OUT/refseq_gpcr_cds.fna"

mkdir -p "$OUT" "${LOGS_DIR}"

# --- Pre-flight: required inputs + the reusable id/cov calibration artifacts. ---
check_file "$BERGHIA_GENOME_PROTEINS" "$GENOME" "$GENOME_CDS"
check_dir "$CAL/gmap_index"
check_file "$CAL/minimap2.paf" "$CAL/gmap.gff3" "$CAL/busco_single_copy.ids"

# --- (a) Select the RefSeq GPCRs and subset their CDS by protein_id. ---
# HMM layer only, no TMbed — design S3a wants paralog-rich GPCR genes with
# known loci, not the strict >=6-TM chemoreceptor definition.
log "Selecting RefSeq GPCRs via identify_gpcr_candidates (HMM layer only)."
identify_gpcr_candidates "$BERGHIA_GENOME_PROTEINS" "$GPCR_FA"
if [ ! -s "$GPCR_FA" ]; then
    log --level=ERROR "identify_gpcr_candidates produced no GPCR-positive RefSeq proteins; nothing to calibrate on."
    exit 1
fi
grep '^>' "$GPCR_FA" | sed 's/^>//; s/ .*//' > "$OUT/gpcr_prot_ids.txt"
log "RefSeq GPCR proteins selected: $(wc -l < "$OUT/gpcr_prot_ids.txt")"

# CDS-subset join key: the [protein_id=XP_...] attribute — verified 1:1 against
# the real cds_from_genomic.fna (every CDS record carries exactly one
# protein_id, no duplicates across the file), so a direct bracket-field match
# is both exact and simpler than stripping the trailing _<n> CDS-instance
# counter off the lcl|<scaffold>_cds_<protid>_<n> id prefix.
awk '
    NR == FNR { keep_id[$1] = 1; next }
    /^>/ {
        keep = 0
        if (match($0, /\[protein_id=[^]]+\]/)) {
            tag = substr($0, RSTART, RLENGTH)
            sub(/^\[protein_id=/, "", tag)
            sub(/\]$/, "", tag)
            if (tag in keep_id) keep = 1
        }
    }
    keep { print }
' "$OUT/gpcr_prot_ids.txt" "$GENOME_CDS" > "$GPCR_CDS"

if [ ! -s "$GPCR_CDS" ]; then
    log --level=ERROR "CDS-subset join matched zero records between ${GPCR_FA} and ${GENOME_CDS}"
    log --level=ERROR "(a protein_id header-format drift?). Aborting before mapping."
    exit 1
fi
log "RefSeq GPCR CDS subset: $(grep -c '^>' "$GPCR_CDS") records -> ${GPCR_CDS}"

# --- (b) Map ONLY the GPCR CDS subset, reusing the id/cov calibration index. ---
# Same minimap2/gmap invocation as 02c_genome_reconcile.sh; no gmap_build here.
log "Mapping the GPCR CDS subset to the genome (minimap2 -x splice + gmap)."
run_command "margin_minimap2" --stdout="$OUT/gpcr_cds.paf" \
    "$MINIMAP2" -x splice -t "$THREADS" "$GENOME" "$GPCR_CDS"
run_command "margin_gmap" --stdout="$OUT/gpcr_cds.gff3" \
    "$GMAP" -D "$CAL/gmap_index" -d "$GMAP_DB_NAME" -f gff3_gene -t "$THREADS" "$GPCR_CDS"

# --- (c) refseq-mode (on-target) + busco-mode (cross-check, on-disk). ---
# busco-mode reuses the id/cov calibration's on-disk alignments verbatim.
log "Running the margin calibration — refseq mode (on-target)."
python3 "${SCRIPTS_DIR}/calibrate_genome_track_margin.py" --mode refseq \
    --minimap2-paf "$OUT/gpcr_cds.paf" --gmap-gff "$OUT/gpcr_cds.gff3" \
    --refseq-cds "$GPCR_CDS" \
    --out-recommendation "$OUT/margin_recommendation_refseq.sh" \
    --out-figure "$OUT/margin_roc_refseq.png" | tee "$OUT/margin_table_refseq.txt"

log "Running the margin calibration — busco mode (cross-check, on-disk)."
python3 "${SCRIPTS_DIR}/calibrate_genome_track_margin.py" --mode busco \
    --minimap2-paf "$CAL/minimap2.paf" --gmap-gff "$CAL/gmap.gff3" \
    --busco-ids "$CAL/busco_single_copy.ids" \
    --out-recommendation "$OUT/margin_recommendation_busco.sh" \
    --out-figure "$OUT/margin_roc_busco.png" | tee "$OUT/margin_table_busco.txt"

log "Margin calibration complete -> ${OUT}"
