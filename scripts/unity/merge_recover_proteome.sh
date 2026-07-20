#!/bin/bash
# merge_recover_proteome.sh — standard recovery procedure for samples that fail
# primary BRAKER4 annotation: TSEBRA-merge a pure-ES (ab-initio) gene set with
# the EP-fallback (protein-evidence) gene set on the SAME genome, then pick the
# best proteome empirically by BUSCO.
#
# Rationale + design: docs/plans/2026-06-19-failed-annotation-merge-recovery.md
# Bead: berghia-chemogpcrs-u25.
#
# Parameterized (env vars; pass via `sbatch --export`):
#   SAMPLE        sample id, e.g. 34576_Nautilus_macromphalus
#   ES_DIR        pure-ES run output dir (has braker.gtf, hintsfile.gff)
#   EP_DIR        EP-fallback run output dir (has braker.gtf, hintsfile_iter2.gff,
#                 prothint_hints_iter2.gff)
#   GENOME        genome FASTA used by BOTH runs (asserted identical)
#   OUTDIR        where to write merged gene sets, proteomes, BUSCO, and report
#   CACHE_WINNER  (optional) es|ep|merge_keep|merge_filter — if set, cache that
#                 candidate's ALL-TRANSCRIPT proteome (+CDS) into proteomes_braker4/
#                 as the standardized recovery endpoint. Unset = stop at the report.
#
# Default flow STOPS at the report so the user picks the winner from the BUSCO
# table; re-run with CACHE_WINNER=<candidate> to cache (feeds the all-557 run).
#
# Nautilus test invocation:
#   R=/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
#   B=$R/species_tree_data/braker4_run
#   sbatch --export=ALL,\
#   SAMPLE=34576_Nautilus_macromphalus,\
#   ES_DIR=$B/array_runs/nautilus_es/output/34576_Nautilus_macromphalus,\
#   EP_DIR=$B/array_runs/nautilus_ep_mollusca/output/34576_Nautilus_macromphalus,\
#   GENOME=$B/../braker4_genomes/34576_Nautilus_macromphalus.fasta,\
#   OUTDIR=$B/array_runs/nautilus_merged_recover/34576_Nautilus_macromphalus \
#     scripts/unity/merge_recover_proteome.sh
#
#SBATCH --job-name=merge_recover
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --output=logs/merge_recover-%j.out
#SBATCH --error=logs/merge_recover-%j.err
#SBATCH --account=pi_pkatz_umass_edu

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
cd "${REPO_ROOT}"

SAMPLE="${SAMPLE:?set SAMPLE}"
ES_DIR="${ES_DIR:?set ES_DIR}"
EP_DIR="${EP_DIR:?set EP_DIR}"
GENOME="${GENOME:?set GENOME}"
OUTDIR="${OUTDIR:?set OUTDIR}"
CPUS="${SLURM_CPUS_PER_TASK:-16}"
MAX_AA="${MAX_AA:-50000}"        # drop proteins longer than this (HMMER-abort guard)
TTABLE="${TTABLE:-1}"            # nuclear genetic code (BRAKER4 default)

SIF="${REPO_ROOT}/species_tree_data/braker4_run/.singularity_cache"
BRAKER_IMG="${SIF}/dbdad69110452498a7d3d6879e51d31a.simg"   # teambraker/braker3:v3.0.10
BUSCO_IMG="${SIF}/5fe4eb5133b24329c9d34e853dcaa155.simg"    # ezlabgva/busco:v6.0.0_cv1
BUSCO_DL="${REPO_ROOT}/external/braker4/shared_data/busco_downloads"
LINEAGE="metazoa_odb12"
GETANNO="${REPO_ROOT}/external/braker4/scripts/getAnnoFastaFromJoingenes.py"

# singularity exec wrapper (matches the pipeline's bind/env convention)
SING() { singularity exec --no-mount resolv-conf -B "${REPO_ROOT}:${REPO_ROOT}" \
            --env PREPEND_PATH=/opt/conda/bin "$@"; }

ES_GTF="${ES_DIR}/braker.gtf"
EP_GTF="${EP_DIR}/braker.gtf"
log() { echo "[$(date +%H:%M:%S)] merge_recover: $*"; }

# --- 0. validate inputs ----------------------------------------------------
for f in "$ES_GTF" "$EP_GTF" "$GENOME" "$BRAKER_IMG" "$BUSCO_IMG" "$GETANNO"; do
    [ -s "$f" ] || { echo "ERROR: missing/empty: $f" >&2; exit 1; }
done
[ -f "${BUSCO_DL}/lineages/${LINEAGE}/dataset.cfg" ] || \
    { echo "ERROR: BUSCO lineage not staged offline: ${BUSCO_DL}/lineages/${LINEAGE}" >&2; exit 1; }
# Assert both runs annotated the SAME genome (coordinate merge requires it).
for run in "$ES_DIR" "$EP_DIR"; do
    sc="$(dirname "$(dirname "$run")")/samples.csv"
    if [ -f "$sc" ]; then
        g=$(awk -F, -v s="$SAMPLE" '$1==s{print $2}' "$sc" | head -1)
        # Compare RESOLVED paths so equivalent spellings that differ only in
        # representation (e.g. .../braker4_run/../braker4_genomes vs the canonical
        # .../braker4_genomes) don't trip the guard.
        if [ -n "$g" ] && [ "$(readlink -f "$g" 2>/dev/null)" != "$(readlink -f "$GENOME" 2>/dev/null)" ]; then
            echo "ERROR: $run samples.csv genome ($g) != GENOME ($GENOME)" >&2; exit 1
        fi
    fi
done
mkdir -p "${OUTDIR}" logs
log "sample=${SAMPLE} ES=$(grep -c $'\tgene\t' "$ES_GTF" 2>/dev/null || echo ?) EP=$(grep -c $'\tgene\t' "$EP_GTF" 2>/dev/null || echo ?) genes"

# --- 1. merge hints (cat, NOT join_mult_hints.pl — preserves grp= linkage) --
MERGED_HINTS="${OUTDIR}/merged_hints.gff"
: > "${MERGED_HINTS}"
for h in "${ES_DIR}/hintsfile.gff" "${EP_DIR}/hintsfile_iter2.gff" "${EP_DIR}/prothint_hints_iter2.gff"; do
    [ -s "$h" ] && cat "$h" >> "${MERGED_HINTS}"
done
[ -s "${MERGED_HINTS}" ] || { echo "ERROR: no hint files found to merge" >&2; exit 1; }
log "merged hints -> ${MERGED_HINTS} ($(wc -l < "${MERGED_HINTS}") lines)"

# --- 2. prefix gene/transcript IDs (avoid collisions between runs) ----------
ES_PFX="${OUTDIR}/es_prefixed.gtf"; EP_PFX="${OUTDIR}/ep_prefixed.gtf"
sed 's/gene_id "\([^"]*\)"/gene_id "es_\1"/g; s/transcript_id "\([^"]*\)"/transcript_id "es_\1"/g' "$ES_GTF" > "$ES_PFX"
sed 's/gene_id "\([^"]*\)"/gene_id "ep_\1"/g; s/transcript_id "\([^"]*\)"/transcript_id "ep_\1"/g' "$EP_GTF" > "$EP_PFX"

# --- 3. TSEBRA, both modes -------------------------------------------------
KEEP_GTF="${OUTDIR}/merged_keep.gtf"; FILTER_GTF="${OUTDIR}/merged_filter.gtf"
log "TSEBRA keep-mode (union-preserving)"
SING "$BRAKER_IMG" bash -c "export PATH=/opt/conda/bin:/opt/TSEBRA/bin:\$PATH; \
    tsebra.py --keep_gtf '${ES_PFX}','${EP_PFX}' --hintfiles '${MERGED_HINTS}' \
      --cfg /opt/TSEBRA/config/default.cfg --out '${KEEP_GTF}'"
log "TSEBRA filter-mode (evidence-filtering)"
SING "$BRAKER_IMG" bash -c "export PATH=/opt/conda/bin:/opt/TSEBRA/bin:\$PATH; \
    tsebra.py --gtf '${ES_PFX}','${EP_PFX}' --hintfiles '${MERGED_HINTS}' \
      --cfg /opt/TSEBRA/config/default.cfg --out '${FILTER_GTF}'"

# --- 4. longest isoform + proteins for all 4 candidates, drop >MAX_AA -------
# Candidate name -> source GTF
declare -A CAND=( [es]="$ES_GTF" [ep]="$EP_GTF" [merge_keep]="$KEEP_GTF" [merge_filter]="$FILTER_GTF" )
for name in es ep merge_keep merge_filter; do
    gtf="${CAND[$name]}"
    long_gtf="${OUTDIR}/${name}.longest.gtf"
    aa_raw="${OUTDIR}/${name}.longest.raw.aa"
    aa="${OUTDIR}/${name}.longest.aa"
    log "extract proteins: ${name}"
    SING "$BRAKER_IMG" bash -c "export PATH=/opt/conda/bin:/opt/TSEBRA/bin:\$PATH; \
        /opt/TSEBRA/bin/get_longest_isoform.py -g '${gtf}' -o '${long_gtf}' -q && \
        python3 '${GETANNO}' -g '${GENOME}' -f '${long_gtf}' -o '${OUTDIR}/${name}.longest' -t ${TTABLE}"
    [ -s "$aa" ] || { echo "ERROR: ${aa} not produced" >&2; exit 1; }
    # drop sequences > MAX_AA (collapses to single-line fasta; fine for BUSCO)
    mv "$aa" "$aa_raw"
    awk -v m="$MAX_AA" '/^>/{if(h!=""&&l<=m){print h; print s} h=$0; s=""; l=0; next}
                        {s=s $0; l+=length($0)} END{if(h!=""&&l<=m){print h; print s}}' "$aa_raw" > "$aa"
    n_raw=$(grep -c "^>" "$aa_raw"); n_keep=$(grep -c "^>" "$aa")
    log "${name}: ${n_keep} proteins (dropped $((n_raw-n_keep)) > ${MAX_AA} aa)"
done

# --- 5. BUSCO proteins mode (offline) on all 4 -----------------------------
BUSCO_ROOT="${OUTDIR}/busco"; mkdir -p "${BUSCO_ROOT}"
for name in es ep merge_keep merge_filter; do
    aa="${OUTDIR}/${name}.longest.aa"
    log "BUSCO proteins: ${name}"
    rm -rf "${BUSCO_ROOT}/${name}"
    SING "$BUSCO_IMG" bash -c "export PATH=/opt/conda/bin:\$PATH; \
        busco -i '${aa}' -o '${name}' --out_path '${BUSCO_ROOT}' \
          -l ${LINEAGE} -m proteins -c ${CPUS} \
          --download_path '${BUSCO_DL}' --offline" || \
        echo "WARN: BUSCO failed for ${name} (continuing)"
done

# --- 6. report TSV ---------------------------------------------------------
REPORT="${OUTDIR}/merge_recover_report.tsv"
printf 'candidate\tn_genes\tn_proteins\tC\tS\tD\tF\tM\tn\n' > "${REPORT}"
for name in es ep merge_keep merge_filter; do
    gtf="${CAND[$name]}"
    ng=$(grep -c $'\tgene\t' "$gtf" 2>/dev/null || true)
    ng=${ng:-0}
    np=$(grep -c "^>" "${OUTDIR}/${name}.longest.aa" 2>/dev/null || true)
    np=${np:-0}
    summ=$(find "${BUSCO_ROOT}/${name}" -name "short_summary*.txt" 2>/dev/null | head -1)
    line=$(grep -oE 'C:[0-9.]+%\[S:[0-9.]+%,D:[0-9.]+%\],F:[0-9.]+%,M:[0-9.]+%,n:[0-9]+' "${summ}" 2>/dev/null | head -1)
    if [ -n "$line" ]; then
        C=$(sed -E 's/C:([0-9.]+)%.*/\1/' <<<"$line"); S=$(sed -E 's/.*S:([0-9.]+)%.*/\1/' <<<"$line")
        D=$(sed -E 's/.*D:([0-9.]+)%.*/\1/' <<<"$line"); F=$(sed -E 's/.*F:([0-9.]+)%.*/\1/' <<<"$line")
        M=$(sed -E 's/.*M:([0-9.]+)%.*/\1/' <<<"$line"); N=$(sed -E 's/.*n:([0-9]+).*/\1/' <<<"$line")
    else C=NA; S=NA; D=NA; F=NA; M=NA; N=NA; fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$ng" "$np" "$C" "$S" "$D" "$F" "$M" "$N" >> "${REPORT}"
done
log "===== ${REPORT} ====="
cat "${REPORT}"

# --- 7. cache the chosen winner (gated by CACHE_WINNER) --------------------
# Standardized recovery endpoint. Set CACHE_WINNER=es|ep|merge_keep|merge_filter
# to write that candidate's ALL-TRANSCRIPT proteome (+CDS) into the Phase-1f
# proteome cache as ${SAMPLE}.{aa,cds}.fna — matching the braker.aa convention
# (all isoforms, no MAX_AA drop) used by every other cached species, so the
# recovered sample is not the lone deflated proteome in the all-557 orthology run.
# Leave CACHE_WINNER unset to stop at the report so the user picks the winner from
# the BUSCO table first.
if [ -n "${CACHE_WINNER:-}" ]; then
    case "${CACHE_WINNER}" in
        es|ep|merge_keep|merge_filter) ;;
        *) echo "ERROR: CACHE_WINNER must be es|ep|merge_keep|merge_filter (got '${CACHE_WINNER}')" >&2; exit 1 ;;
    esac
    win_gtf="${CAND[$CACHE_WINNER]}"
    CACHE_DIR="${CACHE_DIR:-${REPO_ROOT}/references/species_tree/cache/proteomes_braker4}"
    mkdir -p "${CACHE_DIR}"
    log "caching winner '${CACHE_WINNER}' (ALL transcripts) -> ${CACHE_DIR}/${SAMPLE}.{aa,cds}.fna"
    SING "$BRAKER_IMG" bash -c "export PATH=/opt/conda/bin:/opt/TSEBRA/bin:\$PATH; \
        python3 '${GETANNO}' -g '${GENOME}' -f '${win_gtf}' -o '${OUTDIR}/${CACHE_WINNER}.alltx' -t ${TTABLE}"
    [ -s "${OUTDIR}/${CACHE_WINNER}.alltx.aa" ] || { echo "ERROR: all-transcript proteome not produced" >&2; exit 1; }
    cp "${OUTDIR}/${CACHE_WINNER}.alltx.aa" "${CACHE_DIR}/${SAMPLE}.aa.fna"
    if [ -s "${OUTDIR}/${CACHE_WINNER}.alltx.codingseq" ]; then
        cp "${OUTDIR}/${CACHE_WINNER}.alltx.codingseq" "${CACHE_DIR}/${SAMPLE}.cds.fna"
    fi
    log "cached: ${CACHE_DIR}/${SAMPLE}.aa.fna ($(grep -c '^>' "${CACHE_DIR}/${SAMPLE}.aa.fna") proteins) + .cds.fna"
else
    log "DONE — review the table, then re-run with CACHE_WINNER=<candidate> to cache into proteomes_braker4/ (user-gated)."
fi
