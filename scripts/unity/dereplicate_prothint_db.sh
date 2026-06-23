#!/bin/bash
# dereplicate_prothint_db.sh — de-redundify the combined ProtHint evidence DB
# (OrthoDB v12 Metazoa + MolluscaGenes mollusca_aa) with mmseqs2 linclust.
#
# WHY: mollusca_aa is RNAseq/transcriptome-derived (EvidentialGene tr2aacds on
# assembled transcripts) — expression-supported real coding sequences, so the
# sequences are GOOD evidence. The only problem is REDUNDANCY (transcript isoforms
# + near-identical orthologs across 200+ species), which saturates ProtHint's
# DIAMOND -k25 cap: for a mollusc genome the top-25 hits fill with duplicate
# near-identical copies, crowding out diverse independent orthologs and collapsing
# the cross-species agreement ProtHint uses to call high-confidence hints (this is
# why the raw combined DB HURT on Nautilus: 98 -> ~12 high-conf hints).
#
# Clustering at --min-seq-id 0.9 collapses isoforms + near-identical copies into
# one representative each, so the top-25 per seed are DIVERSE (independent), while
# keeping genuine in-clade diversity (orthologs >10% divergent survive — incl. the
# divergent mollusc-specific / lineage-expanded genes we care about).
#
# Output: <OUTDIR>/<TAG>_rep_seq.fasta = the de-redundified DB -> new PROTEIN_DB_ABS.
#
# Submit (defaults shown; override via --export):
#   sbatch scripts/unity/dereplicate_prothint_db.sh
#   sbatch --export=ALL,MINID=0.95 scripts/unity/dereplicate_prothint_db.sh   # conservative
#
#SBATCH --job-name=derep_prothint_db
#SBATCH --partition=cpu
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --output=logs/derep_prothint_db-%j.out
#SBATCH --error=logs/derep_prothint_db-%j.err
#SBATCH --account=pi_pkatz_umass_edu

set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
set -u

REPO_ROOT="${REPO_ROOT:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
cd "${REPO_ROOT}"

IN="${IN:-${REPO_ROOT}/species_tree_data/orthodb/odb12_metazoa_plus_mollusca_aa.fa}"
MINID="${MINID:-0.9}"
COV="${COV:-0.8}"
OUTDIR="${OUTDIR:-${REPO_ROOT}/species_tree_data/orthodb}"
TAG="${TAG:-odb12_metazoa_plus_mollusca_aa.linclust${MINID}}"
# Temp dir for mmseqs scratch. Use a script-local name rooted under OUTDIR and do
# NOT honor the environment's $TMP: on Unity compute nodes $TMP=/tmp, which both
# puts mmseqs scratch on the small shared node disk AND turned the cleanup below
# into `rm -rf /tmp` (it only failed thanks to root-owned /tmp permissions).
MMTMP="${OUTDIR}/mmseqs_tmp_${SLURM_JOB_ID:-derep}"
CPUS="${SLURM_CPUS_PER_TASK:-32}"

[ -s "$IN" ] || { echo "ERROR: input DB missing: $IN" >&2; exit 1; }
[ -n "$OUTDIR" ] || { echo "ERROR: OUTDIR is empty — refusing (would root mmseqs temp at /)" >&2; exit 1; }
mkdir -p "$OUTDIR" "$MMTMP" logs

echo "[$(date +%T)] linclust: $IN"
echo "  -> ${OUTDIR}/${TAG}  (min-seq-id=$MINID cov=$COV cov-mode=1 cpus=$CPUS)"
mmseqs easy-linclust "$IN" "${OUTDIR}/${TAG}" "$MMTMP" \
    --min-seq-id "$MINID" -c "$COV" --cov-mode 1 \
    --split-memory-limit 200G --threads "$CPUS"

REP="${OUTDIR}/${TAG}_rep_seq.fasta"
[ -s "$REP" ] || { echo "ERROR: rep fasta not produced: $REP" >&2; exit 1; }
n_out=$(grep -c '^>' "$REP")
echo "[$(date +%T)] DONE. representatives=${n_out}  -> ${REP}"
echo "Set as ProtHint DB:  PROTEIN_DB_ABS=${REP}"
# Clean up scratch, but never rm outside OUTDIR (guards against the env-$TMP bug).
case "${MMTMP}" in
    "${OUTDIR}"/*) rm -rf -- "${MMTMP}" ;;
    *) echo "WARN: refusing to remove unexpected temp dir '${MMTMP}'" >&2 ;;
esac
