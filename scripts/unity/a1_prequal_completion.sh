#!/bin/bash
#SBATCH --job-name=a1_prequal_completion
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=256G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# A1 filter-stack pilot — COMPLETION of the PREQUAL arm (bead: A1 dedicated tree).
#
# WHY THIS EXISTS. The pilot (61970841) runs both E-INS-i alignments concurrently
# under ONE allocation. Its own #SBATCH asks for 128G, but the job was submitted
# with a `--mem=96G` command-line override, which wins. The two arms shared 96G
# and the PREQUAL arm's MAFFT was OOM-killed in iterative refinement:
#
#   mafft: line 2842: 352366 Killed   "$prefix/dvtditr" ...
#
# leaving prequal.fa at 0 bytes. The pilot's :90 guard then SKIPS the arm (cleanly
# -- it does not feed an empty alignment downstream), so filter_comparison.json
# covers only {canonical, canonical_taper}: it can answer "does TAPER help?" but
# carries ZERO evidence on "does PREQUAL help?" -- and PREQUAL is default-ON in
# production (RUN_PREQUAL=1).
#
# WHAT THIS DOES. Runs the PREQUAL arm ALONE (no concurrent second alignment) with
# the full allocation, then rebuilds the comparison over all FOUR conditions so the
# 2x2 {PREQUAL off,on} x {TAPER off,on} is complete. PREQUAL itself already
# succeeded in the pilot (prequal_masked.fa, 877K); it is deterministic, so this
# reuses that output rather than spending another ~2.5h regenerating it.
#
# Memory: 256G on a ~1TB node. The OOM was caused by sharing, not by an inherent
# 96G ceiling, but E-INS-i peak usage on 1743 seqs was never measured -- so this
# requests generously rather than guessing at the true peak.
#
# Submit AFTER the pilot ends (it writes canonical.nwk / canonical_taper.nwk only
# after `wait`, i.e. once the surviving canonical arm finishes):
#   sbatch --dependency=afterany:61970841 scripts/unity/a1_prequal_completion.sh
#
# Env (all have defaults): REPO, OUTDIR, ANCHOR_FASTA, THREADS, PILOT_JOBID.
set -eo pipefail

REPO="${REPO:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
OUTDIR="${OUTDIR:-${REPO}/results/phylogenies/protein/class_A_a1}"
ANCHOR_FASTA="${ANCHOR_FASTA:-${REPO}/references/anchors/derived/anchor_set_PROD_classA.fasta}"
THREADS="${THREADS:-48}"          # the whole allocation; this arm runs alone
PILOT_JOBID="${PILOT_JOBID:-61970841}"

# Activate conda BEFORE set -u (CONDA_BACKUP_* unbound-var trap on deactivate hooks)
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
set -u
cd "$REPO"

SCRIPTS_DIR="${REPO}/scripts"
MASKED="${OUTDIR}/prequal_masked.fa"

say() { echo "[a1_prequal] $*"; }
die() { echo "[a1_prequal] ERROR: $*" >&2; exit 1; }

# n_seqs: count FASTA records without the doubled-zero idiom (grep -c prints 0 AND
# exits 1 on no-match, so `|| echo 0` yields "0\n0" and breaks every comparison).
n_seqs() { local n; n=$(grep -c '^>' "$1" 2>/dev/null || true); echo "${n:-0}"; }

# --- Preconditions: fail loudly rather than produce a plausible-but-empty verdict
[ -s "$MASKED" ] || die "$MASKED missing/empty; PREQUAL output from the pilot is required"
N_IN=$(n_seqs "$MASKED")
[ "$N_IN" -gt 0 ] || die "$MASKED contains no FASTA records"
say "PREQUAL-masked input reused from the pilot: ${N_IN} seqs"

for t in canonical canonical_taper; do
    [ -s "${OUTDIR}/${t}.nwk" ] \
        || die "${OUTDIR}/${t}.nwk missing; the pilot (${PILOT_JOBID}) has not produced its canonical trees yet"
done
say "pilot canonical trees present: canonical.nwk, canonical_taper.nwk"

[ -s "$ANCHOR_FASTA" ] || die "$ANCHOR_FASTA missing (needed for the reference-id metric)"

# --- The PREQUAL arm, alone, with the full allocation -----------------------
MAFFT_EINSI_DASH=(mafft --dash --originalseqonly --genafpair --maxiterate 1000)

if [ -s "${OUTDIR}/prequal.fa" ]; then
    say "prequal.fa already present and non-empty; skipping the realignment"
else
    say "E-INS-i on PREQUAL-masked input (${THREADS} threads, alone) ..."
    "${MAFFT_EINSI_DASH[@]}" --thread "$THREADS" "$MASKED" \
        > "${OUTDIR}/prequal.fa" 2> "${OUTDIR}/prequal.mafft.log" \
        || die "prequal realignment FAILED; tail of the log:
$(tail -5 "${OUTDIR}/prequal.mafft.log" 2>/dev/null | tr '\r' '\n')"
    say "prequal alignment DONE"
fi

# MAFFT must return every input sequence. A short alignment here means a silent
# truncation, which would bias the whole PREQUAL condition.
N_ALN=$(n_seqs "${OUTDIR}/prequal.fa")
[ "$N_ALN" -eq "$N_IN" ] \
    || die "prequal.fa has ${N_ALN} seqs but the input had ${N_IN}; refusing to compare a truncated alignment"
say "prequal.fa verified: ${N_ALN}/${N_IN} seqs retained"

# --- TAPER branch ----------------------------------------------------------
# NOTE: the checkout's run_taper.sh predates the doubled-zero sweep, but that
# defect is confined to its provenance block (after the failure guard) and cannot
# affect the masking. The counts are asserted here instead of trusting its log.
declare -A COND
COND[prequal]="${OUTDIR}/prequal.fa"

if bash "${SCRIPTS_DIR}/run_taper.sh" --input="${OUTDIR}/prequal.fa" \
        --output="${OUTDIR}/prequal_taper.fa" 2> "${OUTDIR}/prequal_taper.log"; then
    N_TAP=$(n_seqs "${OUTDIR}/prequal_taper.fa")
    if [ "$N_TAP" -eq "$N_IN" ]; then
        COND[prequal_taper]="${OUTDIR}/prequal_taper.fa"
        say "TAPER on prequal DONE: ${N_TAP} seqs"
    else
        say "WARN: prequal_taper.fa has ${N_TAP} seqs (expected ${N_IN}); dropping the condition"
    fi
else
    say "WARN: TAPER on prequal FAILED; continuing with the TAPER-off condition only"
fi

# --- ClipKit + FastTree per new condition ----------------------------------
for name in "${!COND[@]}"; do
    aln="${COND[$name]}"
    trim="${OUTDIR}/${name}_trim.fa"
    clipkit "$aln" -m smart-gap -o "$trim" 2> "${OUTDIR}/${name}_clipkit.log" || cp "$aln" "$trim"
    # Sum the first record's residue lines: taking only the first line reports
    # the FASTA wrap width (60), not the alignment width.
    ncol=$(awk '/^>/{if(seen) exit; seen=1; next} seen{n+=length($0)} END{print n+0}' "$trim")
    say "${name}: $(n_seqs "$trim") seqs, ~${ncol} cols -> FastTree"
    FastTree -lg "$trim" > "${OUTDIR}/${name}.nwk" 2> "${OUTDIR}/${name}.fasttree.log" \
        || say "WARN: FastTree ${name} FAILED"
done

# --- Rebuild the verdict over the complete 2x2 -----------------------------
TREE_ARGS=()
for name in canonical canonical_taper prequal prequal_taper; do
    [ -s "${OUTDIR}/${name}.nwk" ] && TREE_ARGS+=(--tree "${name}=${OUTDIR}/${name}.nwk")
done
say "conditions entering the comparison: ${#TREE_ARGS[@]} of 4 (2 args each -> $((${#TREE_ARGS[@]} / 2)) trees)"
[ "${#TREE_ARGS[@]}" -ge 4 ] || die "fewer than 2 trees available; nothing to compare"

grep '^>' "$ANCHOR_FASTA" | sed 's/^>//; s/ .*//' > "${OUTDIR}/anchor_ids.txt"

# Preserve the pilot's TAPER-only verdict as evidence rather than overwriting it.
PARTIAL="${OUTDIR}/filter_comparison.json"
if [ -s "$PARTIAL" ]; then
    cp -n "$PARTIAL" "${OUTDIR}/filter_comparison_partial_${PILOT_JOBID}.json" || true
    say "pilot's partial verdict preserved as filter_comparison_partial_${PILOT_JOBID}.json"
fi

python3 "${SCRIPTS_DIR}/a1_compare_filter_trees.py" \
    "${TREE_ARGS[@]}" \
    --ref-ids "${OUTDIR}/anchor_ids.txt" \
    --baseline canonical \
    --out "${OUTDIR}/filter_comparison.tmp.json"
mv "${OUTDIR}/filter_comparison.tmp.json" "$PARTIAL"

say "DONE. Complete 2x2 verdict -> ${PARTIAL}"
