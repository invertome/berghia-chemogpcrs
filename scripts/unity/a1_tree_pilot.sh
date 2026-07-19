#!/bin/bash
#SBATCH --job-name=a1_tree_pilot
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# A1 dedicated class-A tree — filter-stack PILOT (bead: A1 dedicated tree).
#
# Builds the reference frame for A1's phylogeny-residualized novelty: the 790
# class-A candidates + the 953 class-A anchor_set_PROD anchors (== the frame the
# proteinclip3b+protrek novelty is scored against). Aligns with MAFFT-DASH
# E-INS-i, original seqs only (DASH structural prior as a GUIDE; homologs are not
# retained). Then EVALUATES whether PREQUAL / TAPER help:
#
#   Conditions (2x2): PREQUAL {off,on} x TAPER {off,on}; ClipKit(smart-gap) held
#   constant. CLOAK is intentionally EXCLUDED — it is a consensus-across-many-MSAs
#   masker (needs the K=5 ensemble) and does not apply to a single E-INS-i
#   alignment; it also masks the alignment-uncertain divergent columns that A1's
#   novelty signal lives in (the same reason the pipeline keeps the un-CLOAKed
#   canonical alignment for 04b ECL analysis).
#
# FastTree per condition (fast) -> a1_compare_filter_trees.py quantifies the
# candidate->nearest-anchor distance-RANK stability (the A1-relevant quantity)
# and mean support. The FINAL tree (IQ-TREE 3) is built by a1_tree_final.sh on
# whichever alignment the pilot selects.
#
# Env (all have defaults): REPO, CAND_FASTA, ANCHOR_FASTA, OUTDIR, THREADS.
set -eo pipefail

REPO="${REPO:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
CAND_FASTA="${CAND_FASTA:-${REPO}/results/chemogpcrs/chemogpcrs_berghia_classA.fa}"
ANCHOR_FASTA="${ANCHOR_FASTA:-${REPO}/references/anchors/derived/anchor_set_PROD_classA.fasta}"
OUTDIR="${OUTDIR:-${REPO}/results/phylogenies/protein/class_A_a1}"
THREADS="${THREADS:-24}"          # per-alignment; two E-INS-i run in parallel

# Activate conda BEFORE set -u (CONDA_BACKUP_* unbound-var trap on deactivate hooks)
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
cd "$REPO"
source config.sh
source functions.sh

mkdir -p "$OUTDIR"
INPUT="${OUTDIR}/a1_input.fa"

echo "[a1_pilot] candidates=$(grep -c '^>' "$CAND_FASTA")  anchors=$(grep -c '^>' "$ANCHOR_FASTA")"
cat "$CAND_FASTA" "$ANCHOR_FASTA" > "$INPUT"
# Self-contained duplicate-id guard (the checkout's functions.sh may predate
# assert_no_duplicate_fasta_ids): candidate ids and ANCHOR_ ids are disjoint by
# construction, so any collision is a real error.
_dupes=$(grep '^>' "$INPUT" | sed 's/^>//; s/[[:space:]].*//' | sort | uniq -d)
if [ -n "$_dupes" ]; then
    echo "[a1_pilot] ERROR: duplicate sequence ids in combined input:"; echo "$_dupes" | head
    exit 1
fi
echo "[a1_pilot] combined input = $(grep -c '^>' "$INPUT") seqs -> $INPUT"

MAFFT_EINSI_DASH=(mafft --dash --originalseqonly --genafpair --maxiterate 1000)

# --- Two E-INS-i alignments in parallel: canonical (raw) and PREQUAL-masked ---
(
  echo "[a1_pilot] E-INS-i canonical (raw) ..."
  "${MAFFT_EINSI_DASH[@]}" --thread "$THREADS" "$INPUT" > "${OUTDIR}/canonical.fa" \
    2> "${OUTDIR}/canonical.mafft.log" && echo "[a1_pilot] canonical DONE" \
    || echo "[a1_pilot] canonical FAILED"
) &
(
  echo "[a1_pilot] PREQUAL -> E-INS-i ..."
  if bash "${SCRIPTS_DIR}/run_prequal.sh" --input="$INPUT" \
        --output="${OUTDIR}/prequal_masked.fa" 2> "${OUTDIR}/prequal.log"; then
    "${MAFFT_EINSI_DASH[@]}" --thread "$THREADS" "${OUTDIR}/prequal_masked.fa" \
      > "${OUTDIR}/prequal.fa" 2> "${OUTDIR}/prequal.mafft.log" \
      && echo "[a1_pilot] prequal DONE" || echo "[a1_pilot] prequal-align FAILED"
  else
    echo "[a1_pilot] PREQUAL FAILED (skipping the prequal arm)"
  fi
) &
wait
echo "[a1_pilot] alignments complete"

# --- Build the 4 condition alignments: {canonical,prequal} x {noTAPER,TAPER} ---
declare -A COND
for base in canonical prequal; do
  aln="${OUTDIR}/${base}.fa"
  [ -s "$aln" ] || { echo "[a1_pilot] $base alignment missing; skipping"; continue; }
  # TAPER-off branch
  COND["${base}"]="$aln"
  # TAPER-on branch
  if bash "${SCRIPTS_DIR}/run_taper.sh" --input="$aln" \
        --output="${OUTDIR}/${base}_taper.fa" 2> "${OUTDIR}/${base}_taper.log"; then
    COND["${base}_taper"]="${OUTDIR}/${base}_taper.fa"
  else
    echo "[a1_pilot] TAPER on $base FAILED"
  fi
done

# --- ClipKit (smart-gap) + FastTree per condition ---
TREE_ARGS=()
for name in "${!COND[@]}"; do
  aln="${COND[$name]}"
  trim="${OUTDIR}/${name}_trim.fa"
  clipkit "$aln" -m smart-gap -o "$trim" 2> "${OUTDIR}/${name}_clipkit.log" || cp "$aln" "$trim"
  ncol=$(awk '/^>/{next}{print length($0); exit}' "$trim")
  echo "[a1_pilot] $name: $(grep -c '^>' "$trim") seqs, ~${ncol} cols -> FastTree"
  FastTree -lg "$trim" > "${OUTDIR}/${name}.nwk" 2> "${OUTDIR}/${name}.fasttree.log" \
    && TREE_ARGS+=(--tree "${name}=${OUTDIR}/${name}.nwk") \
    || echo "[a1_pilot] FastTree $name FAILED"
done

# --- Reference (anchor) ids for the distance metric + comparison verdict ---
grep '^>' "$ANCHOR_FASTA" | sed 's/^>//; s/ .*//' > "${OUTDIR}/anchor_ids.txt"
python3 "${SCRIPTS_DIR}/a1_compare_filter_trees.py" \
  "${TREE_ARGS[@]}" \
  --ref-ids "${OUTDIR}/anchor_ids.txt" \
  --baseline canonical \
  --out "${OUTDIR}/filter_comparison.json"

echo "[a1_pilot] DONE. Verdict -> ${OUTDIR}/filter_comparison.json"
echo "[a1_pilot] alignments kept: canonical.fa, prequal.fa (+ _taper, _trim variants)"
