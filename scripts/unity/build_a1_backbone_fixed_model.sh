#!/bin/bash
#SBATCH --job-name=a1_backbone_qpfam
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# A1 backbone tree search with the model ModelFinder ALREADY SELECTED.
#
# WHY THIS EXISTS. The full-sweep run (61979997) explored 444 of 451 candidate
# models over ~30 hours and was killed OUT_OF_MEMORY at exactly its 128G ceiling
# while fitting C50 -- a 50-component profile mixture over 953 taxa. The seven it
# never reached were C50, C60, EX2, EX3, EHO, LG4M and LG4X: ALL of them profile
# or mixture models from --madd. Every standard empirical matrix, with every
# rate-heterogeneity variant, was scored.
#
# THE SWEEP HAD ALREADY ANSWERED. Best five by BIC across the 444:
#
#     Q.PFAM+R10      1372893.550   <- selected
#     Q.PFAM+I+R10    1372894.018   (+0.47, an invariant-sites term buying nothing)
#     VT+R10          1376329.643   (+3436, next distinct matrix)
#     VT+I+R10        1376338.275
#     JTT+I+R10       1378001.496
#
# and the profile mixtures were the WORST models tested, not merely unselected:
#
#     C30             1420092.466
#     C40             1422061.144   <- WORSE than C30; the series is degrading
#     C20             1433443.907
#     C10             1444181.723
#
# C40 trails the winner by ~49,000 BIC units, and the C-series turns over between
# C30 and C40. For C50/C60 to win they would have to reverse a degrading trend AND
# close a gap forty times larger than the entire spread across the top five. So
# the untested seven cannot change the selection, and re-running the sweep with
# more memory would spend another day confirming a monotone trend.
#
# This does NOT abandon the mixture requirement (bead 1kvh). Mixtures were
# explored -- 444 models, the whole C-series through C40 -- and lost on the
# evidence. Q.PFAM is itself a profile-derived matrix estimated from Pfam
# alignments, so this is a matrix trained on protein-family diversity beating
# generic site-heterogeneity classes on a single-family alignment, which is a
# coherent result rather than an artifact of truncation.
#
# MEMORY. 192G against the 128G that died. The mixture was the thing that did not
# fit -- a fixed matrix with free rates instantiates no 50-component profile -- so
# this should sit far below the old ceiling. Requested generously anyway, because
# a second OOM costs a day and the headroom costs nothing on a ~1TB node.
#
# PREFIX. Deliberately NOT the old one. The previous run's a1_backbone.model.gz
# holds the scored results for all 444 models and is the evidence for this
# selection; writing over it would destroy the justification for the choice.
set -eo pipefail

REPO="${REPO:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
BACKBONE_DIR="${BACKBONE_DIR:-${REPO}/results/phylogenies/protein/a1_backbone}"
THREADS="${THREADS:-32}"
SEED="${SEED:-20260718}"
MODEL="${A1_BACKBONE_MODEL:-Q.PFAM+R10}"

# Activate conda BEFORE set -u (CONDA_BACKUP_* unbound-var trap on deactivate hooks)
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
set -u
cd "$REPO"

REF_ALN="${BACKBONE_DIR}/backbone_reference_alignment.fa"
PREFIX="${BACKBONE_DIR}/a1_backbone_qpfam"

say() { echo "[a1_backbone_fixed] $*"; }
die() { echo "[a1_backbone_fixed] ERROR: $*" >&2; exit 1; }

[ -s "$REF_ALN" ] || die "alignment missing or empty: $REF_ALN"
NSEQ=$(grep -c '^>' "$REF_ALN" 2>/dev/null || true)
[ "${NSEQ:-0}" -gt 0 ] || die "no FASTA records in $REF_ALN"
say "alignment: ${NSEQ} sequences -> $REF_ALN"

# Refuse to clobber the sweep evidence. The model cache is why we can skip
# ModelFinder at all; if the prefix ever collides, stop rather than overwrite.
for ext in .model.gz .log .treefile; do
    [ -e "${PREFIX}${ext}" ] && die "refusing to overwrite existing ${PREFIX}${ext}; \
choose a new prefix or move the previous run aside"
done
say "model FIXED to '${MODEL}' -- ModelFinder is skipped entirely"

iqtree3 -s "$REF_ALN" --prefix "$PREFIX" \
    -m "$MODEL" \
    -B 1000 -alrt 1000 --tbe \
    -seed "$SEED" -T "$THREADS"

[ -s "${PREFIX}.treefile" ] || die "IQ-TREE exited 0 but produced no treefile at ${PREFIX}.treefile"
NTIP=$(grep -o ',' "${PREFIX}.treefile" | wc -l)
say "DONE -> ${PREFIX}.treefile (~$((NTIP + 1)) tips)"
say "model actually used:"
grep -m1 "Best-fit model\|Model of substitution" "${PREFIX}.log" 2>/dev/null || true
