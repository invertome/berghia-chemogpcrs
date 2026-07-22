#!/bin/bash
#SBATCH --job-name=taar_prep
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --account=pi_pkatz_umass_edu
#
# TAAR embedding job PREP (bead h0y0): build the combined 2395-seq FASTA +
# membership.tsv from the four canonical source sets. The two GPU embed jobs
# (esm2_3b, protrek) depend on this via afterok. Pure-stdlib file assembly, but
# it is compute per the login-node guard, so it runs here.
set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"; conda activate berghia-gpcr; set -u

W=/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
cd "$W"
T=results/ranking/taar_embedding
mkdir -p "$T" logs

python3 scripts/analysis/build_taar_embedding_fasta.py \
    --root . --out-fasta "$T/combined.fa" --out-membership "$T/membership.tsv"

N=$(grep -c '^>' "$T/combined.fa")
[ "$N" -eq 2395 ] || { echo "FATAL: combined.fa has $N seqs, expected 2395" >&2; exit 2; }
U=$(grep '^>' "$T/combined.fa" | sort -u | wc -l)
[ "$U" -eq 2395 ] || { echo "FATAL: $U unique ids, expected 2395 (collision)" >&2; exit 2; }
echo "[prep] combined.fa = $N seqs, $U unique. membership:"
tail -n +2 "$T/membership.tsv" | cut -f2 | sort | uniq -c
echo "[prep] DONE $(date -Is)"
