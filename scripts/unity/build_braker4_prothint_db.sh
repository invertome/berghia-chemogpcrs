#!/bin/bash
# build_braker4_prothint_db.sh — build the combined ProtHint protein-evidence DB
# for BRAKER4: OrthoDB v12 Metazoa (40.3M proteins) + MolluscaGenes mollusca_aa
# (17.1M proteins / 200+ mollusc species incl. cephalopods, EvidentialGene).
#
# Rationale (bead: BRAKER4 evidence standardization): the current evidence is the
# 40.3M-protein OrthoDB Metazoa set, yet Nautilus got only ~98-306 ProtHint hits
# from it -> GeneMark-EP parse_ET.pl div-by-zero. The problem is not DB size but
# the lack of in-clade close homologs; mollusca_aa supplies 200+ mollusc species'
# proteins. Combined DB (decided with user): ~57M proteins / ~25 GB.
#
# Headers: mollusca_aa EvidentialGene deflines all start with 'type=protein;'
# (non-unique first token). Rewrite to a guaranteed-unique '>mgN' ID while
# keeping the original defline (organism=/oid=) as the description for provenance.
#
# Submit:  sbatch scripts/unity/build_braker4_prothint_db.sh
#
#SBATCH --job-name=build_prothint_db
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/build_prothint_db-%j.out
#SBATCH --error=logs/build_prothint_db-%j.err
#SBATCH --account=pi_pkatz_umass_edu

set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate ncbi_processing_env   # blastdbcmd 2.16.0+ reads the BLAST v5 mollusca_aa DB
set -u

SC=/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
MOLLUSCA_DB=/work/pi_pkatz_umass_edu/jpm/blast/db/20250320/mollusca_aa
ORTHODB="${SC}/species_tree_data/orthodb/odb12_metazoa.fa"
OUTDIR="${SC}/species_tree_data/orthodb"
MOLLUSCA_FA="${OUTDIR}/mollusca_aa.from_blastdb.fa"
COMBINED="${OUTDIR}/odb12_metazoa_plus_mollusca_aa.fa"

[ -s "$ORTHODB" ] || { echo "[db] ERROR: missing OrthoDB $ORTHODB"; exit 1; }
df -h "$OUTDIR" | tail -1

# 1. Extract mollusca_aa -> FASTA with unique '>mgN <original defline>' headers.
echo "[db] $(date) extracting mollusca_aa (17.1M proteins) ..."
blastdbcmd -db "$MOLLUSCA_DB" -dbtype prot -entry all -outfmt "%f" \
  | awk '/^>/{print ">mg"(++n)" "substr($0,2); next}{print}' > "$MOLLUSCA_FA"
echo "[db] mollusca_aa extracted: $(grep -c '^>' "$MOLLUSCA_FA") proteins, $(du -h "$MOLLUSCA_FA" | cut -f1)"

# 2. Concatenate OrthoDB Metazoa + mollusca_aa -> combined ProtHint evidence DB.
echo "[db] $(date) concatenating combined DB ..."
cat "$ORTHODB" "$MOLLUSCA_FA" > "$COMBINED"
echo "[db] combined: $(grep -c '^>' "$COMBINED") proteins, $(du -h "$COMBINED" | cut -f1)"
echo "[db] $(date) DONE -> $COMBINED"
