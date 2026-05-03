# HPC Installation Guide

This document covers tools that cannot be installed via `environment.yml` (compiled binaries, Java jars, Perl scripts) and must be set up manually on the HPC node before running the corresponding pipeline stages.

For local-only development the conda env (`environment.yml`) is sufficient — the HPC tools below are only needed for full end-to-end runs at scale.

---

## Tools installed by `environment.yml` (no manual step needed)

These are pip/conda-installable and pinned in `environment.yml`. Run `conda env create -f environment.yml` to get them.

| Tool | Purpose | Bead |
|---|---|---|
| `tmbed` | Primary TM-helix predictor (Bernhofer & Rost 2022, +9pp recall over DeepTMHMM) | -6nh |
| `gffutils>=0.13` | GFF parsing for tandem-cluster + JCVI pipelines | -ar8 |
| `dendropy>=5.0` | Yule birth-death simulation for null-calibrated depth thresholds | -m6k |
| `statsmodels>=0.13` | BH-FDR correction (replaces buggy hand-rolled implementation) | -wux |
| `jcvi>=1.4.15` | JCVI MCscan for synteny (Tang 2024) — replaces MCScanX | -e59 |
| `treeshrink>=1.3.9` | Rogue-taxon / outlier-long-branch removal | -iof |

System-level prereqs for JCVI MCscan: `last` aligner and `texlive-latex-base` (for figures).

---

## Tools that require manual install

### MACSE v2 — codon-aware MSA (bead -i61)

Required by `scripts/run_macse.sh` (invoked from `05_selective_pressure_and_asr.sh`).

```bash
# Download to project tools/ directory
mkdir -p $BASE_DIR/tools
cd $BASE_DIR/tools
wget -q https://www.agap-ge2pop.org/wp-content/uploads/macse/releases/macse_v2.07.jar
sha256sum macse_v2.07.jar  # verify against the value posted on the upstream page

# Set in your shell or in config.sh:
export MACSE_JAR="$BASE_DIR/tools/macse_v2.07.jar"
export JAVA_BIN="java"  # or path to your Java 11+ binary
```

Pipeline integration: stage 05 runs MACSE before pal2nal when `RUN_MACSE=1` (default). Falls back to pal2nal when `MACSE_JAR` is unset or the file is missing.

### HmmCleaner.pl — per-sequence segment cleaning (bead -i61)

Required by `scripts/run_hmmcleaner.sh` (invoked from `04_phylogenetic_analysis.sh`).

```bash
# Install via CPAN (requires Perl + cpanminus)
cpanm Bio::MUST::Apps::HmmCleaner

# Verify
which HmmCleaner.pl
HmmCleaner.pl --version
```

Pipeline integration: stage 04 runs HmmCleaner after MAFFT and before ClipKit when `RUN_HMMCLEANER=1` (default). Failures are non-fatal — continue with un-cleaned alignment.

### GeneRax — phylogenetic + reconciliation (bead -30g)

Required by `scripts/hpc/run_generax.sh`. Compile from source:

```bash
git clone --recursive https://github.com/BenoitMorel/GeneRax
cd GeneRax
./install.sh
# Adds generax to ./build/bin
export GENERAX="$PWD/build/bin/generax"
```

Pipeline integration: not auto-invoked. Submit `scripts/hpc/run_generax.sh <families.txt>` manually after orthogroup trees are built.

### GPCRtm matrix file — substitution-model retest (bead -5b0)

Optional. Used by `scripts/hpc/test_models_gpcrtm_c20.sh`.

```bash
# Source: Rios et al. 2015 BMC Bioinf 16:206 supplementary, or contact authors.
# Save the matrix file as:
mv GPCRtm.txt $BASE_DIR/references/GPCRtm.txt
```

Without this file the model retest still runs — just without GPCRtm in the candidate set (LG+C20+R10, LG+C60+R10, LG4X+R10 still tested).

### Foldseek + GPCRdb structures — structural augmentation (bead -s6v, P3)

Optional. Used by `scripts/hpc/foldseek_against_gpcrdb.sh`.

```bash
# Foldseek
wget https://github.com/steineggerlab/foldseek/releases/download/9-427df8a/foldseek-linux-aarch64.tar.gz
tar xzf foldseek-linux-aarch64.tar.gz
export PATH="$PWD/foldseek/bin:$PATH"

# GPCRdb 2025 AlphaFold-Multistate models for 814 human ORs:
# https://gpcrdb.org/structure_models — bulk download required
foldseek createdb gpcrdb_pdbs/ $BASE_DIR/references/foldseek/gpcrdb_2025
```

---

## SLURM submission examples

### Substitution-model retest

```bash
sbatch scripts/hpc/test_models_gpcrtm_c20.sh \
    $RESULTS_DIR/phylogenies/protein/all_berghia_refs_trimmed.fa
```

Expected: 24–72 hours on 16 CPUs / 64 GB.

### Selection stack (GARD → BUSTED-S → BUSTED-MH → aBSREL → MEME)

```bash
# Per orthogroup (run as a SLURM array)
sbatch --array=0-N%50 scripts/hpc/run_selection_stack.sh \
    "$codon_align" "$tree" "$og_base"
```

Expected: ~hours per OG depending on size.

### GeneRax reconciliation

```bash
# Build the families.txt config first (one entry per OG)
sbatch scripts/hpc/run_generax.sh families.txt
```

### Foldseek vs GPCRdb (P3)

```bash
sbatch scripts/hpc/foldseek_against_gpcrdb.sh \
    $RESULTS_DIR/structural_analysis/alphafold/
```

---

## Re-validating with the Aplysia genome (bead -bdu)

Retrospective recall/precision evaluation on a labeled reference set:

```bash
# 1. Download Aplysia californica genome
datasets download genome accession GCF_000002075.2 \
    --include genome,gff3,protein,cds \
    --filename /tmp/aplysia.zip
unzip /tmp/aplysia.zip -d $GENOME_DIR/aplysia/
# 2. Edit references/cummins2009_aplysia_chemoreceptors.csv to add
#    real RefSeq protein IDs (the seed CSV has placeholder gene names).
# 3. Run the pipeline on Aplysia (override BERGHIA_FILE_PREFIX in env)
# 4. Compute recall@N
python3 scripts/validate_against_aplysia.py \
    --ranked-csv results/aplysia_run/ranking/ranked_candidates_sorted.csv \
    --labeled-csv references/cummins2009_aplysia_chemoreceptors.csv \
    --out results/validation/aplysia_recall.tsv
```
