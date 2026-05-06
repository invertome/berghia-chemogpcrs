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

### FAMSA 2 — large-family MSA (bead -align)

Required by `scripts/run_aligner.sh` when N ≥ 1000 sequences. Auto-falls-back
to MAFFT --auto if missing. Available via bioconda:

```bash
conda install -c bioconda famsa>=2.2
# or build from source:
git clone https://github.com/refresh-bio/FAMSA && cd FAMSA && make -j
export FAMSA="$PWD/famsa"
```

Verify: `famsa --help | head -3` should show "FAMSA (Fast and Accurate
Multiple Sequence Alignment) v2.x.x".

Selection rationale (per `run_aligner.sh`):
- N < 200 → MAFFT L-INS-i (gold-standard accuracy small)
- 200 ≤ N < 1000 → MAFFT --auto (MAFFT picks mode)
- N ≥ 1000 → FAMSA 2 (5-10× faster than MAFFT --retree 2; HomFam SP ≈ MAFFT --auto +3pp)

Override via `ALIGNER_BACKEND` env var or `--force-backend` flag.

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

### Alignment-cleanup stack: PREQUAL + CLOAK + TAPER (bead -i61, May 2026 v2)

Replaces HmmCleaner. The legacy HmmCleaner path was abandoned because its
Perl XS dependency tree (`Class::XSAccessor`, `Moose`, the `Bio::MUST::*`
distribution) cannot be built under the conda env's Perl ABI on Unity HPC
without sudo / system headers — every CPAN attempt failed inside
`Bio::MUST::Drivers` tests, no bioconda recipe exists, and no published
container image exists either (see `scripts/unity/README_hmmcleaner_install.md`
for the full post-mortem). HmmCleaner is now deprecated; the gating
variable `RUN_HMMCLEANER` defaults to `0` and the binary is no longer
required.

The replacement is a four-tool stack chained inside
`run_alignment_filter_stack` in `functions.sh` and invoked from
`04_phylogenetic_analysis.sh` at three sites (global concat, per-LSE-level,
per-OG):

1. **PREQUAL** — pre-alignment residue-level mask
   (Whelan, Irisarri & Burki 2018, *Bioinformatics* **34**:3929,
   <https://doi.org/10.1093/bioinformatics/bty448>;
   source <https://github.com/simonwhelan/prequal>). Bioconda binary.
2. **Alignment ensemble** (`scripts/run_alignment_ensemble.sh`) — six
   alignments per dataset feeding CLOAK: canonical MAFFT (regime-based
   via `run_aligner.sh`) + L-INS-i + G-INS-i + E-INS-i + FFT-NS-1 + FAMSA
   (K = 6).
3. **CLOAK** — alignment-uncertainty consensus mask across the ensemble
   (Chatur, Wheeler lab, *bioRxiv* 691663, December 2025,
   <https://doi.org/10.1101/2025.12.06.691663>;
   source <https://github.com/phylowheeler/CLOAK>). Single-file Python
   script; no pip deps beyond stdlib.
4. **TAPER** — post-alignment per-sequence residue-outlier mask
   (Zhang, Zhang, Stamatakis & Mirarab 2021, *Methods Ecol Evol* **13**:91,
   <https://doi.org/10.1111/2041-210X.13696>;
   source <https://github.com/chaoszhang/TAPER>). Julia script; requires
   `julia` (installed from conda-forge by the installer).

Then ClipKit handles column-level trimming downstream.

**One-shot install (Unity SLURM):**

```bash
sbatch scripts/unity/install_alignment_filters.sh
```

Submitted from the repo root, the job:

- `mamba install -c bioconda prequal` into the `berghia-gpcr` env;
- `git clone https://github.com/chaoszhang/TAPER.git` into
  `${WORKDIR}/tools/TAPER` and `mamba install -c conda-forge julia` if
  needed;
- `git clone https://github.com/phylowheeler/CLOAK.git` into
  `${WORKDIR}/tools/CLOAK` (no pip deps).

End-of-job summary prints OK/FAILED per tool plus the `export` lines to
paste into `config.local.sh`. Partial installs are tolerated — a failure
in one tool does not abort the others. End-to-end smoke test
(`scripts/unity/smoke_test_filter_stack.sh`) verified on Unity in 32 s
(job 56830329) and 29 s (job 56831279).

**Gating env vars (config.sh defaults):**

| Var | Default | Effect |
|---|---|---|
| `RUN_PREQUAL`     | `1` | residue mask before MAFFT |
| `RUN_CLOAK`       | `1` | consensus mask across the 6-alignment ensemble |
| `RUN_TAPER`       | `1` | per-sequence residue-outlier mask after CLOAK |
| `RUN_HMMCLEANER`  | `0` | DEPRECATED; if `=1` *and* `HmmCleaner.pl` is on `PATH`, runs after CLOAK as a 4th legacy filter pass (back-compat only) |

Each gate is independent so the stages can be ablated for sensitivity
analysis. Tool paths (`PREQUAL`, `CLOAK`, `TAPER`, `JULIA`) are exported
in `config.sh` and overridable via `config.local.sh`.

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
