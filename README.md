# Berghia Chemoreceptive GPCR Discovery Pipeline

A comprehensive bioinformatics pipeline for identifying, classifying, and characterizing chemoreceptive G protein-coupled receptors (GPCRs) in *Berghia stephanieae* and related gastropod taxa.

**Author:** Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

---

## Table of Contents

1. [Overview](#overview)
2. [Pipeline Architecture](#pipeline-architecture)
3. [Prerequisites](#prerequisites)
4. [Installation](#installation)
5. [Directory Structure](#directory-structure)
6. [Configuration](#configuration)
7. [Pipeline Steps](#pipeline-steps)
8. [Python Utilities](#python-utilities)
9. [Output Files](#output-files)
10. [Running the Pipeline](#running-the-pipeline)
11. [Troubleshooting](#troubleshooting)
12. [Citation](#citation)

---

## Overview

This pipeline integrates multiple computational approaches to discover and characterize chemoreceptive GPCRs:

- **Sequence-based identification** using HMM profiles and HHblits
- **Transmembrane domain filtering** with DeepTMHMM
- **Orthology clustering** via OrthoFinder
- **Lineage-specific expansion (LSE) detection** with NCBI Taxonomy integration
- **Phylogenetic analysis** using IQ-TREE, Phyloformer, and optional MrBayes
- **Selective pressure analysis** with HyPhy aBSREL
- **Ancestral sequence reconstruction** using FastML
- **Synteny conservation** via MCScanX
- **Structural prediction** with AlphaFold
- **Integrated candidate ranking** combining all evidence types

---

## Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           BERGHIA CHEMOGPCR PIPELINE                        │
└─────────────────────────────────────────────────────────────────────────────┘

Step 01: Reference Processing
         │
         ├── Build HMM profiles (hmmbuild)
         ├── Update FASTA headers (Biopython)
         └── Generate ID mapping
         │
         ▼
Step 02: ChemoGPCR Identification
         │
         ├── DeepTMHMM (≥6 TM regions)
         ├── HHblits search
         └── HMMSEARCH (optional)
         │
         ▼
Step 03: Orthology & LSE Classification
         │
         ├── 03a: BUSCO species tree
         ├── 03b: LSE classification (NCBI Taxonomy)
         └── OrthoFinder clustering
         │
         ▼
Step 04: Phylogenetic Analysis
         │
         ├── MAFFT alignment
         ├── ClipKit/TrimAl trimming
         ├── FastTree (seed tree)
         ├── IQ-TREE (ML inference)
         ├── Phyloformer (deep learning)
         └── MrBayes (optional Bayesian)
         │
         ▼
Step 05: Selective Pressure & ASR
         │
         ├── pal2nal (codon alignment)
         ├── HyPhy aBSREL (dN/dS)
         └── FastML (ancestral reconstruction)
         │
         ▼
Step 06: Synteny & Mapping
         │
         ├── minimap2 (transcript mapping)
         ├── BLASTP (all-vs-all)
         └── MCScanX (synteny blocks)
         │
         ▼
Step 07: Candidate Ranking
         │
         ├── Phylogenetic proximity score
         ├── dN/dS score (log-transformed)
         ├── Synteny conservation score
         ├── Expression score
         └── LSE depth score
         │
         ▼
Step 08: Structural Analysis
         │
         ├── AlphaFold predictions
         ├── GPCRdb reference structures
         └── FoldTree (structural phylogeny)
         │
         ▼
Step 09: Report Generation
         │
         └── LaTeX PDF report
```

---

## Prerequisites

### Required Software

| Software | Version | Purpose |
|----------|---------|---------|
| Python | ≥3.8 | Pipeline scripts |
| HMMER | ≥3.3 | HMM building and searching |
| HH-suite | ≥3.3 | HHblits/HHsearch |
| MAFFT | ≥7.0 | Multiple sequence alignment |
| IQ-TREE | ≥2.0 | Maximum likelihood phylogeny |
| FastTree | ≥2.1 | Approximate ML trees |
| TrimAl | ≥1.4 | Alignment trimming |
| ClipKit | ≥1.3 | Smart alignment trimming |
| OrthoFinder | ≥2.5 | Orthology inference |
| HyPhy | ≥2.5 | Selective pressure analysis |
| FastML | ≥3.1 | Ancestral sequence reconstruction |
| pal2nal | ≥14.0 | Codon alignment |
| minimap2 | ≥2.24 | Transcript mapping |
| samtools | ≥1.15 | BAM processing |
| BLAST+ | ≥2.12 | Sequence similarity search |
| MCScanX | ≥1.0 | Synteny detection |
| DeepTMHMM | ≥1.0 | Transmembrane prediction |
| AlphaFold | ≥2.3 | Structure prediction |
| FoldTree | ≥1.0 | Structural phylogeny |
| TMalign | ≥20190822 | Structure alignment |
| seqtk | ≥1.3 | Sequence manipulation |
| BUSCO | ≥5.4 | Ortholog assessment |
| Phyloformer | ≥1.0 | Deep learning phylogeny |
| ASTRAL | ≥5.7 | Species tree estimation |
| MrBayes | ≥3.2 | Bayesian phylogeny (optional) |
| pdfLaTeX | any | Report generation |

### Python Dependencies

```
biopython>=1.79
ete3>=3.1.2
pandas>=1.3.0
numpy>=1.21.0
matplotlib>=3.4.0
scipy>=1.7.0
requests>=2.26.0
```

---

## Installation

### Option 1: Conda Environment (Recommended)

Create a complete conda environment with all dependencies:

```bash
# Create the environment
conda create -n berghia-gpcr python=3.10 -y
conda activate berghia-gpcr

# Install bioconda and conda-forge channels
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Install core bioinformatics tools
conda install -y \
    hmmer=3.3.2 \
    hhsuite=3.3.0 \
    mafft=7.520 \
    iqtree=2.2.2 \
    fasttree=2.1.11 \
    trimal=1.4.1 \
    clipkit=1.3.0 \
    orthofinder=2.5.5 \
    hyphy=2.5.51 \
    minimap2=2.26 \
    samtools=1.17 \
    blast=2.14.0 \
    seqtk=1.4 \
    busco=5.4.7

# Install Python packages
pip install biopython ete3 pandas numpy matplotlib scipy requests

# Install MCScanX (from source)
git clone https://github.com/wyp1125/MCScanX.git
cd MCScanX && make && cd ..
export PATH=$PATH:$(pwd)/MCScanX

# Install DeepTMHMM (requires registration at DTU)
# Download from: https://services.healthtech.dtu.dk/software.php
# Follow installation instructions in the downloaded package

# Install pal2nal
wget http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz
tar -xzf pal2nal.v14.tar.gz
export PATH=$PATH:$(pwd)/pal2nal.v14

# Install FastML
conda install -c bioconda fastml

# Optional: Install AlphaFold (requires significant setup)
# See: https://github.com/deepmind/alphafold
# Or use ColabFold for easier setup:
pip install colabfold

# Optional: Install Phyloformer
pip install phyloformer

# Optional: Install MrBayes
conda install -y mrbayes=3.2.7

# Install FoldTree and TMalign
pip install foldtree
conda install -y tmalign

# Verify installations
which hmmbuild hhblits mafft iqtree2 FastTree trimal orthofinder hyphy
python -c "from Bio import SeqIO; from ete3 import Tree; print('Python packages OK')"
```

### Option 2: Manual Installation

If you prefer manual installation or need specific versions:

```bash
# Create base directory
mkdir -p ~/tools/berghia-gpcr
cd ~/tools/berghia-gpcr

# Download and install each tool following their respective documentation
# Set PATH variables accordingly in your ~/.bashrc or ~/.zshrc
```

### Verify Installation

Run the verification script:

```bash
./verify_installation.sh
```

Or manually check:

```bash
# Check all required tools
for tool in hmmbuild hmmsearch hhmake hhblits mafft iqtree2 FastTree trimal \
            clipkit orthofinder hyphy minimap2 samtools blastp seqtk busco; do
    which $tool && echo "✓ $tool found" || echo "✗ $tool NOT FOUND"
done

# Check Python dependencies
python3 -c "
import sys
packages = ['Bio', 'ete3', 'pandas', 'numpy', 'matplotlib', 'scipy', 'requests']
for pkg in packages:
    try:
        __import__(pkg)
        print(f'✓ {pkg} installed')
    except ImportError:
        print(f'✗ {pkg} NOT FOUND')
"
```

---

## Directory Structure

```
berghia-chemogpcrs/
├── config.sh                    # Global configuration
├── functions.sh                 # Shared bash functions
│
├── 01_reference_processing.sh   # Step 01: HMM building
├── 02_chemogpcrs_identification.sh # Step 02: GPCR identification
├── 03_orthology_clustering.sh   # Step 03: OrthoFinder
├── 03a_busco_species_tree.sh    # Step 03a: Species tree
├── 03b_lse_classification.sh    # Step 03b: LSE detection
├── 04_phylogenetic_analysis.sh  # Step 04: Phylogenetics
├── 05_selective_pressure_and_asr.sh # Step 05: dN/dS & ASR
├── 06_synteny_and_mapping.sh    # Step 06: Synteny
├── 07_candidate_ranking.sh      # Step 07: Ranking
├── 08_structural_analysis.sh    # Step 08: AlphaFold
├── 09_report_generation.sh      # Step 09: Report
│
├── *.py                         # Python utility scripts
│
├── references/                  # Reference sequences
│   ├── taxid1_conserved_refs.aa
│   ├── taxid1_lse_refs.aa
│   └── ...
│
├── transcriptomes/              # Input transcriptomes
│   ├── taxid_berghia_berghia.aa
│   └── ...
│
├── genomes/                     # Optional genome assemblies
│   ├── taxid_berghia_berghia.fasta
│   └── ...
│
├── custom_hmms/                 # Optional custom HMMs
│   ├── conserved.hmm
│   └── lse.hmm
│
└── results/                     # Pipeline outputs
    ├── reference_sequences/
    ├── hmms/
    ├── chemogpcrs/
    ├── orthogroups/
    ├── lse_classification/
    ├── busco/
    ├── phylogenies/
    ├── selective_pressure/
    ├── asr/
    ├── synteny/
    ├── mapping/
    ├── ranking/
    ├── structural_analysis/
    ├── report/
    └── logs/
```

---

## Configuration

Edit `config.sh` to customize the pipeline:

### Base Directories

```bash
export BASE_DIR=$(realpath "$(dirname "$0")")
export RESULTS_DIR="${BASE_DIR}/results"
export SCRIPTS_DIR="${BASE_DIR}/scripts"
export REFERENCE_DIR="${BASE_DIR}/references"
export TRANSCRIPTOME_DIR="${BASE_DIR}/transcriptomes"
export GENOME_DIR="${BASE_DIR}/genomes"
export LOGS_DIR="${RESULTS_DIR}/logs"
```

### Input Files

```bash
export TRANSCRIPTOME="${TRANSCRIPTOME_DIR}/taxid_berghia_berghia.aa"
export CONSERVED_HMM="${BASE_DIR}/custom_hmms/conserved.hmm"  # Optional
export LSE_HMM="${BASE_DIR}/custom_hmms/lse.hmm"              # Optional
export EXPRESSION_DATA="${BASE_DIR}/expression_data.csv"
export GENOME="${GENOME_DIR}/taxid_berghia_berghia.fasta"
```

### Pipeline Parameters

```bash
# Taxa to include in analysis
export TAXA=("taxid1" "taxid2" "taxid_berghia")
export BERGHIA_TAXID="taxid_berghia"

# Computational resources
export CPUS=16
export DEFAULT_TIME="24:00:00"
export DEFAULT_MEM="64G"
export GPU_ENABLED=false

# Tool-specific parameters
export HMM_EVALUE="1e-5"           # HMMSEARCH E-value threshold
export HHBLITS_EVALUE="1e-5"       # HHblits E-value threshold
export MIN_TM_REGIONS=6            # Minimum transmembrane domains
export MIN_SEQ_LENGTH=100          # Minimum sequence length
export MAX_GAP_PERCENT=50          # Maximum alignment gap percentage
export IQTREE_MODEL="TEST"         # IQ-TREE model selection
export IQTREE_BOOTSTRAP=1000       # Bootstrap replicates
export ORTHOFINDER_INFLATION=1.5   # MCL inflation parameter
export MIN_ASR_DISTANCE=0.5        # Minimum branch length for ASR
export NUM_STRUCTURAL_CANDIDATES=10 # Candidates for AlphaFold
```

### Ranking Weights

```bash
# Adjust weights to prioritize different evidence types
export PHYLO_WEIGHT=2        # Phylogenetic proximity to references
export DNDS_WEIGHT=1         # Selective pressure (|log(omega)|)
export SYNTENY_WEIGHT=3      # Synteny conservation
export EXPR_WEIGHT=1         # Expression level
export LSE_DEPTH_WEIGHT=1    # LSE clade depth
```

### NCBI Taxonomy IDs for LSE Classification

```bash
# Customize for your taxonomic groups
export LSE_AEOLID_TAXID=54397      # Aeolidida
export LSE_NUDIBRANCH_TAXID=13843  # Nudibranchia
export LSE_GASTROPOD_TAXID=644     # Gastropoda
```

### GPCRdb Parameters

```bash
export GPCRDB_SEARCH_TERMS="chemoreceptor,invertebrate"
export GPCRDB_SPECIES="Aplysia,Lottia"
export GPCRDB_FAMILIES="all"  # Or: "Class_A,Class_B1,Class_C,Adhesion,Frizzled"
```

---

## Pipeline Steps

### Step 01: Reference Processing (`01_reference_processing.sh`)

**Purpose:** Process reference FASTA files and build HMM profiles.

**Inputs:**
- Reference sequences in `${REFERENCE_DIR}/`: `{taxid}_conserved_refs.aa`, `{taxid}_lse_refs.aa`
- Optional custom HMMs in `${BASE_DIR}/custom_hmms/`

**Outputs:**
- `${RESULTS_DIR}/hmms/conserved.hmm` - Combined conserved GPCR HMM
- `${RESULTS_DIR}/hmms/lse.hmm` - Combined LSE GPCR HMM
- `${RESULTS_DIR}/reference_sequences/all_references.fa` - Combined references
- `${ID_MAP}` - ID mapping CSV (original_id, short_id, taxid)

**Process:**
1. Checks for custom HMMs; if absent, builds HMMs from reference sequences using `hmmbuild`
2. Combines all reference sequences into a single FASTA
3. Updates headers with short IDs (`ref_TAXID_N` format) using Biopython
4. Generates ID mapping file for downstream analysis

---

### Step 02: ChemoGPCR Identification (`02_chemogpcrs_identification.sh`)

**Purpose:** Identify chemoreceptive GPCRs using sequence homology and transmembrane filtering.

**Inputs:**
- Transcriptome files: `${TRANSCRIPTOME_DIR}/*.aa`
- HMMs from Step 01
- Reference sequences from Step 01

**Outputs:**
- `${RESULTS_DIR}/chemogpcrs/chemogpcrs_{taxid}.fa` - GPCR candidates per taxon
- `${RESULTS_DIR}/chemogpcrs/deeptmhmm_{taxid}/` - TM predictions

**Process:**
1. **DeepTMHMM filtering:** Predicts transmembrane topology; retains sequences with ≥6 TM regions
2. **HHblits search:** Searches candidates against reference HMM database
3. **HMMSEARCH (optional):** Searches with conserved and LSE HMMs if provided
4. **Combination:** Merges hits from all methods, removes duplicates

---

### Step 03: Orthology Clustering (`03_orthology_clustering.sh`)

**Purpose:** Cluster orthologous groups across species.

**Inputs:**
- GPCR FASTAs from Step 02
- Reference sequences from Step 01

**Outputs:**
- `${RESULTS_DIR}/orthogroups/OrthoFinder/Results*/` - OrthoFinder results
- Orthogroup FASTA files

**Process:**
1. Prepares one FASTA per species
2. Runs OrthoFinder with DIAMOND for fast homology search
3. Uses MAFFT for alignment and FastTree for gene trees

---

### Step 03a: BUSCO Species Tree (`03a_busco_species_tree.sh`)

**Purpose:** Generate species tree from conserved orthologs.

**Inputs:**
- Transcriptomes

**Outputs:**
- `${RESULTS_DIR}/busco/busco_species_tree.tre`

**Process:**
1. Runs BUSCO to identify conserved single-copy orthologs
2. Concatenates orthologs and builds species tree

---

### Step 03b: LSE Classification (`03b_lse_classification.sh`)

**Purpose:** Classify GPCRs into lineage-specific expansion categories.

**Inputs:**
- GPCR candidates from Step 02
- Reference sequences from Step 01

**Outputs:**
- `${RESULTS_DIR}/lse_classification/lse_{level}.fa` - LSEs at each taxonomic level
- Classification summary

**Process:**
1. Uses NCBI Taxonomy API to determine taxonomic lineages
2. Caches lineage queries for performance
3. Classifies expansions as Aeolid-specific, Nudibranch-specific, or Gastropod-wide

---

### Step 04: Phylogenetic Analysis (`04_phylogenetic_analysis.sh`)

**Purpose:** Construct phylogenetic trees using multiple methods.

**Inputs:**
- GPCR sequences from Step 02
- Reference sequences from Step 01
- LSE FASTAs from Step 03b

**Outputs:**
- `${RESULTS_DIR}/phylogenies/protein/*.treefile` - IQ-TREE results
- `${RESULTS_DIR}/phylogenies/visualizations/` - Tree images (PNG, PDF, SVG)

**Process:**
1. **Alignment:** MAFFT with `--auto` mode
2. **Quality check:** Validates alignment length and gap percentage
3. **Trimming:** ClipKit smart-gap or TrimAl automated1
4. **FastTree seed:** Generates approximate ML tree to avoid local optima
5. **IQ-TREE:** Maximum likelihood with ModelFinder and 1000 ultrafast bootstrap
6. **Phyloformer:** Deep learning-based tree inference
7. **MrBayes (optional):** Bayesian phylogenetic inference
8. **Visualization:** Multi-format output (PNG, PDF, SVG) with taxonomic coloring

**Key features:**
- Uses tree `format=1` for IQ-TREE output with support values
- FastTree seed strategy improves convergence on large, divergent gene families

---

### Step 05: Selective Pressure & ASR (`05_selective_pressure_and_asr.sh`)

**Purpose:** Analyze selective pressure and reconstruct ancestral sequences.

**Inputs:**
- Protein alignments from Step 04
- Nucleotide sequences (`.mrna`, `.cds`, `.fna`, `.fa`)
- Phylogenetic trees from Step 04

**Outputs:**
- `${RESULTS_DIR}/selective_pressure/absrel_results.csv` - dN/dS results
- `${RESULTS_DIR}/asr/*_asr.fa` - Ancestral sequences

**Process:**
1. **Codon alignment:** pal2nal converts protein alignment + nucleotides to codon alignment
2. **aBSREL analysis:** Detects episodic diversifying selection on branches
3. **Deep node selection:** Identifies ancestral nodes at top 10% tree depth
4. **FastML reconstruction:** Infers ancestral sequences at selected nodes

---

### Step 06: Synteny & Mapping (`06_synteny_and_mapping.sh`)

**Purpose:** Analyze synteny conservation across genomes.

**Inputs:**
- Genome assemblies: `${GENOME_DIR}/*.fasta`
- Nucleotide transcriptomes
- Protein annotations

**Outputs:**
- `${RESULTS_DIR}/synteny/synteny_ids.txt` - Synteny-conserved gene IDs
- `${RESULTS_DIR}/synteny/*.collinearity` - MCScanX collinear blocks
- `${RESULTS_DIR}/mapping/*.bam` - Transcript mappings

**Process:**
1. **Transcript mapping:** minimap2 splice-aware alignment
2. **BLAST databases:** Creates protein databases per genome
3. **All-vs-all BLASTP:** Pairwise homology detection
4. **MCScanX:** Identifies collinear synteny blocks
5. **Visualization:** Synteny dot plots and karyotype views

---

### Step 07: Candidate Ranking (`07_candidate_ranking.sh`)

**Purpose:** Integrate all evidence to rank GPCR candidates.

**Inputs:**
- Candidate IDs from Step 02
- Expression data
- Phylogenetic trees from Step 04
- dN/dS results from Step 05
- Synteny IDs from Step 06

**Outputs:**
- `${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv`
- Ranking visualization plots

**Scoring components:**

| Score | Description | Calculation |
|-------|-------------|-------------|
| `phylo_score` | Phylogenetic proximity | `1 / (min_distance_to_refs + ε)` |
| `dnds_score` | Selective pressure | `|log(omega)|` (symmetric for purifying/positive) |
| `synteny_score` | Synteny conservation | Binary (0 or 1) |
| `expression_score` | Expression level | From expression data file |
| `lse_depth_score` | LSE clade depth | Tree depth if > 75th percentile |

**Final rank score:**
```
rank_score = (phylo_score × PHYLO_WEIGHT) +
             (dnds_score × DNDS_WEIGHT) +
             (synteny_score × SYNTENY_WEIGHT) +
             (expression_score × EXPR_WEIGHT) +
             (lse_depth_score × LSE_DEPTH_WEIGHT)
```

---

### Step 08: Structural Analysis (`08_structural_analysis.sh`)

**Purpose:** Predict structures and build structural phylogenies.

**Inputs:**
- Top-ranked candidates from Step 07
- ASR sequences from Step 05

**Outputs:**
- `${RESULTS_DIR}/structural_analysis/alphafold/*.pdb` - Predicted structures
- `${RESULTS_DIR}/structural_analysis/foldtree.tre` - Structural phylogeny
- Structural vs. sequence tree comparison plots

**Process:**
1. **Diverse selection:** Selects phylogenetically diverse top candidates
2. **AlphaFold:** Predicts 3D structures (GPU recommended)
3. **GPCRdb references:** Fetches reference structures with retry logic
4. **FoldTree:** Builds structural phylogeny using TMalign distances
5. **Comparison:** Tanglegram of structural vs. sequence trees

---

### Step 09: Report Generation (`09_report_generation.sh`)

**Purpose:** Generate comprehensive PDF report.

**Inputs:**
- All results from previous steps
- Visualization images

**Outputs:**
- `${RESULTS_DIR}/pipeline_report.pdf`

**Process:**
1. Dynamically finds available visualizations
2. Generates LaTeX document with figures and methods summary
3. Compiles to PDF with pdfLaTeX

---

## Python Utilities

### Core Scripts

| Script | Purpose | Usage |
|--------|---------|-------|
| `rank_candidates.py` | Integrate evidence and rank candidates | `python3 rank_candidates.py <candidates> <expression> <phylo_dir> <selective_dir> <synteny> <output>` |
| `lse_refine.py` | Classify LSEs with NCBI Taxonomy | `python3 lse_refine.py <fasta> <tree> <output_prefix>` |
| `select_deep_nodes.py` | Select deep ancestral nodes for ASR | `python3 select_deep_nodes.py <tree> <taxid> <min_distance>` |
| `update_headers.py` | Standardize FASTA headers (Biopython) | `python3 update_headers.py <fasta> <id_map_output>` |
| `parse_hhr.py` | Parse HHblits/HHsearch results | `python3 parse_hhr.py <input.hhr> <evalue> [output.txt]` |
| `parse_absrel.py` | Parse HyPhy aBSREL JSON output | `python3 parse_absrel.py <input.json> <output.csv>` |
| `fetch_ligands.py` | Fetch GPCRdb structures with retry | `python3 fetch_ligands.py <output_dir> <ligand_csv> <terms> <species>` |

### Visualization Scripts

| Script | Purpose | Outputs |
|--------|---------|---------|
| `visualize_tree.py` | Tree visualization with coloring | PNG, PDF, SVG (basic, colored, circular) |
| `plot_asr.py` | ASR sequences on tree | PNG with sequence length annotations |
| `plot_synteny.py` | Synteny conservation plots | Dot plots, karyotype views |
| `plot_ranking.py` | Candidate ranking barplots | Score breakdown visualization |
| `plot_selective_pressure.py` | dN/dS distribution plots | Branch-wise selection visualization |
| `plot_struct_vs_seq.py` | Tree comparison tanglegrams | Structural vs. sequence phylogeny |
| `plot_phyloformer_iqtree.py` | Method comparison | Phyloformer vs. IQ-TREE trees |
| `plot_heatmap.py` | Expression heatmaps | Clustered expression patterns |
| `plot_pca.py` | PCA of sequence features | Dimensionality reduction plots |

### Analysis Scripts

| Script | Purpose |
|--------|---------|
| `select_diverse_candidates.py` | Select phylogenetically diverse candidates for structural analysis |
| `cluster_structures.py` | Cluster predicted structures by similarity |
| `test_phyloformer_models.py` | Evaluate Phyloformer model performance |
| `compute_lrt.py` | Likelihood ratio tests for model comparison |
| `prune_alignment.py` | Remove problematic sequences from alignments |

---

## Output Files

### Key Results

| File | Description |
|------|-------------|
| `results/ranking/ranked_candidates_sorted.csv` | **Primary output:** Ranked GPCR candidates |
| `results/phylogenies/protein/all_berghia_refs.treefile` | Main phylogenetic tree |
| `results/selective_pressure/absrel_results.csv` | dN/dS analysis results |
| `results/synteny/synteny_ids.txt` | Synteny-conserved genes |
| `results/asr/*_asr.fa` | Reconstructed ancestral sequences |
| `results/structural_analysis/alphafold/*.pdb` | Predicted 3D structures |
| `results/pipeline_report.pdf` | Comprehensive analysis report |

### Intermediate Files

| Directory | Contents |
|-----------|----------|
| `results/reference_sequences/` | Processed reference sequences and ID mappings |
| `results/hmms/` | HMM profiles |
| `results/chemogpcrs/` | Identified GPCR sequences and DeepTMHMM results |
| `results/orthogroups/` | OrthoFinder clustering results |
| `results/lse_classification/` | LSE categories |
| `results/phylogenies/` | Trees, alignments, visualizations |
| `results/selective_pressure/` | aBSREL results, codon alignments |
| `results/mapping/` | BAM files from transcript mapping |
| `results/logs/` | Execution logs and error files |

---

## Running the Pipeline

### Local Execution

Run steps sequentially:

```bash
# Source configuration
source config.sh
source functions.sh

# Run each step
./01_reference_processing.sh
./02_chemogpcrs_identification.sh
./03_orthology_clustering.sh
./03a_busco_species_tree.sh
./03b_lse_classification.sh
./04_phylogenetic_analysis.sh
./05_selective_pressure_and_asr.sh
./06_synteny_and_mapping.sh
./07_candidate_ranking.sh
./08_structural_analysis.sh
./09_report_generation.sh
```

### SLURM Cluster Execution

Submit jobs with dependencies:

```bash
# Submit step 01
jid1=$(sbatch --parsable 01_reference_processing.sh)

# Submit step 02 (depends on 01)
jid2=$(sbatch --parsable --dependency=afterok:$jid1 02_chemogpcrs_identification.sh)

# Submit step 03 (depends on 02)
jid3=$(sbatch --parsable --dependency=afterok:$jid2 03_orthology_clustering.sh)

# Continue with remaining steps...
```

### Partial Execution

Resume from a specific step using completion flags:

```bash
# Check completed steps
ls results/step_completed_*.txt

# The pipeline automatically skips completed steps
./04_phylogenetic_analysis.sh  # Will skip if already done
```

---

## Troubleshooting

### Common Issues

**1. DeepTMHMM not found**
```
Error: deeptmhmm: command not found
```
Solution: Download from DTU and add to PATH, or use Docker:
```bash
docker run --rm -v $(pwd):/data dtubioinformatics/deeptmhmm -f /data/input.fa
```

**2. Tree format errors**
```
Error: Cannot parse tree file
```
Solution: Ensure IQ-TREE output with `format=1`:
```python
t = Tree(tree_file, format=1)  # For IQ-TREE with support values
```

**3. NCBI Taxonomy API rate limiting**
```
Warning: HTTPError 429 Too Many Requests
```
Solution: The pipeline uses lineage caching. For large datasets, consider pre-downloading taxonomy database.

**4. AlphaFold memory issues**
```
Error: CUDA out of memory
```
Solution: Reduce batch size or use CPU mode:
```bash
export GPU_ENABLED=false
```

**5. MCScanX no collinear blocks**
```
Warning: No synteny blocks found
```
Solution: Ensure GFF files have correct format and genomes are from related species.

### Validation Checks

Run the validation suite:

```bash
# Check Python dependencies
python3 -c "
from Bio import SeqIO
from ete3 import Tree
import pandas as pd
import numpy as np
print('All Python dependencies OK')
"

# Verify tree loading
python3 -c "
from ete3 import Tree
t = Tree('(A:0.1,B:0.2):0.3;')
print(f'Tree has {len(t)} leaves')
"

# Test Biopython FASTA parsing
python3 -c "
from Bio import SeqIO
import io
fasta = '>test\nACGT'
for r in SeqIO.parse(io.StringIO(fasta), 'fasta'):
    print(f'Parsed: {r.id}, length {len(r.seq)}')
"
```

### Log Files

Check logs for detailed error messages:

```bash
# View main pipeline log
tail -f results/logs/pipeline.log

# Check specific step errors
cat results/logs/hmmsearch_conserved_berghia.err

# Find all errors
grep -r "Error" results/logs/
```

---

## Citation

If you use this pipeline, please cite:

```
Perez-Moreno, J.L. (2024). Berghia Chemoreceptive GPCR Discovery Pipeline.
Katz Lab, University of Massachusetts, Amherst.
```

And the underlying tools:
- HMMER: Eddy, S.R. (2011) PLoS Comput Biol 7:e1002195
- HH-suite: Steinegger et al. (2019) BMC Bioinformatics 20:473
- IQ-TREE: Minh et al. (2020) Mol Biol Evol 37:1530-1534
- OrthoFinder: Emms & Kelly (2019) Genome Biol 20:238
- HyPhy: Kosakovsky Pond et al. (2020) Mol Biol Evol 37:295-299
- AlphaFold: Jumper et al. (2021) Nature 596:583-589

---

## License

This pipeline is provided for academic research purposes. See individual tool licenses for commercial use restrictions.
