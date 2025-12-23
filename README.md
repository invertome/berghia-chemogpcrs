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
- **Redundancy reduction** via CD-HIT clustering
- **GPCR classification** with InterProScan domain annotation
- **Orthology clustering** via OrthoFinder
- **Lineage-specific expansion (LSE) detection** with CAFE5 birth-death models
- **Gene tree reconciliation** using NOTUNG
- **Phylogenetic analysis** using IQ-TREE, Phyloformer, and optional MrBayes
- **Selective pressure analysis** with HyPhy aBSREL
- **Ancestral sequence reconstruction** using FastML
- **Divergence dating** with molecular clock methods
- **Convergent evolution detection** across lineages
- **Synteny conservation** via MCScanX
- **Structural prediction** with AlphaFold
- **Binding site prediction** from 3D structures
- **Conservation mapping** with per-residue analysis
- **Integrated candidate ranking** with sensitivity analysis

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
Step 02a: Sequence Clustering [NEW]
         │
         └── CD-HIT (98% identity)
         │
         ▼
Step 02b: GPCR Classification [NEW]
         │
         ├── InterProScan domain annotation
         └── Sub-family classification
         │
         ▼
Step 03: Orthology & LSE Classification
         │
         ├── 03a: BUSCO species tree
         ├── 03b: LSE classification (NCBI Taxonomy)
         ├── 03c: CAFE5 analysis (birth-death models) [NEW]
         ├── 03d: NOTUNG reconciliation [NEW]
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
         ├── Weighted phylogenetic proximity
         ├── Separate purifying/positive selection scores
         ├── Quantitative synteny scoring
         ├── Expression analysis (Salmon TPM) [NEW]
         ├── Sensitivity analysis [NEW]
         ├── Cross-validation [NEW]
         └── Confidence tier assignment
         │
         ▼
Step 08: Structural Analysis
         │
         ├── AlphaFold predictions
         ├── Binding site prediction [NEW]
         ├── Conservation mapping [NEW]
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
| R | ≥4.0 | Ultrametric tree generation |
| HMMER | ≥3.3 | HMM building and searching |
| HH-suite | ≥3.3 | HHblits/HHsearch |
| MAFFT | ≥7.0 | Multiple sequence alignment |
| IQ-TREE | ≥2.0 | Maximum likelihood phylogeny |
| FastTree | ≥2.1 | Approximate ML trees |
| TrimAl | ≥1.4 | Alignment trimming |
| ClipKit | ≥1.3 | Smart alignment trimming |
| CD-HIT | ≥4.8 | Sequence clustering |
| OrthoFinder | ≥2.5 | Orthology inference |
| CAFE5 | ≥5.0 | Gene family evolution |
| HyPhy | ≥2.5 | Selective pressure analysis |
| FastML | ≥3.1 | Ancestral sequence reconstruction |
| pal2nal | ≥14.0 | Codon alignment |
| minimap2 | ≥2.24 | Transcript mapping |
| samtools | ≥1.15 | BAM processing |
| BLAST+ | ≥2.12 | Sequence similarity search |
| MCScanX | ≥1.0 | Synteny detection |
| DeepTMHMM | ≥1.0 | Transmembrane prediction |
| InterProScan | ≥5.0 | Domain annotation (optional) |
| AlphaFold | ≥2.3 | Structure prediction |
| FoldTree | ≥1.0 | Structural phylogeny |
| TMalign | ≥20190822 | Structure alignment |
| seqtk | ≥1.3 | Sequence manipulation |
| BUSCO | ≥5.4 | Ortholog assessment |
| Phyloformer | ≥1.0 | Deep learning phylogeny |
| ASTRAL | ≥5.7 | Species tree estimation |
| MrBayes | ≥3.2 | Bayesian phylogeny (optional) |
| NOTUNG | ≥2.9 | Gene tree reconciliation (optional) |
| pdfLaTeX | any | Report generation |

### Python Dependencies

```
biopython>=1.79
ete3>=3.1.2
pandas>=1.3.0
numpy>=1.21.0
matplotlib>=3.4.0
seaborn>=0.11.0
scipy>=1.7.0
requests>=2.26.0
```

### R Dependencies

```r
ape>=5.7
```

---

## Installation

### Option 1: Automated Setup (Recommended)

Use the provided setup script to create a complete conda environment:

```bash
# Run the setup script (full installation)
./setup_conda_env.sh

# Activate the environment
conda activate berghia-gpcr

# Verify installation
./setup_conda_env.sh --verify
```

#### Setup Script Options

| Option | Description |
|--------|-------------|
| `./setup_conda_env.sh` | Full installation with all tools (default) |
| `./setup_conda_env.sh --minimal` | Core tools only (faster, smaller) |
| `./setup_conda_env.sh --full` | Complete installation with optional tools |
| `./setup_conda_env.sh --verify` | Check all tools are correctly installed |
| `./setup_conda_env.sh --external` | Instructions for external tools (DeepTMHMM, NOTUNG, etc.) |
| `./setup_conda_env.sh --help` | Show all options |

#### What Gets Installed

**Core Tools (always installed):**
- HMMER, HH-suite, MAFFT, IQ-TREE, FastTree, TrimAl, ClipKit
- CD-HIT, OrthoFinder, HyPhy, minimap2, samtools, BLAST+
- BUSCO, seqtk, R with ape package

**Python Packages:**
- biopython, ete3, pandas, numpy, matplotlib, seaborn, scipy, requests

**Optional Tools (--full only):**
- CAFE5, MrBayes

**External Tools (manual installation):**
The script provides instructions for tools that require manual setup:
- DeepTMHMM (requires DTU registration)
- NOTUNG (Java-based, download JAR)
- MCScanX (compiled from source automatically)
- pal2nal (downloaded automatically to `tools/`)
- InterProScan (large download)
- AlphaFold (complex setup)

#### Verification Output

The `--verify` option checks each tool and shows status:

```
✓ HMMER (hmmbuild)
✓ HH-suite (hhblits)
✓ MAFFT (mafft)
✓ IQ-TREE (iqtree2)
...
✓ biopython
✓ ete3
...
⚠ InterProScan (interproscan.sh) - not installed [optional]
```

### Option 2: Manual Conda Environment

```bash
# Create environment from YAML file
conda env create -f environment.yml

# Activate
conda activate berghia-gpcr
```

### Option 3: Step-by-Step Installation

```bash
# Create the environment
conda create -n berghia-gpcr python=3.10 -y
conda activate berghia-gpcr

# Add channels
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
    cd-hit=4.8.1 \
    orthofinder=2.5.5 \
    hyphy=2.5.51 \
    minimap2=2.26 \
    samtools=1.17 \
    blast=2.14.0 \
    seqtk=1.4 \
    busco=5.4.7 \
    r-ape=5.7

# Install Python packages
pip install biopython ete3 pandas numpy matplotlib seaborn scipy requests python-dotenv

# Install CAFE5
conda install -c bioconda cafe=5.0.0
```

### Environment Variables (.env file)

Create a `.env` file for API credentials (gitignored by default):

```bash
# Create .env file with your NCBI API key
cat > .env << 'EOF'
# NCBI Entrez API Key
# Get your key at: https://www.ncbi.nlm.nih.gov/account/settings/
# Increases rate limit from 3 to 10 requests/second
NCBI_API_KEY=your_api_key_here

# Email for NCBI Entrez (required)
NCBI_EMAIL=your.email@example.com
EOF
```

The `.env` file is automatically loaded by Python scripts and used for:
- CDS retrieval from NCBI (`fetch_reference_cds.py`)
- Taxonomy lookups (`rename_references_with_taxid.py`)
- Any Biopython Entrez operations

**Security Note:** Never commit your `.env` file to version control.

### API-Independent Mode (Recommended)

Pre-download required databases for offline operation:

```bash
# Download NCBI Taxonomy and GPCRdb data locally
python3 setup_databases.py --db-dir databases

# Verify setup
python3 setup_databases.py --db-dir databases --verify-only
```

---

## Directory Structure

```
berghia-chemogpcrs/
├── config.sh                    # Global configuration
├── functions.sh                 # Shared bash functions (with checkpointing)
├── setup_conda_env.sh           # Conda environment setup script
├── validate_config.sh           # Configuration validator [NEW]
├── environment.yml              # Conda environment specification
│
├── 01_reference_processing.sh   # Step 01: HMM building
├── 02_chemogpcrs_identification.sh # Step 02: GPCR identification
├── 02a_cluster_sequences.sh     # Step 02a: CD-HIT clustering [NEW]
├── 02b_classify_gpcrs.sh        # Step 02b: InterPro classification [NEW]
├── 03_orthology_clustering.sh   # Step 03: OrthoFinder
├── 03a_busco_species_tree.sh    # Step 03a: Species tree
├── 03b_lse_classification.sh    # Step 03b: LSE detection
├── 03c_cafe_analysis.sh         # Step 03c: CAFE5 analysis [NEW]
├── 03d_notung_reconciliation.sh # Step 03d: Gene tree reconciliation [NEW]
├── 04_phylogenetic_analysis.sh  # Step 04: Phylogenetics
├── 05_selective_pressure_and_asr.sh # Step 05: dN/dS & ASR
├── 06_synteny_and_mapping.sh    # Step 06: Synteny
├── 07_candidate_ranking.sh      # Step 07: Ranking
├── 08_structural_analysis.sh    # Step 08: AlphaFold
├── 09_report_generation.sh      # Step 09: Report
│
├── scripts/                     # Python/R utility scripts
│   ├── rename_references_with_taxid.py  # Add taxid prefix to reference files [NEW]
│   ├── fetch_reference_cds.py   # Fetch CDS for reference proteins [NEW]
│   ├── get_metadata.py          # Centralized metadata lookups [NEW]
│   ├── expression_analysis.py   # Salmon TPM parsing [NEW]
│   ├── divergence_dating.py     # Molecular clock dating [NEW]
│   ├── convergent_evolution.py  # Convergent substitution detection [NEW]
│   ├── binding_site_prediction.py # Binding pocket prediction [NEW]
│   ├── conservation_mapping.py  # Per-residue conservation [NEW]
│   └── generate_ultrametric.R   # Ultrametric tree conversion [NEW]
│
├── *.py                         # Core Python scripts
│
├── references/                  # Reference sequences
├── transcriptomes/              # Input transcriptomes
├── genomes/                     # Optional genome assemblies
├── databases/                   # Local databases (API-independent mode)
│
└── results/                     # Pipeline outputs
    ├── clustering/              # CD-HIT results [NEW]
    ├── classification/          # GPCR classification [NEW]
    ├── cafe/                    # CAFE5 results [NEW]
    ├── notung/                  # Gene tree reconciliation [NEW]
    ├── phylogenies/
    ├── selective_pressure/
    ├── ranking/
    │   ├── sensitivity_analysis.csv   # Ranking stability [NEW]
    │   └── weight_importance.json     # Weight analysis [NEW]
    ├── structural_analysis/
    └── logs/
```

---

## Input Data Setup

Before running the pipeline, you must prepare your input files in the correct format and directory structure. This section provides detailed instructions for each required and optional input type.

### Directory Preparation

Create the required input directories:

```bash
mkdir -p references transcriptomes genomes expression_data custom_hmms
```

### 1. Reference Sequences (Required)

Reference sequences are known GPCRs used to build HMM profiles and identify candidates in your transcriptomes.

#### File Location
```
references/
├── taxid1_conserved_refs.aa     # Conserved GPCR references for taxon 1
├── taxid1_lse_refs.aa           # LSE-specific references for taxon 1
├── taxid2_conserved_refs.aa     # Conserved references for taxon 2
├── taxid2_lse_refs.aa           # LSE-specific references for taxon 2
└── ...
```

#### Naming Convention
Files must follow the pattern: `{TAXID}_conserved_refs.aa` and `{TAXID}_lse_refs.aa`

- `TAXID` must match entries in the `TAXA` array in `config.sh`
- Use NCBI taxonomy IDs (e.g., `9606` for human) or custom identifiers
- The `.aa` extension indicates amino acid (protein) sequences

#### File Format
Standard FASTA format with informative headers:

```fasta
>sp|P30872|OR1A1_HUMAN Olfactory receptor 1A1 OS=Homo sapiens OX=9606 GN=OR1A1
MEFRNLTVITDFLMSLPFPFQTNLSCVFILVFVATHTVIALVTLAGNLLIILAVTSDPRLH
TPMYFFLSNLSFVDLCFSSVTVPKMLANFIHSGRPISYAGCITQMFFCVLFGVDSGMLFAF
...
>sp|Q9UBG0|OR2W1_HUMAN Olfactory receptor 2W1 OS=Homo sapiens OX=9606 GN=OR2W1
MSGENGTVSLEFFLLGLSSQPQQQQLLFVFFLSMYLATVLGNLLIMLSVLSDSHLHTPMYF
...
```

**Recommended header fields (for automatic metadata extraction):**
- `OS=` - Organism/species name
- `OX=` - NCBI taxonomy ID
- `GN=` - Gene name

The pipeline will parse these fields to populate the metadata CSV.

#### Sourcing Reference Sequences

| Source | Description | Recommended For |
|--------|-------------|-----------------|
| [GPCRdb](https://gpcrdb.org) | Curated GPCR database | Class A GPCRs, mammalian |
| [UniProt](https://uniprot.org) | Reviewed protein sequences | Broad taxonomic coverage |
| NCBI Protein | Comprehensive, includes predictions | Invertebrate GPCRs |
| Published studies | Species-specific validated sequences | Closely related species |

**Example download from UniProt:**
```bash
# Download all reviewed olfactory receptors
wget "https://rest.uniprot.org/uniprotkb/stream?query=family:olfactory+AND+reviewed:true&format=fasta" \
    -O references/vertebrate_conserved_refs.aa
```

#### Conserved vs. LSE References

- **Conserved (`_conserved_refs.aa`)**: Broadly conserved GPCRs found across taxa (e.g., rhodopsin, adrenergic receptors). Used to identify canonical GPCR domains.

- **LSE (`_lse_refs.aa`)**: Lineage-specific expansions unique to certain clades (e.g., mollusk-specific chemoreceptors). Used to identify novel, clade-specific expansions.

#### Alternative: nath_et_al Directory Structure

The pipeline also supports an alternative reference structure organized by taxonomic groups, useful for large reference sets from published datasets:

```
references/
└── nath_et_al/
    ├── lse/                          # Lineage-specific expansion references
    │   ├── annelida/
    │   │   ├── Capitella_teleta.faa
    │   │   └── ...
    │   ├── gastropoda/
    │   │   ├── Aplysia_californica.faa
    │   │   ├── Biomphalaria_glabrata.faa
    │   │   └── ...
    │   └── bivalvia/
    │       └── ...
    └── one_to_one_ortholog/          # Conserved (1:1 ortholog) references
        ├── annelida/
        ├── gastropoda/
        └── ...
```

**File Naming:**
- Files named as `Species_name.faa` (genus and species with underscore)
- Pipeline automatically extracts species names for taxonomy lookups

**Header Format:**
Headers should contain accession IDs between underscores:
```fasta
>alvmar_JABMCL010001063.1_3302
MKTIIALSYIFCLVFAQDASGRTIFDAVANTSIYTDNPSSSTSREYFTPVVEGINFDVCCF...
>aplcal_XP_005089002.1
MVTELQNNPYPVFQNLMYGIFIYALVFLLGNSVVIAISIDPHLHTPMYFFLANLSLADACF...
```

The accession ID (middle portion) is used for:
- CDS retrieval from NCBI/ENA/UniProt/DDBJ
- Taxonomy lookups when species name fails

**Processing nath_et_al References:**

1. **Rename files with taxonomy IDs:**
   ```bash
   python scripts/rename_references_with_taxid.py references/nath_et_al \
       --output-map results/reference_sequences/taxid_rename_map.csv
   ```

2. **Fetch CDS sequences for dN/dS analysis:**
   ```bash
   python scripts/fetch_reference_cds.py references/nath_et_al \
       -o results/reference_sequences/cds \
       --log-failed results/reference_sequences/failed_cds_accessions.csv
   ```

3. **Run 01_reference_processing.sh** - It will automatically detect the nath_et_al structure and process accordingly.

**CDS Fetching Notes:**
- Requires NCBI API key in `.env` file for faster downloads (10 req/sec vs 3 req/sec)
- Some accessions may not have CDS available (logged to failed_cds_accessions.csv)
- Supports NCBI RefSeq (XP_, NP_), TSA contigs, GenBank, ENA/DDBJ, and UniProt accessions

**Reference Weighting:**

When using general GPCR references (not just chemoreceptors), adjust the ranking weights in `config.sh`:

```bash
# For chemoreceptor-focused analysis (default):
export CHEMORECEPTOR_REF_WEIGHT=2.0
export OTHER_GPCR_REF_WEIGHT=1.0

# For unbiased analysis of all GPCRs:
export CHEMORECEPTOR_REF_WEIGHT=1.0
export OTHER_GPCR_REF_WEIGHT=1.0
```

### 2. Transcriptomes (Required)

Predicted protein sequences from your target species and outgroups.

#### File Location
```
transcriptomes/
├── taxid_berghia.aa             # Target species (Berghia stephanieae)
├── taxid1.aa                    # Outgroup/comparative species 1
├── taxid2.aa                    # Outgroup/comparative species 2
├── taxid_berghia.mrna           # (Optional) Nucleotide CDS for dN/dS
├── taxid_berghia.cds            # (Alternative extension)
└── ...
```

#### Naming Convention
- Protein files: `{TAXID}.aa` (amino acid FASTA)
- Nucleotide files: `{TAXID}.mrna`, `{TAXID}.cds`, or `{TAXID}.fna`
- `TAXID` must match entries in the `TAXA` array in `config.sh`

#### File Format

**Protein sequences (.aa):**
```fasta
>transcript_001 gene=chemoreceptor1
MKTIIALSYIFCLVFAQDASGRTIFDAVANTSIYTDNPSSSTSREYFTPVVEGINFDVCCF
FVVFLHLLGIGIALLSTLVAGNYVVLKIFHFKSLHTPMNFIISMLAVADLLLGFLVMPFAL
...
>transcript_002 gene=or42a
MVTELQNNPYPVFQNLMYGIFIYALVFLLGNSVVIAISIDPHLHTPMYFFLANLSLADACF
...
```

**Nucleotide CDS (.mrna/.cds) - Optional but required for dN/dS analysis:**
```fasta
>transcript_001
ATGAAAACCATTATCGCTCTTTCCTATATATTCTGCCTGGTATTTGCTCAAGATGCATCT
GGACGAACCATATTTGATGCTGTTGCTAATACAAGTATTTACACAGATAATCCATCATCA
...
```

**Requirements:**
- Headers must match between `.aa` and `.mrna`/`.cds` files for dN/dS analysis
- Nucleotide sequences must be in-frame CDS (start with ATG, length divisible by 3)
- No internal stop codons (except terminal)

#### Preparing Transcriptomes

From Trinity assembly:
```bash
# Predict ORFs with TransDecoder
TransDecoder.LongOrfs -t trinity_assembly.fasta
TransDecoder.Predict -t trinity_assembly.fasta

# Rename outputs
mv trinity_assembly.fasta.transdecoder.pep transcriptomes/taxid_berghia.aa
mv trinity_assembly.fasta.transdecoder.cds transcriptomes/taxid_berghia.cds
```

From genome annotation (GFF + genome):
```bash
# Extract protein sequences
gffread annotation.gff -g genome.fa -y transcriptomes/taxid_berghia.aa

# Extract CDS
gffread annotation.gff -g genome.fa -x transcriptomes/taxid_berghia.cds
```

### 3. Genomes (Optional)

Genome assemblies for synteny analysis and structural context.

#### File Location
```
genomes/
├── taxid_berghia.fasta          # Target species genome
├── taxid_berghia.gff            # Gene annotations (GFF3 format)
├── taxid1.fasta                 # Comparative species genome
└── taxid1.gff
```

#### File Format

**Genome FASTA:**
```fasta
>scaffold_1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG...
>scaffold_2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA...
```

**GFF3 annotations:**
```
##gff-version 3
scaffold_1	maker	gene	1000	5000	.	+	.	ID=gene001;Name=chemoreceptor1
scaffold_1	maker	mRNA	1000	5000	.	+	.	ID=mRNA001;Parent=gene001
scaffold_1	maker	CDS	1000	2000	.	+	0	ID=cds001;Parent=mRNA001
scaffold_1	maker	CDS	3000	5000	.	+	2	ID=cds002;Parent=mRNA001
```

### 4. Expression Data (Optional)

Gene expression data for tissue-specific analysis.

#### Option A: Salmon Quantification (Recommended)

```
expression_data/
├── rhinophore/
│   └── quant.sf
├── oral_veil/
│   └── quant.sf
├── tentacle/
│   └── quant.sf
├── foot/
│   └── quant.sf
└── tissue_annotations.tsv       # Optional: map sample names to tissue types
```

**quant.sf format (Salmon output):**
```
Name	Length	EffectiveLength	TPM	NumReads
transcript_001	1245	1095.2	15.234	1523
transcript_002	987	837.4	8.721	892
...
```

**tissue_annotations.tsv (optional):**
```
sample	tissue_type	is_chemosensory
rhinophore	rhinophore	yes
oral_veil	oral_veil	yes
foot	foot	no
mantle	mantle	no
```

#### Option B: Pre-computed Expression Matrix

```
expression_data.csv
```

**Format:**
```csv
gene_id,rhinophore,oral_veil,tentacle,foot,mantle
transcript_001,15.2,12.3,18.7,2.1,0.5
transcript_002,0.3,0.1,0.2,45.6,38.2
...
```

### 5. Custom HMMs (Optional)

Pre-built HMM profiles to skip HMM building in Step 01.

#### File Location
```
custom_hmms/
├── conserved.hmm                # Combined conserved GPCR HMM
└── lse.hmm                      # Combined LSE-specific HMM
```

These should be HMMER3 format profiles built with `hmmbuild`.

### 6. Configuring config.sh

After preparing your input files, update `config.sh` to match your data:

```bash
# --- Taxa Configuration ---
# List all taxa IDs (must match filenames in references/ and transcriptomes/)
export TAXA=("9606" "7227" "taxid_berghia")
export BERGHIA_TAXID="taxid_berghia"

# --- Input Files ---
# Primary transcriptome for target species
export TRANSCRIPTOME="${TRANSCRIPTOME_DIR}/taxid_berghia.aa"

# Genome (optional - set to empty string if not available)
export GENOME="${GENOME_DIR}/taxid_berghia.fasta"

# Expression data (optional)
export EXPRESSION_DATA="${BASE_DIR}/expression_data.csv"
export SALMON_QUANT_DIR="${BASE_DIR}/expression_data"

# Custom HMMs (optional - comment out to build from references)
# export CONSERVED_HMM="${BASE_DIR}/custom_hmms/conserved.hmm"
# export LSE_HMM="${BASE_DIR}/custom_hmms/lse.hmm"

# --- LSE Levels ---
# Define taxonomic groupings for LSE classification
# Format: "LevelName:taxid1,taxid2,taxid3"
export LSE_LEVELS=(
    "Aeolids:taxid_aeolid1,taxid_aeolid2,taxid_berghia"
    "Nudibranchs:taxid_nudi1,taxid_aeolid1,taxid_aeolid2,taxid_berghia"
    "Gastropods:taxid_lottia,taxid_aplysia,taxid_nudi1,taxid_berghia"
)

# Chemosensory tissues for expression analysis
export CHEMOSENSORY_TISSUES="rhinophore,oral_veil,tentacle,cephalic"
```

### 7. Validation Checklist

Before running the pipeline, verify your setup:

```bash
# Run the configuration validator
./validate_config.sh --verbose

# Check for common issues:
# ✓ All TAXA entries have matching reference files
# ✓ All TAXA entries have matching transcriptome files
# ✓ File naming conventions are correct
# ✓ FASTA files are valid (no empty sequences, proper headers)
# ✓ Nucleotide files match protein files (for dN/dS)
```

**Common issues and fixes:**

| Issue | Cause | Fix |
|-------|-------|-----|
| "No reference file for taxid X" | Filename mismatch | Ensure `{taxid}_conserved_refs.aa` exists |
| "Taxid not in TAXA array" | Missing config entry | Add taxid to `TAXA` array in config.sh |
| "Case mismatch" | Linux case sensitivity | Ensure exact case match in filenames |
| "No nucleotide sequences" | Missing .mrna/.cds file | Provide CDS file or dN/dS will be skipped |
| "Empty reference file" | Download failed | Re-download reference sequences |

### 8. Example Setup

Complete example for analyzing *Berghia stephanieae* with two outgroups:

```bash
# 1. Create directories
mkdir -p references transcriptomes genomes expression_data

# 2. Download reference GPCRs
# From GPCRdb - Class A olfactory receptors
wget "https://gpcrdb.org/services/receptorlist/?family=olfactory" -O refs_raw.json
# Process and save as FASTA...

# 3. Prepare references (example structure)
cat > references/9606_conserved_refs.aa << 'EOF'
>sp|P30872|OR1A1_HUMAN OS=Homo sapiens OX=9606 GN=OR1A1
MEFRNLTVIT...
EOF

cat > references/9606_lse_refs.aa << 'EOF'
>sp|Q8NGJ6|OR51E1_HUMAN OS=Homo sapiens OX=9606 GN=OR51E1
MDTQNFSSVT...
EOF

# 4. Copy transcriptomes
cp /path/to/berghia_proteins.fa transcriptomes/taxid_berghia.aa
cp /path/to/berghia_cds.fa transcriptomes/taxid_berghia.cds
cp /path/to/aplysia_proteins.fa transcriptomes/6500.aa
cp /path/to/lottia_proteins.fa transcriptomes/225164.aa

# 5. Copy expression data (if available)
cp -r /path/to/salmon_quants/* expression_data/

# 6. Update config.sh
sed -i 's/TAXA=.*/TAXA=("9606" "6500" "225164" "taxid_berghia")/' config.sh
sed -i 's/BERGHIA_TAXID=.*/BERGHIA_TAXID="taxid_berghia"/' config.sh

# 7. Validate
./validate_config.sh --verbose --fix
```

---

## Configuration

Edit `config.sh` to customize the pipeline. Key new parameters:

### CD-HIT Clustering Parameters

```bash
export CDHIT_IDENTITY=0.98       # 98% identity threshold
export CDHIT_WORDSIZE=5          # Word size for comparison
export CDHIT_MEMORY=16000        # Memory limit in MB
```

### CAFE5 Parameters

```bash
export CAFE_LAMBDA_SEARCH=true   # Search for optimal lambda
export CAFE_PVALUE_THRESHOLD=0.05 # Significance threshold
```

### Sensitivity Analysis Parameters

```bash
export RUN_SENSITIVITY=true              # Run Monte Carlo analysis
export SENSITIVITY_PERTURBATION=0.5      # Weight variation (+/- 50%)
export SENSITIVITY_ITERATIONS=100        # Number of iterations
```

### Cross-Validation Parameters

```bash
export RUN_CROSSVAL=false                # Run k-fold CV
export CROSSVAL_FOLDS=5                  # Number of folds
```

### Expression Analysis Parameters

```bash
export SALMON_QUANT_DIR="${BASE_DIR}/expression_data"
export MIN_TPM_THRESHOLD=1.0             # Minimum TPM for "expressed"
export TAU_THRESHOLD=0.8                 # Tissue specificity threshold
export CHEMOSENSORY_TISSUES="rhinophore,oral_veil,tentacle,cephalic"
```

### Memory Estimation Profiles

The pipeline includes data-aware resource estimation to prevent out-of-memory (OOM) failures on large datasets. These multipliers are used to estimate memory requirements based on dataset size:

```bash
export MEM_MULT_IQTREE=100       # bytes per site per taxon (ML tree inference)
export MEM_MULT_ORTHOFINDER=50   # bytes per sequence pair (all-vs-all DIAMOND)
export MEM_MULT_MAFFT=20         # bytes per residue pair (multiple alignment)
export MEM_MULT_HYPHY=200        # bytes per site per branch (aBSREL)
export MEM_MULT_FASTTREE=50      # bytes per site per taxon (approximate ML)
export MEM_BASE_GB=4             # base memory overhead in GB
export MEM_MAX_GB=128            # maximum memory to request via SLURM
```

**How it works:**
- Before running memory-intensive tools (OrthoFinder, IQ-TREE, HyPhy), the pipeline estimates required memory
- If estimated memory exceeds available resources, a warning is logged with suggestions:
  - Reduce input size with CD-HIT clustering
  - Request more memory via SLURM
  - Use a high-memory node
- Tune the multipliers based on your empirical observations

### Sequence Metadata System

The pipeline uses a centralized metadata CSV (`id_map.csv`) instead of parsing FASTA headers directly. This eliminates fragility when reference database formats change.

**Generated columns:**
- `original_id` - Original sequence identifier
- `short_id` - Standardized short ID (e.g., `ref_9606_1`)
- `taxid` - Taxonomic identifier
- `source_type` - reference/target/outgroup
- `species_name` - Species name (if extractable)
- `gene_name` - Gene name (if extractable)
- `seq_length` - Sequence length

**Usage in custom scripts:**
```python
from get_metadata import MetadataLookup
lookup = MetadataLookup('results/reference_sequences/id_map.csv')
taxid = lookup.get_taxid('ref_9606_1')
source = lookup.get_source_type('ref_9606_1')
```

---

## Pipeline Steps

### New Steps

#### Step 02a: Sequence Clustering (`02a_cluster_sequences.sh`)

**Purpose:** Remove redundant sequences (splice variants, alleles) using CD-HIT.

**Outputs:**
- `results/clustering/candidates_nr098.fa` - Non-redundant sequences
- `results/clustering/cluster_mapping.tsv` - Cluster membership
- `results/clustering/cluster_summary.tsv` - Cluster statistics

#### Step 02b: GPCR Classification (`02b_classify_gpcrs.sh`)

**Purpose:** Classify GPCRs using InterProScan domains and sub-family assignment.

**Outputs:**
- `results/classification/gpcr_classifications.tsv` - Class and sub-family assignments
- `results/classification/putative_chemoreceptors.txt` - Likely chemosensory candidates
- `results/classification/classification_summary.json` - Distribution statistics

#### Step 03c: CAFE5 Analysis (`03c_cafe_analysis.sh`)

**Purpose:** Formal LSE detection using birth-death models.

**Process:**
1. Generates ultrametric species tree using chronos (R/ape)
2. Prepares gene family count matrix from OrthoFinder
3. Runs CAFE5 with lambda estimation
4. Identifies significant expansions/contractions

**Outputs:**
- `results/cafe/species_tree_ultrametric.tre` - Dated species tree
- `results/cafe/significant_expansions.tsv` - Significant gene family changes
- `results/cafe/cafe_summary.json` - Analysis summary

#### Step 03d: NOTUNG Reconciliation (`03d_notung_reconciliation.sh`)

**Purpose:** Reconcile gene trees with species tree to infer duplication/loss events.

**Outputs:**
- `results/notung/reconciliation_summary.tsv` - Duplication/loss counts
- `results/notung/expanded_orthogroups.txt` - Orthogroups with expansions
- `results/notung/events_summary.json` - Per-species event counts

---

## Python Utilities

### New Analysis Scripts

| Script | Purpose | Usage |
|--------|---------|-------|
| `scripts/get_metadata.py` | Centralized sequence metadata lookup | `python3 get_metadata.py -m id_map.csv --seq-id ref_9606_1` |
| `scripts/expression_analysis.py` | Parse Salmon quant.sf, calculate tau index | `python3 expression_analysis.py <quant_dir> <candidates> <output_prefix>` |
| `scripts/divergence_dating.py` | Estimate divergence times with calibrations | `python3 divergence_dating.py <tree> <output_prefix> [calibrations]` |
| `scripts/convergent_evolution.py` | Detect convergent amino acid substitutions | `python3 convergent_evolution.py <alignment> <tree> <output_prefix>` |
| `scripts/binding_site_prediction.py` | Predict ligand binding residues from structures | `python3 binding_site_prediction.py <pdb_dir> <sequences> <output_prefix>` |
| `scripts/conservation_mapping.py` | Per-residue conservation with motif detection | `python3 conservation_mapping.py <alignment> <output_prefix> [pdb_dir]` |
| `scripts/generate_ultrametric.R` | Convert ML tree to ultrametric using chronos | `Rscript generate_ultrametric.R <tree> <output> [model]` |

### Enhanced Scripts

| Script | New Features |
|--------|--------------|
| `rank_candidates.py` | Sensitivity analysis, cross-validation, rank stability metrics |
| `update_headers.py` | Extended metadata CSV with source_type, species_name, gene_name columns |
| `functions.sh` | Checkpointing, provenance tracking, resource estimation, metadata lookup helpers |

---

## Output Files

### New Key Results

| File | Description |
|------|-------------|
| `results/ranking/sensitivity_analysis.csv` | Rank stability across weight perturbations |
| `results/ranking/weight_importance.json` | Relative importance of each scoring component |
| `results/classification/gpcr_classifications.tsv` | GPCR class and sub-family assignments |
| `results/cafe/significant_expansions.tsv` | Statistically significant gene family changes |
| `results/clustering/cluster_mapping.tsv` | Representative-to-cluster member mapping |

---

## Running the Pipeline

### Validate Configuration First

Before running the pipeline, validate your configuration to catch errors early:

```bash
# Basic validation
./validate_config.sh

# With suggested fixes for any issues
./validate_config.sh --fix

# Verbose output showing all checks
./validate_config.sh --verbose
```

The validator checks:
- **Directory existence** - All required directories exist
- **Taxa consistency** - TAXA array matches actual transcriptome files
- **LSE levels** - LSE_LEVELS taxids are all in TAXA array
- **Case sensitivity** - No case mismatches between config and files
- **Tool availability** - Required tools are in PATH
- **Parameter ranges** - Settings are within expected ranges
- **Local databases** - LOCAL_DB_DIR is properly set up (if used)

### Full Pipeline with New Steps

```bash
# Validate configuration first
./validate_config.sh

# Source configuration
source config.sh
source functions.sh

# Core identification
./01_reference_processing.sh
./02_chemogpcrs_identification.sh
./02a_cluster_sequences.sh        # NEW: Remove redundancy
./02b_classify_gpcrs.sh           # NEW: Classify GPCRs

# Orthology and evolution
./03_orthology_clustering.sh
./03a_busco_species_tree.sh
./03b_lse_classification.sh
./03c_cafe_analysis.sh            # NEW: CAFE5 analysis
./03d_notung_reconciliation.sh    # NEW: Gene tree reconciliation

# Phylogenetics and selection
./04_phylogenetic_analysis.sh
./05_selective_pressure_and_asr.sh
./06_synteny_and_mapping.sh

# Ranking and structure
./07_candidate_ranking.sh
./08_structural_analysis.sh
./09_report_generation.sh
```

### Resume from Checkpoint

The pipeline now includes checkpointing for automatic resume:

```bash
# Check completed steps
ls results/checkpoints/

# Resume automatically - completed steps are skipped
./04_phylogenetic_analysis.sh

# Force re-run of a step
./04_phylogenetic_analysis.sh --force
```

### Run Mode Selection

```bash
# Run locally (no SLURM)
export RUN_MODE=local
./07_candidate_ranking.sh

# Run with SLURM
export RUN_MODE=slurm
sbatch 07_candidate_ranking.sh
```

### SLURM Configuration

The pipeline includes a wrapper script (`run_pipeline.sh`) for easy SLURM job submission with proper configuration.

#### Configuration Variables (config.sh)

```bash
# --- SLURM Configuration ---
export DEFAULT_TIME="24:00:00"          # Wall time limit
export DEFAULT_MEM="64G"                # Memory per job
export SLURM_EMAIL="you@example.com"    # Notification email

# Cluster-specific options (leave empty to omit)
export SLURM_PARTITION="cpu"            # Partition/queue name
export SLURM_ACCOUNT="katzlab"          # Account/allocation
export SLURM_QOS="normal"               # Quality of service
export SLURM_CONSTRAINT=""              # Node constraints (e.g., "skylake")
export SLURM_RESERVATION=""             # Reservation name
export SLURM_EXTRA_ARGS=""              # Additional sbatch flags

# Array job settings
export SLURM_ARRAY_LIMIT=50             # Max concurrent array tasks
```

#### Using the Wrapper Script

```bash
# Submit a single step
./run_pipeline.sh 02

# Submit multiple steps
./run_pipeline.sh 02 03 04

# Submit a range of steps
./run_pipeline.sh 02-05

# Submit all steps sequentially
./run_pipeline.sh all

# Wait for each job to complete before submitting next
./run_pipeline.sh --wait 04 05

# Run locally (no SLURM)
./run_pipeline.sh --local 07

# Dry run (show commands without executing)
./run_pipeline.sh --dry-run all

# Verbose output showing SLURM settings
./run_pipeline.sh --verbose 02
```

#### Direct sbatch Submission

You can also submit scripts directly with custom options:

```bash
# Override partition and account
sbatch --partition=gpu --account=mylab 08_structural_analysis.sh

# Add job dependency
sbatch --dependency=afterok:12345 05_selective_pressure_and_asr.sh

# Request GPU resources
sbatch --gres=gpu:1 --partition=gpu 08_structural_analysis.sh
```

#### Monitoring Jobs

```bash
# Check your jobs
squeue -u $USER

# Check specific job
scontrol show job <job_id>

# Cancel a job
scancel <job_id>

# View job output in real-time
tail -f results/logs/02_chemogpcrs_id_<job_id>.out
```

---

## Troubleshooting

### New Common Issues

**1. CD-HIT memory issues**
```
Error: Not enough memory
```
Solution: Adjust memory limit:
```bash
export CDHIT_MEMORY=32000  # Increase to 32GB
```

**2. CAFE5 ultrametric tree error**
```
Error: Tree is not ultrametric
```
Solution: The pipeline uses chronos to convert trees. Check R installation:
```bash
Rscript -e "library(ape); print('OK')"
```

**3. InterProScan not found**
```
Warning: InterProScan not available
```
Solution: Install InterProScan or skip classification:
```bash
conda install -c bioconda interproscan
# Or: pipeline falls back to domain-based classification
```

**4. Sensitivity analysis slow**
```
Sensitivity analysis taking too long
```
Solution: Reduce iterations:
```bash
export SENSITIVITY_ITERATIONS=50
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
- CAFE5: Mendes et al. (2020) Bioinformatics 36:5516-5518
- HyPhy: Kosakovsky Pond et al. (2020) Mol Biol Evol 37:295-299
- CD-HIT: Li & Godzik (2006) Bioinformatics 22:1658-1659
- AlphaFold: Jumper et al. (2021) Nature 596:583-589

---

## License

This pipeline is provided for academic research purposes. See individual tool licenses for commercial use restrictions.
