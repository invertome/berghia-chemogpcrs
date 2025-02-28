#!/bin/bash
# config.sh
# Purpose: Define pipeline-wide variables, paths, and customizable settings for GPCR analysis in Berghia stephanieae.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

# --- Directories ---
# Base directory for the pipeline
PIPELINE_DIR="/path/to/pipeline"
# Input data directory
DATA_DIR="${PIPELINE_DIR}/data"
# Results output directory
RESULTS_DIR="${PIPELINE_DIR}/results"
# Logging directory
LOGS_DIR="${PIPELINE_DIR}/logs"
# Transcriptome input directory
TRANSCRIPTOME_DIR="${PIPELINE_DIR}/transcriptomes"
# Genome input directory
GENOME_DIR="${PIPELINE_DIR}/genomes"
# Reference sequences directory
REFERENCE_DIR="${PIPELINE_DIR}/gpcr_references"
# Quantification data directory
QUANT_DIR="${PIPELINE_DIR}/quantification"
# Scripts directory
SCRIPTS_DIR="${PIPELINE_DIR}/scripts"

# --- Input Files ---
# Genome FASTA file for Berghia
GENOME="${DATA_DIR}/berghia_genome.fa"
# Transcriptome FASTA file for Berghia
TRANSCRIPTOME="${DATA_DIR}/berghia_transcriptome.fa"

# --- Taxonomic Settings ---
# NCBI TaxIDs for species included in the analysis
TAXA=("9606" "10090" "66289" "6447")  # Human, Mouse, Aplysia, Lottia
# Berghia stephanieae TaxID
BERGHIA_TAXID="1263399"

# --- Custom HMMs (Optional) ---
# Path to custom conserved HMM (leave empty to build from references)
CONSERVED_HMM="${REFERENCE_DIR}/conserved.hmm"
# Path to custom LSE HMM (leave empty to build from references)
LSE_HMM="${REFERENCE_DIR}/lse.hmm"

# --- Resource Limits ---
# Default CPU count (global cap, overridden where specified)
DEFAULT_CPUS=16
# Default memory allocation
DEFAULT_MEM="16G"
# Default job time limit
DEFAULT_TIME="24:00:00"
# Enable GPU usage for AlphaFold
GPU_ENABLED=true

# --- Tool Paths ---
# Paths to bioinformatics tools (assumed in PATH or Conda environment)
HMMBUILD="hmmbuild"          # HMMER tool for building HMMs
HMMSEARCH="hmmsearch"        # HMMER tool for searching with HMMs
HHBLITS="hhblits"            # HHSuite tool for profile-profile alignment
HHMAKE="hhmake"              # HHSuite tool for creating HMMs
MAFFT="mafft"                # Multiple sequence alignment tool
IQTREE="iqtree2"             # Phylogenetic tree inference tool
FASTTREE="FastTree"          # Fast phylogenetic tree inference
TRIMAL="trimal"              # Alignment trimming tool
CODEML="codeml"              # PAML tool for codon-based selection analysis
FASTML="fastml"              # Ancestral sequence reconstruction tool
ETE3="ete3"                  # Python library for tree manipulation
MINIMAP2="minimap2"          # Read mapping tool
SAMTOOLS="samtools"          # SAM/BAM file processing tool
MCSCANX="MCScanX"            # Synteny analysis tool
IADHORE="/path/to/i-adhore"  # Synteny analysis tool (custom path)
ALPHAFOLD="/path/to/alphafold"  # Protein structure prediction tool
FOLDTREE="/path/to/foldtree"    # Structural phylogeny tool
TMALIGN="/path/to/TMalign"      # Structural alignment tool
SEQTK="seqtk"                   # Sequence manipulation tool
ORTHOFINDER="orthofinder"       # Orthology inference tool
PHYLOFORMER="/path/to/phyloformer"  # Deep learning phylogeny tool
PDFLATEX="pdflatex"             # LaTeX compiler for report generation
BUSCO="busco"                   # BUSCO gene identification tool
DEEPTMHMM="deeptmhmm"           # Transmembrane prediction tool
ASTRAL="astral"                 # Coalescent-based species tree inference
MRBAYES="mb"                    # Bayesian phylogenetic inference tool
BMGE="bmge"                     # Alignment trimming tool
CLIPKIT="clipkit"               # Phylogenetic-informed alignment trimming

# --- Analysis Settings ---
# E-value threshold for HMM searches
HMM_EVALUE="1e-5"
# E-value threshold for HHblits searches
HHBLITS_EVALUE="1e-5"
# E-value threshold for BLAST in synteny analysis
BLAST_EVALUE="1e-5"
# Model for IQ-TREE (ModelFinder Plus)
IQTREE_MODEL="MFP"
# Number of bootstrap replicates for IQ-TREE
IQTREE_BOOTSTRAP=1000
# Inflation parameter for OrthoFinder clustering
ORTHOFINDER_INFLATION=1.2
# TPM threshold for expression filtering
TPM_THRESHOLD=10
# Weight for rhinophore expression
TPM_WEIGHT_RHINOPHORE=2
# Weight for oral-tentacle expression
TPM_WEIGHT_ORAL=1
# Minimum number of transmembrane regions for GPCR filtering
MIN_TM_REGIONS=6
# Minimum sequence length for quality control
MIN_SEQ_LENGTH=100
# Maximum gap percentage allowed in alignments
MAX_GAP_PERCENT=50
# Toggle Bayesian inference with MrBayes
USE_MRBAYES=false
# Toggle ClipKIT trimming (alternative to BMGE)
USE_CLIPKIT=false
# FoldTree structural alignment method (options: tmalign, rmsd)
FOLDTREE_METHOD="tmalign"

# --- Output Files ---
# Mapping file for reference sequence IDs
ID_MAP="${RESULTS_DIR}/reference_sequences/id_map.csv"
# Species tree file from BUSCO/ASTRAL
SPECIES_TREE="${RESULTS_DIR}/busco/busco_species_tree.tre"

# --- SLURM Settings ---
# Partition for SLURM jobs
SLURM_PARTITION="compute"
# Number of nodes per job
SLURM_NODES=1
# Number of tasks per job
SLURM_NTASKS=1
# Email for SLURM notifications
SLURM_EMAIL="jperezmoreno@umass.edu"
