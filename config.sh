#!/bin/bash
# config.sh
# Purpose: Define pipeline-wide variables, paths, and customizable settings for GPCR analysis in Berghia stephanieae.
# Notes: This script is sourced by all other scripts to ensure consistent configuration.

# --- Directories ---
PIPELINE_DIR="/path/to/pipeline"                # Root directory for the pipeline
DATA_DIR="${PIPELINE_DIR}/data"                 # Input data storage
RESULTS_DIR="${PIPELINE_DIR}/results"           # Output results storage
LOGS_DIR="${PIPELINE_DIR}/logs"                 # Log files storage
TRANSCRIPTOME_DIR="${PIPELINE_DIR}/transcriptomes"  # Directory for transcriptome files
GENOME_DIR="${PIPELINE_DIR}/genomes"            # Directory for genome files
REFERENCE_DIR="${PIPELINE_DIR}/gpcr_references" # Directory for reference sequences and HMMs
QUANT_DIR="${PIPELINE_DIR}/quantification"      # Directory for expression data
SCRIPTS_DIR="${PIPELINE_DIR}/scripts"           # Directory for supporting scripts

# --- Input Files ---
GENOME="${DATA_DIR}/berghia_genome.fa"          # Berghia genome FASTA file
TRANSCRIPTOME="${DATA_DIR}/berghia_transcriptome.fa"  # Berghia transcriptome FASTA file

# --- Taxonomic Settings ---
TAXA=("9606" "10090" "66289" "6447")           # TaxIDs: Human, Mouse, Aplysia, Lottia
BERGHIA_TAXID="1263399"                        # TaxID for Berghia stephanieae

# --- Custom HMMs (Optional) ---
CONSERVED_HMM="${REFERENCE_DIR}/conserved.hmm"  # Optional custom HMM for conserved GPCRs
LSE_HMM="${REFERENCE_DIR}/lse.hmm"              # Optional custom HMM for LSE GPCRs

# --- Resource Limits ---
DEFAULT_CPUS=16                                 # Default number of CPUs for parallel tasks
DEFAULT_MEM="16G"                               # Default memory allocation
DEFAULT_TIME="24:00:00"                         # Default SLURM job time limit
GPU_ENABLED=true                                # Enable GPU usage if available

# --- Tool Paths ---
# Note: Ensure these tools are installed and accessible in your environment
HMMBUILD="hmmbuild"                            # Path to HMMER hmmbuild
HMMSEARCH="hmmsearch"                          # Path to HMMER hmmsearch
HHBLITS="hhblits"                              # Path to HHblits
HHMAKE="hhmake"                                # Path to HHmake
MAFFT="mafft"                                  # Path to MAFFT aligner
IQTREE="iqtree2"                               # Path to IQ-TREE
FASTTREE="FastTree"                            # Path to FastTree
TRIMAL="trimal"                                # Path to TrimAl
CODEML="codeml"                                # Path to PAML codeml
FASTML="fastml"                                # Path to FastML
ETE3="ete3"                                    # Path to ETE3 (used via Python)
MINIMAP2="minimap2"                            # Path to minimap2
SAMTOOLS="samtools"                            # Path to samtools
MCSCANX="MCScanX"                              # Path to MCScanX
ALPHAFOLD="/path/to/alphafold"                 # Path to AlphaFold
FOLDTREE="/path/to/foldtree"                   # Path to FoldTree
TMALIGN="/path/to/TMalign"                     # Path to TMalign
SEQTK="seqtk"                                  # Path to seqtk
ORTHOFINDER="orthofinder"                      # Path to OrthoFinder
PDFLATEX="pdflatex"                            # Path to pdflatex
BUSCO="busco"                                  # Path to BUSCO
DEEPTMHMM="deeptmhmm"                          # Path to DeepTMHMM
ASTRAL="astral"                                # Path to ASTRAL

# --- Analysis Settings ---
HMM_EVALUE="1e-5"                              # E-value threshold for HMMsearch
HHBLITS_EVALUE="1e-5"                          # E-value threshold for HHblits
IQTREE_MODEL="MFP"                             # Model for IQ-TREE (ModelFinder Plus)
IQTREE_BOOTSTRAP=1000                          # Number of bootstrap replicates for IQ-TREE
ORTHOFINDER_INFLATION=1.2                      # Inflation parameter for OrthoFinder
TPM_THRESHOLD=10                               # Minimum TPM for expression filtering
TPM_WEIGHT_RHINOPHORE=2                        # Weight for rhinophore TPM in ranking
TPM_WEIGHT_ORAL=1                              # Weight for oral TPM in ranking
MIN_TM_REGIONS=6                               # Minimum transmembrane regions for GPCR filtering
MIN_SEQ_LENGTH=100                             # Minimum sequence length for alignments
MAX_GAP_PERCENT=50                             # Maximum gap percentage in alignments
FOLDTREE_METHOD="tmalign"                      # Method for FoldTree structural phylogeny

# --- Output Files ---
ID_MAP="${RESULTS_DIR}/reference_sequences/id_map.csv"  # Mapping of short IDs to original IDs
SPECIES_TREE="${RESULTS_DIR}/busco/busco_species_tree.tre"  # Species tree from BUSCO/ASTRAL

# --- SLURM Settings ---
SLURM_PARTITION="compute"                      # SLURM partition name
SLURM_NODES=1                                  # Number of SLURM nodes
SLURM_NTASKS=1                                 # Number of SLURM tasks
SLURM_EMAIL="jperezmoreno@umass.edu"           # Email for SLURM notifications
