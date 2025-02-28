#!/bin/bash
# config.sh
# Purpose: Define pipeline-wide variables, paths, and customizable settings.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

# Base directories
PIPELINE_DIR="/path/to/pipeline"
DATA_DIR="${PIPELINE_DIR}/data"
RESULTS_DIR="${PIPELINE_DIR}/results"
LOGS_DIR="${PIPELINE_DIR}/logs"
TRANSCRIPTOME_DIR="${PIPELINE_DIR}/transcriptomes"
GENOME_DIR="${PIPELINE_DIR}/genomes"
REFERENCE_DIR="${PIPELINE_DIR}/gpcr_references"
QUANT_DIR="${PIPELINE_DIR}/quantification"
SCRIPTS_DIR="${PIPELINE_DIR}/scripts"

# Input data
GENOME="${DATA_DIR}/berghia_genome.fa"
TRANSCRIPTOME="${DATA_DIR}/berghia_transcriptome.fa"

# Specifiable taxa (NCBI TaxIDs)
TAXA=("9606" "10090" "66289" "6447")  # Human, Mouse, Aplysia, Lottia
BERGHIA_TAXID="1263399"              # Berghia stephanieae

# Custom HMMs (optional; leave empty to build from references)
CONSERVED_HMM="${REFERENCE_DIR}/conserved.hmm"
LSE_HMM="${REFERENCE_DIR}/lse.hmm"

# Resource limits
DEFAULT_CPUS=16  # Global cap, overridden where needed
DEFAULT_MEM="16G"
DEFAULT_TIME="24:00:00"
GPU_ENABLED=true

# Tool paths
HMMBUILD="hmmbuild"
HMMSEARCH="hmmsearch"
HHBLITS="hhblits"
HHMAKE="hhmake"
MAFFT="mafft"
IQTREE="iqtree2"
FASTTREE="FastTree"
TRIMAL="trimal"
CODEML="codeml"
FASTML="fastml"
ETE3="ete3"
MINIMAP2="minimap2"
SAMTOOLS="samtools"
MCSCANX="MCScanX"
IADHORE="/path/to/i-adhore"
ALPHAFOLD="/path/to/alphafold"
FOLDTREE="/path/to/foldtree"
TMALIGN="/path/to/TMalign"
SEQTK="seqtk"
ORTHOFINDER="orthofinder"
PHYLOFORMER="/path/to/phyloformer"
PDFLATEX="pdflatex"
BUSCO="busco"
DEEPTMHMM="deeptmhmm"
ASTRAL="astral"
MRBAYES="mb"
BMGE="bmge"

# Analysis settings
HMM_EVALUE="1e-5"
HHBLITS_EVALUE="1e-5"
BLAST_EVALUE="1e-5"
IQTREE_MODEL="MFP"
IQTREE_BOOTSTRAP=1000
ORTHOFINDER_INFLATION=1.2
TPM_THRESHOLD=10
TPM_WEIGHT_RHINOPHORE=2
TPM_WEIGHT_ORAL=1
MIN_TM_REGIONS=6
MIN_SEQ_LENGTH=100
MAX_GAP_PERCENT=50
USE_MRBAYES=false  # Toggle Bayesian inference

# Output files
ID_MAP="${RESULTS_DIR}/reference_sequences/id_map.csv"
SPECIES_TREE="${RESULTS_DIR}/busco/busco_species_tree.tre"

# SLURM settings
SLURM_PARTITION="compute"
SLURM_NODES=1
SLURM_NTASKS=1
SLURM_EMAIL="jperezmoreno@umass.edu"
