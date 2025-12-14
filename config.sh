#!/bin/bash
# config.sh
# Purpose: Define global variables and paths for the pipeline.
# Usage: Source this file in all pipeline scripts (source config.sh).

# --- Base Directories ---
export BASE_DIR=$(realpath "$(dirname "$0")")
export RESULTS_DIR="${BASE_DIR}/results"
export SCRIPTS_DIR="${BASE_DIR}/scripts"
export REFERENCE_DIR="${BASE_DIR}/references"
export TRANSCRIPTOME_DIR="${BASE_DIR}/transcriptomes"
export GENOME_DIR="${BASE_DIR}/genomes"
export LOGS_DIR="${RESULTS_DIR}/logs"

# --- Input Files ---
export TRANSCRIPTOME="${TRANSCRIPTOME_DIR}/taxid_berghia_berghia.aa"
export CONSERVED_HMM="${BASE_DIR}/custom_hmms/conserved.hmm"  # Optional
export LSE_HMM="${BASE_DIR}/custom_hmms/lse.hmm"              # Optional
export ID_MAP="${RESULTS_DIR}/reference_sequences/id_map.csv"
export EXPRESSION_DATA="${BASE_DIR}/expression_data.csv"
export GENOME="${GENOME_DIR}/taxid_berghia_berghia.fasta"

# --- Tool Paths (update based on Conda environment) ---
export HMMBUILD="hmmbuild"
export HMMSEARCH="hmmsearch"
export HHMAKE="hhmake"
export HHBLITS="hhblits"
export MAFFT="mafft"
export IQTREE="iqtree2"
export FASTTREE="FastTree"
export TRIMAL="trimal"
export FASTML="fastml"
export MINIMAP2="minimap2"
export SAMTOOLS="samtools"
export MCSCANX="MCScanX"
export ALPHAFOLD="run_alphafold.sh"
export FOLDTREE="foldtree"
export TMALIGN="TMalign"
export SEQTK="seqtk"
export ORTHOFINDER="orthofinder"
export PHYLOFORMER="phyloformer"
export PDFLATEX="pdflatex"
export BUSCO="busco"
export DEEPTMHMM="deeptmhmm"
export ASTRAL="astral"
export MRBAYES="mb"
export CLIPKIT="clipkit"

# --- Pipeline Parameters ---
export TAXA=("taxid1" "taxid2" "taxid_berghia")
export BERGHIA_TAXID="taxid_berghia"
export CPUS=16
export DEFAULT_TIME="24:00:00"
export DEFAULT_MEM="64G"
export SLURM_EMAIL="your.email@example.com"
export GPU_ENABLED=false

# --- Tool-Specific Parameters ---
export HMM_EVALUE="1e-5"
export HHBLITS_EVALUE="1e-5"
export MIN_TM_REGIONS=6
export MIN_SEQ_LENGTH=100
export MAX_GAP_PERCENT=50
export IQTREE_MODEL="TEST"
export IQTREE_BOOTSTRAP=1000
export ORTHOFINDER_INFLATION=1.5
export FOLDTREE_METHOD="upgma"
export USE_MRBAYES=false
export NUM_STRUCTURAL_CANDIDATES=10
export MIN_ASR_DISTANCE=0.5  # Minimum distance for deep node selection in ASR

# --- Ranking Weights ---
export PHYLO_WEIGHT=2
export DNDS_WEIGHT=1
export SYNTENY_WEIGHT=3
export EXPR_WEIGHT=1
export LSE_DEPTH_WEIGHT=1

# --- GPCRdb Fetch Parameters ---
export GPCRDB_SEARCH_TERMS="chemoreceptor,invertebrate"
export GPCRDB_SPECIES="Aplysia,Lottia"
# GPCR families: "all" for all families, or comma-separated: "Class_A,Class_B1,Class_C,Adhesion,Frizzled"
export GPCRDB_FAMILIES="all"

# --- Taxonomic Levels for LSE ---
export LSE_LEVELS=("Aeolids:taxid_aeolid1,taxid_aeolid2" "Nudibranchs:taxid_nudi1,taxid_nudi2" "Gastropods:taxid_gastro1,taxid_gastro2")

# --- NCBI Taxonomy IDs for LSE Classification ---
# These are used by lse_refine.py to determine taxonomic levels
export LSE_AEOLID_TAXID=54397      # Aeolidida (aeolid nudibranchs)
export LSE_NUDIBRANCH_TAXID=13843  # Nudibranchia (all nudibranchs)
export LSE_GASTROPOD_TAXID=644     # Gastropoda (all gastropods)

# --- Species Tree ---
export SPECIES_TREE="${RESULTS_DIR}/busco/busco_species_tree.tre"
