# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Bioinformatics pipeline for identifying and characterizing chemoreceptor GPCR genes in *Berghia stephanieae* (sea slug). Integrates protein homology detection (HHblits, HMMSEARCH), phylogenetic inference (IQ-TREE, Phyloformer), orthology clustering (OrthoFinder), selective pressure analysis (HyPhy aBSREL), synteny mapping (MCScanX), ancestral sequence reconstruction (FastML), and structural prediction (AlphaFold).

## Running the Pipeline

Each numbered shell script is an independent SLURM job representing one pipeline stage:

```bash
# Submit to SLURM
sbatch 01_reference_processing.sh
sbatch 02_chemogpcrs_identification.sh
# ... through 09_report_generation.sh

# Run locally
bash 07_candidate_ranking.sh
```

Scripts must run sequentially (01 → 09). Each checks for prerequisite `step_completed_*.txt` markers before executing.

## Configuration

All paths, tool locations, and parameters are in `config.sh`. Scripts source this plus `functions.sh` at startup.

Key variables:
- `BASE_DIR`, `RESULTS_DIR`, `SCRIPTS_DIR`, `LOGS_DIR` - directory paths
- Tool paths (28 tools): `IQTREE`, `HMMSEARCH`, `MAFFT`, `ALPHAFOLD`, etc.
- `PHYLO_WEIGHT`, `DNDS_WEIGHT`, `SYNTENY_WEIGHT`, `EXPR_WEIGHT` - ranking weights
- `LSE_LEVELS` - taxonomic groupings for lineage-specific expansion classification

## Architecture

### Pipeline Stages
1. **01_reference_processing.sh** - Build HMMs, prepare reference sequences
2. **02_chemogpcrs_identification.sh** - GPCR detection via HHblits + DeepTMHMM (6+ TM regions)
3. **03_orthology_clustering.sh** - OrthoFinder clustering
4. **03a_busco_species_tree.sh** - Build species tree from BUSCO genes
5. **03b_lse_classification.sh** - Classify orthogroups into LSE taxonomic levels
6. **04_phylogenetic_analysis.sh** - IQ-TREE/Phyloformer trees, visualizations (SLURM array job)
7. **05_selective_pressure_and_asr.sh** - HyPhy aBSREL dN/dS, FastML ASR (SLURM array job)
8. **06_synteny_and_mapping.sh** - Minimap2 mapping, MCScanX synteny
9. **07_candidate_ranking.sh** - Multi-criteria scoring and ranking
10. **08_structural_analysis.sh** - AlphaFold prediction, FoldTree comparison
11. **09_report_generation.sh** - LaTeX PDF report

### Python Utilities (20 scripts)
- `rank_candidates.py` - Multi-criteria weighted scoring using phylogeny, dN/dS, expression, synteny
- `parse_absrel.py` - Parse HyPhy aBSREL JSON output to CSV (branch_id, omega, p_value)
- `visualize_tree.py` - Tree visualization (basic, colored by taxon, circular layouts)
- `plot_*.py` - Various visualization scripts (ASR, synteny, heatmaps, PCA)
- `select_diverse_candidates.py` - Hierarchical clustering selection
- `select_deep_nodes.py` - ASR deep node selection (90th percentile threshold)
- `update_headers.py` - FASTA header standardization
- `test_phyloformer_models.py` - Model selection for Phyloformer

### Utility Functions (functions.sh)
- `log()` - Timestamped logging
- `run_command()` - Execute with logging, error handling, completion markers
- `check_file()` - Verify input files exist (fail-fast)

## Data Flow

Results stored in `${RESULTS_DIR}/` with subdirectories: `phylogenies/`, `selective_pressure/`, `synteny/`, `ranking/`, `structural_analysis/`

ID mapping uses standardized short headers (`ref_1`, `ref_2`, etc.) throughout pipeline.

Logs in `${LOGS_DIR}/pipeline.log` and per-command `.log`/`.err` files.

## Key Dependencies

Python: ete3, pandas, scipy, matplotlib, seaborn, Biopython, requests

External tools: HMMBUILD, HMMSEARCH, HHBLITS, MAFFT, TRIMAL, ClipKit, IQ-TREE, FastTree, Phyloformer, HyPhy, FastML, OrthoFinder, MCScanX, BUSCO, Minimap2, Samtools, AlphaFold, FoldTree, TMalign

## Key Implementation Details

- **FastTree Seed Strategy**: Step 04 generates approximate ML trees with FastTree before running IQ-TREE, avoiding local optima in large GPCR gene families
- **dN/dS Scoring**: `rank_candidates.py` parses aBSREL results from `absrel_results.csv`, scoring both strong purifying selection (conserved function) and positive selection (novel function) as interesting candidates
- **Deep Node Selection**: `select_deep_nodes.py` uses 90th percentile of branch distances as threshold rather than arbitrary cutoffs
- **Synteny Parsing**: `plot_synteny.py` properly parses MCScanX `.collinearity` format including alignment headers and gene pair lines
- **Tree Visualization**: `visualize_tree.py` colors leaves by taxon type (reference, Berghia, other) with legends
