#!/bin/bash
# 09_report_generation.sh
# Purpose: Generate a comprehensive LaTeX report summarizing the pipeline results.
# Inputs: Results from all previous steps
# Outputs: PDF report in ${RESULTS_DIR}/pipeline_report.pdf
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=report_generation
#SBATCH --output=${LOGS_DIR}/09_report_generation_%j.out
#SBATCH --error=${LOGS_DIR}/09_report_generation_%j.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH $(scale_resources)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Create report directory
mkdir -p "${RESULTS_DIR}/report" "${LOGS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

# Check dependency
check_file "${RESULTS_DIR}/step_completed_foldtree.txt"

log "Starting report generation."

# Generate LaTeX report
cat > "${RESULTS_DIR}/report/report.tex" <<'EOF'
\documentclass{article}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{geometry}
\geometry{a4paper, margin=1in}
\title{GPCR Chemoreceptor Analysis in Berghia stephanieae}
\author{Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst}
\date{\today}

\begin{document}

\maketitle

\section{Introduction}
This report presents the results of a comprehensive analysis of chemoreceptive GPCRs in \textit{Berghia stephanieae}, integrating phylogeny, selective pressure, synteny, and structural data.

\subsection{Phylogenetic Analysis}
Trees were constructed using IQ-TREE and Phyloformer, with optional MrBayes analysis.

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{../phylogenies/visualizations/all_berghia_refs_basic.png}
\caption{IQ-TREE phylogeny of all Berghia candidates and references.}
\label{fig:iqtree}
\end{figure}

\subsection{Selective Pressure and ASR}
Selective pressure was analyzed with HyPhy aBSREL, and ancestral sequences were reconstructed with FastML.

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{../asr/OG0000001_asr_plot_circular.png}
\caption{ASR visualization for a representative orthogroup.}
\label{fig:asr}
\end{figure}

\subsection{Synteny Analysis}
Synteny blocks were identified using MCScanX across available genomes.

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{../synteny/synteny_plot_Aeolids.png}
\caption{Synteny at Aeolids level.}
\label{fig:synteny}
\end{figure}

\subsection{Candidate Ranking}
Candidates were ranked based on integrated data.

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{../ranking/ranking_plot.png}
\caption{Top-ranked GPCR candidates.}
\label{fig:ranking}
\end{figure}

\subsection{Structural Analysis}
Structural phylogenies included AlphaFold predictions and GPCRdb references.

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{../structural_analysis/struct_vs_seq_plot.png}
\caption{Comparison of structural and sequence phylogenies.}
\label{fig:struct_vs_seq}
\end{figure}

\end{document}
EOF

# Compile LaTeX
run_command "latex_compile" ${PDFLATEX} -output-directory="${RESULTS_DIR}/report" "${RESULTS_DIR}/report/report.tex"

# Move PDF
mv "${RESULTS_DIR}/report/report.pdf" "${RESULTS_DIR}/pipeline_report.pdf" || { log "Error: Failed to move PDF"; exit 1; }

log "Report generation completed."
