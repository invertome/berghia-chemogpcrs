#!/bin/bash
# 09_report_generation.sh
# Purpose: Generate a comprehensive LaTeX report summarizing the pipeline results.
# Inputs: All analysis results from previous steps
# Outputs: PDF report (${RESULTS_DIR}/pipeline_report.pdf)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

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
mkdir -p "${RESULTS_DIR}/report" "${LOGS_DIR}"

# Check dependency from step 08
if [ ! -f "${RESULTS_DIR}/step_completed_ranking.txt" ]; then
    log "Error: Structural analysis step not completed."
    exit 1
fi

log "Starting report generation."

# Generate LaTeX report template
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
This report summarizes the identification and characterization of chemoreceptive GPCRs in \textit{Berghia stephanieae}, focusing on phylogenetic, evolutionary, and structural analyses.

\section{Phylogenetic Analysis}
Phylogenetic trees were constructed using IQ-TREE and Phyloformer (Fig. \ref{fig:phylo_comparison}). Bayesian inference with MrBayes was optionally applied.

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{../phylogenies/visualizations/all_berghia_refs_comparison.png}
\caption{Phyloformer vs. IQ-TREE Tree Comparison}
\label{fig:phylo_comparison}
\end{figure}

\section{Selective Pressure and ASR}
Selective pressure was analyzed with PAML (Fig. \ref{fig:selective_pressure}). ASR was performed for key nodes (Fig. \ref{fig:asr}).

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{../selective_pressure/lrt_plot.png}
\caption{Selective Pressure Analysis}
\label{fig:selective_pressure}
\end{figure}

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{../asr/asr_plot.png}
\caption{Ancestral Sequence Reconstruction}
\label{fig:asr}
\end{figure}

\section{Synteny Analysis}
Synteny blocks were identified across genomes (Fig. \ref{fig:synteny}).

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{../synteny/synteny_plot.png}
\caption{Synteny Blocks}
\label{fig:synteny}
\end{figure}

\section{Candidate Ranking}
Top GPCR candidates were ranked (Fig. \ref{fig:ranking}).

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{../ranking/ranking_plot.png}
\caption{Top Ranked Candidates}
\label{fig:ranking}
\end{figure}

\section{Structural Analysis}
Structural phylogenies were compared with sequence-based trees (Fig. \ref{fig:struct_vs_seq}, \ref{fig:clustered_vs_seq}). TM-align scores were analyzed (Fig. \ref{fig:tmalign_heatmap}, \ref{fig:tmalign_pca}).

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{../structural_analysis/foldtree_vs_seq_plot.png}
\caption{FoldTree vs. Sequence Tree Comparison}
\label{fig:struct_vs_seq}
\end{figure}

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{../structural_analysis/clustered_vs_seq_plot.png}
\caption{UPGMA vs. Sequence Tree Comparison}
\label{fig:clustered_vs_seq}
\end{figure}

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{../structural_analysis/tmalign_heatmap.png}
\caption{TM-align Score Heatmap}
\label{fig:tmalign_heatmap}
\end{figure}

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{../structural_analysis/tmalign_pca.png}
\caption{PCA of TM-align Scores}
\label{fig:tmalign_pca}
\end{figure}

\end{document}
EOF

# Compile LaTeX report
run_command "latex_compile" ${PDFLATEX} -output-directory="${RESULTS_DIR}/report" "${RESULTS_DIR}/report/report.tex"

# Move compiled PDF to results directory
mv "${RESULTS_DIR}/report/report.pdf" "${RESULTS_DIR}/pipeline_report.pdf"

log "Report generation completed."
