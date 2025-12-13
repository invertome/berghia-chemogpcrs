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

# --- Helper function to find first available image ---
find_image() {
    local dir="$1"
    local pattern="$2"
    local fallback="$3"

    # Find first matching file
    for f in "${dir}"/${pattern}; do
        if [ -f "$f" ]; then
            echo "$f"
            return 0
        fi
    done

    # Return fallback or empty
    echo "$fallback"
    return 1
}

# --- Find available images dynamically ---
PHYLO_VIZ_DIR="${RESULTS_DIR}/phylogenies/visualizations"
ASR_DIR="${RESULTS_DIR}/asr"
SYNTENY_DIR="${RESULTS_DIR}/synteny"
RANKING_DIR="${RESULTS_DIR}/ranking"
STRUCT_DIR="${RESULTS_DIR}/structural_analysis"

# Find phylogenetic tree image (prefer colored, fall back to basic)
PHYLO_IMG=$(find_image "$PHYLO_VIZ_DIR" "all_berghia_refs_colored.png" "")
[ -z "$PHYLO_IMG" ] && PHYLO_IMG=$(find_image "$PHYLO_VIZ_DIR" "all_berghia_refs_basic.png" "")
[ -z "$PHYLO_IMG" ] && PHYLO_IMG=$(find_image "$PHYLO_VIZ_DIR" "*_basic.png" "")

# Find ASR plot
ASR_IMG=$(find_image "$ASR_DIR" "*_asr_plot_circular.png" "")
[ -z "$ASR_IMG" ] && ASR_IMG=$(find_image "$ASR_DIR" "*_asr_plot.png" "")

# Find synteny plot
SYNTENY_IMG=$(find_image "$SYNTENY_DIR" "synteny_plot_*.png" "")

# Find ranking plot
RANKING_IMG=$(find_image "$RANKING_DIR" "ranking_plot.png" "")
[ -z "$RANKING_IMG" ] && RANKING_IMG=$(find_image "$RANKING_DIR" "*.png" "")

# Find structural comparison plot
STRUCT_IMG=$(find_image "$STRUCT_DIR" "struct_vs_seq_plot.png" "")
[ -z "$STRUCT_IMG" ] && STRUCT_IMG=$(find_image "$STRUCT_DIR" "*.png" "")

# --- Generate LaTeX report ---
cat > "${RESULTS_DIR}/report/report.tex" <<'LATEX_HEADER'
\documentclass{article}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{geometry}
\usepackage{hyperref}
\geometry{a4paper, margin=1in}
\title{GPCR Chemoreceptor Analysis in Berghia stephanieae}
\author{Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst}
\date{\today}

\begin{document}

\maketitle

\section{Introduction}
This report presents the results of a comprehensive analysis of chemoreceptive GPCRs in \textit{Berghia stephanieae}, integrating phylogeny, selective pressure, synteny, and structural data.

LATEX_HEADER

# Add phylogenetic section if image exists
if [ -n "$PHYLO_IMG" ] && [ -f "$PHYLO_IMG" ]; then
    # Copy image to report directory
    cp "$PHYLO_IMG" "${RESULTS_DIR}/report/"
    PHYLO_BASENAME=$(basename "$PHYLO_IMG")
    cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\subsection{Phylogenetic Analysis}
Trees were constructed using IQ-TREE (with FastTree seed strategy) and Phyloformer, with optional MrBayes analysis.

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{${PHYLO_BASENAME}}
\caption{Phylogeny of Berghia candidates and references.}
\label{fig:phylo}
\end{figure}

EOF
else
    cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\subsection{Phylogenetic Analysis}
Trees were constructed using IQ-TREE (with FastTree seed strategy) and Phyloformer.
\textit{Note: Phylogenetic visualization not available.}

EOF
fi

# Add ASR section if image exists
if [ -n "$ASR_IMG" ] && [ -f "$ASR_IMG" ]; then
    cp "$ASR_IMG" "${RESULTS_DIR}/report/"
    ASR_BASENAME=$(basename "$ASR_IMG")
    cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\subsection{Selective Pressure and ASR}
Selective pressure was analyzed with HyPhy aBSREL, and ancestral sequences were reconstructed with FastML.

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{${ASR_BASENAME}}
\caption{Ancestral sequence reconstruction visualization.}
\label{fig:asr}
\end{figure}

EOF
else
    cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\subsection{Selective Pressure and ASR}
Selective pressure was analyzed with HyPhy aBSREL, and ancestral sequences were reconstructed with FastML.
\textit{Note: ASR visualization not available.}

EOF
fi

# Add synteny section if image exists
if [ -n "$SYNTENY_IMG" ] && [ -f "$SYNTENY_IMG" ]; then
    cp "$SYNTENY_IMG" "${RESULTS_DIR}/report/"
    SYNTENY_BASENAME=$(basename "$SYNTENY_IMG")
    cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\subsection{Synteny Analysis}
Synteny blocks were identified using MCScanX across available genomes.

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{${SYNTENY_BASENAME}}
\caption{Synteny analysis.}
\label{fig:synteny}
\end{figure}

EOF
else
    cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\subsection{Synteny Analysis}
Synteny blocks were identified using MCScanX across available genomes.
\textit{Note: Synteny visualization not available.}

EOF
fi

# Add ranking section if image exists
if [ -n "$RANKING_IMG" ] && [ -f "$RANKING_IMG" ]; then
    cp "$RANKING_IMG" "${RESULTS_DIR}/report/"
    RANKING_BASENAME=$(basename "$RANKING_IMG")
    cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\subsection{Candidate Ranking}
Candidates were ranked based on integrated data including phylogenetic distance, dN/dS ratios, expression levels, and synteny conservation.

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{${RANKING_BASENAME}}
\caption{Top-ranked GPCR candidates.}
\label{fig:ranking}
\end{figure}

EOF
else
    cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\subsection{Candidate Ranking}
Candidates were ranked based on integrated data including phylogenetic distance, dN/dS ratios, expression levels, and synteny conservation.
\textit{Note: Ranking visualization not available.}

EOF
fi

# Add structural section if image exists
if [ -n "$STRUCT_IMG" ] && [ -f "$STRUCT_IMG" ]; then
    cp "$STRUCT_IMG" "${RESULTS_DIR}/report/"
    STRUCT_BASENAME=$(basename "$STRUCT_IMG")
    cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\subsection{Structural Analysis}
Structural phylogenies included AlphaFold predictions and GPCRdb references.

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{${STRUCT_BASENAME}}
\caption{Comparison of structural and sequence phylogenies.}
\label{fig:struct}
\end{figure}

EOF
else
    cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\subsection{Structural Analysis}
Structural phylogenies included AlphaFold predictions and GPCRdb references.
\textit{Note: Structural visualization not available.}

EOF
fi

# Close LaTeX document
cat >> "${RESULTS_DIR}/report/report.tex" <<'LATEX_FOOTER'

\section{Methods Summary}
\begin{itemize}
\item \textbf{GPCR Identification}: HHblits + DeepTMHMM (6+ TM regions filter)
\item \textbf{Orthology}: OrthoFinder clustering
\item \textbf{Phylogenetics}: FastTree seed + IQ-TREE with 1000 bootstrap replicates
\item \textbf{Selective Pressure}: HyPhy aBSREL
\item \textbf{ASR}: FastML on deep nodes (top 10\% by distance)
\item \textbf{Synteny}: MCScanX
\item \textbf{Structure}: AlphaFold predictions
\end{itemize}

\end{document}
LATEX_FOOTER

# Compile LaTeX (run twice for references)
cd "${RESULTS_DIR}/report" || exit 1
run_command "latex_compile_1" ${PDFLATEX} -interaction=nonstopmode report.tex || log "Warning: First LaTeX pass had issues"
run_command "latex_compile_2" ${PDFLATEX} -interaction=nonstopmode report.tex || { log "Error: LaTeX compilation failed"; exit 1; }

# Move PDF
if [ -f "${RESULTS_DIR}/report/report.pdf" ]; then
    mv "${RESULTS_DIR}/report/report.pdf" "${RESULTS_DIR}/pipeline_report.pdf" || { log "Error: Failed to move PDF"; exit 1; }
    log "Report generated: ${RESULTS_DIR}/pipeline_report.pdf"
else
    log "Error: PDF was not generated"
    exit 1
fi

touch "${RESULTS_DIR}/step_completed_report.txt"
log "Report generation completed."
