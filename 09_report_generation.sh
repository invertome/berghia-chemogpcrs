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
# shellcheck source=scripts/_report_summary_lib.sh
source "${SCRIPTS_DIR}/_report_summary_lib.sh"

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

    for f in "${dir}"/${pattern}; do
        if [ -f "$f" ]; then
            echo "$f"
            return 0
        fi
    done
    echo "$fallback"
    return 1
}

# --- Helper function to count lines in file ---
count_lines() {
    local file="$1"
    if [ -f "$file" ]; then
        wc -l < "$file" | tr -d ' '
    else
        echo "0"
    fi
}

# Header-aware CSV helpers (csv_col_idx, count_by_value, sum_column, read_csv_value)
# come from scripts/_report_summary_lib.sh sourced above.

# --- Find directories ---
PHYLO_VIZ_DIR="${RESULTS_DIR}/phylogenies/visualizations"
ASR_DIR="${RESULTS_DIR}/asr"
SYNTENY_DIR="${RESULTS_DIR}/synteny"
RANKING_DIR="${RESULTS_DIR}/ranking"
STRUCT_DIR="${RESULTS_DIR}/structural_analysis"
SEL_PRESS_DIR="${RESULTS_DIR}/selective_pressure"

# --- Find available images dynamically ---
PHYLO_IMG=$(find_image "$PHYLO_VIZ_DIR" "all_berghia_refs_colored.png" "")
[ -z "$PHYLO_IMG" ] && PHYLO_IMG=$(find_image "$PHYLO_VIZ_DIR" "all_berghia_refs_basic.png" "")
[ -z "$PHYLO_IMG" ] && PHYLO_IMG=$(find_image "$PHYLO_VIZ_DIR" "*_basic.png" "")

ASR_IMG=$(find_image "$ASR_DIR" "*_asr_plot_circular.png" "")
[ -z "$ASR_IMG" ] && ASR_IMG=$(find_image "$ASR_DIR" "*_asr_plot.png" "")

SYNTENY_IMG=$(find_image "$SYNTENY_DIR" "synteny_plot_*.png" "")

RANKING_IMG=$(find_image "$RANKING_DIR" "ranking_plot.png" "")
[ -z "$RANKING_IMG" ] && RANKING_IMG=$(find_image "$RANKING_DIR" "*.png" "")

STRUCT_IMG=$(find_image "$STRUCT_DIR" "struct_vs_seq_plot.png" "")
[ -z "$STRUCT_IMG" ] && STRUCT_IMG=$(find_image "$STRUCT_DIR" "*.png" "")

SEL_PRESS_IMG=$(find_image "$SEL_PRESS_DIR" "selective_pressure_plot.png" "")

HEATMAP_IMG=$(find_image "$STRUCT_DIR" "tmalign_heatmap.png" "")
PCA_IMG=$(find_image "$STRUCT_DIR" "tmalign_pca.png" "")

# Tandem-cluster plot (May-2026 stack; bead -ar8). Emitted by run_tandem_detection.sh.
TANDEM_IMG=$(find_image "$SYNTENY_DIR" "tandem_clusters_plot.png" "")

# JCVI synteny dotplot (May-2026 stack; bead -e59). One PDF per Berghia-vs-target run
# under ${SYNTENY_DIR}/jcvi/<run>/<a>.<b>.pdf. The pattern *.*.pdf (a basename
# with at least one dot) matches JCVI's <a>.<b>.pdf dotplots while excluding
# auxiliary single-token PDFs (karyotype.pdf, density.pdf, etc.) the catalog
# step also emits.
JCVI_IMG=$(find_image "$SYNTENY_DIR/jcvi" "*/*.*.pdf" "")

# --- Collect statistics (header-aware; tolerates schema additions) ---
TOTAL_CANDIDATES=$(count_lines "${RESULTS_DIR}/candidates/chemogpcr_candidates.txt")
RANKED_FILE="${RANKING_DIR}/ranked_candidates_sorted.csv"

HIGH_CONF=""; MED_CONF=""; LOW_CONF=""; SIG_SELECTION=""
BUSTED_S_SIG=""; BUSTED_MH_SIG=""; MEME_TOTAL=""
HCR_FRIENDLY=""; TANDEM_MEMBERS=""
CDS_NATIVE=""; CDS_MINIPROT=""
TOP_N_RANKED=""

if [ -f "$RANKED_FILE" ]; then
    HIGH_CONF=$(count_by_value      "$RANKED_FILE" confidence_tier        High)
    MED_CONF=$(count_by_value       "$RANKED_FILE" confidence_tier        Medium)
    LOW_CONF=$(count_by_value       "$RANKED_FILE" confidence_tier        Low)
    SIG_SELECTION=$(count_by_value  "$RANKED_FILE" selection_significant  True)
    BUSTED_S_SIG=$(count_by_value   "$RANKED_FILE" busted_s_significant   True)
    BUSTED_MH_SIG=$(count_by_value  "$RANKED_FILE" busted_mh_significant  True)
    MEME_TOTAL=$(sum_column         "$RANKED_FILE" meme_n_episodic_sites)
    HCR_FRIENDLY=$(count_by_value   "$RANKED_FILE" hcr_probe_friendly     True)
    TANDEM_MEMBERS=$(count_by_value "$RANKED_FILE" has_tandem_cluster_data True)
    CDS_NATIVE=$(count_by_value     "$RANKED_FILE" cds_source             native)
    CDS_MINIPROT=$(count_by_value   "$RANKED_FILE" cds_source             miniprot)
    TOP_N_RANKED=$(awk -F',' 'NR>1{n++}END{print n+0}' "$RANKED_FILE")
fi

# --- Generate LaTeX report ---
cat > "${RESULTS_DIR}/report/report.tex" <<'LATEX_HEADER'
\documentclass[11pt]{article}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{geometry}
\usepackage{hyperref}
\usepackage{xcolor}
\usepackage{float}
\usepackage{caption}
\usepackage{subcaption}
\geometry{a4paper, margin=1in}

% Define colors for confidence tiers
\definecolor{highconf}{RGB}{46, 160, 67}
\definecolor{medconf}{RGB}{255, 127, 14}
\definecolor{lowconf}{RGB}{214, 39, 40}

\title{GPCR Chemoreceptor Analysis in \textit{Berghia stephanieae}\\
\large Comprehensive Pipeline Report}
\author{Generated by ChemoGPCR Discovery Pipeline\\
Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst}
\date{\today}

\begin{document}

\maketitle

\begin{abstract}
This report presents results from a comprehensive computational pipeline for identifying and characterizing chemoreceptive G protein-coupled receptors (GPCRs) in \textit{Berghia stephanieae}. The analysis integrates phylogenetic inference, selective pressure analysis, synteny mapping, and structural prediction to rank candidate chemoreceptors by their likelihood of functional relevance.
\end{abstract}

\tableofcontents
\newpage

\section{Executive Summary}

LATEX_HEADER

# Add summary statistics
cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\subsection{Key Findings}

\begin{table}[H]
\centering
\caption{Pipeline Summary Statistics (May 2026 schema)}
\begin{tabular}{lr}
\toprule
\textbf{Metric} & \textbf{Value} \\\\
\midrule
Total GPCR Candidates Identified & ${TOTAL_CANDIDATES:-N/A} \\\\
Ranked Candidates & ${TOP_N_RANKED:-N/A} \\\\
\midrule
High Confidence & \textcolor{highconf}{${HIGH_CONF:-0}} \\\\
Medium Confidence & \textcolor{medconf}{${MED_CONF:-0}} \\\\
Low Confidence & \textcolor{lowconf}{${LOW_CONF:-0}} \\\\
\midrule
aBSREL branch-significant & ${SIG_SELECTION:-0} \\\\
BUSTED-S gene-significant & ${BUSTED_S_SIG:-0} \\\\
BUSTED-MH gene-significant & ${BUSTED_MH_SIG:-0} \\\\
MEME episodic sites (sum) & ${MEME_TOTAL:-0} \\\\
\midrule
HCR probe-friendly & ${HCR_FRIENDLY:-0} \\\\
In tandem cluster & ${TANDEM_MEMBERS:-0} \\\\
\midrule
CDS source: native & ${CDS_NATIVE:-0} \\\\
CDS source: miniprot-recovered & ${CDS_MINIPROT:-0} \\\\
\bottomrule
\end{tabular}
\end{table}

EOF

# Add interpretation guide
cat >> "${RESULTS_DIR}/report/report.tex" <<'EOF'

\subsection{Confidence Tier Interpretation}

\begin{itemize}
    \item \textbf{\textcolor{highconf}{High Confidence}}: Strong evidence across multiple data sources (phylogenetic proximity to known chemoreceptors, significant selective pressure, synteny conservation, expression in sensory tissues).
    \item \textbf{\textcolor{medconf}{Medium Confidence}}: Good evidence but incomplete data coverage or moderate scores.
    \item \textbf{\textcolor{lowconf}{Low Confidence}}: Limited evidence or low scores; requires additional validation.
\end{itemize}

\newpage
\section{Methods Overview}

\subsection{Pipeline Architecture}

The analysis pipeline consists of these sequential stages (May 2026 stack):

\begin{enumerate}
    \item \textbf{Reference Processing}: Build per-orthogroup HMMs from Nath et al. references.
    \item \textbf{GPCR Identification}: HHblits homology + TMbed transmembrane prediction (DeepTMHMM / Phobius / Kyte--Doolittle as graded fallbacks; 6+ TM filter). Predictor provenance recorded per run.
    \item \textbf{Orthology Clustering}: OrthoFinder (Diamond mode).
    \item \textbf{Phylogenetic Analysis}: regime-based MAFFT/FAMSA alignment $\to$ HmmCleaner (segment) $\to$ ClipKit (column) $\to$ FastTree seed $\to$ IQ-TREE 3 (ModelFinder; UFBoot+SH-aLRT+TBE; deterministic seed) $\to$ TreeShrink rogue-taxon removal.
    \item \textbf{Selective Pressure}: HyPhy stack — GARD recombination breakpoints $\to$ BUSTED-S (gene-wide) $\to$ BUSTED-MH (multi-hit-aware) $\to$ aBSREL (branch-specific) $\to$ MEME (episodic site-level). MACSE v2 for codon alignment.
    \item \textbf{Ancestral Reconstruction}: IQ-TREE 3 \texttt{-{}-ancestral} (FastML fallback).
    \item \textbf{Reconciliation}: GeneRax (Notung fallback).
    \item \textbf{Synteny}: JCVI MCscan (MCScanX fallback); per-gene anchor counts feed the ranking.
    \item \textbf{Tandem-cluster Detection}: sliding-window clustering of sibling paralogs --- the field's signature chemoreceptor signal.
    \item \textbf{Candidate Ranking}: Multi-criteria weighted scoring (see Section~\ref{sec:weights}).
    \item \textbf{HCR Augmentation}: per-candidate \texttt{cds\_length\_bp}, closest-paralog identity, and \texttt{hcr\_probe\_friendly} boolean for wet-lab probe design.
    \item \textbf{Structural Analysis}: AlphaFold prediction + Foldseek vs.\ GPCRdb (P3).
\end{enumerate}

\subsection{Ranking Algorithm}

Candidates are scored using a weighted combination of:

\begin{itemize}
    \item \textbf{Phylogenetic proximity} to known chemoreceptors (inverse weighted distance).
    \item \textbf{Positive selection}: \texttt{omega\_max} from aBSREL (was weighted-mean prior to May 2026 --- fixed to preserve the episodic-positive-selection signal).
    \item \textbf{Purifying selection}: \emph{disabled by default} (\texttt{PURIFYING\_WEIGHT=0} in \texttt{config.sh}) --- chemoreceptor identification rewards diversifying, not housekeeping conservation.
    \item \textbf{Synteny}: JCVI per-gene anchor count (replaces the broken ID-count input).
    \item \textbf{Tandem cluster}: sibling-paralog clustering signal (chemoreceptor-defining).
    \item \textbf{LSE depth}: position in lineage-specific expansion clades (75\textsuperscript{th}-percentile threshold).
    \item \textbf{Chemosensory expression}: tissue-specific signal where available.
    \item \textbf{G-protein coexpression}: G$\alpha$ family coexpression with the candidate (orthogroup-level).
    \item \textbf{ECL divergence}: extracellular-loop sequence divergence vs.\ paralogs.
    \item \textbf{Family expansion} and \textbf{orthogroup confidence} (auxiliary).
\end{itemize}

P-values from BUSTED-S, BUSTED-MH, aBSREL, and MEME are FDR-corrected via \texttt{statsmodels.multipletests} (Benjamini--Hochberg, $\alpha=0.05$). The composite score is multiplied by an \texttt{evidence\_completeness} factor (floor 0.4) so sparse-data candidates do not co-rank with complete-data candidates.

EOF

# Add phylogenetic section
cat >> "${RESULTS_DIR}/report/report.tex" <<'EOF'

\newpage
\section{Results}

\subsection{Phylogenetic Analysis}

Phylogenetic trees were constructed using a two-stage approach: FastTree for initial topology estimation followed by IQ-TREE for maximum likelihood optimization with 1000 ultrafast bootstrap replicates.

EOF

if [ -n "$PHYLO_IMG" ] && [ -f "$PHYLO_IMG" ]; then
    cp "$PHYLO_IMG" "${RESULTS_DIR}/report/"
    PHYLO_BASENAME=$(basename "$PHYLO_IMG")
    cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\begin{figure}[H]
\centering
\includegraphics[width=0.9\textwidth]{${PHYLO_BASENAME}}
\caption{Phylogenetic tree of \textit{Berghia} GPCR candidates and reference sequences. Colors indicate sequence type: blue = reference GPCRs, green = \textit{Berghia} candidates, orange = other species.}
\label{fig:phylo}
\end{figure}

EOF
else
    cat >> "${RESULTS_DIR}/report/report.tex" <<'EOF'

\textit{Note: Phylogenetic visualization not available.}

EOF
fi

# Add selective pressure section
cat >> "${RESULTS_DIR}/report/report.tex" <<'EOF'

\subsection{Selective Pressure Analysis}

Branch-specific selection was analyzed using HyPhy aBSREL to identify:
\begin{itemize}
    \item \textbf{Purifying selection} ($\omega < 1$): Indicates conserved function under evolutionary constraint
    \item \textbf{Positive selection} ($\omega > 1$): Suggests adaptive evolution or functional diversification
    \item \textbf{Neutral evolution} ($\omega \approx 1$): May indicate relaxed constraint or pseudogenization
\end{itemize}

P-values were corrected using Benjamini-Hochberg FDR at $\alpha = 0.05$.

EOF

if [ -n "$SEL_PRESS_IMG" ] && [ -f "$SEL_PRESS_IMG" ]; then
    cp "$SEL_PRESS_IMG" "${RESULTS_DIR}/report/"
    SEL_BASENAME=$(basename "$SEL_PRESS_IMG")
    cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\begin{figure}[H]
\centering
\includegraphics[width=0.9\textwidth]{${SEL_BASENAME}}
\caption{Selective pressure analysis. A) Volcano plot showing log($\omega$) vs. statistical significance. B) Distribution of selection coefficients. C) Classification summary. D) P-value calibration.}
\label{fig:selection}
\end{figure}

EOF
fi

# Add ASR section
if [ -n "$ASR_IMG" ] && [ -f "$ASR_IMG" ]; then
    cp "$ASR_IMG" "${RESULTS_DIR}/report/"
    ASR_BASENAME=$(basename "$ASR_IMG")
    cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\subsection{Ancestral Sequence Reconstruction}

Deep ancestral nodes (top 10\% by tree distance) were reconstructed using FastML to identify conserved sequence features.

\begin{figure}[H]
\centering
\includegraphics[width=0.8\textwidth]{${ASR_BASENAME}}
\caption{Ancestral sequence reconstruction visualization. Red nodes indicate reconstructed ancestral sequences.}
\label{fig:asr}
\end{figure}

EOF
fi

# Add synteny section (legacy MCScanX visualization)
if [ -n "$SYNTENY_IMG" ] && [ -f "$SYNTENY_IMG" ]; then
    cp "$SYNTENY_IMG" "${RESULTS_DIR}/report/"
    SYNTENY_BASENAME=$(basename "$SYNTENY_IMG")
    cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\subsection{Synteny Analysis}

Conserved gene arrangements were identified using MCScanX (legacy visualization) and JCVI MCscan (May 2026 default) to provide additional evidence for orthology relationships.

\begin{figure}[H]
\centering
\includegraphics[width=0.8\textwidth]{${SYNTENY_BASENAME}}
\caption{Synteny analysis showing conserved gene blocks across species (MCScanX collinearity visualization).}
\label{fig:synteny}
\end{figure}

EOF
fi

# JCVI MCscan dotplot (May-2026 default backend; bead -e59).
if [ -n "$JCVI_IMG" ] && [ -f "$JCVI_IMG" ]; then
    cp "$JCVI_IMG" "${RESULTS_DIR}/report/"
    JCVI_BASENAME=$(basename "$JCVI_IMG")
    cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\begin{figure}[H]
\centering
\includegraphics[width=0.7\textwidth]{${JCVI_BASENAME}}
\caption{JCVI MCscan dotplot (Berghia vs.\ a representative target genome). Diagonals indicate conserved synteny blocks; broken / scattered points indicate rearrangements. Per-gene anchor counts derived from these alignments feed the synteny axis of the ranking algorithm.}
\label{fig:jcvi_dotplot}
\end{figure}

EOF
fi

# Tandem-cluster figure (May-2026 stack; bead -ar8).
if [ -n "$TANDEM_IMG" ] && [ -f "$TANDEM_IMG" ]; then
    cp "$TANDEM_IMG" "${RESULTS_DIR}/report/"
    TANDEM_BASENAME=$(basename "$TANDEM_IMG")
    cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\subsection{Tandem-Cluster Detection}

Intra-genome tandem-cluster detection scans the Berghia GFF for sliding-window clusters (default 100\,kb gap, $\geq$3 members) of paralogous chemoreceptor candidates. Tandem clustering is the field's signature chemoreceptor signal.

\begin{figure}[H]
\centering
\includegraphics[width=0.95\textwidth]{${TANDEM_BASENAME}}
\caption{Berghia GPCR tandem-cluster summary. Panel A: distribution of tandem-cluster sizes across all candidates (singletons through largest). Panel B: top distinct tandem clusters by size.}
\label{fig:tandem}
\end{figure}

EOF
fi

# Add ranking section
cat >> "${RESULTS_DIR}/report/report.tex" <<'EOF'

\subsection{Candidate Ranking}

Candidates were ranked using an integrated multi-criteria scoring system that combines phylogenetic, selective pressure, synteny, and expression evidence.

EOF

if [ -n "$RANKING_IMG" ] && [ -f "$RANKING_IMG" ]; then
    cp "$RANKING_IMG" "${RESULTS_DIR}/report/"
    RANKING_BASENAME=$(basename "$RANKING_IMG")
    cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\begin{figure}[H]
\centering
\includegraphics[width=0.95\textwidth]{${RANKING_BASENAME}}
\caption{Candidate ranking summary. A) Top candidates with score breakdown by component. B) Score distributions. C) Confidence tier distribution. D) Evidence completeness correlation. E) Selection pressure summary.}
\label{fig:ranking}
\end{figure}

EOF
fi

# --- HCR probe-friendliness explainer (always emit; columns may be 0 if not run) ---
cat >> "${RESULTS_DIR}/report/report.tex" <<'EOF'

\subsubsection{HCR Probe-Friendliness}

For each candidate, \texttt{hcr\_probe\_friendly} is a boolean derived from CDS length and the closest-paralog identity (windowed). True iff:
\begin{itemize}
  \item \texttt{cds\_length\_bp} $\geq$ 600 bp (HCR\_DEFAULT\_MIN\_CDS\_LENGTH\_BP)
  \item \texttt{paralog\_min\_identity} $\leq$ 0.80 (HCR\_DEFAULT\_MAX\_PARALOG\_IDENTITY)
\end{itemize}
Candidates marked True are recommended for first-round HCR design without further filtering. False candidates are not disqualified, but require either probe-design effort against a divergent region or a discriminating split-probe strategy.

EOF

# --- Top-N ranked candidates table (header-aware) ---
if [ -f "$RANKED_FILE" ]; then
    cat >> "${RESULTS_DIR}/report/report.tex" <<'EOF'

\subsubsection{Top Ranked Candidates}

\begin{longtable}{rlrlcccr}
\caption{Top 20 ranked GPCR candidates. \textbf{Sel}=aBSREL branch-significant; \textbf{B-MH}=BUSTED-MH gene-significant; \textbf{HCR}=probe-friendly; \textbf{TC size}=tandem-cluster size.} \\
\toprule
\textbf{\#} & \textbf{ID} & \textbf{Score} & \textbf{Conf} & \textbf{Sel} & \textbf{B-MH} & \textbf{HCR} & \textbf{TC size} \\
\midrule
\endfirsthead
\multicolumn{8}{c}{\tablename\ \thetable\ -- continued} \\
\toprule
\textbf{\#} & \textbf{ID} & \textbf{Score} & \textbf{Conf} & \textbf{Sel} & \textbf{B-MH} & \textbf{HCR} & \textbf{TC size} \\
\midrule
\endhead
\bottomrule
\endfoot
EOF

    python3 - "$RANKED_FILE" >> "${RESULTS_DIR}/report/report.tex" <<'PY'
import csv, re, sys

# Escape LaTeX-special characters in interpolated string fields. & % $ # _ { }
# need a backslash; \ ^ ~ need their package-supplied math-mode forms. Without
# this, an ID containing any of these (NCBI provisional names occasionally
# include "&"; gene_name fields routinely contain "_") aborts pdflatex.
_LATEX_SPECIALS = re.compile(r'([&%$#_{}])')
def tex_esc(s):
    s = "" if s is None else str(s)
    s = _LATEX_SPECIALS.sub(r'\\\1', s)
    s = s.replace('\\', r'\textbackslash{}').replace('^', r'\^{}').replace('~', r'\~{}')
    return s

def yn(v):
    return "Y" if str(v).strip().lower() in ("true", "1", "yes", "y") else "N"

with open(sys.argv[1]) as fh:
    r = csv.DictReader(fh)
    for i, row in enumerate(r, start=1):
        if i > 20:
            break
        conf = tex_esc(row.get("confidence_tier", ""))
        color = {"High": "highconf", "Medium": "medconf", "Low": "lowconf"}.get(
            row.get("confidence_tier", ""), "lowconf")
        try:
            score = float(row.get("rank_score", "0") or 0)
        except ValueError:
            score = 0.0
        tc = tex_esc(row.get("tandem_cluster_size", "0") or "0")
        rid = tex_esc(row.get("id", ""))
        print(rf"{i} & {rid} & {score:.3f} & "
              rf"\textcolor{{{color}}}{{{conf}}} & "
              rf"{yn(row.get('selection_significant'))} & "
              rf"{yn(row.get('busted_mh_significant'))} & "
              rf"{yn(row.get('hcr_probe_friendly'))} & {tc} \\")
PY

    cat >> "${RESULTS_DIR}/report/report.tex" <<'EOF'
\end{longtable}

EOF
fi

# --- Positive-control sanity check (HCR-validated markers e.g. Galpha_olf) ---
POSCTRL_TSV="${RANKING_DIR}/positive_controls_check.tsv"
if [ -f "$POSCTRL_TSV" ]; then
    cat >> "${RESULTS_DIR}/report/report.tex" <<'EOF'

\subsection{Positive-Control Sanity Check}

HCR-validated chemoreceptor controls (e.g.\ G$\alpha_{\text{olf}}$) are scored by the same pipeline. A control falling below the 50\textsuperscript{th} percentile is a red flag for the ranking weights, not the control. \textcolor{lowconf}{ALERT} rows indicate either a missing control or one below the alert percentile.

\begin{table}[H]
\centering
\caption{Positive controls — pipeline rank, percentile, and alert status.}
\begin{tabular}{llrrl}
\toprule
\textbf{Control} & \textbf{Found} & \textbf{Rank} & \textbf{Percentile} & \textbf{Status} \\
\midrule
EOF
    python3 - "$POSCTRL_TSV" >> "${RESULTS_DIR}/report/report.tex" <<'PY'
import csv, re, sys

_LATEX_SPECIALS = re.compile(r'([&%$#_{}])')
def tex_esc(s):
    s = "" if s is None else str(s)
    s = _LATEX_SPECIALS.sub(r'\\\1', s)
    s = s.replace('\\', r'\textbackslash{}').replace('^', r'\^{}').replace('~', r'\~{}')
    return s

def truthy(v):
    return str(v).strip().lower() in ("true", "1", "yes", "y")

with open(sys.argv[1], newline="") as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        name = tex_esc(row.get("gene_name") or "")
        found = "Y" if truthy(row.get("found")) else "N"
        rank = tex_esc(row.get("rank", "") or "-")
        pct_raw = row.get("percentile", "")
        try:
            pct = f"{float(pct_raw):.1f}" if pct_raw not in ("", None) else "-"
        except ValueError:
            pct = "-"
        alert = truthy(row.get("alert"))
        status = (r"\textcolor{lowconf}{ALERT}" if alert
                  else r"\textcolor{highconf}{OK}")
        print(rf"{name} & {found} & {rank} & {pct} & {status} \\")
PY
    cat >> "${RESULTS_DIR}/report/report.tex" <<'EOF'
\bottomrule
\end{tabular}
\end{table}

EOF
fi

# Add structural analysis section
cat >> "${RESULTS_DIR}/report/report.tex" <<'EOF'

\subsection{Structural Analysis}

Protein structures were predicted using AlphaFold and compared using TM-align to assess structural similarity to known GPCRs.

EOF

if [ -n "$STRUCT_IMG" ] && [ -f "$STRUCT_IMG" ]; then
    cp "$STRUCT_IMG" "${RESULTS_DIR}/report/"
    STRUCT_BASENAME=$(basename "$STRUCT_IMG")
    cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\begin{figure}[H]
\centering
\includegraphics[width=0.8\textwidth]{${STRUCT_BASENAME}}
\caption{Structural vs. sequence phylogeny comparison.}
\label{fig:struct}
\end{figure}

EOF
fi

if [ -n "$HEATMAP_IMG" ] && [ -f "$HEATMAP_IMG" ]; then
    cp "$HEATMAP_IMG" "${RESULTS_DIR}/report/"
    HEATMAP_BASENAME=$(basename "$HEATMAP_IMG")
    cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\begin{figure}[H]
\centering
\includegraphics[width=0.9\textwidth]{${HEATMAP_BASENAME}}
\caption{Hierarchically clustered heatmap of TM-align structural similarity scores. TM-score $>$ 0.5 indicates same fold; $>$ 0.17 indicates similar topology.}
\label{fig:heatmap}
\end{figure}

EOF
fi

if [ -n "$PCA_IMG" ] && [ -f "$PCA_IMG" ]; then
    cp "$PCA_IMG" "${RESULTS_DIR}/report/"
    PCA_BASENAME=$(basename "$PCA_IMG")
    cat >> "${RESULTS_DIR}/report/report.tex" <<EOF

\begin{figure}[H]
\centering
\includegraphics[width=0.9\textwidth]{${PCA_BASENAME}}
\caption{PCA of structural similarity showing clustering of GPCR candidates.}
\label{fig:pca}
\end{figure}

EOF
fi

# Close document
cat >> "${RESULTS_DIR}/report/report.tex" <<'LATEX_FOOTER'

\newpage
\section{Methods Details}

\subsection{Software Versions}

\begin{table}[H]
\centering
\caption{Key Software Tools (May 2026)}
\begin{tabular}{ll}
\toprule
\textbf{Tool} & \textbf{Purpose} \\
\midrule
HHblits / HMMER & Remote homology / HMM search \\
TMbed (Bernhofer 2022) & Transmembrane topology (primary) \\
DeepTMHMM / Phobius & TM topology (fallbacks) \\
OrthoFinder & Orthology inference \\
MAFFT / FAMSA & Multiple sequence alignment \\
HmmCleaner & Per-sequence segment cleaning \\
ClipKit & Column trimming \\
MACSE v2 & Frameshift-aware codon alignment \\
FastTree & Approximate-ML seed tree \\
IQ-TREE 3 & ML phylogeny + ancestral reconstruction \\
TreeShrink & Rogue-taxon removal \\
HyPhy: GARD / BUSTED-S / BUSTED-MH / aBSREL / MEME & Selection-stack analysis \\
GeneRax & Gene-tree / species-tree reconciliation \\
JCVI MCscan & Synteny / collinearity \\
AlphaFold & Protein structure prediction \\
Foldseek & Structural search vs.\ GPCRdb \\
TM-align & Pairwise structural alignment \\
\bottomrule
\end{tabular}
\end{table}

\subsection{Scoring Weights}\label{sec:weights}

The ranking algorithm uses the following default weights (overridable via environment variables; values shown reflect \texttt{config.sh} defaults at the time of report generation):

\begin{table}[H]
\centering
\caption{Ranking Score Weights}
\begin{tabular}{lcc}
\toprule
\textbf{Component} & \textbf{Default Weight} & \textbf{Variable} \\
\midrule
Phylogenetic proximity & 2.0 & PHYLO\_WEIGHT \\
Positive selection & 1.0 & POSITIVE\_WEIGHT \\
Purifying selection & 0.0 & PURIFYING\_WEIGHT \\
Synteny conservation & 3.0 & SYNTENY\_WEIGHT \\
Tandem cluster & 2.5 & TANDEM\_CLUSTER\_WEIGHT \\
LSE depth & 1.0 & LSE\_DEPTH\_WEIGHT \\
Chemosensory expression & 3.0 & CHEMOSENSORY\_EXPR\_WEIGHT \\
G-protein coexpression & 2.0 & GPROTEIN\_COEXPR\_WEIGHT \\
ECL divergence & 1.5 & ECL\_DIVERGENCE\_WEIGHT \\
Family expansion & 1.5 & EXPANSION\_WEIGHT \\
Orthogroup confidence & 1.0 & OG\_CONFIDENCE\_WEIGHT \\
Expression (legacy) & 1.0 & EXPR\_WEIGHT \\
\bottomrule
\end{tabular}
\end{table}

\section{Recommendations}

\subsection{High-Priority Wet-Lab Shortlist}

For HCR (in situ hybridization chain reaction) probe design, prioritize candidates satisfying \emph{all three}:
\begin{itemize}
    \item \texttt{confidence\_tier = High}
    \item \texttt{hcr\_probe\_friendly = True}
    \item Either \texttt{tandem\_cluster\_size > 1} \emph{or} \texttt{selection\_significant = True} (i.e.\ at least one chemoreceptor-signature corroborator).
\end{itemize}

Secondary candidates (Medium confidence, HCR-friendly, with one signature corroborator) are good targets for second-round screens.

\subsection{Follow-up Analyses}

\begin{itemize}
    \item Expression profiling across head / foot / digestive / sensory tissues --- the current expression-based scores carry the starvation+depth confound noted in the project memory.
    \item Comparative genomics with chromosome-level mollusc assemblies (the synteny weight is highest at 3.0 for a reason).
    \item AlphaFold + Foldseek structural placement against GPCRdb for top-ranked candidates.
    \item Deorphanization / ligand-binding prediction once HCR localization narrows the functional candidate set.
\end{itemize}

\section{Data Availability}

All intermediate files and results are available in the pipeline output directories:
\begin{itemize}
    \item Phylogenies: \texttt{results/phylogenies/}
    \item Selection analysis: \texttt{results/selective\_pressure/}
    \item Rankings: \texttt{results/ranking/}
    \item Structures: \texttt{results/structural\_analysis/}
\end{itemize}

\end{document}
LATEX_FOOTER

# Compile LaTeX (run twice for references and TOC)
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
