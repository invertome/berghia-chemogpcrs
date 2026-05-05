# Stage 09 Report Refresh — Surface May 2026 Pipeline Outputs

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Update `09_report_generation.sh` so the auto-generated PDF reflects the May 2026 pipeline (BUSTED-S/MH, MEME, TMbed, IQ-TREE 3, JCVI MCscan, GeneRax, MACSE, ClipKit, TreeShrink, HCR-friendliness, tandem clusters, CDS provenance, positive-control sanity check).

**Architecture:** The current script is a single bash file that interpolates statistics into a templated LaTeX document via heredocs. The change is a focused in-place edit — no new files, no architectural shift. We add (a) header-aware CSV reads (replace brittle positional `awk $3 / $8` with column-name lookups using a small `csv_col_idx` helper), (b) new aggregate counts (BUSTED-S/MH significant, MEME episodic-site total, HCR-friendly, tandem-clustered, CDS-source breakdown), (c) refreshed Methods / Tools / Weights tables, and (d) a new "HCR Probe-Friendliness" subsection that surfaces the per-candidate HCR columns in the Top-N table and adds an explanatory paragraph. We also wire in the optional positive-control sanity-check CSV if present.

**Tech Stack:** bash (heredoc-driven LaTeX templating), awk (CSV summarization), pdflatex.

**Testing:** Stage 09 is a bash script that builds a LaTeX file. We unit-test the helper logic (column-index lookup + summarization) with a small bash harness in `tests/unit/test_report_summary.sh` against a synthetic CSV that mirrors the May-2026 schema. We do **not** invoke `pdflatex` in tests (heavy/install-dependent). End-to-end PDF compilation is verified manually on Unity after the next chain run.

---

## Task 0: Snapshot the May-2026 ranked-CSV schema as a fixture

**Files:**
- Create: `tests/unit/fixtures/ranked_candidates_may2026_schema.csv` (synthetic, 4 rows + header — covers High/Medium/Low confidence, BUSTED significant, HCR-friendly true/false, native+miniprot CDS sources, one tandem cluster member)

**Why:** Locks the column ordering and naming the report consumes. If `rank_candidates.py` changes the schema, the test fails loud.

**Step 1: Read the canonical schema once**

The full output_cols list lives in `scripts/rank_candidates.py:1903-1928` and `2117-2142` (two writes, identical schema; HCR augmentation in `add_hcr_columns.py` appends `cds_length_bp,paralog_min_identity,closest_paralog_id,hcr_probe_friendly`). Final column order:

```
id,rank_score,confidence_tier,evidence_completeness,
phylo_score,purifying_score,positive_score,selection_significant,
busted_s_p,busted_s_significant,busted_mh_p,busted_mh_significant,
meme_n_episodic_sites,meme_fraction_episodic_sites,
synteny_score,has_synteny_data,expression_score,has_expression_data,
lse_depth_score,raw_tree_depth,
chemosensory_expr_score,has_chemosensory_expr_data,
gprotein_coexpr_score,has_gprotein_data,
ecl_divergence_score,has_ecl_data,
expansion_score,has_expansion_data,
og_confidence_score,has_og_confidence_data,
tandem_cluster_score,tandem_cluster_size,tandem_cluster_id,has_tandem_cluster_data,
cds_source,
cds_length_bp,paralog_min_identity,closest_paralog_id,hcr_probe_friendly
```

**Step 2: Write the fixture CSV** — 4 candidates illustrating: (1) High-conf, BUSTED-S+MH sig, HCR-friendly, native CDS, in tandem cluster; (2) Medium-conf, no BUSTED sig, MEME 3 episodic sites, HCR-friendly, miniprot CDS; (3) Low-conf, sparse data, not HCR-friendly (paralog identity too high); (4) Medium, native CDS, BUSTED-MH only.

**Step 3: Commit**

```bash
git add tests/unit/fixtures/ranked_candidates_may2026_schema.csv
git commit -m "test(09): fixture for ranked CSV May-2026 schema"
```

---

## Task 1: Add bash helper `csv_col_idx` + header-aware `read_csv_count`

**Files:**
- Modify: `09_report_generation.sh:43-61` (replace the existing `count_lines`/`read_csv_value` block with header-aware helpers)
- Create: `tests/unit/test_report_summary.sh`

**Why:** The current report uses positional awk (`$3=="High"`, `$8=="True"`). The new schema is wider; positional reads will break the moment any column shifts. Header-aware lookup is a one-time fix.

**Step 1: Write the failing test**

Create `tests/unit/test_report_summary.sh`:

```bash
#!/usr/bin/env bash
# Unit tests for stage 09 report-summary helpers.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
FIXTURE="$PROJECT_ROOT/tests/unit/fixtures/ranked_candidates_may2026_schema.csv"

# Source the helpers (extracted from 09_report_generation.sh)
# shellcheck source=/dev/null
source "$PROJECT_ROOT/scripts/_report_summary_lib.sh"

fail=0
assert_eq() {
    if [ "$1" != "$2" ]; then
        echo "FAIL: $3 — expected '$2', got '$1'" >&2
        fail=1
    else
        echo "PASS: $3"
    fi
}

assert_eq "$(csv_col_idx "$FIXTURE" id)" "1" "csv_col_idx id == 1"
assert_eq "$(csv_col_idx "$FIXTURE" confidence_tier)" "3" "csv_col_idx confidence_tier"
assert_eq "$(csv_col_idx "$FIXTURE" hcr_probe_friendly)" "39" "csv_col_idx hcr_probe_friendly"

assert_eq "$(count_by_value "$FIXTURE" confidence_tier High)" "1" "count High"
assert_eq "$(count_by_value "$FIXTURE" confidence_tier Medium)" "2" "count Medium"
assert_eq "$(count_by_value "$FIXTURE" confidence_tier Low)" "1" "count Low"
assert_eq "$(count_by_value "$FIXTURE" selection_significant True)" "1" "count selection_significant"
assert_eq "$(count_by_value "$FIXTURE" busted_s_significant True)" "1" "count busted_s_significant"
assert_eq "$(count_by_value "$FIXTURE" busted_mh_significant True)" "2" "count busted_mh_significant"
assert_eq "$(count_by_value "$FIXTURE" hcr_probe_friendly True)" "2" "count hcr_probe_friendly"
assert_eq "$(count_by_value "$FIXTURE" has_tandem_cluster_data True)" "1" "count tandem cluster"
assert_eq "$(count_by_value "$FIXTURE" cds_source native)" "2" "count native CDS"
assert_eq "$(count_by_value "$FIXTURE" cds_source miniprot)" "1" "count miniprot CDS"

assert_eq "$(sum_column "$FIXTURE" meme_n_episodic_sites)" "5" "sum meme episodic sites"

exit $fail
```

**Step 2: Run test to verify it fails**

```bash
bash tests/unit/test_report_summary.sh
```
Expected: FAIL — `_report_summary_lib.sh` does not exist.

**Step 3: Write minimal implementation**

Create `scripts/_report_summary_lib.sh`:

```bash
#!/usr/bin/env bash
# Header-aware CSV summarization helpers used by 09_report_generation.sh.
# All functions echo a numeric or single-token value; never mutate state.

# csv_col_idx <csv> <col_name> -> 1-based index, or empty string if not found.
csv_col_idx() {
    local file="$1" col="$2"
    [ -f "$file" ] || { echo ""; return 0; }
    awk -F',' -v col="$col" 'NR==1{for(i=1;i<=NF;i++)if($i==col){print i;exit}}' "$file"
}

# count_by_value <csv> <col_name> <value> -> count of data rows where col==value.
count_by_value() {
    local file="$1" col="$2" val="$3"
    local idx
    idx="$(csv_col_idx "$file" "$col")"
    [ -z "$idx" ] && { echo "0"; return 0; }
    awk -F',' -v idx="$idx" -v v="$val" 'NR>1 && $idx==v {n++} END{print n+0}' "$file"
}

# sum_column <csv> <col_name> -> numeric sum of column over data rows. Treats blanks as 0.
sum_column() {
    local file="$1" col="$2"
    local idx
    idx="$(csv_col_idx "$file" "$col")"
    [ -z "$idx" ] && { echo "0"; return 0; }
    awk -F',' -v idx="$idx" 'NR>1 {s += ($idx + 0)} END{printf "%g", s+0}' "$file"
}

# read_csv_value <csv> <col_name> [row=2] -> single value at (row, col). Default row=2.
read_csv_value() {
    local file="$1" col="$2" row="${3:-2}"
    local idx
    idx="$(csv_col_idx "$file" "$col")"
    [ -z "$idx" ] && { echo ""; return 0; }
    awk -F',' -v idx="$idx" -v row="$row" 'NR==row {print $idx}' "$file"
}
```

**Step 4: Run test to verify it passes**

```bash
bash tests/unit/test_report_summary.sh
```
Expected: PASS for all 12 assertions.

**Step 5: Commit**

```bash
git add scripts/_report_summary_lib.sh tests/unit/test_report_summary.sh
git commit -m "feat(09): header-aware CSV summary helpers + unit tests"
```

---

## Task 2: Wire helpers into 09_report_generation.sh

**Files:**
- Modify: `09_report_generation.sh:43-106`

**Why:** Replace positional awk with the new header-aware helpers. Source `_report_summary_lib.sh` after `functions.sh`. Recompute existing stats (HIGH_CONF / MED_CONF / LOW_CONF / SIG_SELECTION) using the helpers; add the new aggregates the report will reference in later tasks.

**Step 1: Modify the script**

After `source functions.sh`, add:
```bash
source "${SCRIPTS_DIR}/_report_summary_lib.sh"
```

Replace the existing `count_lines` / `read_csv_value` definitions and the `if [ -f "$RANKED_FILE" ]; then ... fi` block (~ lines 96-106) with:

```bash
TOTAL_CANDIDATES=$(count_lines "${RESULTS_DIR}/candidates/chemogpcr_candidates.txt")
RANKED_FILE="${RANKING_DIR}/ranked_candidates_sorted.csv"

HIGH_CONF=""; MED_CONF=""; LOW_CONF=""; SIG_SELECTION=""
BUSTED_S_SIG=""; BUSTED_MH_SIG=""; MEME_TOTAL=""
HCR_FRIENDLY=""; TANDEM_MEMBERS=""
CDS_NATIVE=""; CDS_MINIPROT=""
TOP_N_RANKED=""

if [ -f "$RANKED_FILE" ]; then
    HIGH_CONF=$(count_by_value "$RANKED_FILE" confidence_tier High)
    MED_CONF=$(count_by_value "$RANKED_FILE" confidence_tier Medium)
    LOW_CONF=$(count_by_value "$RANKED_FILE" confidence_tier Low)
    SIG_SELECTION=$(count_by_value "$RANKED_FILE" selection_significant True)
    BUSTED_S_SIG=$(count_by_value "$RANKED_FILE" busted_s_significant True)
    BUSTED_MH_SIG=$(count_by_value "$RANKED_FILE" busted_mh_significant True)
    MEME_TOTAL=$(sum_column "$RANKED_FILE" meme_n_episodic_sites)
    HCR_FRIENDLY=$(count_by_value "$RANKED_FILE" hcr_probe_friendly True)
    TANDEM_MEMBERS=$(count_by_value "$RANKED_FILE" has_tandem_cluster_data True)
    CDS_NATIVE=$(count_by_value "$RANKED_FILE" cds_source native)
    CDS_MINIPROT=$(count_by_value "$RANKED_FILE" cds_source miniprot)
    TOP_N_RANKED=$(awk -F',' 'NR>1{n++}END{print n+0}' "$RANKED_FILE")
fi
```

**Step 2: Verify the script still parses**

```bash
bash -n 09_report_generation.sh
```
Expected: no output (clean syntax).

**Step 3: Commit**

```bash
git add 09_report_generation.sh
git commit -m "refactor(09): swap positional awk for header-aware helpers"
```

---

## Task 3: Refresh Executive-Summary table with new aggregates

**Files:**
- Modify: `09_report_generation.sh:149-169`

**Why:** Surface BUSTED/MEME/HCR/tandem/CDS-provenance counts in the headline stats table. Without this, the headline numbers say nothing about the May-2026 signals.

**Step 1: Replace the Pipeline Summary Statistics table**

Change the `cat >> ... <<EOF` block currently writing `Total / High / Medium / Low / Sig. Selection` to:

```latex
\begin{table}[H]
\centering
\caption{Pipeline Summary Statistics (May 2026 schema)}
\begin{tabular}{lr}
\toprule
\textbf{Metric} & \textbf{Value} \\\\
\midrule
Total GPCR Candidates Identified & ${TOTAL_CANDIDATES:-N/A} \\\\
Ranked candidates & ${TOP_N_RANKED:-N/A} \\\\
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
HCR probe-friendly candidates & ${HCR_FRIENDLY:-0} \\\\
In tandem cluster & ${TANDEM_MEMBERS:-0} \\\\
\midrule
CDS source: native & ${CDS_NATIVE:-0} \\\\
CDS source: miniprot-recovered & ${CDS_MINIPROT:-0} \\\\
\bottomrule
\end{tabular}
\end{table}
```

**Step 2: bash -n check**

```bash
bash -n 09_report_generation.sh
```

**Step 3: Commit**

```bash
git add 09_report_generation.sh
git commit -m "feat(09): expand executive-summary table for May-2026 signals"
```

---

## Task 4: Refresh Methods section to match current pipeline

**Files:**
- Modify: `09_report_generation.sh:184-213, 380-462` (Methods Overview, Pipeline Architecture, Ranking Algorithm, Software Versions, Scoring Weights)

**Why:** The Methods section currently says DeepTMHMM / FastTree / aBSREL / FastML / MCScanX. The May-2026 stack is TMbed primary, IQ-TREE 3, GARD→BUSTED-S→BUSTED-MH→aBSREL→MEME, IQ-TREE 3 ancestral, JCVI MCscan, GeneRax. Methods narrative must reflect what the pipeline actually ran.

**Step 1: Pipeline Architecture list (~lines 187-199)** — replace with:

```latex
\begin{enumerate}
    \item \textbf{Reference Processing}: Build per-orthogroup HMMs from Nath et al. references.
    \item \textbf{GPCR Identification}: HHblits homology + TMbed transmembrane prediction (DeepTMHMM / Phobius / KD as graded fallbacks; 6+ TM filter). Provenance recorded per run.
    \item \textbf{Orthology Clustering}: OrthoFinder (Diamond mode).
    \item \textbf{Phylogenetic Analysis}: FAMSA/MAFFT regime-based alignment $\to$ HmmCleaner (segment) $\to$ ClipKit (column) $\to$ FastTree seed $\to$ IQ-TREE 3 (ModelFinder; UFBoot+SH-aLRT+TBE; deterministic seed) $\to$ TreeShrink rogue-taxon removal.
    \item \textbf{Selective Pressure}: HyPhy stack — GARD recombination breakpoints $\to$ BUSTED-S (gene-wide) $\to$ BUSTED-MH (multi-hit-aware) $\to$ aBSREL (branch-specific) $\to$ MEME (episodic site-level). MACSE v2 codon alignment.
    \item \textbf{Ancestral Reconstruction}: IQ-TREE 3 \texttt{--ancestral} (FastML fallback).
    \item \textbf{Reconciliation}: GeneRax (Notung fallback).
    \item \textbf{Synteny}: JCVI MCscan (MCScanX fallback); per-gene anchor counts feed ranking.
    \item \textbf{Tandem-cluster detection}: cluster sibling paralogs within sliding genomic windows; the field's signature chemoreceptor signal.
    \item \textbf{Candidate Ranking}: Multi-criteria weighted scoring (see Section~\ref{sec:weights}).
    \item \textbf{HCR augmentation}: per-candidate \texttt{cds\_length\_bp}, closest-paralog identity, and \texttt{hcr\_probe\_friendly} boolean for wet-lab probe design.
    \item \textbf{Structural Analysis}: AlphaFold prediction + Foldseek vs GPCRdb (P3).
\end{enumerate}
```

**Step 2: Ranking Algorithm list (~lines 203-211)** — replace with:

```latex
\begin{itemize}
    \item \textbf{Phylogenetic proximity} to known chemoreceptors (inverse weighted distance).
    \item \textbf{Positive selection}: \texttt{omega\_max} from aBSREL (was weighted-mean — fixed May 2026); episodic-positive selection signal preserved.
    \item \textbf{Purifying selection}: \emph{disabled by default} (\texttt{PURIFYING\_WEIGHT=0}) — chemoreceptor identification rewards diversifying, not housekeeping conservation.
    \item \textbf{Synteny}: JCVI per-gene anchor count (replaces broken count-of-IDs input).
    \item \textbf{Tandem cluster}: sibling paralog clustering signal (chemoreceptor-defining).
    \item \textbf{LSE depth}: position in lineage-specific expansion clades (75th-percentile threshold).
    \item \textbf{Chemosensory expression}: tissue-specific signal where available.
    \item \textbf{G-protein coexpression}: G$\alpha$ family coexpression with candidate (orthogroup-level).
    \item \textbf{ECL divergence}: extracellular-loop sequence divergence vs paralogs.
    \item \textbf{Family expansion} \& \textbf{orthogroup confidence} (auxiliary).
\end{itemize}

P-values from BUSTED-S, BUSTED-MH, aBSREL, and MEME are FDR-corrected via \texttt{statsmodels.multipletests} (Benjamini-Hochberg, $\alpha=0.05$). The composite score is multiplied by an \texttt{evidence\_completeness} factor (floor 0.4) so sparse-data candidates do not co-rank with complete-data candidates.
```

**Step 3: Software Versions table (~lines 442-460)** — replace with the May-2026 toolset:

```latex
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
ClipKit & Column-trimming \\
MACSE v2 & Frameshift-aware codon alignment \\
FastTree & Approximate-ML seed tree \\
IQ-TREE 3 & ML phylogeny + ancestral reconstruction \\
TreeShrink & Rogue-taxon removal \\
HyPhy: GARD / BUSTED-S / BUSTED-MH / aBSREL / MEME & Selection-stack analysis \\
GeneRax & Gene-tree/species-tree reconciliation \\
JCVI MCscan & Synteny / collinearity \\
AlphaFold & Protein structure prediction \\
Foldseek & Structural search vs GPCRdb \\
TM-align & Pairwise structural alignment \\
\bottomrule
\end{tabular}
```

**Step 4: Scoring Weights table (~lines 467-482)** — replace with:

```latex
\label{sec:weights}
\begin{tabular}{lcc}
\toprule
\textbf{Component} & \textbf{Default Weight} & \textbf{Variable} \\
\midrule
Phylogenetic proximity & 2.0 & \texttt{PHYLO\_WEIGHT} \\
Positive selection & 1.0 & \texttt{POSITIVE\_WEIGHT} \\
Purifying selection & 0.0 & \texttt{PURIFYING\_WEIGHT} \\
Synteny conservation & 3.0 & \texttt{SYNTENY\_WEIGHT} \\
Tandem cluster & 2.5 & \texttt{TANDEM\_CLUSTER\_WEIGHT} \\
LSE depth & 1.0 & \texttt{LSE\_DEPTH\_WEIGHT} \\
Chemosensory expression & 1.0 & \texttt{CHEMO\_EXPR\_WEIGHT} \\
G-protein coexpression & 1.0 & \texttt{GPROTEIN\_WEIGHT} \\
ECL divergence & 1.0 & \texttt{ECL\_WEIGHT} \\
Family expansion & 0.5 & \texttt{EXPANSION\_WEIGHT} \\
Orthogroup confidence & 0.5 & \texttt{OG\_CONFIDENCE\_WEIGHT} \\
\bottomrule
\end{tabular}
```

(Weights confirmed from `scripts/rank_candidates.py` env-var defaults; if any weight differs, fix to the actual default before committing.)

**Step 5: bash -n + commit**

```bash
bash -n 09_report_generation.sh
git add 09_report_generation.sh
git commit -m "docs(09): refresh Methods/Tools/Weights for May-2026 stack"
```

---

## Task 5: Add HCR Probe-Friendliness section + extend Top-N table

**Files:**
- Modify: `09_report_generation.sh:343-378` (Top-N candidates table) and a new subsection inserted before Structural Analysis

**Why:** HCR probe-friendliness is the deliverable for the wet-lab. The current Top-15 table shows only `(rank, id, score, confidence, sig_selection)`. The biologists looking at the report need `cds_length_bp`, `paralog_min_identity`, `hcr_probe_friendly`, `tandem_cluster_id` at a glance.

**Step 1: New subsection (insert before `\subsection{Structural Analysis}`)**:

```latex
\subsection{HCR Probe-Friendliness}

For each candidate, HCR probe-friendliness is a boolean derived from CDS length and the closest-paralog identity (windowed). True iff:
\begin{itemize}
  \item \texttt{cds\_length\_bp} $\geq$ minimum probe-design length (default 600 bp)
  \item \texttt{paralog\_min\_identity} $\leq$ maximum cross-reactivity threshold (default 75\%)
\end{itemize}
Candidates marked \texttt{hcr\_probe\_friendly = True} are recommended for first-round HCR design without further filtering. False candidates are not disqualified but require either probe-design effort against a divergent region or a discriminating split-probe strategy.
```

**Step 2: Replace the Top-15 table (~lines 343-378)** with a header-aware writer that includes new columns. Replace the heredoc + awk block with:

```bash
if [ -f "$RANKED_FILE" ]; then
    cat >> "${RESULTS_DIR}/report/report.tex" <<'EOF'

\subsubsection{Top Ranked Candidates}

\begin{longtable}{rlrlcccr}
\caption{Top 20 ranked GPCR candidates. \textbf{Sel}=aBSREL branch-significant; \textbf{B-MH}=BUSTED-MH gene-significant; \textbf{HCR}=probe-friendly; \textbf{TC}=tandem-cluster member; \textbf{CDS}=N(ative) / M(iniprot)} \\
\toprule
\textbf{\#} & \textbf{ID} & \textbf{Score} & \textbf{Conf} & \textbf{Sel} & \textbf{B-MH} & \textbf{HCR} & \textbf{TC size} \\
\midrule
\endfirsthead
\toprule
\textbf{\#} & \textbf{ID} & \textbf{Score} & \textbf{Conf} & \textbf{Sel} & \textbf{B-MH} & \textbf{HCR} & \textbf{TC size} \\
\midrule
\endhead
\bottomrule
\endfoot
EOF

    # Header-aware emit; 20 rows.
    python3 - "$RANKED_FILE" >> "${RESULTS_DIR}/report/report.tex" <<'PY'
import csv, sys
fn = sys.argv[1]
def yn(v): return "Y" if str(v).strip().lower() in ("true","1","yes","y") else "N"
def src1(v):
    s = str(v).strip().lower()
    return "N" if s.startswith("nat") else ("M" if s.startswith("min") else "?")
with open(fn) as fh:
    r = csv.DictReader(fh)
    for i, row in enumerate(r, start=1):
        if i > 20: break
        conf = row.get("confidence_tier", "")
        color = {"High":"highconf","Medium":"medconf","Low":"lowconf"}.get(conf, "lowconf")
        try: score = float(row.get("rank_score","0") or 0)
        except ValueError: score = 0.0
        tc = row.get("tandem_cluster_size", "0") or "0"
        # LaTeX-escape underscores in IDs
        rid = row.get("id","").replace("_", r"\_")
        print(rf"{i} & {rid} & {score:.3f} & \textcolor{{{color}}}{{{conf}}} & "
              rf"{yn(row.get('selection_significant'))} & "
              rf"{yn(row.get('busted_mh_significant'))} & "
              rf"{yn(row.get('hcr_probe_friendly'))} & {tc} \\")
PY

    cat >> "${RESULTS_DIR}/report/report.tex" <<'EOF'
\end{longtable}

EOF
fi
```

**Step 3: bash -n check**

```bash
bash -n 09_report_generation.sh
```

**Step 4: Smoke-test the python writer with the fixture**

```bash
python3 -c "
import csv, sys
fn = 'tests/unit/fixtures/ranked_candidates_may2026_schema.csv'
with open(fn) as fh:
    r = csv.DictReader(fh)
    rows = list(r)
assert len(rows) >= 4, 'fixture too short'
assert {'hcr_probe_friendly','tandem_cluster_size','busted_mh_significant'} <= set(r.fieldnames), 'fixture missing columns'
print('fixture OK,', len(rows), 'rows')
"
```
Expected: `fixture OK, 4 rows`.

**Step 5: Commit**

```bash
git add 09_report_generation.sh
git commit -m "feat(09): HCR probe-friendliness section + extended Top-20 table"
```

---

## Task 6: Wire positive-control sanity check into the report

**Files:**
- Modify: `09_report_generation.sh` (insert after the ranking section)

**Why:** Stage 07 runs `check_positive_controls.py` which writes a CSV of HCR-validated controls (Gα_olf etc.) and their ranking percentile. If a control falls below the 50th percentile we want a glaring warning in the PDF.

**Step 1: Add a section reader**

Locate the actual output path of `check_positive_controls.py`:

```bash
grep -n "check_positive_controls\|positive_control" /home/workspace/Desktop/projects/umass/berghia-chemogpcrs/scripts/check_positive_controls.py /home/workspace/Desktop/projects/umass/berghia-chemogpcrs/07_candidate_ranking.sh | head -20
```

(Read both the script and the call-site to confirm the output filename and columns. Update Step 2 to match.)

**Step 2: Insert subsection** (after the Top-N table, before Structural Analysis):

```bash
POSCTRL_CSV="${RANKING_DIR}/positive_control_check.csv"  # confirm filename in step 1
if [ -f "$POSCTRL_CSV" ]; then
    cat >> "${RESULTS_DIR}/report/report.tex" <<'EOF'

\subsection{Positive-Control Sanity Check}

HCR-validated chemoreceptor controls (e.g.\ G$\alpha_{\text{olf}}$) are scored by the same pipeline. A control falling below the 50th percentile is a red flag for the ranking weights, not the control.

\begin{table}[H]
\centering
\begin{tabular}{lrrr}
\toprule
\textbf{Control} & \textbf{Rank} & \textbf{Percentile} & \textbf{Status} \\
\midrule
EOF
    # Emit one row per positive control.
    python3 - "$POSCTRL_CSV" >> "${RESULTS_DIR}/report/report.tex" <<'PY'
import csv, sys
fn = sys.argv[1]
with open(fn) as fh:
    r = csv.DictReader(fh)
    for row in r:
        name = row.get("control_id","").replace("_", r"\_")
        rank = row.get("rank","-")
        pct  = row.get("percentile","-")
        ok   = (str(row.get("passes_50th","")).strip().lower() in ("true","1","yes","y"))
        status = r"\textcolor{highconf}{PASS}" if ok else r"\textcolor{lowconf}{FAIL}"
        print(rf"{name} & {rank} & {pct} & {status} \\")
PY
    cat >> "${RESULTS_DIR}/report/report.tex" <<'EOF'
\bottomrule
\end{tabular}
\end{table}

EOF
fi
```

(If the column names in `positive_control_check.csv` differ, adjust the `row.get(...)` keys to match the actual file. This is the only step where the test fixture in Task 0 doesn't help — verify against a real run output before merging.)

**Step 3: bash -n + commit**

```bash
bash -n 09_report_generation.sh
git add 09_report_generation.sh
git commit -m "feat(09): surface positive-control sanity check in report"
```

---

## Task 7: Update Recommendations section to mention HCR pathway explicitly

**Files:**
- Modify: `09_report_generation.sh:484-502` (Recommendations section)

**Why:** Make explicit that the deliverable is a ranked list of HCR-prioritized candidates and that the High-confidence + HCR-probe-friendly intersection is the wet-lab shortlist.

**Step 1: Replace the Recommendations subsection content**

```latex
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
    \item Expression profiling across head/foot/digestive/sensory tissues — current expression-based scores carry the starvation+depth confound noted in the project memory.
    \item Comparative genomics with chromosome-level mollusc assemblies (synteny weight is highest at 3.0 for a reason).
    \item AlphaFold + Foldseek structural placement against GPCRdb for top-ranked candidates.
    \item Deorphanization / ligand-binding prediction once HCR localization narrows the functional candidate set.
\end{itemize}
```

**Step 2: bash -n + commit**

```bash
bash -n 09_report_generation.sh
git add 09_report_generation.sh
git commit -m "docs(09): refocus Recommendations on HCR wet-lab shortlist"
```

---

## Task 8: Final verification

**Files:** none (read-only)

**Step 1: Re-run unit tests**

```bash
bash tests/unit/test_report_summary.sh
```
Expected: all PASS.

**Step 2: Confirm no regressions in existing unit suite**

```bash
cd /home/workspace/Desktop/projects/umass/berghia-chemogpcrs
python3 -m pytest tests/unit/ -q 2>&1 | tail -20
```
Expected: previous test count holds (104+).

**Step 3: bash -n on all touched scripts**

```bash
bash -n 09_report_generation.sh
bash -n scripts/_report_summary_lib.sh
bash -n tests/unit/test_report_summary.sh
```
Expected: silent.

**Step 4: Confirm git state is clean and pushed**

```bash
git status
git log --oneline -8
```

**Step 5: Update bead and sync**

```bash
bd update berghia-chemogpcrs-07g --note "Stage 09 templates refreshed for May-2026 schema. PDF compilation will be verified on Unity once stage-08 completes and the chain reaches 09. Engineering is done; verification is runtime."
bd sync
```

(Do **not** close the bead yet — final close happens after Unity produces a PDF that visually checks out.)
