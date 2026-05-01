# Pipeline Review Fixes — Master Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use `superpowers-extended-cc:subagent-driven-development` to implement this plan task-by-task. Heavy-compute tasks (re-run HPC tree, run AlphaFold) are documented but **deferred** to user-side HPC; this plan only delivers the local-runnable code/config/test changes.

**Goal:** Convert the 17 prioritized beads issues from the May 2026 review into a logically-sound, statistically-rigorous, fully-functional ranking pipeline whose deliverable is a top-N list of *Berghia stephanieae* GPCR candidates for prospective HCR validation.

**Architecture:** Keep the 9-stage SLURM/bash + Python architecture. Patch correctness defects in-place. Add the missing tandem-cluster feature once the public Berghia genome is integrated. Modernize tools by config swap (TMbed, JCVI MCscan, GeneRax, BUSTED-MH) without rewriting the orchestration layer. Add a unit-test suite under `tests/unit/` that exercises the scoring/parsing/statistical functions with synthetic inputs, so future regressions surface fast.

**Tech Stack:** Python 3.10/3.13 (numpy/pandas/scipy/statsmodels/ete3/Biopython) · pytest · bash · IQ-TREE 3 / FastTree / MAFFT / ClipKit / HyPhy / OrthoFinder · external new tools (TMbed, MACSE v2, GeneRax, JCVI MCscan, GARD, MEME, TreeShrink, HmmCleaner) · NCBI Datasets CLI · Beads (project tracking).

**Constraints (do not violate):**
- No heavy compute on local 15 GB ARM machine — HPC reruns are deferred to user submission, with the necessary scripts/configs delivered.
- Expression is a *soft* signal, never a gate (rhinophore RNA-seq has starvation + low-depth confound; Gα_olf example).
- Validation method is per-candidate HCR; the ranker output is the deliverable.
- Berghia genome `GCA_034508935.3` (Goodheart 2024) is now public and must be integrated.
- **Use best available tools; do not reinvent.** See Tool Selection Rationale below.

**Beads issues addressed:** -ea9, -wux, -ryr, -6nh, -ce4, -4nu, -ar8, -hg1, -0ku, -cio, -qu9, -m6k, -mqt, -i61, -5b0, -urk, -30g, -o88→-e59, -iof, -j44, -j6f, -bdu, -xqz, -edx, -s6v.

---

## Tool Selection Rationale (audit against "don't reinvent")

For every place where the plan adds non-trivial code, the table below documents whether an existing, maintained library/tool is being used instead of a hand-rolled implementation. The "Why" column states the test that justifies the choice.

| Task | Need | Chosen | Why (vs. alternatives) |
|---|---|---|---|
| 1.1 | Multiple-test correction | `statsmodels.stats.multitest.multipletests` (BH method) | Reference implementation, tested across thousands of papers; reinventing a 25-line BH function caused the bug we are now fixing. |
| 1.2 | Keyword classification of reference names | Standard library `re` with word-boundary patterns | The keyword set is project-specific; using a generic NLP classifier is overkill. The fix is the regex itself. |
| 1.3 / 1.4 | aBSREL JSON parsing | Standard library `json` + `pandas`; no custom parser library exists for HyPhy output | HyPhy doesn't ship a parser; pandas IS the right tool. |
| 1.6 | Multi-criteria score combination | Custom (no library substitute) | `scikit-mcdm` exists for MCDM but is overkill for weighted-sum-with-completeness and adds a dependency for one function. |
| 1.7 | Synteny anchor counting | `pandas` aggregations on JCVI `.anchors` files (Phase 5.5) | Replaces previous hand-rolled anchor counter. |
| 1.8 | Provenance JSON append | `jq` (preferred) or Python `json` stdlib | Idiomatic shell tool; do not reinvent. |
| 1.9 | TM-helix prediction | **TMbed** (Bernhofer & Rost 2022) primary, **DeepTMHMM** fallback, **Phobius** optional consensus | TMbed is current SOTA (+9 pp over DeepTMHMM); both are pip-installable. The deprecated Kyte-Doolittle fallback is being removed. |
| 1.10 (a–h) | Various small fixes | All use existing libs (`re`, `Bio.SeqIO`, ete3, pandas) | No new code introduced where stdlib suffices. |
| 2.1 | Genome download from NCBI | **NCBI Datasets CLI** (`datasets`) — already installed | Official NCBI tool; replaces ad-hoc curl/wget. |
| 2.3 | GFF parsing for CDS lookup | **`gffutils`** (Python) | Battle-tested GFF/GTF library; reads into SQLite for fast lookups; far better than ad-hoc `awk`. Add to `environment.yml`. |
| 4.1 | Tandem-cluster detection | **`gffutils`** for parsing + custom 30-line sliding-window algorithm | The algorithm itself ("≥N paralogs within W kb on same scaffold") is so simple (≤30 lines) that a library import is more overhead than the implementation; gffutils handles the GFF parsing where it matters. Considered MCScanX `dup_classifier` and JCVI tandem — both require pre-computed all-vs-all BLAST and are designed for *whole-genome* tandem analysis, not restricted-to-candidate-set sliding windows. |
| 4.2 | Null-distribution birth-death tree simulation | **`dendropy.simulate.treesim.birth_death_tree`** | Already installed locally; canonical Python phylo simulator. Replaces hand-rolled Yule simulator. |
| 4.2 | TBE bootstrap | **IQ-TREE 3 `--tbe`** built-in | Native; do not reinvent. |
| 5.1 | Per-sequence alignment cleaning | **HmmCleaner.pl** (Di Franco 2019) wrapper | Established tool; do not reinvent. |
| 5.1 | Codon-aware MSA | **MACSE v2** (Ranwez 2018) wrapper | Established tool; replaces naive pal2nal back-translation. |
| 5.3 | Selection analysis | **HyPhy** built-ins: `gard`, `busted` (with `--srv Yes`), `busted` (with `--multiple-hits Double+Triple`), `absrel`, `meme` | All native to HyPhy; only the orchestration script is custom. |
| 5.4 | Reconciliation | **GeneRax** (Morel 2020) | Established tool; replaces Notung. |
| 5.5 | Synteny | **JCVI MCscan** (`python -m jcvi.compara.catalog ortholog` + `jcvi.graphics.synteny`) | Current SOTA per Tang 2024; replaces MCScanX. |
| 5.6 | Rogue-taxon detection | **TreeShrink** (Mai & Mirarab 2018) | Established tool; do not reinvent. |
| 5.7 | Ancestral sequence reconstruction | **IQ-TREE 3 `-asr`** (model-consistent with the inference tree); **GRASP** for indel ancestors on selected deep nodes | Both established; replaces FastML. |
| 6.1 | Aplysia retrospective validation | NCBI Datasets CLI + the existing pipeline; **no new code beyond a labels CSV** | Pipeline IS the validator. |
| 6.2 | Probe-design-friendliness | Pairwise-distance from existing IQ-TREE/MAFFT alignments via Biopython | The actual HCR probe design is done by Molecular Instruments' service; the pipeline only outputs candidates ranked + a `paralog_minimum_distance_aa` diagnostic column. |
| 6.3 | Positive-control check | `pandas` lookup; ~20 lines | No library needed. |
| 6.4 | Structural comparison (deferred P3) | **Foldseek** + **FoldMason** + **GPCRdb 2025 structural set** | All established; do not reinvent. |

**General rule applied:** whenever there is a maintained tool that is canonical in the field, use it; whenever the algorithm is < 50 lines and the library would add a significant dependency surface, write it carefully and cover with unit tests. The `compute_tandem_clusters` and `calculate_fair_rank_score` are the two genuinely-custom pieces in the new code; both are < 50 lines, both have unit tests in this plan.

---

## Phase 0 — Preserve in-progress work and bootstrap test scaffold

### Task 0.1: Commit existing in-progress bugfixes

**Files:**
- Modify (already): `07_candidate_ranking.sh`, `config.sh`, `scripts/rank_candidates.py`
- Add (already): `scripts/generate_all_tree_figures.sh`, `scripts/recover_cds_from_assemblies.py`, `scripts/visualize_*.py`

**Steps:**
1. `git status` to verify three modified + five untracked files (the visualization scripts and CDS recovery script are user-authored and should be kept).
2. Stage all eight files via explicit `git add` per file (no `-A`) so secrets cannot sneak in.
3. Commit with message: `chore: preserve in-progress visualization scripts and rank_candidates env-var support before review-fix sweep`. Do NOT add Co-Authored-By per user preferences (memory).

### Task 0.2: Create unit-test scaffold

**Files:**
- Create: `tests/unit/__init__.py` (empty)
- Create: `tests/unit/conftest.py`
- Create: `tests/unit/README.md` — explains the layout

**Step 1:** Write minimal `conftest.py` with shared fixtures:

```python
# tests/unit/conftest.py
import sys, os
from pathlib import Path
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

@pytest.fixture
def tmpdir_p(tmp_path):
    return tmp_path
```

**Step 2:** Verify pytest discovery: `cd /home/workspace/Desktop/projects/umass/berghia-chemogpcrs && pytest tests/unit/ --collect-only`. Expected: 0 tests collected (no error).

**Step 3:** Commit: `test: add unit-test scaffold under tests/unit/`.

---

## Phase 1 — P0 correctness fixes (no genome dependency)

These five fixes change *which candidates rank highly*. Must be done before any next end-to-end run.

### Task 1.1 — Fix BH-FDR implementation (bead -wux part 1)

**File:** `scripts/rank_candidates.py:283-308` (`benjamini_hochberg`)

**Step 1 — Failing test:** `tests/unit/test_fdr.py`:

```python
import numpy as np
import pytest
from rank_candidates import benjamini_hochberg

def test_bh_known_example():
    # Benjamini-Hochberg 1995 worked example: pvals sorted -> q
    pvals = [0.001, 0.008, 0.039, 0.041, 0.042, 0.060, 0.074, 0.205, 0.212, 0.216,
             0.222, 0.251, 0.269, 0.275, 0.34, 0.456, 0.648, 0.759, 0.901, 0.96]
    q = benjamini_hochberg(pvals)
    # statsmodels reference (computed with multipletests(method='fdr_bh'))
    expected = [0.0200, 0.0800, 0.1640, 0.1640, 0.1680, 0.2000, 0.2114, 0.4860, 0.4860, 0.4860,
                0.4860, 0.4860, 0.4860, 0.4860, 0.5320, 0.6731, 0.9012, 0.9706, 0.96, 0.96]
    np.testing.assert_allclose(q, expected, atol=1e-3)

def test_bh_monotone():
    # q-values must be monotone non-decreasing in p-value order
    pvals = [0.5, 0.01, 0.3, 0.001, 0.04, 0.9]
    q = benjamini_hochberg(pvals)
    sorted_q = sorted(zip(pvals, q))
    qs = [qq for _, qq in sorted_q]
    assert all(b >= a - 1e-9 for a, b in zip(qs, qs[1:]))

def test_bh_empty():
    assert benjamini_hochberg([]) == []
```

**Step 2:** Run test: `pytest tests/unit/test_fdr.py -v` → expected fail.

**Step 3 — Implementation:** Replace lines 283-308 with `statsmodels.stats.multitest.multipletests`-based wrapper:

```python
def benjamini_hochberg(pvals):
    """BH-FDR q-values, preserving input order. Returns list aligned to pvals."""
    if not len(pvals):
        return []
    from statsmodels.stats.multitest import multipletests
    pvals = np.asarray(pvals, dtype=float)
    nan_mask = np.isnan(pvals)
    if nan_mask.all():
        return [float('nan')] * len(pvals)
    q = np.full_like(pvals, np.nan)
    ok = ~nan_mask
    _, q_ok, _, _ = multipletests(pvals[ok], method='fdr_bh')
    q[ok] = q_ok
    return q.tolist()
```

**Step 4:** Run test → expected PASS. Add `statsmodels>=0.13` to `environment.yml` pip section if not present (it imports successfully in current env, so likely already installed via dependency).

**Step 5:** Commit: `fix(rank): replace buggy BH-FDR with statsmodels.multipletests (#wux)`.

### Task 1.2 — Fix reference-categorization keyword regex (bead -wux part 2)

**File:** `scripts/rank_candidates.py:213-219`

**Step 1 — Failing test:** `tests/unit/test_ref_categorization.py`:

```python
import pytest
from rank_candidates import categorize_reference  # to be exposed

def test_chemoreceptor_keywords_match_word_boundary():
    assert categorize_reference("Olfactory receptor 7A1") == "chemoreceptor"
    assert categorize_reference("OR1A1_HUMAN") == "chemoreceptor"
    assert categorize_reference("Vomeronasal receptor V1R") == "chemoreceptor"
    assert categorize_reference("Taste receptor type 2") == "chemoreceptor"

def test_non_chemoreceptor_not_matched():
    # 'or' substring should NOT match these
    assert categorize_reference("Orexin receptor type 1") == "other_gpcr"
    assert categorize_reference("Orphanin FQ receptor") == "other_gpcr"
    assert categorize_reference("Origin recognition complex protein") == "other_gpcr"
    assert categorize_reference("Adrenergic receptor alpha-2A") == "other_gpcr"
```

**Step 2:** Run → fail (function not exposed yet).

**Step 3 — Implementation:** Refactor lines 213-219 into a top-level `categorize_reference(name: str) -> str` using word-boundary regex:

```python
import re
_CHEMO_PATTERNS = re.compile(
    r'\b(?:olfr|or\d+[a-z]?\d*|odorant|taste\s*r|tasr|t[12]r|vomero(?:nasal)?|v[12]r|formyl[\s_-]?peptide|trace[\s_-]?amine|taar)\b',
    re.IGNORECASE,
)
def categorize_reference(name: str) -> str:
    return "chemoreceptor" if _CHEMO_PATTERNS.search(name or "") else "other_gpcr"
```

Wire the call sites at the original location to use `categorize_reference`.

**Step 4:** Run tests → PASS.

**Step 5:** Commit: `fix(rank): word-boundary regex for chemoreceptor keyword matching (#wux)`.

### Task 1.3 — Fix dN/dS asymmetry + remove `+0.2` expression bias (beads -ea9, -mqt part)

**Files:** `config.sh:179-180`, `scripts/rank_candidates.py:553-593` (`get_selection_scores`), `scripts/expression_analysis.py:288`

**Rationale:** Chemoreceptor identification rewards diversifying selection on extracellular loops, not whole-gene purifying selection. Default `PURIFYING_WEIGHT` to 0; keep it user-configurable. Remove the bare `+0.2` bias term in expression scoring (no rationale, compresses score range).

**Step 1 — Failing test:** `tests/unit/test_selection_scoring.py`:

```python
import os
import pytest
from rank_candidates import get_selection_scores  # to be refactored to take weights as args

def test_purifying_does_not_outscore_neutral_when_weight_zero(monkeypatch):
    # With PURIFYING_WEIGHT=0, a strongly purifying ω=0.05 must score 0 on the purifying axis
    s = get_selection_scores(omega=0.05, p_corrected=0.01, purifying_weight=0.0, positive_weight=1.0)
    assert s["purifying_score"] == 0.0
    assert s["positive_score"] == 0.0   # ω<1 => no positive

def test_strong_positive_selection_is_rewarded():
    s = get_selection_scores(omega=20.0, p_corrected=0.001, purifying_weight=0.0, positive_weight=1.0)
    # |log10(20)| ≈ 1.30; with significance 1.5x boost
    assert s["positive_score"] > 1.5
    assert s["purifying_score"] == 0.0

def test_nan_omega_yields_zero():
    s = get_selection_scores(omega=float('nan'), p_corrected=float('nan'),
                             purifying_weight=0.0, positive_weight=1.0)
    assert s["purifying_score"] == 0.0
    assert s["positive_score"] == 0.0
```

**Step 2:** Refactor `get_selection_scores` so weights are explicit arguments (defaults read from env), and the function returns a dict with `{purifying_score, positive_score, is_significant}`. Use `numpy.log10` (the audit found `|log ω|` — clarify the base).

**Step 3 — Config:**
```diff
-PURIFYING_WEIGHT=1
+# PURIFYING_WEIGHT=0 because chemoreceptor identification rewards diversifying
+# selection on extracellular loops, not whole-gene purifying selection.
+# Set >0 only when looking for conserved-function GPCRs.
+PURIFYING_WEIGHT=0
 POSITIVE_WEIGHT=1
```

**Step 4 — `expression_analysis.py:288`:** Remove the `+ 0.2` bias term:
```python
score = 0.4 * min(expr_component, 1.0) + 0.4 * tau_component + chemo_bonus
```

**Step 5:** Run tests → PASS. Commit: `fix(rank): zero out purifying-selection weight + drop expr +0.2 bias (#ea9 #mqt)`.

### Task 1.4 — Fix aBSREL omega reporting (bead -ea9 part 2)

**File:** `scripts/parse_absrel.py:38-49`

**Rationale:** aBSREL fits a mixture of rate classes per branch. Reporting the weighted mean ω across rate classes destroys the episodic-positive-selection signal that aBSREL is designed to detect. Report the **maximum-ω rate class** alongside the mean, and use the LRT p-value as the primary significance signal.

**Step 1 — Failing test:** `tests/unit/test_parse_absrel.py`:

```python
import json
from pathlib import Path
import pytest
from parse_absrel import extract_branch_omega  # to be exposed

def make_branch(rate_classes):
    return {"Rate Distributions": [list(rc) for rc in rate_classes],
            "Corrected P-value": 0.01, "Uncorrected P-value": 0.01,
            "Number of rate classes": len(rate_classes), "LRT": 12.3}

def test_episodic_positive_selection_preserved():
    # 90% sites at ω=0.1, 10% sites at ω=20 — mean ω ≈ 2.09; max ω = 20
    branch = make_branch([(0.1, 0.9), (20.0, 0.1)])
    out = extract_branch_omega(branch)
    assert out["omega_max"] == pytest.approx(20.0)
    assert out["omega_mean"] == pytest.approx(0.1*0.9 + 20.0*0.1)
    assert out["weight_at_max"] == pytest.approx(0.1)

def test_purifying_only():
    branch = make_branch([(0.05, 1.0)])
    out = extract_branch_omega(branch)
    assert out["omega_max"] == pytest.approx(0.05)
    assert out["omega_mean"] == pytest.approx(0.05)
```

**Step 2 — Implementation:** Replace the parsing block with:

```python
def extract_branch_omega(branch_data: dict) -> dict:
    """Return omega_max, omega_mean, weight_at_max, lrt_p, n_rate_classes."""
    rate_classes = branch_data.get("Rate Distributions", []) or []
    if not rate_classes:
        return {"omega_max": float("nan"), "omega_mean": float("nan"),
                "weight_at_max": float("nan"), "n_rate_classes": 0}
    pairs = [(float(rc[0]), float(rc[1])) for rc in rate_classes]
    omegas, weights = zip(*pairs)
    omega_max = max(omegas)
    weight_at_max = weights[omegas.index(omega_max)]
    omega_mean = sum(o*w for o, w in pairs)
    return {"omega_max": omega_max, "omega_mean": omega_mean,
            "weight_at_max": weight_at_max, "n_rate_classes": len(pairs)}
```

Add `omega_max` and `weight_at_max` to the CSV header. Update `rank_candidates.get_selection_scores` to use `omega_max` (not the mean) for the positive-selection score.

**Step 3:** Run tests → PASS. Commit: `fix(absrel): report max-ω rate class to preserve episodic-selection signal (#ea9)`.

### Task 1.5 — Fix cross-validation no-op (bead -wux part 3)

**File:** `scripts/rank_candidates.py:1314-1326` (`run_leave_one_out_crossval`)

**Rationale:** The current CV holds out *references* and looks them up in the *candidate list* — disjoint sets, so the function silently returns ~0 recovery. Replace with a proper labeled-holdout: hold out a *known chemoreceptor reference*, score it as if it were a candidate (using only the remaining references), and check rank percentile.

**Step 1 — Failing test:** `tests/unit/test_crossval.py`:

```python
import pandas as pd
import pytest
from rank_candidates import run_leave_one_out_crossval

def test_known_chemoreceptor_recovers_in_top_quartile(monkeypatch, tmp_path):
    # Synthetic: 100 references, 10 marked as chemoreceptor; phylo + dN/dS scored.
    # Expect mean rank-percentile >= 0.75 for chemoreceptor holdouts.
    ...  # full fixture; details in implementation step
```

**Step 2 — Implementation:** Refactor `run_leave_one_out_crossval` to:
1. Filter references to those labeled `chemoreceptor` in `ref_categories_final.csv`.
2. For each held-out chemoreceptor R: remove R from the reference set, treat R as a candidate, recompute its phylo+dN/dS scores, find its rank in the joint candidate+ref pool, record rank percentile.
3. Report mean rank-percentile, fraction in top-25%, fraction in top-50%.
4. Document that this is a *retrospective* method-validation, distinct from prospective HCR validation (link to bead -bdu).

**Step 3:** Commit: `fix(rank): proper leave-one-out CV on labeled chemoreceptor refs (#wux)`.

### Task 1.6 — Fix composite-score evidence-completeness multiplier (bead -ce4)

**File:** `scripts/rank_candidates.py:1541-1610` (`calculate_fair_rank_score`)

**Rationale:** Currently sparse-data candidates with one strong signal can outrank well-supported ones. Multiply the final score by `evidence_completeness ∈ [0.4, 1.0]` (clipped at 0.4 floor) so missing-data candidates are explicitly down-weighted but not zeroed.

**Step 1 — Failing test:** `tests/unit/test_composite_score.py`:

```python
def test_complete_outranks_sparse_at_same_per_axis_score():
    sparse = calculate_fair_rank_score(scores={"phylo": 1.0}, weights={"phylo": 2, "expr": 1, "dnds": 1})
    full = calculate_fair_rank_score(scores={"phylo": 1.0, "expr": 1.0, "dnds": 1.0}, weights={"phylo": 2, "expr": 1, "dnds": 1})
    assert full > sparse

def test_no_data_floors_at_zero():
    assert calculate_fair_rank_score(scores={}, weights={"phylo": 2}) == 0.0

def test_evidence_completeness_multiplier_in_range():
    out = calculate_fair_rank_score(scores={"phylo": 1.0, "expr": 0.0}, weights={"phylo": 2, "expr": 1, "dnds": 1}, return_diagnostics=True)
    assert 0.4 <= out["evidence_completeness"] <= 1.0
```

**Step 2 — Implementation:**
```python
def calculate_fair_rank_score(scores, weights, return_diagnostics=False, completeness_floor=0.4):
    available = {k: w for k, w in weights.items() if k in scores and scores[k] is not None}
    if not available:
        score = 0.0
        completeness = 0.0
    else:
        total_weight = sum(weights.values())
        avail_weight = sum(available.values())
        completeness = max(completeness_floor, avail_weight / total_weight)
        weighted_sum = sum(available[k] * scores[k] for k in available)
        # Normalize by available weight (so all-1.0 candidate scores 1.0 max), then scale by completeness
        score = (weighted_sum / avail_weight) * completeness
    return ({"score": score, "evidence_completeness": completeness} if return_diagnostics else score)
```

**Step 3:** Commit: `fix(rank): multiply composite score by evidence_completeness (#ce4)`.

### Task 1.7 — Synteny score: skip degenerate normalization (bead -ce4 part 2 / -mqt)

**File:** `scripts/rank_candidates.py:314-354` (`load_synteny_scores`)

**Step 1 — Failing test:**
```python
def test_synteny_skip_when_max_anchors_low():
    # If only one candidate has 1 anchor, do NOT score it 1.0 — return None for all
    out = load_synteny_scores_from_dict({"a": 1, "b": 0})
    assert all(v is None or v == 0 for v in out.values())

def test_synteny_log_scale_when_real_data():
    out = load_synteny_scores_from_dict({"a": 1, "b": 5, "c": 50})
    assert out["c"] > out["b"] > out["a"]
    # All in (0, 1]
    assert all(0 < v <= 1 for v in out.values())
```

**Step 2:** Replace `count / max_anchors` with `np.log1p(count) / np.log1p(max_anchors)` and skip the entire block when `max_anchors < 5` (return `None` per candidate, which the composite-score function then treats as missing data — correctly down-weighting via evidence_completeness).

**Step 3:** Commit: `fix(rank): log-scale synteny score, skip when degenerate (#ce4 #mqt)`.

### Task 1.8 — Reproducibility: seeds, IQTREE_MODEL split, provenance writes (bead -ryr)

**Files:** `04_phylogenetic_analysis.sh`, `config.sh:118`, `functions.sh:238-251` (`record_provenance`)

**Step 1:** Split `IQTREE_MODEL` into two vars:
```diff
-IQTREE_MODEL="MFP -mset LG,WAG,JTT,Dayhoff,mtREV,cpREV"
+IQTREE_MODEL_FIND="MFP"
+IQTREE_MODEL_SET="LG,VT,WAG,JTT,Dayhoff,mtREV,cpREV"
+IQTREE_SEED=12345          # set explicit seed for reproducibility
+FASTTREE_SEED=12345
```
Update all `iqtree2`/`iqtree3` invocations: `iqtree3 -m "$IQTREE_MODEL_FIND" -mset "$IQTREE_MODEL_SET" -seed "$IQTREE_SEED" ...`. Same for FastTree (`-seed $FASTTREE_SEED`) and MAFFT (no seed needed; deterministic given same input).

**Step 2:** Fix `record_provenance` in `functions.sh:238-251` to actually append `step_record` JSON into `$PROVENANCE_FILE` via `jq` (or in-memory Python if `jq` unavailable). Add a unit test that runs `record_provenance` in a tmpdir and asserts the JSON contains the new step.

**Step 3:** Commit: `fix(repro): pin IQ-TREE/FastTree seeds, split model flags, fix provenance write (#ryr)`.

### Task 1.9 — TM-prediction consensus + drop Kyte-Doolittle fallback (bead -6nh part 1)

**Files:** `scripts/run_deeptmhmm.sh`, `scripts/predict_tm_helices.py`, `02_chemogpcrs_identification.sh`, `config.sh`

**Step 1:** Add a config gate `TM_PREDICTOR_PRIMARY=tmbed` (with `deeptmhmm` legacy default), `TM_PREDICTOR_CONSENSUS=true`. When TMbed is not installed, log a clear warning, fall back to DeepTMHMM, and **never** to Kyte-Doolittle.

**Step 2:** Refactor `run_deeptmhmm.sh` so the fallback chain is: TMbed (preferred) → DeepTMHMM container → DeepTMHMM venv → DeepTMHMM biolib CLI → **fail** (do not silently fall back to KD). Move `predict_tm_helices.py` to `scripts/legacy/predict_tm_helices.py` with a header comment marking it deprecated.

**Step 3:** Add a TMbed installer stanza to `environment.yml` (pip: `tmbed`) and a HOWTO in `docs/` for the user to install on HPC.

**Step 4:** Add unit test for the consensus logic: given two TM predictions for the same gene, the consensus accepts iff `|count_A - count_B| ≤ 1` AND `min(count_A, count_B) ≥ MIN_TM_REGIONS`. This implements the Nath 2025 ±1 tolerance.

**Step 5:** Commit: `feat(tm): TMbed-primary TM consensus, drop Kyte-Doolittle fallback (#6nh)`.

### Task 1.10 — Minor bug bundle from audit (bead -mqt)

Sequential, each ≤5 minutes:

| # | File:line | Fix |
|---|---|---|
| a | `scripts/rank_candidates.py:686` | Substring → exact: `[l for l in leaves if l == cand or l.startswith(cand + '_')]` |
| b | `scripts/select_deep_nodes.py:29` | Replace split[0] with `extract_taxid_from_header` from `lse_refine.py`; compute percentile only over `nodes_with_taxid` |
| c | `scripts/select_diverse_candidates.py:24,34` | Validate `tree.search_nodes(name=id)` returns 1 match per ID; assert distance matrix fully populated |
| d | `scripts/parse_absrel.py:108-117` | Read existing CSV, validate columns vs. new row; write per-OG CSVs and concatenate at end |
| e | `scripts/process_expression.py:88-99` | Require `--sample-tissue-tsv` argument; fail loudly if not provided |
| f | `02a_cluster_sequences.sh` (config.sh:127) | Add Trinity-gene-id pre-cluster pass before CD-HIT (group by `TRINITY_DN\d+_c\d+_g\d+`, keep longest isoform per gene) |
| g | `scripts/conservation_mapping.py:99-118` | Multiply Shannon entropy by `(1 - gap_freq)`; require `n_seq_at_col ≥ 0.5 * alignment_depth` |
| h | `scripts/compute_lrt.py:36-37,62` | Sign-optional regex; `df` as CLI argument; document branch-site `0.5*chi2(0) + 0.5*chi2(1)` mode |

Each gets its own test under `tests/unit/test_minor_fixes.py` and its own commit `fix(audit): <one-liner> (#mqt)`.

---

## Phase 2 — P0 Berghia genome integration (bead -4nu)

### Task 2.1 — Acquire genome assembly + GFF + protein FASTA + CDS FASTA

**Tools:** NCBI Datasets CLI (`datasets`).

**Steps:**
1. `datasets download genome accession GCA_034508935.3 --include genome,gff3,protein,cds --filename /tmp/berghia_genome.zip`
2. Unzip into `genomes/` (create dir): expected files `ncbi_dataset/data/GCA_034508935.3/{genomic.fna,genomic.gff,protein.faa,cds_from_genomic.fna}`
3. Move/symlink to canonical paths: `genomes/taxid_berghia_berghia.fasta` (genome), `genomes/taxid_berghia_berghia.gff3`, `genomes/taxid_berghia_berghia.proteins.fa`, `genomes/taxid_berghia_berghia.cds.fna`.
4. Smoke-check: `samtools faidx genomes/taxid_berghia_berghia.fasta` produces `.fai` without errors; `gffread -y test.faa -g genomes/taxid_berghia_berghia.fasta genomes/taxid_berghia_berghia.gff3 | head -2` shows protein.

### Task 2.2 — Update config.sh to point at real genome

**File:** `config.sh`

Add:
```bash
export GENOME_ACCESSION="GCA_034508935.3"
export GENOME_VERSION="UCSD_Bste_1.2"
export GENOME_GFF="${GENOME_DIR}/taxid_berghia_berghia.gff3"
export GENOME_PROTEIN="${GENOME_DIR}/taxid_berghia_berghia.proteins.fa"
export GENOME_CDS="${GENOME_DIR}/taxid_berghia_berghia.cds.fna"
```

Verify in `06_synteny_and_mapping.sh` that `$GENOME` and `$GENOME_GFF` are now read instead of hardcoded paths.

### Task 2.3 — Berghia CDS substitution in step 05

**File:** `05_selective_pressure_and_asr.sh`, `scripts/recover_cds_from_assemblies.py`

When the candidate is Berghia, use `$GENOME_CDS` directly (lookup by gene ID in the GFF) rather than miniprot recovery. Update `find_nucleotide_sequences` helper. Add a unit test that constructs a tiny GFF + FASTA pair and asserts the helper returns the right CDS.

Commit: `feat(genome): integrate Berghia GCA_034508935.3 + use native CDS for dN/dS (#4nu)`.

---

## Phase 3 — P1 quarantines (no further compute on broken scripts)

### Task 3.1 — Quarantine `binding_site_prediction.py` (bead -hg1)

1. `git mv scripts/binding_site_prediction.py scripts/quarantine/binding_site_prediction.py.deprecated`
2. Add header comment explaining quarantine reason (BW dict dead code, median-pocket tautology) and link to bead -hg1.
3. In `08_structural_analysis.sh` and `rank_candidates.py`, remove all references; replace any consumer columns with `NaN` so missing-data multiplier handles it.
4. Add a placeholder `scripts/binding_site_prediction.py` that emits a clear `NotImplementedError` with link to the bead.
5. Commit: `chore: quarantine binding_site_prediction.py until rebuild on GPCRdb structural alignment (#hg1)`.

### Task 3.2 — Quarantine `convergent_evolution.py` (bead -0ku)

Same pattern as 3.1. Note that `08_structural_analysis.sh` does not currently consume this script's output (verified during audit) so removal is safe. Commit: `chore: quarantine convergent_evolution.py — invalid parsimony ASR + arbitrary clade split (#0ku)`.

### Task 3.3 — Fix or quarantine `compute_lrt.py` (bead -cio)

Already covered in Task 1.10 minor-bug item (h): sign-optional regex + df CLI arg + boundary-test mode. After Task 1.10 the script is *fixed* not quarantined. Mark bead -cio resolved on the commit from Task 1.10.

---

## Phase 4 — P1 methodological gaps (some depend on genome from Phase 2)

### Task 4.1 — Tandem-cluster feature (bead -ar8) — REQUIRES Phase 2 done

**Files:** `scripts/compute_tandem_clusters.py` (new), `scripts/rank_candidates.py`, `config.sh`

**Step 1:** Create `scripts/compute_tandem_clusters.py` using **`gffutils`** for GFF parsing (battle-tested, SQLite-backed) plus a small sliding-window cluster identifier (well below the threshold where a library substitution would be net-positive):

```python
#!/usr/bin/env python3
"""Compute intra-genome tandem clusters of GPCR paralogs.

Uses gffutils for GFF parsing (https://daler.github.io/gffutils/).
The sliding-window cluster algorithm itself is ~20 lines and well-tested.
"""
import argparse, sys
from collections import defaultdict
import pandas as pd
import gffutils  # GFF parsing — established library, do not reinvent

def load_gene_index(gff_path, db_path=None):
    """Build (or reuse) a gffutils SQLite DB and yield (gene_id, scaffold, start, end)."""
    db_path = db_path or gff_path + ".db"
    try:
        db = gffutils.FeatureDB(db_path, keep_order=True)
    except Exception:
        db = gffutils.create_db(gff_path, dbfn=db_path, force=True,
                                merge_strategy="merge", keep_order=True)
    for gene in db.features_of_type("gene"):
        yield gene.id, gene.seqid, int(gene.start), int(gene.end)

def find_clusters(genes, candidates, window_kb=100, min_size=3):
    """genes: iterable of (gid, scaffold, start, end). candidates: set of gids.
    Returns dict gid -> cluster_size (1 if isolated)."""
    by_scaff = defaultdict(list)
    for g in genes:
        if g[0] in candidates:
            by_scaff[g[1]].append(g)
    cluster_size = {gid: 1 for gid in candidates}
    for scaff, gs in by_scaff.items():
        if len(gs) < min_size:
            continue
        gs.sort(key=lambda x: x[2])
        i = 0
        while i < len(gs):
            j = i
            while j + 1 < len(gs) and (gs[j+1][2] - gs[i][2]) <= window_kb * 1000:
                j += 1
            n = j - i + 1
            if n >= min_size:
                for k in range(i, j+1):
                    cluster_size[gs[k][0]] = max(cluster_size[gs[k][0]], n)
                i = j + 1
            else:
                i += 1
    return cluster_size

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--window-kb", type=int, default=100)
    ap.add_argument("--min-size", type=int, default=3)
    ap.add_argument("--db-cache", default=None,
                    help="Path for gffutils SQLite cache (default: <gff>.db)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cands = {line.strip() for line in open(a.candidates) if line.strip()}
    genes = list(load_gene_index(a.gff, a.db_cache))
    sizes = find_clusters(genes, cands, a.window_kb, a.min_size)
    pd.DataFrame(sorted(sizes.items()),
                 columns=["candidate_id","tandem_cluster_size"]).to_csv(a.out, index=False)
    n_clustered = sum(1 for s in sizes.values() if s >= a.min_size)
    print(f"Wrote {len(sizes)} candidates; {n_clustered} in clusters (max size={max(sizes.values())})",
          file=sys.stderr)
```

Add `gffutils>=0.13` to `environment.yml` pip section.

**Step 2:** Unit test `tests/unit/test_tandem_clusters.py` with synthetic GFF.

**Step 3:** Integrate in `rank_candidates.py`:
- Load tandem_cluster_size CSV (path from `$TANDEM_CLUSTER_FILE` env).
- Add `tandem_score = log1p(cluster_size) / log1p(20.0)` (cap at cluster size 20).
- Wire into `weights` dict: `TANDEM_CLUSTER_WEIGHT=2.5` in config (top-tier weight given expression unreliability).

**Step 4:** Add invocation to `06_synteny_and_mapping.sh` after the genome is in place.

**Step 5:** Commit: `feat(rank): intra-genome tandem-cluster feature for chemoreceptor signal (#ar8)`.

### Task 4.2 — Stats-thresholds reframe (bead -m6k)

**Files:** `scripts/rank_candidates.py:90` (LSE percentile), `scripts/select_deep_nodes.py:52` (deep node percentile), `04_phylogenetic_analysis.sh` (UFBoot framing), `09_report_generation.sh` (support summary).

**Step 1:** Replace within-dataset 75/90 percentiles with **dual-mode**: (a) within-dataset percentile (current behavior, kept as `--legacy`), (b) **null-calibrated** threshold from a Yule birth-death simulation (default). Add `scripts/calibrate_depth_threshold.py` that uses **`dendropy.simulate.treesim.birth_death_tree`** (canonical Python phylo simulator, already installed locally) to generate 1000 trees of comparable taxon count, computes branch-distance distribution per tree, takes 95th-percentile across simulations, and outputs `results/calibration/depth_thresholds.json`. Skeleton:

```python
import dendropy
from dendropy.simulate import treesim
import json, numpy as np
def calibrate(n_taxa, n_sims=1000, birth_rate=1.0, death_rate=0.0):
    depths = []
    for _ in range(n_sims):
        t = treesim.birth_death_tree(birth_rate=birth_rate, death_rate=death_rate,
                                     ntax=n_taxa, repeat_until_success=True)
        depths.extend(node.distance_from_root() for node in t.preorder_internal_node_iter())
    return float(np.percentile(depths, 95))
```

**Step 2:** In `04_phylogenetic_analysis.sh`, run IQ-TREE with `-B 1000 -alrt 1000 --tbe` (TBE alongside UFBoot+SH-aLRT). In the report (stage 09), show the support-summary as: `% nodes UFBoot ≥95 AND SH-aLRT ≥80` plus `% nodes TBE ≥0.7`. Keep the `≥70` legacy figure but label it `(legacy threshold for comparison)`.

**Step 3:** Commit: `feat(stats): null-calibrated depth thresholds + UFBoot≥95/TBE reporting (#m6k)`.

---

## Phase 5 — P2 tooling modernization (config + scripts; HPC reruns deferred)

For each, deliver: (a) the config flag/path additions, (b) a wrapper script in `scripts/`, (c) install instructions in `docs/`, (d) a smoke test that the wrapper script's `--help` output looks right. **Do not actually run the heavy tool here** — that's HPC.

### Task 5.1 — Alignment & TM stack: TMbed + HmmCleaner + MACSE v2 (bead -i61)

- Already partially covered in Task 1.9 (TMbed primary).
- Add `scripts/run_hmmcleaner.sh` invoked between MAFFT and ClipKit in `04_phylogenetic_analysis.sh`.
- Add `scripts/run_macse.sh` invoked before pal2nal in `05_selective_pressure_and_asr.sh`.
- Update `environment.yml` and `docs/installation.md`.
- Commit: `feat(align): HmmCleaner segment-clean + MACSE codon-aware align (#i61)`.

### Task 5.2 — Substitution-model retest scaffold (bead -5b0)

Add a one-shot HPC script `scripts/hpc/test_models_gpcrtm_c20.sh` that runs `iqtree3 -m MFP -mset LG,VT,WAG,GPCRtm -madd LG+C20+R10,LG+C60+R10,LG4X+R10` on a stratified subsample of the global alignment and reports BIC. Document expected runtime; do not auto-run. Commit: `feat(model): GPCRtm + LG+C20 ModelFinder retest scaffold (#5b0)`.

### Task 5.3 — Selection stack: GARD → BUSTED-S/MH → aBSREL → MEME (bead -urk)

- Add `scripts/hpc/run_selection_stack.sh` that orchestrates the four HyPhy methods per OG.
- Update `parse_absrel.py` with parsers for BUSTED-S, BUSTED-MH, MEME outputs (separate functions; keep existing aBSREL parser).
- Use BH-FDR (Task 1.1) on the *pre-filtered* branch set: Berghia branches + LSE-internal branches.
- Commit: `feat(selection): GARD pre-screen + BUSTED-S/MH stack with FDR on filtered branches (#urk)`.

### Task 5.4 — Reconciliation: Notung → GeneRax (bead -30g)

- Add `scripts/hpc/run_generax.sh` wrapper.
- Keep `03d_notung_reconciliation.sh` as legacy fallback; new default is `03d_generax_reconciliation.sh`.
- Update `lse_refine.py` to consume GeneRax's reconciled-event annotations (duplications/transfers/losses on tree nodes) for cleaner LSE classification.
- Commit: `feat(recon): GeneRax replaces Notung as default reconciliation method (#30g)`.

### Task 5.5 — Synteny: JCVI MCscan implementation (bead -e59) — REQUIRES Phase 2 done

- New `scripts/run_jcvi_synteny.sh` invoking `python -m jcvi.compara.catalog ortholog` and `python -m jcvi.graphics.synteny`.
- Compare Berghia vs. *Aplysia californica*, *Lottia gigantea*, *Biomphalaria glabrata* (download any missing via NCBI Datasets in genome-integration step).
- `rank_candidates.py:load_synteny_scores` reads JCVI `.anchors` files, not `synteny_ids.txt`.
- Commit: `feat(synteny): JCVI MCscan against Aplysia/Lottia/Biomphalaria (#e59)`.

### Task 5.6 — TreeShrink (bead -iof)

`scripts/run_treeshrink.sh` invoked between IQ-TREE final inference and downstream selection. Default `q=0.05` per-clade. Commit: `feat(phylo): TreeShrink rogue-taxon cleaning (#iof)`.

### Task 5.7 — ASR switch: FastML → IQ-TREE --ancestral; GRASP for indels (bead -j44)

- `05_selective_pressure_and_asr.sh` calls `iqtree3 -s aln.fa -m <model_from_step04> -te tree.nwk --ancestral --prefix asr` instead of FastML.
- Optional GRASP wrapper for selected deep nodes.
- Commit: `feat(asr): IQ-TREE --ancestral default, GRASP for indels on deep nodes (#j44)`.

### Task 5.8 — Sensitivity-analysis fix (bead -j6f)

Two fixes:
1. Pass complete `base_weights` dict to `run_leave_one_out_crossval` so all 11 components are exercised.
2. Add Latin-Hypercube weight permutation (100 schemes) and report top-50 Jaccard stability across schemes.

Commit: `fix(sens): full-weight CV + LHS permutation Jaccard stability (#j6f)`.

---

## Phase 6 — P2 validation + deliverable reframing

### Task 6.1 — Aplysia retrospective validation (bead -bdu)

- New `scripts/validate_against_aplysia.sh`: runs the pipeline on the *Aplysia californica* genome (download via NCBI Datasets, accession GCF_000002075.1) using the same config; reports recall@N for known *Aplysia* rhinophore-subfamily members from Cummins 2009 (curated list in `references/cummins2009_aplysia_chemoreceptors.csv`, hand-built from Supp Tables).
- Output: `results/validation/aplysia_recall.tsv` + a one-page methods note in `docs/`.
- Commit: `feat(validate): Aplysia retrospective recall/precision evaluation (#bdu)`.

### Task 6.2 — HCR-ready deliverable reframing (bead -xqz)

- Add per-candidate **probe-design-friendliness** column to the ranked CSV: `paralog_minimum_distance_aa` (smallest pairwise distance to any paralog within ±300 aa region — proxy for HCR probe specificity), `cds_length`, `tandem_cluster_size_warning` (boolean: TRUE if cluster size ≥5 → probe specificity hard).
- Update `09_report_generation.sh` to render the top-N table with these columns.
- Add `scripts/design_hcr_targets.py` skeleton that selects top-K candidates that pass probe-friendliness filters and emits a probe-design-ready FASTA.
- Commit: `feat(hcr): probe-design-friendliness columns + HCR-target selector (#xqz)`.

### Task 6.3 — Positive-control HCR-validated genes (bead -edx)

- Curate `references/hcr_positive_controls.csv`: starting with Gα_olf (the user's example), to be appended as wet-lab data lands.
- Add `scripts/check_positive_controls.py` invoked at the end of `07_candidate_ranking.sh`: looks up each positive-control gene in the ranked CSV, reports its rank+score; raises an alert if rank-percentile drops below 50%.
- Commit: `feat(controls): positive-control sanity check on each pipeline run (#edx)`.

### Task 6.4 — Structural augmentations (bead -s6v) — DEFERRED

This is P3 and depends on AlphaFold runs. Document a script skeleton `scripts/foldseek_against_gpcrdb.sh` and `scripts/check_pocket_integrity.py`; do not implement. Commit a stub + README. (`feat(struct): scaffolds for FoldMason/GPCRdb/pocket-integrity (#s6v)`).

---

## Phase 7 — Final integration + verification

### Task 7.1 — Update environment.yml with all new pip dependencies

Add to conda deps: `gffread`, `dendropy>=5.0`. Add to pip deps: `statsmodels>=0.13`, `gffutils>=0.13`, `tmbed`, `jcvi`, `treeshrink`. Java/compiled tools (MACSE v2, GeneRax, HmmCleaner.pl, Foldseek, FoldMason) — document install in `docs/installation_hpc.md` with version pins, install commands, and SHA256 checksums where applicable.

### Task 7.2 — Run unit-test suite end-to-end

Run: `cd /home/workspace/Desktop/projects/umass/berghia-chemogpcrs && pytest tests/unit/ -v`. Expected: all green. Document any skipped tests (HPC-only).

### Task 7.3 — Run pipeline-validation harness

Run: `bash tests/validate_pipeline.sh`. Confirm step-completion flags appear for every stage.

### Task 7.4 — Update CLAUDE.md and methods_log.md

Append a section: *"May 2026 review fixes — methodology update."* List the eight high-impact changes (dN/dS asymmetry, TMbed, MACSE, GeneRax, JCVI MCscan, tandem-cluster, calibrated thresholds, Aplysia validation) with one-line rationale each.

### Task 7.5 — Close all corresponding beads issues

`bd close berghia-chemogpcrs-{ea9,wux,ryr,6nh,ce4,4nu,ar8,hg1,0ku,cio,qu9,m6k,mqt,i61,5b0,urk,30g,e59,iof,j44,j6f,bdu,xqz,edx} --reason="Implemented in 2026-05-01 review-fix sweep"`.

`bd sync` and final `git push`.

---

## Deferred to user-side HPC (out of scope for this implementation)

These cannot be run on the local 15 GB ARM machine and are flagged for the user:

1. Re-run global Berghia+refs IQ-TREE with `GPCRtm` and `LG+C20+R10` to compare against `VT+I+R10` (Task 5.2 scaffold ready).
2. Re-run aBSREL/BUSTED-S/MH on the new tree once Phase 5 selection stack is wired.
3. AlphaFold runs for top candidates (Task 6.4 stub).
4. Validation pipeline run on *Aplysia* genome (Task 6.1 — script ready, ~24 h on HPC).

For each, the user submits the SLURM script under `scripts/hpc/` and the pipeline picks up results when the artifacts land in `results/`.

---

## Acceptance criteria (whole plan)

- [ ] All Phase 1 unit tests pass.
- [ ] `bash tests/validate_pipeline.sh` completes without error using the new code paths.
- [ ] `git log --oneline` shows ~30 small, focused commits with bead references.
- [ ] `bd ready` shows zero open P0 or P1 issues from the May 2026 review.
- [ ] `methods_log.md` documents the methodology changes with citations.
- [ ] HCR-ready ranked CSV (Task 6.2) is produced from re-running stage 07 on existing preliminary data.
