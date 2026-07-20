"""Reviewed inventory of seam probes with no producer inside the repository.

Every entry here was checked BY HAND against the actual producer. An entry is
a statement that the canary's "no writer found" verdict is understood, not a
way to silence it. Two tests keep the inventory honest:

  * a probe that is unresolved and NOT listed here fails the canary;
  * a listed probe that becomes resolvable fails the canary as a STALE entry,
    so the inventory shrinks as the pipeline is wired up and can never quietly
    accumulate cover for real breakage.

Categories
----------
EXTERNAL_TOOL      A third-party binary creates the path (OrthoFinder, JCVI,
                   MCScanX, AlphaFold, MMseqs, pdflatex). No repo code writes
                   it, and no repo change can.
ONE_TIME_PREP      Produced by a documented one-shot prep step that is run out
                   of band, not by a numbered stage (see the project setup docs).
OFF_REPO           Produced by scripts that live on the cluster, not here.
OPTIONAL_OVERRIDE  A hand-staged file the pipeline deliberately never writes;
                   the consumer probes it and degrades to a documented fallback.
EXTRACTOR_LIMIT    A real producer exists in this repo but static extraction
                   cannot see it (positional CLI args, indirect ``${!VAR}``
                   expansion). The seam is fine; the extractor is not.
SUSPECTED_DEFECT   No producer found anywhere. These are live findings of the
                   exact class this canary exists to catch.

The RESULTS placeholder stands in for the configured RESULTS_DIR so the
inventory is independent of where the repository is checked out.
"""

from __future__ import annotations

EXTERNAL_TOOL = "EXTERNAL_TOOL"
ONE_TIME_PREP = "ONE_TIME_PREP"
OFF_REPO = "OFF_REPO"
EXTRACTOR_LIMIT = "EXTRACTOR_LIMIT"
OPTIONAL_OVERRIDE = "OPTIONAL_OVERRIDE"
SUSPECTED_DEFECT = "SUSPECTED_DEFECT"

# (consuming file, probe pattern) -> (category, why)
KNOWN_GAPS: dict[tuple[str, str], tuple[str, str]] = {
    # --- one-shot preprocessing (the documented pre-pipeline setup steps)
    ("01_reference_processing.sh", "${RESULTS}/reference_sequences/cds/lse"): (
        ONE_TIME_PREP,
        "scripts/run_cds_preprocess.sh recovers reference CDS via miniprot "
        "before the numbered stages run.",
    ),
    (
        "01_reference_processing.sh",
        "${RESULTS}/reference_sequences/cds/one_to_one_ortholog",
    ): (
        ONE_TIME_PREP,
        "Same one-shot CDS recovery step; conserved-ortholog half of the split.",
    ),
    ("06c_classify_non_chemoreceptors.sh", "${RESULTS}/classification/hmms/pfam_fallback"): (
        ONE_TIME_PREP,
        "scripts/build_classification_hmms.py, run once during classification-reference setup.",
    ),
    ("06c_classify_non_chemoreceptors.sh", "${RESULTS}/classification/loo/loo_metrics.tsv"): (
        ONE_TIME_PREP,
        "scripts/validate_classification_hmms.py leave-one-out validation, "
        "run once during classification-reference setup.",
    ),
    ("06c_classify_non_chemoreceptors.sh", "${RESULTS}/classification/trees/backbone.treefile"): (
        ONE_TIME_PREP,
        "scripts/unity/build_classification_reference_trees.sh, run once during classification-reference setup.",
    ),
    ("04_phylogenetic_analysis.sh", "${RESULTS}/p5_phase1a_validation/pools/pool_members_class_*.tsv"): (
        ONE_TIME_PREP,
        "scripts/build_per_class_reference_pools.py via "
        "scripts/unity/run_per_class_pools.sh (per-class pool build).",
    ),

    # --- external binaries own the layout
    # NOTE: the 03_orthology_clustering.sh entry for this same path was REMOVED
    # once the run-selection rule converged on a shared resolver -- the probe now
    # has a producer, so the exemption was stale. The stale-entry guard caught it.
    ("03c_cafe_analysis.sh", "${RESULTS}/orthogroups/input/OrthoFinder"): (
        EXTERNAL_TOOL,
        "Same OrthoFinder-owned directory, probed by the CAFE stage.",
    ),
    (
        "03c_cafe_analysis.sh",
        "${RESULTS}/orthogroups/input/OrthoFinder/Results_*/Orthogroups/Orthogroups.GeneCount.tsv",
    ): (
        EXTERNAL_TOOL,
        "OrthoFinder gene-count table. NOTE: this probe is a three-way "
        "fallback and the canary treats the alternatives as an OR-group, so "
        "ALL THREE alternatives must be listed for the group to be exempt.",
    ),
    (
        "03c_cafe_analysis.sh",
        "${RESULTS}/orthogroups/input/OrthoFinder/Orthogroups.GeneCount.tsv",
    ): (
        EXTERNAL_TOOL,
        "Second alternative of the same OrthoFinder gene-count fallback "
        "chain (older OrthoFinder layout, no Results_* level).",
    ),
    (
        "03c_cafe_analysis.sh",
        "${RESULTS}/orthogroups/input/OrthoFinder/Results_*/Orthogroups.GeneCount.tsv",
    ): (
        EXTERNAL_TOOL,
        "Third alternative of the same fallback chain (gene-count table at "
        "the Results_* root rather than under Orthogroups/).",
    ),
    ("scripts/rank_candidates.py", "${RESULTS}/orthogroups/input/OrthoFinder/Results_*/Resolved_Gene_Trees"): (
        EXTERNAL_TOOL,
        "OrthoFinder writes Resolved_Gene_Trees itself. The ROOT is correct "
        "here -- this is the post-fix path; the pre-fix probe used "
        "results/orthofinder/ and is what the canary catches.",
    ),
    ("06_synteny_and_mapping.sh", "${RESULTS}/synteny/jcvi/berghia_vs_*/berghia.*.lifted.anchors"): (
        EXTERNAL_TOOL,
        "JCVI MCscan writes .lifted.anchors into its own comparison dirs.",
    ),
    ("06_synteny_and_mapping.sh", "${RESULTS}/synteny/*.collinearity"): (
        EXTERNAL_TOOL,
        "MCScanX writes .collinearity next to its input.",
    ),
    ("08_structural_analysis.sh", "${RESULTS}/structural_analysis/alphafold/*/json"): (
        EXTERNAL_TOOL,
        "AlphaFold3 output directory layout is created by the predictor.",
    ),
    ("09_report_generation.sh", "${RESULTS}/report/report.pdf"): (
        EXTERNAL_TOOL,
        "pdflatex writes the PDF from a working directory, so no shell-level "
        "path to results/report/report.pdf appears.",
    ),
    ("01_reference_processing.sh", "${RESULTS}/hmms/*_orthogroups/OG_*.fa"): (
        EXTERNAL_TOOL,
        "Clustering tool writes OG_<category>* members into the HMM dir using "
        "a bare --prefix, so the full path is never spelled out.",
    ),

    # --- produced on the cluster, outside this checkout
    ("07_candidate_ranking.sh", "${RESULTS}/ranking/diagnostics/novelty_esm2_3b_PROD.tsv"): (
        OFF_REPO,
        "PLM novelty tables are produced by scratch scripts on Unity "
        "(see the embedding-channel notes); the stage consumes them if present.",
    ),
    ("07_candidate_ranking.sh", "${RESULTS}/ranking/diagnostics/novelty_protrek_PROD.tsv"): (
        OFF_REPO,
        "Same PLM novelty channel, ProTrek model.",
    ),

    # --- real producer, invisible to static extraction
    ("03b_lse_classification.sh", "${RESULTS}/reference_sequences/id_map.csv"): (
        EXTRACTOR_LIMIT,
        "scripts/update_headers.py writes it via a POSITIONAL argv[2] "
        "argument, not a named --output flag, so the two-hop resolver cannot "
        "bind the path.",
    ),
    ("04_phylogenetic_analysis.sh", "${RESULTS}/p5_phase1a_validation/pools/refs_class_OUTGROUP_SOURCE_CLASS_*.fa"): (
        EXTRACTOR_LIMIT,
        "Path uses bash INDIRECT expansion ${!OUTGROUP_SOURCE_VAR}, which "
        "cannot be resolved without running the stage.",
    ),
    ("07_candidate_ranking.sh", "${RESULTS}/ranking/signal_independence_groups.json"): (
        EXTRACTOR_LIMIT,
        "rank_candidates.py writes it as a sibling of its --output argument "
        "(Path(output_path).parent / ...), and the stage builds that output "
        "path through a variable the extractor cannot follow.",
    ),
    ("06c_classify_non_chemoreceptors.sh", "${RESULTS}/classification/nath_classifications.tsv"): (
        EXTRACTOR_LIMIT,
        "Written by the OG-vote classifier block inside 06c through a shell "
        "variable that is assigned by command substitution.",
    ),
    ("06c_classify_non_chemoreceptors.sh", "${RESULTS}/classification/candidate_placement.tsv"): (
        EXTRACTOR_LIMIT,
        "scripts/classify_via_placement.py --output-tsv writes it; the flag "
        "value is a shell variable assigned outside the invocation line.",
    ),

    # --- optional hand-staged inputs, deliberately unwritten
    ("scripts/classify_gpcr_by_class.py", "${RESULTS}/classification/hmms/pfam_fallback"): (
        ONE_TIME_PREP,
        "Pfam fallback HMM directory built by scripts/build_classification_hmms.py "
        "during one-time prep; same directory 06c probes.",
    ),
    ("scripts/rank_candidates.py", "${RESULTS}/reference_sequences/ref_categories_final.csv"): (
        OPTIONAL_OVERRIDE,
        "Optional explicit reference-category map. rank_candidates.py probes it "
        "and falls back to loaded JSON categories, then a keyword heuristic. "
        "Deliberately hand-staged; nothing in the pipeline writes it.",
    ),
    ("scripts/rank_candidates.py", "${RESULTS}/orthology/gene_orthogroups.csv"): (
        OPTIONAL_OVERRIDE,
        "Legacy hand-supplied gene->orthogroup mapping. The source comment says "
        "outright that nothing in the pipeline writes it; the real path comes "
        "from find_orthogroups_tsv(), which resolves.",
    ),
    (
        "07_candidate_ranking.sh",
        "${RESULTS}/ranking/embeddings/candidate_ref_identity_PROD.tsv",
    ): (
        OPTIONAL_OVERRIDE,
        "candidate->nearest-reference %identity, from an mmseqs/diamond search "
        "run out of band on Unity (see embedding_candidate_diagnostics."
        "_read_identity). Nothing in this repo writes it, and it has never been "
        "generated in this checkout (observed on job 61993481). It is the "
        "OPTIONAL 4th field of a fusion_consensus --models spec, feeding only "
        "the identity_to_nearest CONFOUND -- never the novelty score. The probe "
        "is the documented degradation: present -> spec carries the confound; "
        "absent -> the spec is built without it, every model still scores, and "
        "the omission is logged. Before that guard existed the path was appended "
        "unconditionally, _read_identity raised FileNotFoundError on the missing "
        "file, and stage 07's `|| log --level=WARN` swallowed it -- taking the "
        "whole embedding channel dormant. The same guard already exists in "
        "scripts/unity/rebuild_embedding_channel.sh.",
    ),

    # --- LIVE FINDINGS -----------------------------------------------------
    ("05_selective_pressure_and_asr.sh", "${RESULTS}/calibration/depth_thresholds.json"): (
        SUSPECTED_DEFECT,
        "scripts/calibrate_depth_threshold.py is the ONLY writer of this file "
        "and is invoked by NO stage script. Stage 05's `[ -f ]` guard can "
        "therefore never be true, so the null-calibrated ASR depth threshold "
        "is silently never used and MIN_ASR_DISTANCE is always taken instead. "
        "This is the same orphaned-producer shape as the interpret_cafe.py "
        "defect. Reported, not fixed: stage files are out of scope here.",
    ),
}


def gap_key(stage: str, pattern: str, results_dir: str) -> tuple[str, str]:
    """Inventory key for a probe, with RESULTS_DIR abstracted away."""
    return (stage, pattern.replace(results_dir, "${RESULTS}"))
