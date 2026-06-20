# Genome inventory — single source of truth + how to add species

The Phase-1f BRAKER4 genome inventory lives in **one** manifest:

```
references/species_tree/genome_inventory.tsv
```

14-column superset schema (one row per species):

```
taxid  binomial  clade  policy_class  source  accession  assembly_level
annotation_status  est_protein_count  submission_date  contig_n50
total_length_bp  drop_reason  source_batch
```

`source_batch` records provenance for every row: `nath_phase1e` /
`extension_phase1d` (the two original sourcing phases, merged once by the
migration below) or `datasets_YYYYMMDD` (a later append run).

## Adding species (the one forward path)

Sourcing is driven by an **editable clade-policy config**:

```
references/species_tree/clade_policies.tsv
   clade_name   policy_class   min_assembly_level   max_count
```

**Case 1 — pick up newly-deposited genomes in clades we already cover**
(NCBI has new assemblies since the last run): just re-run the builder. It
queries each clade, dedups against the existing manifest, and **appends**
only the new species (tagged `source_batch=datasets_<today>`). Existing rows
are never modified; re-running with no new deposits changes nothing.

```bash
python scripts/build_genome_inventory.py            # appends in place
# (override the tag with --batch-tag, or the file with --manifest)
```

**Case 2 — add a new clade.** Append a row to `clade_policies.tsv`
(`clade_name`, `policy_class`, `min_assembly_level` — `Chromosome`/`Scaffold`/
empty for any — and optional `max_count` cap), then re-run the builder as in
Case 1.

After either case, download the new genomes and feed the pipeline:

```bash
python scripts/download_species_tree_phase1f_genomes.py   # reads genome_inventory.tsv
python scripts/build_braker4_samples_csv.py --manifest references/species_tree/genome_inventory.tsv \
    --genome-cache species_tree_data/braker4_genomes \
    --protein-db <evidence DB> --out <samples_full.csv>
# then the Phase-1f BRAKER4 array (one task per species)
```

## How the single manifest was created (one-time migration)

`genome_inventory.tsv` was produced by losslessly merging the two original
manifests — `genome_inventory_unannotated.tsv` (M1, Phase 1e per-taxid) and
`extension_inventory.tsv` (M2, Phase 1d per-clade) — with:

```bash
python scripts/migrate_genome_inventory.py    # M1 + M2 -> genome_inventory.tsv
```

"Lossless" is enforced by a test: the unified manifest yields the identical
species set (`build_braker4_samples_csv.read_targets`) as the old two-manifest
build. Verified on the real manifests (450 == 450 identical targets).

## Superseded artifacts (kept for provenance, not used by the pipeline)

- `genome_inventory_unannotated.tsv` (M1) and `extension_inventory.tsv` (M2)
  are the original source manifests. They are inputs to the one-time migration
  only — the pipeline reads `genome_inventory.tsv`. They may be moved under
  `references/species_tree/archive/` at deployment.
- `scripts/build_species_tree_phase1d_extension_inventory.py` is retained as
  the **selection-core library** (`ClaadePolicy`, `query_datasets_for_clade`,
  `select_for_clade`, `build_extension_inventory`, `load_existing_taxids`,
  `ExtensionEntry`) that `build_genome_inventory.py` imports. Its standalone
  `main()` (which wrote the separate M2 file) is superseded by
  `build_genome_inventory.py`.
- The old `sbatch_run_species_tree_phase1d.sh` / `phase1e.sh` wrappers are
  superseded by `sbatch_run_species_tree_genome_inventory.sh` (Unity).

Design + rationale: `docs/plans/2026-06-19-unified-genome-inventory-sourcing-design.md` (local).
