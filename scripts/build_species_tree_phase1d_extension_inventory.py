"""Phase 1d extension inventory: widen the species tree beyond Nath 239.

Bead -pcc (Audit C of -v1c epic, under -dnk umbrella).

Decision 2026-05-21 (tight-scope policy, locked):
  - Heterobranchia ALL (Berghia's subclass) — any assembly quality
  - rare basals (Polyplacophora / Scaphopoda / Caudofoveata /
    Monoplacophora) — any quality
  - other Mollusca subclasses (Bivalvia / Cephalopoda /
    Caenogastropoda / Vetigastropoda / Neritimorpha /
    Patellogastropoda) — CHROMOSOME-LEVEL only
  - outgroup phyla (Annelida / Platyhelminthes / Nemertea /
    Brachiopoda / Phoronida / Bryozoa) — CHROMOSOME-LEVEL with
    per-phylum cap

The script does NOT download anything; it queries NCBI Datasets for
each clade and writes a manifest TSV. Genome download + BRAKER4
annotation happen later (Phase 1f reuses the same plumbing).

Dedup target: union of taxids in the existing Nath references
(.faa filenames) + Phase 1a manifest + Phase 1e manifest.
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Union

from build_species_tree_phase1a_inventory import (
    _ASSEMBLY_LEVEL_RANK,
    AssemblyChoice,
    parse_reference_filename,
    pick_best_assembly,
)


# ----------------------------------------------------------------------
# Existing-taxid dedup source
# ----------------------------------------------------------------------

def load_existing_taxids(
    nath_refs_root: Union[str, Path],
    *manifest_paths: Union[str, Path],
) -> set[int]:
    """Return the union of taxids already covered by:
      - nath_et_al/<clade>/*.faa filename prefix
      - any number of Phase-1a/1e manifest TSVs (`taxid` column)

    Missing refs_root or missing manifests are skipped silently so the
    function works with partial inputs (e.g., when invoked locally
    against repo-tracked manifests only).
    """
    taxids: set[int] = set()

    refs = Path(nath_refs_root)
    if refs.exists():
        for clade_dir in refs.iterdir():
            if not clade_dir.is_dir():
                continue
            for faa in clade_dir.glob("*.faa"):
                try:
                    taxid, _ = parse_reference_filename(faa.name)
                except ValueError:
                    continue
                taxids.add(taxid)

    for mpath in manifest_paths:
        p = Path(mpath)
        if not p.exists():
            continue
        with p.open() as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                try:
                    taxids.add(int(row["taxid"]))
                except (KeyError, TypeError, ValueError):
                    continue

    return taxids


# ----------------------------------------------------------------------
# NCBI Datasets clade query
# ----------------------------------------------------------------------

# Same set of "no hit" stderr phrases Phase 1a tolerates. Some clades
# (Caudofoveata, Monoplacophora) genuinely have 0 chromosome-level
# assemblies; treat that as empty result, not an error.
_NO_HIT_PHRASES = (
    "no assemblies",
    "no results",
    "no genome data",
    "no genomes",
    "no annotation",
    "no proteins",
    "did not match any genomes",
)


def query_datasets_for_clade(
    clade_name: str,
    datasets_bin: str = "datasets",
    timeout: int = 600,
) -> list[dict]:
    """Call `datasets summary genome taxon "<clade_name>" --as-json-lines`
    and return parsed JSON records. Empty list when the clade has no hits.

    Raises RuntimeError on hard failures (network, auth, malformed CLI).
    """
    cmd = [
        datasets_bin, "summary", "genome", "taxon", clade_name,
        "--as-json-lines",
    ]
    result = subprocess.run(
        cmd, capture_output=True, text=True, timeout=timeout,
    )
    stderr_lower = (result.stderr or "").lower()

    if result.returncode != 0:
        if any(p in stderr_lower for p in _NO_HIT_PHRASES):
            return []
        raise RuntimeError(
            f"datasets failed (returncode={result.returncode}) "
            f"for clade={clade_name!r}: {result.stderr.strip()}"
        )

    if not result.stdout.strip() and any(p in stderr_lower for p in _NO_HIT_PHRASES):
        return []

    records: list[dict] = []
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            records.append(json.loads(line))
        except json.JSONDecodeError as e:
            print(
                f"warning: skipping malformed JSON line for clade={clade_name!r}: {e}",
                file=sys.stderr,
            )
    return records


# ----------------------------------------------------------------------
# Per-record helpers
# ----------------------------------------------------------------------

def group_records_by_taxid(records: list[dict]) -> dict[int, list[dict]]:
    """Group flat NCBI Datasets records by organism.tax_id. Drop records
    missing a numeric tax_id."""
    groups: dict[int, list[dict]] = {}
    for r in records:
        org = r.get("organism") or {}
        try:
            taxid = int(org.get("tax_id"))
        except (TypeError, ValueError):
            continue
        groups.setdefault(taxid, []).append(r)
    return groups


def extract_binomial(records: list[dict]) -> str:
    """Best-effort organism_name from a per-taxid record list. Returns
    empty string when no record carries a name."""
    for r in records:
        org = r.get("organism") or {}
        name = (org.get("organism_name") or "").strip()
        if name:
            return name
    return ""


def apply_assembly_level_filter(
    choice: AssemblyChoice | None,
    min_level: str,
) -> AssemblyChoice | None:
    """Drop choices below min_level. Empty min_level => pass-through."""
    if choice is None:
        return None
    if not min_level:
        return choice
    threshold = _ASSEMBLY_LEVEL_RANK.get(min_level, 0)
    have = _ASSEMBLY_LEVEL_RANK.get(choice.assembly_level, 0)
    return choice if have >= threshold else None


# ----------------------------------------------------------------------
# Policy
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class ClaadePolicy:
    """Selection policy for one NCBI taxonomy clade query.

    Fields:
      clade_name:          NCBI taxonomy name (e.g. "Heterobranchia").
      policy_class:        Short label propagated into the manifest's
                           policy_class column (e.g. "heterobranchia").
      min_assembly_level:  "Chromosome" enforces chromosome-or-better;
                           "" lets any level through.
      require_annotation:  False for Phase 1d — we'll BRAKER4-annotate
                           later, so we don't filter out unannotated
                           GenBank assemblies (they're our main target).
      max_count:           Per-clade cap (used for outgroup sampling).
                           None = no cap.
    """
    clade_name: str
    policy_class: str
    min_assembly_level: str
    require_annotation: bool
    max_count: int | None


@dataclass(frozen=True)
class ExtensionEntry:
    taxid: int
    binomial: str
    policy_class: str
    clade_name: str
    choice: AssemblyChoice


# ----------------------------------------------------------------------
# Per-clade selection
# ----------------------------------------------------------------------

def select_for_clade(
    records: list[dict],
    policy: ClaadePolicy,
    exclude_taxids: set[int],
) -> list[ExtensionEntry]:
    """Group records by taxid, pick the best assembly per species,
    apply the policy's quality filter, drop already-covered taxids,
    apply the per-policy cap.

    Cap ordering: by assembly_level rank descending, then taxid
    ascending for a stable tiebreak.
    """
    groups = group_records_by_taxid(records)
    scored: list[tuple[ExtensionEntry, int]] = []
    for taxid, group_records in groups.items():
        if taxid in exclude_taxids:
            continue
        choice = pick_best_assembly(
            group_records, require_annotation=policy.require_annotation,
        )
        choice = apply_assembly_level_filter(choice, policy.min_assembly_level)
        if choice is None:
            continue
        entry = ExtensionEntry(
            taxid=taxid,
            binomial=extract_binomial(group_records),
            policy_class=policy.policy_class,
            clade_name=policy.clade_name,
            choice=choice,
        )
        quality_rank = _ASSEMBLY_LEVEL_RANK.get(choice.assembly_level, 0)
        scored.append((entry, quality_rank))

    scored.sort(key=lambda x: (-x[1], x[0].taxid))
    if policy.max_count is not None:
        scored = scored[: policy.max_count]
    return [e for e, _ in scored]


# ----------------------------------------------------------------------
# Orchestrator
# ----------------------------------------------------------------------

def build_extension_inventory(
    policies: list[ClaadePolicy],
    exclude_taxids: set[int],
    query_fn: Callable[[str], list[dict]] = query_datasets_for_clade,
    progress_fh=sys.stderr,
) -> list[ExtensionEntry]:
    """Iterate policies in order, query each clade once, apply
    `select_for_clade`. A taxid picked up under an earlier policy is
    excluded from subsequent policies — first-policy-wins on labelling.
    """
    seen: set[int] = set(exclude_taxids)
    all_entries: list[ExtensionEntry] = []
    for i, policy in enumerate(policies, start=1):
        try:
            records = query_fn(policy.clade_name)
        except RuntimeError as e:
            print(
                f"  [{i}/{len(policies)}] {policy.clade_name}: query failed "
                f"({e}) — skipping clade",
                file=progress_fh,
            )
            continue
        entries = select_for_clade(records, policy, seen)
        all_entries.extend(entries)
        seen.update(e.taxid for e in entries)
        print(
            f"  [{i}/{len(policies)}] {policy.clade_name:25s} "
            f"policy={policy.policy_class:20s} "
            f"records={len(records):4d}  selected={len(entries):4d}",
            file=progress_fh,
        )
    return all_entries


# ----------------------------------------------------------------------
# Manifest writer
# ----------------------------------------------------------------------

EXTENSION_COLUMNS = (
    "taxid",
    "binomial",
    "policy_class",
    "clade_name",
    "source",
    "accession",
    "assembly_level",
    "annotation_status",
    "est_protein_count",
    "submission_date",
    "contig_n50",
    "total_length_bp",
)


def write_extension_manifest(
    path: Union[str, Path],
    entries: list[ExtensionEntry],
) -> None:
    """Write a TSV with EXTENSION_COLUMNS header. Rows sorted by
    (policy_class, taxid) for stable reviewing."""
    rows = sorted(entries, key=lambda e: (e.policy_class, e.taxid))
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        f.write("\t".join(EXTENSION_COLUMNS) + "\n")
        for e in rows:
            c = e.choice
            f.write("\t".join([
                str(e.taxid),
                e.binomial,
                e.policy_class,
                e.clade_name,
                c.source,
                c.accession,
                c.assembly_level,
                c.annotation_status,
                str(c.est_protein_count),
                c.submission_date,
                str(c.contig_n50),
                str(c.total_length_bp),
            ]) + "\n")


# ----------------------------------------------------------------------
# Hardcoded tight-scope policy list
# ----------------------------------------------------------------------

# Per 2026-05-21 species-set expansion policy decision (memory:
# feedback_species_set_policy). Order matters: earlier entries claim
# taxids first, so list mollusca subclasses before "outgroup" phyla
# and Heterobranchia before generic Gastropoda subclasses.
POLICIES: list[ClaadePolicy] = [
    # ----- Mollusca subclasses (comprehensive) -----
    ClaadePolicy("Heterobranchia",     "heterobranchia",  "",            False, None),
    ClaadePolicy("Polyplacophora",     "rare_basal",      "",            False, None),
    ClaadePolicy("Scaphopoda",         "rare_basal",      "",            False, None),
    ClaadePolicy("Caudofoveata",       "rare_basal",      "",            False, None),
    ClaadePolicy("Monoplacophora",     "rare_basal",      "",            False, None),
    # ----- Other mollusca subclasses, chromosome-only -----
    ClaadePolicy("Bivalvia",           "other_mollusca",  "Chromosome",  False, None),
    ClaadePolicy("Cephalopoda",        "other_mollusca",  "Chromosome",  False, None),
    ClaadePolicy("Caenogastropoda",    "other_mollusca",  "Chromosome",  False, None),
    ClaadePolicy("Vetigastropoda",     "other_mollusca",  "Chromosome",  False, None),
    ClaadePolicy("Neritimorpha",       "other_mollusca",  "Chromosome",  False, None),
    ClaadePolicy("Patellogastropoda",  "other_mollusca",  "Chromosome",  False, None),
    # ----- Outgroup phyla, chromosome-only with per-phylum cap -----
    ClaadePolicy("Annelida",           "outgroup_annelida",         "Chromosome", False, 3),
    ClaadePolicy("Platyhelminthes",    "outgroup_platyhelminthes",  "Chromosome", False, 3),
    ClaadePolicy("Nemertea",           "outgroup_nemertea",         "Chromosome", False, 3),
    ClaadePolicy("Brachiopoda",        "outgroup_brachiopoda",      "Chromosome", False, 3),
    ClaadePolicy("Phoronida",          "outgroup_phoronida",        "Chromosome", False, 3),
    ClaadePolicy("Bryozoa",            "outgroup_bryozoa",          "Chromosome", False, 3),
]


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------

def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Phase 1d extension inventory. Queries NCBI Datasets per "
            "clade under the tight-scope policy and writes a manifest "
            "TSV of NEW species to add to the species tree."
        )
    )
    p.add_argument(
        "--refs-root",
        type=Path,
        default=Path("references/nath_et_al/one_to_one_ortholog"),
        help="Existing Nath references dir (for dedup)",
    )
    p.add_argument(
        "--phase1a-manifest",
        type=Path,
        default=Path("references/species_tree/proteome_manifest.tsv"),
        help="Phase 1a proteome manifest (for dedup)",
    )
    p.add_argument(
        "--phase1e-manifest",
        type=Path,
        default=Path("references/species_tree/genome_inventory_unannotated.tsv"),
        help="Phase 1e unannotated-genome manifest (for dedup)",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=Path("references/species_tree/extension_inventory.tsv"),
        help="Output TSV path",
    )
    p.add_argument(
        "--datasets-bin",
        default="datasets",
        help="Path to NCBI Datasets CLI binary (default: 'datasets' on PATH)",
    )
    p.add_argument(
        "--timeout",
        type=int,
        default=600,
        help="Per-clade query timeout in seconds (clade queries return more "
             "records than per-taxon queries; default 10 minutes)",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_argparser().parse_args(argv)

    exclude = load_existing_taxids(
        args.refs_root, args.phase1a_manifest, args.phase1e_manifest,
    )
    print(
        f"Dedup set: {len(exclude)} existing taxids loaded from refs + manifests",
        file=sys.stderr,
    )

    def query_fn(clade_name: str) -> list[dict]:
        return query_datasets_for_clade(
            clade_name, datasets_bin=args.datasets_bin, timeout=args.timeout,
        )

    entries = build_extension_inventory(
        policies=POLICIES, exclude_taxids=exclude, query_fn=query_fn,
    )
    write_extension_manifest(args.out, entries)

    print(
        f"\nExtension inventory: {len(entries)} NEW species written to {args.out}",
        file=sys.stderr,
    )
    # Per-policy summary
    from collections import Counter
    per_class = Counter(e.policy_class for e in entries)
    for pc, n in sorted(per_class.items()):
        print(f"  {pc:30s}  {n:4d}", file=sys.stderr)

    return 0 if entries else 1


if __name__ == "__main__":
    sys.exit(main())
