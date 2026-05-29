#!/usr/bin/env python3
"""
build_per_class_reference_pools.py

P2 of the per-class refactor.

Consumes:
  - Per-species chemo_candidates.fa files from a P0 scan (glob pattern)
  - A P1 class TSV mapping seq_id → {A,B,C,F,unclassified}
  - (optional) Berghia candidate FASTA + Berghia P1 class TSV

Produces, under --out-dir:
  - refs_class_A.fa  (≤ max-class-A sequences, default 2000)
  - refs_class_B.fa  (≤ max-class-B sequences, default 2900)
  - refs_class_C.fa  (≤ max-class-C sequences, default 2900)
  - refs_class_F.fa  (≤ max-class-F sequences, default 2900)
  - unclassified_log.tsv
  - pool_build_report.json

Design constraints (locked, 2026-05-28):
  - Berghia candidates are classified by P1 into per-class labels.
  - Each per-class pool includes Berghia sequences of that class as
    MUST_INCLUDE (always retained, not subsampled away).
  - MUST_INCLUDE taxids are always retained (their candidates fill first).
  - Remaining slots are filled by Berghia-proximity-weighted sampling
    of the non-must-include candidates.
  - CD-HIT 0.7 deduplication before subsampling.
  - Class A budget reduced (~2000) to account for ~800 Berghia Class A seqs.
  - Class B/C/F budget at ~2900 (no Berghia tax for those classes).

Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import random
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DEFAULT_MUST_INCLUDE_TAXIDS: frozenset[int] = frozenset({
    # Mollusc proximity (highest priority)
    6500,    # Aplysia californica (sea slug, chemoreception model)
    225164,  # Lottia gigantea (limpet, mollusc reference)
    29159,   # Crassostrea gigas (oyster)
    37653,   # Octopus bimaculoides
    6645,    # Octopus vulgaris
    6526,    # Biomphalaria glabrata (snail)
    6523,    # Lymnaea stagnalis (pond snail)
    # Chemoreceptor model organisms
    7227,    # Drosophila melanogaster
    7165,    # Anopheles gambiae
    7460,    # Apis mellifera
    6239,    # Caenorhabditis elegans
    7091,    # Bombyx mori
    # General vertebrate models
    10090,   # Mus musculus
    9606,    # Homo sapiens
    7955,    # Danio rerio
    # Deep outgroups
    7668,    # Strongylocentrotus purpuratus (sea urchin)
    45351,   # Nematostella vectensis (cnidarian)
    10228,   # Trichoplax adhaerens (placozoan)
})

BERGHIA_TAXID_DEFAULT = 1287507

# Default total budget per class: (Berghia count + reference count) per class
DEFAULT_TOTAL_BUDGET_PER_CLASS = 3000

# Legacy default per-class size caps (deprecated; kept for back-compat only).
# Class A is lower because Berghia contributes ~800 Class A candidates.
# Class B/C/F have no Berghia tax, so they get the full 2900-slot budget.
# Note: these are ONLY used if --max-class-A/B/C/F are explicitly passed,
# which triggers a deprecation warning.
DEFAULT_MAX_CLASS_A = 2000
DEFAULT_MAX_CLASS_B = 2900
DEFAULT_MAX_CLASS_C = 2900
DEFAULT_MAX_CLASS_F = 2900

# Berghia lineage (NCBI): from root toward Berghia stephanieae.
# Used to compute proximity scores for other taxa.
# Computed lazily via get_lineage() on first call.
_BERGHIA_LINEAGE_CACHE: Optional[list[int]] = None

GPCR_CLASSES = ("A", "B", "C", "F")


# ---------------------------------------------------------------------------
# Taxonomy helpers
# ---------------------------------------------------------------------------

def _init_ncbi_taxa():
    """Lazily initialise ete3.NCBITaxa (downloads ~10 MB DB on first run)."""
    try:
        from ete3 import NCBITaxa
        return NCBITaxa()
    except ImportError:
        print(
            "WARNING: ete3 not available; proximity scores will all be 0.",
            file=sys.stderr,
        )
        return None


_NCBI_TAXA = None  # initialised on first get_lineage() call


def get_lineage(taxid: int) -> list[int]:
    """Return the NCBI lineage (list of taxids, root to taxid) for a given taxid.

    In production: uses ete3.NCBITaxa (lazy-init).
    In tests: replaced via unittest.mock.patch.
    """
    global _NCBI_TAXA
    if _NCBI_TAXA is None:
        _NCBI_TAXA = _init_ncbi_taxa()
    if _NCBI_TAXA is None:
        return [1, taxid]
    try:
        lineage = _NCBI_TAXA.get_lineage(taxid)
        return lineage if lineage else [1, taxid]
    except Exception:
        return [1, taxid]


def proximity_score(taxid: int, lineage_fn, berghia_taxid: int) -> int:
    """Depth of the lowest common ancestor between *taxid* and *berghia_taxid*.

    Higher = closer to Berghia.  Two taxa that share a large fraction of
    their lineage with Berghia score higher than distant ones.

    Args:
        taxid: taxid to score.
        lineage_fn: callable(taxid) → list[int] — returns lineage root-to-leaf.
        berghia_taxid: the focal taxid (excluded from pools; used as reference).

    Returns:
        int — LCA depth (number of shared ancestors from root).
    """
    berghia_lineage = lineage_fn(berghia_taxid)
    other_lineage = lineage_fn(taxid)

    berghia_set = set(berghia_lineage)
    depth = 0
    for node in other_lineage:
        if node in berghia_set:
            depth += 1
    return depth


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def taxid_from_filename(filename: str) -> Optional[int]:
    """Extract the taxid integer from a scan FASTA filename.

    Expected patterns:
        <taxid>_<Genus>_<species>.chemo_candidates.fa
        <taxid>.chemo_candidates.fa
        <taxid>_<anything>.fa

    Returns None if the leading token is not numeric.
    """
    stem = Path(filename).name
    # Strip known suffixes iteratively
    for suffix in (".chemo_candidates.fa", ".fa", ".faa", ".fasta"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    first_token = stem.split("_")[0]
    if first_token.isdigit():
        return int(first_token)
    return None


def load_class_tsv(path: str) -> dict[str, str]:
    """Parse the P1 class TSV into {seq_id → class_label}.

    Expected columns: seq_id, class, evidence_pfam, top_evalue.
    Extra columns are ignored.
    """
    mapping: dict[str, str] = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            mapping[row["seq_id"]] = row["class"]
    return mapping


def load_must_include_taxids(tsv_path: str) -> frozenset[int]:
    """Load a custom must-include taxid list from a TSV with a 'taxid' column."""
    taxids: set[int] = set()
    with open(tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                taxids.add(int(row["taxid"]))
            except (KeyError, ValueError):
                pass
    return frozenset(taxids)


def load_berghia_candidates(
    berghia_fasta: str,
    berghia_class_tsv: str,
    berghia_taxid: int,
) -> dict[str, list[SeqRecord]]:
    """Load Berghia candidate sequences classified by P1 into per-class bins.

    Args:
        berghia_fasta: path to FASTA of Berghia chemoreceptor candidates.
        berghia_class_tsv: path to P1 class TSV for those candidates.
        berghia_taxid: taxid used for identity tracking (informational only).

    Returns:
        dict mapping class label (A/B/C/F) → list of SeqRecord.
        Unclassified sequences are silently dropped (they add no tree value).
    """
    class_map = load_class_tsv(berghia_class_tsv)
    per_class: dict[str, list[SeqRecord]] = {cls: [] for cls in GPCR_CLASSES}

    for rec in SeqIO.parse(berghia_fasta, "fasta"):
        seq_class = class_map.get(rec.id)
        if seq_class in per_class:
            per_class[seq_class].append(rec)

    total = sum(len(v) for v in per_class.values())
    print(
        f"[build_per_class_reference_pools] Loaded {total} Berghia candidates "
        f"(taxid {berghia_taxid}) from {berghia_fasta}",
        file=sys.stderr,
    )
    return per_class


# ---------------------------------------------------------------------------
# CD-HIT deduplication
# ---------------------------------------------------------------------------

def cdhit_dedup(
    records: list[tuple[int, SeqRecord]],
    identity: float = 0.7,
    threads: int = 4,
    cdhit_path: str = "cd-hit",
) -> list[tuple[int, SeqRecord]]:
    """Run CD-HIT on *records* and return representatives.

    Args:
        records: list of (taxid, SeqRecord) pairs.
        identity: clustering threshold.
        threads: CPU threads for CD-HIT.
        cdhit_path: path to cd-hit binary.

    Returns:
        Filtered list of (taxid, SeqRecord) — one representative per cluster.
    """
    if not records:
        return []

    word_size = 5 if identity >= 0.7 else (4 if identity >= 0.6 else 3)

    with tempfile.TemporaryDirectory() as tmpdir:
        input_fa = os.path.join(tmpdir, "input.fa")
        output_fa = os.path.join(tmpdir, "output.fa")

        # Write input (build id → (taxid, record) map for fast lookup)
        id_to_pair: dict[str, tuple[int, SeqRecord]] = {}
        with open(input_fa, "w") as fh:
            for taxid, rec in records:
                fh.write(f">{rec.id}\n{str(rec.seq)}\n")
                id_to_pair[rec.id] = (taxid, rec)

        cmd = [
            cdhit_path,
            "-i", input_fa,
            "-o", output_fa,
            "-c", str(identity),
            "-n", str(word_size),
            "-T", str(threads),
            "-M", "0",   # no memory limit (let OS decide)
            "-d", "0",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(
                f"CD-HIT failed (return code {result.returncode}), "
                f"keeping all sequences.\nstderr: {result.stderr}",
                file=sys.stderr,
            )
            return records

        # Parse .clstr to find representative IDs
        representatives: set[str] = set()
        clstr_file = output_fa + ".clstr"
        if os.path.exists(clstr_file):
            with open(clstr_file) as fh:
                for line in fh:
                    if line.strip().endswith("*"):
                        start = line.index(">") + 1
                        end = line.index("...", start)
                        representatives.add(line[start:end])
        else:
            # Fallback: parse output FASTA directly
            for rec in SeqIO.parse(output_fa, "fasta"):
                representatives.add(rec.id)

    return [pair for rec_id, pair in id_to_pair.items() if rec_id in representatives]


# ---------------------------------------------------------------------------
# Pool building
# ---------------------------------------------------------------------------

def build_pool_for_class(
    records: list[tuple[int, SeqRecord]],
    must_include_taxids: frozenset[int],
    berghia_taxid: int,
    max_size: int,
    lineage_fn,
    berghia_records: Optional[list[SeqRecord]] = None,
) -> tuple[list[tuple[int, SeqRecord]], dict]:
    """Build a single per-class reference pool.

    Steps:
    1. Prepend Berghia class-specific records as MUST_INCLUDE (always kept).
    2. Separate must-include taxid records from ordinary candidates.
    3. Fill Berghia + must-include sequences first (up to max_size).
    4. Fill remaining slots with Berghia-proximity-weighted sampling
       from the ordinary candidates.

    Args:
        records: (taxid, SeqRecord) pairs for this class (already CD-HIT'd).
        must_include_taxids: taxids whose sequences are always kept.
        berghia_taxid: focal taxid — used for proximity scoring (no longer excluded).
        max_size: per-class sequence cap (refs only; Berghia seqs are on top).
        lineage_fn: callable(taxid) → list[int].
        berghia_records: SeqRecord list for Berghia sequences of this class;
            prepended before any cap logic so they are never subsampled away.

    Returns:
        (selected_pairs, stats_dict)
    """
    # Berghia seqs are stored as (berghia_taxid, rec) to unify the pair format
    berghia_pairs: list[tuple[int, SeqRecord]] = [
        (berghia_taxid, rec) for rec in (berghia_records or [])
    ]

    must_records: list[tuple[int, SeqRecord]] = []
    ordinary_records: list[tuple[int, SeqRecord]] = []

    must_taxids_with_hits: set[int] = set()

    for taxid, rec in records:
        if taxid in must_include_taxids:
            must_records.append((taxid, rec))
            must_taxids_with_hits.add(taxid)
        else:
            ordinary_records.append((taxid, rec))

    missing_must = sorted(must_include_taxids - must_taxids_with_hits)

    # Berghia seqs (class-specific) are always retained first, then must-include
    n_berghia = len(berghia_pairs)
    selected_berghia = berghia_pairs  # all Berghia seqs for this class are kept

    # Must-include refs fill up to max_size after Berghia occupies its slots
    selected_must = must_records[:max(0, max_size - n_berghia)]
    remaining_slots = max_size - n_berghia - len(selected_must)

    selected_ordinary: list[tuple[int, SeqRecord]] = []
    if remaining_slots > 0 and ordinary_records:
        if len(ordinary_records) <= remaining_slots:
            selected_ordinary = ordinary_records
        else:
            # Berghia-proximity weighted sampling
            weights = [
                proximity_score(taxid, lineage_fn, berghia_taxid)
                for taxid, _ in ordinary_records
            ]
            # Normalise; default to 1 if all zero
            total_w = sum(weights)
            if total_w == 0:
                weights = [1.0] * len(ordinary_records)
                total_w = float(len(ordinary_records))
            norm_weights = [w / total_w for w in weights]

            indices = list(range(len(ordinary_records)))
            chosen_indices = random.choices(
                indices,
                weights=norm_weights,
                k=remaining_slots * 4,  # over-sample then deduplicate
            )
            seen: set[int] = set()
            selected_indices: list[int] = []
            for idx in chosen_indices:
                if idx not in seen:
                    seen.add(idx)
                    selected_indices.append(idx)
                if len(selected_indices) == remaining_slots:
                    break
            # If weighted sampling came up short, pad sequentially
            if len(selected_indices) < remaining_slots:
                for idx in range(len(ordinary_records)):
                    if idx not in seen:
                        selected_indices.append(idx)
                    if len(selected_indices) == remaining_slots:
                        break
            selected_ordinary = [ordinary_records[i] for i in selected_indices]

    selected = selected_berghia + selected_must + selected_ordinary
    species_contributing = len({taxid for taxid, _ in selected})

    stats = {
        "n_total_candidates": len(records),
        "n_after_cdhit": len(records),   # caller may update this field
        "n_berghia_included": n_berghia,
        "n_must_include": len(selected_must),
        "n_subsampled": len(selected_ordinary),
        "n_output": len(selected),
        "species_contributing": species_contributing,
        "must_include_taxids_with_hits": sorted(must_taxids_with_hits),
        "must_include_taxids_missing": missing_must,
    }
    return selected, stats


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def build_all_pools(
    scan_fasta_glob: str,
    class_tsv: str,
    out_dir: str,
    max_per_class: int = 2000,
    max_class_caps: Optional[dict[str, int]] = None,
    total_budget_per_class: int = DEFAULT_TOTAL_BUDGET_PER_CLASS,
    cluster_identity: float = 0.7,
    cdhit_path: str = "cd-hit",
    threads: int = 4,
    must_include_taxids: frozenset[int] = DEFAULT_MUST_INCLUDE_TAXIDS,
    berghia_taxid: int = BERGHIA_TAXID_DEFAULT,
    berghia_fasta: Optional[str] = None,
    berghia_class_tsv: Optional[str] = None,
    force: bool = False,
) -> None:
    """Build four per-class reference pools.

    Args:
        scan_fasta_glob: glob pattern matching per-species scan FASTA files.
        class_tsv: path to P1 classifier output TSV.
        out_dir: directory for output files.
        max_per_class: (deprecated) fallback per-class sequence cap. Use total_budget_per_class.
        max_class_caps: (deprecated) per-class overrides, e.g. {"A": 2000, "B": 2900}.
            If provided, these override total_budget_per_class and emit a deprecation warning.
        total_budget_per_class: per-class total budget (refs + Berghia count). Default 3000.
            Per-class ref cap is computed as: total_budget_per_class - count_of_berghia_in_class.
        cluster_identity: CD-HIT identity threshold (default 0.7).
        cdhit_path: path to cd-hit binary.
        threads: threads for CD-HIT.
        must_include_taxids: taxids whose candidates are always kept.
        berghia_taxid: focal species taxid (used for proximity scoring).
        berghia_fasta: optional path to Berghia candidate FASTA.
        berghia_class_tsv: optional path to Berghia P1 class TSV.
            Required when berghia_fasta is provided.
        force: re-run even if outputs already exist.
    """
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    output_fastas = [out_path / f"refs_class_{cls}.fa" for cls in GPCR_CLASSES]
    report_path = out_path / "pool_build_report.json"

    # Idempotency check
    if not force and all(f.exists() for f in output_fastas) and report_path.exists():
        print(
            f"[build_per_class_reference_pools] Outputs exist in {out_dir}; "
            "skipping (use --force to rerun).",
            file=sys.stderr,
        )
        return

    # Resolve per-class caps: explicit overrides > max_per_class fallback
    _caps: dict[str, int] = max_class_caps or {}
    class_caps: dict[str, int] = {
        cls: _caps.get(cls, max_per_class) for cls in GPCR_CLASSES
    }

    # --- 1. Load class map -----------------------------------------------
    class_map = load_class_tsv(class_tsv)

    # --- 2. Load Berghia candidates (optional) ----------------------------
    berghia_per_class: dict[str, list[SeqRecord]] = {cls: [] for cls in GPCR_CLASSES}
    n_berghia_per_class: dict[str, int] = {cls: 0 for cls in GPCR_CLASSES}
    if berghia_fasta and berghia_class_tsv:
        berghia_per_class = load_berghia_candidates(
            berghia_fasta, berghia_class_tsv, berghia_taxid
        )
        # Count Berghia seqs per class for dynamic cap computation
        for cls in GPCR_CLASSES:
            n_berghia_per_class[cls] = len(berghia_per_class[cls])
    elif berghia_fasta:
        print(
            "  WARNING: --berghia-fasta provided but --berghia-class-tsv missing; "
            "Berghia candidates will not be added to pools.",
            file=sys.stderr,
        )

    # --- 3. Read scan FASTAs → per-class (taxid, SeqRecord) lists ---------
    per_class: dict[str, list[tuple[int, SeqRecord]]] = {
        cls: [] for cls in GPCR_CLASSES
    }
    unclassified_records: list[tuple[int, str, str]] = []  # (taxid, seq_id, reason)
    missing_from_class_map: set[str] = set()

    scan_files = sorted(glob.glob(scan_fasta_glob))
    print(
        f"[build_per_class_reference_pools] Found {len(scan_files)} scan FASTA(s)",
        file=sys.stderr,
    )

    for fa_path in scan_files:
        taxid = taxid_from_filename(os.path.basename(fa_path))
        if taxid is None:
            print(
                f"  WARNING: cannot parse taxid from {fa_path}; skipping",
                file=sys.stderr,
            )
            continue

        for rec in SeqIO.parse(fa_path, "fasta"):
            seq_class = class_map.get(rec.id)
            if seq_class is None:
                missing_from_class_map.add(rec.id)
                unclassified_records.append((taxid, rec.id, "not_in_class_tsv"))
            elif seq_class == "unclassified":
                unclassified_records.append((taxid, rec.id, "unclassified"))
            elif seq_class in per_class:
                per_class[seq_class].append((taxid, rec))
            else:
                unclassified_records.append((taxid, rec.id, f"unknown_class:{seq_class}"))

    if missing_from_class_map:
        print(
            f"  WARNING: {len(missing_from_class_map)} seq IDs not found in class TSV "
            "(logged as unclassified)",
            file=sys.stderr,
        )

    # --- 4. Compute dynamic per-class caps + CD-HIT dedup + build pool -----
    # Emit deprecation warning if max_class_caps was explicitly provided
    if max_class_caps:
        print(
            "  WARNING: --max-class-A/B/C/F are deprecated; use --total-budget-per-class instead. "
            "Overriding computed caps with explicit values.",
            file=sys.stderr,
        )
        # Fill in missing classes from max_class_caps with defaults
        class_caps_computed: dict[str, int] = {}
        for cls in GPCR_CLASSES:
            if cls in max_class_caps:
                class_caps_computed[cls] = max_class_caps[cls]
            else:
                # Fill missing classes with the default for that class
                class_caps_computed[cls] = {
                    "A": DEFAULT_MAX_CLASS_A,
                    "B": DEFAULT_MAX_CLASS_B,
                    "C": DEFAULT_MAX_CLASS_C,
                    "F": DEFAULT_MAX_CLASS_F,
                }[cls]
    else:
        # Dynamic computation: max_refs(class) = total_budget - berghia_count(class)
        class_caps_computed = {
            cls: max(0, total_budget_per_class - n_berghia_per_class[cls])
            for cls in GPCR_CLASSES
        }
        print(
            f"  [P2] Computing per-class ref caps from total_budget={total_budget_per_class}, "
            f"berghia_counts={n_berghia_per_class}: {class_caps_computed}",
            file=sys.stderr,
        )

    report: dict = {
        "total_budget_per_class": total_budget_per_class,
        "n_berghia_per_class": n_berghia_per_class,
    }
    lineage_fn = get_lineage  # can be replaced by tests via patch

    for cls in GPCR_CLASSES:
        candidates = per_class[cls]
        n_total = len(candidates)
        cap = class_caps_computed[cls]
        print(
            f"  Class {cls}: {n_total} candidates before CD-HIT (ref_cap={cap}, "
            f"total_budget={total_budget_per_class}, berghia={n_berghia_per_class[cls]})",
            file=sys.stderr,
        )

        after_cdhit = cdhit_dedup(
            candidates,
            identity=cluster_identity,
            threads=threads,
            cdhit_path=cdhit_path,
        )
        n_after = len(after_cdhit)
        print(
            f"  Class {cls}: {n_after} representatives after CD-HIT",
            file=sys.stderr,
        )

        selected, stats = build_pool_for_class(
            records=after_cdhit,
            must_include_taxids=must_include_taxids,
            berghia_taxid=berghia_taxid,
            max_size=cap,
            lineage_fn=lineage_fn,
            berghia_records=berghia_per_class.get(cls),
        )
        stats["n_total_candidates"] = n_total
        stats["n_after_cdhit"] = n_after
        # Add budget info to the per-class stats
        stats["n_refs_target"] = cap
        stats["total_budget_for_class"] = total_budget_per_class

        # Write FASTA
        out_fa = out_path / f"refs_class_{cls}.fa"
        with open(out_fa, "w") as fh:
            for _, rec in selected:
                fh.write(f">{rec.id}\n{str(rec.seq)}\n")

        report[f"class_{cls}"] = stats
        print(
            f"  Class {cls}: wrote {stats['n_output']} sequences to {out_fa} "
            f"({stats['n_berghia_included']} Berghia, {cap} refs targeted)",
            file=sys.stderr,
        )

    # --- 5. Unclassified log ---------------------------------------------
    unclass_path = out_path / "unclassified_log.tsv"
    with open(unclass_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["taxid", "seq_id", "reason"])
        for row in unclassified_records:
            writer.writerow(row)

    report["unclassified"] = {
        "n": len(unclassified_records),
        "log_path": str(unclass_path),
    }
    berghia_totals = {
        cls: len(berghia_per_class[cls]) for cls in GPCR_CLASSES
    }
    report["berghia_included"] = {
        "taxid": berghia_taxid,
        "n_per_class": berghia_totals,
        "n_total": sum(berghia_totals.values()),
    }

    # --- 6. JSON report ---------------------------------------------------
    with open(report_path, "w") as fh:
        json.dump(report, fh, indent=2)
    print(
        f"[build_per_class_reference_pools] Report written to {report_path}",
        file=sys.stderr,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_args_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build per-class GPCR reference pools for phylogenetic tree inference. "
            "Consumes P0 scan FASTAs + P1 class TSV. "
            "Outputs refs_class_{A,B,C,F}.fa (≤MAX_PHYLO_REFS each) + report."
        )
    )
    parser.add_argument(
        "--scan-fasta-glob",
        required=True,
        help=(
            "Glob pattern matching per-species chemo_candidates FASTA files, "
            "e.g. 'scan_output/*.chemo_candidates.fa'"
        ),
    )
    parser.add_argument(
        "--class-tsv",
        required=True,
        help="P1 classifier output TSV (columns: seq_id, class, evidence_pfam, top_evalue)",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        help="Output directory for per-class FASTAs, unclassified log, and JSON report",
    )
    parser.add_argument(
        "--max-per-class",
        type=int,
        default=2000,
        metavar="N",
        help="(deprecated) Maximum sequences per class FASTA. Use --total-budget-per-class.",
    )
    parser.add_argument(
        "--total-budget-per-class",
        type=int,
        default=int(os.environ.get("TOTAL_BUDGET_PER_CLASS", str(DEFAULT_TOTAL_BUDGET_PER_CLASS))),
        metavar="N",
        help=(
            "Total budget per class (refs + Berghia count). "
            "Per-class ref cap = total_budget - berghia_count(class). "
            f"(default: {DEFAULT_TOTAL_BUDGET_PER_CLASS}; env: TOTAL_BUDGET_PER_CLASS)"
        ),
    )
    parser.add_argument(
        "--cluster-identity",
        type=float,
        default=0.7,
        metavar="FLOAT",
        help="CD-HIT identity threshold for deduplication (default: 0.7)",
    )
    parser.add_argument(
        "--cdhit-path",
        default="cd-hit",
        metavar="PATH",
        help="Path to cd-hit binary (default: cd-hit)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=int(os.environ.get("CPUS", "4")),
        metavar="N",
        help="Threads for CD-HIT (default: $CPUS or 4)",
    )
    parser.add_argument(
        "--must-include-taxids",
        default=None,
        metavar="TSV",
        help=(
            "TSV file with a 'taxid' column; overrides built-in MUST_INCLUDE list. "
            "Sequences from these taxa are always retained in their respective pools."
        ),
    )
    parser.add_argument(
        "--berghia-taxid",
        type=int,
        default=BERGHIA_TAXID_DEFAULT,
        metavar="TAXID",
        help=f"Taxid of the focal species (used for proximity scoring; default: {BERGHIA_TAXID_DEFAULT})",
    )
    parser.add_argument(
        "--berghia-fasta",
        default=None,
        metavar="FASTA",
        help=(
            "FASTA of Berghia chemoreceptor candidates from stage-02. "
            "When provided together with --berghia-class-tsv, Berghia sequences "
            "of each class are added as MUST_INCLUDE to that class's pool."
        ),
    )
    parser.add_argument(
        "--berghia-class-tsv",
        default=None,
        metavar="TSV",
        help="P1 class TSV for Berghia candidates (required when --berghia-fasta is used).",
    )
    parser.add_argument(
        "--max-class-A",
        type=int,
        default=None,
        metavar="N",
        help="(deprecated) Max refs for Class A pool. Use --total-budget-per-class.",
    )
    parser.add_argument(
        "--max-class-B",
        type=int,
        default=None,
        metavar="N",
        help="(deprecated) Max refs for Class B pool. Use --total-budget-per-class.",
    )
    parser.add_argument(
        "--max-class-C",
        type=int,
        default=None,
        metavar="N",
        help="(deprecated) Max refs for Class C pool. Use --total-budget-per-class.",
    )
    parser.add_argument(
        "--max-class-F",
        type=int,
        default=None,
        metavar="N",
        help="(deprecated) Max refs for Class F pool. Use --total-budget-per-class.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Re-run even if all output files already exist",
    )
    return parser


def main(argv=None) -> None:
    parser = build_args_parser()
    args = parser.parse_args(argv)

    must_include: frozenset[int]
    if args.must_include_taxids:
        must_include = load_must_include_taxids(args.must_include_taxids)
        print(
            f"[build_per_class_reference_pools] Loaded {len(must_include)} "
            f"must-include taxids from {args.must_include_taxids}",
            file=sys.stderr,
        )
    else:
        must_include = DEFAULT_MUST_INCLUDE_TAXIDS

    # Build per-class cap overrides from CLI (deprecated; emit warning if used)
    max_class_caps: Optional[dict[str, int]] = None
    if any(x is not None for x in [args.max_class_A, args.max_class_B, args.max_class_C, args.max_class_F]):
        max_class_caps = {
            "A": args.max_class_A or DEFAULT_MAX_CLASS_A,
            "B": args.max_class_B or DEFAULT_MAX_CLASS_B,
            "C": args.max_class_C or DEFAULT_MAX_CLASS_C,
            "F": args.max_class_F or DEFAULT_MAX_CLASS_F,
        }

    build_all_pools(
        scan_fasta_glob=args.scan_fasta_glob,
        class_tsv=args.class_tsv,
        out_dir=args.out_dir,
        max_per_class=args.max_per_class,
        max_class_caps=max_class_caps,
        total_budget_per_class=args.total_budget_per_class,
        cluster_identity=args.cluster_identity,
        cdhit_path=args.cdhit_path,
        threads=args.threads,
        must_include_taxids=must_include,
        berghia_taxid=args.berghia_taxid,
        berghia_fasta=args.berghia_fasta,
        berghia_class_tsv=args.berghia_class_tsv,
        force=args.force,
    )


if __name__ == "__main__":
    main()
