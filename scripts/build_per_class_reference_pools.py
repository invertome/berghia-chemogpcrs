#!/usr/bin/env python3
"""
build_per_class_reference_pools.py

P2 of the per-class refactor.

Consumes:
  - Per-species chemo_candidates.fa files from a P0 scan (glob pattern)
  - A P1 class TSV mapping seq_id → {A,B,C,F,unclassified}
  - (optional) Berghia candidate FASTA + Berghia P1 class TSV

Produces, under --out-dir:
  - refs_class_A.fa
  - refs_class_B.fa
  - refs_class_C.fa
  - refs_class_F.fa
  - unclassified_log.tsv
  - pool_build_report.json

Each class pool sized to total_budget_per_class - n_berghia(class) - outgroup_budget_per_class refs
(default 3000 - n_berghia - 10), so total alignment per class ≤ 3000 (mafft-dash cap).

Design constraints (locked, 2026-05-28):
  - Berghia candidates are classified by P1 into per-class labels.
  - Each per-class pool includes Berghia sequences of that class as
    MUST_INCLUDE (always retained, not subsampled away).
  - MUST_INCLUDE taxids are always retained (their candidates fill first).
  - Remaining slots are filled by Berghia-proximity-weighted sampling
    of the non-must-include candidates.
  - CD-HIT 0.9 deduplication before subsampling (0.9, not 0.7: preserves the
    cross-lineage-informative older paralog structure that polarizes shared vs
    Berghia-specific expansions; only lineage-private near-identical recent
    duplicates collapse. Berghia/anchors are never deduped; balanced sampling,
    not CD-HIT, enforces the budget).
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

# Default outgroup budget per class (reserves slots for outgroup sequences)
DEFAULT_OUTGROUP_BUDGET_PER_CLASS = 10

# Balanced reference sampling (option 4b). Every species with candidates gets at
# least FLOOR reps; no taxon exceeds CAP reps; the rest is Berghia-proximity-
# weighted. Replaces the old uncapped-must-include flooding. Tunable; calibrated
# on the Phase-1a pilot. SEED makes selection reproducible.
DEFAULT_REF_FLOOR_PER_SPECIES = 5
DEFAULT_REF_CAP_PER_TAXON = 150
DEFAULT_SELECTION_SEED = 12345

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
        berghia_taxid: the focal taxid (used for proximity scoring).

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


def load_anchor_set(
    anchor_fasta: str,
    anchor_tsv: str,
) -> dict[str, list[tuple[int, SeqRecord]]]:
    """Load the curated anchor set into per-class (taxid, SeqRecord) bins.

    Anchors are characterized GPCR landmarks (build_anchor_set.py). FASTA
    headers are ``ANCHOR_<class>_<tier>_<accession>``; the class is read from
    the header and the taxid from the TSV (keyed by accession). Each anchor's
    SeqRecord.id keeps the full ANCHOR_ header so anchors are traceable in the
    final tree.

    Args:
        anchor_fasta: path to anchor_set.fasta.
        anchor_tsv: path to anchor_set.tsv (cols incl. accession, taxid).

    Returns:
        dict mapping class label (A/B/C/F) → list of (taxid, SeqRecord).
    """
    acc_taxid: dict[str, int] = {}
    with open(anchor_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                acc_taxid[row["accession"]] = int(row["taxid"])
            except (KeyError, ValueError):
                pass

    per_class: dict[str, list[tuple[int, SeqRecord]]] = {
        cls: [] for cls in GPCR_CLASSES
    }
    for rec in SeqIO.parse(anchor_fasta, "fasta"):
        parts = rec.id.split("_")
        if len(parts) < 4 or parts[0] != "ANCHOR":
            continue
        klass = parts[1]
        accession = "_".join(parts[3:])
        taxid = acc_taxid.get(accession, 0)
        if klass in per_class:
            per_class[klass].append((taxid, rec))

    total = sum(len(v) for v in per_class.values())
    print(
        f"[build_per_class_reference_pools] Loaded {total} anchors from "
        f"{anchor_fasta} (per-class: "
        f"{ {c: len(per_class[c]) for c in GPCR_CLASSES} })",
        file=sys.stderr,
    )
    return per_class


# ---------------------------------------------------------------------------
# CD-HIT deduplication
# ---------------------------------------------------------------------------

def _parse_clstr(clstr_path: str) -> list[dict]:
    """Parse a CD-HIT .clstr file into a list of clusters.

    Returns one dict per cluster: {"rep": <representative id or None>,
    "members": [<member ids>]}.
    """
    clusters: list[dict] = []
    cur: Optional[dict] = None
    with open(clstr_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    clusters.append(cur)
                cur = {"rep": None, "members": []}
            elif ">" in line and "..." in line:
                start = line.index(">") + 1
                end = line.index("...", start)
                mid = line[start:end]
                if cur is None:
                    cur = {"rep": None, "members": []}
                cur["members"].append(mid)
                if line.rstrip().endswith("*"):
                    cur["rep"] = mid
    if cur is not None:
        clusters.append(cur)
    return clusters


def _cdhit_survivor_ids(clusters: list[dict],
                        must_keep_ids: frozenset[str]) -> set[str]:
    """Choose surviving ids per cluster.

    A cluster that contains one or more must-keep ids (anchors) keeps ALL of
    them and drops every other member — so a scanned candidate never displaces
    an annotated anchor (keep-annotated-on-collision). A cluster with no
    must-keep id keeps its CD-HIT representative, as usual.
    """
    survivors: set[str] = set()
    for c in clusters:
        keep = [m for m in c["members"] if m in must_keep_ids]
        if keep:
            survivors.update(keep)
        else:
            rep = c["rep"]
            if rep is None and c["members"]:
                rep = c["members"][0]
            if rep is not None:
                survivors.add(rep)
    return survivors


def cdhit_dedup(
    records: list[tuple[int, SeqRecord]],
    identity: float = 0.9,
    threads: int = 4,
    cdhit_path: str = "cd-hit",
    must_keep_ids: frozenset[str] = frozenset(),
) -> list[tuple[int, SeqRecord]]:
    """Run CD-HIT on *records* and return representatives.

    Args:
        records: list of (taxid, SeqRecord) pairs.
        identity: clustering threshold.
        threads: CPU threads for CD-HIT.
        cdhit_path: path to cd-hit binary.
        must_keep_ids: record ids (anchors) that win any cluster they fall in —
            kept as the survivor while colliding candidates are dropped.

    Returns:
        Filtered list of (taxid, SeqRecord) — one representative per cluster,
        except clusters touching a must-keep id, which keep the anchor(s).
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

        # Parse .clstr → clusters → survivors (anchors win their clusters)
        clstr_file = output_fa + ".clstr"
        if os.path.exists(clstr_file):
            clusters = _parse_clstr(clstr_file)
            representatives = _cdhit_survivor_ids(clusters, must_keep_ids)
        else:
            # Fallback: parse output FASTA directly; anchors always kept.
            representatives = {rec.id for rec in SeqIO.parse(output_fa, "fasta")}
            representatives |= {rid for rid in id_to_pair if rid in must_keep_ids}

    return [pair for rec_id, pair in id_to_pair.items() if rec_id in representatives]


# ---------------------------------------------------------------------------
# Pool building
# ---------------------------------------------------------------------------

def balanced_reference_sample(
    records: list[tuple[int, SeqRecord]],
    budget: int,
    *,
    floor: int,
    cap: int,
    must_include_taxids: frozenset[int],
    lineage_fn,
    berghia_taxid: int,
    seed: int,
) -> list[tuple[int, SeqRecord]]:
    """Select <= ``budget`` reference records, balanced across taxa.

    Spends a fixed budget so every species is represented and no taxon floods
    the pool (the uncapped-must-include flooding biased LSE inference):
      - per-species FLOOR: each taxon gets min(floor, available);
      - must-include taxa filled to min(cap, available) (kept, but bounded);
      - per-taxon CAP: no taxon exceeds min(cap, available);
      - remaining budget filled by Berghia-proximity-weighted draw, capped.

    Deterministic given ``seed`` — replaces the old unseeded proximity sampling,
    so reference pools are reproducible (and the anchor with-vs-without eval
    compares against a stable base pool).
    """
    rng = random.Random(seed)
    by_taxid: dict[int, list[SeqRecord]] = {}
    for taxid, rec in records:
        by_taxid.setdefault(taxid, []).append(rec)
    if not by_taxid or budget <= 0:
        return []

    ceil = {t: min(cap, len(recs)) for t, recs in by_taxid.items()}
    alloc = {t: min(floor, ceil[t]) for t in by_taxid}            # per-species floor
    for t in by_taxid:                                            # must-include -> cap
        if t in must_include_taxids:
            alloc[t] = ceil[t]

    total = sum(alloc.values())
    if total > budget:
        # Floor + must-include already exceed the budget (more species than the
        # budget can floor). Trim the largest allocations down until within budget.
        order = sorted(by_taxid, key=lambda t: (alloc[t], t))
        i = len(order) - 1
        while total > budget:
            t = order[i]
            if alloc[t] > 0:
                alloc[t] -= 1
                total -= 1
            i = i - 1 if i > 0 else len(order) - 1
    else:
        # Berghia-proximity-weighted fill of the remaining budget (proximity is
        # constant per taxon → memoise). Each pick respects the per-taxon cap.
        prox = {t: max(proximity_score(t, lineage_fn, berghia_taxid), 0)
                for t in by_taxid}
        remaining = budget - total
        growable = [t for t in by_taxid if alloc[t] < ceil[t]]
        while remaining > 0 and growable:
            weights = [prox[t] for t in growable]
            if sum(weights) == 0:
                weights = [1] * len(growable)
            t = rng.choices(growable, weights=weights, k=1)[0]
            alloc[t] += 1
            remaining -= 1
            if alloc[t] >= ceil[t]:
                growable.remove(t)

    selected: list[tuple[int, SeqRecord]] = []
    for t in sorted(by_taxid):                                    # sorted → deterministic
        recs = by_taxid[t]
        n = alloc[t]
        if n >= len(recs):
            idx = list(range(len(recs)))
        else:
            idx = sorted(rng.sample(range(len(recs)), n))
        selected.extend((t, recs[i]) for i in idx)
    return selected


def build_pool_for_class(
    records: list[tuple[int, SeqRecord]],
    must_include_taxids: frozenset[int],
    berghia_taxid: int,
    max_size: int,
    lineage_fn,
    berghia_records: Optional[list[SeqRecord]] = None,
    *,
    floor: int = DEFAULT_REF_FLOOR_PER_SPECIES,
    cap: int = DEFAULT_REF_CAP_PER_TAXON,
    seed: int = DEFAULT_SELECTION_SEED,
    anchor_records: Optional[list[tuple[int, SeqRecord]]] = None,
) -> tuple[list[tuple[int, SeqRecord]], dict]:
    """Build a single per-class reference pool (balanced — option 4b).

    1. All Berghia class-specific records are kept (guaranteed, on top of the
       refs budget ``max_size``).
    2. Anchors are kept (must-keep, functional landmarks) and charged to the
       refs budget — so they come out of ``max_size`` rather than on top, the
       balanced reference draw being shrunk by the anchor count.
    3. The remaining references are chosen by :func:`balanced_reference_sample`:
       every species with candidates gets a FLOOR, must-include model species
       are filled to CAP, no taxon exceeds CAP, and the remainder is Berghia-
       proximity-weighted — so all species are represented and no taxon floods
       the pool. Deterministic via ``seed``.

    Args:
        records: (taxid, SeqRecord) pairs for this class (already CD-HIT'd).
        must_include_taxids: model-species taxa filled to CAP (kept, but bounded).
        berghia_taxid: focal taxid — used for proximity scoring.
        max_size: refs-only budget; Berghia seqs are guaranteed on top, anchors
            are charged against it.
        lineage_fn: callable(taxid) → list[int].
        berghia_records: Berghia SeqRecords for this class (always kept).
        floor: minimum reps per species (default DEFAULT_REF_FLOOR_PER_SPECIES).
        cap: maximum reps per taxon (default DEFAULT_REF_CAP_PER_TAXON).
        seed: selection seed for reproducibility.
        anchor_records: (taxid, SeqRecord) anchor pairs for this class — kept
            (must-keep) and charged to the refs budget.

    Returns:
        (selected_pairs, stats_dict)
    """
    berghia_pairs: list[tuple[int, SeqRecord]] = [
        (berghia_taxid, rec) for rec in (berghia_records or [])
    ]
    n_berghia = len(berghia_pairs)

    anchor_pairs: list[tuple[int, SeqRecord]] = list(anchor_records or [])
    n_anchors = len(anchor_pairs)

    present_taxids = {taxid for taxid, _ in records}
    must_taxids_with_hits = sorted(must_include_taxids & present_taxids)
    missing_must = sorted(must_include_taxids - present_taxids)

    # Anchors are charged to the refs budget: shrink the balanced-draw budget by
    # the anchor count (anchors themselves are always kept on top of the draw).
    refs_budget = max(0, max_size - n_anchors)
    selected_refs = balanced_reference_sample(
        records,
        budget=refs_budget,
        floor=floor,
        cap=cap,
        must_include_taxids=must_include_taxids,
        lineage_fn=lineage_fn,
        berghia_taxid=berghia_taxid,
        seed=seed,
    )
    selected = berghia_pairs + anchor_pairs + selected_refs

    # n_must_include / n_subsampled retained for report back-compat: refs drawn
    # from must-include taxa vs the rest.
    n_must = sum(1 for taxid, _ in selected_refs if taxid in must_include_taxids)
    n_other = len(selected_refs) - n_must
    species_contributing = len({taxid for taxid, _ in selected})

    stats = {
        "n_total_candidates": len(records),
        "n_after_cdhit": len(records),   # caller may update this field
        "n_berghia_included": n_berghia,
        "n_anchors": n_anchors,
        "n_must_include": n_must,
        "n_subsampled": n_other,
        "n_output": len(selected),
        "species_contributing": species_contributing,
        "must_include_taxids_with_hits": must_taxids_with_hits,
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
    total_budget_per_class: int = DEFAULT_TOTAL_BUDGET_PER_CLASS,
    outgroup_budget_per_class: int = DEFAULT_OUTGROUP_BUDGET_PER_CLASS,
    cluster_identity: float = 0.9,
    cdhit_path: str = "cd-hit",
    threads: int = 4,
    must_include_taxids: frozenset[int] = DEFAULT_MUST_INCLUDE_TAXIDS,
    berghia_taxid: int = BERGHIA_TAXID_DEFAULT,
    berghia_fasta: Optional[str] = None,
    berghia_class_tsv: Optional[str] = None,
    anchor_fasta: Optional[str] = None,
    anchor_tsv: Optional[str] = None,
    ref_floor_per_species: int = DEFAULT_REF_FLOOR_PER_SPECIES,
    ref_cap_per_taxon: int = DEFAULT_REF_CAP_PER_TAXON,
    selection_seed: int = DEFAULT_SELECTION_SEED,
    force: bool = False,
) -> None:
    """Build four per-class reference pools.

    Args:
        scan_fasta_glob: glob pattern matching per-species scan FASTA files.
        class_tsv: path to P1 classifier output TSV.
        out_dir: directory for output files.
        total_budget_per_class: per-class total budget (refs + Berghia + outgroup count). Default 3000.
            Per-class ref cap is computed as: total_budget_per_class - count_of_berghia_in_class - outgroup_budget_per_class.
        outgroup_budget_per_class: reserved slots for outgroup sequences per class. Default 10.
        cluster_identity: CD-HIT identity threshold (default 0.9).
        cdhit_path: path to cd-hit binary.
        threads: threads for CD-HIT.
        must_include_taxids: taxids whose candidates are always kept.
        berghia_taxid: focal species taxid (used for proximity scoring).
        berghia_fasta: optional path to Berghia candidate FASTA.
        berghia_class_tsv: optional path to Berghia P1 class TSV.
            Required when berghia_fasta is provided.
        anchor_fasta: optional path to anchor_set.fasta (characterized GPCR
            landmarks). When provided with anchor_tsv, anchors of each class are
            injected must-keep, through CD-HIT (keep-annotated-on-collision),
            charged to that class's refs budget.
        anchor_tsv: optional path to anchor_set.tsv (required with anchor_fasta).
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

    # --- 2b. Load anchor set (optional) -----------------------------------
    anchor_per_class: dict[str, list[tuple[int, SeqRecord]]] = {
        cls: [] for cls in GPCR_CLASSES
    }
    if anchor_fasta and anchor_tsv:
        if os.path.exists(anchor_fasta) and os.path.exists(anchor_tsv):
            anchor_per_class = load_anchor_set(anchor_fasta, anchor_tsv)
        else:
            print(
                f"  WARNING: anchor set not found ({anchor_fasta} / {anchor_tsv}); "
                "anchors will not be injected. Run build_anchor_set.py first.",
                file=sys.stderr,
            )
    elif anchor_fasta:
        print(
            "  WARNING: --anchor-fasta provided but --anchor-tsv missing; "
            "anchors will not be injected.",
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
    # Dynamic computation: max_refs(class) = total_budget - berghia_count(class) - outgroup_budget
    class_caps_computed = {
        cls: max(0, total_budget_per_class - n_berghia_per_class[cls] - outgroup_budget_per_class)
        for cls in GPCR_CLASSES
    }
    print(
        f"  [P2] Computing per-class ref caps from total_budget={total_budget_per_class}, "
        f"berghia_counts={n_berghia_per_class}, outgroup_budget={outgroup_budget_per_class}: "
        f"{class_caps_computed}",
        file=sys.stderr,
    )

    report: dict = {
        "total_budget_per_class": total_budget_per_class,
        "outgroup_budget_per_class": outgroup_budget_per_class,
        "n_berghia_per_class": n_berghia_per_class,
    }
    lineage_fn = get_lineage  # can be replaced by tests via patch

    for cls in GPCR_CLASSES:
        candidates = per_class[cls]
        n_total = len(candidates)
        cap = class_caps_computed[cls]
        anchor_pairs = anchor_per_class.get(cls, [])
        anchor_ids = frozenset(rec.id for _, rec in anchor_pairs)
        print(
            f"  Class {cls}: {n_total} candidates before CD-HIT (ref_cap={cap}, "
            f"total_budget={total_budget_per_class}, berghia={n_berghia_per_class[cls]}, "
            f"anchors={len(anchor_pairs)})",
            file=sys.stderr,
        )

        # Anchors go THROUGH CD-HIT with the candidates and win any collision
        # (must-keep): a scanned candidate redundant with an annotated anchor is
        # dropped, the anchor is kept. Berghia is never in this set, so it is
        # never absorbed.
        after_cdhit = cdhit_dedup(
            candidates + anchor_pairs,
            identity=cluster_identity,
            threads=threads,
            cdhit_path=cdhit_path,
            must_keep_ids=anchor_ids,
        )
        anchor_survivors = [(t, r) for t, r in after_cdhit if r.id in anchor_ids]
        candidate_survivors = [(t, r) for t, r in after_cdhit if r.id not in anchor_ids]
        n_after = len(candidate_survivors)
        print(
            f"  Class {cls}: {n_after} candidate representatives after CD-HIT "
            f"(+{len(anchor_survivors)} anchors kept)",
            file=sys.stderr,
        )

        selected, stats = build_pool_for_class(
            records=candidate_survivors,
            must_include_taxids=must_include_taxids,
            berghia_taxid=berghia_taxid,
            max_size=cap,
            lineage_fn=lineage_fn,
            berghia_records=berghia_per_class.get(cls),
            floor=ref_floor_per_species,
            cap=ref_cap_per_taxon,
            seed=selection_seed,
            anchor_records=anchor_survivors,
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
        "--total-budget-per-class",
        type=int,
        default=int(os.environ.get("TOTAL_BUDGET_PER_CLASS", str(DEFAULT_TOTAL_BUDGET_PER_CLASS))),
        metavar="N",
        help=(
            "Total budget per class (refs + Berghia + outgroup count). "
            "Per-class ref cap = total_budget - berghia_count(class) - outgroup_budget. "
            f"(default: {DEFAULT_TOTAL_BUDGET_PER_CLASS}; env: TOTAL_BUDGET_PER_CLASS)"
        ),
    )
    parser.add_argument(
        "--outgroup-budget-per-class",
        type=int,
        default=int(os.environ.get("OUTGROUP_BUDGET_PER_CLASS", str(DEFAULT_OUTGROUP_BUDGET_PER_CLASS))),
        metavar="N",
        help=(
            "Reserved slots for outgroup sequences per class. "
            f"(default: {DEFAULT_OUTGROUP_BUDGET_PER_CLASS}; env: OUTGROUP_BUDGET_PER_CLASS)"
        ),
    )
    parser.add_argument(
        "--cluster-identity",
        type=float,
        default=float(os.environ.get("REF_CLUSTER_IDENTITY", "0.9")),
        metavar="FLOAT",
        help=("CD-HIT identity threshold for reference deduplication "
              "(default: 0.9; env: REF_CLUSTER_IDENTITY). 0.9 preserves the "
              "cross-lineage-informative (older, <90%) paralog structure that "
              "polarizes shared vs Berghia-specific expansions, collapsing only "
              "lineage-private near-identical recent duplicates. Berghia and "
              "anchors are never deduped; budget is enforced by balanced "
              "sampling, so CD-HIT need not thin aggressively."),
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
        "--anchor-fasta",
        default=os.environ.get("ANCHOR_FASTA", "references/anchors/anchor_set.fasta"),
        metavar="FASTA",
        help=(
            "anchor_set.fasta of characterized GPCR landmarks (build_anchor_set.py). "
            "Anchors of each class are injected must-keep, through CD-HIT "
            "(keep-annotated-on-collision), charged to the refs budget. Skipped if "
            "the file is absent. (default: references/anchors/anchor_set.fasta; "
            "env: ANCHOR_FASTA)"
        ),
    )
    parser.add_argument(
        "--anchor-tsv",
        default=os.environ.get("ANCHOR_TSV", "references/anchors/anchor_set.tsv"),
        metavar="TSV",
        help=(
            "anchor_set.tsv provenance (required with --anchor-fasta). "
            "(default: references/anchors/anchor_set.tsv; env: ANCHOR_TSV)"
        ),
    )
    parser.add_argument(
        "--ref-floor-per-species",
        type=int,
        default=int(os.environ.get("REF_FLOOR_PER_SPECIES", str(DEFAULT_REF_FLOOR_PER_SPECIES))),
        metavar="N",
        help=("Balanced sampling: minimum reference reps per species "
              f"(default: {DEFAULT_REF_FLOOR_PER_SPECIES}; env: REF_FLOOR_PER_SPECIES)"),
    )
    parser.add_argument(
        "--ref-cap-per-taxon",
        type=int,
        default=int(os.environ.get("REF_CAP_PER_TAXON", str(DEFAULT_REF_CAP_PER_TAXON))),
        metavar="N",
        help=("Balanced sampling: maximum reference reps per taxon "
              f"(default: {DEFAULT_REF_CAP_PER_TAXON}; env: REF_CAP_PER_TAXON)"),
    )
    parser.add_argument(
        "--selection-seed",
        type=int,
        default=int(os.environ.get("SELECTION_SEED", str(DEFAULT_SELECTION_SEED))),
        metavar="N",
        help=("Seed for reproducible reference selection "
              f"(default: {DEFAULT_SELECTION_SEED}; env: SELECTION_SEED)"),
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

    build_all_pools(
        scan_fasta_glob=args.scan_fasta_glob,
        class_tsv=args.class_tsv,
        out_dir=args.out_dir,
        total_budget_per_class=args.total_budget_per_class,
        outgroup_budget_per_class=args.outgroup_budget_per_class,
        cluster_identity=args.cluster_identity,
        cdhit_path=args.cdhit_path,
        threads=args.threads,
        must_include_taxids=must_include,
        berghia_taxid=args.berghia_taxid,
        berghia_fasta=args.berghia_fasta,
        berghia_class_tsv=args.berghia_class_tsv,
        anchor_fasta=args.anchor_fasta,
        anchor_tsv=args.anchor_tsv,
        ref_floor_per_species=args.ref_floor_per_species,
        ref_cap_per_taxon=args.ref_cap_per_taxon,
        selection_seed=args.selection_seed,
        force=args.force,
    )


if __name__ == "__main__":
    main()
