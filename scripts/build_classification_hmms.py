#!/usr/bin/env python3
"""build_classification_hmms.py — Build per-family HMMs from the curated
non-chemoreceptor reference set.

Phase 2 Task 2.1 of the non-chemoreceptor classification pipeline.

For each coarse family (10) AND each medium subfamily (aminergic +
peptide drill-downs), select the curated reference sequences, align
them with MAFFT (--auto), and run hmmbuild. Outputs:

    results/classification/hmms/<family>.hmm                 (coarse)
    results/classification/hmms/<family>.aln                  (alignment)
    results/classification/hmms/<family>_<subfamily>.hmm     (medium)
    results/classification/hmms/<family>_<subfamily>.aln
    results/classification/hmms/manifest.tsv                  (HMM -> family + n_train)

The HMM library is consumed by:
    - scripts/classify_via_hmm.py (Task 4.1) — runs hmmscan against the
      candidate set and assigns best-hit family + subfamily
    - scripts/validate_classification_hmms.py (Task 2.2) — leave-one-out
      cross-validation of the HMM library

Design notes:
    - We use plain MAFFT (--auto) here, NOT the full filter stack
      (PREQUAL/CLOAK/TAPER). HMMs are profile-based and robust to
      alignment uncertainty; aggressive filtering of training input
      would over-trim signal columns. The filter stack matters more for
      tree inference (where positional homology certainty IS critical).
    - Min training-set size is 10 sequences (tunable). Below that,
      hmmbuild produces unreliable HMMs (high E-value variance).
    - For aminergic and peptide, we build BOTH the coarse-family HMM
      AND per-subfamily HMMs. Classification picks the best hit across
      all HMMs (Task 4.1). Hierarchical: a candidate matching both
      'aminergic' and 'aminergic_5HT' gets the more specific label.

Usage:
    python3 build_classification_hmms.py \\
        --reference-fasta references/non_chemo_gpcr/all_references.fasta \\
        --reference-tsv references/non_chemo_gpcr/all_references.tsv \\
        --output-dir results/classification/hmms/ \\
        --threads 4
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

DEFAULT_MIN_TRAIN_SIZE = 10
SUBFAMILY_FAMILIES = {"aminergic", "peptide"}  # medium drill-down only here


# ---- Reference parsing -------------------------------------------------

def parse_reference_tsv(path: str) -> list[dict[str, str]]:
    """Read references/non_chemo_gpcr/all_references.tsv into a list of dicts."""
    records: list[dict[str, str]] = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            records.append({
                "accession": row.get("accession", "").strip(),
                "family": row.get("family", "").strip(),
                "subfamily": row.get("subfamily", "").strip(),
                "species": row.get("species", "").strip(),
                "gene": row.get("gene", "").strip(),
                "source": row.get("source", "").strip(),
                "length": row.get("length", "").strip(),
            })
    return records


def parse_reference_fasta(path: str) -> dict[str, str]:
    """Read references/non_chemo_gpcr/all_references.fasta — return dict
    of {accession: sequence} (accession is the first '|' field of the
    header)."""
    seqs: dict[str, str] = {}
    current: str | None = None
    parts: list[str] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    seqs[current] = "".join(parts)
                # accession is first '|'-separated field after '>'
                current = line[1:].split("|", 1)[0].strip()
                parts = []
            else:
                parts.append(line)
    if current is not None:
        seqs[current] = "".join(parts)
    return seqs


# ---- Grouping ----------------------------------------------------------

def group_by_family(records: list[dict[str, str]]) -> dict[str, list[str]]:
    """Group accession lists by coarse family. Empty/missing family is
    skipped (those are unclassified entries that shouldn't be in the
    consolidated set anyway)."""
    groups: dict[str, list[str]] = defaultdict(list)
    for r in records:
        fam = r.get("family", "")
        acc = r.get("accession", "")
        if fam and acc and fam != "unclassified-gpcr":
            groups[fam].append(acc)
    return dict(groups)


def group_by_subfamily(records: list[dict[str, str]]) -> dict[str, list[str]]:
    """Build subfamily groups ONLY for aminergic and peptide families.
    Records with empty subfamily are NOT included (they only contribute
    to the coarse-family HMM)."""
    groups: dict[str, list[str]] = defaultdict(list)
    for r in records:
        fam = r.get("family", "")
        sub = r.get("subfamily", "")
        acc = r.get("accession", "")
        if fam in SUBFAMILY_FAMILIES and sub and acc:
            key = f"{fam}/{sub}"
            groups[key].append(acc)
    return dict(groups)


def filter_by_min_size(groups: dict[str, list[str]],
                       min_size: int = DEFAULT_MIN_TRAIN_SIZE
                       ) -> dict[str, list[str]]:
    """Drop groups with fewer than min_size members."""
    return {k: v for k, v in groups.items() if len(v) >= min_size}


# ---- Output -----------------------------------------------------------

def write_family_fasta(accessions: list[str], seqs: dict[str, str],
                       path: str) -> int:
    """Write per-family input FASTA preserving the order in `accessions`.
    Returns the number of sequences written."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(path, "w") as f:
        for acc in accessions:
            seq = seqs.get(acc)
            if seq is None:
                print(f"WARN: no sequence for {acc}; skipping", file=sys.stderr)
                continue
            f.write(f">{acc}\n{seq}\n")
            n += 1
    return n


def manifest_row(family_label: str, subfamily_label: str, n_train: int,
                 hmm_path: str, align_path: str) -> dict[str, object]:
    """Canonical manifest row format."""
    return {
        "family": family_label,
        "subfamily": subfamily_label,
        "n_train": n_train,
        "hmm_path": hmm_path,
        "alignment_path": align_path,
    }


def write_manifest(rows: list[dict[str, object]], path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    cols = ["family", "subfamily", "n_train", "hmm_path", "alignment_path"]
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in rows:
            w.writerow([r[c] for c in cols])


# ---- Pipeline ---------------------------------------------------------

def run_mafft(input_fasta: str, output_alignment: str, threads: int = 4) -> None:
    """MAFFT --auto wrapper. Robust default for variable-size inputs."""
    cmd = ["mafft", "--auto", "--thread", str(threads), input_fasta]
    with open(output_alignment, "w") as out:
        result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE,
                                text=True, check=False)
    if result.returncode != 0 or not Path(output_alignment).stat().st_size:
        raise RuntimeError(
            f"mafft failed (rc={result.returncode}) for {input_fasta}: "
            f"{result.stderr[:500]}")


def run_hmmbuild(alignment: str, hmm_out: str, name: str) -> None:
    """hmmbuild wrapper."""
    cmd = ["hmmbuild", "-n", name, hmm_out, alignment]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0 or not Path(hmm_out).stat().st_size:
        raise RuntimeError(
            f"hmmbuild failed (rc={result.returncode}) for {alignment}: "
            f"{result.stderr[:500]}")


def build_hmm_for_group(group_key: str, accessions: list[str],
                       seqs: dict[str, str], output_dir: Path,
                       threads: int = 4) -> dict[str, object] | None:
    """Build alignment + HMM for a single family/subfamily group.
    Returns manifest row dict, or None on failure."""
    safe_key = group_key.replace("/", "_")
    fasta = output_dir / f"{safe_key}.fasta"
    aln = output_dir / f"{safe_key}.aln"
    hmm = output_dir / f"{safe_key}.hmm"

    n = write_family_fasta(accessions, seqs, str(fasta))
    if n < 2:
        print(f"WARN: only {n} sequences for {group_key}; skipping", file=sys.stderr)
        return None

    try:
        run_mafft(str(fasta), str(aln), threads=threads)
    except RuntimeError as e:
        print(f"ERROR ({group_key}): {e}", file=sys.stderr)
        return None

    try:
        run_hmmbuild(str(aln), str(hmm), name=safe_key)
    except RuntimeError as e:
        print(f"ERROR ({group_key}): {e}", file=sys.stderr)
        return None

    # Strip the working FASTA — alignment + HMM are the keep artifacts.
    fasta.unlink(missing_ok=True)

    if "/" in group_key:
        family, subfamily = group_key.split("/", 1)
    else:
        family, subfamily = group_key, ""
    return manifest_row(family, subfamily, n, str(hmm), str(aln))


# ---- Main -------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--reference-fasta", required=True)
    ap.add_argument("--reference-tsv", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--min-train-size", type=int, default=DEFAULT_MIN_TRAIN_SIZE,
                    help=f"Minimum sequences per HMM (default {DEFAULT_MIN_TRAIN_SIZE})")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    if shutil.which("mafft") is None:
        print("ERROR: mafft not on PATH", file=sys.stderr)
        return 2
    if shutil.which("hmmbuild") is None:
        print("ERROR: hmmbuild not on PATH", file=sys.stderr)
        return 2

    print("[hmm-build] Reading references...", file=sys.stderr)
    records = parse_reference_tsv(args.reference_tsv)
    seqs = parse_reference_fasta(args.reference_fasta)
    print(f"[hmm-build]   {len(records)} records, {len(seqs)} sequences", file=sys.stderr)

    fam_groups = group_by_family(records)
    sub_groups = group_by_subfamily(records)
    fam_groups = filter_by_min_size(fam_groups, args.min_train_size)
    sub_groups = filter_by_min_size(sub_groups, args.min_train_size)

    print(f"[hmm-build] Building HMMs for {len(fam_groups)} families "
          f"+ {len(sub_groups)} medium subfamilies...", file=sys.stderr)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest: list[dict[str, object]] = []

    for fam_key, accs in sorted(fam_groups.items()):
        print(f"[hmm-build] {fam_key} (n={len(accs)})...", file=sys.stderr)
        row = build_hmm_for_group(fam_key, accs, seqs, output_dir,
                                  threads=args.threads)
        if row:
            manifest.append(row)

    for sub_key, accs in sorted(sub_groups.items()):
        print(f"[hmm-build] {sub_key} (n={len(accs)})...", file=sys.stderr)
        row = build_hmm_for_group(sub_key, accs, seqs, output_dir,
                                  threads=args.threads)
        if row:
            manifest.append(row)

    write_manifest(manifest, str(output_dir / "manifest.tsv"))
    print(f"\n[hmm-build] DONE: {len(manifest)} HMMs in {output_dir}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
