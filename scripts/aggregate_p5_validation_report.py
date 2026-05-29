"""Aggregate P5 Phase 1a validation results into a structured report.

Reads per-species scan records, per-source class TSVs, and the pool-builder
JSON report; emits both a machine-readable JSON and a human-readable Markdown
report.

Partial tolerance: if classify/ or pools/ outputs are absent (e.g. the
classify-and-pool job hasn't run yet), the report still emits with explicit
``null`` entries for the missing sections.  Callers should check the
``summary.missing_sections`` list.

Idempotent: skips both outputs if they already exist unless --force.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Optional

# ---------------------------------------------------------------------------
# Scan record parsing
# ---------------------------------------------------------------------------

def parse_scan_records(scan_dir: Path) -> list[dict]:
    """Parse all ``*.scan_record.tsv`` files in *scan_dir*.

    Returns a list of per-species dicts:
        taxid, binomial, n_input, n_hmm_positive, n_passed_gate

    Files are expected to be named ``<taxid>_<sanitized_binomial>.scan_record.tsv``.
    Missing or empty ``scan_dir`` returns an empty list.
    """
    if not scan_dir.is_dir():
        return []

    records: list[dict] = []
    for tsv_path in sorted(scan_dir.glob("*.scan_record.tsv")):
        stem = tsv_path.stem  # e.g. "6183_Schistosoma_mansoni.scan_record"
        # Remove trailing ".scan_record" from stem
        if stem.endswith(".scan_record"):
            stem = stem[: -len(".scan_record")]
        # Split on first underscore to get taxid vs binomial
        parts = stem.split("_", 1)
        taxid = parts[0] if parts else stem
        binomial = parts[1] if len(parts) > 1 else ""

        n_input = 0
        n_hmm_positive = 0
        n_passed_gate = 0
        try:
            with tsv_path.open() as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    n_input += 1
                    if str(row.get("gpcr_positive", "0")).strip() == "1":
                        n_hmm_positive += 1
                    if str(row.get("passed_gate", "0")).strip() == "1":
                        n_passed_gate += 1
        except OSError:
            pass

        records.append({
            "taxid": taxid,
            "binomial": binomial,
            "n_input": n_input,
            "n_hmm_positive": n_hmm_positive,
            "n_passed_gate": n_passed_gate,
        })

    return records


# ---------------------------------------------------------------------------
# Class distribution parsing
# ---------------------------------------------------------------------------

def parse_class_distribution(class_tsv: Path) -> Optional[dict[str, int]]:
    """Parse a ``classify_gpcr_by_class.py`` output TSV.

    Returns a dict mapping class label → count, or ``None`` if the file
    does not exist.
    """
    if not class_tsv.exists():
        return None

    dist: dict[str, int] = {}
    try:
        with class_tsv.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                cls = (row.get("class") or "unclassified").strip()
                dist[cls] = dist.get(cls, 0) + 1
    except OSError:
        return None

    return dist


# ---------------------------------------------------------------------------
# Pool report loading
# ---------------------------------------------------------------------------

def load_pool_report(report_path: Path) -> Optional[dict]:
    """Load the pool_build_report.json emitted by build_per_class_reference_pools.py.

    Returns the parsed dict, or ``None`` if the file does not exist.
    """
    if not report_path.exists():
        return None
    try:
        with report_path.open() as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError):
        return None


# ---------------------------------------------------------------------------
# Combined class distribution helper
# ---------------------------------------------------------------------------

def _merge_distributions(
    a: Optional[dict[str, int]],
    b: Optional[dict[str, int]],
) -> Optional[dict[str, int]]:
    """Merge two class-distribution dicts (element-wise sum).

    Returns ``None`` if both inputs are ``None``.
    """
    if a is None and b is None:
        return None
    combined: dict[str, int] = {}
    for d in (a or {}, b or {}):
        for k, v in d.items():
            combined[k] = combined.get(k, 0) + v
    return combined


# ---------------------------------------------------------------------------
# Main report assembly
# ---------------------------------------------------------------------------

def build_validation_report(
    *,
    scan_dir: Path,
    classify_dir: Path,
    pool_report_path: Path,
) -> dict:
    """Assemble the full validation report dict.

    Works even if classify/ or pool/ outputs are absent — those sections
    will be ``None`` and ``summary.missing_sections`` will list them.
    """
    per_species = parse_scan_records(scan_dir)

    phase1a_dist = parse_class_distribution(classify_dir / "class_phase1a.tsv")
    berghia_dist = parse_class_distribution(classify_dir / "class_berghia.tsv")
    combined_dist = _merge_distributions(phase1a_dist, berghia_dist)

    pool_data = load_pool_report(pool_report_path)

    missing_sections: list[str] = []
    if phase1a_dist is None:
        missing_sections.append("class_distribution.phase1a")
    if berghia_dist is None:
        missing_sections.append("class_distribution.berghia")
    if pool_data is None:
        missing_sections.append("pool_composition")

    # Pool composition: extract per-class stats from pool report
    pool_composition: Optional[dict] = None
    if pool_data is not None:
        pool_composition = {}
        for key, val in pool_data.items():
            if key.startswith("class_"):
                pool_composition[key] = val
        # Also pass through unclassified and berghia_included summaries
        if "unclassified" in pool_data:
            pool_composition["unclassified"] = pool_data["unclassified"]
        if "berghia_included" in pool_data:
            pool_composition["berghia_included"] = pool_data["berghia_included"]

    total_candidates = sum(r["n_passed_gate"] for r in per_species)

    return {
        "per_species": per_species,
        "class_distribution": {
            "phase1a": phase1a_dist,
            "berghia": berghia_dist,
            "combined": combined_dist,
        },
        "pool_composition": pool_composition,
        "summary": {
            "n_species": len(per_species),
            "total_candidates": total_candidates,
            "missing_sections": missing_sections,
        },
    }


# ---------------------------------------------------------------------------
# Report writers
# ---------------------------------------------------------------------------

def write_json_report(report: dict, out_path: Path) -> None:
    """Write the report dict as indented JSON."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        json.dump(report, fh, indent=2)


def write_md_report(report: dict, out_path: Path) -> None:
    """Write a human-readable Markdown report."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []

    # Header
    lines.append("# P5 Phase 1a Validation Report\n")

    # Summary
    summary = report.get("summary", {})
    lines.append("## Summary\n")
    lines.append(f"- Species scanned: {summary.get('n_species', 0)}")
    lines.append(f"- Total gate-passed candidates: {summary.get('total_candidates', 0)}")
    missing = summary.get("missing_sections", [])
    if missing:
        lines.append(f"- Missing sections: {', '.join(missing)}")
    lines.append("")

    # Per-species scan stats
    lines.append("## Per-species Scan Stats\n")
    per_species = report.get("per_species", [])
    if per_species:
        lines.append("| taxid | binomial | n_input | n_hmm_positive | n_passed_gate |")
        lines.append("|-------|----------|---------|----------------|---------------|")
        for sp in sorted(per_species, key=lambda x: x.get("taxid", "")):
            lines.append(
                f"| {sp['taxid']} | {sp['binomial']} "
                f"| {sp['n_input']} | {sp['n_hmm_positive']} "
                f"| {sp['n_passed_gate']} |"
            )
    else:
        lines.append("*No scan records found.*")
    lines.append("")

    # Class distribution
    lines.append("## Class Distribution per Source\n")
    class_dist = report.get("class_distribution", {})
    for source_label, source_key in [
        ("Phase 1a candidates", "phase1a"),
        ("Berghia candidates", "berghia"),
        ("Combined", "combined"),
    ]:
        dist = class_dist.get(source_key)
        if dist is None:
            lines.append(f"**{source_label}**: *missing*\n")
            continue
        total = sum(dist.values()) or 1  # avoid division by zero
        lines.append(f"**{source_label}** (total={sum(dist.values())}):\n")
        lines.append("| class | count | % |")
        lines.append("|-------|-------|---|")
        for cls in sorted(dist.keys()):
            pct = 100.0 * dist[cls] / total
            lines.append(f"| {cls} | {dist[cls]} | {pct:.1f} |")
        lines.append("")

    # Per-class pool composition
    lines.append("## Per-class Pool Composition\n")
    pool = report.get("pool_composition")
    if pool is None:
        lines.append("*Pool build report not available.*\n")
    else:
        lines.append("| class | n_berghia | n_must_include | n_refs_target "
                     "| n_output | n_subsampled |")
        lines.append("|-------|-----------|----------------|---------------|"
                     "----------|--------------|")
        for key in sorted(pool.keys()):
            if not key.startswith("class_"):
                continue
            cls = key[len("class_"):]
            stats = pool[key]
            lines.append(
                f"| {cls} "
                f"| {stats.get('n_berghia_included', '')} "
                f"| {stats.get('n_must_include', '')} "
                f"| {stats.get('n_refs_target', '')} "
                f"| {stats.get('n_output', '')} "
                f"| {stats.get('n_subsampled', '')} |"
            )
    lines.append("")

    out_path.write_text("\n".join(lines))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Aggregate P5 Phase 1a validation results from scan records, "
            "class TSVs, and pool-builder JSON into a structured report. "
            "Emits both JSON and Markdown outputs."
        ),
    )
    p.add_argument(
        "--scan-dir", type=Path,
        default=Path("results/p5_phase1a_validation/scan"),
        help="Directory containing *.scan_record.tsv files.",
    )
    p.add_argument(
        "--classify-dir", type=Path,
        default=Path("results/p5_phase1a_validation/classify"),
        help="Directory containing class_phase1a.tsv and class_berghia.tsv.",
    )
    p.add_argument(
        "--pool-report", type=Path,
        default=Path("results/p5_phase1a_validation/pools/pool_build_report.json"),
        help="Path to pool_build_report.json from build_per_class_reference_pools.py.",
    )
    p.add_argument(
        "--out-json", type=Path,
        required=True,
        help="Output path for machine-readable JSON report.",
    )
    p.add_argument(
        "--out-md", type=Path,
        required=True,
        help="Output path for human-readable Markdown report.",
    )
    p.add_argument(
        "--force", action="store_true",
        help="Re-run even if output files already exist.",
    )
    return p


def main(argv: Optional[list[str]] = None) -> int:
    args = _build_argparser().parse_args(argv)

    # Idempotency: skip if BOTH outputs already exist and --force not set
    if (
        args.out_json.exists()
        and args.out_md.exists()
        and not args.force
    ):
        print(
            f"[aggregate_p5_validation_report] Outputs already exist at "
            f"{args.out_json} and {args.out_md}; skipping (use --force to overwrite).",
            file=sys.stderr,
        )
        return 0

    report = build_validation_report(
        scan_dir=args.scan_dir,
        classify_dir=args.classify_dir,
        pool_report_path=args.pool_report,
    )

    write_json_report(report, args.out_json)
    write_md_report(report, args.out_md)

    summary = report["summary"]
    print(
        f"[aggregate_p5_validation_report] Report written: "
        f"{summary['n_species']} species, "
        f"{summary['total_candidates']} candidates. "
        f"JSON: {args.out_json}  MD: {args.out_md}",
        file=sys.stderr,
    )
    if summary["missing_sections"]:
        print(
            f"[aggregate_p5_validation_report] WARNING: missing sections: "
            f"{', '.join(summary['missing_sections'])}",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
