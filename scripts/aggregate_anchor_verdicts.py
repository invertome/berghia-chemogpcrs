#!/usr/bin/env python3
"""aggregate_anchor_verdicts.py — combine the per-class C3 anchor-divergence
verdicts (evaluate_anchor_divergence.py output) into a single, report-only
summary for the user to review.

Deliberately does NOT mutate config: the per-class keep/exclude recommendation
is surfaced (markdown + JSON) so the trees can be eyeballed before the operator
sets ANCHOR_OUTGROUP_CLASSES.

Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts
"""
from __future__ import annotations

import argparse
import glob
import json
import os
from pathlib import Path


def aggregate_verdicts(verdicts: list[dict]) -> dict:
    """Combine per-class verdict dicts into a report (report-only)."""
    per_class = sorted(verdicts, key=lambda v: v.get("class", ""))
    keep = sorted(v["class"] for v in verdicts if v.get("verdict") == "include")
    exclude = sorted(v["class"] for v in verdicts if v.get("verdict") == "exclude")
    return {
        "n_classes": len(verdicts),
        "keep_outgroup_anchors": keep,
        "exclude_outgroup_anchors": exclude,
        "per_class": per_class,
    }


def aggregate_dir(in_dir: str) -> dict:
    """Load all verdict_class_*.json under *in_dir* and aggregate them."""
    verdicts = []
    for path in sorted(glob.glob(os.path.join(in_dir, "verdict_class_*.json"))):
        with open(path) as fh:
            verdicts.append(json.load(fh))
    return aggregate_verdicts(verdicts)


def render_markdown(report: dict) -> str:
    """Render the report as a reviewable markdown summary."""
    lines = [
        "# Anchor divergence calibration (C3) — per-class verdicts",
        "",
        "Report-only. Review the trees, then set `ANCHOR_OUTGROUP_CLASSES` to the",
        "**keep** classes (out-group tier-2/3 anchors retained); excluded classes",
        "should be rebuilt with `--anchor-tiers 1` (in-group anchors only).",
        "",
        "| class | verdict | infiltrations | RF | support_drop | reasons |",
        "|-------|---------|---------------|----|--------------|---------|",
    ]
    for v in report["per_class"]:
        reasons = "; ".join(v.get("reasons", [])) or "—"
        lines.append(
            f"| {v.get('class','?')} | {v.get('verdict','?')} | "
            f"{v.get('n_infiltrations','?')} | {v.get('rf',0):.3f} | "
            f"{v.get('support_drop',0):.1f} | {reasons} |"
        )
    # Flagged placements — out-group anchors nested in a supported focal clade.
    flagged = []
    for v in report["per_class"]:
        for p in v.get("anchor_placements", []):
            if p.get("in_berghia_clade"):
                flagged.append((v.get("class", "?"), p))
    if flagged:
        lines += ["", "## Flagged placements (anchor nested in a supported focal clade)", ""]
        for cls, p in flagged:
            lines.append(
                f"- **{cls}** `{p['anchor']}` — sister to {p['sister_berghia']} "
                f"focal tip(s) (clade size {p['sister_size']}), support "
                f"{p['parent_support']:.0f}"
            )

    lines += [
        "",
        f"**Keep out-group anchors:** {', '.join(report['keep_outgroup_anchors']) or '(none)'}",
        f"**Exclude out-group anchors:** {', '.join(report['exclude_outgroup_anchors']) or '(none)'}",
        "",
        f"Recommended config: `ANCHOR_OUTGROUP_CLASSES=\"{','.join(report['keep_outgroup_anchors'])}\"`",
    ]
    return "\n".join(lines) + "\n"


def main(argv=None) -> None:
    p = argparse.ArgumentParser(
        description="Aggregate per-class anchor-divergence verdicts (report-only).")
    p.add_argument("--in-dir", required=True,
                   help="Directory containing verdict_class_*.json")
    p.add_argument("--out-json", required=True)
    p.add_argument("--out-md", required=True)
    args = p.parse_args(argv)

    report = aggregate_dir(args.in_dir)
    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_json, "w") as fh:
        json.dump(report, fh, indent=2)
    with open(args.out_md, "w") as fh:
        fh.write(render_markdown(report))
    print(f"[aggregate_anchor_verdicts] keep={report['keep_outgroup_anchors']} "
          f"exclude={report['exclude_outgroup_anchors']} -> {args.out_md}")


if __name__ == "__main__":
    main()
