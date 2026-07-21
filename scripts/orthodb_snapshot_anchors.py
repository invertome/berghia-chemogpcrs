#!/usr/bin/env python3
"""Freeze the class-A anchor snapshot used for the OrthoDB orthogroup validation.

The anchor set is being repaired concurrently (class-B/C receptors are being
removed from the class-A pool), so every downstream statistic MUST be pinned to
one identifiable snapshot. This script writes that snapshot plus a provenance
record carrying the sha256 and mtime of every input file, so the numbers in the
report can be re-derived or invalidated later.

Two label-quality filters are applied and both are recorded per row rather than
silently dropped:

  * ``pass_class_a`` from the refexp2 anchor audit. Rows the audit marks False
    are known-bad labels (no rhodopsin-clan domain). Scoring OrthoDB against
    them would understate agreement, so they are excluded from the primary
    statistic and reported separately.
  * ``family == "orphan"`` is not a curated family. It is a grab-bag of
    receptors with no assigned family, so it cannot conflate or fragment in any
    meaningful sense. Excluded from the primary statistic, reported separately.

Outputs (all under RESULTS/ranking/diagnostics/orthodb/):
  anchor_snapshot.tsv       one row per class-A anchor with usability flags
  anchor_accessions.txt     bare accession list, the join key for OrthoDB
  snapshot_provenance.json  sha256/mtime/row-counts of every input

Read-only with respect to references/; writes only under results/.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

# Families that carry no curated identity and therefore cannot be scored.
NON_FAMILY_LABELS = {"orphan", "", "unclassified", "unknown"}


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_audit_flags(path: Path) -> dict[str, bool]:
    """accession -> pass_class_a. Missing file is fatal, not a shrug.

    A missing audit would silently turn known-bad labels back into ground
    truth, which is precisely the failure this snapshot exists to prevent.
    """
    if not path.exists():
        raise SystemExit(
            f"ERROR: anchor audit not found at {path}. Refusing to build a "
            "snapshot without the class-A label audit, because that would "
            "score OrthoDB against labels already known to be wrong."
        )
    flags: dict[str, bool] = {}
    for row in read_tsv(path):
        acc = (row.get("accession") or "").strip()
        if not acc:
            continue
        raw = (row.get("pass_class_a") or "").strip().lower()
        flags[acc] = raw not in {"false", "0", "no"}
    return flags


def build_snapshot(anchors: list[dict], audit: dict[str, bool]) -> list[dict]:
    out = []
    for row in anchors:
        if (row.get("class") or "").strip() != "A":
            continue
        acc = (row.get("accession") or "").strip()
        family = (row.get("family") or "").strip()
        audited = acc in audit
        passes = audit.get(acc, True)
        scoreable_family = family.lower() not in NON_FAMILY_LABELS
        out.append(
            {
                "accession": acc,
                "family": family,
                "class": "A",
                "taxid": (row.get("taxid") or "").strip(),
                "species": (row.get("species") or "").strip(),
                "tier": (row.get("tier") or "").strip(),
                "evidence": (row.get("evidence") or "").strip(),
                "in_audit": str(audited),
                "pass_class_a": str(passes),
                "scoreable_family": str(scoreable_family),
                # The primary-analysis membership flag. One column, so the
                # analysis never has to re-derive the inclusion rule.
                "use_primary": str(passes and scoreable_family),
            }
        )
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--anchors", default="references/anchors/anchor_set_PROD.tsv")
    ap.add_argument(
        "--audit", default="references/non_chemo_gpcr/refexp2_anchor_audit.tsv"
    )
    ap.add_argument("--outdir", default="results/ranking/diagnostics/orthodb")
    args = ap.parse_args(argv)

    anchors_path = Path(args.anchors)
    audit_path = Path(args.audit)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    anchors = read_tsv(anchors_path)
    audit = load_audit_flags(audit_path)
    snapshot = build_snapshot(anchors, audit)

    if not snapshot:
        raise SystemExit("ERROR: no class-A anchors found; refusing to write an empty snapshot.")

    accs = [r["accession"] for r in snapshot]
    if len(set(accs)) != len(accs):
        dupes = [a for a, n in Counter(accs).items() if n > 1]
        raise SystemExit(f"ERROR: duplicate accessions in class-A anchors: {dupes[:10]}")

    snap_path = outdir / "anchor_snapshot.tsv"
    with open(snap_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(snapshot[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(snapshot)

    acc_path = outdir / "anchor_accessions.txt"
    acc_path.write_text("\n".join(accs) + "\n")

    primary = [r for r in snapshot if r["use_primary"] == "True"]
    prov = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            str(p): {
                "sha256": sha256_of(p),
                "mtime_utc": datetime.fromtimestamp(
                    os.path.getmtime(p), timezone.utc
                ).isoformat(),
                "bytes": p.stat().st_size,
            }
            for p in (anchors_path, audit_path)
        },
        "counts": {
            "anchor_rows_total": len(anchors),
            "class_a_total": len(snapshot),
            "class_a_in_audit": sum(1 for r in snapshot if r["in_audit"] == "True"),
            "class_a_failing_audit": sum(1 for r in snapshot if r["pass_class_a"] == "False"),
            "class_a_unscoreable_family": sum(
                1 for r in snapshot if r["scoreable_family"] == "False"
            ),
            "class_a_primary": len(primary),
        },
        "family_counts_primary": dict(Counter(r["family"] for r in primary).most_common()),
        "family_counts_all_class_a": dict(
            Counter(r["family"] for r in snapshot).most_common()
        ),
    }
    (outdir / "snapshot_provenance.json").write_text(json.dumps(prov, indent=2) + "\n")

    print(json.dumps(prov["counts"], indent=2))
    print(f"wrote {snap_path}")
    print(f"wrote {acc_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
