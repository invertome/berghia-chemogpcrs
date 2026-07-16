#!/usr/bin/env python3
"""Parse scratch_lofo_bakeoff.py stdout logs into a tidy bake-off metrics table.

The proven harness PRINTS a fixed-width table per (model, refset); it is not
touched here. Parsing its logs (rather than patching it to emit TSV) works on
runs already completed and keeps the validated metric code untouched.

One row per (model, refset, scoring config):
  model refset n_refs n_families dim score calib cent lofo_auroc invariance loo_famacc

Usage:
  embedding_bakeoff_metrics.py logs/harness_PROD-*.out --out results/ranking/diagnostics/bakeoff_metrics.tsv
"""
from __future__ import annotations

import argparse
import glob
import re
from typing import Dict, List

BAKEOFF_COLUMNS = (
    "model", "refset", "n_refs", "n_families", "dim",
    "score", "calib", "cent", "lofo_auroc", "invariance", "loo_famacc",
)

# "########## protrek :: FULL-clean ##########"  (refset label is FULL / VERIFIED,
# with an optional -clean / -nr / -723 style suffix the harness variants append).
_BLOCK = re.compile(
    r"#{3,}\s+(?P<model>\S+)\s+::\s+(?P<refset>FULL|VERIFIED)\S*\s+#{3,}"
    r"(?P<body>.*?)(?=#{3,}\s+\S+\s+::|\Z)",
    re.S,
)
_META = re.compile(r"(?P<refs>\d+)\s+refs,\s+(?P<fam>\d+)\s+families,.*?dim\s+(?P<dim>\d+)")
_ROW = re.compile(
    r"^(?P<score>cos|knn|maha|rel_maha)\s+(?P<calib>\w+)\s+(?P<cent>\w+)\s+"
    r"(?P<auroc>[\d.]+)\s+(?P<inv>[\d.]+|nan)\s+(?P<famacc>[\d.]+)\s*$",
    re.M,
)


def _f(x: str) -> float:
    return float("nan") if x == "nan" else float(x)


def parse_harness_log(text: str) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for blk in _BLOCK.finditer(text):
        body = blk.group("body")
        meta = _META.search(body)
        n_refs = int(meta.group("refs")) if meta else 0
        n_fam = int(meta.group("fam")) if meta else 0
        dim = int(meta.group("dim")) if meta else 0
        for m in _ROW.finditer(body):
            rows.append({
                "model": blk.group("model"), "refset": blk.group("refset"),
                "n_refs": n_refs, "n_families": n_fam, "dim": dim,
                "score": m.group("score"), "calib": m.group("calib"), "cent": m.group("cent"),
                "lofo_auroc": _f(m.group("auroc")), "invariance": _f(m.group("inv")),
                "loo_famacc": _f(m.group("famacc")),
            })
    return rows


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("logs", nargs="+", help="harness .out log files or globs")
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)

    paths: List[str] = []
    for g in a.logs:
        paths.extend(sorted(glob.glob(g)) or [g])
    rows: List[Dict[str, object]] = []
    for path in paths:
        with open(path) as fh:
            rows.extend(parse_harness_log(fh.read()))

    with open(a.out, "w") as fh:
        fh.write("\t".join(BAKEOFF_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in BAKEOFF_COLUMNS) + "\n")
    models = sorted({r["model"] for r in rows})
    print(f"parsed {len(rows)} rows from {len(paths)} logs, {len(models)} models -> {a.out}")


if __name__ == "__main__":
    main()
