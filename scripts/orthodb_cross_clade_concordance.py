#!/usr/bin/env python3
"""Does an OrthoDB orthogroup connect a characterized vertebrate receptor to a
characterized invertebrate receptor of the SAME curated family?

That is the exact inference the expansion route depends on. The harvest works
like this: seed from a characterized anchor (characterization lives
overwhelmingly in vertebrates), take its orthogroup, and harvest the
non-vertebrate members as candidate family assignments. The step that can
silently fail is the one in the middle - the orthogroup has to actually span
~550 Myr of divergence AND stay family-specific while doing it.

We have ground truth for exactly that step: the anchor set contains
invertebrate entries with curated families. So for every invertebrate anchor
that maps into OrthoDB we ask which of three things happened at a given level:

  AGREE        it shares an orthogroup with >=1 vertebrate anchor, and the
               modal curated family of those vertebrate anchors matches its
               own. The route would have made the right call.
  DISAGREE     it shares an orthogroup with vertebrate anchors of a DIFFERENT
               modal family. The route would have made a confident wrong call
               - the dangerous outcome.
  NOT_SPANNED  no vertebrate anchor shares its orthogroup at all. The route
               returns nothing for it. Not a wrong answer, but no answer, and
               for a family that genuinely exists in molluscs it means the
               orthogroup structure does not reach across the distance.

NOT_SPANNED is not automatically an indictment of OrthoDB: for families that
molluscs simply lack, it is the correct result. It is only damning for families
known to have genuine non-vertebrate representation.

A structural constraint dominates the level sweep and is worth stating up
front: in odb12v2 there is NO level between Metazoa and either Mollusca or
Vertebrata - Bilateria, Protostomia and Deuterostomia are not OrthoDB levels.
So the only levels at which a mollusc and a vertebrate can possibly share an
orthogroup are Metazoa (33208) and Eukaryota (2759). At Mollusca, Vertebrata or
Arthropoda level the test is vacuous by construction, because those levels
contain only their own clade's organisms. The sweep reports every level so this
is visible rather than assumed.

Usage:
    python3 scripts/orthodb_cross_clade_concordance.py \
        --indir results/ranking/diagnostics/orthodb
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

VERTEBRATA = "Vertebrata"
NON_METAZOA = "non_Metazoa"
UNKNOWN = "unknown"


def modal_family(families: list[str]) -> tuple[str, int, bool]:
    """Most common family, its count, and whether the mode is unique.

    A tie is reported rather than silently resolved by dict order, because
    'the vertebrate anchors disagree among themselves' is a different result
    from 'they agree on family X'.
    """
    if not families:
        raise ValueError("no families to take a mode of")
    counts = Counter(families)
    top = counts.most_common()
    best_n = top[0][1]
    tied = [f for f, n in top if n == best_n]
    return top[0][0], best_n, len(tied) == 1


def classify_anchor(
    inv_family: str, vert_families: list[str]
) -> tuple[str, str]:
    """-> (bucket, modal vertebrate family or '')."""
    if not vert_families:
        return "NOT_SPANNED", ""
    mode, _, unique = modal_family(vert_families)
    if not unique:
        return "AMBIGUOUS_VERTEBRATE_MODE", mode
    return ("AGREE" if mode == inv_family else "DISAGREE"), mode


def load(indir: Path) -> dict:
    need = ["mapping_audit.tsv", "anchor_gene_map.tsv", "anchor_og2genes.tsv",
            "anchor_og_meta.tsv"]
    missing = [f for f in need if not (indir / f).exists()]
    if missing:
        raise SystemExit(
            f"ERROR: missing {missing} in {indir}. Run the Unity join and then "
            "scripts/orthodb_mapping_audit.py first."
        )
    with open(indir / "mapping_audit.tsv", newline="") as fh:
        audit = {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}

    gene_to_accs: dict[str, set] = defaultdict(set)
    with open(indir / "anchor_gene_map.tsv", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2:
                gene_to_accs[row[1]].add(row[0])

    og_level: dict[str, str] = {}
    og_name: dict[str, str] = {}
    with open(indir / "anchor_og_meta.tsv", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2:
                og_level[row[0]] = row[1]
                og_name[row[0]] = row[2] if len(row) > 2 else ""

    # (level, accession) -> orthogroups
    assign: dict[tuple[str, str], set] = defaultdict(set)
    with open(indir / "anchor_og2genes.tsv", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 2:
                continue
            lvl = og_level.get(row[0])
            if lvl is None:
                continue
            for acc in gene_to_accs.get(row[1], ()):
                assign[(lvl, acc)].add(row[0])

    levels: dict[str, str] = {}
    lv = indir / "odb_levels.tsv"
    if lv.exists():
        with open(lv, newline="") as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if len(row) >= 2:
                    levels[row[0]] = row[1]
    return {"audit": audit, "assign": assign, "og_level": og_level,
            "og_name": og_name, "levels": levels}


def analyse_level(
    lvl: str, data: dict
) -> tuple[list[dict], dict, dict]:
    audit = data["audit"]
    assign = data["assign"]

    # orthogroup -> primary anchors in it at this level
    og_members: dict[str, list[dict]] = defaultdict(list)
    for (l, acc), ogs in assign.items():
        if l != lvl:
            continue
        rec = audit.get(acc)
        if rec is None or rec["use_primary"] != "True":
            continue
        for og in ogs:
            og_members[og].append(rec)

    per_anchor = []
    for (l, acc), ogs in assign.items():
        if l != lvl:
            continue
        rec = audit.get(acc)
        if rec is None or rec["use_primary"] != "True":
            continue
        clade = rec["clade"]
        if clade in (VERTEBRATA, NON_METAZOA, UNKNOWN):
            continue
        vert_fams: list[str] = []
        vert_accs: list[str] = []
        for og in ogs:
            for m in og_members[og]:
                if m["clade"] == VERTEBRATA and m["accession"] != acc:
                    vert_fams.append(m["family"])
                    vert_accs.append(m["accession"])
        bucket, mode = classify_anchor(rec["family"], vert_fams)
        per_anchor.append(
            {
                "level_taxid": lvl,
                "accession": acc,
                "clade": clade,
                "curated_family": rec["family"],
                "orthogroups": ";".join(sorted(ogs)),
                "n_vertebrate_anchors_cogrouped": len(vert_fams),
                "modal_vertebrate_family": mode,
                "bucket": bucket,
                "vertebrate_anchor_accessions": ";".join(sorted(set(vert_accs))[:12]),
            }
        )

    by_family: dict[str, Counter] = defaultdict(Counter)
    by_clade: dict[str, Counter] = defaultdict(Counter)
    confusion: dict[tuple[str, str], int] = Counter()
    for r in per_anchor:
        by_family[r["curated_family"]][r["bucket"]] += 1
        by_clade[r["clade"]][r["bucket"]] += 1
        if r["bucket"] in ("AGREE", "DISAGREE"):
            confusion[(r["curated_family"], r["modal_vertebrate_family"])] += 1

    summary = {
        "level_taxid": lvl,
        "level_name": data["levels"].get(lvl, "?"),
        "invertebrate_anchors_tested": len(per_anchor),
        "buckets": dict(Counter(r["bucket"] for r in per_anchor).most_common()),
        "by_family": {k: dict(v.most_common()) for k, v in by_family.items()},
        "by_clade": {k: dict(v.most_common()) for k, v in by_clade.items()},
        "confusion_matrix": {f"{a} -> {b}": n for (a, b), n in
                             sorted(confusion.items(), key=lambda kv: -kv[1])},
    }
    return per_anchor, summary, confusion


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--indir", default="results/ranking/diagnostics/orthodb")
    args = ap.parse_args(argv)
    indir = Path(args.indir)
    data = load(indir)

    levels = sorted({l for l, _ in data["assign"]})
    all_rows: list[dict] = []
    summaries: list[dict] = []
    for lvl in levels:
        rows, summ, _ = analyse_level(lvl, data)
        if summ["invertebrate_anchors_tested"] == 0:
            continue
        all_rows.extend(rows)
        summaries.append(summ)

    summaries.sort(key=lambda s: -s["invertebrate_anchors_tested"])

    if all_rows:
        out = indir / "cross_clade_concordance_per_anchor.tsv"
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(all_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(all_rows)
    (indir / "cross_clade_concordance_summary.json").write_text(
        json.dumps(summaries, indent=2) + "\n"
    )

    print("=== CROSS-CLADE FAMILY CONCORDANCE, PER LEVEL ===")
    print("(invertebrate anchors sharing an orthogroup with vertebrate anchors)\n")
    hdr = (f"{'level':<24}{'tested':>7}{'AGREE':>8}{'DISAGR':>8}"
           f"{'NOTSPAN':>9}{'AMBIG':>7}{'%agree_of_spanned':>19}")
    print(hdr); print("-" * len(hdr))
    for s in summaries:
        b = s["buckets"]
        a, d = b.get("AGREE", 0), b.get("DISAGREE", 0)
        ns, am = b.get("NOT_SPANNED", 0), b.get("AMBIGUOUS_VERTEBRATE_MODE", 0)
        spanned = a + d + am
        pct = f"{a/spanned*100:.1f}%" if spanned else "n/a"
        name = f"{s['level_name']}({s['level_taxid']})"
        print(f"{name:<24}{s['invertebrate_anchors_tested']:>7}{a:>8}{d:>8}{ns:>9}{am:>7}{pct:>19}")

    for s in summaries:
        if s["invertebrate_anchors_tested"] < 5:
            continue
        print(f"\n--- {s['level_name']} ({s['level_taxid']}) per curated family ---")
        print(f"{'family':<24}{'AGREE':>7}{'DISAGR':>8}{'NOTSPAN':>9}{'AMBIG':>7}")
        for fam, b in sorted(s["by_family"].items(), key=lambda kv: -sum(kv[1].values())):
            print(f"{fam:<24}{b.get('AGREE',0):>7}{b.get('DISAGREE',0):>8}"
                  f"{b.get('NOT_SPANNED',0):>9}{b.get('AMBIGUOUS_VERTEBRATE_MODE',0):>7}")
        print(f"{'-- by clade --':<24}")
        for cl, b in sorted(s["by_clade"].items(), key=lambda kv: -sum(kv[1].values())):
            print(f"{cl:<24}{b.get('AGREE',0):>7}{b.get('DISAGREE',0):>8}"
                  f"{b.get('NOT_SPANNED',0):>9}{b.get('AMBIGUOUS_VERTEBRATE_MODE',0):>7}")
        if s["confusion_matrix"]:
            print("  confusion (invertebrate family -> modal vertebrate family):")
            for k, n in list(s["confusion_matrix"].items())[:15]:
                print(f"    {k:<50} {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
