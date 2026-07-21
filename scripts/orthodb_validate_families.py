#!/usr/bin/env python3
"""Test whether OrthoDB orthogroups recover curated GPCR family structure.

The question this answers is narrow and empirical: if a molluscan sequence
shares an OrthoDB orthogroup with several experimentally characterized
vertebrate receptors of a known family, is that structured evidence of family
membership? It only is if orthogroups track curated families in the first
place. So we take anchors whose family is known from curation, map them into
OrthoDB, and measure agreement at every taxonomic level OrthoDB exposes.

Two failure modes are measured separately because they have opposite cures:

  CONFLATION    several curated families sharing one orthogroup. A level that
                conflates is too coarse: membership in the group tells you
                nothing about which family a sequence belongs to. Measured as
                the fraction of anchors sitting in a mixed-family orthogroup,
                and summarised by *homogeneity*.

  FRAGMENTATION one curated family scattered across many orthogroups. Some is
                expected (families are large paralogous expansions and OrthoDB
                splits by orthology, not by family), so this is only fatal in
                the extreme. Measured as orthogroups-per-anchor within a
                family, and summarised by *completeness*.

Headline statistic is the ADJUSTED RAND INDEX. It is chance-corrected, which
matters here because the anchor families are wildly unbalanced (peptide alone
is ~45% of the set) and an uncorrected index would score high for a trivial
partition. Homogeneity and completeness are reported alongside because they
decompose ARI's single number into exactly the conflation and fragmentation
axes the question asks about; V-measure is their harmonic mean.

All three are computed here rather than imported from scikit-learn so the
implementation is unit-testable in isolation and the script has no heavyweight
dependency. The unit tests cross-check against sklearn when it is available.

Usage:
    python3 scripts/orthodb_validate_families.py \
        --indir results/ranking/diagnostics/orthodb
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

# Minimum anchors at a level before its statistics mean anything. Below this,
# conflation and fragmentation are dominated by sampling noise.
MIN_ANCHORS_PER_LEVEL = 20

# OrthoDB level taxids used to bin orthogroup members into clades. Verified
# present in odb12v2_levels.tab: Bilateria, Protostomia, Deuterostomia,
# Chordata and Ecdysozoa are NOT OrthoDB levels, which is itself a finding.
CLADE_LEVELS = [
    ("Vertebrata", 7742),
    ("Mollusca", 6447),
    ("other_Lophotrochozoa", 1206795),
    ("Arthropoda", 6656),
    ("Nematoda", 6231),
    ("Cnidaria", 6073),
]
METAZOA = 33208


# --------------------------------------------------------------- statistics --

def _comb2(n: int) -> int:
    return n * (n - 1) // 2


def contingency(labels_a: list, labels_b: list) -> dict:
    """Joint counts of two labelings over the same items."""
    if len(labels_a) != len(labels_b):
        raise ValueError("labelings must be the same length")
    table: dict = defaultdict(int)
    for a, b in zip(labels_a, labels_b):
        table[(a, b)] += 1
    return dict(table)


def adjusted_rand_index(labels_true: list, labels_pred: list) -> float:
    """Chance-corrected agreement between two partitions.

    Returns 1.0 for identical partitions, ~0.0 for independent ones, and can
    go negative when agreement is worse than chance.
    """
    n = len(labels_true)
    if n != len(labels_pred):
        raise ValueError("labelings must be the same length")
    if n < 2:
        raise ValueError("adjusted Rand index is undefined for fewer than 2 items")

    table = contingency(labels_true, labels_pred)
    a_counts = Counter(labels_true)
    b_counts = Counter(labels_pred)

    sum_ij = sum(_comb2(v) for v in table.values())
    sum_a = sum(_comb2(v) for v in a_counts.values())
    sum_b = sum(_comb2(v) for v in b_counts.values())
    total = _comb2(n)

    expected = (sum_a * sum_b) / total if total else 0.0
    maximum = 0.5 * (sum_a + sum_b)
    denom = maximum - expected
    if denom == 0:
        # Both partitions are trivial (all-singletons or all-one-cluster).
        # Convention matches sklearn: perfectly matching trivial partitions
        # score 1.0.
        return 1.0
    return (sum_ij - expected) / denom


def _entropy(counts: list[int], n: int) -> float:
    if n <= 0:
        return 0.0
    h = 0.0
    for c in counts:
        if c > 0:
            p = c / n
            h -= p * math.log(p)
    return h


def homogeneity_completeness_vmeasure(
    labels_true: list, labels_pred: list
) -> tuple[float, float, float]:
    """Rosenberg & Hirschberg (2007) cluster-evaluation triple.

    homogeneity  1.0 iff every cluster contains members of a single class.
                 This is the direct inverse of CONFLATION.
    completeness 1.0 iff every class is contained in a single cluster.
                 This is the direct inverse of FRAGMENTATION.
    v_measure    harmonic mean of the two.
    """
    n = len(labels_true)
    if n != len(labels_pred):
        raise ValueError("labelings must be the same length")
    if n == 0:
        raise ValueError("cannot evaluate an empty labeling")

    table = contingency(labels_true, labels_pred)
    true_counts = Counter(labels_true)
    pred_counts = Counter(labels_pred)

    h_c = _entropy(list(true_counts.values()), n)
    h_k = _entropy(list(pred_counts.values()), n)

    # H(C|K) and H(K|C)
    h_c_given_k = 0.0
    h_k_given_c = 0.0
    for (c, k), nck in table.items():
        p = nck / n
        h_c_given_k -= p * math.log(nck / pred_counts[k])
        h_k_given_c -= p * math.log(nck / true_counts[c])

    homogeneity = 1.0 if h_c == 0 else 1.0 - h_c_given_k / h_c
    completeness = 1.0 if h_k == 0 else 1.0 - h_k_given_c / h_k
    if homogeneity + completeness == 0:
        v = 0.0
    else:
        v = 2 * homogeneity * completeness / (homogeneity + completeness)
    return homogeneity, completeness, v


def conflation_stats(pairs: list[tuple[str, str]]) -> dict:
    """How much do distinct curated families share one orthogroup?

    pairs: (family, orthogroup) per anchor.
    """
    by_og: dict[str, Counter] = defaultdict(Counter)
    for fam, og in pairs:
        by_og[og][fam] += 1

    mixed = {og: fams for og, fams in by_og.items() if len(fams) > 1}
    anchors_in_mixed = sum(sum(f.values()) for f in mixed.values())
    total = len(pairs)

    return {
        "orthogroups_total": len(by_og),
        "orthogroups_mixed_family": len(mixed),
        "frac_orthogroups_mixed": len(mixed) / len(by_og) if by_og else 0.0,
        "anchors_in_mixed_orthogroup": anchors_in_mixed,
        "frac_anchors_in_mixed_orthogroup": anchors_in_mixed / total if total else 0.0,
        "max_families_in_one_orthogroup": max(
            (len(f) for f in by_og.values()), default=0
        ),
        "mean_families_per_orthogroup": (
            sum(len(f) for f in by_og.values()) / len(by_og) if by_og else 0.0
        ),
    }


def fragmentation_stats(pairs: list[tuple[str, str]]) -> dict:
    """How many orthogroups does one curated family span?

    Raw orthogroup counts scale with family size, so the normalised measure
    (orthogroups per anchor) is the comparable one: 1.0 means every anchor sits
    in its own orthogroup, 1/n means the family is perfectly contained.
    """
    by_family: dict[str, set] = defaultdict(set)
    size: Counter = Counter()
    for fam, og in pairs:
        by_family[fam].add(og)
        size[fam] += 1

    per_family = {
        fam: {
            "anchors": size[fam],
            "orthogroups": len(ogs),
            "orthogroups_per_anchor": len(ogs) / size[fam],
        }
        for fam, ogs in by_family.items()
    }
    ratios = [v["orthogroups_per_anchor"] for v in per_family.values()]
    ratios.sort()
    median = (
        0.0
        if not ratios
        else (
            ratios[len(ratios) // 2]
            if len(ratios) % 2
            else (ratios[len(ratios) // 2 - 1] + ratios[len(ratios) // 2]) / 2
        )
    )
    return {
        "families_total": len(by_family),
        "median_orthogroups_per_anchor": median,
        "max_orthogroups_per_anchor": max(ratios, default=0.0),
        "per_family": per_family,
    }


# --------------------------------------------------------------------- I/O --

def read_tsv(path: Path, fieldnames=None) -> list[dict]:
    with open(path, newline="") as fh:
        if fieldnames:
            return list(csv.DictReader(fh, delimiter="\t", fieldnames=fieldnames))
        return list(csv.DictReader(fh, delimiter="\t"))


def load_inputs(indir: Path) -> dict:
    need = [
        "anchor_snapshot.tsv",
        "anchor_gene_map.tsv",
        "anchor_og2genes.tsv",
        "anchor_og_meta.tsv",
        "odb_levels.tsv",
    ]
    missing = [f for f in need if not (indir / f).exists()]
    if missing:
        raise SystemExit(
            f"ERROR: missing inputs in {indir}: {missing}. "
            "Run scripts/unity/orthodb_fetch_and_join.sh first."
        )

    snapshot = read_tsv(indir / "anchor_snapshot.tsv")
    gene_map = read_tsv(
        indir / "anchor_gene_map.tsv", ["accession", "odb_gene_id", "raw_xref"]
    )
    og2genes = read_tsv(indir / "anchor_og2genes.tsv", ["og_id", "odb_gene_id"])
    og_meta = read_tsv(indir / "anchor_og_meta.tsv", ["og_id", "level_taxid", "og_name"])
    levels = read_tsv(
        indir / "odb_levels.tsv",
        ["taxid", "name", "n_genes", "n_ogs", "n_species"],
    )
    return {
        "snapshot": snapshot,
        "gene_map": gene_map,
        "og2genes": og2genes,
        "og_meta": og_meta,
        "levels": {r["taxid"]: r for r in levels},
    }


def assert_join(snapshot: list[dict], gene_map: list[dict]) -> dict:
    """Cross-namespace join UniProt -> OrthoDB gene. Never assume it worked.

    Reports the mapping rate as a first-class finding, because every downstream
    statistic is computed on whatever fraction actually mapped.
    """
    anchors = {r["accession"] for r in snapshot}
    primary = {r["accession"] for r in snapshot if r["use_primary"] == "True"}
    mapped_by_acc: dict[str, set] = defaultdict(set)
    for r in gene_map:
        if r["accession"] in anchors:
            mapped_by_acc[r["accession"]].add(r["odb_gene_id"])

    stray = {r["accession"] for r in gene_map} - anchors
    multi = {a: g for a, g in mapped_by_acc.items() if len(g) > 1}

    report = {
        "anchors_class_a": len(anchors),
        "anchors_primary": len(primary),
        "anchors_mapped_any": len(mapped_by_acc),
        "anchors_primary_mapped": len(primary & set(mapped_by_acc)),
        "mapping_rate_class_a": len(mapped_by_acc) / len(anchors) if anchors else 0.0,
        "mapping_rate_primary": (
            len(primary & set(mapped_by_acc)) / len(primary) if primary else 0.0
        ),
        "anchors_mapping_to_multiple_odb_genes": len(multi),
        "accessions_in_gene_map_not_in_snapshot": len(stray),
    }
    if report["anchors_mapped_any"] == 0:
        raise SystemExit(
            "ERROR: zero anchors mapped to OrthoDB gene ids. The join key is "
            "wrong. Refusing to emit statistics over an empty intersection."
        )
    if stray:
        raise SystemExit(
            f"ERROR: {len(stray)} accessions in the gene map are absent from the "
            "anchor snapshot; the two files are out of sync."
        )
    return report


def per_level_analysis(data: dict) -> tuple[list[dict], dict]:
    snapshot = {r["accession"]: r for r in data["snapshot"]}
    gene_to_accs: dict[str, set] = defaultdict(set)
    for r in data["gene_map"]:
        gene_to_accs[r["odb_gene_id"]].add(r["accession"])

    og_level = {r["og_id"]: r["level_taxid"] for r in data["og_meta"]}
    og_name = {r["og_id"]: r.get("og_name", "") for r in data["og_meta"]}

    # (level, accession) -> set of orthogroups
    assign: dict[tuple[str, str], set] = defaultdict(set)
    for r in data["og2genes"]:
        og = r["og_id"]
        lvl = og_level.get(og)
        if lvl is None:
            continue
        for acc in gene_to_accs.get(r["odb_gene_id"], ()):
            assign[(lvl, acc)].add(og)

    levels_seen = sorted({lvl for lvl, _ in assign})
    rows = []
    details: dict = {}

    for lvl in levels_seen:
        accs = [a for (l, a) in assign if l == lvl]
        primary_pairs = []
        ambiguous = 0
        for acc in accs:
            rec = snapshot.get(acc)
            if rec is None or rec["use_primary"] != "True":
                continue
            ogs = assign[(lvl, acc)]
            if len(ogs) != 1:
                ambiguous += 1
                continue
            primary_pairs.append((rec["family"], next(iter(ogs))))

        n = len(primary_pairs)
        lvl_meta = data["levels"].get(lvl, {})
        row = {
            "level_taxid": lvl,
            "level_name": lvl_meta.get("name", "?"),
            "level_species": lvl_meta.get("n_species", ""),
            "anchors_assigned": len(accs),
            "anchors_scored": n,
            "anchors_ambiguous_multi_og": ambiguous,
        }
        if n < MIN_ANCHORS_PER_LEVEL:
            row.update(
                {
                    "families": len(set(f for f, _ in primary_pairs)),
                    "orthogroups": len(set(o for _, o in primary_pairs)),
                    "note": f"below MIN_ANCHORS_PER_LEVEL={MIN_ANCHORS_PER_LEVEL}",
                }
            )
            rows.append(row)
            continue

        fams = [f for f, _ in primary_pairs]
        ogs = [o for _, o in primary_pairs]
        conf = conflation_stats(primary_pairs)
        frag = fragmentation_stats(primary_pairs)
        ari = adjusted_rand_index(fams, ogs)
        hom, comp, v = homogeneity_completeness_vmeasure(fams, ogs)

        row.update(
            {
                "families": len(set(fams)),
                "orthogroups": conf["orthogroups_total"],
                "orthogroups_mixed_family": conf["orthogroups_mixed_family"],
                "frac_anchors_in_mixed_og": round(
                    conf["frac_anchors_in_mixed_orthogroup"], 4
                ),
                "max_families_in_one_og": conf["max_families_in_one_orthogroup"],
                "median_ogs_per_anchor": round(frag["median_orthogroups_per_anchor"], 4),
                "adjusted_rand_index": round(ari, 4),
                "homogeneity": round(hom, 4),
                "completeness": round(comp, 4),
                "v_measure": round(v, 4),
                "note": "",
            }
        )
        rows.append(row)

        by_og: dict[str, Counter] = defaultdict(Counter)
        for f, o in primary_pairs:
            by_og[o][f] += 1
        worst = sorted(by_og.items(), key=lambda kv: -len(kv[1]))[:10]
        details[lvl] = {
            "level_name": lvl_meta.get("name", "?"),
            "conflation": conf,
            "fragmentation": frag,
            "worst_conflated_orthogroups": [
                {
                    "og_id": og,
                    "og_name": og_name.get(og, ""),
                    "families": dict(c.most_common()),
                    "anchors": sum(c.values()),
                }
                for og, c in worst
                if len(c) > 1
            ],
        }

    rows.sort(key=lambda r: -r["anchors_scored"])
    return rows, details


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--indir", default="results/ranking/diagnostics/orthodb")
    args = ap.parse_args(argv)
    indir = Path(args.indir)

    data = load_inputs(indir)
    join = assert_join(data["snapshot"], data["gene_map"])
    print("=== JOIN ASSERTION ===")
    print(json.dumps(join, indent=2))

    rows, details = per_level_analysis(data)

    cols = [
        "level_taxid", "level_name", "level_species", "anchors_assigned",
        "anchors_scored", "anchors_ambiguous_multi_og", "families", "orthogroups",
        "orthogroups_mixed_family", "frac_anchors_in_mixed_og",
        "max_families_in_one_og", "median_ogs_per_anchor", "adjusted_rand_index",
        "homogeneity", "completeness", "v_measure", "note",
    ]
    out_tsv = indir / "per_level_agreement.tsv"
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})

    (indir / "per_level_detail.json").write_text(
        json.dumps({"join": join, "levels": details}, indent=2) + "\n"
    )

    print("\n=== PER-LEVEL AGREEMENT ===")
    hdr = f"{'level':<22}{'scored':>7}{'fams':>6}{'OGs':>6}{'mixOG':>7}{'%mixed':>8}{'ARI':>8}{'homog':>8}{'compl':>8}"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        if r.get("anchors_scored", 0) < MIN_ANCHORS_PER_LEVEL:
            continue
        name = f"{r['level_name']}({r['level_taxid']})"
        print(
            f"{name:<22}{r['anchors_scored']:>7}{r['families']:>6}{r['orthogroups']:>6}"
            f"{r['orthogroups_mixed_family']:>7}{r['frac_anchors_in_mixed_og']*100:>7.1f}%"
            f"{r['adjusted_rand_index']:>8.3f}{r['homogeneity']:>8.3f}{r['completeness']:>8.3f}"
        )
    print(f"\nwrote {out_tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
