#!/usr/bin/env python3
"""compute_tandem_clusters.py — Find intra-genome tandem clusters of GPCR paralogs.

Bead -ar8. The intra-genome tandem-cluster signal is the field's signature
evidence type for chemoreceptor claims (Cummins 2009 *Aplysia*, Vogt 2023
gastropod survey, Robertson 2015 insect superfamily, Nath 2025 lophotrochozoa).
For each candidate gene, this script reports:
  - tandem_cluster_size: how many candidates fall within W kb on the same
    scaffold (1 if isolated)
  - tandem_cluster_id:   group label (scaffold:cluster_index) so paralogs in
    the same cluster can be tracked downstream

Uses gffutils (https://daler.github.io/gffutils/) for GFF parsing. The
sliding-window cluster algorithm itself is small enough (~30 lines) that
adding a heavy dependency like JCVI just for this is not worth it.

Usage:
    python compute_tandem_clusters.py \\
        --gff genomes/${BERGHIA_FILE_PREFIX}.gff3 \\
        --candidates results/chemogpcrs/chemogpcrs_berghia.fa \\
        --window-kb 100 --min-size 3 \\
        --out results/synteny/tandem_clusters.csv

The candidates file may be either a FASTA (gene IDs are the headers) or a
plain text file with one gene ID per line.
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Iterator

import gffutils  # battle-tested GFF parser; do not roll our own
import pandas as pd


def load_candidate_ids(path: str) -> set[str]:
    """Read candidate IDs from a FASTA (>header) or one-per-line text file.

    Detects FASTA vs plain text from the first non-empty line. In FASTA mode,
    only header lines contribute IDs (sequence lines are skipped). In plain
    text mode, every non-empty line is an ID; a leading '>' is tolerated.
    """
    with open(path) as f:
        lines = f.readlines()

    first_content = next(
        (line.strip() for line in lines if line.strip()), ""
    )
    is_fasta = first_content.startswith(">")
    ids: set[str] = set()
    if is_fasta:
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    else:
        for line in lines:
            line = line.strip()
            if not line:
                continue
            ids.add(line.lstrip(">").split()[0])
    return ids


def build_or_load_gff_db(gff_path: str, db_path: str | None = None) -> gffutils.FeatureDB:
    """Build (or reuse) a gffutils SQLite DB for fast GFF lookups."""
    db_path = db_path or f"{gff_path}.db"
    try:
        return gffutils.FeatureDB(db_path, keep_order=True)
    except Exception:
        return gffutils.create_db(
            gff_path, dbfn=db_path, force=True,
            merge_strategy="merge", keep_order=True,
            disable_infer_genes=True, disable_infer_transcripts=True,
        )


def iter_genes(db: gffutils.FeatureDB) -> Iterator[tuple[str, str, int, int]]:
    """Yield (gene_id, scaffold, start, end) for every gene in the DB."""
    for gene in db.features_of_type("gene"):
        yield (gene.id, gene.seqid, int(gene.start), int(gene.end))


def find_tandem_clusters(
    candidate_genes: Iterable[tuple[str, str, int, int]],
    *,
    window_kb: int = 100,
    min_size: int = 3,
) -> dict[str, tuple[int, str | None]]:
    """For each candidate gene, return (cluster_size, cluster_id).

    Cluster_id is "<scaffold>:cluster_<i>" or None when the gene is isolated.
    Two candidates belong to the same cluster iff they are on the same
    scaffold and within ``window_kb`` of consecutive cluster members.
    """
    by_scaffold: dict[str, list[tuple[str, str, int, int]]] = defaultdict(list)
    all_ids: list[str] = []
    for gene in candidate_genes:
        by_scaffold[gene[1]].append(gene)
        all_ids.append(gene[0])

    out: dict[str, tuple[int, str | None]] = {gid: (1, None) for gid in all_ids}
    cluster_counter: dict[str, int] = defaultdict(int)

    for scaff, genes in by_scaffold.items():
        if len(genes) < min_size:
            continue
        genes.sort(key=lambda x: x[2])
        i = 0
        while i < len(genes):
            j = i
            # Greedy extend the cluster while next gene is within window of
            # the previous one (use end-to-start gap, not span from i).
            while j + 1 < len(genes) and (genes[j + 1][2] - genes[j][3]) <= window_kb * 1000:
                j += 1
            n = j - i + 1
            if n >= min_size:
                cluster_counter[scaff] += 1
                cid = f"{scaff}:cluster_{cluster_counter[scaff]}"
                for k in range(i, j + 1):
                    out[genes[k][0]] = (n, cid)
                i = j + 1
            else:
                i += 1
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--gff", required=True, help="Path to genome GFF3")
    ap.add_argument("--candidates", required=True,
                    help="FASTA or text file of candidate gene IDs")
    ap.add_argument("--window-kb", type=int, default=100,
                    help="Max gap between consecutive cluster members in kb (default 100)")
    ap.add_argument("--min-size", type=int, default=3,
                    help="Minimum cluster size (default 3)")
    ap.add_argument("--db-cache", default=None,
                    help="Path for gffutils SQLite cache (default: <gff>.db)")
    ap.add_argument("--out", required=True, help="Output CSV path")
    args = ap.parse_args()

    cand_ids = load_candidate_ids(args.candidates)
    if not cand_ids:
        print(f"WARN: no candidate IDs found in {args.candidates}", file=sys.stderr)

    db = build_or_load_gff_db(args.gff, args.db_cache)
    candidate_genes = [g for g in iter_genes(db) if g[0] in cand_ids]
    matched = {g[0] for g in candidate_genes}
    missing = cand_ids - matched
    if missing:
        print(f"WARN: {len(missing)} candidate IDs not found in GFF "
              f"(first 5: {list(missing)[:5]})", file=sys.stderr)

    cluster_info = find_tandem_clusters(
        candidate_genes, window_kb=args.window_kb, min_size=args.min_size,
    )

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    rows = [
        {"candidate_id": gid, "tandem_cluster_size": size, "tandem_cluster_id": cid or ""}
        for gid, (size, cid) in cluster_info.items()
    ]
    # Also include candidates that were missing from the GFF (size=NaN, cid=NaN)
    for missing_id in missing:
        rows.append({"candidate_id": missing_id, "tandem_cluster_size": None,
                     "tandem_cluster_id": ""})
    df = pd.DataFrame(rows).sort_values("candidate_id")
    df.to_csv(args.out, index=False)

    n_in_clusters = sum(1 for v in cluster_info.values() if v[0] >= args.min_size)
    print(f"Wrote {len(df)} rows to {args.out}; {n_in_clusters} candidates "
          f"in tandem clusters (max size: "
          f"{max((v[0] for v in cluster_info.values()), default=0)})",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
