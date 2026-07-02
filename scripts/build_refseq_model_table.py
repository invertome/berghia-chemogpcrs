#!/usr/bin/env python3
"""build_refseq_model_table.py — RefSeq gene-model + loci table builder.

Extracted verbatim (Task-6 review follow-up) from the inline RefSeq GFF-walk
heredoc in ``02c_genome_reconcile.sh`` so its logic — especially the Task-6
``da7142e`` fail-loud-on-unresolvable never-drop guard — can be unit-tested.

Given the RefSeq chemoreceptor candidate ids (one versioned ``XP_<ver>`` per
line), the RefSeq GFF3, the candidate protein FASTA, and the candidate
DeepTMHMM prediction, it writes the two inputs ``reconcile_candidates.py``
consumes:

  --refseq-models  TSV ``query chrom start end strand complete length n_tm``
                   (``#``-commented header; the ``complete`` cell is the literal
                   ``complete``/``partial``). One row per candidate; coords come
                   from the candidate's *gene* feature. ``complete`` = NOT
                   partial, where partial = ``partial=true`` / ``start_range=`` /
                   ``end_range=`` on the gene, mRNA, OR CDS. ``length`` = residue
                   count from the protein FASTA; ``n_tm`` from the prediction.
  --refseq-loci    Headerless 4-column ``chrom start end strand`` of ALL gene
                   intervals in the GFF (not just candidates) — so a candidate
                   overlapping a NON-candidate RefSeq gene is still
                   chimeric-detectable.

Each candidate ``XP_<ver>`` is resolved to its gene via the CDS ``gene=LOC…``
attribute, or failing that the ``Parent=rna-… -> mRNA -> gene`` chain. The join
key is the VERSIONED accession (consistent across the ids file, the GFF
``protein_id=``, and the FASTA header's first token) — versions are NOT stripped.

FAIL LOUD: any candidate whose CDS cannot resolve to a gene is collected and,
if the list is non-empty, printed to stderr with a nonzero exit — a genome
candidate is never silently dropped (the ``da7142e`` contract, kept as strict
as ``reconcile_candidates.py``'s readers).

Stdlib-only, no network access.
"""
from __future__ import annotations

import argparse
import os
import re
import sys

MODELS_HEADER = "#query\tchrom\tstart\tend\tstrand\tcomplete\tlength\tn_tm\n"


def attr(attrs: str, key: str):
    """Return GFF3 attribute ``key``'s value in the col-9 string, or None.

    Anchored on ``^`` / ``;`` so ``gene=`` is not matched inside a longer key
    such as ``gene_biotype=`` (a common RefSeq decoy)."""
    m = re.search(r'(?:^|;)' + re.escape(key) + r'=([^;]+)', attrs)
    return m.group(1) if m else None


def is_partial(attrs: str) -> bool:
    """True if the attribute string marks the feature incomplete
    (``partial=true`` / ``start_range=`` / ``end_range=``)."""
    a = attrs.lower()
    return ("partial=true" in a) or ("start_range=" in a) or ("end_range=" in a)


def read_candidate_ids(ids_file: str):
    """Read the candidate id file (one versioned ``XP_<ver>`` per line)."""
    cand = set()
    with open(ids_file) as fh:
        for line in fh:
            s = line.strip()
            if s:
                cand.add(s)
    return cand


def parse_gff(gff: str, cand):
    """Walk the RefSeq GFF3 once, returning
    ``(genes_by_id, genes_by_loc, mrna_parent, all_loci, cds_for_prot)``:

    * ``genes_by_id``  gene feature id -> ``(chrom, start, end, strand, partial)``
    * ``genes_by_loc`` ``LOC…`` symbol -> the same tuple
    * ``mrna_parent``  rna id -> ``(gene feature id, mrna_partial)``
    * ``all_loci``     ``(chrom, start, end, strand)`` for EVERY gene (GFF order)
    * ``cds_for_prot`` candidate protein_id -> ``{loc, parent, partial}``
    """
    genes_by_id, genes_by_loc = {}, {}   # -> (chrom,start,end,strand,partial)
    mrna_parent = {}                     # rna id -> (gene feature id, mrna_partial)
    all_loci = []                        # (chrom,start,end,strand) for every gene
    cds_for_prot = {}                    # protein_id -> {loc, parent, partial}

    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, ftype, start, end, strand, a = f[0], f[2], f[3], f[4], f[6], f[8]
            if ftype == "gene":
                gid = attr(a, "ID")
                loc = attr(a, "gene") or (gid[len("gene-"):] if gid and gid.startswith("gene-") else None)
                tup = (chrom, int(start), int(end), strand, is_partial(a))
                if gid:
                    genes_by_id[gid] = tup
                if loc:
                    genes_by_loc[loc] = tup
                all_loci.append((chrom, int(start), int(end), strand))
            elif ftype in ("mRNA", "transcript"):
                rid = attr(a, "ID")
                if rid:
                    mrna_parent[rid] = (attr(a, "Parent"), is_partial(a))
            elif ftype == "CDS":
                pid = attr(a, "protein_id")
                if pid and pid in cand:
                    d = cds_for_prot.setdefault(pid, {"loc": None, "parent": None, "partial": False})
                    d["loc"] = d["loc"] or attr(a, "gene")
                    d["parent"] = d["parent"] or attr(a, "Parent")
                    d["partial"] = d["partial"] or is_partial(a)
    return genes_by_id, genes_by_loc, mrna_parent, all_loci, cds_for_prot


def read_protein_lengths(prot_fa: str):
    """Residue-count length per protein, keyed by the header's first token."""
    plen, name, ln = {}, None, 0
    if os.path.exists(prot_fa):
        with open(prot_fa) as fh:
            for line in fh:
                if line.startswith(">"):
                    if name is not None:
                        plen[name] = ln
                    name, ln = line[1:].split()[0], 0
                else:
                    ln += len(line.strip())
            if name is not None:
                plen[name] = ln
    return plen


def read_n_tm(pred: str):
    """TM-region count per id from a DeepTMHMM prediction
    (col1 = id, col5 = TM count). Short/unparseable rows are skipped."""
    ntm = {}
    if os.path.exists(pred):
        with open(pred) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 5:
                    continue
                try:
                    ntm[f[0]] = int(float(f[4]))
                except ValueError:
                    continue
    return ntm


def resolve(pid, genes_by_id, genes_by_loc, mrna_parent, cds_for_prot):
    """Resolve a candidate protein id to ``(chrom, start, end, strand, partial)``
    or ``None`` if it cannot reach a gene locus.

    Route A: the CDS ``gene=LOC…`` symbol resolves via ``genes_by_loc``.
    Route B: else the CDS ``Parent=rna-… -> mRNA -> gene`` chain via
    ``genes_by_id``. ``partial`` ORs the gene, mRNA, and CDS partial flags."""
    d = cds_for_prot.get(pid)
    if not d:
        return None
    if d["loc"] and d["loc"] in genes_by_loc:
        g, part = genes_by_loc[d["loc"]], d["partial"]
    else:
        gidf = None
        if d["parent"] and d["parent"] in mrna_parent:
            gidf = mrna_parent[d["parent"]][0]
        if not (gidf and gidf in genes_by_id):
            return None
        g, part = genes_by_id[gidf], d["partial"]
    mrna_part = mrna_parent[d["parent"]][1] if d["parent"] in mrna_parent else False
    partial = g[4] or mrna_part or part
    return g[0], g[1], g[2], g[3], partial


def build_tables(ids_file, gff, prot_fa, pred, models_out, loci_out):
    """Write the ``--refseq-models`` and ``--refseq-loci`` TSVs.

    Returns ``(n_written, all_loci, unresolved)``. Never silently drops a
    candidate: ids that fail to resolve to a gene are collected into
    ``unresolved`` for the caller to fail loud on (both output files are still
    written first, exactly as the original heredoc did)."""
    cand = read_candidate_ids(ids_file)
    genes_by_id, genes_by_loc, mrna_parent, all_loci, cds_for_prot = parse_gff(gff, cand)
    plen = read_protein_lengths(prot_fa)
    ntm = read_n_tm(pred)

    n_written = 0
    unresolved = []
    with open(models_out, "w") as mf:
        mf.write(MODELS_HEADER)
        for pid in sorted(cand):
            r = resolve(pid, genes_by_id, genes_by_loc, mrna_parent, cds_for_prot)
            if r is None:
                # Never silently drop a genome candidate (the module's never-drop
                # contract): collect and fail loud after the loop. Unreachable on
                # the RefSeq GFF (every protein-coding CDS carries gene= AND
                # Parent=), but this keeps the stage as strict as
                # reconcile_candidates.py's readers.
                unresolved.append(pid)
                continue
            chrom, gs, ge, strand, partial = r
            mf.write(f"{pid}\t{chrom}\t{gs}\t{ge}\t{strand}\t"
                     f"{'partial' if partial else 'complete'}\t{plen.get(pid, 0)}\t{ntm.get(pid, 0)}\n")
            n_written += 1

    with open(loci_out, "w") as lf:
        for chrom, s, e, strand in all_loci:
            lf.write(f"{chrom}\t{s}\t{e}\t{strand}\n")

    return n_written, all_loci, unresolved


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="build_refseq_model_table.py",
        description="Build the reconcile module's --refseq-models and "
                    "--refseq-loci TSV inputs by walking the Berghia RefSeq "
                    "GFF3 for each candidate XP_ accession.")
    p.add_argument("--ids", required=True,
                   help="RefSeq candidate ids file (one versioned XP_<ver> per line)")
    p.add_argument("--gff", required=True, help="RefSeq GFF3 (coords / loci)")
    p.add_argument("--proteins", required=True,
                   help="RefSeq candidate protein FASTA (residue-count length)")
    p.add_argument("--prediction", required=True,
                   help="RefSeq DeepTMHMM prediction (col1=id, col5=n_tm)")
    p.add_argument("--models-out", required=True,
                   help="output --refseq-models TSV")
    p.add_argument("--loci-out", required=True,
                   help="output --refseq-loci TSV")
    return p


def main(argv=None) -> int:
    """CLI: build both tables, then fail loud (nonzero) if any candidate did
    not resolve to a gene locus — never silently dropping a genome candidate."""
    args = _build_arg_parser().parse_args(argv)
    n_written, all_loci, unresolved = build_tables(
        args.ids, args.gff, args.proteins, args.prediction,
        args.models_out, args.loci_out)
    if unresolved:
        sys.stderr.write(
            f"ERROR: {len(unresolved)} RefSeq candidate(s) did not resolve to a gene "
            f"locus in {args.gff} (a genome candidate would otherwise vanish silently): "
            f"{', '.join(unresolved)}\n")
        return 1
    sys.stderr.write(f"refseq models: {n_written} written; refseq loci: {len(all_loci)}\n")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
