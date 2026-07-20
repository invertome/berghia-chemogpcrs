#!/usr/bin/env python3
"""molluscan_calibration_frame.py — build the sampling frame for the molluscan
background null, and emit the FASTAs to embed.

WHY THIS EXISTS
---------------
The embedding channel's "is this candidate inside a known non-chemoreceptor
family envelope?" test compares a Berghia candidate's distance to family F
against the distances F's OWN reference members have to F. The reference set is
overwhelmingly vertebrate (1094 anchors; ~61 molluscan; four of the seven
families have ZERO molluscan members). Every within-family distance is
therefore vertebrate-to-vertebrate while every query distance is
mollusc-to-vertebrate across ~550 My. That is a systematic, directional
inflation applied to every candidate regardless of true family membership.

This script builds the population needed to measure that inflation: a
background of REAL molluscan class-A GPCRs drawn from the 90-species phase-1a
proteome scan, so a Berghia candidate's distance can be expressed relative to
"a typical molluscan class-A receptor" instead of "a vertebrate family member".

WHAT IT READS
-------------
  <scan-dir>/<sample>.scan_record.tsv   per-sequence record from
        scan_proteome_for_chemoreceptors.sh. Rows are the HMM-GPCR-positive
        sequences ONLY (every row has gpcr_positive=1); `passed_gate` records
        the >=6TM + confidence outcome. So the file DOES retain the pre-gate
        set relative to the TM gate -- which is what deliverable 3 needs.
  <scan-dir>/<sample>.chemo_candidates.fa   sequences that passed the gate.
  <classify>/class_phase1a.tsv          GPCR class (A/B/C/F/unclassified) for
        the 40,026 post-gate sequences ONLY.
  <proteome-dir>/<sample>.faa           source proteome; the ONLY place the
        pre-gate-only (TM-failed) sequences can be recovered from.
  <taxonomy-tsv>                        sample -> phylum/class/clade flags,
        obtained programmatically from NCBI taxonomy (not typed by hand).

DEDUPLICATION (measured, not assumed)
-------------------------------------
Three sample pairs share a taxid because the same assembly is filed under two
binomial synonyms; their candidate FASTAs are byte-identical. Left in, each of
those species would carry double weight in the null. One sample per taxid is
kept (lexicographically first) and the drop is recorded in the frame.

SAMPLING
--------
The default is a self-weighting simple random sample from the pooled molluscan
class-A post-gate set: an unbiased draw of "a random molluscan class-A
receptor", with species composition preserved as it actually is rather than
flattened by a per-species cap. Per-species and per-taxon-class composition is
written out so the draw can be audited, and nested phylogenetic scopes
(Mollusca / Gastropoda / Heterobranchia) are tagged per row so the same
embedded sample supports the phylogenetic-distance confound test without a
second embedding run.

The pre-gate-only (TM-failed) sequences are taken IN FULL rather than sampled:
they are the smaller stratum and they carry the whole of deliverable 3's
enrichment measurement, so subsampling them would be the one place a sampling
error would matter most.

NOTHING HERE IS A PRODUCTION RANKING. Outputs are calibration/validation
artifacts under results/ranking/diagnostics/.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import random
import sys
from collections import Counter, defaultdict


def read_fasta(path):
    """{first header token: sequence}. Preserves insertion order."""
    seqs, sid, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if sid is not None:
                    seqs[sid] = "".join(chunks).rstrip("*")
                sid, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if sid is not None:
        seqs[sid] = "".join(chunks).rstrip("*")
    return seqs


def write_fasta(path, records, width=60):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path + ".tmp", "w") as fh:
        for sid, seq in records:
            fh.write(f">{sid}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")
    os.replace(path + ".tmp", path)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scan-dir", required=True)
    p.add_argument("--class-tsv", required=True,
                   help="classify/class_phase1a.tsv (post-gate sequences only)")
    p.add_argument("--proteome-dir", required=True,
                   help="dir of <sample>.faa source proteomes; needed to "
                        "recover the pre-gate-only (TM-failed) sequences")
    p.add_argument("--taxonomy-tsv", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--n-postgate", type=int, default=6000,
                   help="size of the post-gate class-A molluscan sample")
    p.add_argument("--seed", type=int, default=20260720)
    a = p.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    rng = random.Random(a.seed)

    # ---- taxonomy -------------------------------------------------------
    tax = {}
    with open(a.taxonomy_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tax[r["sample"]] = r
    print(f"[frame] taxonomy rows: {len(tax)}")

    # ---- dedup samples that share a taxid --------------------------------
    by_taxid = defaultdict(list)
    for s, r in tax.items():
        by_taxid[r["taxid"]].append(s)
    keep, dropped = set(), []
    for tid, ss in by_taxid.items():
        ss = sorted(ss)
        keep.add(ss[0])
        for s in ss[1:]:
            dropped.append((s, ss[0], tid))
    if dropped:
        print(f"[frame] dropped {len(dropped)} duplicate-taxid samples:")
        for s, kept, tid in sorted(dropped):
            print(f"        {s}  (taxid {tid}, kept {kept})")
    print(f"[frame] samples kept: {len(keep)} of {len(tax)}")

    # ---- class assignment for the post-gate set --------------------------
    gpcr_class = {}
    with open(a.class_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gpcr_class[r["seq_id"]] = r["class"]
    print(f"[frame] class calls loaded: {len(gpcr_class)}")

    # ---- walk the scan records ------------------------------------------
    frame_rows = []
    postgate_seq = {}      # sid -> seq  (from chemo_candidates.fa)
    pregate_ids = defaultdict(list)   # sample -> [sid] needing proteome lookup
    n_dupskip = 0

    for sample in sorted(tax):
        rec = os.path.join(a.scan_dir, f"{sample}.scan_record.tsv")
        fa = os.path.join(a.scan_dir, f"{sample}.chemo_candidates.fa")
        if not os.path.exists(rec):
            print(f"[frame] FATAL: missing scan record for {sample}", file=sys.stderr)
            return 1
        if sample not in keep:
            n_dupskip += 1
            continue
        t = tax[sample]
        seqs = read_fasta(fa) if os.path.exists(fa) else {}
        with open(rec) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                sid = r["seq_id"]
                passed = r["passed_gate"] == "1"
                row = {
                    "seq_id": sid,
                    "sample": sample,
                    "taxid": t["taxid"],
                    "phylum": t["phylum"],
                    "tax_class": t["tax_class"],
                    "is_mollusca": t["is_mollusca"],
                    "is_gastropoda": t["is_gastropoda"],
                    "is_heterobranchia": t["is_heterobranchia"],
                    "passed_tm_gate": int(passed),
                    "tm_count": r["tm_count"],
                    "gpcr_class": gpcr_class.get(sid, "") if passed else "",
                    "seq_len": "",
                }
                if passed:
                    s = seqs.get(sid)
                    if s is None:
                        print(f"[frame] FATAL: {sid} passed the gate but is "
                              f"absent from {fa}", file=sys.stderr)
                        return 1
                    row["seq_len"] = len(s)
                    postgate_seq[sid] = s
                elif t["is_mollusca"] == "1":
                    # only molluscan pre-gate sequences feed the null, so only
                    # those proteomes are worth parsing
                    pregate_ids[sample].append(sid)
                frame_rows.append(row)

    print(f"[frame] skipped {n_dupskip} duplicate samples; "
          f"frame rows: {len(frame_rows)}")

    # ---- recover pre-gate-only sequences from the source proteomes -------
    pregate_seq = {}
    n_missing_proteome = 0
    for sample, ids in sorted(pregate_ids.items()):
        pfa = os.path.join(a.proteome_dir, f"{sample}.faa")
        if not os.path.exists(pfa):
            n_missing_proteome += 1
            print(f"[frame] WARN: no proteome for {sample}; "
                  f"{len(ids)} pre-gate-only sequences unrecoverable")
            continue
        pool = read_fasta(pfa)
        miss = [i for i in ids if i not in pool]
        if miss:
            print(f"[frame] WARN: {len(miss)}/{len(ids)} pre-gate ids absent "
                  f"from {sample}.faa (e.g. {miss[:2]})")
        for i in ids:
            if i in pool:
                pregate_seq[i] = pool[i]
    print(f"[frame] pre-gate-only sequences recovered: {len(pregate_seq)}; "
          f"{n_missing_proteome} samples had no proteome")

    for row in frame_rows:
        if not row["seq_len"] and row["seq_id"] in pregate_seq:
            row["seq_len"] = len(pregate_seq[row["seq_id"]])

    # ---- write the frame -------------------------------------------------
    frame_path = os.path.join(a.out_dir, "molluscan_calibration_frame.tsv")
    with open(frame_path + ".tmp", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(frame_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(frame_rows)
    os.replace(frame_path + ".tmp", frame_path)
    print(f"[frame] wrote {frame_path} ({len(frame_rows)} rows)")

    # ---- universes -------------------------------------------------------
    moll = [r for r in frame_rows if r["is_mollusca"] == "1"]
    post_A = [r for r in moll if r["passed_tm_gate"] == 1 and r["gpcr_class"] == "A"]
    pre_only = [r for r in moll if r["passed_tm_gate"] == 0
                and r["seq_id"] in pregate_seq]
    print(f"\n[frame] molluscan rows                       : {len(moll)}")
    print(f"[frame] molluscan post-gate                  : "
          f"{sum(1 for r in moll if r['passed_tm_gate'] == 1)}")
    print(f"[frame] molluscan post-gate class A          : {len(post_A)}")
    print(f"[frame] molluscan pre-gate-only (TM-failed)  : {len(pre_only)}")
    print(f"[frame]   of which gastropod                 : "
          f"{sum(1 for r in post_A if r['is_gastropoda'] == '1')} post-gate A")
    print(f"[frame]   of which heterobranch              : "
          f"{sum(1 for r in post_A if r['is_heterobranchia'] == '1')} post-gate A")

    # ---- draw the post-gate sample --------------------------------------
    n = min(a.n_postgate, len(post_A))
    sample_rows = rng.sample(post_A, n) if n < len(post_A) else list(post_A)
    sample_ids = {r["seq_id"] for r in sample_rows}
    print(f"\n[frame] post-gate class-A sample: {n} of {len(post_A)} "
          f"(seed {a.seed}, simple random, self-weighting)")
    comp = Counter(r["tax_class"] for r in sample_rows)
    universe = Counter(r["tax_class"] for r in post_A)
    for k in sorted(universe):
        print(f"        {k:<14} sample {comp[k]:>5}  universe {universe[k]:>6} "
              f"({100.0 * comp[k] / max(n, 1):5.1f}% vs "
              f"{100.0 * universe[k] / len(post_A):5.1f}%)")

    write_fasta(os.path.join(a.out_dir, "molluscan_null_postgate.fa"),
                [(r["seq_id"], postgate_seq[r["seq_id"]]) for r in sample_rows])
    write_fasta(os.path.join(a.out_dir, "molluscan_null_pregate_only.fa"),
                [(r["seq_id"], pregate_seq[r["seq_id"]]) for r in pre_only])

    with open(os.path.join(a.out_dir, "molluscan_calibration_sample.tsv"),
              "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["seq_id", "stratum", "sample", "taxid", "tax_class",
                    "is_gastropoda", "is_heterobranchia", "seq_len"])
        for r in sample_rows:
            w.writerow([r["seq_id"], "postgate_classA", r["sample"], r["taxid"],
                        r["tax_class"], r["is_gastropoda"],
                        r["is_heterobranchia"], r["seq_len"]])
        for r in pre_only:
            w.writerow([r["seq_id"], "pregate_only", r["sample"], r["taxid"],
                        r["tax_class"], r["is_gastropoda"],
                        r["is_heterobranchia"], r["seq_len"]])

    meta = {
        "seed": a.seed,
        "n_samples_total": len(tax),
        "n_samples_kept": len(keep),
        "duplicate_taxid_drops": [{"dropped": s, "kept": k, "taxid": t}
                                  for s, k, t in sorted(dropped)],
        "n_frame_rows": len(frame_rows),
        "n_molluscan_postgate_classA": len(post_A),
        "n_molluscan_pregate_only_recovered": len(pre_only),
        "n_postgate_sampled": n,
        "sample_taxon_class_composition": dict(comp),
        "universe_taxon_class_composition": dict(universe),
    }
    with open(os.path.join(a.out_dir, "molluscan_calibration_frame.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"\n[frame] wrote FASTAs + sample TSV + frame.json to {a.out_dir}")
    assert len(sample_ids) == n, "sampled ids are not unique"
    return 0


if __name__ == "__main__":
    sys.exit(main())
