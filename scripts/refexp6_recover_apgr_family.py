#!/usr/bin/env python3
"""Recover the Aplysia californica Apgr chemoreceptor family by HMM search.

WHY THIS EXISTS
---------------
The chemoreceptor probe had three members: the Aplysia californica receptors
Cummins et al. deposited (C5H877, C5H675, C5H674 -- "Chemosensory receptor
A/B/C", PMID 19493360). As a falsification instrument three is underpowered.
At the exclusion rates this pipeline calls (1.0-3.0%), three probes have only a
3.0-8.7% chance of producing a single exclusion under the null, so surviving is
the expected outcome whether the machinery is right or wrong.

Cummins reported roughly 90 genes and deposited three. The rest are in the
genome and unreachable by name query -- the RefSeq records carry generic
"G-protein coupled receptor"-style product names, not "chemosensory receptor".
An HMM search is the only route to them.

THE SEARCH IS ANCHORED ON THE SEEDS, NOT ON AN INVENTED CUTOFF
--------------------------------------------------------------
Every threshold here is derived from the three seeds' own measured behaviour.
A hit threshold chosen independently of the family is not calibrated to it.

  1. WHAT THE SEEDS HIT. hmmscan of the three seeds against all 27,481 Pfam-A
     profiles (HMMER 3.4). Above Pfam's curated gathering thresholds they hit
     exactly three profiles -- PF00001 (7tm_1), PF10324 (7TM_GPCR_Srw) and
     PF10328 (7TM_GPCR_Srx). PF10328 was NOT named in the task; it was found
     by asking the seeds rather than assuming. Measured scores:

         profile            C5H877   C5H675   C5H674   Pfam GA
         PF10324 Srw          58.1     52.2     24.6      22.2
         PF00001 7tm_1        56.2     41.7     51.3      30.5
         PF10328 Srx          28.6       --     33.2      23.8

     All three seeds clear PF10324 and PF00001, so those two are the family
     boundary the seeds themselves certify. Srx is hit by only two of three and
     is not used as a criterion.

  2. HOW RELATED THE SEEDS ARE TO EACH OTHER. phmmer all-vs-all gives the six
     off-diagonal bitscores 95.5, 103.4, 114.6, 122.1, 231.1, 231.8. The
     MINIMUM, 95.5, is the bar: "at least as related to a seed as the two most
     distant seeds are to each other". The family defines its own yardstick.

  3. SRW IS NOT PURE. InterPro describes Srw as "a solo family amongst the
     superfamilies of chemoreceptors", and it holds BOTH the nematode
     chemoreceptor-class Srw genes AND the deorphanized peptide receptors
     (Drosophila SPR, MIP-R). A bare Srw-GA cut therefore recovers peptide
     receptors alongside chemoreceptors -- it returns 266 of the proteome's
     26,656 proteins, far above Cummins's ~90. Each survivor is required to
     score higher against a SEED than against any of the 460 characterized
     peptide anchors in the reference set. This is the step that separates the
     target family from the reference family it is entangled with.

WHAT THE RECOVERED MEMBERS ARE, AND ARE NOT
-------------------------------------------
They are homologs of three published chemosensory receptors, recovered from the
same species and the same family those three represent. They are NOT
individually published chemoreceptors. The original probe roster admitted only
entries with a primary publication proposing them as chemosensory receptors,
on the reasoning that a probe whose own chemoreceptor claim is weak cannot
falsify anything -- if the machinery excludes it, the honest reading is that it
was never a chemoreceptor.

That reasoning still holds, and it is why every member carries a `stratum`:

    published-deposited   the 3 seeds. GenBank-deposited, named in the paper.
    published-independent C6FGJ8, a different species and a different
                          publication. The only member not derived from the
                          seeds by homology.
    hmm-recovered         everything this module finds. Evidence is homology
                          to the seeds plus family membership, NOT an
                          individual publication.

Exclusion of a `published-deposited` member is a strong red flag. Exclusion of
an `hmm-recovered` member is weaker: it is consistent with the machinery being
wrong AND with that particular gene model not being a real chemoreceptor.
Never pool the strata into a single pass/fail number.

Usage:
    python3 scripts/refexp6_recover_apgr_family.py --dry-run
    python3 scripts/refexp6_recover_apgr_family.py
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import os
import shutil
import sys
import urllib.parse
import urllib.request
from pathlib import Path

# The assembly the proteome came from. Resolved via `datasets summary genome
# taxon "Aplysia californica" --reference` on 2026-07-21, never recalled:
# AplCal3.0, RefSeq GCF_000002075.1, tax_id 6500, NCBI Annotation Release 102,
# 19,405 protein-coding genes, 26,656 proteins in protein.faa.
ASSEMBLY = "GCF_000002075.1"
ASSEMBLY_NAME = "AplCal3.0"
APLYSIA_CALIFORNICA_TAXID = "6500"
PROTEOME_MD5 = "80f91fdbcb5b2eb66b9a023423a51667"

# Pfam gathering thresholds, read from the Pfam-A.hmm actually searched.
PFAM_GA = {"PF00001": 30.5, "PF10324": 22.2, "PF10328": 23.8}

# The measured seed-vs-seed phmmer bitscores (all six off-diagonal values).
# The MINIMUM is the calibration bar. Pinned so a drift in the search cannot
# silently move the threshold.
SEED_VS_SEED_BITSCORES = (95.5, 103.4, 114.6, 122.1, 231.1, 231.8)
CALIBRATION_BAR = min(SEED_VS_SEED_BITSCORES)

RECOVERED_COLUMNS = [
    "refseq_accession", "organism", "taxid", "protein_name", "sequence_length",
    "best_seed", "best_seed_bitscore", "best_peptide_bitscore",
    "srw_bitscore", "tm7_bitscore", "assembly", "added_utc",
]

NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def read_fasta_titles(path: Path) -> dict[str, str]:
    """accession -> full defline description (the RefSeq product name)."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                head = line[1:].rstrip()
                acc = head.split()[0]
                rest = head[len(acc):].strip()
                # strip the trailing "[Organism]" RefSeq appends
                if rest.endswith("]") and "[" in rest:
                    rest = rest[:rest.rindex("[")].strip()
                out[acc] = rest
    return out


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def parse_tblout(path: Path) -> dict[str, tuple[str, float]]:
    """target -> (best query, best bitscore) from an hmmer tblout."""
    best: dict[str, tuple[str, float]] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            tgt, query, score = f[0], f[2], float(f[5])
            if tgt not in best or score > best[tgt][1]:
                best[tgt] = (query, score)
    return best


def parse_profile_hits(path: Path, profile: str) -> dict[str, float]:
    """target -> bitscore for one profile, from an hmmsearch tblout."""
    out: dict[str, float] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if f[2] == profile:
                sc = float(f[5])
                if f[0] not in out or sc > out[f[0]]:
                    out[f[0]] = sc
    return out


def verify_accessions_ncbi(accessions: list[str]) -> dict[str, dict]:
    """Confirm every recovered accession really is an A. californica protein.

    Accessions are never trusted because they came out of a file we wrote --
    they are re-resolved against NCBI and the organism is checked. A recovered
    "chemoreceptor" that is actually a different species' protein would corrupt
    the probe silently.
    """
    out: dict[str, dict] = {}
    for i in range(0, len(accessions), 200):
        chunk = accessions[i:i + 200]
        params = urllib.parse.urlencode({
            "db": "protein", "id": ",".join(chunk), "retmode": "json"})
        with urllib.request.urlopen(f"{NCBI_ESUMMARY}?{params}", timeout=180) as r:
            payload = json.load(r)
        res = payload.get("result", {})
        for uid in res.get("uids", []):
            rec = res[uid]
            out[rec.get("accessionversion", "")] = {
                "organism": rec.get("organism", ""),
                "taxid": str(rec.get("taxid", "")),
                "title": rec.get("title", ""),
                "slen": rec.get("slen"),
            }
    return out


def select_members(res_dir: Path, probe_seqs: dict[str, str]) -> dict:
    """Apply the seed-calibrated criteria. Returns a result dict."""
    proteome = read_fasta(res_dir / "protein.faa")
    srw_ga = {l.strip() for l in open(res_dir / "srw_ga_ids.txt") if l.strip()}
    tm7 = parse_profile_hits(res_dir / "proteome_hits_ga.tbl", "7tm_1")
    srw = parse_profile_hits(res_dir / "proteome_hits_ga.tbl", "7TM_GPCR_Srw")
    seed_best = parse_tblout(res_dir / "seeds_vs_proteome.tbl")
    pep_best = parse_tblout(res_dir / "pep_vs_srw.tbl")

    steps = {"srw_ga": len(srw_ga)}
    s1 = {p for p in srw_ga if seed_best.get(p, ("", 0.0))[1] >= CALIBRATION_BAR}
    steps["ge_calibration_bar"] = len(s1)
    s2 = {p for p in s1 if p in tm7}
    steps["class_a_confirmed"] = len(s2)

    kept, peptide_like = set(), []
    for p in s2:
        s = seed_best.get(p, ("", 0.0))[1]
        q = pep_best.get(p, ("", 0.0))[1]
        if s > q:
            kept.add(p)
        else:
            peptide_like.append((p, s, q))
    steps["closer_to_seed_than_peptide"] = len(kept)

    # Sequence-level dedup. Two recovered RefSeq entries ARE seeds under
    # different accessions (NP_001191494.1 == C5H675, NP_001191495.1 ==
    # C5H674); an accession-only check would admit the same protein twice.
    by_seq = {s: pid for pid, s in probe_seqs.items() if s}
    dupes, seen, final = [], {}, []
    for p in sorted(kept):
        s = proteome[p]
        if s in by_seq:
            dupes.append((p, by_seq[s]))
            continue
        if s in seen:
            dupes.append((p, f"INTRASET:{seen[s]}"))
            continue
        seen[s] = p
        final.append(p)
    steps["after_sequence_dedup"] = len(final)

    return {
        "final": final, "dupes": dupes, "peptide_like": peptide_like,
        "steps": steps, "proteome": proteome, "seed_best": seed_best,
        "pep_best": pep_best, "srw": srw, "tm7": tm7,
    }


def snapshot(paths: list[Path], outdir: Path) -> Path:
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    snap = outdir / f".apgr_snapshot_{stamp}"
    snap.mkdir(parents=True, exist_ok=True)
    for p in paths:
        if p.exists():
            shutil.copy2(p, snap / p.name)
    return snap


def write_atomic(path: Path, rows: list[dict], columns: list[str]) -> None:
    if not rows:
        raise SystemExit(f"ERROR: refusing to write an empty table to {path}")
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)


def write_fasta_atomic(path: Path, records: list[tuple[str, str]]) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")
    os.replace(tmp, path)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--res-dir", required=True,
                    help="directory holding the HMMER outputs and protein.faa")
    ap.add_argument("--probe-fasta",
                    default="references/anchors/chemoreceptor_probe_set.fasta")
    ap.add_argument("--outdir", default="references/anchors")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    res_dir, outdir = Path(args.res_dir), Path(args.outdir)
    now = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    probe_seqs = read_fasta(Path(args.probe_fasta))

    r = select_members(res_dir, probe_seqs)
    final, steps = r["final"], r["steps"]

    print("=== SEED-CALIBRATED RECOVERY ===")
    print(f"  calibration bar (min seed-seed bitscore) {CALIBRATION_BAR}")
    for k, v in steps.items():
        print(f"  {k:<32}{v:>6}")
    print(f"\n  dropped as peptide-receptor-like: {len(r['peptide_like'])}")
    print(f"  dropped as byte-identical duplicates: {len(r['dupes'])}")
    for p, w in r["dupes"]:
        print(f"      {p} == {w}")

    if not final:
        raise SystemExit("ERROR: recovered nothing -- refusing to write")

    # --- verify every accession against NCBI ---------------------------
    print(f"\n=== VERIFYING {len(final)} ACCESSIONS AGAINST NCBI ===")
    ncbi = verify_accessions_ncbi(final)
    titles = read_fasta_titles(res_dir / "protein.faa")
    rows, seqs, errors = [], [], []
    for p in final:
        rec = ncbi.get(p)
        seq = r["proteome"][p]
        if rec is None:
            errors.append(f"{p}: NCBI returned no record")
            continue
        if rec["taxid"] != APLYSIA_CALIFORNICA_TAXID:
            errors.append(f"{p}: taxid {rec['taxid']} != "
                          f"{APLYSIA_CALIFORNICA_TAXID} ({rec['organism']}) -- "
                          "the accession does not resolve to A. californica")
            continue
        if rec["slen"] is not None and int(rec["slen"]) != len(seq):
            errors.append(f"{p}: NCBI length {rec['slen']} != proteome length "
                          f"{len(seq)}")
            continue
        rows.append({
            "refseq_accession": p, "organism": rec["organism"],
            "taxid": rec["taxid"],
            "protein_name": titles.get(p, "") or rec["title"],
            "sequence_length": str(len(seq)),
            "best_seed": r["seed_best"][p][0],
            "best_seed_bitscore": f"{r['seed_best'][p][1]:.1f}",
            "best_peptide_bitscore": f"{r['pep_best'].get(p, ('', 0.0))[1]:.1f}",
            "srw_bitscore": f"{r['srw'].get(p, 0.0):.1f}",
            "tm7_bitscore": f"{r['tm7'].get(p, 0.0):.1f}",
            "assembly": ASSEMBLY, "added_utc": now,
        })
        seqs.append((p, seq))

    for e in errors:
        print(f"  FAIL {e}", file=sys.stderr)
    if errors:
        raise SystemExit(f"refusing to write: {len(errors)} accession(s) failed")
    print(f"  all {len(rows)} resolve to Aplysia californica "
          f"(taxid {APLYSIA_CALIFORNICA_TAXID}) with matching lengths")

    if args.dry_run:
        print("\n[dry-run] nothing written")
        return 0

    tsv_path = outdir / "apgr_recovered_members.tsv"
    fa_path = outdir / "apgr_recovered_members.fasta"
    snap = snapshot([tsv_path, fa_path], outdir)
    print(f"\nsnapshot -> {snap}")
    write_atomic(tsv_path, rows, RECOVERED_COLUMNS)
    write_fasta_atomic(fa_path, seqs)

    # --- re-read and verify --------------------------------------------
    reread = read_tsv(tsv_path)
    refa = read_fasta(fa_path)
    failures = []
    if len(reread) != len(rows):
        failures.append(f"rows {len(reread)} != {len(rows)}")
    accs = [x["refseq_accession"] for x in reread]
    if len(set(accs)) != len(accs):
        failures.append("recovered accessions are not unique")
    if len({refa[a] for a in refa}) != len(refa):
        failures.append("recovered sequences are not unique")
    for x in reread:
        a = x["refseq_accession"]
        if a not in refa:
            failures.append(f"{a} has no sequence")
        elif len(refa[a]) != int(x["sequence_length"]):
            failures.append(f"{a} length {len(refa[a])} != {x['sequence_length']}")

    print("\n=== POST-WRITE VERIFICATION ===")
    print(f"  recovered rows      {len(reread)}")
    print(f"  accessions unique   {len(set(accs)) == len(accs)}")
    print(f"  sequences unique    {len({refa[a] for a in refa}) == len(refa)}")
    print(f"  failures            {len(failures)}")
    if failures:
        for f in failures:
            print(f"  * {f}", file=sys.stderr)
        return 1
    print("\nALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
