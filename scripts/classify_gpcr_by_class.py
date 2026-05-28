#!/usr/bin/env python3
"""classify_gpcr_by_class.py — Sequence-level GPCR class classifier.

P1 of the per-class refactor.

Input : protein FASTA
Output: per-sequence TSV
    seq_id  class  evidence_pfam  evidence_family_hmm  top_evalue

Class call hierarchy:
    1. hmmscan against Pfam 7TM-class HMMs (17 accessions covering A/B/C/F).
       Best E-value below --evalue wins → class label.
    2. hmmscan against 06c family HMMs (aminergic, opsin, etc.) independently.
       Best passing hit → evidence_family_hmm (subfamily annotation).
       This NEVER overrides the Pfam-derived class label.
    3. Special case: PF10324 (insect 7tm_6) → class=A but
       evidence_family_hmm=insect_OR_atypical even without a 06c hit.

TIAMMAT preference: on startup, glob family_hmm_dir for tiammat_*.hmm.
If found, prefer those over pfam_fallback/. Today this is a no-op (TIAMMAT
HMMs don't exist on disk yet — molluscagenes Stage 9 hasn't been run).

Bootstrap: if --bootstrap-pfam (default ON), download any missing Pfam HMMs
from InterPro + hmmpress them on first run. Idempotent.

Usage:
    python3 classify_gpcr_by_class.py \\
        --input candidates.fa \\
        --out candidate_classes.tsv \\
        [--pfam-dir results/classification/hmms/pfam_fallback] \\
        [--family-hmm-dir results/classification/hmms] \\
        [--tiammat-glob "tiammat_*.hmm"] \\
        [--evalue 1e-5] \\
        [--threads 4] \\
        [--no-bootstrap] \\
        [--force]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path
from typing import Callable, Optional

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# 17 Pfam accessions -> GPCR class (A/B/C/F).
# Class A covers rhodopsin-like (PF00001) plus all the expanded sr / 7tm_N
# families used in invertebrate chemoreceptor research.
_PFAM_TO_CLASS: dict[str, str] = {
    "PF00001": "A",   # 7tm_1      rhodopsin-like (canonical Class A)
    "PF00002": "B",   # 7tm_2      secretin/adhesion (Class B)
    "PF00003": "C",   # 7tm_3      glutamate-like (Class C)
    "PF01534": "F",   # Frizzled   frizzled/smoothened (Class F)
    "PF02949": "A",   # 7tm_4      vertebrate olfactory receptors
    "PF05296": "A",   # TAS2R      taste type-2 receptors
    "PF10324": "A",   # 7tm_6      insect olfactory receptors (inverted topology)
    "PF08395": "A",   # 7tm_5      C. elegans sr family
    "PF03402": "A",   # 7tm_7      C. elegans sr family
    "PF12022": "A",   # 7tm_8      C. elegans
    "PF11399": "A",   # 7tm_9      C. elegans
    "PF13853": "A",   # Srx        C. elegans nematode chemoreceptors
    "PF10326": "A",   # Srsx       C. elegans
    "PF13863": "A",   # Srbc       C. elegans
    "PF13886": "A",   # Srt        C. elegans
    "PF13887": "A",   # Sri        C. elegans
    "PF13889": "A",   # Srh        C. elegans
}

# Pfam accessions that carry an intrinsic subfamily annotation even
# without a 06c family-HMM hit.  Currently only PF10324 (atypical inverted
# topology — must NOT be grouped with mollusc Class A in per-class trees).
_PFAM_SUBFAMILY: dict[str, str] = {
    "PF10324": "insect_OR_atypical",
}

# 06c family-HMM prefixes that represent GPCR classes, not subfamilies.
# These are excluded from the refine_subfamily scan so they can't
# accidentally override the Pfam-derived class call via the subfamily field.
_CLASS_HMM_PREFIXES: tuple[str, ...] = (
    "class-B",
    "class-C",
    "class-F",
)

# InterPro HMM download endpoint (gzip-compressed .hmm)
_INTERPRO_HMM_URL = (
    "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{accession}?annotation=hmm"
)

# ---------------------------------------------------------------------------
# tblout parser (shared format with classify_via_hmm.py)
# ---------------------------------------------------------------------------

def parse_hmmscan_tblout(
    path: str,
) -> dict[str, list[tuple[str, str, float]]]:
    """Parse HMMER --tblout format.

    Returns dict: query_id -> list of (target_name, accession, evalue)
    sorted ascending by E-value (best hit first).

    Column layout (0-indexed):
        0  target name
        1  target accession
        2  query name
        3  query accession
        4  E-value (sequence, full)
        5  score
        ...
    """
    if not os.path.exists(path):
        return {}
    by_query: dict[str, list[tuple[str, str, float]]] = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            target = parts[0]
            accession = parts[1]   # '-' when not present
            query = parts[2]
            try:
                evalue = float(parts[4])
            except ValueError:
                continue
            by_query.setdefault(query, []).append((target, accession, evalue))
    for q in by_query:
        by_query[q].sort(key=lambda t: t[2])
    return by_query


# ---------------------------------------------------------------------------
# Pfam class call
# ---------------------------------------------------------------------------

def call_class(
    hits: list[tuple[str, str, float]],
    pfam_to_class: dict[str, str],
    evalue_threshold: float,
) -> tuple[str, str, str | float]:
    """Determine GPCR class from a list of Pfam hits for one query.

    Parameters
    ----------
    hits:
        List of (target_name, accession, evalue) tuples sorted best first.
        target_name is the HMM NAME field (e.g. "PF00001" or "7tm_1" after
        hmmpress); accession is the HMM accession field ("-" when absent).
    pfam_to_class:
        Mapping from Pfam accession (or HMM name) to class letter.
    evalue_threshold:
        Only hits with evalue < threshold are considered.

    Returns
    -------
    (class_label, evidence_pfam, top_evalue)
        class_label ∈ {A, B, C, F, unclassified}
        evidence_pfam = Pfam accession of the best hit (or "")
        top_evalue = float evalue or "" when unclassified
    """
    for target, accession, evalue in hits:
        if evalue >= evalue_threshold:
            continue
        # Try target name first (may be "PF00001"), then accession
        cls = pfam_to_class.get(target) or pfam_to_class.get(accession)
        if cls is not None:
            # Use whichever key matched as the evidence label
            evidence = target if target in pfam_to_class else accession
            return (cls, evidence, evalue)
    return ("unclassified", "", "")


# ---------------------------------------------------------------------------
# 06c subfamily refinement
# ---------------------------------------------------------------------------

def refine_subfamily(
    hits_06c: list[tuple[str, str, float]],
    query_id: str,
    *,
    pfam_accession: str = "",
    evalue_threshold: float = 1e-5,
) -> str:
    """Determine subfamily annotation from 06c family HMM hits.

    The class-B/C/F HMMs are EXCLUDED — they don't refine subfamily.
    Returns the name of the best-passing hit, or the _PFAM_SUBFAMILY
    annotation when pfam_accession is supplied (covers the PF10324 case
    where there may be no 06c hit at all but the topology is distinctive).

    Parameters
    ----------
    hits_06c:
        List of (target_name, accession, evalue) tuples for this query.
    query_id:
        Only hits whose query field matches query_id are used.
        (The list may already be pre-filtered per-query by the caller.)
    pfam_accession:
        If set and in _PFAM_SUBFAMILY, that value wins unconditionally
        (before the 06c scan).
    evalue_threshold:
        Hits at or above this threshold are ignored.
    """
    # Pfam intrinsic subfamily (e.g. insect_OR_atypical from PF10324)
    if pfam_accession and pfam_accession in _PFAM_SUBFAMILY:
        return _PFAM_SUBFAMILY[pfam_accession]

    for target, _acc, evalue in hits_06c:
        if evalue >= evalue_threshold:
            continue
        # Skip class-level HMMs
        if any(target.startswith(p) for p in _CLASS_HMM_PREFIXES):
            continue
        return target

    return ""


# ---------------------------------------------------------------------------
# TIAMMAT / pfam_fallback library selection
# ---------------------------------------------------------------------------

def select_hmm_library(
    pfam_dir: str,
    family_hmm_dir: str,
    tiammat_glob: str = "tiammat_*.hmm",
) -> dict:
    """Determine which Pfam HMM library to use.

    Prefers TIAMMAT HMMs (glob ``tiammat_glob`` under ``family_hmm_dir``)
    when present; falls back to ``pfam_dir`` otherwise.

    Returns a dict with keys:
        mode:   "tiammat" | "pfam_fallback"
        paths:  list[Path]  — the HMM files selected
        dir:    Path        — parent directory of the selected set
    """
    tiammat_paths = sorted(Path(family_hmm_dir).glob(tiammat_glob))
    if tiammat_paths:
        print(
            f"[classify_gpcr_by_class] TIAMMAT HMMs found ({len(tiammat_paths)}); "
            "using TIAMMAT library.",
            file=sys.stderr,
        )
        return {"mode": "tiammat", "paths": tiammat_paths, "dir": Path(family_hmm_dir)}

    pfam_paths = sorted(Path(pfam_dir).glob("*.hmm")) if pfam_dir else []
    print(
        "[classify_gpcr_by_class] No TIAMMAT HMMs found; using pfam_fallback.",
        file=sys.stderr,
    )
    return {"mode": "pfam_fallback", "paths": pfam_paths, "dir": Path(pfam_dir)}


# ---------------------------------------------------------------------------
# Pfam bootstrap download
# ---------------------------------------------------------------------------

def _download_pfam_hmm(accession: str, out_path: str) -> None:
    """Download a single Pfam HMM from InterPro (gzip -> decompress)."""
    url = _INTERPRO_HMM_URL.format(accession=accession)
    print(f"[classify_gpcr_by_class] Downloading {accession} from {url} ...",
          file=sys.stderr)
    req = urllib.request.Request(url, headers={"Accept": "application/octet-stream"})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            raw = resp.read()
    except Exception as exc:
        raise RuntimeError(f"Failed to download {accession}: {exc}") from exc

    # Response may be gzipped (content-type application/x-gzip or gzip magic)
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    Path(out_path).write_bytes(raw)

    # hmmpress the downloaded HMM for fast hmmscan
    _hmmpress(out_path)


def _hmmpress(hmm_path: str) -> None:
    """Run hmmpress on a single HMM file (removes stale index files first)."""
    for ext in ("h3f", "h3i", "h3m", "h3p"):
        stale = f"{hmm_path}.{ext}"
        if os.path.exists(stale):
            os.remove(stale)
    result = subprocess.run(
        ["hmmpress", hmm_path],
        stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True,
    )
    if result.returncode != 0:
        print(
            f"[classify_gpcr_by_class] WARN: hmmpress failed for {hmm_path}: "
            f"{result.stderr[:300]}",
            file=sys.stderr,
        )


def bootstrap_pfams(
    pfam_dir: str,
    required_pfams: list[str],
    *,
    _download_fn: Callable[[str, str], None] = _download_pfam_hmm,
    force: bool = False,
) -> None:
    """Ensure all required Pfam HMMs are present under pfam_dir.

    Downloads (via _download_fn) any that are missing.  Idempotent:
    already-present files are skipped unless force=True.

    Parameters
    ----------
    pfam_dir:
        Target directory (will be created if absent).
    required_pfams:
        List of Pfam accessions to ensure are present (e.g. ["PF00001", ...]).
    _download_fn:
        Override for tests — defaults to _download_pfam_hmm.
    force:
        Re-download even when file already exists.
    """
    Path(pfam_dir).mkdir(parents=True, exist_ok=True)
    for accession in required_pfams:
        out_path = str(Path(pfam_dir) / f"{accession}.hmm")
        if not force and Path(out_path).exists() and Path(out_path).stat().st_size > 0:
            continue
        _download_fn(accession, out_path)


# ---------------------------------------------------------------------------
# hmmscan runner
# ---------------------------------------------------------------------------

def run_hmmscan(
    query_fa: str,
    hmm_db: str,
    evalue: float,
    threads: int,
    *,
    tblout_path: Optional[str] = None,
) -> dict[str, list[tuple[str, str, float]]]:
    """Run hmmscan and return per-query hit lists.

    Parameters
    ----------
    query_fa:   Input FASTA
    hmm_db:     hmmpress-indexed HMM database
    evalue:     Report threshold passed to hmmscan -E
    threads:    --cpu value
    tblout_path: Where to write the tblout. If None, uses a tempfile.

    Returns the parse_hmmscan_tblout dict.
    """
    _own_tmp = tblout_path is None
    if _own_tmp:
        tmp = tempfile.NamedTemporaryFile(suffix=".tbl", delete=False)
        tblout_path = tmp.name
        tmp.close()
    try:
        result = subprocess.run(
            [
                "hmmscan",
                "--cpu", str(threads),
                "--tblout", tblout_path,
                "-E", str(evalue),
                hmm_db,
                query_fa,
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"hmmscan failed (rc={result.returncode}): {result.stderr[:500]}"
            )
        return parse_hmmscan_tblout(tblout_path)
    finally:
        if _own_tmp and os.path.exists(tblout_path):
            os.remove(tblout_path)


# ---------------------------------------------------------------------------
# TSV output
# ---------------------------------------------------------------------------

_TSV_COLUMNS = ["seq_id", "class", "evidence_pfam", "evidence_family_hmm", "top_evalue"]


def write_classification_tsv(rows: list[dict], out_path: str) -> None:
    """Write classification rows to out_path as a tab-separated file.

    Column order is fixed: seq_id, class, evidence_pfam,
    evidence_family_hmm, top_evalue.
    """
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(_TSV_COLUMNS)
        for row in rows:
            w.writerow([row.get(c, "") for c in _TSV_COLUMNS])


# ---------------------------------------------------------------------------
# Main classification logic
# ---------------------------------------------------------------------------

def _read_seq_ids_from_fasta(fasta_path: str) -> list[str]:
    ids: list[str] = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                tok = line[1:].split()[0]
                if tok:
                    ids.append(tok)
    return ids


def _build_concatenated_db(
    hmm_paths: list[Path],
    out_db: str,
) -> None:
    """Concatenate HMM files and hmmpress into a single searchable database."""
    with open(out_db, "w") as out:
        for hmm_path in hmm_paths:
            out.write(hmm_path.read_text())
    _hmmpress(out_db)


def classify_fasta(
    input_fa: str,
    pfam_dir: str,
    family_hmm_dir: str,
    *,
    tiammat_glob: str = "tiammat_*.hmm",
    evalue: float = 1e-5,
    threads: int = 4,
    bootstrap_pfam: bool = True,
) -> list[dict]:
    """Classify every sequence in input_fa.

    Returns a list of dicts with keys matching _TSV_COLUMNS.
    """
    required_pfams = list(_PFAM_TO_CLASS.keys())

    # Bootstrap missing Pfam HMMs
    if bootstrap_pfam:
        bootstrap_pfams(pfam_dir, required_pfams)

    # Select Pfam library (TIAMMAT or pfam_fallback)
    lib_info = select_hmm_library(pfam_dir, family_hmm_dir, tiammat_glob)
    pfam_hmm_paths = lib_info["paths"]

    if not pfam_hmm_paths:
        print(
            "[classify_gpcr_by_class] WARN: no Pfam HMMs found; "
            "all sequences will be unclassified.",
            file=sys.stderr,
        )

    # Family (06c) HMMs — exclude class-B/C/F-prefixed ones from refinement
    family_hmm_paths = [
        p for p in sorted(Path(family_hmm_dir).glob("*.hmm"))
        if not any(p.stem.startswith(pfx) for pfx in _CLASS_HMM_PREFIXES)
        and not p.stem.startswith("tiammat_")  # tiammat handled separately
    ]

    seq_ids = _read_seq_ids_from_fasta(input_fa)

    with tempfile.TemporaryDirectory(prefix="cgc_") as work:
        # Build Pfam db and scan
        pfam_hits: dict[str, list] = {}
        if pfam_hmm_paths:
            pfam_db = os.path.join(work, "pfam.hmm")
            _build_concatenated_db(pfam_hmm_paths, pfam_db)
            pfam_hits = run_hmmscan(input_fa, pfam_db, evalue, threads)

        # Build 06c family db and scan
        family_hits: dict[str, list] = {}
        if family_hmm_paths:
            fam_db = os.path.join(work, "family.hmm")
            _build_concatenated_db(family_hmm_paths, fam_db)
            family_hits = run_hmmscan(input_fa, fam_db, evalue, threads)

    rows: list[dict] = []
    for seq_id in seq_ids:
        p_hits = pfam_hits.get(seq_id, [])
        f_hits = family_hits.get(seq_id, [])

        cls, evidence_pfam, top_evalue = call_class(p_hits, _PFAM_TO_CLASS, evalue)

        # Determine subfamily: Pfam intrinsic first, then 06c refinement
        evidence_family_hmm = refine_subfamily(
            f_hits,
            seq_id,
            pfam_accession=evidence_pfam,
            evalue_threshold=evalue,
        )

        top_evalue_str = (
            f"{top_evalue:.3e}" if isinstance(top_evalue, float) else ""
        )

        rows.append({
            "seq_id": seq_id,
            "class": cls,
            "evidence_pfam": evidence_pfam,
            "evidence_family_hmm": evidence_family_hmm,
            "top_evalue": top_evalue_str,
        })

    return rows


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    default_results = os.environ.get("RESULTS_DIR", "results")
    default_pfam_dir = os.path.join(
        default_results, "classification", "hmms", "pfam_fallback"
    )
    default_family_dir = os.path.join(default_results, "classification", "hmms")
    default_threads = int(os.environ.get("CPUS", "4"))

    ap = argparse.ArgumentParser(
        description=(
            "Classify GPCR protein sequences into classes A/B/C/F "
            "using Pfam 7TM HMMs (class call) + 06c family HMMs (subfamily)."
        )
    )
    ap.add_argument("--input", required=True, metavar="FA",
                    help="Input protein FASTA")
    ap.add_argument("--out", required=True, metavar="TSV",
                    help="Output TSV path")
    ap.add_argument("--pfam-dir", default=default_pfam_dir, metavar="DIR",
                    help=f"Pfam fallback HMM directory (default: {default_pfam_dir})")
    ap.add_argument("--family-hmm-dir", default=default_family_dir, metavar="DIR",
                    help=f"06c family HMM directory (default: {default_family_dir})")
    ap.add_argument("--tiammat-glob", default="tiammat_*.hmm", metavar="GLOB",
                    help="Glob pattern for TIAMMAT HMMs (default: tiammat_*.hmm)")
    ap.add_argument("--evalue", type=float, default=1e-5, metavar="E",
                    help="E-value threshold (default: 1e-5)")
    ap.add_argument("--threads", type=int, default=default_threads, metavar="N",
                    help=f"CPU threads for hmmscan (default: {default_threads})")

    # Bootstrap flags (mutually exclusive pair)
    boot_grp = ap.add_mutually_exclusive_group()
    boot_grp.add_argument(
        "--bootstrap-pfam",
        dest="bootstrap_pfam",
        action="store_true",
        default=True,
        help="Download missing Pfam HMMs from InterPro on first run (default ON)",
    )
    boot_grp.add_argument(
        "--no-bootstrap",
        dest="bootstrap_pfam",
        action="store_false",
        help="Disable automatic Pfam bootstrap download",
    )

    ap.add_argument("--force", action="store_true", default=False,
                    help="Re-run even if --out already exists")
    return ap


def main() -> int:
    ap = build_arg_parser()
    args = ap.parse_args()

    if not args.force and Path(args.out).exists():
        print(
            f"[classify_gpcr_by_class] Output exists: {args.out}  "
            "(use --force to re-run)",
            file=sys.stderr,
        )
        return 0

    if shutil.which("hmmscan") is None:
        print("[classify_gpcr_by_class] ERROR: hmmscan not on PATH", file=sys.stderr)
        return 2
    if shutil.which("hmmpress") is None:
        print("[classify_gpcr_by_class] ERROR: hmmpress not on PATH", file=sys.stderr)
        return 2

    rows = classify_fasta(
        input_fa=args.input,
        pfam_dir=args.pfam_dir,
        family_hmm_dir=args.family_hmm_dir,
        tiammat_glob=args.tiammat_glob,
        evalue=args.evalue,
        threads=args.threads,
        bootstrap_pfam=args.bootstrap_pfam,
    )

    write_classification_tsv(rows, args.out)
    n_classified = sum(1 for r in rows if r["class"] != "unclassified")
    print(
        f"[classify_gpcr_by_class] DONE: {n_classified}/{len(rows)} "
        f"sequences classified -> {args.out}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
