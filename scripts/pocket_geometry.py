#!/usr/bin/env python3
# pocket_geometry.py
# Purpose: Receptor-intrinsic small-molecule-vs-peptide discriminant for
#   class-A GPCR candidates, following Foster et al. 2019 (Cell 177:1933,
#   "Discovery of Human Signaling Systems: Pairing Peptides to G
#   Protein-Coupled Receptors"): peptide-binding class-A GPCRs have a LONG
#   extracellular loop 2 (ECL2, > ~20 residues) and a LARGE orthosteric
#   cavity, whereas small-molecule / odorant receptors (the targets of this
#   pipeline) have a SHORT ECL2 and a small, enclosed pocket.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst
#
# Scope / status: DESCRIPTIVE ONLY. This extractor emits per-candidate
#   columns (ecl2_length, pocket_volume, sm_vs_peptide_score, method_used)
#   intended as a 4th source for the non-chemoreceptor classifier (to flag
#   neuropeptide-receptor contaminants). It MUST NOT influence ranking until
#   its thresholds are calibrated against known controls (see
#   small_molecule_vs_peptide_score) and the change is consult-gated.
#
# Inputs:
#   * TMbed per-residue topology (3-line format: >id / sequence / labels,
#     labels B/b/H/h = TM strand/helix, S = signal, i/o = inside/outside
#     loops). ECL2 is parsed positionally as the loop between the 4th and
#     5th TM helices.
#   * AlphaFold3 structure (mmCIF or PDB). Pocket size is estimated with
#     fpocket when available, else a documented convex-hull geometric proxy.
#
# Usage:
#   python pocket_geometry.py --topology <file_or_dir> \
#                             --structures <file_or_dir> --output <tsv>

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# TMbed helix labels: 'H' (IN->OUT) and 'h' (OUT->IN). Both are TM alpha
# helices; orientation is not needed here because ECL2 is located
# positionally (the loop between the 4th and 5th helix), not by label.
HELIX_LABELS = "Hh"

# Structure file extensions recognized for the two Bio.PDB parsers.
CIF_SUFFIXES = {".cif", ".mmcif"}
PDB_SUFFIXES = {".pdb", ".ent"}
STRUCTURE_SUFFIXES = CIF_SUFFIXES | PDB_SUFFIXES

# ---------------------------------------------------------------------------
# PROVISIONAL thresholds for small_molecule_vs_peptide_score.
# These are placeholders grounded in Foster et al. 2019's ECL2 observation
# (~20-residue split) and typical fpocket cavity volumes. They MUST be
# calibrated against a labeled control set (known small-molecule/odorant vs
# known peptide receptors) before the score is used for anything beyond
# description. In particular the POCKET_* volumes are on the fpocket Angstrom^3
# scale; the geometric proxy below is a DIFFERENT scale and needs its own
# calibration.
# ---------------------------------------------------------------------------
ECL2_SM_MAX_LEN = 10.0        # <= this ECL2 length reads as small-molecule-like
ECL2_PEPTIDE_MIN_LEN = 20.0   # >= this ECL2 length reads as peptide-like
POCKET_SM_MAX_VOL = 600.0     # fpocket A^3: <= this reads as small-molecule-like
POCKET_PEPTIDE_MIN_VOL = 1200.0  # fpocket A^3: >= this reads as peptide-like


# ===========================================================================
# ECL2 length from a TMbed topology string
# ===========================================================================
def parse_tm_segments(
    topology_labels: str, helix_chars: str = HELIX_LABELS
) -> List[Tuple[int, int]]:
    """Return the (start, end) inclusive spans of each contiguous helix run.

    A "TM segment" is a maximal contiguous run of helix labels (default
    'H'/'h'). Consecutive TM helices in a TMbed prediction are always
    separated by at least one non-helix loop residue, so each run
    corresponds to one transmembrane helix. Positions are 0-indexed.
    """
    if not topology_labels:
        return []
    helix = set(helix_chars)
    segments: List[Tuple[int, int]] = []
    start: Optional[int] = None
    for i, ch in enumerate(topology_labels):
        if ch in helix:
            if start is None:
                start = i
        else:
            if start is not None:
                segments.append((start, i - 1))
                start = None
    if start is not None:
        segments.append((start, len(topology_labels) - 1))
    return segments


def ecl2_length(
    topology_labels: Optional[str], helix_chars: str = HELIX_LABELS
) -> Optional[int]:
    """Length (residue count) of ECL2, the loop between TM4 and TM5.

    ECL2 is located positionally: parse the contiguous helix runs, then take
    the gap between the 4th and 5th of them. Returns None when the input is
    empty / None or when fewer than 5 TM helices are present (ECL2
    undefined). Extra helices beyond the 7th are ignored.
    """
    if not topology_labels:
        return None
    segments = parse_tm_segments(topology_labels.strip(), helix_chars)
    if len(segments) < 5:
        return None
    tm4_end = segments[3][1]
    tm5_start = segments[4][0]
    return tm5_start - tm4_end - 1


# ===========================================================================
# Pocket volume: fpocket if available, else a geometric proxy
# ===========================================================================
def parse_fpocket_info(text: str) -> Optional[float]:
    """Return the top-ranked (Pocket 1) volume from an fpocket _info.txt.

    fpocket ranks pockets by score, so "Pocket 1" is the best-scored cavity.
    Returns None if no pocket / volume line is found.
    """
    lines = text.splitlines()
    in_pocket1 = False
    for line in lines:
        stripped = line.strip()
        m = re.match(r"^Pocket\s+(\d+)\s*:", stripped)
        if m:
            in_pocket1 = (m.group(1) == "1")
            continue
        if in_pocket1 and stripped.lower().startswith("volume"):
            vm = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", stripped.split(":", 1)[-1])
            if vm:
                return float(vm.group(0))
    return None


def _read_ca_coords(structure_path: str) -> Optional[np.ndarray]:
    """Read CA coordinates (fall back to all atoms) as an (N,3) array.

    Dispatches on file suffix: .cif/.mmcif -> MMCIFParser, else PDBParser.
    Returns None if the file cannot be parsed into any atoms.
    """
    from Bio.PDB import MMCIFParser, PDBParser

    suffix = Path(structure_path).suffix.lower()
    parser = MMCIFParser(QUIET=True) if suffix in CIF_SUFFIXES else PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("model", structure_path)
    except Exception as exc:  # malformed structure file
        print(f"Warning: could not parse structure {structure_path}: {exc}",
              file=sys.stderr)
        return None

    ca_coords: List[np.ndarray] = []
    all_coords: List[np.ndarray] = []
    for atom in structure.get_atoms():
        all_coords.append(atom.get_coord())
        if atom.get_id() == "CA":
            ca_coords.append(atom.get_coord())

    coords = ca_coords if len(ca_coords) >= 4 else all_coords
    if not coords:
        return None
    return np.asarray(coords, dtype=float)


def _bbox_volume(coords: np.ndarray) -> float:
    """Axis-aligned bounding-box volume (fallback when a hull is degenerate)."""
    spans = coords.max(axis=0) - coords.min(axis=0)
    return float(np.prod(spans))


def _convex_hull_volume(coords: np.ndarray) -> Tuple[float, str]:
    """Convex-hull volume of `coords`, robust to degenerate inputs.

    Returns (volume, method_suffix). Tries an exact hull, then a joggled
    hull for near-coplanar inputs, then the bounding box as a last resort.
    """
    from scipy.spatial import ConvexHull, QhullError

    try:
        return float(ConvexHull(coords).volume), "hull_A3"
    except QhullError:
        try:
            return float(ConvexHull(coords, qhull_options="QJ").volume), "hull_qj_A3"
        except QhullError:
            return _bbox_volume(coords), "bbox_A3"


def _proxy_pocket_volume(coords: np.ndarray) -> Tuple[float, str]:
    """Geometric proxy for orthosteric cavity size (no fpocket available).

    Method ("ca_upper_third"): the receptor's long (membrane-normal) axis is
    the first principal axis of the CA cloud. Residues are projected onto it
    and the third with the largest projection -- i.e. one membrane face,
    where the class-A orthosteric pocket sits -- are taken; the convex-hull
    volume of that sub-cloud is the proxy. This is a RELATIVE, provisional
    measure on its own Angstrom^3 scale (not comparable to fpocket) and must
    be calibrated before use.

    Note: the principal-axis sign is fixed deterministically (largest-|.|
    component made positive) but which membrane face is "extracellular" is
    not known from coordinates alone; pair with TMbed topology at
    calibration time to orient it.
    """
    n = len(coords)
    if n < 4:
        return _bbox_volume(coords), "proxy:bbox_A3"

    center = coords.mean(axis=0)
    centered = coords - center
    # first principal axis via SVD
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    axis = vt[0]
    # deterministic sign fix
    if axis[int(np.argmax(np.abs(axis)))] < 0:
        axis = -axis

    proj = centered @ axis
    k = max(4, int(math.ceil(n / 3.0)))
    k = min(k, n)
    top_idx = np.argsort(proj)[-k:]
    sub = coords[top_idx]

    vol, suffix = _convex_hull_volume(sub)
    return vol, f"proxy:ca_upper_third_{suffix}"


def _run_fpocket(structure_path: str) -> Optional[float]:
    """Run fpocket on a copy of the structure and return Pocket 1's volume.

    fpocket writes results next to its input, so we copy into a temp dir to
    avoid polluting the source tree. Returns None on any failure.
    """
    fpocket = shutil.which("fpocket")
    if not fpocket:
        return None
    src = Path(structure_path)
    with tempfile.TemporaryDirectory() as tmp:
        local = Path(tmp) / src.name
        shutil.copyfile(src, local)
        try:
            subprocess.run(
                [fpocket, "-f", str(local)],
                check=True, capture_output=True, text=True,
            )
        except (subprocess.CalledProcessError, OSError) as exc:
            print(f"Warning: fpocket failed on {structure_path}: {exc}",
                  file=sys.stderr)
            return None
        out_dir = Path(tmp) / (local.stem + "_out")
        info = out_dir / (local.stem + "_info.txt")
        if not info.exists():
            return None
        return parse_fpocket_info(info.read_text())


def pocket_volume(structure_path: str) -> Tuple[Optional[float], str]:
    """Estimate orthosteric pocket size from an AF3 structure.

    Returns (volume, method_used). method_used is:
      * "fpocket"                       -- fpocket was on PATH and succeeded.
      * "proxy:ca_upper_third_*"        -- documented geometric fallback.
      * "none"                          -- structure missing / unparseable.
    Handles both mmCIF (.cif/.mmcif) and PDB (.pdb/.ent).
    """
    if not structure_path or not os.path.isfile(structure_path):
        return None, "none"

    if shutil.which("fpocket"):
        vol = _run_fpocket(structure_path)
        if vol is not None:
            return vol, "fpocket"
        # fpocket present but produced nothing usable -> fall through to proxy

    coords = _read_ca_coords(structure_path)
    if coords is None:
        return None, "none"
    vol, method = _proxy_pocket_volume(coords)
    return vol, method


# ===========================================================================
# small-molecule-vs-peptide score
# ===========================================================================
def _linear_ramp(x: float, sm_bound: float, pep_bound: float) -> float:
    """1.0 at x<=sm_bound, 0.0 at x>=pep_bound, linear between (sm<pep)."""
    if x <= sm_bound:
        return 1.0
    if x >= pep_bound:
        return 0.0
    return (pep_bound - x) / (pep_bound - sm_bound)


def small_molecule_vs_peptide_score(
    ecl2_len: Optional[float], pocket_vol: Optional[float]
) -> Optional[float]:
    """Provisional discriminant in [0, 1]; HIGHER = more small-molecule-like.

    Short ECL2 + small pocket -> high (odorant / small-molecule-like);
    long ECL2 + large pocket -> low (peptide-like). Computed as the mean of
    the available sub-scores (ECL2 ramp, pocket ramp). Returns None only when
    BOTH inputs are missing.

    THRESHOLDS ARE PROVISIONAL (see ECL2_* / POCKET_* module constants) and
    must be calibrated against known-ligand-class controls before this score
    is used to influence anything. The pocket ramp assumes fpocket-scale
    Angstrom^3 volumes; the geometric proxy needs separate calibration.
    """
    components: List[float] = []
    if ecl2_len is not None:
        components.append(_linear_ramp(float(ecl2_len), ECL2_SM_MAX_LEN, ECL2_PEPTIDE_MIN_LEN))
    if pocket_vol is not None:
        components.append(_linear_ramp(float(pocket_vol), POCKET_SM_MAX_VOL, POCKET_PEPTIDE_MIN_VOL))
    if not components:
        return None
    return float(sum(components) / len(components))


# ===========================================================================
# I/O + CLI
# ===========================================================================
def parse_topology_file(path: str) -> Dict[str, str]:
    """Parse a TMbed/DeepTMHMM 3-line file into {candidate_id: label_string}.

    Each record is: header line (>id [description]) / sequence / topology
    labels. The candidate id is the first whitespace-delimited token of the
    header.
    """
    result: Dict[str, str] = {}
    with open(path) as fh:
        lines = fh.readlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith(">"):
            cid = lines[i][1:].strip().split()[0]
            if i + 2 < len(lines):
                result[cid] = lines[i + 2].strip()
            i += 3
        else:
            i += 1
    return result


def _resolve_structures(path: str, candidate_ids: List[str]) -> Dict[str, str]:
    """Map candidate_id -> structure path.

    A directory is scanned for structure files keyed by filename stem. A
    single file is bound to the sole candidate when there is exactly one,
    else keyed by its stem.
    """
    p = Path(path)
    mapping: Dict[str, str] = {}
    if p.is_dir():
        for f in sorted(p.iterdir()):
            if f.suffix.lower() in STRUCTURE_SUFFIXES:
                mapping[f.stem] = str(f)
    elif p.is_file():
        if len(candidate_ids) == 1:
            mapping[candidate_ids[0]] = str(p)
        else:
            mapping[p.stem] = str(p)
    return mapping


def _load_topologies(path: str) -> Dict[str, str]:
    """Load topologies from a single 3-line file or a directory of them."""
    p = Path(path)
    topos: Dict[str, str] = {}
    if p.is_dir():
        patterns = ("*.3line", "*.3lines", "*.tmbed", "*.txt")
        files: List[Path] = []
        for pat in patterns:
            files.extend(sorted(p.glob(pat)))
        for f in files:
            topos.update(parse_topology_file(str(f)))
    elif p.is_file():
        topos.update(parse_topology_file(str(p)))
    return topos


def _fmt(value: Optional[float]) -> str:
    """Format a possibly-None numeric for TSV (None -> blank)."""
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Extract ECL2 length + orthosteric pocket size and a provisional "
            "small-molecule-vs-peptide score per class-A GPCR candidate "
            "(Foster et al. 2019). DESCRIPTIVE only -- do not rank on it "
            "until calibrated."
        )
    )
    parser.add_argument(
        "--topology", required=True,
        help="TMbed 3-line topology file, or a directory of them.",
    )
    parser.add_argument(
        "--structures", required=True,
        help="AF3 structure file (mmCIF/PDB), or a directory of them.",
    )
    parser.add_argument("--output", required=True, help="Output TSV path.")
    args = parser.parse_args(argv)

    topologies = _load_topologies(args.topology)
    if not topologies:
        print(f"Warning: no topology records found in {args.topology}",
              file=sys.stderr)

    candidate_ids = sorted(topologies)
    structures = _resolve_structures(args.structures, candidate_ids)

    out_path = Path(args.output)
    if out_path.parent != Path(""):
        out_path.parent.mkdir(parents=True, exist_ok=True)

    columns = [
        "candidate_id", "ecl2_length", "pocket_volume",
        "sm_vs_peptide_score", "method_used",
    ]
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(columns)
        for cid in candidate_ids:
            ecl2 = ecl2_length(topologies[cid])
            struct_path = structures.get(cid)
            if struct_path:
                pvol, method = pocket_volume(struct_path)
            else:
                pvol, method = None, "none"
            score = small_molecule_vs_peptide_score(ecl2, pvol)
            writer.writerow([
                cid, _fmt(ecl2), _fmt(pvol), _fmt(score), method,
            ])

    print(f"Wrote {len(candidate_ids)} candidate rows to {out_path}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
