#!/usr/bin/env python3
"""composition_distance.py — per-candidate amino-acid composition confound (A2).

A2 (epic v4bs): the embedding-novelty axis partly rewards amino-acid-composition
oddity (compositional signal lives in early PLM layers; Tule et al.), a confound
separable from both length and phylogeny. This producer computes each
candidate's composition distance — the Euclidean distance between its 20-dim
amino-acid frequency (one-hot-averaged) vector and the characterized-reference
centroid composition — for use as a ``--residual-confound`` in
``fusion_consensus.py`` (mirroring the seq_len / tree_distance confounds).

Pure/numpy (no torch/GPU); the FASTA IO lives in ``main``. It is a DORMANT
confound source: it adds a residual column/voter, never a veto.
"""
from __future__ import annotations

from typing import Dict, Mapping, Optional

import numpy as np

# Standard 20 proteinogenic amino acids; non-standard chars (X, B, Z, U, gaps)
# are ignored so composition is over resolved residues only.
_AA = "ACDEFGHIKLMNPQRSTVWY"
_AA_INDEX = {a: i for i, a in enumerate(_AA)}


def aa_composition(seq: str) -> np.ndarray:
    """20-dim normalized amino-acid frequency vector (one-hot averaged).

    Non-standard residues are ignored; a sequence with no standard residues
    returns the zero vector. Length-normalized, so composition is scale-free.
    """
    counts = np.zeros(20, dtype=float)
    for ch in seq.upper():
        idx = _AA_INDEX.get(ch)
        if idx is not None:
            counts[idx] += 1.0
    total = counts.sum()
    return counts / total if total > 0 else counts


def composition_distance(
    candidate_seqs: Mapping[str, str],
    reference_seqs: Optional[Mapping[str, str]] = None,
) -> Dict[str, float]:
    """Euclidean distance of each candidate's composition to a centroid.

    The centroid is the mean composition of ``reference_seqs`` (the
    characterized set) when given, else the mean composition of the candidates
    themselves (self-centered). Returns ``{candidate_id: distance}``; a larger
    distance means a more compositionally atypical candidate.
    """
    cand_vecs = {cid: aa_composition(s) for cid, s in candidate_seqs.items()}
    if reference_seqs:
        ref_vecs = [aa_composition(s) for s in reference_seqs.values()]
    else:
        ref_vecs = list(cand_vecs.values())
    centroid = np.mean(ref_vecs, axis=0) if ref_vecs else np.zeros(20)
    return {cid: float(np.linalg.norm(v - centroid)) for cid, v in cand_vecs.items()}


def main(argv=None) -> None:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--candidate-fasta", required=True, help="candidate protein FASTA")
    ap.add_argument(
        "--reference-fasta", default="",
        help="characterized-reference FASTA for the centroid "
             "(default: self-centered on the candidates)",
    )
    ap.add_argument("--out", required=True, help="output TSV (id\\tcomposition_distance)")
    args = ap.parse_args(argv)

    from Bio import SeqIO  # lazy: keep pure functions Biopython-free

    def _read(path: str) -> Dict[str, str]:
        return {rec.id: str(rec.seq) for rec in SeqIO.parse(path, "fasta")}

    cands = _read(args.candidate_fasta)
    refs = _read(args.reference_fasta) if args.reference_fasta else None
    dist = composition_distance(cands, refs)
    with open(args.out, "w") as fh:
        fh.write("id\tcomposition_distance\n")
        for cid in sorted(dist):
            fh.write(f"{cid}\t{dist[cid]:.6g}\n")
    print(f"[composition_distance] wrote {len(dist)} candidates -> {args.out}")


if __name__ == "__main__":
    main()
