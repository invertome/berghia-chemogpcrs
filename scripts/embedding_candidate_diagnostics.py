#!/usr/bin/env python3
"""Per-candidate embedding diagnostics + scoring-confound report (bead cw3.11).

For each candidate, in a given model's embedding space, scored the SAME way the
production channel scores (tied-covariance Mahalanobis vs CHARACTERIZED class-A
family prototypes — the class restriction is reused from build_embedding_channel,
not re-derived):
    novelty          min squared-Mahalanobis distance to the nearest family
    nearest_family   which characterized class-A family it most resembles
    margin           2nd-nearest-family distance minus nearest (assignment confidence)
plus two nuisance variables carried per candidate (seq_len, identity_to_nearest).

`confound_report` then measures whether novelty is really tracking a nuisance
variable rather than biology: Spearman(novelty, identity_to_nearest) is the
Dijkhof low-identity-generalization check (our orphan LSEs are low-identity to
every reference, so a strong negative correlation means the top-novelty ranks are
an IDENTITY ARTIFACT, not learned novelty). The two axes are read differently
(see CONFOUND_ROLE / confound_report): length is a pure artifact to minimize,
whereas identity is a redundancy/added-value axis — a MODERATE correlation is
expected (divergent sequences genuinely are more novel), and the decisive test is
whether novelty still discriminates WITHIN the low-identity tail rather than the
global correlation alone. The pathological case remains the ESM-C length confound
(rho -0.53), a pure artifact.
"""
from __future__ import annotations

import argparse
import sys
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, "scripts")
from embedding_evidence import family_prototypes, shrinkage_precision  # noqa: E402
from build_embedding_channel import (  # noqa: E402
    load_embeddings, load_ref_labels, novelty_reference_labels,
)

DIAG_COLUMNS = ["candidate_id", "novelty", "nearest_family", "margin",
                "seq_len", "identity_to_nearest"]
CONFOUNDS = ("seq_len", "identity_to_nearest")

# The two nuisance axes are NOT the same kind of thing:
#   seq_len            — a pure ARTIFACT (length has no legitimate bearing on
#                        biological novelty), so we want ~zero correlation.
#   identity_to_nearest— NOT a pure nuisance: it is a free, alignment-based proxy
#                        for the same divergence novelty is meant to capture, so a
#                        MODERATE correlation is expected and healthy. Here the
#                        check is REDUNDANCY / added-value: a near-perfect global
#                        correlation means the PLM just re-derives identity (we get
#                        it for free from mmseqs), and — decisively — novelty must
#                        still discriminate WITHIN the low-identity tail, where
#                        alignment is uninformative and the divergent orphans live.
CONFOUND_ROLE = {"seq_len": "artifact", "identity_to_nearest": "redundancy"}


def _tied_precision(emb: Dict[str, np.ndarray], labels: Dict[str, str]) -> np.ndarray:
    """Within-family-pooled, Ledoit-Wolf-shrunk precision (matches the producer)."""
    by_fam: Dict[str, list] = {}
    for k, fam in labels.items():
        if k in emb:
            by_fam.setdefault(fam, []).append(np.asarray(emb[k], float))
    centered = np.vstack([np.array(v) - np.array(v).mean(0) for v in by_fam.values()])
    return shrinkage_precision(centered)


def _per_family_dist(vec: np.ndarray, protos: Dict[str, np.ndarray],
                     precision: np.ndarray) -> Dict[str, float]:
    """Min squared-Mahalanobis distance from vec to each family's prototypes."""
    out = {}
    for fam, M in protos.items():
        M = np.atleast_2d(M)
        diff = vec[None, :] - M
        d = np.einsum("ni,ij,nj->n", diff, precision, diff)
        out[fam] = float(d.min())
    return out


def candidate_diagnostics(
    candidate_embeddings: Dict[str, np.ndarray],
    reference_embeddings: Dict[str, np.ndarray],
    reference_labels: Dict[str, str],
    seq_len: Dict[str, int],
    identity_to_nearest: Dict[str, float],
    k: int = 3,
) -> pd.DataFrame:
    labels = novelty_reference_labels(reference_labels)      # characterized class-A only
    protos = family_prototypes(reference_embeddings, labels, k)
    precision = _tied_precision(reference_embeddings, labels)

    rows: List[dict] = []
    for cid, vec in candidate_embeddings.items():
        vec = np.asarray(vec, float)
        d = _per_family_dist(vec, protos, precision)
        ordered = sorted(d.values())
        nearest_family = min(d, key=d.get)
        margin = (ordered[1] - ordered[0]) if len(ordered) > 1 else float("nan")
        rows.append({
            "candidate_id": cid,
            "novelty": ordered[0],
            "nearest_family": nearest_family,
            "margin": margin,
            "seq_len": seq_len.get(cid, np.nan),
            "identity_to_nearest": identity_to_nearest.get(cid, np.nan),
        })
    return pd.DataFrame(rows, columns=DIAG_COLUMNS)


def confound_report(df: pd.DataFrame, confounds=CONFOUNDS,
                    low_identity_quantile: float = 0.33,
                    artifact_rho: float = 0.4, redundant_rho: float = 0.85,
                    generalization_rho: float = 0.4) -> Dict[str, dict]:
    """Role-aware confound report (see CONFOUND_ROLE for why the axes differ).

    Every entry carries `spearman`/`p`/`n` (global Spearman(novelty, confound)),
    a `role`, and a `verdict`:

    - artifact axes (length): `verdict` = "artifact-confound" when |rho| >
      `artifact_rho`, else "clean". Here ~zero is the goal.
    - redundancy axes (identity): a MODERATE global rho is expected. The entry
      adds a `low_identity` sub-report — Spearman(novelty, identity) restricted to
      the bottom `low_identity_quantile` of identity (a scale-free cut, since the
      subset boundary is a quantile of whatever identity scale is supplied). The
      verdict is:
        "redundant-with-identity" if |global rho| > `redundant_rho` (the PLM just
            re-derives the free alignment identity → no added value);
        "low-identity-artifact"   if within the low-identity tail novelty still
            tracks residual identity (|low rho| > `generalization_rho`) → an
            identity proxy in the regime that matters (the Dijkhof failure);
        "adds-value"              otherwise (novelty is orthogonal where it counts).
    """
    rep: Dict[str, dict] = {}
    for c in confounds:
        if c not in df.columns:
            continue
        role = CONFOUND_ROLE.get(c, "artifact")
        sub = df[["novelty", c]].dropna()
        if len(sub) < 3:
            rep[c] = {"spearman": float("nan"), "p": float("nan"), "n": len(sub),
                      "role": role, "verdict": "insufficient-data"}
            continue
        rho, p = spearmanr(sub["novelty"], sub[c])
        entry = {"spearman": float(rho), "p": float(p), "n": int(len(sub)), "role": role}
        if role == "redundancy":
            thr = float(sub[c].quantile(low_identity_quantile))
            low = sub[sub[c] <= thr]
            if len(low) >= 3 and low[c].nunique() > 1:
                rho_low, p_low = spearmanr(low["novelty"], low[c])
            else:
                rho_low, p_low = float("nan"), float("nan")
            entry["low_identity"] = {"spearman": float(rho_low), "p": float(p_low),
                                     "n": int(len(low)), "threshold": thr,
                                     "quantile": low_identity_quantile}
            if abs(rho) > redundant_rho:
                entry["verdict"] = "redundant-with-identity"
            elif not np.isnan(rho_low) and abs(rho_low) > generalization_rho:
                entry["verdict"] = "low-identity-artifact"
            else:
                entry["verdict"] = "adds-value"
        else:
            entry["verdict"] = "artifact-confound" if abs(rho) > artifact_rho else "clean"
        rep[c] = entry
    return rep


def _read_lengths(fasta: str) -> Dict[str, int]:
    lens, sid, n = {}, None, 0
    for line in open(fasta):
        if line.startswith(">"):
            if sid is not None:
                lens[sid] = n
            sid = line[1:].split()[0]; n = 0
        else:
            n += len(line.strip())
    if sid is not None:
        lens[sid] = n
    return lens


def _read_identity(path: str) -> Dict[str, float]:
    """candidate_id \t max_pct_identity_to_any_reference (from an mmseqs/diamond
    search, computed on Unity). Missing candidates default to 0.0 identity."""
    if not path:
        return {}
    d = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            try:
                d[p[0]] = float(p[1])
            except ValueError:
                pass
    return d


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidate-npz", required=True)
    p.add_argument("--ref-npz", required=True)
    p.add_argument("--ref-labels", required=True)
    p.add_argument("--candidate-fasta", required=True, help="for seq_len")
    p.add_argument("--identity-tsv", default="", help="candidate_id<TAB>max_pct_id (mmseqs); optional")
    p.add_argument("--model", required=True, help="model tag, recorded in output")
    p.add_argument("--out", required=True)
    p.add_argument("-k", type=int, default=3)
    a = p.parse_args(argv)

    cand = load_embeddings(a.candidate_npz)
    ref = load_embeddings(a.ref_npz)
    labels = load_ref_labels(a.ref_labels)
    if not cand or not ref or not labels:
        sys.exit(f"[candidate_diagnostics] empty input(s): cand={len(cand)} ref={len(ref)} labels={len(labels)}")
    seq_len = _read_lengths(a.candidate_fasta)
    identity = _read_identity(a.identity_tsv)

    df = candidate_diagnostics(cand, ref, labels, seq_len, identity, k=a.k)
    df.insert(0, "model", a.model)
    import os
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    df.to_csv(a.out, sep="\t", index=False)

    rep = confound_report(df)
    print(f"[candidate_diagnostics] {a.model}: {len(df)} candidates -> {a.out}")
    for c, r in rep.items():
        line = (f"  [{r['role']}] novelty~{c}: Spearman={r['spearman']:+.3f} "
                f"(p={r['p']:.1e}, n={r['n']}) -> {r['verdict']}")
        li = r.get("low_identity")
        if li:
            line += (f"; low-identity(bottom {li['quantile']:.0%}) "
                     f"Spearman={li['spearman']:+.3f} (n={li['n']})")
        print(line)


if __name__ == "__main__":
    main()
