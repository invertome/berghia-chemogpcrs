"""The chemosensory tissue vocabulary must have exactly ONE definition.

config.sh defines it twice, and the two disagree:

    :513  CHEMOSENSORY_TISSUES="rhinophore,oral_veil,tentacle,cephalic"
    :743  CHEMOSENSORY_TISSUE_WEIGHTS="rhinophore:2.0,oral-tentacle:1.0"

`oral-tentacle` names nothing the producer emits. `tissue_tpms` keys come
from expression_summary.csv's ``<tissue>_tpm`` columns and are looked up by
EXACT key (`tissue_tpms.get(tissue, 0.0)`), so the weights map's second entry
is a phantom tissue. Two consequences in
get_chemosensory_expression_score() (rank_candidates.py:1112,1136), feeding
CHEMOSENSORY_EXPR_WEIGHT=3:

1. A candidate with 400 TPM in `tentacle` and 4 TPM in `rhinophore`
   contributes the 400 TPM as exactly ZERO, while the phantom tissue's
   weight still counts in the divisor -- so the mean is diluted by a tissue
   that can never contribute.

2. `other_chemo_tpms` is always `[0.0]`, so `rhino_tpm >= max_other_chemo` is
   ALWAYS True whenever rhino_tpm > 0. The 1.75x "rhinophore-dominant" bonus
   is therefore granted to every chemosensory-specific candidate regardless
   of its tentacle expression, and the 1.5x branch is unreachable.

The fix makes CHEMOSENSORY_TISSUES the single controlled vocabulary and
demotes CHEMOSENSORY_TISSUE_WEIGHTS to weights-only over that vocabulary: a
tissue in the vocabulary with no explicit weight defaults to 1.0, and a
weight naming a tissue OUTSIDE the vocabulary is a hard error rather than a
silent zero.

rank_candidates.py is not import-safe (top-level ``sys.argv``); these tests
exec single top-level function fragments, the pattern established in
tests/unit/test_rank_candidates_dnds_reliability.py.
"""
from __future__ import annotations

import ast
import re
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
CONFIG = PROJECT_ROOT / "config.sh"
RANK = PROJECT_ROOT / "scripts" / "rank_candidates.py"


def _exec_function(name: str) -> Any:
    """Compile just the named top-level function out of rank_candidates.py.

    Uses ast to take the function's exact line span rather than slicing to the
    next ``\\n\\ndef ``, so a module-level statement sitting between two defs
    can't be dragged into the fragment.
    """
    src = RANK.read_text()
    for node in ast.parse(src).body:
        if isinstance(node, ast.FunctionDef) and node.name == name:
            ns: dict[str, Any] = {"pd": pd, "np": np, "os": __import__("os")}
            exec(ast.unparse(node), ns)
            return ns[name]
    raise AssertionError(f"{name} not found as a top-level def in {RANK}")


def _config_value(name: str) -> str:
    """Read an exported scalar from config.sh without sourcing it."""
    m = re.search(rf'^export {name}="([^"]*)"', CONFIG.read_text(), re.MULTILINE)
    assert m, f"{name} not found in config.sh"
    return m.group(1)


# --------------------------------------------------------------------------
# the two definitions must agree
# --------------------------------------------------------------------------

def test_config_is_syntactically_valid() -> None:
    assert subprocess.run(["bash", "-n", str(CONFIG)]).returncode == 0


def test_weight_keys_are_a_subset_of_the_tissue_vocabulary() -> None:
    """Every weighted tissue must be a real member of CHEMOSENSORY_TISSUES."""
    vocab = {t.strip() for t in _config_value("CHEMOSENSORY_TISSUES").split(",")}
    weighted = {p.split(":")[0].strip()
                for p in _config_value("CHEMOSENSORY_TISSUE_WEIGHTS").split(",")
                if p.strip()}
    phantom = weighted - vocab
    assert not phantom, (
        f"CHEMOSENSORY_TISSUE_WEIGHTS names {sorted(phantom)}, which is not in "
        f"CHEMOSENSORY_TISSUES {sorted(vocab)}; the lookup is by exact key so "
        f"those weights match nothing the producer emits."
    )


# --------------------------------------------------------------------------
# resolution: vocabulary + weights -> effective weight map
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def resolve():
    return _exec_function("resolve_tissue_weights")


def test_unweighted_vocabulary_members_default_to_one(resolve) -> None:
    """A tissue in the vocabulary but absent from the weights still counts."""
    assert resolve("rhinophore,oral_veil,tentacle,cephalic",
                   "rhinophore:2.0") == {
        "rhinophore": 2.0, "oral_veil": 1.0,
        "tentacle": 1.0, "cephalic": 1.0}


def test_weight_outside_the_vocabulary_is_a_hard_error(resolve) -> None:
    """A phantom tissue must fail loudly, not silently contribute zero."""
    with pytest.raises(ValueError, match="oral-tentacle"):
        resolve("rhinophore,oral_veil,tentacle,cephalic",
                "rhinophore:2.0,oral-tentacle:1.0")


def test_explicit_weights_are_preserved(resolve) -> None:
    assert resolve("rhinophore,tentacle",
                   "rhinophore:2.0,tentacle:1.5") == {
        "rhinophore": 2.0, "tentacle": 1.5}


# --------------------------------------------------------------------------
# the two scoring consequences
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def score_fn():
    return _exec_function("get_chemosensory_expression_score")


def _score(score_fn, tissue_tpms: dict, weights: dict, specific: bool = False,
           tau: float = 0.0, enrichment: float = 0.0) -> float:
    data = {
        "cand": {
            "has_data": True,
            "tissue_tpms": tissue_tpms,
            "tau_index": tau,
            "chemosensory_enrichment": enrichment,
            "is_chemosensory_specific": specific,
        }
    }
    score, has = score_fn("cand", data, tissue_weights=weights)
    assert has is True
    return score


def test_tentacle_expression_actually_contributes(score_fn, resolve) -> None:
    """400 TPM tentacle / 4 TPM rhinophore must not score as rhinophore-only."""
    weights = resolve("rhinophore,oral_veil,tentacle,cephalic", "rhinophore:2.0")
    with_tentacle = _score(score_fn, {"rhinophore": 4.0, "tentacle": 400.0}, weights)
    rhino_only = _score(score_fn, {"rhinophore": 4.0}, weights)
    assert with_tentacle > rhino_only, (
        "the 400 TPM tentacle signal contributed nothing to the score"
    )


def test_rhinophore_dominance_bonus_is_not_automatic(score_fn, resolve) -> None:
    """The 1.5x branch must be reachable: tentacle-dominant gets the lower bonus.

    With the phantom-tissue bug, other_chemo_tpms is always [0.0], so
    rhino_tpm >= max_other_chemo always holds and every chemosensory-specific
    candidate collects the 1.75x rhinophore-dominant multiplier.
    """
    weights = resolve("rhinophore,oral_veil,tentacle,cephalic", "rhinophore:2.0")
    tpms = {"rhinophore": 4.0, "tentacle": 400.0}
    boosted = _score(score_fn, tpms, weights, specific=True, tau=0.9)
    plain = _score(score_fn, tpms, weights, specific=False, tau=0.9)
    ratio = boosted / plain
    assert ratio == pytest.approx(1.5, abs=1e-9), (
        f"tentacle-dominant candidate got a {ratio:.2f}x bonus; the 1.75x "
        f"rhinophore-dominant branch fired on a candidate whose rhinophore "
        f"TPM (4) is 100x BELOW its tentacle TPM (400)."
    )


def test_rhinophore_dominant_still_earns_the_higher_bonus(score_fn, resolve) -> None:
    weights = resolve("rhinophore,oral_veil,tentacle,cephalic", "rhinophore:2.0")
    tpms = {"rhinophore": 400.0, "tentacle": 4.0}
    boosted = _score(score_fn, tpms, weights, specific=True, tau=0.9)
    plain = _score(score_fn, tpms, weights, specific=False, tau=0.9)
    assert boosted / plain == pytest.approx(1.75, abs=1e-9)
