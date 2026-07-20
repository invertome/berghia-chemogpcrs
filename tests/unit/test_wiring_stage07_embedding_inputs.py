"""Stage 07's embedding-channel inputs: right reference set, degrading optionals.

Two independent defects on the same producer call, both of which end in the
embedding channel going dormant or scoring against the wrong thing, and neither
of which surfaces as a failure because the whole call is wrapped in
``|| log --level=WARN "... (channel stays dormant)"``.

1. HARDCODED REFERENCE LABELS. The legacy ESM-C branch passed
   ``--ref-labels ${REFERENCE_DIR}/anchors/anchor_set.tsv`` literally, bypassing
   ``${EMB_REF_LABELS}`` whose default was corrected to the 1094-row
   ``anchor_set_PROD.tsv``. anchor_set.tsv is the stale 206-row June set --
   scoring novelty against it is scoring against a different reference universe.
   Latent today (the locked default is ``EMB_SCORER=consensus``, so that branch
   is unreachable), which is precisely why it needs pinning: it is a landmine
   that arms itself the moment anyone flips EMB_SCORER for a comparison run.
   ``scripts/unity/rebuild_embedding_channel.sh`` already treats pointing at
   anchor_set.tsv as FATAL; stage 07 must not quietly do it.

2. AN OPTIONAL INPUT THAT RAISES. Stage 07 appended
   ``${EMB_CANDIDATE_IDENTITY_TSV}`` to EVERY model spec unconditionally, but
   ``_read_identity`` degraded gracefully only on an EMPTY path -- on a MISSING
   one it went straight into ``open()`` and raised FileNotFoundError. That file
   has never been generated in this checkout. In production the exception is
   swallowed by the ``|| log --level=WARN`` and the ENTIRE embedding channel --
   the locked proteinclip3b + protrek consensus, not just the optional identity
   confound -- silently goes dormant. A real job died on exactly this.
   ``scripts/unity/rebuild_embedding_channel.sh`` already guards it with
   ``[ -f ]`` before appending; stage 07 must do the same, and the library
   function must degrade rather than raise regardless of its caller.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
STAGE_07 = REPO_ROOT / "07_candidate_ranking.sh"
CONFIG_SH = REPO_ROOT / "config.sh"

# conftest.py adds scripts/ to sys.path
import embedding_candidate_diagnostics as ecd


# --------------------------------------------------------------------------- #
# 1. --ref-labels must come from ${EMB_REF_LABELS}
# --------------------------------------------------------------------------- #

def test_config_default_ref_labels_is_the_prod_anchor_set():
    """Pin the locked production reference set (1094-row PROD)."""
    text = CONFIG_SH.read_text()
    m = re.search(r'^export EMB_REF_LABELS=.*$', text, re.M)
    assert m, "EMB_REF_LABELS not exported from config.sh"
    assert "anchor_set_PROD.tsv" in m.group(0), m.group(0)


def test_stage07_never_hardcodes_the_stale_anchor_set_for_ref_labels():
    """No --ref-labels may point at a literal anchor_set.tsv path."""
    text = STAGE_07.read_text()
    offenders = [
        line.strip() for line in text.splitlines()
        if "--ref-labels" in line and "anchors/anchor_set.tsv" in line
    ]
    assert not offenders, (
        "stage 07 hardcodes the stale 206-row anchor_set.tsv instead of "
        f"${{EMB_REF_LABELS}} (default anchor_set_PROD.tsv): {offenders}"
    )


def test_every_ref_labels_argument_in_stage07_uses_the_variable():
    text = STAGE_07.read_text()
    args = re.findall(r'--ref-labels\s+("?[^\s\\]+"?)', text)
    assert args, "no --ref-labels arguments found in stage 07"
    for arg in args:
        assert "EMB_REF_LABELS" in arg, (
            f"--ref-labels {arg} does not use ${{EMB_REF_LABELS}}"
        )


def test_ref_labels_is_read_via_the_variable_in_the_helper_invocation():
    """The A1 ref-id extraction reads the same variable, so the tree-distance
    confound and the scorer cannot disagree about the reference universe."""
    text = STAGE_07.read_text()
    assert "load_ref_labels('${EMB_REF_LABELS}')" in text


# --------------------------------------------------------------------------- #
# 2a. _read_identity must DEGRADE on a missing optional input, never raise
# --------------------------------------------------------------------------- #

def test_read_identity_empty_path_returns_empty(tmp_path):
    """Pre-existing behaviour, pinned so the fix does not regress it."""
    assert ecd._read_identity("") == {}
    assert ecd._read_identity(None) == {}


def test_read_identity_missing_path_degrades_instead_of_raising(tmp_path):
    """The defect: a MISSING path raised FileNotFoundError into a caller whose
    exception handler is `|| log --level=WARN`, voiding the whole channel."""
    missing = str(tmp_path / "candidate_ref_identity_PROD.tsv")

    result = ecd._read_identity(missing)  # must not raise

    assert result == {}


def test_read_identity_missing_path_is_visible_on_stderr(tmp_path, capsys):
    """Degrading must be EXPLICIT: a silent {} is how an input goes missing for
    weeks without anyone noticing."""
    missing = str(tmp_path / "candidate_ref_identity_PROD.tsv")

    ecd._read_identity(missing)

    err = capsys.readouterr().err
    assert missing in err, err
    assert "identity" in err.lower(), err


def test_read_identity_directory_path_degrades(tmp_path, capsys):
    """A path that exists but is not a readable file (e.g. a directory left by
    a half-finished producer) must degrade the same way, not raise IsADirectory."""
    d = tmp_path / "not_a_file"
    d.mkdir()

    assert ecd._read_identity(str(d)) == {}
    assert str(d) in capsys.readouterr().err


def test_read_identity_reads_a_real_file(tmp_path):
    """The happy path is untouched: candidate_id <TAB> max_pct_identity."""
    p = tmp_path / "identity.tsv"
    p.write_text(
        "BersteEVm000001t1\t42.5\n"
        "BersteEVm000002t1\t93.1\n"
        "malformed_row_without_a_value\n"
        "BersteEVm000003t1\tnot_a_number\n"
    )

    result = ecd._read_identity(str(p))

    assert result == {"BersteEVm000001t1": 42.5, "BersteEVm000002t1": 93.1}


def test_read_identity_empty_file_returns_empty(tmp_path):
    p = tmp_path / "identity.tsv"
    p.write_text("")
    assert ecd._read_identity(str(p)) == {}


# --------------------------------------------------------------------------- #
# 2b. Stage 07 must not append a missing optional input to the model spec
# --------------------------------------------------------------------------- #

def _spec_block(text: str) -> str:
    """The consensus-branch loop that builds the tag:cand:ref[:identity] specs."""
    start = text.index("_emb_specs=()")
    end = text.index("_emb_residual_args=()", start)
    return text[start:end]


def test_stage07_guards_the_optional_identity_tsv_before_appending():
    """The spec may only carry the identity field when the file actually exists."""
    block = _spec_block(STAGE_07.read_text())
    assert "EMB_CANDIDATE_IDENTITY_TSV" in block, block

    spec_lines = [ln.strip() for ln in block.splitlines()
                  if "_emb_specs+=(" in ln]
    assert spec_lines, block
    # There must be a spec form WITHOUT the identity field (the degraded one).
    assert any("EMB_CANDIDATE_IDENTITY_TSV" not in ln for ln in spec_lines), (
        "every _emb_specs append carries ${EMB_CANDIDATE_IDENTITY_TSV} "
        f"unconditionally; a missing optional input still voids the channel: {spec_lines}"
    )
    # And the identity-bearing form must be guarded by an existence test.
    assert re.search(r'\[\s*-f\s+"\$\{EMB_CANDIDATE_IDENTITY_TSV\}"\s*\]', block), (
        "no [ -f \"${EMB_CANDIDATE_IDENTITY_TSV}\" ] guard in the spec block"
    )


def test_stage07_logs_when_the_optional_identity_tsv_is_absent():
    """Degrading must be visible in the pipeline log, not inferred from silence."""
    block = _spec_block(STAGE_07.read_text())
    guard_idx = block.index("EMB_CANDIDATE_IDENTITY_TSV")
    assert "log " in block[guard_idx:], (
        "the identity-TSV guard does not log its degraded path"
    )
    assert "identity" in block.lower()


def test_stage07_does_not_let_a_missing_identity_tsv_reduce_the_model_count():
    """The channel needs >=2 models. Degrading one model's optional confound
    must never drop the model itself -- that is what turns a missing optional
    input into a dormant channel."""
    block = _spec_block(STAGE_07.read_text())
    # The only condition that may set _emb_missing / skip a model is the npz
    # existence check, not the identity TSV.
    for line in block.splitlines():
        if "_emb_missing=1" in line or "-- excluded" in line:
            assert "IDENTITY" not in line.upper(), line
