"""Tests for ancestral-sequence record-id namespacing (scripts/extract_iqtree_asr.py).

The defect
----------
``extract_iqtree_asr.py`` wrote its FASTA header as the bare IQ-TREE internal
node label::

    f.write(f">{args.node_name}\\n{seq}\\n")            # line 65

IQ-TREE numbers internal nodes *per tree* (``Node1``, ``Node2``, ...), so two
different orthogroups near-certainly both emit ``>Node5`` denoting completely
different ancestral sequences. The output FILENAME was already namespaced
(``${base}_${node}_asr.fa``) but the record id inside it was not.

That blocks the stage-08 ancestral-folding path: stage 08 concatenates every
``*_asr.fa`` into one extraction source and guards it with
``assert_no_duplicate_fasta_ids``, which correctly fails loud on the
collision. Folding an arbitrary one of two different ancestral sequences under
the same id would corrupt the structural phylogeny, so failing is right --
but the ids have to be fixed for the path to work at all.

The fix
-------
Namespace at the PRODUCER: an ancestral record's id becomes
``{orthogroup}_{node}`` (e.g. ``OG0000001_Node5``), matching the existing
filename convention, so the producer owns a globally unique id and every
consumer gets unambiguous records.

The design constraint these tests pin
-------------------------------------
``node_name`` was doing two different jobs: it was the lookup key into the
IQ-TREE ``.state`` file AND the emitted FASTA header. Those are now separate.
The ``.state`` lookup must keep using the bare node label (``Node5``) because
that is what IQ-TREE writes; only the emitted record id changes. The record id
is a distinct OPTIONAL flag that defaults to the bare node name, so the
three-positional-arg CLI stays backward compatible.

Note that these ids are not newly minted: the orthogroup base and the node
label both already exist, and the record id is their composition. No record's
identity is re-issued.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
EXTRACTOR = REPO_ROOT / "scripts" / "extract_iqtree_asr.py"
STAGE05 = REPO_ROOT / "05_selective_pressure_and_asr.sh"
STAGE08 = REPO_ROOT / "08_structural_analysis.sh"
FUNCTIONS_SH = REPO_ROOT / "functions.sh"


# --------------------------------------------------------------------------
# Fixtures: a synthetic IQ-TREE .state file
# --------------------------------------------------------------------------

def _write_state_file(path: Path, node_sequences: dict[str, str]) -> None:
    """Write a minimal IQ-TREE ``--ancestral`` .state file.

    Real format: ``#``-prefixed comment banner, then a header row
    (``Node Site State p_A p_R ...``), then one tab-separated row per
    (node, site). The extractor takes column 2 (``State``) as the marginal ML
    state, so a single probability column is enough to keep rows well-formed.
    """
    lines = [
        "# Ancestral state reconstruction",
        "# This file is best viewed in a spreadsheet program",
        "Node\tSite\tState\tp_A\tp_R\tp_N",
    ]
    for node, seq in node_sequences.items():
        for offset, residue in enumerate(seq):
            site = offset + 1
            lines.append(f"{node}\t{site}\t{residue}\t0.9\t0.05\t0.05")
    path.write_text("\n".join(lines) + "\n")


def _run_extractor(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(EXTRACTOR), *args],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )


def _records(fasta: Path) -> list[tuple[str, str]]:
    """(bare id, sequence) pairs in file order."""
    out: list[tuple[str, str]] = []
    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            out.append((line[1:].split()[0], ""))
        elif line.strip() and out:
            ident, seq = out[-1]
            out[-1] = (ident, seq + line.strip())
    return out


# --------------------------------------------------------------------------
# The emitted record id
# --------------------------------------------------------------------------

def test_emitted_header_is_the_namespaced_record_id(tmp_path):
    """THE fix: the header carries the orthogroup-namespaced id, not `Node5`."""
    state = tmp_path / "OG0000001_asr.state"
    _write_state_file(state, {"Node5": "MKFAA"})
    out = tmp_path / "OG0000001_Node5_asr.fa"

    result = _run_extractor(
        str(state), "Node5", str(out), "--record-id", "OG0000001_Node5"
    )
    assert result.returncode == 0, result.stderr

    assert _records(out) == [("OG0000001_Node5", "MKFAA")]


def test_state_lookup_still_resolves_the_bare_node_label(tmp_path):
    """The .state key stays the bare label -- IQ-TREE writes `Node5`, not the
    namespaced id. A namespaced record id must not leak into the lookup, or the
    extractor would find zero sites and fail."""
    state = tmp_path / "OG0000001_asr.state"
    _write_state_file(state, {"Node5": "MKFAA", "Node9": "WWWWW"})
    out = tmp_path / "asr.fa"

    result = _run_extractor(
        str(state), "Node5", str(out), "--record-id", "OG0000001_Node5"
    )
    assert result.returncode == 0, result.stderr

    # Node5's sequence, not Node9's -- the lookup keyed on the bare label.
    assert _records(out) == [("OG0000001_Node5", "MKFAA")]


def test_namespaced_id_is_not_used_as_the_state_lookup_key(tmp_path):
    """Guard against the naive implementation that reuses one value for both
    jobs: if the record id were used as the lookup key it would match nothing
    and the extractor would exit non-zero."""
    state = tmp_path / "asr.state"
    _write_state_file(state, {"Node5": "MKFAA"})
    out = tmp_path / "asr.fa"

    result = _run_extractor(
        str(state), "Node5", str(out), "--record-id", "OG0000001_Node5"
    )
    assert result.returncode == 0, result.stderr
    assert "no sites found" not in result.stderr


def test_default_record_id_preserves_the_bare_node_name(tmp_path):
    """Backward compatibility: with no --record-id the old behaviour stands, so
    the three-positional CLI keeps working for any existing caller."""
    state = tmp_path / "asr.state"
    _write_state_file(state, {"Node5": "MKFAA"})
    out = tmp_path / "asr.fa"

    result = _run_extractor(str(state), "Node5", str(out))
    assert result.returncode == 0, result.stderr

    assert _records(out) == [("Node5", "MKFAA")]


def test_two_orthogroups_at_the_same_node_number_get_distinct_ids(tmp_path):
    """The collision that blocked stage 08: both trees have a `Node5`."""
    first_state = tmp_path / "OG0000001_asr.state"
    second_state = tmp_path / "OG0000002_asr.state"
    _write_state_file(first_state, {"Node5": "MKFAA"})
    _write_state_file(second_state, {"Node5": "WWCCY"})

    first_out = tmp_path / "OG0000001_Node5_asr.fa"
    second_out = tmp_path / "OG0000002_Node5_asr.fa"

    assert _run_extractor(
        str(first_state), "Node5", str(first_out), "--record-id", "OG0000001_Node5"
    ).returncode == 0
    assert _run_extractor(
        str(second_state), "Node5", str(second_out), "--record-id", "OG0000002_Node5"
    ).returncode == 0

    first_ids = [i for i, _ in _records(first_out)]
    second_ids = [i for i, _ in _records(second_out)]
    assert first_ids == ["OG0000001_Node5"]
    assert second_ids == ["OG0000002_Node5"]
    assert set(first_ids).isdisjoint(second_ids)

    # ...and each id still carries its OWN orthogroup's ancestral sequence.
    assert _records(first_out)[0][1] == "MKFAA"
    assert _records(second_out)[0][1] == "WWCCY"


def test_missing_node_still_fails_loud(tmp_path):
    """A record id must not paper over an absent node: no sites is still an
    error, and the message names the bare label that was actually looked up."""
    state = tmp_path / "asr.state"
    _write_state_file(state, {"Node5": "MKFAA"})
    out = tmp_path / "asr.fa"

    result = _run_extractor(
        str(state), "Node404", str(out), "--record-id", "OG0000001_Node404"
    )
    assert result.returncode != 0
    assert "Node404" in result.stderr
    assert not out.exists()


# --------------------------------------------------------------------------
# The stage-05 call site
# --------------------------------------------------------------------------

def _iqtree_asr_call(src: str) -> str:
    """The extract_iqtree_asr.py invocation, with line continuations folded."""
    flat = src.replace("\\\n", " ")
    match = re.search(r"[^\n]*extract_iqtree_asr\.py[^\n]*", flat)
    assert match, "05_selective_pressure_and_asr.sh no longer calls extract_iqtree_asr.py"
    return match.group(0)


def test_stage05_passes_the_namespaced_record_id():
    """Wiring: the producer is only fixed if the call site actually uses it.

    Stage 05 already has both halves in scope as ``${base}`` and ``${node}``.
    """
    call = _iqtree_asr_call(STAGE05.read_text())
    assert "--record-id" in call, (
        "stage 05 still calls the extractor without --record-id, so every "
        "orthogroup would emit bare `NodeN` ids and collide again"
    )
    assert "${base}_${node}" in call, (
        "the record id must be the orthogroup-namespaced form matching the "
        f"filename convention; got: {call}"
    )


def test_stage05_still_passes_the_bare_node_as_the_state_lookup_key():
    """The positional node argument must stay bare -- it keys the .state file."""
    call = _iqtree_asr_call(STAGE05.read_text())
    assert re.search(r'"\$\{asr_prefix\}\.state"\s+"\$node"', call), (
        f"the .state lookup key must remain the bare ${{node}}; got: {call}"
    )


# --------------------------------------------------------------------------
# End-to-end: the stage-08 duplicate guard now passes
# --------------------------------------------------------------------------

def _extract_shell_function(name: str, src: str) -> str:
    """Exact source text of one top-level ``name() { ... }``.

    Same extraction pattern as tests/unit/test_stage08_asr_folding.py: the
    logic under test is the shipping code, not a copy of it.
    """
    match = re.search(
        rf"^{re.escape(name)}\(\)\s*\{{$.*?^\}}$",
        src,
        flags=re.MULTILINE | re.DOTALL,
    )
    assert match, f"{STAGE08.name} does not define a shell function '{name}()'"
    return match.group(0)


def _run_stage08_guard(asr_dir: Path, tmp_path: Path) -> subprocess.CompletedProcess:
    """Run stage 08's real collect + build-source path over ``asr_dir``.

    ``build_structural_seq_source`` ends in ``assert_no_duplicate_fasta_ids``,
    so its exit status IS the duplicate-guard verdict.
    """
    src = STAGE08.read_text()
    fns = "\n".join(
        _extract_shell_function(n, src)
        for n in ("collect_asr_sequences", "build_structural_seq_source")
    )
    extant = tmp_path / "chemogpcrs_berghia.fa"
    extant.write_text(">cand1\nMKFAA\n>cand2\nMKVCC\n")
    asr_seqs = tmp_path / "asr_seqs.fa"
    combined = tmp_path / "extraction_source.fa"

    env = os.environ.copy()
    env["LOGS_DIR"] = str(tmp_path)
    body = (
        f'collect_asr_sequences "{asr_dir}" "{asr_seqs}"\n'
        f'build_structural_seq_source "{extant}" "{asr_seqs}" "{combined}"\n'
    )
    script = f'source "{FUNCTIONS_SH}"\ntrap - EXIT\n{fns}\n{body}'
    result = subprocess.run(
        ["bash", "-c", script], capture_output=True, text=True, env=env
    )
    result.combined_path = combined  # type: ignore[attr-defined]
    return result


def test_stage08_duplicate_guard_passes_on_a_two_orthogroup_asr_set(tmp_path):
    """THE point of the change, demonstrated end-to-end.

    Two orthogroups both reconstruct their own `Node5`. Run the REAL extractor
    on each, then feed the results through stage 08's REAL collect +
    build-source path, whose final act is `assert_no_duplicate_fasta_ids`.
    """
    asr_dir = tmp_path / "asr"
    asr_dir.mkdir()

    first_state = tmp_path / "OG0000001_asr.state"
    second_state = tmp_path / "OG0000002_asr.state"
    _write_state_file(first_state, {"Node5": "MKFAA"})
    _write_state_file(second_state, {"Node5": "WWCCY"})

    assert _run_extractor(
        str(first_state), "Node5", str(asr_dir / "OG0000001_Node5_asr.fa"),
        "--record-id", "OG0000001_Node5",
    ).returncode == 0
    assert _run_extractor(
        str(second_state), "Node5", str(asr_dir / "OG0000002_Node5_asr.fa"),
        "--record-id", "OG0000002_Node5",
    ).returncode == 0

    result = _run_stage08_guard(asr_dir, tmp_path)
    assert result.returncode == 0, (
        "stage 08's duplicate-id guard rejected the namespaced ASR set:\n"
        f"{result.stdout}\n{result.stderr}"
    )

    combined = result.combined_path  # type: ignore[attr-defined]
    ids = [i for i, _ in _records(combined)]
    assert ids == ["cand1", "cand2", "OG0000001_Node5", "OG0000002_Node5"]
    # Each ancestral id resolves to its OWN orthogroup's sequence.
    seqs = dict(_records(combined))
    assert seqs["OG0000001_Node5"] == "MKFAA"
    assert seqs["OG0000002_Node5"] == "WWCCY"


def test_stage08_duplicate_guard_would_have_failed_on_the_old_bare_ids(tmp_path):
    """The collision was real, not hypothetical.

    Same two orthogroups via the DEFAULT (un-namespaced) id -- i.e. the old
    behaviour -- must still trip the guard. This is what pins that the passing
    test above is caused by the namespacing rather than by some unrelated
    property of the fixture.
    """
    asr_dir = tmp_path / "asr"
    asr_dir.mkdir()

    first_state = tmp_path / "OG0000001_asr.state"
    second_state = tmp_path / "OG0000002_asr.state"
    _write_state_file(first_state, {"Node5": "MKFAA"})
    _write_state_file(second_state, {"Node5": "WWCCY"})

    assert _run_extractor(
        str(first_state), "Node5", str(asr_dir / "OG0000001_Node5_asr.fa")
    ).returncode == 0
    assert _run_extractor(
        str(second_state), "Node5", str(asr_dir / "OG0000002_Node5_asr.fa")
    ).returncode == 0

    result = _run_stage08_guard(asr_dir, tmp_path)
    assert result.returncode != 0, (
        "the guard should have rejected two different ancestral sequences "
        "sharing the id `Node5`"
    )
    assert "Duplicate FASTA IDs" in (result.stdout + result.stderr)
