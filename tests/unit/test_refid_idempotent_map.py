"""Tests pinning WRITE-ONCE identity semantics for scripts/update_headers.py.

Why these tests exist
---------------------
``update_headers.py`` is the pipeline's ID minter: it turns source accessions
into the ``ref_<taxid>_<N>`` short IDs that every downstream stage joins on.
It had three defects that all share one root cause -- the minter assigned IDs
*positionally* and had no path that reuses an existing map.

1. **protein/CDS misalignment.** Stage 01 runs the script twice against one
   map: proteins, then CDS with ``--append``. ``main()`` never passed
   ``existing_counters`` (the parameter existed and was honored, but was
   discarded at the call site), and the two runs are separate processes, so
   the module-level per-taxid counters restarted at 1 on the CDS pass. CDS
   recovery is partial by design (``CDS_PROTEIN_IDENTITY_MIN=0.95`` drops
   hits), so the CDS file is a same-order SUBSET and the counters shift:
   protein ``ACC2`` is ``ref_phau_2`` while the CDS written as ``ref_phau_2``
   is really ``ACC3``. ``find_nucleotide_sequences`` looks up by short ID,
   SUCCEEDS, and returns a different gene's nucleotides, which flow through
   MACSE/pal2nal into HyPhy and the dN/dS ranking axis. No crash, no warning.

2. **Stage 01 was not idempotent.** Stage 01 moves ``*_updated.fa`` over its
   own input, so a second run feeds minted short IDs back into the minter.
   ``extract_taxid('ref_phau_2')`` returns ``'ref'`` (``parts[0]``, alnum, in
   the 2-10 window), the double-prefix guard tests ``startswith('ref_')`` and
   ``'ref'`` has no trailing underscore so it never fires, and every taxon
   collapses into one ``ref`` namespace renumbered by position. The legacy
   ``preliminary/results/reference_sequences/id_map.csv`` is the receipt:
   108,895 of 109,334 rows have ``ref_ref_`` short IDs and 372 of its 373
   distinct "taxids" are exactly 10 characters (a fallback), i.e. a re-run
   already happened and destroyed the taxonomy.

3. **The minter re-issued IDs over a pruned set.** Prune one record and every
   subsequent short ID shifted to denote a different sequence, while both the
   old and new files contained the same ID -- so every join SUCCEEDED and
   returned confidently wrong results.

The invariant these tests pin
-----------------------------
IDs are WRITE-ONCE. Assignment must be idempotent: look for an existing map
first; if one exists, REUSE it, assert every retained record resolves, and
ABORT LOUDLY on any unmapped record. Never silently fall back to minting.
Genuine new minting into an existing map is gated behind ``ALLOW_NEW_IDS=1``
and refuses by default. Gaps left by pruning are CORRECT and are never closed
by renumbering.

No test here mints an ID as an expected outcome except via the explicit
opt-in; the rest assert that prior identity is preserved unchanged.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
UPDATE_HEADERS = REPO_ROOT / "scripts" / "update_headers.py"


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def write_fasta(path: Path, records) -> None:
    """Write ``records`` -- an iterable of (header, sequence) -- as FASTA."""
    path.write_text("".join(f">{h}\n{s}\n" for h, s in records))


def read_fasta_ids(path: Path):
    return [
        line[1:].split()[0]
        for line in path.read_text().splitlines()
        if line.startswith(">")
    ]


def read_map(path: Path):
    """Return the map as {original_id: short_id}, preserving nothing else."""
    rows = path.read_text().strip().splitlines()
    out = {}
    for row in rows[1:]:
        fields = row.split(",")
        out[fields[0]] = fields[1]
    return out


def run_update(fasta: Path, id_map: Path, *extra, env_extra=None):
    """Invoke the minter as stage 01 does, returning the CompletedProcess."""
    env = dict(os.environ)
    if env_extra:
        env.update(env_extra)
    return subprocess.run(
        [sys.executable, str(UPDATE_HEADERS), str(fasta), str(id_map),
         "--source-type", "reference", *extra],
        capture_output=True, text=True, env=env,
    )


def apply_stage01_move(fasta: Path) -> None:
    """Stage 01 moves the updated FASTA over its own input; replicate that."""
    updated = Path(f"{fasta}_updated.fa")
    assert updated.exists(), "minter did not produce an updated FASTA"
    fasta.write_text(updated.read_text())
    updated.unlink()


PROTEINS = [
    ("phau_ACC1", "MAAA"),
    ("phau_ACC2", "MBBB"),
    ("phau_ACC3", "MCCC"),
    ("phau_ACC4", "MDDD"),
]


@pytest.fixture()
def refs(tmp_path):
    """A fresh, first-ever run: proteins minted, map created."""
    fasta = tmp_path / "all_references.fa"
    id_map = tmp_path / "id_map.csv"
    write_fasta(fasta, PROTEINS)
    proc = run_update(fasta, id_map)
    assert proc.returncode == 0, proc.stderr
    apply_stage01_move(fasta)
    return fasta, id_map


# --------------------------------------------------------------------------
# Bootstrap: with no prior map there is nothing to re-issue, so minting is the
# correct behavior and must not require the opt-in.
# --------------------------------------------------------------------------
def test_first_run_without_existing_map_assigns_ids(refs):
    fasta, id_map = refs
    mapping = read_map(id_map)
    assert mapping == {
        "phau_ACC1": "ref_phau_1",
        "phau_ACC2": "ref_phau_2",
        "phau_ACC3": "ref_phau_3",
        "phau_ACC4": "ref_phau_4",
    }
    assert read_fasta_ids(fasta) == [
        "ref_phau_1", "ref_phau_2", "ref_phau_3", "ref_phau_4",
    ]


# --------------------------------------------------------------------------
# Finding #1 -- the dN/dS corruption.
# --------------------------------------------------------------------------
def test_append_on_partial_cds_subset_agrees_with_protein_ids(refs, tmp_path):
    """CDS recovery is partial; the surviving CDS must keep the protein's ID.

    ACC2 fails the identity filter and is absent from the CDS file. A
    positional minter restarts at 1 and writes ACC3's nucleotides under
    ``ref_phau_2`` -- the protein ``ref_phau_2`` is ACC2. That silent swap is
    what reaches HyPhy.
    """
    _, id_map = refs
    cds = tmp_path / "all_references_cds.fna"
    write_fasta(cds, [("phau_ACC1", "ATGAAA"),
                      ("phau_ACC3", "ATGCCC"),
                      ("phau_ACC4", "ATGGGG")])

    proc = run_update(cds, id_map, "--append")
    assert proc.returncode == 0, proc.stderr

    cds_ids = read_fasta_ids(Path(f"{cds}_updated.fa"))
    assert cds_ids == ["ref_phau_1", "ref_phau_3", "ref_phau_4"]


def test_append_never_reassigns_a_short_id_to_a_second_original(refs, tmp_path):
    """No short ID may denote two different source records across passes."""
    _, id_map = refs
    before = read_map(id_map)
    cds = tmp_path / "all_references_cds.fna"
    write_fasta(cds, [("phau_ACC1", "ATGAAA"), ("phau_ACC3", "ATGCCC")])

    proc = run_update(cds, id_map, "--append")
    assert proc.returncode == 0, proc.stderr

    after = read_map(id_map)
    for original, short in before.items():
        assert after.get(original) == short, (
            f"{original} was re-issued: {short} -> {after.get(original)}"
        )
    # A short ID must resolve to exactly one original record.
    shorts = list(after.values())
    assert len(shorts) == len(set(shorts)), f"duplicate short IDs: {shorts}"


# --------------------------------------------------------------------------
# Finding #2 -- stage 01 must be idempotent.
# --------------------------------------------------------------------------
def test_rerunning_stage01_is_a_noop(refs):
    """Re-running over already-mapped input changes neither map nor FASTA."""
    fasta, id_map = refs
    map_before = id_map.read_text()
    ids_before = read_fasta_ids(fasta)

    proc = run_update(fasta, id_map)
    assert proc.returncode == 0, proc.stderr
    apply_stage01_move(fasta)

    assert id_map.read_text() == map_before
    assert read_fasta_ids(fasta) == ids_before


def test_rerun_never_double_prefixes_or_collapses_taxonomy(refs):
    """The ``ref_ref_*`` / taxid=='ref' collapse must not reappear."""
    fasta, id_map = refs
    proc = run_update(fasta, id_map)
    assert proc.returncode == 0, proc.stderr
    apply_stage01_move(fasta)

    assert "ref_ref_" not in id_map.read_text()
    assert not any(i.startswith("ref_ref_") for i in read_fasta_ids(fasta))
    taxids = {row.split(",")[2] for row in id_map.read_text().splitlines()[1:]}
    assert taxids == {"phau"}, f"taxonomy collapsed: {taxids}"


# --------------------------------------------------------------------------
# Finding #3 -- a pruned re-run must reuse, not re-issue.
# --------------------------------------------------------------------------
def test_pruned_rerun_reuses_prior_map_and_preserves_gaps(refs):
    """Dropping ACC2 must leave a GAP, not shift ACC3/ACC4 down."""
    fasta, id_map = refs
    map_before = id_map.read_text()

    kept = [(h, s) for h, s in
            zip(read_fasta_ids(fasta), ["MAAA", "MBBB", "MCCC", "MDDD"])
            if h != "ref_phau_2"]
    write_fasta(fasta, kept)

    proc = run_update(fasta, id_map)
    assert proc.returncode == 0, proc.stderr
    apply_stage01_move(fasta)

    assert read_fasta_ids(fasta) == ["ref_phau_1", "ref_phau_3", "ref_phau_4"]
    assert id_map.read_text() == map_before, "pruning must not rewrite the map"


def test_unmapped_record_aborts_loudly_and_leaves_map_untouched(refs):
    """An unresolvable record must fail hard, never fall back to minting."""
    fasta, id_map = refs
    map_before = id_map.read_text()
    records = list(zip(read_fasta_ids(fasta), ["MAAA", "MBBB", "MCCC", "MDDD"]))
    write_fasta(fasta, records + [("phau_ACC9", "MEEE")])

    proc = run_update(fasta, id_map)

    assert proc.returncode != 0, "unmapped record was silently accepted"
    combined = proc.stdout + proc.stderr
    assert "phau_ACC9" in combined, combined
    assert "ALLOW_NEW_IDS" in combined, "abort must name the opt-in gate"
    assert id_map.read_text() == map_before
    assert not Path(f"{fasta}_updated.fa").exists(), (
        "aborting run must not leave a partial FASTA"
    )


def test_allow_new_ids_opt_in_extends_without_reissuing(refs):
    """The opt-in mints ABOVE the high-water mark; prior IDs are untouched."""
    fasta, id_map = refs
    before = read_map(id_map)
    records = list(zip(read_fasta_ids(fasta), ["MAAA", "MBBB", "MCCC", "MDDD"]))
    write_fasta(fasta, records + [("phau_ACC9", "MEEE")])

    proc = run_update(fasta, id_map, env_extra={"ALLOW_NEW_IDS": "1"})
    assert proc.returncode == 0, proc.stderr
    apply_stage01_move(fasta)

    after = read_map(id_map)
    for original, short in before.items():
        assert after[original] == short
    assert after["phau_ACC9"] == "ref_phau_5"


def test_stage01_propagates_the_cds_refusal(refs):
    """Stage 01's CDS pass must fail the stage, not fall through to the mv.

    The bare ``python3 update_headers.py ... --append`` ignored its exit
    status, so a refusal surfaced only as a confusing "failed to move updated
    CDS file" two lines later.
    """
    stage01 = (REPO_ROOT / "01_reference_processing.sh").read_text()
    cds_block = stage01.split("Update CDS headers to match protein headers")[1]
    invocation = cds_block.split("--append")[0]
    assert "if ! python3" in invocation, (
        "stage 01 must check update_headers.py's exit status on the CDS pass:\n"
        + invocation
    )


def test_opt_in_after_pruning_does_not_recycle_a_retired_id(refs):
    """A retired ID must never be handed to a different record."""
    fasta, id_map = refs
    records = list(zip(read_fasta_ids(fasta), ["MAAA", "MBBB", "MCCC", "MDDD"]))
    pruned = [r for r in records if r[0] != "ref_phau_4"]
    write_fasta(fasta, pruned + [("phau_NEW", "MEEE")])

    proc = run_update(fasta, id_map, env_extra={"ALLOW_NEW_IDS": "1"})
    assert proc.returncode == 0, proc.stderr

    after = read_map(id_map)
    assert after["phau_ACC4"] == "ref_phau_4", "retired mapping was dropped"
    assert after["phau_NEW"] != "ref_phau_4", "retired ID was recycled"
    assert after["phau_NEW"] == "ref_phau_5"
