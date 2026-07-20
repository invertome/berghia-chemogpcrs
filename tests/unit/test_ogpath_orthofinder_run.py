"""Every stage must resolve the SAME OrthoFinder run, and the per-OG FASTA
must be the intended sequences -- not whichever file ``find`` happened to
return first.

Two coupled defects.

**A. The per-OG FASTA was resolved by an unsorted ``find | head -1``.**
Stage 04 and stage 05 both did::

    og=$(find "${RESULTS_DIR}/orthogroups" -name "${base}.fa" -type f | head -1)

On a real tree ``OG0000000.fa`` exists under SEVERAL directories of the same
OrthoFinder run with genuinely different contents:

  * ``Results_*/Orthogroup_Sequences/OG0000000.fa``  -- the intended unaligned
    protein sequences, original headers. THE one every stage means.
  * ``Results_*/MultipleSequenceAlignments/OG0000000.fa`` -- the already-aligned
    MSA. Re-aligning an alignment produces a garbage tree and exits 0.
  * ``Results_*/WorkingDirectory/...`` -- sequences keyed by OrthoFinder's
    internal integer ids, so ``get_taxids_from_fasta`` yields nothing and
    stage 05's ``[ "$taxa_count" -gt 1 ]`` gate silently skips selection
    analysis for that orthogroup.
  * zero-byte leftovers.

``find`` order is filesystem-dependent, so which of those a stage got was not
even stable between runs on the same tree. Resolution is now anchored to the
authoritative run's ``Orthogroup_Sequences/`` and requires a non-empty file.

**B. Five call sites disagreed about which OrthoFinder run is authoritative.**
``Results_<Mon><DD>`` carries no year and is not lexicographically ordered:
``sort -r`` puts ``Results_Sep05`` above ``Results_Dec01`` (S > D) although
December is later, and it cannot order ``Results_Dec01`` (last year) against
``Results_Jan05`` (this year) at all. Same-day reruns add ``_1``..``_10``,
which lexicographic sorting also gets wrong (``_10`` < ``_2``). Meanwhile
stages 04/06c/07/03b used an *unsorted* ``find ... | head -1``, so the same
orthogroup id could denote different gene sets in different stages while
every join still succeeded.

The adopted rule is: **the run whose ``Orthogroups/Orthogroups.tsv`` has the
newest mtime, ties broken by reverse path order.** mtime is the only
chronologically correct signal available -- the directory name is missing the
year entirely -- and requiring ``Orthogroups.tsv`` to exist means a crashed
or half-written run can never win.
"""
from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
HELPER = PROJECT_ROOT / "scripts" / "orthofinder_paths.sh"

OWNED_STAGES = [
    "03_orthology_clustering.sh",
    "03b_lse_classification.sh",
    "04_phylogenetic_analysis.sh",
    "05_selective_pressure_and_asr.sh",
    "06c_classify_non_chemoreceptors.sh",
    "07_candidate_ranking.sh",
]

OG_TSV_TEXT = (
    "Orthogroup\tberghia\tref_species\n"
    "OG0000000\tbste_cand_a, bste_cand_b\tref_1, ref_2\n"
)


# --------------------------------------------------------------------------- #
# harness
# --------------------------------------------------------------------------- #

def _bash(snippet: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["bash", "-c", f'set -o pipefail; source "{HELPER}"\n{snippet}'],
        capture_output=True, text=True,
    )


def _call(func: str, *args: str) -> subprocess.CompletedProcess:
    quoted = " ".join(f'"{a}"' for a in args)
    return _bash(f"{func} {quoted}")


def _make_run(root: Path, stamp: str, mtime: float,
              og_sequences=("OG0000000",), with_tsv: bool = True) -> Path:
    """Build one realistic OrthoFinder Results_<stamp> tree.

    Mirrors the layout ``orthofinder -f ${RESULTS_DIR}/orthogroups/input``
    produces, including the decoy copies of ``<OG>.fa`` that made the
    unsorted ``find`` non-deterministic.
    """
    run = root / "orthogroups" / "input" / "OrthoFinder" / stamp
    (run / "Orthogroups").mkdir(parents=True, exist_ok=True)
    (run / "Orthogroup_Sequences").mkdir(parents=True, exist_ok=True)
    (run / "MultipleSequenceAlignments").mkdir(parents=True, exist_ok=True)
    (run / "WorkingDirectory").mkdir(parents=True, exist_ok=True)

    for og in og_sequences:
        # 1. the intended unaligned sequences
        (run / "Orthogroup_Sequences" / f"{og}.fa").write_text(
            f">bste_cand_a\nMEALTKVFG\n>ref_1\nMEALSKVFG\n")
        # 2. the already-aligned MSA -- same basename, different content
        (run / "MultipleSequenceAlignments" / f"{og}.fa").write_text(
            f">bste_cand_a\nMEA-LTKVFG\n>ref_1\nMEA-LSKVFG\n")
        # 3. internal-integer-id sequences
        (run / "WorkingDirectory" / f"{og}.fa").write_text(
            ">0_1\nMEALTKVFG\n>1_7\nMEALSKVFG\n")
        # 4. a zero-byte leftover
        (run / "WorkingDirectory" / "decoy").mkdir(exist_ok=True)
        (run / "WorkingDirectory" / "decoy" / f"{og}.fa").write_text("")

    if with_tsv:
        tsv = run / "Orthogroups" / "Orthogroups.tsv"
        tsv.write_text(OG_TSV_TEXT)
        os.utime(tsv, (mtime, mtime))
    os.utime(run, (mtime, mtime))
    return run


@pytest.fixture
def results(tmp_path: Path) -> Path:
    d = tmp_path / "results"
    (d / "orthogroups").mkdir(parents=True)
    return d


def _root(results: Path) -> str:
    return str(results / "orthogroups")


# --------------------------------------------------------------------------- #
# the helper exists and is sourceable
# --------------------------------------------------------------------------- #

def test_helper_script_exists_and_parses():
    assert HELPER.exists(), f"missing shared helper: {HELPER}"
    proc = subprocess.run(["bash", "-n", str(HELPER)], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr


# --------------------------------------------------------------------------- #
# B. resolve_orthofinder_run -- one deterministic, chronological rule
# --------------------------------------------------------------------------- #

def test_single_run_resolves_to_itself(results: Path):
    run = _make_run(results, "Results_Jul19", mtime=1_700_000_000)

    proc = _call("resolve_orthofinder_run", _root(results))

    assert proc.returncode == 0, proc.stderr
    assert Path(proc.stdout.strip()).resolve() == run.resolve()


def test_newest_run_wins_when_the_name_sorts_the_wrong_way(results: Path):
    """THE finding: `sort -r` ranks Results_Sep05 above Results_Dec01 because
    'S' > 'D', but December is the later run."""
    _make_run(results, "Results_Sep05", mtime=1_700_000_000)
    december = _make_run(results, "Results_Dec01", mtime=1_707_000_000)

    proc = _call("resolve_orthofinder_run", _root(results))

    assert Path(proc.stdout.strip()).name == "Results_Dec01", proc.stdout
    assert Path(proc.stdout.strip()).resolve() == december.resolve()


def test_year_rollover_is_resolved_correctly(results: Path):
    """Results_<Mon><DD> has no year at all, so NO name-based rule can order
    a December run against the following January's."""
    _make_run(results, "Results_Dec01", mtime=1_700_000_000)   # last year
    january = _make_run(results, "Results_Jan05", mtime=1_710_000_000)  # this year

    proc = _call("resolve_orthofinder_run", _root(results))

    assert Path(proc.stdout.strip()).resolve() == january.resolve()


def test_same_day_rerun_suffix_is_ordered_numerically_not_lexically(results: Path):
    """Same-day reruns append _1.._10; lexicographically '_10' < '_2'."""
    _make_run(results, "Results_Jul19", mtime=1_700_000_000)
    _make_run(results, "Results_Jul19_2", mtime=1_700_001_000)
    tenth = _make_run(results, "Results_Jul19_10", mtime=1_700_002_000)

    proc = _call("resolve_orthofinder_run", _root(results))

    assert Path(proc.stdout.strip()).resolve() == tenth.resolve()


def test_run_without_orthogroups_tsv_never_wins(results: Path):
    """A crashed / half-written run is the NEWEST directory on disk. It must
    not become authoritative."""
    good = _make_run(results, "Results_Jul19", mtime=1_700_000_000)
    _make_run(results, "Results_Jul20", mtime=1_800_000_000, with_tsv=False)

    proc = _call("resolve_orthofinder_run", _root(results))

    assert Path(proc.stdout.strip()).resolve() == good.resolve()


def test_no_run_at_all_fails_nonzero_and_prints_nothing(results: Path):
    proc = _call("resolve_orthofinder_run", _root(results))

    assert proc.returncode != 0
    assert proc.stdout.strip() == ""


def test_missing_root_fails_nonzero(tmp_path: Path):
    proc = _call("resolve_orthofinder_run", str(tmp_path / "nope"))
    assert proc.returncode != 0


def test_identical_mtimes_break_ties_deterministically(results: Path):
    """Two runs written in the same second must still resolve to the same
    answer on every invocation, or stages drift apart again."""
    _make_run(results, "Results_Jul19", mtime=1_700_000_000)
    _make_run(results, "Results_Jul20", mtime=1_700_000_000)

    answers = {_call("resolve_orthofinder_run", _root(results)).stdout.strip()
               for _ in range(5)}

    assert len(answers) == 1, f"non-deterministic tie-break: {answers}"


def test_resolution_is_stable_across_repeated_calls(results: Path):
    _make_run(results, "Results_Sep05", mtime=1_700_000_000)
    _make_run(results, "Results_Dec01", mtime=1_707_000_000)

    answers = {_call("resolve_orthofinder_run", _root(results)).stdout.strip()
               for _ in range(5)}

    assert len(answers) == 1


# --------------------------------------------------------------------------- #
# derived resolvers -- every stage's need served off the SAME run
# --------------------------------------------------------------------------- #

def test_resolve_orthogroups_tsv_points_into_the_authoritative_run(results: Path):
    _make_run(results, "Results_Sep05", mtime=1_700_000_000)
    december = _make_run(results, "Results_Dec01", mtime=1_707_000_000)

    proc = _call("resolve_orthogroups_tsv", _root(results))

    assert proc.returncode == 0, proc.stderr
    tsv = Path(proc.stdout.strip())
    assert tsv.name == "Orthogroups.tsv"
    assert tsv.resolve() == (december / "Orthogroups" / "Orthogroups.tsv").resolve()


def test_resolve_orthogroup_sequences_dir_points_into_the_same_run(results: Path):
    _make_run(results, "Results_Sep05", mtime=1_700_000_000)
    december = _make_run(results, "Results_Dec01", mtime=1_707_000_000)

    proc = _call("resolve_orthogroup_sequences_dir", _root(results))

    assert proc.returncode == 0, proc.stderr
    assert Path(proc.stdout.strip()).resolve() == \
        (december / "Orthogroup_Sequences").resolve()


def test_all_resolvers_agree_on_one_run(results: Path):
    """The join-safety property: the OG id a stage reads from Orthogroups.tsv
    and the sequences another stage aligns must come from ONE run."""
    _make_run(results, "Results_Sep05", mtime=1_700_000_000)
    _make_run(results, "Results_Dec01", mtime=1_707_000_000)
    root = _root(results)

    tsv = Path(_call("resolve_orthogroups_tsv", root).stdout.strip())
    seqs = Path(_call("resolve_orthogroup_sequences_dir", root).stdout.strip())
    fasta = Path(_call("resolve_orthogroup_fasta", "OG0000000", root).stdout.strip())

    run = tsv.parent.parent
    assert seqs.parent.resolve() == run.resolve()
    assert fasta.parent.parent.resolve() == run.resolve()


# --------------------------------------------------------------------------- #
# A. resolve_orthogroup_fasta -- the intended sequences, never a decoy
# --------------------------------------------------------------------------- #

def test_orthogroup_fasta_is_the_unaligned_orthogroup_sequences_copy(results: Path):
    run = _make_run(results, "Results_Jul19", mtime=1_700_000_000)

    proc = _call("resolve_orthogroup_fasta", "OG0000000", _root(results))

    assert proc.returncode == 0, proc.stderr
    assert Path(proc.stdout.strip()).resolve() == \
        (run / "Orthogroup_Sequences" / "OG0000000.fa").resolve()


def test_orthogroup_fasta_is_never_the_multiple_sequence_alignment(results: Path):
    """Re-aligning an MSA yields a garbage tree and exits 0 -- silent."""
    _make_run(results, "Results_Jul19", mtime=1_700_000_000)

    resolved = _call("resolve_orthogroup_fasta", "OG0000000",
                     _root(results)).stdout.strip()

    assert "MultipleSequenceAlignments" not in resolved
    assert "-" not in Path(resolved).read_text(), "resolved an aligned file"


def test_orthogroup_fasta_is_never_the_internal_integer_id_copy(results: Path):
    """The WorkingDirectory copy is keyed on OrthoFinder's internal integer
    ids, so taxid extraction yields nothing and stage 05's `-gt 1` gate
    silently skips selection analysis."""
    _make_run(results, "Results_Jul19", mtime=1_700_000_000)

    resolved = _call("resolve_orthogroup_fasta", "OG0000000",
                     _root(results)).stdout.strip()

    assert "WorkingDirectory" not in resolved
    assert ">bste_cand_a" in Path(resolved).read_text()


def test_orthogroup_fasta_rejects_a_zero_byte_file(results: Path):
    run = _make_run(results, "Results_Jul19", mtime=1_700_000_000)
    (run / "Orthogroup_Sequences" / "OG0000000.fa").write_text("")

    proc = _call("resolve_orthogroup_fasta", "OG0000000", _root(results))

    assert proc.returncode != 0
    assert proc.stdout.strip() == ""


def test_orthogroup_fasta_absent_og_fails_nonzero(results: Path):
    _make_run(results, "Results_Jul19", mtime=1_700_000_000)

    proc = _call("resolve_orthogroup_fasta", "OG9999999", _root(results))

    assert proc.returncode != 0
    assert proc.stdout.strip() == ""


def test_orthogroup_fasta_comes_from_the_newest_run_not_a_stale_one(results: Path):
    """The identity hazard: the same OG id in an older run denotes a
    different gene set."""
    stale = _make_run(results, "Results_Sep05", mtime=1_700_000_000)
    (stale / "Orthogroup_Sequences" / "OG0000000.fa").write_text(
        ">stale_gene\nMKKKKK\n")
    fresh = _make_run(results, "Results_Dec01", mtime=1_707_000_000)

    resolved = Path(_call("resolve_orthogroup_fasta", "OG0000000",
                          _root(results)).stdout.strip())

    assert resolved.resolve() == \
        (fresh / "Orthogroup_Sequences" / "OG0000000.fa").resolve()
    assert "stale_gene" not in resolved.read_text()


# --------------------------------------------------------------------------- #
# the owned stages actually route through the helper
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("stage", OWNED_STAGES)
def test_owned_stage_sources_the_shared_helper(stage):
    text = (PROJECT_ROOT / stage).read_text()
    assert "orthofinder_paths.sh" in text, \
        f"{stage} does not source the shared OrthoFinder path helper"


@pytest.mark.parametrize("stage", OWNED_STAGES)
def test_owned_stage_has_no_adhoc_orthofinder_run_resolution(stage):
    """No owned stage may re-derive the authoritative run for itself -- that
    divergence is the bug."""
    text = (PROJECT_ROOT / stage).read_text()
    offenders = [
        line.strip() for line in text.splitlines()
        if not line.lstrip().startswith("#")
        and re.search(r"find\s+.*(Orthogroups\.tsv|Orthogroup_Sequences|"
                      r'-name\s+"Results_\*"|-name\s+.Results_\*.)', line)
    ]
    assert not offenders, f"{stage} still resolves OrthoFinder paths ad hoc: {offenders}"


@pytest.mark.parametrize("stage", ["04_phylogenetic_analysis.sh",
                                    "05_selective_pressure_and_asr.sh"])
def test_stage_no_longer_resolves_per_og_fasta_by_unsorted_find(stage):
    text = (PROJECT_ROOT / stage).read_text()
    offenders = [
        line.strip() for line in text.splitlines()
        if not line.lstrip().startswith("#")
        and "find " in line and '.fa"' in line and "head -1" in line
    ]
    assert not offenders, f"{stage} still uses `find | head -1` for the OG FASTA: {offenders}"


@pytest.mark.parametrize("stage", ["04_phylogenetic_analysis.sh",
                                    "05_selective_pressure_and_asr.sh"])
def test_no_manifest_fallback_glob_points_at_a_directory_that_exists(stage):
    """Same defect class, two lines up: the no-manifest fallback globbed
    ``orthogroups/OrthoFinder/Results*/Orthogroups/OG*.fa`` -- wrong on BOTH
    components. The real root is ``orthogroups/input/OrthoFinder`` and the
    per-OG FASTAs live in ``Orthogroup_Sequences/``, not ``Orthogroups/``
    (which holds only the TSV/TXT tables). The glob therefore never matched,
    bash left it literal, and ``basename`` turned the unexpanded pattern into
    the orthogroup name.
    """
    text = (PROJECT_ROOT / stage).read_text()
    offenders = [
        line.strip() for line in text.splitlines()
        if not line.lstrip().startswith("#")
        and "OrthoFinder/Results" in line and "OG*.fa" in line
    ]
    assert not offenders, f"{stage} still uses the non-matching fallback glob: {offenders}"


@pytest.mark.parametrize("stage", OWNED_STAGES)
def test_owned_stage_still_parses(stage):
    proc = subprocess.run(["bash", "-n", str(PROJECT_ROOT / stage)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr


def test_stage_03_no_longer_sorts_results_dirs_lexicographically():
    """`sort -r` on Results_<Mon><DD> is neither chronological nor numeric."""
    text = (PROJECT_ROOT / "03_orthology_clustering.sh").read_text()
    offenders = [
        line.strip() for line in text.splitlines()
        if not line.lstrip().startswith("#")
        and "Results_" in line and "sort -r" in line
    ]
    assert not offenders, f"stage 03 still lexicographically sorts runs: {offenders}"
