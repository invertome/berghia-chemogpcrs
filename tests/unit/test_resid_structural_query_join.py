"""RESID REPAIR 3 — structural_channel() must join on candidate ids.

``structural_evidence.structural_channel()`` keyed its returned channel on
the Foldseek ``query`` field VERBATIM -- the same defect
``build_structural_channel.py`` was already fixed for. A raw search-tool
output field is not a join key.

MEASURED foldseek behaviour (recorded in build_structural_channel.py from
foldseek 10.941cd33 in this project's berghia-gpcr env on Unity, srun
61999969, with a flat query dir staged ``<candidate_id>.<ext>`` exactly as
``scripts/unity/run_foldseek_candidates.sh`` stages it):

    BersteEVm000001t1.cif  (1 chain)  -> query "BersteEVm000001t1"
    BersteEVm000002t1.pdb  (1 chain)  -> query "BersteEVm000002t1"
    BersteEVm000003t1.cif  (2 chains) -> queries "BersteEVm000003t1_A"
                                                 "BersteEVm000003t1_B"

Two consequences these tests pin down:

  * foldseek STRIPS the extension, so "<id>.cif" is a query form the tool
    NEVER emits. A fixture that writes one is testing an imagined tool.
  * a multi-chain model (AF3 receptor + peptide/G-protein, or any complex)
    emits ONE ROW PER CHAIN, and neither "<id>_A" nor "<id>_B" equals the
    bare candidate id -- so one candidate silently becomes two orphan rows
    that left-join to nothing.

Chain stripping cannot be blind: a candidate id may legitimately end in
"_<alnum>". It is applied only when corroborated against the real candidate
id universe -- resolution by identity, never by pattern-guess.

TARGET-side fixtures here use real foldseek target formats
("AF-P31356-F1-model_v4", "1f88_A"), never a bare accession: a bare
accession is a form the tool does not emit, and a fixture that uses one
passes while production matches nothing.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from structural_evidence import (  # noqa: E402
    canonical_query_id,
    query_keys,
    structural_channel,
)

# The real candidate id universe shape: Berghia RefSeq-style model ids.
CANDIDATE_IDS = {
    "BersteEVm000001t1",
    "BersteEVm000002t1",
    "BersteEVm000003t1",
}

# A real AFDB target for an opsin-family anchor, in the form foldseek prints.
AFDB_OPSIN_TARGET = "AF-P31356-F1-model_v4"
FAMILY_MAP = {"P31356": "opsin"}       # anchor_set.tsv is keyed on accession


def _write_hits(tmp_path: Path, lines: list[str]) -> str:
    p = tmp_path / "foldseek_hits.tsv"
    p.write_text("\n".join(lines) + "\n")
    return str(p)


# --------------------------------------------------------------------------
# the query-key helpers
# --------------------------------------------------------------------------
def test_query_keys_orders_raw_first_then_chain_stripped():
    keys = query_keys("BersteEVm000003t1_A")
    assert keys[0] == "BersteEVm000003t1_A", "the literal query must win first"
    assert "BersteEVm000003t1" in keys, "the chain-stripped form must be offered"


def test_query_keys_empty_input_yields_nothing():
    assert query_keys(None) == []
    assert query_keys("  ") == []


def test_canonical_query_id_resolves_a_chain_row_to_the_candidate():
    assert canonical_query_id("BersteEVm000003t1_A", CANDIDATE_IDS) == \
        "BersteEVm000003t1"


def test_canonical_query_id_prefers_an_exact_candidate_over_stripping():
    """A candidate whose id legitimately ends in '_<alnum>' must not be
    truncated into a different, shorter id that also happens to exist."""
    ids = {"Berghia_scaffold_12", "Berghia_scaffold"}
    assert canonical_query_id("Berghia_scaffold_12", ids) == "Berghia_scaffold_12"


def test_canonical_query_id_refuses_to_invent_a_key():
    """An unresolvable query is reported, never rewritten into something
    that merely looks joinable."""
    assert canonical_query_id("AF3_job_00017_model_0", CANDIDATE_IDS) is None


def test_canonical_query_id_without_a_universe_leaves_the_chain_alone():
    """Nothing to corroborate against -> no guessing."""
    assert canonical_query_id("Berghia_scaffold_12", None) == "Berghia_scaffold_12"


# --------------------------------------------------------------------------
# structural_channel — the actual repair
# --------------------------------------------------------------------------
def test_multichain_query_collapses_to_one_candidate_row(tmp_path):
    """THE defect. Two chain rows for one model must become ONE channel row
    keyed on the candidate id, not two orphans keyed on '<id>_A' / '<id>_B'.
    """
    hits = _write_hits(tmp_path, [
        f"BersteEVm000003t1_A\t{AFDB_OPSIN_TARGET}\t0.62\t0.71\t1e-30",
        f"BersteEVm000003t1_B\t{AFDB_OPSIN_TARGET}\t0.20\t0.31\t1e-3",
    ])

    channel = structural_channel(hits, FAMILY_MAP, candidate_ids=CANDIDATE_IDS)

    assert set(channel) == {"BersteEVm000003t1"}, (
        f"expected one row keyed on the candidate id, got {sorted(channel)} -- "
        f"these would left-join to nothing in rank_aggregation"
    )
    assert channel["BersteEVm000003t1"]["has_struct_data"] is True


def test_chain_collapse_keeps_the_strongest_hit(tmp_path):
    """Chain A here is the confident opsin match; chain B is noise. The
    collapsed row must reflect A, or a real exclusion signal is lost."""
    hits = _write_hits(tmp_path, [
        f"BersteEVm000003t1_B\tAF-Q9Y5N1-F1-model_v4\t0.10\t0.18\t1e-1",
        f"BersteEVm000003t1_A\t{AFDB_OPSIN_TARGET}\t0.62\t0.88\t1e-40",
    ])
    channel = structural_channel(hits, FAMILY_MAP, candidate_ids=CANDIDATE_IDS)

    row = channel["BersteEVm000003t1"]
    assert row["struct_state"] == "known_non_chemoreceptor"
    assert row["struct_nonchemo_corrob"] == 1
    assert row["struct_novelty"] == 0


def test_single_chain_query_is_unchanged(tmp_path):
    """The 1-chain case foldseek emits bare must keep working untouched."""
    hits = _write_hits(tmp_path, [
        f"BersteEVm000001t1\t{AFDB_OPSIN_TARGET}\t0.62\t0.88\t1e-40",
    ])
    channel = structural_channel(hits, FAMILY_MAP, candidate_ids=CANDIDATE_IDS)
    assert set(channel) == {"BersteEVm000001t1"}


def test_zero_key_overlap_raises_instead_of_going_dormant(tmp_path):
    """The cross-cutting invariant: a namespace mapping must assert non-zero
    key overlap on real data.

    AF3 job names are the realistic wrong-scheme case (bead 5ubd). With a
    verbatim key this returns a full, plausible-looking channel that
    left-joins to nothing -- the structural voter silently stops voting
    while every exit code stays 0.
    """
    hits = _write_hits(tmp_path, [
        f"AF3_job_00017_model_0\t{AFDB_OPSIN_TARGET}\t0.62\t0.88\t1e-40",
        f"AF3_job_00018_model_0\t{AFDB_OPSIN_TARGET}\t0.55\t0.80\t1e-35",
    ])
    with pytest.raises(ValueError, match="ZERO|zero"):
        structural_channel(hits, FAMILY_MAP, candidate_ids=CANDIDATE_IDS)


def test_empty_candidate_universe_raises(tmp_path):
    """An empty universe means the candidate FASTA could not be read; every
    query would fail to resolve for a reason that has nothing to do with
    the data."""
    hits = _write_hits(tmp_path, [
        f"BersteEVm000001t1\t{AFDB_OPSIN_TARGET}\t0.62\t0.88\t1e-40",
    ])
    with pytest.raises(ValueError, match="empty"):
        structural_channel(hits, FAMILY_MAP, candidate_ids=set())


def test_partially_unresolved_queries_are_dropped_and_named(tmp_path, capsys):
    hits = _write_hits(tmp_path, [
        f"BersteEVm000001t1\t{AFDB_OPSIN_TARGET}\t0.62\t0.88\t1e-40",
        f"AF3_job_00017_model_0\t{AFDB_OPSIN_TARGET}\t0.62\t0.88\t1e-40",
    ])
    channel = structural_channel(hits, FAMILY_MAP, candidate_ids=CANDIDATE_IDS)

    assert set(channel) == {"BersteEVm000001t1"}
    err = capsys.readouterr().err
    assert "AF3_job_00017_model_0" in err, (
        "a dropped query must be named, not silently discarded"
    )


def test_no_candidate_ids_preserves_the_legacy_behaviour(tmp_path):
    """Backward compatibility: existing callers pass no universe. Nothing to
    corroborate against means no chain stripping and no assertion."""
    hits = _write_hits(tmp_path, [
        f"cand_novel\ttargetX\t0.10\t0.20\t1e-2",
        f"BersteEVm000003t1_A\t{AFDB_OPSIN_TARGET}\t0.62\t0.88\t1e-40",
    ])
    channel = structural_channel(hits, FAMILY_MAP)
    assert set(channel) == {"cand_novel", "BersteEVm000003t1_A"}


def test_a_missing_file_still_returns_empty_even_with_a_universe(tmp_path):
    """No hits is not zero-overlap; there is nothing to overlap with."""
    missing = str(tmp_path / "nope.tsv")
    assert structural_channel(missing, FAMILY_MAP,
                              candidate_ids=CANDIDATE_IDS) == {}


# --------------------------------------------------------------------------
# fixture-realism guard
# --------------------------------------------------------------------------
def test_fixtures_use_target_ids_foldseek_actually_emits():
    """A test here once wrote a bare accession as a foldseek target -- a form
    the real tool never emits -- and passed while production matched nothing.

    Pin the fixture to the real AFDB format and prove the family lookup has
    to traverse target_keys() to succeed, exactly as it does in production.
    """
    from structural_evidence import family_for_target

    assert AFDB_OPSIN_TARGET not in FAMILY_MAP, (
        "the fixture must exercise target-id parsing, not sidestep it by "
        "keying the family map on the literal target"
    )
    assert family_for_target(AFDB_OPSIN_TARGET, FAMILY_MAP) == "opsin"
