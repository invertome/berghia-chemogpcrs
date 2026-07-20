"""The Foldseek `query` field must be normalised into the candidate id namespace.

build_structural_channel.py used the raw Foldseek `query` VERBATIM as the
channel join key while TARGET ids went through structural_evidence.target_keys()
/ _strip_structure_suffixes(). That asymmetry is the whole defect: the channel
TSV is left-joined onto the ranking table by bare candidate id, so any decoration
Foldseek adds to the query makes the row join to nothing -- and an all-miss join
is silent, because merge_evidence_channels() simply reads has_struct_data=0 for
every candidate and the structural voter goes dormant while every job exits 0.

GROUND TRUTH (this is the part the previous fixtures got wrong -- they encoded
the ASSUMED format). Measured on Unity, foldseek 10.941cd33, the binary in this
project's own `berghia-gpcr` env, on a compute node (srun job 61999969). Three
real RCSB structures were copied into a flat query dir under this repo's staging
convention (<candidate_id>.<ext>, as scripts/unity/run_foldseek_candidates.sh
symlinks them) and searched against themselves with
`--format-output query,target,fident,alntmscore,evalue`:

    query dir                          -> foldseek query/target id
    BersteEVm000001t1.cif   (1UBQ, 1 chain)  -> BersteEVm000001t1
    BersteEVm000002t1.pdb   (1UBQ, 1 chain)  -> BersteEVm000002t1
    BersteEVm000003t1.cif   (1F88, 2 chains) -> BersteEVm000003t1_A
                                                BersteEVm000003t1_B

So, concretely:
  1. Foldseek STRIPS the structure extension. It does NOT emit "<id>.cif", and
     therefore never emits "<id>.cif_A" either -- the form the audit predicted
     does not exist. (structcreatedb.cpp calls Util::remove_extension, twice for
     .gz/.zst.)
  2. Foldseek APPENDS "_<chain>" when the file contains more than one chain
     (CHAIN_MODE_AUTO), and emits one query row PER CHAIN. This is the real
     break: a multi-chain AF3 model (receptor + peptide/ligand/G-protein chain,
     or any complex) yields <cand_id>_A / <cand_id>_B, neither of which joins to
     the bare candidate id, and one candidate silently becomes two channel rows.

Note that target_keys() ALONE does not repair this: its _PDB_RE requires a
4-character id beginning with a digit, so "BersteEVm000003t1_A" matches nothing
and normalises to itself. The chain suffix has to be handled explicitly, and --
because a candidate id may legitimately contain underscores -- it may only be
stripped when the result is confirmed against the real candidate id universe.
Hence the join assertion these tests require.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

# conftest.py adds scripts/ to sys.path
import build_structural_channel as bsc
import structural_evidence as se

ANCHOR_FIELDS = ["accession", "tier", "taxid", "species", "family", "class", "evidence"]


def _write_foldseek(tmp_path: Path, name: str, lines: list) -> str:
    p = tmp_path / name
    p.write_text("\n".join(lines) + ("\n" if lines else ""))
    return str(p)


def _write_anchor_set(tmp_path: Path, rows: list) -> str:
    p = tmp_path / "anchor_set.tsv"
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ANCHOR_FIELDS, delimiter="\t")
        w.writeheader()
        for acc, fam in rows:
            w.writerow({"accession": acc, "tier": "1", "taxid": "6637",
                        "species": "Todarodes pacificus", "family": fam,
                        "class": "A", "evidence": "reviewed"})
    return str(p)


# The candidate id universe: bare ids exactly as they key the ranking table.
CANDIDATE_IDS = {"BersteEVm000001t1", "BersteEVm000002t1", "BersteEVm000003t1"}


# --------------------------------------------------------------------------- #
# query normalisation
# --------------------------------------------------------------------------- #

def test_bare_query_id_is_unchanged():
    """The single-chain case foldseek actually produces today."""
    assert bsc.canonical_query_id("BersteEVm000001t1", CANDIDATE_IDS) == "BersteEVm000001t1"


def test_chain_suffixed_query_resolves_to_the_candidate_id():
    """The measured multi-chain form: <cand_id>_A / <cand_id>_B -> <cand_id>."""
    assert bsc.canonical_query_id("BersteEVm000003t1_A", CANDIDATE_IDS) == "BersteEVm000003t1"
    assert bsc.canonical_query_id("BersteEVm000003t1_B", CANDIDATE_IDS) == "BersteEVm000003t1"


def test_extension_suffixed_query_resolves():
    """Defensive: a hits file from a foldseek build/config that kept the
    extension must still join. Uses the SAME helper the target side uses."""
    assert bsc.canonical_query_id("BersteEVm000001t1.cif", CANDIDATE_IDS) == "BersteEVm000001t1"
    assert bsc.canonical_query_id("BersteEVm000001t1.pdb", CANDIDATE_IDS) == "BersteEVm000001t1"
    assert bsc.canonical_query_id("BersteEVm000001t1.cif.gz", CANDIDATE_IDS) == "BersteEVm000001t1"


def test_query_normalisation_reuses_the_target_side_helper():
    """Symmetry is the point: the two sides must not drift apart."""
    for raw in ("BersteEVm000001t1.cif", "BersteEVm000001t1.cif.gz",
                "BersteEVm000001t1.pdb"):
        assert se._strip_structure_suffixes(raw) in bsc.query_keys(raw)


def test_path_prefixed_query_resolves():
    assert bsc.canonical_query_id(
        "/scratch/tmp/query_structures/BersteEVm000001t1.cif", CANDIDATE_IDS
    ) == "BersteEVm000001t1"


def test_unresolvable_query_returns_none_rather_than_a_guess():
    """A query that is in no sense a known candidate must NOT be silently
    rewritten into something that looks joinable."""
    assert bsc.canonical_query_id("SomeOtherThing_A", CANDIDATE_IDS) is None
    assert bsc.canonical_query_id("", CANDIDATE_IDS) is None


def test_underscore_bearing_candidate_id_is_not_truncated():
    """A candidate id ending in _<alnum> is a real id, not a chain suffix.

    This is why the chain strip may only be applied when corroborated by the
    candidate universe: blind regex stripping would map the real candidate
    'Berghia_scaffold_12' onto a nonexistent 'Berghia_scaffold'.
    """
    ids = {"Berghia_scaffold_12", "Berghia_scaffold"}
    assert bsc.canonical_query_id("Berghia_scaffold_12", ids) == "Berghia_scaffold_12"


def test_without_a_candidate_universe_the_ambiguous_chain_strip_is_not_applied():
    """No universe -> only the unambiguous extension strip. Never guess."""
    assert bsc.canonical_query_id("BersteEVm000001t1.cif", None) == "BersteEVm000001t1"
    assert bsc.canonical_query_id("BersteEVm000003t1_A", None) == "BersteEVm000003t1_A"


# --------------------------------------------------------------------------- #
# end-to-end: the channel joins, and multi-chain collapses to ONE row
# --------------------------------------------------------------------------- #

def test_multichain_queries_collapse_to_one_row_keeping_the_best_hit(tmp_path):
    """Two chains of one AF3 model must not become two candidate rows."""
    hits = _write_foldseek(tmp_path, "PDB.tsv", [
        "BersteEVm000003t1_A\tAF-P31356-F1-model_v4\t0.30\t0.40\t1e-5",
        "BersteEVm000003t1_B\tAF-P31356-F1-model_v4\t0.85\t0.92\t1e-30",
    ])
    fam = bsc.build_family_map(_write_anchor_set(tmp_path, [("P31356", "opsin")]))

    channel = bsc.build_structural_channel([hits], fam, candidate_ids=CANDIDATE_IDS)

    assert set(channel) == {"BersteEVm000003t1"}, channel
    # The better chain (alntmscore 0.92 > 0.40) wins, so the exclusion signal
    # survives instead of being diluted by the weaker chain.
    assert channel["BersteEVm000003t1"]["struct_state"] == "known_non_chemoreceptor"


def test_channel_ids_are_all_in_the_candidate_namespace(tmp_path):
    hits = _write_foldseek(tmp_path, "PDB.tsv", [
        "BersteEVm000001t1\tAF-P31356-F1-model_v4\t0.85\t0.92\t1e-30",
        "BersteEVm000003t1_A\tAF-P31356-F1-model_v4\t0.85\t0.92\t1e-30",
    ])
    fam = bsc.build_family_map(_write_anchor_set(tmp_path, [("P31356", "opsin")]))

    channel = bsc.build_structural_channel([hits], fam, candidate_ids=CANDIDATE_IDS)

    assert set(channel) <= CANDIDATE_IDS, set(channel) - CANDIDATE_IDS


# --------------------------------------------------------------------------- #
# the cross-cutting invariant: assert non-zero, in-range key overlap
# --------------------------------------------------------------------------- #

def test_zero_overlap_raises_instead_of_writing_an_orphan_channel(tmp_path):
    """Every query resolving to nothing is the exact silent-dormancy failure.

    Had the AF3-job-name bug (bead 5ubd) recurred, or had the query dir been
    staged under the wrong id scheme, this is the assertion that catches it
    instead of emitting a channel that left-joins to zero rows.
    """
    hits = _write_foldseek(tmp_path, "PDB.tsv", [
        "af3_job_00017_model\tAF-P31356-F1-model_v4\t0.85\t0.92\t1e-30",
        "af3_job_00042_model\tAF-P31356-F1-model_v4\t0.85\t0.92\t1e-30",
    ])
    fam = bsc.build_family_map(_write_anchor_set(tmp_path, [("P31356", "opsin")]))

    with pytest.raises(ValueError) as excinfo:
        bsc.build_structural_channel([hits], fam, candidate_ids=CANDIDATE_IDS)

    msg = str(excinfo.value)
    assert "af3_job_00017_model" in msg, msg          # what we saw
    assert "BersteEVm000001t1" in msg or "BersteEVm" in msg, msg  # what we expected


def test_partial_overlap_is_allowed_but_reported(tmp_path, capsys):
    """A partial foldseek run is legitimate (some candidates have no model);
    it must NOT raise, but the resolved/unresolved counts must be visible."""
    hits = _write_foldseek(tmp_path, "PDB.tsv", [
        "BersteEVm000001t1\tAF-P31356-F1-model_v4\t0.85\t0.92\t1e-30",
        "not_a_candidate\tAF-P31356-F1-model_v4\t0.85\t0.92\t1e-30",
    ])
    fam = bsc.build_family_map(_write_anchor_set(tmp_path, [("P31356", "opsin")]))

    channel = bsc.build_structural_channel([hits], fam, candidate_ids=CANDIDATE_IDS)

    assert set(channel) == {"BersteEVm000001t1"}
    err = capsys.readouterr().err
    assert "not_a_candidate" in err or "unresolved" in err.lower(), err


def test_no_hits_at_all_does_not_raise(tmp_path):
    """Zero hits is a different condition from zero OVERLAP; the guard must not
    turn an empty (but honest) foldseek result into a crash."""
    empty = _write_foldseek(tmp_path, "PDB.tsv", [])
    assert bsc.build_structural_channel([empty], {}, candidate_ids=CANDIDATE_IDS) == {}


def test_empty_candidate_universe_is_rejected(tmp_path):
    """An empty id set means the caller's own candidate source was empty --
    silently joining nothing would look like 'no structural evidence'."""
    hits = _write_foldseek(tmp_path, "PDB.tsv", [
        "BersteEVm000001t1\tAF-P31356-F1-model_v4\t0.85\t0.92\t1e-30",
    ])
    with pytest.raises(ValueError):
        bsc.build_structural_channel([hits], {}, candidate_ids=set())


# --------------------------------------------------------------------------- #
# backwards compatibility: candidate_ids is optional
# --------------------------------------------------------------------------- #

def test_without_candidate_ids_behaviour_is_the_previous_one(tmp_path):
    hits = _write_foldseek(tmp_path, "PDB.tsv", [
        "cand_nonchemo\tAF-P31356-F1-model_v4\t0.85\t0.92\t1e-30",
    ])
    fam = bsc.build_family_map(_write_anchor_set(tmp_path, [("P31356", "opsin")]))

    channel = bsc.build_structural_channel([hits], fam)

    assert channel["cand_nonchemo"]["struct_state"] == "known_non_chemoreceptor"


# --------------------------------------------------------------------------- #
# CLI: the candidate universe must be reachable from the command line
# --------------------------------------------------------------------------- #

def test_cli_accepts_a_candidate_fasta(tmp_path):
    fasta = tmp_path / "candidates.fa"
    fasta.write_text(
        ">BersteEVm000003t1 some description here\nMEEPMEEP\n"
        ">BersteEVm000001t1\nMEEP\n"
    )
    hits = _write_foldseek(tmp_path, "PDB.tsv", [
        "BersteEVm000003t1_A\tAF-P31356-F1-model_v4\t0.85\t0.92\t1e-30",
    ])
    anchors = _write_anchor_set(tmp_path, [("P31356", "opsin")])
    out = tmp_path / "structural_channel.tsv"

    rc = bsc.main([
        "--foldseek-tsvs", hits, "--anchor-set", anchors,
        "--candidate-fasta", str(fasta), "--out", str(out),
    ])

    assert rc == 0
    rows = list(csv.DictReader(open(out), delimiter="\t"))
    assert [r["id"] for r in rows] == ["BersteEVm000003t1"], rows


def test_read_candidate_ids_takes_the_first_whitespace_token(tmp_path):
    """FASTA headers carry descriptions; the id is the first token only."""
    fasta = tmp_path / "c.fa"
    fasta.write_text(">BersteEVm000001t1 len=350 class=A\nMEEP\n")
    assert bsc.read_candidate_ids(str(fasta)) == {"BersteEVm000001t1"}


# --------------------------------------------------------------------------- #
# stage 07 wiring
# --------------------------------------------------------------------------- #

STAGE_07 = Path(__file__).resolve().parents[2] / "07_candidate_ranking.sh"


def test_stage07_passes_the_candidate_universe_to_the_producer():
    text = STAGE_07.read_text()
    assert "--candidate-fasta" in text and "EMB_CANDIDATE_FASTA" in text, (
        "stage 07 does not hand build_structural_channel.py the candidate id "
        "universe, so query ids cannot be resolved and the join is unverified"
    )


def test_stage07_degrades_rather_than_voiding_when_the_candidate_fasta_is_absent():
    """An absent optional input must cost the join ASSERTION, not the channel."""
    text = STAGE_07.read_text()
    # Anchor on the INVOCATION, not the several comment mentions of the script.
    start = text.index('"${SCRIPTS_DIR}/build_structural_channel.py"')
    block = text[max(0, start - 2000):start]
    assert '[ -f "${EMB_CANDIDATE_FASTA}" ]' in block, block[-800:]
    assert "log --level=WARN" in block
