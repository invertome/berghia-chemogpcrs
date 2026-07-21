"""Tests that the chemoreceptor probe can never reach the reference side.

The probe set holds TARGET-class sequences: published molluscan chemosensory
receptors. The reference set holds things a candidate is NOT. Novelty is a MIN
over family prototypes, so a single chemoreceptor admitted as a reference
becomes the prototype nearest to every chemoreceptor-like candidate and
suppresses exactly the signal the pipeline hunts -- silently, with no error and
no count changing.

These tests therefore assert the separation on the REAL on-disk artifacts, not
only on fixtures. A fixture encodes the separation it was written to assume; it
cannot catch the case where the production files actually overlap.
"""
import csv
from pathlib import Path

import pytest

from refexp5_build_chemoreceptor_probe import (
    IDENTITY_CLUSTERS,
    IDENTITY_MEDIAN_NEAREST_NEIGHBOUR,
    IDENTITY_MEDIAN_PAIRWISE,
    KNOWN_PROBE_TWINS,
    N_EFF_BRACKET,
    PROBE_ID_PREFIX,
    REFERENCE_ID_PREFIX,
    assert_no_probe_in_reference,
    assert_no_probe_sequence_in_reference,
    effective_exclusion_power,
    effective_sample_size,
    exclusion_power,
    probe_id,
    verify_separation,
)

REPO = Path(__file__).resolve().parents[2]
PROBE_TSV = REPO / "references/anchors/chemoreceptor_probe_set.tsv"
PROBE_FASTA = REPO / "references/anchors/chemoreceptor_probe_set.fasta"
ANCHORS_TSV = REPO / "references/anchors/anchor_set_PROD.tsv"

# references/ is gitignored, so a clean checkout has no artifacts to check.
# Skipping is correct there; silently passing would not be.
needs_artifacts = pytest.mark.skipif(
    not (PROBE_TSV.exists() and ANCHORS_TSV.exists()),
    reason="probe set / anchor table absent (references/ is gitignored)",
)


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# --- the guard itself --------------------------------------------------------

def test_guard_passes_on_disjoint_sets():
    assert_no_probe_in_reference(
        ["ANCHOR_A_10_P12345", "ANCHOR_A_1_Q99999"],
        ["PROBE_A_C5H877", "PROBE_A_C5H675"],
    )


def test_guard_catches_a_probe_copied_into_the_reference():
    with pytest.raises(AssertionError, match="reached the"):
        assert_no_probe_in_reference(
            ["ANCHOR_A_10_P12345", "PROBE_A_C5H877"],
            ["PROBE_A_C5H877"],
        )


def test_guard_catches_a_probe_smuggled_in_under_an_unlisted_id():
    """The prefix check is not redundant with the overlap check.

    A probe re-keyed to an id the caller did not pass in `probe_ids` defeats the
    overlap check entirely; only the prefix check sees it.
    """
    with pytest.raises(AssertionError, match="probe prefix"):
        assert_no_probe_in_reference(
            ["ANCHOR_A_10_P12345", "PROBE_A_SOMETHING_ELSE"],
            ["PROBE_A_C5H877"],
        )


def test_probe_id_cannot_collide_with_an_anchor_composite():
    """Even for the same accession, the two id schemes are disjoint."""
    acc = "C5H877"
    assert probe_id(acc) != f"{REFERENCE_ID_PREFIX}A_10_{acc}"
    assert probe_id(acc).startswith(PROBE_ID_PREFIX)
    assert not probe_id(acc).startswith(REFERENCE_ID_PREFIX)


# --- the real artifacts ------------------------------------------------------

@needs_artifacts
def test_no_probe_id_appears_in_the_real_anchor_table():
    """The load-bearing assertion, on production data."""
    assert not verify_separation(PROBE_TSV, ANCHORS_TSV)


@needs_artifacts
def test_no_probe_accession_is_a_real_anchor():
    probe_acc = {r["accession"] for r in read_tsv(PROBE_TSV)}
    anchor_acc = {r["accession"] for r in read_tsv(ANCHORS_TSV)}
    assert probe_acc, "probe set is empty"
    assert not probe_acc & anchor_acc


@needs_artifacts
def test_probes_are_absent_from_the_novelty_prototype_labels():
    """Check the ACTUAL producer, not a re-implementation of its rule.

    `load_ref_labels` mints the keys that become prototype members, and
    `novelty_reference_labels` restricts them to the characterized class-A
    families that define "known" space. Neither may ever yield a probe.
    """
    from build_embedding_channel import load_ref_labels, novelty_reference_labels

    ref_labels = load_ref_labels(str(ANCHORS_TSV))
    assert ref_labels, "reference labels are empty"
    novelty_labels = novelty_reference_labels(ref_labels)
    assert novelty_labels, "novelty labels are empty"

    probe_ids = [r["probe_id"] for r in read_tsv(PROBE_TSV)]
    assert_no_probe_in_reference(ref_labels, probe_ids, context="ref labels")
    assert_no_probe_in_reference(novelty_labels, probe_ids,
                                 context="novelty prototype labels")


@needs_artifacts
def test_probe_tsv_has_no_family_column():
    """Layer 3 of the separation: pointing a reference loader here must FAIL.

    `load_ref_labels` reads df["family"], df["tier"] and df["class"]. Withholding
    `family` and `tier` turns a misdirected loader into a loud KeyError instead
    of a silent chemoreceptor prototype.
    """
    rows = read_tsv(PROBE_TSV)
    assert rows
    assert "family" not in rows[0]
    assert "tier" not in rows[0]


@needs_artifacts
def test_misdirecting_the_reference_loader_at_the_probe_set_raises():
    from build_embedding_channel import load_ref_labels

    with pytest.raises(KeyError):
        load_ref_labels(str(PROBE_TSV))


@needs_artifacts
def test_probe_sequence_lengths_match_the_recorded_lengths():
    seqs, name = {}, None
    with open(PROBE_FASTA) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.strip())
    seqs = {k: "".join(v) for k, v in seqs.items()}

    rows = read_tsv(PROBE_TSV)
    assert len(seqs) == len(rows)
    for r in rows:
        assert len(seqs[r["probe_id"]]) == int(r["sequence_length"])


# --- the duplicate-under-another-accession leak ------------------------------

def test_sequence_guard_passes_when_nothing_is_shared():
    assert_no_probe_sequence_in_reference(
        {"ANCHOR_A_10_P12345": "MKTAYIAK"}, {"PROBE_A_C5H877": "MGGWSSLL"})


def test_sequence_guard_catches_the_same_protein_under_another_accession():
    """The leak every key-based check passes through.

    Accessions differ, ids differ, counts are right -- and a chemoreceptor
    prototype sits at distance zero from the probe it duplicates.
    """
    with pytest.raises(AssertionError, match="byte-identical"):
        assert_no_probe_sequence_in_reference(
            {"ANCHOR_A_9_A0ABM0GH42": "MGGWSSLL"},
            {"PROBE_A_C5H675": "MGGWSSLL"})


def test_known_probe_twins_are_recorded():
    """Pinned 2026-07-21; all four really are byte-identical to a probe.

    The UniProt pair was found by sequence comparison during the reference
    audit. The RefSeq pair was found the same way while recovering the family:
    AplCal3.0 carries two of the three seeds under RefSeq accessions.
    """
    assert KNOWN_PROBE_TWINS == {"A0ABM0GH42": "C5H675",
                                 "A0ABM0GH43": "C5H674",
                                 "NP_001191494.1": "C5H675",
                                 "NP_001191495.1": "C5H674"}


@needs_artifacts
def test_no_known_probe_twin_is_in_the_real_anchor_set():
    anchor_acc = {r["accession"] for r in read_tsv(ANCHORS_TSV)}
    assert not set(KNOWN_PROBE_TWINS) & anchor_acc


@needs_artifacts
def test_no_real_anchor_sequence_duplicates_a_probe():
    """The load-bearing sequence check, on production data.

    Anchor sequences come from the UniProt dump plus the additions FASTA; the
    test asserts it actually resolved sequences before concluding anything, so
    an empty join cannot masquerade as a clean result.
    """
    uniprot_tsv = REPO / "references/anchors/anchor_set_PROD_uniprot.tsv"
    additions = REPO / "references/anchors/anchor_set_PROD_additions.fasta"
    if not (uniprot_tsv.exists() and additions.exists()):
        pytest.skip("anchor sequence sources absent")

    by_acc = {r["queried_accession"]: (r.get("Sequence") or "").strip()
              for r in read_tsv(uniprot_tsv)}
    add = {}
    name = None
    with open(additions) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                add[name] = []
            elif name:
                add[name].append(line.strip())
    add = {k: "".join(v) for k, v in add.items()}

    anchor_seqs = {}
    for r in read_tsv(ANCHORS_TSV):
        cid = f"ANCHOR_{r['class']}_{r['tier']}_{r['accession']}"
        seq = by_acc.get(r["accession"]) or add.get(cid, "")
        if seq:
            anchor_seqs[cid] = seq
    assert len(anchor_seqs) > 1000, (
        f"only {len(anchor_seqs)} anchor sequences resolved -- an empty join "
        "would make this test vacuously pass")

    probe_seqs, name = {}, None
    with open(PROBE_FASTA) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                probe_seqs[name] = []
            elif name:
                probe_seqs[name].append(line.strip())
    probe_seqs = {k: "".join(v) for k, v in probe_seqs.items()}
    assert probe_seqs

    assert_no_probe_sequence_in_reference(anchor_seqs, probe_seqs,
                                          context="anchor set")


# --- the honesty of the readout ---------------------------------------------

def test_exclusion_readout_was_underpowered_at_three_probes():
    """Pin why the set had to be expanded.

    At the exclusion rates this pipeline calls, a 3-probe set had under a 10%
    chance of producing even one exclusion if the probes behave like average
    candidates. "Not excluded" was therefore the expected outcome under BOTH
    hypotheses -- which is what made expansion necessary rather than cosmetic.
    """
    for rate in (0.010, 0.0177, 0.030):
        assert exclusion_power(3, rate) < 0.10


def test_exclusion_power_is_monotone_and_bounded():
    assert exclusion_power(3, 0.0) == 0.0
    assert exclusion_power(3, 1.0) == 1.0
    assert exclusion_power(1, 0.05) < exclusion_power(3, 0.05)


# --- correlation: the expanded set is NOT n independent probes ---------------

def test_effective_sample_size_collapses_as_correlation_rises():
    """The limits are the whole point.

    rho=0 -> n independent probes. rho=1 -> the set behaves as ONE probe,
    which is correct for a group of near-identical sequences: they either all
    trip the exclusion or none do.
    """
    assert effective_sample_size(155, 0.0) == 155
    assert effective_sample_size(155, 1.0) == pytest.approx(1.0)
    assert effective_sample_size(1, 0.5) == 1


def test_effective_sample_size_is_monotone_in_correlation():
    prev = None
    for rho in (0.0, 0.1, 0.25, 0.5, 0.9, 1.0):
        n_eff = effective_sample_size(159, rho)
        if prev is not None:
            assert n_eff < prev
        prev = n_eff


def test_effective_sample_size_rejects_impossible_correlation():
    with pytest.raises(ValueError):
        effective_sample_size(10, 1.5)
    with pytest.raises(ValueError):
        effective_sample_size(10, -0.1)


def test_correlated_power_never_exceeds_the_independent_bound():
    """The naive number is an upper bound and must behave like one."""
    for rate in (0.010, 0.0177, 0.030):
        for rho in (0.0, 0.2, 0.5, 0.9):
            assert (effective_exclusion_power(159, rate, rho)
                    <= exclusion_power(159, rate) + 1e-12)


def test_expansion_materially_improves_falsification_power():
    """The honest claim: at the LOW end of the n_eff bracket the expanded set
    still beats 3 probes by a wide margin. If this ever fails, the expansion
    stopped buying what it was built to buy."""
    lo, hi = N_EFF_BRACKET
    for rate in (0.010, 0.0177, 0.030):
        assert exclusion_power(lo, rate) > 4 * exclusion_power(3, rate)
        assert exclusion_power(hi, rate) > exclusion_power(lo, rate)


def test_n_eff_bracket_comes_from_the_measured_clusters():
    assert N_EFF_BRACKET == (IDENTITY_CLUSTERS[50], IDENTITY_CLUSTERS[90])
    assert IDENTITY_CLUSTERS[50] < IDENTITY_CLUSTERS[90] < 159


def test_cluster_counts_are_monotone_in_identity_threshold():
    """Looser threshold -> fewer, larger clusters. A violation means the
    numbers were transcribed rather than computed."""
    ths = sorted(IDENTITY_CLUSTERS)
    counts = [IDENTITY_CLUSTERS[t] for t in ths]
    assert counts == sorted(counts)


def test_nearest_neighbour_identity_exceeds_median_pairwise():
    """The two identity numbers must both be reported.

    Median pairwise (31.6%) alone would suggest a diffuse set; median
    nearest-neighbour (75.9%) reveals it is clustered into close relatives.
    Quoting only the first understates the redundancy.
    """
    assert IDENTITY_MEDIAN_NEAREST_NEIGHBOUR > IDENTITY_MEDIAN_PAIRWISE


# --- the expanded set: strata and separation ---------------------------------

@needs_artifacts
def test_expanded_probe_set_carries_all_three_strata():
    rows = read_tsv(PROBE_TSV)
    strata = {r["stratum"] for r in rows}
    assert strata == {"published-deposited", "published-independent",
                      "hmm-recovered"}
    assert sum(1 for r in rows if r["stratum"] == "published-deposited") == 3
    assert sum(1 for r in rows if r["stratum"] == "published-independent") == 1


@needs_artifacts
def test_every_expanded_member_still_carries_the_probe_prefix():
    """Layer 2 must hold for all 159, not only the original three."""
    for r in read_tsv(PROBE_TSV):
        assert r["probe_id"].startswith(PROBE_ID_PREFIX)
        assert not r["probe_id"].startswith(REFERENCE_ID_PREFIX)


@needs_artifacts
def test_expanded_probe_tsv_still_withholds_family_and_tier():
    """Layer 3 -- adding columns must not have added the two forbidden ones."""
    rows = read_tsv(PROBE_TSV)
    assert rows
    assert "family" not in rows[0]
    assert "tier" not in rows[0]


@needs_artifacts
def test_the_fragment_member_is_flagged():
    """C6FGJ8 is 246 aa against 354 for the seeds. Scoring it as though it were
    full-length would report partly on length rather than on identity."""
    rows = {r["accession"]: r for r in read_tsv(PROBE_TSV)}
    assert rows["C6FGJ8"]["is_fragment"] == "1"
    assert rows["C6FGJ8"]["stratum"] == "published-independent"
    assert int(rows["C6FGJ8"]["sequence_length"]) == 246
    # and it is the ONLY fragment, so a pooled statistic has exactly one
    # member whose scores are not length-comparable
    frags = [a for a, r in rows.items() if r["is_fragment"] == "1"]
    assert frags == ["C6FGJ8"]


@needs_artifacts
def test_no_expanded_probe_sequence_duplicates_another():
    """Sequence-level uniqueness across the whole expanded set."""
    seqs = {}
    name = None
    with open(PROBE_FASTA) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.strip())
    seqs = {k: "".join(v) for k, v in seqs.items()}
    assert len(seqs) == len(read_tsv(PROBE_TSV))
    assert len(set(seqs.values())) == len(seqs), "duplicate probe sequences"
