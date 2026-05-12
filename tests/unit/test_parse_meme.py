"""Unit tests for scripts/parse_meme.py (bead -urk + bead -7cy dual-mode)."""
import json
from pathlib import Path

import pytest

from parse_meme import OG_FIELDS, SITE_FIELDS, aggregate, extract_sites, parse_one


def _meme_json(rows):
    """Build a minimal MEME-style JSON dict.
    rows is a list of [alpha, beta_plus, prop_plus, p_value, n_branches] sublists.
    """
    return {
        "MLE": {
            "headers": [
                ["alpha"], ["beta-"], ["p-"], ["beta+"], ["p+"],
                ["LRT"], ["p-value"], ["# branches"], ["post-mean"]
            ],
            "content": {
                "0": [
                    [r[0], 0.0, 0.0, r[1], r[2], 0.0, r[3], r[4], 0.0]
                    for r in rows
                ]
            }
        }
    }


class TestExtractSites:
    def test_basic_three_sites(self):
        d = _meme_json([
            [0.5, 5.0, 0.1, 0.001, 3],
            [0.5, 1.0, 0.0, 0.5,   0],
            [0.5, 8.0, 0.2, 0.04,  2],
        ])
        sites = extract_sites(d)
        assert len(sites) == 3
        assert sites[0]["beta_plus"] == 5.0
        assert sites[0]["p_value"] == pytest.approx(0.001)
        assert sites[2]["n_branches_episodic"] == 2

    def test_empty(self):
        assert extract_sites({}) == []


class TestAggregate:
    def test_counts_episodic_sites_correctly(self):
        sites = [
            {"p_value": 0.001, "beta_plus": 10.0},
            {"p_value": 0.04,  "beta_plus": 5.0},
            {"p_value": 0.5,   "beta_plus": 1.0},
            {"p_value": 0.02,  "beta_plus": 0.5},
        ]
        out = aggregate(sites, alpha=0.05)
        assert out["n_sites_total"] == 4
        assert out["n_sites_episodic"] == 2
        assert out["fraction_episodic"] == pytest.approx(0.5)

    def test_no_episodic(self):
        sites = [{"p_value": 0.5, "beta_plus": 0.5}]
        out = aggregate(sites, alpha=0.05)
        assert out["n_sites_episodic"] == 0
        assert out["fraction_episodic"] == 0.0


class TestParseOne:
    def test_schema_fields(self, tmp_path):
        f = tmp_path / "OGxx_meme.json"
        f.write_text(json.dumps(_meme_json([[0.5, 8.0, 0.1, 0.01, 2]])))
        sites, og = parse_one(str(f), "OGxx", alpha=0.05)
        assert og["og_name"] == "OGxx"
        for k in OG_FIELDS:
            assert k in og, k
        for k in ("og_name", "site", "is_episodic"):
            assert k in sites[0], k

    def test_episodic_flag_per_site(self, tmp_path):
        f = tmp_path / "x_meme.json"
        f.write_text(json.dumps(_meme_json([
            [0.5, 8.0, 0.1, 0.001, 2],
            [0.5, 0.5, 0.0, 0.001, 0],
        ])))
        sites, og = parse_one(str(f), "OGyy", alpha=0.05)
        assert sites[0]["is_episodic"] == 1
        assert sites[1]["is_episodic"] == 0
        assert og["n_sites_episodic"] == 1


# ---- bead -7cy: dual-mode concordance stratification --------------------

from parse_meme import (
    map_strict_to_lenient,
    parse_one_dual,
    OG_FIELDS_DUAL,
    SITE_FIELDS_DUAL,
)


def _write_fa(path: Path, seqs: dict[str, str]) -> None:
    """Write a minimal FASTA file from a dict of {name: sequence}."""
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")


class TestColumnVectorMapping:
    """The trimmed-to-trimmed column mapping is the load-bearing primitive
    for concordance: strict's columns are a subset of lenient's columns,
    and we identify which-is-which by residue-vector equality across all
    sequences in the OG. Tests cover the typical case + edge cases that
    could mis-align the mapping."""

    def test_strict_is_subset_of_lenient(self, tmp_path):
        """Strict drops columns 2 and 5 of the (implicit) original. Lenient
        keeps all 6 columns. Each strict column must map to its lenient
        counterpart by residue-vector match."""
        strict = tmp_path / "s.fa"
        lenient = tmp_path / "l.fa"
        # 4 sequences, columns chosen so each column-vector is unique
        _write_fa(strict, {
            # cols:    1   2   3   4    (kept from original 1,3,4,6)
            "A":      "M" "L" "Y" "K",
            "B":      "M" "F" "Y" "R",
            "C":      "M" "L" "W" "K",
            "Berghia":"M" "I" "Y" "K",
        })
        _write_fa(lenient, {
            # cols:    1   2   3   4   5   6   (original kept entirely)
            "A":      "M" "A" "L" "Y" "G" "K",
            "B":      "M" "S" "F" "Y" "G" "R",
            "C":      "M" "A" "L" "W" "T" "K",
            "Berghia":"M" "A" "I" "Y" "G" "K",
        })
        mapping = map_strict_to_lenient(str(strict), str(lenient))
        # Strict's 4 columns must map to lenient's columns 1, 3, 4, 6
        assert mapping == [1, 3, 4, 6], (
            f"expected strict→lenient mapping [1,3,4,6], got {mapping}"
        )

    def test_identical_alignments_map_one_to_one(self, tmp_path):
        """If strict == lenient (nothing was trimmed in either pass), every
        column maps to itself."""
        strict = tmp_path / "s.fa"
        lenient = tmp_path / "l.fa"
        seqs = {"A": "MKLY", "B": "MFLY", "C": "MKWY"}
        _write_fa(strict, seqs)
        _write_fa(lenient, seqs)
        assert map_strict_to_lenient(str(strict), str(lenient)) == [1, 2, 3, 4]

    def test_strict_empty_returns_empty_mapping(self, tmp_path):
        """Degenerate: strict trimmed to zero columns. Mapping is []."""
        strict = tmp_path / "s.fa"
        lenient = tmp_path / "l.fa"
        _write_fa(strict, {"A": "", "B": ""})
        _write_fa(lenient, {"A": "MK", "B": "MF"})
        assert map_strict_to_lenient(str(strict), str(lenient)) == []


class TestConcordanceStratification:
    """OG-level concordance counts derived from strict/lenient MEME positives
    cross-referenced via the column-vector mapping."""

    def test_concordance_tiers_counted_correctly(self, tmp_path):
        """Strict drops cols 2 & 5; lenient keeps all 6. MEME-positives:
            strict sites {1, 3}     → in lenient coords (via map): {1, 4}
            lenient sites {1, 2, 4} → as-is in lenient coords:     {1, 2, 4}
        Expected:
            high_confidence = {1, 4}     (sites positive in both passes)
            lenient_only    = {2}        (only positive in lenient)
            strict_only     = {}         (none — both strict positives also positive in lenient)
            alignment_robustness_index = high_confidence_n / lenient_n
                                       = 2 / 3 ≈ 0.6667
        """
        strict_fa = tmp_path / "trimmed.fa"
        lenient_fa = tmp_path / "trimmed_lenient.fa"
        _write_fa(strict_fa, {
            "A":      "M" "L" "Y" "K",   # cols 1,3,4,6 from original
            "B":      "M" "F" "Y" "R",
            "C":      "M" "L" "W" "K",
            "Berghia":"M" "I" "Y" "K",
        })
        _write_fa(lenient_fa, {
            "A":      "M" "A" "L" "Y" "G" "K",
            "B":      "M" "S" "F" "Y" "G" "R",
            "C":      "M" "A" "L" "W" "T" "K",
            "Berghia":"M" "A" "I" "Y" "G" "K",
        })

        # Strict MEME: sites 1 and 3 are episodic (p<0.05, beta+ > 1)
        strict_json = tmp_path / "strict_meme.json"
        strict_json.write_text(json.dumps(_meme_json([
            [0.5, 8.0, 0.1, 0.001, 2],   # site 1: positive
            [0.5, 0.5, 0.0, 0.50,  0],   # site 2: negative
            [0.5, 6.0, 0.1, 0.01,  3],   # site 3: positive
            [0.5, 1.0, 0.0, 0.30,  0],   # site 4: negative
        ])))

        # Lenient MEME: sites 1, 2, 4 are episodic
        lenient_json = tmp_path / "lenient_meme.json"
        lenient_json.write_text(json.dumps(_meme_json([
            [0.5, 8.0, 0.1, 0.001, 2],   # lenient site 1 (→ strict site 1 via map): positive
            [0.5, 5.0, 0.1, 0.02,  2],   # lenient site 2: positive (NOT in strict)
            [0.5, 0.5, 0.0, 0.50,  0],   # lenient site 3: negative (→ strict site 2)
            [0.5, 7.0, 0.1, 0.005, 2],   # lenient site 4 (→ strict site 3 via map): positive
            [0.5, 0.5, 0.0, 0.50,  0],   # lenient site 5: negative
            [0.5, 0.5, 0.0, 0.30,  0],   # lenient site 6 (→ strict site 4): negative
        ])))

        sites, og = parse_one_dual(
            strict_json=str(strict_json), strict_fa=str(strict_fa),
            lenient_json=str(lenient_json), lenient_fa=str(lenient_fa),
            og_name="OGtest", alpha=0.05,
        )

        # OG-level counts
        assert og["og_name"] == "OGtest"
        assert og["n_strict_positive_sites"] == 2
        assert og["n_lenient_positive_sites"] == 3
        assert og["high_confidence_sites_n"] == 2, (
            f"expected high_confidence=2, got {og['high_confidence_sites_n']}"
        )
        assert og["lenient_only_sites_n"] == 1, (
            f"expected lenient_only=1, got {og['lenient_only_sites_n']}"
        )
        assert og["strict_only_sites_n"] == 0, (
            f"expected strict_only=0, got {og['strict_only_sites_n']}"
        )
        assert og["alignment_robustness_index"] == pytest.approx(2.0 / 3.0, abs=1e-4)

    def test_strict_only_sites_when_strict_finds_extra(self, tmp_path):
        """Pathological case: strict declares a site positive that lenient
        doesn't. (Rare in practice — usually lenient is more permissive —
        but the stratifier must handle it.)"""
        strict_fa = tmp_path / "trimmed.fa"
        lenient_fa = tmp_path / "trimmed_lenient.fa"
        # Identical alignments (no actual trim difference) so mapping is 1:1
        seqs = {"A": "MKL", "B": "MFL", "C": "MKW"}
        _write_fa(strict_fa, seqs)
        _write_fa(lenient_fa, seqs)

        strict_json = tmp_path / "strict.json"
        strict_json.write_text(json.dumps(_meme_json([
            [0.5, 9.0, 0.1, 0.001, 2],   # positive
            [0.5, 0.5, 0.0, 0.5,   0],
            [0.5, 0.5, 0.0, 0.5,   0],
        ])))
        lenient_json = tmp_path / "lenient.json"
        lenient_json.write_text(json.dumps(_meme_json([
            [0.5, 0.5, 0.0, 0.5,   0],   # NOT positive in lenient
            [0.5, 0.5, 0.0, 0.5,   0],
            [0.5, 0.5, 0.0, 0.5,   0],
        ])))
        _, og = parse_one_dual(
            strict_json=str(strict_json), strict_fa=str(strict_fa),
            lenient_json=str(lenient_json), lenient_fa=str(lenient_fa),
            og_name="OGedge", alpha=0.05,
        )
        assert og["high_confidence_sites_n"] == 0
        assert og["lenient_only_sites_n"] == 0
        assert og["strict_only_sites_n"] == 1
        # robustness = high_conf / lenient — with lenient=0, define as 0.0 (not NaN/error)
        assert og["alignment_robustness_index"] == pytest.approx(0.0)

    def test_dual_mode_per_site_csv_carries_concordance_tier(self, tmp_path):
        """Per-site CSV: each MEME-positive site under either pass gets one
        row, with concordance_tier ∈ {high_confidence, lenient_only,
        strict_only}. Negatives in both are excluded (per design — we don't
        report non-positive sites)."""
        strict_fa = tmp_path / "trimmed.fa"
        lenient_fa = tmp_path / "trimmed_lenient.fa"
        seqs = {"A": "MK", "B": "MF"}
        _write_fa(strict_fa, seqs)
        _write_fa(lenient_fa, seqs)
        strict_json = tmp_path / "strict.json"
        strict_json.write_text(json.dumps(_meme_json([
            [0.5, 9.0, 0.1, 0.001, 2],   # site 1: positive (concordant)
            [0.5, 0.5, 0.0, 0.5,   0],   # site 2: negative
        ])))
        lenient_json = tmp_path / "lenient.json"
        lenient_json.write_text(json.dumps(_meme_json([
            [0.5, 9.0, 0.1, 0.001, 2],   # site 1: positive (concordant)
            [0.5, 6.0, 0.1, 0.01,  2],   # site 2: positive (lenient-only)
        ])))
        sites, _ = parse_one_dual(
            strict_json=str(strict_json), strict_fa=str(strict_fa),
            lenient_json=str(lenient_json), lenient_fa=str(lenient_fa),
            og_name="OGsites", alpha=0.05,
        )
        tiers = {(s["lenient_site"], s["concordance_tier"]) for s in sites}
        assert tiers == {(1, "high_confidence"), (2, "lenient_only")}
        # All per-site rows must have the required dual-mode fields
        for s in sites:
            for k in SITE_FIELDS_DUAL:
                assert k in s, f"per-site row missing field: {k}"


class TestCodonAwareMapping:
    """Codon-aware mode groups input FASTA columns in triplets so column-
    vector matching operates at codon granularity. Critical for stage 05's
    use case: ClipKit -co trims codons; MEME indexes by codon site; the
    mapping must be codon-to-codon, not nucleotide-to-nucleotide."""

    def test_codon_aware_groups_triplets(self, tmp_path):
        """Strict drops codon 2 of an input with 3 codons. Lenient keeps
        all 3. Each codon = 3 nucleotides; the column-vector at codon i
        is the tuple of 3-char codons at that position across all seqs.
        Strict→lenient map under codon_aware=True: [1, 3]."""
        strict = tmp_path / "s.fa"
        lenient = tmp_path / "l.fa"
        _write_fa(strict, {
            # codons:  1     2 (dropped)  3
            "A":      "AAA"               "CCC",  # 6 nt = 2 codons
            "B":      "AAG"               "CCT",
            "C":      "AAA"               "CCC",
        })
        _write_fa(lenient, {
            # codons:  1     2     3
            "A":      "AAA" "GGG" "CCC",          # 9 nt = 3 codons
            "B":      "AAG" "GTA" "CCT",
            "C":      "AAA" "GGG" "CCC",
        })
        mapping = map_strict_to_lenient(str(strict), str(lenient),
                                         codon_aware=True)
        assert mapping == [1, 3], (
            f"expected codon mapping [1,3], got {mapping}"
        )

    def test_codon_aware_concordance_via_parse_one_dual(self, tmp_path):
        """End-to-end: strict-trimmed codon (2 codons) + lenient codon
        (3 codons), MEME positives in strict={1, 2} (both kept codons)
        and lenient={1, 2, 3}. Expected:
            strict positives in lenient coords: {1, 3}
            lenient positives: {1, 2, 3}
            high_confidence = {1, 3}      (concordant across passes)
            lenient_only    = {2}         (codon dropped by strict)
            strict_only     = {}
        """
        strict_fa = tmp_path / "codon_strict.fa"
        lenient_fa = tmp_path / "codon_lenient.fa"
        _write_fa(strict_fa, {
            "A": "AAACCC", "B": "AAGCCT", "C": "AAACCC",
        })
        _write_fa(lenient_fa, {
            "A": "AAAGGGCCC", "B": "AAGGTACCT", "C": "AAAGGGCCC",
        })
        strict_json = tmp_path / "strict.json"
        strict_json.write_text(json.dumps(_meme_json([
            [0.5, 9.0, 0.1, 0.001, 2],   # site 1 (codon 1): positive
            [0.5, 7.0, 0.1, 0.005, 2],   # site 2 (codon 2 in strict = codon 3 in lenient): positive
        ])))
        lenient_json = tmp_path / "lenient.json"
        lenient_json.write_text(json.dumps(_meme_json([
            [0.5, 9.0, 0.1, 0.001, 2],   # site 1: positive
            [0.5, 6.0, 0.1, 0.02,  2],   # site 2 (only in lenient): positive
            [0.5, 7.0, 0.1, 0.005, 2],   # site 3: positive
        ])))
        from parse_meme import parse_one_dual
        sites, og = parse_one_dual(
            strict_json=str(strict_json), strict_fa=str(strict_fa),
            lenient_json=str(lenient_json), lenient_fa=str(lenient_fa),
            og_name="OGcodon", alpha=0.05, codon_aware=True,
        )
        assert og["high_confidence_sites_n"] == 2
        assert og["lenient_only_sites_n"] == 1
        assert og["strict_only_sites_n"] == 0
        assert og["alignment_robustness_index"] == pytest.approx(2 / 3, abs=1e-4)


class TestDualModeSchema:
    def test_OG_FIELDS_DUAL_contains_concordance_columns(self):
        for required in ("og_name", "n_strict_positive_sites",
                         "n_lenient_positive_sites", "high_confidence_sites_n",
                         "lenient_only_sites_n", "strict_only_sites_n",
                         "alignment_robustness_index"):
            assert required in OG_FIELDS_DUAL, required

    def test_SITE_FIELDS_DUAL_contains_concordance_tier(self):
        for required in ("og_name", "lenient_site", "strict_site",
                         "strict_pvalue", "lenient_pvalue",
                         "concordance_tier"):
            assert required in SITE_FIELDS_DUAL, required
