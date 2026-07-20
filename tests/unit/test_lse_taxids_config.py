"""Tests pinning the LSE clade taxids in config.sh and scripts/lse_refine.py.

Background: all three LSE clade taxids were wrong. 644 resolved to
*Aeromonas hydrophila* (a bacterium), 54397 to a *Lamellibrachia* sp.
endosymbiont (also a bacterium), and 13843 did not resolve at all. Because
`lse_refine.py` classifies via `if LSE_AEOLID_TAXID in common_lineage`, and a
bacterium's taxid can never occur in a mollusc's lineage, every orthogroup fell
through every branch: LSE classification silently produced nothing.

The corrected values are ancestors of *Berghia stephanieae* (taxid 1287507),
verified against the NCBI taxonomy API:

    Gastropoda    6448   (class)
    Nudibranchia  70849  (order)
    Aeolidioidea  71481  (superfamily)

These tests are deliberately offline: they pin the literal values recorded in
the two source files and that the two agree with each other. Live resolution of
those values against NCBI lives in test_lse_taxids_ncbi.py.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
CONFIG_SH = REPO_ROOT / "config.sh"
LSE_REFINE = REPO_ROOT / "scripts" / "lse_refine.py"

# The corrected, API-verified clade taxids. Every assertion below keys off this
# one table so the expected values live in exactly one place.
EXPECTED_TAXIDS = {
    "LSE_AEOLID_TAXID": 71481,      # Aeolidioidea (superfamily)
    "LSE_NUDIBRANCH_TAXID": 70849,  # Nudibranchia (order)
    "LSE_GASTROPOD_TAXID": 6448,    # Gastropoda (class)
}

# The wrong values that shipped previously. None may reappear in either file.
HISTORICAL_WRONG_TAXIDS = {
    54397: "Lamellibrachia sp. endosymbiont (bacterium)",
    13843: "does not resolve",
    644: "Aeromonas hydrophila (bacterium)",
    # The old config comment's "Example with real taxids" was wrong too.
    6524: "Planorbidae (freshwater snail family, not Nudibranchia)",
    1514845: "Bartonella sp. (bacterium)",
    285658: "Parasesarma (a crab genus)",
}


def _sourced_config_value(name: str) -> str:
    """Source config.sh in a subshell and echo one variable.

    config.sh is not `set -u`-safe, so it is sourced without nounset.
    """
    proc = subprocess.run(
        ["bash", "-c", f'source "{CONFIG_SH}" >/dev/null 2>&1; printf "%s" "${{{name}}}"'],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, f"sourcing config.sh failed: {proc.stderr}"
    return proc.stdout.strip()


def _sourced_lse_levels() -> list[str]:
    """Source config.sh and emit the LSE_LEVELS array, one element per line."""
    proc = subprocess.run(
        ["bash", "-c",
         f'source "{CONFIG_SH}" >/dev/null 2>&1; printf "%s\\n" "${{LSE_LEVELS[@]}}"'],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, f"sourcing config.sh failed: {proc.stderr}"
    return [line for line in proc.stdout.splitlines() if line.strip()]


def _lse_refine_default(name: str) -> int:
    """Extract the os.getenv default for a taxid constant in lse_refine.py."""
    text = LSE_REFINE.read_text()
    m = re.search(
        rf"^{re.escape(name)}\s*=\s*int\(\s*os\.getenv\(\s*['\"]{re.escape(name)}['\"]\s*,\s*(\d+)\s*\)",
        text,
        re.MULTILINE,
    )
    assert m, f"could not find an os.getenv default for {name} in {LSE_REFINE}"
    return int(m.group(1))


# --- config.sh: the three exported taxids ------------------------------------

@pytest.mark.parametrize("name,expected", sorted(EXPECTED_TAXIDS.items()))
def test_config_exports_corrected_taxid(name, expected):
    assert _sourced_config_value(name) == str(expected)


@pytest.mark.parametrize("wrong,what", sorted(HISTORICAL_WRONG_TAXIDS.items()))
def test_config_no_longer_carries_wrong_taxid(wrong, what):
    """No historical wrong value may survive anywhere in config.sh.

    This scans the raw text, not just the exports, so a stale value left behind
    in a comment or an example is caught too.
    """
    hits = [
        (n, line) for n, line in enumerate(CONFIG_SH.read_text().splitlines(), 1)
        if re.search(rf"(?<!\d){wrong}(?!\d)", line)
    ]
    assert not hits, (
        f"config.sh still references taxid {wrong} ({what}) at "
        + "; ".join(f"line {n}: {line.strip()}" for n, line in hits)
    )


# --- config.sh: LSE_LEVELS ---------------------------------------------------

def test_lse_levels_has_no_placeholder_text():
    """LSE_LEVELS shipped with literal unsubstituted placeholders."""
    levels = _sourced_lse_levels()
    assert levels, "LSE_LEVELS is empty"
    for level in levels:
        assert not re.search(r"taxid_[a-z]+\d*", level), (
            f"LSE_LEVELS still contains placeholder text: {level!r}"
        )


def test_lse_levels_taxids_are_all_numeric():
    for level in _sourced_lse_levels():
        _, _, taxid_field = level.partition(":")
        assert taxid_field, f"LSE_LEVELS entry has no taxid field: {level!r}"
        for taxid in taxid_field.split(","):
            assert taxid.strip().isdigit(), (
                f"LSE_LEVELS entry {level!r} has a non-numeric taxid {taxid!r}"
            )


def test_lse_levels_agrees_with_the_three_exported_taxids():
    """Each level must carry exactly the clade taxid of its matching variable."""
    levels = {}
    for level in _sourced_lse_levels():
        name, _, taxid_field = level.partition(":")
        levels[name] = [int(t) for t in taxid_field.split(",") if t.strip()]

    assert levels.get("Aeolids") == [EXPECTED_TAXIDS["LSE_AEOLID_TAXID"]]
    assert levels.get("Nudibranchs") == [EXPECTED_TAXIDS["LSE_NUDIBRANCH_TAXID"]]
    assert levels.get("Gastropods") == [EXPECTED_TAXIDS["LSE_GASTROPOD_TAXID"]]


def test_lse_levels_ordered_most_specific_to_most_general():
    """The documented contract, and what lse_refine.py's if/elif chain assumes."""
    names = [level.partition(":")[0] for level in _sourced_lse_levels()]
    assert names == ["Aeolids", "Nudibranchs", "Gastropods"]


# --- scripts/lse_refine.py: the os.getenv defaults ---------------------------

@pytest.mark.parametrize("name,expected", sorted(EXPECTED_TAXIDS.items()))
def test_lse_refine_default_matches_corrected_taxid(name, expected):
    """The defaults are live whenever the env is not exported."""
    assert _lse_refine_default(name) == expected


@pytest.mark.parametrize("name,expected", sorted(EXPECTED_TAXIDS.items()))
def test_lse_refine_default_matches_config(name, expected):
    """config.sh and lse_refine.py must never drift apart again."""
    assert _lse_refine_default(name) == int(_sourced_config_value(name))


# --- the defect itself: does classification actually fire now? ---------------

def test_classify_lse_level_fires_with_the_corrected_taxids(monkeypatch):
    """The end-to-end proof, at the level the bug actually lived.

    Correct constants are necessary but not sufficient: what matters is that
    `classify_lse_level` stops returning None for real molluscan input. Real
    NCBI lineages are supplied directly so this test needs neither the network
    nor a local ete3 taxonomy database.
    """
    import lse_refine

    # Genuine NCBI lineages. Berghia stephanieae and a second aeolid
    # (Aeolidiella, in Aeolidiidae) share ancestry up through Aeolidioidea.
    berghia = [131567, 2759, 33154, 33208, 6072, 33213, 33317, 2697495,
               1206795, 6447, 6448, 216305, 216307, 680346, 70849, 1707744,
               71481, 195871, 929455, 1287507]
    other_aeolid = berghia[:18] + [123456]          # diverges below Aeolidiidae
    other_nudibranch = berghia[:16] + [222222]      # shares up to Cladobranchia
    other_gastropod = berghia[:12] + [333333]       # shares up to Heterobranchia

    monkeypatch.setattr(
        lse_refine, "get_cached_lineage",
        lambda taxid: {
            1287507: berghia, 123456: other_aeolid,
            222222: other_nudibranch, 333333: other_gastropod,
        }.get(taxid),
    )

    assert lse_refine.classify_lse_level([1287507, 123456]) == "aeolids"
    assert lse_refine.classify_lse_level([1287507, 222222]) == "nudibranchs"
    assert lse_refine.classify_lse_level([1287507, 333333]) == "gastropods"


def test_classify_lse_level_returned_none_for_everything_before_the_fix(monkeypatch):
    """Pins WHY the bug was invisible: the old constants matched nothing.

    With the historical values restored, every molluscan input falls through
    every branch and yields None -- no error, no warning, just silence. This is
    the regression the guard exists to prevent.
    """
    import lse_refine

    berghia = [131567, 2759, 33154, 33208, 6072, 33213, 33317, 2697495,
               1206795, 6447, 6448, 216305, 216307, 680346, 70849, 1707744,
               71481, 195871, 929455, 1287507]
    monkeypatch.setattr(lse_refine, "get_cached_lineage", lambda taxid: berghia)
    monkeypatch.setattr(lse_refine, "LSE_AEOLID_TAXID", 54397)
    monkeypatch.setattr(lse_refine, "LSE_NUDIBRANCH_TAXID", 13843)
    monkeypatch.setattr(lse_refine, "LSE_GASTROPOD_TAXID", 644)

    assert lse_refine.classify_lse_level([1287507, 1287507]) is None
