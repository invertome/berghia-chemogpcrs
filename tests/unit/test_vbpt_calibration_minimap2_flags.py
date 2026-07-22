"""Bead vbpt / kfqz: every spliced minimap2 alignment that feeds
reconcile_candidates.parse_minimap2_paf must emit the intron-aware identity
source, i.e. carry `-c --cs`.

Under `-x splice` alone, PAF `alen` spans introns, so identity read as
nmatch/alen collapses on multi-exon transcripts; the fix reads the `de:f`/`cs`
tag instead and FAILS LOUD on a tag-less PAF. So stage 02c AND both genome-track
calibration wrappers (which produce their PAF "exactly as 02c does" and then feed
it to the same parser) must generate the PAF the same way. This guards against any
of the three drifting back to a bare `-x splice`.
"""
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

# scripts whose `$MINIMAP2 -x splice` output is parsed by parse_minimap2_paf
SPLICE_ALIGNERS = [
    PROJECT_ROOT / "02c_genome_reconcile.sh",
    PROJECT_ROOT / "scripts" / "unity" / "calibrate_genome_track.sh",
    PROJECT_ROOT / "scripts" / "unity" / "calibrate_genome_track_margin.sh",
]


def _minimap2_command_lines(path: Path):
    """Real `$MINIMAP2 -x splice ...` invocations (comments/prose stripped)."""
    cmds = []
    for raw in path.read_text().splitlines():
        code = raw.split("#", 1)[0]  # drop trailing comments; comment-only lines -> ""
        if "MINIMAP2" in code and "-x splice" in code:
            cmds.append(code)
    return cmds


def test_splice_aligners_emit_intron_aware_identity_tag():
    missing = []
    for f in SPLICE_ALIGNERS:
        assert f.exists(), f"aligner script not found: {f}"
        cmds = _minimap2_command_lines(f)
        assert cmds, f"no `$MINIMAP2 -x splice` command found in {f.name}"
        for cmd in cmds:
            if not ((" -c " in f" {cmd} ") and ("--cs" in cmd)):
                missing.append((f.name, cmd.strip()))
    assert not missing, (
        "spliced minimap2 commands lacking `-c --cs` (needed so the PAF carries the "
        "de:f/cs identity tag reconcile_candidates.py requires):\n"
        + "\n".join(f"  {name}: {cmd}" for name, cmd in missing)
    )
