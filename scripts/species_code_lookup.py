#!/usr/bin/env python3
"""Resolve reference-leaf codes (``ref_<code>_<locus>``) to species names.

Leaf codes in the trees are the header-prefix abbreviations used by the
reference proteomes (e.g. ``ref_seph_CAHIK_3`` → code ``seph``). This module
turns those codes into full ``Genus species`` names by merging two authoritative,
in-repo sources — nothing is hand-typed:

  1. ``SPECIES_MAP`` in ``recover_cds_from_assemblies.py`` (the miniprot-recovered
     reference assemblies: code → taxid, species, accession).
  2. The reference proteome filenames ``references/**/<taxid>_<Genus_species>.faa``,
     keyed by the header-prefix code found inside each file.

Where a single filename disagrees with the curated ``SPECIES_MAP``, the curated
entry wins and a warning is emitted. Where TWO reference proteomes claim the
same code, the code is AMBIGUOUS and lookups raise
:class:`AmbiguousSpeciesCodeError` rather than returning whichever species was
scanned first.

Known live ambiguity (2026-07): ``phau`` is the header prefix of both
``lse/gastropoda/109671_Physella_acuta.faa`` (848 seqs, Mollusca) and
``lse/other_lophotrochozoan_phyla/115415_Phoronis_australis.faa`` (251 seqs,
Phoronida). ``recover_cds_from_assemblies.SPECIES_MAP`` documents the one-time
``sed 's/^>phau_/>phaust_/'`` rename that fixes it; it has not been applied, so
phoronid leaves raise until the data is remediated.

Run as a script to dump the full code→species table as TSV.
"""
from __future__ import annotations

import glob
import os
import re
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPTS_DIR.parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from recover_cds_from_assemblies import SPECIES_MAP  # code -> (taxid, Species_name, accession)

# <taxid>_<Genus_species[...]>.faa  (filename is the authoritative species label)
FAA_RE = re.compile(r'(\d+)_([A-Z][A-Za-z]+_[a-z0-9_]+)\.(?:faa|fasta)$')


def _first_header_code(faa_path):
    """Header-prefix code of a proteome file (token after '>' up to first '_')."""
    with open(faa_path) as fh:
        for line in fh:
            if line.startswith('>'):
                return line[1:].strip().split()[0].split('_')[0]
    return None


class AmbiguousSpeciesCodeError(KeyError):
    """A header-prefix code is claimed by more than one reference proteome.

    Raised instead of returning whichever species happened to be seen first —
    a silent wrong answer here mislabels whole clades (see the ``phau``
    Physella/Phoronis collision documented in the module docstring).
    """


class CodeSpeciesMap(dict):
    """``{code: (taxid, 'Genus species')}`` plus the ambiguity register.

    Subclasses ``dict`` so existing callers keep working unchanged; the extra
    ``ambiguous`` attribute maps an ambiguous code to the sorted list of
    species that claim it.
    """

    def __init__(self, *args, ambiguous=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.ambiguous = dict(ambiguous or {})


def build_code_species_map(references_dir=None, warn=True):
    """Return a :class:`CodeSpeciesMap` merged from both sources.

    A code claimed by two different proteome FILES is recorded as ambiguous.
    The curated ``SPECIES_MAP`` still wins ordinary single-file disagreements,
    but it deliberately does NOT clear an ambiguity flag: that override is
    exactly what silently turned 251 phoronid sequences into a gastropod.
    """
    references_dir = Path(references_dir) if references_dir else PROJECT_ROOT / 'references'
    code_map = {}
    claims = {}   # code -> {species: first source file}

    files = sorted(glob.glob(str(references_dir / '**' / '*.faa'), recursive=True) +
                   glob.glob(str(references_dir / '**' / '*.fasta'), recursive=True))
    for faa in files:
        m = FAA_RE.search(os.path.basename(faa))
        if not m:
            continue
        taxid, species = m.group(1), m.group(2).replace('_', ' ')
        code = _first_header_code(faa)
        if not code:
            continue
        claims.setdefault(code, {}).setdefault(species, faa)
        prev = code_map.get(code)
        if prev and prev[1] != species and warn:
            print(f"[species_code_lookup] WARN: code '{code}' claimed by both "
                  f"'{prev[1]}' and '{species}'; marking AMBIGUOUS", file=sys.stderr)
        code_map.setdefault(code, (taxid, species))

    ambiguous = {code: sorted(sp) for code, sp in claims.items() if len(sp) > 1}

    # SPECIES_MAP is the curated source → it wins any single-file disagreement.
    for code, (taxid, species, _acc) in SPECIES_MAP.items():
        species = species.replace('_', ' ')
        prev = code_map.get(code)
        if prev and prev[1] != species and warn and code not in ambiguous:
            print(f"[species_code_lookup] WARN: code '{code}' filename='{prev[1]}' "
                  f"overridden by SPECIES_MAP='{species}'", file=sys.stderr)
        code_map[code] = (taxid, species)

    if ambiguous and warn:
        for code, species in sorted(ambiguous.items()):
            print(f"[species_code_lookup] AMBIGUOUS code '{code}': "
                  f"{', '.join(species)} — resolve by renaming the header prefix "
                  f"in the offending proteome; lookups will raise until then.",
                  file=sys.stderr)

    return CodeSpeciesMap(code_map, ambiguous=ambiguous)


def is_berghia(name):
    return bool(name) and (name.startswith('Berste') or name.startswith('TRINITY_'))


def parse_leaf(name):
    """``name`` → ``(kind, code, locus)``; kind ∈ {'berghia','ref','other'}."""
    if is_berghia(name):
        return 'berghia', None, name
    if name and name.startswith('ref_'):
        code, _sep, locus = name[4:].partition('_')
        return 'ref', code, locus
    return 'other', None, name or ''


def species_for(name, code_map, strict=True):
    """``name`` → ``(species_or_None, locus)`` for a reference leaf.

    Raises :class:`AmbiguousSpeciesCodeError` when the leaf's code is claimed
    by more than one reference proteome. Pass ``strict=False`` to degrade to
    ``(None, locus)`` — 'unknown' is an acceptable answer, a confident wrong
    species is not.
    """
    kind, code, locus = parse_leaf(name)
    if kind == 'ref':
        ambiguous = getattr(code_map, 'ambiguous', {})
        if code in ambiguous:
            if strict:
                raise AmbiguousSpeciesCodeError(
                    f"species code '{code}' is claimed by "
                    f"{len(ambiguous[code])} proteomes ({', '.join(ambiguous[code])}); "
                    f"refusing to guess for leaf '{name}'")
            return None, locus
        entry = code_map.get(code)
        return (entry[1] if entry else None), locus
    return None, locus


def berghia_display_name(name):
    """Berghia transcript id → readable species label.

    'BersteEVm009528t6' -> 'Berghia stephanieae EVm009528t6'. The 'Berste'
    species prefix is expanded to the binomial; the transcript suffix is kept
    to distinguish paralogs.
    """
    if name and name.startswith('Berste'):
        return 'Berghia stephanieae ' + name[len('Berste'):]
    return 'Berghia stephanieae ' + (name or '')


def main():
    code_map = build_code_species_map()
    print("code\ttaxid\tspecies")
    for code in sorted(code_map):
        taxid, species = code_map[code]
        print(f"{code}\t{taxid}\t{species}")
    print(f"# {len(code_map)} codes", file=sys.stderr)
    if code_map.ambiguous:
        print(f"# {len(code_map.ambiguous)} AMBIGUOUS code(s): "
              f"{', '.join(sorted(code_map.ambiguous))}", file=sys.stderr)


if __name__ == '__main__':
    main()
