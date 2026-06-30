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

On a code collision the SPECIES_MAP entry wins (it is the curated source) and a
warning is emitted. Run as a script to dump the full code→species table as TSV.
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


def build_code_species_map(references_dir=None, warn=True):
    """Return ``{code: (taxid, 'Genus species')}`` merged from both sources."""
    references_dir = Path(references_dir) if references_dir else PROJECT_ROOT / 'references'
    code_map = {}

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
        prev = code_map.get(code)
        if prev and prev[1] != species and warn:
            print(f"[species_code_lookup] WARN: code '{code}' seen for both "
                  f"'{prev[1]}' and '{species}'; keeping '{prev[1]}'", file=sys.stderr)
        code_map.setdefault(code, (taxid, species))

    # SPECIES_MAP is the curated source → it wins any disagreement.
    for code, (taxid, species, _acc) in SPECIES_MAP.items():
        species = species.replace('_', ' ')
        prev = code_map.get(code)
        if prev and prev[1] != species and warn:
            print(f"[species_code_lookup] WARN: code '{code}' filename='{prev[1]}' "
                  f"overridden by SPECIES_MAP='{species}'", file=sys.stderr)
        code_map[code] = (taxid, species)

    return code_map


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


def species_for(name, code_map):
    """``name`` → ``(species_or_None, locus)`` for a reference leaf."""
    kind, code, locus = parse_leaf(name)
    if kind == 'ref':
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


if __name__ == '__main__':
    main()
