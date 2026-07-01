"""Drift guard: the coarse non-chemoreceptor GPCR family list is hand-mirrored
across three places (structural_evidence.py's own docstring: "an
authoritative mirror of COARSE_FAMILIES in
scripts/validate_classification_hmms.py (identical to the family list in
scripts/curate_gpcr_references.py). Keep the three copies in sync."). This
test binds the two importable Python constants together so an edit to one
copy without the other fails CI immediately, rather than silently letting a
Foldseek EXCLUSION call (structural_evidence.classify_hit) and the HMM LOO
validation (validate_classification_hmms) disagree about which coarse
families count as "known non-chemoreceptor".

validate_classification_hmms.py only imports stdlib modules at module scope
(argparse/csv/os/subprocess/sys/tempfile/collections/pathlib/typing) and has
no import-time side effects, so importing the module directly to reach its
constant is cheap and safe (matches the existing
tests/unit/test_validate_classification_hmms.py import style).
"""
from structural_evidence import NON_CHEMORECEPTOR_FAMILIES
from validate_classification_hmms import COARSE_FAMILIES


def test_structural_evidence_and_validate_classification_hmms_family_lists_match():
    assert NON_CHEMORECEPTOR_FAMILIES == COARSE_FAMILIES
