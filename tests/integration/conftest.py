"""Session fixtures for the pipeline seam canary.

The seam graph is built once per session: sourcing config.sh and scanning
every stage costs well under a second, but there is no reason to pay it per
test.
"""

from __future__ import annotations

import collections
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import seamlib  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return REPO_ROOT


@pytest.fixture(scope="session")
def config_vars() -> dict[str, str]:
    """Every variable the REAL config.sh defines, obtained by sourcing it."""
    return seamlib.load_config_vars(REPO_ROOT)


@pytest.fixture(scope="session")
def results_dir(config_vars: dict[str, str]) -> str:
    return config_vars["RESULTS_DIR"]


@pytest.fixture(scope="session")
def stage_scripts() -> list[Path]:
    """The numbered stage scripts, in registry order on disk."""
    return sorted(REPO_ROOT.glob("[0-9]*_*.sh"))


@pytest.fixture(scope="session")
def seam_graph(
    config_vars: dict[str, str], results_dir: str, stage_scripts: list[Path]
) -> seamlib.SeamGraph:
    """Every probe and every write the pipeline expresses.

    Producers deliberately include scripts/ (one-shot prep and helpers) and
    functions.sh, because "some other part of this repository does write it"
    is a legitimate answer to "who produces this path?".
    """
    scripts_dir = REPO_ROOT / "scripts"
    producers = (
        stage_scripts
        + sorted(scripts_dir.rglob("*.sh"))
        + [REPO_ROOT / "functions.sh"]
    )
    helpers = seamlib.infer_helper_output_positions(REPO_ROOT / "functions.sh")

    writes: list[seamlib.Write] = []
    for script in producers:
        writes += seamlib.extract_writes(script, config_vars, results_dir)
        writes += seamlib.extract_heredoc_writes(script, config_vars, results_dir)
        writes += seamlib.extract_heredoc_var_writes(script, config_vars, results_dir)
        writes += seamlib.extract_two_hop_writes(
            script, config_vars, results_dir, scripts_dir
        )
        writes += seamlib.extract_helper_call_writes(
            script, config_vars, results_dir, helpers
        )

    probes: list[seamlib.Probe] = []
    for stage in stage_scripts:
        probes += seamlib.extract_probes(stage, config_vars, results_dir)
    for script in sorted(scripts_dir.glob("*.py")):
        probes += seamlib.extract_python_probes(script, results_dir)

    return seamlib.SeamGraph(probes=probes, writes=writes)


@pytest.fixture(scope="session")
def probe_groups(seam_graph: seamlib.SeamGraph) -> dict[str, list[seamlib.Probe]]:
    """Probes bucketed into OR-groups.

    A fallback chain (`for f in a b c; do [ -f "$f" ]`) is satisfied when ANY
    alternative has a producer. Stage 02a's five-way candidate lookup is that
    shape, and the historical defect was that all of its alternatives were
    unwritable at once.
    """
    groups: dict[str, list[seamlib.Probe]] = collections.defaultdict(list)
    for probe in seam_graph.probes:
        groups[probe.group].append(probe)
    return dict(groups)


def pytest_configure(config: pytest.Config) -> None:
    """Register the marker locally so the canary needs no root pytest.ini.

    Adding a root config file would change how `pytest tests/unit -q` resolves
    its rootdir, and this canary must not disturb the existing suite.
    """
    config.addinivalue_line(
        "markers", "seam: pipeline stage-to-stage data-flow seam checks"
    )
