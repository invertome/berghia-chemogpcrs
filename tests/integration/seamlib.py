"""Seam extraction library for the pipeline integration canary.

WHY THIS EXISTS
---------------
The pipeline's ~2057 unit tests all passed while an entire stage (02c) never
produced output, four scored ranking axes were dormant, and two producer
scripts were invoked by nothing. Every *unit* was correct. The defects lived
in the *seams*: the path stage B probes versus the path stage A writes.

A unit test cannot reach a seam, because a seam is not inside any unit. This
module extracts both sides of every seam from the real stage scripts and lets
the canary assert the graph connects.

THE CONTRACT WE CHECK
---------------------
    Every existence-guard probe of a path under RESULTS_DIR must be
    satisfiable by some write performed elsewhere in the pipeline.

"Existence guard" is deliberately narrow -- ``[ -f X ]``, ``[ -d X ]``,
``check_file X``, ``Path(results).glob(...)``. That narrowness is the point:
every historical defect manifested at a guard that was never satisfied, and
a guard is the exact construct that lets a stage skip its work *silently*
while still touching its completion flag.

Paths are compared as PATTERNS, never as fixed filenames, so that a stage may
freely rename or re-derive a file as long as producer and consumer still agree.
"""

from __future__ import annotations

import fnmatch
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

# A rare token substituted for a wildcard when testing pattern intersection
# from the other direction. Must not occur in any real path.
_SENTINEL = "ZqSeamWildcardZq"


# --------------------------------------------------------------------------
# Config resolution: we SOURCE the real config.sh rather than reimplementing
# its variable expansion. Anything else would be a second source of truth,
# and a second source of truth is how seams rot in the first place.
# --------------------------------------------------------------------------

def load_config_vars(repo_root: Path = REPO_ROOT) -> dict[str, str]:
    """Return every scalar variable defined by the real config.sh."""
    script = (
        "set +u\n"
        'source "%s/config.sh" >/dev/null 2>&1\n'
        "declare -p\n" % repo_root
    )
    proc = subprocess.run(
        ["bash", "-c", script], capture_output=True, text=True, timeout=60
    )
    out: dict[str, str] = {}
    pattern = re.compile(r'^declare -\S+ ([A-Za-z_][A-Za-z0-9_]*)="(.*)"$')
    for line in proc.stdout.splitlines():
        match = pattern.match(line)
        if match:
            out[match.group(1)] = match.group(2)
    if "RESULTS_DIR" not in out:
        raise RuntimeError(
            "could not source config.sh; canary cannot resolve any path.\n"
            f"stderr: {proc.stderr[:2000]}"
        )
    return out


# --------------------------------------------------------------------------
# Shell reading
# --------------------------------------------------------------------------

def logical_lines(text: str) -> list[tuple[int, str]]:
    """Join backslash continuations, keeping the first physical line number.

    Multi-line command invocations are the norm in these stages, and both
    sides of a seam are routinely split across continuations.
    """
    joined: list[tuple[int, str]] = []
    buffer = ""
    start: int | None = None
    for lineno, raw in enumerate(text.splitlines(), 1):
        if start is None:
            start = lineno
        stripped = raw.rstrip()
        if stripped.endswith("\\"):
            buffer += stripped[:-1] + " "
            continue
        buffer += stripped
        joined.append((start, buffer))
        buffer, start = "", None
    if buffer and start is not None:
        joined.append((start, buffer))
    return joined


_ASSIGN = re.compile(
    r"^\s*(?:local\s+|export\s+|declare\s+-\S+\s+)?"
    r"([A-Za-z_][A-Za-z0-9_]*)=(.+)$"
)
_VAR = re.compile(
    r"\$\{!?([A-Za-z_][A-Za-z0-9_]*)(?:(?::-|:=|##|%%|#|%)[^}]*)?\}"
    r"|\$([A-Za-z_][A-Za-z0-9_]*)"
)


def build_symbol_table(
    lines: list[tuple[int, str]], base: dict[str, str]
) -> dict[str, str]:
    """Resolve in-script assignments on top of the config.sh variables.

    Stages routinely alias a directory once (``RECON_DIR="${RESULTS_DIR}/..."``)
    and then use the alias everywhere, so without this pass most seams would
    be invisible.
    """
    table = dict(base)
    for _lineno, line in lines:
        head = unquote(line).split(";")[0]
        match = _ASSIGN.match(head)
        if not match:
            continue
        name, rhs = match.group(1), match.group(2).strip()
        # Command substitution is genuinely dynamic; leave the name unknown so
        # it degrades to a wildcard rather than to a confident wrong value.
        if "$(" in rhs or "`" in rhs:
            continue
        if len(rhs) >= 2 and rhs[0] == "'" and rhs[-1] == "'":
            rhs = rhs[1:-1]
        if "/" in rhs or "$" in rhs:
            table[name] = rhs
    return table


_POSITIONAL_REF = re.compile(r"\$\{?([1-9])\}?")


def expand(text: str, table: dict[str, str], depth: int = 0) -> str:
    """Expand shell variables; unknown ones collapse to a ``*`` wildcard."""
    if depth > 10:
        return text

    # Positional parameters inside helper functions are caller-supplied, so
    # they are wildcards rather than literals. Left as "$1" they would read as
    # a filename in failure messages and confuse whoever is debugging.
    text = _POSITIONAL_REF.sub("*", text)

    def replace(match: re.Match) -> str:
        name = match.group(1) or match.group(2)
        if name in table:
            return expand(table[name], table, depth + 1)
        return "*"

    return _VAR.sub(replace, text)


# --------------------------------------------------------------------------
# Seam records
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Probe:
    """A place where some stage asks 'does this path exist?'."""

    stage: str
    lineno: int
    raw: str
    pattern: str
    kind: str          # "file" | "dir" | "tree"
    group: str         # probes sharing a group are alternatives (OR)

    @property
    def location(self) -> str:
        return f"{self.stage}:{self.lineno}"


@dataclass(frozen=True)
class Write:
    """A place where some stage claims to produce a path."""

    stage: str
    lineno: int
    pattern: str
    how: str

    @property
    def location(self) -> str:
        return f"{self.stage}:{self.lineno}"


@dataclass
class SeamGraph:
    probes: list[Probe] = field(default_factory=list)
    writes: list[Write] = field(default_factory=list)


# --------------------------------------------------------------------------
# Pattern intersection
#
# Both sides carry wildcards, so this is glob-vs-glob, not glob-vs-string.
# Matching is SEGMENT-WISE (split on "/") because the historical 04->05 defect
# was purely a directory-depth mismatch: a matcher whose "*" crossed "/" would
# have called those two paths compatible and missed the bug entirely.
# --------------------------------------------------------------------------

def _segments_intersect(left: str, right: str) -> bool:
    if left == right:
        return True
    if fnmatch.fnmatchcase(left.replace("*", _SENTINEL), right):
        return True
    if fnmatch.fnmatchcase(right.replace("*", _SENTINEL), left):
        return True
    return False


def patterns_intersect(left: str, right: str) -> bool:
    """True if some concrete path could match both patterns."""
    lparts, rparts = left.split("/"), right.split("/")
    if len(lparts) != len(rparts):
        return False
    return all(_segments_intersect(a, b) for a, b in zip(lparts, rparts))


def write_is_under(write_pattern: str, dir_pattern: str) -> bool:
    """True if the write lands inside the probed directory (or is it)."""
    dparts = dir_pattern.split("/")
    wparts = write_pattern.split("/")
    if len(wparts) < len(dparts):
        return False
    return all(_segments_intersect(a, b) for a, b in zip(wparts[: len(dparts)], dparts))


def write_satisfies(write: Write, probe: Probe) -> bool:
    if probe.kind in ("dir", "tree"):
        return write_is_under(write.pattern, probe.pattern)
    return patterns_intersect(write.pattern, probe.pattern)


# --------------------------------------------------------------------------
# Probe extraction
# --------------------------------------------------------------------------

_GUARD = re.compile(r"\[\[?\s+-([fsdeLr])\s+\"?([^\"\]\s]+)\"?")
_CHECK_FILE = re.compile(r"\bcheck_file\s+\"?([^\"\s;)]+)\"?")
_FOR_IN = re.compile(r"^\s*for\s+([A-Za-z_][A-Za-z0-9_]*)\s+in\s+(.+?)(?:;\s*do|\s+do)?\s*$")
_FIND_ROOT = re.compile(r"\bfind\s+\"?([^\"\s;)]+)\"?")

_DIR_FLAGS = {"d"}


def unquote(line: str) -> str:
    """Drop double quotes so adjacent quoted/unquoted fragments stay ONE word.

    ``"${RESULTS_DIR}/busco/busco_"*`` is a single glob, but any tokenizer that
    treats the closing quote as a word boundary silently truncates it to
    ``.../busco_`` and then reports a false unresolved seam. Paths in this
    repository never contain spaces, so removing the quote character is both
    safe and much more faithful than trying to model shell word splitting.

    One case must be protected first: a QUOTED shell metacharacter is an
    argument, not an operator. ``grep -c ">" "$fa"`` counts FASTA headers; if
    the quotes are stripped naively it reads as a redirect into ``$fa`` and the
    file is misreported as a write. Neutralise those before unquoting.
    """
    line = re.sub(r"\"[<>|&]+\"|'[<>|&]+'", "SHELLMETA", line)
    return line.replace('"', "")


def _tokenise(items: str) -> list[str]:
    """Split an (already unquoted) shell word list."""
    return items.split()


def extract_probes(path: Path, config: dict[str, str], results_dir: str) -> list[Probe]:
    """Pull every existence-guard probe of a RESULTS_DIR path out of a shell script."""
    text = path.read_text(errors="replace")
    lines = logical_lines(text)
    table = build_symbol_table(lines, config)
    stage = path.name
    probes: list[Probe] = []

    for lineno, raw_line in lines:
        if raw_line.strip().startswith("#"):
            continue
        line = unquote(raw_line)

        # --- `for VAR in a b c` : the items are alternatives, not requirements.
        # Stage 02a's five-way candidate fallback is exactly this shape, and the
        # historical defect was that ALL of its alternatives were unwritable.
        for_match = _FOR_IN.match(line)
        if for_match:
            group = f"{stage}:{lineno}:for"
            emitted = False
            for token in _tokenise(for_match.group(2)):
                if "$" not in token:
                    continue
                pattern = expand(token, table)
                if pattern.startswith(results_dir):
                    probes.append(
                        Probe(stage, lineno, token, pattern, "file", group)
                    )
                    emitted = True
            if emitted:
                continue

        # --- plain existence guards. Guards joined by `||` on one line are
        # alternatives; anything else stands alone and must resolve.
        or_joined = "||" in line
        group = f"{stage}:{lineno}:or" if or_joined else None
        for match in _GUARD.finditer(line):
            flag, token = match.group(1), match.group(2)
            if "$" not in token:
                continue
            pattern = expand(token, table)
            if not pattern.startswith(results_dir):
                continue
            kind = "dir" if flag in _DIR_FLAGS else "file"
            probes.append(
                Probe(
                    stage, lineno, token, pattern, kind,
                    group or f"{stage}:{lineno}:{token}",
                )
            )

        for match in _CHECK_FILE.finditer(line):
            token = match.group(1)
            if "$" not in token:
                continue
            pattern = expand(token, table)
            if pattern.startswith(results_dir):
                probes.append(
                    Probe(stage, lineno, token, pattern, "file",
                          f"{stage}:{lineno}:check_file")
                )

        # --- `find <root>` : the search root is a read of a whole subtree.
        for match in _FIND_ROOT.finditer(line):
            token = match.group(1)
            if "$" not in token:
                continue
            pattern = expand(token, table).rstrip("/")
            if pattern.startswith(results_dir) and pattern != results_dir:
                probes.append(
                    Probe(stage, lineno, token, pattern, "tree",
                          f"{stage}:{lineno}:find")
                )

    return probes


_PY_GLOB = re.compile(
    r"Path\(\s*(?:self\.)?(\w*results\w*)\s*\)\s*\.\s*glob\(\s*[\"']([^\"']+)[\"']"
)
# Multi-segment joins, e.g. os.path.join(results_dir, 'cafe', 'x.csv').
# The ranking script reads its scored axes this way, and every one of the four
# axes that were found dormant had this exact shape -- so a probe extractor
# that only understood .glob() would miss the whole class.
_PY_JOIN_EXISTS = re.compile(
    r"os\.path\.join\(\s*(\w*results\w*)\s*,\s*"
    r"((?:[\"'][^\"']+[\"']\s*,\s*)*[\"'][^\"']+[\"'])\s*\)"
)
_PY_SEGMENT = re.compile(r"[\"']([^\"']+)[\"']")


def extract_python_probes(
    path: Path, results_dir: str, only_exists: bool = True
) -> list[Probe]:
    """Probes expressed in Python against a results-directory root.

    Only ``.glob()`` is treated as a probe: it is the construct that silently
    yields nothing when the path is wrong, which is precisely how the ranking
    script's dormant axes failed. ``os.path.join`` results are collected only
    when the same line tests for existence.
    """
    text = path.read_text(errors="replace")
    probes: list[Probe] = []
    stage = f"scripts/{path.name}"

    for lineno, line in enumerate(text.splitlines(), 1):
        if line.strip().startswith("#"):
            continue
        for match in _PY_GLOB.finditer(line):
            rel = match.group(2)
            probes.append(
                Probe(stage, lineno, rel, f"{results_dir}/{rel}", "file",
                      f"{stage}:{lineno}:glob")
            )

    for match in _PY_JOIN_EXISTS.finditer(text):
        segments = _PY_SEGMENT.findall(match.group(2))
        if not segments:
            continue
        rel = "/".join(segment.strip("/") for segment in segments)
        lineno = text[: match.start()].count("\n") + 1
        # A joined path whose last segment carries no extension names a
        # DIRECTORY (results_dir/'orthogroups'), and a directory probe is
        # satisfied by any write landing inside it. Treating it as a file
        # would demand an exact depth match and invent a broken seam.
        kind = "file" if "." in segments[-1] else "dir"
        probes.append(
            Probe(stage, lineno, rel, f"{results_dir}/{rel}", kind,
                  f"{stage}:{lineno}:join")
        )

    # Multi-line `.glob(` calls: the ranking script wraps them, so a purely
    # line-based scan would miss the very defect this canary exists for.
    flat = re.sub(r"\s*\n\s*", " ", text)
    for match in _PY_GLOB.finditer(flat):
        rel = match.group(2)
        pattern = f"{results_dir}/{rel}"
        if not any(p.pattern == pattern for p in probes):
            lineno = text[: text.find(rel)].count("\n") + 1 if rel in text else 0
            probes.append(
                Probe(stage, lineno, rel, pattern, "file", f"{stage}:{lineno}:glob")
            )
    return probes


# --------------------------------------------------------------------------
# Write extraction
# --------------------------------------------------------------------------

_REDIRECT = re.compile(r"(?<![0-9&])>>?\s*\"?([^\"\s;|&)]+)\"?")
_OUT_FLAG = re.compile(
    r"(?:^|\s)(?:-o|-out|--out|--output|--output-dir|--outdir|--output-file|--out-dir)"
    r"[=\s]+\"?([^\"\s;|&)]+)\"?"
)
_PREFIX_FLAG = re.compile(r"(?:^|\s)(?:--prefix|-pre)[=\s]+\"?([^\"\s;|&)]+)\"?")
_MKDIR = re.compile(r"\bmkdir\s+(?:-p\s+)?(.+?)(?:\s*(?:\|\||&&|;|$))")
_TOUCH = re.compile(r"\btouch\s+(.+?)(?:\s*(?:\|\||&&|;|$))")
_STDOUT_FLAG = re.compile(r"--std(?:out|err)=\"?([^\"\s;|&)]+)\"?")
_TEE = re.compile(r"\btee\s+(?:-a\s+)?\"?([^\"\s;|&)]+)\"?")
_COPYLIKE = re.compile(r"\b(?:cp|mv|rsync|ln)\s+(.+?)(?:\s*(?:\|\||&&|;|$))")
_VALIDATE = re.compile(r"\bvalidate_outputs\s+(.+?)$")


def extract_writes(path: Path, config: dict[str, str], results_dir: str) -> list[Write]:
    """Pull every path a shell script claims to produce."""
    text = path.read_text(errors="replace")
    lines = logical_lines(text)
    table = build_symbol_table(lines, config)
    stage = path.name
    writes: list[Write] = []

    def add(token: str, lineno: int, how: str, also_dotstar: bool = False) -> None:
        if "$" not in token:
            return
        pattern = expand(token, table).rstrip("/")
        if not pattern.startswith(results_dir):
            return
        writes.append(Write(stage, lineno, pattern, how))
        if also_dotstar:
            writes.append(Write(stage, lineno, pattern + ".*", how))

    for lineno, raw_line in lines:
        if raw_line.strip().startswith("#"):
            continue
        line = unquote(raw_line)

        for match in _REDIRECT.finditer(line):
            target = match.group(1)
            if target.startswith("/dev/") or target.startswith("&"):
                continue
            add(target, lineno, "redirect")

        for match in _OUT_FLAG.finditer(line):
            add(match.group(1), lineno, "output-flag", also_dotstar=True)

        for match in _STDOUT_FLAG.finditer(line):
            add(match.group(1), lineno, "run_command-stdout")

        for match in _PREFIX_FLAG.finditer(line):
            # IQ-TREE/HyPhy prefixes name a family of siblings, not one file.
            add(match.group(1), lineno, "prefix", also_dotstar=True)

        for match in _MKDIR.finditer(line):
            for token in _tokenise(match.group(1)):
                if token.startswith("-"):
                    continue
                add(token, lineno, "mkdir")

        for match in _TOUCH.finditer(line):
            for token in _tokenise(match.group(1)):
                if token.startswith("-"):
                    continue
                add(token, lineno, "touch")

        for match in _TEE.finditer(line):
            add(match.group(1), lineno, "tee")

        for match in _COPYLIKE.finditer(line):
            tokens = [t for t in _tokenise(match.group(1)) if not t.startswith("-")]
            if tokens:
                add(tokens[-1], lineno, "copy-dest")

        for match in _VALIDATE.finditer(line):
            for token in _tokenise(match.group(1)):
                if token.startswith("-"):
                    continue
                add(token, lineno, "validate_outputs")

    return writes


# --------------------------------------------------------------------------
# Producers that are not shell redirects
#
# Two mechanisms in this pipeline create files without any shell-visible write,
# and both must be modelled or the canary drowns in false alarms:
#
#   (a) inline `python3 << 'HEREDOC'` blocks embedded in a stage, which take
#       their output directory from an exported shell variable;
#   (b) scripts/*.py invoked with an output flag, which then join their own
#       relative filenames onto whatever directory the stage handed them.
# --------------------------------------------------------------------------

_PY_ENV_BIND = re.compile(
    r"([A-Za-z_]\w*)\s*=\s*os\.environ\.get\(\s*[\"'](\w+)[\"']"
)
_PY_FSTRING_WRITE = re.compile(r"f[\"']\{([A-Za-z_]\w*)\}/([^\"'{}]+)[\"']")
_PY_JOIN_WRITE = re.compile(
    r"os\.path\.join\(\s*([A-Za-z_]\w*)\s*,\s*[\"']([^\"']+)[\"']\s*\)"
)
_PY_PATH_DIV = re.compile(
    r"Path\(\s*([A-Za-z_]\w*)\s*\)(?:\.parent)?\s*/\s*[\"']([^\"']+)[\"']"
)


def extract_heredoc_writes(
    path: Path, config: dict[str, str], results_dir: str
) -> list[Write]:
    """Writes performed by Python embedded inside a shell stage.

    Stage 03c builds the CAFE input this way: the heredoc reads ``CAFE_DIR``
    out of the environment and writes ``{cafe_dir}/gene_families.tsv``, which
    no shell-level redirect ever mentions.
    """
    text = path.read_text(errors="replace")
    lines = logical_lines(text)
    table = build_symbol_table(lines, config)
    stage = path.name

    bindings: dict[str, str] = {}
    for match in _PY_ENV_BIND.finditer(text):
        pyvar, shellvar = match.group(1), match.group(2)
        if shellvar in table:
            bindings[pyvar] = expand(table[shellvar], table)

    writes: list[Write] = []
    for regex in (_PY_FSTRING_WRITE, _PY_JOIN_WRITE, _PY_PATH_DIV):
        for match in regex.finditer(text):
            pyvar, leaf = match.group(1), match.group(2)
            root = bindings.get(pyvar)
            if not root or not root.startswith(results_dir):
                continue
            lineno = text[: match.start()].count("\n") + 1
            writes.append(
                Write(stage, lineno, f"{root.rstrip('/')}/{leaf}", "heredoc-python")
            )
    return writes


_OUTPUT_FLAG_NAME = re.compile(
    r"^(?:-o|--out|--output.*|--.*-out|--.*-dir|--.*-tsv|--.*-csv|--.*-json|--.*-fasta)$"
)
_PY_INVOKE = re.compile(r"(?:\$\{SCRIPTS_DIR\}|scripts)/([A-Za-z0-9_]+\.py)")


def script_relative_writes(script: Path) -> set[str]:
    """Relative filenames a Python script joins onto a caller-supplied dir."""
    text = script.read_text(errors="replace")
    leaves: set[str] = set()
    for regex in (_PY_JOIN_WRITE, _PY_FSTRING_WRITE, _PY_PATH_DIV):
        for match in regex.finditer(text):
            leaf = match.group(2).strip("/")
            if leaf and not leaf.startswith("."):
                leaves.add(leaf)
    return leaves


def extract_two_hop_writes(
    path: Path, config: dict[str, str], results_dir: str, scripts_dir: Path
) -> list[Write]:
    """Writes a stage delegates to a Python script it invokes.

    The stage supplies the directory (``--output-dir results/classification/loo``)
    and the script supplies the filename (``loo_metrics.tsv``). Neither half is
    a complete path on its own, so a single-file scan sees no producer at all.
    """
    text = path.read_text(errors="replace")
    lines = logical_lines(text)
    table = build_symbol_table(lines, config)
    stage = path.name
    writes: list[Write] = []

    for lineno, raw_line in lines:
        if raw_line.strip().startswith("#"):
            continue
        line = unquote(raw_line)
        invoked = _PY_INVOKE.search(line)
        if not invoked:
            continue
        script = scripts_dir / invoked.group(1)
        if not script.is_file():
            continue
        leaves = script_relative_writes(script)
        if not leaves:
            continue

        tokens = line.split()
        for index, token in enumerate(tokens[:-1]):
            flag, value = token, tokens[index + 1]
            if "=" in token and token.startswith("-"):
                flag, value = token.split("=", 1)
            elif not _OUTPUT_FLAG_NAME.match(flag):
                continue
            if not _OUTPUT_FLAG_NAME.match(flag):
                continue
            root = expand(value, table).rstrip("/")
            if not root.startswith(results_dir):
                continue
            for leaf in leaves:
                writes.append(
                    Write(stage, lineno, f"{root}/{leaf}", f"via {invoked.group(1)}")
                )
    return writes


# --------------------------------------------------------------------------
# Shell helper functions whose OUTPUTS are positional arguments
#
# `identify_gpcr_candidates "$in" "$out" "$census"` produces two files, but a
# regex scan of the call site sees only three bare words. Rather than hardcode
# a table (which would rot), infer each helper's output positions from the body
# of the helper itself in functions.sh.
# --------------------------------------------------------------------------

_FUNC_START = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*\(\)\s*\{")
_POSITIONAL = re.compile(
    r"^\s*local\s+([A-Za-z_]\w*)\s*=\s*\"?\$\{?(\d+)(?::-[^}]*)?\}?\"?"
)


def infer_helper_output_positions(functions_sh: Path) -> dict[str, set[int]]:
    """Map helper name -> positional argument indices that it writes to."""
    text = functions_sh.read_text(errors="replace")
    lines = text.splitlines()
    helpers: dict[str, set[int]] = {}

    current: str | None = None
    depth = 0
    params: dict[str, int] = {}
    body: list[str] = []

    def finish() -> None:
        if current is None:
            return
        joined = "\n".join(body)
        written: set[int] = set()
        for name, position in params.items():
            reference = re.compile(r"\$\{?" + re.escape(name) + r"\}?\b")
            for _lineno, line in logical_lines(joined):
                clean = unquote(line)
                if not reference.search(clean):
                    continue
                if (
                    _REDIRECT.search(clean)
                    or _OUT_FLAG.search(clean)
                    or _TOUCH.search(clean)
                    or _TEE.search(clean)
                    or _STDOUT_FLAG.search(clean)
                ):
                    # Only count it when the reference is on the write side.
                    for regex in (_REDIRECT, _OUT_FLAG, _TOUCH, _TEE, _STDOUT_FLAG):
                        for match in regex.finditer(clean):
                            if reference.search(match.group(0)):
                                written.add(position)
        if written:
            helpers[current] = written

    for raw in lines:
        start = _FUNC_START.match(raw)
        if start and depth == 0:
            finish()
            current, params, body, depth = start.group(1), {}, [], 1
            continue
        if current is None:
            continue
        depth += raw.count("{") - raw.count("}")
        body.append(raw)
        positional = _POSITIONAL.match(raw)
        if positional:
            params[positional.group(1)] = int(positional.group(2))
        if depth <= 0:
            finish()
            current, params, body, depth = None, {}, [], 0
    finish()
    return helpers


def extract_helper_call_writes(
    path: Path,
    config: dict[str, str],
    results_dir: str,
    helpers: dict[str, set[int]],
) -> list[Write]:
    """Writes made by helper calls whose output paths are positional args."""
    text = path.read_text(errors="replace")
    lines = logical_lines(text)
    table = build_symbol_table(lines, config)
    stage = path.name
    writes: list[Write] = []

    for lineno, raw_line in lines:
        if raw_line.strip().startswith("#"):
            continue
        line = unquote(raw_line)
        tokens = line.split()
        if not tokens:
            continue
        for index, token in enumerate(tokens):
            positions = helpers.get(token)
            if not positions:
                continue
            for position in positions:
                target = index + position
                if target >= len(tokens):
                    continue
                pattern = expand(tokens[target], table).rstrip("/")
                if pattern.startswith(results_dir):
                    writes.append(
                        Write(stage, lineno, pattern, f"{token} arg {position}")
                    )
    return writes


# --------------------------------------------------------------------------
# Heredoc Python that takes its paths by shell interpolation
#
# An UNQUOTED heredoc (`python3 << PYTHON_SCRIPT`) lets the shell expand
# `${VAR}` inside the Python source, so stage 03d writes its NOTUNG species
# tree through a Python variable that never appears as a literal path.
# --------------------------------------------------------------------------

_PY_SHELL_BIND = re.compile(r"([A-Za-z_]\w*)\s*=\s*[\"']\$\{([A-Za-z_]\w*)\}[\"']")
_PY_VAR_WRITE = re.compile(
    r"(?:to_csv|savefig|write_text|write_bytes)\(\s*([A-Za-z_]\w*)"
    r"|\bopen\(\s*([A-Za-z_]\w*)\s*,\s*[\"'][wa]"
    r"|outfile\s*=\s*([A-Za-z_]\w*)"
    r"|\bto_csv\(\s*([A-Za-z_]\w*)"
)


def extract_heredoc_var_writes(
    path: Path, config: dict[str, str], results_dir: str
) -> list[Write]:
    """Heredoc writes routed through a Python variable bound to a shell path."""
    text = path.read_text(errors="replace")
    lines = logical_lines(text)
    table = build_symbol_table(lines, config)
    stage = path.name

    bindings: dict[str, str] = {}
    for match in _PY_SHELL_BIND.finditer(text):
        pyvar, shellvar = match.group(1), match.group(2)
        if shellvar in table:
            bindings[pyvar] = expand(table[shellvar], table)

    writes: list[Write] = []
    for match in _PY_VAR_WRITE.finditer(text):
        pyvar = next((g for g in match.groups() if g), None)
        if not pyvar:
            continue
        value = bindings.get(pyvar)
        if value and value.startswith(results_dir):
            lineno = text[: match.start()].count("\n") + 1
            writes.append(Write(stage, lineno, value.rstrip("/"), "heredoc-var"))
    return writes
