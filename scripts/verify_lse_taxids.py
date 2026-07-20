#!/usr/bin/env python3
# verify_lse_taxids.py
# Purpose: Verify the configured LSE clade taxids against the NCBI taxonomy API,
#          so a wrong taxid can never again pass silently.
# Inputs:  LSE_AEOLID_TAXID / LSE_NUDIBRANCH_TAXID / LSE_GASTROPOD_TAXID (env)
# Outputs: A per-variable report on stdout; exit 0 (ok) / 1 (wrong) / 2 (unchecked)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.
#
# WHY THIS EXISTS
# ---------------
# All three LSE clade taxids shipped wrong. 644 was *Aeromonas hydrophila* and
# 54397 a *Lamellibrachia* sp. endosymbiont -- both bacteria -- and 13843 did
# not resolve at all. lse_refine.py classifies with
# `if LSE_AEOLID_TAXID in common_lineage`, and a bacterium's taxid can never
# appear in a mollusc's lineage, so every orthogroup fell through every branch.
# LSE classification, the scientific core of the project, silently produced
# nothing and nothing complained. validate_config.sh only checked `taxid > 0`.
#
# WHAT IS CHECKED
# ---------------
# Two things per taxid, because either alone is too weak:
#   1. The scientific name NCBI returns is the expected clade.
#   2. The taxid is a genuine ancestor of Berghia stephanieae (1287507).
# (2) is what actually catches the original defect: a bacterial taxid resolves
# perfectly well, so a name check alone is only accidentally protective. (1)
# catches the subtler case of a real Berghia ancestor at the wrong rank, e.g.
# Aeolidiidae (family, 195871) where Aeolidioidea (superfamily, 71481) is meant.
#
# OFFLINE BEHAVIOUR
# -----------------
# The pipeline runs on compute nodes that frequently have no outbound network,
# so this guard must never turn a connectivity problem into a false accusation
# or a blocked run. Outcomes are therefore three-valued, not two:
#
#   exit 0  OK       - checked, correct
#   exit 1  MISMATCH - checked, WRONG. A real error. Blocks.
#   exit 2  UNKNOWN  - could NOT be checked. A warning. Does not block.
#
# "NCBI says this taxid has no record" is a MISMATCH, not an UNKNOWN: that is a
# definitive answer from the API, and it is exactly what 13843 returns. Only
# transport-level failures (DNS, timeout, connection refused, HTTP 5xx/429,
# unparseable payload) produce UNKNOWN. Being unable to check is reported as
# being unable to check, and never as being wrong.

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

BERGHIA_TAXID = 1287507

STATUS_OK = "ok"
STATUS_MISMATCH = "mismatch"
STATUS_UNKNOWN = "unknown"

# The authoritative expectation: variable -> (taxid, scientific name, rank).
# Verified against the NCBI taxonomy API; all three are ancestors of
# Berghia stephanieae (1287507). Do not edit these from memory -- re-resolve
# them against NCBI first.
EXPECTED_CLADES = {
    "LSE_AEOLID_TAXID": (71481, "Aeolidioidea", "superfamily"),
    "LSE_NUDIBRANCH_TAXID": (70849, "Nudibranchia", "order"),
    "LSE_GASTROPOD_TAXID": (6448, "Gastropoda", "class"),
}

# Convenience view: variable -> taxid, shaped like a configuration mapping.
EXPECTED_CLADES_AS_CONFIG = {k: v[0] for k, v in EXPECTED_CLADES.items()}


class NetworkUnavailable(Exception):
    """The API could not be reached or gave an unusable answer.

    Raised ONLY for transport-level problems. A well-formed reply saying a
    taxid does not exist is a real answer and must not raise this.
    """


@dataclass
class Result:
    variable: str
    taxid: int
    expected_name: str
    status: str
    detail: str


class NCBIFetcher:
    """Minimal, rate-limited NCBI E-utilities client.

    NCBI permits 3 requests/second without an API key and 10 with one, so
    requests are spaced by `min_interval` and every call carries a timeout.
    The whole guard costs two requests: one batched esummary for all the
    taxids, one efetch for the Berghia lineage.
    """

    def __init__(self, timeout: float = 15.0, min_interval: Optional[float] = None,
                 api_key: Optional[str] = None, tool: str = "berghia-chemogpcrs",
                 email: Optional[str] = None):
        self.timeout = float(timeout)
        self.api_key = api_key or os.getenv("NCBI_API_KEY") or None
        # 10 req/s with a key, 3 without; stay a little under either ceiling.
        default_interval = 0.11 if self.api_key else 0.34
        self.min_interval = default_interval if min_interval is None else float(min_interval)
        self.tool = tool
        self.email = email or os.getenv("NCBI_EMAIL") or None
        self._last_request = 0.0

    # -- transport ----------------------------------------------------------

    def _throttle(self) -> None:
        elapsed = time.monotonic() - self._last_request
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)

    def _params(self, **kwargs) -> str:
        params = {"db": "taxonomy", "tool": self.tool}
        if self.email:
            params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key
        params.update(kwargs)
        return urllib.parse.urlencode(params)

    def _get(self, endpoint: str, retries: int = 2, **params) -> bytes:
        url = f"{EUTILS}/{endpoint}?{self._params(**params)}"
        last_exc: Optional[Exception] = None
        for attempt in range(retries + 1):
            self._throttle()
            try:
                with urllib.request.urlopen(url, timeout=self.timeout) as resp:
                    self._last_request = time.monotonic()
                    return resp.read()
            except urllib.error.HTTPError as exc:
                self._last_request = time.monotonic()
                # 429/5xx are transient; 4xx otherwise is not worth retrying.
                if exc.code == 429 or exc.code >= 500:
                    last_exc = exc
                else:
                    raise NetworkUnavailable(f"HTTP {exc.code} from {endpoint}") from exc
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                self._last_request = time.monotonic()
                last_exc = exc
            if attempt < retries:
                time.sleep(self.min_interval * (2 ** attempt))
        raise NetworkUnavailable(f"{endpoint} unreachable: {last_exc}")

    # -- queries ------------------------------------------------------------

    def summaries(self, taxids: Sequence[int]) -> Dict[int, Optional[dict]]:
        """Batch-resolve taxids. A value of None means NCBI has no such record."""
        if os.getenv("LSE_VERIFY_FORCE_OFFLINE"):
            raise NetworkUnavailable("forced offline via LSE_VERIFY_FORCE_OFFLINE")
        taxids = list(dict.fromkeys(int(t) for t in taxids))
        raw = self._get("esummary.fcgi", id=",".join(str(t) for t in taxids),
                        retmode="json")
        try:
            payload = json.loads(raw)
        except (ValueError, UnicodeDecodeError) as exc:
            raise NetworkUnavailable(f"unparseable esummary payload: {exc}") from exc
        result = payload.get("result")
        if not isinstance(result, dict):
            raise NetworkUnavailable("esummary payload had no result block")

        out: Dict[int, Optional[dict]] = {}
        for taxid in taxids:
            record = result.get(str(taxid))
            # A record carrying an "error" key is NCBI definitively saying the
            # id is unusable -- a real answer, so None rather than an exception.
            if not isinstance(record, dict) or "error" in record or not record.get("scientificname"):
                out[taxid] = None
            else:
                out[taxid] = record
        return out

    def lineage(self, taxid: int) -> List[int]:
        """Full ancestor taxid list for `taxid`, root first."""
        if os.getenv("LSE_VERIFY_FORCE_OFFLINE"):
            raise NetworkUnavailable("forced offline via LSE_VERIFY_FORCE_OFFLINE")
        raw = self._get("efetch.fcgi", id=str(int(taxid)), retmode="xml")
        try:
            root = ET.fromstring(raw)
        except ET.ParseError as exc:
            raise NetworkUnavailable(f"unparseable efetch payload: {exc}") from exc
        node = root.find("Taxon")
        if node is None:
            raise NetworkUnavailable(f"efetch returned no Taxon for {taxid}")
        lineage_ex = node.find("LineageEx")
        if lineage_ex is None:
            raise NetworkUnavailable(f"efetch returned no LineageEx for {taxid}")
        return [int(t.findtext("TaxId")) for t in lineage_ex.findall("Taxon")
                if t.findtext("TaxId")]


def verify_taxids(configured: Dict[str, int], fetcher) -> List[Result]:
    """Check each configured taxid for name and Berghia-ancestry.

    Never raises on network trouble: unreachable checks come back as UNKNOWN.
    """
    variables = [v for v in EXPECTED_CLADES if v in configured]
    variables += [v for v in configured if v not in EXPECTED_CLADES]

    # One batched summary request for every taxid we care about.
    try:
        records = fetcher.summaries([configured[v] for v in variables])
        summary_error = None
    except NetworkUnavailable as exc:
        records, summary_error = {}, exc

    # One lineage request, reused across all three checks.
    lineage: Optional[List[int]] = None
    lineage_error: Optional[Exception] = None
    if summary_error is None:
        try:
            lineage = fetcher.lineage(BERGHIA_TAXID)
        except NetworkUnavailable as exc:
            lineage_error = exc

    results: List[Result] = []
    for variable in variables:
        taxid = int(configured[variable])
        expected = EXPECTED_CLADES.get(variable)
        expected_name = expected[1] if expected else "<unknown clade>"

        if summary_error is not None:
            results.append(Result(
                variable, taxid, expected_name, STATUS_UNKNOWN,
                f"could not check taxid {taxid}: NCBI taxonomy API unreachable "
                f"({summary_error}). This is NOT a mismatch -- the value was "
                f"not verified either way.",
            ))
            continue

        if expected is None:
            results.append(Result(
                variable, taxid, expected_name, STATUS_UNKNOWN,
                f"could not check taxid {taxid}: no expected clade is defined "
                f"for {variable}.",
            ))
            continue

        _, expected_name, expected_rank = expected
        record = records.get(taxid)

        if record is None:
            results.append(Result(
                variable, taxid, expected_name, STATUS_MISMATCH,
                f"{variable}={taxid} has no NCBI taxonomy record; expected "
                f"{expected_name} ({expected_rank}, taxid {expected[0]}).",
            ))
            continue

        actual_name = record.get("scientificname", "")
        actual_rank = record.get("rank", "")

        if actual_name != expected_name:
            results.append(Result(
                variable, taxid, expected_name, STATUS_MISMATCH,
                f"{variable}={taxid} resolves to {actual_name!r} "
                f"(rank {actual_rank or 'unknown'}), but {expected_name} "
                f"({expected_rank}, taxid {expected[0]}) was expected.",
            ))
            continue

        if actual_rank and expected_rank and actual_rank != expected_rank:
            results.append(Result(
                variable, taxid, expected_name, STATUS_MISMATCH,
                f"{variable}={taxid} resolves to {actual_name!r} but at rank "
                f"{actual_rank!r}, not the expected {expected_rank!r}.",
            ))
            continue

        # Name is right. Now the check that actually matters.
        if lineage is None:
            results.append(Result(
                variable, taxid, expected_name, STATUS_UNKNOWN,
                f"could not check whether taxid {taxid} is an ancestor of "
                f"Berghia stephanieae: NCBI unreachable ({lineage_error}). "
                f"The name resolved correctly to {actual_name!r}.",
            ))
            continue

        if taxid not in lineage:
            results.append(Result(
                variable, taxid, expected_name, STATUS_MISMATCH,
                f"{variable}={taxid} ({actual_name}) is NOT an ancestor of "
                f"Berghia stephanieae ({BERGHIA_TAXID}). lse_refine.py tests "
                f"`taxid in common_lineage`, so this level could never match.",
            ))
            continue

        results.append(Result(
            variable, taxid, expected_name, STATUS_OK,
            f"{variable}={taxid} is {actual_name} ({actual_rank}), a verified "
            f"ancestor of Berghia stephanieae.",
        ))

    return results


def exit_code_for(results: Sequence[Result]) -> int:
    """1 if anything was checked and wrong; 2 if anything could not be checked."""
    if any(r.status == STATUS_MISMATCH for r in results):
        return 1
    if any(r.status == STATUS_UNKNOWN for r in results):
        return 2
    return 0


def configured_from_env(env=None) -> Dict[str, int]:
    env = os.environ if env is None else env
    configured: Dict[str, int] = {}
    for variable, (default, _, _) in EXPECTED_CLADES.items():
        raw = env.get(variable, "").strip()
        if not raw:
            configured[variable] = default
            continue
        try:
            configured[variable] = int(raw)
        except ValueError:
            # Non-numeric config is a definite error, but parsing is not this
            # function's job to report; -1 can never be a valid taxid and will
            # be reported as a mismatch downstream.
            configured[variable] = -1
    return configured


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Verify the configured LSE clade taxids against NCBI taxonomy.",
    )
    parser.add_argument("--timeout", type=float, default=15.0,
                        help="Per-request timeout in seconds (default: 15).")
    parser.add_argument("--min-interval", type=float, default=None,
                        help="Minimum seconds between requests (default: 0.34, "
                             "or 0.11 when NCBI_API_KEY is set).")
    parser.add_argument("--print-configured", action="store_true",
                        help="Print the taxids that would be checked, then exit 0.")
    parser.add_argument("--quiet", action="store_true",
                        help="Only report problems.")
    args = parser.parse_args(argv)

    configured = configured_from_env()

    if args.print_configured:
        for variable in EXPECTED_CLADES:
            print(f"{variable}={configured[variable]}")
        return 0

    fetcher = NCBIFetcher(timeout=args.timeout, min_interval=args.min_interval)
    results = verify_taxids(configured, fetcher)

    for r in results:
        if r.status == STATUS_OK:
            if not args.quiet:
                print(f"[ok]       {r.detail}")
        elif r.status == STATUS_MISMATCH:
            print(f"[MISMATCH] {r.detail}", file=sys.stderr)
        else:
            print(f"[unchecked] {r.detail}", file=sys.stderr)

    code = exit_code_for(results)
    if code == 1:
        print(
            "\nLSE taxid verification FAILED. Fix the values in config.sh "
            "(and the matching os.getenv defaults in scripts/lse_refine.py) "
            "before running the pipeline: an LSE level whose taxid is not a "
            "Berghia ancestor silently classifies nothing.",
            file=sys.stderr,
        )
    elif code == 2 and not args.quiet:
        print(
            "\nLSE taxids could NOT be verified (no network). This is not a "
            "failure and does not block the run; the values were simply not "
            "checked. Re-run this guard from a networked host.",
            file=sys.stderr,
        )
    return code


if __name__ == "__main__":
    sys.exit(main())
