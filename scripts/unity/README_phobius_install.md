# Phobius install on Unity HPC (tertiary TM-predictor fallback)

## TL;DR — what to run

Phobius is **academic-licensed** software and the upstream distributors
(Stockholm University SBC) require manual acceptance of the license before
they release the tarball. There is no way to fully automate this. The
two-step procedure is:

```bash
# Step 1 (one-time, on a workstation with a browser):
#   - Visit https://phobius.sbc.su.se/data.html
#   - Follow the link to the request form
#         http://software.sbc.su.se/cgi-bin/request.cgi?project=phobius
#   - Accept the academic-use license (required by the authors)
#   - Download phobius101_linux.tar.gz (~2 MB)
#   - scp it to Unity, renaming as expected by the install script:
scp phobius101_linux.tar.gz \
    unity:/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05/tools/phobius/phobius_dist.tar.gz

# Step 2 (on Unity, from the repo root):
sbatch scripts/unity/install_phobius.sh
```

After SLURM completes, the wrapper `phobius.pl` is on `$PATH` whenever the
`berghia-gpcr` conda env is activated, and
`${WORKDIR}/.phobius_installed` exists as a success flag. The pipeline's
TM-prediction chain auto-detects `phobius.pl` and uses it as the
**third-tier fallback** behind TMbed (primary) and DeepTMHMM (secondary).

## Why this is needed

The pipeline's TM-prediction chain is:

| Tier | Tool | Status as of May 2026 |
|---|---|---|
| 1 (primary)   | TMbed (Bernhofer 2022)                       | installed via `git clone BernhoferM/TMbed`, works |
| 2 (secondary) | DeepTMHMM (three sub-strategies: Apptainer / venv / system) | **not installed**, not blocking |
| 3 (tertiary)  | **Phobius (Käll et al. 2004)** — this script | **not installed**; bioconda attempt in job 56714199 used wrong package |

Phobius is a classical, very-robust HMM-based predictor that combines TM
topology + signal-peptide prediction in a single model. It is the safety
net for the rare case where both ML predictors fail or are unavailable.
For chemoreceptor GPCRs (the targets of this pipeline) we expect 7 TM
helices; Phobius produces interpretable per-residue topology output that
matches what TMbed / DeepTMHMM emit.

## Why bioconda doesn't have a real Phobius package

| Channel | Package | Reality |
|---|---|---|
| `bioconda` (main)            | none                            | **does not exist** — Stockholm Univ. has not signed off on redistribution |
| `bioconda` (`perl-bio-must-...` packages)  | none                | **wrong package** — these are alignment-cleaning Perl modules, not Phobius. This is what failed in job 56714199. |
| `predector` (anaconda.org)   | `predector::phobius`            | placeholder; install completes only after the user runs `phobius-register phobius101_linux.tar.gz` against a manually-downloaded tarball. Not a one-shot automated install. |
| `quay.io/biocontainers`      | none                            | same license restriction blocks container redistribution |
| InterProScan bundle          | "Phobius support"               | activated only after the user drops the same tarball into `interproscan/data/phobius/` (see ebi-pf-team/interproscan-docs `ActivatingLicensedAnalyses.rst`) |

So every published "Phobius install path" reduces to: user accepts license
→ downloads tarball → registers it locally. This script automates
everything *after* the manual download.

## What the install script does

1. **Idempotency check** — if `${WORKDIR}/.phobius_installed` exists and
   `phobius.pl -h` works, exit 0.
2. **Activates `berghia-gpcr`** and prints perl/conda diagnostics.
3. **Tarball presence check** — if the tarball is missing, prints the exact
   manual download instructions and exits 2 (does **not** attempt to
   auto-download — that would require bypassing the license form).
4. **Untars** to `${WORKDIR}/tools/phobius/dist/` and locates the unpacked
   `phobius/` directory robustly (handles both `phobius/` and
   `phobius101_linux/` layouts seen in different upstream versions).
5. **64-bit normalisation** — the tarball ships a legacy 32-bit `decodeanhmm`
   plus a `decodeanhmm.64bit`. On x86_64 hosts (Unity is x86_64), the 32-bit
   binary segfaults / silently returns the famous `"could not read provided
   fasta sequence at phobius.pl line 408"` error. The fix, documented at
   biostars.org/p/238642, predector#82, and funannotate#696, is to copy
   `decodeanhmm.64bit` over `decodeanhmm`. The script does this
   automatically.
6. **Sanity-checks** that `phobius.model` and `phobius.options` are present
   beside `phobius.pl` (Phobius bails without them).
7. **Wrapper** — writes `~/.local/bin/phobius.pl` as a tiny `exec`-shim
   pointing at the unpacked `phobius.pl`. Phobius's own `phobius.pl` uses
   `FindBin::Bin` to locate `decodeanhmm` and the model files, so the
   wrapper preserves that behavior by `exec`-ing the original (no symlink
   needed; symlink would also have worked but a shell shim is grep-friendly
   in logs).
8. **Smoke test** — runs `phobius.pl -short` on a bovine rhodopsin (P02699)
   FASTA shipped beside the script (`phobius_smoke_test.fasta`). Rhodopsin
   is the canonical 7TM GPCR; Phobius should report TM ≥ 1 (typically TM=7).
   Script bails with `exit 7` if TM count is below 1.
9. **Success flag** — touches `${WORKDIR}/.phobius_installed`.

## Citations

- **Phobius original paper.** Käll L, Krogh A, Sonnhammer ELL. "A combined
  transmembrane topology and signal peptide prediction method." *J Mol
  Biol* 338:1027-1036 (2004). doi:[10.1016/j.jmb.2004.03.016](https://doi.org/10.1016/j.jmb.2004.03.016).
- **Phobius web service.** Käll L, Krogh A, Sonnhammer ELL. "Advantages of
  combined transmembrane topology and signal peptide prediction—the Phobius
  web server." *Nucleic Acids Res* 35:W429-W432 (2007).
  doi:[10.1093/nar/gkm256](https://doi.org/10.1093/nar/gkm256).
- **Phobius distribution + license page** —
  https://phobius.sbc.su.se/data.html (license terms);
  http://software.sbc.su.se/cgi-bin/request.cgi?project=phobius (download form).
- **Phobius readme** — https://phobius.sbc.su.se/readme.html.

## Sources for the install workarounds

- `decodeanhmm.64bit` rename fix on x86_64 — Biostars
  [#238642](https://www.biostars.org/p/238642/), funannotate
  [#82](https://github.com/nextgenusfs/funannotate/issues/82) and
  [#696](https://github.com/nextgenusfs/funannotate/issues/696),
  predector wiki.
- License-bound install pattern (the same pattern used by InterProScan and
  predector) — ebi-pf-team/interproscan-docs
  [`ActivatingLicensedAnalyses.rst`](https://github.com/ebi-pf-team/interproscan-docs/blob/v5/docs/ActivatingLicensedAnalyses.rst).
- predector's `phobius-register` shim (anaconda.org `predector::phobius`)
  confirms the manual-tarball-then-register workflow is canonical.

## Residual risks and what to do if this fails

1. **`~/.local/bin` not on PATH inside the conda env.** Some conda envs
   prepend `${CONDA_PREFIX}/bin` and nothing else. The script exits 5 with
   a diagnostic if `command -v phobius.pl` fails right after writing the
   wrapper. Fix by adding `~/.local/bin` to PATH in `~/.bashrc` (or use
   conda's `env_vars` mechanism: `conda env config vars set
   PATH="$HOME/.local/bin:$PATH" -n berghia-gpcr`), then rerun the script.
   Alternatively, change the wrapper path in the script to
   `${CONDA_PREFIX}/bin/phobius.pl` — that's always on PATH but persists
   only as long as the env exists.
2. **Different tarball filename.** If Stockholm releases a future version
   (e.g. `phobius102_linux.tar.gz`), the script's robust `*/phobius.pl`
   glob still finds the binary. Just rename your local download to
   `phobius_dist.tar.gz` before scp'ing — the install script does not
   inspect the tarball name.
3. **Smoke test reports TM=0 on rhodopsin.** That would indicate a deeper
   problem (broken decodeanhmm, missing model file, FASTA parser bug from
   pre-64-bit-fix decodeanhmm). Inspect `${LOGDIR}/phobius_smoke_test.out`
   for raw stderr; verify
   `${WORKDIR}/tools/phobius/dist/*/phobius.model` exists and is
   non-empty; confirm `file ${WORKDIR}/tools/phobius/dist/*/decodeanhmm`
   reports `ELF 64-bit`.
4. **Tarball changes layout.** If the upstream tarball ever stops shipping
   `decodeanhmm.64bit` and only ships a single 64-bit `decodeanhmm`, the
   script falls through cleanly (the `elif` branch handles that — see
   inline comment in `install_phobius.sh`).
5. **Pipeline integration.** The pipeline's TM-prediction wrapper resolves
   `phobius.pl` via `command -v`. No config edit is required — installing
   Phobius is sufficient to enable the tertiary fallback. If you want to
   verify the chain at runtime, set `TM_PREDICTOR_PRIMARY=phobius` in
   `config.local.sh` (forces Phobius first, useful for debugging the
   wrapper without disabling TMbed).

## Files in this directory

| File | Purpose |
|---|---|
| `install_phobius.sh`            | sbatch script — idempotent install + verify |
| `README_phobius_install.md`     | this file |
| `phobius_smoke_test.fasta`      | bovine rhodopsin (UniProt P02699), used by step 8 of the install script |
