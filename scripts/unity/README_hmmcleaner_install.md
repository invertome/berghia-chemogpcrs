# HmmCleaner.pl install on Unity HPC

> **DEPRECATED (May 2026).** HmmCleaner is no longer in the active
> Berghia chemoGPCRs pipeline. The Perl XS dependency stack
> (`Class::XSAccessor`, `Moose`, `Bio::MUST::*`) refused to build under the
> `berghia-gpcr` conda env's Perl ABI on Unity HPC (no sudo, no system
> headers, no bioconda recipe, no published container — see the
> post-mortem table below). The alignment-cleanup step has been replaced
> by a PREQUAL + alignment-ensemble + CLOAK + TAPER stack; see
> `scripts/unity/install_alignment_filters.sh` and
> `docs/installation_hpc.md` for the replacement. `RUN_HMMCLEANER`
> defaults to `0` going forward and the install script below is retained
> only for back-compat / reproducibility.

## TL;DR — what to run

```bash
# from the repo root on Unity:
sbatch scripts/unity/install_hmmcleaner.sh
```

After SLURM completes, the binary `HmmCleaner.pl` is on `$PATH` whenever the
`berghia-gpcr` conda env is activated, and `${WORKDIR}/.hmmcleaner_installed`
exists as a success flag. The main pipeline auto-detects the binary; set
`RUN_HMMCLEANER=1` in `config.local.sh` if you want to be explicit.

## What was tried before and why it failed

| Attempt | Outcome | Root cause |
|---|---|---|
| `cpanm Bio::MUST::Apps::HmmCleaner` | fails at Bio::MUST::Drivers test phase | `t/exo_aligned.t` hangs / fails on exonerate sanity test; cpan testers show **zero** successes ([PhyloSuite #67](https://github.com/dongzhang0725/PhyloSuite/issues/67)) |
| `cpanm --installdeps Bio::MUST::Apps::HmmCleaner` | same failure | dep resolver still descends into Bio::MUST::Drivers tests |
| `conda install bioconda::perl-bio-must-apps-hmmcleaner` (job 56714199) | "completed" but binary missing | **the package does not exist on bioconda** — no `perl-bio-must-*` recipe has ever been published |
| `apptainer pull <somewhere>/hmmcleaner` | N/A — no image to pull | searched quay.io/biocontainers, depot.galaxyproject.org/singularity, DockerHub — **no published container** for HmmCleaner / Bio-MUST anywhere |

## Chosen approach: ordered cpanm with `--notest`, system bio-tools via conda

The install script does, in order:

1. Activates `berghia-gpcr` and confirms `perl` / `cpanm` paths.
2. `conda install -c bioconda exonerate hmmer blast` — these are the system
   tools `Bio::MUST::Drivers` shells out to at runtime. Installing them via
   conda is far more reliable than letting CPAN search for them, and it
   eliminates the runtime "exonerate not found" failures the CPAN tests
   were trying (poorly) to detect.
3. Walks `cpanm_dep_order.txt` topologically:
   - Tier 1: low-level Perl infra (Test::Most, Test::Files, Try::Tiny, …)
   - Tier 2: Moose stack (Moose, MooseX::*, namespace::autoclean, …)
   - Tier 3: utility / templating (Modern::Perl, Path::Class, Tie::IxHash,
     Template, JSON, YAML::XS, LWP, …)
   - Tier 4: BioPerl + Bio::LITE::Taxonomy + Bio::Phylo (Bio::MUST::Core
     uses these at runtime; **not** "BioPerl" the umbrella —
     [Bio::MUST::Core POD explicitly says so](https://metacpan.org/pod/Bio::MUST::Core))
   - Tier 5: Bio::FastParsers → Bio::MUST::Core → Bio::MUST::Drivers →
     Bio::MUST::Apps::HmmCleaner
4. Every cpanm invocation uses `--notest --no-man-pages` (set via
   `PERL_CPANM_OPT`). This is what unblocks the `exo_aligned.t` failure:
   the test was the only thing failing, and skipping it is safe because
   we substitute our own end-to-end verification (`HmmCleaner.pl --help`)
   plus runtime use in stage 04 of the pipeline.
5. Verification: `command -v HmmCleaner.pl` AND `HmmCleaner.pl --help`
   must both succeed before the success flag is touched.

The script is **idempotent** — re-running it after success is a no-op
(checks the flag and the working binary).

## Sources / citations

- **Bio::MUST::Apps::HmmCleaner distribution + Requires list**
  https://metacpan.org/dist/Bio-MUST-Apps-HmmCleaner
- **Bio::MUST::Core distribution + Requires list (latest 0.250380)**
  https://metacpan.org/release/DBAURAIN/Bio-MUST-Core-0.250380
- **Bio::MUST::Core POD — "depends on Bio::LITE::Taxonomy and Bio::Phylo, not BioPerl per se"**
  https://metacpan.org/pod/Bio::MUST::Core
- **PhyloSuite issue #67 — documented "zero successes" on cpan testers**
  https://github.com/dongzhang0725/PhyloSuite/issues/67
- **PhyloSuite plugin install docs — recommends Perl 5.18 over 5.26+**
  https://dongzhang0725.github.io/PhyloSuite-demo/how-to-configure-plugins/
- **HmmCleaner upstream description**
  https://bioinformaticshome.com/tools/msa/descriptions/HmmCleaner.html
- **cpanm `--notest` / `--force` semantics**
  https://metacpan.org/dist/App-cpanminus/view/bin/cpanm

## Residual risks and what to do if this fails

1. **`Bio::Perl` install hangs.** If tier-4 wedges on Bio::Perl's network
   tests, run `cpanm --notest Bio::Perl` interactively first (the
   `PERL_CPANM_OPT` should already do this, but BioPerl's `Test::Database`
   tries DNS lookups).
2. **`Bio::MUST::Drivers` final install still fails with `--notest`.**
   The script then retries `cpanm --force Bio::MUST::Apps::HmmCleaner`.
   If even `--force` fails, the binary may still install (Makefile install
   target is independent of test target); check
   `$CONDA_PREFIX/bin/HmmCleaner.pl` directly.
3. **Perl ABI mismatch.** The `berghia-gpcr` env uses bioconda's perl
   (which is 5.32.x as of May 2026). Bio::MUST::Drivers tests were filed
   when the upstream author developed against 5.18 — runtime works on
   5.32, only tests are flaky. We sidestep this with `--notest`.
4. **If we ever DO get a published Apptainer image** for Bio-MUST-Apps,
   replace this whole script with a 3-line `apptainer pull` + wrapper.
   None exists today (May 2026).
