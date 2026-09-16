# CRAN submission comments — admixr2 0.4.1

## Summary

This is a major update from the current CRAN version (0.2.0); versions 0.3.0
and 0.4.0 were released on GitHub only. The most relevant changes for
reviewers:

- **New default: `covMethod = "r,s"`.** An Asymptotically Distribution-Free
  (ADF) sandwich covariance correction (Browne, 1984) that scores the model's
  summary statistics against their own asymptotic sampling law instead of
  assuming multivariate normality. This is now the default for all four
  estimators, so **reported standard errors change** for any existing script
  that does not name `covMethod` explicitly; point estimates and objective
  function values are unaffected (the covariance is computed after
  convergence and never feeds back into the fit). `covMethod = "r"` restores
  the previous behaviour.
- **Transform-both-sides endpoints (`boxCox`, `yeoJohnson`, `logitNorm`,
  `probitNorm`) are now composed exactly** rather than via a second-order
  approximation that did not converge with node count. This **changes fit
  results** (point estimates and OFV) for these four residual families on the
  `adgh`/`admc` estimators; `adfo` is unaffected.
- New estimator features: multi-compartment (multi-output) fitting, parallel
  restarts on `mirai` daemons, an Iterative Reweighting Monte Carlo estimator
  (`est = "adirmc"`), Student-t and ordinal-categorical residual error models,
  analytical structural-theta gradients (including second-order `linCmt()`
  sensitivities), and finite-difference steps measured per parameter (Shi
  2021) in place of a fixed step.
- Numerous correctness fixes across the residual-error, gradient, and
  covariance-reporting code, none of which involve a package API change.

See `NEWS.md` for the complete, itemised list of changes across all three
releases.

## Test environments

- Local: Debian 13 (trixie), R 4.5.0; Windows 11, R 4.5.3
- GitHub Actions (PR #125 and #126, all green): ubuntu-latest (R release,
  R devel), windows-latest (R release), macOS-latest (R release), plus a
  separate integration-test job against the current CRAN dependency stack
- [MAINTAINER TODO before submission: run win-builder (R-devel and R-release)
  and/or R-hub, and update this section with the results]

## R CMD check results

0 errors | 0 warnings | 1 note

- `Skipping checking math rendering: package 'V8' unavailable`: the local
  check environment does not have the 'V8' R package installed; this is an
  environment gap of the check machine, not a package issue, and 'V8' is not
  a dependency of admixr2.

The first submission of 0.4.1 was rejected on three notes and a warning from
the CRAN incoming checks, all caused by local working files that had been
picked up by the build rather than by the package itself. All are fixed:

- hidden directories (`.scratch/`, `.codex/`) and non-standard top-level
  files (`AGENTS.md`, `dead.txt`, `meanlog`, `zrep.out`) are removed from the
  package directory, and the editor/agent configuration that remains there
  locally is now listed in `.Rbuildignore` next to the existing `.claude`
  and `CLAUDE.md` entries;
- the undeclared executable `.scratch/.../admixr2.dll` went with `.scratch/`;
- the tarball is 0.7Mb, down from 13.3Mb (`.scratch/` was the whole of it);
- the five "possibly invalid file URIs" in `inst/doc/admixr2.html` were
  relative links from the shipped vignette to vignettes that are not shipped
  (see "Vignettes" below). They now point at the package website.

Earlier, a local `R CMD check --as-cran` found two further real issues, both
fixed at the time: a stray top-level `index.md` (pkgdown homepage source, now
in `.Rbuildignore`) and a `URL:` in DESCRIPTION that redirected (trailing
slash added). Everything substantive -- installation, examples (including
`--run-donttest`), the full test suite, R/Rd consistency, and the vignette
build/re-build -- was independently confirmed clean.

## Reverse dependencies

`nlmixr2` lists `admixr2` in `Suggests` (a soft, ecosystem-convention link --
admixr2 in turn Suggests `nlmixr2` for its own examples/tests). No exported
admixr2 function signature changed in this release, only internal behaviour
and defaults, so no impact on `nlmixr2` is expected; a full `revdepcheck` run
was not performed as part of this submission.

## Notes on dependencies

- `nlmixr2est` and `rxode2` are on CRAN and provide the model specification
  and ODE-solving infrastructure this package integrates with. Minimum versions
  are declared in `Imports`.
- `nlmixr2` (Suggests) is used in examples and tests; a minimum version is
  declared.
- `mirai` (Suggests) is used only when `workers > 1` is set in the control
  object; it is loaded conditionally via `requireNamespace()`.
- `memuse` (Suggests) is never called by this package. It is declared because
  `rxode2::rxSolve()` estimates free RAM on every call, and its fallback path
  shells out to the macOS-only `vm_stat` when `memuse` is absent -- a spawned
  process per solve, on every platform where that command does not exist. With
  `memuse` installed the fallback is never reached. See `?adfoControl`.
- `patchwork` (Suggests) is loaded conditionally via `requireNamespace()` for
  optional 2x2 diagnostic plot layouts.
- `knitr` (Suggests) is used for a `knit_print.admFit` S3 method registered
  at package load time.
- `expm` (Suggests) is a fallback for Hessian inversion when Cholesky and
  `solve()` both fail; loaded conditionally via `requireNamespace()`.

## Vignettes

One vignette ("Getting started with admixr2") is included and fully executed
during package building. The remaining vignettes (diagnostic plots, multiple
studies, multi-compartment models, PK/PD, aggregate-data theory, estimator
comparison, advanced usage, residual error models, data generation) are
excluded from the tarball via `.Rbuildignore` because they each require
long-running fits. They are available on the package website at
<https://leidenpharmacology.github.io/admixr2>.

## Integration tests

Integration tests (those requiring `rxode2` model compilation and ODE solving)
are skipped on CRAN via `skip_on_cran()`. All remaining tests pass with 0
failures.
