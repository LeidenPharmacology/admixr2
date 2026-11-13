# CRAN submission comments -- admixr2 0.4.2

## Summary

Update from the current CRAN version (0.4.1). Main changes (see `NEWS.md`):

- **RcppEigen 5 compatibility.** An importance-weight computation depended on
  Eigen 3.4's handling of non-finite values; it is now explicit and gives the
  same result under Eigen 3.4 and Eigen 5.
- New: inter-occasion variability for `adgh`/`admc` (needs nlmixr2est >= 7.2.1,
  skipped otherwise), fixed omega is held rather than refused, simulated trials
  and covariate-dependent dosing in `admStudy()`, `datagen()` on study objects.
- Several silent-wrong-answer bugs fixed (unsorted observation times, lognormal
  covariate margins, random number stream reset); bad study input is now an error.

## Test environments

- Local: Windows 11, R 4.5.3, `R CMD check --as-cran` on a tarball built from a
  clean checkout
- GitHub Actions: ubuntu-latest (R release, R devel), windows-latest (R release),
  macOS-latest (R release), plus a job against the current CRAN dependency stack

## R CMD check results

0 errors | 0 warnings | 2 notes

- `unable to verify current time`: the check machine could not reach a time
  server.
- `Skipping checking math rendering: package 'V8' unavailable`: environment gap
  of the check machine; V8 is not a dependency.

## Reverse dependencies

`nlmixr2` lists `admixr2` in `Suggests`. No exported function signature was
removed or reordered (new arguments were added last), so no impact is expected.
