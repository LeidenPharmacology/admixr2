# Probe B -- the SUITE curve, per file, on whichever R this runs under.
#
# Reproduces the profiling that characterised the devel failure (RSS 229 -> 1457
# MB, native 112 -> 917 MB, diffuse over ~19 files) but on BOTH R versions from
# one commit, so the devel-vs-release delta can finally be read off directly --
# the original profiling only ever ran on devel, so there was never a release
# curve to compare against.
#
# Carries the same anonymous/file-backed/n_so split as Probe A, so if the delta
# is real this says in which column it lives.
#
# Emits one flushed row per test file, so an OOM kill still leaves the curve up
# to the kill -- which is the point, since devel is expected to die.
# Run from the repository root.

source(".github/mem-probe/meminfo.R")

Sys.setenv(NOT_CRAN = "true")
.memenv()
suppressMessages(library(admixr2))
suppressMessages(library(testthat))
cat(sprintf("# admixr2        : %s\n", as.character(utils::packageVersion("admixr2"))))

.memhdr("PROBE_B", c("step", "secs", "fail"))

owd <- setwd("tests/testthat")
on.exit(setwd(owd), add = TRUE)

for (h in list.files(".", "^helper.*\\.R$")) sys.source(h, envir = globalenv())
.memrow("PROBE_B", c("after_helpers", "0.0", "0"), .meminfo())

# Alphabetical -- matches how the suite actually runs, so the curve is
# comparable to the R CMD check one rather than a re-ordering of it.
for (f in sort(list.files(".", "^test-.*\\.R$"))) {
  nfail <- 0L
  t <- system.time(
    tryCatch({
      res <- testthat::test_file(f, reporter = "silent")
      nfail <- sum(vapply(res, function(x)
        sum(vapply(x$results, function(r)
          inherits(r, "expectation_failure") || inherits(r, "expectation_error"),
          logical(1))), integer(1)))
    }, error = function(e) {
      cat(sprintf("# ERROR in %s: %s\n", f, conditionMessage(e)))
      nfail <<- -1L
    })
  )[["elapsed"]]
  .memrow("PROBE_B", c(f, sprintf("%.1f", t), as.character(nfail)), .meminfo())
}
