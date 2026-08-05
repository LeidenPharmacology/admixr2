# Probe D -- WHAT is retained?
#
# Probe B settled where the devel/release difference lives: it is the R HEAP,
# not native memory. After a full gc(), devel holds ~5.3 GB of LIVE R objects
# 27 files in, where release holds 0.57 GB after all 58. Same package versions,
# same machine. Per file, the worst offenders are stark:
#
#   test-integration-cache-hit.R   release 5.0s / +40 MB   devel 117.7s / +2014 MB
#   test-integration-gill.R        release small           devel 238.2s / +3020 MB
#   test-integration-edge-cases.R  release small           devel 133.5s / +1010 MB
#
# The slowness follows from the heap: a multi-GB live heap makes every
# subsequent gc() expensive, which is why devel runs the suite ~6x slower.
#
# So: run the worst files, then walk every reachable environment and report the
# biggest holders. Whatever differs between the two R versions has to show up
# as a binding that is large on devel and small on release.
#
# Run from the repository root.

source(".github/mem-probe/meminfo.R")

Sys.setenv(NOT_CRAN = "true")
.memenv()
suppressMessages(library(admixr2))
suppressMessages(library(testthat))

.sz <- function(x) tryCatch(as.numeric(utils::object.size(x)) / 1024^2,
                            error = function(e) NA_real_)

# Sum the bindings of an environment WITHOUT recursing into child environments
# (those are reported on their own line, so nothing is double counted).
.envsz <- function(e) {
  nms <- tryCatch(ls(e, all.names = TRUE), error = function(x) character(0))
  if (!length(nms)) return(c(mb = 0, n = 0))
  tot <- 0
  for (n in nms) {
    o <- tryCatch(get(n, envir = e, inherits = FALSE), error = function(x) NULL)
    if (is.environment(o)) next
    v <- .sz(o)
    if (!is.na(v)) tot <- tot + v
  }
  c(mb = tot, n = length(nms))
}


# base::memory.profile() counts LIVE SEXPs by type. This is the instrument
# object.size() cannot be: it does not care who references what, so nothing can
# hide in a closure environment. If devel retains something release does not, a
# type count has to diverge -- and WHICH type names the culprit (LANGSXP/EXPRSXP
# => parsed code and srcrefs; ENVSXP => environments; PROMSXP => unforced
# promises; REALSXP => numeric data).
.prof0 <- NULL
.profdiff <- function() {
  p <- tryCatch(base::memory.profile(), error = function(e) NULL)
  if (is.null(p)) return(invisible())
  if (is.null(.prof0)) { .prof0 <<- p; return(invisible()) }
  d <- p - .prof0[names(p)]
  d <- d[order(-abs(d))]
  d <- d[abs(d) > 1000][seq_len(min(8L, sum(abs(d) > 1000)))]
  if (!length(d)) return(invisible())
  cat("  SEXP delta vs baseline: ",
      paste(sprintf("%s %+d", names(d), as.integer(d)), collapse = "  "), "\n")
}

.report <- function(label) {
  gc(full = TRUE)
  m <- .meminfo()
  cat(sprintf("\n===== %s : heap %.1f MB (ncells %.1f / vcells %.1f) | rss %.1f MB =====\n",
              label, m[["heap"]], m[["ncells"]], m[["vcells"]], m[["rss"]]))
  .profdiff()

  rows <- list()
  add <- function(nm, e) {
    s <- .envsz(e)
    if (!is.na(s[["mb"]]) && s[["mb"]] > 5)
      rows[[length(rows) + 1L]] <<- data.frame(where = nm, mb = s[["mb"]],
                                               n = s[["n"]])
  }

  add("globalenv", globalenv())
  for (p in loadedNamespaces()) {
    ns <- tryCatch(asNamespace(p), error = function(e) NULL)
    if (is.null(ns)) next
    add(paste0("ns:", p), ns)
    # Package-level environments are where registries and memo caches live --
    # rxode2's .rxModels, admixr2's .adm_sens_env, and anything similar.
    for (n in tryCatch(ls(ns, all.names = TRUE), error = function(e) character(0))) {
      o <- tryCatch(get(n, envir = ns, inherits = FALSE), error = function(e) NULL)
      if (is.environment(o)) add(paste0("ns:", p, "$", n), o)
    }
  }

  if (!length(rows)) { cat("  (nothing over 5 MB)\n"); return(invisible()) }
  d <- do.call(rbind, rows)
  d <- d[order(-d$mb), ][seq_len(min(20L, nrow(d))), ]
  for (i in seq_len(nrow(d)))
    cat(sprintf("  %-52s %9.1f MB  (%d bindings)\n", d$where[i], d$mb[i], d$n[i]))
  flush(stdout())
}

owd <- setwd("tests/testthat")
on.exit(setwd(owd), add = TRUE)
for (h in list.files(".", "^helper.*\\.R$")) sys.source(h, envir = globalenv())
.report("baseline (helpers only)")

# The files that carried the growth on devel in the run-3 curve, in suite order.
# Three files from a cold session only reached 1.3 GB; the 5.3 GB figure is
# CUMULATIVE, so the sequence has to be long enough to show the accumulation.
for (f in c("test-adgh-nodes.R", "test-datagen.R", "test-integration-beta.R",
            "test-integration-cache-hit.R", "test-integration-cold-session.R",
            "test-integration-datagen-roundtrip.R",
            "test-integration-edge-cases.R", "test-integration-gill.R")) {
  if (!file.exists(f)) { cat(sprintf("# missing: %s\n", f)); next }
  t <- system.time(tryCatch(testthat::test_file(f, reporter = "silent"),
                            error = function(e)
                              cat(sprintf("# ERROR %s: %s\n", f,
                                          conditionMessage(e)))))[["elapsed"]]
  .report(sprintf("%s (%.1fs)", f, t))
}
