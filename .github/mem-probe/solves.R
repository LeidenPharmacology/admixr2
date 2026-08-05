# Probe C -- is the leak per MODEL, or per SOLVE?
#
# Probe A measured ~1.7 MB of RSS per newly compiled model, of which ~0.9 MB is
# anonymous non-heap that survives gc(), rxUnloadAll() and rxClean(). But the
# suite grows ~1.2 GB while compiling only ~57 distinct models -- roughly 100 MB
# on Probe A's rate. So model compilation cannot be the dominant term, and the
# obvious other candidate is the thing the suite does thousands of times:
# rxSolve() on an ALREADY compiled model.
#
# This holds the model fixed and varies only the number of solves, which is the
# one comparison that separates the two. Reported per 100 solves, at two subject
# counts, since admixr2 solves n_sim subjects at a time (n_sim defaults to
# 20000 for admc) and a per-subject term would look very different from a
# per-call one.
#
# Run from the repository root.

source(".github/mem-probe/meminfo.R")

args   <- commandArgs(trailingOnly = TRUE)
tag    <- if (length(args)) args[[1L]] else "default"
NSOLVE <- as.integer(Sys.getenv("PROBE_SOLVES", "300"))

.memenv()
suppressMessages(library(rxode2))
.memhdr("PROBE_C", c("tag", "step", "nsub"))
.memrow("PROBE_C", c(tag, "baseline", "0"), .meminfo())

m <- rxode2::rxode2(paste0(
  "ka = 1.0;\ncl = 1.0;\nv = 20;\n",
  "d/dt(depot) = -ka*depot;\n",
  "d/dt(central) = ka*depot - (cl/v)*central;\n",
  "cp = central/v;\n"))
ev <- rxode2::eventTable()
ev$add.dosing(dose = 100)
ev$add.sampling(seq(0, 24, by = 4))
.memrow("PROBE_C", c(tag, "after_compile", "0"), .meminfo())

# Two subject counts: a per-CALL leak is flat in nsub, a per-SUBJECT leak is not.
for (nsub in c(1L, 200L)) {
  pars <- data.frame(ka = rep(1, nsub), cl = rep(1, nsub), v = rep(20, nsub))
  for (i in seq_len(NSOLVE)) {
    invisible(rxode2::rxSolve(m, ev, params = pars,
                              nDisplayProgress = .Machine$integer.max))
    if (i %% 25L == 0L)
      .memrow("PROBE_C", c(tag, sprintf("solve_%04d", i), as.character(nsub)),
              .meminfo())
  }
  gc(full = TRUE)
  .memrow("PROBE_C", c(tag, sprintf("after_gc_nsub%d", nsub), as.character(nsub)),
          .meminfo())
}

tryCatch({ rxode2::rxUnloadAll(); gc(full = TRUE) }, error = function(e) NULL)
.memrow("PROBE_C", c(tag, "after_rxUnloadAll", "0"), .meminfo())
