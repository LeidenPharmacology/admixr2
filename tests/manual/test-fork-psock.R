# Manual verification script for issue #45:
# Fork (Unix/WSL) vs PSOCK (Windows) parallel-restart paths.
#
# Run from WSL:
#   cd /mnt/c/package/admixr2
#   Rscript tests/manual/test-fork-psock.R
#
# The script covers:
#   A. Correctness: parallel (fork) NLL == sequential NLL (same seed/Sobol draws)
#   B. Benchmark:   wall-clock time for n_restarts = 4 at multiple worker counts

# Install/update admixr2 from local source so the test always reflects the
# current working tree (pak skips the reinstall when nothing has changed).
if (!requireNamespace("pak", quietly = TRUE))
  install.packages("pak")
cat("Installing admixr2 from local source...\n")
pak::local_install(".", upgrade = FALSE, ask = FALSE)

# Install CRAN packages needed by this script that may not be present.
needed <- c("future", "furrr")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  cat(sprintf("Installing missing packages: %s\n", paste(missing, collapse = ", ")))
  pak::pak(missing, ask = FALSE)
}

suppressPackageStartupMessages({
  library(admixr2)
  library(rxode2)
  library(nlmixr2est)
  library(future)
})

# ---- Model + study (mirrors helper-integration.R) ----------------------------

one_cmt_fn <- function() {
  ini({
    tcl     <- log(5)  ; label("Log CL")
    tv      <- log(20) ; label("Log V")
    add.err <- 0.1     ; label("Additive SD")
    eta.cl  ~ 0.09
    eta.v   ~ 0.04
  })
  model({
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    d/dt(central) <- -(cl / v) * central
    cp <- central / v
    cp ~ add(add.err)
  })
}

times  <- c(0.5, 1, 2, 4)
E_true <- (100 / 20) * exp(-(5 / 20) * times)
V_true <- diag((0.3 * E_true)^2)

study_spec <- list(
  s1 = list(E = E_true, V = V_true, n = 200L,
            times = times, ev = rxode2::et(amt = 100))
)

# ---- Section A: correctness (n_restarts = 2, workers = 2) --------------------

cat("\n=== A. CORRECTNESS: parallel vs sequential (n_restarts = 2) ===\n\n")

fork_supported <- future::supportsMulticore()
cat(sprintf("future::supportsMulticore(): %s\n", fork_supported))
if (!fork_supported) {
  cat("  WARNING: fork not supported on this platform — PSOCK will be used instead.\n")
  cat("  To test the fork path run this script outside RStudio on Linux/WSL.\n\n")
}

ctl_base <- admControl(
  studies  = study_spec,
  n_sim    = 300L,
  maxeval  = 10L,
  seed     = 1L,
  grad     = "sens"
)

cat("Running sequential fit (workers = 1)...\n")
t_seq <- system.time(
  fit_seq <- suppressMessages(
    nlmixr2(one_cmt_fn, admData(), est = "admc",
            control = modifyList(ctl_base, list(workers = 1L, n_restarts = 2L)))
  )
)
cat(sprintf("  Sequential NLL: %.6f  (%.1f s)\n", fit_seq$objective, t_seq["elapsed"]))

cat("Running parallel fit (workers = 2)...\n")
t_par2 <- system.time(
  fit_par2 <- suppressMessages(
    nlmixr2(one_cmt_fn, admData(), est = "admc",
            control = modifyList(ctl_base, list(workers = 2L, n_restarts = 2L)))
  )
)
cat(sprintf("  Parallel NLL:   %.6f  (%.1f s)\n", fit_par2$objective, t_par2["elapsed"]))

diff_A <- abs(fit_seq$objective - fit_par2$objective)
pass_A <- diff_A < 1e-4
cat(sprintf("\n  |NLL_par - NLL_seq| = %.2e  [%s]\n",
            diff_A, if (pass_A) "PASS" else "FAIL — objectives differ"))

# ---- Section B: benchmark (n_restarts = 4) -----------------------------------

cat("\n=== B. BENCHMARK: wall-clock time (n_restarts = 4) ===\n\n")

ctl_bench <- admControl(
  studies    = study_spec,
  n_sim      = 300L,
  maxeval    = 20L,
  seed       = 1L,
  grad       = "sens",
  n_restarts = 4L
)

worker_counts <- c(1L, 2L, 4L)
results <- vector("list", length(worker_counts))

for (i in seq_along(worker_counts)) {
  w <- worker_counts[[i]]
  cat(sprintf("  workers = %d ...", w))
  t <- system.time(
    fit <- suppressMessages(
      nlmixr2(one_cmt_fn, admData(), est = "admc",
              control = modifyList(ctl_bench, list(workers = w)))
    )
  )
  results[[i]] <- list(workers = w, elapsed = t["elapsed"], nll = fit$objective)
  cat(sprintf("  %.1f s  (NLL = %.4f)\n", t["elapsed"], fit$objective))
}

cat("\n  Summary:\n")
cat(sprintf("  %-10s  %-12s  %-10s  %-10s\n",
            "workers", "elapsed (s)", "NLL", "speedup"))
baseline <- results[[1L]]$elapsed
for (r in results) {
  cat(sprintf("  %-10d  %-12.1f  %-10.4f  %-10.2fx\n",
              r$workers, r$elapsed, r$nll, baseline / r$elapsed))
}

cat(sprintf("\nPlatform: %s\nfork supported: %s\n\n",
            .Platform$OS.type, fork_supported))
