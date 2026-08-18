# ---------------------------------------------------------------------------
# validation/estimator-solve-count.R
#
# Benchmark harness for the admc / adfo / adirmc estimators.
#
# Counts rxSolve CALLS (not wall time).  Per CLAUDE.md an rxSolve call costs
# ~11 ms before it integrates anything and is flat in the number of subjects,
# so the CALL COUNT is the robust cost metric on a noisy machine.  Every call is
# attributed to the deepest admixr2 frame on the stack, so the table says WHERE
# the solves are spent, not merely how many there are.
#
# It also writes a GOLDEN file (objective, gradient, fitted parameters, and a
# grid of NLL probes) so any claimed speed-up can be checked for BIT identity
# with identical(), not all.equal().
#
# Run from the package root:
#   Rscript validation/estimator-solve-count.R baseline
#   ...make a change, recompile with pkgload...
#   Rscript validation/estimator-solve-count.R after
#   Rscript validation/estimator-solve-count.R compare baseline after
# ---------------------------------------------------------------------------

args  <- commandArgs(trailingOnly = TRUE)
mode  <- if (length(args)) args[[1]] else "baseline"
outdir <- "validation/solve-count"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- compare mode: no model work, just diff two runs -----------------------
if (identical(mode, "compare")) {
  a <- readRDS(file.path(outdir, paste0(args[[2]], ".rds")))
  b <- readRDS(file.path(outdir, paste0(args[[3]], ".rds")))
  cat("\n=== rxSolve CALL COUNTS ===\n")
  ests <- union(names(a$counts), names(b$counts))
  for (e in ests) {
    ca <- a$counts[[e]]; cb <- b$counts[[e]]
    cat(sprintf("\n-- %s : %s -> %s calls\n", e,
                if (is.null(ca)) "NA" else sum(ca),
                if (is.null(cb)) "NA" else sum(cb)))
    who <- union(names(ca), names(cb))
    for (w in who)
      cat(sprintf("     %-28s %6s -> %6s\n", w,
                  if (is.null(ca[[w]])) "." else ca[[w]],
                  if (is.null(cb[[w]])) "." else cb[[w]]))
  }
  cat("\n=== BIT IDENTITY ===\n")
  for (k in names(a$golden)) {
    ok <- identical(a$golden[[k]], b$golden[[k]])
    cat(sprintf("  %-34s %s\n", k, if (ok) "IDENTICAL" else "*** DIFFERS ***"))
    if (!ok) {
      x <- a$golden[[k]]; y <- b$golden[[k]]
      if (is.numeric(x) && is.numeric(y) && length(x) == length(y)) {
        d <- max(abs(x - y) / pmax(abs(x), 1e-300))
        cat(sprintf("      max rel diff %.3e\n", d))
        cat(sprintf("      a: %s\n      b: %s\n",
                    paste(sprintf("%.17g", x), collapse = " "),
                    paste(sprintf("%.17g", y), collapse = " ")))
      }
    }
  }
  cat("\n=== WALL TIME (indicative only) ===\n")
  for (e in ests)
    cat(sprintf("  %-10s %8.2fs -> %8.2fs\n", e,
                a$times[[e]] %||% NA_real_, b$times[[e]] %||% NA_real_))
  quit(save = "no")
}

suppressMessages({
  library(rxode2)
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
})
`%||%` <- function(a, b) if (is.null(a)) b else a

stopifnot(requireNamespace("memuse", quietly = TRUE))  # else every rxSolve
                                                       # spawns a doomed process

# ---------------------------------------------------------------------------
# Instrumentation: count rxSolve calls, attributed to the deepest admixr2 frame
# ---------------------------------------------------------------------------
.SC <- new.env(parent = emptyenv())
.SC$on <- FALSE
.SC$tab <- integer(0)

.sc_frame <- function() {
  cl  <- sys.calls()
  nms <- vapply(cl, function(x) {
    f <- x[[1L]]
    if (is.name(f)) return(as.character(f))
    if (is.call(f) && length(f) == 3L &&
        as.character(f[[1L]]) %in% c("::", ":::")) return(as.character(f[[3L]]))
    "<anon>"
  }, character(1))
  hit <- grep("^\\.(adm|adfo|adgh|adirmc)", nms, value = TRUE)
  if (length(hit)) hit[[length(hit)]] else "<other>"
}

sc_start <- function() { .SC$tab <- integer(0); .SC$depth <- 0L; .SC$on <- TRUE }
sc_stop  <- function() { .SC$on <- FALSE; sort(.SC$tab, decreasing = TRUE) }

# The tracer body is evaluated in rxSolve's own frame; lookup then walks the
# rxode2 namespace -> imports -> base -> globalenv, which is where .SC lives.
#
# rxSolve re-enters itself (the generic normalises its arguments and calls back
# in), so a naive tracer counts one user-visible solve twice.  A depth counter,
# incremented on entry and decremented by the `exit` tracer, counts only the
# OUTERMOST entry -- that is the call whose ~11 ms fixed cost is paid.
sc_install <- function() {
  suppressMessages(trace(
    "rxSolve", where = asNamespace("rxode2"), print = FALSE,
    tracer = quote({
      if (isTRUE(.SC$on)) {
        if (.SC$depth == 0L) {
          .w <- .sc_frame()
          .SC$tab[.w] <- (if (is.na(.SC$tab[.w])) 0L else .SC$tab[.w]) + 1L
        }
        .SC$depth <- .SC$depth + 1L
      }
    }),
    exit = quote({ if (isTRUE(.SC$on)) .SC$depth <- .SC$depth - 1L })))
}

sc_install()

# self-check: the instrumentation must actually fire, or every number is a lie
local({
  m <- rxode2::rxode2({ d/dt(x) <- -k * x })
  sc_start()
  invisible(rxode2::rxSolve(m, params = c(k = 1), events = rxode2::et(1:3),
                            inits = c(x = 1), cores = 1L))
  n <- sum(sc_stop())
  if (n != 1L) stop("solve counter did not fire (got ", n, ")")
  cat("instrumentation self-check OK\n")
})

# ---------------------------------------------------------------------------
# Representative problem: 2-cmt ODE, 2 mu-referenced etas, 2 unpaired thetas,
# proportional + additive residual, TWO studies with different dosing.
# ---------------------------------------------------------------------------
mod_fn <- function() {
  ini({
    tcl  <- log(4)
    tv   <- log(30)
    tq   <- log(2)     # unpaired: no eta -> exercises THETA_j_ sens columns
    tvp  <- log(40)    # unpaired
    add.err  <- 0.05
    prop.err <- 0.1
    eta.cl ~ 0.09
    eta.v  ~ 0.04
  })
  model({
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    q  <- exp(tq)
    vp <- exp(tvp)
    d/dt(central)  <-  q / vp * periph - q / v * central - cl / v * central
    d/dt(periph)   <- -q / vp * periph + q / v * central
    cp <- central / v
    cp ~ add(add.err) + prop(prop.err)
  })
}

set.seed(11)
times1 <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
times2 <- c(0.5, 1, 3, 6, 12, 24)

# synthetic aggregate data from datagen-free closed form is not available for a
# 2-cmt model, so simulate it once with rxode2 (outside the counter).
sim_studies <- function() {
  ui <- suppressMessages(rxode2::rxode2(mod_fn))
  m  <- admixr2:::.admLoadModel(ui)
  mk <- function(times, dose, n) {
    ev <- rxode2::et(amt = dose) |> rxode2::et(times)
    om <- lotri::lotri(eta.cl ~ 0.09, eta.v ~ 0.04)
    s  <- rxode2::rxSolve(m, params = c(tcl = log(4), tv = log(30), tq = log(2),
                                        tvp = log(40)),
                          omega = om, events = ev, nSub = n,
                          returnType = "data.frame", cores = 1L)
    mat <- matrix(s$cp, nrow = n, byrow = TRUE)
    mat <- mat * (1 + 0.1 * matrix(rnorm(length(mat)), nrow = n)) +
      0.05 * matrix(rnorm(length(mat)), nrow = n)
    list(E = colMeans(mat), V = stats::cov.wt(mat, method = "ML")$cov,
         n = as.integer(n), times = times, ev = rxode2::et(amt = dose))
  }
  list(low = mk(times1, 100, 200L), high = mk(times2, 300, 150L))
}
studies <- sim_studies()

dat <- admData()

# ---------------------------------------------------------------------------
# Runs
# ---------------------------------------------------------------------------
counts <- list(); times <- list(); golden <- list()

run <- function(name, est, control) {
  cat(sprintf("\n### %s\n", name)); flush.console()
  sc_start()
  t0 <- proc.time()[["elapsed"]]
  fit <- suppressMessages(suppressWarnings(
    nlmixr2est::nlmixr2(mod_fn, dat, est = est, control = control)))
  el <- proc.time()[["elapsed"]] - t0
  counts[[name]] <<- sc_stop()
  times[[name]]  <<- el
  golden[[paste0(name, ".objective")]] <<- fit$objective
  golden[[paste0(name, ".theta")]]     <<- unname(fit$theta)
  cv <- tryCatch(unname(as.vector(fit$cov)), error = function(e) NULL)
  golden[[paste0(name, ".cov")]] <<- cv
  cat(sprintf("  %d rxSolve calls, %.1fs, objective %.17g\n",
              sum(counts[[name]]), el, fit$objective))
  print(counts[[name]])
  invisible(fit)
}

run("admc", "admc",
    admControl(studies = studies, n_sim = 300L, maxeval = 25L, print = 0L,
               cores = 1L, covMethod = "r", cov_n_sim = 300L))

run("adfo", "adfo",
    adfoControl(studies = studies, maxeval = 25L, print = 0L,
                cores = 1L, covMethod = "r"))

run("adirmc", "adirmc",
    adirmcControl(studies = studies, n_sim = 300L, outer_iter = 3L,
                  maxeval = 40L, print = 0L, cores = 1L, covMethod = "r",
                  cov_n_sim = 300L))

# ---------------------------------------------------------------------------
# Fine-grained probes: NLL and gradient at a fixed parameter grid.
# These are the bit-identity witnesses -- a full fit's objective can hide a
# change inside a path the optimiser never reached.
# ---------------------------------------------------------------------------
probe <- local({
  ui    <- suppressMessages(rxode2::rxode2(mod_fn))
  pinfo <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  sens  <- tryCatch(admixr2:::.admLoadSensModel(ui), error = function(e) NULL)
  rxMod <- admixr2:::.admLoadModel(ui)
  rxode2::rxLoad(rxMod)
  st <- lapply(seq_along(studies), function(i)
    admixr2:::.admNormaliseStudy(studies[[i]], names(studies)[i]))
  names(st) <- names(studies)
  st <- lapply(st, function(s) { s$ev_full <- s$ev |> rxode2::et(s$times); s })
  pinfo$sim_cache_file <- admixr2:::.admModelCacheFile(ui)
  list(ui = ui, pinfo = pinfo, sens = sens, rxMod = rxMod, studies = st)
})

pv <- admixr2:::.admBuildOptVec(probe$pinfo)
p0 <- pv$p0
set.seed(4)
pgrid <- list(p0, p0 + 0.05, p0 - 0.03, p0 * 1.02)

nsim <- 300L
zl <- admixr2:::.admMakeZ(nsim, probe$pinfo, length(probe$studies), "sobol")
pl <- admixr2:::.admMakeParamsList(nsim, probe$pinfo, length(probe$studies))

sc_start()
golden$probe.admNLL <- vapply(pgrid, function(p)
  admixr2:::.admNLL(p, probe$pinfo, probe$studies, zl, probe$rxMod, "cp", pl,
                    cores = 1L), double(1))
counts$probe.admNLL <- sc_stop()

sc_start()
golden$probe.admGrad <- unlist(lapply(pgrid, function(p)
  unname(admixr2:::.admGrad(p, probe$pinfo, probe$studies, zl, probe$rxMod, "cp",
                            pl, cores = 1L, h = 1e-4, sensModel = probe$sens,
                            use_central = TRUE))))
counts$probe.admGrad <- sc_stop()

sc_start()
golden$probe.admNLLBatch <- admixr2:::.admNLLBatch(
  pgrid, probe$pinfo, probe$studies, zl, probe$rxMod, "cp", pl, cores = 1L)
counts$probe.admNLLBatch <- sc_stop()

sc_start()
golden$probe.admGradBatch <- unname(unlist(admixr2:::.admGradBatch(
  pgrid, probe$pinfo, probe$studies, zl, probe$rxMod, "cp", pl, cores = 1L,
  h = 1e-4, sensModel = probe$sens)))
counts$probe.admGradBatch <- sc_stop()

sens2 <- tryCatch(admixr2:::.admLoadSensModel(probe$ui, order = 2L),
                  error = function(e) NULL)
rxode2::rxLoad(probe$rxMod)

pl1 <- admixr2:::.admMakeParamsList(1L, probe$pinfo, length(probe$studies))

sc_start()
golden$probe.adfoNLL <- vapply(pgrid, function(p)
  admixr2:::.adfoNLL(p, probe$pinfo, probe$studies, sens2, probe$rxMod, "cp",
                     pl1, cores = 1L), double(1))
counts$probe.adfoNLL <- sc_stop()

sc_start()
golden$probe.adfoGrad <- unlist(lapply(pgrid, function(p)
  unname(admixr2:::.adfoGrad(p, probe$pinfo, probe$studies, sens2, probe$rxMod,
                             "cp", pl1, cores = 1L, grad_h = 1e-4))))
counts$probe.adfoGrad <- sc_stop()

sc_start()
golden$probe.adfoFDGrad <- unlist(lapply(pgrid, function(p)
  unname(admixr2:::.adfoFDGrad(p, probe$pinfo, probe$studies, sens2, probe$rxMod,
                               "cp", pl1, cores = 1L, grad_h = 1e-4))))
counts$probe.adfoFDGrad <- sc_stop()

sc_start()
golden$probe.adfoCalcCov.grad <- unname(as.vector(admixr2:::.adfoCalcCov(
  p0, probe$pinfo, probe$studies, sens2, probe$rxMod, "cp", pl1, cores = 1L,
  use_grad = TRUE)))
counts$probe.adfoCalcCov.grad <- sc_stop()

sc_start()
golden$probe.adfoCalcCov.nll <- unname(as.vector(admixr2:::.adfoCalcCov(
  p0, probe$pinfo, probe$studies, sens2, probe$rxMod, "cp", pl1, cores = 1L,
  use_grad = FALSE)))
counts$probe.adfoCalcCov.nll <- sc_stop()

sc_start()
golden$probe.admCalcCov <- unname(as.vector(admixr2:::.admCalcCov(
  p0, probe$pinfo, probe$studies, zl, probe$rxMod, "cp", pl, cores = 1L,
  cov_n_sim = 300L, use_grad = TRUE, sensModel = probe$sens)))
counts$probe.admCalcCov <- sc_stop()

for (nm in grep("^probe", names(counts), value = TRUE))
  cat(sprintf("%-24s %4d rxSolve  %s\n", nm, sum(counts[[nm]]),
              paste(sprintf("%s=%d", names(counts[[nm]]), counts[[nm]]),
                    collapse = " ")))

saveRDS(list(counts = counts, times = times, golden = golden),
        file.path(outdir, paste0(mode, ".rds")))
cat(sprintf("\nwrote %s\n", file.path(outdir, paste0(mode, ".rds"))))
