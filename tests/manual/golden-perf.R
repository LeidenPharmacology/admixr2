# Golden-reference harness for the performance work (items 1-7).
#
# Captures every NLL / gradient / covariance / fit result that the optimisation
# work could plausibly perturb, across all four estimators, every gradient mode,
# and the single-output / multi-output / joint model shapes.
#
#   Rscript tests/manual/golden-perf.R record   # write the baseline
#   Rscript tests/manual/golden-perf.R check    # compare against the baseline
#
# The transformations in items 3-6 are meant to be EXACT (they remove redundant
# solves; they do not change any arithmetic), so the default tolerance is 0.
# Item 7's C++ changes alter floating-point summation order, so those are
# compared at 1e-10 relative.

suppressPackageStartupMessages({
  library(admixr2); library(rxode2); library(nlmixr2est)
})
A <- asNamespace("admixr2")

mode  <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(mode)) mode <- "record"
store <- file.path("tests", "manual", "golden-perf.rds")

# ---- models ------------------------------------------------------------------

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

# Unpaired struct theta (ka has no eta) -> exercises the CRN-FD / unpaired paths
unpaired_fn <- function() {
  ini({
    tcl     <- log(5)  ; label("Log CL")
    tv      <- log(20) ; label("Log V")
    tka     <- log(1)  ; label("Log ka")
    add.err <- 0.1     ; label("Additive SD")
    eta.cl  ~ 0.09
    eta.v   ~ 0.04
  })
  model({
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    ka <- exp(tka)
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - (cl / v) * central
    cp <- central / v
    cp ~ add(add.err)
  })
}

# Two observed outputs -> exercises multi-output + joint branches
two_out_fn <- function() {
  ini({
    tcl      <- log(5)  ; label("Log CL")
    tv       <- log(20) ; label("Log V")
    tq       <- log(2)  ; label("Log Q")
    add.err  <- 0.1     ; label("Plasma additive SD")
    add.csf  <- 0.05    ; label("CSF additive SD")
    eta.cl   ~ 0.09
    eta.v    ~ 0.04
  })
  model({
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    q  <- exp(tq)
    d/dt(central) <- -(cl / v) * central - q * central + q * csfc
    d/dt(csfc)    <-  q * central - q * csfc
    cp   <- central / v
    cCSF <- csfc / v
    cp   ~ add(add.err)
    cCSF ~ add(add.csf)
  })
}

# ---- study builders ----------------------------------------------------------

times  <- c(0.5, 1, 2, 4, 8, 12)
E_true <- (100 / 20) * exp(-(5 / 20) * times)

study_1 <- list(s1 = list(E = E_true, V = diag((0.3 * E_true)^2), n = 200L,
                          times = times, ev = rxode2::et(amt = 100)))

# two studies, DIFFERENT dosing (blocks the same-ev grouping of item 6)
study_2 <- list(
  s1 = list(E = E_true, V = diag((0.3 * E_true)^2), n = 200L,
            times = times, ev = rxode2::et(amt = 100)),
  s2 = list(E = 2 * E_true, V = diag((0.3 * 2 * E_true)^2), n = 150L,
            times = times, ev = rxode2::et(amt = 200))
)

# two studies, SAME dosing, different times (item 6 SHOULD merge these)
t_a <- c(0.5, 1, 2, 4)
t_b <- c(1, 3, 8, 12)
study_sameev <- list(
  sa = list(E = (100/20) * exp(-(5/20) * t_a), V = diag((0.3 * (100/20) * exp(-(5/20)*t_a))^2),
            n = 200L, times = t_a, ev = rxode2::et(amt = 100)),
  sb = list(E = (100/20) * exp(-(5/20) * t_b), V = diag((0.3 * (100/20) * exp(-(5/20)*t_b))^2),
            n = 120L, times = t_b, ev = rxode2::et(amt = 100))
)

# a "var" method study (diagonal V) -> exercises the var branches
study_var <- list(sv = list(E = E_true, V = (0.3 * E_true)^2, n = 200L,
                            times = times, ev = rxode2::et(amt = 100)))

mk_multi <- function(joint) {
  tp <- c(0.5, 1, 2, 4); tc <- c(1, 4, 8)
  Ep <- (100/20) * exp(-(5/20) * tp)
  Ec <- 0.3 * (100/20) * exp(-(5/20) * tc)
  obs <- list(
    plasma = list(output = "cp",   ev = rxode2::et(amt = 100), times = tp,
                  E = Ep, V = diag((0.3 * Ep)^2), n = 40L),
    csf    = list(output = "cCSF", ev = rxode2::et(amt = 100), times = tc,
                  E = Ec, V = diag((0.3 * Ec)^2), n = 40L)
  )
  s <- list(observations = obs)
  if (joint) {
    # same subjects -> joint unit with cross-covariance
    s$n     <- 40L
    s$ev    <- rxode2::et(amt = 100)
    s$cross <- list(list(a = "plasma", b = "csf",
                         V = matrix(0.002, nrow = length(tp), ncol = length(tc))))
  }
  list(m1 = s)
}

# ---- internal setup ----------------------------------------------------------

setup <- function(fn, studies_raw, out_var = "cp", tag = FALSE) {
  ui        <- rxode2::rxode2(fn)
  sensModel <- A$.admLoadSensModel(ui)
  rxMod     <- A$.admLoadModel(ui)
  pinfo     <- A$.admParseIniDf(ui$iniDf, ui)
  ov        <- A$.admBuildOptVec(pinfo)
  st <- lapply(names(studies_raw), function(nm)
    A$.admNormaliseStudy(studies_raw[[nm]], nm, out_var))
  names(st) <- names(studies_raw)
  studies <- A$.admBuildEvFull(A$.admFlattenStudies(st), tag_cmt = tag)
  list(ui = ui, sensModel = sensModel, rxMod = rxMod, pinfo = pinfo, ov = ov,
       studies = studies, out_var = out_var)
}

# a couple of parameter points, not just the initial one
pts <- function(ov) {
  p0 <- ov$p0
  set.seed(7)
  list(p0 = p0,
       p1 = p0 + rnorm(length(p0), sd = 0.15),
       p2 = p0 - rnorm(length(p0), sd = 0.10))
}

res <- list()
add <- function(key, val) { res[[key]] <<- val; invisible(NULL) }

# ---- capture: estimator internals -------------------------------------------

capture_internals <- function(label, fn, studies_raw, out_var = "cp", tag = FALSE,
                              do_adirmc = TRUE) {
  s  <- setup(fn, studies_raw, out_var, tag)
  P  <- pts(s$ov)
  ns <- length(s$studies)
  pl1 <- A$.admMakeParamsList(1L, s$pinfo, ns)
  n_sim <- 400L
  z    <- A$.admMakeZ(n_sim, s$pinfo, ns, "sobol")
  plN  <- A$.admMakeParamsList(n_sim, s$pinfo, ns)
  grid <- A$.adghNodeGrid(3L, s$pinfo$n_eta)
  plQ  <- A$.admMakeParamsList(nrow(grid$X), s$pinfo, ns)

  for (pn in names(P)) {
    p <- P[[pn]]
    k <- function(w) paste(label, pn, w, sep = ".")

    add(k("adfoNLL"),  A$.adfoNLL(p, s$pinfo, s$studies, s$sensModel, s$rxMod,
                                  s$out_var, pl1, 1L))
    add(k("adfoGrad"), A$.adfoGrad(p, s$pinfo, s$studies, s$sensModel, s$rxMod,
                                   s$out_var, pl1, 1L))
    add(k("adfoFDGradF"), A$.adfoFDGrad(p, s$pinfo, s$studies, s$sensModel, s$rxMod,
                                        s$out_var, pl1, 1L, use_central = FALSE))
    add(k("adfoFDGradC"), A$.adfoFDGrad(p, s$pinfo, s$studies, s$sensModel, s$rxMod,
                                        s$out_var, pl1, 1L, use_central = TRUE))

    add(k("adghNLL"),  A$.adghNLL(p, s$pinfo, s$studies, s$rxMod, s$out_var, grid, 1L))
    add(k("adghGrad"), A$.adghGrad(p, s$pinfo, s$studies, s$sensModel, s$rxMod,
                                   s$out_var, grid, 1L))
    add(k("adghFDGrad"), A$.adghFDGrad(p, s$pinfo, s$studies, s$rxMod, s$out_var,
                                       grid, 1L))

    add(k("admNLL"),  A$.admNLL(p, s$pinfo, s$studies, z, s$rxMod, s$out_var, plN, 1L))
    add(k("admGradS"), A$.admGrad(p, s$pinfo, s$studies, z, s$rxMod, s$out_var,
                                  plN, 1L, 1e-4, s$sensModel))
    add(k("admGradF"), A$.admGrad(p, s$pinfo, s$studies, z, s$rxMod, s$out_var,
                                  plN, 1L, 1e-4, NULL))
    add(k("admNLLGradFD"), A$.admNLLGradFD(p, s$pinfo, s$studies, z, s$rxMod,
                                           s$out_var, plN, 1L, 1e-4))
    # batch paths (used by the Hessian)
    add(k("admNLLBatch"), A$.admNLLBatch(list(P$p0, P$p1, P$p2), s$pinfo, s$studies,
                                         z, s$rxMod, s$out_var, plN, 1L))
  }
  invisible(NULL)
}

# ---- capture: end-to-end fits ------------------------------------------------

capture_fit <- function(label, fn, studies_raw, est, ctl) {
  f <- tryCatch(suppressMessages(nlmixr2(fn, admData(A$.admOutputVars(rxode2::rxode2(fn))),
                                         est = est, control = ctl)),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(f)) { add(paste(label, est, "fit", sep = "."), f); return(invisible(NULL)) }
  add(paste(label, est, "obj", sep = "."), f$objective)
  add(paste(label, est, "est", sep = "."), f$parFixedDf$Estimate)
  add(paste(label, est, "se",  sep = "."), f$parFixedDf$SE)
  invisible(NULL)
}

cat("capturing internals: single-output, 1 study\n")
capture_internals("s1", one_cmt_fn, study_1)
cat("capturing internals: single-output, 2 studies (different ev)\n")
capture_internals("s2", one_cmt_fn, study_2)
cat("capturing internals: single-output, 2 studies (SAME ev, different times)\n")
capture_internals("sameev", one_cmt_fn, study_sameev)
cat("capturing internals: var-method study\n")
capture_internals("var", one_cmt_fn, study_var)
cat("capturing internals: unpaired struct theta (ka)\n")
capture_internals("unp", unpaired_fn, study_1)
cat("capturing internals: multi-output (independent blocks)\n")
capture_internals("multi", two_out_fn, mk_multi(FALSE), out_var = "cp", tag = TRUE)
cat("capturing internals: joint (same-subject, cross-cov)\n")
capture_internals("joint", two_out_fn, mk_multi(TRUE), out_var = "cp", tag = TRUE)

cat("capturing end-to-end fits\n")
for (est in c("adfo", "adgh", "admc", "adirmc")) {
  ctl <- switch(est,
    adfo   = adfoControl(studies = study_2, maxeval = 30L, seed = 1L,
                         grad = "analytical", cores = 2L, covMethod = "r", print = 0L),
    adgh   = adghControl(studies = study_2, maxeval = 30L, seed = 1L,
                         cores = 2L, covMethod = "r", print = 0L),
    admc   = admControl(studies = study_2, n_sim = 500L, maxeval = 20L, seed = 1L,
                        grad = "sens", cores = 2L, covMethod = "r", print = 0L),
    adirmc = adirmcControl(studies = study_2, n_sim = 500L, seed = 1L,
                           grad = "analytical", phases = c(1, 0.5), outer_iter = 3L,
                           cores = 2L, covMethod = "none", print = 0L))
  capture_fit("fit2", one_cmt_fn, study_2, est, ctl)
}
# multi-output fits (adirmc rejects these by design)
for (est in c("adfo", "adgh", "admc")) {
  ctl <- switch(est,
    adfo = adfoControl(studies = mk_multi(FALSE), maxeval = 20L, seed = 1L,
                       grad = "analytical", cores = 2L, covMethod = "none", print = 0L),
    adgh = adghControl(studies = mk_multi(FALSE), maxeval = 20L, seed = 1L,
                       cores = 2L, covMethod = "none", print = 0L),
    admc = admControl(studies = mk_multi(FALSE), n_sim = 400L, maxeval = 15L, seed = 1L,
                      grad = "sens", cores = 2L, covMethod = "none", print = 0L))
  capture_fit("fitmulti", two_out_fn, mk_multi(FALSE), est, ctl)
}
for (est in c("adfo", "adgh", "admc")) {
  ctl <- switch(est,
    adfo = adfoControl(studies = mk_multi(TRUE), maxeval = 20L, seed = 1L,
                       grad = "analytical", cores = 2L, covMethod = "none", print = 0L),
    adgh = adghControl(studies = mk_multi(TRUE), maxeval = 20L, seed = 1L,
                       cores = 2L, covMethod = "none", print = 0L),
    admc = admControl(studies = mk_multi(TRUE), n_sim = 400L, maxeval = 15L, seed = 1L,
                      grad = "sens", cores = 2L, covMethod = "none", print = 0L))
  capture_fit("fitjoint", two_out_fn, mk_multi(TRUE), est, ctl)
}

# ---- record / check ----------------------------------------------------------

if (identical(mode, "record")) {
  saveRDS(res, store)
  cat(sprintf("\nrecorded %d values -> %s\n", length(res), store))
} else {
  old <- readRDS(store)
  tol <- as.numeric(Sys.getenv("GOLDEN_TOL", "0"))
  cat(sprintf("\nchecking %d values against %s (tol = %g)\n", length(res), store, tol))
  keys <- union(names(old), names(res))
  bad <- character(0)
  for (k in keys) {
    if (!k %in% names(res)) { bad <- c(bad, sprintf("%-40s MISSING in new", k)); next }
    if (!k %in% names(old)) { bad <- c(bad, sprintf("%-40s NEW (not in baseline)", k)); next }
    a <- old[[k]]; b <- res[[k]]
    if (is.character(a) || is.character(b)) {
      if (!identical(a, b)) bad <- c(bad, sprintf("%-40s %s -> %s", k, a, b))
      next
    }
    if (length(a) != length(b)) { bad <- c(bad, sprintf("%-40s length %d -> %d", k, length(a), length(b))); next }
    fa <- is.finite(a); fb <- is.finite(b)
    if (!identical(fa, fb)) { bad <- c(bad, sprintf("%-40s finiteness changed", k)); next }
    if (!any(fa)) next
    d <- max(abs(a[fa] - b[fb]))
    r <- d / max(1e-12, max(abs(a[fa])))
    if (tol == 0) {
      if (d != 0) bad <- c(bad, sprintf("%-40s max|abs diff| = %.3e (rel %.2e)", k, d, r))
    } else if (r > tol) {
      bad <- c(bad, sprintf("%-40s rel diff %.3e > tol %.1e", k, r, tol))
    }
  }
  if (length(bad)) {
    cat(sprintf("\n*** %d MISMATCH(ES) ***\n", length(bad)))
    cat(paste(bad, collapse = "\n"), "\n")
    quit(status = 1)
  }
  cat("\nALL MATCH\n")
}
