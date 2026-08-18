# Bit-identity harness for the residual "moment tail" refactor -----------------
#
# The moment tail is .admResidApply -> .admApplyResidTail -> .admResidDeriv ->
# .admResidVChain, hand-assembled at six consumers (.admGrad, .adghGrad /
# .adghMomentsFromCp, .admNLLBatch, .admGradBatch, datagen() and .admAggData).
# Collapsing it into shared helpers must leave the OBJECTIVE bit-for-bit
# unchanged; CLAUDE.md records the previous errmodel refactor as "0 NLL
# mismatches, 62 gradient mismatches at <=1e-17 relative".
#
# This script evaluates, per error model and per study method (var and cov):
#   .adfoNLL/.adfoGrad, .adghNLL/.adghGrad, .admNLL/.admGrad,
#   .admNLLBatch/.admGradBatch, .adghMoments, the raw .admResidApply/
#   .admApplyResidTail assembly plot.R's .admAggData() repeats, the gradient
#   chain (.admResidDeriv/.admResidVChain/.admSigmaGrad/.admResidMuCoupling),
#   and datagen()'s mc/gh/fo branches.
# adirmc needs a proposal object and has its own probe:
#   validation/resid-tail-adirmc-probe.R
#
# Usage (four run modes plus compare):
#   Rscript validation/resid-moment-tail-bitidentity.R before
#   Rscript validation/resid-moment-tail-bitidentity.R before-dg
#   Rscript validation/resid-moment-tail-bitidentity.R after
#   Rscript validation/resid-moment-tail-bitidentity.R after-dg
#   Rscript validation/resid-moment-tail-bitidentity.R compare

suppressMessages(devtools::load_all(".", quiet = TRUE))
mode <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(mode)) mode <- "before"
# datagen() calls rxode2::rxUnloadAll() on exit, which pulls every compiled model
# out from under the rest of the process -- so it gets its own run/file.
dg_only <- grepl("-dg$", mode)
stem    <- sub("-dg$", "", mode)
outf <- file.path("validation", sprintf("resid-tail-%s.rds", mode))

if (identical(mode, "compare")) {
  ld <- function(w) {
    f1 <- file.path("validation", sprintf("resid-tail-%s.rds", w))
    f2 <- file.path("validation", sprintf("resid-tail-%s-dg.rds", w))
    c(if (file.exists(f1)) readRDS(f1) else list(),
      if (file.exists(f2)) readRDS(f2) else list())
  }
  b <- ld("before"); a <- ld("after")
  nm <- union(names(b), names(a))
  rows <- lapply(nm, function(k) {
    x <- b[[k]]; y <- a[[k]]
    if (is.null(x) || is.null(y))
      return(data.frame(key = k, n = NA_integer_, bit_eq = NA_integer_,
                        max_rel = NA_real_, note = "missing"))
    x <- as.numeric(x); y <- as.numeric(y)
    if (length(x) != length(y))
      return(data.frame(key = k, n = length(x), bit_eq = NA_integer_,
                        max_rel = NA_real_, note = "length"))
    # bitwise on the raw IEEE-754 payload
    bx <- vapply(x, function(v) paste(writeBin(v, raw()), collapse = ""), "")
    by <- vapply(y, function(v) paste(writeBin(v, raw()), collapse = ""), "")
    eq <- bx == by
    rel <- ifelse(eq, 0, abs(x - y) / pmax(abs(x), 1e-300))
    data.frame(key = k, n = length(x), bit_eq = sum(eq),
               max_rel = if (length(rel)) max(rel, na.rm = TRUE) else 0,
               note = "")
  })
  tab <- do.call(rbind, rows)
  print(tab, row.names = FALSE)
  cat("\n== TOTALS ==\n")
  cat(sprintf("values compared : %d\n", sum(tab$n, na.rm = TRUE)))
  cat(sprintf("bit-identical   : %d\n", sum(tab$bit_eq, na.rm = TRUE)))
  cat(sprintf("mismatches      : %d\n",
              sum(tab$n, na.rm = TRUE) - sum(tab$bit_eq, na.rm = TRUE)))
  cat(sprintf("max rel diff    : %.3e\n", max(tab$max_rel, na.rm = TRUE)))
  nll <- tab[grepl("nll", tab$key), ]
  cat(sprintf("NLL mismatches  : %d (of %d)\n",
              sum(nll$n) - sum(nll$bit_eq), sum(nll$n)))
  gr <- tab[grepl("grad", tab$key), ]
  cat(sprintf("grad mismatches : %d (of %d), max rel %.3e\n",
              sum(gr$n) - sum(gr$bit_eq), sum(gr$n), max(gr$max_rel, na.rm = TRUE)))
  quit(save = "no")
}

# -- model factory (mirrors test-integration-resid-moments.R) ------------------
.mkf <- function(ini, err) eval(parse(text = sprintf(
  'function() { ini({ tcl <- log(5); tv <- log(20); %s
                      eta.cl ~ 0.09; eta.v ~ 0.04 })
     model({ cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
             d/dt(central) <- -(cl/v)*central; cp <- central/v
             %s }) }', ini, err)))
.mk <- function(ini, err) suppressMessages(rxode2::rxode2(.mkf(ini, err)))

cases <- list(
  add    = c("a <- 0.5",                       "cp ~ add(a)"),
  prop   = c("b <- 0.15",                      "cp ~ prop(b)"),
  comb2  = c("a <- 0.5; b <- 0.15",            "cp ~ add(a) + prop(b)"),
  comb1  = c("a <- 0.5; b <- 0.15",            "cp ~ add(a) + prop(b) + combined1()"),
  lnorm  = c("a <- 0.3",                       "cp ~ lnorm(a)"),
  pow    = c("b <- 0.15; c1 <- 1.3",           "cp ~ pow(b, c1)"),
  propt  = c("b <- 0.15; nu <- fix(5)",        "cp ~ prop(b) + t(nu)"),
  bc     = c("a <- 0.3; lam <- 0.4",           "cp ~ add(a) + boxCox(lam)"),
  arm    = c("a <- 0.5; rho <- 0.4",           "cp ~ add(a) + ar(rho)"),
  pois   = c("",                               "y ~ pois(cp)")
)

times <- c(0.5, 1, 2, 4)
E0    <- 100 / 20 * exp(-(5 / 20) * times)
# a genuinely NON-diagonal V -> method "cov"; a diagonal one -> method "var"
Vcov  <- outer(0.3 * E0, 0.3 * E0) * (0.4 ^ abs(outer(times, times, "-")))
Vvar  <- diag((0.3 * E0) ^ 2)

res <- list()
put <- function(key, val) {
  v <- suppressWarnings(tryCatch(as.numeric(val), error = function(e) NA_real_))
  res[[key]] <<- v
}

for (nm in names(cases)) {
  cs <- cases[[nm]]
  if (dg_only) {
    # datagen() takes the model FUNCTION, and unloads every compiled model on exit
    for (mth in c("mc", "gh", "fo")) {
      put(sprintf("%s/datagen.%s", nm, mth), tryCatch({
        set.seed(42L)
        dg <- datagen(list(d = list(times = times, n = 100L,
                                    ev = rxode2::et(amt = 100))),
                      model = .mkf(cs[1], cs[2]),
                      control = datagenControl(method = mth, n_sim = 200L,
                                               cores = 1L))
        st1 <- dg[["d"]]
        c(as.numeric(st1$E), as.numeric(st1$V))
      }, error = function(e) NA_real_))
    }
    message("done (datagen): ", nm)
    next
  }
  ui <- tryCatch(.mk(cs[1], cs[2]), error = function(e) NULL)
  if (is.null(ui)) { message("skip (parse): ", nm); next }
  pinfo <- suppressWarnings(.admParseIniDf(ui$iniDf, ui))
  sens  <- suppressMessages(tryCatch(.admLoadSensModel(ui), error = function(e) NULL))
  rxMod <- .admLoadModel(ui)
  ovar  <- .admOutputVar(ui)

  for (meth in c("var", "cov")) {
    # ar() is refused on a var-method study (rho is unidentifiable from a diagonal V)
    if (identical(nm, "arm") && identical(meth, "var")) next
    Vs <- if (meth == "cov") Vcov else Vvar
    s <- tryCatch(.admNormaliseStudy(list(E = E0, V = Vs, n = 200L, times = times,
                                          ev = rxode2::et(amt = 100)), "s"),
                  error = function(e) NULL)
    if (is.null(s)) next
    s$ev_full <- rxode2::et(s$ev, s$times)
    st <- list(s = s)
    z  <- .admMakeZ(300L, pinfo, 1L, "sobol")
    pm <- .admMakeParamsList(300L, pinfo, 1L)
    pm1 <- .admMakeParamsList(1L, pinfo, 1L)
    grid <- .adghNodeGrid(5L, pinfo$n_eta)
    p0 <- .admBuildOptVec(pinfo)$p0 * 1.03 + 0.01
    k <- function(w) sprintf("%s/%s/%s", nm, meth, w)

    put(k("admc.nll"),  tryCatch(.admNLL(p0, pinfo, st, z, rxMod, ovar, pm, 1L),
                                 error = function(e) NA_real_))
    put(k("admc.grad"), tryCatch(.admGrad(p0, pinfo, st, z, rxMod, ovar, pm, 1L,
                                          1e-5, sens), error = function(e) NA_real_))
    put(k("adgh.nll"),  tryCatch(.adghNLL(p0, pinfo, st, rxMod, ovar, grid, 1L),
                                 error = function(e) NA_real_))
    put(k("adgh.grad"), tryCatch(.adghGrad(p0, pinfo, st, sens, rxMod, ovar, grid,
                                           1L, 1e-5), error = function(e) NA_real_))
    put(k("adfo.nll"),  tryCatch(.adfoNLL(p0, pinfo, st, sens, rxMod, ovar, pm1, 1L),
                                 error = function(e) NA_real_))
    put(k("adfo.grad"), tryCatch(.adfoGrad(p0, pinfo, st, sens, rxMod, ovar, pm1,
                                           1L, 1e-5), error = function(e) NA_real_))

    # batched paths (the post-fit Hessian evaluators)
    pl <- list(p0, p0 * 1.01, p0 * 0.99)
    put(k("batch.nll"),  tryCatch(.admNLLBatch(pl, pinfo, st, z, rxMod, ovar, pm, 1L),
                                  error = function(e) NA_real_))
    put(k("batch.grad"), tryCatch(unlist(.admGradBatch(pl, pinfo, st, z, rxMod, ovar,
                                                       pm, 1L, 1e-5, sens)),
                                  error = function(e) NA_real_))

    # adgh's shared moment assembly, called directly
    put(k("adghmom"), tryCatch({
      pars <- .admUnpack(p0, pinfo)
      m <- .adghMoments(pars, pinfo, s, rxMod, ovar, grid, 1L)
      c(as.numeric(m$E), as.numeric(m$V))
    }, error = function(e) NA_real_))

    # .admAggData's assembly, exercised through .admJointResidual's sibling path:
    # plot.R rebuilds arr + .admResidApply + .admApplyResidTail on a cp_mat.
    put(k("aggdata"), tryCatch({
      pars <- .admUnpack(p0, pinfo)
      eta  <- z[[1L]] %*% t(pars$L); colnames(eta) <- pinfo$eta_col_names
      cp   <- .admSimulate(rxMod, pars$struct, pinfo$sigma_names, eta, s, ovar,
                           pm[[1L]], 1L)
      arr  <- .admUnitResidRows(pinfo, ovar, pars$sigma_var, length(times),
                                phi = attr(cp, "phi"))
      mu   <- colMeans(cp); cpc <- sweep(cp, 2L, mu)
      V    <- crossprod(cpc) / nrow(cp)
      ap   <- .admResidApply(mu, diag(V), arr, times, V)
      c(as.numeric(ap$mu), as.numeric(.admApplyResidTail(V, ap)))
    }, error = function(e) NA_real_))

    # the gradient chain tail, evaluated directly
    put(k("chain"), tryCatch({
      pars <- .admUnpack(p0, pinfo)
      eta  <- z[[1L]] %*% t(pars$L); colnames(eta) <- pinfo$eta_col_names
      cp   <- .admSimulate(rxMod, pars$struct, pinfo$sigma_names, eta, s, ovar,
                           pm[[1L]], 1L)
      arr  <- .admUnitResidRows(pinfo, ovar, pars$sigma_var, length(times),
                                phi = attr(cp, "phi"))
      mu   <- colMeans(cp); cpc <- sweep(cp, 2L, mu)
      cov_f <- crossprod(cpc) / nrow(cp); var_f <- diag(cov_f)
      dres <- .admResidDeriv(mu, var_f, arr, pinfo)
      vch  <- .admResidVChain(mu, var_f, arr, pinfo, times, deriv = dres)
      dV   <- diag(length(times)) + 0.1
      dmu  <- rep(0.3, length(times))
      sg   <- .admSigmaGrad(mu, arr, pinfo, diag(dV), dmu, var_f, dV, times,
                            cov_f, deriv = dres)
      mc   <- .admResidMuCoupling(mu, arr, pinfo, diag(dV), dmu, var_f, dV,
                                  cov_f, times, deriv = dres)
      c(as.numeric(vch), as.numeric(attr(vch, "dmu_dv0")), as.numeric(sg),
        as.numeric(mc))
    }, error = function(e) NA_real_))
  }

  message("done: ", nm)
}

saveRDS(res, outf)
cat(sprintf("wrote %s  (%d keys, %d doubles)\n", outf, length(res),
            sum(vapply(res, length, integer(1)))))
