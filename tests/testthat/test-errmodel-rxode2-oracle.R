# The residual variance, checked against rxode2's OWN definition of it.
#
# rxode2 exports `.rxGetVarianceForErrorType(env, pred1)`, which returns the exact
# `rx_r_` LANGUAGE OBJECT it would emit for an endpoint -- the authoritative
# statement of what a residual error model means:
#
#   add(a)                     (a)^2
#   prop(b)                    (rx_pred_f_ * b)^2
#   pow(b, c)                  ((rx_pred_f_)^(c) * b)^2
#   add(a) + prop(b)           (a)^2 + (rx_pred_f_)^2 * (b)^2       [combined2]
#   ... + combined1()          ((a) + (rx_pred_f_) * (b))^2         [combined1]
#   ... + propT(b)             (a)^2 + (rx_pred_)^2 * (b)^2         [TRANSFORMED pred]
#
# Testing against this rather than against a formula retyped here is the point:
# both historical residual bugs in this package -- pow()'s exponent being treated
# as a variance, and combined1() being computed as combined2 -- are exactly the
# kind of thing that survives a hand-written expectation and dies against the
# generator. It also pins the recently-fixed behaviour that a prop() term on a
# TRANSFORMED endpoint contributes at all, and that propT() scales by the
# transformed prediction.

.orc_model <- function(ini, err) {
  eval(parse(text = sprintf(
    'function() { ini({ tcl <- log(5); tv <- log(20); %s
                        eta.cl ~ 0.09 })
       model({ cl <- exp(tcl + eta.cl); v <- exp(tv)
               d/dt(central) <- -(cl/v)*central; cp <- central/v
               %s }) }', ini, err)))
}

# rx_r_ evaluated at a given structural prediction f.
.orc_rx_r <- function(ui, f, pars) {
  ex <- rxode2::.rxGetVarianceForErrorType(ui, as.data.frame(ui$predDf)[1L, ])
  e  <- list2env(pars, parent = baseenv())
  e$rx_pred_f_ <- f                                    # untransformed prediction
  tr  <- as.character(ui$predDf$transform[1L])
  yj  <- if (tr %in% names(admixr2:::.ADM_TBS_YJ)) unname(admixr2:::.ADM_TBS_YJ[[tr]]) else 2L
  lam <- pars$lam %||% 1
  e$rx_pred_ <- if (identical(tr, "lnorm")) log(f)
                else if (yj == 2L) f
                else admixr2:::.admTBS(f, lam, yj,
                                       suppressWarnings(as.numeric(ui$predDf$trLow[1L] %||% 0)),
                                       suppressWarnings(as.numeric(ui$predDf$trHi[1L]  %||% 1)))
  eval(ex, envir = e)
}

test_that("UNTRANSFORMED residual variance equals rxode2's rx_r_ exactly", {
  skip_if_not_installed("rxode2")
  skip_if(is.null(tryCatch(rxode2::.rxGetVarianceForErrorType, error = function(e) NULL)),
          "this rxode2 does not expose .rxGetVarianceForErrorType")
  cases <- list(
    add       = list("a <- 0.5",             "cp ~ add(a)",                          list(a = 0.5)),
    prop      = list("b <- 0.2",             "cp ~ prop(b)",                         list(b = 0.2)),
    pow       = list("b <- 0.2; c1 <- 0.75", "cp ~ pow(b, c1)",                      list(b = 0.2, c1 = 0.75)),
    pow_hi    = list("b <- 0.2; c1 <- 1.30", "cp ~ pow(b, c1)",                      list(b = 0.2, c1 = 1.30)),
    combined2 = list("a <- 0.5; b <- 0.2",   "cp ~ add(a) + prop(b)",                list(a = 0.5, b = 0.2)),
    combined1 = list("a <- 0.5; b <- 0.2",   "cp ~ add(a) + prop(b) + combined1()",  list(a = 0.5, b = 0.2))
  )
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    ui <- suppressMessages(rxode2::rxode2(.orc_model(cs[[1L]], cs[[2L]])))
    p  <- suppressWarnings(admixr2:::.admParseIniDf(ui$iniDf, ui))
    for (f in c(0.3, 1, 4.7)) {
      arr <- admixr2:::.admResidRows(p, "cp", admixr2:::.admSigmaNat(p$sigma_init, p), 1L)
      # var_f = 0 isolates the CONDITIONAL variance, which is what rx_r_ states.
      dv  <- admixr2:::.admResidApply(f, 0, arr)$dv
      expect_equal(dv, .orc_rx_r(ui, f, cs[[3L]]), tolerance = 1e-12,
                   info = sprintf("%s at f = %g", nm, f))
    }
  }
})

test_that("Student-t scales rxode2's rx_r_ by nu/(nu-2)", {
  skip_if_not_installed("rxode2")
  skip_if(is.null(tryCatch(rxode2::.rxGetVarianceForErrorType, error = function(e) NULL)), "")
  ui <- suppressMessages(rxode2::rxode2(
    .orc_model("a <- 0.5; nu <- fix(5)", "cp ~ add(a) + t(nu)")))
  p   <- suppressWarnings(admixr2:::.admParseIniDf(ui$iniDf, ui))
  arr <- admixr2:::.admResidRows(p, "cp", admixr2:::.admSigmaNat(p$sigma_init, p), 1L)
  # nlmixr2 writes t() as a SCALE family (sim = rx_pred_ + sqrt(rx_r_)*rxt(nu)),
  # so rx_r_ is the SCALE^2 and the variance is that times nu/(nu-2).
  expect_equal(admixr2:::.admResidApply(2, 0, arr)$dv,
               .orc_rx_r(ui, 2, list(a = 0.5)) * (5 / 3), tolerance = 1e-12)
})

test_that("a prop()/propT() term on a TRANSFORMED endpoint uses rxode2's rx_r_", {
  skip_if_not_installed("rxode2")
  skip_if(is.null(tryCatch(rxode2::.rxGetVarianceForErrorType, error = function(e) NULL)), "")
  # For a transform-both-sides endpoint rx_r_ is the residual variance on the
  # TRANSFORMED scale, i.e. y = g(h(f) + sqrt(rx_r_) * eps). admixr2 integrates
  # that by quadrature, so the check is against a direct simulation of exactly
  # that generative statement -- with rxode2 supplying the variance.
  cases <- list(
    bc_add   = list("a <- 0.3; lam <- fix(0.5)",             "cp ~ add(a) + boxCox(lam)",             list(a = 0.3, lam = 0.5)),
    bc_prop  = list("a <- 0.3; b <- 0.2; lam <- fix(0.5)",   "cp ~ add(a) + prop(b) + boxCox(lam)",   list(a = 0.3, b = 0.2, lam = 0.5)),
    bc_propT = list("a <- 0.3; b <- 0.2; lam <- fix(0.5)",   "cp ~ add(a) + propT(b) + boxCox(lam)",  list(a = 0.3, b = 0.2, lam = 0.5)),
    logit    = list("a <- 0.3",                              "cp ~ logitNorm(a, 0, 40)",              list(a = 0.3))
  )
  set.seed(7); M <- 400000L; eps <- stats::rnorm(M)
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    ui <- suppressMessages(rxode2::rxode2(.orc_model(cs[[1L]], cs[[2L]])))
    p   <- suppressWarnings(admixr2:::.admParseIniDf(ui$iniDf, ui))
    arr <- admixr2:::.admResidRows(p, "cp", admixr2:::.admSigmaNat(p$sigma_init, p), 1L)
    f   <- 3.5
    lam <- cs[[3L]]$lam %||% 1
    yj  <- unname(admixr2:::.ADM_TBS_YJ[[as.character(ui$predDf$transform[1L])]])
    lo  <- suppressWarnings(as.numeric(ui$predDf$trLow[1L] %||% 0))
    hi  <- suppressWarnings(as.numeric(ui$predDf$trHi[1L]  %||% 1))
    if (!is.finite(lo)) lo <- 0
    if (!is.finite(hi)) hi <- 1

    sd_t <- sqrt(.orc_rx_r(ui, f, cs[[3L]]))            # rxode2's transformed-scale SD
    y    <- admixr2:::.admTBSi(admixr2:::.admTBS(f, lam, yj, lo, hi) + sd_t * eps, lam, yj, lo, hi)
    ap   <- admixr2:::.admResidApply(f, 0, arr)
    mc   <- stats::sd(y) / sqrt(M)
    expect_equal(ap$mu, mean(y), tolerance = 20 * mc / max(abs(mean(y)), 1e-8),
                 info = paste(nm, "mean"))
    expect_lt(abs(ap$dv - stats::var(y)), 40 * mc * max(stats::sd(y), 1))
  }
})

test_that("a prop() term on a transformed endpoint is not silently ignored", {
  skip_if_not_installed("rxode2")
  # Regression: the TBS branch once read only arr$a2, so b had no effect at all on
  # the objective and an exactly-zero gradient -- it was reported at its start value.
  ui1 <- suppressMessages(rxode2::rxode2(.orc_model(
    "a <- 0.3; lam <- fix(0.5)", "cp ~ add(a) + boxCox(lam)")))
  ui2 <- suppressMessages(rxode2::rxode2(.orc_model(
    "a <- 0.3; b <- 0.2; lam <- fix(0.5)", "cp ~ add(a) + prop(b) + boxCox(lam)")))
  dv <- function(ui) {
    p <- suppressWarnings(admixr2:::.admParseIniDf(ui$iniDf, ui))
    admixr2:::.admResidApply(3.5, 0.25,
                             admixr2:::.admResidRows(p, "cp", admixr2:::.admSigmaNat(p$sigma_init, p), 1L))$dv
  }
  expect_gt(abs(dv(ui2) - dv(ui1)), 1e-6)
})

test_that("ar() is refused with a prediction-dependent residual variance", {
  skip_if_not_installed("rxode2")
  # rxode2 emits
  #   rx.arRes = phi*lag0(rx.arRes,1) + sqrt(rx_r_*(1 - phi^2))*rxerr
  # whose marginal variance telescopes to rx_r_ ONLY when rx_r_ is constant. With
  # add() it does, and admixr2's sqrt(v_i v_j)*rho^|dt| form matches rxode2's own
  # simulation to within Monte-Carlo noise. With prop()/pow() it does not: the
  # process is non-stationary and admixr2's variance was measured 2.4-12x too
  # high (and its correlations 3-11x too high) against 1e6 simulated subjects --
  # silently, because nothing refused the combination.
  ok <- function(ini, err) {
    ui <- suppressMessages(rxode2::rxode2(.orc_model(ini, err)))
    suppressWarnings(admixr2:::.admParseIniDf(ui$iniDf, ui))
  }
  # add() + ar() is the case that DOES match the simulator -- it must keep working.
  expect_type(ok("a <- 0.5; rho <- 0.5",      "cp ~ add(a) + ar(rho)"), "list")
  expect_type(ok("a <- 0.5; rho <- fix(0.5)", "cp ~ add(a) + ar(rho)"), "list")
  for (m in list(c("b <- 0.2; rho <- 0.5",                 "cp ~ prop(b) + ar(rho)"),
                 c("a <- .5; b <- .2; rho <- 0.5",         "cp ~ add(a) + prop(b) + ar(rho)"),
                 c("b <- .2; c1 <- .8; rho <- 0.5",        "cp ~ pow(b, c1) + ar(rho)"),
                 # a FIXED prop coefficient still makes the variance f-dependent
                 c("a <- .5; b <- fix(.2); rho <- 0.5",    "cp ~ add(a) + prop(b) + ar(rho)")))
    expect_error(ok(m[1L], m[2L]), "proportional or power residual")
})
