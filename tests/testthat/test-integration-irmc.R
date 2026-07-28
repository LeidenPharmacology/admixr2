skip_if_not_installed("rxode2")
skip_on_cran()

# Reuses .int_grad_setup() — proposals already cached, no extra rxSolve.
# Covers paths not tested elsewhere: .adirmcInnerNLL() (both branches) and
# .adirmcNLL() monotonicity.

.int_irmc_setup <- function() {
  if (!is.null(.int_irmc_cache)) return(.int_irmc_cache)

  env <- .int_grad_setup()
  if (!env$proposals_ok) return(NULL)

  p0    <- env$vec$p0
  pars0 <- admixr2:::.admUnpack(p0, env$pinfo)
  prop  <- env$irmc_proposals[[1]]

  study_var <- env$studies[[1]]   # method = "var" (auto-detected, diagonal V)
  study_cov <- study_var
  study_cov$method <- "cov"

  # .adirmcInnerNLL now takes pinfo: the IS-weighted mean is computed inside the
  # C++ kernel, so the residual arrays must be built from the parameter info here
  # rather than pre-baked into the proposal.
  nll_var <- tryCatch(admixr2:::.adirmcInnerNLL(pars0, prop, study_var, env$pinfo),
                      error = function(e) NA_real_)
  nll_cov <- tryCatch(admixr2:::.adirmcInnerNLL(pars0, prop, study_cov, env$pinfo),
                      error = function(e) NA_real_)

  p_bad   <- p0; p_bad["tcl"] <- p_bad["tcl"] + 1.0
  nll_p0  <- admixr2:::.adirmcNLL(p0,    env$pinfo, env$studies, env$irmc_proposals)
  nll_bad <- admixr2:::.adirmcNLL(p_bad, env$pinfo, env$studies, env$irmc_proposals)

  .int_irmc_cache <<- list(
    nll_var = nll_var,
    nll_cov = nll_cov,
    nll_p0  = nll_p0,
    nll_bad = nll_bad
  )
  .int_irmc_cache
}

test_that("irmcInnerNLL: use_var branch finite and positive", {
  env <- .int_grad_setup()
  if (!env$proposals_ok) skip("proposal draw failed")
  irmc <- .int_irmc_setup()
  expect_true(is.finite(irmc$nll_var))
  expect_gt(irmc$nll_var, 0)
})

test_that("irmcInnerNLL: use_cov branch finite and positive", {
  env <- .int_grad_setup()
  if (!env$proposals_ok) skip("proposal draw failed")
  irmc <- .int_irmc_setup()
  expect_true(is.finite(irmc$nll_cov))
  expect_gt(irmc$nll_cov, 0)
})

test_that("irmcNLL: NLL at true params < NLL at substantially perturbed params", {
  env <- .int_grad_setup()
  if (!env$proposals_ok) skip("proposal draw failed")
  irmc <- .int_irmc_setup()
  expect_lt(irmc$nll_p0, irmc$nll_bad)
})

test_that("irmcInnerGrad with linearized kappa: ratio vs FD within 5%", {
  env <- .int_grad_lin_kappa_setup()
  if (!env$proposals_ok) skip("proposal draw failed")

  ratio    <- env$g_irmc_ana / env$g_irmc_fd
  ok       <- is.finite(env$g_irmc_fd) & abs(env$g_irmc_fd) > 1e-6
  if (sum(ok) == 0L) skip("All IRMC FD gradients near-zero or non-finite at p0")
  ratio_ok <- ratio[ok]
  bad      <- names(ratio_ok)[abs(ratio_ok - 1) > 0.05]
  expect_equal(length(bad), 0L,
    info = paste("Params with |ratio - 1| > 0.05:",
                 paste(sprintf("%s=%.4f", bad, ratio_ok[bad]), collapse = ", ")))
})

test_that("irmcInnerGrad with exact kappa: ratio vs FD within 5%", {
  env <- .int_irmc_exact_kappa_setup()
  if (!env$proposals_ok) skip("proposal draw failed")

  ratio    <- env$g_irmc_ana / env$g_irmc_fd
  ok       <- is.finite(env$g_irmc_fd) & abs(env$g_irmc_fd) > 1e-6
  if (sum(ok) == 0L) skip("All IRMC FD gradients near-zero or non-finite at p0")
  ratio_ok <- ratio[ok]
  bad      <- names(ratio_ok)[abs(ratio_ok - 1) > 0.05]
  expect_equal(length(bad), 0L,
    info = paste("Params with |ratio - 1| > 0.05:",
                 paste(sprintf("%s=%.4f", bad, ratio_ok[bad]), collapse = ", ")))
})

# The IRMC importance-sampling mean-shift for a mu-referenced paired theta is
# p_new - p_orig for EVERY transform: eta and theta share the argument of
# `param <- h(theta + eta)`, so shifting theta by Delta shifts eta's target mean by
# Delta regardless of h. The code only had this right for exp (log(exp(p)) = p);
# additive fell to log(p) and expit/probit to a natural-scale-log form -- all
# wrong. Measured against a direct adgh objective the expit shift was ~140 -2LL
# units off a few tenths from the proposal point. Pin the shift itself (exact, no
# importance-sampling variance): paired_type 0 and mean_new == the raw p-shift.
test_that("a paired theta's IS shift is the identity p-shift for every transform", {
  cases <- list(
    additive = list(
      fn = function() { ini({ tcl <- log(5); tv <- 20; a <- 0.5; eta.v ~ 4 })
        model({ cl <- exp(tcl); v <- tv + eta.v
                d/dt(central) <- -(cl/v)*central; cp <- central/v; cp ~ add(a) }) },
      theta = "tv", curEval = "", delta = 1.5),
    expit = list(
      fn = function() { ini({ tcl <- logit(5, 0, 20); tv <- log(20); a <- 0.4; eta.cl ~ 0.09 })
        model({ cl <- expit(tcl + eta.cl, 0, 20); v <- exp(tv)
                d/dt(central) <- -(cl/v)*central; cp <- central/v; cp ~ add(a) }) },
      theta = "tcl", curEval = "expit", delta = 0.3))

  for (nm in names(cases)) {
    cc    <- cases[[nm]]
    ui    <- suppressMessages(rxode2::rxode2(cc$fn))
    pinfo <- admixr2:::.admParseIniDf(ui$iniDf, ui)
    expect_identical(pinfo$struct_transforms[[cc$theta]]$curEval, cc$curEval,
                     label = sprintf("%s curEval", nm))
    tt <- c(1, 2, 4, 8); ev <- rxode2::et(amt = 100)
    E  <- 100/20*exp(-(5/20)*tt)
    s  <- admixr2:::.admNormaliseStudy(list(E = E, V = diag((0.2*E)^2), n = 100L,
                                            times = tt, ev = ev), "s")
    s$ev_full <- rxode2::et(s$ev, s$times)
    rxMod <- admixr2:::.admLoadModel(ui)
    z  <- admixr2:::.admMakeZ(1000L, pinfo, 1L, "sobol")
    pm <- admixr2:::.admMakeParamsList(1000L, pinfo, 1L)
    p0    <- admixr2:::.admBuildOptVec(pinfo)$p0
    pars0 <- admixr2:::.admUnpack(p0, pinfo)
    prop <- admixr2:::.adirmcProposal(
      rxMod, pars0$struct, pinfo$sigma_names, pars0$omega, omega_expansion = 2,
      s, z[[1]], "cp", pm[[1]], cores = 1L, pinfo$eta_col_names,
      has_kappa = pinfo$has_kappa, struct_transforms = pinfo$struct_transforms,
      struct_eta_idx = pinfo$struct_eta_idx, use_grad = TRUE)
    skip_if(is.null(prop), "proposal draw failed")

    k <- match(cc$theta, names(prop$log_origbeta))
    expect_false(is.na(k), label = sprintf("%s paired", nm))
    expect_identical(unname(prop$paired_types[k]), 0L,       # identity, not 1/2/3
                     label = sprintf("%s paired_type", nm))
    struct_new <- pars0$struct[names(prop$log_origbeta)]
    struct_new[cc$theta] <- struct_new[cc$theta] + cc$delta
    mn <- admixr2:::compute_mean_new_cpp(struct_new, prop$log_origbeta,
            prop$paired_types, prop$paired_lows, prop$paired_his)
    # exact p-shift, NOT log(back(p_new)) - log(back(p_orig))
    expect_equal(mn[k], cc$delta, tolerance = 1e-9,
                 label = sprintf("%s mean_new", nm))
  }
})
