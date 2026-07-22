skip_if_not_installed("rxode2")
skip_on_cran()

# Setup in helper-integration.R. .admCalcCov() uses .admNLLBatch() (use_grad=FALSE)
# or .admGradBatch() (use_grad=TRUE) to compute a numerical Hessian, then
# returns 2*H^-1 restricted to struct + sigma parameters (omega Cholesky excluded).
# All rxSolve calls are batched: one call per study for all perturbed configs.
#
# .int_cov_setup() evaluates at p_cov where sigma_sd = 1 (not the true 0.1).
# At the true params, sigma contributes only 0.01 variance vs ~1.7 from IIV —
# making H[sigma,sigma] near-zero and non-PD regardless of step size or n_sim.
# Shifting sigma to sigma_sd = 1 makes it identifiable for structural tests.

# ---- NLL-FD Hessian (use_grad = FALSE, default) ------------------------------

test_that("admCalcCov NLL-FD: result is matrix of correct dimensions", {
  env_cov <- .int_cov_setup()
  result  <- env_cov$result_nll
  expect_false(is.null(result), info = "NLL-FD Hessian should be PD at sigma_sd=1")
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(env_cov$n_cov, env_cov$n_cov))
  expect_equal(rownames(result), env_cov$cov_names)
})

test_that("admCalcCov NLL-FD: result is symmetric", {
  env_cov <- .int_cov_setup()
  result  <- env_cov$result_nll
  expect_false(is.null(result), info = "NLL-FD Hessian should be PD at sigma_sd=1")
  expect_equal(result, t(result), tolerance = 1e-10)
})

test_that("admCalcCov NLL-FD: result is positive definite", {
  env_cov <- .int_cov_setup()
  result  <- env_cov$result_nll
  expect_false(is.null(result), info = "NLL-FD Hessian should be PD at sigma_sd=1")
  eigs <- eigen(result, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})

test_that("admCalcCov NLL-FD: omega is in the Hessian but NOT in the returned matrix", {
  env_cov <- .int_cov_setup()
  result  <- env_cov$result_nll
  expect_false(is.null(result), info = "NLL-FD Hessian should be PD at sigma_sd=1")
  # Omega used to be excluded. That did not merely forgo omega's own SEs -- it made
  # the STRUCTURAL ones too small, because a theta carrying an eta is correlated
  # with that eta's variance and profiling it out is not the same as fixing it.
  # Measured against the empirical sampling SD over simulated datasets, reported
  # SE / empirical SD for the eta-carrying theta went from 0.67 to 1.17 (prop) and
  # 0.67 to 1.06 (lnorm) once omega was included; a purely additive model barely
  # moved. Under-stated SEs give over-confident intervals, so that was the
  # dangerous direction of error.
  #
  # Omega is NOT returned, though. .admFullTheta() reports omega as the
  # variance/covariance entries while the optimizer holds the log-Cholesky, and
  # that map is not diagonal once omega is correlated -- so a correct omega SE
  # needs a full Jacobian, not a per-row factor. A mixed-scale matrix would be
  # worse than none, and nlmixr2est prints `Estimate +- 1.96*SE` from whatever
  # this returns.
  expect_false(any(env_cov$env$pinfo$omega_par_names %in% rownames(result)))
})

# ---- Grad-FD Hessian (use_grad = TRUE) ---------------------------------------

test_that("admCalcCov grad-FD: result dimensions match NLL-FD path", {
  env_cov  <- .int_cov_setup()
  res_nll  <- env_cov$result_nll
  res_grad <- env_cov$result_grad
  expect_false(is.null(res_nll),  info = "NLL-FD Hessian should be PD at sigma_sd=1")
  expect_false(is.null(res_grad), info = "grad-FD Hessian should be PD at sigma_sd=1")
  expect_equal(dim(res_grad), dim(res_nll))
})

# ---- Error paths: non-finite NLL and non-finite Hessian ----------------------

test_that("admCalcCov: non-finite Hessian entries warns and returns NULL", {
  env_cov <- .int_cov_setup()
  # cov_h_outer = 1e6 → h_gill[sigma] = max(|0|, 0.1) * 1e6 = 1e5.
  # sigma perturbed to exp(1e5/2) → Inf → Inf NLL → H[sigma,sigma] = Inf.
  result <- NULL
  expect_warning(
    { result <- admixr2:::.admCalcCov(
        env_cov$p_cov, env_cov$env$pinfo, env_cov$env$studies, env_cov$z_cov,
        env_cov$env$rxMod, env_cov$env$output_var, env_cov$params_cov, 1L,
        cov_n_sim = NULL, use_grad = FALSE, cov_h_outer = 1e6
      ) },
    "Hessian has non-finite entries"
  )
  expect_null(result)
})

test_that("admCalcCov: non-finite NLL at p_hat warns and returns NULL", {
  env <- .int_grad_setup()

  p_nonpd <- env$vec$p0
  p_nonpd[env$pinfo$omega_par_names[env$pinfo$chol_diag][1L]] <- -1e10

  result <- NULL
  expect_warning(
    { result <- admixr2:::.admCalcCov(
        p_nonpd, env$pinfo, env$studies, env$z_list,
        env$rxMod, env$output_var, env$params_list, 1L,
        cov_n_sim = NULL, use_grad = FALSE
      ) },
    "NLL not finite"
  )
  expect_null(result)
})

# ---- adfoCalcCov: NLL-FD Hessian (use_grad = FALSE, default) -----------------

test_that("adfoCalcCov NLL-FD: result is matrix of correct dimensions", {
  env_cov <- .int_adfo_cov_setup()
  result  <- env_cov$result_nll
  expect_false(is.null(result), info = "NLL-FD Hessian should be PD at sigma_sd=1")
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(env_cov$n_cov, env_cov$n_cov))
  expect_equal(rownames(result), env_cov$cov_names)
})

test_that("adfoCalcCov NLL-FD: result is symmetric", {
  env_cov <- .int_adfo_cov_setup()
  result  <- env_cov$result_nll
  expect_false(is.null(result), info = "NLL-FD Hessian should be PD at sigma_sd=1")
  expect_equal(result, t(result), tolerance = 1e-10)
})

test_that("adfoCalcCov NLL-FD: result is positive definite", {
  env_cov <- .int_adfo_cov_setup()
  result  <- env_cov$result_nll
  expect_false(is.null(result), info = "NLL-FD Hessian should be PD at sigma_sd=1")
  eigs <- eigen(result, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})

test_that("adfoCalcCov NLL-FD: omega is in the Hessian but NOT in the returned matrix", {
  env_cov <- .int_adfo_cov_setup()
  result  <- env_cov$result_nll
  expect_false(is.null(result), info = "NLL-FD Hessian should be PD at sigma_sd=1")
  # Omega is in the Hessian but not returned -- see the admCalcCov test above.
  expect_false(any(env_cov$env$pinfo$omega_par_names %in% rownames(result)))
})

# ---- adfoCalcCov: Grad-FD Hessian (use_grad = TRUE) --------------------------

test_that("adfoCalcCov grad-FD: result dimensions match NLL-FD path", {
  env_cov  <- .int_adfo_cov_setup()
  res_nll  <- env_cov$result_nll
  res_grad <- env_cov$result_grad
  expect_false(is.null(res_nll),  info = "NLL-FD Hessian should be PD at sigma_sd=1")
  expect_false(is.null(res_grad), info = "grad-FD Hessian should be PD at sigma_sd=1")
  expect_equal(dim(res_grad), dim(res_nll))
})

# ---- adfoCalcCov: Error paths -------------------------------------------------

test_that("adfoCalcCov: non-finite Hessian entries warns and returns NULL", {
  env_cov <- .int_adfo_cov_setup()
  # cov_h_outer = 1e6 → h_gill[sigma] = max(|0|, 0.1) * 1e6 = 1e5.
  # sigma perturbed to exp(1e5/2) → Inf → Inf NLL → H[sigma,sigma] = Inf.
  result <- NULL
  expect_warning(
    { result <- admixr2:::.adfoCalcCov(
        env_cov$p_cov, env_cov$env$pinfo, env_cov$env$studies,
        env_cov$env$sensModel, env_cov$env$rxMod, env_cov$env$output_var,
        env_cov$env$params_list, 1L,
        use_grad = FALSE, cov_h_outer = 1e6
      ) },
    "Hessian has non-finite entries"
  )
  expect_null(result)
})

test_that("adfoCalcCov: non-finite NLL at p_hat warns and returns NULL", {
  env <- .int_adfo_setup()

  p_nonfinite <- env$p0
  p_nonfinite["tcl"] <- 1e10   # exp(1e10) → Inf CL → rxSolve returns 0/NaN → Inf NLL

  nll_check <- suppressWarnings(admixr2:::.adfoNLL(
    p_nonfinite, env$pinfo, env$studies,
    env$sensModel, env$rxMod, env$output_var,
    env$params_list, 1L
  ))
  skip_if(is.finite(nll_check),
          "rxSolve returned finite NLL at tcl=1e10; error-path assumption invalid")

  result <- NULL
  expect_warning(
    { result <- admixr2:::.adfoCalcCov(
        p_nonfinite, env$pinfo, env$studies,
        env$sensModel, env$rxMod, env$output_var,
        env$params_list, 1L,
        use_grad = FALSE
      ) },
    "NLL not finite"
  )
  expect_null(result)
})

# ---- The returned covariance must be on the REPORTED scale -------------------

test_that("a sigma SE is on the same scale as the sigma estimate it is printed with", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  skip_if_not_installed("lotri")
  # nlmixr2est prints `Estimate +- 1.96*SE` straight from parFixedDf. .admFullTheta()
  # reports a "var"-role sigma back-transformed to an SD, while the Hessian is taken
  # w.r.t. log(sigma^2) -- so returning the raw optimizer-scale block made the printed
  # SE too large by a factor of 2/a. At a = 0.1 that is 20x, and it looked plausible
  # rather than obviously broken, which is why it needs an INDEPENDENT check rather
  # than one phrased in terms of the transform the code applies.
  #
  # The independent statement: with 2000 subjects, a residual SD is estimated
  # precisely, so its %RSE is a few percent. Under the scale bug it would be ~2/a in
  # relative terms, i.e. of order 2000%. The assertion below has a ~50x margin, so it
  # cannot pass by luck and does not restate the implementation.
  fn <- function() {
    ini({ tcl <- log(3); tv <- log(30); a <- 0.1; eta.cl ~ 0.09 })
    model({ cl <- exp(tcl + eta.cl); v <- exp(tv)
            d/dt(central) <- -(cl / v) * central; cp <- central / v
            cp ~ add(a) })
  }
  ui <- suppressMessages(rxode2::rxode2(fn))
  tt <- c(0.5, 1, 2, 4, 8); N <- 2000L
  set.seed(11)
  sim <- suppressWarnings(rxode2::rxSolve(
    ui$simulationModel, rxode2::et(rxode2::et(amt = 100), tt), nSub = N,
    omega = lotri::lotri(eta.cl ~ 0.09), sigma = lotri::lotri(rxerr.cp ~ 1),
    returnType = "data.frame", nDisplayProgress = .Machine$integer.max))
  Y <- matrix(sim$sim[sim$time %in% tt], N, length(tt), byrow = TRUE)
  d <- list(E = colMeans(Y), V = stats::cov.wt(Y, method = "ML")$cov,
            n = N, times = tt, ev = rxode2::et(amt = 100))
  fit <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    fn, admData(), est = "adgh",
    control = adghControl(studies = list(s = d), maxeval = 200L, covMethod = "r"))))
  skip_if(is.null(fit$cov), "covariance not computed")

  pf  <- fit$parFixedDf
  est <- unname(pf["a", "Estimate"])
  se  <- unname(pf["a", "SE"])
  expect_true(is.finite(se) && se > 0)
  expect_equal(est, 0.1, tolerance = 0.15)          # the estimate itself is sane
  expect_lt(se / est, 0.4)                          # scale bug would give ~20
  # ... and sigma is reported at all (it used to be uninitialised memory).
  expect_true("a" %in% rownames(fit$cov))
  expect_false(any(grepl("^logchol_", rownames(fit$cov))))
})
