# FO integration tests — Tier 2 (requires rxode2, slow)
# Skip on CRAN; each test guards with skip_on_cran() + skip_if_not_installed().
#
# Tests .adfoNLL() and .adfoGrad() against FD with the 1-cmt ODE model
# from helper-integration.R.

# ---- FO setup ----------------------------------------------------------------

.int_adfo_kappa_cache <- NULL

.int_adfo_setup <- function() {
  if (!is.null(.int_adfo_cache)) return(.int_adfo_cache)

  skip_on_cran()
  skip_if_not_installed("rxode2")

  # Reuse sens+rxMod from grad setup (ordering invariant already satisfied).
  env <- .int_grad_setup()
  if (is.null(env$rxMod)) skip("rxMod unavailable from grad setup")

  pinfo      <- env$pinfo
  studies    <- env$studies
  output_var <- env$output_var
  p0         <- env$vec$p0

  params_list <- admixr2:::.admMakeParamsList(1L, pinfo, length(studies))

  # Compute FO NLL at p0.
  nll_p0 <- admixr2:::.adfoNLL(p0, pinfo, studies, env$sensModel, env$rxMod,
                                 output_var, params_list, cores = 1L)

  # Perturb tcl away from truth.
  p_bad <- p0; p_bad["tcl"] <- p_bad["tcl"] + 0.5
  nll_p_bad <- admixr2:::.adfoNLL(p_bad, pinfo, studies, env$sensModel, env$rxMod,
                                    output_var, params_list, cores = 1L)

  # Analytical gradient at p0.
  h_fd <- 1e-4
  g_ana <- admixr2:::.adfoGrad(p0, pinfo, studies, env$sensModel, env$rxMod,
                                 output_var, params_list, cores = 1L, grad_h = h_fd)

  # FD gradient at p0 for comparison.
  g_fd <- vapply(seq_along(p0), function(k) {
    ph <- p0; ph[k] <- ph[k] + h_fd
    pl <- p0; pl[k] <- pl[k] - h_fd
    nh <- admixr2:::.adfoNLL(ph, pinfo, studies, env$sensModel, env$rxMod,
                               output_var, params_list, 1L)
    nl <- admixr2:::.adfoNLL(pl, pinfo, studies, env$sensModel, env$rxMod,
                               output_var, params_list, 1L)
    (nh - nl) / (2 * h_fd)
  }, double(1))
  names(g_fd) <- names(p0)

  .int_adfo_cache <<- list(
    pinfo = pinfo, studies = studies,
    rxMod = env$rxMod, sensModel = env$sensModel,
    output_var = output_var, params_list = params_list,
    p0 = p0, p_bad = p_bad,
    nll_p0 = nll_p0, nll_p_bad = nll_p_bad,
    g_ana = g_ana, g_fd = g_fd, h_fd = h_fd
  )
  .int_adfo_cache
}

.int_adfo_kappa_setup <- function() {
  if (!is.null(.int_adfo_kappa_cache)) return(.int_adfo_kappa_cache)

  skip_on_cran()
  skip_if_not_installed("rxode2")

  ui <- suppressMessages(tryCatch(
    rxode2::rxode2(one_cmt_kappa_fn),
    error = function(e) NULL
  ))
  if (is.null(ui)) skip("rxode2 model parse failed")

  pinfo      <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  output_var <- "cp"

  sensModel <- tryCatch(admixr2:::.admLoadSensModel(ui), error = function(e) NULL)
  rxMod     <- tryCatch(admixr2:::.admLoadModel(ui), error = function(e) NULL)
  if (is.null(rxMod)) skip("Model compilation failed")

  times  <- c(0.5, 1, 2, 4)
  E_true <- .one_cmt_mean(5, 20, 100, times)
  V_true <- diag((0.3 * E_true)^2)

  study <- list(E = E_true, V = V_true, n = 200L, times = times,
                ev = rxode2::et(amt = 100))
  study <- admixr2:::.admNormaliseStudy(study, "s")
  study$ev_full <- study$ev |> rxode2::et(study$times)
  studies <- list(s = study)

  params_list <- admixr2:::.admMakeParamsList(1L, pinfo, length(studies))
  p0          <- admixr2:::.admBuildOptVec(pinfo)$p0
  h_fd        <- 1e-4

  g_ana <- admixr2:::.adfoGrad(p0, pinfo, studies, sensModel, rxMod,
                               output_var, params_list, cores = 1L, grad_h = h_fd)

  g_fd <- vapply(seq_along(p0), function(k) {
    ph <- p0; ph[k] <- ph[k] + h_fd
    pl <- p0; pl[k] <- pl[k] - h_fd
    nh <- admixr2:::.adfoNLL(ph, pinfo, studies, sensModel, rxMod,
                               output_var, params_list, 1L)
    nl <- admixr2:::.adfoNLL(pl, pinfo, studies, sensModel, rxMod,
                               output_var, params_list, 1L)
    (nh - nl) / (2 * h_fd)
  }, double(1))
  names(g_fd) <- names(p0)

  .int_adfo_kappa_cache <<- list(
    pinfo = pinfo, studies = studies,
    rxMod = rxMod, sensModel = sensModel,
    output_var = output_var, params_list = params_list,
    p0 = p0, g_ana = g_ana, g_fd = g_fd, h_fd = h_fd
  )
  .int_adfo_kappa_cache
}

# ==============================================================================

test_that("FO NLL is finite and positive at truth", {
  skip_on_cran()
  env <- .int_adfo_setup()
  expect_true(is.finite(env$nll_p0))
  expect_gt(env$nll_p0, 0)
})

test_that("FO NLL is finite and positive at perturbed params", {
  skip_on_cran()
  env <- .int_adfo_setup()
  expect_true(is.finite(env$nll_p_bad))
  expect_gt(env$nll_p_bad, 0)
})

test_that("FO NLL increases when struct theta perturbed away from truth", {
  skip_on_cran()
  env <- .int_adfo_setup()
  # NLL at true theta should be lower than at perturbed theta
  expect_lt(env$nll_p0, env$nll_p_bad)
})

test_that("FO analytical gradient is finite at p0", {
  skip_on_cran()
  env <- .int_adfo_setup()
  expect_true(all(is.finite(env$g_ana)))
})

test_that("FO analytical gradient direction agrees with FD (ratio within 5% for all params)", {
  skip_on_cran()
  env <- .int_adfo_setup()

  ratio <- env$g_ana / env$g_fd
  # Exclude near-zero FD gradients to avoid 0/0 issues
  ok <- abs(env$g_fd) > 1e-6
  if (sum(ok) == 0) skip("All FD gradients near-zero at p0")

  ratio_ok <- ratio[ok]
  expect_true(all(abs(ratio_ok - 1) < 0.05),
    info = sprintf("max ratio deviation: %.4f (param: %s)",
                   max(abs(ratio_ok - 1)),
                   names(ratio_ok)[which.max(abs(ratio_ok - 1))]))
})

test_that("FO analytical diagonal omega gradients agree with FD (ratio within 5%)", {
  skip_on_cran()
  env <- .int_adfo_setup()

  n_s <- length(env$pinfo$struct_names)
  n_e <- length(env$pinfo$sigma_names)
  diag_idx <- which(env$pinfo$chol_diag) + n_s + n_e

  ratio <- env$g_ana[diag_idx] / env$g_fd[diag_idx]
  ok <- abs(env$g_fd[diag_idx]) > 1e-6
  if (sum(ok) == 0) skip("All diagonal omega FD gradients near-zero at p0")

  ratio_ok <- ratio[ok]
  expect_true(all(abs(ratio_ok - 1) < 0.05),
    info = sprintf("max diagonal omega ratio deviation: %.4f (param: %s)",
                   max(abs(ratio_ok - 1)),
                   names(ratio_ok)[which.max(abs(ratio_ok - 1))]))
})

test_that("FO analytical gradient matches FD for unpaired struct theta", {
  skip_on_cran()
  env <- .int_adfo_kappa_setup()

  expect_true(isFALSE(env$pinfo$struct_has_eta[["tsc"]]))

  ratio <- env$g_ana["tsc"] / env$g_fd["tsc"]
  if (!is.finite(env$g_fd["tsc"]) || abs(env$g_fd["tsc"]) <= 1e-6)
    skip("Unpaired struct theta FD gradient near-zero at p0")

  expect_lt(abs(ratio - 1), 0.05,
    info = sprintf("unpaired struct theta ratio: %.4f", ratio))
})


test_that("FO .adfoNLL() produces no R-level warnings at true params", {
  skip_on_cran()
  env <- .int_adfo_setup()

  warns <- character(0)
  withCallingHandlers(
    admixr2:::.adfoNLL(env$p0, env$pinfo, env$studies, env$sensModel,
                        env$rxMod, env$output_var, env$params_list, 1L),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warns, 0)
})
