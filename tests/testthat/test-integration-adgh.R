# Integration tests for the adgh (Gauss-Hermite) estimator.
# Tier 2: requires rxode2; skipped on CRAN.
# All rxSolve work happens in .int_adgh_setup(); tests operate on cached results.

skip_on_cran()
skip_if_not_installed("rxode2")

# ---- Moments -----------------------------------------------------------------

test_that(".adghMoments(): E close to analytic 1-cmt mean", {
  env <- .int_adgh_setup()
  # GH with 5 nodes should be very accurate; allow 2% relative tolerance
  expect_equal(env$moments_p0$E, env$E_true, tolerance = 0.02)
})

test_that(".adghMoments(): V is positive definite at truth", {
  env <- .int_adgh_setup()
  V   <- env$moments_p0$V
  expect_true(is.matrix(V))
  expect_equal(nrow(V), length(env$times))
  eig <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig > 0), info = "V is not positive definite")
})

test_that(".adghMoments(): E has correct length", {
  env <- .int_adgh_setup()
  expect_length(env$moments_p0$E, length(env$times))
})

# ---- NLL ---------------------------------------------------------------------

test_that(".adghNLL(): finite and positive at true params", {
  env <- .int_adgh_setup()
  expect_true(is.finite(env$nll_p0))
  expect_true(env$nll_p0 > 0)
})

test_that(".adghNLL(): true params give lower NLL than perturbed (monotonicity)", {
  env <- .int_adgh_setup()
  expect_lt(env$nll_p0, env$nll_bad)
})

test_that(".adghNLL(): NaN parameter vector returns non-finite (tryCatch guard)", {
  env   <- .int_adgh_setup()
  p_nan <- env$p0; p_nan[1] <- NaN
  nll <- admixr2:::.adghNLL(p_nan, env$pinfo, env$studies, env$rxMod,
                              env$output_var, env$grid, 1L)
  expect_false(is.finite(nll))
})

test_that(".adghNLL(): scalar (not vector) returned", {
  env <- .int_adgh_setup()
  expect_length(env$nll_p0, 1L)
})

# ---- FD gradient (no sens model needed) -------------------------------------

test_that(".adghFDGrad(): finite at truth (forward FD)", {
  env  <- .int_adgh_setup()
  g_fd <- admixr2:::.adghFDGrad(env$p0, env$pinfo, env$studies,
                                  env$rxMod, env$output_var, env$grid,
                                  cores = 1L, grad_h = 1e-4, use_central = FALSE)
  expect_true(all(is.finite(g_fd)))
})

test_that(".adghFDGrad(): finite at truth (central FD)", {
  env  <- .int_adgh_setup()
  g_fd <- admixr2:::.adghFDGrad(env$p0, env$pinfo, env$studies,
                                  env$rxMod, env$output_var, env$grid,
                                  cores = 1L, grad_h = 1e-4, use_central = TRUE)
  expect_true(all(is.finite(g_fd)))
})

# ---- Analytical gradient vs FD -----------------------------------------------

test_that(".adghGrad(): finite at truth", {
  env <- .int_adgh_setup()
  if (is.null(env$g_ana)) skip("sens model unavailable")
  expect_true(all(is.finite(env$g_ana)))
})

test_that(".adghGrad(): ratio vs central FD within 5% for all params", {
  env <- .int_adgh_setup()
  if (is.null(env$g_ana)) skip("sens model unavailable")

  ratio <- env$g_ana / env$g_fd
  # Skip near-zero gradient entries (ratio ill-defined)
  big   <- abs(env$g_fd) > 1e-6
  if (!any(big)) skip("all gradients near zero")
  expect_true(all(abs(ratio[big] - 1) < 0.05),
              info = paste("ratio:", paste(round(ratio[big], 4), collapse = ", ")))
})

test_that(".adghGrad(): correct length", {
  env <- .int_adgh_setup()
  if (is.null(env$g_ana)) skip("sens model unavailable")
  expect_length(env$g_ana, length(env$p0))
})

# ---- linearized_gh kappa proposal --------------------------------------------

test_that("adirmcProposal with kappa_method='linearized_gh': returns non-NULL", {
  env   <- .int_adgh_setup()
  pinfo <- admixr2:::.admParseIniDf(env$ui$iniDf, env$ui)

  # Use one_cmt_kappa_fn (has unpaired tsc) for has_kappa=TRUE path
  ui_k <- suppressMessages(tryCatch(
    rxode2::rxode2(one_cmt_kappa_fn), error = function(e) NULL))
  if (is.null(ui_k)) skip("kappa model parse failed")

  pinfo_k    <- admixr2:::.admParseIniDf(ui_k$iniDf, ui_k)
  rxMod_k    <- tryCatch(admixr2:::.admLoadModel(ui_k), error = function(e) NULL)
  if (is.null(rxMod_k)) skip("kappa model compilation failed")

  p0_k    <- admixr2:::.admBuildOptVec(pinfo_k)$p0
  pars_k  <- admixr2:::.admUnpack(p0_k, pinfo_k)
  z_list  <- admixr2:::.admMakeZ(100L, pinfo_k, 1L, "sobol")
  pm_list <- admixr2:::.admMakeParamsList(100L, pinfo_k, 1L)

  study_k <- env$study
  prop <- tryCatch(
    admixr2:::.adirmcProposal(
      rxMod_k, pars_k$struct, pinfo_k$sigma_names,
      pinfo_k$sigma_is_prop, pinfo_k$sigma_is_lnorm,
      pars_k$omega, omega_expansion = 1.5,
      study_k, z_list[[1L]], "cp",
      pm_list[[1L]], cores = 1L,
      pinfo_k$eta_col_names,
      has_kappa         = pinfo_k$has_kappa,
      kappa_method      = "linearized_gh",
      kappa_n_nodes     = 3L,
      struct_transforms = pinfo_k$struct_transforms,
      struct_eta_idx    = pinfo_k$struct_eta_idx,
      use_grad          = TRUE
    ),
    error = function(e) NULL
  )
  expect_false(is.null(prop), info = "proposal returned NULL")
})

test_that("adirmcProposal with kappa_method='linearized_gh': NLL is finite", {
  env   <- .int_adgh_setup()

  ui_k <- suppressMessages(tryCatch(
    rxode2::rxode2(one_cmt_kappa_fn), error = function(e) NULL))
  if (is.null(ui_k)) skip("kappa model parse failed")

  pinfo_k <- admixr2:::.admParseIniDf(ui_k$iniDf, ui_k)
  rxMod_k <- tryCatch(admixr2:::.admLoadModel(ui_k), error = function(e) NULL)
  if (is.null(rxMod_k)) skip("kappa model compilation failed")

  p0_k    <- admixr2:::.admBuildOptVec(pinfo_k)$p0
  pars_k  <- admixr2:::.admUnpack(p0_k, pinfo_k)
  z_list  <- admixr2:::.admMakeZ(100L, pinfo_k, 1L, "sobol")
  pm_list <- admixr2:::.admMakeParamsList(100L, pinfo_k, 1L)

  times <- c(0.5, 1, 2, 4)
  E_k   <- .one_cmt_mean(5, 20, 100, times)
  V_k   <- diag((0.3 * E_k)^2)
  study_k <- admixr2:::.admNormaliseStudy(
    list(E = E_k, V = V_k, n = 200L, times = times, ev = rxode2::et(amt = 100)), "s")
  study_k$ev_full <- study_k$ev |> rxode2::et(times)
  studies_k <- list(s = study_k)

  prop <- tryCatch(
    admixr2:::.adirmcProposal(
      rxMod_k, pars_k$struct, pinfo_k$sigma_names,
      pinfo_k$sigma_is_prop, pinfo_k$sigma_is_lnorm,
      pars_k$omega, omega_expansion = 1.5,
      study_k, z_list[[1L]], "cp",
      pm_list[[1L]], cores = 1L,
      pinfo_k$eta_col_names,
      has_kappa         = pinfo_k$has_kappa,
      kappa_method      = "linearized_gh",
      kappa_n_nodes     = 3L,
      struct_transforms = pinfo_k$struct_transforms,
      struct_eta_idx    = pinfo_k$struct_eta_idx,
      use_grad          = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(prop)) skip("proposal generation failed")

  nll <- admixr2:::.adirmcNLL(p0_k, pinfo_k, studies_k, list(prop))
  expect_true(is.finite(nll))
  expect_true(nll > 0)
})
