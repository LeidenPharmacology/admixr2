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

  nll_var <- tryCatch(admixr2:::.adirmcInnerNLL(pars0, prop, study_var),
                      error = function(e) NA_real_)
  nll_cov <- tryCatch(admixr2:::.adirmcInnerNLL(pars0, prop, study_cov),
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
