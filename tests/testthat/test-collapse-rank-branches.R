test_that("collapse rank is stable over fitted coefficients", {
  skip_if_not_installed("randtoolbox")
  cd <- list(W1 = list(meanlog = log(70), sdlog = 0.3),
             W2 = list(meanlog = log(70), sdlog = 0.3))
  ui <- list(lstExpr = list(
    quote(cl <- (W1 / 70)^0.6),
    quote(v <- (W1 / 70)^0.6 * (W2 / 70)^b)),
    allCovs = c("W1", "W2"))
  pin <- list(eta_col_names = character(), n_eta = 0L,
              struct_names = "b", struct_init = c(b = 0), cov_nodes = 7L)

  # At b = 0 the two loading columns are collinear; at b != 0 they are not.
  # There is no structural reduction from two covariates to one direction.
  expect_null(admixr2:::.admCovCollapse(ui, pin, cd, 7L))

  # Refresh is the backstop when admission lacks parameter metadata.
  pin$struct_names <- character()
  co <- admixr2:::.admCovCollapse(ui, pin, cd, 7L)
  expect_false(is.null(co))
  expect_true(isTRUE(admixr2:::.admCovRefresh(co, list(b = 1))$stale))

  # The joint space still reduces from three latent axes to two readers.
  pin$n_eta <- 1L
  pin$eta_col_names <- "eta.cl"
  pin$struct_names <- "b"
  ui$lstExpr <- list(
    quote(cl <- exp(eta.cl) * (W1 / 70)^0.6),
    quote(v <- exp(eta.cl) * (W1 / 70)^0.6 * (W2 / 70)^b))
  jc <- admixr2:::.admJointCollapse(ui, pin, cd, 7L, NULL, NULL)
  jc <- admixr2:::.admJointAdmit(jc, list(b = 0), matrix(0.3))
  expect_false(is.null(jc))
  expect_equal(jc$r, 2L)
  expect_false(is.null(admixr2:::.admJointDesign(jc, list(b = 1), matrix(0.3))))
})

test_that("joint collapse refuses random effects inside branches", {
  skip_if_not_installed("randtoolbox")
  cd <- list(WT = list(meanlog = log(70), sdlog = 0.3))
  pin <- list(eta_col_names = c("eta.cl", "eta.v"), n_eta = 2L,
              struct_names = character(), cov_nodes = 7L)
  ui <- list(lstExpr = list(
    quote(cl <- exp(eta.cl) * (WT / 70)^0.6),
    quote(if (eta.v > 0) v <- 20 else v <- 10)), allCovs = "WT")

  expect_null(admixr2:::.admJointCollapse(ui, pin, cd, 7L, NULL, NULL))
})
