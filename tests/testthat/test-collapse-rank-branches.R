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

test_that("a stale re-aim is an unsolvable point, not an error", {
  skip_if_not_installed("randtoolbox")
  cd <- list(W1 = list(meanlog = log(70), sdlog = 0.3),
             W2 = list(meanlog = log(70), sdlog = 0.3))
  ui <- list(lstExpr = list(
    quote(cl <- (W1 / 70)^0.6),
    quote(v <- (W1 / 70)^0.6 * (W2 / 70)^b)),
    allCovs = c("W1", "W2"))
  pin <- list(eta_col_names = character(), n_eta = 0L,
              struct_names = character(), struct_init = c(b = 0),
              cov_nodes = 7L)
  co <- admixr2:::.admCovCollapse(ui, pin, cd, 7L)
  expect_false(is.null(co))

  st   <- list(cov_dist = cd, .adm_cov_collapse = co)
  pr   <- list(struct = c(b = 1), L = matrix(0, 0, 0))
  grid <- list(X = matrix(0, 1L, 0L), W = 1)
  g    <- admixr2:::.adghGrid(pr, pin, grid, st)
  expect_true(isTRUE(g$failed))
  expect_null(g$eta)
  # rxMod = NULL: a marked grid must be reported BEFORE the solve. This errored
  # with "non-numeric matrix extent" on nrow(NULL), aborting the whole fit.
  expect_equal(admixr2:::.adghMoments(pr, pin, st, NULL, "cp", grid, 1L),
               list(failed = TRUE))
})

test_that("a degenerate point margin is not a dimension", {
  skip_if_not_installed("randtoolbox")
  pt <- list(.point = TRUE,
             quantile = function(u) rep(70, length.out = length(u)))
  cd <- list(W1 = list(meanlog = log(70), sdlog = 0.3), W2 = pt)
  ui <- list(lstExpr = list(quote(cl <- (W1 / 70)^0.6 * (W2 / 70)^0.4)),
             allCovs = c("W1", "W2"))
  pin <- list(eta_col_names = character(), n_eta = 0L,
              struct_names = character(), struct_init = numeric(),
              cov_nodes = 7L)
  # One real axis plus a point: the product grid already costs cov_nodes rows,
  # and the collapse used to price itself at 2 * cov_nodes against it.
  expect_null(admixr2:::.admCovCollapse(ui, pin, cd, 7L))
})

test_that("a reader of a discrete covariate alone does not refuse the study", {
  skip_if_not_installed("randtoolbox")
  cd <- list(WT = list(meanlog = log(70), sdlog = 0.3),
             SEX = list(values = c(0, 1), probs = c(0.5, 0.5)))
  ui <- list(lstExpr = list(
    quote(cl <- exp(eta.cl) * (WT / 70)^0.75),
    quote(f <- exp(0.2 * SEX))), allCovs = c("WT", "SEX"))
  pin <- list(eta_col_names = "eta.cl", n_eta = 1L,
              struct_names = character(), struct_init = numeric(),
              cov_nodes = 7L)
  jc <- admixr2:::.admJointCollapse(ui, pin, cd, 7L, NULL, NULL)
  expect_false(is.null(jc))
  # `f` is constant within a cell -- carried by the discrete cross, not by the
  # rotation -- so it contributes no direction and must not refuse the probe.
  expect_false(is.null(admixr2:::.admJointAdmit(jc, list(), matrix(0.3))))
})

test_that("coefficient probes keep the full single-index certificate", {
  skip_if_not_installed("randtoolbox")
  cd <- list(W1 = list(mu = 0, sd = 1), W2 = list(mu = 0, sd = 1))
  ui <- list(lstExpr = list(quote(cl <- exp(0.1 * W1 + b * W2^2))),
             allCovs = c("W1", "W2"))
  pin <- list(eta_col_names = character(), n_eta = 0L,
              struct_names = "b", struct_init = c(b = 0), cov_nodes = 7L)

  # At b = 0 this looks one-dimensional. Once b moves, W2^2 is not a linear
  # index; probing it at only one latent point admitted a wrong rank-1 law.
  expect_null(admixr2:::.admCovCollapse(ui, pin, cd, 7L))
})

test_that("collapse refuses correlated discrete cell probabilities", {
  cd <- list(W1 = list(mu = 0, sd = 1), W2 = list(mu = 0, sd = 1),
             S1 = list(values = c(0, 1), probs = c(0.5, 0.5)),
             S2 = list(values = c(0, 1), probs = c(0.5, 0.5)))
  R <- diag(4); R[3L, 4L] <- R[4L, 3L] <- 0.8
  expect_null(admixr2:::.admCovLatentBlock(
    cd, names(cd), c("W1", "W2"), c("S1", "S2"), R))
})

test_that("eta invariance is checked one random effect at a time", {
  skip_if_not_installed("randtoolbox")
  cd <- list(W1 = list(mu = 0, sd = 1), W2 = list(mu = 0, sd = 1))
  ui <- list(lstExpr = list(quote(
    cl <- exp(0.1 * W1 + (0.1 + eta.cl - eta.v) * W2))),
    allCovs = c("W1", "W2"))
  pin <- list(eta_col_names = c("eta.cl", "eta.v"), n_eta = 2L,
              struct_names = character(), struct_init = numeric(), cov_nodes = 7L)

  # Moving both etas by 0.5 cancels in eta.cl - eta.v; independent probes do not.
  expect_null(admixr2:::.admCovCollapse(ui, pin, cd, 7L))
})

test_that("fitted-parameter probes also certify eta and strata", {
  skip_if_not_installed("randtoolbox")
  cd <- list(W1 = list(mu = 0, sd = 1), W2 = list(mu = 0, sd = 1),
             SEX = list(values = c(0, 1), probs = c(.5, .5)))
  pin <- list(eta_col_names = "eta.cl", n_eta = 1L, struct_names = "b",
              struct_init = c(b = 0), cov_nodes = 7L)

  ui <- list(lstExpr = list(quote(
    cl <- exp(eta.cl + .1 * W1 + b * SEX * W2))), allCovs = names(cd))
  co <- admixr2:::.admCovCollapse(ui, pin, cd, 7L)
  expect_length(co$by_cell, 2L)
  expect_false(isTRUE(admixr2:::.admCovRefresh(co, list(b = 1))$stale))
  jc <- admixr2:::.admJointCollapse(ui, pin, cd, 7L, NULL, NULL)
  expect_null(admixr2:::.admJointAdmit(jc, list(b = 0), matrix(.3)))

  cd$SEX <- NULL
  ui$allCovs <- names(cd)
  ui$lstExpr <- list(quote(
    cl <- exp(eta.cl + .1 * W1 + b * eta.cl * W2)))
  expect_null(admixr2:::.admCovCollapse(ui, pin, cd, 7L))
})

test_that("deterministic verification keeps strong exact collapses", {
  skip_if_not_installed("randtoolbox")
  b <- 3 * c(.2, .35, -.15)
  cd <- stats::setNames(rep(list(list(mu = 0, sd = 1)), 3L), paste0("W", 1:3))
  ui <- list(lstExpr = list(bquote(
    cl <- exp(.(b[1]) * W1 + .(b[2]) * W2 + .(b[3]) * W3))),
    allCovs = names(cd))
  pin <- list(eta_col_names = character(), n_eta = 0L,
              struct_names = character(), struct_init = numeric(), cov_nodes = 7L)

  co <- admixr2:::.admCovCollapse(ui, pin, cd, 7L)
  expect_false(is.null(co))
  y <- exp(drop(as.matrix(co$X) %*% b))
  expect_equal(sum(co$W * y^2), exp(2 * sum(b^2)), tolerance = 1e-10)
})

test_that("ordinal cells keep their own continuous collapse", {
  skip_if_not_installed("randtoolbox")
  pr <- c(.2, .3, .5)
  cd <- list(W1 = list(mu = 0, sd = 1), W2 = list(mu = 0, sd = 1),
             GROUP = list(values = 0:2, probs = pr))
  ui <- list(lstExpr = list(quote(
    cl <- exp((.1 + .2 * GROUP) * W1 + (.4 - .1 * GROUP) * W2))),
    allCovs = names(cd))
  pin <- list(eta_col_names = character(), n_eta = 0L,
              struct_names = character(), struct_init = numeric(), cov_nodes = 7L)

  co <- admixr2:::.admCovCollapse(ui, pin, cd, 7L)
  expect_length(co$by_cell, 3L)
  expect_equal(nrow(co$X), 42L)
  y <- with(as.data.frame(co$X),
            exp((.1 + .2 * GROUP) * W1 + (.4 - .1 * GROUP) * W2))
  v <- vapply(0:2, function(k) (.1 + .2 * k)^2 + (.4 - .1 * k)^2, numeric(1))
  expect_equal(sum(co$W * y), sum(pr * exp(v / 2)), tolerance = 1e-10)
  expect_equal(sum(co$W * y^2), sum(pr * exp(2 * v)), tolerance = 1e-10)
})

test_that("eta invariance is certified in every discrete cell", {
  skip_if_not_installed("randtoolbox")
  cd <- list(W1 = list(mu = 0, sd = 1), W2 = list(mu = 0, sd = 1),
             SEX = list(values = c(0, 1), probs = c(.5, .5)))
  pin <- list(eta_col_names = "eta.cl", n_eta = 1L, struct_names = "b",
              struct_init = c(b = 1), cov_nodes = 7L)
  ui <- list(lstExpr = list(quote(
    cl <- exp(eta.cl + .1 * W1 + b * SEX * eta.cl * W2))),
    allCovs = names(cd))

  expect_null(admixr2:::.admCovCollapse(ui, pin, cd, 7L))
})

test_that("refresh retains the certificate points around a stationary link", {
  skip_if_not_installed("randtoolbox")
  cd <- list(W1 = list(mu = 0, sd = 1), W2 = list(mu = 0, sd = 1))
  pin <- list(eta_col_names = character(), n_eta = 0L,
              struct_names = "b", struct_init = c(b = 0), cov_nodes = 7L)
  ui <- list(lstExpr = list(quote(cl <- exp((.1 * W1 + .2 * W2 - b)^2))),
             allCovs = names(cd))
  co <- admixr2:::.admCovCollapse(ui, pin, cd, 7L)
  expect_false(is.null(co))
  expect_gt(nrow(co$z0), 1L)

  for (i in seq_len(nrow(co$z0))) {
    b <- sum(co$z0[i, ] * c(.1, .2))
    cr <- admixr2:::.admCovRefresh(co, list(b = b))
    expect_false(isTRUE(cr$stale))
    got <- sum(cr$W * exp((.1 * cr$X[, "W1"] + .2 * cr$X[, "W2"] - b)^2))
    expect_equal(got, exp(b^2 / .9) / sqrt(.9), tolerance = 1e-8)
  }
})

test_that("joint sizing includes the requested eta resolution", {
  skip_if_not_installed("randtoolbox")
  cd <- list(WT = list(mu = 0, sd = 1))
  ui <- list(lstExpr = list(quote(cl <- exp(eta.cl + b * WT))), allCovs = "WT")
  pin <- list(eta_col_names = "eta.cl", n_eta = 1L, struct_names = "b",
              struct_init = c(b = .1), cov_nodes = 7L)
  jc <- admixr2:::.admJointCollapse(ui, pin, cd, 7L, NULL, NULL,
                                    eta_nodes = 31L)
  jc <- admixr2:::.admJointAdmit(jc, list(b = .1), matrix(.7))
  expect_false(is.null(jc))
  expect_equal(jc$m, 38L)
  expect_equal(nrow(admixr2:::.admJointDesign(jc, list(b = .1), matrix(.7))$eta),
               38L)

  jc <- admixr2:::.admJointCollapse(ui, pin, cd, 7L, NULL, NULL,
                                    eta_nodes = 151L)
  jc <- admixr2:::.admJointAdmit(jc, list(b = .1), matrix(.7))
  expect_false(is.null(jc))
  expect_equal(jc$m, 158L)
})

test_that("an admitted reduction does not retain its product grid", {
  skip_if_not_installed("randtoolbox")
  cd <- stats::setNames(rep(list(list(mu = 0, sd = 1)), 5L), paste0("W", 1:5))
  ui <- list(lstExpr = list(quote(cl <- exp(.02 * (W1 + W2 + W3 + W4 + W5)))),
             allCovs = names(cd))
  pin <- list(eta_col_names = character(), n_eta = 0L,
              struct_names = character(), struct_init = numeric(),
              cov_nodes = 7L, cov_integration = "on")
  s <- admixr2:::.admCheckCovariates(ui, pin, list(s = list(cov_dist = cd)),
                                     "adgh")[[1L]]
  expect_false(is.null(s$.adm_cov_collapse))
  expect_null(s$.adm_cov_grid)
})

test_that("the base-point ORDER is frozen, not re-ranked when a point goes flat", {
  z0 <- rbind(c(0, 0), c(1, 0), c(0, 1))
  # one reader through a stationary link: log p = (z1 + z2 - a)^2, whose gradient
  # direction is constant and whose MAGNITUDE vanishes at z1 + z2 = a -- so the
  # coefficient decides which candidate point carries no signal.
  f_at <- function(a) function(Z) matrix(exp((Z[, 1] + Z[, 2] - a)^2), ncol = 1L)

  ord <- attr(admixr2:::.admCovGradB(f_at(0.3), z0), "at")
  expect_setequal(ord, seq_len(nrow(z0)))

  # replayed whole, at any parameter: the ranking is never redone
  for (a in c(0.1, 0.5, 0.9))
    expect_identical(attr(admixr2:::.admCovGradB(f_at(a), z0, i0 = ord), "at"),
                     ord)

  # with the leading point exactly stationary the loading is read at the NEXT
  # candidate in the frozen order -- not at an argmax re-ranked at these
  # parameters, which is what would rotate the SVD basis mid-fit.
  a  <- sum(z0[ord[1L], ])
  B  <- admixr2:::.admCovGradB(f_at(a), z0, i0 = ord)
  Bn <- admixr2:::.admCovGradB(f_at(a), z0, i0 = ord[-1L])
  expect_false(is.null(B))
  expect_equal(as.numeric(B), as.numeric(Bn), tolerance = 1e-8)
})
