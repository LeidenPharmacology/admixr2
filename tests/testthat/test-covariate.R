# Tier 1 -- covariate marginalisation helpers (no rxode2).
# Quadrature builders, covariance-of-covariate spec, per-study combine.

test_that("admBuildQuadrature: 1-d GL nodes + normalised weights", {
  skip_if_not_installed("statmod")
  q <- admBuildQuadrature(list(wt = list(mu = 70, sd = 10)), "gl", n_nodes = 9L)
  expect_equal(q$method, "gl")
  expect_equal(q$d, 1L)
  expect_length(q$wt_nodes, 9L)
  expect_length(q$weights, 9L)
  expect_equal(sum(q$weights), 1, tolerance = 1e-8)
  # nodes lie within the truncation interval and are centred on mu
  expect_true(all(q$wt_nodes >= 70 - 3.5 * 10 & q$wt_nodes <= 70 + 3.5 * 10))
  expect_equal(mean(range(q$wt_nodes)), 70, tolerance = 1e-8)
})

test_that("admBuildQuadrature: 1-d GH weights sum to 1, nodes exact-normal", {
  skip_if_not_installed("statmod")
  q <- admBuildQuadrature(list(wt = list(mu = 70, sd = 10)), "gh", n_nodes = 7L)
  expect_length(q$wt_nodes, 7L)
  expect_equal(sum(q$weights), 1, tolerance = 1e-8)
  expect_equal(mean(range(q$wt_nodes)), 70, tolerance = 1e-8)
})

test_that("admBuildQuadrature: 1-d Taylor has 3 evaluation nodes, no weights", {
  q <- admBuildQuadrature(list(wt = list(mu = 70, sd = 10)), "taylor", h = 2)
  expect_equal(q$wt_nodes, c(68, 70, 72))
  expect_null(q$weights)
  expect_equal(q$node_signs, matrix(c(-1L, 0L, 1L), 3L, 1L))
  expect_equal(q$Sigma_cov, matrix(100, 1, 1))
  expect_false(q$is_correlated)
  expect_equal(q$h, 2)  # default is 2.0
})

test_that("admBuildQuadrature: multi-d Taylor builds axis + cross nodes", {
  cd <- list(wt = list(mu = 70, sd = 10), age = list(mu = 40, sd = 8), rho = 0.5)
  q  <- admBuildQuadrature(cd, "taylor", h = c(2, 1.5))
  expect_equal(q$d, 2L)
  expect_true(q$is_correlated)
  sq <- rowSums(q$node_signs^2)
  expect_equal(sum(sq == 0L), 1L)   # one central node
  expect_equal(sum(sq == 1L), 4L)   # +/- per covariate
  expect_equal(sum(sq == 2L), 4L)   # one correlated pair -> 4 cross nodes
  # off-diagonal of Sigma_cov = rho * sd_wt * sd_age
  expect_equal(q$Sigma_cov[1, 2], 0.5 * 10 * 8, tolerance = 1e-8)
})

test_that(".admGetSigmaCov: rho scalar and explicit Sigma", {
  cd <- list(a = list(mu = 1, sd = 2), b = list(mu = 3, sd = 4), rho = 0.25)
  S  <- admixr2:::.admGetSigmaCov(cd)
  expect_equal(diag(S), c(4, 16))
  expect_equal(S[1, 2], 0.25 * 2 * 4)
  Sig <- matrix(c(1, 0.3, 0.3, 1), 2)
  cd2 <- list(a = list(mu = 1, sd = 2), b = list(mu = 3, sd = 4), Sigma = Sig)
  expect_equal(admixr2:::.admGetSigmaCov(cd2), Sig)
})

test_that(".admCovNodeVals pulls covariate values by name, 0 otherwise", {
  f <- admixr2:::.admCovNodeVals(list(wt = 70, age = 40), c("wt", "age", "other"))
  expect_equal(unname(f), c(70, 40, 0))
})

test_that(".admCovCols adds ONLY declared covariates, never a blanket fill", {
  m <- matrix(1, 3L, 1L, dimnames = list(NULL, "tcl"))
  # `vb` is a hard-coded model constant and `lam` an estimated TBS lambda:
  # both are model params, neither is a covariate, so neither may be added.
  got <- admixr2:::.admCovCols(m, c("tcl", "wt", "vb", "lam"), list(wt = 70))
  expect_equal(colnames(got), c("tcl", "wt"))
  expect_equal(unname(got[, "wt"]), rep(70, 3L))

  # a covariate the model does not read is not added either
  expect_equal(colnames(admixr2:::.admCovCols(m, c("tcl", "vb"), list(wt = 70))),
               "tcl")
  # no covariates at all -> untouched
  expect_identical(admixr2:::.admCovCols(m, c("tcl", "wt"), NULL), m)
})

test_that(".adm_combine_nll: weighted sum (no quadrature / GL)", {
  studies <- list(list(weight = 0.2), list(weight = 0.5), list(weight = 0.3))
  expect_equal(admixr2:::.adm_combine_nll(c(10, 20, 30), studies, NULL),
               0.2 * 10 + 0.5 * 20 + 0.3 * 30)
  # missing weights default to 1 (plain sum)
  expect_equal(admixr2:::.adm_combine_nll(c(1, 2, 3), list(list(), list(), list()), NULL), 6)
  q_gl <- list(method = "gl")
  expect_equal(admixr2:::.adm_combine_nll(c(10, 20, 30), studies, q_gl),
               0.2 * 10 + 0.5 * 20 + 0.3 * 30)
})

test_that(".adm_combine_nll: Taylor Hessian combine (1D)", {
  # nodes: minus(1), center(2), plus(3); sigma2=100, h=2
  q <- admBuildQuadrature(list(wt = list(mu = 70, sd = 10)), "taylor", h = 2)
  studies <- list(list(weight = 1), list(weight = 1), list(weight = 1))
  nll_m <- 10; nll_0 <- 8; nll_p <- 12
  # expected: nll_0 + 0.5 * (100/4) * (nll_p - 2*nll_0 + nll_m)
  expected <- nll_0 + 0.5 * (100 / 4) * (nll_p - 2 * nll_0 + nll_m)
  expect_equal(admixr2:::.adm_combine_nll(c(nll_m, nll_0, nll_p), studies, q), expected)
})

test_that(".adm_combine_grad: weighted sum (GL)", {
  studies <- list(list(weight = 0.5), list(weight = 0.5))
  g <- admixr2:::.adm_combine_grad(list(c(2, 4), c(6, 8)), c(1, 1), studies, NULL)
  expect_equal(g, c(0.5 * 2 + 0.5 * 6, 0.5 * 4 + 0.5 * 8))
})

test_that(".adm_combine_grad: Taylor Hessian combine (1D)", {
  q <- admBuildQuadrature(list(wt = list(mu = 70, sd = 10)), "taylor", h = 2)
  studies <- list(list(weight = 1), list(weight = 1), list(weight = 1))
  gm <- c(1, 2); g0 <- c(3, 4); gp <- c(5, 6)
  expected <- g0 + 0.5 * (100 / 4) * (gp - 2 * g0 + gm)
  result <- admixr2:::.adm_combine_grad(list(gm, g0, gp), c(10, 8, 12), studies, q)
  expect_equal(result, expected)
})

test_that("admBuildCovStudies pairs nodes with aggregate data (GL)", {
  skip_if_not_installed("statmod")
  qg  <- admBuildQuadrature(list(wt = list(mu = 70, sd = 10)), "gl", n_nodes = 5L)
  agg5 <- lapply(seq_len(5), function(k) list(E = c(1, 2), V = diag(c(0.1, 0.2))))
  stg <- admBuildCovStudies(agg5, qg, ev = NULL, times = c(1, 2), n = 100L)
  expect_length(stg, 5L)
  expect_equal(sum(vapply(stg, function(s) s$weight, numeric(1))), 1, tolerance = 1e-8)
  expect_equal(stg[[1]]$n, 100L)
})

test_that("admBuildCovStudies Taylor: 3-element agg creates 3 per-stratum studies", {
  q   <- admBuildQuadrature(list(wt = list(mu = 70, sd = 10)), "taylor", h = 2)
  agg <- lapply(seq_len(3), function(k) list(E = c(1, 2), V = diag(c(0.1, 0.2))))
  st  <- admBuildCovStudies(agg, q, ev = NULL, times = c(1, 2), n = 100L)
  expect_length(st, 3L)
  # covariate values at each node: 68, 70, 72
  expect_equal(unname(st[[1L]]$cov), 68)
  expect_equal(unname(st[[2L]]$cov), 70)
  expect_equal(unname(st[[3L]]$cov), 72)
  expect_equal(st[[1L]]$n, 100L)
})

test_that("admBuildCovStudies errors on node/data length mismatch", {
  skip_if_not_installed("statmod")
  qg  <- admBuildQuadrature(list(wt = list(mu = 70, sd = 10)), "gl", n_nodes = 5L)
  agg <- list(list(E = c(1, 2), V = diag(c(0.1, 0.2))))  # 1 != 5 nodes
  expect_error(admBuildCovStudies(agg, qg, ev = NULL, times = c(1, 2), n = 100L),
               "quadrature nodes")
  qt  <- admBuildQuadrature(list(wt = list(mu = 70, sd = 10)), "taylor", h = 2)
  agg1 <- list(list(E = c(1, 2), V = diag(c(0.1, 0.2))))  # 1 != 3 nodes
  expect_error(admBuildCovStudies(agg1, qt, ev = NULL, times = c(1, 2), n = 100L),
               "quadrature nodes")
})

test_that("admBuildCovStudies + admBuildQuadrature: quadrature attr attached manually", {
  q   <- admBuildQuadrature(list(wt = list(mu = 70, sd = 10)), "taylor", h = 2)
  agg <- lapply(seq_len(3), function(k) list(E = c(1, 2), V = diag(c(0.1, 0.2))))
  st  <- admBuildCovStudies(agg, q, times = c(1, 2), ev = NULL, n = 100L)
  attr(st, "quadrature") <- q
  expect_length(st, 3L)
  expect_false(is.null(attr(st, "quadrature")))
  expect_equal(attr(st, "quadrature")$method, "taylor")
  expect_equal(unname(st[[2L]]$cov), 70)
  expect_equal(st[[1L]]$n, 100L)
})

test_that("admBuildQuadrature: order=4 Taylor has 5 nodes (1D)", {
  q <- admBuildQuadrature(list(wt = list(mu = 70, sd = 10)), "taylor", h = 2, order = 4L)
  expect_equal(q$wt_nodes, c(66, 68, 70, 72, 74))
  expect_equal(q$order, 4L)
  expect_equal(q$node_signs, matrix(c(-2L, -1L, 0L, 1L, 2L), 5L, 1L))
  expect_null(q$weights)
})

test_that("admBuildQuadrature: order=4 multi-d Taylor has axis ±1 and ±2 nodes", {
  cd <- list(wt = list(mu = 70, sd = 10), age = list(mu = 40, sd = 8))
  q  <- admBuildQuadrature(cd, "taylor", h = c(2, 1.5), order = 4L)
  expect_equal(q$d, 2L)
  expect_equal(q$order, 4L)
  sq <- rowSums(q$node_signs^2)
  expect_equal(sum(sq == 0L), 1L)   # center
  expect_equal(sum(sq == 1L), 4L)   # ±1 per dim
  expect_equal(sum(sq == 4L), 4L)   # ±2 per dim (sign^2=4)
  expect_equal(nrow(q$node_signs), 9L)  # 1 + 4 + 4
})

test_that("admBuildQuadrature: invalid order errors", {
  expect_error(
    admBuildQuadrature(list(wt = list(mu = 70, sd = 10)), "taylor", order = 3L),
    "order must be 2 or 4")
})

test_that(".adm_combine_nll: Taylor order=4 Hessian combine (1D)", {
  # nodes order: -2, -1, 0, +1, +2; sigma2=100, h=2
  q <- admBuildQuadrature(list(wt = list(mu = 70, sd = 10)), "taylor", h = 2, order = 4L)
  studies <- replicate(5, list(weight = 1), simplify = FALSE)
  nll_m2 <- 9; nll_m1 <- 10; nll_0 <- 8; nll_p1 <- 12; nll_p2 <- 11
  sig2 <- 100; h <- 2
  # 2nd-order term: sig2 * (nll_p1 - 2*nll_0 + nll_m1) / h^2
  term2 <- sig2 * (nll_p1 - 2*nll_0 + nll_m1) / h^2
  # 4th-order term: sig2^2/8 * (nll_p2 - 4*nll_p1 + 6*nll_0 - 4*nll_m1 + nll_m2) / h^4
  term4 <- sig2^2 / 8 * (nll_p2 - 4*nll_p1 + 6*nll_0 - 4*nll_m1 + nll_m2) / h^4
  expected <- nll_0 + 0.5 * term2 + term4
  # node_signs column: -2,-1,0,+1,+2 -> nll_vec order matches wt_nodes order
  nll_vec <- c(nll_m2, nll_m1, nll_0, nll_p1, nll_p2)
  expect_equal(admixr2:::.adm_combine_nll(nll_vec, studies, q), expected)
})

test_that(".adm_combine_grad: Taylor order=4 combine (1D)", {
  q <- admBuildQuadrature(list(wt = list(mu = 70, sd = 10)), "taylor", h = 2, order = 4L)
  studies <- replicate(5, list(weight = 1), simplify = FALSE)
  gm2 <- c(1, 2); gm1 <- c(3, 4); g0 <- c(5, 6); gp1 <- c(7, 8); gp2 <- c(9, 10)
  sig2 <- 100; h <- 2
  term2_g <- sig2 * (gp1 - 2*g0 + gm1) / h^2
  term4_g <- sig2^2 / 8 * (gp2 - 4*gp1 + 6*g0 - 4*gm1 + gm2) / h^4
  expected <- g0 + 0.5 * term2_g + term4_g
  result <- admixr2:::.adm_combine_grad(
    list(gm2, gm1, g0, gp1, gp2), rep(8, 5), studies, q)
  expect_equal(result, expected)
})

test_that("datagen rejects multiple/mixed covariate studies (validation)", {
  f  <- function() NULL                 # validation never calls the model
  cs <- list(wt = list(mu = 70, sd = 10))
  expect_error(
    datagen(list(a = list(times = 1, ev = 1, covariate = cs),
                 b = list(times = 1, ev = 1, covariate = cs)), model = f),
    "at most one")
  expect_error(
    datagen(list(a = list(times = 1, ev = 1, covariate = cs),
                 b = list(times = 1, ev = 1)), model = f),
    "cannot be mixed")
})

