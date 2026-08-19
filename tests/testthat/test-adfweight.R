# The ADF weight: the sampling law of (ybar, vech V) implied by the model.

test_that("the fast weight equals the reference expansion", {
  # .admAdfWeight is the readable Wick expansion; .admAdfWeightFast hoists the
  # node contraction out of the q x q loop and is what runs. They must not drift:
  # every term is a weighted sum over nodes of products of C and Dv columns, and
  # the fast form only reorders when that sum is taken.
  set.seed(4)
  for (m in c(3L, 5L, 8L)) {
    Q  <- 15L
    C  <- matrix(stats::rnorm(Q * m), Q, m)
    C  <- sweep(C, 2L, colMeans(C))
    w  <- stats::runif(Q); w <- w / sum(w)
    Dv <- matrix(stats::runif(Q * m, 0.01, 0.3), Q, m)
    a  <- admixr2:::.admAdfWeight(C, w, Dv, 200)
    b  <- admixr2:::.admAdfWeightFast(C, w, Dv, 200)
    expect_equal(b, a, tolerance = 1e-12)
    expect_true(isSymmetric(round(b, 12)))
  }
})

test_that("an exactly normal ensemble collapses onto the normal-theory weight", {
  # The correction IS the difference between the two, so where the marginal is
  # genuinely multivariate normal they must agree -- otherwise the estimator is
  # not a generalisation of the one it replaces. A single node carries no mixture
  # spread, which is that case exactly.
  set.seed(7)
  m  <- 4L
  C  <- matrix(0, 1L, m)                    # one node: no between-node spread
  w  <- 1
  Dv <- matrix(stats::runif(m, 0.05, 0.4), 1L, m)
  S  <- diag(as.numeric(Dv))
  expect_equal(admixr2:::.admAdfWeightFast(C, w, Dv, 300),
               admixr2:::.admAdfWeightNormal(S, 300), tolerance = 1e-12)
})

test_that("the mean block is Vt/N under either weight", {
  # Var(ybar) = Vt/N holds without normality -- it is the CLT, not the
  # assumption -- so this block is the one the correction must NOT move.
  set.seed(9)
  m <- 3L; Q <- 11L
  C  <- matrix(stats::rnorm(Q * m), Q, m); C <- sweep(C, 2L, colMeans(C))
  w  <- stats::runif(Q); w <- w / sum(w)
  Dv <- matrix(stats::runif(Q * m, 0.02, 0.2), Q, m)
  S  <- crossprod(C, w * C); diag(S) <- diag(S) + colSums(w * Dv)
  W  <- admixr2:::.admAdfWeightFast(C, w, Dv, 250)
  expect_equal(W[seq_len(m), seq_len(m)], S / 250, tolerance = 1e-12)
})
