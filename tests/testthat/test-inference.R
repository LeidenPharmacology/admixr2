# Inference under misspecification: G at tau, the LRT correction, TIC.
# See algorithm/adf/HANDOFF-INFERENCE.md.

test_that("G is evaluated at tau, so it does not move with the observed mean", {
  # J is DEFINED as Var(S), and the expansion of the score about t = tau leaves
  # G evaluated AT tau. Building it from the realised residual instead adds
  # E[K delta Omega delta' K'] -- a quadratic form, so non-negative -- biasing J
  # upward by O(1/N). Measured +0.33% at omega = 0.2 rising to +2.9% with a
  # proportional residual, over 200,000 paired replicates.
  #
  # The whole of that fix is that `r` is zero, and this is the one-line assertion
  # that says so: perturb the study mean and G must not budge.
  m  <- 4L; N <- 120
  set.seed(5)
  E  <- 10 * exp(-0.1 * seq_len(m))
  V  <- diag((0.3 * E)^2) + 0.02
  dE <- matrix(stats::rnorm(m * 3L), m, 3L)
  dV <- lapply(seq_len(3L), function(k) {
    A <- matrix(stats::rnorm(m * m), m, m); (A + t(A)) / 20
  })
  for (meth in c("cov", "var")) {
    s1 <- list(E = E,               V = V, n = N, method = meth)
    s2 <- list(E = E * 1.25 + 0.7,  V = V, n = N, method = meth)  # wildly off
    G1 <- admixr2:::.admScoreCross(E, V, dE, dV, s1, N)
    G2 <- admixr2:::.admScoreCross(E, V, dE, dV, s2, N)
    expect_equal(G1, G2, info = meth)
  }
})

test_that("Ruben's series is exact when the weights are equal", {
  # Self-checking by construction: with all lambda equal every g_r is zero, so
  # every correction term vanishes and the series collapses to a single
  # chi-squared. If this is not EXACT, the recursion is wrong.
  for (q in 1:4) for (lam in c(0.5, 1, 3.7)) {
    l <- rep(lam, q)
    for (x in c(0.1, 1, 4.3, 12)) {
      expect_equal(admixr2:::.admRubenP(x, l), stats::pchisq(x / lam, q),
                   tolerance = 1e-12, info = paste(q, lam, x))
    }
  }
})

test_that("Ruben's series matches Monte Carlo for spread weights", {
  # The q = 2 weights from the multiparameter study were 2.21 and 0.90 -- the
  # coefficient on the parameter carrying the random effect needs the full
  # correction, the one without it needs almost none. A single averaged weight
  # is not good enough, which is why the eigenvalues are used rather than a
  # trace whenever the block is well conditioned.
  set.seed(11)
  lam <- c(2.21, 0.90)
  Z   <- matrix(stats::rnorm(2L * 4e5L)^2, 2L)
  q   <- colSums(lam * Z)
  for (x in c(1.5, 4.0, 9.0))
    expect_equal(admixr2:::.admRubenP(x, lam), mean(q < x), tolerance = 3e-3,
                 info = x)
})

# A minimal object with the shape .admLRT reads: env$admExtra$sandwich and
# env$admExtra$par_names, plus an objective.
.mock_fit <- function(H, J, nms, objective) {
  structure(list(
    objective = objective,
    env = list(admExtra = list(
      sandwich  = list(H = H, J = J, par_names = nms),
      par_names = nms))), class = "admFit")
}

test_that("the tested block is subset AFTER inverting, not before", {
  # The inverse's gamma block is the Schur complement, which PROFILES OUT the
  # nuisance parameters; inverting H[gamma, gamma] conditions on them instead.
  # Measured, the wrong order inflated sum(lambda) 12x (4.66 -> 58.4) -- in the
  # CONSERVATIVE direction, so the test rejected at 0.0002 instead of 0.05 and
  # no covariate would ever have been selected. Silent either way.
  set.seed(3)
  p <- 4L; nms <- c("tcl", "tv", "tcov", "err")
  A0 <- matrix(stats::rnorm(p * p), p, p)
  H  <- crossprod(A0) + diag(p)                     # PD, and NOT diagonal
  B0 <- matrix(stats::rnorm(p * p), p, p)
  J  <- crossprod(B0) + diag(p)
  k  <- 3L                                          # test "tcov"
  Hi <- solve(H)
  right <- (2 * Hi)[k, k, drop = FALSE]
  wrong <- 2 * solve(H[k, k, drop = FALSE])
  # the two orders genuinely differ, or this test proves nothing
  expect_false(isTRUE(all.equal(as.numeric(right), as.numeric(wrong))))

  full <- .mock_fit(H, J, nms, 100)
  red  <- .mock_fit(H[-k, -k], J[-k, -k], nms[-k], 110)
  got  <- admixr2:::.admLRT(full, red)
  B    <- (Hi %*% J %*% Hi)[k, k, drop = FALSE]
  expect_equal(got$lambda, as.numeric(B / right))
  expect_false(isTRUE(all.equal(got$lambda, as.numeric(B / wrong))))
})

test_that("correct specification reduces the test to the ordinary chi-squared", {
  # J = 2H is the information equality. Then sandwich = 2H^-1 = naive, every
  # eigenvalue is 1, and the corrected p must equal pchisq(dOFV, q) to machine
  # precision -- not approximately. Holds at q = 1 (scaled path) and q > 1
  # (Ruben path), which are different branches.
  set.seed(9)
  p <- 5L; nms <- paste0("p", seq_len(p))
  A0 <- matrix(stats::rnorm(p * p), p, p)
  H  <- crossprod(A0) + diag(p)
  full <- .mock_fit(H, 2 * H, nms, 100)
  for (drop in list(5L, c(4L, 5L), c(3L, 4L, 5L))) {
    red <- .mock_fit(H[-drop, -drop], 2 * H[-drop, -drop], nms[-drop], 100 + 6.5)
    got <- admixr2:::.admLRT(full, red)
    expect_equal(got$lambda, rep(1, length(drop)), tolerance = 1e-10)
    expect_equal(got$p, stats::pchisq(6.5, length(drop), lower.tail = FALSE),
                 tolerance = 1e-10)
    expect_equal(got$p, got$naive_p, tolerance = 1e-10)
  }
})

test_that("anova() refuses what it cannot test rather than returning a p-value", {
  set.seed(4)
  p <- 3L; nms <- c("a", "b", "c")
  A0 <- matrix(stats::rnorm(p * p), p, p); H <- crossprod(A0) + diag(p)
  full <- .mock_fit(H, 2 * H, nms, 100)
  # non-nested: the smaller model has a parameter the larger one does not
  odd <- .mock_fit(H[1:2, 1:2], 2 * H[1:2, 1:2], c("a", "z"), 105)
  expect_error(admixr2:::.admLRT(full, odd), "not nested")
  # identical parameter sets
  expect_error(admixr2:::.admLRT(full, .mock_fit(H, 2 * H, nms, 100)),
               "nothing to test")
  # dOFV < 0 is impossible under nesting and means the LARGER fit failed. It
  # must not be clamped to zero, which would turn a failed fit into a silent
  # non-rejection.
  bad <- .mock_fit(H[1:2, 1:2], 2 * H[1:2, 1:2], c("a", "b"), 90)
  expect_error(admixr2:::.admLRT(full, bad), "nesting forbids")
  # and the correction needs the larger fit's H and J
  nosw <- structure(list(objective = 100,
                         env = list(admExtra = list(par_names = nms))),
                    class = "admFit")
  expect_error(admixr2:::.admLRT(nosw, .mock_fit(H[1:2, 1:2], 2 * H[1:2, 1:2],
                                                 nms[1:2], 105)),
               "covMethod")
  # anova() itself needs a pair
  expect_error(anova(full), "nothing to compare")
})

test_that("TIC reduces to AIC exactly under correct specification", {
  # THE convention check, and the one that catches a factor of two: J = 2H gives
  # tr(H^-1 J) = 2p, so p_eff = p and TIC = objective + 2p = AIC.
  set.seed(2)
  for (p in c(1L, 3L, 7L)) {
    A0 <- matrix(stats::rnorm(p * p), p, p); H <- crossprod(A0) + diag(p)
    st <- admixr2:::.admTICStats(list(H = H, J = 2 * H), 250)
    expect_equal(st$p_eff, p, tolerance = 1e-9, info = p)
    expect_equal(st$TIC, 250 + 2 * p, tolerance = 1e-9, info = p)
  }
  # and a misspecified J moves the penalty UP -- the direction that matters,
  # because it means AIC systematically OVER-selects here.
  p  <- 4L
  A0 <- matrix(stats::rnorm(p * p), p, p); H <- crossprod(A0) + diag(p)
  st <- admixr2:::.admTICStats(list(H = H, J = 5 * H), 250)
  expect_gt(st$p_eff, p)
  expect_null(admixr2:::.admTICStats(NULL, 250))
})

# BIC and BIC_h are pinned in test-utils.R, beside .admCalcObjStats itself.

test_that("a boundary restriction is refused, not rescaled", {
  # Fixing a VARIANCE to zero puts the null on the edge of the parameter space,
  # where the limit is a chi-bar-squared mixture. The eigenvalue weights assume
  # an interior null, so rescaling dOFV there returns a finite, plausible
  # p-value from the wrong reference distribution. The roxygen claimed this was
  # refused before the code did it.
  set.seed(6)
  p <- 4L
  nms <- c("tcl", "tv", "logchol_eta.cl", "chol_eta.v_eta.cl")
  A0 <- matrix(stats::rnorm(p * p), p, p); H <- crossprod(A0) + diag(p)
  full <- .mock_fit(H, 2 * H, nms, 100)
  # dropping the VARIANCE is a boundary restriction -> refuse
  red_v <- .mock_fit(H[-3L, -3L], 2 * H[-3L, -3L], nms[-3L], 106)
  expect_error(admixr2:::.admLRT(full, red_v), "BOUNDARY")
  expect_error(anova(full, red_v), "BOUNDARY")
  # dropping the COVARIANCE is interior -- it may take either sign -- and is
  # tested normally
  red_c <- .mock_fit(H[-4L, -4L], 2 * H[-4L, -4L], nms[-4L], 106)
  got <- admixr2:::.admLRT(full, red_c)
  expect_equal(got$df, 1L)
  expect_true(is.finite(got$p))
})
