# admixr2's transform functions, pinned to rxode2's EXPORTED ones.
#
# .admTBS/.admTBSi are now thin wrappers over `rxode2::.rxTransform()` -- the same
# entry point rxode2's own boxCox()/yeoJohnson()/logit()/probit() use, bottoming
# out in `.Call(_rxode2_powerD, ...)`. They used to be a line-by-line R port of
# that C code, and this file existed to catch the port DRIFTING from it. That
# specific risk is gone; the file is kept, and still earns its place, because the
# wrappers add two things of their own that can break:
#
#   * dim() restoration -- `.rxTransform()` drops it, and a dropped dim turned
#     cp_mat into a flat vector so every downstream colMeans()/sweep() gave NA;
#   * the yj code mapping (0 boxCox, 1 yeoJohnson, 2 untransformed, 4 logit,
#     6 probit) -- if rxode2 ever renumbers those, admixr2 would silently apply
#     the WRONG transform, and these comparisons against the named functions
#     (boxCox/yeoJohnson/logit/probit and their inverses) are what would catch it.
#
# .admTBSid -- the derivative of the INVERSE -- remains admixr2's own: rxode2
# exposes no equivalent (`_powerDD` is the derivative of the FORWARD transform,
# and composing 1/_powerDD(g(z)) blows up to 6.7e7 exactly where the 81-node grid's
# +-12 SD tails land), and rxode2's `_powerDD` additionally has a sign error on the
# Yeo-Johnson negative branch which admixr2 deliberately does not reproduce. It is
# checked below against a central difference of rxode2's OWN inverse.

test_that("forward transforms match rxode2's exported ones", {
  skip_if_not_installed("rxode2")
  y <- c(0.05, 0.3, 1, 2.5, 7, 19.5)                  # inside every support used
  for (lam in c(-0.5, 0, 0.3, 0.5, 1, 1.5, 2, 2.5)) {
    expect_equal(.admTBS(y, lam, 0L, 0, 1), rxode2::boxCox(y, lam),
                 tolerance = 1e-12, info = paste("boxCox lambda =", lam))
    expect_equal(.admTBS(y, lam, 1L, 0, 1), rxode2::yeoJohnson(y, lam),
                 tolerance = 1e-12, info = paste("yeoJohnson lambda =", lam))
  }
  # Yeo-Johnson is defined for negative inputs too -- the branch boxCox lacks.
  yn <- c(-9, -2.5, -0.4, 0, 0.4, 2.5, 9)
  for (lam in c(-0.5, 0, 0.5, 1, 2, 2.5))
    expect_equal(.admTBS(yn, lam, 1L, 0, 1), rxode2::yeoJohnson(yn, lam),
                 tolerance = 1e-12, info = paste("yeoJohnson neg lambda =", lam))
  # logit/probit: bounds come from predDf$trLow/trHi, so vary them.
  for (b in list(c(0, 20), c(-1, 1), c(2, 40))) {
    expect_equal(.admTBS(seq(b[1] + 0.01, b[2] - 0.01, length.out = 9), 1, 4L, b[1], b[2]),
                 rxode2::logit(seq(b[1] + 0.01, b[2] - 0.01, length.out = 9), b[1], b[2]),
                 tolerance = 1e-10, info = paste("logit", b[1], b[2]))
    expect_equal(.admTBS(seq(b[1] + 0.01, b[2] - 0.01, length.out = 9), 1, 6L, b[1], b[2]),
                 rxode2::probit(seq(b[1] + 0.01, b[2] - 0.01, length.out = 9), b[1], b[2]),
                 tolerance = 1e-10, info = paste("probit", b[1], b[2]))
  }
})

test_that("inverse transforms match rxode2's exported ones", {
  skip_if_not_installed("rxode2")
  z <- c(-3.2, -1, -0.2, 0, 0.4, 1.7, 4)
  for (lam in c(0, 0.3, 0.5, 1, 1.5, 2, 2.5)) {
    # boxCox's inverse has bounded support (lam*z + 1 > 0); rxode2 clamps the
    # RESULT to sqrt(DBL_EPSILON) there and so does admixr2, so compare on the
    # whole grid including the out-of-support tail.
    expect_equal(.admTBSi(z, lam, 0L, 0, 1), rxode2::boxCoxInv(z, lam),
                 tolerance = 1e-10, info = paste("boxCoxInv lambda =", lam))
    expect_equal(.admTBSi(z, lam, 1L, 0, 1), rxode2::yeoJohnsonInv(z, lam),
                 tolerance = 1e-10, info = paste("yeoJohnsonInv lambda =", lam))
  }
  for (b in list(c(0, 20), c(-1, 1), c(2, 40))) {
    expect_equal(.admTBSi(z, 1, 4L, b[1], b[2]), rxode2::expit(z, b[1], b[2]),
                 tolerance = 1e-12, info = paste("expit", b[1], b[2]))
    expect_equal(.admTBSi(z, 1, 6L, b[1], b[2]), rxode2::probitInv(z, b[1], b[2]),
                 tolerance = 1e-12, info = paste("probitInv", b[1], b[2]))
  }
})

test_that("the inverse derivative matches a finite difference of rxode2's inverse", {
  skip_if_not_installed("rxode2")
  # .admTBSid is g'(z). rxode2 has no g', so the reference is a central difference
  # of rxode2's OWN inverse -- an independent check that does not reuse any
  # admixr2 code.
  fd <- function(f, z, h = 1e-6) (f(z + h) - f(z - h)) / (2 * h)
  z  <- c(-1.4, -0.5, 0.2, 0.9, 2.1)
  for (lam in c(0.3, 0.5, 1, 1.5)) {
    expect_equal(.admTBSid(z, lam, 0L, 0, 1),
                 fd(function(q) rxode2::boxCoxInv(q, lam), z),
                 tolerance = 1e-5, info = paste("d boxCoxInv lambda =", lam))
    expect_equal(.admTBSid(z, lam, 1L, 0, 1),
                 fd(function(q) rxode2::yeoJohnsonInv(q, lam), z),
                 tolerance = 1e-5, info = paste("d yeoJohnsonInv lambda =", lam))
  }
  expect_equal(.admTBSid(z, 1, 4L, 0, 20), fd(function(q) rxode2::expit(q, 0, 20), z),
               tolerance = 1e-5)
  expect_equal(.admTBSid(z, 1, 6L, 0, 20), fd(function(q) rxode2::probitInv(q, 0, 20), z),
               tolerance = 1e-5)
})

test_that("transforms preserve dim() and handle non-finite input like rxode2", {
  skip_if_not_installed("rxode2")
  # A dropped dim turned cp_mat into a flat vector and produced NA downstream;
  # an NA input used to raise "NAs are not allowed in subscripted assignments"
  # where rxode2's `if (!R_finite(x)) return NA_REAL` gives NA.
  m <- matrix(c(0.5, 1, 2, 3, 4, 5), 3, 2)
  for (yj in c(0L, 1L, 2L, 4L, 6L)) for (lam in c(0, 0.5, 1, 2)) {
    for (f in list(.admTBS, .admTBSi, .admTBSid)) {
      out <- f(m, lam, yj, 0, 40)
      expect_identical(dim(out), c(3L, 2L),
                       info = sprintf("yj=%d lambda=%g", yj, lam))
    }
  }
  for (yj in c(0L, 1L, 4L)) {
    expect_true(is.na(.admTBSi(c(1, NA), 0.5, yj, 0, 40)[2L]))
    expect_true(is.na(.admTBS(c(1, Inf), 0.5, yj, 0, 40)[2L]))
  }
})

test_that("_eps matches rxode2's, and is not confused with DBL_EPSILON", {
  skip_if_not_installed("rxode2")
  # rxode2.h: `#define _eps sqrt(DBL_EPSILON)` for the transform clamps, but
  # safeLog/safeZero/safePow use DBL_EPSILON itself. Conflating them would move
  # every out-of-support value by eight orders of magnitude.
  expect_identical(.ADM_EPS, sqrt(.Machine$double.eps))
  # Compare the RATIO, not the values: all.equal() falls back to an absolute
  # comparison for tiny numbers, so all.equal(sqrt(eps), eps) is TRUE even though
  # they differ by eight orders of magnitude.
  expect_gt(.ADM_EPS / .Machine$double.eps, 1e7)
})
