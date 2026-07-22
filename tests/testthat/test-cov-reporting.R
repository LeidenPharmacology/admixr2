# The REPORTED covariance: scale, names, and surviving nlmixr2est.
#
# The Hessian is built on the optimizer's parameterisation (log thetas,
# log(sigma^2), log-Cholesky omega) but `fit$cov` sits next to `Estimate` in
# nlmixr2est's parFixed table, where `Estimate +- 1.96*SE` has to mean something.
# So the block is delta-transformed on the way out. Sigma is a per-row factor;
# omega is NOT -- Omega = L L' is dense in the Cholesky once omega is correlated,
# and using a per-row factor there would be silently wrong only for correlated
# models, i.e. exactly where nobody would notice.
#
# Tier 1: no rxode2. .admParseIniDf() takes a plain iniDf, so the whole Jacobian
# can be checked against a finite difference of .admUnpack() itself.

test_that(".admOmegaJacobian matches a finite difference of Omega(p), 1 eta", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_1eta(omega = 0.09))
  p0    <- admixr2:::.admBuildOptVec(pinfo)$p0
  oi    <- length(pinfo$struct_names) + length(pinfo$sigma_names) +
             seq_along(pinfo$omega_par)

  omvec <- function(pv) {
    L <- admixr2:::.admUnpack(pv, pinfo)$L
    Om <- L %*% t(L)
    vapply(seq_along(pinfo$omega_par),
           function(r) Om[pinfo$chol_i[r], pinfo$chol_j[r]], double(1))
  }
  J   <- admixr2:::.admOmegaJacobian(pinfo, admixr2:::.admUnpack(p0, pinfo)$L)
  h   <- 1e-6
  Jfd <- vapply(seq_along(oi), function(k) {
    pp <- p0; pm <- p0
    pp[oi[k]] <- pp[oi[k]] + h; pm[oi[k]] <- pm[oi[k]] - h
    (omvec(pp) - omvec(pm)) / (2 * h)
  }, double(length(oi)))
  dim(Jfd) <- dim(J)

  expect_equal(J, Jfd, tolerance = 1e-6)
  # p = log(Omega_11), so d(Omega_11)/dp = Omega_11 exactly. A d(L_11)/dp chain
  # rule of 1 instead of L_11/2 would give twice this.
  expect_equal(J[1, 1], 0.09, tolerance = 1e-10)
})

test_that(".admOmegaJacobian matches a finite difference for a CORRELATED omega", {
  # The case a per-row delta factor cannot represent: J is lower-triangular here,
  # not diagonal, because moving L_21 moves both Omega_21 and Omega_22.
  pinfo <- admixr2:::.admParseIniDf(make_inidf_2eta())
  p0    <- admixr2:::.admBuildOptVec(pinfo)$p0
  oi    <- length(pinfo$struct_names) + length(pinfo$sigma_names) +
             seq_along(pinfo$omega_par)

  omvec <- function(pv) {
    L <- admixr2:::.admUnpack(pv, pinfo)$L
    Om <- L %*% t(L)
    vapply(seq_along(pinfo$omega_par),
           function(r) Om[pinfo$chol_i[r], pinfo$chol_j[r]], double(1))
  }
  J   <- admixr2:::.admOmegaJacobian(pinfo, admixr2:::.admUnpack(p0, pinfo)$L)
  h   <- 1e-6
  Jfd <- vapply(seq_along(oi), function(k) {
    pp <- p0; pm <- p0
    pp[oi[k]] <- pp[oi[k]] + h; pm[oi[k]] <- pm[oi[k]] - h
    (omvec(pp) - omvec(pm)) / (2 * h)
  }, double(length(oi)))
  dim(Jfd) <- dim(J)

  expect_equal(J, Jfd, tolerance = 1e-6)
  expect_true(any(abs(J[upper.tri(J)]) + abs(J[lower.tri(J)]) > 0),
              info = "a correlated omega Jacobian must not be diagonal")
})

test_that(".admOmegaJacobian returns NULL when there is no omega", {
  pinfo <- admixr2:::.admParseIniDf(make_inidf_0eta())
  expect_null(admixr2:::.admOmegaJacobian(pinfo, matrix(0, 0, 0)))
})

test_that("omega report names follow nlmixr2est's convention, not the optimizer's", {
  # nlmixr2est's .foceiOmegaCovNames(): `om.<eta>` on the diagonal,
  # `cov.<eta_i>.<eta_j>` off it with i the LATER eta. The optimizer-scale names
  # (logchol_*) describe the Cholesky and must never reach the user.
  p1 <- admixr2:::.admParseIniDf(make_inidf_1eta())
  expect_identical(admixr2:::.admOmegaReportNames(p1), "om.eta.cl")

  p2 <- admixr2:::.admParseIniDf(make_inidf_2eta())
  nm <- admixr2:::.admOmegaReportNames(p2)
  expect_setequal(nm, c("om.eta.cl", "om.eta.v", "cov.eta.v.eta.cl"))
  # Same order as omega_par / the Hessian rows -- a permutation here would label
  # every SE with the wrong parameter.
  expect_identical(length(nm), length(p2$omega_par))
  for (r in seq_along(nm)) {
    i <- p2$chol_i[r]; j <- p2$chol_j[r]
    expect_identical(nm[r], if (i == j) paste0("om.", p2$eta_names[i])
                            else paste0("cov.", p2$eta_names[i], ".", p2$eta_names[j]))
  }
  expect_false(any(nm %in% p2$omega_par_names))
})

test_that(".admRestoreCovNames puts names back, and refuses on a shape mismatch", {
  # nlmixr2est's C++ foceiFitCpp_ re-dimnames the covariance from its own theta
  # name vector and blanks the omega rows -- IN PLACE, so the driver's own copy is
  # blanked too and the names have to be snapshot beforehand. This is that repair.
  mk <- function(n, nms) {
    m <- diag(n); dimnames(m) <- list(nms, nms); m
  }
  good <- c("tcl", "tv", "a", "om.eta.cl")
  fit  <- list(env = new.env(parent = emptyenv()))
  fit$env$cov <- mk(4L, c("tcl", "tv", "a", ""))

  admixr2:::.admRestoreCovNames(fit, good)
  expect_identical(rownames(fit$env$cov), good)
  expect_identical(colnames(fit$env$cov), good)

  # Wrong size: leave it alone. A wrongly labelled SE is worse than an
  # unlabelled one.
  fit$env$cov <- mk(3L, c("tcl", "tv", "a"))
  admixr2:::.admRestoreCovNames(fit, good)
  expect_identical(rownames(fit$env$cov), c("tcl", "tv", "a"))

  # Nothing to do / nothing to do it to.
  expect_silent(admixr2:::.admRestoreCovNames(fit, NULL))
  expect_silent(admixr2:::.admRestoreCovNames(list(env = new.env()), good))
})

test_that("the omega entry ORDER matches rxode2's own (row, col) enumeration", {
  skip_if_not_installed("rxode2")
  # admixr2 cannot reuse nlmixr2est's covariance machinery for omega -- see the
  # note in .admOmegaJacobian() -- but it MUST agree with it on which matrix entry
  # row k of the block refers to, or every omega SE is attached to the wrong
  # parameter and nothing about the numbers looks wrong.
  #
  # rxode2::rxOmegaVarCovDeriv() is EXPORTED and enumerates the free entries of a
  # symmetric Omega in `$elements` as (row, col); nlmixr2est's own
  # .foceiOmegaPairs()/.foceiOmegaCovNames() build the `om.`/`cov.` names off that
  # same enumeration. Pinning to the exported function keeps the convention tied
  # to upstream without a ::: dependency (admixr2 has none, deliberately).
  pinfo <- admixr2:::.admParseIniDf(make_inidf_2eta())
  Om    <- admixr2:::.admUnpack(admixr2:::.admBuildOptVec(pinfo)$p0, pinfo)$omega
  el    <- rxode2::rxOmegaVarCovDeriv(unname(Om), order = 1L)$elements

  # The SET of entries must match. The ORDER need not, and does not in general:
  # rxode2 enumerates the lower triangle column-major ((1,1), (2,1), (2,2)) while
  # admixr2 follows iniDf ROW order, which is what .admFullTheta() reports in. A
  # model written the usual way (`eta.cl + eta.v ~ c(...)`) comes out of rxode2 in
  # rxode2's order anyway; a hand-built iniDf that lists both diagonals before the
  # off-diagonal does not. That is safe ONLY because every row of the block is
  # named from its OWN (chol_i, chol_j) -- which is what the rest of this test
  # pins down. Never index the omega block positionally against rxode2's list.
  ours <- paste(pinfo$chol_i, pinfo$chol_j, sep = "-")
  theirs <- paste(el[, "row"], el[, "col"], sep = "-")
  expect_setequal(ours, theirs)
  expect_equal(length(ours), length(theirs))
  # Lower triangle in both conventions, so `cov.<eta_i>.<eta_j>` has i the LATER
  # eta exactly as .foceiOmegaCovNames() produces it.
  expect_true(all(pinfo$chol_i >= pinfo$chol_j))
  expect_true(all(el[, "row"] >= el[, "col"]))
})

test_that(".admSigmaReportJac is the derivative of the map .admFullTheta reports", {
  # These factors used to be a switch() copied into each of the three CalcCov
  # functions, checked by nothing. They convert an optimizer-scale SE onto the
  # scale the Estimate is printed on, so getting one wrong does not fail loudly --
  # it prints a confident interval of the wrong width. For a "var" sigma the
  # factor is 2/a, so an SD of 0.1 with a missing factor prints an SE 20x too big.
  #
  # The reference is a finite difference of .admSigmaNat() itself (post-composed
  # with sqrt for the "var" role, since .admFullTheta reports that one as an SD),
  # so the test does not reuse the formulas it is checking.
  reported <- function(p, pinfo) {
    nat  <- admixr2:::.admSigmaNat(p, pinfo)
    role <- admixr2:::.admSigmaRole(pinfo)
    ifelse(role == "var", sqrt(nat), nat)
  }
  h <- 1e-6
  for (role in c("var", "t_df", "ar_cor", "nb_size", "pow_exp", "tbs_lam")) {
    pinfo <- list(sigma_names = "s", sigma_role = role)
    # A value in the sensible interior of each role's range.
    p <- switch(role, var = 2 * log(0.35), t_df = log(6 - 2),
                ar_cor = log(0.6 / 0.4), nb_size = log(4), 0.7)
    ana <- admixr2:::.admSigmaReportJac(p, pinfo)
    fd  <- (reported(p + h, pinfo) - reported(p - h, pinfo)) / (2 * h)
    expect_equal(unname(ana), unname(fd), tolerance = 1e-6, info = role)
  }
  # Absent sigma_role (legacy/hand-built pinfo) must behave as "var", not as 1.
  legacy <- list(sigma_names = "s")
  expect_equal(admixr2:::.admSigmaReportJac(2 * log(0.35), legacy),
               0.35 / 2, tolerance = 1e-10)
  expect_identical(admixr2:::.admSigmaReportJac(numeric(0), legacy), numeric(0))
})

test_that(".admCovNames snapshots row names, and tolerates a NULL covariance", {
  m <- diag(2); dimnames(m) <- list(c("a", "b"), c("a", "b"))
  expect_identical(admixr2:::.admCovNames(m), c("a", "b"))
  expect_null(admixr2:::.admCovNames(NULL))
  expect_null(admixr2:::.admCovNames("not a matrix"))
})
