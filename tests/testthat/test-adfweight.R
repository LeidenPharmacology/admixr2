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

test_that("the information equality holds: J(working weight) == 2H", {
  # THE test for the sandwich. J = sum_s G_s W_s G_s' with W the weight the
  # objective actually uses must equal 2H exactly at a well-specified fixture.
  # It fails loudly if either half is wrong, and it caught three real bugs:
  #
  #   * a GLS-surrogate bread and dtau/dPsi as G. Eq. (1) is GLS on t with the
  #     normal-theory weight only ASYMPTOTICALLY -- F is linear in V and
  #     quadratic in ybar, a GLS criterion is quadratic in both -- so the two
  #     agree at t = tau and nowhere else.
  #   * full-matrix derivatives on the `var` branch, whose objective is the
  #     DIAGONAL one, so dF/dV_ii = N/Vt_ii rather than N (Vt^-1)_ii.
  #   * the marginal of the normal-theory weight as the var baseline, where the
  #     objective in fact assumes WORKING INDEPENDENCE across the m variances.
  skip_if_not_installed("rxode2")
  TIMES <- c(4, 8, 12, 16); DOSE <- 100; NQ <- 9L
  .mod <- function() {
    ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16; prop.err <- 0.15 })
    model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt()
            cp ~ prop(prop.err) })
  }
  ui <- suppressMessages(rxode2::rxode2(.mod))
  ov <- admixr2:::.admOutputVar(ui); rx <- admixr2:::.admLoadModel(ui)
  E0 <- DOSE / 10 * exp(-0.1 * TIMES)
  mk <- function(E, V, N, meth)
    list(s = list(E = E, V = if (meth == "var") diag(V) else V, n = N,
                  times = TIMES, ev = rxode2::et(amt = DOSE)))
  for (meth in c("var", "cov")) {
    N   <- 100L
    st  <- mk(E0, diag((0.3 * E0)^2), N, meth)
    ctl <- adghControl(studies = st, grad = "none", n_nodes = NQ, print = 0L,
                       covMethod = "none")
    pin <- admixr2:::.admDriverPinfo(ui, ctl)
    u   <- admixr2:::.admDriverUnits(st, ui, ov)
    g   <- admixr2:::.adghNodeGrid(NQ, pin$n_eta)
    p   <- admixr2:::.admBuildOptVec(pin)$p0
    pars <- admixr2:::.admUnpack(p, pin)
    # a WELL-SPECIFIED fixture: the summary IS what the model predicts, t = tau
    pt  <- admixr2:::.admAdfParts(pars, pin, u$studies[[1L]], rx, ov, g, 1L)
    u2  <- admixr2:::.admDriverUnits(mk(as.numeric(pt$E), pt$V, N, meth), ui, ov)
    s1  <- u2$studies[[1L]]
    f   <- function(q) admixr2:::.adghNLL(q, pin, u2$studies, rx, ov, g, 1L)
    H   <- numDeriv::hessian(f, p)
    md  <- admixr2:::.admMomentDeriv(p, pin, u2$studies, rx, ov, g, 1L)
    p2  <- admixr2:::.admAdfParts(pars, pin, s1, rx, ov, g, 1L)
    G   <- admixr2:::.admScoreCross(md$E[[1L]], md$V[[1L]], md$dE[[1L]],
                                    md$dV[[1L]], s1, N)
    W   <- admixr2:::.admWorkingWeight(p2$V, N, s1$method)
    J   <- G %*% W %*% t(G)
    ev  <- sort(Re(eigen(solve(2 * H) %*% J)$values))
    expect_equal(ev, rep(1, length(ev)), tolerance = 1e-4, info = meth)
  }
})

test_that("the true weight moves the sandwich away from 2H", {
  # ... and the correction must not be trivial, or the machinery is inert.
  set.seed(12)
  m <- 4L; Q <- 21L
  C  <- matrix(stats::rnorm(Q * m), Q, m); C <- sweep(C, 2L, colMeans(C))
  w  <- stats::runif(Q); w <- w / sum(w)
  Dv <- matrix(stats::runif(Q * m, 0.02, 0.2), Q, m)
  S  <- crossprod(C, w * C); diag(S) <- diag(S) + colSums(w * Dv)
  Om <- admixr2:::.admAdfWeightFast(C, w, Dv, 150)
  Wn <- admixr2:::.admAdfWeightNormal(S, 150)
  # the mean block is the CLT and must be untouched by the correction
  expect_equal(Om[seq_len(m), seq_len(m)], Wn[seq_len(m), seq_len(m)],
               tolerance = 1e-12)
  # the covariance block is not
  expect_gt(max(abs(Om[-seq_len(m), -seq_len(m)] /
                    Wn[-seq_len(m), -seq_len(m)] - 1)), 0.05)
  # and the cross block, zero under normality, is not zero here
  expect_gt(max(abs(Om[seq_len(m), -seq_len(m)])), 0)
})
