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
  skip_on_cran()
  skip_if_not_installed("rxode2")
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
  skip_if_not_installed("numDeriv")   # used for the reference Hessian
  TIMES <- c(4, 8, 12, 16); DOSE <- 100; NQ <- 9L
  .mod <- function() {
    ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16; prop.err <- 0.15 })
    model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt()
            cp ~ prop(prop.err) })
  }
  ui <- suppressMessages(rxode2::rxode2(.mod))
  ov <- admixr2:::.admOutputVar(ui)
  # SENS BEFORE SIM, always: .admLoadModel's cache-miss path caches
  # foceiModel$inner = NULL for a linCmt() model, after which the sensitivity
  # equations can never be retrieved -- and the cache is user-wide, so this
  # poisons other sessions too. CLAUDE.md records that there is no recovery.
  sens <- admixr2:::.admLoadSensModel(ui)
  rx <- admixr2:::.admLoadModel(ui)
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
    p2  <- admixr2:::.admAdfParts(pars, pin, s1, rx, ov, g, 1L)
    W   <- admixr2:::.admWorkingWeight(p2$V, N, s1$method)
    # BOTH routes to the moment Jacobian must satisfy it -- the identity is what
    # says the analytic one describes the same objective the FD one does.
    mds <- list(
      fd       = admixr2:::.admMomentDeriv(p, pin, u2$studies, rx, ov, g, 1L),
      analytic = admixr2:::.admMomentJac(p, pin, u2$studies, sens, rx, ov, g, 1L))
    for (route in names(mds)) {
      md <- mds[[route]]
      expect_false(is.null(md), info = paste(meth, route))
      G  <- admixr2:::.admScoreCross(md$E[[1L]], md$V[[1L]], md$dE[[1L]],
                                     md$dV[[1L]], s1, N)
      J  <- G %*% W %*% t(G)
      ev <- sort(Re(eigen(solve(2 * H) %*% J)$values))
      expect_equal(ev, rep(1, length(ev)), tolerance = 1e-4,
                   info = paste(meth, route))
    }
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

test_that("a missing or unusable H is an error, not a silent NULL", {
  # .admSandwich's tryCatch exists for a singular H and cannot distinguish that
  # from an H the caller never passed. Before this guard, a study script calling
  # .admSandwichCov() positionally with seven arguments got NULL for every
  # replicate, skipped them all, and still printed a table -- of the replicates
  # it had not skipped in an earlier, differently-shaped version of the call.
  # Fail fast and say which argument, at the top, before anything is solved.
  expect_error(
    admixr2:::.admSandwichCov(1, NULL, list(), NULL, "cp", NULL, 1L),
    "`H` must be a finite square Hessian")
  expect_error(
    admixr2:::.admSandwichCov(1, NULL, list(), NULL, "cp", NULL, 1L,
                              H = matrix(c(1, NA, NA, 1), 2L)),
    "`H` must be a finite square Hessian")
  expect_error(
    admixr2:::.admSandwichCov(1, NULL, list(), NULL, "cp", NULL, 1L,
                              H = matrix(1, 2L, 3L)),
    "`H` must be a finite square Hessian")
})

test_that("the expansion is exact for a SKEWED conditional residual", {
  # The Wick expansion is often described as needing a conditionally NORMAL
  # residual. It does not: it needs one that is conditionally INDEPENDENT across
  # timepoints, plus its 3rd and 4th central moments. Supplying only the variance
  # and letting the odd pairings collapse is what assumes symmetry, and for a
  # skewed residual that wrecks the cross block Cov(ybar, vech V) -- the block
  # probe 08 found load-bearing -- while leaving Cov(ybar) untouched and so
  # looking healthy from the mean side.
  #
  # Measured against 40k simulated studies: lognormal residual, cross block
  # 0.861 relative error with variance only, 0.027 with the moments; Poisson
  # 0.686 -> 0.047. This is the cheap version of that check.
  set.seed(77)
  TIMES <- c(1, 3, 6); DOSE <- 100; V0 <- 10; OM <- 0.4
  N <- 80L; R <- 6000L; sv <- 0.09
  gh <- admixr2:::.adghNodes1(31L); x <- gh$x; w <- gh$w / sum(gh$w)
  fmat <- function(eta) t(vapply(eta, function(e)
    DOSE / V0 * exp(-exp(e) / V0 * TIMES), numeric(length(TIMES))))
  F  <- fmat(x * OM); m <- ncol(F)
  ij <- which(lower.tri(diag(m), diag = TRUE), arr.ind = TRUE)
  ms <- exp(sv / 2); wl <- exp(sv); M1 <- ms * F
  C  <- sweep(M1, 2L, colSums(w * M1))
  Dv <- M1^2 * (wl - 1)
  T3 <- M1^3 * (wl - 1)^2 * (wl + 2)
  Q4 <- M1^4 * (wl - 1)^2 * (wl^4 + 2 * wl^3 + 3 * wl^2 - 3)

  tv <- matrix(0, R, m + nrow(ij))
  for (r in seq_len(R)) {
    f <- fmat(stats::rnorm(N, 0, OM))
    y <- matrix(stats::rlnorm(length(f), log(f), sqrt(sv)), nrow(f))
    V <- stats::cov.wt(y, method = "ML")$cov
    tv[r, ] <- c(colMeans(y), V[cbind(ij[, 1L], ij[, 2L])])
  }
  ref <- stats::cov(tv)
  rel <- function(A, r, c)
    sqrt(sum((A[r, c] - ref[r, c])^2)) / sqrt(sum(ref[r, c]^2))
  gen <- admixr2:::.admAdfWeightFast(C, w, Dv, N, T3, Q4)
  sym <- admixr2:::.admAdfWeightFast(C, w, Dv, N)
  # with the moments: at the Monte Carlo floor, in every block
  expect_lt(rel(gen, seq_len(m), seq_len(m)),   0.15)
  expect_lt(rel(gen, seq_len(m), -seq_len(m)),  0.20)
  expect_lt(rel(gen, -seq_len(m), -seq_len(m)), 0.20)
  # without them: the mean block is still right and the cross block is not,
  # which is exactly why this cannot be caught from the mean side
  expect_lt(rel(sym, seq_len(m), seq_len(m)),   0.15)
  expect_gt(rel(sym, seq_len(m), -seq_len(m)),  0.50)
})

test_that("a residual the expansion cannot reach is refused, not approximated", {
  skip_if_not_installed("rxode2")
  # ar() correlates the residual ACROSS timepoints, so conditional independence
  # fails and every cross term the expansion drops is real. Returning a weight
  # built as though it held would be a plausible wrong answer.
  arr <- list(form = rep(0L, 3L), a2 = rep(0.1, 3L), b2 = rep(0, 3L),
              cc = rep(1, 3L), vmul = rep(1, 3L), csz = rep(NA_real_, 3L),
              phi = rep(NA_real_, 3L), rho = rep(0.4, 3L))
  expect_null(admixr2:::.admAdfCondMom(matrix(1, 5L, 3L), arr))
  arr$rho <- rep(NA_real_, 3L)
  expect_false(is.null(admixr2:::.admAdfCondMom(matrix(1, 5L, 3L), arr)))
  # a t residual with nu <= 4 has no finite kurtosis: the sampling law of the
  # reported V does not have the variance the weight is made of.
  arr$vmul <- rep(4 / 2, 3L)                      # nu/(nu-2) at nu = 4
  expect_null(admixr2:::.admAdfCondMom(matrix(1, 5L, 3L), arr))
  arr$vmul <- rep(8 / 6, 3L)                      # nu = 8, fine
  expect_false(is.null(admixr2:::.admAdfCondMom(matrix(1, 5L, 3L), arr)))
})

test_that("the weight's own S is the V_pred the objective scores against", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  # .admAdfWeightFast rebuilds S from (C, Dv) via the law of total variance. If
  # that disagrees with V_pred, the weight and the objective describe different
  # laws and every block downstream is scaled wrong -- silently, since both are
  # finite and plausible. On lnorm the conditional mean is f exp(s/2), so an
  # unscaled C broke this while add/prop looked fine.
  skip_if_not_installed("rxode2")
  TIMES <- c(2, 5, 9, 14); DOSE <- 100; NQ <- 11L; N <- 100L
  mods <- list(
    add   = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16; e <- 0.5 })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt(); cp ~ add(e) }) },
    prop  = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16; e <- 0.15 })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt(); cp ~ prop(e) }) },
    lnorm = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16; e <- 0.15 })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt(); cp ~ lnorm(e) }) })
  for (nm in names(mods)) {
    ui <- suppressMessages(rxode2::rxode2(mods[[nm]]))
    ov <- admixr2:::.admOutputVar(ui)
    invisible(admixr2:::.admLoadSensModel(ui))  # sens BEFORE sim; see above
    rx <- admixr2:::.admLoadModel(ui)
    E0 <- DOSE / 10 * exp(-0.1 * TIMES)
    st <- list(s = list(E = E0, V = diag((0.3 * E0)^2), n = N, times = TIMES,
                        ev = rxode2::et(amt = DOSE)))
    ctl <- adghControl(studies = st, grad = "none", n_nodes = NQ, print = 0L,
                       covMethod = "none")
    pin <- admixr2:::.admDriverPinfo(ui, ctl)
    u   <- admixr2:::.admDriverUnits(st, ui, ov)
    g   <- admixr2:::.adghNodeGrid(NQ, pin$n_eta)
    pars <- admixr2:::.admUnpack(admixr2:::.admBuildOptVec(pin)$p0, pin)
    pt  <- admixr2:::.admAdfParts(pars, pin, u$studies[[1L]], rx, ov, g, 1L)
    m   <- length(pt$E)
    W   <- admixr2:::.admAdfWeightFast(pt$C, pt$w, pt$Dv, N, pt$T3, pt$Q4)
    expect_equal(N * W[seq_len(m), seq_len(m)], pt$V, tolerance = 1e-10,
                 ignore_attr = TRUE, info = nm)
  }
})

test_that("the weight's S tracks V_pred to the objective's own accuracy on TBS", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  # For add/prop/lnorm the identity above is exact to ~1e-15. TBS is the one
  # family where it is not, and the gap is the OBJECTIVE's approximation rather
  # than the weight's: .admTBSRow builds the predicted variance as a
  # second-order Taylor expansion in var_f (ev = q$v + 0.5 d2v var_f), while the
  # weight sums the conditional variance exactly over the node ensemble.
  #
  # Measured 2.4e-04 at omega^2 = 0.16 rising to 6.3e-04 at 0.49, on the
  # DIAGONAL only -- the off-diagonal stays at 1e-16 because ms enters it
  # multiplicatively and cancels. That propagates to ~0.06% on a standard
  # error, two orders below anything that matters, but it must not grow
  # silently: a larger gap would mean the weight had stopped describing the law
  # the objective scores against.
  skip_if_not_installed("rxode2")
  TT <- c(2, 5, 9, 14); DOSE <- 100; NQ <- 11L; N <- 200L
  .bc <- function() {
    ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16; e <- 0.5
          lam <- fix(0.4) })
    model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt()
            cp ~ add(e) + boxCox(lam) })
  }
  ui  <- suppressMessages(rxode2::rxode2(.bc))
  ov  <- admixr2:::.admOutputVar(ui)
  invisible(admixr2:::.admLoadSensModel(ui))   # sens BEFORE sim; see above
  rx  <- admixr2:::.admLoadModel(ui)
  pin <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  pin$nDisplayProgress <- .Machine$integer.max; pin$resid_nodes <- 81L
  E0 <- DOSE / 10 * exp(-0.1 * TT)
  st <- list(E = E0, V = diag((0.3 * E0)^2), n = N, times = TT,
             ev = rxode2::et(amt = DOSE))
  u  <- admixr2:::.admFlattenStudies(
          list(s1 = admixr2:::.admNormaliseStudy(st, "s1", "cp")))
  u  <- admixr2:::.admBuildEvFull(u)
  g  <- admixr2:::.adghNodeGrid(NQ, pin$n_eta)
  for (om in c(0.16, 0.49)) {
    p <- admixr2:::.admBuildOptVec(pin)$p0
    p[grep("^logchol", names(p))[1L]] <- log(om)
    pt <- admixr2:::.admAdfParts(admixr2:::.admUnpack(p, pin), pin, u[[1L]],
                                 rx, ov, g, 1L)
    m  <- length(pt$E)
    W  <- admixr2:::.admAdfWeightFast(pt$C, pt$w, pt$Dv, N, pt$T3, pt$Q4)
    S  <- N * W[seq_len(m), seq_len(m)]
    expect_lt(max(abs(S - pt$V) / abs(pt$V)), 2e-3, label = paste("om", om))
    # the off-diagonal is exact: ms enters it multiplicatively and cancels
    expect_lt(max(abs((S - pt$V)[upper.tri(S)]) /
                  abs(pt$V[upper.tri(S)])), 1e-12, label = paste("offdiag", om))
  }
})

test_that("the analytic moment Jacobian matches the finite-difference oracle", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  # G = d2F/(dPsi dt') is closed form in (dE/dPsi, dV/dPsi), so these two are the
  # only derivatives the sandwich takes. .admMomentJac forms them from one
  # sensitivity solve; .admMomentDeriv central-differences the moments and is
  # kept as the oracle precisely because the forward residual composition is a
  # seventh consumer of the residual row arrays -- the place CLAUDE.md says the
  # misses happen.
  #
  # linCmt() throughout on purpose: an ODE model puts the two routes on
  # SEPARATELY COMPILED models (sensitivity vs plain simulation), which take
  # different adaptive steps and disagree by ~7e-5 no matter the FD step. That
  # is a model-consistency difference, not a derivative error, and it would make
  # this test measure the solver instead of the algebra.
  skip_if_not_installed("rxode2")
  TT <- c(2, 5, 9, 14); DOSE <- 100; NQ <- 9L; N <- 100L
  mods <- list(
    add     = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16; e <- 0.5 })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt(); cp ~ add(e) }) },
    prop    = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16; e <- 0.15 })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt(); cp ~ prop(e) }) },
    lnorm   = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16; e <- 0.15 })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt(); cp ~ lnorm(e) }) },
    addprop = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16
                                 a <- 0.2; b <- 0.1 })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt()
              cp ~ add(a) + prop(b) }) },
    eta2    = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16
                                 eta.v ~ 0.04; e <- 0.15 })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v); cp <- linCmt()
              cp ~ prop(e) }) })
  E0 <- DOSE / 10 * exp(-0.1 * TT); Vv <- diag((0.25 * E0)^2)
  for (nm in names(mods)) {
    ui   <- suppressMessages(rxode2::rxode2(mods[[nm]]))
    ov   <- admixr2:::.admOutputVar(ui)
    sens <- admixr2:::.admLoadSensModel(ui); rx <- admixr2:::.admLoadModel(ui)
    for (meth in c("cov", "var")) {
      st <- list(
        s1 = list(E = E0, V = if (meth == "cov") Vv else diag(Vv), n = N,
                  times = TT, ev = rxode2::et(amt = DOSE)),
        s2 = list(E = 1.1 * E0, V = if (meth == "cov") 1.2 * Vv else diag(1.2 * Vv),
                  n = 60L, times = TT, ev = rxode2::et(amt = DOSE)))
      ctl <- adghControl(studies = st, grad = "none", n_nodes = NQ, print = 0L,
                         covMethod = "none")
      pin <- admixr2:::.admDriverPinfo(ui, ctl)
      u   <- admixr2:::.admDriverUnits(st, ui, ov)
      g   <- admixr2:::.adghNodeGrid(NQ, pin$n_eta)
      p   <- admixr2:::.admBuildOptVec(pin)$p0 + 0.05      # off the initial point
      fd  <- admixr2:::.admMomentDeriv(p, pin, u$studies, rx, ov, g, 1L)
      an  <- admixr2:::.admMomentJac(p, pin, u$studies, sens, rx, ov, g, 1L)
      expect_false(is.null(an), info = paste(nm, meth))
      for (i in seq_along(fd$dE)) {
        expect_equal(an$dE[[i]], fd$dE[[i]], tolerance = 1e-6,
                     info = paste(nm, meth, "dE", i))
        for (k in seq_along(fd$dV[[i]]))
          expect_equal(an$dV[[i]][[k]], fd$dV[[i]][[k]], tolerance = 1e-6,
                       info = paste(nm, meth, "dV", i, k))
      }
    }
  }
})

test_that("the analytic Jacobian covers the covariate paths", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  # The quadrature and Taylor covariate grids change the node set the moments are
  # taken over, and a structural theta then reaches f through the shift nodes as
  # well as its own sensitivity column. Both moments are LINEAR in the raw
  # sensitivity column, which is why summing the columns before applying the
  # chain is sufficient -- and this is the test that says so.
  skip_if_not_installed("rxode2")
  TT <- c(2, 5, 9, 14); DOSE <- 100; NQ <- 7L; N <- 200L
  fn <- function() { ini({ tcl <- log(1); tv <- log(10); tcov <- 0.4
                           eta.cl ~ 0.09; add.err <- 0.3 })
    model({ cl <- exp(tcl + tcov * WT + eta.cl); v <- exp(tv)
            cp <- linCmt(); cp ~ add(add.err) }) }
  ui   <- suppressMessages(rxode2::rxode2(fn))
  ov   <- admixr2:::.admOutputVar(ui)
  sens <- admixr2:::.admLoadSensModel(ui); rx <- admixr2:::.admLoadModel(ui)
  E0 <- DOSE / 10 * exp(-0.1 * TT)
  st <- list(s1 = list(E = E0, V = diag((0.25 * E0)^2), n = N, times = TT,
                       ev = rxode2::et(amt = DOSE), cov = list(WT = 0),
                       cov_dist = list(WT = list(mu = 0, sd = 0.3))))
  for (ci in c("quadrature", "taylor")) {
    ctl <- adghControl(studies = st, grad = "none", n_nodes = NQ, print = 0L,
                       covMethod = "none", cov_integration = ci)
    pin <- admixr2:::.admDriverPinfo(ui, ctl)
    u   <- admixr2:::.admDriverUnits(st, ui, ov)
    u$studies <- admixr2:::.admCheckCovariates(ui, pin, u$studies)
    g   <- admixr2:::.adghNodeGrid(NQ, pin$n_eta)
    p   <- admixr2:::.admBuildOptVec(pin)$p0 + 0.03
    fd  <- admixr2:::.admMomentDeriv(p, pin, u$studies, rx, ov, g, 1L)
    an  <- admixr2:::.admMomentJac(p, pin, u$studies, sens, rx, ov, g, 1L)
    expect_false(is.null(an), info = ci)
    expect_equal(an$dE[[1L]], fd$dE[[1L]], tolerance = 1e-6, info = ci)
    for (k in seq_along(fd$dV[[1L]]))
      expect_equal(an$dV[[1L]][[k]], fd$dV[[1L]][[k]], tolerance = 1e-6,
                   info = paste(ci, k))
  }
})

test_that(".admMomentJac refuses rather than approximates what it cannot reach", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  skip_if_not_installed("rxode2")
  TT <- c(2, 5, 9); DOSE <- 100
  fn <- function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16; e <- 0.5 })
    model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt(); cp ~ add(e) }) }
  ui  <- suppressMessages(rxode2::rxode2(fn)); ov <- admixr2:::.admOutputVar(ui)
  invisible(admixr2:::.admLoadSensModel(ui))   # sens BEFORE sim; see above
  rx  <- admixr2:::.admLoadModel(ui)
  E0  <- DOSE / 10 * exp(-0.1 * TT)
  st  <- list(s1 = list(E = E0, V = diag((0.25 * E0)^2), n = 100L, times = TT,
                        ev = rxode2::et(amt = DOSE)))
  ctl <- adghControl(studies = st, grad = "none", n_nodes = 7L, print = 0L,
                     covMethod = "none")
  pin <- admixr2:::.admDriverPinfo(ui, ctl)
  u   <- admixr2:::.admDriverUnits(st, ui, ov)
  g   <- admixr2:::.adghNodeGrid(7L, pin$n_eta)
  p   <- admixr2:::.admBuildOptVec(pin)$p0
  # no sensitivity model -> no analytic route, and it must say so rather than
  # silently reaching for the plain model's predictions
  expect_null(admixr2:::.admMomentJac(p, pin, u$studies, NULL, rx, ov, g, 1L))
})

test_that("TBS conditional moments match the pushforward they describe", {
  # For a transform-both-sides endpoint the residual is normal on the
  # TRANSFORMED scale, so the conditional law given a node is a pushforward of a
  # normal through m^-1. Its moments are the same Gauss-Hermite quadrature the
  # NLL already runs for the mean and variance, carried to third and fourth
  # order -- which is what lets boxCox/yeoJohnson/logitNorm/probitNorm reach the
  # sandwich instead of degrading to covMethod = "r".
  skip_if_not_installed("rxode2")
  YJ <- admixr2:::.ADM_TBS_YJ
  set.seed(8); n <- 4e5L
  cases <- list(
    list("boxCox lam=0.4",  0.4, YJ[["boxCox"]],     0, 1, 3.0, 0.25),
    list("boxCox lam=0",    0.0, YJ[["boxCox"]],     0, 1, 3.0, 0.25),
    list("yeoJohnson",      0.6, YJ[["yeoJohnson"]], 0, 1, 1.5, 0.30),
    list("logitNorm",       1.0, YJ[["logit"]],      0, 1, 0.4, 0.35),
    list("probitNorm",      1.0, YJ[["probit"]],     0, 1, 0.4, 0.35))
  for (cs in cases) {
    lbl <- cs[[1L]]; lam <- cs[[2L]]; yj <- cs[[3L]]
    lo <- cs[[4L]]; hi <- cs[[5L]]; f <- cs[[6L]]; sdv <- cs[[7L]]
    z <- admixr2:::.admTBS(f, lam, yj, lo, hi) + sdv * stats::rnorm(n)
    y <- admixr2:::.admTBSi(z, lam, yj, lo, hi)
    y <- y[is.finite(y)]
    yc <- y - mean(y)
    q  <- admixr2:::.admTBSCondMom(f, sdv^2, 0, 1, lam, yj, lo, hi,
                                   FALSE, FALSE, 81L)
    # relative, because the four transforms span three orders of magnitude
    expect_equal(q$d,  mean(yc^2), tolerance = 0.02, info = paste(lbl, "var"))
    expect_equal(q$t3, mean(yc^3), tolerance = 0.05, info = paste(lbl, "m3"))
    expect_equal(q$q4, mean(yc^4), tolerance = 0.05, info = paste(lbl, "m4"))
  }
})

test_that(".admHeatSE reads Omega's vech block, and refuses a shape mismatch", {
  # plot.admFit standardises the observed-minus-predicted covariance. Omega's
  # vech block holds Var(V_ij) directly, so the SE is the square root of its
  # diagonal -- no closed form, no normality assumption. Everything is guarded
  # on shape: the plot re-derives its own moments, and a mismatch must fall back
  # rather than index into the wrong matrix.
  m <- 3L
  ij <- which(lower.tri(diag(m), diag = TRUE), arr.ind = TRUE)
  q  <- nrow(ij)                                   # 6
  Om <- diag(c(rep(1, m), (seq_len(q) / 10)^2))    # known vech variances
  sw <- list(Om = list(Om), study_names = "s1")
  se <- admixr2:::.admHeatSE(sw, "s1", m, "cov")
  expect_false(is.null(se))
  for (b in seq_len(q)) {
    expect_equal(se[ij[b, 1L], ij[b, 2L]], b / 10)
    expect_equal(se[ij[b, 2L], ij[b, 1L]], b / 10)  # symmetric
  }
  # no Omega, unknown study, or the wrong shape -> NULL, so the caller keeps
  # the closed form
  expect_null(admixr2:::.admHeatSE(NULL, "s1", m, "cov"))
  expect_null(admixr2:::.admHeatSE(list(study_names = "s1"), "s1", m, "cov"))
  expect_null(admixr2:::.admHeatSE(sw, "other", m, "cov"))
  expect_null(admixr2:::.admHeatSE(sw, "s1", m + 1L, "cov"))
  # a `var` study reports only the diagonal, so .admWeightSel kept only those
  # rows and the matrix is 2m -- the off-diagonals come back NA
  Omv <- diag(c(rep(1, m), c(0.04, 0.09, 0.16)))
  sev <- admixr2:::.admHeatSE(list(Om = list(Omv), study_names = "s1"),
                              "s1", m, "var")
  expect_equal(diag(sev), c(0.2, 0.3, 0.4))
  expect_true(all(is.na(sev[upper.tri(sev)])))
})

# ---- the model-source split in the meat -------------------------------------

test_that("skip and extra PARTITION the meat exactly", {
  # The failure mode the split invites: a block that is skipped but not covered
  # by an `extra` term -- or covered twice -- is silent, finite and plausible.
  # Asked for by review of HANDOFF-model-source-n-IMPLEMENTED.md.
  set.seed(4)
  p <- 3L
  mk <- function(m) matrix(stats::rnorm(p * m), p, m)
  psd <- function(m) { A <- matrix(stats::rnorm(m * m), m, m); crossprod(A) + diag(m) }
  H  <- psd(p)
  G  <- list(mk(4L), mk(5L))
  Om <- list(psd(4L), psd(5L))
  M  <- psd(p)                                   # one model source's whole term

  # every block is data: the classic sum, unchanged
  a <- admixr2:::.admSandwich(H, G, Om)
  expect_equal(a$J, G[[1L]] %*% Om[[1L]] %*% t(G[[1L]]) +
                    G[[2L]] %*% Om[[2L]] %*% t(G[[2L]]), tolerance = 1e-12)

  # block 1 is a MODEL source: it leaves the data sum and rejoins as `extra`
  b <- admixr2:::.admSandwich(H, G, Om, extra = list(M), skip = 1L)
  expect_equal(b$J, M + G[[2L]] %*% Om[[2L]] %*% t(G[[2L]]), tolerance = 1e-12)

  # ... so adding the data source increases the meat by EXACTLY its own term
  c1 <- admixr2:::.admSandwich(H, G[1L], Om[1L], extra = list(M), skip = 1L)
  expect_equal(c1$J, M, tolerance = 1e-12)
  expect_equal(b$J - c1$J, G[[2L]] %*% Om[[2L]] %*% t(G[[2L]]), tolerance = 1e-12)

  # a skipped block must not also be counted: skipping BOTH leaves only `extra`
  d <- admixr2:::.admSandwich(H, G, Om, extra = list(M), skip = c(1L, 2L))
  expect_equal(d$J, M, tolerance = 1e-12)
  # and the sandwich itself is H^-1 J H^-1 throughout
  Hi <- solve(H)
  expect_equal(b$cov, Hi %*% b$J %*% Hi, tolerance = 1e-12)
})

test_that(".admCondCheck names the parameters a flat direction is carried by", {
  # Var = H^-1 (meat) H^-1, so a near-singular H poisons the answer however
  # cleanly the meat was formed -- and it does not fail loudly: the inverse is
  # finite, large and plausible.
  nms <- c("a", "b", "c")
  expect_null(admixr2:::.admCondCheck(diag(3), nms))
  # a exactly-flat direction in (a, b): H is rank 2
  V <- cbind(c(1, -1, 0) / sqrt(2), c(1, 1, 0) / sqrt(2), c(0, 0, 1))
  H <- V %*% diag(c(1e-18, 2, 3)) %*% t(V)
  chk <- admixr2:::.admCondCheck(H, nms)
  expect_false(is.null(chk))
  expect_identical(chk$level, "undetermined")
  expect_setequal(chk$pars, c("a", "b"))
  expect_equal(chk$ndir, 1L)
  # ... and those rows come back NA rather than as a large finite number
  cv <- diag(3); dimnames(cv) <- list(nms, nms)
  bl <- admixr2:::.admCondBlank(cv, chk)
  expect_true(all(is.na(diag(bl)[c("a", "b")])))
  expect_false(is.na(bl["c", "c"]))
})

test_that("the yardstick fires only when OUR SE beats the source's own", {
  # A summary of a published model cannot carry more information than the
  # analyst who had every patient.
  C <- diag(c(0.08, 0.06)^2); dimnames(C) <- list(c("tcl", "tv"), c("tcl", "tv"))
  st <- list(s1 = list(.adm_src = list(id = "t1", cov = C)))
  mkcov <- function(se) { m <- diag(se^2)
    dimnames(m) <- list(c("tcl", "tv"), c("tcl", "tv")); m }
  expect_null(admixr2:::.admSrcYardstick(mkcov(c(0.08, 0.06)), st))   # equal
  expect_null(admixr2:::.admSrcYardstick(mkcov(c(0.09, 0.06)), st))   # wider
  y <- admixr2:::.admSrcYardstick(mkcov(c(0.068, 0.06)), st)          # adfo's case
  expect_equal(y$pars, "tcl")
  expect_equal(unname(y$ratio), 0.85, tolerance = 1e-6)
  # a few tenths of a percent is quadrature/MC noise, not a finding
  expect_null(admixr2:::.admSrcYardstick(mkcov(c(0.0797, 0.06)), st))
  # with SEVERAL sources ADM legitimately beats any one of them
  st2 <- c(st, list(s2 = list(.adm_src = list(id = "t2", cov = C))))
  expect_null(admixr2:::.admSrcYardstick(mkcov(c(0.05, 0.04)), st2))
})

test_that("a WEAKLY determined Hessian is flagged but NOT blanked", {
  # Two different failures with two different remedies. Below sqrt(eps) the
  # inverse is set by rounding error, so NA is the honest report. Between that
  # and 1e-4 the inverse is still meaningful and the sandwich stays VALID --
  # measured coverage at or above nominal, because conservative is not broken --
  # so the numbers are kept and the fit says the interval is uninformative.
  nms <- c("a", "b", "c")
  V <- cbind(c(1, -1, 0) / sqrt(2), c(1, 1, 0) / sqrt(2), c(0, 0, 1))
  weak <- V %*% diag(c(1e-6, 2, 3)) %*% t(V)
  chk <- admixr2:::.admCondCheck(weak, nms)
  expect_identical(chk$level, "weak")
  expect_setequal(chk$pars, c("a", "b"))
  cv <- diag(3); dimnames(cv) <- list(nms, nms)
  # kept, not blanked -- a wide valid interval beats no interval
  expect_equal(admixr2:::.admCondBlank(cv, chk), cv)
  # ... and a well-conditioned Hessian is neither
  expect_null(admixr2:::.admCondCheck(
    V %*% diag(c(0.5, 2, 3)) %*% t(V), nms))
})

test_that("an independent between-value block enters the weight exactly", {
  # cov_integration = "taylor" splits V into the WITHIN-covariate-value spread
  # the node ensemble carries and a derivative term no node carries. The weight
  # used to get only the first, so Om described a narrower law than the
  # objective scored and every sandwich SE came out too small -- while
  # .admMomentJac already carried the term, so G and Om described different
  # models. Verified against a brute-force simulation of the SAME law rather
  # than against a second derivation.
  set.seed(11)
  m <- 3L; nq <- 6L; N <- 400
  C  <- matrix(stats::rnorm(nq * m, 0, 0.6), nq, m)
  w  <- stats::runif(nq, 0.5, 1.5); w <- w / sum(w)
  C  <- sweep(C, 2L, colSums(w * C))
  Dv <- matrix(stats::runif(nq * m, 0.05, 0.25), nq, m)
  A  <- matrix(stats::rnorm(m * m, 0, 0.3), m, m)
  Vx <- crossprod(A) + diag(0.05, m)

  W  <- admixr2:::.admAdfWeightFast(C, w, Dv, N, Vx = Vx)
  # NULL is bit-identical, which is what keeps every non-taylor fit unchanged
  expect_identical(admixr2:::.admAdfWeightFast(C, w, Dv, N),
                   admixr2:::.admAdfWeightFast(C, w, Dv, N, NULL, NULL, NULL))

  R  <- 8000L
  ij <- which(lower.tri(diag(m), diag = TRUE), arr.ind = TRUE)
  Lx <- t(chol(Vx))
  tau <- matrix(0, R, m + nrow(ij))
  for (r in seq_len(R)) {
    q  <- sample.int(nq, N, replace = TRUE, prob = w)
    y  <- C[q, , drop = FALSE] +
      matrix(stats::rnorm(N * m), N, m) * sqrt(Dv[q, , drop = FALSE]) +
      matrix(stats::rnorm(N * m), N, m) %*% t(Lx)
    yb <- colMeans(y); V <- crossprod(sweep(y, 2L, yb)) / N
    tau[r, ] <- c(yb, V[cbind(ij[, 1L], ij[, 2L])])
  }
  Ob <- stats::cov(tau)
  expect_equal(max(abs(W - Ob)) / max(abs(Ob)), 0, tolerance = 0.06)
  # and it is strictly LARGER than the version that omitted it -- the omission
  # was in the over-confident direction
  Wn <- admixr2:::.admAdfWeightFast(C, w, Dv, N)
  expect_true(all(diag(W) >= diag(Wn) - 1e-12))
  expect_lt(min(sqrt(diag(Wn)) / sqrt(diag(W))), 0.9)
})
