# The ADF weight: the sampling law of (ybar, vech V) implied by the model.

test_that(".admSandwichGrid never exceeds max_nodes, even for many etas", {
  # The decrement loop used to floor at nq = 3 and stop, so it never re-checked
  # whether 3^n_eta was still over the cap -- for n_eta >= 8, 3^n_eta > 5000 and
  # the "capped" grid silently blew past it (3^8 = 6561, 3^15 ~= 14.3M).
  for (n_eta in c(1L, 5L, 8L, 12L)) {
    grid <- admixr2:::.admSandwichGrid(list(n_eta = n_eta), max_nodes = 5000L)
    expect_lte(nrow(grid$X), 5000L)
    expect_length(grid$W, nrow(grid$X))
  }
})

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
  # numDeriv is Suggests, so this has to degrade rather than error where it is
  # absent -- the companion TBS test below already guarded it and this one did not.
  skip_if_not_installed("numDeriv")
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
    p2  <- admixr2:::.admAdfParts(pars, pin, s1, rx, ov, g, 1L)
    W   <- admixr2:::.admWorkingWeight(p2$V, N, s1$method)
    # BOTH routes to the moment Jacobian must satisfy it -- the identity is what
    # says the analytic one describes the same objective the FD one does.
    sens <- admixr2:::.admLoadSensModel(ui)
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
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt(); cp ~ lnorm(e) }) },
    # TBS belongs in this test, not in one of its own. Once the objective
    # composes at the nodes the weight and the objective are the SAME
    # aggregation, so the identity holds to machine precision here too -- it was
    # only the delta expansion that made these families a special case.
    boxCox = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16
                                e <- 0.3; lam <- fix(0.5) })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt()
              cp ~ add(e) + boxCox(lam) }) },
    bc_prop = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16
                                 e <- 0.3; b <- 0.2; lam <- fix(0.5) })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt()
              cp ~ add(e) + prop(b) + boxCox(lam) }) },
    yeoJohnson = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16
                                    e <- 0.3; lam <- fix(0.5) })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt()
              cp ~ add(e) + yeoJohnson(lam) }) },
    logitNorm = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16
                                   e <- 0.25 })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt()
              cp ~ logitNorm(e, 0, 40) }) },
    probitNorm = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16
                                    e <- 0.25 })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt()
              cp ~ probitNorm(e, 0, 40) }) })
  for (nm in names(mods)) {
    ui <- suppressMessages(rxode2::rxode2(mods[[nm]]))
    ov <- admixr2:::.admOutputVar(ui); rx <- admixr2:::.admLoadModel(ui)
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

test_that("the analytic moment Jacobian matches the finite-difference oracle", {
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
    # The TBS family: its forward composition goes through dmu_dv0/dms_df, the
    # mean-from-covariance path no closed-form endpoint exercises, so the oracle
    # matters more here than anywhere else in this list.
    boxCox  = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16
                                 e <- 0.3; lam <- fix(0.5) })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt()
              cp ~ add(e) + boxCox(lam) }) },
    logitN  = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16; e <- 0.25 })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt()
              cp ~ logitNorm(e, 0, 40) }) },
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

test_that(".admMomentJac refuses rather than approximates what it cannot reach", {
  skip_if_not_installed("rxode2")
  TT <- c(2, 5, 9); DOSE <- 100
  fn <- function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16; e <- 0.5 })
    model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt(); cp ~ add(e) }) }
  ui  <- suppressMessages(rxode2::rxode2(fn)); ov <- admixr2:::.admOutputVar(ui)
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

test_that("TBS conditional moments match a simulation of the same law", {
  # boxCox / yeoJohnson / logitNorm / probitNorm are the one residual family the
  # sandwich used to refuse. Nothing about them breaks the Isserlis expansion --
  # they are conditionally independent across timepoints like every other family
  # here -- they simply had no closed form for the third and fourth moments. The
  # quadrature that already computes the mean and variance supplies them.
  #
  # The oracle is a direct simulation of y | f = g(h(f) + sd * eps), NOT another
  # admixr2 path, and the moments are compared on a scale-free footing (each
  # divided by the matching power of the SD) so that a mu3 which is truly zero is
  # judged absolutely rather than against Monte Carlo noise.
  skip_on_cran()
  set.seed(303)
  n <- 4e5
  cases <- list(
    list("boxCox l=0.5",     0.5, 0L, c(2, 5, 10),      0.30),
    list("boxCox l=0",       0.0, 0L, c(2, 5, 10),      0.25),
    list("yeoJohnson l=0.5", 0.5, 1L, c(-2, 2, 8),      0.30),
    list("logitNorm",        1.0, 4L, c(0.2, 0.5, 0.8), 0.40),
    list("probitNorm",       1.0, 6L, c(0.2, 0.5, 0.8), 0.40))
  for (cs in cases) {
    lbl <- cs[[1L]]; lam <- cs[[2L]]; yj <- cs[[3L]]; f <- cs[[4L]]; sd <- cs[[5L]]
    tm <- admixr2:::.admTBSCentral(f, rep(sd, length(f)), lam, yj, 0, 1)
    for (i in seq_along(f)) {
      z <- admixr2:::.admTBS(f[i], lam, yj, 0, 1) + sd * stats::rnorm(n)
      y <- admixr2:::.admTBSi(z, lam, yj, 0, 1)
      y <- y[is.finite(y)]
      yc <- y - mean(y)
      mc <- c(mean(yc^2), mean(yc^3), mean(yc^4))
      qd <- c(tm$v[i], tm$mu3[i], tm$mu4[i])
      sc <- sqrt(mc[1L])^c(2, 3, 4)
      expect_lt(max(abs(qd - mc) / sc), 0.05)
    }
  }
  # An untransformed endpoint is exactly normal, so this is a closed-form check
  # rather than a simulated one: v = sd^2, mu3 = 0, mu4 = 3 sd^4.
  ex <- admixr2:::.admTBSCentral(c(3, 7), c(0.25, 0.25), 1, 2L, 0, 1)
  expect_equal(ex$v,   rep(0.0625, 2))
  expect_equal(ex$mu3, rep(0, 2), tolerance = 1e-10)
  expect_equal(ex$mu4, rep(3 * 0.0625^2, 2))
})

test_that("the TBS weight reproduces the sampling law of a simulated study", {
  # The same check da4fddb ran for lognormal and Poisson, on the family this
  # commit adds. `Omega` is scored against the empirical covariance of
  # (ybar, vech V) over simulated studies -- the definition of what it is meant
  # to be -- and the third/fourth moments are shown to be load-bearing by
  # dropping them: on boxCox at lambda = 0 (a log transform, so a strongly
  # skewed residual) the cross block goes from ~0.04 to ~0.85.
  skip_on_cran()
  set.seed(909)
  TIMES <- c(1, 3, 6); DOSE <- 100; V0 <- 10; OM <- 0.4; N <- 80L; R <- 3000L
  gh <- admixr2:::.adghNodes1(31L); x <- gh$x; w <- gh$w / sum(gh$w)
  fmat <- function(eta) t(vapply(eta, function(e)
    DOSE / V0 * exp(-exp(e) / V0 * TIMES), numeric(length(TIMES))))
  F <- fmat(x * OM); m <- ncol(F)
  ij <- which(lower.tri(diag(m), diag = TRUE), arr.ind = TRUE)
  lam <- 0.0; yj <- 0L; sd <- 0.25                       # boxCox, lambda 0

  tm <- admixr2:::.admTBSCentral(as.numeric(F), rep(sd, length(F)), lam, yj, 0, 1)
  M1 <- matrix(tm$m, nrow(F)); Dv <- matrix(tm$v, nrow(F))
  T3 <- matrix(tm$mu3, nrow(F)); Q4 <- matrix(tm$mu4, nrow(F))
  # The conditional mean is NONLINEAR in f here, so C is the centred conditional
  # means themselves -- the linearisation every closed-form family uses is only
  # exact because their means are linear in f.
  C <- sweep(M1, 2L, as.numeric(crossprod(w, M1)))

  tv <- matrix(0, R, m + nrow(ij))
  for (r in seq_len(R)) {
    f <- fmat(stats::rnorm(N, 0, OM))
    z <- admixr2:::.admTBS(f, lam, yj, 0, 1) + sd * stats::rnorm(length(f))
    y <- matrix(admixr2:::.admTBSi(z, lam, yj, 0, 1), nrow(f))
    V <- stats::cov.wt(y, method = "ML")$cov
    tv[r, ] <- c(colMeans(y), V[cbind(ij[, 1L], ij[, 2L])])
  }
  ref <- stats::cov(tv)
  rel <- function(A, r, c)
    sqrt(sum((A[r, c] - ref[r, c])^2)) / sqrt(sum(ref[r, c]^2))
  gen <- admixr2:::.admAdfWeightFast(C, w, Dv, N, T3, Q4)
  sym <- admixr2:::.admAdfWeightFast(C, w, Dv, N)
  expect_lt(rel(gen, seq_len(m), seq_len(m)),   0.15)
  expect_lt(rel(gen, seq_len(m), -seq_len(m)),  0.20)
  expect_lt(rel(gen, -seq_len(m), -seq_len(m)), 0.20)
  # dropping the moments wrecks the cross block and leaves the mean block intact,
  # which is why this class of error cannot be seen from the mean side
  expect_lt(rel(sym, seq_len(m), seq_len(m)),   0.15)
  expect_gt(rel(sym, seq_len(m), -seq_len(m)),  0.50)
})

test_that("the information equality holds on a TBS endpoint too", {
  # Separate from the equality test above rather than folded into it: that one
  # caught three real bugs and is left structurally alone. This asks the same
  # question of the family added later, where the residual composition runs
  # through the mean-from-covariance path (dmu_dv0, dms_df) that no closed-form
  # endpoint exercises.
  #
  # It also settles what the S != V_pred gap on TBS is and is not. G here comes
  # from the OBJECTIVE's moment map and H from the objective's Hessian, so this
  # holding exactly says the two halves the sandwich shares are consistent for
  # TBS. The gap lives entirely in the objective's own second-order expansion,
  # not in the sandwich's reading of it.
  skip_if_not_installed("rxode2")
  skip_if_not_installed("numDeriv")
  TIMES <- c(4, 8, 12, 16); DOSE <- 100; NQ <- 9L; N <- 100L
  mods <- list(
    boxCox = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16
                                aa <- 0.3; lam <- fix(0.5) })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt()
              cp ~ add(aa) + boxCox(lam) }) },
    logitNorm = function() { ini({ tcl <- log(1); tv <- log(10); eta.cl ~ 0.16
                                   aa <- 0.25 })
      model({ cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt()
              cp ~ logitNorm(aa, 0, 40) }) })
  for (nm in names(mods)) {
    ui <- suppressMessages(rxode2::rxode2(mods[[nm]]))
    ov <- admixr2:::.admOutputVar(ui); rx <- admixr2:::.admLoadModel(ui)
    sens <- admixr2:::.admLoadSensModel(ui)
    E0 <- DOSE / 10 * exp(-0.1 * TIMES)
    mk <- function(E, V, meth)
      list(s = list(E = E, V = if (meth == "var") diag(V) else V, n = N,
                    times = TIMES, ev = rxode2::et(amt = DOSE)))
    for (meth in c("var", "cov")) {
      st  <- mk(E0, diag((0.3 * E0)^2), meth)
      ctl <- adghControl(studies = st, grad = "none", n_nodes = NQ, print = 0L,
                         covMethod = "none")
      pin <- admixr2:::.admDriverPinfo(ui, ctl)
      u   <- admixr2:::.admDriverUnits(st, ui, ov)
      g   <- admixr2:::.adghNodeGrid(NQ, pin$n_eta)
      p   <- admixr2:::.admBuildOptVec(pin)$p0
      pars <- admixr2:::.admUnpack(p, pin)
      pt  <- admixr2:::.admAdfParts(pars, pin, u$studies[[1L]], rx, ov, g, 1L)
      # well specified: the summary IS what the model predicts, so t = tau
      u2  <- admixr2:::.admDriverUnits(mk(as.numeric(pt$E), pt$V, meth), ui, ov)
      s1  <- u2$studies[[1L]]
      f   <- function(q) admixr2:::.adghNLL(q, pin, u2$studies, rx, ov, g, 1L)
      H   <- numDeriv::hessian(f, p)
      p2  <- admixr2:::.admAdfParts(pars, pin, s1, rx, ov, g, 1L)
      W   <- admixr2:::.admWorkingWeight(p2$V, N, s1$method)
      mds <- list(
        fd       = admixr2:::.admMomentDeriv(p, pin, u2$studies, rx, ov, g, 1L),
        analytic = admixr2:::.admMomentJac(p, pin, u2$studies, sens, rx, ov, g, 1L))
      for (route in names(mds)) {
        md <- mds[[route]]
        expect_false(is.null(md), info = paste(nm, meth, route))
        G  <- admixr2:::.admScoreCross(md$E[[1L]], md$V[[1L]], md$dE[[1L]],
                                       md$dV[[1L]], s1, N)
        ev <- sort(Re(eigen(solve(2 * H) %*% (G %*% W %*% t(G)))$values))
        expect_equal(ev, rep(1, length(ev)), tolerance = 1e-4,
                     info = paste(nm, meth, route))
      }
    }
  }
})
