# The loading, and the contract that makes it safe.
#
# The collapse rests on ONE claim: each covariate-reading assignment depends on
# the latent normal only through a single linear combination, p = G(b'xi). The
# loading is then the relative gradient d log p / d xi, whose DIRECTION is b at
# every xi. This file is that claim, checked against CLOSED FORM rather than
# against another of the package's own routes -- an agreement test between two
# internal implementations passes just as happily when both are wrong.
#
# Replaces test-collapse-syntactic.R. That file tested a second, mu-referenced
# route which no longer exists: the gradient IS the arithmetic that route
# computed symbolically, so a model now collapses on what it DOES rather than on
# how it was spelled. Its three traps are the three tests below.

.cl_mk <- function(clline, pars) {
  f <- function() {}
  body(f) <- as.call(c(list(quote(`{`), as.call(list(quote(ini), pars)),
    as.call(list(quote(model), as.call(c(list(quote(`{`), clline),
      as.list(quote({ v <- exp(tv)
                      d/dt(central) <- -(cl / v) * central
                      cp <- central / v; cp ~ prop(prop.sd) }))[-1L])))))))
  f
}
.cl_pars  <- quote({ tcl <- log(5); tv <- log(10); b1 <- 0.75
                     prop.sd <- c(0, 0.2); eta.cl ~ 0.09 })
.cl_pars2 <- quote({ tcl <- log(5); tv <- log(10); b1 <- 0.75; b2 <- 0.4
                     prop.sd <- c(0, 0.2); eta.cl ~ 0.09 })
.cl_small <- quote({ tcl <- log(5); tv <- log(10); b1 <- 0.004
                     prop.sd <- c(0, 0.2); eta.cl ~ 0.09 })

# Returns list(jc, st, L, B) for one model, or NULL if the collapse declined.
.cl_case <- function(f, cd, om = matrix(0.09, 1, 1)) {
  ui    <- suppressMessages(rxode2::rxode2(f))
  pinfo <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  pinfo$cov_nodes <- 7L
  s  <- list(times = c(1, 4, 12), ev = rxode2::et(amt = 100), n = 200L,
             cov_dist = cd)
  jc <- admixr2:::.admJointCollapse(ui, pinfo, cd, 7L, s,
                                    admixr2:::.admOutputVar(ui))
  if (is.null(jc)) return(NULL)
  st <- admixr2:::.admShiftStruct(
          pinfo, admixr2:::.admUnpack(
                   admixr2:::.admBuildOptVec(pinfo)$p0, pinfo)$struct)
  L <- t(chol(as.matrix(om)))
  list(jc = jc, st = st, L = L,
       B = admixr2:::.admJointB(jc, st, L, jc$Xi))
}
.cl_lnWT  <- function() covDist(WT = c(mean = 78, sd = 78 * 0.30), dist = "lnorm")
.cl_nrmWT <- function() covDist(WT = c(mean = 78, sd = 8), dist = "normal")
.cl_sdlog <- function(m, cv) sqrt(log(1 + cv^2))

test_that("the loading equals the derivative in CLOSED FORM", {
  skip_if_not_installed("rxode2")
  # cl = exp(tcl + eta.cl + b1*log(WT/70)) with WT lognormal, so
  #   log cl = tcl + eta.cl + b1*(meanlog + sdlog*z - log 70)
  # giving d log cl/d(eta axis) = L[1,1] and d log cl/dz = b1 * sdlog, both
  # exactly. Nothing here is a tolerance on a proxy.
  r <- .cl_case(.cl_mk(quote(cl <- exp(tcl + eta.cl + b1 * log(WT / 70))),
                       .cl_pars), .cl_lnWT())
  expect_false(is.null(r))
  expect_false(is.null(r$B))
  expect_equal(dim(r$B), c(r$jc$nl, 1L))
  expect_equal(r$B[1L, 1L], r$L[1L, 1L], tolerance = 1e-6)
  expect_equal(r$B[r$jc$ne + 1L, 1L],
               0.75 * .cl_sdlog(78, 0.30) * r$jc$Lc[1L, 1L], tolerance = 1e-6)
})

test_that("correlated covariates load through the right ROW of Lc", {
  skip_if_not_installed("rxode2")
  # Zc = xi %*% Lc, so d z / d xi_k is ROW k, not column k. The two coincide at
  # one covariate, so a correlated PAIR is the only thing that can catch the
  # swap -- it was wrong by ~33% and every single-covariate test still passed.
  cd <- covDist(WT = c(mean = 78, sd = 78 * 0.30),
                CRCL = c(mean = 90, sd = 90 * 0.35), dist = "lnorm",
                cor = c(WT.CRCL = 0.45))
  r <- .cl_case(.cl_mk(quote(cl <- exp(tcl + eta.cl + b1 * log(WT / 70) +
                                         b2 * log(CRCL / 90))), .cl_pars2), cd)
  expect_false(is.null(r))
  expect_false(is.null(r$B))
  b   <- c(0.75 * .cl_sdlog(78, 0.30), 0.4 * .cl_sdlog(90, 0.35))
  got <- r$B[r$jc$ne + seq_len(r$jc$pc), 1L]
  expect_equal(as.numeric(got), as.numeric(r$jc$Lc %*% b), tolerance = 1e-6)
  # and the transpose is genuinely different, so the test can fail
  expect_false(isTRUE(all.equal(as.numeric(r$jc$Lc %*% b),
                                as.numeric(t(r$jc$Lc) %*% b), tolerance = 1e-3)))
})

test_that("the certificate admits a NORMAL margin and refuses a LOGNORMAL one", {
  skip_if_not_installed("rxode2")
  # `b1*WT` depends on xi through one linear combination when WT is NORMAL
  # (d/dxi = b1*sd, constant) and does NOT when WT is lognormal
  # (d/dxi = b1*WT*sdlog, which moves with z). BOTH halves matter: without the
  # first, a certificate that declined everything would look like a pass.
  ok <- .cl_case(.cl_mk(quote(cl <- exp(tcl + eta.cl + b1 * WT)), .cl_small),
                 .cl_nrmWT())
  expect_false(is.null(ok))
  expect_false(is.null(ok$B))
  expect_equal(ok$B[ok$jc$ne + 1L, 1L],
               0.004 * 8 * ok$jc$Lc[1L, 1L], tolerance = 1e-5)

  no <- .cl_case(.cl_mk(quote(cl <- exp(tcl + eta.cl + b1 * WT)), .cl_small),
                 .cl_lnWT())
  expect_true(is.null(no) || is.null(no$B))
})

test_that("the loading does not depend on HOW the model was SPELLED", {
  skip_if_not_installed("rxode2")
  # exp(tcl+eta.cl)*(WT/70)^b1 and exp(tcl+eta.cl+b1*log(WT/70)) are the same
  # model. One is mu-referenced and the other is not, and that used to decide
  # which route computed the loading. It no longer decides anything.
  a <- .cl_case(.cl_mk(quote(cl <- exp(tcl + eta.cl + b1 * log(WT / 70))),
                       .cl_pars), .cl_lnWT())
  b <- .cl_case(.cl_mk(quote(cl <- exp(tcl + eta.cl) * (WT / 70)^b1),
                       .cl_pars), .cl_lnWT())
  expect_false(is.null(a$B)); expect_false(is.null(b$B))
  expect_equal(as.numeric(a$B), as.numeric(b$B), tolerance = 1e-6)
})

test_that("a SINGLE-INDEX link survives the random effect", {
  skip_if_not_installed("randtoolbox")
  # .admCovCollapse re-probes at eta = 0.5 and requires the SAME loadings,
  # because a covariate-by-eta interaction would otherwise collapse onto the
  # wrong subspace silently. The RAW slope moves with eta -- it scales by
  # exp(eta), 6.5e-01 at 0.5 -- so every model whose covariates entered through
  # a nonlinear LINK failed that check and fell back to the product grid, even
  # though its INDEX was perfectly affine. The relative gradient drops eta
  # exactly, since log p = theta + eta + log G(b'z). See
  # algorithm/collapse-derivation/relgrad_eta.txt.
  cd  <- covDist(WT = c(mean = 78, sd = 78 * .3),
                 CRCL = c(mean = 90, sd = 90 * .35),
                 AGE = c(mean = 50, sd = 50 * .25), dist = "lnorm")
  pin <- list(eta_col_names = "eta.cl", struct_names = character(0),
              cov_nodes = 7L)
  rk  <- function(expr) {
    co <- admixr2:::.admCovCollapse(
            list(lstExpr = expr, allCovs = c("WT", "CRCL", "AGE")), pin, cd, 7L)
    if (is.null(co)) NA_integer_ else co$r
  }
  idx <- quote(s <- 0.6 * log(WT / 70) + 0.4 * log(CRCL / 90) +
                    0.3 * log(AGE / 50))
  # three covariates through one index, an eta on the same parameter: rank 1
  expect_equal(rk(list(idx, quote(cl <- exp(0.5 + eta.cl) * sqrt(1.5 + s)))), 1L)
  expect_equal(rk(list(idx, quote(cl <- exp(0.5 + eta.cl) * (1 + s) / (2 + s)))),
               1L)
  # a GENUINE covariate-by-eta interaction must still be refused, or the check
  # above only proves the certificate stopped looking
  expect_true(is.na(rk(list(quote(
    cl <- exp(0.5 + eta.cl + 0.6 * eta.cl * log(WT / 70) +
              0.4 * log(CRCL / 90) + 0.3 * log(AGE / 50)))))))
})
