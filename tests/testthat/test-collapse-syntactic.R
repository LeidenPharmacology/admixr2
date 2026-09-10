# The syntactic loading route, and the contract that makes it safe.
#
# Compared as NUMBERS throughout: the two routes carry different `routes`
# attributes by design -- the syntactic one has none to replay -- so stripping
# the class is not enough to compare them.
#
# .admJointB() tries mu-referencing first and falls back to .admCovLoadings().
# That is only sound if the two AGREE wherever both apply -- otherwise which
# answer a fit gets depends on how the model was typed. Everything here is that
# contract. Each case below corresponds to a defect found while building it.

.syn_mk <- function(clline, pars) {
  f <- function() {}
  body(f) <- as.call(c(list(quote(`{`), as.call(list(quote(ini), pars)),
    as.call(list(quote(model), as.call(c(list(quote(`{`), clline),
      as.list(quote({ v <- exp(tv)
                      d/dt(central) <- -(cl / v) * central
                      cp <- central / v; cp ~ prop(prop.sd) }))[-1L])))))))
  f
}
.syn_pars  <- quote({ tcl <- log(5); tv <- log(10); b1 <- 0.75
                      prop.sd <- c(0, 0.2); eta.cl ~ 0.09 })
.syn_pars2 <- quote({ tcl <- log(5); tv <- log(10); b1 <- 0.75; b2 <- 0.4
                      prop.sd <- c(0, 0.2); eta.cl ~ 0.09 })
.syn_small <- quote({ tcl <- log(5); tv <- log(10); b1 <- 0.004
                      prop.sd <- c(0, 0.2); eta.cl ~ 0.09 })

# Returns list(syn, num, used) for one model, or NULL if the collapse declined.
.syn_case <- function(f, cd, om) {
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
  L   <- t(chol(as.matrix(om)))
  jc0 <- jc; jc0$syn <- NULL                       # force the numeric route
  list(syn  = admixr2:::.admSynB(jc, st, L),
       num  = admixr2:::.admJointB(jc0, st, L, jc$Xi, NULL),
       used = admixr2:::.admJointB(jc, st, L, jc$Xi, NULL))
}

.syn_lnWT  <- function() covDist(WT = c(mean = 78, sd = 78 * 0.30), dist = "lnorm")
.syn_nrmWT <- function() covDist(WT = c(mean = 78, sd = 8), dist = "normal")

test_that("a mu-referenced covariate gives the same loadings as probing does", {
  skip_if_not_installed("rxode2")
  r <- .syn_case(.syn_mk(quote(cl <- exp(tcl + eta.cl + b1 * log(WT / 70))),
                         .syn_pars), .syn_lnWT(), matrix(0.09, 1, 1))
  expect_false(is.null(r$syn))
  expect_equal(as.numeric(r$syn), as.numeric(r$num), tolerance = 1e-8)
  # ... and it is the one the fit actually uses
  expect_equal(as.numeric(r$used), as.numeric(r$syn), tolerance = 1e-14)
})

test_that("correlated covariates load through the right row of Lc", {
  skip_if_not_installed("rxode2")
  # Zc = xi %*% Lc, so d z / d xi_k is ROW k, not column k. The two coincide at
  # one covariate, so this pair is the only thing that can catch the swap -- it
  # was wrong by ~33% and every single-covariate test still passed.
  cd <- covDist(WT = c(mean = 78, sd = 78 * 0.30),
                CRCL = c(mean = 90, sd = 90 * 0.35), dist = "lnorm",
                cor = c(WT.CRCL = 0.45))
  r <- .syn_case(.syn_mk(quote(cl <- exp(tcl + eta.cl + b1 * log(WT / 70) +
                                           b2 * log(CRCL / 90))), .syn_pars2),
                 cd, matrix(0.09, 1, 1))
  expect_false(is.null(r$syn))
  expect_equal(as.numeric(r$syn), as.numeric(r$num), tolerance = 1e-8)
})

test_that("the affine certificate admits a NORMAL margin and refuses a LOGNORMAL one", {
  skip_if_not_installed("rxode2")
  # `b1*WT` is affine in the latent when WT is NORMAL (dT/dxi = b1*sd) and is
  # NOT when WT is lognormal (dT/dxi = b1*WT*sdlog, which moves with z). Both
  # halves matter: without the first, a certificate that declined everything
  # would look like a pass.
  ok <- .syn_case(.syn_mk(quote(cl <- exp(tcl + eta.cl + b1 * WT)), .syn_small),
                  .syn_nrmWT(), matrix(0.09, 1, 1))
  expect_false(is.null(ok$syn))
  expect_equal(as.numeric(ok$syn), as.numeric(ok$num), tolerance = 1e-8)

  no <- .syn_case(.syn_mk(quote(cl <- exp(tcl + eta.cl + b1 * WT)), .syn_small),
                  .syn_lnWT(), matrix(0.09, 1, 1))
  expect_null(no$syn)
  expect_equal(as.numeric(no$used), as.numeric(no$num), tolerance = 1e-12)
})

test_that("a model that is not mu-referenced falls back, it does not fail", {
  skip_if_not_installed("rxode2")
  # `(WT/70)^b1` outside the exp is caught by no mu-ref frame at any edition,
  # and rxode2 is right -- that model is not mu-referenced. The numeric route
  # must still produce the loadings, unchanged.
  r <- .syn_case(.syn_mk(quote(cl <- exp(tcl + eta.cl) * (WT / 70)^b1),
                         .syn_pars), .syn_lnWT(), matrix(0.09, 1, 1))
  expect_null(r$syn)
  expect_false(is.null(r$num))
  expect_equal(as.numeric(r$used), as.numeric(r$num), tolerance = 1e-12)
})
