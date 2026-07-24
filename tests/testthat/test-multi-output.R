# Tier-1 unit tests for the multi-compartment (multiple observed output) data
# model: observation-list normalisation, flattening to units, sigma->output
# selection, and multi-output detection. No rxode2 required.

# ---- .admNormaliseObs / .admNormaliseStudy: observations list ----------------

test_that("observations list expands into per-compartment units", {
  study <- list(
    observations = list(
      plasma = list(output = "cp",   times = c(1, 2, 4), E = c(3, 2, 1),
                    V = c(0.1, 0.1, 0.1), n = 40L, ev = "EVp"),
      csf    = list(output = "cCSF", times = c(2, 4),    E = c(0.5, 0.4),
                    V = c(0.02, 0.02),    n = 12L, ev = "EVc")
    )
  )
  ns <- admixr2:::.admNormaliseStudy(study, "rat", default_output = "cp")

  expect_true(ns$multi)
  expect_named(ns$observations, c("plasma", "csf"))
  expect_equal(ns$observations$plasma$output, "cp")
  expect_equal(ns$observations$csf$output,    "cCSF")
  expect_equal(ns$observations$plasma$label,  "rat.plasma")
  expect_equal(ns$observations$csf$label,     "rat.csf")
  # vector V -> diag + method 'var'
  expect_equal(ns$observations$plasma$method, "var")
  expect_equal(ns$observations$plasma$v_diag, c(0.1, 0.1, 0.1))
  expect_equal(ns$observations$csf$n, 12L)
})

test_that("observations inherit study-level n / ev / output when omitted", {
  study <- list(
    n = 25L, ev = "sharedEV", output = "cp",
    observations = list(
      a = list(times = c(1, 2), E = c(1, 2), V = c(0.1, 0.1)),
      b = list(times = 3,       E = 4,       V = 0.2, output = "cCSF")
    )
  )
  ns <- admixr2:::.admNormaliseStudy(study, "s", default_output = "cp")
  expect_equal(ns$observations$a$n, 25L)
  expect_equal(ns$observations$a$ev, "sharedEV")
  expect_equal(ns$observations$a$output, "cp")     # from study-level output
  expect_equal(ns$observations$b$output, "cCSF")   # per-obs override
  expect_equal(ns$observations$b$ev, "sharedEV")
})

test_that("unnamed observations get obs<i> names", {
  study <- list(observations = list(
    list(output = "cp",   times = 1, E = 1, V = 0.1, n = 5L),
    list(output = "cCSF", times = 1, E = 1, V = 0.1, n = 5L)
  ))
  ns <- admixr2:::.admNormaliseStudy(study, "s")
  expect_named(ns$observations, c("obs1", "obs2"))
})

test_that("dimension mismatches error with the observation label", {
  bad_E <- list(observations = list(
    plasma = list(output = "cp", times = c(1, 2, 3), E = c(1, 2),
                  V = diag(2), n = 5L)))
  expect_error(admixr2:::.admNormaliseStudy(bad_E, "rat"),
               regexp = "rat.plasma.*length\\(E\\)")

  bad_V <- list(observations = list(
    csf = list(output = "cCSF", times = c(1, 2), E = c(1, 2),
               V = matrix(0, 3, 3), n = 5L)))
  expect_error(admixr2:::.admNormaliseStudy(bad_V, "rat"),
               regexp = "rat.csf.*V must be")
})

test_that("empty observations list errors", {
  expect_error(
    admixr2:::.admNormaliseStudy(list(observations = list()), "s"),
    regexp = "non-empty")
})

# ---- legacy single-output studies still normalise as before ------------------

test_that("legacy flat study keeps top-level fields and gains one unit", {
  s  <- list(E = c(1, 2), V = diag(c(0.1, 0.2)), n = 50L, times = c(1, 2))
  ns <- admixr2:::.admNormaliseStudy(s, "leg", default_output = "cp")
  expect_false(ns$multi)
  expect_equal(ns$method, "var")             # top-level preserved
  expect_equal(ns$v_diag, c(0.1, 0.2))
  expect_length(ns$observations, 1L)
  expect_equal(ns$observations[[1]]$output, "cp")
  expect_equal(ns$observations[[1]]$label, "leg")
})

# ---- .admFlattenStudies ------------------------------------------------------

test_that("flatten produces one unit per observation across studies", {
  studies <- list(
    ratA = admixr2:::.admNormaliseStudy(list(observations = list(
      plasma = list(output = "cp",   times = 1, E = 1, V = 0.1, n = 5L),
      csf    = list(output = "cCSF", times = 1, E = 1, V = 0.1, n = 5L))), "ratA"),
    ratB = admixr2:::.admNormaliseStudy(
      list(E = 1, V = 0.1, n = 8L, times = 1), "ratB", default_output = "cp")
  )
  units <- admixr2:::.admFlattenStudies(studies)
  expect_length(units, 3L)
  expect_equal(names(units), c("ratA.plasma", "ratA.csf", "ratB"))
  expect_equal(vapply(units, function(u) u$output, character(1)),
               c(ratA.plasma = "cp", ratA.csf = "cCSF", ratB = "cp"))
})

# ---- .admIsMultiOutput -------------------------------------------------------

test_that("multi-output detection keys on distinct output vars", {
  one <- list(list(output = "cp"), list(output = "cp"))
  two <- list(list(output = "cp"), list(output = "cCSF"))
  expect_false(admixr2:::.admIsMultiOutput(one, "cp"))
  expect_true(admixr2:::.admIsMultiOutput(two, "cp"))
})

# ---- .admSigmaSel ------------------------------------------------------------

test_that("sigma selector maps sigmas to their output, all-TRUE when unknown", {
  # Known mapping
  pinfo <- list(sigma_names  = c("prop.cp", "add.cCSF"),
                sigma_output = c("cp", "cCSF"))
  expect_equal(admixr2:::.admSigmaSel(pinfo, "cp"),   c(TRUE,  FALSE))
  expect_equal(admixr2:::.admSigmaSel(pinfo, "cCSF"), c(FALSE, TRUE))

  # Unknown mapping (NA) -> every sigma belongs to the single output
  pinfo_na <- list(sigma_names  = c("a", "b"),
                   sigma_output = c(NA_character_, NA_character_))
  expect_equal(admixr2:::.admSigmaSel(pinfo_na, "cp"), c(TRUE, TRUE))

  # No sigmas
  pinfo0 <- list(sigma_names = character(0), sigma_output = character(0))
  expect_equal(admixr2:::.admSigmaSel(pinfo0, "cp"), logical(0))
})

test_that("sigma selector maps linCmt endpoint names to the ipredSim column", {
  # sigma_output holds the endpoint name (iniDf condition = "rxLinCmt") while the
  # rxSolve output column is "ipredSim"; the selector must still match so linCmt
  # residual error is not silently dropped.
  pinfo <- list(sigma_names = "prop.err", sigma_output = "rxLinCmt")
  expect_equal(admixr2:::.admSigmaSel(pinfo, "ipredSim"), TRUE)
})

# ---- .admOutputVars ----------------------------------------------------------

test_that("output vars lists every predDf endpoint; linCmt maps to ipredSim", {
  ui_two <- list(predDf = data.frame(var = c("cp", "cCSF"),
                                     stringsAsFactors = FALSE))
  expect_equal(admixr2:::.admOutputVars(ui_two), c("cp", "cCSF"))

  ui_lin <- list(predDf = data.frame(var = "linCmtB", stringsAsFactors = FALSE))
  expect_equal(admixr2:::.admOutputVars(ui_lin), "ipredSim")
})

# ---- .admParseIniDf: sigma -> output mapping ---------------------------------

test_that("pinfo$sigma_output comes from the iniDf condition column", {
  iniDf <- data.frame(
    name  = c("tcl", "prop.cp", "add.cCSF", "eta.cl"),
    est   = c(1,      0.2,       0.1,        0.09),
    lower = c(-Inf,   0,         0,          NA),
    upper = c(Inf,    Inf,       Inf,        NA),
    fix   = c(FALSE,  FALSE,     FALSE,      FALSE),
    neta1 = c(NA,     NA,        NA,         1L),
    neta2 = c(NA,     NA,        NA,         1L),
    err   = c(NA,     "prop",    "add",      NA),
    condition = c(NA,  "cp",     "cCSF",     NA),
    stringsAsFactors = FALSE
  )
  pinfo <- admixr2:::.admParseIniDf(iniDf)
  expect_equal(pinfo$sigma_names,  c("prop.cp", "add.cCSF"))
  expect_equal(pinfo$sigma_output, c("cp", "cCSF"))
  expect_equal(admixr2:::.admSigmaSel(pinfo, "cCSF"), c(FALSE, TRUE))
})

test_that("pinfo$sigma_output is all-NA when iniDf has no condition column", {
  iniDf <- data.frame(
    name  = c("tcl", "prop.err", "eta.cl"),
    est   = c(1,      0.2,        0.09),
    lower = c(-Inf,   0,          NA),
    upper = c(Inf,    Inf,        NA),
    fix   = c(FALSE,  FALSE,      FALSE),
    neta1 = c(NA,     NA,         1L),
    neta2 = c(NA,     NA,         1L),
    err   = c(NA,     "prop",     NA),
    stringsAsFactors = FALSE
  )
  pinfo <- admixr2:::.admParseIniDf(iniDf)
  expect_true(all(is.na(pinfo$sigma_output)))
  expect_equal(admixr2:::.admSigmaSel(pinfo, "cp"), TRUE)
})

# ---- Joint (same-subject) data model: .admBuildJointUnit ---------------------

test_that("cross pairs assemble a symmetric joint V with correct blocks", {
  Vp <- matrix(c(0.2, 0.05, 0.05, 0.15), 2, 2)
  Vc <- matrix(0.1, 1, 1)
  Cpc <- matrix(c(0.03, 0.02), 2, 1)   # 2 (plasma times) x 1 (csf time)
  study <- list(n = 40L, ev = "EV", observations = list(
    plasma = list(output = "cp",   times = c(1, 2), E = c(3, 2), V = Vp),
    csf    = list(output = "cCSF", times = 4,       E = 0.5,     V = Vc)),
    cross = list("plasma:csf" = Cpc))
  ns <- admixr2:::.admNormaliseStudy(study, "rat", "cp")
  u  <- ns$observations[["rat"]]
  expect_true(isTRUE(u$is_joint))
  expect_equal(u$n_total, 3L)
  expect_equal(u$E, c(3, 2, 0.5))
  expect_equal(u$row_output, c(1L, 1L, 2L))
  # diagonal blocks
  expect_equal(u$V[1:2, 1:2], Vp)
  expect_equal(u$V[3, 3, drop = FALSE], Vc)
  # cross block + symmetry
  expect_equal(u$V[1:2, 3, drop = FALSE], Cpc)
  expect_equal(u$V[3, 1:2], as.numeric(Cpc))
  expect_equal(u$V, t(u$V))
  # block metadata
  expect_equal(u$blocks[[1]]$output, "cp")
  expect_equal(u$blocks[[2]]$output, "cCSF")
  expect_equal(u$blocks[[1]]$rows, 1:2)
  expect_equal(u$blocks[[2]]$rows, 3L)
})

test_that("study-level full V is used directly for a joint study", {
  Vj <- diag(c(0.2, 0.15, 0.1)); Vj[1, 3] <- Vj[3, 1] <- 0.02
  study <- list(n = 30L, ev = "EV", observations = list(
    plasma = list(output = "cp",   times = c(1, 2), E = c(3, 2)),
    csf    = list(output = "cCSF", times = 4,       E = 0.5)),
    V = Vj)
  u <- admixr2:::.admNormaliseStudy(study, "rat", "cp")$observations[["rat"]]
  expect_true(u$is_joint)
  expect_equal(u$V, Vj)
  expect_equal(u$E, c(3, 2, 0.5))
})

test_that("joint times are sorted and E reordered to match", {
  study <- list(n = 20L, ev = "EV", observations = list(
    plasma = list(output = "cp", times = c(4, 1, 2), E = c(30, 10, 20),
                  V = diag(c(0.3, 0.1, 0.2))),
    csf    = list(output = "cCSF", times = 8, E = 1, V = matrix(0.05, 1, 1))),
    cross = NULL)
  # cross=NULL but study-level neither -> still joint? Needs cross or V. Use V:
  study$cross <- NULL
  study$V <- diag(c(0.1, 0.2, 0.3, 0.05))   # order = sorted plasma times 1,2,4 then csf 8
  u <- admixr2:::.admNormaliseStudy(study, "rat", "cp")$observations[["rat"]]
  expect_equal(u$blocks[[1]]$times, c(1, 2, 4))
  expect_equal(u$E[1:3], c(10, 20, 30))       # reordered to sorted times
})

test_that("non-positive-definite assembled joint V errors clearly", {
  Vp <- diag(c(0.2, 0.15)); Vc <- matrix(0.1, 1, 1)
  Cpc <- matrix(c(5, 5), 2, 1)   # absurd cross -> not PD
  study <- list(n = 40L, ev = "EV", observations = list(
    plasma = list(output = "cp",   times = c(1, 2), E = c(3, 2), V = Vp),
    csf    = list(output = "cCSF", times = 4,       E = 0.5,     V = Vc)),
    cross = list("plasma:csf" = Cpc))
  expect_error(admixr2:::.admNormaliseStudy(study, "rat", "cp"),
               regexp = "positive-definite")
})

test_that("cross validation: bad pair name and wrong block dims error", {
  base <- list(n = 40L, ev = "EV", observations = list(
    plasma = list(output = "cp",   times = c(1, 2), E = c(3, 2), V = diag(2)),
    csf    = list(output = "cCSF", times = 4,       E = 0.5,     V = matrix(0.1,1,1))))
  bad_name <- base; bad_name$cross <- list("plasma:xxx" = matrix(0, 2, 1))
  expect_error(admixr2:::.admNormaliseStudy(bad_name, "rat", "cp"),
               regexp = "must name two observations")
  bad_dim <- base; bad_dim$cross <- list("plasma:csf" = matrix(0, 3, 1))
  expect_error(admixr2:::.admNormaliseStudy(bad_dim, "rat", "cp"),
               regexp = "cross block")
})

test_that("joint residual applies each output's sigma to its own rows", {
  # 3 rows: rows 1-2 = output cp (prop), row 3 = output cCSF (add)
  unit <- list(blocks = list(
    list(output = "cp",   rows = 1:2),
    list(output = "cCSF", rows = 3L)))
  pinfo <- list(sigma_names = c("prop.cp", "add.cCSF"),
                sigma_output = c("cp", "cCSF"),
                sigma_is_prop  = c(TRUE, FALSE),
                sigma_is_lnorm = c(FALSE, FALSE))
  mu <- c(4, 2, 1); V <- diag(3)
  jr <- admixr2:::.admJointResidual(mu, V, unit, pinfo, c(0.09, 0.05))
  # Law of total variance: a PROPORTIONAL residual contributes 0.09*E[f^2] =
  # 0.09*(var_f + mu^2), not 0.09*mu^2 -- var_f is 1 here (V = diag(3)). The
  # ADDITIVE row is unaffected, since its variance does not depend on f.
  expect_equal(diag(jr$V), c(1 + 0.09*(1 + 16), 1 + 0.09*(1 + 4), 1 + 0.05))
})

# ---- endpoint NAMES vs solve COLUMNS -----------------------------------------

test_that(".admEndpointNames gives nlmixr2's endpoint names, not the solve column", {
  skip_if_not_installed("rxode2")
  # The two coincide for an ordinary endpoint and diverge for exactly the ones
  # .admEndpointVar() exists for: `y ~ pois(lam)` is SOLVED through `lam` while
  # nlmixr2 knows the endpoint as `y`. These names go into the DVID column of the
  # dummy frame handed to nlmixr2CreateOutputFromUi(), whose dvid->cmt translation
  # rejects anything that is not an endpoint -- so a converged fit died at the
  # output-building step with "'dvid'->'cmt' ... undefined compartment".
  fn <- function() {
    ini({ tcl <- log(5); tv <- log(20); a <- 0.5; eta.cl ~ 0.09 })
    model({ cl <- exp(tcl + eta.cl); v <- exp(tv)
            d/dt(central) <- -(cl / v) * central
            cp  <- central / v
            lam <- cp * 2 + 0.5
            cp ~ add(a)
            y  ~ pois(lam) })
  }
  ui <- suppressMessages(rxode2::rxode2(fn))
  expect_equal(admixr2:::.admOutputVars(ui),   c("cp", "lam"))   # solve columns
  expect_equal(admixr2:::.admEndpointNames(ui), c("cp", "y"))    # endpoint names

  # ... and this combination is refused rather than silently scoring Inf: the
  # multi-endpoint solve routes observations by compartment, and `lam` is not one.
  expect_error(admixr2:::.admCheckMixedEndpoints(ui),
               "cannot be combined with other endpoints")

  # a single-endpoint count model is the supported case and must stay silent
  fn1 <- function() {
    ini({ tcl <- log(5); tv <- log(20); eta.cl ~ 0.09 })
    model({ cl <- exp(tcl + eta.cl); v <- exp(tv)
            d/dt(central) <- -(cl / v) * central
            cp <- central / v
            y ~ pois(cp) })
  }
  ui1 <- suppressMessages(rxode2::rxode2(fn1))
  expect_silent(admixr2:::.admCheckMixedEndpoints(ui1))
  expect_equal(admixr2:::.admEndpointNames(ui1), "y")
  expect_equal(admixr2:::.admOutputVars(ui1),    "cp")
})
