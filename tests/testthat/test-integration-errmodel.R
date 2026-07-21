# End-to-end: error models that admixr2 could not previously represent.
#
# The unit tests in test-errmodel.R pin the maths. These check the whole chain --
# rxode2 ui -> predDf -> pinfo -> NLL/gradient -- because the original defect was
# not in the maths but in the PARSE: admixr2 read only iniDf$err and never
# ui$predDf, so it never saw errType/addProp at all.

.em_model <- function(err_line, extra_ini = "") {
  txt <- sprintf('function() {
    ini({
      tcl <- log(5)
      tv  <- log(20)
      %s
      eta.cl ~ 0.09
    })
    model({
      cl <- exp(tcl + eta.cl)
      v  <- exp(tv)
      d/dt(central) <- -(cl / v) * central
      cp <- central / v
      %s
    })
  }', extra_ini, err_line)
  rxode2::rxode2(eval(parse(text = txt)))
}

.em_studies <- function(ui) {
  tp <- c(0.5, 1, 2, 4, 8)
  E  <- (100 / 20) * exp(-(5 / 20) * tp)
  raw <- list(observations = list(
    plasma = list(output = "cp", ev = rxode2::et(amt = 100), times = tp,
                  E = E, V = diag((0.25 * E)^2), n = 40L)))
  st <- list(m1 = admixr2:::.admNormaliseStudy(raw, "m1", "cp"))
  admixr2:::.admBuildEvFull(admixr2:::.admFlattenStudies(st), tag_cmt = FALSE)
}

# Central-FD gradient of the FO NLL: the independent reference for the analytical
# sigma gradient. If the analytical residual derivatives are wrong, this catches it.
.em_fd_grad <- function(p, pinfo, studies, sensModel, rxMod, pl, h = 1e-5) {
  vapply(seq_along(p), function(k) {
    hi <- p; hi[k] <- p[k] + h
    lo <- p; lo[k] <- p[k] - h
    (admixr2:::.adfoNLL(hi, pinfo, studies, sensModel, rxMod, "cp", pl, 1L) -
     admixr2:::.adfoNLL(lo, pinfo, studies, sensModel, rxMod, "cp", pl, 1L)) / (2 * h)
  }, numeric(1))
}

.em_check_grad <- function(ui, tol = 1e-3) {
  sensModel <- admixr2:::.admLoadSensModel(ui)
  rxMod     <- admixr2:::.admLoadModel(ui)
  pinfo     <- admixr2:::.admParseIniDf(ui$iniDf, ui)
  ov        <- admixr2:::.admBuildOptVec(pinfo)
  studies   <- .em_studies(ui)
  pl        <- admixr2:::.admMakeParamsList(1L, pinfo, length(studies))

  g_ana <- admixr2:::.adfoGrad(ov$p0, pinfo, studies, sensModel, rxMod, "cp",
                               pl, 1L, 1e-4)
  g_fd  <- .em_fd_grad(ov$p0, pinfo, studies, sensModel, rxMod, pl)

  # Compare only the residual parameters: the struct thetas are themselves FD in
  # adfo, so comparing them against FD is circular.
  n_s   <- length(pinfo$struct_names)
  n_e   <- length(pinfo$sigma_names)
  slots <- n_s + seq_len(n_e)
  expect_equal(unname(g_ana[slots]), unname(g_fd[slots]),
               tolerance = tol,
               info = paste("sigma gradient:", paste(pinfo$sigma_names, collapse = ", ")))
  list(pinfo = pinfo, grad = g_ana)
}

test_that("pow(): the exponent is parsed as an exponent, and the fit is not silently additive", {
  skip_on_cran()
  skip_if_not_installed("rxode2")

  ui    <- .em_model("cp ~ pow(pow.sd, pow.c)", "pow.sd <- 0.2; pow.c <- 1.1")
  pinfo <- admixr2:::.admParseIniDf(ui$iniDf, ui)

  # rxode2 emits err = "pow" (coefficient) and err = "pow2" (exponent).
  expect_setequal(pinfo$sigma_names, c("pow.sd", "pow.c"))
  k <- match("pow.c", pinfo$sigma_names)
  expect_identical(unname(admixr2:::.admSigmaRole(pinfo)[k]), "pow_exp")

  # Before the fix the exponent was treated as an additive VARIANCE: it entered
  # the optimizer as 2*log(1.1) and contributed exp(that) = 1.21 to diag(V).
  expect_equal(unname(pinfo$sigma_init[k]), 1.1)

  spec <- admixr2:::.admResidSpecFor(pinfo, "cp")
  expect_equal(spec$form, 0L)                       # combined2
  expect_true(is.na(spec$k_add))                    # pow() alone has no additive term
  expect_equal(spec$k_pow, k)
})

test_that("pow(): analytical sigma gradient matches central FD of the NLL", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  .em_check_grad(.em_model("cp ~ pow(pow.sd, pow.c)", "pow.sd <- 0.2; pow.c <- 1.1"))
})

test_that("add + prop defaults to combined2 and its gradient matches FD", {
  skip_on_cran()
  skip_if_not_installed("rxode2")

  ui  <- .em_model("cp ~ add(add.sd) + prop(prop.sd)", "add.sd <- 0.1; prop.sd <- 0.2")
  res <- .em_check_grad(ui)
  spec <- admixr2:::.admResidSpecFor(res$pinfo, "cp")
  expect_equal(spec$form, 0L)   # nlmixr2's addProp default resolves to combined2
})

test_that("combined1() is honoured, not silently computed as combined2", {
  skip_on_cran()
  skip_if_not_installed("rxode2")

  ui <- .em_model("cp ~ add(add.sd) + prop(prop.sd) + combined1()",
                  "add.sd <- 0.1; prop.sd <- 0.2")
  # predDf$addProp is the authoritative selector; admixr2 used to ignore it.
  expect_identical(as.character(as.data.frame(ui$predDf)$addProp[1]), "combined1")

  res  <- .em_check_grad(ui)
  spec <- admixr2:::.admResidSpecFor(res$pinfo, "cp")
  expect_equal(spec$form, 1L)

  # ... and it really is a different variance from combined2.
  ui2  <- .em_model("cp ~ add(add.sd) + prop(prop.sd) + combined2()",
                    "add.sd <- 0.1; prop.sd <- 0.2")
  p2   <- admixr2:::.admParseIniDf(ui2$iniDf, ui2)
  expect_equal(admixr2:::.admResidSpecFor(p2, "cp")$form, 0L)

  f   <- c(1, 3, 6)
  sv  <- admixr2:::.admSigmaNat(admixr2:::.admBuildOptVec(res$pinfo)$p0[
           length(res$pinfo$struct_names) + seq_along(res$pinfo$sigma_names)], res$pinfo)
  a1  <- admixr2:::.admResidRows(res$pinfo, "cp", sv, length(f))
  a2  <- admixr2:::.admResidRows(p2,        "cp", sv, length(f))
  v1  <- admixr2:::.admResidApply(f, rep(0, 3), a1)$dv
  v2  <- admixr2:::.admResidApply(f, rep(0, 3), a2)$dv
  expect_true(all(v1 > v2))     # the 2ab*f cross term
})

test_that("add + pow: gradient matches FD across all three residual parameters", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  .em_check_grad(.em_model("cp ~ add(add.sd) + pow(pow.sd, pow.c)",
                           "add.sd <- 0.1; pow.sd <- 0.2; pow.c <- 1.1"))
})

test_that("an unrepresentable residual model is refused, not silently approximated", {
  skip_on_cran()
  skip_if_not_installed("rxode2")

  # logitNorm is now SUPPORTED (transform-both-sides, by quadrature), so the
  # refusal example is cauchy -- which is refused on mathematics rather than
  # missing machinery: it is the nu = 1 Student-t, has no finite mean or
  # variance, and averaging cannot rescue it (the Cauchy is a STABLE
  # distribution, so a study mean of n subjects is distributed exactly like one
  # subject and never concentrates).
  ui <- .em_model("cp ~ add(add.sd) + cauchy()", "add.sd <- 0.2")
  expect_error(admixr2:::.admParseIniDf(ui$iniDf, ui),
               regexp = "Unsupported residual error model")

  # ... and the ones that ARE representable now parse.
  for (spec in list(
    c("cp ~ logitNorm(lg.sd, 0, 10)",           "lg.sd <- 0.2"),
    c("cp ~ probitNorm(lg.sd, 0, 10)",          "lg.sd <- 0.2"),
    c("cp ~ add(add.sd) + boxCox(lam)",         "add.sd <- 0.2; lam <- 0.5"),
    c("cp ~ add(add.sd) + yeoJohnson(lam)",     "add.sd <- 0.2; lam <- 0.5"),
    c("cp ~ add(add.sd) + ar(rho)",             "add.sd <- 0.2; rho <- 0.5"))) {
    u <- .em_model(spec[1], spec[2])
    expect_no_error(suppressWarnings(admixr2:::.admParseIniDf(u$iniDf, u)))
  }
})

test_that("refusal messages name the model, the endpoint, the reason and the alternatives", {
  skip_on_cran()
  skip_if_not_installed("rxode2")

  # The refusal replaces a silent wrong fit, so the message has to carry its
  # weight: what was asked for, on which endpoint, WHY it cannot be done, and what
  # to use instead. A bare "unsupported" would just move the confusion.
  msg <- tryCatch(
    admixr2:::.admParseIniDf(
      .em_model("cp ~ add(add.sd) + cauchy()", "add.sd <- 0.2")$iniDf,
      .em_model("cp ~ add(add.sd) + cauchy()", "add.sd <- 0.2")),
    error = conditionMessage)

  expect_match(msg, "cauchy")              # what was asked for
  expect_match(msg, "endpoint 'cp'")       # where
  expect_match(msg, "Why:")                # why not
  expect_match(msg, "STABLE")              # ... specifically: why averaging fails
  expect_match(msg, "aggregate", ignore.case = TRUE)
  expect_match(msg, "add\\(a\\) \\+ prop\\(b\\)")   # what IS supported
  expect_match(msg, "ADDITIVE")            # warns that old results are suspect
  # no line runs off the terminal
  expect_true(all(nchar(strsplit(msg, "\n")[[1]]) <= 80))
})

test_that("propF() gets its own targeted message, not the generic one", {
  skip_on_cran()
  skip_if_not_installed("rxode2")

  # propF()'s second argument is a model variable; `v` is defined in .em_model().
  ui <- .em_model("cp ~ propF(prop.sd, v)", "prop.sd <- 0.2")
  msg <- tryCatch(admixr2:::.admParseIniDf(ui$iniDf, ui), error = conditionMessage)
  expect_match(msg, "propF")
  expect_match(msg, "Fix:")
  expect_match(msg, "prop\\(\\) or pow\\(\\)")
})
