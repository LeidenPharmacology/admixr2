test_that("optimizer xtol_rel is exposed and plumbed to nloptr", {
  controls <- list(admControl, adfoControl, adghControl, adirmcControl)
  for (control in controls) {
    expect_identical(tail(setdiff(names(formals(control)), "..."), 1L), "xtol_rel")
    expect_identical(control(xtol_rel = 1e-9)$xtol_rel, 1e-9)
  }

  shared <- paste(deparse(body(admixr2:::.admScaledOptimize)), collapse = "\n")
  irmc <- paste(deparse(body(admixr2:::.adirmcPhaseLoop)), collapse = "\n")
  expect_match(shared, "xtol_rel = xtol_rel", fixed = TRUE)
  expect_match(irmc, "xtol_rel = xtol_rel", fixed = TRUE)

  for (driver in list(admixr2:::nlmixr2Est.admc, admixr2:::nlmixr2Est.adfo,
                      admixr2:::nlmixr2Est.adgh, admixr2:::nlmixr2Est.adirmc))
    expect_match(paste(deparse(body(driver)), collapse = "\n"),
                 ".ctl$xtol_rel", fixed = TRUE)

  restarts <- paste(deparse(body(admixr2:::.admRunRestarts)), collapse = "\n")
  expect_match(restarts, "pinfo$.xtol_rel <- .ctl$xtol_rel", fixed = TRUE)
})
