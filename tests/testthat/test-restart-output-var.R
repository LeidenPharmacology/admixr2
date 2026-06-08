test_that("restart workers expose output_var argument with cp default", {
  expect_true("output_var" %in% names(formals(admixr2:::.adfoRestartWorker)))
  expect_true("output_var" %in% names(formals(admixr2:::.admRestartWorker)))
  expect_true("output_var" %in% names(formals(admixr2:::.adirmcRestartWorker)))

  expect_identical(eval(formals(admixr2:::.adfoRestartWorker)$output_var), "cp")
  expect_identical(eval(formals(admixr2:::.admRestartWorker)$output_var), "cp")
  expect_identical(eval(formals(admixr2:::.adirmcRestartWorker)$output_var), "cp")
})

test_that("restart estimators pass detected output_var into restart workers", {
  adfo_txt <- paste(deparse(body(admixr2:::nlmixr2Est.adfo)), collapse = "\n")
  admc_txt <- paste(deparse(body(admixr2:::nlmixr2Est.admc)), collapse = "\n")
  irmc_txt <- paste(deparse(body(admixr2:::nlmixr2Est.adirmc)), collapse = "\n")

  pat <- "extra_args\\s*=\\s*list\\([\\s\\S]*output_var\\s*=\\s*output_var"

  expect_true(grepl(pat, adfo_txt, perl = TRUE))
  expect_true(grepl(pat, admc_txt, perl = TRUE))
  expect_true(grepl(pat, irmc_txt, perl = TRUE))
})
