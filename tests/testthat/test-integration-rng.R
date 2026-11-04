test_that("the RNG snapshot restores a stream and an unseeded session", {
  set.seed(1); expected <- runif(1)
  set.seed(1)
  restore <- admixr2:::.admRngSnapshot()
  set.seed(999); runif(5)
  restore()
  expect_identical(runif(1), expected)

  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    rm(".Random.seed", envir = globalenv())
  restore <- admixr2:::.admRngSnapshot()
  set.seed(999)
  restore()
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  set.seed(1)
})

test_that("fits and datagen leave the caller's RNG stream alone", {
  skip_on_cran()
  skip_if_not_installed("rxode2")

  mod <- function() {
    ini({ tcl <- log(2); tv <- log(10); eta.cl ~ 0.09; add.err <- 0.3 })
    model({
      cl <- exp(tcl + eta.cl); v <- exp(tv)
      d/dt(central) <- -cl / v * central
      cp <- central / v
      cp ~ add(add.err)
    })
  }
  times <- c(1, 2, 4, 8)
  E <- 10 * exp(-0.2 * times)
  st <- list(s1 = list(E = E, V = diag(0.3 + 0.02 * E^2), n = 100L, times = times,
                       ev = rxode2::et(amt = 100)))
  ctl <- list(
    adfo   = admixr2::adfoControl(studies = st, covMethod = "none", cores = 2L, maxeval = 5L),
    adgh   = admixr2::adghControl(studies = st, covMethod = "none", cores = 2L, maxeval = 5L),
    admc   = admixr2::admControl(studies = st, n_sim = 100L, covMethod = "none", cores = 2L, maxeval = 5L),
    adirmc = admixr2::adirmcControl(studies = st, n_sim = 100L, covMethod = "none", cores = 2L, outer_iter = 2L)
  )
  untouched <- function(code) {
    set.seed(101); expected <- runif(1)
    set.seed(101)
    suppressMessages(suppressWarnings(force(code)))
    expect_identical(runif(1), expected)
  }
  for (est in names(ctl))
    untouched(nlmixr2est::nlmixr2(mod, admixr2::admData(), est = est, control = ctl[[est]]))

  for (s in c("rnorm", "lhs"))
    untouched(admixr2::datagen(st, mod, admixr2::datagenControl(sampling = s, n_sim = 100L)))
})
