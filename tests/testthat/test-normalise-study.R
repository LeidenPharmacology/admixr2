test_that("Vector V expanded to diagonal matrix with method = 'var'", {
  s  <- list(E = c(1.0, 2.0), V = c(0.1, 0.2), n = 50L, times = c(1, 2))
  ns <- admixr2:::.admNormaliseStudy(s, "s1")

  expect_equal(ns$method, "var")
  expect_equal(ns$V, diag(c(0.1, 0.2)))
  expect_equal(ns$v_diag, c(0.1, 0.2))
})

test_that("Diagonal matrix auto-detected as method = 'var'", {
  s  <- list(E = c(1.0, 2.0), V = diag(c(0.1, 0.2)), n = 50L, times = c(1, 2))
  ns <- admixr2:::.admNormaliseStudy(s, "s2")
  expect_equal(ns$method, "var")
  expect_equal(ns$v_diag, c(0.1, 0.2))
})

test_that("Full matrix auto-detected as method = 'cov'", {
  V  <- matrix(c(0.1, 0.02, 0.02, 0.2), 2, 2)
  s  <- list(E = c(1.0, 2.0), V = V, n = 50L, times = c(1, 2))
  ns <- admixr2:::.admNormaliseStudy(s, "s3")
  expect_equal(ns$method, "cov")
  expect_null(ns$v_diag)
})

test_that("v_diag only set for 'var' studies", {
  V  <- matrix(c(0.1, 0.02, 0.02, 0.2), 2, 2)
  s  <- list(E = c(1.0, 2.0), V = V, n = 50L, times = c(1, 2))
  ns <- admixr2:::.admNormaliseStudy(s, "s4")
  expect_null(ns$v_diag)
})

test_that("Vector V + method='cov' warns and coerces to 'var'", {
  s <- list(E = c(1.0, 2.0), V = c(0.1, 0.2), n = 50L, times = c(1, 2),
            method = "cov")
  expect_warning(
    ns <- admixr2:::.admNormaliseStudy(s, "s5"),
    regexp = "method='var'"
  )
  expect_equal(ns$method, "var")
})

test_that("Non-diagonal V + method='var' warns about off-diagonal entries", {
  V <- matrix(c(0.1, 0.02, 0.02, 0.2), 2, 2)
  s <- list(E = c(1.0, 2.0), V = V, n = 50L, times = c(1, 2), method = "var")
  expect_warning(
    admixr2:::.admNormaliseStudy(s, "s6"),
    regexp = "off-diagonal"
  )
})

test_that("Missing E stops with informative message", {
  s <- list(V = diag(2), n = 50L, times = c(1, 2))
  expect_error(admixr2:::.admNormaliseStudy(s, "study_x"), regexp = "missing 'E'")
})

test_that("Missing V stops with informative message", {
  s <- list(E = c(1, 2), n = 50L, times = c(1, 2))
  expect_error(admixr2:::.admNormaliseStudy(s, "study_x"), regexp = "missing 'V'")
})

test_that("Missing n stops with informative message", {
  s <- list(E = c(1, 2), V = diag(2), times = c(1, 2))
  expect_error(admixr2:::.admNormaliseStudy(s, "study_x"), regexp = "missing 'n'")
})

test_that("Missing times stops with informative message", {
  s <- list(E = c(1, 2), V = diag(2), n = 50L)
  expect_error(admixr2:::.admNormaliseStudy(s, "study_x"), regexp = "missing 'times'")
})

test_that("Explicit method='cov' on full matrix is respected", {
  V  <- matrix(c(0.1, 0.02, 0.02, 0.2), 2, 2)
  s  <- list(E = c(1.0, 2.0), V = V, n = 50L, times = c(1, 2), method = "cov")
  ns <- admixr2:::.admNormaliseStudy(s, "s7")
  expect_equal(ns$method, "cov")
})

test_that("Study name appears in error message", {
  s <- list(V = diag(2), n = 50L, times = c(1, 2))
  expect_error(admixr2:::.admNormaliseStudy(s, "my_study"), regexp = "my_study")
})

# -- long-format studies (data frame keyed by DVID / CMT) ----------------------

# One row per observed endpoint/time: the nlmixr2 way of keying observations.
.long_df <- function()
  data.frame(DVID = c("cp", "cp", "cp", "cb", "cb"),
             TIME = c(0.5, 1, 2, 1, 2),
             E    = c(8.8, 7.8, 6.9, 3.0, 3.4),
             V    = c(1.21, 0.81, 0.64, 0.09, 0.09))

test_that("Long format without a joint V gives one independent unit per endpoint", {
  ns <- admixr2:::.admNormaliseStudy(
    list(n = 60L, ev = "EV", data = .long_df()), "lit")

  expect_named(ns$observations, c("cp", "cb"))
  expect_false(isTRUE(ns$joint))
  expect_equal(ns$observations$cp$times, c(0.5, 1, 2))
  expect_equal(ns$observations$cb$E, c(3.0, 3.4))
  expect_equal(ns$observations$cb$output, "cb")
  # variance column -> diagonal V, method "var", n/ev inherited from the study
  expect_equal(ns$observations$cp$method, "var")
  expect_equal(ns$observations$cp$v_diag, c(1.21, 0.81, 0.64))
  expect_equal(ns$observations$cb$n, 60L)
  expect_equal(ns$observations$cb$ev, "EV")
})

test_that("Long format matches the equivalent observations spec", {
  long <- admixr2:::.admNormaliseStudy(
    list(n = 60L, ev = "EV", data = .long_df()), "lit")
  obs <- admixr2:::.admNormaliseStudy(list(n = 60L, ev = "EV", observations = list(
    cp = list(output = "cp", times = c(0.5, 1, 2), E = c(8.8, 7.8, 6.9),
              V = c(1.21, 0.81, 0.64)),
    cb = list(output = "cb", times = c(1, 2), E = c(3.0, 3.4), V = c(0.09, 0.09)))),
    "lit")

  for (k in c("cp", "cb"))
    expect_equal(long$observations[[k]][c("output", "times", "E", "V", "n", "method", "v_diag")],
                 obs$observations[[k]][c("output", "times", "E", "V", "n", "method", "v_diag")])
})

test_that("Long format + study-level V builds one joint (same-subject) unit", {
  V  <- diag(c(1.21, 0.81, 0.64, 0.09, 0.09))
  V[1, 4] <- V[4, 1] <- 0.15                    # plasma(t=0.5) <-> brain(t=1)
  ns <- admixr2:::.admNormaliseStudy(
    list(n = 60L, ev = "EV", data = .long_df(), V = V), "lit")

  expect_length(ns$observations, 1L)
  u <- ns$observations[[1L]]
  expect_true(isTRUE(u$is_joint))
  expect_equal(u$E, c(8.8, 7.8, 6.9, 3.0, 3.4))   # stacked, blocks in data order
  expect_equal(u$V, unname(V))                    # rows aligned to the data rows
  expect_equal(u$row_output, c(1L, 1L, 1L, 2L, 2L))
  expect_equal(vapply(u$blocks, `[[`, character(1), "output"), c("cp", "cb"))
  expect_equal(u$n, 60L)
})

test_that("Long format keeps E and V aligned when data rows are unordered", {
  d    <- .long_df()
  V    <- diag(c(1.21, 0.81, 0.64, 0.09, 0.09)); V[1, 4] <- V[4, 1] <- 0.15
  ord  <- c(4, 1, 5, 3, 2)                        # scramble the rows
  ns   <- admixr2:::.admNormaliseStudy(
    list(n = 60L, ev = "EV", data = d[ord, ], V = V[ord, ord]), "lit")
  u    <- ns$observations[[1L]]

  # blocks follow first appearance (cb now leads); E and V permute together
  stack <- c(4, 5, 1, 2, 3)
  expect_equal(vapply(u$blocks, `[[`, character(1), "output"), c("cb", "cp"))
  expect_equal(u$E, d$E[stack])
  expect_equal(u$V, unname(V[stack, stack]))
})

test_that("Long format accepts an SD column and per-endpoint n / ev", {
  d  <- data.frame(CMT = c("cp", "cp", "cb"), time = c(1, 2, 1),
                   mean = c(8.8, 7.8, 3.0), SD = c(1.1, 0.9, 0.3),
                   n = c(60, 60, 12))
  ns <- admixr2:::.admNormaliseStudy(
    list(ev = list(cp = "EVp", cb = "EVb"), data = d), "lit")

  expect_equal(ns$observations$cp$v_diag, c(1.21, 0.81))   # SD^2
  expect_equal(ns$observations$cp$n, 60)
  expect_equal(ns$observations$cb$n, 12)                   # separate experiment
  expect_equal(ns$observations$cb$ev, "EVb")
})

test_that("Long-format input errors are informative", {
  d <- .long_df()
  expect_error(admixr2:::.admNormaliseStudy(
    list(n = 1L, ev = "EV", data = d[, c("DVID", "TIME")]), "lit"),
    regexp = "mean column")
  expect_error(admixr2:::.admNormaliseStudy(
    list(n = 1L, ev = "EV", data = d[, c("DVID", "TIME", "E")]), "lit"),
    regexp = "variance column")
  expect_error(admixr2:::.admNormaliseStudy(
    list(n = 1L, ev = "EV", data = rbind(d, d[1, ])), "lit"),
    regexp = "duplicate endpoint/time")
  expect_error(admixr2:::.admNormaliseStudy(
    list(n = 1L, ev = "EV", data = d, V = diag(4)), "lit"),
    regexp = "must be 5 x 5")
  # a joint study is one experiment: it cannot carry per-endpoint n or ev
  expect_error(admixr2:::.admNormaliseStudy(
    list(ev = "EV", data = transform(d, n = ifelse(DVID == "cp", 60, 12)),
         V = diag(5)), "lit"),
    regexp = "one shared `n`")
  expect_error(admixr2:::.admNormaliseStudy(
    list(n = 1L, ev = list(cp = "a", cb = "b"), data = d, V = diag(5)), "lit"),
    regexp = "shares one `ev`")
})


test_that("a second normalisation fills a JOINT unit's per-block output", {
  skip_if_not_installed("rxode2")
  # The fixtures (and anything that prepares studies before a model exists)
  # normalise without a default_output, so every blk$output comes out NULL --
  # .admBuildJointUnit() takes each from `ob$output %||% default_output`. The
  # driver's later pass is the ONLY chance to fill them, and the short-circuit
  # for an already-normalised study used to skip joint units entirely, making the
  # NULL permanent. .admBuildEvFull() then calls et(blk$times, cmt = NULL) per
  # block: the joint sens solve either errors (dropping to FD) or reads an
  # untagged compartment -- a finite but wrong joint objective, silently.
  ev <- rxode2::et(amt = 100)
  raw <- list(n = 40L, ev = ev, joint = TRUE,
              observations = list(
                a = list(times = c(1, 2), E = c(5, 4), V = diag(c(1, 1))),
                b = list(times = c(3),    E = 3,       V = matrix(1))))

  once <- admixr2:::.admNormaliseStudy(raw, "s")               # no model yet
  u1   <- once$observations[[1L]]
  expect_true(isTRUE(u1$is_joint))
  expect_null(u1$blocks[[1L]]$output)                          # the NULL in question

  twice <- admixr2:::.admNormaliseStudy(once, "s", default_output = "cp")
  u2    <- twice$observations[[1L]]
  expect_true(isTRUE(u2$is_joint))
  expect_identical(vapply(u2$blocks, `[[`, character(1), "output"),
                   c("cp", "cp"))
  # The unit's own `output` exists only for cmt-tagging and is taken from
  # blocks[[1]], not stamped independently.
  expect_identical(u2$output, "cp")
  # row_output holds block INDICES, so it must be untouched by any of this.
  expect_identical(u2$row_output, once$observations[[1L]]$row_output)
})

test_that("every unit records the STUDY it came from", {
  # `label` identifies a unit, `study` identifies the trial, and they are not
  # the same thing as soon as a study has several outputs. Study-level effects
  # -- a per-study baseline, a random study effect, within- vs between-study
  # covariate coefficients -- all group by trial, and once units are flattened
  # the trial is unrecoverable unless it was recorded here.
  E <- c(1, 2); V <- diag(2); tm <- c(1, 2)
  one <- admixr2:::.admNormaliseStudy(
    list(E = E, V = V, n = 10L, times = tm), "trialA")
  expect_identical(one$observations[[1L]]$study, "trialA")

  multi <- admixr2:::.admNormaliseStudy(
    list(n = 10L, ev = rxode2::et(amt = 1), observations = list(
      plasma = list(output = "cp", times = tm, E = E, V = V),
      csf    = list(output = "cc", times = tm, E = E, V = V))), "trialB")
  us <- multi$observations
  expect_true(all(vapply(us, function(u) u$study, character(1)) == "trialB"))
  # ... while the labels stay distinct, which is the point
  expect_equal(length(unique(vapply(us, function(u) u$label, character(1)))), 2L)

  # grouping survives flattening, which is where `label` alone would lose it
  flat <- admixr2:::.admFlattenStudies(list(trialA = one, trialB = multi))
  g <- admixr2:::.admStudyGroups(flat)
  expect_identical(names(g), c("trialA", "trialB"))
  expect_identical(lengths(g), c(trialA = 1L, trialB = 2L))

  # a hand-built unit that never declared a study reads as its own study,
  # which is the right default for one unit per study
  expect_identical(admixr2:::.admUnitStudy(list(label = "solo")), "solo")

  # idempotent: normalising twice must not disturb it
  expect_identical(
    admixr2:::.admNormaliseStudy(one, "trialA")$observations[[1L]]$study,
    "trialA")
})

test_that("v_denom converts a published (n-1) covariance to the ML one", {
  # The two input types disagree about what V is: a digitised figure gives
  # V = SD^2 with SD the UNBIASED sample SD, while cov.wt(method = "ML") and
  # datagen() give the n covariance the likelihood is exact for. Declaring it
  # per study is what lets one meta-analysis mix both.
  N <- 20L; f <- (N - 1) / N
  base <- list(E = c(10, 8, 6), n = N, times = c(1, 2, 4),
               ev = rxode2::et(amt = 100))
  a <- admixr2:::.admNormaliseStudy(c(base, list(V = c(4, 2.25, 1))), "s")
  b <- admixr2:::.admNormaliseStudy(
    c(base, list(V = c(4, 2.25, 1), v_denom = "unbiased")), "s")
  expect_equal(diag(b$V), diag(a$V) * f)
  expect_equal(b$v_diag, diag(b$V))          # the cached diagonal tracks it
  expect_identical(b$v_denom, "ml")          # records that it has been applied

  # idempotent: normalising twice must not convert twice
  expect_equal(admixr2:::.admNormaliseStudy(b, "s")$V, b$V)

  # a full covariance converts off-diagonals too
  Vf <- diag(c(4, 2.25, 1)); Vf[1, 2] <- Vf[2, 1] <- 0.5
  cf <- admixr2:::.admNormaliseStudy(c(base, list(V = Vf, v_denom = "unbiased")), "s")
  expect_equal(cf$V, Vf * f)

  # the default is unchanged behaviour
  expect_equal(a$V, diag(c(4, 2.25, 1)))
})

test_that("v_denom reaches the paths that bypass .admNormaliseObs", {
  # A joint study assembles its covariance in .admBuildJointUnit from the RAW
  # blocks, so converting inside the per-observation normaliser would miss it.
  N <- 25L; f <- (N - 1) / N; ev <- rxode2::et(amt = 100)
  Vj <- diag(c(4, 2, 1, 0.5)); Vj[1, 3] <- Vj[3, 1] <- 0.4
  jt <- function(vd) admixr2:::.admNormaliseStudy(list(
    n = N, ev = ev, V = Vj, joint = TRUE, v_denom = vd,
    observations = list(
      plasma = list(output = "cp", times = c(1, 2), E = c(9, 7)),
      csf    = list(output = "cc", times = c(2, 8), E = c(3, 1)))), "s")
  ja <- jt("ml")$observations[[1L]]; jb <- jt("unbiased")$observations[[1L]]
  expect_true(isTRUE(jb$is_joint))
  expect_equal(jb$V, ja$V * f)               # diagonal AND cross blocks

  # multi-output uses each observation's OWN n, not the study's
  mo <- function(vd) admixr2:::.admNormaliseStudy(list(
    n = N, ev = ev, v_denom = vd,
    observations = list(
      plasma = list(output = "cp", times = c(1, 2), E = c(9, 7),
                    V = c(4, 2), n = N),
      csf    = list(output = "cc", times = c(2, 8), E = c(3, 1),
                    V = c(1, 0.5), n = 10L))), "s")
  ma <- mo("ml"); mb <- mo("unbiased")
  expect_equal(diag(mb$observations$plasma$V), diag(ma$observations$plasma$V) * f)
  expect_equal(diag(mb$observations$csf$V),
               diag(ma$observations$csf$V) * (10 - 1) / 10)
})

test_that("v_denom refuses what it cannot convert", {
  base <- list(E = c(10, 8), times = c(1, 2), ev = rxode2::et(amt = 100))
  expect_error(admixr2:::.admNormaliseStudy(
    c(base, list(V = c(4, 2), v_denom = "unbiased")), "s"), "n")
  expect_error(admixr2:::.admNormaliseStudy(
    c(base, list(V = c(4, 2), n = 1L, v_denom = "unbiased")), "s"), "n > 1")
  expect_error(admixr2:::.admNormaliseStudy(
    c(base, list(V = c(4, 2), n = 20L, v_denom = "n-1")), "s"), "unbiased")
})
