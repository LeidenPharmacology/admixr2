## ============================================================================
## Covariate marginalisation, end to end, through the real estimators
## ============================================================================
##
##   Rscript validation/covariate-collapse-endtoend.R
##
## Three studies of a meta-analysis, each enrolling a DIFFERENT population:
## subjects within a study span a covariate distribution, and only the study's
## aggregate mean vector + covariance matrix are observed.
##
## Covariate columns in the post-fit dummy frame are added by the driver
## (.admDummyData), so admData() is passed unchanged.
##
## The reference data are built from an ANALYTIC 1-cmt solution in plain R, by
## brute-force simulation over (covariate, random effect). Nothing in the data
## generation touches the code under test, so a bug in the collapse cannot hide
## itself by also corrupting the reference.
##
## Fits compared, both with adgh and admc:
##   collapse   study carries `cov_dist`  -> Omega inflated to Omega + J Sig J'
##   mean-only  study carries `cov` only  -> covariate at its mean, spread ignored
##                                           (what you get today without nodes)

suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(42)

## ---- truth -----------------------------------------------------------------
TCL <- log(1.0); TV <- log(10); TCOV <- 0.75      # covariate coefficient
OM  <- 0.30                                        # IIV SD on log CL
ADD <- 0.30                                        # additive residual SD
DOSE <- 100; TIMES <- c(0.5, 1, 2, 4, 8)

## Studies differ in BOTH covariate mean and spread. That matters: from marginal
## moments alone theta_cov and omega enter V only through theta_cov^2*sd_a^2 +
## omega^2, so a single population would not separate them. Different means give
## a mean shift that pins theta_cov; different spreads separate the variances.
POPS <- list(low  = list(mu = -0.40, sd = 0.25, n = 240L),
             mid  = list(mu =  0.00, sd = 0.60, n = 300L),
             high = list(mu =  0.45, sd = 0.30, n = 260L))

conc <- function(cl) DOSE / exp(TV) * exp(outer(-cl / exp(TV), TIMES))

## Reference aggregate moments by EXACT nested 2-D Gauss-Hermite over (a, b),
## on the analytic solution. Deliberately NOT the collapse identity and not Monte
## Carlo: an MC reference at 2e5 draws is noise-dominated when the covariate
## effect is small (it made tcov = 0.15 read as 0.095), and using the collapse to
## build the data the collapse is then tested on would prove nothing.
.gh <- function(k) {
  i <- seq_len(k - 1L)
  J <- matrix(0, k, k); J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2)
}
QQ <- .gh(40L)
ref_study <- function(p, tc = TCOV) {
  m1 <- 0; M2 <- 0
  for (ia in seq_along(QQ$x)) for (ib in seq_along(QQ$x)) {
    w <- QQ$w[ia] * QQ$w[ib]
    Y <- as.numeric(conc(exp(TCL + tc * (p$mu + p$sd * QQ$x[ia]) + OM * QQ$x[ib])))
    m1 <- m1 + w * Y; M2 <- M2 + w * outer(Y, Y)
  }
  V <- M2 - outer(m1, m1)
  diag(V) <- diag(V) + ADD^2
  list(E = m1, V = V, n = p$n, times = TIMES, ev = rxode2::et(amt = DOSE))
}
ref <- lapply(POPS, ref_study)

## ---- model -----------------------------------------------------------------
mod <- function() {
  ini({ tcl <- log(1.0); tv <- log(10); tcov <- 0.5
        eta.cl ~ 0.09
        add.err <- 0.3 })
  model({ cl <- exp(tcl + tcov * WT + eta.cl)
          v  <- exp(tv)
          d/dt(centr) <- -cl / v * centr
          cp <- centr / v
          cp ~ add(add.err) })
}

## ---- study lists -----------------------------------------------------------
## mean-only: covariate supplied at its mean, spread discarded
st_mean <- Map(function(r, p) c(r, list(cov = list(WT = p$mu))), ref, POPS)
## collapse: same, plus the distribution those subjects span
st_coll <- Map(function(r, p) c(r, list(cov      = list(WT = p$mu),
                                        cov_dist = list(WT = list(mu = p$mu, sd = p$sd)))),
               ref, POPS)

fit_one <- function(est, studies) {
  ctl <- switch(est,
    adgh = adghControl(studies = studies, grad = "none", n_nodes = 9L,
                       print = 0L, covMethod = "none"),
    admc = admControl(studies = studies, grad = "none", n_sim = 4000L,
                      print = 0L, covMethod = "none"))
  t0 <- Sys.time()
  f  <- try(suppressWarnings(suppressMessages(
          nlmixr2est::nlmixr2(mod, admData(), est = est, control = ctl))), silent = TRUE)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(f, "try-error")) return(list(err = conditionMessage(attr(f, "condition"))))
  p <- setNames(f$parFixedDf$Estimate, rownames(f$parFixedDf))
  list(tcl = unname(p["tcl"]), tv = unname(p["tv"]), tcov = unname(p["tcov"]),
       om = sqrt(as.numeric(f$omega[1, 1])), obj = f$objective, secs = el)
}

row <- function(lbl, r) {
  if (!is.null(r$err)) { cat(sprintf("%-18s  ERROR: %s\n", lbl, r$err)); return(invisible()) }
  cat(sprintf("%-18s %9.4f %9.4f %9.4f %9.4f %11.1f %8.1f\n",
              lbl, r$tcl, r$tv, r$tcov, r$om, r$obj, r$secs))
}

cat("\n=========================================================================\n")
cat("Recovery across three populations (n_draw = 4e5 reference draws)\n")
cat("=========================================================================\n")
cat(sprintf("%-18s %9s %9s %9s %9s %11s %8s\n",
            "", "tcl", "tv", "tcov", "omega", "objective", "secs"))
cat(sprintf("%-18s %9.4f %9.4f %9.4f %9.4f\n", "TRUTH", TCL, TV, TCOV, OM))
cat(strrep("-", 78), "\n")
for (est in c("adgh", "admc")) {
  row(paste0(est, " collapse"),  fit_one(est, st_coll))
  row(paste0(est, " mean-only"), fit_one(est, st_mean))
}

## ---- how wrong is mean-only, as the covariate effect grows? ----------------
cat("\n=========================================================================\n")
cat("Sensitivity: covariate effect vs IIV  (adgh, tcov estimated)\n")
cat("=========================================================================\n")
cat("ratio = tcov*sd_a/omega, using the mid population's spread.\n\n")
cat(sprintf("%7s %9s | %19s | %19s\n", "ratio", "tcov true",
            "collapse", "mean-only"))
cat(sprintf("%7s %9s | %9s %9s | %9s %9s\n", "", "", "tcov", "omega", "tcov", "omega"))
for (tc in c(0.15, 0.40, 0.75, 1.10)) {
  refs <- lapply(POPS, ref_study, tc = tc)
  sc <- Map(function(r, p) c(r, list(cov = list(WT = p$mu),
                                     cov_dist = list(WT = list(mu = p$mu, sd = p$sd)))),
            refs, POPS)
  sm <- Map(function(r, p) c(r, list(cov = list(WT = p$mu))), refs, POPS)
  a1 <- fit_one("adgh", sc); a2 <- fit_one("adgh", sm)
  f <- function(x, d) if (is.null(x$err)) sprintf("%9.4f", x[[d]]) else "    ERROR"
  cat(sprintf("%7.2f %9.2f | %s %s | %s %s\n", tc * POPS$mid$sd / OM, tc,
              f(a1, "tcov"), f(a1, "om"), f(a2, "tcov"), f(a2, "om")))
}
cat("\n")
