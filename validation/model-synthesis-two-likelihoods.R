## TWO objectives for synthesising published models, side by side.
##
##   AGGREGATE   -- what admixr2 fits. Re-simulate each source over its own
##                  population/design, treat the result as data with weight n_s,
##                  score the mean-and-covariance likelihood.
##                  Data = predicted profiles.   Weight = n_s (a CHOICE).
##
##   PARAMETER   -- a formal likelihood on the PUBLISHED NUMBERS:
##                     theta_hat_s ~ N( g_s(phi), Sigma_s )
##                  g_s(phi) = the INDUCED map: what study s would have estimated,
##                  under its own model form, design and population, if the truth
##                  were phi. Computed by projecting the unified model ONTO the
##                  source's structure -- the inner problem run backwards.
##                  Data = published estimates.  Weight = Sigma_s (REPORTED).
##
## Key structural fact exploited below: when a source's model FORM equals the
## unified form, g_s is the IDENTITY and costs nothing. Only structurally reduced
## sources need an inner projection. Verified numerically, not assumed.

setwd("C:/package/admixr2/.claude/worktrees/feature-covariate-quadrature")
suppressMessages(devtools::load_all(".", quiet = TRUE)); suppressMessages(library(rxode2))
say <- function(...) { cat(..., "\n"); utils::flush.console() }
T0 <- proc.time()[["elapsed"]]
GEN <- datagenControl(n_sim = 8000L, seed = 11L)
KK  <- 4L

## ---------------------------------------------------------------- the sources
lnorm_q <- function(mean, sd) {
  ml <- log(mean^2/sqrt(sd^2 + mean^2)); sl <- sqrt(log(1 + sd^2/mean^2))
  function(u) stats::qlnorm(u, ml, sl)
}
## unified / source-A form (identical structure); one eta
mk_U <- function(p) eval(parse(text = sprintf('function() {
  ini({ lcl <- %.8f; lv <- %.8f; clwt <- %.8f; vwt <- %.8f
        eta.cl ~ %.8f; prop.err <- %.8f })
  model({ cl <- exp(lcl + eta.cl) * (WT/70)^clwt
          v  <- exp(lv) * (WT/70)^vwt
          d/dt(centr) <- -cl/v * centr
          cp <- centr/v
          cp ~ prop(prop.err) }) }',
  p[["lcl"]], p[["lv"]], p[["clwt"]], p[["vwt"]], p[["om.cl"]], p[["prop.err"]])))
## source-B form: NO volume-weight term
mk_B <- function(p) eval(parse(text = sprintf('function() {
  ini({ lcl <- %.8f; lv <- %.8f; clwt <- %.8f
        eta.cl ~ %.8f; prop.err <- %.8f })
  model({ cl <- exp(lcl + eta.cl) * (WT/70)^clwt
          v  <- exp(lv)
          d/dt(centr) <- -cl/v * centr
          cp <- centr/v
          cp ~ prop(prop.err) }) }',
  p[["lcl"]], p[["lv"]], p[["clwt"]], p[["om.cl"]], p[["prop.err"]])))

PHI <- c(lcl = log(4.2), lv = log(48), clwt = 0.75, vwt = 1.0,
         om.cl = 0.08, prop.err = 0.13)                      # A's published set
THB <- c(lcl = log(5.1), lv = log(55), clwt = 0.62,
         om.cl = 0.05, prop.err = 0.16)                      # B's published set
## published uncertainty. A: n=40 but RICH sampling -> tight exponents.
## B: n=150 but TROUGH-ONLY -> loose exponents. n and precision deliberately
## point in opposite directions; that is where the two objectives separate.
SEA <- c(lcl = 0.05, lv = 0.06, clwt = 0.05, vwt = 0.07, om.cl = 0.020, prop.err = 0.012)
SEB <- c(lcl = 0.04, lv = 0.05, clwt = 0.16,             om.cl = 0.011, prop.err = 0.009)

design <- list(
  A = list(n = 40,  times = c(0.5,1,2,4,6,8,12), ev = rxode2::et(amt = 500),
           cov_dist = list(WT = list(quantile = lnorm_q(68, 12)))),
  B = list(n = 150, times = c(11.5, 13),
           ev = rxode2::et(amt = 750, ii = 12, addl = 5),
           cov_dist = list(WT = list(quantile = lnorm_q(84, 19)))))

stratify <- function(study, K, covariate = "WT") {
  q <- study$cov_dist[[covariate]]$quantile
  stats::setNames(lapply(seq_len(K), function(k) {
    force(k); s <- study; s$n <- study$n/K
    s$cov_dist[[covariate]] <- list(quantile = function(u) q(((k-1)+u)/K)); s
  }), paste0(covariate, seq_len(K)))
}
gen_from <- function(mod, d, K) datagen(stratify(utils::modifyList(d, list(model = mod)), K),
                                        control = GEN)
fit_form <- function(startmod, studies)
  suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    startmod, admData(), est = "adgh",
    control = adghControl(studies = studies, n_nodes = 5L, grad = "analytical",
                          print = 0L, covMethod = "none"))))
grab <- function(f, nms) {
  fx <- nlmixr2est::fixef(f); om <- f$omega
  stats::setNames(vapply(nms, function(n)
    if (n == "om.cl") om[["eta.cl", "eta.cl"]] else fx[[n]], 0), nms)
}

## ============================== 1. AGGREGATE ================================
agg_blocks <- function(K) {
  a <- stats::setNames(stratify(utils::modifyList(design$A, list(model = mk_U(PHI))), K),
                       paste0("A_", seq_len(K)))
  b <- stats::setNames(stratify(utils::modifyList(design$B, list(model = mk_B(THB))), K),
                       paste0("B_", seq_len(K)))
  datagen(c(a, b), control = GEN)
}
t1 <- proc.time()[["elapsed"]]
START <- c(lcl = log(4.5), lv = log(50), clwt = 0.5, vwt = 0.5,
           om.cl = 0.07, prop.err = 0.15)
fit_agg <- fit_form(mk_U(START), agg_blocks(KK))
PN  <- names(PHI)
agg <- grab(fit_agg, PN)
say(sprintf("[%4.0fs] AGGREGATE  : %s", proc.time()[["elapsed"]] - T0,
            paste(sprintf("%s %.3f", PN, agg), collapse = "  ")))
say(sprintf("         one aggregate fit = %.0f s", proc.time()[["elapsed"]] - t1))

## ============================== 2. INDUCED MAP ==============================
## g_A: A's form == unified form, so this SHOULD be the identity. Verify.
t1 <- proc.time()[["elapsed"]]
gA_check <- grab(fit_form(mk_U(START), gen_from(mk_U(PHI), design$A, KK)), PN)
say(sprintf("[%4.0fs] g_A(phi) vs phi (identity check): max |diff| = %.2e   (%.0f s)",
            proc.time()[["elapsed"]] - T0, max(abs(gA_check - PHI)),
            proc.time()[["elapsed"]] - t1))

BN <- names(THB)
gB <- function(phi) grab(fit_form(mk_B(THB), gen_from(mk_U(phi), design$B, KK)), BN)
t1 <- proc.time()[["elapsed"]]
gB0 <- gB(PHI)
say(sprintf("[%4.0fs] g_B(phi_A)  : %s   (%.0f s per inner fit)",
            proc.time()[["elapsed"]] - T0,
            paste(sprintf("%s %.3f", BN, gB0), collapse = "  "),
            proc.time()[["elapsed"]] - t1))

## ============================== 3. PARAMETER LL =============================
## scaled residual vector: both sources stacked, each divided by its own SE
resid <- function(phi) c((PHI[PN] - phi[PN]) / SEA[PN],      # g_A = identity
                         (THB[BN] - gB(phi))  / SEB[BN])
FD <- c(lcl = 4e-3, lv = 4e-3, clwt = 8e-3, vwt = 8e-3, om.cl = 4e-3, prop.err = 3e-3)
phi <- PHI
for (it in 1:3) {
  r0 <- resid(phi)
  J  <- vapply(PN, function(p) {
    ph <- phi; ph[[p]] <- ph[[p]] + FD[[p]]
    (resid(ph) - r0) / FD[[p]] }, numeric(length(r0)))
  step <- tryCatch(solve(crossprod(J) + diag(1e-8, ncol(J)), crossprod(J, r0)),
                   error = function(e) rep(0, ncol(J)))
  phi <- phi - as.numeric(step)
  say(sprintf("[%4.0fs]   GN %d  SSR %.4f -> %s", proc.time()[["elapsed"]] - T0,
              it, sum(r0^2), paste(sprintf("%s %.3f", PN, phi), collapse = "  ")))
}
rF   <- resid(phi)
Vphi <- tryCatch(solve(crossprod(J)), error = function(e) matrix(NA, 6, 6))
Q    <- sum(rF^2); dfQ <- length(rF) - length(PN)

## ---- between-source heterogeneity (random-effects meta-analysis) ------------
## Q rejects => one common structure does not explain the published numbers to
## within their reported precision. The meta-analytic response is a between-source
## variance tau^2 added to every Sigma_s, estimated so that the weighted residual
## sum equals its degrees of freedom (Paule-Mandel). All of this reuses the
## Jacobian already computed -- no further model fits.
SEs  <- c(SEA[PN], SEB[BN])
D    <- J * SEs                    # d(raw residual)/d(phi)
rraw <- rF * SEs
solve_at <- function(tau2) {
  w <- 1/(SEs^2 + tau2)
  A <- crossprod(D, w * D) + diag(1e-10, ncol(D))
  d <- -solve(A, crossprod(D, w * rraw))
  list(d = d, Q = sum(w * (rraw + D %*% d)^2), V = solve(A))
}
tau2 <- tryCatch(stats::uniroot(function(t2) solve_at(t2)$Q - dfQ,
                                c(0, 5), tol = 1e-10)$root,
                 error = function(e) 0)
re      <- solve_at(tau2)
phi_re  <- phi + as.numeric(re$d)
se_fe   <- sqrt(pmax(diag(Vphi), 0))
se_re   <- sqrt(pmax(diag(re$V), 0))

say("\n================================ RESULT ================================")
say(sprintf("%-10s %10s %10s %10s %10s", "param", "A pub", "B pub", "AGGREGATE", "PARAM-LL"))
for (p in PN)
  say(sprintf("%-10s %10.3f %10s %10.3f %10.3f", p, PHI[[p]],
              if (p %in% BN) sprintf("%.3f", THB[[p]]) else "-", agg[[p]], phi[[p]]))
say(sprintf("\ngoodness of fit  Q = %.2f on %d df,  p = %.4f", Q, dfQ,
            stats::pchisq(Q, dfQ, lower.tail = FALSE)))
say(sprintf("between-source SD tau = %.3f\n", sqrt(tau2)))
say(sprintf("%-10s %10s %8s %10s %8s", "param", "PARAM-LL", "SE(fe)", "+ tau^2", "SE(re)"))
for (i in seq_along(PN))
  say(sprintf("%-10s %10.3f %8.3f %10.3f %8.3f", PN[i], phi[[i]], se_fe[i],
              phi_re[i], se_re[i]))
say(sprintf("\ntotal %.0f s", proc.time()[["elapsed"]] - T0))
saveRDS(list(agg = agg, phi = phi, V = Vphi, Q = Q, df = dfQ, gB0 = gB0,
             gA_check = gA_check, tau2 = tau2, phi_re = phi_re,
             se_fe = se_fe, se_re = se_re),
        "C:/Users/hidde/AppData/Local/Temp/claude/C--package-admixr2/3ff305c7-64a5-4cb0-ba61-1436e2e9b16e/scratchpad/two-lik.rds")
