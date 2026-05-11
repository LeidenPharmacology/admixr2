# Estimator comparison: admc vs adirmc

## Two estimators, one interface

Both estimators target the same aggregate-data log-likelihood and accept
the same model and study specification. They differ only in how they
approximate the population distribution:

**`admc` — Monte Carlo:** Draws `n_sim` eta samples from the current
Omega, simulates each trajectory via rxSolve, and forms a sample
covariance to enter the NLL. The analytical gradient (sensitivity
equations or CRN finite-differences) makes individual NLL evaluations
fast. Works well for standard PK models with good starting values.

**`adirmc` — Importance Resampling MC:** Draws proposals from an
inflated Omega, then reweights them by their likelihood under the
current parameters (importance sampling). The inner optimisation given
fixed proposals is deterministic and fast; proposals are refreshed at
each outer phase. IRMC is more robust to importance weight degeneracy in
higher-dimensional models and can recover from poor starting values more
reliably.

## Common setup

``` r

library(admixr2)
library(rxode2)
library(nlmixr2)
library(ggplot2)

data("examplomycin")
obs   <- examplomycin[examplomycin$EVID == 0, ]
obs   <- obs[order(obs$ID, obs$TIME), ]
times <- sort(unique(obs$TIME))
ids   <- unique(obs$ID)
n     <- length(ids)

dv_mat <- matrix(NA_real_, nrow = n, ncol = length(times))
for (i in seq_along(ids)) {
  sub         <- obs[obs$ID == ids[i], ]
  dv_mat[i, ] <- sub$DV[order(sub$TIME)]
}
E <- colMeans(dv_mat)
V <- cov(dv_mat)

pk_model <- function() {
  ini({
    tcl     <- log(5)  ; label("Log clearance (L/hr)")
    tv1     <- log(10) ; label("Log central volume (L)")
    tv2     <- log(30) ; label("Log peripheral volume (L)")
    tq      <- log(10) ; label("Log inter-compartmental CL (L/hr)")
    tka     <- log(1)  ; label("Log absorption rate constant (1/hr)")
    prop.sd <- c(0, 0.2); label("Proportional residual error SD")
    eta.cl ~ 0.09
    eta.v1 ~ 0.09
    eta.v2 ~ 0.09
    eta.q  ~ 0.09
    eta.ka ~ 0.09
  })
  model({
    cl <- exp(tcl + eta.cl)
    v1 <- exp(tv1 + eta.v1)
    v2 <- exp(tv2 + eta.v2)
    q  <- exp(tq  + eta.q)
    ka <- exp(tka + eta.ka)
    d/dt(depot)      <- -ka * depot
    d/dt(central)    <- ka * depot - (cl/v1 + q/v1) * central + (q/v2) * peripheral
    d/dt(peripheral) <- (q/v1) * central - (q/v2) * peripheral
    cp <- central / v1
    cp ~ prop(prop.sd)
  })
}

study <- list(E = E, V = V, n = n, times = times, ev = rxode2::et(amt = 100))
```

The model uses mu-referenced parameterisation
(`cl <- exp(tcl + eta.cl)`), which enables analytical gradient
computation via sensitivity equations in both estimators. See the
[Advanced
usage](https://leidenpharmacology.github.io/admixr2/articles/advanced.html#mu-referencing-and-sensitivity-equations)
vignette for details.

## Fitting with admc

``` r

fit_mc <- nlmixr2(
  pk_model, admData(), est = "admc",
  control = admControl(
    studies   = list(examplomycin = study),
    n_sim     = 5000L,
    cov_n_sim = 10000L,
    maxeval   = 300L,
    seed      = 1L
  )
)
#> [====|====|====|====|====|====|====|====|====|====] 0:00:08
```

## Fitting with adirmc

`adirmc` is slower per evaluation but more robust to poor starting
values and high-dimensional Omega. Run it interactively with larger
`n_sim` for a production fit; the settings below are tuned for a fast
vignette build.

``` r

fit_irmc <- nlmixr2(
  pk_model, admData(), est = "adirmc",
  control = adirmcControl(
    studies         = list(examplomycin = study),
    n_sim           = 2000L,
    phases          = c(2, 1, 0.5, 0.1),
    cov_n_sim       = 10000L,
    omega_expansion = 1.5,
    seed            = 1L
  )
)
```

## Comparing parameter estimates

Once both fits are available, compare parameters in a table. Both
estimators should recover values close to the true parameters (CL=5,
V1=10, V2=30, Q=10, ka=1):

``` r

get_pars <- function(fit) {
  s  <- fit$env$admExtra$struct
  om <- diag(fit$env$admExtra$omega)
  sg <- sqrt(fit$env$admExtra$sigma_var)
  c(exp(s), om, sg)
}

tbl <- data.frame(
  Parameter = c(
    paste0("exp(", names(fit_mc$env$admExtra$struct), ")"),
    paste0("V(", fit_mc$env$admExtra$eta_col_names, ")"),
    "prop.sd"
  ),
  True  = c(5, 10, 30, 10, 1, rep(0.09, 5), 0.2),
  admc  = round(get_pars(fit_mc),   4),
  adirmc = round(get_pars(fit_irmc), 4)
)

knitr::kable(tbl, caption = "Parameter estimates vs true values")
```

## Comparing objectives

Both estimators target the same likelihood; comparable -2LL values
confirm both have converged:

``` r

cat(sprintf("admc  -2LL = %.2f   AIC = %.2f\n",
            fit_mc$objective,   AIC(fit_mc)))
cat(sprintf("adirmc -2LL = %.2f  AIC = %.2f\n",
            fit_irmc$objective, AIC(fit_irmc)))
```

A large gap in the final objective suggests convergence failure in one
of the estimators — try more restarts or a better starting point.

## NLL convergence traces

``` r

plot(fit_mc, which = "nll")
```

![admc NLL convergence
trace.](estimator-comparison_files/figure-html/nll-trace-1.png)

admc NLL convergence trace.

For `adirmc`, the NLL trace reflects outer phase iterations rather than
individual LBFGS steps. Jumps between phases (as box constraints
tighten) are normal and expected.

## When to use each

| Situation | Recommendation |
|----|----|
| Standard 1–2 compartment PK, good starting values | `admc` with `grad = "sens"` |
| Initial exploration or poor starting values | `admc` with `n_restarts >= 3` |
| High-dimensional Omega (≥ 5 etas) | `adirmc` — IS weights degrade more slowly |
| Non-Gaussian or bounded IIV | `adirmc` with `omega_expansion > 1` |
| Fastest possible runtime | `admc` with `grad = "none"` (BOBYQA) |
| Maximum reproducibility | Both estimators accept `seed` |

## Control objects at a glance

``` r

# admc key arguments
admControl(
  studies    = list(...),
  n_sim      = 5000L,       # MC sample count
  grad       = "sens",      # gradient mode: "sens", "fd", "cfd", "none"
  n_restarts = 1L,          # number of optimizer restarts
  workers    = 1L,          # parallel workers for restarts
  covMethod  = "r",         # "r" = numerical Hessian, "none" = skip
  seed       = 1L
)

# adirmc key arguments
adirmcControl(
  studies         = list(...),
  n_sim           = 5000L,
  phases          = c(2, 1, 0.5, 0.1),   # box constraint half-widths per phase
  omega_expansion = 1.5,                   # inflate proposal Omega
  grad            = "analytical",          # "analytical", "none", "fd"
  kappa_method    = "first-order",         # "first-order" or "second-order"
  seed            = 1L
)
```
