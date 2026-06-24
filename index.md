# admixr2

`admixr2` fits pharmacometric PK/PD models directly to **aggregate-level
data** — the observed mean vector **E** and covariance matrix **V**
reported per clinical study — rather than requiring individual patient
records. It integrates with the [nlmixr2](https://nlmixr2.org/) /
[rxode2](https://cran.r-project.org/package=rxode2) ecosystem and
provides four estimation backends:

| Estimator | `est =` | Control |
|----|----|----|
| First-Order | `"adfo"` | [`adfoControl()`](https://leidenpharmacology.github.io/admixr2/reference/adfoControl.md) |
| Monte Carlo | `"admc"` | [`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md) |
| Gauss-Hermite | `"adgh"` | [`adghControl()`](https://leidenpharmacology.github.io/admixr2/reference/adghControl.md) |
| Iterative Reweighting MC | `"adirmc"` | [`adirmcControl()`](https://leidenpharmacology.github.io/admixr2/reference/adirmcControl.md) |

## Model-Based Meta-Analysis

**Model-Based Meta-Analysis (MBMA)** is a pharmacometric framework for
synthesising evidence across multiple clinical studies by fitting a
shared mechanistic PK/PD model to the aggregate outcomes (means,
variances) reported in each study. Unlike classical meta-analysis, which
pools effect estimates, MBMA preserves the full pharmacometric model
structure — including nonlinear dose-response, inter-individual
variability, and residual error — enabling principled extrapolation and
dose optimisation across the evidence base.

`admixr2` fits a single population model jointly to all studies.
Between-study differences in outcomes are accounted for through the
population model structure: inter-individual variability captures
subject-level spread within each study, and residual error absorbs
remaining discrepancies. Each study contributes its own dosing regimen,
observation times, and sample size, but shares the same structural and
variance parameters.

## When to use admixr2

**Individual patient data are unavailable** — the most common scenario
in MBMA. Published papers report means and standard deviations;
regulatory submissions and competitive reasons prevent individual
patient data sharing across companies or institutions. `admixr2`
extracts the maximum information from what is publicly available.

**Leveraging the literature for trial design** — fit a mechanistic PK/PD
model to aggregated results from existing trials, then simulate new
dosing regimens or patient populations before committing to a costly
study.

**Combining evidence across heterogeneous trials** — studies differ in
dose, formulation, population, or observation schedule. `admixr2`
handles multi-study fits with per-study dosing events and time grids
under a single shared population model.

**Reproducing and extending published models** — digitised mean
concentration– time profiles from figures are sufficient input. No
individual patient data required.

## Installation

``` r

# From CRAN
install.packages("admixr2")

# Development version from GitHub
pak::pak("LeidenPharmacology/admixr2")
```

## Quick start

``` r

library(admixr2)
library(rxode2)
library(nlmixr2)

# 1. Compute aggregate statistics from individual data (or digitise from paper)
data("examplomycin")
obs    <- examplomycin[examplomycin$EVID == 0, ]
obs    <- obs[order(obs$ID, obs$TIME), ]
times  <- sort(unique(obs$TIME))
ids    <- unique(obs$ID)
dv_mat <- matrix(NA_real_, nrow = length(ids), ncol = length(times))
for (i in seq_along(ids)) {
  sub         <- obs[obs$ID == ids[i], ]
  dv_mat[i, ] <- sub$DV[order(sub$TIME)]
}
E <- colMeans(dv_mat)
V <- cov.wt(dv_mat, method = "ML")$cov

# 2. Define the model (standard nlmixr2 syntax)
pk_model <- function() {
  ini({
    tcl     <- log(5);  label("Log clearance (L/hr)")
    tv1     <- log(10); label("Log central volume (L)")
    tv2     <- log(30); label("Log peripheral volume (L)")
    tq      <- log(10); label("Log inter-compartmental CL (L/hr)")
    tka     <- log(1);  label("Log absorption rate constant (1/hr)")
    prop.sd <- c(0, 0.2)
    eta.cl ~ 0.09; eta.v1 ~ 0.09; eta.v2 ~ 0.09
    eta.q  ~ 0.09; eta.ka ~ 0.09
  })
  model({
    cl <- exp(tcl + eta.cl); v1 <- exp(tv1 + eta.v1)
    v2 <- exp(tv2 + eta.v2); q  <- exp(tq  + eta.q)
    ka <- exp(tka + eta.ka)
    d/dt(depot)      <- -ka * depot
    d/dt(central)    <- ka * depot - (cl/v1 + q/v1) * central + (q/v2) * peripheral
    d/dt(peripheral) <- (q/v1) * central - (q/v2) * peripheral
    cp <- central / v1
    cp ~ prop(prop.sd)
  })
}

# 3. Fit
fit <- nlmixr2(
  pk_model, admData(), est = "admc",
  control = admControl(
    studies = list(examplomycin = list(
      E = E, V = V, n = length(ids),
      times = times, ev = et(amt = 100)
    )),
    n_sim = 5000L, seed = 1L
  )
)

print(fit)
plot(fit)
```

## Vignettes

| Vignette | Topic |
|----|----|
| [Getting started](https://leidenpharmacology.github.io/admixr2/articles/admixr2.html) | Core workflow: data prep, model, fit, diagnostics |
| [Diagnostic plots](https://leidenpharmacology.github.io/admixr2/articles/diagnostic-plots.html) | All four plot panels explained; IIV heatmap |
| [Multiple studies](https://leidenpharmacology.github.io/admixr2/articles/multiple-studies.html) | Joint fitting across studies with different designs |
| [Estimator comparison](https://leidenpharmacology.github.io/admixr2/articles/estimator-comparison.html) | adfo, admc, adgh and adirmc: mathematical foundations and when to use each |
| [Advanced usage](https://leidenpharmacology.github.io/admixr2/articles/advanced.html) | Gradient modes, parallel restarts, AIC/BIC model comparison |

## Citation

If you use `admixr2` in your work, please cite the software paper, which
introduces the Iterative Reweighting Monte Carlo estimator:

> van de Beek H., Välitalo P.A.J., van Hasselt J.G.C., Zwep L.B. (2025).
> Aggregate data modelling: A fast implementation for fitting
> pharmacometrics models to summary-level data in R. *Journal of
> Pharmacokinetics and Pharmacodynamics*, 53(1), 3.
> <https://doi.org/10.1007/s10928-025-10011-w>

The aggregate data modelling methodology is introduced in:

> Välitalo P.A.J. (2021). Pharmacometric estimation methods for
> aggregate data, including data simulated from other pharmacometric
> models. *Journal of Pharmacokinetics and Pharmacodynamics*, 48(5),
> 623–638. <https://doi.org/10.1007/s10928-021-09760-1>
