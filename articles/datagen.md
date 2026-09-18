# Generating aggregate data with datagen

## When would you use `datagen`?

Two situations call for generating aggregate data programmatically:

**Simulation studies.** To check that an estimator recovers the truth:
define a data-generating model with known parameters, call
[`datagen()`](https://leidenpharmacology.github.io/admixr2/reference/datagen.md)
for E and V, fit your analysis model, compare.

**Published models as input.** Each published study was analysed with
its own structural model — a different number of compartments, a
different parameterisation, only some of the IIV terms.
[`datagen()`](https://leidenpharmacology.github.io/admixr2/reference/datagen.md)
turns such a model into the aggregate data it implies, making it a
direct input alongside digitised summaries. Several of them in
`admControl(studies = ...)` is a meta-analysis across the literature
(see
[`vignette("multiple-studies")`](https://leidenpharmacology.github.io/admixr2/articles/multiple-studies.md)).

Either way the output is a named list of `(E, V, n, times, ev)` objects
that plugs straight into `admControl(studies = ...)`.

## Simulation study: same model across studies

A one-compartment oral model generates two synthetic studies — a
low-dose and a high-dose arm — which the estimator should recover the
truth from.

### Data-generating model

``` r

library(admixr2)
library(rxode2)
library(nlmixr2)
library(ggplot2)

true_model <- function() {
  ini({
    tcl     <- log(5)  ; label("Log clearance (L/h)")
    tv      <- log(10) ; label("Log volume (L)")
    tka     <- log(1)  ; label("Log absorption rate (1/h)")
    prop.sd <- c(0, 0.2) ; label("Proportional error SD")
    eta.cl  ~ 0.09
    eta.v   ~ 0.04
    eta.ka  ~ 0.04
  })
  model({
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    ka <- exp(tka + eta.ka)
    d/dt(depot)   <- -ka * depot
    d/dt(central) <- ka * depot - (cl/v) * central
    cp <- central / v
    cp ~ prop(prop.sd)
  })
}
```

True parameter values: CL = 5 L/h, V = 10 L, ka = 1 h⁻¹; IIV of 0.3 (SD
on log scale) for CL and 0.2 for V and ka; proportional residual error
SD = 0.2.

### Generating the data

A top-level `model` is the default every study inherits, so each spec
needs only `times`, `ev` and `n`:

``` r

times <- c(0.5, 1, 2, 4, 8, 12, 24)

study_data <- datagen(
  studies = list(
    low_dose  = list(times = times, ev = rxode2::et(amt =  50), n = 250L),
    high_dose = list(times = times, ev = rxode2::et(amt = 100), n = 250L)
  ),
  model   = true_model,
  control = datagenControl(n_sim = 10000L, seed = 1L)
)

# Each study returns E, V, n, times, ev
names(study_data$low_dose)
#> [1] "E"        "V"        "n"        "times"    "ev"       "v_denom"  ".adm_src"
#> [8] "output"
round(study_data$low_dose$E, 2)   # population mean at each time
#>  0.5    1    2    4    8   12   24 
#> 1.75 2.38 2.27 1.18 0.24 0.06 0.00
```

`V` carries the IIV-driven between-time correlation and, on the
diagonal, the residual error:

``` r

knitr::kable(round(study_data$low_dose$V, 3),
             caption = "Generated covariance matrix V for the low-dose study")
```

|     |    0.5 |      1 |     2 |      4 |      8 |     12 |    24 |
|:----|-------:|-------:|------:|-------:|-------:|-------:|------:|
| 0.5 |  0.297 |  0.187 | 0.106 | -0.005 | -0.018 | -0.007 | 0.000 |
| 1   |  0.187 |  0.461 | 0.171 |  0.051 |  0.001 | -0.001 | 0.000 |
| 2   |  0.106 |  0.171 | 0.428 |  0.161 |  0.051 |  0.015 | 0.001 |
| 4   | -0.005 |  0.051 | 0.161 |  0.263 |  0.084 |  0.028 | 0.001 |
| 8   | -0.018 |  0.001 | 0.051 |  0.084 |  0.047 |  0.016 | 0.001 |
| 12  | -0.007 | -0.001 | 0.015 |  0.028 |  0.016 |  0.007 | 0.000 |
| 24  |  0.000 |  0.000 | 0.001 |  0.001 |  0.001 |  0.000 | 0.000 |

Generated covariance matrix V for the low-dose study {.table}

### Inspecting the generated profiles

``` r

df_gen <- do.call(rbind, lapply(names(study_data), function(nm) {
  s  <- study_data[[nm]]
  sd <- sqrt(diag(s$V))
  data.frame(study = nm, time = s$times,
             mean = s$E, lo = s$E - sd, hi = s$E + sd)
}))

ggplot(df_gen, aes(x = time, y = mean, colour = study, fill = study)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_x_log10(breaks = times, labels = times) +
  scale_colour_manual(values = c(low_dose = "#0072B2", high_dose = "#D55E00")) +
  scale_fill_manual(  values = c(low_dose = "#0072B2", high_dose = "#D55E00")) +
  labs(title  = "Generated aggregate data — mean ± 1 SD",
       x = "Time (h, log scale)", y = "Concentration (mg/L)",
       colour = NULL, fill = NULL) +
  theme_bw()
```

![Generated population mean ± 1 SD for both dose
levels.](datagen_files/figure-html/plot-generated-1.png)

Generated population mean ± 1 SD for both dose levels.

The ribbon is the square root of the diagonal of V — IIV spread plus
residual error. The off-diagonal within-subject correlation is not
drawn, but it reaches the estimator and enters the likelihood.

### Parameter recovery

The generated data plug straight into `admControl(studies = ...)`. Using
the same structural model for the analysis is what makes this a
consistency check:

``` r

analysis_model <- function() {
  ini({
    tcl     <- log(4)  ; label("Log clearance (L/h)")
    tv      <- log(12) ; label("Log volume (L)")
    tka     <- log(1.5); label("Log absorption rate (1/h)")
    prop.sd <- c(0, 0.3); label("Proportional error SD")
    eta.cl  ~ 0.12
    eta.v   ~ 0.05
    eta.ka  ~ 0.05
  })
  model({
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    ka <- exp(tka + eta.ka)
    d/dt(depot)   <- -ka * depot
    d/dt(central) <- ka * depot - (cl/v) * central
    cp <- central / v
    cp ~ prop(prop.sd)
  })
}
```

The starting values are deliberately off the truth (CL 4 against 5, V 12
against 10), so finding it is the estimator’s doing rather than the
initialisation’s.

``` r
fit_sim <- nlmixr2(
  analysis_model, admData(), est = "admc",
  control = admControl(
    studies   = study_data,
      maxeval = 300L
  )
)

print(fit_sim)
── nlmixr² admc ──

          OBJF       AIC       BIC Log-likelihood
admc -7359.742 -7345.742 -7302.618       3679.871

── Time (sec fit_sim$time): ──

  optimize covariance other elapsed
1   46.325          0     0  46.325

── Population Parameters (fit_sim$parFixed or fit_sim$parFixedDf): ──

                        Parameter     Est. SE %RSE Back-transformed(95%CI)
tcl           Log clearance (L/h)    1.608                           4.991
tv                 Log volume (L)    2.303                           10.00
tka     Log absorption rate (1/h) 0.001072                           1.001
prop.sd     Proportional error SD   0.2001                          0.2001
        BSV(CV%) Shrink(SD)%
tcl        30.41         NaN
tv         19.60         NaN
tka        20.88         NaN
prop.sd                     
 
  No correlations in between subject variability (BSV) matrix
  Full BSV covariance (fit_sim$omega) 
    or correlation (fit_sim$omegaR; diagonals=SDs)
  Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
  Censoring (fit_sim$censInformation): No censoring
  Minimization message (fit_sim$message):  
    NLOPT_XTOL_REACHED: Optimization stopped because xtol_rel or xtol_abs (above) was reached. 
```

**No standard errors, and the blank `SE` / `%RSE` / `95%CI` columns are
why.** A study generated from a model is not a sample: `n` is the source
study’s true sample size, but the parameter uncertainty a sampling law
would need is not available, so admixr2 sets `covMethod = "none"` and
reports no uncertainty at all. Naming a `covMethod` explicitly is an
error rather than an override. The columns are nlmixr2’s table; the
emptiness is the answer.

Structural parameter estimates and the truth:

``` r

est   <- fit_sim$env$admExtra$struct
truth <- c(tcl = log(5), tv = log(10), tka = log(1))

knitr::kable(
  data.frame(
    parameter = names(truth),
    truth     = round(exp(truth), 2),
    estimate  = round(exp(est[names(truth)]), 2)
  ),
  row.names = FALSE, caption = "Structural parameter estimates vs true values")
```

| parameter | truth | estimate |
|:----------|------:|---------:|
| tcl       |     5 |     4.99 |
| tv        |    10 |    10.00 |
| tka       |     1 |     1.00 |

Structural parameter estimates vs true values {.table}

## Literature workflow: per-study models

Each publication comes with its own model — a 2016 paper with one
compartment, a 2022 follow-up with two. A `model` inside a study spec
overrides the top-level default:

``` r

# One-compartment oral model (as published in a simpler earlier study)
model_1cmt <- function() {
  ini({
    tcl     <- log(5)
    tv      <- log(40)    # apparent volume, peripheral compartment lumped in
    tka     <- log(1)
    prop.sd <- c(0, 0.25)
    eta.cl  ~ 0.04
    eta.v   ~ 0.01
    eta.ka  ~ 0.01
  })
  model({
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    ka <- exp(tka + eta.ka)
    d/dt(depot)   <- -ka * depot
    d/dt(central) <- ka * depot - (cl/v) * central
    cp <- central / v
    cp ~ prop(prop.sd)
  })
}

# Two-compartment oral model (as published in a more detailed later study)
model_2cmt <- function() {
  ini({
    tcl     <- log(5)
    tv1     <- log(10)
    tv2     <- log(30)
    tq      <- log(10)
    tka     <- log(1)
    prop.sd <- c(0, 0.2)
    eta.cl  ~ 0.09
    eta.v1  ~ 0.04
  })
  model({
    cl <- exp(tcl + eta.cl)
    v1 <- exp(tv1 + eta.v1)
    v2 <- exp(tv2)
    q  <- exp(tq)
    ka <- exp(tka)
    d/dt(depot)      <- -ka * depot
    d/dt(central)    <- ka * depot - (cl/v1 + q/v1) * central + (q/v2) * peripheral
    d/dt(peripheral) <- (q/v1) * central - (q/v2) * peripheral
    cp <- central / v1
    cp ~ prop(prop.sd)
  })
}
```

Both share a CL of 5 L/h and differ only in how they describe
distribution: the earlier study lumped the periphery into one larger
apparent volume (V = 40 L), the later one resolved it (V₁ = 10 L, V₂ =
30 L, Q = 10 L/h).

``` r

lit_data <- datagen(
  studies = list(
    early_study = list(
      model = model_1cmt,
      times = c(0.5, 1, 2, 4, 8, 12),
      ev    = rxode2::et(amt = 100),
      n     = 120L
    ),
    later_study = list(
      model = model_2cmt,
      times = c(2, 4, 8, 12, 24),
      ev    = rxode2::et(amt = 200),
      n     = 180L
    )
  ),
  control = datagenControl(n_sim = 10000L, seed = 2L)
)
```

Same drug, same CL, different E profiles — the dosing and the
distributional assumptions differ:

``` r

df_lit <- do.call(rbind, lapply(names(lit_data), function(nm) {
  s  <- lit_data[[nm]]
  sd <- sqrt(diag(s$V))
  data.frame(study = nm, time = s$times,
             mean = s$E, lo = s$E - sd, hi = s$E + sd)
}))

ggplot(df_lit, aes(x = time, y = mean, colour = study, fill = study)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c(early_study = "#009E73", later_study = "#CC79A7")) +
  scale_fill_manual(  values = c(early_study = "#009E73", later_study = "#CC79A7")) +
  labs(title  = "Literature-derived aggregate profiles",
       x = "Time (h)", y = "Concentration (mg/L)",
       colour = NULL, fill = NULL) +
  theme_bw()
```

![Simulated aggregate profiles from two publications with different
structural models.](datagen_files/figure-html/lit-plot-1.png)

Simulated aggregate profiles from two publications with different
structural models.

### Fitting an analysis model to literature data

The data-generating models encode each publication’s assumptions; the
analysis model is your unified structural hypothesis, fitted to both
studies at once:

``` r

lit_analysis <- function() {
  ini({
    tcl     <- log(5)
    tv1     <- log(10)
    tv2     <- log(30)
    tq      <- log(10)
    tka     <- log(1)
    prop.sd <- c(0, 0.25)
    eta.cl  ~ 0.09
    eta.v1  ~ 0.04
  })
  model({
    cl <- exp(tcl + eta.cl)
    v1 <- exp(tv1 + eta.v1)
    v2 <- exp(tv2)
    q  <- exp(tq)
    ka <- exp(tka)
    d/dt(depot)      <- -ka * depot
    d/dt(central)    <- ka * depot - (cl/v1 + q/v1) * central + (q/v2) * peripheral
    d/dt(peripheral) <- (q/v1) * central - (q/v2) * peripheral
    cp <- central / v1
    cp ~ prop(prop.sd)
  })
}

fit_lit <- nlmixr2(
  lit_analysis, admData(), est = "admc",
  control = admControl(
    studies   = lit_data,
    n_sim     = 5000L,
    cov_n_sim = 10000L,
    seed      = 3L
  )
)
```

Because the two used genuinely different structural models, the analysis
model lands on a weighted compromise between the two parameterisations,
minimising the joint discrepancy rather than matching either exactly.
That is the point: aggregate-data modelling takes each study as
published and looks for one set of parameters consistent with all of
them.

## Examining individual simulated samples

`return_samples = TRUE` adds the raw `n_sim × n_times` prediction matrix
to each study’s output — useful for inspecting the shape of the
simulated distribution, or for custom summaries:

``` r

study_with_samples <- datagen(
  studies = list(
    single_study = list(times = times, ev = rxode2::et(amt = 100), n = 250L)
  ),
  model   = true_model,
  control = datagenControl(n_sim = 2000L, seed = 4L, return_samples = TRUE)
)

cp_mat <- study_with_samples$single_study$samples   # 2000 x 7
mu     <- study_with_samples$single_study$E

df_samp <- data.frame(
  time  = rep(times, each = 200),
  conc  = as.vector(cp_mat[1:200, ]),
  id    = rep(seq_len(200), times = length(times))
)

ggplot(df_samp, aes(x = time, y = conc, group = id)) +
  geom_line(alpha = 0.06, colour = "grey30") +
  geom_line(data = data.frame(time = times, conc = mu),
            aes(group = NULL), colour = "#0072B2", linewidth = 1.5) +
  labs(title = "Simulated individual trajectories (n = 200 shown)",
       x = "Time (h)", y = "Concentration (mg/L)") +
  theme_bw()
```

![First 200 individual simulated trajectories (grey) with the population
mean (blue).](datagen_files/figure-html/samples-1.png)

First 200 individual simulated trajectories (grey) with the population
mean (blue).

Blue is the population mean E, the grey spaghetti the IIV spread.
`samples` holds concentrations **before** residual error; only the
diagonal of V carries that.

## Choosing the moment method

[`datagen()`](https://leidenpharmacology.github.io/admixr2/reference/datagen.md)
integrates over the IIV by Monte Carlo, the same engine as
`est = "admc"`. Two deterministic alternatives match the other
estimators, and both ignore `n_sim`, `sampling` and `seed`, so their
output is exactly reproducible:

``` r

one <- list(single_study = list(times = times, ev = rxode2::et(amt = 100),
                                n = 250L))
fo_data <- datagen(studies = one, model = true_model,
                   control = datagenControl(method = "fo"))
gh_data <- datagen(studies = one, model = true_model,
                   control = datagenControl(method = "gh", n_nodes = 5L))
round(rbind(fo = fo_data$single_study$E, gh = gh_data$single_study$E), 3)
#>      0.5     1     2     4     8    12    24
#> fo 3.445 4.773 4.651 2.340 0.360 0.049 0.000
#> gh 3.495 4.753 4.549 2.365 0.485 0.113 0.003
```

`method = "fo"` is the first-order expansion `est = "adfo"` uses,
`E = f(theta, 0)` and `V = J Omega J' + Sigma`; `method = "gh"` is the
Gauss-Hermite quadrature `est = "adgh"` uses, unbiased at any IIV.
`n_nodes` (per eta dimension) trades accuracy against cost: 3 is fast, 5
is near-exact to IIV SD ~0.5, 7 reaches ~0.7.

For design work, **generate and analyse under the same approximation**.
The expected information matrix is the Hessian of the log-likelihood at
the generating parameters, which is only a genuine maximum when the two
agree – generate with Monte Carlo, analyse under FO, and those
parameters are not an FO maximum-likelihood estimate, so the FIM is
biased.

## See also

- [Multiple
  studies](https://leidenpharmacology.github.io/admixr2/articles/multiple-studies.md)
  — feed per-study models into a meta-analysis
- [From a published figure to E, V and
  n](https://leidenpharmacology.github.io/admixr2/articles/aggregate-data.md)
  — the other input type
- [Estimator
  comparison](https://leidenpharmacology.github.io/admixr2/articles/estimator-comparison.md)
