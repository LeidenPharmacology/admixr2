# Package index

## Fitting

Main entry point and data helper.

- [`admData()`](https://leidenpharmacology.github.io/admixr2/reference/admData.md)
  : Dummy data frame for nlmixr2 dispatch
- [`nlmixr2Est(`*`<adfo>`*`)`](https://leidenpharmacology.github.io/admixr2/reference/nlmixr2Est.adfo.md)
  : Fit an aggregate data model via First-Order (FO) approximation
- [`nlmixr2Est(`*`<admc>`*`)`](https://leidenpharmacology.github.io/admixr2/reference/nlmixr2Est.admc.md)
  : Fit an aggregate data model via Monte Carlo (admc estimator)
- [`nlmixr2Est(`*`<adgh>`*`)`](https://leidenpharmacology.github.io/admixr2/reference/nlmixr2Est.adgh.md)
  : Fit an aggregate data model via Gauss-Hermite quadrature
- [`nlmixr2Est(`*`<adirmc>`*`)`](https://leidenpharmacology.github.io/admixr2/reference/nlmixr2Est.adirmc.md)
  : Fit an aggregate data model via Iterative Reweighting MC (adirmc
  estimator)

## Writing down the studies

Describe a source the way its paper does — the model it published, the
cohort it enrolled — and let admixr2 assemble the summary it implies.

- [`admStudy()`](https://leidenpharmacology.github.io/admixr2/reference/admStudy.md)
  : Write down a published study
- [`admStudies()`](https://leidenpharmacology.github.io/admixr2/reference/admStudies.md)
  : Collect the studies a meta-analysis draws on
- [`admPopulation()`](https://leidenpharmacology.github.io/admixr2/reference/admPopulation.md)
  : Describe the population a study enrolled

## Covariates

Declare the covariate distribution a source enrolled, and how it enters
the fit.

- [`covDist()`](https://leidenpharmacology.github.io/admixr2/reference/covDist.md)
  : Describe the covariate distribution a study's subjects span

- [`covStrata()`](https://leidenpharmacology.github.io/admixr2/reference/covStrata.md)
  : Cut a covariate distribution into strata

- [`covDraw()`](https://leidenpharmacology.github.io/admixr2/reference/covDraw.md)
  :

  Draw covariate values from a `cov_dist` specification

## Simulation

Generate aggregate data from models for simulation studies.

- [`datagen()`](https://leidenpharmacology.github.io/admixr2/reference/datagen.md)
  : Generate aggregate study data from (possibly different)
  pharmacometric models

- [`datagenControl()`](https://leidenpharmacology.github.io/admixr2/reference/datagenControl.md)
  :

  Control parameters for
  [`datagen()`](https://leidenpharmacology.github.io/admixr2/reference/datagen.md)

## Control objects

Configuration for each estimator.

- [`adfoControl()`](https://leidenpharmacology.github.io/admixr2/reference/adfoControl.md)
  : Control settings for the FO (First-Order) estimator
- [`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md)
  : Control settings for the ADM estimator
- [`adghControl()`](https://leidenpharmacology.github.io/admixr2/reference/adghControl.md)
  : Control settings for the Gauss-Hermite (GH) quadrature estimator
- [`adirmcControl()`](https://leidenpharmacology.github.io/admixr2/reference/adirmcControl.md)
  : Control settings for the IRMC estimator
- [`admStopWorkers()`](https://leidenpharmacology.github.io/admixr2/reference/admStopWorkers.md)
  : Stop parallel workers

## Methods

S3 methods for fit objects.

- [`print(`*`<admFit>`*`)`](https://leidenpharmacology.github.io/admixr2/reference/print.admFit.md)
  : Print method for admFit objects
- [`plot(`*`<admFit>`*`)`](https://leidenpharmacology.github.io/admixr2/reference/plot.admFit.md)
  : Diagnostic plots for an admixr2 fit
- [`anova(`*`<admFit>`*`)`](https://leidenpharmacology.github.io/admixr2/reference/anova.admFit.md)
  : Compare nested admixr2 fits by a likelihood-ratio test

## Data

Example datasets.

- [`examplomycin`](https://leidenpharmacology.github.io/admixr2/reference/examplomycin.md)
  : Examplomycin dataset
