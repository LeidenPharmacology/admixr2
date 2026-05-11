# Package index

## Fitting

Main entry point and data helper.

- [`admData()`](https://leidenpharmacology.github.io/admixr2/reference/admData.md)
  : Dummy data frame for nlmixr2 dispatch
- [`nlmixr2Est(`*`<admc>`*`)`](https://leidenpharmacology.github.io/admixr2/reference/nlmixr2Est.admc.md)
  : Fit an aggregate data model via Monte Carlo (admc estimator)
- [`nlmixr2Est(`*`<adirmc>`*`)`](https://leidenpharmacology.github.io/admixr2/reference/nlmixr2Est.adirmc.md)
  : Fit an aggregate data model via Importance Resampling MC (adirmc
  estimator)

## Control objects

Configuration for each estimator.

- [`admControl()`](https://leidenpharmacology.github.io/admixr2/reference/admControl.md)
  : Control settings for the ADM estimator
- [`adirmcControl()`](https://leidenpharmacology.github.io/admixr2/reference/adirmcControl.md)
  : Control settings for the IRMC estimator
- [`admStopWorkers()`](https://leidenpharmacology.github.io/admixr2/reference/admStopWorkers.md)
  : Stop PSOCK workers

## Methods

S3 methods for fit objects.

- [`print(`*`<admFit>`*`)`](https://leidenpharmacology.github.io/admixr2/reference/print.admFit.md)
  : Print method for admFit objects
- [`plot(`*`<admFit>`*`)`](https://leidenpharmacology.github.io/admixr2/reference/plot.admFit.md)
  : Diagnostic plots for an admixr2 fit

## Data

Example datasets.

- [`examplomycin`](https://leidenpharmacology.github.io/admixr2/reference/examplomycin.md)
  : Examplomycin dataset
