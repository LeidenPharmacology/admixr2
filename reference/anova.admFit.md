# Compare nested admixr2 fits by a likelihood-ratio test

The ordinary LRT: the objective difference against a chi-squared
reference with `Df` equal to the number of parameters the larger model
adds.

## Usage

``` r
# S3 method for class 'admFit'
anova(object, ...)
```

## Arguments

- object:

  An `admFit`.

- ...:

  Further `admFit`s to compare it with.

## Value

A data frame of class `anova.admFit`, smallest model first.

## Details

Both fits must come from the same estimator and, for the quadrature
estimators, the same node count (`n_nodes`) or, for the Monte Carlo ones
(`admc`, `adirmc`), the same sample size (`n_sim`). Each scores its own
approximation to the likelihood, so objectives from different ones are
not comparable and the comparison is refused rather than reported.

The **stratum** resolution is checked per covariate, and only where the
two fits overlap. A covariate a model does not read cannot put the two
on different scales: if its prediction does not move across a source's
nodes, the mixture those nodes collapse to is a sufficient statistic for
it, so its objective is the same at either resolution. That is what
makes the nested pair of a covariate test comparable — the null model
drops the term, so its sources are not cut along it. Two fits that both
read a covariate and cut it differently, including one cutting it and
the other integrating over it whole, are refused.

Testing a variance AT ZERO puts the null on the boundary of the
parameter space, where the exact reference is a chi-bar-squared mixture
rather than a chi-squared. The p-value reported there is CONSERVATIVE –
too large – so a significant result stays significant, but treat a
borderline one with care.
