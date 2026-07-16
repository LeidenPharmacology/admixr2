# Estimation-method classification tags.  Each aggregate-data estimator is
# tagged with a `type` (its category) and a `description` so it joins the
# category-grouped estimation-method list nlmixr2est prints for an unsupported
# `est=` (or a bare `nlmixr2()` call).  All admixr2 estimators are Model Based
# Meta Analysis (aggregate-data) methods.  Sourced after the ad*.R method
# definitions (alphabetical collation), so the functions already exist here.

attr(nlmixr2Est.adfo, "type") <- "Model Based Meta Analysis"
attr(nlmixr2Est.adfo, "description") <- "Aggregate data, First-Order approximation"

attr(nlmixr2Est.adgh, "type") <- "Model Based Meta Analysis"
attr(nlmixr2Est.adgh, "description") <- "Aggregate data, Gauss-Hermite quadrature"

attr(nlmixr2Est.adirmc, "type") <- "Model Based Meta Analysis"
attr(nlmixr2Est.adirmc, "description") <- "Aggregate data, Iterative Reweighting Monte Carlo"

attr(nlmixr2Est.admc, "type") <- "Model Based Meta Analysis"
attr(nlmixr2Est.admc, "description") <- "Aggregate data, Monte Carlo"
