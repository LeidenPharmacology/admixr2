# `type`/`description` tags so each estimator joins nlmixr2est's category-grouped
# estimation-method list for an unsupported `est=`. Sourced after the ad*.R
# method definitions (alphabetical collation), so the functions already exist.

attr(nlmixr2Est.adfo, "type") <- "Model Based Meta Analysis"
attr(nlmixr2Est.adfo, "description") <- "Aggregate data, First-Order approximation"

attr(nlmixr2Est.adgh, "type") <- "Model Based Meta Analysis"
attr(nlmixr2Est.adgh, "description") <- "Aggregate data, Gauss-Hermite quadrature"

attr(nlmixr2Est.adirmc, "type") <- "Model Based Meta Analysis"
attr(nlmixr2Est.adirmc, "description") <- "Aggregate data, Iterative Reweighting Monte Carlo"

attr(nlmixr2Est.admc, "type") <- "Model Based Meta Analysis"
attr(nlmixr2Est.admc, "description") <- "Aggregate data, Monte Carlo"
