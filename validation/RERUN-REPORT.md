# RERUN REPORT — the mismatched-targets bug in the covariate-quadrature scripts

Date: 2026-08-16. Scope: re-derivation of every conclusion that rested on an arm
scoring **marginal observations against conditional predictions**.

Nothing in `NUMBERS-REGISTRY.md` was edited. This file states, for each affected
registry label, the OLD value, the NEW value, and whether the conclusion is
**unchanged**, **changed**, or **withdrawn**.

---

## 0. VERDICT ON THE CRITICAL QUESTION — the headline result STANDS

**The ~68% attenuation headline is NOT affected by this bug.**

`aggregate-marginal.R` and `overnight-simulation.R` both use **construction 2**:
observed and predicted blocks are conditional on the **same** covariate node.

Code evidence — `aggregate-marginal.R:106` builds the observed side as
`mall(r$thA, GA, cl$f)` and `fitb` (lines 88–91) builds the predicted side as
`mall(p, GL[[s]], FM[s])`, then scores

```r
nl1(BL[[s]][[k]], pr[[k]], GL[[s]]$wk[k] * GL[[s]]$n)
```

`BL[[s]][[k]]` and `pr[[k]]` share **both** the block index `k` and the grid
object `GL[[s]]`. `overnight-simulation.R:152` is the identical call with
`mkgrid(DA)`. There is no pooled `(E, V)` anywhere in either objective.

Numeric confirmation — `headline-construction-check.R` (new, in `validation/`),
which re-implements `mkgrid`/`mall`/`nl1` verbatim:

| test | result |
|---|---|
| obs block *k* vs pred block *k*, same parameters, both grids | `max abs diff` in E and V = **0.000e+00** (K = 16 and K = 4) |
| the same predictions against a POOLED `(E, V)` from the same grid | `max abs diff` in E = **1.18e+01** / **6.83e+00** — plainly not zero |
| objective, matched vs mismatched (BOTH stratified) | **803.41** vs 1083.57 — numerically distinguishable |
| profile `bCL2`, everything else at truth, **matched** | **0.4500** (truth 0.450, **+0.0%**) |
| profile `bCL2`, same grid, **mismatched** (hypothetical) | 0.3515 (**−21.9%**) |

The matched construction recovers `bCL2` **exactly**; a hypothetical mismatched
version of the same code would give −21.9%, which is neither the reported −68%
nor anything the scripts print. The headline attenuation therefore comes from
what the registry already documents in §0.1/§4.3 — **stratifying a source on a
covariate its own model never conditioned on, so every node returns the same
flat answer** — and not from the target mismatch fixed here. Those are two
different defects; only the first is the headline's subject.

**`AGGMARG.*`, `OVN.*` and the §7 conflict-3 resolution stand as written.**

---

## 1. What the bug was

Let `f(a, η)` be the individual model, `g(a) = E_η[f|a]`, `Vc(a) = Cov_η(f|a)`.

| # | construction | definition |
|---|---|---|
| **1** | **marginal / correct GH** | `E = E_a[g(a)]`, `V = E_a[Vc(a)] + Cov_a(g(a))`, then **one** `nll2` |
| **2** | **stratified / matched** | `Σ_k w_k · n · nll2(obs_k, pred_k)`, both conditional on node *k* |
| **3** | **Taylor (corrected)** | 2nd-order expansion of the **moments**, then **one** `nll2` |
| **4** | **mismatched targets** | `Σ_k w_k · nll2(obs_MARGINAL, pred_COND(a_k), n)` — a category error |

The arms named `gh` and `taylor` in `covariate-threeway.R`,
`covariate-matched-conditional.R`, `covariate-node-retest.R` (rows A/B/C) and
target **C** in `covariate-integration-comparison.R` are all construction 4.
The arm named `marginal` in those scripts was **already** the correct
Gauss-Hermite implementation — it integrates the moments over 21 GH nodes and
evaluates the likelihood once. It is renamed **`gh_marginal`** throughout the
`-v2` scripts so nothing reads it as a non-quadrature alternative.

Construction 4 is not a quadrature method with a truncation error. Its
displacement is a difference of *targets*: the observed `(E, V)` already has the
covariate integrated out, while each predicted block pretends the whole study
sat at one covariate value. Averaging −2LLs over nodes adds a term
`tr(V⁻¹ Cov_a(μ))` that carries no data and is minimised by driving the
covariate coefficient to zero — the decomposition reproduced in
`covariate-node-retest-v2.R`.

---

## 2. Corrected scripts

Originals are untouched. Each `-v2` file is headed by a comment naming the wrong
arm and why.

| new file | replaces | console output |
|---|---|---|
| `covariate-threeway-v2.R` | `covariate-threeway.R` | `out-threeway-v2.txt` |
| `covariate-matched-conditional-v2.R` | `covariate-matched-conditional.R` | `out-matched-v2.txt` |
| `covariate-node-retest-v2.R` | `covariate-node-retest.R` | `out-node-v2.txt` |
| `covariate-integration-comparison-v2.R` | `covariate-integration-comparison.R` | `out-integ-v2.txt` |
| `headline-construction-check.R` | (new) verifies §0 | run inline, see §0 |

All exit 0. Data are exact/noiseless throughout, so every displaced argmin is a
property of the objective.

---

## 3. REGISTRY CORRECTIONS, label by label

### 3.1 §2.4 `NODE.argmin` — `covariate-node-retest.R`

**The `−59.4%` is construction 4. It is an implementation artefact of a category
error, not a property of the stratified construction.**

Row-by-row classification and the re-run (`covariate-node-retest-v2.R`):

| row | construction | OLD argmin | NEW argmin | status |
|---|---|---|---|---|
| **M** marginal (shipped) | **1** correct | 0.7500 (+0.0%) | **0.7500 (+0.0%)** | **unchanged** |
| **D** gh, NODE data + `V_cond` | **2** stratified | 0.7500 (+0.0%) | **0.7500 (+0.0%)** | **unchanged** |
| **T** taylor on the MOMENTS | **3** approx | — (did not exist) | **0.7354 (−1.9%)** | **new** |
| **A** gh, pooled data, `V_cond` | **4 MISMATCHED** | 0.3045 (**−59.4%**) | 0.3045 (−59.4%) | **value reproduces; label must change** |
| **B** gh, pooled data, `V_marginal` | **4 MISMATCHED** | 0.4567 (−39.1%) | 0.4567 (−39.1%) | same |
| **C** gh, cond obs V + `V_cond` | **4 MISMATCHED** | 0.3248 (−56.7%) | 0.3248 (−56.7%) | same |

`NODE.decomp` (0.852572 / 0.000000 / 0.852572, identity matches) **reproduces
bit-for-bit and is unchanged**, but it should be labelled as what it is: the
decomposition of **construction 4's** mean term, i.e. the mechanism by which the
category error drives `b → 0`. It is not a property of the shipped objective.

**Required registry edits.** Relabel rows A/B/C as *mismatched targets*, not as
"the node route". Add row T. Rewrite the "Reading" paragraph: the correct
statement is that rows A–C are not a method at all, and that both the marginal
(row M, construction 1) and the stratified (row D, construction 2) constructions
are exact. **Delete the framing "is the node objective biased, or merely a
different objective?"** — it presupposes A is a candidate objective.

### 3.2 §2.5 `MATCHED.onepop` / `MATCHED.threepop` — `covariate-matched-conditional.R`

**`MATCHED.threepop`** — truth 0.7500:

| arm | construction | OLD | NEW | status |
|---|---|---|---|---|
| marginal → **`gh_marginal`** | 1 | 0.7500 (+0.0000) | **0.7500 (+0.0000)** | unchanged (rename only) |
| — → **`taylor_moments`** | 3 | — | **0.7353 (−0.0147)** | new |
| `gh_matched` → **`stratified`** | 2 | 0.7501 (+0.0001) | **0.7501 (+0.0001)** | unchanged (rename only) |
| `gh` | **4** | 0.3045 (−0.4455) | 0.3045 (−0.4455) | **relabel `mismatched`; withdraw as a "gh" result** |
| `taylor` | **4** | 0.3128 (−0.4372) | 0.3128 (−0.4372) | **relabel; withdraw as a "taylor" result** |

**`MATCHED.onepop`** — truth 0.7500:

| arm | construction | OLD | NEW | status |
|---|---|---|---|---|
| marginal → `gh_marginal` | 1 | 0.8956 (+0.1456, +19.4%) | **0.8956 (+0.1456)** | **value unchanged; interpretation must be sharpened — see below** |
| `taylor_moments` | 3 | — | **0.0637 (boundary)** | new; boundary, do not quote |
| `gh_matched` → `stratified` | 2 | 0.7500 (−0.0000) | **0.7500 (−0.0000)** | unchanged |
| `gh` | **4** | 0.0502 (boundary) | 0.0502 (boundary) | relabel; already uncitable |
| `taylor` | **4** | 0.0502 (boundary) | 0.0502 (boundary) | relabel; already uncitable |

> **Sharpening required on the "+19.4%".** The re-run adds the profile spread
> over `tcov ∈ 0.45…1.05`. For one population, `gh_marginal` gives
> **`0.000 0.000 0.000 0.000 5.912`** — the objective is **exactly flat** from
> 0.45 to at least 0.90. So `0.8956` is not a bias of magnitude 19.4%; it is
> where an optimiser lands inside a flat region. The registry's framing ("an
> identification statement about the design, not a defect of the method") is
> **correct and should be kept**, but the figure **+19.4% must not be quoted as
> a bias magnitude**. Quote the flat profile instead. Construction 3 landing at
> the boundary in the same cell is the same flatness, not a separate defect.
>
> For three populations the profiles are strongly curved (`gh_marginal`
> `43.903 12.416 0.000 14.797 60.511`) and all of constructions 1, 2, 3 recover
> the truth.

### 3.3 §2.6 `THREEWAY.bias` / `THREEWAY.perN` — `covariate-threeway.R`

**`THREEWAY.bias`** — invariant to n (10/100/1000), confirmed again:

| design | `gh_marginal` [1] | `taylor_moments` [3] | `mismatched` [4] | old `taylor` [4t] |
|---|---|---|---|---|
| rich | **0.7500** (+0.0000) | **0.7353** (−0.0147) | 0.3045 (−0.4455) | 0.3128 (−0.4372) |
| sparse | **0.7501** (+0.0001) | **0.7401** (−0.0099) | 0.3477 (−0.4023) | 0.3624 (−0.3876) |

- the `marginal` column is **unchanged**, rename to `gh_marginal`;
- the `gh` column is **construction 4** — keep the number only if it is
  explicitly labelled *mismatched targets*; **withdraw it as "gh"**;
- the `taylor` column is **withdrawn and replaced**: the old −0.4372 / −0.3876
  measured a 2nd-order expansion of the category error. The corrected Taylor
  (construction 3) is **−0.0147 / −0.0099**, i.e. **a 30× smaller error**, and
  it is a usable approximation.

**`THREEWAY.perN`** — per-subject objective at the true parameters:

| design | OLD label | OLD value | NEW label | NEW value |
|---|---|---|---|---|
| rich | marginal/n | −22.86914 | `gh_marginal`/n | **−22.86914** (unchanged) |
| rich | taylor/n | −15.55568 | `taylor_moments`/n | **−22.76895** |
| rich | gh/n | −14.39631 | `mismatched`/n | −14.39631 (relabel) |
| rich | "taylor − gh = truncation" | **−1.15938** | **`[3]−[1]` truncation** | **+0.10019** |
| rich | "gh − marginal = target" | +8.47283 | `[4]−[1]` target mismatch | +8.47283 (relabel) |
| sparse | marginal / taylor / gh | 2.86045 / 7.36017 / 8.07093 | `[1]` / `[3]` / `[4]` | **2.86045 / 2.86445** / 8.07093 |
| sparse | truncation / target gaps | −0.71076 / +5.21048 | `[3]−[1]` / `[4]−[1]` | **+0.00400** / +5.21048 |

> **The "truncation error" figures −1.15938 and −0.71076 are WITHDRAWN.** They
> were the truncation error of an expansion of the wrong functional. The
> corrected truncation errors are **+0.10019 (rich)** and **+0.00400 (sparse)**,
> i.e. an order of magnitude smaller and one to two orders smaller respectively.
> Any prose saying "Taylor's truncation error is large" must be reversed.

New, worth adding — moment accuracy of construction 3 against construction 1 at
the true parameters (rich design):

| population | `mu_a` | `sd_a` | rel err E | rel err V |
|---|---|---|---|---|
| pop1 | −0.45 | 0.30 | 1.853e−04 | 2.167e−02 |
| pop2 | 0.00 | 0.55 | 2.252e−03 | 1.096e−01 |
| pop3 | 0.50 | 0.35 | 4.057e−04 | 2.737e−02 |

The expansion's weak point is `V`, and it degrades with `sd_a` — which is the
regime statement §3.4 makes quantitative.

### 3.4 §2.7 `INTEG.*` — `covariate-integration-comparison.R`

**`INTEG.targets`** — the "three targets" table is restructured, because (B) and
(C) are not targets a published aggregate dataset can have. Truth 0.7500,
omega truth 0.3000:

| row | construction | OLD minimiser | NEW minimiser | omega OLD | omega NEW | status |
|---|---|---|---|---|---|---|
| (A) marginal → **`gh_marginal`** | **1** | 0.7500 (+0.0000) | **0.7500 (+0.0000)** | 0.3000 | **0.3000** | **unchanged** |
| **`stratified`** | **2** | — | **0.7500 (+0.0000)** | — | **0.3000** | **new** |
| **`taylor_moments`** | **3** | — | **1.0864 (+44.9%)** | — | **0.5970 (+0.2970)** | **new — see warning** |
| (C) `E[NLL]` → **`mismatched`** | **4** | 0.2000 (−73.3%) | 0.2000 (−73.3%) | 1.2000 | 1.2000 | **value reproduces; relabel as the category error** |
| (B) mixture | 4-family | 0.8390 (+11.9%) | **0.8390 (+11.9%)** | — | — | unchanged; it is a valid likelihood for a *different* data model (all subjects share one unknown `a`) — keep, but do not present as an option for summary data |

> **New finding — construction 3 fails in this scenario, and that is real.** This
> script uses `sd_a = 1` with `omega = 0.30`, i.e. `tcov·sd_a/omega = 2.5` at
> truth: the covariate spread is 2.5× the IIV. The moment expansion is a
> second-order approximation in `sd_a` and is simply out of its regime there
> (+44.9% on `tcov`, +99% on `omega`). Contrast `THREEWAY.bias`, where
> `sd_a ≤ 0.55` and construction 3 errs by −1.9%. **The corrected Taylor is not
> unconditionally safe** — §3.4's regime table below is the boundary.

**`INTEG.taylor.h` — WITHDRAWN.** The old sweep
(`h=0.05 → −18758.5 … h=3.00 → +14249.0`) is the differencing step of
construction 4. It measures the step-size sensitivity of the wrong functional and
supports no claim about Taylor expansion as a technique. The *separate* defect it
also illustrated — that `h` was in **absolute covariate units**, so a default
`h = 2.0` lands arbitrarily for a covariate on any other scale — **remains valid
and should be kept**, but must be cited from the corrected sweep. Replacement, with
`h` as a fraction of `sd_a`, errors against construction 1:

| `hfrac` | rel err E | rel err V | NLL err vs [1] |
|---|---|---|---|
| 0.05 | 1.072e−01 | 3.821e−01 | 2142.33 |
| 0.10 | 1.067e−01 | 3.805e−01 | 2119.89 |
| 0.25 | 1.039e−01 | 3.696e−01 | 1970.86 |
| 0.50 | 9.404e−02 | 3.316e−01 | 1534.79 |
| 0.75 | 7.845e−02 | 2.721e−01 | 1043.16 |
| 1.00 | 5.830e−02 | 1.961e−01 | 650.20 |
| 2.00 | 4.056e−02 | 1.758e−01 | 812.96 |

(These are large because `sd_a = 1` puts the expansion out of regime; the point of
the table is that the error is now *monotone and interpretable in `sd_a` units*
rather than sign-flipping across the sweep as the old one did.)

**`INTEG.ess` — UNCHANGED, and explicitly re-verified.** `ESS/N = exp(−(D/omega)²)`
reproduces exactly: `0.5 → 7.79e−01` (15576 effective), `1.0 → 3.68e−01` (7358),
`2.0 → 1.83e−02` (366), `3.0 → 1.23e−04` (2.5), `4.0 → 1.13e−07` (0.0). The ESS
mechanics do not depend on what is combined afterwards.

**But the IS accuracy figure is replaced.** OLD: "at `tcov = 0.75` the IS absolute
error is **4308.06**" — that was error against construction 4. The corrected arm
reweights the pool to per-node **moments**, combines them into marginal moments,
and evaluates one likelihood:

| shift/omega | `tcov` | [1] exact | IS-moments | abs err | min ESS/N |
|---|---|---|---|---|---|
| 0.33 | 0.05 | 8640.471 | 8650.174 | 9.70 | 5.69e−01 |
| 0.67 | 0.10 | 7746.849 | 7767.370 | 20.52 | 1.13e−01 |
| 1.33 | 0.20 | 5176.798 | 5297.598 | 120.80 | 2.27e−03 |
| 2.33 | 0.35 | 1751.091 | 2397.202 | 646.11 | 1.69e−04 |
| 3.33 | 0.50 | −56.846 | 1045.176 | 1102.02 | 8.13e−05 |
| 5.00 | 0.75 | −694.365 | 285.301 | **979.67** | 5.74e−05 |

**Conclusion unchanged**: IS reweighting is unusable at realistic covariate effect
sizes. Only the number changes, **4308.06 → 979.67**.

**`INTEG.scollapse` — UNCHANGED**, reproduced verbatim: `E 9.99e−16, V 5.28e−14`.
This section never involved construction 4.

**`INTEG.regime` — REPLACED.** OLD: relative errors against target C, with `gl-9`
`≤ 8.24e−04` and `taylor2` reaching `7.52e−01`, "target gap" `2.63e−03 → 2.35e+02`.
NEW — relative errors against construction 1:

| ratio | `tcov` | `gh_mom-9` | `taylor_moments` | IS-moments | `stratified` | mismatch size |
|---|---|---|---|---|---|---|
| 0.10 | 0.030 | 6.50e−16 | 2.72e−10 | 4.30e−08 | 2.56e−03 | 2.63e−03 |
| 0.25 | 0.075 | 3.73e−15 | 3.70e−07 | 1.40e−07 | 1.59e−02 | 1.87e−02 |
| 0.50 | 0.150 | 1.95e−15 | 6.73e−05 | 3.05e−06 | 6.27e−02 | 1.09e−01 |
| 1.00 | 0.300 | 2.81e−15 | 7.78e−03 | 8.99e−04 | 2.46e−01 | 1.12e+00 |
| 2.00 | 0.600 | 1.42e−07 | 5.93e−01 | 2.73e−01 | 1.14e+00 | 2.41e+01 |
| 3.00 | 0.900 | 6.26e−04 | 6.85e+00 | 5.70e+00 | 4.20e+00 | 2.35e+02 |

Readings: **9-node GH over the moments is essentially exact** to ratio ≈ 2
(`1.4e−07`) and still `6.3e−04` at ratio 3 — better than the old `gl-9` column and
now against the right reference. The **corrected Taylor is good to ~1e−2 up to
ratio 1 and breaks beyond ratio 2** — this is the usable statement about
construction 3 and it should replace any blanket claim that Taylor is inadequate.
The last column (formerly "target gap") is the **size of the category error** and
is the number to quote when explaining why construction 4 must not be used:
`2.6e−03 → 2.35e+02` as the covariate effect grows.

Also new, §2 of the corrected script — approximation error against construction 1
at `tcov = 0.75`: `gh_mom-5` **3.74e+00**, `gh_mom-9` **1.25e−02**,
`taylor_moments` 1.54e+03, `mismatched` 5.42e+04.

### 3.5 §4.1 — "0.196 against a truth of 0.75 (−74%), omega inflated ~40%"

The registry already corrected the numbers to `0.2000 / −73.3%` and
`omega 1.2000 / +300%`. **Those corrected numbers reproduce exactly.** What must
change is the *label*: they are **construction 4**, i.e. the size of a category
error, **not** "how biased the aggregate method is". The paragraph currently
recommends `NODE.argmin` row A (−59.4%) as "the better-conditioned measurement of
the same phenomenon" — that is still true (both are construction 4) but the
phenomenon must be renamed. Suggested wording: *"scoring a marginal observation
against conditional predictions displaces the coefficient by −59% to −73% and
inflates omega 4×; it is a mismatch of targets, not a quadrature error."*

### 3.6 §7 conflict 1 — "How biased is the mismatched node route?"

**Conclusion changed in framing, unchanged in arithmetic.** All five quoted
figures (−59.4%, −56.7%, −43.7%, −73.3%, −93.3%) are construction 4 or an
expansion of it. The registry's advice — believe `NODE.argmin` and quote −59.4%
— is still the right choice **provided the quantity is renamed** from "the node
route" to "the mismatched-targets construction". The genuine node/stratified
route is `NODE.argmin` row **D**, which is **exact (0.0%)**.

The one figure that must be *removed* from that list is `THREEWAY.bias` taylor
**−43.7%**: it is now superseded by the corrected Taylor's **−1.9%**, and leaving
it in the list implies Taylor expansion is a 44%-biased technique, which the
re-run refutes.

### 3.7 §7 conflict 2 — "Is the marginal/aggregate objective biased at all?"

**Unchanged in substance.** Three populations: 0.0% under constructions 1, 2 and
−1.9% under 3. One population: flat objective, argmin lands at 0.8956. Keep the
identification framing; drop "+19.4%" as a bias magnitude (see §3.2).

### 3.8 §7 conflict 3 and §2.1/§2.2 — **unaffected**

See §0. `AGGMARG.*` and `OVN.*` are construction 2 with matched blocks.

---

## 4. Other scripts in `validation/` containing a construction-4 arm

All 36 remaining `.R` files in `validation/` were read and classified. **Exactly
one carries a construction-4 arm, and it was not re-run here.**

### 4.1 NOT RE-RUN — contains construction 4

| script | arm | lines | why not re-run |
|---|---|---|---|
| **`page-abstract-taylor.R`** | `obj()`, fed by `nll_node()` | objective **111–123**, node construction **78–90** | Needs `rxode2` + `nlmixr2est` + a full `datagen()` + two BOBYQA fits (minutes to tens of minutes). Out of the pure-base-R scope of this task. |

Both construction-4 tells are present at once. Lines 83–87 take the **marginal**
generated study (`gen[[nm]]`, built at 67–72 *with* `cov_dist`, i.e. the covariate
already integrated out of `E`/`V`), strip `cov_dist`, and pin only the
**prediction** to a covariate value:

```r
s <- gen[[nm]]; s$cov_dist <- NULL
s$cov <- list(WT = p$mu + off * h)
```

Lines 116–120 then finite-difference the **NLL scalars**:

```r
f0 <- nll_node(p, paste0(nm, "_0")); fp <- ...("_p"); fm <- ...("_m")
tot <- tot + f0 + 0.5 * sd^2 * (fp - 2*f0 + fm) / h^2
```

That is `Σ_k w_k · nll2(obs_marginal, pred_cond(a_k))` expanded to second order —
the same object as the old `covariate-threeway.R` `taylor` arm, now applied to the
PAGE vancomycin case and to a real `admixr2` `.adghNLL`.

**Its script header (lines 10–17) says this is deliberate**: it is pricing
`.adm_combine_nll(method = "taylor")` — what the *package* does today — against the
marginal fit at lines 144–149. So the script is not itself confused; it is a
measurement of a package code path. **But the finding it produces must now be
stated as "the package's taylor combine implements construction 4", not as "the
Taylor approximation is inaccurate".** The corrected Taylor (construction 3) errs by
**−1.9%**, not −44%, in the regime where it applies (§3.3/§3.4).

**Registry exposure: none.** §6 lists all `page-abstract-*.R` as "NOT RUN TODAY, no
`.rds`", so no registry label depends on this script's output. It should be re-run
with a construction-3 arm before any PAGE Taylor claim is made, and it belongs on the
§8 re-run list.

### 4.2 Covariate-node work, but NOT construction 4 — no action

Constructions 1 (marginal moments) and/or 2 (matched per-node blocks), verified by
the obs side being a per-node **list** indexed with the same `k` and the same grid
object as the prediction:

`aggregate-marginal-taylor.R` (1+2+**3** — its `mallT`, lines 123–134, is a correct
Taylor on the *moments*, matching `taylor-corrected.R`), `mc-averaging-harness.R`,
`mc-averaging-study.R`, `mc-averaging-nocontrast.R`, `mc-averaging-validate.R`,
`three-way-replicates.R`, `gma-vs-ipd-stage2.R`, `gma-minimal.R`,
`model-synthesis-omitted-covariate.R`, `model-synthesis-two-likelihoods.R`,
`model-synthesis-standalone.R`, `theory-tests.R`, `covariate-collapse-endtoend.R`,
`covariate-uq-endtoend.R`, `covariate-varonly-robustness.R`,
`framework-simulation-estimation.R`, `covariate-u-distribution.R`,
`covariate-dependent-copula.R`, `page-abstract-mbma.R`, `page-abstract-debug.R`,
`page-abstract-varonly.R`, `page-framework-test.R`, `page-weights.R`,
`page-param-likelihood.R`, `param-likelihood-1cmt-2cmt.R`,
`param-likelihood-peak-trough.R`, `param-likelihood-two-publications.R`.

Two of these deserve a note because they *look* like construction 4 and are not:

- **`theory-tests.R` T2 (lines 95–101)** computes the mismatched quadratic form
  `Σ_k w_k r_k' V⁻¹ r_k` with a marginal `obs` and conditional `MU` — **on
  purpose**, as an algebraic identity, with no fit and no optimiser. It is the file's
  proof that construction 4 decomposes into `rbar'V⁻¹rbar + tr(V⁻¹Cov_a(μ))` and
  therefore attenuates `b`. `THEORY.T2` is **unchanged and correctly labelled**;
  it is a *diagnosis* of the bug, not an instance of it.
- **`gma-minimal.R`'s `gb_pool` (lines 81–85)** is pooled-marginal on **both** sides
  in a single `nll1` call with no node loop — a construction-1 control deliberately
  printed as "the earlier error". Not construction 4.

### 4.3 No covariate-node integration at all

`model-synthesis-formal-likelihood.R`, `covariate-rvine-sampling-check.R`,
`gma-vs-ipd.R`, `regression-check.R`, `overnight-summary.R`,
`page-abstract-panelA.R`, `page-abstract-weightplot.R`, `page-cross-dosing.R`.

---

## 5. Summary of what to correct or delete in `NUMBERS-REGISTRY.md`

**Delete / withdraw**

- `INTEG.taylor.h` (the whole six-value sweep) — step size of the wrong functional.
- `THREEWAY.perN` "taylor truncation error" **−1.15938 / −0.71076**.
- `THREEWAY.bias` `taylor` column **−0.4372 / −0.3876** as a Taylor result.
- `MATCHED.*` `gh` and `taylor` rows as "gh"/"taylor" results.
- "+19.4%" as a *bias magnitude* in `MATCHED.onepop` (keep the identification point).
- §7 conflict 1's inclusion of **−43.7%** in the list of mismatch measurements.

**Relabel (numbers reproduce, names are wrong)**

- `NODE.argmin` rows A/B/C → *mismatched targets*, not "the node route".
- `NODE.decomp` → the decomposition of construction 4's mean term.
- `THREEWAY.bias` `gh` column, `MATCHED.*` `gh` rows, `INTEG.targets` target C →
  all *mismatched targets*.
- `marginal` → `gh_marginal` everywhere (it **is** the GH implementation).
- `gh_matched` → `stratified` (construction 2).
- §4.1's `−73.3%` / `+300%` → size of a category error.

**Replace with a new value**

| label | old | new |
|---|---|---|
| `THREEWAY.bias` taylor, rich / sparse | −0.4372 / −0.3876 | **−0.0147 / −0.0099** |
| `THREEWAY.perN` taylor/n, rich / sparse | −15.55568 / 7.36017 | **−22.76895 / 2.86445** |
| `THREEWAY.perN` truncation, rich / sparse | −1.15938 / −0.71076 | **+0.10019 / +0.00400** |
| `MATCHED.threepop` taylor | 0.3128 (−0.4372) | **0.7353 (−0.0147)** |
| `INTEG` IS absolute error at `tcov = 0.75` | 4308.06 | **979.67** |
| `INTEG.regime` method columns | vs target C | **vs construction 1** (table in §3.4) |

**Add**

- `NODE.argmin` row **T** (construction 3): **0.7354 (−1.9%)**.
- `INTEG.targets` rows `stratified` (**0.7500**, omega **0.3000**) and
  `taylor_moments` (**1.0864**, omega **0.5970**).
- A regime boundary for construction 3: exact to ~1e−2 up to
  `tcov·sd_a/omega ≈ 1`, unusable beyond ≈ 2.
- `gh_mom-9` accuracy against construction 1: **1.42e−07** at ratio 2,
  **6.26e−04** at ratio 3.

**Unchanged and re-verified**

`NODE.argmin` rows M and D; `MATCHED.threepop` marginal and `gh_matched`;
`THREEWAY.bias` marginal; `INTEG.targets` target A; `INTEG.scollapse`;
`INTEG.ess` decay law; `AGGMARG.*`; `OVN.*`; `THEORY.T2`/`T3`.

**Add to §8 (what to re-run before the paper is finalised)**

- `page-abstract-taylor.R` with a construction-3 arm beside its current
  construction-4 one. It is the only remaining script with the bug, and it is the
  only place the bug is measured **inside a real `admixr2` code path**
  (`.adm_combine_nll(method = "taylor")`), so the corrected version is the one that
  decides whether that package option should exist in its present form.
