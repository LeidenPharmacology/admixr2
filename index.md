<div class="adm-hero">
<div class="adm-hero-text">
<p class="adm-eyebrow">Aggregate data modelling · nlmixr2 / rxode2</p>
<h1 class="adm-h1">A meta-analysis framework for population PK/PD</h1>
<p class="adm-tag">Fit one unified population model to what the literature reports — aggregate summary data or previously published PK/PD models — without individual patient data.</p>
<p class="adm-badges"><a href="https://github.com/LeidenPharmacology/admixr2/actions/workflows/R-CMD-check.yaml"><img src="https://img.shields.io/github/actions/workflow/status/LeidenPharmacology/admixr2/R-CMD-check.yaml?style=flat-square&amp;label=R-CMD-check" alt="R-CMD-check"></a> <a href="https://cran.r-project.org/package=admixr2"><img src="https://img.shields.io/cran/v/admixr2?style=flat-square" alt="CRAN"></a> <a href="https://app.codecov.io/gh/LeidenPharmacology/admixr2"><img src="https://img.shields.io/codecov/c/github/LeidenPharmacology/admixr2?style=flat-square&amp;label=coverage" alt="Coverage"></a> <a href="https://lifecycle.r-lib.org/articles/stages.html#stable"><img src="https://img.shields.io/badge/lifecycle-stable-4A6B5D?style=flat-square" alt="Lifecycle: stable"></a> <a href="https://doi.org/10.1007/s10928-025-10011-w"><img src="https://img.shields.io/badge/DOI-10.1007%2Fs10928--025--10011--w-1F4E79?style=flat-square" alt="DOI"></a></p>
<p class="adm-cta"><a class="btn btn-primary btn-lg" href="articles/admixr2.html">Get started</a> <a class="btn btn-outline-secondary btn-lg" href="https://github.com/LeidenPharmacology/admixr2">GitHub</a></p>
<p class="adm-install"><span class="adm-prompt">&gt;</span> <code>install.packages(&quot;admixr2&quot;)</code></p>
</div>
<div class="adm-card">
<div class="adm-card-top"><i></i><i></i><i></i><span>meta-analysis.R</span></div>
<pre class="adm-code"><code># Fit a model to published means + covariances
fit &lt;- nlmixr2(model, admData(), est = "admc",
  control = admControl(studies = list(
    trial_A = list(E = Ea, V = Va, n = 120L,
                   times = c(1, 2, 4, 8, 24),
                   ev = et(amt = 100)))))

plot(fit)   # observed mean vs prediction, ±1 SD</code></pre>
<div class="adm-plot">
<svg viewBox="0 0 460 150" width="100%" style="display:block"><line x1="42" y1="12" x2="42" y2="122" stroke="var(--bs-border-color)" stroke-width="1"/><line x1="42" y1="122" x2="440" y2="122" stroke="var(--bs-border-color)" stroke-width="1"/><path d="M42,30 C120,44 210,74 300,96 C350,108 400,114 440,117 L440,127 C400,124 350,118 300,106 C210,84 120,54 42,40 Z" fill="var(--adm-navy)" opacity="0.12"/><path d="M42,35 C120,49 210,79 300,101 C360,114 410,118 440,122" fill="none" stroke="var(--adm-navy)" stroke-width="2"/><g stroke="var(--bs-body-color)" stroke-width="1.3"><line x1="90" y1="40" x2="90" y2="58"/><line x1="175" y1="62" x2="175" y2="82"/><line x1="260" y1="82" x2="260" y2="100"/><line x1="345" y1="100" x2="345" y2="116"/></g><g fill="var(--adm-orange)"><circle cx="90" cy="49" r="3.6"/><circle cx="175" cy="72" r="3.6"/><circle cx="260" cy="91" r="3.6"/><circle cx="345" cy="108" r="3.6"/></g></svg>
<div class="cap">observed study mean ± SD &nbsp;vs&nbsp; model prediction</div>
</div>
</div>
</div>

<div class="adm-steps">
<div class="adm-step">
<div class="adm-inpair">
<div class="adm-inp"><div class="ico-sm"><svg width="19" height="19" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><rect x="4" y="3" width="16" height="18" rx="2"/><path d="M8 8h8M8 12h8M8 16h5" stroke-linecap="round"/></svg></div><div><b>Aggregate data</b><span>published means, error bars &amp; covariances</span></div></div>
<div class="adm-inp"><div class="ico-sm"><svg width="19" height="19" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><rect x="3" y="4" width="18" height="16" rx="2"/><path d="M6 15c3-6 8-6 12 0" stroke-linecap="round"/><circle cx="9" cy="10.5" r="1.2" fill="currentColor" stroke="none"/></svg></div><div><b>PK/PD models</b><span>previously published population models</span></div></div>
</div>
</div>
<div class="adm-arrow"><svg width="28" height="28" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.6"><path d="M4 12h15M14 6l6 6-6 6" stroke-linecap="round" stroke-linejoin="round"/></svg></div>
<div class="adm-step">
<div class="ico"><svg width="26" height="26" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"><path d="M4 4v16h16"/><path d="M7.5 16.5c2.5-6 7-9 11.5-11"/><circle cx="10" cy="13.2" r="1.35" fill="currentColor" stroke="none"/><circle cx="14.5" cy="9.6" r="1.35" fill="currentColor" stroke="none"/></svg></div>
<h4>Aggregate data modelling</h4>
<p>One likelihood over the summary-level information — a meta-analysis across studies.</p>
</div>
<div class="adm-arrow"><svg width="28" height="28" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.6"><path d="M4 12h15M14 6l6 6-6 6" stroke-linecap="round" stroke-linejoin="round"/></svg></div>
<div class="adm-step">
<div class="ico"><svg width="26" height="26" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round"><path d="M3 5c6.5 0 6 7 11.5 7"/><path d="M3 12h11.5"/><path d="M3 19c6.5 0 6-7 11.5-7"/><circle cx="18.6" cy="12" r="2.4"/></svg></div>
<h4>Unified population model</h4>
<p>One fit with interpretable fixed, random and covariate effects.</p>
</div>
</div>

<div class="adm-sec">
<div class="adm-kick">Four estimators, one interface</div>
<h2>Choose the method that fits your problem</h2>
<p class="adm-lead">Every backend plugs into <code>nlmixr2(..., est=)</code> with analytical gradients. Compare them side by side in the <a href="articles/estimator-comparison.html">estimator vignette</a>.</p>
</div>
<div class="adm-cards">
<div class="adm-est"><span class="tag">adfo</span><h3>First-order</h3><p>Analytical first-order linearisation — deterministic, no Monte-Carlo draws. Typically the fastest backend.</p></div>
<div class="adm-est"><span class="tag">admc</span><h3>Monte Carlo</h3><p>Simulation of the aggregate likelihood, with analytical common-random-number gradients. Accuracy improves with the number of draws.</p></div>
<div class="adm-est"><span class="tag">adgh</span><h3>Gauss–Hermite</h3><p>Deterministic quadrature over the random effects — noise-free, with an exact gradient. Most efficient when the random-effect dimension is small.</p></div>
<div class="adm-est"><span class="tag">adirmc</span><h3>Iterative RMC</h3><p>Iterative Reweighting Monte Carlo — reweights simulated draws across iterations, with a kappa correction for non-linear models.</p></div>
</div>

<div class="adm-sec">
<div class="adm-kick">Documentation</div>
<h2>Guides, grouped by what you're doing</h2>
<p class="adm-lead">Nine worked vignettes, from your first fit to publication-grade diagnostics.</p>
</div>
<div class="adm-vcards">
<div class="adm-vcard"><h3>Getting started</h3><ul><li><a href="articles/admixr2.html">Introduction to admixr2</a></li></ul></div>
<div class="adm-vcard"><h3>Inputs — data &amp; models</h3><ul><li><a href="articles/aggregate-data.html">From a published figure to E, V and n</a></li><li><a href="articles/datagen.html">Simulating data &amp; using published models</a></li></ul></div>
<div class="adm-vcard"><h3>Meta-analysis across studies</h3><ul><li><a href="articles/multiple-studies.html">Multiple studies</a></li><li><a href="articles/multi-compartment.html">Several outputs (plasma + brain)</a></li><li><a href="articles/pkpd.html">PD and PK/PD data</a></li></ul></div>
<div class="adm-vcard"><h3>Estimation &amp; diagnostics</h3><ul><li><a href="articles/estimator-comparison.html">Comparing the estimators</a></li><li><a href="articles/advanced.html">Advanced usage</a></li><li><a href="articles/diagnostic-plots.html">Diagnostic plots</a></li></ul></div>
</div>

<div class="adm-close">
<div><h5>How to cite</h5><p class="cite">van de Beek H., Välitalo P.A.J., Zwep L.B., van Hasselt J.G.C. (2025). <em>admixr2: Aggregate Data Modelling.</em> Journal of Pharmacokinetics and Pharmacodynamics. <a href="https://doi.org/10.1007/s10928-025-10011-w">doi:10.1007/s10928-025-10011-w</a></p><p class="leiden">Developed at Leiden University</p></div>
<div><h5>Package</h5><a href="articles/admixr2.html">Get started</a><a href="reference/index.html">Reference</a><a href="news/index.html">Changelog</a></div>
<div><h5>Ecosystem</h5><a href="https://nlmixr2.org/">nlmixr2</a><a href="https://cran.r-project.org/package=rxode2">rxode2</a><a href="https://github.com/LeidenPharmacology/admixr2">GitHub</a></div>
</div>
