# fdid 1.0.2

* Added `fdid_contrast()` for continuous-G level comparisons and
  `fdid_derivative()` for fixed-G derivative extraction from `kernel` and
  `dml_flex` curve objects.
* Stored kernel bootstrap curve replicates and `dml_flex` multiplier curve
  replicates for reporting-level contrast and derivative inference when those
  resampling routines run.
* Clarified continuous-G curve-grid semantics in package markdown docs:
  `eval_g` is the estimator grid for `kernel` and `dml_flex`, while
  `plot(type = "curve")` displays stored curve estimates.
* Aligned continuous-curve support overlays with the displayed curve window by
  using observed `G` values restricted to the plotted range.
* Improved continuous-G curve plotting by drawing `curve = "both"` level and
  derivative panels side-by-side with embedded support displays.
* Updated the rendered tutorial so DML examples use `learner = "grf"` and
  kernel/DML estimator chapters use a consistent histogram support display.
* Refined continuous-G tutorial figures: kernel bootstrap and contrast examples
  now use the mortality application throughout, DML mortality derivative and
  incremental displays emphasize point estimates, and a synthetic DML example
  illustrates curve uncertainty.
* Polished continuous-G plotting defaults with compact curve legend labels,
  automatic legend placement, suppressible legends via `legend.position`, and
  cleaner contrast-plot axis labels.
* Updated the continuous-G tutorial to clarify that `eval_g = NULL` uses
  `trim = 0.05` by default; the kernel chapter now deliberately sets
  `trim = 0` to estimate the full observed mortality `G` support, while DML
  curve examples keep the default 50-point grid.
* Consolidated the synthetic DML curve workflow so one example carries the
  data-generating process, flexible-curve fit, curve inference display, and
  reporting workflow; restored learner choice and high-dimensional expansions
  to top-level tutorial sections.
* Updated the kernel tutorial bootstrap example to render with `boot = 200`
  and revised the DML continuous-G chapter to use one synthetic DGP throughout
  the scalar, curve, inference, reporting, incremental, learner, and
  high-dimensional examples.
* Clarified why kernel bootstrap uniform-band effective draw counts can be
  smaller than `boot`: finite-grid uniform envelopes use only bootstrap curves
  that are finite at every stored `eval_g` point.
* Updated the kernel tutorial reporting examples so
  `fdid_contrast()`/`fdid_derivative()` use the analytical robust kernel fit
  with pointwise inference, while the bootstrap section remains the source for
  bootstrap pointwise and uniform-band displays.
* Hid manual table-construction code in the Chapter 4 and Chapter 5 rendered
  tutorial while keeping substantive `fdid` calls visible.
* Moved automatic plot legends into a boxed area below the plot so continuous-G
  curve legends no longer cover fitted curves or interval layers.
* Replaced `options(warn = -1)` with `suppressWarnings()` in internal helper
  `silent_ebalance()` for CRAN compliance (avoids global state modification).

# fdid 0.2.0

* Initial public release.
* Implements six estimation methods: `ols1`, `ols2`, `did`, `ebal`, `ipw`, `aipw`.
* Supports robust, bootstrap, and jackknife variance estimation.
* Parallel computation via `foreach`/`doParallel` for bootstrap and jackknife.
* `fdid_prepare()` for reshaping long-format panel data.
* `fdid_list()` for comparing multiple estimates.
* S3 methods: `print`, `summary`, `plot` for `fdid` and `fdid_list` objects.
* Built-in `mortality` dataset (Cao et al. 2022, China Great Famine).
* Reference: Xu, Zhao & Ding (2026) <doi:10.1080/01621459.2026.2628343>.
