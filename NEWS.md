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
