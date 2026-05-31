# Consolidated TODO

Timestamp: 2026-05-19, America/Chicago

This is the single active TODO document for the continuous-G extension of the
`fdid` package. Historical markdowns were archived under:

```text
historical_mds_2026-05-19/
```

The other two active project markdowns are:

- `CONTINUOUS_G_PROJECT_LOG.md`: detailed log, validation results, and caveats.
- `CONTINUOUS_G_ESTIMATOR_REFERENCE.md`: current estimator and inference
  reference.
- `HANDOFF.md`: compact restart entry point that points back to these active
  docs; it is excluded from package builds through `.Rbuildignore`.

Package-standard files `NEWS.md` and `LICENSE.md` remain in place.

## Tutorial Cleanup Status As Of 2026-05-31

The P0 tutorial reorganization has been completed in the public tutorial
source. The canonical source is now `tutorial/`, with rendered GitHub Pages
output in `docs/`.

Final chapter files:

```text
01-start.Rmd
02-main.Rmd
03-continuous-overview.Rmd
04-kernel-continuous-g.Rmd
05-dml-continuous-g.Rmd
06-visualization.Rmd
07-reporting-inference.Rmd
08-sensitivity.Rmd
```

The obsolete source copies `tutorial-full-continuous-g/` and
`tutorial-next-2026-05-30/` were removed after migration. Remaining tutorial
follow-ups are future plotting-layer modernization documentation if that API is
implemented, and any future theory-note sync if the inference semantics change.
The tutorial cleanup validation is recorded in
`TUTORIAL_CLEANUP_PLAN_2026-05-31.md` and the StatsClaw run `audit.md`.

## Current Status

The package now supports the original binary-G methods and the new
continuous-G/DML methods through the existing `fdid(..., method = ...)` API:

```text
ols1
ols2
did
ebal
ipw
aipw
kernel
dml_binary
dml_plr
dml_flex
dml_incremental
```

Completed continuous-G work includes:

- continuous-G kernel estimator;
- DML binary, PLR, flexible curve, and scalar incremental estimators;
- learner options including linear, lasso, ridge, elasticnet, ranger, grf, gam,
  nnet, and xgboost;
- high-dimensional covariate expansion through `fdid_prepare()` and
  `fdid_expand_covariates()`;
- conditional-density options `residual_kde`, `location_scale`, and
  `local_kde`;
- signal-to-G mapping options `local_poly`, `kernel`, `spline`,
  `blp_spline`, and `gam`;
- cross-fitting with user-controlled `K` folds and repeated-split median
  aggregation through `S`;
- scalar DML robust SEs with cluster support when `cluster_label` is supplied;
- scalar DML score-multiplier intervals through `dml_inference`,
  `dml_boot`, and `dml_multiplier`;
- pointwise curve intervals, practical multiplier bands, and BLP-spline
  covariance/Gaussian max-t bands for `dml_flex`;
- kernel pointwise intervals, bootstrap pointwise intervals, level and
  derivative quantile-envelope bands, and bootstrap curve covariance matrices;
- DML diagnostics for signals, density, derivative nuisance components,
  correction terms, low-density shares, and cluster counts;
- updated plotting with `type = "curve"`, explicit interval controls,
  `curve = "level"`, `"derivative"`, or `"both"`, support histogram/density/rug
  displays, evaluation-grid ticks, and optional marked `G` values;
- fixed-reference continuous-G contrast coefficient plots through
  `plot(type = "contrast")` for kernel and `dml_flex`, including a default
  algebraic-zero reference row, character `ref.g` rules, and quantile-based
  target selection;
- review-driven package-impact cleanup: DML Gaussian max-t no longer uses a
  local variable named `T`, `HANDOFF.md` is build-ignored, and the latest
  tutorial copy avoids stale local-install and dependency-check examples;
- first-class reporting helpers `fdid_contrast()` and `fdid_derivative()` for
  level comparisons and fixed-`g` derivative values;
- optional storage of kernel bootstrap curve replicates, kernel bootstrap curve
  covariance matrices, `dml_flex` multiplier curve replicates, and
  `dml_flex` BLP-spline covariance matrices for reporting inference;
- validation scripts for oracle-signal tests, stress tests, and
  nuisance-oracle ablations;
- tutorial updates and FDID Notes additions.

## P0: Reporting Notes-Style Estimands

### 1. Add first-class reporting for level comparisons

Problem:

The FDID Notes estimands are often not the whole curve itself, but a contrast
such as:

```text
theta(g1) - theta(g0)
```

Current state:

- `kernel` and `dml_flex` store level curves.
- `fdid_contrast()` now reports level comparisons from stored curve values.
- `kernel` stores bootstrap replicate curves when fit with
  `vartype = "bootstrap"`.
- `kernel` bootstrap fits store `mu_boot_vcov` and `delta_boot_vcov`, so exact
  grid contrasts can use `V11 + V22 - 2 V12`.
- `dml_flex` stores practical multiplier replicate curves when its multiplier
  band routine runs.
- When replicate curves are unavailable, pointwise endpoint-independent
  intervals remain available but are labeled approximate; conservative
  band-implied intervals are available when stored bands exist.
- Full curve covariance remains open for non-kernel curve paths.
- Plotting supports vertical curve markers through `diff.values` and
  `contrast.values`. It also now supports selected fixed-reference contrast
  coefficient plots through `plot(type = "contrast")`, with estimates and
  intervals delegated to `fdid_contrast()`.

Needed implementation:

- Complete as of 2026-05-29 for the practical reporting layer, with the
  first-class contrast coefficient plot added on 2026-05-31.
- Remaining research-level work: formalize whether the stored `dml_flex`
  multiplier process is sufficient for theorem-level contrast inference, and
  decide whether to add full curve covariance or influence-function storage for
  non-kernel curve paths.

Temporary reporting rule:

- Point estimate is acceptable.
- Independence-style SE should be labeled approximate and should not be the
  recommended default.
- A conservative band-implied CI can be constructed as:

```text
lower = lower_band(g1) - upper_band(g0)
upper = upper_band(g1) - lower_band(g0)
```

This is conservative and not efficient.

### 2. Add first-class reporting for derivative values

Problem:

Notes may target a derivative value:

```text
theta_prime(g0)
```

Current state:

- `kernel` stores `curve_event$delta_hat` and `curve_event$se_delta`.
- `dml_flex` stores `curve_event$delta_hat` and `curve_event$se_delta`.
- `dml_flex` also stores practical multiplier bands for the derivative curve.
- `fdid_derivative()` now reports fixed-`g` derivative values with replicate,
  band, or pointwise inference depending on the fitted object.
- `dml_incremental` targets a different scalar:

```text
E[partial_g mu(G, X)]
```

Needed implementation:

- Complete as of 2026-05-29 for the extraction helper.
- Remaining work: add derivative-step stability diagnostics and validation
  examples for interpolation under sparse support.

### 3. Decide final inference language for `dml_flex` bands

Current state:

- `dml_flex` stores practical multiplier simultaneous bands over
  signal-mapping residuals for non-BLP signal maps.
- `dml_flex` with `signal_map = "blp_spline"` stores grid covariance matrices
  and uses Gaussian max-t bands from those covariance matrices.
- These are useful diagnostics and are plotted.
- Non-BLP practical multiplier bands are not yet a finalized theorem-level
  uniform band for the full DML estimator.

Needed work:

- Keep the current language conservative.
- Decide whether to develop a formal multiplier process using cross-fitted
  orthogonal scores and mapping influence functions.
- If theory is not finalized, keep bands in the package but label them
  "practical signal-mapping bands."

## P0: Tutorial Reorganization

### 4. Rewrite the Continuous G tutorial into a main path plus appendices

Status as of 2026-05-31: complete for the public tutorial cleanup. The final
structure is target-first rather than appendix-heavy, with focused chapters for
continuous overview, kernel, DML, visualization, reporting/inference, and
sensitivity. Historical details below are retained for provenance.

The current tutorial remains too dense. The recommended organization from the
2026-05-18 plan is:

```text
3 Continuous G
  3.1 Continuous G in one page
  3.2 Load the local package and mortality data
  3.3 Inspect treatment support
  3.4 Scalar benchmarks
  3.5 Kernel level and derivative curves
  3.6 Flexible DML level and derivative curves
  3.7 Scalar incremental DML
  3.8 Learner and HD expansion choices
  3.9 One known-truth synthetic sandbox
  3.10 Method choice cheatsheet

4 Continuous G advanced examples
  4.1 Optional learner availability
  4.2 HD expansion design
  4.3 Learner grid for dml_flex
  4.4 Learner grid for dml_incremental
  4.5 Interaction-heavy synthetic DGP
  4.6 Nonsmooth synthetic DGP
  4.7 Heteroskedastic and skewed G | X
  4.8 What these stress tests teach

5 Continuous G technical appendix
  5.1 Estimands
  5.2 Orthogonal signals
  5.3 Cross-fitting
  5.4 Nuisance learners
  5.5 Density estimation
  5.6 Inference
  5.7 Validation tests and saved results
  5.8 Caveats and open TODOs
```

Rewrite rules:

- Start with the mortality example, not local installation details.
- Show support diagnostics before any curve estimator.
- Use one table for scalar benchmarks.
- Teach kernel curves before DML curves.
- Treat `dml_flex` as curve-first.
- Treat `dml_incremental` as scalar-average-derivative first.
- Move full learner grids, oracle validations, and stress tests out of the main
  rendered path.
- Keep default RStudio render fast.

### 5. Add a concise estimator-reporting section to the tutorial

Status as of 2026-05-31: complete in `tutorial/07-reporting-inference.Rmd`.

The tutorial should teach the difference between:

- reporting a whole curve;
- reporting a level value `theta(g0)`;
- reporting a level comparison `theta(g1) - theta(g0)`;
- reporting a derivative value `theta_prime(g0)`;
- reporting a scalar observed-population average derivative from
  `dml_incremental`.

Current package support and caveats should be explicit.

## P1: Inference And Diagnostics

### 6. Add curve covariance or replicate storage

Needed for formal curve contrasts:

- covariance matrix for `theta_hat(eval_g)`;
- covariance matrix for `delta_hat(eval_g)`;
- covariance between level and derivative curves if needed;
- multiplier or bootstrap replicate curves for flexible contrast inference.

Storage should be optional or lightweight enough not to bloat ordinary objects.

Status as of 2026-05-29:

- optional replicate storage is implemented for `kernel` bootstrap fits and
  `dml_flex` multiplier fits;
- full covariance matrix and influence-function storage remain open.

Status as of 2026-05-30:

- `kernel` bootstrap fits store `mu_boot_vcov` and `delta_boot_vcov`;
- kernel scalar summaries and exact-grid contrasts can use covariance-aware
  endpoint differences;
- `dml_flex` with `signal_map = "blp_spline"` stores `theta_vcov` and
  `delta_vcov`, so exact-grid and linearly interpolated level/derivative
  reporting can use covariance-aware inference;
- full covariance matrix and influence-function storage remain open for
  non-BLP `dml_flex` curve paths.

### 7. Add scalar DML multiplier/bootstrap inference

Current DML variance support:

- DML methods support `vartype = "robust"` only.
- Scalar DML robust SEs use cross-fitted orthogonal score variance and
  repeated-split dispersion.
- Scalar DML robust SEs are cluster-aware when clusters are supplied.
- As of 2026-05-30, scalar DML methods store influence values and
  score-multiplier diagnostics in `scalar_event`.
- `dml_inference = "score_multiplier"` uses the multiplier critical value for
  scalar CIs when available. `dml_boot` controls multiplier draws and
  `dml_multiplier` chooses `normal`, `wild`, or `exponential` weights.

Remaining work:

- Consider nonparametric DML refit bootstrap only if computationally
  acceptable.
- Validate new inference against oracle-signal simulations before documenting
  it as recommended.

### 8. Improve derivative-nuisance stability

Derivative targets are fragile because they require smooth nuisance estimates.

Needed work:

- Add user control over finite-difference step sizes.
- Store derivative-step diagnostics.
- Warn when derivative estimates change materially across step sizes.
- Recommend smooth learners for derivative-sensitive targets.
- Avoid presenting discontinuous learners as default for derivative targets.

### 9. Strengthen density diagnostics and warnings

Current diagnostics store density quantiles, minimum density, and low-density
shares.

Needed work:

- Warn when `s_hat(G | X)` is very small for many observations.
- Warn when `partial_g log s_hat(G | X)` has extreme tails.
- Warn when `eval_g` extends into sparse support.
- Recommend sensitivity checks across `density_method`.
- Add a compact density diagnostic plot for `dml_flex` and `dml_incremental`.

### 10. Keep conditional-density research active

Existing options:

- `residual_kde`;
- `location_scale`;
- `local_kde`.

Further options to consider:

- series conditional density;
- flexible distribution regression over `G`;
- quantile-based density approximation;
- local likelihood density;
- external density learners if stable and lightweight.

Any new density method must be tested through the oracle-signal, stress-test,
and nuisance-oracle ablation battery.

## P1: Validation

### 11. Turn validation scripts into a clean test battery

Current validation scripts:

```text
inst/oracle_signal_validation.R
inst/oracle_signal_stress_tests.R
inst/nuisance_oracle_ablation.R
inst/plot_validation_results.R
```

Needed organization:

- fast CI validation with low replication counts;
- local research validation with higher replication counts;
- frozen result files for Notes and tutorial figures;
- standardized result schemas;
- fixed seeds;
- plot generation that can be rerun without modifying package tests.

### 12. Add external benchmark checks

Recommended comparisons:

- compare `dml_plr` with a standard DoubleML PLR implementation on simple
  DGPs;
- compare `dml_binary` with standard IRM/AIPW implementations on binary-G DGPs;
- retain oracle validation as the primary validation for `dml_flex` and
  `dml_incremental`, because exact package equivalents may not exist.

## P1: Build And Release Hygiene

### 13. Keep generated artifacts out of package builds

Current `.Rbuildignore` was expanded, but this should remain checked after
every validation run.

Generated outputs that should not accidentally ship unless intentionally
frozen:

```text
tests/testthat/Rplots.pdf
inst/oracle_signal_validation_results/
inst/oracle_signal_stress_results/
inst/nuisance_oracle_ablation_results/
inst/validation_figures/
tutorial/_book/
tutorial/*_cache/
tutorial/*_files/
fdid.Rcheck/
```

### 14. Make command-line tutorial rendering reproducible

Current state:

- RStudio's bundled Quarto rendered the tutorial.
- Command-line package builds can fail if Pandoc is unavailable on `PATH`.

Needed work:

- document RStudio rendering;
- document command-line rendering;
- decide whether tutorial rendering belongs inside package checks or outside
  the package build;
- keep the default render free of long Monte Carlo jobs.

## P2: API And Helpers

### 15. Add method-comparison helpers

Useful applied helper:

```r
fdid_compare_continuous(
  s,
  tr_period,
  ref_period,
  methods = c("dml_plr", "kernel", "dml_flex", "dml_incremental"),
  ...
)
```

Output should include:

- scalar table;
- curve list;
- learner/density diagnostics;
- support diagnostics;
- warnings when estimates target different estimands.

### 16. Add support and diagnostic plotting helpers

Potential helpers:

```r
plot_g_support()
plot_dml_diagnostics()
plot_density_sensitivity()
plot_curve_sensitivity()
```

These should avoid requiring users to manually extract internal fields.

### 17. Consider refactoring the monolithic `fdid.R`

Current architecture keeps original methods inside `fdid.R` and continuous-G
methods in `kernel.R` and `dml.R`.

Future cleanup:

- split original method estimators into separate internal functions;
- reduce closure-heavy logic in `fdid.R`;
- add unit tests for internal helpers where feasible;
- keep the public API unchanged.

## P2: Documentation

### 18. Keep Notes, tutorial, and estimator reference synchronized

Every future estimator change should update:

- FDID Notes theory section;
- tutorial examples;
- `CONTINUOUS_G_ESTIMATOR_REFERENCE.md`;
- `CONTINUOUS_G_PROJECT_LOG.md`;
- package roxygen documentation if user-facing.

### 19. Add a method-choice cheatsheet to the tutorial

Recommended table:

| Need | Use first | Then check | Avoid saying |
|---|---|---|---|
| Quick scalar benchmark | `dml_plr` | OLS sensitivity | "This is the whole dose-response curve" |
| Smooth level curve | `kernel` | bandwidth/support | "Kernel removes all confounding concerns" |
| Flexible adjusted level curve | `dml_flex` | learner, HD expansion, density diagnostics | "The scalar summary is the main result" |
| Scalar marginal effect | `dml_incremental` | density diagnostics and fold sensitivity | "It is just the average plotted slope" |
| Rich nuisance structure | `dml_flex` plus HD expansion | optional learners and support | "More flexible is always better" |

## Completed And Archived Work

The following work is complete and preserved in the consolidated log:

- strict continuous-G input validation for `alpha`, `trim`, `eval_g`, `h0`,
  `boot`, `K_folds`, `K`, and `S`;
- repeated sample-size-aware validation after missing-data removal;
- cleanup of local install metadata from `DESCRIPTION`;
- `.Rbuildignore` expansion;
- cluster-aware scalar DML robust SEs;
- DML diagnostics;
- practical `dml_flex` multiplier bands;
- curve plotting with explicit level, derivative, and two-panel options;
- kernel curve plotting with pointwise versus simultaneous interval controls,
  interflex-style `show.uniform.CI`, support rug/separate-panel displays, and
  robust safeguards against labeling pointwise normal intervals as uniform
  bands;
- summary changes for DML diagnostics and curve-band method;
- the earlier Mapping Extensions, Density Extensions, and Summary Change issue
  list from `DML_ESTIMATOR_REVIEW_VERDICT.md`;
- xgboost 3.x compatibility through `xgb.train()`;
- continuous-G regression tests;
- mortality examples;
- HD expansion and learner sensitivity examples;
- known-truth synthetic examples;
- oracle-signal validation;
- stress tests varying `G | X` density and nonsmooth/interacting `mu(g, x)`;
- nuisance-oracle ablation tests;
- validation visualizations;
- FDID Notes updates;
- markdown consolidation into three project documents.

## Recommended Immediate Next Step

The reporting layer for Notes-style estimands was implemented on 2026-05-29.
The next priority is to validate and harden it:

1. Run oracle-signal checks that explicitly cover `fdid_contrast()` and
   `fdid_derivative()`.
2. Stress-test interpolation and sparse-support warnings.
3. Decide whether to add full curve covariance or influence-function storage
   for theorem-level contrast inference.
4. Keep the tutorial, FDID Notes, and estimator reference synchronized with
   any inference-language changes.

## Tutorial Audit Update Added 2026-05-31

The items below were first completed in the temporary tutorial copy
`tutorial-next-2026-05-30/` and were then migrated into the canonical
`tutorial/` source during the target-first cleanup:

- Added a focused kernel chapter covering support diagnostics, estimation-grid
  choices, level and derivative curves, pointwise versus bootstrap
  quantile-envelope bands, `show.uniform.CI`, `fdid_contrast()`,
  `fdid_derivative()`, and a user option checklist.
- Added a focused DML chapter covering target mapping, scalar PLR, flexible
  DML curves, scalar score versus score-multiplier inference, practical
  `dml_flex` bands, BLP-spline Gaussian max-\(t\) covariance bands,
  incremental average derivatives, and learner/design sensitivity.
- Expanded the DML chapter's conditional-density section so the three
  `density_method` choices (`residual_kde`, `location_scale`, `local_kde`) are
  described separately from `signal_map`; `blp_spline` is now described as a
  signal-mapping and covariance option, not a density method.
- Expanded Chapter 6, Visualization Options, with the continuous-G plotting
  controls for `curve`, `interval`, `show.uniform.CI`, `Xdistr`,
  `support.panel`, `show.eval_g`, marker values, support display styles, and
  base-graphics styling.
- Replaced the earlier manual fixed-reference contrast example in Chapter 6
  with the formal `plot(type = "contrast")` coefficient plot, while keeping
  selected inferential contrasts directed to `fdid_contrast()`.
- Added `PACKAGE_PLOTTING_LAYER_UPDATE_PLAN_2026-05-31.md`. The revised
  plotting sequence now puts a formal fixed-reference continuous-G contrast
  coefficient plot before package-wide `fdid_autoplot()` / ggplot work.
- Implemented the recommended fixed-reference contrast coefficient plot:
  `plot(type = "contrast")` now displays selected
  `mu(g_j) - mu(g_ref)` or `theta(g_j) - theta(g_ref)` estimates with 95%
  intervals computed through `fdid_contrast()`.
- Refined the contrast plot so it displays the reference row by default,
  supports `ref.g = "median"`, `"mean"`, `"min"`, `"max"`, or `"first"`,
  and supports quantile target selection through
  `target.type = "quantile"`.
- Rendered the full canonical tutorial book successfully into `docs/`.

Open tutorial tasks:

- Reviewer should inspect the completed fallback audit before ship.
- Keep the new `plot(type = "contrast")` examples synchronized with any future
  changes to `fdid_contrast()` inference labels or package-wide plotting APIs.
- After the package-wide plotting-layer decision is made, update the kernel
  and DML chapters if a `ggplot2`/`autoplot()` path is added.

## Full TODO Crosswalk Added 2026-05-19

This section was added after re-reading the active docs, archived markdowns,
and current R code. It preserves TODOs that were mentioned in chat or older
markdowns even when they overlap with the shorter prioritized list above.

### Code-Audit TODOs Added 2026-05-19

These come directly from the current R implementation audit.

- Decide whether `fdid()` should keep checking `estimatr` and `car` before
  continuous-G and DML dispatch. The current code requires them even when a DML
  path would not otherwise use them.
- Decide whether DML continuous-G dispatchers should include cluster column `c`
  in their initial complete-case filter. Kernel now includes the cluster column
  for bootstrap resampling when available.
- Add or document cluster-aware curve inference for `dml_flex`. Kernel
  bootstrap resampling is cluster-aware; analytical kernel sandwiches and
  `dml_flex` curve inference are not.
- Decide whether `dml_flex` diagnostics should report cluster counts when a
  cluster column exists. Current `dml_flex` diagnostics pass `cluster = NULL`.
- Keep an eye on runtime for large `boot` values in `dml_flex`; `boot` now
  controls both practical residual-multiplier draws and Gaussian max-t
  simulation draws.
- Optional future plot extension: add a `ggplot2`/`autoplot.fdid()` path that
  returns a reusable plot object while keeping the base graphics method
  dependency-free.
- DML plot Phase A is implemented for the base curve path: `dml_flex` plot data
  now carries DML-specific labels and the legend distinguishes practical
  signal-residual multiplier bands from BLP-spline Gaussian max-t covariance
  bands. Remaining DML plotting work is optional and broader: scalar DML
  object-returning plots, validated `dml_binary` overlap plot data, opt-in DML
  diagnostics, and any package-wide ggplot/autoplot path.
- Add a clearer label or warning that continuous-G `raw_means` are median-split
  descriptive displays only.
- Keep `historical_mds_2026-05-19/` ignored from package builds unless a
  release intentionally includes the archive.
- Add a quick documentation check after future code edits: every new object
  field in `R/kernel.R`, `R/dml.R`, `R/plot.R`, `R/summary.R`, and
  `R/fdid_prepare.R` should be reflected in
  `CONTINUOUS_G_ESTIMATOR_REFERENCE.md`.

### Historical Inference TODOs Preserved

- Keep current inference support explicit everywhere:
  original methods support robust/bootstrap/jackknife, kernel supports
  robust/bootstrap, and DML supports robust only.
- First-class practical contrast reporting for `theta(g1) - theta(g0)` was
  added through `fdid_contrast()`; full theorem-level efficient covariance or
  influence-function storage remains open for non-BLP maps and all nuisance
  stages. The narrower `signal_map = "blp_spline"` path now stores finite-grid
  level and derivative covariance matrices.
- Derivative-value extraction for `theta_prime(g0)` was added through
  `fdid_derivative()` with explicit interpolation control.
- Optional multiplier/bootstrap replicate curves are now stored for
  `dml_flex` multiplier fits and kernel bootstrap fits; kernel bootstrap curve
  covariance matrices are stored, and `dml_flex` BLP-spline fits store level
  and derivative grid covariance matrices.
- Scalar score multiplier inference is implemented for `dml_binary`,
  `dml_plr`, and `dml_incremental`; nonparametric DML refit bootstrap remains
  open.
- Decide whether nonparametric bootstrap for DML is computationally acceptable.
- Develop or clearly reject theorem-level uniform bands for `dml_flex`.
- Keep practical `dml_flex` bands labeled as signal-mapping bands until the
  full theory is finalized.
- Avoid reporting endpoint-difference SEs as formal contrast inference unless a
  stored curve covariance matrix or replicate-curve contrast is available.

### Historical Diagnostics TODOs Preserved

- Add automatic warnings when estimated `s_hat(G | X)` is near zero for many
  observations.
- Add automatic warnings when `partial_g log s_hat(G | X)` has extreme tails.
- Add warnings when `eval_g` extends into sparse treatment support.
- Add density-method sensitivity helpers and diagnostic plots.
- Add nuisance-quality diagnostics: outcome nuisance loss, treatment nuisance
  loss, density diagnostics, fold count, split count, split/fold failures, and
  orthogonal-score summaries.
- Store finite-difference step sizes for derivative nuisances.
- Compare derivative estimates across multiple finite-difference step sizes and
  warn when they are unstable.
- Warn when `local_kde` is requested after high-dimensional expansion.
- Expose conditional-density bandwidth controls separately from signal-to-`G`
  mapping bandwidth controls.
- Consider series density, distribution regression, quantile-based density,
  local likelihood density, and eventually density ensembles after diagnostics
  are stable.

### Historical Learner TODOs Preserved

- Keep documenting that derivative targets need smooth nuisance estimates.
- Recommend smooth learners for derivative-sensitive targets: `gam`,
  spline-expanded linear models, and regularized smooth-basis regressions.
- Treat tree, boosting, and neural learners as sensitivity checks for
  derivative targets unless derivative stability diagnostics are added.
- Consider smoothing fitted nuisance predictions before finite-difference
  derivatives for tree-based learners.
- Keep optional learner availability tables in the tutorial so skipped
  `ranger`, `grf`, `nnet`, or `xgboost` rows are not interpreted as theoretical
  method failures.

### Historical Validation TODOs Preserved

- Turn the validation scripts into a documented test battery with three tiers:
  fast CI validation, local research validation, and frozen evidence for Notes
  and tutorial figures.
- Standardize result schemas across oracle validation, stress tests, nuisance
  ablation, and validation figures.
- Keep fixed seeds and broad pass/fail tolerances that do not overfit to Monte
  Carlo noise.
- Add external benchmark checks where possible:
  compare `dml_plr` with a standard DoubleML PLR implementation and
  `dml_binary` with a standard IRM/AIPW implementation.
- Keep oracle validation as the main validation for `dml_flex` and
  `dml_incremental`, because direct package equivalents may not exist.
- Preserve the existing oracle-signal, stress-test, and nuisance-oracle
  ablation scripts and saved CSV/PNG outputs.

### Historical Tutorial TODOs Preserved

- Keep the Continuous G tutorial organized by user question and estimand, not
  by an exhaustive sequence of every method and option.
- Start the main tutorial with mortality, treatment support, and scalar
  benchmarks before advanced theory.
- Show support diagnostics before curve estimators and explain why the main
  grid avoids extreme support.
- Explain that a dense `eval_g` grid is for curve estimation and a small set of
  quantile landmarks is for display; the estimator is not partitioning `G` into
  five bins unless a binning method is explicitly used.
- Keep `kernel` and `dml_flex` curve-first; avoid leading with dynamic estimate
  tables for `dml_flex`.
- Keep `dml_incremental` scalar-first; do not plot it as a curve.
- Move `dml_binary` out of the main continuous-G flow or label it as a binary
  companion with a different estimand.
- Keep learner/HD/density/signal-map comparisons as sensitivity analysis, not
  the main empirical story.
- Add a "bad example" panel showing density-method sensitivity, especially
  `residual_kde` versus `location_scale` for `dml_incremental`.
- Keep full learner grids and heavy synthetic examples behind a render switch
  such as `FDID_TUTORIAL_HEAVY=true`.
- Include local RStudio workflow instructions: reinstall the local package,
  restart R, call `library(fdid)`, check `packageVersion("fdid")`, and render
  with RStudio's bundled Quarto when needed.
- Keep formulas in tutorial markdown compatible with the renderer and avoid
  TeX display problems.

### Historical Synthetic-Example TODOs Preserved

- Keep a smooth nonlinear known-truth example showing why scalar slopes are
  insufficient.
- Keep a nonlinear-confounding example showing why flexible covariate
  adjustment and HD expansion matter.
- Keep a heteroskedastic `G | X` example showing why density method choice is
  not a minor tuning detail.
- Keep a nonsmooth or kinked outcome example showing why derivative targets are
  fragile even when level curves are useful.
- Use truth overlays for level and derivative curves whenever possible.
- Avoid turning the synthetic sandbox into a second exhaustive tutorial; each
  synthetic example should teach one diagnostic lesson.

### Historical Build And Release TODOs Preserved

- Keep generated outputs out of package builds unless intentionally frozen:
  Rcheck folders, tarballs, tutorial build output, caches, validation result
  folders, validation figures, and `tests/testthat/Rplots.pdf`.
- Decide whether validation outputs under `inst/` are intentionally shipped as
  package evidence or kept only as local research artifacts.
- Keep command-line tutorial rendering reproducible, including Pandoc/Quarto
  availability.
- Keep default tutorial render fast and free of long Monte Carlo jobs.
- Update roxygen/man pages after user-facing API or argument changes.
- Keep `NEWS.md` synchronized with release-ready user-facing changes.

### Historical Documentation Synchronization TODOs Preserved

- Every future continuous-G estimator change should update:
  `CONTINUOUS_G_ESTIMATOR_REFERENCE.md`, `CONTINUOUS_G_PROJECT_LOG.md`,
  `TODO.md`, FDID Notes, tutorial files, and package help files when
  user-facing.
- Keep the memory/log markdown complete enough that the workflow can be
  reconstructed without reading archived markdowns.
- Preserve historical markdowns separately but treat the three active docs as
  the source of truth going forward.
