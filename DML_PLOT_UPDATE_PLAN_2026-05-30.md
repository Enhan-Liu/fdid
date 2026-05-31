# DML Plot Update Plan

Date: 2026-05-30

This plan adapts the kernel plot update plan to the DML estimators in `fdid`.
It is a plotting and documentation plan only. No estimator, inference
calculation, stored object field, tutorial, `NEWS.md`, or `HANDOFF.md` is
changed by this document.

Execution priority: run this DML plot plan before the package-wide ggplot /
plotting-layer plan. The first DML implementation pass should align
`dml_flex` base-graphics curve plotting with the kernel plot behavior that is
already implemented. Package-wide `ggplot2`, `fdid_autoplot()`, raw/dynamic
refactors, and tutorial plotting rewrites should remain deferred until this
DML-kernel alignment pass has been reviewed.

Implementation status, 2026-05-30: Phase A has been implemented in the
base-graphics path. The update extends DML curve plot-data metadata, improves
`dml_flex` band labels, adds focused tests, and updates reference/log/TODO/Rd
documentation. Phase B data-builder expansion and Phase C ggplot/autoplot
remain deferred.

The short answer is yes: most of the kernel plot design can migrate to DML,
but only after separating DML into scalar, curve, overlap, and diagnostic
plotting surfaces. The kernel plot plan is mostly a curve-plan; DML is not.

## Objective

Build a DML-aware plotting layer that:

- keeps `plot.fdid()` backward compatible and base-graphics first;
- reuses the kernel curve plot grammar for `dml_flex` where appropriate;
- does not turn scalar DML methods into fake curve methods;
- makes DML interval and band labels match the implemented inference exactly;
- adds a future object-returning ggplot path for DML scalar, curve, overlap,
  and diagnostic plots;
- makes DML diagnostics visible without cluttering ordinary estimator plots.

## What Transfers From The Kernel Plot Plan

The following kernel plot ideas transfer directly to `dml_flex`:

- canonical plot-data helpers before rendering;
- explicit interval controls:
  `interval = c("auto", "pointwise", "uniform", "both", "none")`;
- `show.uniform.CI` as an interflex-style control for simultaneous bands;
- `curve = c("auto", "level", "derivative", "both")`;
- support displays for observed `G`:
  `Xdistr = c("histogram", "hist", "density", "rug", "none")`;
- optional evaluation-grid ticks through `show.eval_g`;
- `diff.values` and `contrast.values` vertical markers;
- pointwise intervals as ribbons and simultaneous bands as visually separate
  boundaries;
- `xlim` filtering before drawing, not only coordinate clipping;
- optional `ggplot2` path that returns a reusable object.

These ideas should be implemented through shared curve plot data so base
graphics and future ggplot rendering use the same rows and labels.

## What Does Not Transfer Directly

The kernel plot plan cannot be copied wholesale because DML has multiple target
types:

- `dml_binary` is a scalar binary-treatment DML/AIPW-style contrast and stores
  propensity scores for overlap diagnostics;
- `dml_plr` is a scalar partially linear slope and does not have a dose-response
  curve;
- `dml_flex` is the only DML method with level and derivative curves over
  `eval_g`;
- `dml_incremental` is a scalar observed-population average derivative and
  must not be plotted as a fixed-`g` derivative curve.

The kernel bootstrap language also does not transfer:

- scalar DML score multiplier intervals are not row bootstrap intervals;
- `dml_flex` practical bands are residual-multiplier signal-mapping bands;
- `dml_flex` with `signal_map = "blp_spline"` can display finite-grid
  Gaussian max-t BLP-spline covariance bands;
- DML refit bootstrap remains unsupported.

## Current DML Plot State

Current `plot.fdid()` behavior relevant to DML:

- `type = "raw"` works for all DML fits through stored descriptive raw means;
- `type = "dynamic"` works for DML scalar summaries over time;
- `type = "overlap"` can be meaningful for `dml_binary` because it stores
  averaged propensity scores in `result$ps`;
- `type = "curve"` is available only for `dml_flex`, because it stores
  `result$curve_event`;
- `dml_plr` and `dml_incremental` correctly remain scalar methods with no curve
  plot;
- `dml_flex` base curve plots already support level, derivative, both,
  pointwise intervals, stored bands, support displays, evaluation-grid ticks,
  and marked `G` values.

Current limitations:

- the curve plot-data helper is shared with kernel but does not yet store all
  DML-specific labels users may need, such as `signal_map`,
  `dml_metadata$target_estimand`, `dml_metadata$inference_method`,
  `curve_event$practical_band_method`, `gaussian_crit_theta`, or
  `gaussian_crit_delta`;
- scalar DML intervals are visible in tables and summaries but do not have a
  dedicated scalar ggplot/object-returning display;
- `dml_binary` overlap plots still use the older binary histogram path, not a
  reusable overlap plot-data layer;
- `dml_diagnostics` are stored but do not have first-class diagnostic plots;
- plot labels can say "Uniform band" generically instead of distinguishing
  practical signal-mapping bands from BLP-spline Gaussian max-t bands.

## DML Plotting Surfaces

### `dml_binary`

Primary displays:

- scalar event estimate and confidence interval;
- dynamic scalar estimates;
- propensity-score overlap diagnostic using `result$ps` and binary `G`;
- optional score/influence diagnostic after the data layer stabilizes.

Do not add:

- continuous-G curve plot;
- continuous support histogram as if `G` were a dose.

Inference label:

- analytical score interval by default;
- score-multiplier interval only when `scalar_event$scalar_inference` is
  `"score_multiplier"`;
- never label as bootstrap.

### `dml_plr`

Primary displays:

- scalar partially linear slope estimate;
- dynamic scalar slope estimates;
- optional coefficient-style plot for a single object;
- scalar comparison through `plot.fdid_list()` or future
  `fdid_autoplot(fdid_list(...))`.

Do not add:

- curve plot;
- support strip implying a nonlinear dose-response curve.

Inference label:

- analytical score interval by default;
- optional score-multiplier interval if requested and stored;
- no DML refit bootstrap label.

### `dml_flex`

Primary displays:

- level curve `theta(g)`;
- derivative curve `theta'(g)`;
- two-panel level plus derivative display;
- observed `G` support;
- optional `diff.values` and `contrast.values`;
- optional diagnostics for signal, density, correction terms, and sparse
  support.

Inference labels:

- pointwise interval: second-stage mapping pointwise interval;
- non-BLP simultaneous band: practical residual multiplier over cross-fitted
  signal-mapping residuals;
- `signal_map = "blp_spline"` simultaneous band: finite-grid Gaussian max-t
  BLP-spline covariance band when stored covariance is usable;
- never describe practical bands as final theorem-level DML uniform bands;
- never describe these bands as DML refit bootstrap bands.

### `dml_incremental`

Primary displays:

- scalar observed-population average derivative estimate;
- dynamic scalar estimates if produced;
- diagnostic plots for signal, conditional density, density derivative, and
  correction terms;
- scalar comparison through `fdid_list`.

Do not add:

- `plot(type = "curve")`;
- fixed-`g` derivative interpretation.

Inference label:

- analytical score interval by default;
- optional score-multiplier interval if requested and stored;
- target label should say observed-population average derivative.

## Proposed User-Facing API

Keep existing base calls:

```r
plot(fit_flex, type = "curve", curve = "level")
plot(fit_flex, type = "curve", curve = "derivative")
plot(fit_flex, type = "curve", curve = "both", interval = "both")
plot(fit_binary, type = "overlap")
plot(fit_plr, type = "dynamic")
plot(fit_incremental, type = "dynamic")
```

Do not make these calls valid:

```r
plot(fit_plr, type = "curve")
plot(fit_incremental, type = "curve")
```

Add a future object-returning path:

```r
fdid_autoplot(fit, type = "scalar")
fdid_autoplot(fit, type = "curve", curve = "both", interval = "both")
fdid_autoplot(fit, type = "overlap")
fdid_autoplot(fit, type = "diagnostics", diagnostic = "dml_density")
fdid_autoplot(fit, type = "diagnostics", diagnostic = "dml_signal")
```

Suggested semantics:

- `type = "scalar"` displays the event-period scalar estimate and interval for
  scalar methods;
- `type = "curve"` is accepted only when `curve_event` exists;
- `type = "overlap"` is accepted only when valid binary propensity-score data
  exist;
- `type = "diagnostics"` is opt-in and should never replace the main estimate
  plot;
- scalar and curve plot labels are generated from `dml_metadata` when
  available.

## Proposed Internal Data Layer

Create or extend a shared plotting data layer:

```r
.fdid_plot_data_dml_scalar(x, alpha = NULL)
.fdid_plot_data_dml_curve(x, curve, interval, xlim = NULL)
.fdid_plot_data_dml_overlap(x, xlim = NULL, bins = 30L)
.fdid_plot_data_dml_diagnostics(x, diagnostic, xlim = NULL)
```

Scalar DML plot data should include:

```text
method
target_estimand
target_scale
estimate
std_error
conf_low
conf_high
inference_method
scalar_inference
scalar_multiplier_crit
n_multiplier
n_splits
n_clusters
```

Curve DML plot data should extend the current curve data with:

```text
signal_map
density_method
target_estimand
inference_scope
pointwise_method
band_method
practical_band_method
gaussian_crit
simultaneous_crit
curve_vcov_available
band_label
band_caveat
```

Overlap DML plot data should include:

```text
score
group
bin_start
bin_end
count
density
method
```

Diagnostic DML plot data should include a long format:

```text
diagnostic
component
value
statistic
method
density_method
signal_map
```

Recommended first diagnostic components:

- conditional density quantiles;
- `d_log_density` quantiles;
- correction-term quantiles;
- signal quantiles;
- propensity-score overlap for `dml_binary`;
- scalar influence values when stored and safe to display.

## Rendering Plan

### Base Graphics

Base graphics should remain the default:

- preserve `plot.fdid()` return value `invisible(NULL)`;
- preserve `plot.fdid(type = "curve")` only for `kernel` and `dml_flex`;
- improve DML-specific legend labels in the existing curve path;
- keep scalar methods on raw, dynamic, overlap, list, and future scalar
  displays;
- avoid adding new hard dependencies.

Near-term base improvement:

- for `dml_flex`, change generic "Uniform band" legend text into method-aware
  labels:
  - `Practical band (signal residual multiplier)`;
  - `Gaussian max-t band (BLP spline covariance)`.

### ggplot Path

The future ggplot path should:

- live behind `fdid_autoplot()` while `ggplot2` is only suggested;
- return visible `ggplot` objects;
- share the same canonical data used by base graphics;
- draw scalar DML estimates with `geom_pointrange()`;
- draw `dml_binary` overlap with histogram/density layers;
- draw `dml_flex` pointwise intervals and simultaneous bands as separate
  layers;
- draw diagnostics as opt-in plots, not as default estimator plots;
- use facets for `curve = "both"` in the first pass.

## Inference Display Rules

Scalar DML:

- `scalar_inference = "analytical_score"`:
  label as analytical score interval;
- `scalar_inference = "score_multiplier"`:
  label as score-multiplier critical-value interval;
- no scalar DML plot should say bootstrap unless a future refit bootstrap is
  actually implemented.

`dml_flex` curve:

- pointwise ribbon:
  label as pointwise second-stage mapping interval;
- practical non-BLP band:
  label as practical signal-mapping residual multiplier band;
- BLP-spline Gaussian band:
  label as finite-grid Gaussian max-t BLP-spline covariance band;
- if `theta_vcov` or `delta_vcov` is available, plot metadata may note that
  reporting helpers can use covariance-aware contrasts;
- no plot should call non-BLP practical bands theorem-level uniform DML bands.

`dml_incremental`:

- scalar label must say observed-population average derivative;
- diagnostics may show density sensitivity but must not imply a fixed-`g`
  derivative curve.

## Relationship To `interflex`

Borrow:

- canonical plot-data frame before rendering;
- pointwise ribbon plus separate uniform-band boundary;
- support display below the main curve;
- marked moderator/treatment values;
- ggplot object return.

Do not borrow blindly:

- interflex's continuous-treatment DML estimand;
- a monolithic plotting function;
- labels that imply FDID `dml_flex` is a PLR moderator-CME estimator;
- uniform-band language that is stronger than FDID's stored DML inference.

## Documentation Plan

When implementation begins, update:

- `R/plot.R` roxygen and `man/plot.fdid.Rd`;
- new `R/plot-data.R` and `R/plot-ggplot.R` roxygen if created;
- `CONTINUOUS_G_ESTIMATOR_REFERENCE.md`;
- `CONTINUOUS_G_PROJECT_LOG.md`;
- `TODO.md`;
- the package-wide plotting plan.

Do not update in the DML plot planning pass:

- tutorials;
- `NEWS.md`;
- `HANDOFF.md`.

FDID Notes should be updated only if the plot layer changes how uncertainty is
described in the theory notes. Object-returning plotting is package
documentation, not theory.

## Tests

Structural tests:

- scalar DML plot data includes target and inference labels;
- scalar DML score-multiplier fits report multiplier metadata;
- `dml_plr` and `dml_incremental` still reject curve plotting;
- `dml_flex` curve data includes signal-map and band-method labels;
- `dml_flex` BLP-spline curve data labels Gaussian max-t BLP covariance bands;
- non-BLP `dml_flex` curve data labels practical residual-multiplier bands;
- `dml_binary` overlap data validates binary `G` and stored propensity scores;
- diagnostics data handles missing diagnostics gracefully.

Base graphics smoke tests:

- `plot(dml_flex, type = "curve", curve = "level")`;
- `plot(dml_flex, type = "curve", curve = "derivative")`;
- `plot(dml_flex, type = "curve", curve = "both")`;
- `plot(dml_binary, type = "overlap")`;
- scalar DML methods use `type = "dynamic"` without error.

ggplot tests, when implemented:

- skip cleanly if `ggplot2` is unavailable;
- `fdid_autoplot(dml_plr, type = "scalar")` returns a `ggplot`;
- `fdid_autoplot(dml_incremental, type = "scalar")` returns a `ggplot`;
- `fdid_autoplot(dml_flex, type = "curve")` returns a `ggplot`;
- expected layer data include pointwise and band layers separately;
- DML diagnostic plots return ggplot objects and have expected diagnostic
  labels.

## Acceptance Criteria

The immediate DML-kernel alignment pass is complete when:

- `dml_flex` curve plots clearly distinguish pointwise intervals from stored
  simultaneous bands;
- `dml_flex` BLP-spline Gaussian max-t bands and practical residual-multiplier
  bands have different labels;
- `dml_incremental` is never plotted as a fixed-`g` derivative curve;
- `dml_plr` and `dml_binary` are not presented as continuous-G curve methods;
- `dml_flex` supports the same base curve controls as the current kernel plot:
  `curve`, `interval`, `show.uniform.CI`, `Xdistr`, `support.panel`,
  `show.eval_g`, `diff.values`, `contrast.values`, and `xlim` filtering;
- base plot behavior and return semantics remain backward compatible;
- documentation states that scalar score multiplier inference is not bootstrap;
- documentation states that practical `dml_flex` bands are not final
  theorem-level DML uniform bands;
- tests cover DML curve plot-data labels, DML method restrictions, inference
  labels, and base plotting smoke behavior.

The broader DML plotting plan is complete only after later passes add scalar
DML object-returning plots, validated `dml_binary` overlap plot data, opt-in
DML diagnostics, and any ggplot/autoplot path chosen by the package-wide plan.

## Implementation Order

### Phase A: Immediate DML-Kernel Base Plot Alignment

1. Extend the current curve plot data helper with DML-specific metadata labels.
2. Improve base `dml_flex` curve legend labels for practical versus BLP bands.
3. Confirm `dml_flex` uses the same base curve controls as kernel:
   `curve = "level"`, `"derivative"`, `"both"`;
   `interval = "auto"`, `"pointwise"`, `"uniform"`, `"both"`, `"none"`;
   `show.uniform.CI`; `Xdistr`; `support.panel`; `show.eval_g`;
   `diff.values`; `contrast.values`; and `xlim`.
4. Keep `plot.fdid()` return semantics unchanged: base graphics and
   `invisible(NULL)`.
5. Add or update focused tests for DML curve labels, DML curve controls, and
   scalar-DML curve rejection.
6. Update roxygen/Rd, estimator reference, project log, TODO, and this plan.
7. Run continuous-G tests, plotting tests, and then the full package test
   suite.

### Phase B: Optional DML-Specific Data Builders

After Phase A is reviewed, consider:

1. Add internal scalar DML plot-data builder.
2. Add internal binary overlap plot-data builder.
3. Add internal DML diagnostic plot-data builder.
4. Add structural tests for these data builders.

### Phase C: Optional ggplot/autoplot DML Path

Run this only if the package-wide plotting plan is activated:

1. Add `fdid_autoplot()` scalar DML ggplot path.
2. Add `fdid_autoplot()` DML curve ggplot path.
3. Add `fdid_autoplot()` DML overlap path.
4. Add opt-in DML diagnostic plots.

## Open Decisions

1. Should single-object scalar DML plots be `type = "scalar"` or
   `type = "estimate"`?
2. Should DML diagnostics live under `fdid_autoplot(type = "diagnostics")` or
   separate helpers such as `plot_dml_diagnostics()`?
3. Should DML scalar influence values be plotted by default when available, or
   only behind an explicit diagnostic option?
4. Should `plot.fdid()` get a base scalar display for single objects, or should
   scalar displays remain summaries plus `fdid_list`/ggplot paths?
5. Should `plot(engine = "ggplot")` be added after `fdid_autoplot()`, or should
   ggplot use stay in a separate object-returning function?
