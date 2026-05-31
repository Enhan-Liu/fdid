# Package-Wide Plotting Layer Update Plan

Date: 2026-05-30

This document records a package-wide plan for modernizing the plotting layer in
`fdid`. It is a planning artifact only. No R implementation is changed by this
document.

2026-05-30 DML adaptation note: this package-wide plan is now paired with
`DML_PLOT_UPDATE_PLAN_2026-05-30.md`. The global plan should treat DML as a
four-surface plotting family: scalar estimates, binary overlap, `dml_flex`
curves, and opt-in diagnostics. It should not treat all DML methods as
continuous-G curve methods.

Sequencing decision: do not execute this package-wide plotting layer plan
before the DML plot plan. The next plotting code pass should execute
`DML_PLOT_UPDATE_PLAN_2026-05-30.md` Phase A first, with the limited goal of
aligning `dml_flex` base-graphics curve plotting with the already implemented
kernel curve plot. This package-wide ggplot/autoplot plan should be revisited
only after that DML-kernel alignment pass is implemented and reviewed.

Implementation note, 2026-05-30: the DML Phase A base-graphics alignment pass
has been implemented. The package-wide ggplot/autoplot plan remains deferred
pending a later explicit decision.

The plan is based on a full package-level read of the current `fdid-main`
plotting surface, continuous-G estimator outputs, package markdowns, FDID
Notes, tutorials, and the local `interflex-master` plotting implementation.

## StatsClaw Workflow Status

Workflow: Docs Only.

Applied sequence:

1. Planner pass: inventory the existing package plotting surface and identify
   implementation boundaries.
2. Explorer pass: read `interflex-master` plotting code and summarize reusable
   design patterns.
3. Scriber pass: write this dated package-wide plotting plan.
4. DML adaptation pass: split DML plotting into scalar, overlap, curve, and
   diagnostic surfaces.
5. Reviewer pass: review this plan for completeness before implementation.

This pass intentionally does not edit R code, tests, tutorials, `NEWS.md`, or
`HANDOFF.md`.

## Objective

Add a coherent, package-wide plotting layer that:

- keeps the current base-graphics `plot.fdid()` and `plot.fdid_list()` behavior
  backward compatible;
- adds a reusable plot-data layer shared by base graphics and any future
  ggplot path;
- adds optional ggplot/autoplot-style plot objects for downstream customization;
- covers every current plotting target in the package, not only the
  continuous-G kernel curve;
- covers DML-specific plotting targets without misclassifying scalar DML
  estimators as curve estimators;
- makes interval semantics explicit in every plot that displays uncertainty;
- makes DML scalar score, score-multiplier, practical curve-band, and
  BLP-spline Gaussian max-t labels explicit;
- follows the strongest reusable ideas from `interflex` without copying its
  monolithic plotting implementation.

The update should be a rendering-layer modernization. It should not change any
estimator, inference calculation, stored estimate, or default scientific
interpretation.

## Package Inventory

Current exported user workflow:

- `fdid_prepare()` prepares panel data in the package's wide format.
- `fdid_expand_covariates()` expands prepared covariates.
- `fdid()` estimates the selected FDID method.
- `print.fdid()` and `summary.fdid()` summarize results.
- `plot.fdid()` draws a single result object.
- `fdid_list()` bundles multiple `fdid` objects.
- `plot.fdid_list()` compares scalar estimates across bundled objects.
- `fdid_contrast()` and `fdid_derivative()` report continuous-G curve
  contrasts and fixed-G derivative values.

Current registered S3 plotting methods:

- `plot.fdid`
- `plot.fdid_list`

No `autoplot` method is registered. `ggplot2` is not currently listed in
`Imports` or `Suggests`.

Current estimation methods that can feed plotting:

- original binary/discrete methods: `ols1`, `ols2`, `did`, `ebal`, `ipw`,
  `aipw`;
- continuous-G kernel method: `kernel`;
- DML methods: `dml_binary`, `dml_plr`, `dml_flex`, `dml_incremental`.

Current common output fields:

- `est`: scalar pre/event/post summaries;
- `dynamic`: time-indexed dynamic estimates with confidence intervals;
- `raw_means`: descriptive group mean outcomes over time;
- `tr_period`, `ref_period`, `entire_period`, `times`;
- `method`, `vartype`, `alpha`, `target.pop`;
- `G`: observed baseline factor;
- `ps`: propensity scores where applicable;
- `curve_event`: continuous-G curve fields for `kernel` and `dml_flex`;
- `curve_dynamic`: dynamic curve fields for `kernel`;
- `scalar_event`: scalar DML score, influence, split, and multiplier fields;
- `incremental_event`: scalar incremental DML fields;
- `dml_metadata`: DML target, scale, population, and inference labels;
- `dml_diagnostics`: DML signal, density, derivative, correction, and cluster
  diagnostics.

Current documentation-only plotting surfaces:

- `fdid-main/vignettes/fdid.Rmd` uses `plot(result, type = "raw")` and
  `plot(result, type = "dynamic")` in the package vignette;
- `fdid-main/tutorial/04-visualization.Rmd` and the mirrored
  `tutorial-full-continuous-g/04-visualization.Rmd` contain many base
  `plot.fdid()` examples for raw, dynamic, overlap, and curve displays;
- `fdid-main/tutorial/03-continuous.Rmd` and its full tutorial mirror include
  a tutorial-only helper, `plot_curve_overlay()`, for overlaying stored curve
  outputs and truth curves in synthetic examples.

These tutorial-only plotting helpers are not package API, but they should be
reviewed during the tutorial rewrite and either replaced by the new plotting
layer or kept as clearly local tutorial code.

## Current Plotting Surface

### `plot.fdid(type = "raw")`

Source object: `x$raw_means`.

Current behavior:

- uses base graphics;
- plots descriptive raw outcome means over time by group;
- supports connected or point-only displays through `connected`;
- shades treatment periods through `shade_periods`;
- displays confidence intervals if `CI_Lower` and `CI_Upper` exist;
- relies on `RColorBrewer` for palette selection;
- returns `invisible(NULL)`.

Important interpretation:

- for continuous-G methods, `raw_means` are median-split descriptive displays;
- they are not the continuous-G estimator and should not be styled as an
  adjusted curve estimate.

### `plot.fdid(type = "dynamic")`

Source object: `x$dynamic`.

Current behavior:

- uses base graphics;
- plots dynamic estimates indexed by time or event period;
- supports connected or point-only displays;
- displays confidence intervals when available;
- draws a zero reference line;
- shades treatment periods;
- returns `invisible(NULL)`.

Important interpretation:

- dynamic rows are scalar summaries at each time period;
- for continuous-G curve methods, dynamic scalar estimates should remain
  distinct from full curve plots.

### `plot.fdid(type = "overlap")`

Source objects: `x$ps`, `x$G`.

Current behavior:

- uses base graphics;
- available for methods that store propensity scores;
- draws treated and control propensity score histograms;
- currently assumes binary `G` by plotting `G == 0` and `G == 1` groups;
- returns `invisible(NULL)`.

Important interpretation:

- this is a binary-treatment overlap diagnostic;
- it should not be reused as the continuous-G support diagnostic without a
  separate continuous support data builder;
- the new plot-data builder should explicitly validate binary `G` before
  constructing overlap data.

### `plot.fdid(type = "curve")`

Source objects: `x$eval_g`, `x$curve_event`, `x$G`.

Current behavior after the 2026-05-30 kernel plot update:

- uses base graphics;
- available for `method = "kernel"` and `method = "dml_flex"`;
- supports `curve = "auto"`, `"level"`, `"derivative"`, and `"both"`;
- supports `interval = "auto"`, `"pointwise"`, `"uniform"`, `"both"`, and
  `"none"`;
- supports `show.uniform.CI`;
- uses stored `eval_g` rows and does not re-estimate, densify, or interpolate;
- filters curve rows and support overlays by `xlim` before drawing;
- draws pointwise intervals as ribbons when available;
- draws stored simultaneous bands as dashed boundary lines when available;
- distinguishes kernel bootstrap percentile pointwise intervals from
  bootstrap quantile-envelope bands;
- avoids showing robust pointwise normal intervals as uniform bands;
- supports `Xdistr = "histogram"`, `"hist"`, `"density"`, `"rug"`, and
  `"none"`;
- supports `support.panel = "embedded"` or `"separate"`;
- supports `show.eval_g`, `diff.values`, and `contrast.values`;
- exposes base graphics style controls for curve lines, ribbons, bands, and
  support overlays;
- returns `invisible(NULL)`.

Important interpretation:

- the continuous-G curve plot is a display method over stored curve estimates;
- the support layer is descriptive and restricted to the displayed window;
- kernel bootstrap simultaneous bands are finite-grid bootstrap
  quantile-envelope bands, not studentized sup-t bands;
- `dml_flex` simultaneous bands remain practical signal-mapping bands unless
  the theory layer is later upgraded.

### `plot.fdid_list()`

Source objects: each element's `obj$est$event`.

Current behavior:

- uses base graphics;
- compares scalar event estimates across a list of `fdid` objects;
- supports vertical and horizontal orientation through `vertical`;
- includes `vartype` in labels for bootstrap and jackknife results;
- returns `invisible(x)`.

Important interpretation:

- this is a scalar estimate comparison plot, not a curve comparison plot;
- it should be modernized alongside `plot.fdid()` so the package has one
  coherent plotting vocabulary.

## DML-Specific Plotting Surface

DML methods need an explicit plotting split:

- `dml_binary`: scalar event/dynamic estimates plus binary propensity-score
  overlap through stored `ps`;
- `dml_plr`: scalar event/dynamic partially linear slope estimates, with no
  curve plot;
- `dml_flex`: level and derivative curves through `curve_event`, with
  pointwise mapping intervals and method-specific stored bands;
- `dml_incremental`: scalar observed-population average derivative, with no
  fixed-`g` derivative curve plot;
- DML diagnostics: optional views of stored signal, density, density
  derivative, correction, propensity, influence, split, or cluster summaries.

The DML plot layer should read `dml_metadata` when available. Plot titles,
subtitles, legends, and data columns should distinguish:

- analytical score intervals;
- scalar score-multiplier intervals;
- practical `dml_flex` residual-multiplier signal-mapping bands;
- finite-grid Gaussian max-t BLP-spline covariance bands;
- unsupported DML refit bootstrap.

The package-wide ggplot plan should therefore include a single-object scalar
plot path for scalar DML methods, not only `plot.fdid_list()` comparisons.

## Current Limitations

The current plotting layer works, but it has structural limits:

- plot methods are mostly side-effect-only base graphics calls;
- users cannot easily customize, compose, or save returned plot objects;
- raw, dynamic, overlap, curve, and list plots do not share a common plot-data
  contract;
- only continuous-G curve plots currently have a reusable internal data helper;
- aesthetic options are not consistently available across plot types;
- interval semantics are clearer for curve plots than for raw, dynamic, and
  list plots;
- support and density diagnostics are not first-class plot-data objects;
- scalar DML methods do not yet have a single-object scalar plotting data
  contract;
- DML diagnostics are stored but not exposed through opt-in diagnostic plots;
- `dml_flex` curve labels should distinguish practical residual-multiplier
  bands from BLP-spline Gaussian max-t bands;
- `theme.bw` exists only as a reserved compatibility argument in the base
  curve path;
- there is no package-wide convention for returning plot data, ggplot objects,
  or saved files;
- test coverage for plotting is necessarily more smoke-test oriented because
  base graphics is difficult to inspect structurally.

## Lessons From `interflex`

The following `interflex` patterns are worth adopting.

### Canonicalize Plot Data Early

`interflex` normalizes estimator outputs into stable plotting tables before
rendering. FDID should do the same through named columns, not positional column
assumptions.

Recommended FDID canonical columns:

- time/dynamic plots: `time`, `estimate`, `std.error`, `conf.low`,
  `conf.high`, `period_role`, `method`, `vartype`;
- raw plots: `time`, `group`, `mean`, `conf.low`, `conf.high`, `display_role`;
- overlap plots: `score`, `group`, `density`, `bin_start`, `bin_end`, `count`;
- continuous-G curve plots: `g`, `estimate`, `std.error`, `pointwise_lower`,
  `pointwise_upper`, `band_lower`, `band_upper`, `curve`, `interval_method`,
  `band_method`, `target_estimand`, `signal_map`, `density_method`,
  `band_label`, `band_caveat`;
- DML scalar plots: `method`, `target_estimand`, `target_scale`, `estimate`,
  `std.error`, `conf.low`, `conf.high`, `inference_method`,
  `scalar_inference`, `scalar_multiplier_crit`, `n_multiplier`, `n_splits`,
  `n_clusters`;
- DML diagnostic plots: `diagnostic`, `component`, `value`, `statistic`,
  `method`, `density_method`, `signal_map`;
- support plots: `g`, `support_type`, `density`, `bin_start`, `bin_end`,
  `count`;
- list plots: `label`, `method`, `vartype`, `estimate`, `conf.low`,
  `conf.high`, `index`.

### Return Plot Objects

`interflex` returns ggplot-compatible objects, and optionally returns panel
lists. FDID should add this capability without replacing base graphics as the
default S3 `plot()` behavior in the first implementation pass.

Recommended FDID convention:

- `plot.fdid()` keeps base graphics as the default and preserves existing
  return behavior;
- a new object-returning path returns ggplot objects;
- plot-data helpers are separately testable;
- object-returning functions are preferred for new examples and future
  tutorials.

### Separate Pointwise And Uniform Layers

`interflex` treats pointwise intervals and uniform bands as separate geoms.
FDID should use the same conceptual separation:

- pointwise intervals: ribbon or errorbar;
- simultaneous bands: dashed boundary lines, a different ribbon style, or a
  clearly labeled second layer;
- no plot should imply that pointwise normal intervals are uniform bands;
- curve-band labels should track stored inference metadata.

### Put Support Information Below The Estimate Scale

`interflex` rescales histograms or density displays into a lower strip under
the estimate curve. FDID already has this behavior for base curve plots. The
new package-wide layer should generalize it:

- support strips for continuous-G curves;
- overlap density/histogram layers for binary propensity scores;
- optional support diagnostics for DML density-risk plots;
- x-window filtering before rendering, not only coordinate clipping.

### Preserve Marked Values

`interflex` has `diff.values` for highlighting moderator values. FDID already
added `diff.values` and `contrast.values` for continuous-G curves. The new
plotting layer should preserve this idea in the ggplot path and eventually
connect it to `fdid_contrast()` and `fdid_derivative()` examples.

### Use Explicit Layer-Oriented Style Arguments

`interflex` exposes colors, line widths, interval alpha, support alpha, theme
choice, grid controls, labels, and output dimensions. FDID should collect these
into consistent defaults rather than scattering hard-coded values in each plot
branch.

### Filter Data Before Drawing

`interflex` filters estimator tables, density vectors, histogram bins, and
support displays before constructing geoms. FDID should adopt this as a
package-wide rule:

- `xlim` filters plot data first;
- only then should the renderer apply coordinate limits;
- tests should inspect filtered plot data.

## What Not To Copy From `interflex`

FDID should not copy these implementation patterns:

- one large branch-heavy plotting function;
- positional column selection when named columns are available;
- large helper functions hidden inside the main plotting method;
- fragile layout logic that depends on already-rendered layer scales;
- scattered hard-coded colors and alpha values;
- an implicit assumption that every estimator shares the same uncertainty
  columns.

FDID should instead use small internal data builders plus renderer-specific
layer helpers.

## Design Principles

1. Backward compatibility comes first.

   Existing calls such as `plot(fit, type = "dynamic")` and
   `plot(fit, type = "curve")` should continue to work and draw base graphics.

2. Estimation and plotting stay separate.

   Plotting should never re-estimate curves, densify sparse `eval_g` grids, or
   silently interpolate estimates unless a future public helper explicitly says
   so.

3. Plot data is the contract.

   Every plot type should have a testable data builder. Base graphics and
   ggplot renderers should consume the same canonical data.

4. Inference labels follow stored fields.

   The plot layer should display uncertainty exactly as stored by the
   estimator. It should not invent a stronger inference claim.

5. Optional dependencies should be deliberate.

   The first implementation should avoid forcing new plotting dependencies
   unless the public API requires them.

6. Package-wide consistency should not erase plot-specific meaning.

   Raw descriptive means, dynamic scalar estimates, propensity overlap,
   continuous-G curves, and fdid-list scalar comparisons have different
   interpretations. A common style system should not make them look like the
   same estimand.

## Proposed User API

### Phase-1 Public Function

Add an exported object-returning helper:

```r
fdid_autoplot(
  object,
  type = c("raw", "dynamic", "scalar", "overlap", "curve", "diagnostics"),
  ...
)
```

Rationale:

- avoids hard-registering an `autoplot` S3 method before deciding dependency
  policy;
- can require `ggplot2` at runtime with a clear message;
- gives users a stable ggplot-returning entry point;
- can later delegate to `autoplot.fdid()` if `ggplot2` moves to `Imports`.

DML-specific semantics:

- `type = "scalar"` displays a single object's event-period scalar estimate
  and interval, using `dml_metadata` and `scalar_event` labels when available;
- `type = "curve"` is valid for `dml_flex` and other methods with
  `curve_event`, but remains invalid for `dml_plr` and `dml_incremental`;
- `type = "overlap"` is valid for binary propensity-score objects such as
  `dml_binary`;
- `type = "diagnostics"` is opt-in and should never replace the primary
  scalar or curve display.

### Later S3 API

After dependency policy is settled, add:

```r
autoplot.fdid(object, ...)
autoplot.fdid_list(object, ...)
```

Possible registration choices:

- if `ggplot2` is moved to `Imports`, register S3 methods in `NAMESPACE`;
- if `ggplot2` remains in `Suggests`, keep `fdid_autoplot()` as the primary
  exported function and avoid fragile S3 registration.

### Optional `plot()` Engine Argument

In a later compatibility pass, consider:

```r
plot(x, type = "curve", engine = c("base", "ggplot"), ...)
```

Recommended semantics:

- default `engine = "base"` keeps current behavior;
- `engine = "ggplot"` builds the ggplot object, prints it, and returns it
  invisibly;
- `fdid_autoplot()` returns the ggplot object visibly and does not print unless
  the normal R console does so.

This lets base users keep old behavior and lets object-oriented plotting users
choose a modern path explicitly.

### Plot Data Helper

Consider an internal helper first:

```r
.fdid_plot_data(x, type, ...)
```

and export a public helper only after the schema stabilizes:

```r
fdid_plot_data(x, type, ...)
```

Reasons to keep it internal at first:

- schemas may change during the first implementation;
- exported data APIs become long-term maintenance commitments;
- tests can still cover internal helpers.

## Proposed Internal Architecture

### New File: `R/plot-data.R`

Purpose: canonical data extraction and validation.

Proposed helpers:

```r
.fdid_plot_data_raw(x, xlim = NULL, ci = TRUE)
.fdid_plot_data_dynamic(x, xlim = NULL, ci = TRUE)
.fdid_plot_data_overlap(x, xlim = NULL, bins = 30L, density = TRUE)
.fdid_plot_data_dml_scalar(x, alpha = NULL)
.fdid_plot_data_curve(
  x,
  curve = c("auto", "level", "derivative", "both"),
  interval = c("auto", "pointwise", "uniform", "both", "none"),
  xlim = NULL,
  ci = TRUE,
  show.uniform.CI = TRUE
)
.fdid_plot_data_support(x, xlim = NULL, Xdistr = "histogram")
.fdid_plot_data_dml_diagnostics(x, diagnostic, xlim = NULL)
.fdid_plot_data_list(x, vertical = TRUE, show_vartype = TRUE)
```

Notes:

- move or wrap the current `.fdid_curve_plot_data()` into this file;
- preserve the current curve helper's interval resolution and metadata logic;
- extend curve data with DML-specific method labels where available;
- keep scalar DML data separate from curve data so scalar methods are not
  presented as functions of `G`;
- keep named columns stable;
- make all helpers return ordinary data frames or small named lists of data
  frames;
- use these helpers in tests before changing renderers.

### New File: `R/plot-style.R`

Purpose: centralize defaults and normalize style arguments.

Proposed helpers:

```r
.fdid_plot_style(...)
.fdid_plot_palette(n, palette = "Set2")
.fdid_plot_labels(x, type, ...)
```

Style fields should cover:

- estimate line color and size;
- point color and size;
- pointwise interval color and alpha;
- simultaneous band color, linetype, and alpha;
- support histogram/density color and alpha;
- reference line color and linetype;
- treatment-period shade color and alpha;
- grid display;
- base theme choice for ggplot output.

### New File: `R/plot-ggplot.R`

Purpose: ggplot renderers.

Proposed helpers:

```r
fdid_autoplot(object, ...)
.fdid_autoplot_fdid(object, ...)
.fdid_autoplot_fdid_list(object, ...)
.fdid_ggplot_raw(data, style, labels, ...)
.fdid_ggplot_dynamic(data, style, labels, ...)
.fdid_ggplot_scalar(data, style, labels, ...)
.fdid_ggplot_overlap(data, style, labels, ...)
.fdid_ggplot_curve(curve_data, support_data, style, labels, ...)
.fdid_ggplot_dml_diagnostics(data, style, labels, ...)
.fdid_ggplot_fdid_list(data, style, labels, ...)
```

Keep the renderer small by assuming data has already been validated and
canonicalized.

Implementation choice to settle before coding:

- either make `fdid_autoplot()` a single exported function with internal
  dispatch on `inherits(object, "fdid")` and `inherits(object, "fdid_list")`;
- or make `fdid_autoplot()` an S3 generic through
  `UseMethod("fdid_autoplot")`, with methods for `fdid` and `fdid_list`.

The S3-generic version is cleaner if the package expects the helper to grow.
The single-function version is simpler while `ggplot2` remains only suggested.
For `fdid_list`, the default display should be a scalar comparison plot; if a
`type` argument is accepted, use `type = "comparison"` rather than reusing
`"raw"`, `"dynamic"`, or `"curve"`.

### Existing File: `R/plot.R`

Purpose after refactor:

- keep base `plot.fdid()` entry point;
- delegate data extraction to `R/plot-data.R`;
- keep base renderers or move them to `R/plot-base.R` if the file becomes too
  large;
- preserve current default return behavior.

### Existing File: `R/plot.fdid_list.R`

Purpose after refactor:

- delegate scalar extraction to `.fdid_plot_data_list()`;
- preserve base graphics behavior by default;
- allow the ggplot path to share the same data.

## Plot-Type Implementation Plan

### Raw Means Plot

Goal:

- make descriptive group means easier to customize while preserving their
  descriptive interpretation.

Data builder:

- input: `x$raw_means`;
- required columns: `time`, `group`, `meanY`;
- optional columns: `CI_Lower`, `CI_Upper`;
- output columns: `time`, `group`, `mean`, `conf.low`, `conf.high`,
  `has_interval`.

ggplot renderer:

- `geom_line()` when `connected = TRUE`;
- `geom_point()` always or when `connected = FALSE`;
- `geom_errorbar()` or `geom_ribbon()` for intervals;
- treatment period shading via `annotate("rect", ...)`;
- stable group palette.

Tests:

- accepts binary original-method raw means;
- accepts continuous-G median-split raw means;
- handles missing interval columns;
- respects `xlim`;
- does not relabel descriptive median-split means as adjusted estimates.

### Dynamic Plot

Goal:

- make event-time scalar estimates visually consistent across all methods.

Data builder:

- input: `x$dynamic`;
- row names or explicit time index become `time`;
- output columns: `time`, `estimate`, `std.error`, `conf.low`, `conf.high`,
  `period_role`.

ggplot renderer:

- zero line through `geom_hline()`;
- point estimates through `geom_point()`;
- optional line through `geom_line()`;
- confidence intervals through `geom_errorbar()` or `geom_linerange()`;
- treatment-period shading through rectangles;
- optional pre/event/post visual grouping.

Tests:

- handles reference-period zero row;
- respects `xlim`;
- handles missing confidence intervals;
- works for original, kernel, and DML methods.

### DML Scalar Plot

Goal:

- give scalar DML methods a single-object object-returning display without
  pretending they are continuous-G curves.

Applies to:

- `dml_binary`;
- `dml_plr`;
- `dml_incremental`;
- scalar summaries from `dml_flex` only when explicitly requested and clearly
  labeled as summaries, not as the primary curve estimand.

Data builder:

- input: `x$est$event`, `x$scalar_event`, `x$incremental_event`,
  `x$dml_metadata`, and `x$dml_diagnostics`;
- output columns: `method`, `target_estimand`, `target_scale`, `estimate`,
  `std.error`, `conf.low`, `conf.high`, `inference_method`,
  `scalar_inference`, `scalar_multiplier_crit`, `n_multiplier`, `n_splits`,
  `n_clusters`.

ggplot renderer:

- `geom_pointrange()` for the scalar estimate and interval;
- zero reference line when scientifically meaningful;
- subtitle or caption from `dml_metadata$target_estimand` and
  `scalar_event$scalar_inference`;
- optional small annotation for score-multiplier critical value when used.

Tests:

- `dml_plr` scalar plot data reports partially linear slope;
- `dml_incremental` scalar plot data reports observed-population average
  derivative;
- `dml_binary` scalar plot data reports binary AIPW/IRM contrast;
- score-multiplier fits are labeled score multiplier, not bootstrap;
- scalar DML methods still reject `type = "curve"`.

### Binary Overlap Plot

Goal:

- modernize propensity score diagnostics for binary-G methods.

Data builder:

- input: `x$ps`, `x$G`;
- validate binary `G`;
- compute histogram bins or densities by group;
- output named data frames for histogram and optional density.

ggplot renderer:

- overlay or mirrored histogram, matching current interpretability;
- optional density overlay;
- clear x-axis label as propensity score;
- no continuous-G support semantics.

Tests:

- errors clearly when `ps` is unavailable;
- errors clearly for continuous `G`;
- supports `xlim`;
- uses consistent bin breaks across groups.

### Continuous-G Curve Plot

Goal:

- make kernel and `dml_flex` curve plots publication-ready and object-returning
  while preserving current interval semantics.

Data builder:

- reuse current `.fdid_curve_plot_data()` logic;
- output one table per requested curve or a single long table with `curve`;
- output support data from observed `G` restricted to the displayed window;
- preserve fields for pointwise method, band method, zeta, empirical coverage,
  and effective bootstrap/multiplier draws.

ggplot renderer:

- estimate curve: `geom_line()`;
- pointwise interval: `geom_ribbon()`;
- simultaneous band: dashed `geom_line()` boundaries or distinct ribbon
  styling;
- support display: rug, histogram strip, or density strip below estimate
  scale;
- `curve = "both"`: use facets by default to avoid new layout dependencies in
  the first pass;
- optional marked values through `geom_vline()`;
- optional evaluation-grid ticks through `geom_rug()`;
- titles and subtitles should identify `kernel` versus `dml_flex` and the
  displayed interval method.

Interflex parity target:

- pointwise and uniform uncertainty are separate layers;
- support displays are x-window filtered before drawing;
- user can hide uniform bands through `show.uniform.CI = FALSE`;
- marked moderator values are preserved.

Tests:

- kernel robust plot has pointwise intervals but no uniform band unless a true
  band is stored;
- kernel bootstrap plot uses stored percentile pointwise intervals and
  quantile-envelope bands;
- `dml_flex` plot uses stored pointwise and practical multiplier band fields;
- `dml_flex` with `signal_map = "blp_spline"` labels finite-grid Gaussian
  max-t BLP-spline covariance bands separately from practical bands;
- `curve = "level"`, `"derivative"`, and `"both"` work;
- `Xdistr = "histogram"`, `"density"`, `"rug"`, and `"none"` work;
- missing values split line/ribbon segments instead of connecting across gaps;
- ggplot layer data can be inspected with `ggplot_build()` when available.

### `fdid_list` Scalar Comparison Plot

Goal:

- modernize multi-estimator scalar comparison while preserving the current
  vertical/horizontal base plot.

Data builder:

- input: list of `fdid` objects;
- source: `obj$est$event`;
- output columns: `label`, `method`, `vartype`, `estimate`, `conf.low`,
  `conf.high`, `index`.

ggplot renderer:

- `geom_pointrange()` for vertical plots;
- `coord_flip()` for horizontal plots;
- zero reference line;
- optional ordering controls later.

Tests:

- labels match current `show_vartype` behavior;
- vertical and horizontal outputs preserve estimate-axis semantics;
- errors clearly when an object lacks event estimates.

### Diagnostics Plotting

Goal:

- add diagnostics after the main plotting layer is stable.

Candidates:

- continuous-G observed support plot;
- DML conditional-density diagnostics;
- DML signal and correction term tail diagnostics;
- DML score or influence diagnostics for scalar methods when stored;
- `dml_binary` propensity-score overlap diagnostics;
- `dml_incremental` density-derivative diagnostics, clearly labeled as
  diagnostics for an observed-population scalar average derivative;
- learner/split/fold failure summaries;
- density-method sensitivity plots;
- curve sensitivity plots across `eval_g`, `trim`, `h0`, learner, or
  `density_method`.

Recommended timing:

- do not include diagnostics in the first ggplot implementation;
- add them after core raw/dynamic/overlap/curve/list plots share a stable
  data layer.

## Dependency Plan

### Recommended First Step

Add `ggplot2` to `Suggests`, not `Imports`, and implement
`fdid_autoplot()` with runtime checks:

```r
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Install ggplot2 to use fdid_autoplot().")
}
```

Advantages:

- avoids adding a hard dependency to existing users;
- keeps base plotting available everywhere;
- allows tests to use `testthat::skip_if_not_installed("ggplot2")`.

Tradeoff:

- true S3 `autoplot.fdid()` registration is less clean when `ggplot2` is only
  suggested.

### Later Option

Move `ggplot2` to `Imports` if:

- `fdid_autoplot()` becomes a central public API;
- documentation and tutorials rely on ggplot objects;
- maintainers prefer standard `autoplot()` registration.

### Avoid Initially

Avoid adding `gridExtra`, `patchwork`, `ggplotify`, or `scales` in the first
pass unless a specific plot cannot be implemented cleanly with ggplot2 alone.

For `curve = "both"`, use facets first. Consider multi-panel layout packages
only after the basic ggplot path is stable.

## Migration Phases

### Phase 0: Freeze Current Base Behavior

Status: planning baseline.

Actions:

1. Record current `plot.fdid()` and `plot.fdid_list()` behavior.
2. Keep existing smoke tests for base graphics.
3. Avoid changing default return semantics.

Acceptance:

- all existing plot calls still draw without error;
- base graphics path still returns `invisible(NULL)` for `plot.fdid()`;
- `plot.fdid_list()` still returns `invisible(x)`.

### Phase 1: Add Plot-Data Builders

Actions:

1. Create `R/plot-data.R`.
2. Move or wrap `.fdid_curve_plot_data()` into the shared data layer.
3. Add raw, dynamic, overlap, support, and fdid-list data builders.
4. Add structural tests for every data builder.

Acceptance:

- every plot type can be represented as canonical data;
- data helpers are independent of graphics devices;
- no public plotting behavior changes.

### Phase 2: Refactor Base Plots To Use Plot Data

Actions:

1. Make `plot.fdid(type = "raw")` consume raw plot data.
2. Make `plot.fdid(type = "dynamic")` consume dynamic plot data.
3. Make `plot.fdid(type = "overlap")` consume overlap plot data.
4. Keep the current curve base path aligned with shared curve data.
5. Make `plot.fdid_list()` consume list plot data.

Acceptance:

- base visual output remains materially unchanged;
- tests prove that rendering reads the same canonical data that ggplot will
  later consume;
- no estimator output changes.

### Phase 3: Add `fdid_autoplot()` For Core Plots

Actions:

1. Add `ggplot2` to `Suggests`.
2. Add exported `fdid_autoplot()`.
3. Implement ggplot renderers for `type = "raw"`, `"dynamic"`, and
   `fdid_list`.
4. Add `type = "scalar"` for single-object scalar DML estimates.
5. Add object-class and layer-data tests.

Acceptance:

- `fdid_autoplot(fit, type = "raw")` returns a `ggplot` object;
- `fdid_autoplot(fit, type = "dynamic")` returns a `ggplot` object;
- `fdid_autoplot(fit, type = "scalar")` returns a `ggplot` object for scalar
  DML methods;
- `fdid_autoplot(fdid_list(...))` returns a `ggplot` object;
- scalar score-multiplier intervals are not labeled as bootstrap intervals;
- tests skip cleanly when `ggplot2` is not installed.

### Phase 4: Add ggplot Continuous-G Curve Path

Actions:

1. Implement `.fdid_ggplot_curve()`.
2. Support level, derivative, and both-curve displays.
3. Support pointwise ribbons and uniform-band boundaries.
4. Support support displays: histogram, density, rug, none.
5. Support `show.uniform.CI`, `show.eval_g`, `diff.values`, and
   `contrast.values`.
6. Add metadata-aware labels for interval type.

Acceptance:

- kernel and `dml_flex` curve plots return ggplot objects;
- interval layers match stored inference fields;
- support displays respect `xlim`;
- no plot labels bootstrap quantile-envelope bands as studentized sup-t bands;
- no plot labels practical `dml_flex` bands as theorem-level bands.

### Phase 5: Add ggplot Overlap And Diagnostics

Actions:

1. Implement `.fdid_ggplot_overlap()`.
2. Add continuous-G support-only plotting if useful.
3. Add DML diagnostic plot prototypes after the core API stabilizes.
4. Keep DML diagnostics under an explicit opt-in `type = "diagnostics"` or a
   named diagnostic helper.

Acceptance:

- binary propensity overlap has a ggplot path;
- continuous-G support diagnostics are clearly separated from binary
  propensity overlap;
- DML diagnostics are opt-in and do not clutter ordinary `plot.fdid()`.
- `dml_incremental` diagnostics do not imply a fixed-`g` derivative curve.

### Phase 6: Decide `autoplot` S3 And `plot(engine = "ggplot")`

Actions:

1. Decide `ggplot2` dependency status.
2. If `ggplot2` is imported, register `autoplot.fdid()` and
   `autoplot.fdid_list()`.
3. Optionally add `engine = c("base", "ggplot")` to `plot.fdid()` and
   `plot.fdid_list()`.

Acceptance:

- dependency policy is documented;
- base default remains unchanged;
- ggplot path return semantics are explicit.

### Phase 7: Documentation And Tutorial Integration

Actions:

1. Update roxygen/Rd docs for new plotting helpers and arguments.
2. Update `CONTINUOUS_G_ESTIMATOR_REFERENCE.md` for object-returning plots.
3. Update `CONTINUOUS_G_PROJECT_LOG.md` with implementation details.
4. Update `TODO.md` by closing or reframing the current optional ggplot TODO.
5. Update FDID Notes only where plotting semantics or displayed uncertainty
   wording changes.
6. Update the package vignette decision explicitly: either keep its base
   `plot()` examples as the stable default path or add short
   `fdid_autoplot()` examples after the API stabilizes.
7. Leave tutorial rewrite to the separate full tutorial reorganization unless
   the user explicitly requests an interim tutorial patch.
8. During the tutorial rewrite, replace or explicitly preserve local helpers
   such as `plot_curve_overlay()` after checking whether the new plotting
   layer covers their use case.
9. Leave `NEWS.md` and `HANDOFF.md` out of scope unless a release or handoff is
   requested.

Acceptance:

- reference docs identify base versus ggplot paths;
- docs state what each interval layer means;
- docs state that curve plots display stored `eval_g` estimates;
- package vignette strategy is explicit;
- tutorial-only plotting helpers have a replacement or preservation decision;
- tutorial changes happen in the dedicated tutorial rewrite.

## Test Plan

### Structural Data Tests

Add tests for:

- raw plot data columns, row counts, group labels, interval availability;
- dynamic plot data time index, reference row, interval availability;
- DML scalar plot data target labels, inference labels, and multiplier fields;
- overlap plot data group validation and shared breaks;
- curve plot data interval resolution and metadata;
- DML curve plot data signal-map, density-method, and band-method labels;
- DML diagnostic plot data long-format labels and graceful missing-data
  behavior;
- support plot data x-window filtering;
- fdid-list plot data labels and confidence intervals.

### ggplot Tests

Use `testthat::skip_if_not_installed("ggplot2")`.

Add tests for:

- returned object inherits from `ggplot`;
- expected layer count for each plot type;
- expected mapped data columns appear in `ggplot_build()`;
- scalar DML object plots render as pointrange-style layers;
- pointwise and uniform curve layers can be distinguished;
- support layers disappear when `Xdistr = "none"`;
- marked values appear when requested.

### Backward Compatibility Tests

Keep or add smoke tests for:

- `plot(fit, type = "raw")`;
- `plot(fit, type = "dynamic")`;
- `plot(fit, type = "overlap")`;
- `plot(fit, type = "curve", curve = "level")`;
- `plot(fit, type = "curve", curve = "derivative")`;
- `plot(fit, type = "curve", curve = "both")`;
- `plot(fdid_list(...))` vertical and horizontal.

### Inference Display Tests

Add tests that verify:

- kernel robust pointwise intervals are not exposed as uniform bands;
- kernel bootstrap pointwise percentile intervals and quantile-envelope bands
  are labeled separately;
- `show.uniform.CI = FALSE` suppresses uniform-band data/layers;
- `interval = "pointwise"`, `"uniform"`, `"both"`, and `"none"` have the
  expected layer behavior;
- `dml_flex` non-BLP band labels remain practical signal-mapping labels;
- `dml_flex` BLP-spline band labels identify Gaussian max-t BLP covariance;
- scalar DML score multiplier intervals are not labeled as bootstrap;
- `dml_incremental` scalar plots identify the observed-population average
  derivative target.

### Device Hygiene Tests

Because base graphics can leave device state behind, keep tests for:

- graphics parameters restored after multi-panel curve plots;
- no unwanted `Rplots.pdf` dependency in package tests;
- no snapshots that are brittle across graphics devices unless deliberately
  introduced later.

## Documentation Plan

### Roxygen/Rd

Update:

- `plot.fdid`;
- `plot.fdid_list`;
- new `fdid_autoplot`;
- future `autoplot.fdid` and `autoplot.fdid_list` if registered;
- any public `fdid_plot_data` helper if exported.

Must state:

- base graphics path return semantics;
- ggplot path return semantics;
- optional `ggplot2` requirement;
- interval-layer meaning;
- continuous-G curve grid semantics.

### Markdown Reference

Update `CONTINUOUS_G_ESTIMATOR_REFERENCE.md`:

- add a plotting section that separates base plot, object-returning ggplot
  plot, and plot-data semantics;
- state that ggplot output uses the same stored estimates and interval fields
  as base plots;
- keep the current inference wording for kernel quantile-envelope bands and
  `dml_flex` practical bands.
- add a DML plotting subsection that separates scalar DML, `dml_binary`
  overlap, `dml_flex` curves, and opt-in diagnostics.

### Project Log

Update `CONTINUOUS_G_PROJECT_LOG.md` when implementation begins:

- record files changed;
- record new public API;
- record dependency policy;
- record tests run;
- record any deliberate departures from the plan.

### TODO

Update `TODO.md`:

- replace the current optional ggplot TODO with phase-specific remaining work;
- keep diagnostics plotting as a later phase until implemented;
- keep scalar DML plotting and DML diagnostics as explicit plotting-layer
  work items;
- keep tutorial rewrite as a separate item.

### FDID Notes

Update `FDID Notes/main.tex` only if the implementation changes the way plots
are described in the paper notes:

- object-returning plotting is package documentation, not theory;
- interval semantics in Notes should change only if displayed uncertainty
  semantics change;
- no new theory claim should be introduced by a plotting-layer update.

### Tutorials

Do not patch tutorials in the first plotting-layer implementation unless the
user requests it. The tutorial needs a broader rewrite, and plotting examples
should be incorporated there after the API stabilizes.

Tutorial rewrite checklist:

- replace ordinary `plot.fdid()` examples only where object-returning plots add
  real pedagogical value;
- keep base `plot()` examples where backward-compatible API teaching is the
  goal;
- decide whether `plot_curve_overlay()` should become unnecessary once
  curve-comparison or truth-overlay support exists in the new layer;
- keep synthetic truth overlays local to tutorials unless they become a
  supported package feature.

### Package Vignette

Update `fdid-main/vignettes/fdid.Rmd` only after the new API stabilizes:

- keep the current base `plot(result, type = "raw")` and
  `plot(result, type = "dynamic")` examples if the vignette should emphasize
  the default dependency-light path;
- add a compact `fdid_autoplot()` example if maintainers want the vignette to
  advertise object-returning plotting.

### NEWS And HANDOFF

Do not update `NEWS.md` or `HANDOFF.md` for this planning artifact. Update
them only for a release, shipment, or explicit handoff request.

## Risk Register

### Dependency Risk

Risk:

- adding `ggplot2` to `Imports` changes the package dependency footprint.

Mitigation:

- start with `Suggests` and `fdid_autoplot()`;
- move to `Imports` only after maintainers accept the API.

### API Risk

Risk:

- exporting plot-data helpers too early locks the package into unstable column
  schemas.

Mitigation:

- keep plot-data helpers internal first;
- export only after at least one implementation cycle.

### Return-Semantics Risk

Risk:

- changing `plot.fdid()` to return visible objects may surprise existing users
  or tests.

Mitigation:

- keep base `plot()` default behavior unchanged;
- use `fdid_autoplot()` for visible ggplot objects.

### Inference-Label Risk

Risk:

- a nicer plot could make practical or pointwise intervals look stronger than
  they are.

Mitigation:

- build interval labels from stored metadata;
- test layer behavior for pointwise versus uniform intervals;
- keep kernel quantile-envelope and `dml_flex` practical-band caveats in docs.

### Layout Risk

Risk:

- support strips and multi-curve panels can become hard to align in pure
  ggplot2.

Mitigation:

- use facets for `curve = "both"` in the first ggplot pass;
- avoid new layout dependencies until a clear need appears.

### Maintenance Risk

Risk:

- supporting base graphics and ggplot can duplicate code.

Mitigation:

- make canonical plot data the shared layer;
- keep renderer helpers thin.

## Open Decisions

1. Should `ggplot2` remain in `Suggests` with `fdid_autoplot()`, or move to
   `Imports` to support clean S3 `autoplot` registration?
2. Should `fdid_plot_data()` become an exported helper, or stay internal?
3. Should `plot(..., engine = "ggplot")` be added, or should object-returning
   plots remain separate through `fdid_autoplot()`?
4. Should curve `support.panel = "separate"` be replicated in ggplot now, or
   should the first ggplot path use embedded support/facets only?
5. Should single-object scalar DML plots use `type = "scalar"` or
   `type = "estimate"`?
6. Should diagnostic plots be part of `fdid_autoplot(type = "diagnostics")` or
   separate named helpers such as `plot_dml_diagnostics()`?
7. Should DML scalar influence values be plotted only as diagnostics, or should
   scalar displays include optional influence insets?
8. Should file-saving be built into the plotting API, following `interflex`, or
   should users rely on `ggplot2::ggsave()`?

## Recommended First Implementation Slice

This package-wide plan is not the next code pass. The next code pass should be
the DML-specific Phase A in `DML_PLOT_UPDATE_PLAN_2026-05-30.md`:

1. align `dml_flex` base curve labels and controls with kernel;
2. keep scalar DML methods from being plotted as curves;
3. update focused docs and tests for DML plot semantics;
4. avoid adding `ggplot2`, `fdid_autoplot()`, or package-wide raw/dynamic/list
   refactors in that pass.

If the package-wide plan is activated later, its first code pass should be
deliberately small:

1. Add internal plot-data helpers for raw, dynamic, fdid-list, and curve.
2. Add internal scalar DML plot-data helpers.
3. Refactor base raw/dynamic/fdid-list plots to use the helpers.
4. Add `ggplot2` to `Suggests`.
5. Add exported `fdid_autoplot()` for raw, dynamic, scalar DML, and fdid-list
   plots.
6. Add structural tests and ggplot object tests.
7. Update roxygen/Rd, `CONTINUOUS_G_PROJECT_LOG.md`, `TODO.md`, and the
   estimator reference.

Defer to the second code pass:

- continuous-G ggplot curve path;
- overlap ggplot path;
- diagnostics plotting;
- `autoplot` S3 registration;
- `plot(engine = "ggplot")`;
- tutorial rewrite.

Because current user interest is now concentrated on DML plot alignment, the
package-wide ggplot/autoplot plan should remain an explicit later decision
rather than becoming an automatic next step.

## Acceptance Criteria For The Full Plotting-Layer Update

The package-wide update is complete when:

- every existing plot type has a canonical plot-data builder;
- scalar DML methods have a canonical single-object scalar plot-data builder;
- every existing base plot still works with the old default call style;
- users can obtain ggplot objects for raw, dynamic, overlap, curve, and
  fdid-list plots;
- users can obtain ggplot objects for scalar DML estimates;
- continuous-G ggplot curve plots support level, derivative, both, pointwise
  intervals, stored simultaneous bands, support displays, and marked G values;
- DML diagnostic plots are available through an explicit opt-in surface;
- interval labels and docs match the implemented inference methods;
- DML plots never label score-multiplier intervals as bootstrap and never label
  practical `dml_flex` bands as theorem-level uniform DML bands;
- tests cover data builders, ggplot object returns, interval-layer behavior,
  support filtering, and base-graphics smoke behavior;
- package docs explain base versus ggplot return semantics;
- `CONTINUOUS_G_PROJECT_LOG.md`, `CONTINUOUS_G_ESTIMATOR_REFERENCE.md`, and
  `TODO.md` are updated;
- FDID Notes are updated only if plot semantics or displayed uncertainty
  language changes;
- tutorials are updated in the separate tutorial rewrite.
