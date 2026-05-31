# Kernel Plot Update Plan

Date: 2026-05-30

This plan records a proposed update to the `fdid` continuous-`G` kernel curve
plot. It is motivated by the current `fdid-main` plotting code, the kernel
inference update, and the plotting conventions used by `interflex`.

This is a plotting and documentation plan only. Tutorial rewrite, `NEWS.md`,
and `HANDOFF.md` are intentionally out of scope for this pass.

Implementation status: implemented in the base-graphics path on 2026-05-30.
The optional `ggplot2`/`autoplot.fdid()` path remains a future extension.

## Objective

The current kernel estimator now stores the information needed for richer curve
plots: level and derivative point estimates, analytical pointwise intervals,
bootstrap pointwise intervals, bootstrap quantile-envelope simultaneous bands,
and observed `G` support. The plot update should make these objects visible in
a way that is clear to users and consistent with the inference language in the
FDID Notes.

The target is practical parity with `interflex` plotting for kernel moderation
curves:

- show the estimate as a curve over the evaluation grid;
- show pointwise confidence intervals as a ribbon;
- optionally show simultaneous bands as dashed boundary lines;
- show the observed moderator or treatment support below the curve;
- expose enough style and display controls for publication-quality plots;
- clearly label the implemented band as a finite-grid bootstrap
  quantile-envelope band, not as a fully studentized sup-t band.

## Current FDID Plot State

`plot.fdid(type = "curve")` currently:

- supports kernel and `dml_flex` curve objects;
- supports `curve = "auto"`, `"level"`, or `"derivative"`;
- uses stored `eval_g` as the plotting grid and respects `xlim`;
- draws a line for the point estimate;
- draws a pointwise confidence ribbon when interval fields are available;
- draws stored band boundaries when `*_band_lower` and `*_band_upper` exist;
- supports `Xdistr = "histogram"`, `"hist"`, `"density"`, or `"none"` for the
  observed `G` support display;
- uses base R graphics and returns `invisible(NULL)`.

Current limitations:

- interval display is controlled implicitly by `ci` and stored fields, rather
  than by explicit pointwise versus uniform display choices;
- pointwise intervals and simultaneous bands are visually present but not
  always clearly distinguished in the legend or plot annotation;
- the plot does not expose a `show.uniform.CI` style option analogous to
  `interflex`;
- the plot does not support `curve = "both"` for side-by-side level and
  derivative views;
- support displays are embedded in the same plotting region, which can make the
  y-axis harder to read in some empirical examples;
- there is no optional rug display for observed `G` values or evaluation-grid
  ticks;
- missing values or trimmed grid gaps are not handled through a separate
  plotting-data layer;
- the function cannot return a reusable `ggplot` object;
- style controls are limited compared with `interflex`;
- there is no direct plot annotation for contrast values such as `g0` and `g1`.

## Reference: interflex

The relevant `interflex` plot behavior is:

- `plot.interflex()` builds a canonical plotting data frame with the moderator
  grid, point estimates, standard errors, pointwise interval fields, and
  optional uniform interval fields;
- pointwise intervals are rendered as a semi-transparent `geom_ribbon`;
- uniform confidence bands are rendered as dashed gray boundary lines when
  `show.uniform.CI = TRUE`;
- `Xdistr` can display histogram, density, or no support information;
- histogram and density inputs are filtered to `xlim` before drawing;
- support displays are scaled into a lower strip of the plot;
- style options expose line size, line color, confidence-interval color,
  confidence-interval opacity, histogram color, density color, and grid/theme
  choices;
- `diff.values` can mark moderator values for contrast-style displays;
- the function returns a `ggplot` object, making downstream customization easy.

FDID should borrow these plotting ideas, while preserving FDID-specific
inference labels and keeping the base-graphics method backward compatible.

## Proposed User-Facing API

Keep the current API working:

```r
plot(fit, type = "curve", curve = "derivative")
plot(fit, type = "curve", curve = "level")
plot(fit, type = "curve", Xdistr = "density")
```

Add explicit interval controls:

```r
interval = c("auto", "pointwise", "uniform", "both", "none")
show.uniform.CI = TRUE
```

Suggested semantics:

- `interval = "auto"`: show pointwise intervals when available and show uniform
  bands when the object stores them;
- `interval = "pointwise"`: show only pointwise intervals;
- `interval = "uniform"`: show only simultaneous bands;
- `interval = "both"`: show pointwise intervals and simultaneous bands;
- `interval = "none"`: show only the point estimate;
- `show.uniform.CI` is retained as an interflex-style alias for users familiar
  with that API.

Add curve selection:

```r
curve = c("auto", "level", "derivative", "both")
```

Suggested semantics:

- `curve = "level"`: plot `mu_hat(g)` or `theta_hat(g)`;
- `curve = "derivative"`: plot `delta_hat(g)` or `theta_prime_hat(g)`;
- `curve = "both"`: show level and derivative views in two panels;
- `curve = "auto"`: preserve current method-specific defaults.

Extend support display choices:

```r
Xdistr = c("histogram", "hist", "density", "rug", "none")
support.panel = c("embedded", "separate")
show.eval_g = FALSE
```

Suggested semantics:

- `rug`: draw observed `G` values as ticks along the bottom of the curve;
- `support.panel = "embedded"`: preserve the current lower-strip behavior;
- `support.panel = "separate"`: draw support in a small aligned lower panel;
- `show.eval_g = TRUE`: draw a second set of ticks for the stored evaluation
  grid so users can see where the estimator was evaluated.

Add style options modeled on `interflex`:

```r
line.color
line.size
ci.color
ci.alpha
band.color
band.lty
hist.color
density.color
support.alpha
theme.bw
show.grid
```

Add contrast annotation:

```r
diff.values = NULL
contrast.values = NULL
```

Suggested semantics:

- accept one or more `G` values and draw vertical markers;
- if exactly two values are supplied and `fdid_contrast()` can evaluate them,
  optionally annotate the estimated level contrast in the margin or legend;
- keep this annotation optional so the default plot stays uncluttered.

Add an optional ggplot path:

```r
engine = c("base", "ggplot")
autoplot.fdid()
```

Suggested semantics:

- keep `plot.fdid()` base-graphics compatible by default;
- implement `autoplot.fdid()` or `plot(..., engine = "ggplot")` as an optional
  `ggplot2` path when `ggplot2` is installed;
- return the `ggplot` object from the ggplot path for downstream editing.

## Proposed Internal Data Layer

Add a helper that converts fitted curve objects into a plotting data frame:

```r
.fdid_curve_plot_data <- function(x, curve, interval, xlim = NULL) {
  # returns one row per displayed eval_g point
}
```

The returned data frame should include:

```text
g
estimate
pointwise_lower
pointwise_upper
band_lower
band_upper
curve
method
vartype
pointwise_method
band_method
band_zeta
band_coverage
n_boot_eff
```

Advantages:

- one place resolves whether pointwise intervals come from analytical normal
  intervals or bootstrap percentile intervals;
- one place resolves whether uniform bands come from kernel
  `bootstrap_quantile_envelope`, DML multiplier bands, or older pointwise
  normal fallback fields;
- base graphics and future ggplot rendering can share the same data;
- `xlim` filtering and missing-value handling become testable;
- line and ribbon drawing can avoid connecting across invalid or missing grid
  segments.

Add a small segment helper:

```r
.fdid_contiguous_segments <- function(data, cols)
```

This helper should split plotting data into valid contiguous segments before
drawing ribbons or band lines. This prevents a polygon from connecting across
`NA` intervals or intentionally trimmed gaps.

## Rendering Plan

### Phase 1: Base-Graphics Upgrade

Update `R/plot.R` without adding a hard plotting dependency:

- add `interval`, `show.uniform.CI`, and `curve = "both"` support;
- route curve plotting through `.fdid_curve_plot_data()`;
- draw pointwise intervals as the primary translucent ribbon;
- draw simultaneous bands as dashed boundary lines;
- add a compact legend that distinguishes:
  - point estimate;
  - pointwise interval;
  - simultaneous band;
- add optional annotation text for the simultaneous-band method:
  `bootstrap quantile-envelope`;
- add `Xdistr = "rug"` and `show.eval_g`;
- make `support.panel = "separate"` available through a base `layout()` path;
- expose style arguments with conservative defaults matching current FDID
  colors.

### Phase 2: Optional ggplot/autoplot Path

Add a ggplot rendering path in a later pass:

- add `ggplot2` to `Suggests`, not `Imports`;
- implement `autoplot.fdid()` or `plot.fdid(..., engine = "ggplot")`;
- mimic the `interflex` visual grammar:
  - curve line;
  - pointwise ribbon;
  - dashed uniform bands;
  - lower support histogram/density/rug layer;
  - optional `diff.values` markers;
- return the `ggplot` object invisibly or visibly according to normal R plotting
  conventions;
- fall back to base graphics with a clear message if `ggplot2` is unavailable.

## Inference Display Rules

For kernel robust fits:

- pointwise interval: analytical local-linear sandwich interval;
- uniform band: unavailable unless future analytical or multiplier bands are
  added;
- do not label robust pointwise intervals as simultaneous bands.

For kernel bootstrap fits:

- pointwise interval: bootstrap percentile interval when stored;
- pointwise SE: bootstrap standard deviation when stored;
- simultaneous band: finite-grid bootstrap quantile-envelope band;
- display method label: `bootstrap quantile-envelope`;
- optional metadata: `zeta`, realized coverage, effective bootstrap draws, and
  row versus cluster bootstrap.

For `dml_flex` fits:

- pointwise interval: currently stored curve pointwise interval;
- simultaneous band: practical multiplier signal-mapping band when stored;
- display method label should avoid theorem-level language unless the
  underlying inference is upgraded.

Important language rule:

- Use `finite-grid bootstrap quantile-envelope band` for the current kernel
  bootstrap band.
- Use `studentized sup-t band` only if the code implements a max-statistic of
  the form `sup_g |(theta_hat^*(g) - theta_hat(g)) / se_hat^*(g)|`.

## Documentation Plan

Update only the active reference documents:

- `R/plot.R` roxygen and the generated `man/plot.fdid.Rd`;
- `CONTINUOUS_G_ESTIMATOR_REFERENCE.md`;
- `CONTINUOUS_G_PROJECT_LOG.md`;
- `TODO.md`;
- FDID Notes, only if inference-display wording needs to distinguish the plot
  layer from the estimator layer.

Do not update in this pass:

- tutorial files;
- `NEWS.md`;
- `HANDOFF.md`.

## Tests

Add focused tests for:

- `plot(type = "curve", curve = "level")` with kernel robust output;
- `plot(type = "curve", curve = "derivative")` with kernel robust output;
- bootstrap kernel curve plot reads stored pointwise percentile intervals;
- bootstrap kernel curve plot reads stored quantile-envelope bands;
- `interval = "pointwise"`, `"uniform"`, `"both"`, and `"none"`;
- `show.uniform.CI = FALSE` suppresses simultaneous bands;
- `curve = "both"` draws without error;
- `Xdistr = "histogram"`, `"density"`, `"rug"`, and `"none"`;
- `xlim` filters curve rows and support rows consistently;
- plotting-data helper does not connect across missing interval values;
- ggplot path returns a `ggplot` object when implemented and `ggplot2` is
  installed.

Because plot tests are brittle across devices, prefer structural tests for the
plot-data helper and smoke tests that ensure plotting calls complete without
error.

## Acceptance Criteria

The update is complete when:

- users can separately request pointwise intervals, simultaneous bands, both, or
  neither;
- kernel bootstrap plots clearly distinguish bootstrap percentile pointwise
  intervals from bootstrap quantile-envelope simultaneous bands;
- robust kernel plots do not present pointwise normal intervals as uniform
  inference;
- level and derivative plots are both supported, including a two-panel option;
- observed `G` support can be shown as histogram, density, rug, or omitted;
- support displays respect `xlim` and do not obscure the main curve;
- plot documentation explicitly states what each displayed interval means;
- tests cover interval selection, curve selection, support display, and xlim
  filtering.

## Implementation Order

1. Add `.fdid_curve_plot_data()` and unit-test its interval resolution.
2. Add `.fdid_contiguous_segments()` and unit-test missing-value behavior.
3. Update base curve plotting to consume the plotting data helper.
4. Add `interval`, `show.uniform.CI`, style options, and improved legends.
5. Add `curve = "both"`.
6. Add `Xdistr = "rug"` and optional evaluation-grid ticks.
7. Add `support.panel = "separate"` if base layout remains stable.
8. Update roxygen, Rd, reference md, project log, and TODO.
9. Run the continuous-G test subset and then the full package test suite.
10. Consider the optional ggplot/autoplot path as a second implementation pass.
