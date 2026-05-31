# Continuous-G Contrast Plot Refinement Plan

Date: 2026-05-31

This document plans a focused refinement of the already implemented
continuous-G fixed-reference contrast coefficient plot:

```r
plot(fit, type = "contrast")
```

The current contrast plot is a good first implementation: it delegates
estimates and inference to `fdid_contrast()`, works for `method = "kernel"` and
`method = "dml_flex"`, supports explicit reference and target values, and keeps
the plotting layer out of estimator logic. This refinement should keep those
properties while making the default display more interpretable and making target
selection easier.

Execution note, later on 2026-05-31: this refinement has been implemented in
the R plotting layer, tests, help, estimator reference, tutorial copy, TODO,
project log, and the package-wide plotting plan. The broader ggplot/autoplot
modernization remains out of scope.

## Pre-Refinement State

Before implementing this refinement, the package code behaved as follows:

- `ref.g = NULL` chooses the finite fitted `eval_g` value closest to the median
  of the fitted evaluation grid.
- `target.g = NULL` and `contrast.values = NULL` choose up to `n.contrast`
  target points from finite fitted `eval_g`, excluding `ref.g`.
- `n.contrast = 7L` and `max.contrast = 10L`.
- `target.g` and `contrast.values` are numeric level values, not quantile
  instructions.
- The reference level was not drawn unless the user explicitly included it in
  `target.g` or `contrast.values`.
- The plot data are built by calling `fdid_contrast()` once for each displayed
  non-reference target.
- `ci = FALSE` suppresses contrast confidence-interval bars.

## Check On `PACKAGE_PLOTTING_LAYER_UPDATE_PLAN_2026-05-31.md`

`PACKAGE_PLOTTING_LAYER_UPDATE_PLAN_2026-05-31.md` has already been updated to
the first implemented contrast-plot version. It records that
`plot.fdid(type = "contrast")` has been implemented for kernel and `dml_flex`,
that the fixed-reference contrast slice is complete, and that the larger
package-wide plot-data / `fdid_autoplot()` modernization remains deferred.

The package-wide plan has since been patched again so its contrast-plot section
reflects the implemented reference-row and quantile-target semantics.

## User-Facing Goals

### 1. Draw The Reference Row By Default

The contrast plot should display the reference point with coefficient zero:

```text
mu(g_ref) - mu(g_ref) = 0
theta(g_ref) - theta(g_ref) = 0
```

Recommended API:

```r
plot(fit, type = "contrast", include.ref = TRUE)
```

with `include.ref = TRUE` as the new default.

Implementation recommendation:

- Add a reference row to the contrast plot data whenever `include.ref = TRUE`.
- Mark that row with `is_reference = TRUE`.
- Set `estimate = 0`, `std.error = 0`, `conf.low = 0`, and `conf.high = 0`.
- Set `inference = "reference"` and a note such as
  `"Reference row; contrast of fitted curve value with itself."`
- Do not rely on `fdid_contrast(x, g0 = ref.g, g1 = ref.g)` for this row,
  because pointwise fallback logic intentionally ignores endpoint covariance
  and would otherwise give a nonzero approximate interval for a contrast that
  is algebraically zero.
- Keep non-reference rows delegated to `fdid_contrast()`.
- In the plot's top inference label, summarize inference for non-reference
  rows only, so the added `"reference"` row does not make every plot say that
  inference varies.

Visual recommendation:

- Draw the reference point at zero on the same x-axis as the other target
  values.
- Use the same point style or a lightly distinct point style; do not add a
  second legend unless needed.
- Keep the horizontal zero line.
- Order points by numeric `target_g` in automatic and quantile modes.

### 2. Add Quantile-Based Target Selection

Users should be able to request contrast targets by quantile points rather than
typing absolute `G` levels.

Important ambiguity:

- A bare numeric `target.g = 5` currently means the absolute level `G = 5`.
- Reinterpreting bare numeric `5` as "five quantile points" would break or
  confuse level-based usage.

Recommended API:

```r
plot(fit, type = "contrast", target.type = "level",
     target.g = c(0.2, 0.5, 0.8))

plot(fit, type = "contrast", target.type = "quantile",
     target.g = 5)
```

Semantics:

- `target.type = "level"` preserves the current meaning: `target.g` and
  `contrast.values` are absolute `G` levels.
- `target.type = "quantile"` changes the meaning of `target.g`:
  - scalar integer `K > 1`: use `K` quantile points, labeled `Q1`, ..., `QK`;
  - numeric vector in `[0, 1]`: use those quantile probabilities directly;
  - numeric vector in `[1, 100]` may optionally be accepted as percent values,
    but only if this can be documented without ambiguity.
- For the example `target.type = "quantile", target.g = 5`, the plot should
  draw five quantile points including the reference, and the first quantile
  point should be the default reference.

Recommended quantile source:

```r
quantile.source = c("support", "eval_g")
```

Default recommendation:

- `quantile.source = "support"` should use observed continuous `G` values
  stored in the fitted object, restricted to the fitted `eval_g` range.
- If stored `G` support is missing or has too few unique finite values, fall
  back to `eval_g`.
- `quantile.source = "eval_g"` should compute quantile points directly from
  finite fitted `eval_g`; this guarantees all selected values are in fitted
  support.

Validation:

- Quantile-selected values must be finite and inside the fitted `eval_g` range.
- Duplicated quantile values after rounding or ties should be uniqued with a
  clear warning if fewer rows than requested remain.
- If `interpolate = "none"`, quantile-selected off-grid values should either be
  snapped to the nearest fitted grid point with an explicit documented rule, or
  error clearly. The safer default is to keep `interpolate = "linear"` for
  quantile mode and document that choice.

Possible shorthand:

```r
plot(fit, type = "contrast", target.g = "q5")
plot(fit, type = "contrast", target.g = "5")
```

This shorthand is optional. If added, it should be parsed as
`target.type = "quantile", target.g = 5`. The implemented refinement accepts
both `"q5"` and bare character `"5"` as quantile-count shorthand, while leaving
bare numeric `target.g = 5` in level mode unless the user sets
`target.type = "quantile"`.

### 3. Upgrade Reference Selection

Users should be able to request common reference rules without manually
computing a number.

Recommended API:

```r
plot(fit, type = "contrast", ref.g = "median")
plot(fit, type = "contrast", ref.g = "min")
plot(fit, type = "contrast", ref.g = "mean")
plot(fit, type = "contrast", ref.g = 0.5)
```

Recommended semantics:

- `ref.g = NULL`: preserve current default for level mode, equivalent to
  `ref.g = "median"`.
- `ref.g = "median"`: finite fitted `eval_g` value closest to the median of
  finite fitted `eval_g`.
- `ref.g = "min"` or `"lowest"`: minimum finite fitted `eval_g`.
- `ref.g = "max"` or `"highest"`: maximum finite fitted `eval_g`; optional but
  easy to include.
- `ref.g = "mean"`: fitted `eval_g` value closest to the mean of finite fitted
  `eval_g`, so it remains on the fitted grid.
- Numeric scalar: preserve current behavior; if off-grid, it is allowed only
  when `interpolate = "linear"` can evaluate it within the fitted range.
- In quantile target mode, `ref.g = NULL` should default to the first selected
  quantile target, matching the requested "Q1 as reference" behavior.

Optional future extension:

```r
ref.source = c("eval_g", "support")
```

This can be deferred. For the first refinement, compute character reference
rules from `eval_g` to avoid off-range and missing-support edge cases.

## Proposed Function Signature Changes

Update `plot.fdid()` and `.fdid_contrast_plot_data()` with:

```r
include.ref = TRUE,
target.type = c("level", "quantile"),
quantile.source = c("support", "eval_g")
```

Existing arguments remain:

```r
ref.g = NULL,
target.g = NULL,
contrast.values = NULL,
n.contrast = 7L,
max.contrast = 10L,
inference = c("auto", "replicate", "vcov", "band", "pointwise"),
interpolate = c("linear", "none")
```

Recommended count semantics after adding a default reference row:

- Treat `n.contrast` as the desired total number of displayed coefficient rows
  in automatic and quantile modes, including the reference row.
- Treat `max.contrast` as the automatic display cap including the reference
  row.
- For explicit user-supplied level targets, add the reference row if
  `include.ref = TRUE`; if the total row count exceeds `max.contrast`, warn but
  honor the user-supplied display.

This is slightly different from the first implementation, where `n.contrast`
counted non-reference target rows. The documentation and tests must call out
the change.

## Implementation Plan

### R Code

Update `fdid-main/R/plot.R`:

1. Add helper `.fdid_resolve_ref_g()`:
   - inputs: `x`, `evg`, `ref.g`, optional resolved target values,
     `target.type`;
   - handles numeric and character reference rules;
   - for quantile mode with `ref.g = NULL`, uses the first resolved quantile
     target.
2. Add helper `.fdid_resolve_contrast_targets()`:
   - supports level targets, quantile targets, defaults, duplicate removal,
     range validation, and display labels;
   - returns `target_g`, `label`, optional `target_rank`, and
     `target_source`.
3. Add helper `.fdid_reference_contrast_row()`:
   - builds the zero reference row without calling pointwise fallback
     inference.
4. Update `.fdid_contrast_plot_data()`:
   - accept `include.ref`, `target.type`, and `quantile.source`;
   - resolve targets before reference when quantile mode needs first target as
     default reference;
   - insert the reference row when requested;
   - maintain `target_g`, `ref_g`, `reference_g`, `is_reference`,
     `target_estimand`, `label`, and inference metadata columns.
5. Update `.fdid_plot_contrast_base()`:
   - compute inference label from non-reference rows;
   - optionally draw the reference row with a distinct color or point symbol;
   - keep `ci = FALSE` behavior.
6. Update `plot.fdid()`:
   - expose the new arguments;
   - pass them through to `.fdid_contrast_plot_data()`.

### Tests

Update `fdid-main/tests/testthat/test-continuous.R`:

- default `plot(type = "contrast")` includes exactly one `is_reference` row;
- the reference row has estimate, SE, and CI limits equal to zero;
- non-reference rows still inherit inference from `fdid_contrast()`;
- the top-level plot call still returns `invisible(NULL)`;
- `include.ref = FALSE` restores old no-reference-row behavior;
- `ref.g = "median"`, `"mean"`, `"min"`, and numeric values work;
- unsupported character `ref.g` errors clearly;
- `target.type = "quantile", target.g = 5` returns five displayed rows
  including the first quantile as reference;
- quantile mode labels are stable and duplicates are handled;
- `target.type = "level"` preserves absolute level behavior;
- `target.g` and `contrast.values` conflict remains an error for level mode;
- scalar DML methods still reject `type = "contrast"`.

### Documentation

Update:

- `fdid-main/R/plot.R` roxygen;
- `fdid-main/man/plot.fdid.Rd`;
- `fdid-main/CONTINUOUS_G_ESTIMATOR_REFERENCE.md`;
- `fdid-main/tutorial-next-2026-05-30/04-visualization.Rmd`;
- `fdid-main/CONTINUOUS_G_PROJECT_LOG.md`;
- `fdid-main/TODO.md`;
- `fdid-main/PACKAGE_PLOTTING_LAYER_UPDATE_PLAN_2026-05-31.md`.

Documentation should state:

- the current default reference rule;
- how `include.ref` changes the coefficient plot;
- the difference between absolute level targets and quantile targets;
- that quantile mode with `target.g = 5` displays five quantile-selected
  coefficient rows including the reference;
- that non-reference intervals still come from `fdid_contrast()`;
- that the reference row is an algebraic zero anchor, not an independently
  estimated contrast.

`FDID Notes` only needs an update if its plotting or package-interface section
mentions contrast-plot defaults. The statistical theory does not change.

### Tutorial Examples

Add examples such as:

```r
plot(fdid_results_flex_curve, type = "contrast")

plot(fdid_results_flex_curve, type = "contrast",
     ref.g = "min")

plot(fdid_results_flex_curve, type = "contrast",
     target.type = "quantile", target.g = 5)

plot(fdid_results_flex_curve, type = "contrast",
     ref.g = 0.5, target.g = c(0.2, 0.35, 0.65, 0.8))
```

The tutorial should explain that the first displayed quantile is the reference
when `target.type = "quantile"` and `ref.g` is left at `NULL`.

## Inference Implications

This refinement should not change contrast inference for non-reference rows.
All non-reference intervals must still come from `fdid_contrast()`.

The reference row is different:

- its estimand is exactly zero by construction;
- its interval should be exactly zero-width;
- it should be labeled as a reference anchor;
- it should be excluded from inference-source summaries.

This keeps the plot visually useful without weakening the package's existing
inference discipline.

## Validation Plan

Minimum validation after implementation:

```sh
Rscript -e 'r_files <- list.files("fdid-main/R", pattern = "[.]R$", full.names = TRUE); invisible(lapply(r_files, parse)); rd_files <- list.files("fdid-main/man", pattern = "[.]Rd$", full.names = TRUE); invisible(lapply(rd_files, tools::checkRd)); cat("parse/checkRd ok\n")'
R CMD INSTALL -l /private/tmp/fdidlib fdid-main
Rscript -e '.libPaths(c("/private/tmp/fdidlib", .libPaths())); setwd("fdid-main"); testthat::test_file("tests/testthat/test-continuous.R")'
Rscript -e '.libPaths(c("/private/tmp/fdidlib", .libPaths())); devtools::test("fdid-main")'
env R_LIBS=/private/tmp/fdidlib quarto render fdid-main/tutorial-next-2026-05-30
git diff --check -- fdid-main/R/plot.R fdid-main/man/plot.fdid.Rd fdid-main/tests/testthat/test-continuous.R fdid-main/CONTINUOUS_G_ESTIMATOR_REFERENCE.md fdid-main/CONTINUOUS_G_PROJECT_LOG.md fdid-main/TODO.md fdid-main/PACKAGE_PLOTTING_LAYER_UPDATE_PLAN_2026-05-31.md fdid-main/tutorial-next-2026-05-30/04-visualization.Rmd
```

After tests, restore generated plot artifacts such as
`fdid-main/tests/testthat/Rplots.pdf` if they are modified only by validation.

## Recommended First Code Pass

Implement the refinement in this order:

1. Add reference-row support with `include.ref = TRUE`.
2. Add character reference strategies for `ref.g`.
3. Add quantile target mode.
4. Update tests.
5. Update roxygen/Rd, estimator reference, tutorial, TODO, project log, and
   the package-wide plotting plan.

This keeps the first refinement focused on the contrast coefficient plot. It
does not start the broader package-wide ggplot/autoplot modernization.
