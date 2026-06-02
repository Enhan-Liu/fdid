# Continuous-G Estimator Reference

Timestamp: 2026-06-01, America/Chicago

This is the current technical reference for the `fdid` package after the
continuous-G update. It describes the estimators, options, algorithms, output
objects, and current SE/CI support.

For restart context, begin with the local `HANDOFF.md` note, then read
`TODO.md` and `CONTINUOUS_G_PROJECT_LOG.md` alongside this reference.
`HANDOFF.md` is excluded from package builds through `.Rbuildignore`.

## Public Entry Points

Main functions:

```r
fdid_prepare()
fdid_expand_covariates()
fdid()
fdid_contrast()
fdid_derivative()
summary()
plot()
fdid_list()
```

Continuous-G estimators enter through:

```r
fdid(..., method = "kernel")
fdid(..., method = "dml_plr")
fdid(..., method = "dml_flex")
fdid(..., method = "dml_incremental")
```

The DML binary method also enters through:

```r
fdid(..., method = "dml_binary")
```

No new estimator is exposed as a separate top-level estimation function.

Reporting helpers:

```r
fdid_contrast(result, g0, g1)
fdid_derivative(result, g0)
```

These helpers report Notes-style estimands from stored continuous-G curve
objects. `fdid_contrast()` reports level comparisons such as
`theta(g1) - theta(g0)` for `dml_flex` or `mu(g1) - mu(g0)` for `kernel`.
`fdid_derivative()` reports fixed-`g` derivative values such as
`theta_prime(g0)` or `delta(g0)`. Both helpers use stored `eval_g` rows, can
linearly interpolate between fitted grid points, and return a diagnostic note
describing the inference source.

## Continuous-G Inference Reference

This section is the package-level map for continuous-G standard errors,
confidence intervals, simultaneous bands, summaries, and reporting helpers. It
is organized by output object because the same fitted estimator can support
several inferential targets.

### Object 1: `summary()` Scalar Rows

`summary.fdid()` first prints scalar aggregate rows because the original `fdid`
object structure is organized around pre-event, event, post-event, and dynamic
tables. For continuous-G estimators these scalar rows should be read as compact
summaries of the curve object, not as substitutes for the stored curve.

Current scalar targets:

| Method | `summary()` scalar target | Primary output |
|---|---|---|
| `kernel` | interval-average level slope over the fitted `eval_g` range | `curve_event$mu_hat`, `curve_event$delta_hat` |
| `dml_flex` | interval-average level slope over the fitted `eval_g` range | `curve_event$theta_hat`, `curve_event$delta_hat` |
| `dml_plr` | partially linear residualized slope | scalar DML score output |
| `dml_incremental` | observed-population average derivative | scalar DML average-derivative signal |
| `dml_binary` | binary AIPW/IRM contrast | scalar DML AIPW output |

For `kernel` and `dml_flex`, the scalar event row is:

```text
level_slope = {level(max eval_g) - level(min eval_g)}
              / {max eval_g - min eval_g}
```

where `level(g)` is `mu(g)` for `kernel` and `theta(g)` for `dml_flex`.
The row is useful for quick scalar comparison across methods, but reporting
should still lead with the level curve, derivative curve, or targeted
`fdid_contrast()` / `fdid_derivative()` call when those are the research target.

`summary()` prints the scalar event SE and CI method when the fitted row stores
them. These labels are part of the interpretation:

| Label | Meaning |
|---|---|
| `analytical_stacked_sandwich` | kernel analytical cross-grid sandwich covariance |
| `analytical_cluster_stacked_sandwich` | kernel cluster-summed cross-grid sandwich covariance |
| `bootstrap_replicate` | scalar computed in each bootstrap curve replication |
| `bootstrap_curve_vcov` | scalar normal CI from bootstrap curve covariance |
| `curve_vcov` | DML-flex interval-average normal CI from stored `theta_vcov` |
| `endpoint_independent_fallback` | pointwise endpoint SEs combined without covariance |

The endpoint-independent fallback is an approximation and should not be treated
as equivalent to covariance-aware contrast inference.

### Object 2: Kernel Curve Inference

The kernel estimator stores a level curve `mu(g)` and derivative curve
`delta(g)` over `eval_g`.

With `vartype = "robust"`, the package now stores analytical cross-grid curve
covariance matrices:

```r
result$curve_event$mu_vcov
result$curve_event$delta_vcov
result$curve_event$curve_vcov_method
```

The covariance is computed from stacked local-linear estimating equations. For
grid points `g_j` and `g_k`,

```text
Cov(theta_hat(g_j), theta_hat(g_k))
  = A_j^{-1} [sum_i psi_ij psi_ik'] A_k^{-1},

psi_ij = W_ij R_ij e_ij.
```

The level and derivative covariance matrices are then obtained by linear
projection:

```text
Cov(mu_hat(g_j), mu_hat(g_k))     = a' Cov(theta_j, theta_k) a
Cov(delta_hat(g_j), delta_hat(g_k)) = b' Cov(theta_j, theta_k) b
```

If a prepared cluster column `c` is present, the internal sandwich meat is based
on cluster-summed scores. Otherwise it is unit-level HC0-style stacked
sandwich. The pointwise `se_mu` and `se_delta` are the square roots of the
diagonals of these matrices.

With `vartype = "bootstrap"`, the package stores bootstrap replicate curves and
their empirical covariance matrices:

```r
result$curve_event$mu_bootstrap
result$curve_event$delta_bootstrap
result$curve_event$mu_boot_vcov
result$curve_event$delta_boot_vcov
```

Bootstrap inference resamples rows or clusters, reselects the baseline
bandwidth, refits the whole curve on the stored `eval_g` grid, and stores
pointwise percentile intervals plus finite-grid quantile-envelope simultaneous
bands. These bands follow the practical interflex-style quantile-envelope
construction; they are not the studentized sup-t band.

### Object 3: Scalar DML Inference

The scalar DML estimators are `dml_binary`, `dml_plr`, and `dml_incremental`.
They use cross-fitting and median aggregation over `S` random splits.

For each split, the package forms a Neyman-orthogonal score or signal and
computes a score variance. Across splits it reports the median estimate and
adds split variability:

```text
SE = sqrt(median_s { var_s + (estimate_s - median_estimate)^2 })
```

`dml_inference = "score"` is the default analytical score interval.
`dml_inference = "score_multiplier"` keeps the same scalar SE but replaces the
normal critical value with a multiplier critical value computed from stored
score contributions. This is not a nonparametric refit bootstrap.

### Object 4: DML-Flex Curve Inference

`dml_flex` first builds cross-fitted continuous-treatment pseudo-outcomes and
then maps those signals back to `G` over `eval_g`. The mapping choice is
controlled by `signal_map`, not by `density_method`.

| `signal_map` | Mapping role | Pointwise SE source | Stored covariance | Band route |
|---|---|---|---|---|
| `local_poly` | local polynomial smoother of signal on `G` | local polynomial sandwich | no `theta_vcov` | practical signal-residual multiplier |
| `kernel` | local linear smoother alias | local polynomial sandwich | no `theta_vcov` | practical signal-residual multiplier |
| `spline` | B-spline linear model | model-matrix delta method | no grid `theta_vcov` currently | practical signal-residual multiplier |
| `gam` | `mgcv::gam()` smooth | finite-difference design rows | no grid `theta_vcov` currently | practical signal-residual multiplier |
| `blp_spline` | B-spline best linear projection | BLP sandwich covariance | `theta_vcov`, `delta_vcov`, `level_derivative_vcov` | Gaussian max-t from covariance |

For non-BLP maps, `theta_multiplier_reps` and `delta_multiplier_reps` store
practical multiplier replicate curves. These support replicate-based contrast
reporting, but should be described as signal-mapping replicate inference rather
than full DML refit bootstrap.

For `signal_map = "blp_spline"`, the second-stage map is:

```text
theta_hat(eval_g) = B_grid beta_hat
theta_vcov = B_grid V_beta B_grid'
delta_vcov = D_grid V_beta D_grid'
```

The diagonal of the averaged covariance matrix is rescaled to match the
repeated-split pointwise standard errors. Gaussian max-t critical values are
then simulated from the stored covariance matrices for finite-grid bands. This
is the closest current path to the Practical Guide/interflex BLP logic, but it
remains a second-stage grid covariance approximation, not full theorem-level
DML uniform inference over all nuisance stages.

### Object 5: `fdid_contrast()` And `fdid_derivative()`

`fdid_contrast(x, g0, g1)` reports level changes:

```text
mu(g1) - mu(g0)       for kernel
theta(g1) - theta(g0) for dml_flex
```

`fdid_derivative(x, g0)` reports fixed-G derivative values:

```text
delta(g0)       for kernel
theta_prime(g0) for dml_flex
```

Both helpers use `eval_g` and can either require exact grid points
(`interpolate = "none"`) or linearly interpolate between grid points
(`interpolate = "linear"`).

With `inference = "auto"`, the helper tries:

```text
vcov -> replicate -> band -> pointwise
```

The covariance route uses a linear-contrast quadratic form:

```text
SE = sqrt(w' V w)
```

For exact-grid two-point level contrasts this reduces to:

```text
Var(level(g1) - level(g0)) = V11 + V00 - 2 V10
```

The replicate route computes the target in every stored replicate curve:

```text
contrast_b = level_b(g1) - level_b(g0)
SE = sd(contrast_b)
CI = percentile(contrast_b)
```

Thus replicate contrast and `theta_vcov` are two representations of joint curve
uncertainty. Replicates preserve distributional shape and can yield percentile
intervals; covariance matrices give a compact second-moment representation that
supports interpolation, Wald intervals, scalar summaries, and Gaussian max-t
simulation.

### Object 6: Current Limitations And Planned Improvements

The current inference layer intentionally labels practical approximations.

- Kernel robust now has analytical cross-grid covariance for curve values and
  summary/contrast inference.
- Kernel bootstrap has replicate curves, bootstrap curve covariance, pointwise
  percentile intervals, and finite-grid quantile-envelope bands.
- `dml_flex + blp_spline` has stored grid covariance and Gaussian max-t bands.
- Non-BLP `dml_flex` maps rely on practical multiplier replicate curves for
  contrast inference and practical signal-mapping bands for curves.
- A future extension can add replicate-derived covariance matrices for non-BLP
  `dml_flex`, or map-specific analytical covariance for `spline`, `gam`,
  `local_poly`, and `kernel` maps.
- Full DML refit bootstrap and theorem-level uniform bands for all `dml_flex`
  paths remain research-level work.

## Prepared Data

`fdid_prepare()` takes long panel data and creates a wide unit-level data frame.

Important output columns:

| Column pattern | Meaning |
|---|---|
| `unit` | unit identifier |
| `Y_<time>` | outcome at each time |
| `G` | baseline factor or treatment intensity |
| `x1`, `x2`, ... | averaged or expanded covariates |
| `c` | optional cluster identifier |

The common event-period change score used by estimators is:

```text
DeltaY_i = mean(Y_i,t over tr_period) - Y_i,ref
```

Dynamic estimates replace the event-period mean with each period-specific
change from the reference period.

## High-Dimensional Covariate Expansion

Covariate expansion can be requested inside `fdid_prepare()`:

```r
s_hd <- fdid_prepare(
  data = my_panel,
  Y_label = "outcome",
  X_labels = c("z1", "z2", "z3"),
  G_label = "G",
  unit_label = "id",
  time_label = "time",
  basis_type = "polynomial",
  poly_degree = 2,
  include_interactions = TRUE
)
```

or after preparation:

```r
s_hd <- fdid_expand_covariates(
  s,
  basis_type = "polynomial",
  poly_degree = 2,
  include_interactions = TRUE
)
```

Options:

| Argument | Default | Description |
|---|---|---|
| `basis_type` | `none` | `none`, `polynomial`, or `bspline` |
| `include_interactions` | `FALSE` | add all pairwise interactions among expanded columns |
| `poly_degree` | `2` | maximum polynomial degree |
| `spline_df` | `4` | B-spline degrees of freedom |
| `spline_degree` | `2` | B-spline degree |
| `standardize_basis` | `TRUE` | center and scale expanded columns |

Expansion affects only the `x1`, `x2`, ... covariate columns. It never changes
the factor/treatment variable `G`.

Implementation details:

- `fdid_prepare()` pivots outcomes to `Y_<time>` columns;
- time-varying covariates are averaged by unit before estimation;
- `G_label` is renamed to `G`;
- `cluster_label`, when supplied, is renamed to `c`;
- `fdid_expand_covariates()` finds columns named `x1`, `x2`, ... and replaces
  them with expanded `x1`, `x2`, ... columns;
- polynomial expansion uses powers from 1 through `poly_degree`;
- B-spline expansion falls back to the raw covariate if support is too small or
  spline construction fails;
- pairwise interactions are formed after the main-effect basis is built;
- nonfinite or near-constant expanded columns are dropped;
- expanded columns are standardized unless `standardize_basis = FALSE`;
- expansion metadata is stored in `attr(prepared_data, "fdid_basis")`.

## Method Map

| Method | G type | Main target | Main output |
|---|---|---|---|
| `did` | binary or continuous | no-covariate slope/contrast | scalar |
| `ols1` | binary or continuous | additive regression slope/contrast | scalar |
| `ols2` | binary or continuous | interacted regression slope/contrast at target covariate population | scalar |
| `ebal` | binary | entropy-balanced ATT | scalar |
| `ipw` | binary | IPW contrast | scalar |
| `aipw` | binary | AIPW/causal-forest contrast | scalar |
| `kernel` | continuous | local-linear level and derivative curves | curve |
| `dml_binary` | binary | cross-fitted AIPW/IRM contrast | scalar |
| `dml_plr` | continuous | partially linear scalar slope | scalar |
| `dml_flex` | continuous | level curve `theta(g)` and derivative curve `theta'(g)` | curve |
| `dml_incremental` | continuous | scalar `E[partial_g mu(G, X)]` | scalar |

## Common Arguments

| Argument | Methods | Meaning |
|---|---|---|
| `tr_period` | all | treatment/event periods |
| `ref_period` | all | reference period |
| `entire_period` | all | periods used for dynamic output |
| `method` | all | estimator choice |
| `vartype` | non-DML and kernel | variance type; DML currently robust only |
| `missing_data` | all | `listwise` or `available` |
| `nsims` | bootstrap/jackknife original methods | resampling iterations |
| `parallel`, `cores` | bootstrap/jackknife original methods | parallel execution |
| `target.pop` | original binary methods | `all`, `1`, or `0` |
| `alpha` | continuous-G and DML | significance level |

## Continuous-G Arguments

| Argument | Default | Methods | Meaning |
|---|---|---|---|
| `eval_g` | `NULL` | `kernel`, `dml_flex` | G-grid for curve evaluation |
| `h0` | `NULL` | `kernel` | baseline bandwidth |
| `K_folds` | `10` | `kernel` | bandwidth CV folds |
| `boot` | `500` | `kernel`, `dml_flex` | bootstrap or multiplier replications for bands and stored replicate curves |
| `trim` | `0.05` | `kernel`, `dml_flex` | tail trimming for default grid |
| `learner` | `linear` | DML | nuisance learner |
| `K` | `5` | DML | cross-fitting folds |
| `S` | `3` | DML | random splits for median aggregation |
| `signal_map` | `local_poly` | `dml_flex` | maps orthogonal signals to `G` |
| `map_degree` | `2` | `dml_flex` | polynomial/spline degree |
| `map_df` | `NULL` | `dml_flex` | spline/GAM basis size |
| `density_method` | `residual_kde` | `dml_flex`, `dml_incremental` | conditional density model |
| `dml_inference` | `score` | scalar DML | scalar interval rule: analytical score or score multiplier |
| `dml_boot` | `500` | scalar DML | multiplier draws for scalar score multiplier inference |
| `dml_multiplier` | `normal` | scalar DML | multiplier law: `normal`, `wild`, or `exponential` |

## Curve Grids, Estimation, And Plotting

`eval_g` is an estimation grid for the continuous-curve methods, not a
plotting-only option. For `kernel` and `dml_flex`, the package computes and
stores curve quantities only at the values in `eval_g`; `plot(type = "curve")`
then displays those stored rows.

Default FDID behavior:

- `kernel`: when `eval_g = NULL`, the package builds 50 equally spaced points
  from `quantile(G, trim)` to `quantile(G, 1 - trim)`, with `trim = 0.05` by
  default;
- `dml_flex`: uses the same 50-point trimmed default grid when `eval_g = NULL`;
- `dml_plr` and `dml_incremental`: scalar estimators; `eval_g` does not define
  a curve for these methods.

Estimator roles:

- `kernel` runs one local-linear weighted least squares fit at each
  `g0 in eval_g`. Supplying five values means five local fits and five stored
  curve points.
- `dml_flex` first constructs orthogonal pseudo-outcomes for all observations,
  then maps those signals back to `G` only at `eval_g` through `signal_map`.
  Thus all observations contribute to the signal, but the reported
  `theta_hat(g)` and `delta_hat(g)` curves are grid-based.

Plotting roles:

- `plot.fdid(type = "curve")` does not re-estimate, densify, or interpolate a
  new curve. It reads `result$eval_g` and the matching columns in
  `result$curve_event`.
- The lower support display is descriptive only. It uses observed `G` values
  restricted to the displayed curve window so the histogram or density aligns
  with the plotted curve range.

Comparison to `interflex`:

- `interflex` kernel is also grid-based. Its `X.eval` is user-supplied or,
  by default, `seq(min(X), max(X), length.out = neval)` with `neval = 50`
  unless an explicit plotting window is used; it then loops over `X.eval` and
  runs a local WLS fit at each grid point.
- `interflex` DML extracts pseudo-outcomes for all observations, then projects
  them onto a B-spline basis and predicts the CATE on a grid, with
  `.compute_cate_blp(..., n_grid = 50)` by default.
- The practical difference is not "FDID estimates only a few points while
  interflex estimates all data points." Both are grid-based for continuous
  curves by default. The main differences are the default range and user
  choices: FDID defaults to 50 points on trimmed support, while `interflex`
  defaults to 50 points over the observed min-max support unless `X.eval` or
  `xlim` changes that range. A sparse user-supplied `eval_g` in FDID produces
  a sparse estimated curve.

## Variance Support

| Method family | `robust` | `bootstrap` | `jackknife` |
|---|---|---|---|
| `did`, `ols1`, `ols2`, `ebal`, `ipw`, `aipw` | yes | yes | yes |
| `kernel` | yes | yes | no |
| `dml_binary`, `dml_plr`, `dml_flex`, `dml_incremental` | yes | no | no |

DML methods reject unsupported variance requests.

## Package Wiring And Validation

This section records code-level behavior from the current R implementation.

Package status:

- package version in `DESCRIPTION` is `1.0.2`;
- exported functions are `fdid()`, `fdid_prepare()`,
  `fdid_expand_covariates()`, and `fdid_list()`;
- S3 methods are registered for `print.fdid()`, `summary.fdid()`,
  `plot.fdid()`, and `plot.fdid_list()`;
- continuous-G estimators are internal helpers called by `fdid()`;
- optional ML packages are loaded lazily when their learner is requested.

Important dependency behavior:

- `fdid()` currently checks for `estimatr` and `car` before method dispatch,
  even for continuous-G and DML methods;
- `plot.fdid()` checks for `RColorBrewer` when plotting;
- original `ipw` and `aipw` require `grf`;
- `ebal` requires `ebal`;
- clustered `ols2` variance requires `sandwich`;
- optional DML learners require their packages only when requested:
  `glmnet`, `ranger`, `grf`, `mgcv`, `nnet`, or `xgboost`.

Early validation in `fdid()`:

- validates `method` and `vartype`;
- rejects `vartype = "jackknife"` for `kernel`;
- rejects any non-`robust` `vartype` for DML;
- requires a `G` column;
- validates `alpha`, `trim`, `eval_g`, `h0`, `boot`, `K_folds`, `K`, and
  `S`;
- checks that treatment and reference outcome columns exist;
- maps multiple reference periods to `max(ref_period)` before continuous-G
  dispatch.

Continuous-G dispatch:

- `method = "kernel"` dispatches to `run_kernel_method()`;
- DML methods dispatch to `run_dml_method()`;
- original binary/scalar methods continue through the legacy `fdid.R` path;
- `target.pop != "all"` produces a warning for `did` and `ols1`, because it
  does not change those estimates.

Missing-data behavior:

- original methods include `G`, covariates, and cluster column `c` in the
  minimal complete-case set;
- `kernel` and DML dispatchers complete-case only `G` and `x1`, `x2`, ...
  before reading the optional cluster vector;
- for `missing_data = "listwise"`, continuous-G dispatchers also require all
  outcome columns in `all_times`;
- for `missing_data = "available"`, continuous-G dispatchers allow partially
  missing outcomes and drop observations inside each period-specific call;
- dynamic and pre/post DML calls require at least `max(10, K)` nonmissing
  observations; kernel calls require at least 10.

Display-only raw means:

- for `kernel` and DML continuous-G outputs, `raw_means` are based on a median
  split of continuous `G`;
- those raw means are only descriptive display objects and are not estimator
  targets.

## Original Estimators

### `method = "did"`

Target:

```text
DeltaY_i = alpha + tau G_i + error_i
```

Implementation:

- fits `tempY ~ G`;
- uses `estimatr::lm_robust()`;
- cluster-robust version is used if prepared data contain `c`;
- returns scalar estimate, SE, and CI.

For binary `G`, this is a no-covariate difference-in-differences contrast. For
continuous `G`, it is a no-covariate slope summary.

### `method = "ols1"`

Target:

```text
DeltaY_i = alpha + tau G_i + X_i beta + error_i
```

Implementation:

- fits `tempY ~ G + x1 + x2 + ...`;
- uses `estimatr::lm_robust()`;
- cluster-robust version is used if clusters are supplied;
- returns the coefficient on `G`.

### `method = "ols2"`

Target:

```text
DeltaY_i = alpha + tau G_i + X_i beta + G_i X_i gamma + error_i
```

Implementation:

- scales covariates;
- fits `lm(Y ~ Z * X)`;
- estimates the marginal effect of `G` at a target covariate mean;
- `target.pop = "all"` uses zero after scaling;
- `target.pop = "1"` uses treated-group covariate mean;
- `target.pop = "0"` uses control-group covariate mean;
- uses `car::hccm()` for robust variance or `sandwich::vcovCL()` for clustered
  variance;
- adds a Lin-style super-population correction when `target.pop = "all"`.

### `method = "ebal"`

Target:

```text
ATT for binary G
```

Implementation:

- requires `target.pop = "1"`;
- uses `ebal::ebalance()` to weight controls to match treated covariates;
- fits weighted `tempY ~ G` with `estimatr::lm_robust()`;
- supports robust, bootstrap, and jackknife through the original variance
  wrappers.

### `method = "ipw"`

Target:

Binary-G weighted contrast for `target.pop = "all"`, `"1"`, or `"0"`.

Implementation:

- estimates propensity scores using `grf::probability_forest()` when
  covariates are present;
- uses marginal treatment probability when no covariates are present;
- clips probabilities away from 0 and 1;
- constructs ATE, ATT, or ATC weights;
- fits weighted `tempY ~ G` with `estimatr::lm_robust()`.

### `method = "aipw"`

Target:

Binary-G AIPW contrast.

Implementation:

- uses `grf::causal_forest()`;
- calls `grf::average_treatment_effect()`;
- maps `target.pop` to `all`, `treated`, or `control`;
- returns estimate, standard error, and normal CI.

## Kernel Continuous-G Estimator

### Target

The kernel estimator targets a local-linear DID-response curve:

```text
mu(g) = E[DeltaY | G = g, adjusted for X through local-linear regression]
```

and its local derivative:

```text
delta(g) = d mu(g) / dg
```

### Algorithm

For each evaluation point `g0`, fit a weighted local-linear model with Gaussian
kernel weights:

```text
K((G_i - g0) / h_g0) / h_g0
```

The local design contains:

```text
1
X
G - g0
(G - g0) * X
```

The level is evaluated at `Xbar`:

```text
mu_hat(g0) = a(Xbar)' theta_hat(g0)
```

The derivative is:

```text
delta_hat(g0) = b(Xbar)' theta_hat(g0)
```

Bandwidth:

- if `h0` is supplied, use it;
- otherwise choose baseline `h0` by K-fold least-squares cross-validation;
- adapt bandwidths using pilot density of `G`;
- default `K_folds = 10`.

Code-level details:

- at least 10 complete observations are required;
- continuous `G` must have more than two unique values;
- the pilot density uses `density(G, bw = 1.06 * sd(G) * n^(-1/5), n = 1024)`;
- bandwidth search uses 25 log-spaced candidates;
- local fits return `NA` if all kernel weights are negligible, effective local
  sample size is below 5, or the weighted normal matrix is singular;
- cross-validation penalizes bandwidths that fail on more than half of
  held-out observations;
- analytical kernel curve inference uses a stacked sandwich covariance over the
  fitted grid; when a cluster column is available, the sandwich meat is
  cluster-summed;
- when `vartype = "bootstrap"` and a cluster column is available, kernel
  bootstrap resampling is cluster-aware.

Default grid:

- if `eval_g = NULL`, use 50 estimation points over the trimmed support of
  `G`;
- trimming uses `trim = 0.05` by default.

### Output Fields

```r
result$eval_g
result$curve_event$mu_hat
result$curve_event$delta_hat
result$curve_event$se_mu
result$curve_event$se_delta
result$curve_event$h_adaptive
result$curve_event$mu_vcov
result$curve_event$delta_vcov
result$curve_event$curve_vcov_method
result$curve_event$mu_band_lower
result$curve_event$mu_band_upper
result$curve_event$delta_band_lower
result$curve_event$delta_band_upper
result$curve_event$mu_bootstrap        # when vartype = "bootstrap"
result$curve_event$delta_bootstrap     # when vartype = "bootstrap"
result$curve_event$mu_boot_se          # when vartype = "bootstrap"
result$curve_event$delta_boot_se       # when vartype = "bootstrap"
result$curve_event$mu_boot_ci_lower    # when vartype = "bootstrap"
result$curve_event$mu_boot_ci_upper    # when vartype = "bootstrap"
result$curve_event$delta_boot_ci_lower # when vartype = "bootstrap"
result$curve_event$delta_boot_ci_upper # when vartype = "bootstrap"
result$curve_event$mu_boot_vcov        # when vartype = "bootstrap"
result$curve_event$delta_boot_vcov     # when vartype = "bootstrap"
result$band_lower
result$band_upper
result$h0
```

### SE And CI

Pointwise and cross-grid SEs:

- computed by a stacked local-linear sandwich over the fitted `eval_g` grid;
- `se_mu` for level;
- `se_delta` for derivative;
- `mu_vcov` and `delta_vcov` store the analytical cross-grid covariance
  matrices used by covariance-aware robust contrasts and interval-average
  summaries.

Inference fields:

- if `vartype = "robust"`, `band_lower` and `band_upper` for derivatives are
  pointwise normal intervals;
- if `vartype = "bootstrap"`, the event curve stores `mu_bootstrap` and
  `delta_bootstrap` matrices for downstream contrast and derivative reporting;
- bootstrap pointwise uncertainty is stored as `mu_boot_se`,
  `delta_boot_se`, and percentile interval fields;
- bootstrap simultaneous bands are stored separately for the level and
  derivative curves as `mu_band_*` and `delta_band_*`;
- the implemented bootstrap simultaneous band is an interflex-style
  quantile-envelope band over the evaluation grid, not a fully studentized
  sup-t band;
- `mu_boot_vcov` and `delta_boot_vcov` store grid-by-grid bootstrap covariance
  matrices for covariance-aware contrasts.

Implemented quantile-envelope algorithm:

1. For each bootstrap replication, resample rows or clusters, reselect the
   baseline bandwidth, and refit the kernel curve on the stored `eval_g` grid.
2. For each curve target, collect a bootstrap matrix with rows equal to
   bootstrap draws and columns equal to evaluation points.
3. Search over `zeta` in `[alpha / (2 * length(eval_g)), alpha / 2]`.
4. At each grid point, compute the empirical `zeta` and `1 - zeta` quantiles of
   the bootstrap estimates.
5. Choose the largest calibrated `zeta` whose rectangular envelope contains
   whole bootstrap curves with empirical probability at least `1 - alpha`.
6. Store the resulting envelope as `mu_band_*` or `delta_band_*`, along with
   `*_uniform_zeta`, `*_uniform_coverage`, and `*_uniform_n_boot_eff`.

This is the same practical simultaneous-inference target as the current
`interflex` kernel implementation and the CME Practical Guide's described
bootstrap quantile-envelope procedure. It should be reported as a finite-grid
bootstrap quantile-envelope band, not as a studentized max-t/sup-t band.

Kernel bootstrap resampling:

- row bootstrap is used when no cluster labels are available;
- cluster bootstrap is used when the prepared data include a cluster column
  `c`;
- the fitted object records `bootstrap_resampling` and `n_clusters`.

Scalar event summary:

- if the grid has one point, report derivative at that point;
- otherwise report an interval-average slope:

```text
(mu_hat(max eval_g) - mu_hat(min eval_g)) / (max eval_g - min eval_g)
```

Scalar inference priority:

- when bootstrap replicate curves are available, the scalar event summary uses
  bootstrap replicate contrasts and percentile intervals;
- when a bootstrap curve covariance matrix is available, endpoint covariance is
  incorporated through `V22 + V11 - 2 V12`;
- when analytical curve covariance is available, the robust scalar event summary
  uses the same cross-grid covariance formula;
- otherwise the scalar SE falls back to the endpoint-independent pointwise
  approximation and labels the method as a fallback.

## DML Learners

DML nuisance learners are used for outcome regressions, treatment regressions,
propensity models, and density submodels.

| Normalized learner | Regression backend | Classification backend |
|---|---|---|
| `linear` | `lm()` | binomial `glm()` |
| `lasso` | `glmnet::cv.glmnet(alpha = 1)` | binomial glmnet |
| `ridge` | `glmnet::cv.glmnet(alpha = 0)` | binomial glmnet |
| `elasticnet` | `glmnet::cv.glmnet(alpha = 0.5)` | binomial glmnet |
| `ranger` | `ranger::ranger()` | probability forest |
| `grf` | `grf::regression_forest()` | regression forest on 0/1 outcome |
| `gam` | `mgcv::gam()` | binomial `mgcv::gam()` |
| `nnet` | `nnet::nnet()` | entropy neural net |
| `xgboost` | `xgboost::xgb.train()` | binary logistic xgboost |

All optional package dependencies are checked only when requested.

Learner aliases:

- `linear`: `default`, `logistic`, `lm`, `ols`, `l`, `d`;
- `lasso`: `glmnet`, `regularization`, `regularisation`, `r`;
- `elasticnet`: `enet`;
- `ranger`: `rf`, `randomforest`, `forest`;
- `grf`: `regressionforest`, `probabilityforest`;
- `gam`: `mgcv`;
- `nnet`: `nn`, `network`, `neuralnetwork`;
- `xgboost`: `xgb`, `boosting`, `boost`, `gradientboosting`,
  `gradientboost`, `histgradientboosting`, `hgb`.

Backend details:

- `glmnet` uses `cv.glmnet()` and sets its internal folds to
  `max(3, min(5, floor(n / 2)))`;
- one-column `glmnet` designs are padded with a zero column to avoid backend
  failures;
- `ranger` regression uses 200 trees and `min.node.size = 5`;
- `ranger` classification uses probability forests with 200 trees;
- `grf` uses 500-tree regression forests, including for 0/1 probability
  prediction;
- `gam` uses `mgcv::gam(..., method = "REML")`; continuous variables with at
  least six unique values get smooth terms with `k <= 5`;
- `nnet` standardizes features, uses `size = min(5, p)`, `decay = 0.01`,
  `maxit = 500`, and `MaxNWts = 10000`;
- `xgboost` uses `xgb.train()` with `eta = 0.05`, `max_depth = 3`,
  `min_child_weight = 5`, `subsample = 0.9`, `colsample_bytree = 0.9`,
  `nthread = 1`, and 100 rounds;
- classification predictions are clipped to `[0.01, 0.99]`.

## DML Cross-Fitting

All DML estimators use:

```text
K = number of folds
S = number of repeated random splits
```

For each split:

1. Randomly assign observations to `K` folds.
2. Fit nuisance functions on training folds.
3. Predict nuisance functions only on held-out folds.
4. Construct held-out orthogonal signals.
5. Estimate the target for that split.

Across splits:

- point estimates use the median of split estimates;
- SEs combine within-split variance and repeated-split dispersion.

Scalar DML variance helper:

- unit-level variance uses `sum(influence_i^2) / n^2`;
- cluster variance aggregates influence contributions by cluster and uses the
  finite-cluster correction `G_clusters / (G_clusters - 1)`;
- if the Jacobian is nonfinite or below `1e-10`, the split variance is
  reported as infinite;
- cluster variance falls back to unit-level variance if cluster data are
  missing or unusable.

Scalar DML inference options added on 2026-05-30:

- `dml_inference = "score"` is the default analytical score interval:
  `estimate +/- z_(1-alpha/2) * se`;
- `dml_inference = "score_multiplier"` uses the stored multiplier critical
  value for the scalar CI when the score multiplier process is well defined;
- `dml_boot` controls the number of scalar multiplier draws;
- `dml_multiplier` chooses standard normal, Rademacher/wild, or centered
  exponential weights;
- scalar DML fits store `scalar_event$influence`,
  `scalar_event$split_estimates`, `scalar_event$scalar_multiplier_reps`,
  `scalar_event$scalar_multiplier_crit`, and `scalar_event$scalar_inference`;
- nonparametric DML bootstrap and jackknife remain unsupported through
  `vartype`; the multiplier option is a score/influence multiplier path, not a
  refit bootstrap.

## `method = "dml_binary"`

### Target

Binary-G average treatment effect style target:

```text
E[mu_1(X) - mu_0(X)]
```

### Signal

```text
psi_i =
  mu1_hat(X_i) - mu0_hat(X_i)
  + 1(G_i = 1) / pi_hat(X_i) * (DeltaY_i - mu1_hat(X_i))
  - 1(G_i = 0) / (1 - pi_hat(X_i)) * (DeltaY_i - mu0_hat(X_i))
```

### Implementation

- requires `G` exactly in `{0, 1}`;
- cross-fits `mu0`, `mu1`, and `pi`;
- clips `pi_hat` to `[0.01, 0.99]`;
- stores averaged propensity scores in `result$ps`.

Fallbacks:

- if a training fold has fewer than five observations in a treatment cell, the
  corresponding outcome regression uses the global treatment-cell mean;
- if there are no covariates or only one treatment class in the training fold,
  the propensity model uses the global treatment probability.

### SE And CI

For each split:

```text
tau_s = mean(psi_i)
var_s = variance of psi_i / n
```

Final:

```text
tau = median(tau_s)
se = sqrt(median(var_s + (tau_s - tau)^2))
CI = tau +/- z_(1 - alpha/2) se
```

Clustered scalar DML SEs aggregate influence contributions by cluster.

## `method = "dml_plr"`

### Target

Partially linear model:

```text
DeltaY_i = beta G_i + h(X_i) + U_i
```

### Algorithm

Cross-fit:

```text
l_hat(X) = E[DeltaY | X]
m_hat(X) = E[G | X]
```

Residualize:

```text
U_tilde_i = DeltaY_i - l_hat(X_i)
V_tilde_i = G_i - m_hat(X_i)
```

Estimate:

```text
beta_hat = sum(V_tilde_i U_tilde_i) / sum(V_tilde_i^2)
```

If `sum(V_tilde_i^2) < 1e-10` in a split, the implementation sets that split
estimate to zero and its variance to infinity.

### SE And CI

Score:

```text
psi_i = V_tilde_i * (U_tilde_i - beta_hat V_tilde_i)
```

Jacobian:

```text
J = mean(V_tilde_i^2)
```

Split variance divides the score by `J`, then uses unit-level or cluster-level
score variance. Final aggregation uses the same DML median-plus-dispersion
formula.

## Conditional Density For Continuous-G DML

`dml_flex` and `dml_incremental` need:

```text
s(g | x) = conditional density of G given X
```

### `density_method = "residual_kde"`

Algorithm:

1. Fit `m_hat(X) = E[G | X]`.
2. Compute residuals `V = G - m_hat(X)`.
3. Estimate a univariate KDE of `V`.
4. Predict:

```text
s_hat(g | x) = f_hat_V(g - m_hat(x))
```

Interpretation:

- fast;
- assumes location-shift structure;
- can fail under strong heteroskedasticity.

### `density_method = "location_scale"`

Algorithm:

1. Fit `m_hat(X) = E[G | X]`.
2. Compute residuals.
3. Fit a regression for log residual variance.
4. Standardize residuals by estimated conditional scale.
5. KDE standardized residuals.
6. Predict:

```text
s_hat(g | x) = f_hat_Z((g - m_hat(x)) / sigma_hat(x)) / sigma_hat(x)
```

Interpretation:

- allows X-dependent scale;
- useful when `G | X` is heteroskedastic.

### `density_method = "local_kde"`

Algorithm:

- uses product Gaussian kernels over scaled covariates and a Gaussian KDE over
  `G`;
- reduces to marginal KDE when there are no covariates.

Interpretation:

- more nonparametric;
- sensitive to covariate dimension.

Shared density implementation details:

- `.safe_kde_bw()` uses Sheather-Jones bandwidth when possible and falls back
  to Silverman's rule;
- predicted densities are floored at `0.005`;
- `location_scale` uses
  `sigma_floor = max(0.05 * sd(G), 1e-4)`;
- `local_kde` uses a covariate bandwidth of
  `max(n^(-1 / (4 + p_x)), 0.25)`;
- density derivative quantities use finite differences in `G`;
- the finite-difference step is based on
  `max(1e-5 * support_span, .Machine$double.eps^0.25)` and is adjusted near
  the support boundary.

## `method = "dml_flex"`

### Target

Level curve:

```text
theta(g) = E_X[mu(g, X)]
mu(g, x) = E[DeltaY | G = g, X = x]
```

Derivative curve:

```text
theta_prime(g) = d theta(g) / dg
```

This is a fixed-G derivative of the level curve.

### Orthogonal Signal

Nuisance functions:

```text
mu(g, x) = E[DeltaY | G = g, X = x]
s(g | x) = conditional density of G given X
omega(g) = marginal density of G
m(g) = E_X[mu(g, X)]
```

Cross-fitted pseudo-outcome:

```text
xi_i =
  ((DeltaY_i - mu_hat(G_i, X_i)) / s_hat(G_i | X_i)) * omega_hat(G_i)
  + m_hat(G_i)
```

The implementation estimates `m_hat(G_i)` using held-out-fold covariates. For
linear learners it uses an exact leave-one-out formula within the held-out
fold. For nonlinear learners it predicts at `G_i` over other held-out
covariates and averages.

### Signal Mapping

After constructing `xi_i`, the package maps the signal to the `G` grid:

```text
xi_i -> theta_hat(g)
```

This second-stage grid is `eval_g`. The pseudo-outcome `xi_i` is available for
each observation, but the stored curve is evaluated only at the requested grid
points. To obtain a denser curve, fit with a denser `eval_g` or use the default
50-point grid; the plot will not create additional estimated points after the
fact.

Available `signal_map` values:

| `signal_map` | What it does | Derivative SE |
|---|---|---|
| `local_poly` | local polynomial smoother, default degree 2 | local polynomial sandwich |
| `kernel` | local-linear version of local polynomial | local polynomial sandwich |
| `spline` | B-spline linear model | delta method from `vcov()` |
| `blp_spline` | B-spline best linear projection with stored grid covariance | BLP sandwich covariance |
| `gam` | `mgcv::gam()` smooth | delta method using finite-difference design rows |

`map_degree` controls local polynomial or spline degree. `map_df` controls
spline/GAM basis size.

Mapping implementation details:

- `local_poly` and `kernel` choose a second-stage bandwidth by leave-one-out
  cross-validation on the signal-to-`G` mapping;
- `kernel` is the degree-1 version of the local-polynomial mapper;
- `spline` uses `splines::bs()` and a linear model;
- `blp_spline` uses a B-spline linear projection and stores level and
  derivative grid covariance matrices;
- `gam` uses `mgcv::gam(..., method = "REML")`;
- `map_degree` must be a single integer between 1 and 5;
- `map_df`, when supplied, must be a single integer at least 3;
- spline mapping clips the effective degree to 1 through 3 and clips the
  effective basis size to the available number of unique `G` values;
- GAM mapping uses `k = min(10, n_unique - 1)` by default, then clips it to
  feasible support;
- derivative SEs for spline/GAM maps are computed by finite-difference design
  rows;
- the mapper stores `map_info` when available, such as bandwidth, degree, or
  basis size.

### Output Fields

```r
result$eval_g
result$curve_event$theta_hat
result$curve_event$delta_hat
result$curve_event$se_theta
result$curve_event$se_delta
result$curve_event$theta_ci_lower
result$curve_event$theta_ci_upper
result$curve_event$delta_ci_lower
result$curve_event$delta_ci_upper
result$curve_event$theta_band_lower
result$curve_event$theta_band_upper
result$curve_event$delta_band_lower
result$curve_event$delta_band_upper
result$curve_event$theta_multiplier_reps
result$curve_event$delta_multiplier_reps
result$curve_event$theta_vcov
result$curve_event$delta_vcov
result$curve_event$level_derivative_vcov
result$curve_event$simultaneous_band_method
result$curve_event$interval_average
result$curve_event$mean_delta
result$curve_event$metadata
result$curve_event$diagnostics
```

### Pointwise SE And CI

For each split:

- map `xi_i` to `theta_hat_s(g)`;
- compute mapping SEs `se_theta_s(g)` and `se_delta_s(g)`.

Across splits:

```text
theta_hat(g) = median_s theta_hat_s(g)
delta_hat(g) = median_s delta_hat_s(g)
se_theta(g) = sqrt(median_s(se_theta_s(g)^2 + split_deviation_s(g)^2))
se_delta(g) = sqrt(median_s(se_delta_s(g)^2 + split_deviation_s(g)^2))
```

Pointwise CIs:

```text
theta_hat(g) +/- z_(1 - alpha/2) se_theta(g)
delta_hat(g) +/- z_(1 - alpha/2) se_delta(g)
```

### Practical Multiplier Bands

The implementation builds a practical multiplier band over the signal-mapping
step:

1. Map observed signals to `G`.
2. Compute mapping residuals.
3. Multiply residuals by standard normal weights.
4. Remap perturbed signals to the evaluation grid.
5. Use the maximum standardized deviation over the grid.

Stored method label:

```text
practical_signal_residual_multiplier
```

The practical multiplier routine stores `theta_multiplier_reps` and
`delta_multiplier_reps`, with one row per multiplier replication and one column
per `eval_g` grid point.

Caveat:

These are practical signal-mapping bands. They are not yet a final
theorem-level uniform DML band for the whole estimator.

### BLP Spline Covariance And Gaussian Max-t Bands

When `signal_map = "blp_spline"`, the second-stage mapping is a B-spline best
linear projection of the cross-fitted signal on `G`. The package stores:

```r
result$curve_event$theta_vcov
result$curve_event$delta_vcov
result$curve_event$level_derivative_vcov
```

`theta_vcov` and `delta_vcov` are grid covariance matrices after rescaling the
diagonal to match the repeated-split pointwise standard errors. Gaussian
max-t critical values are simulated from the stored covariance matrices and
used for `theta_band_*` and `delta_band_*` when available. The method label is:

```text
gaussian_max_t_blp_spline_covariance
```

This is the closest current DML-flex path to the Practical Guide/interflex BLP
logic: pseudo-outcomes are mapped through a spline basis, covariance is stored
on the evaluation grid, and simultaneous finite-grid bands use a Gaussian
max-t calculation. It is still a second-stage grid covariance approximation,
not a full refit bootstrap for all nuisance stages.

The public `boot` argument controls both the practical residual multiplier
replications and the Gaussian max-t simulation draws for `dml_flex`.

### Scalar Compatibility Output

`summary()` includes a scalar `Event` row because the `fdid` object structure
expects scalar aggregate sections.

For `dml_flex`:

```text
est$event$Estimate =
  (theta_hat(max eval_g) - theta_hat(min eval_g))
  / (max eval_g - min eval_g)

est$event$Std.Error =
  sqrt(V[max,max] + V[min,min] - 2V[max,min])
  / (max eval_g - min eval_g), when theta_vcov is available
```

This is an interval-average level slope compatibility summary. The primary
substantive output is still the curve in `curve_event`.

Current cluster caveat:

The fitted `fdid` object stores `cluster` when present, but `dml_flex` curve
diagnostics currently pass `cluster = NULL` and curve inference is not
cluster-adjusted. Cluster-aware DML variance currently applies to scalar DML
methods, not the flexible curve.

## `method = "dml_incremental"`

### Target

Observed-population average derivative:

```text
bar_delta_obs = E[partial_g mu(G_i, X_i)]
```

This differs from the fixed-G derivative curve `theta_prime(g)` from
`dml_flex`.

### Orthogonal Signal

Nuisances:

```text
mu(g, x)
s(g | x)
partial_g mu(g, x)
partial_g log s(g | x)
```

Signal:

```text
phi_i =
  partial_g mu_hat(G_i, X_i)
  - partial_g log s_hat(G_i | X_i) * (DeltaY_i - mu_hat(G_i, X_i))
```

The implementation estimates derivatives by finite differences in `G`.

### Output Fields

```r
result$est$event
result$incremental_event$signal
result$incremental_event$split_estimates
result$incremental_event$density_method
result$dml_diagnostics
```

It does not store:

```r
result$eval_g
result$curve_event
```

### SE And CI

For each split:

```text
tau_s = mean(phi_i)
var_s = score variance of phi_i around tau_s
```

Final:

```text
tau = median(tau_s)
se = sqrt(median(var_s + (tau_s - tau)^2))
CI = tau +/- z_(1 - alpha/2) se
```

Clustered SEs are supported for scalar DML by aggregating influence
contributions at the cluster level.

## Continuous-G Object Fields

All continuous-G methods return class `"fdid"` objects with the original shape:

```r
result$est$pre
result$est$event
result$est$post
result$dynamic
result$raw_means
result$tr_period
result$ref_period
result$entire_period
result$method
result$vartype
result$times
result$G
result$ps
result$call
result$target.pop
```

Kernel-specific fields:

```r
result$eval_g
result$curve_event
result$curve_dynamic
result$band_lower
result$band_upper
result$h0
result$h_adaptive
```

DML-specific fields:

```r
result$cluster
result$dml_method
result$K
result$S
result$learner
result$signal_map
result$map_degree
result$map_df
result$density_method
result$dml_inference
result$dml_boot
result$dml_multiplier
result$eval_g
result$curve_event
result$scalar_event
result$incremental_event
result$dml_metadata
result$dml_diagnostics
```

For scalar DML methods, `curve_event` is `NULL`. For `dml_binary`, `ps` stores
the split-averaged propensity scores. For `dml_incremental`,
`incremental_event` stores the scalar orthogonal signal and split estimates.
For scalar DML methods, `scalar_event` stores the scalar result, influence
values, split estimates, and multiplier diagnostics. `dml_metadata` records
the target estimand, target scale, population, inference method, inference
scope, and orthogonal-score label.

## Plotting

Continuous-G curve plotting:

```r
plot(result, type = "curve", curve = "level")
plot(result, type = "curve", curve = "derivative")
plot(result, type = "curve", curve = "both")
plot(result, type = "curve", interval = "both", show.uniform.CI = TRUE)
plot(result, type = "contrast", ref.g = 0.5)
plot(result, type = "contrast", ref.g = 0.5,
     target.g = c(0.25, 0.75))
```

Defaults:

- `kernel`: `curve = "auto"` means derivative curve;
- `dml_flex`: `curve = "auto"` means level curve.

Y-axis labels:

- kernel level: `Level DID-response mu(g)`;
- `dml_flex` level: `Level DID-response theta(g)`;
- derivative: `Derivative DID effect delta(g) = d theta(g)/dG`.

`dml_incremental` does not support `type = "curve"` because it is scalar.

Plot implementation notes:

- `plot.fdid()` supports `type = "raw"`, `"dynamic"`, `"overlap"`,
  `"curve"`, and `"contrast"`;
- raw and dynamic plots support connected lines, CI drawing, event-period
  shading, `RColorBrewer` palettes, custom labels, and axis limits;
- curve plots expose `interval = "auto"`, `"pointwise"`, `"uniform"`,
  `"both"`, and `"none"` to separate pointwise intervals from simultaneous
  bands;
- `show.uniform.CI = FALSE` suppresses method-specific simultaneous bands;
- curve plots draw pointwise intervals using the fitted object's `alpha` when
  available, with a 0.05 fallback for older or malformed objects;
- kernel bootstrap pointwise intervals use stored percentile intervals when
  available;
- kernel bootstrap simultaneous bands are displayed as finite-grid bootstrap
  quantile-envelope bands for both level and derivative curves;
- kernel robust curve plots display pointwise normal intervals only; the
  analytical cross-grid covariance is used for scalar summaries and reporting
  helpers, not as a uniform band;
- `dml_flex` curve bands display the stored curve-band method: practical
  signal-residual multiplier bands for most signal maps, or Gaussian max-t
  BLP-spline covariance bands when `signal_map = "blp_spline"` stores valid
  covariance matrices;
- `dml_flex` base-graphics legends use method-aware labels:
  `Practical band (signal residual multiplier)` for non-BLP practical
  signal-mapping bands and `Gaussian max-t band (BLP spline covariance)` for
  BLP-spline covariance bands;
- the internal curve plot-data helper stores DML target, inference, signal-map,
  density-method, band-label, band-caveat, Gaussian-critical-value, and
  covariance-availability fields so DML curve plots use the same base controls
  as kernel plots without overstating the inference;
- curve plots display the already-estimated `eval_g` rows and do not add new
  estimation points;
- curve support overlays use observed `G` values restricted to the displayed
  curve window;
- `Xdistr` supports `"histogram"`, `"hist"`, `"density"`, `"rug"`, and
  `"none"`;
- `support.panel = "embedded"` draws support in a lower strip of the curve
  panel, while `support.panel = "separate"` draws an aligned lower support
  panel for a single curve. When `curve = "both"`, support displays are
  embedded in each panel to avoid nested graphics layouts;
- `show.eval_g = TRUE` adds ticks for the stored evaluation grid;
- `diff.values` and `contrast.values` can mark selected `G` values with
  vertical reference lines;
- `type = "contrast"` draws a fixed-reference coefficient plot for selected
  level contrasts, `mu(g_j) - mu(g_ref)` for kernel and
  `theta(g_j) - theta(g_ref)` for `dml_flex`;
- by default, contrast plots include the reference row as an algebraic zero
  anchor, with `is_reference = TRUE` and a zero-width interval;
- contrast plots call `fdid_contrast()` target-by-target, so their SEs and
  intervals for non-reference rows inherit the same inference rule order:
  stored covariance, replicate curves, stored bands, then labeled pointwise
  fallback;
- `ref.g` controls the reference value and can be numeric or one of
  `"median"`, `"mean"`, `"min"`, `"max"`, or `"first"`;
- `target.g` or `contrast.values` controls displayed level targets; with
  `target.type = "quantile"`, `target.g = 5`, `target.g = "5"`, or
  `target.g = "q5"` requests five quantile-selected rows, using the first
  quantile as the default reference;
- `n.contrast` and `max.contrast` count displayed rows, including the reference
  row when `include.ref = TRUE`;
- the curve plotting path uses internal plot-data and contiguous-segment
  helpers so ribbons and band lines do not connect across missing interval
  values;
- `fdid_list()` accepts multiple `fdid` objects or a single list of them;
- `plot.fdid_list()` plots scalar event estimates and CIs only, with optional
  vertical or horizontal orientation.

## Summary Output

`summary.fdid()` reports:

- method;
- variance type;
- reference period;
- event periods;
- aggregate pre, event, and post estimates;
- dynamic estimates;
- curve summary for `kernel` and `dml_flex`;
- DML incremental signal summary for `dml_incremental`;
- DML diagnostics when available.

For `dml_flex`, summary explicitly says:

```text
Full theta(g) and theta'(g) curves in $curve_event;
scalar $est$event is interval-average theta slope.
```

For `dml_incremental`, summary explicitly says:

```text
Target: observed-population average derivative E[partial_g mu(G, X)]
```

`print.fdid()` remains intentionally brief: it prints the method, variance
type, reference period, event periods, and event-period estimate table.

## Reporting Notes-Style Estimands

This section answers how to report estimates beyond a full curve.

### Level value `theta(g0)`

Current support:

- supported by `dml_flex` as `theta(g0)` through stored `curve_event` rows;
- supported by `kernel` as `mu(g0)` through stored `curve_event` rows;
- helpers can linearly interpolate between stored `eval_g` rows, or require an
  exact grid value with `interpolate = "none"`;
- pointwise SE and CI are available;
- replicate-curve intervals are available when the fitted object stores
  multiplier or bootstrap replicate curves;
- band intervals are available when method-specific curve bands exist.

Recommended current report:

```text
curve_event row at g0:
estimate = theta_hat(g0) or mu_hat(g0)
SE = se_theta(g0) or se_mu(g0)
```

For level changes relative to a reference value, use
`fdid_contrast(result, g0 = reference_g, g1 = target_g)`.

For a compact visual summary of several such selected contrasts, use:

```text
plot(result, type = "contrast", ref.g = reference_g)
```

### Derivative value `theta_prime(g0)`

Current support:

- supported by `dml_flex` through `delta_hat(g0)`;
- supported by `kernel` through `delta_hat(g0)`;
- exposed through `fdid_derivative(result, g0)`;
- pointwise SE and CI are available;
- `dml_flex` practical multiplier derivative bands are available;
- kernel bootstrap derivative bands are available when `vartype = "bootstrap"`;
- replicate-curve intervals are available when the fitted object stores
  multiplier or bootstrap replicate curves.

Recommended current report:

```text
fdid_derivative(result, g0)
```

Important distinction:

This fixed-G derivative value is not the same as the scalar
observed-population average derivative from `dml_incremental`.

### Level comparison `theta(g1) - theta(g0)`

Current support:

- exposed through `fdid_contrast(result, g0, g1)`;
- point estimate is computed from stored curve values, with optional linear
  interpolation between stored `eval_g` rows;
- when a stored curve covariance matrix is available, the helper can report
  covariance-aware normal inference using `V11 + V22 - 2 V12` for exact grid
  endpoints, or the corresponding linear-interpolation weight formula;
- when stored replicate curves are available but covariance is unavailable,
  the helper computes contrast estimates for each replicate and reports
  replicate-based SE and CI;
- when replicate curves are unavailable but curve bands exist, the helper can
  report a conservative band-implied contrast interval;
- when neither replicate curves nor bands are available, the helper can report
  an endpoint-independent pointwise interval labeled approximate.

Efficient covariance target:

```text
Var(theta_hat(g1) - theta_hat(g0))
  = Var(theta_hat(g1)) + Var(theta_hat(g0))
    - 2 Cov(theta_hat(g1), theta_hat(g0))
```

For kernel bootstrap fits, the package stores `mu_boot_vcov` and
`delta_boot_vcov`, so exact-grid level and derivative contrasts can use the
full covariance expression. For `dml_flex` with `signal_map = "blp_spline"`,
the package stores `theta_vcov` and `delta_vcov`, so DML-flex contrasts and
derivative values can also use covariance-aware normal inference. The
endpoint-independent pointwise interval is only a fallback.

Temporary conservative interval:

```text
lower = lower_band(g1) - upper_band(g0)
upper = upper_band(g1) - lower_band(g0)
```

This is conservative and not efficient. It should be labeled as a band-implied
contrast interval, not a final contrast SE.

Recommended next implementation:

- validate replicate-contrast inference in oracle designs;
- extend curve covariance storage beyond kernel bootstrap and
  `dml_flex`/`blp_spline` if a theorem-level path is developed for the other
  signal maps;
- keep non-BLP `dml_flex` multiplier language conservative until the theory is
  final.

### Scalar average derivative

If the desired target is:

```text
E[partial_g mu(G, X)]
```

use:

```r
fdid(..., method = "dml_incremental")
```

This has a direct orthogonal score and scalar robust SE/CI now. It should not be
reported as the average of the `dml_flex` derivative curve.

## Diagnostics

DML diagnostics may include:

- learner;
- density method;
- complete-observation count;
- cluster count;
- signal quantiles;
- conditional-density quantiles;
- minimum estimated density;
- share of densities below `0.001`;
- share of densities below `0.01`;
- `partial_g mu` quantiles;
- `partial_g log s` quantiles;
- residual-correction quantiles.

Low density or extreme derivative-score tails should be treated as warnings
that derivative targets may be unstable.

## Validation And Test Record

The detailed validation and test records were consolidated in
`CONTINUOUS_G_PROJECT_LOG.md`. This section imports the estimator-relevant
record so the estimator reference can be read without opening the project log.

Primary validation scripts and saved outputs:

| Layer | Script | Saved outputs |
|---|---|---|
| Early Monte Carlo | `inst/simulate_mc.R` | script only |
| Smooth oracle-signal validation | `inst/oracle_signal_validation.R` | `inst/oracle_signal_validation_results/*.csv` |
| Stress tests | `inst/oracle_signal_stress_tests.R` | `inst/oracle_signal_stress_results/*.csv` |
| Nuisance-oracle ablation | `inst/nuisance_oracle_ablation.R` | `inst/nuisance_oracle_ablation_results/*.csv` |
| Validation figures | `inst/plot_validation_results.R` | `inst/validation_figures/*.png`, `*.csv` |

The saved validation figure files are:

```text
inst/validation_figures/oracle_rmse.png
inst/validation_figures/stress_diagnostics.png
inst/validation_figures/ablation_scalar_rmse.png
inst/validation_figures/oracle_rmse_plot_data.csv
inst/validation_figures/stress_diagnostics_plot_data.csv
inst/validation_figures/ablation_scalar_plot_data.csv
```

### Early Monte Carlo Validation

Script:

```text
inst/simulate_mc.R
```

Older generated documentation recorded early Monte Carlo checks with `n = 300`
and 100 simulations. These were preliminary regression checks, not the later
oracle-signal validation battery.

DGP 1, linear moderation:

```text
DeltaY = 2G + 0.5X + error
true delta = 2
true beta = 2
theta(g) = 2g
```

| Criterion | Result | Target |
|---|---:|---:|
| Kernel RMSE on derivative curve | 0.277 | < 0.4 |
| PLR absolute bias | 0.023 | < 0.2 |
| PLR 95% coverage | 0.940 | >= 0.85 |
| Flex RMSE on level curve | 0.074 | < 0.5 |

DGP 5, null effect:

| Criterion | Result | Target |
|---|---:|---:|
| Kernel rejection rate | 0.060 | <= 0.12 |
| PLR rejection rate | 0.040 | <= 0.12 |
| DML binary rejection rate | 0.030 | <= 0.12 |

All early targets were recorded as met. Later oracle-signal validation is the
primary evidence.

### Smooth Oracle-Signal Validation

Script:

```text
inst/oracle_signal_validation.R
```

DGP:

```text
mu(g, x) = 0.40 + 0.50 x1 - 0.25 x2 + 0.80 g + 0.25 g^2 + sin(0.70 g)
G | X ~ N(0.40 x1 - 0.30 x2, 0.70^2)
partial_g mu(g, x) = 0.80 + 0.50 g + 0.70 cos(0.70 g)
partial_g log s(g | x) = -(g - 0.40 x1 + 0.30 x2) / 0.70^2
```

Settings:

```text
Master seed: 20260503
Evaluation grid: -0.892, -0.615, -0.390, -0.190, 0.000, 0.190, 0.390, 0.615, 0.892
Estimated-nuisance learner: gam
Signal map: spline
Density method: residual_kde
```

Oracle signal sanity check:

| n | reps | level_rmse | level_coverage | derivative_rmse | derivative_coverage | scalar_bias | scalar_rmse | scalar_coverage |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 500 | 120 | 0.0536 | 0.9528 | 0.1609 | 0.9324 | -0.0009 | 0.0298 | 0.9917 |
| 1500 | 120 | 0.0342 | 0.9454 | 0.1034 | 0.9176 | 0.0024 | 0.0164 | 0.9917 |

Estimated-nuisance comparison:

| n | reps | flex_failures | incremental_failures | oracle_level_rmse | flex_level_rmse | flex_level_vs_oracle_rmse | flex_level_coverage | oracle_derivative_rmse | flex_derivative_rmse | flex_derivative_vs_oracle_rmse | flex_derivative_coverage | oracle_scalar_bias | oracle_scalar_rmse | oracle_scalar_coverage | incremental_scalar_bias | incremental_scalar_rmse | incremental_scalar_coverage | incremental_vs_oracle_bias | incremental_vs_oracle_rmse |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 400 | 40 | 0 | 0 | 0.0588 | 0.0640 | 0.0352 | 0.9722 | 0.1925 | 0.1837 | 0.1318 | 0.9583 | -0.0061 | 0.0300 | 1.0000 | -0.0061 | 0.0326 | 0.9750 | 0.0001 | 0.0146 |
| 800 | 40 | 0 | 0 | 0.0404 | 0.0401 | 0.0133 | 0.9694 | 0.1184 | 0.1131 | 0.0508 | 0.9472 | 0.0007 | 0.0239 | 0.9750 | 0.0005 | 0.0254 | 0.9500 | -0.0001 | 0.0090 |

Interpretation:

- oracle signal RMSE falls as `n` increases;
- estimated-nuisance RMSE is generally larger than oracle RMSE;
- scalar `dml_incremental` targets `E[partial_g mu(G, X)]`;
- derivative-curve comparisons belong to `dml_flex`, not
  `dml_incremental`.

### Stress Tests

Script:

```text
inst/oracle_signal_stress_tests.R
```

Settings:

```text
Master seed: 20260504
Oracle reps per scenario/n: 80
Estimated reps per scenario/n/density: 12
Estimated learner: gam
Signal map: spline
Map df: 7
Density methods: residual_kde, location_scale
```

Scenarios:

- `heteroskedastic_interactions`: normal `G | X` with X-dependent scale and
  smooth G-by-X interactions in `mu(g, x)`.
- `mixture_kink`: two-component normal mixture `G | X` with X-dependent
  weights and kinks in `mu(g, x)`.
- `heavy_tail_interactions`: scaled t(3) residuals in `G | X` with nonlinear
  trigonometric interactions in `mu(g, x)`.

Oracle signal results:

| scenario | n | reps | level_rmse | level_coverage | derivative_rmse | derivative_coverage | scalar_bias | scalar_rmse | scalar_coverage |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| heavy_tail_interactions | 600 | 80 | 0.0400 | 0.9222 | 0.0972 | 0.9028 | 0.0009 | 0.0240 | 0.9875 |
| heavy_tail_interactions | 1200 | 80 | 0.0306 | 0.9056 | 0.0873 | 0.7681 | 0.0010 | 0.0171 | 0.9625 |
| heteroskedastic_interactions | 600 | 80 | 0.0484 | 0.9389 | 0.1578 | 0.9458 | 0.0020 | 0.0286 | 1.0000 |
| heteroskedastic_interactions | 1200 | 80 | 0.0372 | 0.9444 | 0.1239 | 0.9097 | -0.0017 | 0.0200 | 0.9875 |
| mixture_kink | 600 | 80 | 0.0418 | 0.9278 | 0.1865 | 0.8278 | 0.0001 | 0.0268 | 0.9875 |
| mixture_kink | 1200 | 80 | 0.0327 | 0.9319 | 0.1686 | 0.7597 | -0.0007 | 0.0170 | 1.0000 |

Estimated-nuisance results:

| scenario | density_method | n | reps | flex_failures | incremental_failures | oracle_level_rmse | flex_level_rmse | flex_level_vs_oracle_rmse | flex_level_coverage | oracle_derivative_rmse | flex_derivative_rmse | flex_derivative_vs_oracle_rmse | flex_derivative_coverage | incremental_scalar_bias | incremental_scalar_rmse | incremental_scalar_coverage | incremental_vs_oracle_rmse |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| heavy_tail_interactions | location_scale | 500 | 12 | 0 | 0 | 0.0415 | 0.0614 | 0.0361 | 0.9537 | 0.0903 | 0.1489 | 0.1133 | 0.9352 | 0.0157 | 0.0307 | 1.0000 | 0.0210 |
| heavy_tail_interactions | location_scale | 1000 | 12 | 0 | 0 | 0.0289 | 0.0350 | 0.0211 | 0.9537 | 0.0759 | 0.0855 | 0.0655 | 0.9630 | -0.0018 | 0.0320 | 0.9167 | 0.0170 |
| heavy_tail_interactions | residual_kde | 500 | 12 | 0 | 0 | 0.0420 | 0.0446 | 0.0377 | 0.9722 | 0.0864 | 0.1040 | 0.1063 | 0.9630 | 0.0099 | 0.0282 | 1.0000 | 0.0210 |
| heavy_tail_interactions | residual_kde | 1000 | 12 | 0 | 0 | 0.0262 | 0.0304 | 0.0204 | 0.9722 | 0.0827 | 0.0994 | 0.0570 | 0.8796 | -0.0033 | 0.0233 | 1.0000 | 0.0190 |
| heteroskedastic_interactions | location_scale | 500 | 12 | 0 | 0 | 0.1112 | 0.1658 | 0.1088 | 0.9444 | 0.3788 | 0.5509 | 0.4129 | 0.9537 | 0.0241 | 0.0601 | 1.0000 | 0.0494 |
| heteroskedastic_interactions | location_scale | 1000 | 12 | 0 | 0 | 0.0363 | 0.0765 | 0.0625 | 0.9074 | 0.1124 | 0.2307 | 0.2054 | 0.8241 | 0.0332 | 0.0443 | 0.9167 | 0.0343 |
| heteroskedastic_interactions | residual_kde | 500 | 12 | 0 | 0 | 0.0432 | 0.1103 | 0.1027 | 0.9537 | 0.1132 | 0.4211 | 0.4112 | 0.8981 | 0.1439 | 0.1524 | 0.0833 | 0.1528 |
| heteroskedastic_interactions | residual_kde | 1000 | 12 | 0 | 0 | 0.0659 | 0.1061 | 0.1086 | 0.8148 | 0.1784 | 0.3229 | 0.3555 | 0.6667 | 0.1461 | 0.1484 | 0.0000 | 0.1331 |
| mixture_kink | location_scale | 500 | 12 | 0 | 0 | 0.0427 | 0.0783 | 0.0690 | 1.0000 | 0.2008 | 0.3263 | 0.2638 | 0.9074 | 0.0042 | 0.0351 | 1.0000 | 0.0422 |
| mixture_kink | location_scale | 1000 | 12 | 0 | 0 | 0.0384 | 0.0681 | 0.0501 | 0.9444 | 0.1718 | 0.2825 | 0.2031 | 0.8056 | -0.0005 | 0.0334 | 1.0000 | 0.0205 |
| mixture_kink | residual_kde | 500 | 12 | 0 | 0 | 0.0548 | 0.0640 | 0.0444 | 0.8611 | 0.2081 | 0.2574 | 0.1743 | 0.7778 | -0.0715 | 0.0790 | 0.6667 | 0.0823 |
| mixture_kink | residual_kde | 1000 | 12 | 0 | 0 | 0.0352 | 0.0480 | 0.0362 | 0.8519 | 0.1599 | 0.2028 | 0.1322 | 0.7593 | -0.0689 | 0.0759 | 0.3333 | 0.0681 |

Stress-test interpretation:

- `flex_*` columns evaluate `dml_flex` level and derivative curves.
- `incremental_*` columns evaluate scalar `dml_incremental`.
- `*_vs_oracle_*` columns compare estimated-nuisance fits to oracle-signal
  fits on the same DGP.
- density choice matters strongly for scalar incremental targets in
  heteroskedastic and mixture settings.

### Nuisance-Oracle Ablation

Script:

```text
inst/nuisance_oracle_ablation.R
```

Settings:

```text
Master seed: 20260505
Replications per scenario/n: 24
Sample sizes: 600, 1000
Learner for estimated nuisances: gam
Signal map: spline
Map df: 7
Estimated density methods: residual_kde, location_scale
```

Ablation layers:

| Layer | Meaning |
|---|---|
| `oracle_all` | true outcome and true density nuisances |
| `estimated_outcome_oracle_density` | estimated outcome, true density |
| `oracle_outcome_estimated_density` | true outcome, estimated density |
| `estimated_all` | estimated outcome and estimated density |

Representative complete-summary rows:

| scenario | n | layer | density_method | reps | failures | level_rmse | level_coverage | derivative_rmse | derivative_coverage | scalar_bias | scalar_rmse | scalar_coverage | level_vs_oracle_rmse | derivative_vs_oracle_rmse | scalar_vs_oracle_rmse |
|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| heavy_tail_interactions | 600 | oracle_all | oracle | 24 | 0 | 0.0359 | 0.9444 | 0.0897 | 0.9167 | -0.0036 | 0.0218 | 1.0000 | 0.0000 | 0.0000 | 0.0000 |
| heavy_tail_interactions | 600 | estimated_all | residual_kde | 24 | 0 | 0.0550 | 0.9398 | 0.1445 | 0.8750 | -0.0014 | 0.0230 | 1.0000 | 0.0407 | 0.1164 | 0.0221 |
| heteroskedastic_interactions | 600 | oracle_all | oracle | 24 | 0 | 0.0429 | 0.9630 | 0.1469 | 0.9352 | 0.0065 | 0.0314 | 1.0000 | 0.0000 | 0.0000 | 0.0000 |
| heteroskedastic_interactions | 600 | estimated_all | residual_kde | 24 | 0 | 0.1247 | 0.8981 | 0.3976 | 0.8704 | 0.1308 | 0.1368 | 0.0417 | 0.1142 | 0.3731 | 0.1278 |
| heteroskedastic_interactions | 600 | estimated_all | location_scale | 24 | 0 | 0.0946 | 0.8889 | 0.3032 | 0.8889 | 0.0350 | 0.0476 | 0.9583 | 0.0779 | 0.2632 | 0.0379 |
| mixture_kink | 600 | oracle_all | oracle | 24 | 0 | 0.0426 | 0.9074 | 0.1879 | 0.8241 | -0.0002 | 0.0260 | 1.0000 | 0.0000 | 0.0000 | 0.0000 |
| mixture_kink | 600 | estimated_all | residual_kde | 24 | 0 | 0.0612 | 0.9028 | 0.2506 | 0.7546 | -0.0745 | 0.0826 | 0.5000 | 0.0456 | 0.1606 | 0.0808 |
| mixture_kink | 600 | estimated_all | location_scale | 24 | 0 | 0.0946 | 0.9213 | 0.3636 | 0.8611 | -0.0102 | 0.0316 | 1.0000 | 0.0753 | 0.2761 | 0.0215 |

Full ablation records are saved in:

```text
inst/nuisance_oracle_ablation_results/nuisance_oracle_ablation_summary.csv
inst/nuisance_oracle_ablation_results/nuisance_oracle_ablation_replications.csv
```

Ablation interpretation:

- `*_vs_oracle_*` columns compare each layer with `oracle_all` on the same
  replication;
- large `oracle_outcome_estimated_density` errors isolate density and
  density-derivative problems;
- large `estimated_outcome_oracle_density` errors isolate outcome-regression
  and outcome-derivative problems;
- `dml_flex` curve diagnostics are in level and derivative columns;
- `dml_incremental` diagnostics are in scalar columns.

### Package Tests And Build Verification

Recorded verification:

- `testthat::test_file("tests/testthat/test-continuous.R")`: 139 passed,
  3 skipped on CRAN;
- `devtools::test(filter = "continuous")`: 128 passed;
- `devtools::test()`: 199 passed, with one pre-existing warning about
  `target.pop` for `did`;
- RStudio bundled Quarto render of `tutorial/03-continuous.Rmd` produced
  `tutorial/_book/03-continuous.html`;
- `R CMD build --no-build-vignettes fdid-main`: succeeded;
- `R CMD check --no-manual --no-build-vignettes fdid_1.0.2.tar.gz`:
  2 expected no-vignette warnings, 0 errors, 0 notes;
- 2026-05-19 documentation-audit packaging check:
  `R CMD build --no-build-vignettes fdid-main` succeeded, and a tar listing
  confirmed that `historical_mds_2026-05-19/` and validation result folders
  were excluded while active project markdowns remained included.
- 2026-05-30 DML inference implementation check:
  `devtools::test("fdid-main", filter = "continuous")`: 233 passed, 0
  warnings, 0 failures.

### Validation Conclusions

- Oracle signals recover the intended targets.
- Estimated nuisance versions behave close to oracle baselines in smooth DGPs.
- `dml_flex` level curves are more stable than derivative curves.
- `dml_incremental` is sensitive to conditional density and derivative
  nuisance quality.
- Density method choice matters, especially under heteroskedastic or mixture
  treatment density.
- Derivative targets require more caution than level targets.

## Practical Method Choice

| Need | Start with | Primary display | Main caveat |
|---|---|---|---|
| quick scalar slope | `dml_plr` | scalar table | assumes partially linear effect |
| low-dimensional curve | `kernel` | level and derivative curves | weak in high-dimensional covariates |
| flexible adjusted level curve | `dml_flex` | level curve | density and mapping matter |
| fixed-G derivative curve | `dml_flex` or `kernel` | derivative curve | derivative is fragile |
| scalar observed average derivative | `dml_incremental` | scalar table | density derivative is fragile |
| binary companion | `dml_binary` | scalar table | target differs from continuous-G targets |

## Current Caveats

1. `dml_flex` bands are second-stage curve bands. For `signal_map =
   "blp_spline"`, they use Gaussian max-t simulation from the stored BLP-spline
   grid covariance; for other maps, they remain practical signal-mapping
   multiplier bands, not final theorem-level uniform DML bands.
2. Level-comparison and fixed-G derivative reporting are first-class through
   `fdid_contrast()` and `fdid_derivative()`. Kernel bootstrap fits now support
   replicate-based and curve-covariance-based grid contrasts; kernel robust fits
   support analytical stacked-sandwich curve covariance; `dml_flex` with
   `signal_map = "blp_spline"` supports curve-covariance-based contrasts.
   Theorem-level efficient contrast inference remains open for non-BLP
   `dml_flex` paths.
3. Derivative targets rely on finite-difference nuisance derivatives.
4. Conditional density estimation remains a first-order practical issue.
5. DML refit bootstrap and jackknife are unsupported through `vartype`.
   Scalar DML score multiplier intervals are available through
   `dml_inference = "score_multiplier"`.
6. External package benchmark checks remain future validation work.
7. `boot` controls the `dml_flex` practical residual-multiplier draws and
   Gaussian max-t simulation draws.
8. Kernel bootstrap resampling is cluster-aware when a cluster column is
   available. Kernel robust stacked sandwiches are also cluster-summed when a
   cluster column is available. `dml_flex` curve inference is not currently
   cluster-adjusted.
9. Continuous-G raw means are median-split descriptive displays, not
    continuous-G estimators.
