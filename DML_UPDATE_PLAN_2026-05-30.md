# DML Update Plan

Date: 2026-05-30

Status: implementation pass started and partially completed on 2026-05-30.
This document remains the roadmap; implemented items are recorded below and in
`CONTINUOUS_G_PROJECT_LOG.md`.

Workflow: StatsClaw docs-only planning artifact now; future implementation
should use the full StatsClaw code-change workflow because it will touch
estimator code, tests, documentation, and validation evidence.

Implementation note, 2026-05-30: the first DML inference implementation pass
used the StatsClaw staged workflow in this workspace and updated estimator
code, reporting helpers, tests, FDID Notes, the estimator reference, project
log, TODO tracking, and Rd reference files. `NEWS.md`, `HANDOFF.md`, and
tutorial rewrites stayed out of scope, consistent with earlier planning.

## Not Completed In The 2026-05-30 Implementation Pass

The following items remain research/validation work after the first DML
inference implementation pass:

1. DML refit bootstrap:
   - nonparametric refit bootstrap for DML nuisance stages and final targets is
     not implemented;
   - `vartype = "bootstrap"` remains unsupported for DML methods;
   - scalar score multiplier inference is implemented, but it is not a refit
     bootstrap.
2. Non-BLP `dml_flex` theorem-level uniform bands:
   - `signal_map = "blp_spline"` now has stored grid covariance and Gaussian
     max-t finite-grid bands;
   - other `dml_flex` signal maps still use practical residual-multiplier
     signal-mapping bands;
   - these non-BLP bands should not be described as theorem-level uniform DML
     bands for the full estimator.
3. Systematic simulation and external-package benchmarks:
   - additional finite-sample simulations are still needed before recommending
     scalar score multiplier inference as the default;
   - external benchmarks against packages such as DoubleML remain future
     validation work where comparable estimands exist;
   - direct external benchmarks for `dml_flex` and `dml_incremental` may remain
     limited because standard package equivalents are not one-to-one.

## Purpose

This plan records a complete update path for the DML estimators in `fdid`.
It follows the audit of:

- `FDID Notes/main.tex`;
- `fdid-main/R/dml.R`;
- the local `interflex-master` DML implementation;
- the local Practical Guide chapter,
  `A Practical Guide to Estimating Conditional Marginal Effects- Modern Approaches/sections/4_dml.tex`;
- external continuous-treatment DML, derivative-effect, and continuous-DID
  references.

The central recommendation is not to replace the current FDID DML estimands.
The current signal split is conceptually coherent:

- `dml_binary` targets a binary-G AIPW/IRM contrast;
- `dml_plr` targets a scalar partially linear continuous-G slope;
- `dml_flex` targets the fixed-G level curve `theta(g)` and derivative curve
  `theta'(g)`;
- `dml_incremental` targets the observed-population scalar average derivative
  `E[partial_g mu(G, X)]`.

The main update should instead build a more systematic inference layer,
stabilize derivative-sensitive paths, improve diagnostics, and align user-facing
documentation with what the implementation can and cannot prove.

## Current Baseline

### Current DML Methods

| Method | Current target | Output shape | Primary code |
|---|---|---|---|
| `dml_binary` | `E[mu_1(X) - mu_0(X)]` for binary `G` | scalar | `.est_dml_binary()` |
| `dml_plr` | best partially linear scalar slope | scalar | `.est_dml_plr()` |
| `dml_flex` | `theta(g) = E_X[mu(g, X)]` and `theta'(g)` | curve | `.est_dml_flex()` |
| `dml_incremental` | `E[partial_g mu(G, X)]` | scalar | `.est_dml_incremental_scalar()` |

### Current DML Options

Shared DML options:

- `learner`: nuisance learner;
- `K`: cross-fitting folds;
- `S`: repeated random sample splits;
- `alpha`: confidence level control.

Curve-specific options:

- `eval_g`: estimator grid for `dml_flex`;
- `trim`: trimmed support used to build default `eval_g`;
- `boot`: currently used by `dml_flex` practical multiplier bands;
- `signal_map`: `local_poly`, `kernel`, `spline`, or `gam`;
- `map_degree`: local polynomial or spline degree;
- `map_df`: spline or GAM basis size;
- `density_method`: `residual_kde`, `location_scale`, or `local_kde`.

### Current Inference

Scalar DML methods currently use score-based analytical standard errors:

```text
tau_s = mean(score_i)
var_s = split-specific variance of the influence contribution
tau = median_s(tau_s)
se = sqrt(median_s(var_s + (tau_s - tau)^2))
```

For cluster data, scalar DML aggregates influence contributions by cluster and
uses a finite-cluster correction.

`dml_flex` currently uses:

- pointwise second-stage signal-mapping standard errors;
- repeated-split dispersion folded into pointwise SEs;
- practical multiplier bands over signal-mapping residuals.

The `dml_flex` practical bands are useful finite-grid uncertainty diagnostics,
but they should not be described as final theorem-level DML uniform bands for
the full estimator.

## Main Design Principle

Keep the estimand split sharp.

Do not collapse these three continuous-G quantities:

```text
1. PLR scalar:
   beta in DeltaY = beta G + h(X) + U

2. Fixed-G level and derivative curves:
   theta(g) = E_X[mu(g, X)]
   theta'(g)

3. Observed-population scalar average derivative:
   E[partial_g mu(G, X)]
```

The first is a partially linear projection target. The second is a dose-response
curve and its derivative. The third is a scalar average derivative evaluated at
observed `G_i`. They are numerically and conceptually different except under
special restrictions.

## StatsClaw Workflow Classification

The future implementation should use a full code-change workflow, not a
docs-only shortcut.

Suggested workflow:

```text
leader -> planner -> builder -> tester -> scriber -> reviewer
```

If simulation evidence is added in the same run:

```text
leader -> planner -> [builder || simulator] -> tester -> scriber -> reviewer
```

Pipeline separation:

- builder receives the implementation spec only;
- tester receives an independent test spec only;
- simulator receives a simulation spec only;
- scriber updates documentation after code and tests are stable;
- reviewer checks estimator theory, code behavior, test evidence, and docs.

## Phase 0: Pre-Implementation Lockdown

### Goal

Before code changes, lock down the exact DML targets, user-facing language, and
testable acceptance criteria.

### Actions

1. Read and freeze the current DML theory sections in `FDID Notes/main.tex`.
2. Confirm the current source map:
   - `R/dml.R`;
   - `R/fdid.R`;
   - `R/summary.R`;
   - `R/plot.R`;
   - `R/reporting.R`;
   - `tests/testthat/test-continuous.R`;
   - `CONTINUOUS_G_ESTIMATOR_REFERENCE.md`;
   - `CONTINUOUS_G_PROJECT_LOG.md`;
   - `TODO.md`;
   - `man/fdid.Rd`.
3. Decide whether this implementation run will include simulation evidence.
4. Decide whether to update tutorials now or defer them. Current recommendation:
   defer tutorial rewriting until the core inference API stabilizes.
5. Avoid updating `NEWS.md` and `HANDOFF.md` unless this becomes a release or
   handoff task.

### Acceptance Criteria

- There is a written spec stating which estimands each method targets.
- The spec explicitly says that `dml_incremental` is not a curve method.
- The spec explicitly says that `dml_flex` scalar `est$event` is a compatibility
  summary, not the primary estimand.
- The spec distinguishes analytical, multiplier, bootstrap, pointwise, and
  uniform inference.

## Phase 1: Estimand Registry And Object Metadata

### Goal

Make every fitted DML object self-describing enough that summary, plot,
reporting helpers, and documentation do not infer targets from method names
alone.

### Recommended Code Updates

Add internal metadata fields to DML output objects.

Suggested fields:

```r
target_family
target_estimand
target_scale
target_population
inference_type
inference_scope
inference_method
uniform_band_method
orthogonal_score
```

Example metadata:

```r
# dml_plr
target_family = "continuous_G_DML"
target_estimand = "partially_linear_slope"
target_scale = "scalar"
target_population = "all_complete_observations"
inference_type = "analytical_score"
inference_scope = "pointwise_scalar"
orthogonal_score = "partialling_out"

# dml_flex
target_estimand = "fixed_G_level_curve"
target_scale = "curve"
inference_type = "second_stage_mapping"
inference_scope = "pointwise_curve_and_practical_band"
uniform_band_method = "signal_mapping_residual_multiplier"

# dml_incremental
target_estimand = "observed_population_average_derivative"
target_scale = "scalar"
inference_type = "analytical_score"
orthogonal_score = "average_derivative_density_score"
```

### Why This Matters

Current object fields are already usable, but the method name carries too much
meaning. Explicit metadata will reduce future ambiguity in:

- `summary.fdid()`;
- `plot.fdid()`;
- `fdid_contrast()`;
- `fdid_derivative()`;
- tutorials;
- external validation scripts.

### Files

Likely touched files:

- `R/dml.R`;
- `R/summary.R`;
- `R/plot.R`;
- `R/reporting.R`;
- `CONTINUOUS_G_ESTIMATOR_REFERENCE.md`;
- `TODO.md`;
- tests for object fields.

### Tests

Add tests that assert:

- every DML object has metadata fields;
- `dml_flex` reports curve-first metadata;
- `dml_incremental` reports scalar observed-average-derivative metadata;
- plotting and reporting helpers read metadata or remain consistent with it.

## Phase 2: Scalar DML Score Multiplier Inference

### Goal

Add multiplier inference for scalar DML methods:

- `dml_binary`;
- `dml_plr`;
- `dml_incremental`.

This should complement the current analytical score SEs. It should not replace
them as the default until validated.

### Current State

Scalar DML currently supports `vartype = "robust"` only. Non-robust requests are
rejected early in `fdid()`.

The current scalar SE is analytical:

```text
score_i -> influence_i -> var(influence_i) / n
```

For repeated splits:

```text
se = sqrt(median(var_s + split_dispersion_s))
```

### Proposed API

Keep `vartype = "robust"` as the public default for DML.

Add a DML-specific option rather than overloading `vartype` immediately:

```r
dml_inference = c("score", "score_multiplier")
dml_boot = 500L
dml_multiplier = c("normal", "wild", "exponential")
```

Alternative:

```r
vartype = "robust"
dml_band = c("none", "multiplier")
```

Recommendation: use `dml_inference` because it is scalar and curve neutral.

### Algorithm: Non-Clustered Scalar Multiplier

For a split `s`, obtain estimated influence contributions:

```text
IF_i = (score_i - center) / J
```

For each multiplier draw `b`:

```text
T_b = (1 / sqrt(n)) * sum_i xi_i_b * IF_i / sd(IF_i)
```

where `xi_i_b` can be standard normal by default.

Use the empirical quantile of `|T_b|` for scalar robust-multiplier intervals.
For one scalar target this will often be close to a normal interval; the value
is greater for joint scalar families if later extended.

### Algorithm: Clustered Scalar Multiplier

For cluster `c`, aggregate:

```text
IF_c = sum_{i in c} IF_i
```

Then draw cluster-level multipliers:

```text
T_b = scale_factor * sum_c xi_c_b * IF_c / estimated_cluster_sd
```

Use the same finite-cluster correction convention as current scalar DML where
possible. Store the number of clusters and warn when the cluster count is small.

### Split Aggregation

Two implementation choices:

1. Split-wise multiplier:
   - compute critical values per split;
   - combine intervals with median point estimate and split dispersion.

2. Pooled influence approximation:
   - store averaged or median-centered influence contributions across splits;
   - run multiplier once on the pooled influence object.

Recommendation: implement split-wise first because it follows current
`tau_s`/`var_s` architecture. Add pooled influence only after validation.

### Output Fields

For scalar DML event objects, store:

```r
result$scalar_inference
result$scalar_multiplier_reps
result$scalar_multiplier_crit
result$scalar_multiplier_method
result$scalar_influence
```

For `dml_incremental`, also store:

```r
result$incremental_event$influence
result$incremental_event$multiplier_reps
result$incremental_event$multiplier_crit
```

### Reference

Use DoubleML's multiplier bootstrap as the closest software reference. The
DoubleML documentation describes multiplier weights and joint confidence
intervals based on bootstrap t-statistics:

- https://docs.doubleml.org/dev/guide/se_confint.html

### Acceptance Criteria

- Analytical score SEs remain unchanged by default.
- Optional scalar multiplier intervals are available for all scalar DML methods.
- Clustered scalar multiplier inference uses cluster-level multipliers.
- All output labels distinguish analytical-score CI from multiplier CI.
- Tests verify reproducibility under `set.seed()`.
- Tests verify fallback behavior when `dml_boot` is too small or influence
  variance is degenerate.

## Phase 3: Theorem-Oriented `dml_flex` BLP/Sieve Inference

### Goal

Add a second `dml_flex` inference path closer to the Practical Guide and
`interflex` DML approach:

```text
orthogonal signal -> B-spline/sieve projection -> HC/sandwich covariance
                 -> Gaussian max-t uniform band
```

This should coexist with the current local-polynomial/spline/GAM mapping path.

### Current State

`dml_flex` currently constructs cross-fitted pseudo-outcomes:

```text
xi_i =
  ((DeltaY_i - mu_hat(G_i, X_i)) / s_hat(G_i | X_i)) * omega_hat(G_i)
  + m_hat(G_i)
```

Then it maps:

```text
xi_i -> theta(g)
```

using:

- local polynomial;
- kernel/local linear;
- spline lm;
- GAM.

Pointwise SEs come from the second-stage map. Practical uniform bands come from
residual multipliers around the map.

### Proposed New Mapping/Inferences

Add:

```r
signal_map = "blp_spline"
```

or:

```r
dml_curve_inference = c("mapping", "blp_gaussian")
```

Recommendation:

- keep `signal_map` for the fitted curve shape;
- add `curve_inference = c("mapping", "blp_gaussian", "practical_multiplier")`;
- allow `signal_map = "spline"` to use either ordinary mapping SEs or the new
  BLP/Gaussian covariance path.

### Algorithm

Let `B(G_i)` be a spline basis.

Fit:

```text
xi_i = B(G_i)' beta + residual_i
```

Then for grid `g`:

```text
theta_hat(g) = B(g)' beta_hat
```

Use HC0/HC1 sandwich covariance:

```text
V_beta = (B'B)^(-1) B' diag(residual_i^2) B (B'B)^(-1)
V_grid = B_grid V_beta B_grid'
se_theta(g) = sqrt(diag(V_grid))
```

For derivative:

```text
D_B(g) = d B(g) / dg
delta_hat(g) = D_B(g)' beta_hat
V_delta_grid = D_B_grid V_beta D_B_grid'
se_delta(g) = sqrt(diag(V_delta_grid))
```

If analytical spline derivatives are inconvenient initially, use stable finite
differences of `B(g)` with stored step sizes. Later replace with analytical
derivative basis for B-splines.

### Uniform Bands

Use a Gaussian max-t critical value:

```text
Z ~ N(0, V_grid)
T = max_g |Z_g| / se_theta(g)
c_{1-alpha} = quantile(T, 1-alpha)
band(g) = theta_hat(g) +/- c_{1-alpha} se_theta(g)
```

This mirrors the current local `interflex-master/R/uniform.R` function
`calculate_delta_uniformCI()`, which simulates from a multivariate normal with
the estimated grid covariance matrix.

### Contrast Inference

Once `V_grid` is stored, implement exact grid-level contrast variance:

```text
Var(theta(g2) - theta(g1)) =
  V[g2, g2] + V[g1, g1] - 2 V[g2, g1]
```

This directly answers the earlier concern about endpoint covariance. Kernel
bootstrap already stores covariance matrices; `dml_flex` should do the same for
the BLP/sieve path.

### Output Fields

Store:

```r
curve_event$theta_vcov
curve_event$delta_vcov
curve_event$theta_delta_vcov
curve_event$theta_blp_basis
curve_event$theta_blp_df
curve_event$theta_blp_degree
curve_event$theta_uniform_crit
curve_event$delta_uniform_crit
curve_event$uniform_band_method
```

If cross-covariance between level and derivative is available:

```r
curve_event$level_derivative_vcov
```

### Reference

Primary theory reference:

- Semenova and Chernozhukov (2021), "Debiased machine learning of conditional
  average treatment effects and other causal functions",
  The Econometrics Journal.
  https://academic.oup.com/ectj/article-pdf/24/2/264/46748086/grf.pdf

Practical software reference:

- `interflex-master/R/DML.R`, especially `.compute_cate_blp()`;
- `interflex-master/R/uniform.R`, especially `calculate_delta_uniformCI()`.

Guide reference:

- local Practical Guide chapter `sections/4_dml.tex`, especially the DML
  inference subsection that uses a projected influence function and Gaussian
  multiplier process.

### Acceptance Criteria

- `dml_flex` can still run with current local-polynomial default.
- A spline/BLP inference path can return `theta_vcov`.
- `fdid_contrast()` uses `theta_vcov` for grid-level `dml_flex` contrasts when
  available.
- Uniform bands are labeled as Gaussian max-t / BLP covariance bands, not as
  bootstrap bands.
- Practical residual-multiplier bands remain labeled as practical bands.
- Tests compare `Var(g2-g1)` with direct covariance formula.

## Phase 4: `dml_flex` Practical Band Cleanup

### Goal

Keep the current practical multiplier band, but make it more transparent and
less confusing.

### Current State

The current `dml_flex` multiplier process:

1. maps observed signals to observed `G`;
2. computes mapping residuals;
3. multiplies residuals by normal weights;
4. remaps perturbed signals to `eval_g`;
5. records maximum standardized deviations.

This is useful, but it is not the same as the theorem-level projected influence
function process in the Practical Guide.

### Recommended Changes

1. Rename internal method label from:

```text
multiplier over cross-fitted signal-mapping residuals
```

to:

```text
practical residual multiplier over second-stage signal mapping
```

2. Add fields:

```r
curve_event$band_level = "practical"
curve_event$band_scope = "second_stage_mapping"
curve_event$band_warning
```

3. Make `summary()` print:

```text
Simultaneous bands: practical residual multiplier over second-stage mapping
```

not just "multiplier".

4. Reconcile public `boot = 500` with internal cap at 100:
   - either expose `dml_flex_max_multiplier = 100`;
   - or remove the cap and use progress-friendly tests;
   - or document the cap explicitly in `?fdid`.

Recommendation: expose the cap as an internal guardrail for now but document it
clearly in reference docs and summary output.

### Acceptance Criteria

- No user-facing output implies that current practical bands are final
  theorem-level DML uniform bands.
- `boot` cap is documented wherever `boot` is documented.
- Tests assert the stored method label.

## Phase 5: Cluster-Aware Curve Inference

### Goal

Extend cluster-aware inference beyond scalar DML to curve objects where
possible.

### Current State

Scalar DML can use cluster-aware score variance. `dml_flex` stores cluster
labels in the outer `fdid` object, but curve diagnostics currently pass
`cluster = NULL`; curve inference is not cluster-adjusted.

### Recommended Implementation Path

Start with the BLP/sieve path because cluster sandwich covariance is natural
there.

For basis regressors `B_i` and residuals `u_i`, compute cluster meat:

```text
M_cluster = sum_c (sum_{i in c} B_i u_i) (sum_{i in c} B_i u_i)'
V_beta_cluster = (B'B)^(-1) M_cluster (B'B)^(-1)
```

Apply finite-cluster correction:

```text
G_clusters / (G_clusters - 1)
```

optionally with an `n` and `p` adjustment if adopted elsewhere in the package.

Then propagate to grid:

```text
V_grid_cluster = B_grid V_beta_cluster B_grid'
```

### API

Do not add a new public option at first. If `cluster_label` was supplied to
`fdid_prepare()`, use cluster-aware curve covariance for the BLP/sieve path and
record it.

Store:

```r
curve_event$curve_vcov_type = "cluster"
curve_event$n_clusters
```

Fallback:

- if fewer than 2 clusters, use non-cluster covariance and warn;
- if cluster labels are missing after complete-case filtering, use non-cluster
  covariance and record fallback.

### Acceptance Criteria

- Clustered `dml_flex` BLP/sieve path returns `n_clusters`.
- Tests compare clustered and non-clustered paths on a toy dataset.
- Summary reports whether curve covariance is cluster-adjusted.

## Phase 6: Conditional Density And Derivative Stability

### Goal

Reduce silent instability in `dml_flex` derivatives and `dml_incremental`.

### Current State

Continuous-G DML uses estimated conditional densities:

```text
s(g | x)
```

and for `dml_incremental` also uses:

```text
partial_g mu(g, x)
partial_g log s(g | x)
```

Derivatives are currently finite-difference approximations. This is practical
but fragile, especially for nonsmooth learners and boundary points.

### Recommended Changes

1. Store finite-difference step sizes.

Suggested fields:

```r
dml_diagnostics$finite_diff_step_mu
dml_diagnostics$finite_diff_step_density
dml_diagnostics$finite_diff_support_span
```

2. Add derivative sensitivity diagnostics.

For derivative-sensitive methods, compute derivative estimates at multiple step
scales:

```text
eps
2 eps
4 eps
```

Store max relative change:

```r
dml_diagnostics$d_mu_step_sensitivity
dml_diagnostics$d_log_density_step_sensitivity
```

3. Warn when sensitivity is high.

Example warning:

```text
Derivative-sensitive DML diagnostics indicate high finite-difference
sensitivity. Consider smoother learners, a different density_method, or
reporting level curves/interval averages instead of point derivatives.
```

4. Add learner guidance.

For derivative-sensitive targets, docs should recommend smoother learners:

- `linear`;
- `gam`;
- possibly `glmnet` with smooth basis-expanded covariates.

Warn that tree learners can be useful for prediction but may produce unstable
finite-difference derivatives unless smoothed through the final signal map.

5. Add density support diagnostics.

Store and summarize:

```r
min_density
share_density_below_1e_3
share_density_below_1e_2
density_quantiles
d_log_density_quantiles
correction_quantiles
```

Most of these already exist. The update should make them more visible for
`dml_incremental`.

### Acceptance Criteria

- `dml_incremental` stores derivative step diagnostics.
- `dml_flex` stores derivative step diagnostics when derivative curve is
  requested or computed.
- Tests verify diagnostic fields exist.
- Tests verify warnings can be triggered on a constructed unstable case.

## Phase 7: Incremental Estimator Reference Upgrade

### Goal

Ground `dml_incremental` more explicitly in the average derivative literature
and explain what can be borrowed from related continuous-treatment DML work.

### Current Target

The current scalar signal is:

```text
phi_i =
  partial_g mu_hat(G_i, X_i)
  - partial_g log s_hat(G_i | X_i) * (DeltaY_i - mu_hat(G_i, X_i))
```

This targets:

```text
E[partial_g mu(G_i, X_i)]
```

### Recommended Text Clarification

Add language to Notes/reference docs:

```text
This target is an observed-population average derivative. It is not the
fixed-g derivative curve theta'(g), and it is not generally equal to a simple
marginal-density-weighted average of theta'(g) unless the weighting over X is
handled explicitly.
```

### Candidate Method Extensions

Do not replace the current estimator immediately. Instead add a roadmap:

1. Current density-score estimator:
   - useful because it directly matches the Notes signal;
   - fragile because it needs `partial_g log s(g | x)`.

2. Colangelo-Lee style localized DDML:
   - useful for average dose-response and partial effects at fixed treatment
     values;
   - kernel localizes the continuous treatment;
   - can guide future fixed-g derivative and level inference.

3. ADML / Riesz-based approaches:
   - potentially useful for automatic debiasing of average derivatives;
   - more ambitious and should not be mixed into the current implementation
     without a separate design spec.

4. Derivative-effect-specific doubly robust methods:
   - useful for derivative curves and derivative inference;
   - should be reviewed before changing `dml_incremental`.

### References

Primary references:

- Semenova and Chernozhukov (2021), structural functions and derivatives:
  https://academic.oup.com/ectj/article-pdf/24/2/264/46748086/grf.pdf
- Colangelo and Lee, "Double Debiased Machine Learning Nonparametric Inference
  with Continuous Treatments":
  https://arxiv.org/abs/2004.03036
- Klosin, "Automatic Double Machine Learning for Continuous Treatment Effects":
  https://arxiv.org/abs/2104.10334
- Zhang and Chen, "Doubly Robust Inference on Causal Derivative Effects for
  Continuous Treatments":
  https://arxiv.org/abs/2501.06969
- Chernozhukov, Newey, and Singh, "Automatic Debiased Machine Learning of
  Causal and Structural Effects":
  https://www.econometricsociety.org/publications/econometrica/2022/05/01/automatic-debiased-machine-learning-causal-and-structural

### Acceptance Criteria

- Reference docs clearly separate fixed-G derivative curves from observed
  average derivatives.
- `summary()` for `dml_incremental` describes it as scalar, not curve based.
- `plot(type = "curve")` continues to reject `dml_incremental`.
- `fdid_derivative()` continues to target fixed-G derivative values from
  `kernel` or `dml_flex`, not `dml_incremental`.

## Phase 8: Interflex/Guide Alignment

### Goal

Borrow the strongest reusable DML inference patterns from `interflex` without
misstating FDID's estimand.

### What Interflex Does

The local `interflex-master/R/DML.R` DML path:

- uses `DoubleMLIRM` for discrete treatment;
- uses `DoubleMLPLR` for continuous treatment;
- extracts pseudo-outcomes from DoubleML score components;
- projects pseudo-outcomes onto a B-spline basis of the moderator;
- computes HC0 covariance for the projection;
- computes pointwise intervals;
- computes uniform intervals using a Gaussian max-t simulation from the grid
  covariance matrix.

### What FDID Should Copy

Copy the architecture, not the estimand.

Reusable pieces:

1. Pseudo-outcome extraction discipline:
   - store score components;
   - store pseudo-outcomes;
   - document scaling.

2. BLP/spline projection:
   - stable basis;
   - explicit covariance matrix;
   - grid covariance for contrasts and bands.

3. Gaussian max-t uniform band:
   - simulate from estimated grid covariance;
   - standardize by pointwise SE;
   - report a single critical value.

4. First-stage diagnostics:
   - nuisance losses where available;
   - support/range checks;
   - clear warning when first-stage information is unavailable.

### What FDID Should Not Copy Blindly

Do not make `dml_flex` into interflex's continuous-treatment PLR CME estimator.

Interflex's continuous-treatment DML target is a CME over a moderator `X` under
a PLR structure. FDID's `dml_flex` target is a continuous-G level curve
`theta(g) = E_X[mu(g, X)]`, which needs density adjustment in the current
theory. The two methods can share second-stage inference machinery, but not the
same first-stage score.

### Acceptance Criteria

- Docs include a short "Comparison to interflex" subsection.
- It says interflex DML is moderator-CME oriented.
- It says FDID `dml_flex` is treatment-level-curve oriented.
- It identifies the borrowed piece as BLP/Gaussian covariance inference, not
  the PLR estimand.

## Phase 9: External Software Benchmarks

### Goal

Add sanity checks against established packages where targets overlap.

### Benchmark Candidates

#### DoubleML

Use for:

- `dml_plr` scalar PLR slope;
- `dml_binary` IRM/AIPW style checks;
- multiplier bootstrap design ideas.

Reference:

- https://docs.doubleml.org/dev/guide/se_confint.html

Expected limitation:

- DoubleML does not directly benchmark FDID `dml_flex` dose-response curves.

#### grf

Use for:

- continuous-treatment average partial effect comparison under PLR-like DGPs;
- forest-based diagnostic checks.

Reference:

- https://grf-labs.github.io/grf/reference/average_treatment_effect.html

Expected limitation:

- For continuous treatment, `grf` reports
  `E[Cov(W, Y | X) / Var(W | X)]`, which is closer to a PLR/linear projection
  target than to FDID `E[partial_g mu(G, X)]`.

#### contdid

Use for:

- conceptual and plotting comparison for continuous-treatment DID;
- level versus slope terminology;
- multiplier/bootstrap option language.

References:

- https://bcallaway11.github.io/contdid/
- https://rdrr.io/cran/contdid/man/cont_did.html

Expected limitation:

- `contdid` targets continuous-treatment DID with treated/untreated or timing
  structure; FDID has universal exposure and a different identification map.

#### moderndid

Use for:

- Python-side conceptual comparison for continuous-treatment DID;
- cband/bootstrap API design.

Reference:

- https://moderndid.readthedocs.io/en/latest/api/generated/didcont/moderndid.cont_did.html

Expected limitation:

- Same target mismatch as `contdid`.

#### Colangelo-Lee Code

Use for:

- continuous-treatment DDML simulation design;
- density and kernel localization implementation ideas.

Reference:

- https://github.com/KColangelo/Double-ML-Continuous-Treatment

Expected limitation:

- Not an R package API to import directly into `fdid`; use as methodological
  and simulation reference.

### Acceptance Criteria

- Tests or simulation scripts compare `dml_plr` to DoubleML PLR on simple DGPs.
- Tests or simulation scripts compare `dml_binary` to DoubleML IRM on binary-G
  DGPs.
- Documentation states that no direct package equivalent currently validates
  FDID `dml_flex` and `dml_incremental` targets.

## Phase 10: Reporting Helper Upgrades

### Goal

Make `fdid_contrast()` and `fdid_derivative()` use the best available DML curve
inference.

### Current State

Helpers can use:

- replicate curves when stored;
- kernel covariance matrices when available;
- band-implied intervals;
- approximate endpoint-independent pointwise intervals.

### Recommended Changes

1. For `dml_flex`, use `theta_vcov` when available:

```text
Var(theta(g2)-theta(g1)) = V22 + V11 - 2 V12
```

2. For `dml_flex` derivative values, use `delta_vcov` when available.

3. For interpolated points:
   - either compute interpolation weights and apply them to `V_grid`;
   - or state that covariance inference is grid-only and require
     `interpolate = FALSE` for exact covariance.

Recommendation: implement interpolation weights.

If `g0` lies between grid points `a` and `b`, write:

```text
theta(g0) = w theta(a) + (1-w) theta(b)
```

Then use:

```text
Var(theta(g0)) = w_vec' V_grid w_vec
```

4. Add `inference_source` to helper outputs:

```r
inference_source = c(
  "curve_replicates",
  "curve_vcov",
  "uniform_band",
  "pointwise_independent_approximation"
)
```

### Acceptance Criteria

- DML flex contrasts use covariance when available.
- Helper outputs label the source of inference.
- Tests verify exact grid contrast variance.
- Tests verify interpolated contrast variance.

## Phase 11: Simulation And Validation Plan

### Goal

Validate the updated DML inference layer separately from estimator signal
correctness.

### Simulation Families

#### Family A: PLR Truth

DGP:

```text
G = m(X) + V
DeltaY = beta G + h(X) + U
```

Targets:

- `dml_plr` should recover `beta`;
- `dml_flex` should recover approximately linear `theta(g)`;
- `dml_incremental` should recover `beta`.

Use this family for:

- scalar score SE coverage;
- scalar multiplier coverage;
- DoubleML benchmark.

#### Family B: Nonlinear Smooth Dose Response

DGP:

```text
DeltaY = f(G) + h(X) + interaction(G, X) + U
```

Targets:

- `theta(g)`;
- `theta'(g)`;
- observed scalar `E[partial_g mu(G, X)]`.

Use this family for:

- `dml_flex` level coverage;
- derivative coverage;
- BLP/Gaussian band behavior;
- density-method sensitivity.

#### Family C: Heteroskedastic Treatment Density

DGP:

```text
G | X = m(X) + sigma(X) V
```

Use this family to compare:

- `residual_kde`;
- `location_scale`;
- `local_kde`.

Expected result:

- `location_scale` should improve over `residual_kde` when conditional scale
  variation is strong.

#### Family D: Weak Support / Boundary Stress

DGP:

- sparse tails in `G`;
- nonlinear density;
- evaluation grid near boundaries.

Use this family to test:

- density floor diagnostics;
- derivative-step diagnostics;
- warnings;
- trimmed grid defaults.

#### Family E: Clustered Panel

DGP:

- cluster-level random effects;
- repeated observations by unit or region.

Use this family to test:

- scalar cluster SE;
- curve cluster sandwich in the BLP path.

### Metrics

Report:

- bias;
- RMSE;
- pointwise coverage;
- uniform band coverage;
- average interval width;
- failure rate;
- density tail diagnostics;
- derivative-step sensitivity;
- runtime.

### Acceptance Criteria

- Scalar DML analytical SEs retain existing behavior.
- Scalar multiplier intervals have reasonable coverage in PLR truth.
- `dml_flex` BLP/Gaussian pointwise intervals have reasonable pointwise
  coverage in smooth DGPs.
- Uniform bands are conservative enough in finite grid checks.
- Derivative targets are documented as fragile if coverage remains poor.

## Phase 12: Documentation Updates

### Goal

After code and tests stabilize, update documentation to reflect implemented
behavior precisely.

### Files To Update

Required:

- `FDID Notes/main.tex`;
- `fdid-main/CONTINUOUS_G_ESTIMATOR_REFERENCE.md`;
- `fdid-main/CONTINUOUS_G_PROJECT_LOG.md`;
- `fdid-main/TODO.md`;
- `fdid-main/man/fdid.Rd`;
- `fdid-main/man/fdid_contrast.Rd`;
- `fdid-main/man/fdid_derivative.Rd`;
- `fdid-main/man/plot.fdid.Rd` if plot labels change;
- `fdid-main/R/fdid.R` roxygen comments if regenerated.

Deferred unless the user asks:

- tutorial rewrite;
- `NEWS.md`;
- `HANDOFF.md`.

### Notes Text To Add

Add a DML inference subsection distinguishing:

1. scalar analytical-score inference;
2. scalar score multiplier inference if implemented;
3. `dml_flex` pointwise mapping inference;
4. `dml_flex` practical residual-multiplier bands;
5. `dml_flex` BLP/Gaussian max-t bands if implemented;
6. cluster-aware scalar and curve inference.

### Reference Markdown Text To Add

Add a table:

| Method | Pointwise CI | Uniform band | Cluster-aware | Bootstrap |
|---|---|---|---|---|
| `dml_binary` | analytical score | optional score multiplier if implemented | yes | no row bootstrap |
| `dml_plr` | analytical score | optional score multiplier if implemented | yes | no row bootstrap |
| `dml_incremental` | analytical score | optional score multiplier if implemented | yes | no row bootstrap |
| `dml_flex` | mapping or BLP/sieve | practical multiplier; optional Gaussian max-t | BLP path if implemented | no row bootstrap |

### Acceptance Criteria

- Docs never call practical `dml_flex` bands theorem-level bands.
- Docs never present `dml_incremental` as a curve.
- Docs state which inference paths are analytical, multiplier, or bootstrap.
- Docs list all DML options and which methods they affect.

## Phase 13: API Compatibility And Migration Rules

### Goal

Improve inference without breaking current user code.

### Compatibility Rules

1. Do not change default point estimates.
2. Do not change the default `dml_flex` curve unless the user opts into a new
   `signal_map` or `curve_inference`.
3. Keep `vartype = "robust"` as the accepted DML default.
4. Do not suddenly make `vartype = "bootstrap"` work for DML unless semantics
   are fully implemented and documented.
5. Add new fields without removing old fields:
   - keep `theta_band_lower`;
   - keep `theta_multiplier_reps`;
   - add clearer method labels and covariance fields.
6. Keep `plot.fdid()` behavior compatible, but improve labels.

### Deprecation Warnings

No immediate deprecation is recommended.

Possible future warning:

```text
For DML methods, `boot` controls multiplier replications for curve bands, not
row bootstrap resampling.
```

### Acceptance Criteria

- Existing DML tests pass.
- Existing examples still run.
- New fields are additive.

## Phase 14: Test Plan

### Unit Tests

Add tests for:

- metadata fields;
- scalar multiplier reproducibility;
- cluster scalar multiplier;
- `dml_flex` BLP covariance dimensions;
- `dml_flex` Gaussian max-t critical value;
- `fdid_contrast()` covariance formula;
- `fdid_derivative()` covariance source;
- derivative-step diagnostics;
- warning labels.

### Snapshot/Output Tests

Add tests for:

- `summary(dml_flex)` reports exact band method;
- `summary(dml_incremental)` reports scalar observed-average-derivative target;
- `plot(dml_flex, type = "curve")` labels pointwise versus uniform bands
  correctly.

### Error Tests

Add tests for:

- invalid `dml_inference`;
- invalid `dml_boot`;
- requesting covariance-based contrast when no covariance is stored and
  fallback is disabled;
- too few clusters for cluster curve covariance.

### Numerical Tests

Add small deterministic tests:

- BLP covariance matrix is symmetric positive semidefinite up to tolerance;
- contrast variance equals `V22 + V11 - 2 V12`;
- interpolated variance equals `w' V w`;
- scalar influence variance matches current analytical SE when multiplier is
  disabled.

### Simulation Tests

Keep heavy simulations outside CRAN tests. Store scripts under an existing
validation/simulation directory and summarize results in markdown.

CRAN tests should remain fast:

- small `n`;
- small `boot`;
- deterministic seeds;
- optional package skips.

## Phase 15: Reviewer Checklist

The final reviewer should check:

1. Estimand labels match formulas in Notes.
2. No method silently changes target.
3. `dml_incremental` remains scalar.
4. `dml_flex` curve inference labels are accurate.
5. Analytical, multiplier, and bootstrap language is never mixed.
6. Cluster support is correctly scoped.
7. Tests cover new output fields.
8. Documentation matches actual code.
9. Practical bands are not oversold.
10. External references are cited in the correct places.

## Implementation Order

Recommended order:

1. Add metadata fields.
2. Add scalar DML influence storage.
3. Add optional scalar score multiplier inference.
4. Add BLP/sieve covariance path for `dml_flex`.
5. Add Gaussian max-t uniform bands for BLP/sieve path.
6. Add DML covariance-aware `fdid_contrast()` and `fdid_derivative()`.
7. Add derivative-step diagnostics.
8. Add cluster covariance for BLP/sieve curve path.
9. Add tests.
10. Update Notes, reference docs, TODO, and project log.
11. Run validation/simulation if requested.
12. Review and then decide whether tutorials should be rewritten.

## Open Questions

1. Should `dml_flex` default remain local polynomial, or should the BLP/spline
   path become the recommended inference path?

Recommendation: keep local polynomial default for now; introduce BLP/spline as
an opt-in inference path until validation is complete.

2. Should `boot` continue controlling `dml_flex` multiplier replications?

Recommendation: keep for backward compatibility, but add DML-specific aliases
or clearer docs.

3. Should scalar DML multiplier intervals be exposed through `vartype`?

Recommendation: not initially. Use `dml_inference` to avoid implying row
bootstrap.

4. Should `dml_incremental` get a curve plot?

Recommendation: no. It is scalar. Add diagnostics plots later, but not a
treatment-effect curve.

5. Should tree learners be discouraged for derivative targets?

Recommendation: warn through docs and diagnostics, but do not forbid them.

6. Should the package import DoubleML?

Recommendation: no for core implementation. Use DoubleML only in optional
benchmark scripts or tests skipped when unavailable. The current package
already has native DML code; importing DoubleML would add a large dependency
surface.

## Reference List

### Primary Theory

- Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C.,
  Newey, W., and Robins, J. (2018). "Double/debiased machine learning for
  treatment and structural parameters." The Econometrics Journal.
  https://academic.oup.com/ectj/article/21/1/C1/5056401

- Semenova, V. and Chernozhukov, V. (2021). "Debiased machine learning of
  conditional average treatment effects and other causal functions." The
  Econometrics Journal.
  https://academic.oup.com/ectj/article-pdf/24/2/264/46748086/grf.pdf

- Chernozhukov, V., Chetverikov, D., and Kato, K. (2013). "Gaussian
  approximations and multiplier bootstrap for maxima of sums of high-dimensional
  random vectors."
  https://arxiv.org/abs/1212.6906

### Continuous Treatment And Incremental Effects

- Kennedy, E. H., Ma, Z., McHugh, M. D., and Small, D. S. (2017).
  "Non-parametric methods for doubly robust estimation of continuous treatment
  effects."
  https://arxiv.org/abs/1507.00747

- Colangelo, K. and Lee, Y.-Y. "Double Debiased Machine Learning Nonparametric
  Inference with Continuous Treatments."
  https://arxiv.org/abs/2004.03036

- Klosin, S. "Automatic Double Machine Learning for Continuous Treatment
  Effects."
  https://arxiv.org/abs/2104.10334

- Zhang, Y. and Chen, Y.-C. "Doubly Robust Inference on Causal Derivative
  Effects for Continuous Treatments."
  https://arxiv.org/abs/2501.06969

- Chernozhukov, V., Newey, W., and Singh, R. "Automatic Debiased Machine
  Learning of Causal and Structural Effects."
  https://www.econometricsociety.org/publications/econometrica/2022/05/01/automatic-debiased-machine-learning-causal-and-structural

### Continuous DID References

- Callaway, B., Goodman-Bacon, A., and Sant'Anna, P. H. C. "Difference-in-
  Differences with a Continuous Treatment."
  https://bcallaway11.github.io/contdid/

- `contdid` package reference:
  https://rdrr.io/cran/contdid/man/cont_did.html

- `moderndid` continuous-DID documentation:
  https://moderndid.readthedocs.io/en/latest/api/generated/didcont/moderndid.cont_did.html

### Software References

- `interflex` local implementation:
  - `interflex-master/R/DML.R`;
  - `interflex-master/R/uniform.R`;
  - `interflex-master/vignettes/05_dml.qmd`.

- DoubleML variance and multiplier bootstrap documentation:
  https://docs.doubleml.org/dev/guide/se_confint.html

- `grf` average treatment effect / continuous-treatment average partial effect
  documentation:
  https://grf-labs.github.io/grf/reference/average_treatment_effect.html

- Colangelo-Lee continuous-treatment DDML code:
  https://github.com/KColangelo/Double-ML-Continuous-Treatment

## Final Recommendation

The next DML implementation run should prioritize:

```text
1. scalar DML influence storage and optional score multiplier inference;
2. dml_flex BLP/spline covariance path;
3. Gaussian max-t uniform bands for dml_flex;
4. covariance-aware fdid_contrast() and fdid_derivative();
5. derivative and density stability diagnostics;
6. documentation updates that clearly label every inference path.
```

This sequence improves inference in the areas the audit identified as weakest
while preserving the current, theoretically coherent estimand split.
