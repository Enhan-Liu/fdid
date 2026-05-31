# Kernel Inference Update Plan

Date: 2026-05-30

This plan records the proposed update to the `fdid` continuous-`G` kernel
estimator. The update is motivated by the current `fdid-main` implementation,
the FDID Notes, and the kernel inference design in `interflex` and the
Practical Guide to Estimating Conditional Marginal Effects.

## Objective

The point estimator in `R/kernel.R` is already close to the local-linear SVCM
kernel estimator described in the FDID Notes and in the Practical Guide. The
main gap is the inference layer. The update should keep the existing estimator
structure stable and systematically improve:

- analytical pointwise inference;
- bootstrap pointwise inference;
- bootstrap covariance storage for level and derivative curves;
- covariance-aware level and derivative contrasts;
- bootstrap uniform bands;
- cluster-aware bootstrap resampling;
- documentation, tests, and project logs.

The first implementation target is parity with the current `interflex` kernel
practice: bootstrap pointwise intervals, bootstrap curve covariance matrices,
and bootstrap quantile-envelope uniform bands. A stricter studentized sup-t band
should be documented separately and either implemented as an explicit optional
mode or retained as a second-stage theory-aligned extension.

## Current State

The current `fdid` kernel implementation:

- estimates the level curve `mu_hat(g)` and derivative curve `delta_hat(g)` on
  an evaluation grid;
- uses a local-linear WLS design with Gaussian kernel weights and adaptive
  bandwidths;
- selects the baseline bandwidth by LSCV when `h0 = NULL`;
- computes analytical sandwich standard errors for `mu_hat(g)` and
  `delta_hat(g)`;
- when `vartype = "bootstrap"`, resamples rows, reselects bandwidths, refits
  the kernel curves, and stores bootstrap replicate matrices for `mu` and
  `delta`;
- forms a practical derivative band by calibrating empirical bootstrap quantile
  envelopes.

Important limitations:

- bootstrap pointwise standard errors and percentile confidence intervals are
  not systematically exposed;
- the practical uniform band is not a studentized sup-t band;
- only the derivative band is exposed through legacy top-level `band_lower` and
  `band_upper` fields;
- bootstrap covariance matrices for the full `mu` and `delta` curves are not
  stored;
- scalar summaries and level contrasts can still fall back to endpoint
  independence approximations;
- row bootstrap is used even when cluster labels are available;
- dynamic-curve bootstrap inference is not systematically implemented;
- documentation still risks conflating quantile-envelope bands with
  studentized sup-t bands.

## Reference: interflex

The relevant `interflex` behavior is:

- kernel bootstrap stores curve replicate matrices;
- pointwise bootstrap standard errors are computed as row-wise standard
  deviations of the bootstrap curve matrix;
- pointwise bootstrap confidence intervals are percentile intervals;
- uniform bands are produced by `calculate_uniform_quantiles()`, which searches
  for a quantile level `zeta` such that the empirical rectangular envelope
  covers whole bootstrap curves with probability at least `1 - alpha`;
- bootstrap covariance matrices are stored as `cov(t(curve_boot_matrix))`;
- `inter.test()` uses the covariance identity
  `Var(theta(x2) - theta(x1)) = V22 + V11 - 2 V12` for contrasts;
- cluster bootstrap is supported by resampling clusters and keeping all rows
  within selected clusters.

The Practical Guide text describes interflex uniform intervals as a sup-t band
approach introduced by Montiel Olea and Plagborg-Moller, but the procedure
shown in the text and implemented in `interflex` is an empirical bootstrap
quantile-envelope calibration. It is not the same as a fully studentized
max-t statistic of the form:

```text
sup_g |(theta_hat^*(g) - theta_hat(g)) / se_hat^*(g)|
```

Therefore, the FDID implementation and documentation should use precise names:

- `bootstrap_quantile_envelope` for the interflex-style band;
- `studentized_sup_t` only for a true studentized maximum-statistic band.

## R Code Update Plan

### 1. Refactor Kernel Inference Helpers

File: `R/kernel.R`

Introduce internal helpers instead of keeping all bootstrap logic inside one
large routine:

```r
kernel_bootstrap_curves()
kernel_pointwise_ci()
kernel_uniform_quantile_envelope()
kernel_curve_vcov()
kernel_curve_contrast()
kernel_studentized_sup_t_band()
```

The point-estimation helpers should remain separate from inference helpers so
that robust, bootstrap, and future studentized inference paths share the same
curve estimator.

### 2. Replace the Existing Uniform-Band Search

Replace the current coarse-grid `find_zeta_star()` logic with an
interflex-style bisection helper.

Inputs:

- `theta_matrix`: matrix with dimensions `grid_points x bootstrap_draws`;
- `alpha`;
- optional tolerance and maximum iterations.

Behavior:

- drop invalid bootstrap columns;
- set the search interval to `[alpha / (2 * k), alpha / 2]`, where `k` is the
  number of grid points;
- for each candidate `zeta`, compute pointwise empirical quantiles;
- compute the fraction of bootstrap curves lying entirely inside the envelope;
- return the calibrated envelope, `zeta_hat`, realized coverage, and effective
  number of bootstrap draws.

Return structure:

```r
list(
  lower = ...,
  upper = ...,
  zeta_hat = ...,
  coverage = ...,
  n_boot_eff = ...,
  method = "bootstrap_quantile_envelope"
)
```

### 3. Add Bootstrap Pointwise Inference

When `vartype = "bootstrap"`, expose:

```r
curve_event$mu_boot_se
curve_event$delta_boot_se
curve_event$mu_boot_ci_lower
curve_event$mu_boot_ci_upper
curve_event$delta_boot_ci_lower
curve_event$delta_boot_ci_upper
```

Default pointwise bootstrap confidence intervals should be percentile
intervals:

```r
quantile(theta_boot, probs = c(alpha / 2, 1 - alpha / 2))
```

Optionally, a future public option can distinguish:

```r
boot_ci = c("percentile", "se")
```

For the initial update, percentile intervals are preferred because they match
the current `interflex` kernel bootstrap output.

### 4. Add Uniform Bands for Both Level and Derivative Curves

Expose both level and derivative uniform bands:

```r
curve_event$mu_band_lower
curve_event$mu_band_upper
curve_event$delta_band_lower
curve_event$delta_band_upper
curve_event$mu_uniform_method
curve_event$delta_uniform_method
curve_event$mu_uniform_zeta
curve_event$delta_uniform_zeta
curve_event$mu_uniform_coverage
curve_event$delta_uniform_coverage
```

For backward compatibility:

```r
out$band_lower <- curve_event$delta_band_lower
out$band_upper <- curve_event$delta_band_upper
```

If `vartype = "robust"`, these top-level fields should continue to represent
pointwise normal derivative intervals unless the plotting and documentation are
updated to distinguish pointwise intervals from uniform bands.

### 5. Store Bootstrap Curve Covariance Matrices

For bootstrap fits, store:

```r
curve_event$mu_boot_vcov
curve_event$delta_boot_vcov
```

where:

```r
mu_boot_vcov <- cov(mu_boot_matrix, use = "pairwise.complete.obs")
delta_boot_vcov <- cov(delta_boot_matrix, use = "pairwise.complete.obs")
```

The chosen orientation should be documented and tested. If bootstrap matrices
are stored as `B x grid`, then `cov(mu_boot_matrix)` produces a `grid x grid`
covariance matrix.

### 6. Fix Level and Derivative Contrast Inference

For any grid contrast:

```text
theta(g2) - theta(g1)
```

use:

```text
Var(theta(g2) - theta(g1)) = V22 + V11 - 2 V12
```

This should apply to:

- level contrasts using `mu_boot_vcov`;
- derivative contrasts using `delta_boot_vcov`;
- scalar endpoint summaries such as `(mu(g_U) - mu(g_L)) / (g_U - g_L)`.

Recommended priority:

1. If bootstrap replicate curves are available, compute contrast replicates
   directly and derive SE and percentile intervals from those replicates.
2. If only a bootstrap covariance matrix is available, use
   `V22 + V11 - 2 V12`.
3. If only analytical pointwise SEs are available, keep the current independent
   endpoint fallback but label it clearly as an approximation.

Suggested metadata:

```r
est$se_method
est$ci_method
est$contrast_vcov_used
est$contrast_inference_note
```

### 7. Add a Kernel Contrast Helper

Add an internal helper:

```r
kernel_curve_contrast(curve, g1, g2, target = c("mu", "delta"))
```

First version:

- supports grid points exactly;
- optionally uses nearest-grid matching with a warning;
- uses bootstrap replicates first;
- otherwise uses bootstrap covariance;
- otherwise falls back to analytical endpoint approximation.

Avoid copying `interflex::inter.test()` covariance-surface smoothing in the
first implementation. A transparent grid-based contrast is easier to test and
less likely to introduce smoothing artifacts.

### 8. Add Cluster-Aware Bootstrap

Extend `run_kernel_method()` to accept:

```r
cluster = NULL
```

Dispatcher updates:

- pass the cluster labels from the prepared FDID object when available;
- preserve existing behavior when no cluster labels exist.

Bootstrap behavior:

- if `cluster = NULL`, use row bootstrap;
- if `cluster` is supplied, sample clusters with replacement and include all
  rows in each sampled cluster;
- record the bootstrap resampling mode.

Suggested output:

```r
out$bootstrap_resampling <- "row" or "cluster"
out$n_clusters <- ...
curve_event$bootstrap_resampling <- ...
```

### 9. Add a Studentized Sup-t Option Separately

Add, or reserve, a public option:

```r
uniform_band = c("quantile", "studentized", "none")
```

Recommended default:

```r
uniform_band = "quantile"
```

because this is the current `interflex`-style behavior and can be tested
immediately.

A first studentized implementation could use:

```text
T_b = sup_g |(delta_boot[b, g] - delta_hat[g]) / se_delta[g]|
c_alpha = quantile(T_b, 1 - alpha)
band(g) = delta_hat(g) +/- c_alpha * se_delta(g)
```

A stricter version would recompute or store bootstrap-specific standard errors
`se_delta_boot[b, g]` and use them in the denominator. That version is more
expensive and should be documented as a separate theory-aligned extension.

The documentation must not call the quantile-envelope band a studentized sup-t
band.

## Public API Considerations

Minimal API change:

- keep `vartype = c("robust", "bootstrap")`;
- keep `boot`, `alpha`, `trim`, `eval_g`, `h0`, and `K_folds`;
- add optional `uniform_band = c("quantile", "studentized", "none")` only if
  studentized or no-band behavior is implemented;
- add optional `boot_ci = c("percentile", "se")` only if both pointwise CI
  types are implemented;
- add optional `cluster` only internally if the dispatcher can pass cluster
  labels from prepared data.

Backward-compatible output:

- retain `band_lower` and `band_upper`;
- retain `curve_event$mu_bootstrap` and `curve_event$delta_bootstrap`;
- add new fields without removing old fields.

## Tests

File: `tests/testthat/test-continuous.R`

Add or extend tests for:

1. Robust analytical inference:
   - `se_mu` and `se_delta` exist;
   - pointwise normal intervals exist;
   - bootstrap-only fields are absent or `NULL`.

2. Bootstrap pointwise inference:
   - `mu_bootstrap` and `delta_bootstrap` dimensions are correct;
   - bootstrap SE fields exist;
   - pointwise percentile intervals exist.

3. Uniform quantile-envelope inference:
   - level and derivative bands have length `length(eval_g)`;
   - `zeta_hat`, `coverage`, and `n_boot_eff` are stored;
   - lower bands are not above upper bands.

4. Curve covariance matrices:
   - `mu_boot_vcov` and `delta_boot_vcov` are square `grid x grid` matrices;
   - diagonal entries match bootstrap variances up to tolerance.

5. Contrast covariance:
   - with a synthetic bootstrap matrix, verify that
     `V22 + V11 - 2 V12` matches the variance of bootstrap contrast replicates.

6. Scalar endpoint summary:
   - when bootstrap replicates are available, scalar SE uses replicate
     contrasts or covariance-aware endpoint variance;
   - independent-endpoint fallback is only used when covariance is unavailable.

7. Cluster bootstrap:
   - cluster labels trigger cluster resampling;
   - all rows in selected clusters are retained together;
   - bootstrap metadata records `"cluster"`.

8. Backward compatibility:
   - `out$band_lower` and `out$band_upper` remain available;
   - `plot(type = "curve")` and `summary()` continue to work.

## Documentation Updates

### FDID Notes

File: `../FDID Notes/main.tex`

Update the kernel inference section to distinguish:

- analytical pointwise sandwich inference;
- bootstrap pointwise inference;
- bootstrap quantile-envelope uniform bands;
- studentized sup-t uniform bands.

Add the level contrast variance identity:

```text
Var(mu_hat(g2) - mu_hat(g1)) = V(g2, g2) + V(g1, g1) - 2 V(g2, g1).
```

Explain that `V(gj, gk)` can be estimated from bootstrap curve replicates.

Clarify implementation status:

- implemented or planned first: interflex-style quantile-envelope band;
- theory-aligned target: studentized sup-t band;
- these are related but not identical.

### Continuous-G Estimator Reference

File: `CONTINUOUS_G_ESTIMATOR_REFERENCE.md`

Add a kernel inference table:

| Option or field | Meaning |
| --- | --- |
| `vartype = "robust"` | analytical pointwise sandwich inference |
| `vartype = "bootstrap"` | bootstrap curve inference |
| `*_boot_se` | pointwise bootstrap standard errors |
| `*_boot_ci_lower`, `*_boot_ci_upper` | pointwise bootstrap percentile intervals |
| `*_band_lower`, `*_band_upper` | uniform or pointwise bands, depending on method |
| `*_boot_vcov` | bootstrap curve covariance matrix |
| `uniform_method` | e.g. `bootstrap_quantile_envelope` |

Document the contrast rule:

```text
V22 + V11 - 2 V12
```

Document resampling:

- row bootstrap when no cluster labels exist;
- cluster bootstrap when cluster labels are available.

### Project Log and TODO

Files:

- `CONTINUOUS_G_PROJECT_LOG.md`
- `TODO.md`

Record:

- bootstrap pointwise SE and percentile CI added for kernel `mu` and `delta`;
- interflex-style quantile-envelope uniform bands added or clarified;
- bootstrap curve covariance matrices added;
- level and derivative contrasts updated to use covariance-aware inference;
- cluster bootstrap support added if implemented;
- studentized sup-t remains a separate optional/theory-aligned extension if not
  implemented in the first pass.

Suggested log entry:

```text
Kernel inference update:
- aligned bootstrap uniform bands with interflex quantile-envelope calibration;
- added bootstrap pointwise SE/CI for level and derivative curves;
- added bootstrap curve covariance matrices;
- updated endpoint and curve contrasts to use covariance-aware inference;
- added cluster bootstrap support when cluster labels are available;
- clarified documentation so quantile-envelope bands are not described as
  studentized sup-t bands.
```

## Validation

Minimum validation before marking complete:

1. Run focused testthat tests for continuous-`G` kernel.
2. Confirm existing plots still render with legacy `band_lower` and
   `band_upper`.
3. Confirm small-bootstrap tests do not rely on exact numerical coverage.
4. Confirm docs no longer conflate quantile-envelope and studentized sup-t
   bands.

Recommended full validation:

```sh
R CMD check fdid-main
```

or the equivalent package check used in the local workflow.

## Recommended Implementation Order

1. Refactor kernel inference helpers.
2. Add bootstrap pointwise SE and percentile intervals.
3. Replace uniform-band calibration with interflex-style bisection.
4. Store `mu` and `delta` bootstrap covariance matrices.
5. Fix scalar endpoint summaries and grid contrasts using covariance-aware
   formulas.
6. Add cluster bootstrap if cluster labels are available from the dispatcher.
7. Update plot, summary, and reporting code for new fields while keeping legacy
   fields.
8. Add tests.
9. Update FDID Notes, the continuous-G markdown reference, project log, and
   TODO.
10. Consider studentized sup-t as a separate optional mode.

## Practical Guide Sup-t Note

The Practical Guide does mention a sup-t band approach for uniform confidence
intervals in the classic/kernel discussion. However, the procedure it lays out
for `interflex` is the quantile-envelope procedure:

1. compute bootstrap curves on a grid;
2. compute empirical pointwise quantiles at each grid point;
3. choose `zeta` so that whole bootstrap curves fall inside the rectangular
   quantile band with probability at least `1 - alpha`;
4. report the resulting simultaneous band.

This matches the `interflex` code path in `R/uniform.R`. It is related to the
same simultaneous-inference motivation as sup-t bands, but it is not the same
as explicitly computing a studentized maximum statistic over the grid. FDID
documentation should therefore distinguish:

- what the Practical Guide calls or motivates as a sup-t-style uniform band;
- what `interflex` actually implements: bootstrap quantile-envelope
  calibration;
- what FDID Notes define as a studentized sup-t statistic for derivative
  curves.
