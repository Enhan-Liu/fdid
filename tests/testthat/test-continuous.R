library(fdid)

# ===========================================================================
# Continuous-G method tests (R3–R7 from regression_tests.md)
# ===========================================================================

# Shared synthetic panel: continuous G, 60 units × 3 periods
set.seed(42)
n_units <- 60
times   <- 1:3

cran_data_cont <- data.frame(
  id      = rep(seq_len(n_units), each = length(times)),
  time    = rep(times, times = n_units),
  G       = rep(runif(n_units, 0, 1), each = length(times)),
  outcome = rnorm(n_units * length(times)),
  covar1  = rnorm(n_units * length(times))
)
s_cont <- fdid_prepare(
  data       = cran_data_cont,
  Y_label    = "outcome",
  X_labels   = "covar1",
  G_label    = "G",
  unit_label = "id",
  time_label = "time"
)

test_that("[CRAN] fdid_prepare expands covariates with basis functions", {
  d_x <- data.frame(
    id = rep(1:50, each = 3),
    time = rep(1:3, 50),
    G = rep(runif(50), each = 3),
    outcome = rnorm(150),
    covar1 = rep(rnorm(50), each = 3),
    covar2 = rep(runif(50), each = 3)
  )
  s_raw <- fdid_prepare(d_x, "outcome", c("covar1", "covar2"),
                        "G", "id", "time")
  s_poly <- fdid_expand_covariates(
    s_raw,
    basis_type = "polynomial",
    poly_degree = 2L,
    include_interactions = TRUE
  )
  s_bs <- fdid_prepare(
    d_x, "outcome", c("covar1", "covar2"), "G", "id", "time",
    basis_type = "bspline", spline_df = 3L, spline_degree = 2L
  )

  expect_equal(length(grep("^x[0-9]+$", names(s_raw))), 2L)
  expect_gt(length(grep("^x[0-9]+$", names(s_poly))), 2L)
  expect_gt(length(grep("^x[0-9]+$", names(s_bs))), 2L)
  expect_equal(attr(s_poly, "fdid_basis")$basis_type, "polynomial")
  expect_true(all(is.finite(as.matrix(s_poly[, grep("^x[0-9]+$", names(s_poly))]))))
})

test_that("[CRAN] flexible DML works with expanded high-dimensional covariates", {
  d_x <- data.frame(
    id = rep(1:70, each = 3),
    time = rep(1:3, 70),
    G = rep(runif(70), each = 3),
    outcome = rnorm(210),
    covar1 = rep(rnorm(70), each = 3),
    covar2 = rep(rnorm(70), each = 3),
    covar3 = rep(rnorm(70), each = 3)
  )
  s_hd <- fdid_prepare(
    d_x, "outcome", c("covar1", "covar2", "covar3"), "G", "id", "time",
    basis_type = "polynomial", poly_degree = 2L,
    include_interactions = TRUE
  )
  expect_gt(length(grep("^x[0-9]+$", names(s_hd))), 6L)
  evg <- seq(0.2, 0.8, length.out = 3)
  res_f <- fdid(s_hd, tr_period = 2L, ref_period = 1L,
                entire_period = 1L, method = "dml_flex",
                K = 3L, S = 1L, eval_g = evg, boot = 3L)
  res_i <- fdid(s_hd, tr_period = 2L, ref_period = 1L,
                entire_period = 1L, method = "dml_incremental",
                K = 3L, S = 1L)
  expect_s3_class(res_f, "fdid")
  expect_s3_class(res_i, "fdid")
  expect_true(all(is.finite(res_f$curve_event$theta_hat)))
  expect_true(all(is.finite(res_f$curve_event$delta_hat)))
  expect_true(is.finite(res_i$est$event$Estimate))
  expect_true(!is.null(res_i$incremental_event$signal))
})

# Required fields all new methods must have (superset check with %in%)
REQUIRED_FIELDS <- c("est", "dynamic", "raw_means", "tr_period", "ref_period",
                      "entire_period", "method", "vartype", "times", "G", "ps",
                      "call", "target.pop")

# ===========================================================================
# R3 — kernel structural checks
# ===========================================================================

test_that("[CRAN] fdid kernel: returns valid fdid object", {
  res <- fdid(s_cont, tr_period = 2L, ref_period = 1L, method = "kernel")
  expect_s3_class(res, "fdid")
  expect_true(all(REQUIRED_FIELDS %in% names(res)))
  expect_true(all(c("eval_g", "curve_event", "h0", "h_adaptive") %in% names(res)))
  expect_equal(res$method, "kernel")
  expect_true(is.finite(res$est$event$Estimate))
  expect_true(is.finite(res$est$event$Std.Error))
  expect_true(nrow(res$dynamic) > 0)
})

test_that("[CRAN] fdid kernel: curve fields populated", {
  res <- fdid(s_cont, tr_period = 2L, ref_period = 1L, method = "kernel")
  expect_true(!is.null(res$eval_g) && length(res$eval_g) > 0)
  expect_equal(length(res$curve_event$delta_hat), length(res$eval_g))
  expect_equal(length(res$curve_event$se_delta),  length(res$eval_g))
  expect_equal(length(res$curve_event$mu_hat),    length(res$eval_g))
  expect_true(is.matrix(res$curve_event$mu_vcov))
  expect_true(is.matrix(res$curve_event$delta_vcov))
  expect_equal(dim(res$curve_event$mu_vcov), c(length(res$eval_g), length(res$eval_g)))
  expect_equal(dim(res$curve_event$delta_vcov), c(length(res$eval_g), length(res$eval_g)))
})

test_that("[CRAN] fdid kernel: h0 is finite positive", {
  res <- fdid(s_cont, tr_period = 2L, ref_period = 1L, method = "kernel")
  expect_true(is.finite(res$h0) && res$h0 > 0)
})

test_that("[CRAN] fdid kernel: custom eval_g is respected", {
  evg <- seq(0.2, 0.8, length.out = 7)
  res <- fdid(s_cont, tr_period = 2L, ref_period = 1L, method = "kernel", eval_g = evg)
  expect_equal(length(res$eval_g), 7L)
  expect_equal(length(res$curve_event$delta_hat), 7L)
})

test_that("[CRAN] fdid kernel: single eval_g has finite event summary", {
  res <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
              method = "kernel", eval_g = 0.5)
  expect_true(is.finite(res$est$event$Estimate))
  expect_true(is.finite(res$est$event$Std.Error))
})

test_that("[CRAN] fdid kernel: vartype=bootstrap runs without error", {
  evg <- seq(0.2, 0.8, length.out = 5)
  res <- suppressWarnings(
    fdid(s_cont, tr_period = 2L, ref_period = 1L, method = "kernel",
         vartype = "bootstrap", eval_g = evg, boot = 5L)
  )
  expect_s3_class(res, "fdid")
  expect_true(!is.null(res$band_lower))
  expect_equal(length(res$band_lower), length(res$eval_g))
  expect_true(is.matrix(res$curve_event$mu_bootstrap))
  expect_true(is.matrix(res$curve_event$delta_bootstrap))
  expect_equal(ncol(res$curve_event$mu_bootstrap), length(res$eval_g))
  expect_equal(length(res$curve_event$mu_boot_se), length(res$eval_g))
  expect_equal(length(res$curve_event$delta_boot_ci_lower), length(res$eval_g))
  expect_equal(length(res$curve_event$mu_band_lower), length(res$eval_g))
  expect_equal(length(res$curve_event$delta_band_upper), length(res$eval_g))
  expect_true(is.matrix(res$curve_event$mu_boot_vcov))
  expect_equal(dim(res$curve_event$mu_boot_vcov), c(length(res$eval_g), length(res$eval_g)))
  expect_equal(res$curve_event$delta_uniform_method, "bootstrap_quantile_envelope")
  expect_equal(res$bootstrap_resampling, "row")
  expect_equal(res$est$event$SE_Method, "bootstrap_replicate")
})

test_that("[CRAN] fdid kernel: bootstrap curve vcov supports covariance-aware contrasts", {
  evg <- seq(0.2, 0.8, length.out = 4)
  res <- suppressWarnings(
    fdid(s_cont, tr_period = 2L, ref_period = 1L, method = "kernel",
         vartype = "bootstrap", eval_g = evg, boot = 5L)
  )
  V <- res$curve_event$mu_boot_vcov
  expected <- sqrt(pmax(0, V[4, 4] + V[1, 1] - 2 * V[4, 1]))
  con <- fdid_contrast(res, g0 = evg[1], g1 = evg[4], inference = "vcov")
  expect_equal(con$inference, "vcov")
  expect_equal(con$std.error, expected, tolerance = 1e-8)
})

test_that("[CRAN] fdid kernel: robust curve vcov supports covariance-aware contrasts", {
  evg <- seq(0.2, 0.8, length.out = 4)
  res <- fdid(s_cont, tr_period = 2L, ref_period = 1L, method = "kernel",
              vartype = "robust", eval_g = evg)
  V <- res$curve_event$mu_vcov
  span <- evg[4] - evg[1]
  expected_contrast <- sqrt(pmax(0, V[4, 4] + V[1, 1] - 2 * V[4, 1]))
  expected_scalar <- expected_contrast / span
  con <- fdid_contrast(res, g0 = evg[1], g1 = evg[4], inference = "vcov")
  expect_equal(con$inference, "vcov")
  expect_equal(con$std.error, expected_contrast, tolerance = 1e-8)
  expect_equal(res$est$event$Std.Error, expected_scalar, tolerance = 1e-8)
  expect_equal(res$est$event$SE_Method, "analytical_stacked_sandwich")
})

test_that("[CRAN] fdid kernel: cluster bootstrap records cluster resampling", {
  set.seed(225)
  n_cl <- 40
  d_cl <- data.frame(
    id = rep(seq_len(n_cl), each = 3),
    time = rep(1:3, n_cl),
    cluster = rep(rep(seq_len(8), each = 5), each = 3),
    G = rep(stats::runif(n_cl), each = 3),
    outcome = stats::rnorm(n_cl * 3),
    covar1 = rep(stats::rnorm(n_cl), each = 3)
  )
  s_cl <- fdid_prepare(d_cl, "outcome", "covar1", "G", "id", "time",
                       cluster_label = "cluster")
  res <- suppressWarnings(
    fdid(s_cl, tr_period = 2L, ref_period = 1L, method = "kernel",
         vartype = "bootstrap", eval_g = seq(0.25, 0.75, length.out = 3),
         boot = 3L)
  )
  expect_equal(res$bootstrap_resampling, "cluster")
  expect_equal(res$n_clusters, 8L)
  expect_true(is.matrix(res$curve_event$delta_boot_vcov))
})

# ===========================================================================
# R4 — DML structural checks
# ===========================================================================

# Binary G dataset for dml_binary
set.seed(42)
n_bin <- 60
d_bin <- data.frame(
  id      = rep(1:n_bin, each = 3),
  time    = rep(1:3, n_bin),
  G       = rep(rbinom(n_bin, 1, 0.5), each = 3),
  outcome = rnorm(n_bin * 3),
  covar1  = rnorm(n_bin * 3)
)
s_bin <- fdid_prepare(d_bin, "outcome", "covar1", "G", "id", "time")

test_that("[CRAN] fdid dml_binary: returns valid fdid object", {
  res <- fdid(s_bin, tr_period = 2L, ref_period = 1L,
              method = "dml_binary", K = 3L, S = 1L)
  expect_s3_class(res, "fdid")
  expect_true(all(REQUIRED_FIELDS %in% names(res)))
  expect_true(all(c("dml_method", "K", "S", "learner") %in% names(res)))
  expect_equal(res$method, "dml_binary")
  expect_true(is.finite(res$est$event$Estimate))
  expect_true(is.finite(res$est$event$Std.Error))
  expect_true(nrow(res$dynamic) > 0)
})

test_that("[CRAN] fdid dml_plr: returns valid fdid object", {
  res <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
              method = "dml_plr", K = 3L, S = 1L)
  expect_s3_class(res, "fdid")
  expect_true(all(REQUIRED_FIELDS %in% names(res)))
  expect_equal(res$method, "dml_plr")
  expect_true(is.finite(res$est$event$Estimate))
  expect_true(is.finite(res$est$event$Std.Error))
  expect_true(nrow(res$dynamic) > 0)
})

test_that("[CRAN] scalar DML stores metadata and score multiplier inference", {
  res <- fdid(
    s_cont, tr_period = 2L, ref_period = 1L,
    method = "dml_plr", K = 3L, S = 1L,
    dml_inference = "score_multiplier", dml_boot = 5L,
    dml_multiplier = "normal"
  )
  expect_equal(res$dml_inference, "score_multiplier")
  expect_equal(res$dml_boot, 5L)
  expect_equal(res$dml_metadata$target_estimand, "partially_linear_slope")
  expect_equal(res$scalar_event$scalar_inference, "score_multiplier")
  expect_true(is.numeric(res$scalar_event$influence))
  expect_equal(length(res$scalar_event$split_estimates), 1L)
  expect_equal(length(res$scalar_event$scalar_multiplier_reps), 5L)
  expect_true(is.finite(res$scalar_event$scalar_multiplier_crit))
})

test_that("[CRAN] fdid DML accepts expanded learner aliases", {
  skip_if_not_installed("glmnet")
  res_ridge <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
                    method = "dml_plr", K = 3L, S = 1L,
                    learner = "ridge")
  res_enet <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
                   method = "dml_plr", K = 3L, S = 1L,
                   learner = "elasticnet")
  expect_s3_class(res_ridge, "fdid")
  expect_s3_class(res_enet, "fdid")
  expect_true(is.finite(res_ridge$est$event$Estimate))
  expect_true(is.finite(res_enet$est$event$Estimate))
})

test_that("[CRAN] DML learner helpers support optional ML backends", {
  set.seed(99)
  X_tr <- matrix(rnorm(120), ncol = 3)
  X_new <- matrix(rnorm(30), ncol = 3)
  y_tr <- 0.5 * X_tr[, 1] - X_tr[, 2] + rnorm(nrow(X_tr), sd = 0.2)
  g_tr <- as.integer(y_tr > stats::median(y_tr))

  optional <- list(grf = "grf", gam = "mgcv", nnet = "nnet")
  for (learner_name in names(optional)) {
    if (!requireNamespace(optional[[learner_name]], quietly = TRUE)) next
    pred <- fdid:::.fit_reg(y_tr, X_tr, X_new, learner_name)
    ps <- fdid:::.fit_classif(g_tr, X_tr, X_new, learner_name)
    expect_equal(length(pred), nrow(X_new))
    expect_equal(length(ps), nrow(X_new))
    expect_true(all(is.finite(pred)))
    expect_true(all(is.finite(ps)))
    expect_true(all(ps >= 0.01 & ps <= 0.99))
  }

  if (!requireNamespace("xgboost", quietly = TRUE)) {
    expect_error(fdid:::.fit_reg(y_tr, X_tr, X_new, "xgboost"),
                 regexp = "xgboost")
  }
})

test_that("[CRAN] fdid dml_flex: returns curve fields", {
  evg <- seq(0.2, 0.8, length.out = 5)
  res <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
              method = "dml_flex", K = 3L, S = 1L, eval_g = evg, boot = 3L)
  expect_s3_class(res, "fdid")
  expect_true(all(REQUIRED_FIELDS %in% names(res)))
  expect_equal(res$method, "dml_flex")
  expect_true(is.finite(res$est$event$Estimate))
  expect_true(!is.null(res$curve_event$theta_hat))
  expect_equal(length(res$curve_event$theta_hat), 5L)
  expect_true(!is.null(res$curve_event$delta_hat))
  expect_equal(length(res$curve_event$delta_hat), 5L)
  expect_true(all(is.finite(res$curve_event$se_theta)))
  expect_true(all(is.finite(res$curve_event$se_delta)))
  expect_true(nrow(res$dynamic) > 0)
  expect_equal(res$signal_map, "local_poly")
  expect_equal(res$density_method, "residual_kde")
  expect_true(is.matrix(res$curve_event$theta_multiplier_reps))
  expect_true(is.matrix(res$curve_event$delta_multiplier_reps))
  expect_equal(ncol(res$curve_event$theta_multiplier_reps), length(res$eval_g))
})

test_that("[CRAN] fdid dml_flex blp_spline stores covariance-aware inference", {
  evg <- seq(0.25, 0.75, length.out = 4)
  res <- fdid(
    s_cont, tr_period = 2L, ref_period = 1L,
    method = "dml_flex", K = 3L, S = 1L, eval_g = evg,
    signal_map = "blp_spline", map_degree = 2L, map_df = 4L,
    boot = 5L
  )
  expect_equal(res$signal_map, "blp_spline")
  expect_equal(res$dml_metadata$target_estimand, "fixed_G_level_curve")
  expect_true(is.matrix(res$curve_event$theta_vcov))
  expect_true(is.matrix(res$curve_event$delta_vcov))
  expect_equal(dim(res$curve_event$theta_vcov), c(length(evg), length(evg)))
  expect_equal(dim(res$curve_event$delta_vcov), c(length(evg), length(evg)))
  expect_equal(res$curve_event$simultaneous_band_method,
               "gaussian_max_t_blp_spline_covariance")

  con <- fdid_contrast(res, g0 = evg[1], g1 = evg[4], inference = "vcov")
  V <- res$curve_event$theta_vcov
  expected_se <- sqrt(pmax(0, V[4, 4] + V[1, 1] - 2 * V[4, 1]))
  expect_equal(con$std.error, expected_se, tolerance = 1e-8)
  expect_equal(con$inference, "vcov")
  expected_interval_se <- expected_se / (evg[4] - evg[1])
  expect_equal(unname(res$curve_event$interval_average["Std.Error"]),
               expected_interval_se, tolerance = 1e-8)
  expect_equal(unname(res$est$event$Estimate),
               unname(res$curve_event$interval_average["Estimate"]),
               tolerance = 1e-8)
  expect_equal(res$est$event$SE_Method, "curve_vcov")

  der <- fdid_derivative(res, g0 = evg[2], inference = "vcov")
  expect_equal(der$std.error, sqrt(pmax(0, res$curve_event$delta_vcov[2, 2])),
               tolerance = 1e-8)
  expect_equal(der$inference, "vcov")
})

test_that("[CRAN] continuous-G reporting helpers extract contrasts and derivatives", {
  evg <- seq(0.2, 0.8, length.out = 5)
  res_f <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
                method = "dml_flex", K = 3L, S = 1L, eval_g = evg, boot = 3L)
  res_k <- suppressWarnings(
    fdid(s_cont, tr_period = 2L, ref_period = 1L, method = "kernel",
         vartype = "bootstrap", eval_g = evg, boot = 3L)
  )
  con_f <- fdid_contrast(res_f, g0 = 0.2, g1 = 0.8)
  der_f <- fdid_derivative(res_f, g0 = 0.5)
  con_k <- fdid_contrast(res_k, g0 = 0.2, g1 = 0.8)
  der_k <- fdid_derivative(res_k, g0 = 0.5)

  expect_equal(con_f$target, "level contrast")
  expect_equal(der_f$target, "derivative value")
  expect_true(is.finite(con_f$estimate))
  expect_true(is.finite(der_f$estimate))
  expect_true(is.finite(con_k$estimate))
  expect_true(is.finite(der_k$estimate))
  expect_true(con_f$inference %in% c("replicate", "vcov", "band", "pointwise_approx"))
  expect_true(con_k$inference %in% c("replicate", "vcov", "band", "pointwise_approx"))
  expect_error(fdid_contrast(res_f, 0.21, 0.8, interpolate = "none"),
               regexp = "eval_g grid")
})

test_that("[CRAN] fdid dml_incremental: returns scalar average-derivative fields", {
  res <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
              method = "dml_incremental", K = 3L, S = 1L)
  expect_s3_class(res, "fdid")
  expect_true(all(REQUIRED_FIELDS %in% names(res)))
  expect_equal(res$method, "dml_incremental")
  expect_equal(res$dml_method, "incremental")
  expect_true(is.finite(res$est$event$Estimate))
  expect_null(res$curve_event)
  expect_null(res$eval_g)
  expect_true(!is.null(res$incremental_event$signal))
  expect_true(!is.null(res$incremental_event$split_estimates))
  expect_true(!is.null(res$scalar_event$influence))
  expect_equal(res$dml_metadata$target_estimand,
               "observed_population_average_derivative")
  expect_equal(res$density_method, "residual_kde")
  expect_true(nrow(res$dynamic) > 0)
})

test_that("[CRAN] fdid flexible DML supports mapping and density extensions", {
  evg <- seq(0.25, 0.75, length.out = 4)
  res_spline <- fdid(
    s_cont, tr_period = 2L, ref_period = 1L,
    method = "dml_flex", K = 3L, S = 1L, eval_g = evg,
    signal_map = "spline", map_degree = 2L, map_df = 4L,
    density_method = "location_scale", boot = 3L
  )
  res_kernel <- fdid(
    s_cont, tr_period = 2L, ref_period = 1L,
    method = "dml_flex", K = 3L, S = 1L, eval_g = evg,
    signal_map = "kernel",
    density_method = "local_kde", boot = 3L
  )

  expect_equal(res_spline$signal_map, "spline")
  expect_equal(res_spline$density_method, "location_scale")
  expect_equal(res_kernel$signal_map, "kernel")
  expect_equal(res_kernel$density_method, "local_kde")
  expect_true(all(is.finite(res_spline$curve_event$theta_hat)))
  expect_true(all(is.finite(res_kernel$curve_event$theta_hat)))

  if (requireNamespace("mgcv", quietly = TRUE)) {
    res_gam <- fdid(
      s_cont, tr_period = 2L, ref_period = 1L,
      method = "dml_flex", K = 3L, S = 1L, eval_g = evg,
      signal_map = "gam", map_df = 4L,
      density_method = "location_scale", boot = 3L
    )
    expect_equal(res_gam$signal_map, "gam")
    expect_true(all(is.finite(res_gam$curve_event$delta_hat)))
  }

  res_incremental <- fdid(
    s_cont, tr_period = 2L, ref_period = 1L,
    method = "dml_incremental", K = 3L, S = 1L,
    density_method = "location_scale"
  )
  expect_equal(res_incremental$density_method, "location_scale")
  expect_true(is.finite(res_incremental$est$event$Estimate))
})

# ===========================================================================
# R5 — S3 dispatch on new objects
# ===========================================================================

test_that("[CRAN] print.fdid works for kernel and dml objects", {
  res_k <- fdid(s_cont, tr_period = 2L, ref_period = 1L, method = "kernel")
  res_p <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
                method = "dml_plr", K = 3L, S = 1L)
  expect_output(print(res_k))
  expect_output(print(res_p))
})

test_that("[CRAN] summary.fdid works for kernel and dml objects", {
  res_k <- fdid(s_cont, tr_period = 2L, ref_period = 1L, method = "kernel")
  res_p <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
                method = "dml_plr", K = 3L, S = 1L)
  expect_output(summary(res_k))
  expect_output(summary(res_p))
  out_k <- capture.output(summary(res_k))
  expect_true(any(grepl("Scalar event SE", out_k)))
  expect_true(any(grepl("interval average", out_k)))
})

test_that("[CRAN] summary.fdid output is ASCII-safe", {
  res_p <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
                method = "dml_plr", K = 3L, S = 1L)
  out <- capture.output(summary(res_p))
  expect_false(any(grepl("<U\\+", out)))
  expect_false(any(grepl("[\u2500\u2550\u00B7]", out, useBytes = TRUE)))
  expect_true(any(grepl("={10,}", out)))
  expect_true(any(grepl("-{10,}", out)))
})

test_that("[CRAN] plot.fdid 'raw' and 'dynamic' work on new objects", {
  res_k <- fdid(s_cont, tr_period = 2L, ref_period = 1L, method = "kernel")
  expect_invisible(plot(res_k, type = "raw"))
  expect_invisible(plot(res_k, type = "dynamic"))
})

test_that("[CRAN] plot.fdid 'curve' works for kernel object", {
  res_k <- fdid(s_cont, tr_period = 2L, ref_period = 1L, method = "kernel")
  expect_equal(res_k$alpha, 0.05)
  expect_invisible(plot(res_k, type = "curve"))
  expect_invisible(plot(res_k, type = "curve", curve = "level"))
  expect_invisible(plot(res_k, type = "curve", curve = "derivative"))
  expect_invisible(plot(res_k, type = "curve", curve = "level",
                        Xdistr = "density", xlim = c(0.2, 0.8)))
  expect_error(plot(res_k, type = "curve", xlim = c(2, 1)),
               regexp = "xlim")
})

test_that("[CRAN] plot.fdid kernel curve controls use stored inference cleanly", {
  evg <- seq(0.2, 0.8, length.out = 5)
  res_b <- suppressWarnings(
    fdid(s_cont, tr_period = 2L, ref_period = 1L, method = "kernel",
         vartype = "bootstrap", eval_g = evg, boot = 5L)
  )
  res_r <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
                method = "kernel", eval_g = evg)

  expect_invisible(plot(res_b, type = "curve", interval = "auto"))
  expect_invisible(plot(res_b, type = "curve", interval = "pointwise"))
  expect_invisible(plot(res_b, type = "curve", interval = "uniform"))
  expect_invisible(plot(res_b, type = "curve", interval = "both"))
  expect_invisible(plot(res_b, type = "curve", interval = "none"))
  expect_invisible(plot(res_b, type = "curve", show.uniform.CI = FALSE))
  expect_invisible(plot(res_b, type = "curve", Xdistr = "rug",
                        show.eval_g = TRUE, diff.values = c(0.35, 0.65)))
  expect_invisible(plot(res_b, type = "curve", Xdistr = "density",
                        support.panel = "separate"))
  expect_error(plot(res_b, type = "curve", interval = "bad"))
  expect_error(plot(res_b, type = "curve", Xdistr = "bad"))
  expect_error(plot(res_b, type = "curve", support.panel = "bad"))

  old_mfrow <- graphics::par("mfrow")
  expect_invisible(plot(res_b, type = "curve", curve = "both",
                        interval = "both"))
  expect_equal(graphics::par("mfrow"), old_mfrow)

  pd <- fdid:::.fdid_curve_plot_data(res_b, curve = "derivative",
                                     interval = "both",
                                     xlim = c(0.3, 0.7))
  expect_s3_class(pd, "data.frame")
  expect_true(all(c("g", "estimate", "pointwise_lower", "pointwise_upper",
                    "band_lower", "band_upper", "curve", "method",
                    "vartype", "pointwise_method", "band_method",
                    "band_zeta", "band_coverage", "n_boot_eff") %in% names(pd)))
  expect_true(all(pd$g >= 0.3 & pd$g <= 0.7))
  expect_equal(unique(pd$curve), "derivative")
  expect_equal(unique(pd$method), "kernel")
  expect_true(any(is.finite(pd$pointwise_lower)))
  expect_true(any(is.finite(pd$band_lower)))
  expect_equal(unique(stats::na.omit(pd$band_method)),
               "bootstrap_quantile_envelope")

  pd_none <- fdid:::.fdid_curve_plot_data(res_b, curve = "derivative",
                                          interval = "none")
  expect_true(all(is.na(pd_none$pointwise_lower)))
  expect_true(all(is.na(pd_none$band_lower)))

  pd_robust <- fdid:::.fdid_curve_plot_data(res_r, curve = "derivative",
                                            interval = "uniform")
  expect_true(all(is.na(pd_robust$band_lower)))

  pd_gap <- pd
  pd_gap$pointwise_lower[2] <- NA_real_
  seg <- fdid:::.fdid_contiguous_segments(
    pd_gap, cols = c("estimate", "pointwise_lower", "pointwise_upper")
  )
  expect_true(is.list(seg))
  expect_gt(length(seg), 1L)
})

test_that("[CRAN] plot.fdid 'curve' works for dml_flex object", {
  res_f <- fdid(s_cont, tr_period = 2L, ref_period = 1L, method = "dml_flex",
                K = 3L, S = 1L, eval_g = seq(0.2, 0.8, length.out = 5),
                boot = 3L)
  expect_equal(res_f$alpha, 0.05)
  expect_invisible(plot(res_f, type = "curve"))
  expect_invisible(plot(res_f, type = "curve", curve = "level"))
  expect_invisible(plot(res_f, type = "curve", curve = "derivative"))
  expect_invisible(plot(res_f, type = "curve", curve = "level",
                        Xdistr = "none", xlim = c(0.25, 0.75)))
})

test_that("[CRAN] dml_flex curve plots expose DML-specific interval labels", {
  evg <- seq(0.2, 0.8, length.out = 5)
  res_f <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
                method = "dml_flex", K = 3L, S = 1L,
                eval_g = evg, boot = 3L)

  expect_invisible(plot(res_f, type = "curve", interval = "auto"))
  expect_invisible(plot(res_f, type = "curve", interval = "pointwise"))
  expect_invisible(plot(res_f, type = "curve", interval = "uniform"))
  expect_invisible(plot(res_f, type = "curve", interval = "both"))
  expect_invisible(plot(res_f, type = "curve", interval = "none"))
  expect_invisible(plot(res_f, type = "curve", show.uniform.CI = FALSE))
  expect_invisible(plot(res_f, type = "curve", curve = "both",
                        interval = "both", Xdistr = "rug",
                        show.eval_g = TRUE,
                        diff.values = c(0.35, 0.65)))
  expect_invisible(plot(res_f, type = "curve", curve = "level",
                        interval = "both", Xdistr = "density",
                        support.panel = "separate",
                        contrast.values = c(0.3, 0.7),
                        xlim = c(0.25, 0.75)))

  pd <- fdid:::.fdid_curve_plot_data(
    res_f, curve = "level", interval = "both", xlim = c(0.25, 0.75)
  )
  expect_true(all(c("target_estimand", "inference_scope",
                    "inference_method", "signal_map", "density_method",
                    "pointwise_label", "band_label", "band_caveat",
                    "practical_band_method", "gaussian_crit",
                    "simultaneous_crit", "curve_vcov_available") %in%
                    names(pd)))
  expect_equal(unique(pd$target_estimand), "fixed_G_level_curve")
  expect_equal(unique(pd$signal_map), "local_poly")
  expect_equal(unique(stats::na.omit(pd$band_method)),
               "practical_signal_residual_multiplier")
  expect_equal(unique(stats::na.omit(pd$band_label)),
               "Practical DML band")
  expect_true(any(grepl("not a final theorem-level",
                        stats::na.omit(pd$band_caveat),
                        fixed = TRUE)))
  expect_false(any(pd$curve_vcov_available))

  pd_no_uniform <- fdid:::.fdid_curve_plot_data(
    res_f, curve = "level", interval = "both", show.uniform.CI = FALSE
  )
  expect_true(all(is.na(pd_no_uniform$band_lower)))
  expect_true(all(is.na(pd_no_uniform$band_label)))

  res_blp <- fdid(
    s_cont, tr_period = 2L, ref_period = 1L,
    method = "dml_flex", K = 3L, S = 1L, eval_g = evg,
    signal_map = "blp_spline", map_degree = 2L, map_df = 4L,
    boot = 3L
  )
  pd_blp <- fdid:::.fdid_curve_plot_data(
    res_blp, curve = "level", interval = "uniform"
  )
  expect_equal(unique(stats::na.omit(pd_blp$band_method)),
               "gaussian_max_t_blp_spline_covariance")
  expect_equal(unique(stats::na.omit(pd_blp$band_label)),
               "BLP max-t band")
  expect_true(any(pd_blp$curve_vcov_available))
  expect_true(any(grepl("not a DML refit bootstrap",
                        stats::na.omit(pd_blp$band_caveat),
                        fixed = TRUE)))
})

test_that("[CRAN] plot.fdid contrast works for kernel and dml_flex objects", {
  evg <- seq(0.2, 0.8, length.out = 5)
  res_k <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
                method = "kernel", eval_g = evg)
  res_f <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
                method = "dml_flex", K = 3L, S = 1L,
                eval_g = evg, boot = 3L)

  expect_invisible(plot(res_k, type = "contrast"))
  expect_invisible(plot(res_k, type = "contrast", ci = FALSE))
  expect_invisible(plot(res_f, type = "contrast"))
  expect_invisible(plot(res_f, type = "curve", curve = "level",
                        legend.position = "none"))
})

test_that("[CRAN] plot.fdid contrast helper selects and validates targets", {
  evg <- seq(0.1, 0.9, length.out = 11)
  res_k <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
                method = "kernel", eval_g = evg)

  pd_default <- fdid:::.fdid_contrast_plot_data(res_k)
  expect_lte(nrow(pd_default), 10L)
  expect_true(all(c("reference_g", "is_reference", "target_estimand") %in%
                    names(pd_default)))
  expect_equal(unique(pd_default$target_estimand), "mu(g) - mu(g_ref)")
  expect_equal(sum(pd_default$is_reference), 1L)
  ref_row <- pd_default[pd_default$is_reference, ]
  expect_equal(ref_row$estimate, 0)
  expect_equal(ref_row$std.error, 0)
  expect_equal(ref_row$conf.low, 0)
  expect_equal(ref_row$conf.high, 0)
  expect_equal(ref_row$inference, "reference")
  expect_false(any(abs(pd_default$target_g[!pd_default$is_reference] -
                         pd_default$ref_g[1L]) <= sqrt(.Machine$double.eps)))

  pd_no_ref <- fdid:::.fdid_contrast_plot_data(res_k, include.ref = FALSE)
  expect_false(any(pd_no_ref$is_reference))

  pd_contrast <- fdid:::.fdid_contrast_plot_data(
    res_k, ref.g = evg[3], contrast.values = c(evg[1], evg[5], evg[5]),
    include.ref = FALSE
  )
  expect_equal(pd_contrast$ref_g, rep(evg[3], 2L))
  expect_equal(pd_contrast$target_g, c(evg[1], evg[5]))

  pd_contrast_ref <- fdid:::.fdid_contrast_plot_data(
    res_k, ref.g = evg[3], contrast.values = c(evg[3], evg[5]),
    include.ref = FALSE
  )
  expect_equal(sum(pd_contrast_ref$is_reference), 1L)
  expect_equal(pd_contrast_ref$target_g, c(evg[3], evg[5]))

  pd_target <- fdid:::.fdid_contrast_plot_data(
    res_k, ref.g = evg[3], target.g = c(evg[9], evg[2]),
    include.ref = FALSE
  )
  expect_equal(pd_target$target_g, c(evg[2], evg[9]))
  expect_error(
    fdid:::.fdid_contrast_plot_data(
      res_k, target.g = evg[1], contrast.values = evg[2]
    ),
    regexp = "only one"
  )
  expect_warning(
    fdid:::.fdid_contrast_plot_data(
      res_k, ref.g = evg[3], target.g = evg[1:4], max.contrast = 2L
    ),
    regexp = "max.contrast"
  )
})

test_that("[CRAN] plot.fdid contrast supports reference rules and quantile targets", {
  evg <- seq(0.1, 0.9, length.out = 9)
  res_k <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
                method = "kernel", eval_g = evg)

  pd_min <- fdid:::.fdid_contrast_plot_data(res_k, ref.g = "min")
  expect_equal(unique(pd_min$reference_g), min(evg))
  expect_true(any(pd_min$is_reference))

  pd_mean <- fdid:::.fdid_contrast_plot_data(res_k, ref.g = "mean")
  expect_equal(unique(pd_mean$reference_g), evg[which.min(abs(evg - mean(evg)))])

  pd_median <- fdid:::.fdid_contrast_plot_data(res_k, ref.g = "median")
  expect_equal(unique(pd_median$reference_g),
               evg[which.min(abs(evg - stats::median(evg)))])

  expect_error(
    fdid:::.fdid_contrast_plot_data(res_k, ref.g = "bad_rule"),
    regexp = "Unknown ref.g"
  )

  pd_q <- fdid:::.fdid_contrast_plot_data(
    res_k, target.type = "quantile", target.g = 5
  )
  expect_equal(nrow(pd_q), 5L)
  expect_equal(sum(pd_q$is_reference), 1L)
  expect_true(pd_q$is_reference[1L])
  expect_true(grepl("^Q1", pd_q$label[pd_q$is_reference]))
  expect_equal(unique(pd_q$reference_g), pd_q$target_g[pd_q$is_reference])

  pd_q_short <- fdid:::.fdid_contrast_plot_data(res_k, target.g = "q5")
  expect_equal(nrow(pd_q_short), 5L)
  expect_equal(sum(pd_q_short$is_reference), 1L)

  pd_q_bare_short <- fdid:::.fdid_contrast_plot_data(res_k, target.g = "5")
  expect_equal(nrow(pd_q_bare_short), 5L)
  expect_equal(sum(pd_q_bare_short$is_reference), 1L)

  pd_q_probs <- fdid:::.fdid_contrast_plot_data(
    res_k, target.type = "quantile", target.g = c(0.25, 0.5, 0.75)
  )
  expect_equal(nrow(pd_q_probs), 3L)
  expect_true(all(grepl("^Q", pd_q_probs$label)))

  expect_error(
    fdid:::.fdid_contrast_plot_data(
      res_k, target.type = "quantile", contrast.values = evg[1:3]
    ),
    regexp = "contrast.values"
  )
})

test_that("[CRAN] plot.fdid contrast interpolation and inference choices work", {
  evg <- seq(0.2, 0.8, length.out = 5)
  res_k <- suppressWarnings(
    fdid(s_cont, tr_period = 2L, ref_period = 1L, method = "kernel",
         vartype = "bootstrap", eval_g = evg, boot = 5L)
  )

  pd_linear <- fdid:::.fdid_contrast_plot_data(
    res_k, ref.g = 0.45, target.g = 0.8, interpolate = "linear",
    include.ref = FALSE
  )
  expect_equal(nrow(pd_linear), 1L)
  expect_true(is.finite(pd_linear$estimate))
  expect_error(
    fdid:::.fdid_contrast_plot_data(
      res_k, ref.g = 0.45, target.g = 0.8, interpolate = "none",
      include.ref = FALSE
    ),
    regexp = "eval_g grid"
  )

  pd_auto <- fdid:::.fdid_contrast_plot_data(
    res_k, ref.g = evg[1], target.g = evg[5], inference = "auto",
    include.ref = FALSE
  )
  pd_rep <- fdid:::.fdid_contrast_plot_data(
    res_k, ref.g = evg[1], target.g = evg[5], inference = "replicate",
    include.ref = FALSE
  )
  pd_vcov <- fdid:::.fdid_contrast_plot_data(
    res_k, ref.g = evg[1], target.g = evg[5], inference = "vcov",
    include.ref = FALSE
  )
  expect_equal(pd_auto$inference, "vcov")
  expect_equal(pd_rep$inference, "replicate")
  expect_equal(pd_vcov$inference, "vcov")
})

test_that("[CRAN] plot.fdid contrast uses BLP vcov and rejects scalar DML", {
  evg <- seq(0.25, 0.75, length.out = 4)
  res_blp <- fdid(
    s_cont, tr_period = 2L, ref_period = 1L,
    method = "dml_flex", K = 3L, S = 1L, eval_g = evg,
    signal_map = "blp_spline", map_degree = 2L, map_df = 4L,
    boot = 3L
  )
  pd_blp <- fdid:::.fdid_contrast_plot_data(
    res_blp, ref.g = evg[1], target.g = evg[4], inference = "auto",
    include.ref = FALSE
  )
  expect_equal(pd_blp$inference, "vcov")
  expect_equal(pd_blp$target_estimand, "theta(g) - theta(g_ref)")

  res_p <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
                method = "dml_plr", K = 3L, S = 1L)
  expect_error(plot(res_p, type = "contrast"), regexp = "contrast plots")
})

test_that("[CRAN] plot.fdid 'curve' rejects scalar dml_incremental object", {
  res_i <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
                method = "dml_incremental", K = 3L, S = 1L)
  expect_error(plot(res_i, type = "curve"), regexp = "No curve data")
})

test_that("[CRAN] plot.fdid 'curve' rejects scalar DML objects", {
  res_b <- fdid(s_bin, tr_period = 2L, ref_period = 1L,
                method = "dml_binary", K = 3L, S = 1L)
  res_p <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
                method = "dml_plr", K = 3L, S = 1L)
  expect_error(plot(res_b, type = "curve"), regexp = "No curve data")
  expect_error(plot(res_p, type = "curve"), regexp = "No curve data")
})

# ===========================================================================
# R6 — Error handling
# ===========================================================================

test_that("[CRAN] fdid kernel errors on binary G", {
  d_bin2 <- data.frame(
    id = rep(1:40, each = 3), time = rep(1:3, 40),
    G  = rep(rbinom(40, 1, 0.5), each = 3),
    outcome = rnorm(120)
  )
  s_b2 <- fdid_prepare(d_bin2, "outcome", NULL, "G", "id", "time")
  expect_error(fdid(s_b2, tr_period = 2L, ref_period = 1L, method = "kernel"),
               regexp = "continuous")
})

test_that("[CRAN] fdid errors on unknown method", {
  expect_error(fdid(s_cont, tr_period = 2L, ref_period = 1L,
                    method = "dml_nonexistent"))
})

test_that("[CRAN] fdid dml_binary errors on continuous G", {
  expect_error(
    fdid(s_cont, tr_period = 2L, ref_period = 1L, method = "dml_binary"),
    regexp = "binary"
  )
})

test_that("[CRAN] fdid dml_binary requires both binary groups", {
  d_one <- data.frame(
    id = rep(1:40, each = 3), time = rep(1:3, 40),
    G  = rep(1, 120), outcome = rnorm(120)
  )
  s_one <- fdid_prepare(d_one, "outcome", NULL, "G", "id", "time")
  expect_error(
    fdid(s_one, tr_period = 2L, ref_period = 1L, method = "dml_binary"),
    regexp = "binary"
  )
})

test_that("[CRAN] fdid dml_plr, dml_flex, and dml_incremental require continuous G", {
  expect_error(
    fdid(s_bin, tr_period = 2L, ref_period = 1L,
         method = "dml_plr", K = 3L, S = 1L),
    regexp = "continuous"
  )
  expect_error(
    fdid(s_bin, tr_period = 2L, ref_period = 1L,
         method = "dml_flex", K = 3L, S = 1L),
    regexp = "continuous"
  )
  expect_error(
    fdid(s_bin, tr_period = 2L, ref_period = 1L,
         method = "dml_incremental", K = 3L, S = 1L),
    regexp = "continuous"
  )
})

test_that("[CRAN] continuous-G methods reject unsupported variance types", {
  expect_error(
    fdid(s_cont, tr_period = 2L, ref_period = 1L,
         method = "kernel", vartype = "bad"),
    regexp = "vartype"
  )
  expect_error(
    fdid(s_cont, tr_period = 2L, ref_period = 1L,
         method = "kernel", vartype = "jackknife"),
    regexp = "kernel.*supports"
  )
  expect_error(
    fdid(s_cont, tr_period = 2L, ref_period = 1L,
         method = "dml_plr", vartype = "bootstrap"),
    regexp = "DML methods.*robust"
  )
})

test_that("[CRAN] flexible DML rejects unknown mapping and density options", {
  expect_error(
    fdid(s_cont, tr_period = 2L, ref_period = 1L,
         method = "dml_flex", K = 3L, S = 1L,
         signal_map = "bad_map"),
    regexp = "signal_map"
  )
  expect_error(
    fdid(s_cont, tr_period = 2L, ref_period = 1L,
         method = "dml_flex", K = 3L, S = 1L,
         density_method = "bad_density"),
    regexp = "density_method"
  )
})

test_that("[CRAN] scalar DML rejects unknown inference options", {
  expect_error(
    fdid(s_cont, tr_period = 2L, ref_period = 1L,
         method = "dml_plr", dml_inference = "bad_inference"),
    regexp = "dml_inference"
  )
  expect_error(
    fdid(s_cont, tr_period = 2L, ref_period = 1L,
         method = "dml_plr", dml_multiplier = "bad_multiplier"),
    regexp = "dml_multiplier"
  )
  expect_error(
    fdid(s_cont, tr_period = 2L, ref_period = 1L,
         method = "dml_plr", dml_boot = 0),
    regexp = "dml_boot"
  )
})

test_that("[CRAN] continuous-G tuning validation fails clearly", {
  expect_error(
    fdid(s_cont, tr_period = 2L, ref_period = 1L,
         method = "kernel", alpha = 1),
    regexp = "alpha"
  )
  expect_error(
    fdid(s_cont, tr_period = 2L, ref_period = 1L,
         method = "kernel", trim = 0.5),
    regexp = "trim"
  )
  expect_error(
    fdid(s_cont, tr_period = 2L, ref_period = 1L,
         method = "kernel", eval_g = c(0.2, Inf)),
    regexp = "eval_g"
  )
  expect_error(
    fdid(s_cont, tr_period = 2L, ref_period = 1L,
         method = "kernel", h0 = -1),
    regexp = "h0"
  )
  expect_error(
    fdid(s_cont, tr_period = 2L, ref_period = 1L,
         method = "kernel", K_folds = 1),
    regexp = "K_folds"
  )
  expect_error(
    fdid(s_cont, tr_period = 2L, ref_period = 1L,
         method = "dml_plr", K = 1),
    regexp = "K"
  )
  expect_error(
    fdid(s_cont, tr_period = 2L, ref_period = 1L,
         method = "dml_plr", S = 0),
    regexp = "S"
  )
  expect_error(
    fdid(s_cont, tr_period = 2L, ref_period = 1L,
         method = "dml_flex", boot = 0),
    regexp = "boot"
  )
})

test_that("[CRAN] DML stores diagnostics and cluster-aware scalar paths", {
  set.seed(224)
  n_cl <- 50
  d_cl <- data.frame(
    id = rep(seq_len(n_cl), each = 3),
    time = rep(1:3, n_cl),
    cluster = rep(rep(seq_len(10), each = 5), each = 3),
    G = rep(stats::runif(n_cl), each = 3),
    outcome = stats::rnorm(n_cl * 3),
    covar1 = rep(stats::rnorm(n_cl), each = 3)
  )
  s_cl <- fdid_prepare(d_cl, "outcome", "covar1", "G", "id", "time",
                       cluster_label = "cluster")
  res_plr <- fdid(s_cl, tr_period = 2L, ref_period = 1L,
                  method = "dml_plr", K = 3L, S = 1L)
  res_inc <- fdid(s_cl, tr_period = 2L, ref_period = 1L,
                  method = "dml_incremental", K = 3L, S = 1L)
  res_flex <- fdid(s_cl, tr_period = 2L, ref_period = 1L,
                   method = "dml_flex", K = 3L, S = 1L,
                   eval_g = seq(0.2, 0.8, length.out = 4), boot = 3L)

  expect_equal(length(unique(res_plr$cluster)), 10L)
  expect_equal(res_inc$dml_diagnostics$n_clusters, 10L)
  expect_true(is.finite(res_inc$est$event$Std.Error))
  expect_true(!is.null(res_inc$dml_diagnostics$d_mu_quantiles))
  expect_true(!is.null(res_flex$dml_diagnostics$density_quantiles))
  expect_true(!is.null(res_flex$curve_event$theta_band_lower))
  expect_equal(length(res_flex$curve_event$theta_band_lower), length(res_flex$eval_g))
  expect_output(summary(res_flex), regexp = "Simultaneous bands")
  expect_output(summary(res_inc), regexp = "DML Diagnostics")
})

# ===========================================================================
# R7 — entire_period restriction (skip_on_cran for speed)
# ===========================================================================

test_that("fdid kernel: entire_period restricts dynamic output", {
  skip_on_cran()
  res <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
              method = "kernel", entire_period = c(1L, 2L))
  expect_equal(sort(res$entire_period), c(1L, 2L))
  expect_equal(nrow(res$dynamic), 2L)
})

test_that("fdid dml_plr: entire_period restricts dynamic output", {
  skip_on_cran()
  res <- fdid(s_cont, tr_period = 2L, ref_period = 1L,
              method = "dml_plr", K = 3L, S = 1L,
              entire_period = c(1L, 2L))
  expect_equal(sort(res$entire_period), c(1L, 2L))
  expect_equal(nrow(res$dynamic), 2L)
})

# ===========================================================================
# R2 — Existing ols1 estimate unchanged (regression guard)
# ===========================================================================

test_that("R2: ols1 estimate unchanged after update", {
  skip_on_cran()
  data(fdid)
  mortality$uniqueid <- paste(mortality$provid, mortality$countyid, sep = "-")
  mortality$G <- ifelse(mortality$pczupu >= median(mortality$pczupu, na.rm = TRUE), 1, 0)
  s_mort <- fdid_prepare(
    data = mortality, Y_label = "mortality",
    X_labels = c("avggrain", "lnpop"),
    G_label = "G", unit_label = "uniqueid", time_label = "year"
  )
  res <- fdid(s_mort, tr_period = 1958, ref_period = 1957,
              method = "ols1", vartype = "robust")
  expect_equal(res$est$event$Estimate, -1.024863, tolerance = 1e-5)
  expect_equal(res$est$event$Std.Error,  0.3705082, tolerance = 1e-5)
})
