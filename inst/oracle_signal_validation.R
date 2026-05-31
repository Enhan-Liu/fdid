# Oracle-signal validation for continuous-G DML estimators.
#
# Usage from the package root:
#   LC_ALL=C Rscript inst/oracle_signal_validation.R
#
# The script simulates a DGP where mu(g,x), s(g|x), partial_g mu(g,x),
# and partial_g log s(g|x) are known. It first validates the oracle signals,
# then compares the implemented estimated-nuisance fdid estimators against the
# oracle baselines on the same DGP.

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(fdid)
  }
})

`%||%` <- function(x, y) if (is.null(x)) y else x

validation_config <- function() {
  list(
    master_seed = as.integer(Sys.getenv("FDID_ORACLE_SEED", "20260503")),
    oracle_reps = as.integer(Sys.getenv("FDID_ORACLE_REPS", "120")),
    estimated_reps = as.integer(Sys.getenv("FDID_EST_REPS", "40")),
    oracle_n = as.integer(strsplit(Sys.getenv("FDID_ORACLE_N", "500,1500"), ",")[[1]]),
    estimated_n = as.integer(strsplit(Sys.getenv("FDID_EST_N", "400,800"), ",")[[1]]),
    learner = Sys.getenv(
      "FDID_ORACLE_LEARNER",
      if (requireNamespace("mgcv", quietly = TRUE)) "gam" else "linear"
    ),
    signal_map = Sys.getenv("FDID_ORACLE_SIGNAL_MAP", "spline"),
    map_degree = as.integer(Sys.getenv("FDID_ORACLE_MAP_DEGREE", "3")),
    map_df = as.integer(Sys.getenv("FDID_ORACLE_MAP_DF", "6")),
    density_method = Sys.getenv("FDID_ORACLE_DENSITY", "residual_kde"),
    K = as.integer(Sys.getenv("FDID_ORACLE_K", "3")),
    S = as.integer(Sys.getenv("FDID_ORACLE_S", "1")),
    out_dir = Sys.getenv("FDID_ORACLE_OUT", "inst/oracle_signal_validation_results")
  )
}

dgp_constants <- function() {
  list(
    sigma_g = 0.70,
    sigma_y = 0.45,
    sd_g_marginal = sqrt(0.40^2 + (-0.30)^2 + 0.70^2)
  )
}

f_g <- function(g) {
  0.80 * g + 0.25 * g^2 + sin(0.70 * g)
}

df_g <- function(g) {
  0.80 + 0.50 * g + 0.70 * cos(0.70 * g)
}

draw_oracle_dgp <- function(n, seed) {
  cc <- dgp_constants()
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  m_g <- 0.40 * x1 - 0.30 * x2
  g <- m_g + cc$sigma_g * stats::rnorm(n)
  base_x <- 0.40 + 0.50 * x1 - 0.25 * x2
  mu <- base_x + f_g(g)
  eps <- stats::rnorm(n, sd = cc$sigma_y)
  delta_y <- mu + eps
  s_cond <- stats::dnorm((g - m_g) / cc$sigma_g) / cc$sigma_g
  omega <- stats::dnorm(g / cc$sd_g_marginal) / cc$sd_g_marginal
  dlog_s <- -(g - m_g) / cc$sigma_g^2
  dmu <- df_g(g)
  alpha <- stats::rnorm(n)

  list(
    n = n,
    x1 = x1,
    x2 = x2,
    G = g,
    base_x = base_x,
    mu = mu,
    DeltaY = delta_y,
    s_cond = s_cond,
    omega = omega,
    dlog_s = dlog_s,
    dmu = dmu,
    alpha = alpha
  )
}

theta_sample <- function(g, d) {
  mean(d$base_x) + f_g(g)
}

delta_sample <- function(g, d) {
  df_g(g)
}

observed_avg_derivative <- function(d) {
  mean(d$dmu)
}

oracle_eval_grid <- function() {
  cc <- dgp_constants()
  as.numeric(stats::qnorm(seq(0.15, 0.85, length.out = 9),
                          mean = 0, sd = cc$sd_g_marginal))
}

oracle_signals <- function(d, eval_g, signal_map, map_degree, map_df) {
  map_signal <- getFromNamespace(".map_signal_to_g", "fdid")
  xi <- ((d$DeltaY - d$mu) / d$s_cond) * d$omega + theta_sample(d$G, d)
  phi <- d$dmu - d$dlog_s * (d$DeltaY - d$mu)
  mapped <- map_signal(
    xi_hat = xi,
    G = d$G,
    eval_g = eval_g,
    signal_map = signal_map,
    map_degree = map_degree,
    map_df = map_df
  )

  list(
    xi = xi,
    phi = phi,
    theta_hat = mapped$theta,
    delta_hat = mapped$delta,
    se_theta = mapped$se_theta,
    se_delta = mapped$se_delta,
    scalar = mean(phi),
    scalar_se = stats::sd(phi) / sqrt(length(phi))
  )
}

make_panel_s <- function(d) {
  n <- d$n
  panel <- data.frame(
    id = rep(seq_len(n), each = 2L),
    time = rep(1:2, times = n),
    G = rep(d$G, each = 2L),
    x1 = rep(d$x1, each = 2L),
    x2 = rep(d$x2, each = 2L),
    Y = as.vector(t(cbind(d$alpha, d$alpha + d$DeltaY)))
  )
  fdid_prepare(
    data = panel,
    Y_label = "Y",
    X_labels = c("x1", "x2"),
    G_label = "G",
    unit_label = "id",
    time_label = "time"
  )
}

curve_metrics <- function(est, se, truth) {
  z <- stats::qnorm(0.975)
  ok <- is.finite(est) & is.finite(truth)
  if (!any(ok)) {
    return(c(rmse = NA_real_, bias = NA_real_, coverage = NA_real_))
  }
  cover <- if (!is.null(se) && any(is.finite(se))) {
    mean(est[ok] - z * se[ok] <= truth[ok] &
           est[ok] + z * se[ok] >= truth[ok], na.rm = TRUE)
  } else {
    NA_real_
  }
  c(
    rmse = sqrt(mean((est[ok] - truth[ok])^2)),
    bias = mean(est[ok] - truth[ok]),
    coverage = cover
  )
}

scalar_metrics <- function(est, se, truth) {
  z <- stats::qnorm(0.975)
  c(
    error = est - truth,
    abs_error = abs(est - truth),
    covered = as.numeric(est - z * se <= truth && est + z * se >= truth)
  )
}

run_oracle_rep <- function(n, seed, cfg, eval_g) {
  d <- draw_oracle_dgp(n, seed)
  o <- oracle_signals(d, eval_g, cfg$signal_map, cfg$map_degree, cfg$map_df)
  theta_true <- theta_sample(eval_g, d)
  delta_true <- delta_sample(eval_g, d)
  scalar_true <- observed_avg_derivative(d)
  scalar <- scalar_metrics(o$scalar, o$scalar_se, scalar_true)

  data.frame(
    n = n,
    seed = seed,
    level_rmse = curve_metrics(o$theta_hat, o$se_theta, theta_true)[["rmse"]],
    level_coverage = curve_metrics(o$theta_hat, o$se_theta, theta_true)[["coverage"]],
    derivative_rmse = curve_metrics(o$delta_hat, o$se_delta, delta_true)[["rmse"]],
    derivative_coverage = curve_metrics(o$delta_hat, o$se_delta, delta_true)[["coverage"]],
    scalar_error = scalar[["error"]],
    scalar_abs_error = scalar[["abs_error"]],
    scalar_covered = scalar[["covered"]],
    stringsAsFactors = FALSE
  )
}

run_estimated_rep <- function(n, seed, cfg, eval_g) {
  d <- draw_oracle_dgp(n, seed)
  s <- make_panel_s(d)
  o <- oracle_signals(d, eval_g, cfg$signal_map, cfg$map_degree, cfg$map_df)
  theta_true <- theta_sample(eval_g, d)
  delta_true <- delta_sample(eval_g, d)
  scalar_true <- observed_avg_derivative(d)

  flex <- tryCatch(
    fdid(
      s = s,
      tr_period = 2L,
      ref_period = 1L,
      entire_period = 1L,
      method = "dml_flex",
      learner = cfg$learner,
      signal_map = cfg$signal_map,
      map_degree = cfg$map_degree,
      map_df = cfg$map_df,
      density_method = cfg$density_method,
      K = cfg$K,
      S = cfg$S,
      eval_g = eval_g
    ),
    error = function(e) e
  )

  inc <- tryCatch(
    fdid(
      s = s,
      tr_period = 2L,
      ref_period = 1L,
      entire_period = 1L,
      method = "dml_incremental",
      learner = cfg$learner,
      density_method = cfg$density_method,
      K = cfg$K,
      S = cfg$S
    ),
    error = function(e) e
  )

  if (inherits(flex, "error")) {
    flex_level <- c(rmse = NA_real_, bias = NA_real_, coverage = NA_real_)
    flex_deriv <- c(rmse = NA_real_, bias = NA_real_, coverage = NA_real_)
    flex_level_oracle_rmse <- NA_real_
    flex_deriv_oracle_rmse <- NA_real_
    flex_error <- conditionMessage(flex)
  } else {
    flex_level <- curve_metrics(
      flex$curve_event$theta_hat,
      flex$curve_event$se_theta,
      theta_true
    )
    flex_deriv <- curve_metrics(
      flex$curve_event$delta_hat,
      flex$curve_event$se_delta,
      delta_true
    )
    flex_level_oracle_rmse <- sqrt(mean(
      (flex$curve_event$theta_hat - o$theta_hat)^2,
      na.rm = TRUE
    ))
    flex_deriv_oracle_rmse <- sqrt(mean(
      (flex$curve_event$delta_hat - o$delta_hat)^2,
      na.rm = TRUE
    ))
    flex_error <- ""
  }

  o_scalar_metrics <- scalar_metrics(o$scalar, o$scalar_se, scalar_true)

  if (inherits(inc, "error")) {
    inc_metrics <- c(error = NA_real_, abs_error = NA_real_, covered = NA_real_)
    inc_minus_oracle <- NA_real_
    inc_se <- NA_real_
    inc_error <- conditionMessage(inc)
  } else {
    inc_est <- inc$est$event$Estimate
    inc_se <- inc$est$event$Std.Error
    inc_metrics <- scalar_metrics(inc_est, inc_se, scalar_true)
    inc_minus_oracle <- inc_est - o$scalar
    inc_error <- ""
  }

  data.frame(
    n = n,
    seed = seed,
    oracle_level_rmse = curve_metrics(o$theta_hat, o$se_theta, theta_true)[["rmse"]],
    oracle_derivative_rmse = curve_metrics(o$delta_hat, o$se_delta, delta_true)[["rmse"]],
    flex_level_rmse = flex_level[["rmse"]],
    flex_level_coverage = flex_level[["coverage"]],
    flex_level_vs_oracle_rmse = flex_level_oracle_rmse,
    flex_derivative_rmse = flex_deriv[["rmse"]],
    flex_derivative_coverage = flex_deriv[["coverage"]],
    flex_derivative_vs_oracle_rmse = flex_deriv_oracle_rmse,
    oracle_scalar_error = o_scalar_metrics[["error"]],
    oracle_scalar_abs_error = o_scalar_metrics[["abs_error"]],
    oracle_scalar_covered = o_scalar_metrics[["covered"]],
    incremental_scalar_error = inc_metrics[["error"]],
    incremental_scalar_abs_error = inc_metrics[["abs_error"]],
    incremental_scalar_covered = inc_metrics[["covered"]],
    incremental_minus_oracle = inc_minus_oracle,
    incremental_se = inc_se,
    flex_error = flex_error,
    incremental_error = inc_error,
    stringsAsFactors = FALSE
  )
}

summarise_oracle <- function(rows) {
  do.call(rbind, lapply(split(rows, rows$n), function(d) {
    data.frame(
      n = unique(d$n),
      reps = nrow(d),
      level_rmse = mean(d$level_rmse, na.rm = TRUE),
      level_coverage = mean(d$level_coverage, na.rm = TRUE),
      derivative_rmse = mean(d$derivative_rmse, na.rm = TRUE),
      derivative_coverage = mean(d$derivative_coverage, na.rm = TRUE),
      scalar_bias = mean(d$scalar_error, na.rm = TRUE),
      scalar_rmse = sqrt(mean(d$scalar_error^2, na.rm = TRUE)),
      scalar_coverage = mean(d$scalar_covered, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

summarise_estimated <- function(rows) {
  do.call(rbind, lapply(split(rows, rows$n), function(d) {
    data.frame(
      n = unique(d$n),
      reps = nrow(d),
      flex_failures = sum(nzchar(d$flex_error)),
      incremental_failures = sum(nzchar(d$incremental_error)),
      oracle_level_rmse = mean(d$oracle_level_rmse, na.rm = TRUE),
      flex_level_rmse = mean(d$flex_level_rmse, na.rm = TRUE),
      flex_level_vs_oracle_rmse = mean(d$flex_level_vs_oracle_rmse, na.rm = TRUE),
      flex_level_coverage = mean(d$flex_level_coverage, na.rm = TRUE),
      oracle_derivative_rmse = mean(d$oracle_derivative_rmse, na.rm = TRUE),
      flex_derivative_rmse = mean(d$flex_derivative_rmse, na.rm = TRUE),
      flex_derivative_vs_oracle_rmse = mean(d$flex_derivative_vs_oracle_rmse, na.rm = TRUE),
      flex_derivative_coverage = mean(d$flex_derivative_coverage, na.rm = TRUE),
      oracle_scalar_bias = mean(d$oracle_scalar_error, na.rm = TRUE),
      oracle_scalar_rmse = sqrt(mean(d$oracle_scalar_error^2, na.rm = TRUE)),
      oracle_scalar_coverage = mean(d$oracle_scalar_covered, na.rm = TRUE),
      incremental_scalar_bias = mean(d$incremental_scalar_error, na.rm = TRUE),
      incremental_scalar_rmse = sqrt(mean(d$incremental_scalar_error^2, na.rm = TRUE)),
      incremental_scalar_coverage = mean(d$incremental_scalar_covered, na.rm = TRUE),
      incremental_vs_oracle_bias = mean(d$incremental_minus_oracle, na.rm = TRUE),
      incremental_vs_oracle_rmse = sqrt(mean(d$incremental_minus_oracle^2, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }))
}

format_number <- function(x) {
  ifelse(is.na(x), NA_character_, sprintf("%.4f", x))
}

markdown_table <- function(x) {
  x2 <- x
  for (nm in names(x2)) {
    if (is.numeric(x2[[nm]]) && !nm %in% c("n", "reps", "flex_failures",
                                           "incremental_failures")) {
      x2[[nm]] <- format_number(x2[[nm]])
    }
  }
  header <- paste("|", paste(names(x2), collapse = " | "), "|")
  sep <- paste("|", paste(rep("---", ncol(x2)), collapse = " | "), "|")
  body <- apply(x2, 1, function(row) {
    paste("|", paste(row, collapse = " | "), "|")
  })
  paste(c(header, sep, body), collapse = "\n")
}

write_report <- function(cfg, eval_g, oracle_summary, estimated_summary,
                         oracle_rows, estimated_rows) {
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(oracle_rows, file.path(cfg$out_dir, "oracle_replications.csv"),
                   row.names = FALSE)
  utils::write.csv(estimated_rows, file.path(cfg$out_dir, "estimated_replications.csv"),
                   row.names = FALSE)
  utils::write.csv(oracle_summary, file.path(cfg$out_dir, "oracle_summary.csv"),
                   row.names = FALSE)
  utils::write.csv(estimated_summary, file.path(cfg$out_dir, "estimated_summary.csv"),
                   row.names = FALSE)

  lines <- c(
    "# Oracle Signal Validation Results",
    "",
    "This validation uses a DGP where the following nuisance functions are known:",
    "",
    "- `mu(g, x) = 0.40 + 0.50 x1 - 0.25 x2 + 0.80 g + 0.25 g^2 + sin(0.70 g)`;",
    "- `G | X ~ N(0.40 x1 - 0.30 x2, 0.70^2)`;",
    "- `partial_g mu(g, x) = 0.80 + 0.50 g + 0.70 cos(0.70 g)`;",
    "- `partial_g log s(g | x) = -(g - 0.40 x1 + 0.30 x2) / 0.70^2`.",
    "",
    sprintf("Master seed: `%d`.", cfg$master_seed),
    sprintf("Evaluation grid: `%s`.", paste(sprintf("%.3f", eval_g), collapse = ", ")),
    sprintf("Estimated-nuisance learner: `%s`.", cfg$learner),
    sprintf("Signal map: `%s`, density method: `%s`.", cfg$signal_map, cfg$density_method),
    "",
    "## Oracle Signal Sanity Check",
    "",
    "The oracle rows use true `mu`, true `s(g|x)`, true `partial_g mu`, and true `partial_g log s`.",
    "",
    markdown_table(oracle_summary),
    "",
    "## Estimated Nuisance Comparison",
    "",
    "These rows compare the implemented `fdid` DML estimators to the oracle baselines on the same DGP.",
    "",
    markdown_table(estimated_summary),
    "",
    "## Interpretation",
    "",
    "- Oracle signal RMSE should fall as `n` increases.",
    "- Estimated-nuisance RMSE should be larger than oracle RMSE, because nuisance estimation adds error.",
    "- The scalar `dml_incremental` comparison uses the correct target `E[partial_g mu(G, X)]`.",
    "- Derivative-curve comparisons belong to `dml_flex`, not `dml_incremental`."
  )

  writeLines(lines, file.path(cfg$out_dir, "ORACLE_SIGNAL_VALIDATION_RESULTS.md"))
  invisible(file.path(cfg$out_dir, "ORACLE_SIGNAL_VALIDATION_RESULTS.md"))
}

run_oracle_signal_validation <- function(cfg = validation_config()) {
  eval_g <- oracle_eval_grid()
  message("Running oracle signal validation")
  message("  Oracle reps per n: ", cfg$oracle_reps)
  message("  Estimated reps per n: ", cfg$estimated_reps)
  message("  Learner: ", cfg$learner)

  oracle_rows <- do.call(rbind, lapply(cfg$oracle_n, function(n) {
    do.call(rbind, lapply(seq_len(cfg$oracle_reps), function(r) {
      if (r %% 25L == 0L) message("  oracle n=", n, " rep=", r)
      run_oracle_rep(n, cfg$master_seed + n * 10000L + r, cfg, eval_g)
    }))
  }))

  estimated_rows <- do.call(rbind, lapply(cfg$estimated_n, function(n) {
    do.call(rbind, lapply(seq_len(cfg$estimated_reps), function(r) {
      if (r %% 10L == 0L) message("  estimated n=", n, " rep=", r)
      run_estimated_rep(n, cfg$master_seed + n * 20000L + r, cfg, eval_g)
    }))
  }))

  oracle_summary <- summarise_oracle(oracle_rows)
  estimated_summary <- summarise_estimated(estimated_rows)
  report_path <- write_report(cfg, eval_g, oracle_summary, estimated_summary,
                              oracle_rows, estimated_rows)
  print(oracle_summary)
  print(estimated_summary)
  message("Wrote report: ", report_path)
  invisible(list(
    oracle_summary = oracle_summary,
    estimated_summary = estimated_summary,
    report = report_path
  ))
}

if (!interactive()) {
  run_oracle_signal_validation()
}
