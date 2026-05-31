# Stress tests for continuous-G DML oracle signals.
#
# Usage from the package root:
#   LC_ALL=C Rscript inst/oracle_signal_stress_tests.R
#
# This script extends inst/oracle_signal_validation.R by using DGPs with
# heteroskedastic, mixture, and heavy-tailed G | X density shapes plus nonlinear,
# interaction-heavy, and nonsmooth mu(g, x). Each DGP provides known
# mu(g,x), s(g|x), partial_g mu(g,x), and partial_g log s(g|x).

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(fdid)
  }
})

stress_config <- function() {
  list(
    master_seed = as.integer(Sys.getenv("FDID_STRESS_SEED", "20260504")),
    oracle_reps = as.integer(Sys.getenv("FDID_STRESS_ORACLE_REPS", "80")),
    estimated_reps = as.integer(Sys.getenv("FDID_STRESS_EST_REPS", "12")),
    oracle_n = as.integer(strsplit(Sys.getenv("FDID_STRESS_ORACLE_N", "600,1200"), ",")[[1]]),
    estimated_n = as.integer(strsplit(Sys.getenv("FDID_STRESS_EST_N", "500,1000"), ",")[[1]]),
    learner = Sys.getenv(
      "FDID_STRESS_LEARNER",
      if (requireNamespace("mgcv", quietly = TRUE)) "gam" else "linear"
    ),
    signal_map = Sys.getenv("FDID_STRESS_SIGNAL_MAP", "spline"),
    map_degree = as.integer(Sys.getenv("FDID_STRESS_MAP_DEGREE", "3")),
    map_df = as.integer(Sys.getenv("FDID_STRESS_MAP_DF", "7")),
    density_methods = strsplit(
      Sys.getenv("FDID_STRESS_DENSITIES", "residual_kde,location_scale"),
      ","
    )[[1]],
    K = as.integer(Sys.getenv("FDID_STRESS_K", "3")),
    S = as.integer(Sys.getenv("FDID_STRESS_S", "1")),
    out_dir = Sys.getenv("FDID_STRESS_OUT", "inst/oracle_signal_stress_results")
  )
}

sign0 <- function(x) ifelse(x > 0, 1, ifelse(x < 0, -1, 0))

scenario_heteroskedastic <- function() {
  sigma_y <- 0.45
  list(
    name = "heteroskedastic_interactions",
    description = paste(
      "G | X is normal with X-dependent scale;",
      "mu(g,x) has smooth G-by-X interactions."
    ),
    draw = function(n, seed) {
      set.seed(seed)
      x1 <- stats::rnorm(n)
      x2 <- stats::rnorm(n)
      m <- 0.35 * x1 - 0.25 * x2
      sigma <- exp(-0.35 + 0.25 * x1)
      sigma <- pmin(pmax(sigma, 0.30), 1.40)
      g <- m + sigma * stats::rnorm(n)
      base <- 0.30 + 0.45 * x1 - 0.20 * x2
      mu <- base + 0.55 * g + 0.20 * g^2 + sin(0.90 * g) +
        0.35 * g * x1 - 0.20 * g^2 * x2 + 0.18 * sin(g * x1)
      dmu <- 0.55 + 0.40 * g + 0.90 * cos(0.90 * g) +
        0.35 * x1 - 0.40 * g * x2 + 0.18 * x1 * cos(g * x1)
      eps <- stats::rnorm(n, sd = sigma_y)
      alpha <- stats::rnorm(n)
      density <- function(g0, x1_new, x2_new) {
        m0 <- 0.35 * x1_new - 0.25 * x2_new
        s0 <- exp(-0.35 + 0.25 * x1_new)
        s0 <- pmin(pmax(s0, 0.30), 1.40)
        stats::dnorm((g0 - m0) / s0) / s0
      }
      dlog <- function(g0, x1_new, x2_new) {
        m0 <- 0.35 * x1_new - 0.25 * x2_new
        s0 <- exp(-0.35 + 0.25 * x1_new)
        s0 <- pmin(pmax(s0, 0.30), 1.40)
        -(g0 - m0) / s0^2
      }
      mu_fn <- function(g0, x1_new, x2_new) {
        base0 <- 0.30 + 0.45 * x1_new - 0.20 * x2_new
        base0 + 0.55 * g0 + 0.20 * g0^2 + sin(0.90 * g0) +
          0.35 * g0 * x1_new - 0.20 * g0^2 * x2_new +
          0.18 * sin(g0 * x1_new)
      }
      dmu_fn <- function(g0, x1_new, x2_new) {
        0.55 + 0.40 * g0 + 0.90 * cos(0.90 * g0) +
          0.35 * x1_new - 0.40 * g0 * x2_new +
          0.18 * x1_new * cos(g0 * x1_new)
      }
      list(
        scenario = "heteroskedastic_interactions",
        n = n,
        x1 = x1,
        x2 = x2,
        G = g,
        mu = mu,
        DeltaY = mu + eps,
        s_cond = density(g, x1, x2),
        dlog_s = dlog(g, x1, x2),
        dmu = dmu,
        alpha = alpha,
        density = density,
        dlog = dlog,
        mu_fn = mu_fn,
        dmu_fn = dmu_fn
      )
    }
  )
}

scenario_mixture_kink <- function() {
  sigma_y <- 0.45
  list(
    name = "mixture_kink",
    description = paste(
      "G | X is a two-component normal mixture with X-dependent weights;",
      "mu(g,x) has kinks through abs(g) and positive-part terms."
    ),
    draw = function(n, seed) {
      set.seed(seed)
      x1 <- stats::rnorm(n)
      x2 <- stats::rnorm(n)
      w <- stats::plogis(0.80 * x1 - 0.50 * x2)
      m_base <- 0.30 * x1 - 0.20 * x2
      m1 <- m_base - 0.55
      m2 <- m_base + 0.65
      s1 <- 0.45
      s2 <- 0.90
      comp <- stats::rbinom(n, 1L, w)
      g <- ifelse(comp == 1L, m1 + s1 * stats::rnorm(n),
                  m2 + s2 * stats::rnorm(n))
      sp <- sign0(g)
      base <- 0.20 + 0.30 * x1 - 0.25 * x2
      mu <- base + 0.55 * g + 0.35 * abs(g) + 0.25 * pmax(g, 0)^2 +
        0.30 * g * x1 - 0.20 * abs(g) * x2 + 0.20 * sin(1.50 * g * x1)
      dmu <- 0.55 + 0.35 * sp + 0.50 * pmax(g, 0) +
        0.30 * x1 - 0.20 * sp * x2 + 0.30 * x1 * cos(1.50 * g * x1)
      eps <- stats::rnorm(n, sd = sigma_y)
      alpha <- stats::rnorm(n)
      density <- function(g0, x1_new, x2_new) {
        w0 <- stats::plogis(0.80 * x1_new - 0.50 * x2_new)
        mb <- 0.30 * x1_new - 0.20 * x2_new
        d1 <- stats::dnorm((g0 - (mb - 0.55)) / s1) / s1
        d2 <- stats::dnorm((g0 - (mb + 0.65)) / s2) / s2
        w0 * d1 + (1 - w0) * d2
      }
      dlog <- function(g0, x1_new, x2_new) {
        w0 <- stats::plogis(0.80 * x1_new - 0.50 * x2_new)
        mb <- 0.30 * x1_new - 0.20 * x2_new
        mu1 <- mb - 0.55
        mu2 <- mb + 0.65
        d1 <- stats::dnorm((g0 - mu1) / s1) / s1
        d2 <- stats::dnorm((g0 - mu2) / s2) / s2
        s <- w0 * d1 + (1 - w0) * d2
        sp <- w0 * d1 * (-(g0 - mu1) / s1^2) +
          (1 - w0) * d2 * (-(g0 - mu2) / s2^2)
        sp / pmax(s, 1e-12)
      }
      mu_fn <- function(g0, x1_new, x2_new) {
        sp0 <- sign0(g0)
        base0 <- 0.20 + 0.30 * x1_new - 0.25 * x2_new
        base0 + 0.55 * g0 + 0.35 * abs(g0) + 0.25 * pmax(g0, 0)^2 +
          0.30 * g0 * x1_new - 0.20 * abs(g0) * x2_new +
          0.20 * sin(1.50 * g0 * x1_new)
      }
      dmu_fn <- function(g0, x1_new, x2_new) {
        sp0 <- sign0(g0)
        0.55 + 0.35 * sp0 + 0.50 * pmax(g0, 0) +
          0.30 * x1_new - 0.20 * sp0 * x2_new +
          0.30 * x1_new * cos(1.50 * g0 * x1_new)
      }
      list(
        scenario = "mixture_kink",
        n = n,
        x1 = x1,
        x2 = x2,
        G = g,
        mu = mu,
        DeltaY = mu + eps,
        s_cond = density(g, x1, x2),
        dlog_s = dlog(g, x1, x2),
        dmu = dmu,
        alpha = alpha,
        density = density,
        dlog = dlog,
        mu_fn = mu_fn,
        dmu_fn = dmu_fn
      )
    }
  )
}

scenario_heavy_tail <- function() {
  sigma_y <- 0.45
  df <- 3
  list(
    name = "heavy_tail_interactions",
    description = paste(
      "G | X has scaled t(3) residuals and X-dependent scale;",
      "mu(g,x) has nonlinear trigonometric interactions."
    ),
    draw = function(n, seed) {
      set.seed(seed)
      x1 <- stats::rnorm(n)
      x2 <- stats::rnorm(n)
      m <- 0.25 * x1 + 0.25 * x2
      sigma <- 0.55 + 0.25 * stats::plogis(x1 - x2)
      z <- stats::rt(n, df = df)
      g <- m + sigma * z
      base <- -0.10 + 0.35 * x1 + 0.20 * x2
      mu <- base + 0.45 * g + 0.12 * g^3 / (1 + g^2) +
        0.25 * sin(1.20 * g) + 0.25 * g * x1 * x2 +
        0.18 * cos(g * x2)
      dmu <- 0.45 + 0.12 * (3 * g^2 + g^4) / (1 + g^2)^2 +
        0.30 * cos(1.20 * g) + 0.25 * x1 * x2 -
        0.18 * x2 * sin(g * x2)
      eps <- stats::rnorm(n, sd = sigma_y)
      alpha <- stats::rnorm(n)
      density <- function(g0, x1_new, x2_new) {
        m0 <- 0.25 * x1_new + 0.25 * x2_new
        s0 <- 0.55 + 0.25 * stats::plogis(x1_new - x2_new)
        stats::dt((g0 - m0) / s0, df = df) / s0
      }
      dlog <- function(g0, x1_new, x2_new) {
        m0 <- 0.25 * x1_new + 0.25 * x2_new
        s0 <- 0.55 + 0.25 * stats::plogis(x1_new - x2_new)
        z0 <- (g0 - m0) / s0
        -((df + 1) * z0 / (df + z0^2)) / s0
      }
      mu_fn <- function(g0, x1_new, x2_new) {
        base0 <- -0.10 + 0.35 * x1_new + 0.20 * x2_new
        base0 + 0.45 * g0 + 0.12 * g0^3 / (1 + g0^2) +
          0.25 * sin(1.20 * g0) + 0.25 * g0 * x1_new * x2_new +
          0.18 * cos(g0 * x2_new)
      }
      dmu_fn <- function(g0, x1_new, x2_new) {
        0.45 + 0.12 * (3 * g0^2 + g0^4) / (1 + g0^2)^2 +
          0.30 * cos(1.20 * g0) + 0.25 * x1_new * x2_new -
          0.18 * x2_new * sin(g0 * x2_new)
      }
      list(
        scenario = "heavy_tail_interactions",
        n = n,
        x1 = x1,
        x2 = x2,
        G = g,
        mu = mu,
        DeltaY = mu + eps,
        s_cond = density(g, x1, x2),
        dlog_s = dlog(g, x1, x2),
        dmu = dmu,
        alpha = alpha,
        density = density,
        dlog = dlog,
        mu_fn = mu_fn,
        dmu_fn = dmu_fn
      )
    }
  )
}

stress_scenarios <- function() {
  list(
    scenario_heteroskedastic(),
    scenario_mixture_kink(),
    scenario_heavy_tail()
  )
}

theta_sample <- function(g, d) {
  vapply(g, function(g0) mean(d$mu_fn(g0, d$x1, d$x2)), numeric(1L))
}

delta_sample <- function(g, d) {
  vapply(g, function(g0) mean(d$dmu_fn(g0, d$x1, d$x2)), numeric(1L))
}

omega_empirical <- function(g, d) {
  vapply(g, function(g0) mean(d$density(g0, d$x1, d$x2)), numeric(1L))
}

observed_avg_derivative <- function(d) {
  mean(d$dmu)
}

scenario_eval_grid <- function(scenario, seed) {
  pilot <- scenario$draw(25000L, seed)
  as.numeric(stats::quantile(pilot$G, probs = seq(0.15, 0.85, length.out = 9)))
}

oracle_signals <- function(d, eval_g, cfg) {
  map_signal <- getFromNamespace(".map_signal_to_g", "fdid")
  xi <- ((d$DeltaY - d$mu) / d$s_cond) * omega_empirical(d$G, d) +
    theta_sample(d$G, d)
  phi <- d$dmu - d$dlog_s * (d$DeltaY - d$mu)
  mapped <- map_signal(
    xi_hat = xi,
    G = d$G,
    eval_g = eval_g,
    signal_map = cfg$signal_map,
    map_degree = cfg$map_degree,
    map_df = cfg$map_df
  )
  list(
    theta_hat = mapped$theta,
    delta_hat = mapped$delta,
    se_theta = mapped$se_theta,
    se_delta = mapped$se_delta,
    scalar = mean(phi),
    scalar_se = stats::sd(phi) / sqrt(length(phi))
  )
}

make_panel_s <- function(d) {
  panel <- data.frame(
    id = rep(seq_len(d$n), each = 2L),
    time = rep(1:2, times = d$n),
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
  c(
    rmse = sqrt(mean((est[ok] - truth[ok])^2)),
    bias = mean(est[ok] - truth[ok]),
    coverage = mean(est[ok] - z * se[ok] <= truth[ok] &
                      est[ok] + z * se[ok] >= truth[ok],
                    na.rm = TRUE)
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

run_oracle_rep <- function(scenario, n, seed, cfg, eval_g) {
  d <- scenario$draw(n, seed)
  o <- oracle_signals(d, eval_g, cfg)
  theta_true <- theta_sample(eval_g, d)
  delta_true <- delta_sample(eval_g, d)
  scalar_true <- observed_avg_derivative(d)
  scalar <- scalar_metrics(o$scalar, o$scalar_se, scalar_true)
  level <- curve_metrics(o$theta_hat, o$se_theta, theta_true)
  deriv <- curve_metrics(o$delta_hat, o$se_delta, delta_true)
  data.frame(
    scenario = scenario$name,
    n = n,
    seed = seed,
    level_rmse = level[["rmse"]],
    level_coverage = level[["coverage"]],
    derivative_rmse = deriv[["rmse"]],
    derivative_coverage = deriv[["coverage"]],
    scalar_error = scalar[["error"]],
    scalar_abs_error = scalar[["abs_error"]],
    scalar_covered = scalar[["covered"]],
    stringsAsFactors = FALSE
  )
}

run_estimated_rep <- function(scenario, n, seed, density_method, cfg, eval_g) {
  d <- scenario$draw(n, seed)
  s <- make_panel_s(d)
  o <- oracle_signals(d, eval_g, cfg)
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
      density_method = density_method,
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
      density_method = density_method,
      K = cfg$K,
      S = cfg$S
    ),
    error = function(e) e
  )

  if (inherits(flex, "error")) {
    flex_level <- c(rmse = NA_real_, coverage = NA_real_)
    flex_deriv <- c(rmse = NA_real_, coverage = NA_real_)
    level_vs_oracle <- NA_real_
    deriv_vs_oracle <- NA_real_
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
    level_vs_oracle <- sqrt(mean(
      (flex$curve_event$theta_hat - o$theta_hat)^2,
      na.rm = TRUE
    ))
    deriv_vs_oracle <- sqrt(mean(
      (flex$curve_event$delta_hat - o$delta_hat)^2,
      na.rm = TRUE
    ))
    flex_error <- ""
  }

  o_scalar <- scalar_metrics(o$scalar, o$scalar_se, scalar_true)
  if (inherits(inc, "error")) {
    inc_scalar <- c(error = NA_real_, abs_error = NA_real_, covered = NA_real_)
    inc_vs_oracle <- NA_real_
    inc_error <- conditionMessage(inc)
  } else {
    inc_scalar <- scalar_metrics(
      inc$est$event$Estimate,
      inc$est$event$Std.Error,
      scalar_true
    )
    inc_vs_oracle <- inc$est$event$Estimate - o$scalar
    inc_error <- ""
  }

  data.frame(
    scenario = scenario$name,
    density_method = density_method,
    n = n,
    seed = seed,
    oracle_level_rmse = curve_metrics(o$theta_hat, o$se_theta, theta_true)[["rmse"]],
    flex_level_rmse = flex_level[["rmse"]],
    flex_level_vs_oracle_rmse = level_vs_oracle,
    flex_level_coverage = flex_level[["coverage"]],
    oracle_derivative_rmse = curve_metrics(o$delta_hat, o$se_delta, delta_true)[["rmse"]],
    flex_derivative_rmse = flex_deriv[["rmse"]],
    flex_derivative_vs_oracle_rmse = deriv_vs_oracle,
    flex_derivative_coverage = flex_deriv[["coverage"]],
    oracle_scalar_error = o_scalar[["error"]],
    incremental_scalar_error = inc_scalar[["error"]],
    incremental_scalar_abs_error = inc_scalar[["abs_error"]],
    incremental_scalar_covered = inc_scalar[["covered"]],
    incremental_vs_oracle = inc_vs_oracle,
    flex_error = flex_error,
    incremental_error = inc_error,
    stringsAsFactors = FALSE
  )
}

summarise_oracle <- function(rows) {
  key <- interaction(rows$scenario, rows$n, drop = TRUE)
  out <- do.call(rbind, lapply(split(rows, key), function(d) {
    data.frame(
      scenario = unique(d$scenario),
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
  out[order(out$scenario, out$n), ]
}

summarise_estimated <- function(rows) {
  key <- interaction(rows$scenario, rows$density_method, rows$n, drop = TRUE)
  out <- do.call(rbind, lapply(split(rows, key), function(d) {
    data.frame(
      scenario = unique(d$scenario),
      density_method = unique(d$density_method),
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
      incremental_scalar_bias = mean(d$incremental_scalar_error, na.rm = TRUE),
      incremental_scalar_rmse = sqrt(mean(d$incremental_scalar_error^2, na.rm = TRUE)),
      incremental_scalar_coverage = mean(d$incremental_scalar_covered, na.rm = TRUE),
      incremental_vs_oracle_rmse = sqrt(mean(d$incremental_vs_oracle^2, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }))
  out[order(out$scenario, out$density_method, out$n), ]
}

fmt <- function(x) ifelse(is.na(x), NA_character_, sprintf("%.4f", x))

markdown_table <- function(x) {
  x2 <- x
  for (nm in names(x2)) {
    if (is.numeric(x2[[nm]]) &&
        !nm %in% c("n", "reps", "flex_failures", "incremental_failures")) {
      x2[[nm]] <- fmt(x2[[nm]])
    }
  }
  header <- paste("|", paste(names(x2), collapse = " | "), "|")
  sep <- paste("|", paste(rep("---", ncol(x2)), collapse = " | "), "|")
  body <- apply(x2, 1, function(row) paste("|", paste(row, collapse = " | "), "|"))
  paste(c(header, sep, body), collapse = "\n")
}

write_report <- function(cfg, scenario_info, oracle_summary, estimated_summary,
                         oracle_rows, estimated_rows) {
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(oracle_rows, file.path(cfg$out_dir, "stress_oracle_replications.csv"),
                   row.names = FALSE)
  utils::write.csv(estimated_rows, file.path(cfg$out_dir, "stress_estimated_replications.csv"),
                   row.names = FALSE)
  utils::write.csv(oracle_summary, file.path(cfg$out_dir, "stress_oracle_summary.csv"),
                   row.names = FALSE)
  utils::write.csv(estimated_summary, file.path(cfg$out_dir, "stress_estimated_summary.csv"),
                   row.names = FALSE)

  scenario_lines <- unlist(lapply(scenario_info, function(s) {
    c(
      paste0("- `", s$name, "`: ", s$description),
      paste0("  Eval grid: `", paste(sprintf("%.3f", s$eval_g), collapse = ", "), "`.")
    )
  }))

  lines <- c(
    "# Oracle Signal Stress-Test Results",
    "",
    sprintf("Master seed: `%d`.", cfg$master_seed),
    sprintf("Oracle reps per scenario/n: `%d`.", cfg$oracle_reps),
    sprintf("Estimated reps per scenario/n/density: `%d`.", cfg$estimated_reps),
    sprintf("Estimated learner: `%s`.", cfg$learner),
    sprintf("Signal map: `%s`, map df: `%d`.", cfg$signal_map, cfg$map_df),
    sprintf("Density methods: `%s`.", paste(cfg$density_methods, collapse = "`, `")),
    "",
    "## Stress Scenarios",
    "",
    scenario_lines,
    "",
    "## Oracle Signal Results",
    "",
    markdown_table(oracle_summary),
    "",
    "## Estimated-Nuisance Results",
    "",
    markdown_table(estimated_summary),
    "",
    "## Reading The Table",
    "",
    "- `flex_*` columns evaluate `method = \"dml_flex\"` level and derivative curves.",
    "- `incremental_*` columns evaluate the scalar `method = \"dml_incremental\"` target.",
    "- `*_vs_oracle_*` columns compare estimated-nuisance results with oracle-signal results on the same DGP.",
    "- These are stress tests; higher RMSE than the smooth baseline is expected."
  )

  path <- file.path(cfg$out_dir, "ORACLE_SIGNAL_STRESS_RESULTS.md")
  writeLines(lines, path)
  invisible(path)
}

run_oracle_signal_stress_tests <- function(cfg = stress_config()) {
  scenarios <- stress_scenarios()
  scenario_info <- lapply(seq_along(scenarios), function(i) {
    evg <- scenario_eval_grid(scenarios[[i]], cfg$master_seed + i * 1000L)
    list(
      name = scenarios[[i]]$name,
      description = scenarios[[i]]$description,
      eval_g = evg
    )
  })

  message("Running oracle signal stress tests")
  message("  Scenarios: ", paste(vapply(scenarios, `[[`, character(1L), "name"),
                                  collapse = ", "))
  message("  Oracle reps: ", cfg$oracle_reps)
  message("  Estimated reps: ", cfg$estimated_reps)
  message("  Density methods: ", paste(cfg$density_methods, collapse = ", "))

  oracle_rows <- do.call(rbind, lapply(seq_along(scenarios), function(si) {
    scenario <- scenarios[[si]]
    eval_g <- scenario_info[[si]]$eval_g
    do.call(rbind, lapply(cfg$oracle_n, function(n) {
      do.call(rbind, lapply(seq_len(cfg$oracle_reps), function(r) {
        if (r %% 20L == 0L) {
          message("  oracle ", scenario$name, " n=", n, " rep=", r)
        }
        run_oracle_rep(
          scenario, n,
          cfg$master_seed + si * 10000000L + n * 1000L + r,
          cfg, eval_g
        )
      }))
    }))
  }))

  estimated_rows <- do.call(rbind, lapply(seq_along(scenarios), function(si) {
    scenario <- scenarios[[si]]
    eval_g <- scenario_info[[si]]$eval_g
    do.call(rbind, lapply(cfg$density_methods, function(dm) {
      do.call(rbind, lapply(cfg$estimated_n, function(n) {
        do.call(rbind, lapply(seq_len(cfg$estimated_reps), function(r) {
          if (r %% 6L == 0L) {
            message("  estimated ", scenario$name, " density=", dm,
                    " n=", n, " rep=", r)
          }
          run_estimated_rep(
            scenario, n,
            cfg$master_seed + si * 20000000L + n * 2000L +
              match(dm, cfg$density_methods) * 10000L + r,
            dm, cfg, eval_g
          )
        }))
      }))
    }))
  }))

  oracle_summary <- summarise_oracle(oracle_rows)
  estimated_summary <- summarise_estimated(estimated_rows)
  report <- write_report(cfg, scenario_info, oracle_summary, estimated_summary,
                         oracle_rows, estimated_rows)
  print(oracle_summary)
  print(estimated_summary)
  message("Wrote report: ", report)
  invisible(list(
    oracle_summary = oracle_summary,
    estimated_summary = estimated_summary,
    report = report
  ))
}

if (!interactive()) {
  run_oracle_signal_stress_tests()
}
