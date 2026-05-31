# Nuisance-oracle ablation for continuous-G DML signals.
#
# Usage from the package root:
#   LC_ALL=C Rscript inst/nuisance_oracle_ablation.R
#
# The script uses the stress-test DGPs with known mu(g,x), s(g|x),
# partial_g mu(g,x), and partial_g log s(g|x). It keeps cross-fitting for
# estimated nuisance components and replaces outcome and density nuisances one
# block at a time:
#
#   1. oracle_all
#   2. estimated_outcome_oracle_density
#   3. oracle_outcome_estimated_density
#   4. estimated_all

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(fdid)
  }
})

load_stress_environment <- function() {
  stress_file <- file.path("inst", "oracle_signal_stress_tests.R")
  exprs <- parse(stress_file)
  env <- new.env(parent = globalenv())
  if (length(exprs) > 1L) {
    for (ii in seq_len(length(exprs) - 1L)) {
      eval(exprs[[ii]], envir = env)
    }
  }
  env
}

ablation_config <- function() {
  list(
    master_seed = as.integer(Sys.getenv("FDID_ABLATION_SEED", "20260505")),
    reps = as.integer(Sys.getenv("FDID_ABLATION_REPS", "24")),
    n = as.integer(strsplit(Sys.getenv("FDID_ABLATION_N", "600,1000"), ",")[[1]]),
    learner = Sys.getenv(
      "FDID_ABLATION_LEARNER",
      if (requireNamespace("mgcv", quietly = TRUE)) "gam" else "linear"
    ),
    signal_map = Sys.getenv("FDID_ABLATION_SIGNAL_MAP", "spline"),
    map_degree = as.integer(Sys.getenv("FDID_ABLATION_MAP_DEGREE", "3")),
    map_df = as.integer(Sys.getenv("FDID_ABLATION_MAP_DF", "7")),
    density_methods = strsplit(
      Sys.getenv("FDID_ABLATION_DENSITIES", "residual_kde,location_scale"),
      ","
    )[[1]],
    K = as.integer(Sys.getenv("FDID_ABLATION_K", "3")),
    S = as.integer(Sys.getenv("FDID_ABLATION_S", "1")),
    out_dir = Sys.getenv("FDID_ABLATION_OUT", "inst/nuisance_oracle_ablation_results")
  )
}

theta_truth <- function(g, d) {
  vapply(g, function(g0) mean(d$mu_fn(g0, d$x1, d$x2)), numeric(1L))
}

delta_truth <- function(g, d) {
  vapply(g, function(g0) mean(d$dmu_fn(g0, d$x1, d$x2)), numeric(1L))
}

omega_truth <- function(g, d) {
  vapply(g, function(g0) mean(d$density(g0, d$x1, d$x2)), numeric(1L))
}

observed_avg_derivative <- function(d) {
  mean(d$dmu)
}

curve_metrics <- function(est, se, truth) {
  z <- stats::qnorm(0.975)
  ok <- is.finite(est) & is.finite(truth)
  if (!any(ok)) {
    return(c(rmse = NA_real_, bias = NA_real_, coverage = NA_real_))
  }
  coverage <- if (!is.null(se) && any(is.finite(se[ok]))) {
    mean(est[ok] - z * se[ok] <= truth[ok] &
           est[ok] + z * se[ok] >= truth[ok],
         na.rm = TRUE)
  } else {
    NA_real_
  }
  c(
    rmse = sqrt(mean((est[ok] - truth[ok])^2)),
    bias = mean(est[ok] - truth[ok]),
    coverage = coverage
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

finite_kde_fn <- function(x) {
  safe_bw <- getFromNamespace(".safe_kde_bw", "fdid")
  kd <- stats::density(x, bw = safe_bw(x), n = 1024)
  stats::approxfun(kd, rule = 2)
}

predict_mu_average_loo <- function(mu_model, G_ho, X_ho) {
  predict_reg_model <- getFromNamespace(".predict_reg_model", "fdid")
  n_ho <- length(G_ho)
  if (n_ho <= 1L) {
    return(predict_reg_model(mu_model, cbind(G_ho, X_ho)))
  }

  if (identical(mu_model$type, "lm")) {
    co <- stats::coef(mu_model$model)
    b0 <- if (!is.na(co["(Intercept)"])) co["(Intercept)"] else 0
    bG <- if (!is.na(co["v1"])) co["v1"] else 0
    p_x <- if (is.null(mu_model$p_x)) ncol(X_ho) else mu_model$p_x
    bX_names <- if (p_x > 0L) paste0("v", seq.int(2L, p_x + 1L)) else character(0)
    bX <- if (p_x > 0L) co[bX_names] else numeric(0)
    bX[is.na(bX)] <- 0
    X_sum_ho <- if (ncol(X_ho) > 0L) colSums(X_ho) else numeric(0)
    mu_sum <- n_ho * b0 + n_ho * bG * G_ho +
      if (length(bX) > 0L) as.numeric(bX %*% X_sum_ho) else 0
    mu_own <- b0 + bG * G_ho +
      if (ncol(X_ho) > 0L) as.numeric(X_ho %*% bX) else 0
    return((mu_sum - mu_own) / (n_ho - 1L))
  }

  vapply(seq_len(n_ho), function(ii) {
    mate_idx <- setdiff(seq_len(n_ho), ii)
    X_mates <- if (ncol(X_ho) > 0L) {
      X_ho[mate_idx, , drop = FALSE]
    } else {
      matrix(0, length(mate_idx), 0)
    }
    feat_i <- cbind(rep(G_ho[ii], length(mate_idx)), X_mates)
    mean(predict_reg_model(mu_model, feat_i), na.rm = TRUE)
  }, numeric(1L))
}

layer_grid <- function(density_methods) {
  rbind(
    data.frame(
      layer = c("oracle_all", "estimated_outcome_oracle_density"),
      outcome_mode = c("oracle", "estimated"),
      density_mode = "oracle",
      density_method = "oracle",
      stringsAsFactors = FALSE
    ),
    do.call(rbind, lapply(density_methods, function(dm) {
      data.frame(
        layer = c("oracle_outcome_estimated_density", "estimated_all"),
        outcome_mode = c("oracle", "estimated"),
        density_mode = "estimated",
        density_method = dm,
        stringsAsFactors = FALSE
      )
    }))
  )
}

estimate_ablation_layer <- function(d, eval_g, cfg, outcome_mode,
                                    density_mode, density_method) {
  fit_reg_model <- getFromNamespace(".fit_reg_model", "fdid")
  predict_reg_model <- getFromNamespace(".predict_reg_model", "fdid")
  fit_cond_density <- getFromNamespace(".fit_cond_density", "fdid")
  predict_cond_density <- getFromNamespace(".predict_cond_density", "fdid")
  predict_mu_derivative <- getFromNamespace(".predict_mu_g_derivative", "fdid")
  predict_log_density_derivative <- getFromNamespace(".predict_log_density_g_derivative", "fdid")
  map_signal <- getFromNamespace(".map_signal_to_g", "fdid")

  n <- d$n
  X <- cbind(d$x1, d$x2)
  all_theta <- matrix(NA_real_, cfg$S, length(eval_g))
  all_delta <- matrix(NA_real_, cfg$S, length(eval_g))
  all_se_theta <- matrix(NA_real_, cfg$S, length(eval_g))
  all_se_delta <- matrix(NA_real_, cfg$S, length(eval_g))
  scalar_s <- numeric(cfg$S)
  var_s <- numeric(cfg$S)

  for (ss in seq_len(cfg$S)) {
    fold_id <- sample(rep(seq_len(cfg$K), length.out = n))
    xi <- numeric(n)
    phi <- numeric(n)

    for (kk in seq_len(cfg$K)) {
      tr <- which(fold_id != kk)
      ho <- which(fold_id == kk)
      X_tr <- X[tr, , drop = FALSE]
      X_ho <- X[ho, , drop = FALSE]

      if (identical(outcome_mode, "oracle")) {
        mu_hat <- d$mu[ho]
        d_mu_hat <- d$dmu[ho]
        m_hat <- theta_truth(d$G[ho], d)
      } else {
        feat_tr <- cbind(d$G[tr], X_tr)
        mu_model <- fit_reg_model(d$DeltaY[tr], feat_tr, cfg$learner)
        mu_model$p_x <- ncol(X_tr)
        feat_ho <- cbind(d$G[ho], X_ho)
        mu_hat <- predict_reg_model(mu_model, feat_ho)
        d_mu_hat <- predict_mu_derivative(mu_model, d$G[ho], X_ho, d$G[tr])
        m_hat <- predict_mu_average_loo(mu_model, d$G[ho], X_ho)
      }

      if (identical(density_mode, "oracle")) {
        s_hat <- d$s_cond[ho]
        d_log_s_hat <- d$dlog_s[ho]
        omega_hat <- omega_truth(d$G[ho], d)
      } else {
        s_model <- fit_cond_density(d$G[tr], X_tr, cfg$learner, density_method)
        s_hat <- predict_cond_density(s_model, d$G[ho], X_ho)
        d_log_s_hat <- predict_log_density_derivative(s_model, d$G[ho], X_ho, d$G[tr])
        omega_fn <- finite_kde_fn(d$G[tr])
        omega_hat <- as.numeric(omega_fn(d$G[ho]))
      }

      s_hat <- pmax(s_hat, 1e-8)
      omega_hat <- pmax(omega_hat, 1e-8)
      xi[ho] <- ((d$DeltaY[ho] - mu_hat) / s_hat) * omega_hat + m_hat
      phi[ho] <- d_mu_hat - d_log_s_hat * (d$DeltaY[ho] - mu_hat)
    }

    mapped <- map_signal(
      xi_hat = xi,
      G = d$G,
      eval_g = eval_g,
      signal_map = cfg$signal_map,
      map_degree = cfg$map_degree,
      map_df = cfg$map_df
    )
    all_theta[ss, ] <- mapped$theta
    all_delta[ss, ] <- mapped$delta
    all_se_theta[ss, ] <- mapped$se_theta
    all_se_delta[ss, ] <- mapped$se_delta
    scalar_s[ss] <- mean(phi, na.rm = TRUE)
    var_s[ss] <- mean((phi - scalar_s[ss])^2, na.rm = TRUE) / n
  }

  theta <- apply(all_theta, 2, stats::median, na.rm = TRUE)
  delta <- apply(all_delta, 2, stats::median, na.rm = TRUE)
  theta_dev <- sweep(all_theta, 2, theta, "-")
  delta_dev <- sweep(all_delta, 2, delta, "-")
  se_theta <- sqrt(apply(all_se_theta^2 + theta_dev^2, 2, stats::median, na.rm = TRUE))
  se_delta <- sqrt(apply(all_se_delta^2 + delta_dev^2, 2, stats::median, na.rm = TRUE))
  scalar <- stats::median(scalar_s, na.rm = TRUE)
  scalar_se <- sqrt(max(0, stats::median(var_s + (scalar_s - scalar)^2, na.rm = TRUE)))

  list(
    theta = theta,
    delta = delta,
    se_theta = se_theta,
    se_delta = se_delta,
    scalar = scalar,
    scalar_se = scalar_se
  )
}

run_ablation_rep <- function(scenario, n, seed, eval_g, cfg) {
  set.seed(seed)
  d <- scenario$draw(n, seed)
  theta_true <- theta_truth(eval_g, d)
  delta_true <- delta_truth(eval_g, d)
  scalar_true <- observed_avg_derivative(d)
  layers <- layer_grid(cfg$density_methods)
  fits <- vector("list", nrow(layers))

  for (ii in seq_len(nrow(layers))) {
    fits[[ii]] <- tryCatch(
      estimate_ablation_layer(
        d = d,
        eval_g = eval_g,
        cfg = cfg,
        outcome_mode = layers$outcome_mode[ii],
        density_mode = layers$density_mode[ii],
        density_method = layers$density_method[ii]
      ),
      error = function(e) e
    )
  }

  oracle_idx <- which(layers$layer == "oracle_all")[1L]
  oracle_fit <- fits[[oracle_idx]]
  rows <- lapply(seq_len(nrow(layers)), function(ii) {
    fit <- fits[[ii]]
    if (inherits(fit, "error")) {
      return(data.frame(
        scenario = scenario$name,
        n = n,
        seed = seed,
        layer = layers$layer[ii],
        density_method = layers$density_method[ii],
        level_rmse = NA_real_,
        level_coverage = NA_real_,
        derivative_rmse = NA_real_,
        derivative_coverage = NA_real_,
        scalar_error = NA_real_,
        scalar_abs_error = NA_real_,
        scalar_covered = NA_real_,
        level_vs_oracle_rmse = NA_real_,
        derivative_vs_oracle_rmse = NA_real_,
        scalar_vs_oracle_error = NA_real_,
        error = conditionMessage(fit),
        stringsAsFactors = FALSE
      ))
    }

    level <- curve_metrics(fit$theta, fit$se_theta, theta_true)
    deriv <- curve_metrics(fit$delta, fit$se_delta, delta_true)
    scalar <- scalar_metrics(fit$scalar, fit$scalar_se, scalar_true)
    if (inherits(oracle_fit, "error")) {
      level_vs_oracle <- NA_real_
      deriv_vs_oracle <- NA_real_
      scalar_vs_oracle <- NA_real_
    } else {
      level_vs_oracle <- sqrt(mean((fit$theta - oracle_fit$theta)^2, na.rm = TRUE))
      deriv_vs_oracle <- sqrt(mean((fit$delta - oracle_fit$delta)^2, na.rm = TRUE))
      scalar_vs_oracle <- fit$scalar - oracle_fit$scalar
    }

    data.frame(
      scenario = scenario$name,
      n = n,
      seed = seed,
      layer = layers$layer[ii],
      density_method = layers$density_method[ii],
      level_rmse = level[["rmse"]],
      level_coverage = level[["coverage"]],
      derivative_rmse = deriv[["rmse"]],
      derivative_coverage = deriv[["coverage"]],
      scalar_error = scalar[["error"]],
      scalar_abs_error = scalar[["abs_error"]],
      scalar_covered = scalar[["covered"]],
      level_vs_oracle_rmse = level_vs_oracle,
      derivative_vs_oracle_rmse = deriv_vs_oracle,
      scalar_vs_oracle_error = scalar_vs_oracle,
      error = "",
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

summarise_ablation <- function(rows) {
  key <- interaction(rows$scenario, rows$n, rows$layer, rows$density_method, drop = TRUE)
  out <- do.call(rbind, lapply(split(rows, key), function(d) {
    data.frame(
      scenario = unique(d$scenario),
      n = unique(d$n),
      layer = unique(d$layer),
      density_method = unique(d$density_method),
      reps = nrow(d),
      failures = sum(!is.na(d$error) & nzchar(d$error)),
      level_rmse = mean(d$level_rmse, na.rm = TRUE),
      level_coverage = mean(d$level_coverage, na.rm = TRUE),
      derivative_rmse = mean(d$derivative_rmse, na.rm = TRUE),
      derivative_coverage = mean(d$derivative_coverage, na.rm = TRUE),
      scalar_bias = mean(d$scalar_error, na.rm = TRUE),
      scalar_rmse = sqrt(mean(d$scalar_error^2, na.rm = TRUE)),
      scalar_coverage = mean(d$scalar_covered, na.rm = TRUE),
      level_vs_oracle_rmse = mean(d$level_vs_oracle_rmse, na.rm = TRUE),
      derivative_vs_oracle_rmse = mean(d$derivative_vs_oracle_rmse, na.rm = TRUE),
      scalar_vs_oracle_rmse = sqrt(mean(d$scalar_vs_oracle_error^2, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }))
  out[order(out$scenario, out$n, out$layer, out$density_method), ]
}

format_number <- function(x) {
  ifelse(is.na(x), NA_character_, sprintf("%.4f", x))
}

markdown_table <- function(x) {
  x2 <- x
  for (nm in names(x2)) {
    if (is.numeric(x2[[nm]]) &&
        !nm %in% c("n", "reps", "failures")) {
      x2[[nm]] <- format_number(x2[[nm]])
    }
  }
  header <- paste("|", paste(names(x2), collapse = " | "), "|")
  sep <- paste("|", paste(rep("---", ncol(x2)), collapse = " | "), "|")
  body <- apply(x2, 1, function(row) paste("|", paste(row, collapse = " | "), "|"))
  paste(c(header, sep, body), collapse = "\n")
}

write_ablation_report <- function(cfg, scenario_info, rows, summary) {
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(rows, file.path(cfg$out_dir, "nuisance_oracle_ablation_replications.csv"),
                   row.names = FALSE)
  utils::write.csv(summary, file.path(cfg$out_dir, "nuisance_oracle_ablation_summary.csv"),
                   row.names = FALSE)

  scenario_lines <- unlist(lapply(scenario_info, function(s) {
    c(
      paste0("- `", s$name, "`: ", s$description),
      paste0("  Eval grid: `", paste(sprintf("%.3f", s$eval_g), collapse = ", "), "`.")
    )
  }))

  lines <- c(
    "# Nuisance-Oracle Ablation Results",
    "",
    sprintf("Master seed: `%d`.", cfg$master_seed),
    sprintf("Replications per scenario/n: `%d`.", cfg$reps),
    sprintf("Sample sizes: `%s`.", paste(cfg$n, collapse = ", ")),
    sprintf("Learner for estimated nuisances: `%s`.", cfg$learner),
    sprintf("Signal map: `%s`, map df: `%d`.", cfg$signal_map, cfg$map_df),
    sprintf("Estimated density methods: `%s`.", paste(cfg$density_methods, collapse = "`, `")),
    "",
    "## Stress Scenarios",
    "",
    scenario_lines,
    "",
    "## Ablation Layers",
    "",
    "- `oracle_all`: true outcome and true density nuisances.",
    "- `estimated_outcome_oracle_density`: estimated `mu` and `partial_g mu`, true density nuisances.",
    "- `oracle_outcome_estimated_density`: true outcome nuisances, estimated density and density derivative.",
    "- `estimated_all`: estimated outcome and estimated density nuisances.",
    "",
    "## Summary",
    "",
    markdown_table(summary),
    "",
    "## Interpretation",
    "",
    "- `*_vs_oracle_*` columns compare each ablation layer with `oracle_all` on the same replication.",
    "- Large `oracle_outcome_estimated_density` errors isolate density and density-derivative problems.",
    "- Large `estimated_outcome_oracle_density` errors isolate outcome-regression and outcome-derivative problems.",
    "- `dml_flex` curve diagnostics are in the level and derivative columns.",
    "- `dml_incremental` diagnostics are in the scalar columns."
  )

  path <- file.path(cfg$out_dir, "NUISANCE_ORACLE_ABLATION_RESULTS.md")
  writeLines(lines, path)
  invisible(path)
}

run_nuisance_oracle_ablation <- function(cfg = ablation_config()) {
  stress_env <- load_stress_environment()
  scenarios <- get("stress_scenarios", envir = stress_env)()
  scenario_eval_grid <- get("scenario_eval_grid", envir = stress_env)

  scenario_info <- lapply(seq_along(scenarios), function(ii) {
    list(
      name = scenarios[[ii]]$name,
      description = scenarios[[ii]]$description,
      eval_g = scenario_eval_grid(scenarios[[ii]], cfg$master_seed + ii * 1000L)
    )
  })

  message("Running nuisance-oracle ablation")
  message("  Scenarios: ", paste(vapply(scenarios, `[[`, character(1L), "name"),
                                  collapse = ", "))
  message("  Reps per scenario/n: ", cfg$reps)
  message("  Sample sizes: ", paste(cfg$n, collapse = ", "))
  message("  Learner: ", cfg$learner)
  message("  Density methods: ", paste(cfg$density_methods, collapse = ", "))

  rows <- do.call(rbind, lapply(seq_along(scenarios), function(si) {
    scenario <- scenarios[[si]]
    eval_g <- scenario_info[[si]]$eval_g
    do.call(rbind, lapply(cfg$n, function(n) {
      do.call(rbind, lapply(seq_len(cfg$reps), function(rr) {
        if (rr %% 6L == 0L) {
          message("  ", scenario$name, " n=", n, " rep=", rr)
        }
        seed <- cfg$master_seed + si * 10000000L + n * 1000L + rr
        run_ablation_rep(scenario, n, seed, eval_g, cfg)
      }))
    }))
  }))

  summary <- summarise_ablation(rows)
  report <- write_ablation_report(cfg, scenario_info, rows, summary)
  print(summary)
  message("Wrote report: ", report)
  invisible(list(rows = rows, summary = summary, report = report))
}

if (!interactive()) {
  run_nuisance_oracle_ablation()
}
