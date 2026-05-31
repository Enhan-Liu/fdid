# Monte Carlo simulation study for continuous-G FDID estimators
# DGPs from dgp_spec.md; run after devtools::load_all('.')
# Usage: source("inst/simulate_mc.R")

suppressPackageStartupMessages(devtools::load_all("."))

# ── Helper: build wide s from synthetic panel ─────────────────────────────────

make_s <- function(n_units, n_pre, n_post, G_vec, X_mat, beta_fn,
                   ref_period = 0, sigma = 0.5, seed = 42) {
  set.seed(seed)
  times   <- seq(-n_pre, n_post)
  long_df <- do.call(rbind, lapply(seq_len(n_units), function(i) {
    alpha_i <- rnorm(1)
    data.frame(
      id     = i,
      time   = times,
      G      = G_vec[i],
      Y      = alpha_i + ifelse(times > ref_period,
                                 beta_fn(G_vec[i], X_mat[i, , drop = FALSE]) + rnorm(1, 0, sigma),
                                 rnorm(1, 0, sigma)),
      covar1 = X_mat[i, 1],
      stringsAsFactors = FALSE
    )
  }))
  fdid_prepare(data = long_df, Y_label = "Y", X_labels = "covar1",
               G_label = "G", unit_label = "id", time_label = "time")
}

# ── DGP 1: Linear moderation, constant incremental effect ─────────────────────
#   ΔY = 2G + 0.5X + ε   →   δ_DID(g) = 2, θ(g) = 2g, β_PLR = 2

cat("\n══════════════════════════════════════════════════════\n")
cat("DGP 1 — Linear moderation (n=300, 100 simulations)\n")
cat("══════════════════════════════════════════════════════\n")

set.seed(2025); n <- 300
G1 <- runif(n); X1 <- matrix(rnorm(n), n, 1)
beta_fn_1 <- function(g, x) 2 * g + 0.5 * x[1]
TRUE_delta <- 2; TRUE_beta  <- 2
eval_g1    <- seq(0.1, 0.9, length.out = 9)

n_sim <- 100
bias_plr    <- numeric(n_sim)
cov_plr     <- numeric(n_sim)
cov_binary  <- numeric(n_sim)
rmse_kernel <- numeric(n_sim)
cov_kernel  <- numeric(n_sim)
rmse_flex   <- numeric(n_sim)
cov_flex    <- numeric(n_sim)
rmse_incr   <- numeric(n_sim)
cov_incr    <- numeric(n_sim)
G1_bin      <- ifelse(G1 >= stats::median(G1), 1, 0)

for (b in seq_len(n_sim)) {
  s_b <- make_s(n, n_pre = 2, n_post = 2, G1, X1, beta_fn_1, seed = b)

  # Kernel: derivative curve should be flat at 2
  fk <- tryCatch(
    fdid(s_b, tr_period = 1:2, ref_period = 0, method = "kernel",
         eval_g = eval_g1),
    error = function(e) NULL)
  if (!is.null(fk)) {
    rmse_kernel[b] <- sqrt(mean((fk$curve_event$delta_hat - TRUE_delta)^2, na.rm = TRUE))
    cov_kernel[b]  <- mean(
      fk$curve_event$delta_hat - qnorm(0.975) * fk$curve_event$se_delta <= TRUE_delta &
        fk$curve_event$delta_hat + qnorm(0.975) * fk$curve_event$se_delta >= TRUE_delta,
      na.rm = TRUE
    )
  } else {
    rmse_kernel[b] <- NA; cov_kernel[b] <- NA
  }

  # PLR: scalar β ≈ 2
  fp <- tryCatch(
    fdid(s_b, tr_period = 1:2, ref_period = 0, method = "dml_plr", K = 5L, S = 1L),
    error = function(e) NULL)
  if (!is.null(fp)) {
    bias_plr[b] <- fp$est$event$Estimate - TRUE_beta
    cov_plr[b]  <- (fp$est$event$CI_Lower <= TRUE_beta) &
                   (fp$est$event$CI_Upper >= TRUE_beta)
  } else {
    bias_plr[b] <- NA; cov_plr[b] <- NA
  }

  # Binary DML: same linear contrast after median split of G
  fb <- tryCatch(
    fdid(make_s(n, n_pre = 2, n_post = 2, G1_bin, X1, beta_fn_1, seed = b),
         tr_period = 1:2, ref_period = 0, method = "dml_binary", K = 5L, S = 1L),
    error = function(e) NULL)
  if (!is.null(fb)) {
    cov_binary[b] <- (fb$est$event$CI_Lower <= TRUE_beta) &
      (fb$est$event$CI_Upper >= TRUE_beta)
  } else {
    cov_binary[b] <- NA
  }

  # DML flex: level curve ≈ 2g
  ff <- tryCatch(
    fdid(s_b, tr_period = 1:2, ref_period = 0, method = "dml_flex",
         K = 3L, S = 1L, eval_g = eval_g1),
    error = function(e) NULL)
  if (!is.null(ff)) {
    theta_true <- 2 * eval_g1
    rmse_flex[b] <- sqrt(mean((ff$curve_event$theta_hat - theta_true)^2, na.rm = TRUE))
    cov_flex[b]  <- mean(ff$curve_event$ci_lower <= theta_true &
                           ff$curve_event$ci_upper >= theta_true,
                         na.rm = TRUE)
  } else {
    rmse_flex[b] <- NA; cov_flex[b] <- NA
  }

  # DML incremental: derivative curve ≈ 2
  fi <- tryCatch(
    fdid(s_b, tr_period = 1:2, ref_period = 0, method = "dml_incremental",
         K = 3L, S = 1L, eval_g = eval_g1),
    error = function(e) NULL)
  if (!is.null(fi)) {
    rmse_incr[b] <- sqrt(mean((fi$curve_event$delta_hat - TRUE_delta)^2, na.rm = TRUE))
    cov_incr[b]  <- mean(fi$curve_event$delta_ci_lower <= TRUE_delta &
                           fi$curve_event$delta_ci_upper >= TRUE_delta,
                         na.rm = TRUE)
  } else {
    rmse_incr[b] <- NA; cov_incr[b] <- NA
  }
}

cat(sprintf("  Kernel RMSE (target < 0.4):    %.3f  [%s]\n",
            mean(rmse_kernel, na.rm = TRUE),
            ifelse(mean(rmse_kernel, na.rm = TRUE) < 0.4, "PASS", "FAIL")))
cat(sprintf("  Kernel pointwise coverage (target≥0.75): %.3f  [%s]\n",
            mean(cov_kernel, na.rm = TRUE),
            ifelse(mean(cov_kernel, na.rm = TRUE) >= 0.75, "PASS", "FAIL")))
cat(sprintf("  PLR bias  (target |b| < 0.2):  %.3f  [%s]\n",
            abs(mean(bias_plr, na.rm = TRUE)),
            ifelse(abs(mean(bias_plr, na.rm = TRUE)) < 0.2, "PASS", "FAIL")))
cat(sprintf("  PLR 95%% coverage (target≥0.85): %.3f  [%s]\n",
            mean(cov_plr, na.rm = TRUE),
            ifelse(mean(cov_plr, na.rm = TRUE) >= 0.85, "PASS", "FAIL")))
cat(sprintf("  DML binary 95%% coverage (target≥0.85): %.3f  [%s]\n",
            mean(cov_binary, na.rm = TRUE),
            ifelse(mean(cov_binary, na.rm = TRUE) >= 0.85, "PASS", "FAIL")))
cat(sprintf("  Flex RMSE (target < 0.5):       %.3f  [%s]\n",
            mean(rmse_flex, na.rm = TRUE),
            ifelse(mean(rmse_flex, na.rm = TRUE) < 0.5, "PASS", "FAIL")))
cat(sprintf("  Flex pointwise coverage (target≥0.75):  %.3f  [%s]\n",
            mean(cov_flex, na.rm = TRUE),
            ifelse(mean(cov_flex, na.rm = TRUE) >= 0.75, "PASS", "FAIL")))
cat(sprintf("  DML incremental RMSE (target < 0.6): %.3f  [%s]\n",
            mean(rmse_incr, na.rm = TRUE),
            ifelse(mean(rmse_incr, na.rm = TRUE) < 0.6, "PASS", "FAIL")))
cat(sprintf("  DML incremental coverage (target≥0.70): %.3f  [%s]\n",
            mean(cov_incr, na.rm = TRUE),
            ifelse(mean(cov_incr, na.rm = TRUE) >= 0.70, "PASS", "FAIL")))

# ── DGP 5: Null DGP — size control ────────────────────────────────────────────
#   ΔY = 0.5X + ε   →   δ_DID = 0 for all g

cat("\n══════════════════════════════════════════════════════\n")
cat("DGP 5 — Null DGP, size control (n=300, 100 simulations)\n")
cat("══════════════════════════════════════════════════════\n")

set.seed(2025); n5 <- 300
G5 <- runif(n5); X5 <- matrix(rnorm(n5), n5, 1)
beta_fn_5 <- function(g, x) 0.5 * x[1]
n_sim5 <- 100

rej_kernel5 <- numeric(n_sim5)
rej_plr5    <- numeric(n_sim5)
rej_bin5    <- numeric(n_sim5)
rej_incr5   <- numeric(n_sim5)
G_bin5      <- ifelse(G5 >= 0.5, 1, 0)

for (b in seq_len(n_sim5)) {
  # Kernel: pointwise test at g=0.5
  sk <- tryCatch(
    fdid(make_s(n5, 2, 2, G5, X5, beta_fn_5, seed = b + 500),
         tr_period = 1:2, ref_period = 0, method = "kernel", eval_g = c(0.5)),
    error = function(e) NULL)
  rej_kernel5[b] <- if (!is.null(sk) && !is.na(sk$curve_event$delta_hat) &&
                         !is.na(sk$curve_event$se_delta) &&
                         sk$curve_event$se_delta > 0) {
    abs(sk$curve_event$delta_hat / sk$curve_event$se_delta) > qnorm(0.975)
  } else NA

  # PLR: scalar test
  sp <- tryCatch(
    fdid(make_s(n5, 2, 2, G5, X5, beta_fn_5, seed = b + 500),
         tr_period = 1:2, ref_period = 0, method = "dml_plr", K = 5L, S = 1L),
    error = function(e) NULL)
  rej_plr5[b] <- if (!is.null(sp))
    abs(sp$est$event$Estimate / sp$est$event$Std.Error) > qnorm(0.975) else NA

  # DML binary: binary split test
  sb <- tryCatch(
    fdid(make_s(n5, 2, 2, G_bin5, X5, beta_fn_5, seed = b + 500),
         tr_period = 1:2, ref_period = 0, method = "dml_binary", K = 5L, S = 1L),
    error = function(e) NULL)
  rej_bin5[b] <- if (!is.null(sb))
    abs(sb$est$event$Estimate / sb$est$event$Std.Error) > qnorm(0.975) else NA

  # DML incremental: pointwise test at g=0.5
  si <- tryCatch(
    fdid(make_s(n5, 2, 2, G5, X5, beta_fn_5, seed = b + 500),
         tr_period = 1:2, ref_period = 0, method = "dml_incremental",
         eval_g = c(0.5), K = 3L, S = 1L),
    error = function(e) NULL)
  rej_incr5[b] <- if (!is.null(si) && !is.na(si$curve_event$delta_hat) &&
                        !is.na(si$curve_event$se_delta) &&
                        si$curve_event$se_delta > 0) {
    abs(si$curve_event$delta_hat / si$curve_event$se_delta) > qnorm(0.975)
  } else NA
}

cat(sprintf("  Kernel size (target ≤ 0.12):  %.3f  [%s]\n",
            mean(rej_kernel5, na.rm = TRUE),
            ifelse(mean(rej_kernel5, na.rm = TRUE) <= 0.12, "PASS", "FAIL")))
cat(sprintf("  PLR size (target ≤ 0.12):     %.3f  [%s]\n",
            mean(rej_plr5, na.rm = TRUE),
            ifelse(mean(rej_plr5, na.rm = TRUE) <= 0.12, "PASS", "FAIL")))
cat(sprintf("  DML binary size (target≤0.12): %.3f  [%s]\n",
            mean(rej_bin5, na.rm = TRUE),
            ifelse(mean(rej_bin5, na.rm = TRUE) <= 0.12, "PASS", "FAIL")))
cat(sprintf("  DML incremental size (target≤0.15): %.3f  [%s]\n",
            mean(rej_incr5, na.rm = TRUE),
            ifelse(mean(rej_incr5, na.rm = TRUE) <= 0.15, "PASS", "FAIL")))

cat("\nSimulation complete.\n")
