# Kernel estimator for continuous-G FDID
# Theory: FDID Notes §2.4, Liu et al. (2025) §2
# Engineering reference: interflex-master/R/kernel.R (G replaces moderator X)
# Called internally by fdid() when method = "kernel"

run_kernel_method <- function(
  s,
  tr_period,
  ref_period,
  entire_period = NULL,
  eval_g        = NULL,
  h0            = NULL,
  K_folds       = 10,
  boot          = 500,
  alpha         = 0.05,
  trim          = 0.05,
  vartype       = "robust",
  missing_data  = "listwise",
  cluster       = NULL,
  the_call      = NULL
) {

  ycol <- function(t) paste0("Y_", t)

  covar <- grep("^x\\d+$", names(s), value = TRUE)

  all_y_cols    <- grep("^Y_[-]?[0-9]+(\\.[0-9]+)?$", names(s), value = TRUE)
  numeric_times <- sort(unique(suppressWarnings(
    as.numeric(sub("^Y_", "", all_y_cols)))))
  numeric_times <- numeric_times[!is.na(numeric_times)]

  if (length(ref_period) > 1) ref_period <- max(ref_period)

  all_times <- if (!is.null(entire_period)) {
    sort(intersect(numeric_times, entire_period))
  } else {
    numeric_times
  }

  cluster_col <- if (!is.null(cluster) && length(cluster) == 1L &&
                     cluster %in% names(s)) cluster else NULL
  needed_cols <- c("G", covar, cluster_col)
  if (missing_data == "listwise") {
    keep <- complete.cases(s[, c(needed_cols, ycol(all_times)), drop = FALSE])
  } else {
    keep <- complete.cases(s[, needed_cols, drop = FALSE])
  }
  s <- s[keep, , drop = FALSE]

  n     <- nrow(s)
  p     <- length(covar)
  G_vec <- as.numeric(s$G)
  X_mat <- if (p > 0) as.matrix(s[, covar, drop = FALSE]) else matrix(0, n, 0)
  X_bar <- if (p > 0) colMeans(X_mat) else numeric(0)
  cluster_vec <- if (!is.null(cluster_col)) s[[cluster_col]] else NULL
  if (!is.null(cluster_vec) && length(unique(cluster_vec[!is.na(cluster_vec)])) < 2L) {
    cluster_vec <- NULL
  }

  if (n < 10L) {
    stop("method='kernel' needs at least 10 complete observations after missing-data handling.")
  }
  if (!is.numeric(K_folds) || length(K_folds) != 1L || !is.finite(K_folds) ||
      K_folds < 2L || K_folds > n) {
    stop("K_folds must be between 2 and the number of complete observations (", n, ").")
  }
  K_folds <- as.integer(K_folds)
  if (!is.numeric(boot) || length(boot) != 1L || !is.finite(boot) || boot < 1L) {
    stop("boot must be a positive integer-like value.")
  }
  boot <- as.integer(boot)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single finite number strictly between 0 and 1.")
  }
  if (!is.numeric(trim) || length(trim) != 1L || !is.finite(trim) ||
      trim < 0 || trim >= 0.5) {
    stop("trim must be a single finite number in [0, 0.5).")
  }
  if (!is.null(h0) &&
      (!is.numeric(h0) || length(h0) != 1L || !is.finite(h0) || h0 <= 0)) {
    stop("h0 must be a single positive finite number when supplied.")
  }

  if (length(unique(G_vec)) <= 2)
    stop("method='kernel' requires continuous G (more than 2 unique values). ",
         "Use method='aipw' or method='did' for binary G.")

  if (is.null(eval_g)) {
    g_lo   <- as.numeric(quantile(G_vec, trim))
    g_hi   <- as.numeric(quantile(G_vec, 1 - trim))
    eval_g <- seq(g_lo, g_hi, length.out = 50)
  } else {
    eval_g <- sort(as.numeric(eval_g))
  }
  if (length(eval_g) == 0L || any(!is.finite(eval_g))) {
    stop("eval_g must be a nonempty finite numeric vector.")
  }
  m_eval <- length(eval_g)
  np     <- 2L + 2L * p

  # ── Selection vectors (FDID Notes §2.4) ──
  a_vec <- function(x) if (p > 0) c(1, x, rep(0, 1 + p)) else c(1, 0)
  b_vec <- function(x) if (p > 0) c(rep(0, 1 + p), 1, x) else c(0, 1)

  # ── WLS local-linear fit at one evaluation point ──
  fit_local <- function(g0, h_g0, DeltaY, G, X) {
    u      <- (G - g0) / h_g0
    W_diag <- dnorm(u) / h_g0
    max_w  <- max(W_diag)
    if (max_w < .Machine$double.eps) return(rep(NA_real_, 2L + 2L * ncol(X)))
    eff_n  <- sum(W_diag) / max_w
    if (eff_n < 5) return(rep(NA_real_, 2L + 2L * ncol(X)))
    dG <- G - g0
    R  <- if (ncol(X) > 0) cbind(1, X, dG, dG * X) else cbind(1, dG)
    WR <- R * W_diag
    A  <- crossprod(WR, R)
    if (rcond(A) < .Machine$double.eps * 100) return(rep(NA_real_, ncol(R)))
    bv <- crossprod(WR, DeltaY)
    tryCatch(as.numeric(solve(A, bv)), error = function(e) rep(NA_real_, ncol(R)))
  }

  # ── Sandwich SE at one evaluation point (HC0) ──
  sandwich_SE_kernel <- function(g0, h_g0, theta_g0, DeltaY, G, X) {
    k <- 2L + 2L * ncol(X)
    if (any(is.na(theta_g0))) return(matrix(NA_real_, k, k))
    dG     <- G - g0
    W_diag <- dnorm((G - g0) / h_g0) / h_g0
    R      <- if (ncol(X) > 0) cbind(1, X, dG, dG * X) else cbind(1, dG)
    e_hat  <- as.numeric(DeltaY - R %*% theta_g0)
    WR     <- R * W_diag
    A      <- crossprod(WR, R)
    A_inv  <- tryCatch(solve(A), error = function(e) matrix(NA_real_, k, k))
    if (anyNA(A_inv)) return(matrix(NA_real_, k, k))
    meat   <- crossprod(R * (W_diag^2 * e_hat^2), R)
    A_inv %*% meat %*% A_inv
  }

  # ── Adaptive bandwidth at evaluation points ──
  adapt_bw <- function(h0_cand, rho_pts, rho_gm) {
    h0_cand * sqrt(rho_gm / pmax(rho_pts, 1e-8))
  }

  # ── 10-fold LSCV for baseline bandwidth ──
  select_h0_LSCV <- function(DeltaY, G, X, nfolds) {
    n_cv    <- length(DeltaY)
    sg      <- sd(G); if (sg < 1e-10) sg <- 1
    # Floor ensures enough local observations; ceiling avoids over-smoothing
    h_min   <- max(0.05 * sg, (2L + 2L * ncol(X) + 1L) / n_cv)
    h_max   <- 2.0 * sg
    if (h_min >= h_max) h_max <- 4.0 * h_min
    h0_grid <- exp(seq(log(h_min), log(h_max), length.out = 25))
    fold_id <- sample(rep(seq_len(nfolds), length.out = n_cv))

    CV <- vapply(h0_grid, function(h_cand) {
      err2 <- numeric(n_cv)
      for (k in seq_len(nfolds)) {
        tr <- which(fold_id != k)
        ho <- which(fold_id == k)
        if (length(tr) < 10) next
        sg_tr  <- sd(G[tr]); if (sg_tr < 1e-10) sg_tr <- 1
        dens_tr <- density(G[tr], bw = 1.06 * sg_tr * length(tr)^(-1/5), n = 512)
        rho_fn  <- approxfun(dens_tr, rule = 2)
        rho_tr  <- pmax(rho_fn(G[tr]), 1e-8)
        rho_gm_k <- exp(mean(log(rho_tr)))
        for (i in ho) {
          rho_i <- pmax(rho_fn(G[i]), 1e-8)
          h_i   <- h_cand * sqrt(rho_gm_k / rho_i)
          th    <- fit_local(G[i], h_i, DeltaY[tr], G[tr],
                             if (ncol(X) > 0) X[tr, , drop = FALSE] else X)
          xi    <- if (p > 0) X[i, ] else numeric(0)
          err2[i] <- if (!any(is.na(th))) {
            (DeltaY[i] - as.numeric(a_vec(xi) %*% th))^2
          } else {
            NA_real_
          }
        }
      }
      # exclude NA (failed fits) from mean; penalise h that causes mass failures
      n_ok <- sum(!is.na(err2))
      if (n_ok < n_cv / 2) return(Inf)
      mean(err2, na.rm = TRUE)
    }, numeric(1))

    h0_grid[which.min(CV)]
  }

  # ── Estimate curve over evaluation grid ──
  estimate_curve_kernel <- function(DeltaY, G, X, evg, h0_star, rho_ev, rho_gm, Xbar) {
    m_c  <- length(evg)
    np_c <- 2L + 2L * ncol(X)
    theta_mat <- matrix(NA_real_, m_c, np_c)
    mu_hat    <- rep(NA_real_, m_c)
    delta_hat <- rep(NA_real_, m_c)
    se_mu     <- rep(NA_real_, m_c)
    se_delta  <- rep(NA_real_, m_c)
    h_ev      <- adapt_bw(h0_star, rho_ev, rho_gm)

    for (j in seq_len(m_c)) {
      g0   <- evg[j]
      h_g0 <- h_ev[j]
      th   <- fit_local(g0, h_g0, DeltaY, G, X)
      V_j  <- sandwich_SE_kernel(g0, h_g0, th, DeltaY, G, X)
      theta_mat[j, ] <- th
      if (!any(is.na(th))) {
        av <- a_vec(Xbar); bv <- b_vec(Xbar)
        mu_hat[j]    <- as.numeric(av %*% th)
        delta_hat[j] <- as.numeric(bv %*% th)
        if (!any(is.na(V_j))) {
          se_mu[j]    <- sqrt(pmax(0, as.numeric(t(av) %*% V_j %*% av)))
          se_delta[j] <- sqrt(pmax(0, as.numeric(t(bv) %*% V_j %*% bv)))
        }
      }
    }
    list(theta_mat = theta_mat, mu_hat = mu_hat, delta_hat = delta_hat,
         se_mu = se_mu, se_delta = se_delta, h_adaptive = h_ev)
  }

  kernel_pointwise_ci <- function(boot_mat, alpha_val) {
    if (is.null(boot_mat) || !is.matrix(boot_mat)) return(NULL)
    list(
      se = apply(boot_mat, 2L, stats::sd, na.rm = TRUE),
      lower = apply(boot_mat, 2L, stats::quantile, probs = alpha_val / 2,
                    na.rm = TRUE, names = FALSE),
      upper = apply(boot_mat, 2L, stats::quantile, probs = 1 - alpha_val / 2,
                    na.rm = TRUE, names = FALSE),
      method = "bootstrap_percentile"
    )
  }

  kernel_curve_vcov <- function(boot_mat) {
    if (is.null(boot_mat) || !is.matrix(boot_mat) || nrow(boot_mat) < 2L) {
      return(NULL)
    }
    V <- stats::cov(boot_mat, use = "pairwise.complete.obs")
    rownames(V) <- colnames(V) <- NULL
    V
  }

  kernel_uniform_quantile_envelope <- function(boot_mat, alpha_val,
                                               tol = 1e-4, max_iter = 60L) {
    if (is.null(boot_mat) || !is.matrix(boot_mat) || nrow(boot_mat) < 2L) {
      return(list(
        lower = rep(NA_real_, if (is.null(boot_mat)) 0L else ncol(boot_mat)),
        upper = rep(NA_real_, if (is.null(boot_mat)) 0L else ncol(boot_mat)),
        zeta_hat = NA_real_, coverage = NA_real_, n_boot_eff = 0L,
        method = "bootstrap_quantile_envelope"
      ))
    }
    k <- ncol(boot_mat)
    complete <- apply(boot_mat, 1L, function(row) all(is.finite(row)))
    theta <- boot_mat[complete, , drop = FALSE]
    if (nrow(theta) < 2L) {
      theta <- boot_mat
      theta <- theta[rowSums(is.finite(theta)) > 0L, , drop = FALSE]
    }
    n_eff <- nrow(theta)
    if (n_eff < 2L || k < 1L) {
      return(list(
        lower = rep(NA_real_, k), upper = rep(NA_real_, k),
        zeta_hat = NA_real_, coverage = NA_real_, n_boot_eff = n_eff,
        method = "bootstrap_quantile_envelope"
      ))
    }

    coverage_at <- function(zeta) {
      q_lo <- apply(theta, 2L, stats::quantile, probs = zeta,
                    na.rm = TRUE, names = FALSE)
      q_hi <- apply(theta, 2L, stats::quantile, probs = 1 - zeta,
                    na.rm = TRUE, names = FALSE)
      inside <- apply(theta, 1L, function(row) {
        ok <- is.finite(row) & is.finite(q_lo) & is.finite(q_hi)
        any(ok) && all(row[ok] >= q_lo[ok] & row[ok] <= q_hi[ok])
      })
      list(coverage = mean(inside, na.rm = TRUE), lower = q_lo, upper = q_hi)
    }

    lo <- alpha_val / (2 * k)
    hi <- alpha_val / 2
    best <- coverage_at(lo)
    zeta_hat <- lo
    if (!is.finite(best$coverage) || best$coverage < 1 - alpha_val) {
      return(list(
        lower = best$lower, upper = best$upper, zeta_hat = zeta_hat,
        coverage = best$coverage, n_boot_eff = n_eff,
        method = "bootstrap_quantile_envelope"
      ))
    }
    for (iter in seq_len(max_iter)) {
      mid <- (lo + hi) / 2
      cand <- coverage_at(mid)
      if (is.finite(cand$coverage) && cand$coverage >= 1 - alpha_val) {
        zeta_hat <- mid
        best <- cand
        lo <- mid
      } else {
        hi <- mid
      }
      if (abs(hi - lo) < tol) break
    }
    list(
      lower = best$lower, upper = best$upper, zeta_hat = zeta_hat,
      coverage = best$coverage, n_boot_eff = n_eff,
      method = "bootstrap_quantile_envelope"
    )
  }

  kernel_bootstrap_indices <- function(n_obs, cluster_labels = NULL) {
    if (is.null(cluster_labels)) return(sample.int(n_obs, replace = TRUE))
    cl <- unique(cluster_labels[!is.na(cluster_labels)])
    draw <- sample(cl, length(cl), replace = TRUE)
    unlist(lapply(draw, function(z) which(cluster_labels == z)), use.names = FALSE)
  }

  # ── Scalar summary from curve ──
  curve_scalar <- function(cv, evg, z_c) {
    if (length(evg) == 1L) {
      db <- cv$delta_hat[1L]
      if (!is.null(cv$delta_bootstrap) && is.matrix(cv$delta_bootstrap) &&
          nrow(cv$delta_bootstrap) >= 2L) {
        rr <- cv$delta_bootstrap[, 1L]
        rr <- rr[is.finite(rr)]
        sb <- stats::sd(rr, na.rm = TRUE)
        ci <- stats::quantile(rr, probs = c(alpha / 2, 1 - alpha / 2),
                              na.rm = TRUE, names = FALSE)
        return(data.frame(Estimate = db, Std.Error = sb,
                          CI_Lower = ci[1L], CI_Upper = ci[2L],
                          SE_Method = "bootstrap_replicate",
                          CI_Method = "bootstrap_percentile"))
      }
      sb <- cv$se_delta[1L]
      return(data.frame(Estimate = db, Std.Error = sb,
                        CI_Lower = db - z_c * sb, CI_Upper = db + z_c * sb,
                        SE_Method = "analytical_pointwise",
                        CI_Method = "normal_pointwise"))
    }
    span <- evg[length(evg)] - evg[1]
    tau  <- cv$mu_hat[length(evg)] - cv$mu_hat[1]
    db   <- if (!is.na(tau)  && span > 0) tau  / span else NA_real_
    if (!is.null(cv$mu_bootstrap) && is.matrix(cv$mu_bootstrap) &&
        nrow(cv$mu_bootstrap) >= 2L && span > 0) {
      rr <- (cv$mu_bootstrap[, length(evg)] - cv$mu_bootstrap[, 1L]) / span
      rr <- rr[is.finite(rr)]
      if (length(rr) >= 2L) {
        sb <- stats::sd(rr, na.rm = TRUE)
        ci <- stats::quantile(rr, probs = c(alpha / 2, 1 - alpha / 2),
                              na.rm = TRUE, names = FALSE)
        return(data.frame(Estimate = db, Std.Error = sb,
                          CI_Lower = ci[1L], CI_Upper = ci[2L],
                          SE_Method = "bootstrap_replicate",
                          CI_Method = "bootstrap_percentile"))
      }
    }
    if (!is.null(cv$mu_boot_vcov) && is.matrix(cv$mu_boot_vcov) &&
        all(dim(cv$mu_boot_vcov) >= length(evg)) && span > 0) {
      V <- cv$mu_boot_vcov
      v_t <- V[length(evg), length(evg)] + V[1L, 1L] - 2 * V[length(evg), 1L]
      sb <- sqrt(pmax(0, v_t)) / span
      return(data.frame(Estimate = db, Std.Error = sb,
                        CI_Lower = db - z_c * sb, CI_Upper = db + z_c * sb,
                        SE_Method = "bootstrap_curve_vcov",
                        CI_Method = "normal_curve_vcov"))
    }
    se_t <- sqrt(pmax(0, cv$se_mu[length(evg)]^2 + cv$se_mu[1]^2, na.rm = TRUE))
    sb   <- if (!is.na(se_t) && span > 0) se_t / span else NA_real_
    data.frame(Estimate = db, Std.Error = sb,
               CI_Lower = db - z_c * sb, CI_Upper = db + z_c * sb,
               SE_Method = "endpoint_independent_fallback",
               CI_Method = "normal_pointwise_fallback")
  }

  bootstrap_band_kernel <- function(DeltaY, G, X, evg, h0_star, boot_B, Xbar,
                                    cluster_labels = NULL) {
    mu_boot <- matrix(NA_real_, boot_B, length(evg))
    delta_boot <- matrix(NA_real_, boot_B, length(evg))
    for (b in seq_len(boot_B)) {
      idx_b <- kernel_bootstrap_indices(length(DeltaY), cluster_labels)
      G_b   <- G[idx_b]
      DY_b  <- DeltaY[idx_b]
      X_b   <- if (ncol(X) > 0) X[idx_b, , drop = FALSE] else X
      sg_b  <- sd(G_b); if (sg_b < 1e-10) sg_b <- 1
      dens_b <- density(G_b, bw = 1.06 * sg_b * length(G_b)^(-1/5), n = 512)
      rho_fn_b   <- approxfun(dens_b, rule = 2)
      rho_b_ev   <- pmax(rho_fn_b(evg), 1e-8)
      rho_b_pts  <- pmax(rho_fn_b(G_b), 1e-8)
      rho_gm_b   <- exp(mean(log(rho_b_pts)))
      h0_b       <- tryCatch(
        select_h0_LSCV(DY_b, G_b, X_b, min(K_folds, max(2L, length(DY_b) %/% 5L))),
        error = function(e) h0_star
      )
      cv_b <- estimate_curve_kernel(DY_b, G_b, X_b, evg, h0_b, rho_b_ev, rho_gm_b, Xbar)
      mu_boot[b, ] <- cv_b$mu_hat
      delta_boot[b, ] <- cv_b$delta_hat
    }
    colnames(mu_boot) <- colnames(delta_boot) <- paste0("g", seq_along(evg))
    mu_pt <- kernel_pointwise_ci(mu_boot, alpha)
    delta_pt <- kernel_pointwise_ci(delta_boot, alpha)
    mu_uni <- kernel_uniform_quantile_envelope(mu_boot, alpha)
    delta_uni <- kernel_uniform_quantile_envelope(delta_boot, alpha)
    list(
      mu_band_lower = mu_uni$lower,
      mu_band_upper = mu_uni$upper,
      delta_band_lower = delta_uni$lower,
      delta_band_upper = delta_uni$upper,
      mu_bootstrap = mu_boot,
      delta_bootstrap = delta_boot,
      mu_boot_se = mu_pt$se,
      delta_boot_se = delta_pt$se,
      mu_boot_ci_lower = mu_pt$lower,
      mu_boot_ci_upper = mu_pt$upper,
      delta_boot_ci_lower = delta_pt$lower,
      delta_boot_ci_upper = delta_pt$upper,
      mu_boot_vcov = kernel_curve_vcov(mu_boot),
      delta_boot_vcov = kernel_curve_vcov(delta_boot),
      mu_uniform_method = mu_uni$method,
      delta_uniform_method = delta_uni$method,
      mu_uniform_zeta = mu_uni$zeta_hat,
      delta_uniform_zeta = delta_uni$zeta_hat,
      mu_uniform_coverage = mu_uni$coverage,
      delta_uniform_coverage = delta_uni$coverage,
      mu_uniform_n_boot_eff = mu_uni$n_boot_eff,
      delta_uniform_n_boot_eff = delta_uni$n_boot_eff
    )
  }

  # ── Pilot density (global sample) ──
  sg_full   <- sd(G_vec); if (sg_full < 1e-10) sg_full <- 1
  bw_pilot  <- 1.06 * sg_full * n^(-1/5)
  dens_full <- density(G_vec, bw = bw_pilot, n = 1024)
  rho_fn    <- approxfun(dens_full, rule = 2)
  rho_G     <- pmax(rho_fn(G_vec), 1e-8)
  rho_eval  <- pmax(rho_fn(eval_g), 1e-8)
  rho_bar_GM <- exp(mean(log(rho_G)))

  # ── Event-period tempY ──
  Y_ref_col  <- ycol(ref_period)
  Y_tr_cols  <- setdiff(ycol(tr_period), Y_ref_col)
  tempY_event <- rowMeans(s[, Y_tr_cols, drop = FALSE], na.rm = TRUE) - s[[Y_ref_col]]

  # ── Bandwidth selection ──
  h0_star <- if (is.null(h0)) {
    select_h0_LSCV(tempY_event, G_vec, X_mat, K_folds)
  } else {
    as.numeric(h0)
  }

  # ── Event-period curve ──
  cv_event <- estimate_curve_kernel(tempY_event, G_vec, X_mat,
                                    eval_g, h0_star, rho_eval, rho_bar_GM, X_bar)

  # ── Confidence band ──
  z_crit <- qnorm(1 - alpha / 2)
  band <- if (vartype == "bootstrap") {
    bootstrap_band_kernel(tempY_event, G_vec, X_mat, eval_g, h0_star, boot, X_bar,
                          cluster_vec)
  } else {
    list(
      mu_band_lower = cv_event$mu_hat - z_crit * cv_event$se_mu,
      mu_band_upper = cv_event$mu_hat + z_crit * cv_event$se_mu,
      delta_band_lower = cv_event$delta_hat - z_crit * cv_event$se_delta,
      delta_band_upper = cv_event$delta_hat + z_crit * cv_event$se_delta,
      mu_bootstrap = NULL,
      delta_bootstrap = NULL,
      mu_boot_se = NULL,
      delta_boot_se = NULL,
      mu_boot_ci_lower = NULL,
      mu_boot_ci_upper = NULL,
      delta_boot_ci_lower = NULL,
      delta_boot_ci_upper = NULL,
      mu_boot_vcov = NULL,
      delta_boot_vcov = NULL,
      mu_uniform_method = "normal_pointwise",
      delta_uniform_method = "normal_pointwise",
      mu_uniform_zeta = NA_real_,
      delta_uniform_zeta = NA_real_,
      mu_uniform_coverage = NA_real_,
      delta_uniform_coverage = NA_real_,
      mu_uniform_n_boot_eff = 0L,
      delta_uniform_n_boot_eff = 0L
    )
  }
  cv_event$mu_band_lower <- band$mu_band_lower
  cv_event$mu_band_upper <- band$mu_band_upper
  cv_event$delta_band_lower <- band$delta_band_lower
  cv_event$delta_band_upper <- band$delta_band_upper
  for (nm in setdiff(names(band), c("mu_band_lower", "mu_band_upper",
                                    "delta_band_lower", "delta_band_upper"))) {
    cv_event[[nm]] <- band[[nm]]
  }
	  cv_event$simultaneous_band_method <- if (vartype == "bootstrap") {
	    "bootstrap_quantile_envelope"
	  } else {
	    "normal_pointwise"
	  }
	  cv_event$bootstrap_zeta <- cv_event$delta_uniform_zeta

  # ── Scalar event summary ──
  est_event_df <- curve_scalar(cv_event, eval_g, z_crit)

  # ── Dynamic loop ──
  dynamic_df <- data.frame(
    Estimate = numeric(length(all_times)), Std.Error = numeric(length(all_times)),
    CI_Lower = numeric(length(all_times)), CI_Upper  = numeric(length(all_times)),
    row.names = as.character(all_times)
  )
  curve_dynamic <- vector("list", length(all_times))
  names(curve_dynamic) <- as.character(all_times)

  for (i in seq_along(all_times)) {
    t_yr <- all_times[i]
    if (t_yr == ref_period) {
      dynamic_df[i, ] <- 0
      next
    }
    t_col <- ycol(t_yr)
    if (!t_col %in% names(s)) next
	    tempY_t <- s[[t_col]] - s[[Y_ref_col]]
	    ok      <- !is.na(tempY_t)
	    if (sum(ok) < 10) next
	    Xbar_t <- if (p > 0) colMeans(X_mat[ok, , drop = FALSE]) else numeric(0)
	    cv_t   <- estimate_curve_kernel(tempY_t[ok], G_vec[ok],
	                                    if (p > 0) X_mat[ok, , drop = FALSE] else X_mat,
	                                    eval_g, h0_star, rho_eval, rho_bar_GM, Xbar_t)
	    df_t   <- curve_scalar(cv_t, eval_g, z_crit)
	    dynamic_df[i, ]  <- unlist(df_t[names(dynamic_df)])
	    curve_dynamic[[i]] <- cv_t
  }

  # ── Pre/post aggregates ──
  pre_times  <- all_times[all_times < min(tr_period)]
  post_times <- all_times[all_times > max(tr_period)]

  kernel_aggregate <- function(time_vec) {
    if (length(time_vec) == 0)
      return(data.frame(Estimate = NA, Std.Error = NA, CI_Lower = NA, CI_Upper = NA))
    t_cols <- setdiff(unique(ycol(time_vec)), Y_ref_col)
    if (length(t_cols) == 0)
      return(data.frame(Estimate = NA, Std.Error = NA, CI_Lower = NA, CI_Upper = NA))
    tmp <- rowMeans(s[, t_cols, drop = FALSE], na.rm = TRUE)
    tmp[is.nan(tmp)] <- NA_real_
    DY_a <- tmp - s[[Y_ref_col]]
    ok   <- !is.na(DY_a)
    if (sum(ok) < 10)
      return(data.frame(Estimate = NA, Std.Error = NA, CI_Lower = NA, CI_Upper = NA))
    Xbar_a <- if (p > 0) colMeans(X_mat[ok, , drop = FALSE]) else numeric(0)
    cv_a   <- estimate_curve_kernel(DY_a[ok], G_vec[ok],
                                    if (p > 0) X_mat[ok, , drop = FALSE] else X_mat,
                                    eval_g, h0_star, rho_eval, rho_bar_GM, Xbar_a)
    curve_scalar(cv_a, eval_g, z_crit)
  }

  pre_event_result  <- kernel_aggregate(pre_times)
  post_event_result <- kernel_aggregate(post_times)

  # ── raw_means: median-split approximation for continuous G ──
  zval  <- stats::qnorm(0.975)
  G_med <- median(G_vec, na.rm = TRUE)
  rawdf <- do.call(rbind, lapply(all_times, function(yr) {
    col_yr <- ycol(yr)
    if (!col_yr %in% names(s)) return(NULL)
    y1 <- s[[col_yr]][G_vec >= G_med]; y0 <- s[[col_yr]][G_vec < G_med]
    n1 <- sum(!is.na(y1)); m1 <- mean(y1, na.rm = TRUE)
    n0 <- sum(!is.na(y0)); m0 <- mean(y0, na.rm = TRUE)
    se1 <- if (n1 <= 1) NA_real_ else stats::sd(y1, na.rm = TRUE) / sqrt(n1)
    se0 <- if (n0 <= 1) NA_real_ else stats::sd(y0, na.rm = TRUE) / sqrt(n0)
    rbind(
      data.frame(time = yr, group = "Group 1", meanY = m1, n = n1, seY = se1,
                 CI_Lower = m1 - zval * se1, CI_Upper = m1 + zval * se1,
                 stringsAsFactors = FALSE),
      data.frame(time = yr, group = "Group 0", meanY = m0, n = n0, seY = se0,
                 CI_Lower = m0 - zval * se0, CI_Upper = m0 + zval * se0,
                 stringsAsFactors = FALSE)
    )
  }))
  if (is.null(rawdf)) rawdf <- data.frame(
    time = numeric(0), group = character(0), meanY = numeric(0), n = integer(0),
    seY = numeric(0), CI_Lower = numeric(0), CI_Upper = numeric(0),
    stringsAsFactors = FALSE)

  # ── Assemble fdid object ──
  out <- list(
    est = list(
      pre   = pre_event_result,
      event = est_event_df,
      post  = post_event_result
    ),
    dynamic       = dynamic_df,
    raw_means     = rawdf,
    tr_period     = tr_period,
    ref_period    = ref_period,
	    entire_period = all_times,
	    method        = "kernel",
	    vartype       = vartype,
	    alpha         = alpha,
	    times         = numeric_times,
	    G             = G_vec,
	    cluster       = cluster_vec,
	    ps            = NULL,
	    call          = the_call,
	    target.pop    = "all",
	    eval_g        = eval_g,
	    curve_event   = cv_event,
	    curve_dynamic = curve_dynamic,
	    band_lower    = band$delta_band_lower,
	    band_upper    = band$delta_band_upper,
	    bootstrap_resampling = if (is.null(cluster_vec)) "row" else "cluster",
	    n_clusters    = if (is.null(cluster_vec)) NA_integer_ else
	      length(unique(cluster_vec[!is.na(cluster_vec)])),
	    h0            = h0_star,
	    h_adaptive    = cv_event$h_adaptive
  )
  class(out) <- "fdid"
  out
}
