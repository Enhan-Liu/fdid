.fdid_curve_vec <- function(value, n) {
  if (is.null(value) || length(value) != n) {
    return(rep(NA_real_, n))
  }
  as.numeric(value)
}

.fdid_curve_has_values <- function(value, n) {
  !is.null(value) && length(value) == n && any(is.finite(as.numeric(value)))
}

.fdid_scalar_label <- function(value, default = NA_character_) {
  if (is.null(value) || length(value) == 0L || is.na(value[1L])) {
    return(default)
  }
  as.character(value[1L])
}

.fdid_scalar_number <- function(value, default = NA_real_) {
  if (is.null(value) || length(value) == 0L) return(default)
  value <- suppressWarnings(as.numeric(value[1L]))
  if (is.finite(value)) value else default
}

.fdid_curve_band_label <- function(method, band_method,
                                   practical_band_method = NA_character_) {
  method <- .fdid_scalar_label(method)
  band_method <- .fdid_scalar_label(band_method)
  practical_band_method <- .fdid_scalar_label(practical_band_method)
  if (is.na(band_method) || !nzchar(band_method)) {
    return(NA_character_)
  }
  if (identical(method, "kernel") &&
      identical(band_method, "bootstrap_quantile_envelope")) {
    return("Bootstrap uniform band")
  }
  if (identical(method, "dml_flex")) {
    if (identical(band_method, "gaussian_max_t_blp_spline_covariance")) {
      return("BLP max-t band")
    }
    if (identical(band_method, "practical_signal_residual_multiplier") ||
        grepl("practical residual multiplier", practical_band_method,
              fixed = TRUE)) {
      return("Practical DML band")
    }
    return("DML curve band")
  }
  "Uniform band"
}

.fdid_curve_band_caveat <- function(method, band_method) {
  method <- .fdid_scalar_label(method)
  band_method <- .fdid_scalar_label(band_method)
  if (!identical(method, "dml_flex") || is.na(band_method) ||
      !nzchar(band_method)) {
    return(NA_character_)
  }
  if (identical(band_method, "gaussian_max_t_blp_spline_covariance")) {
    return("Finite-grid Gaussian max-t band from stored BLP-spline covariance; not a DML refit bootstrap.")
  }
  "Practical signal-mapping band; not a final theorem-level uniform DML band."
}

.fdid_curve_plot_data <- function(x, curve,
                                  interval = c("auto", "pointwise", "uniform",
                                               "both", "none"),
                                  xlim = NULL,
                                  ci = TRUE,
                                  show.uniform.CI = TRUE) {
  curve <- match.arg(curve, c("auto", "level", "derivative"))
  interval <- match.arg(interval)
  if (identical(curve, "auto")) {
    curve <- if (identical(x$method, "dml_flex")) "level" else "derivative"
  }
  ce <- x$curve_event
  evg <- as.numeric(x$eval_g)
  n <- length(evg)
  alpha_curve <- x$alpha %||% 0.05
  if (!is.numeric(alpha_curve) || length(alpha_curve) != 1L ||
      !is.finite(alpha_curve) || alpha_curve <= 0 || alpha_curve >= 1) {
    alpha_curve <- 0.05
  }
  z <- stats::qnorm(1 - alpha_curve / 2)

  method <- x$method %||% NA_character_
  vartype <- x$vartype %||% NA_character_
  dml_meta <- x$dml_metadata %||% ce$metadata %||% list()
  target_estimand <- .fdid_scalar_label(dml_meta$target_estimand)
  inference_scope <- .fdid_scalar_label(dml_meta$inference_scope)
  inference_method <- .fdid_scalar_label(dml_meta$inference_method)
  signal_map <- .fdid_scalar_label(x$signal_map %||% ce$signal_map)
  density_method <- .fdid_scalar_label(x$density_method %||% ce$density_method)
  practical_band_method <- .fdid_scalar_label(ce$practical_band_method)
  pointwise_method <- "normal_pointwise"
  pointwise_label <- "Pointwise CI"
  band_method <- NA_character_
  band_zeta <- NA_real_
  band_coverage <- NA_real_
  n_boot_eff <- NA_integer_
  gaussian_crit <- NA_real_
  simultaneous_crit <- NA_real_
  curve_vcov_available <- FALSE

  if (identical(curve, "level")) {
    if (identical(method, "kernel")) {
      estimate <- .fdid_curve_vec(ce$mu_hat, n)
      se <- .fdid_curve_vec(ce$mu_boot_se %||% ce$se_mu, n)
      ci_lo_stored <- ce$mu_boot_ci_lower
      ci_hi_stored <- ce$mu_boot_ci_upper
      ci_method_stored <- "bootstrap_percentile"
      band_lo_stored <- ce$mu_band_lower
      band_hi_stored <- ce$mu_band_upper
      band_method <- ce$mu_uniform_method %||% ce$simultaneous_band_method
      band_zeta <- ce$mu_uniform_zeta %||% NA_real_
      band_coverage <- ce$mu_uniform_coverage %||% NA_real_
      n_boot_eff <- ce$mu_uniform_n_boot_eff %||% NA_integer_
      ylab_default <- "Level DID-response mu(g)"
      curve_label <- "Level curve"
    } else {
      estimate <- .fdid_curve_vec(ce$theta_hat, n)
      se <- .fdid_curve_vec(ce$se_theta, n)
      ci_lo_stored <- ce$theta_ci_lower
      ci_hi_stored <- ce$theta_ci_upper
      ci_method_stored <- "second_stage_mapping_pointwise"
      band_lo_stored <- ce$theta_band_lower
      band_hi_stored <- ce$theta_band_upper
      band_method <- ce$simultaneous_band_method %||% "stored_simultaneous_band"
      gaussian_crit <- .fdid_scalar_number(ce$gaussian_crit_theta)
      simultaneous_crit <- .fdid_scalar_number(ce$simultaneous_crit_theta)
      curve_vcov_available <- !is.null(ce$theta_vcov) && is.matrix(ce$theta_vcov)
      ylab_default <- "Level DID-response theta(g)"
      curve_label <- "Level curve"
    }
  } else {
    estimate <- .fdid_curve_vec(ce$delta_hat, n)
    se <- .fdid_curve_vec(ce$delta_boot_se %||% ce$se_delta, n)
    ci_lo_stored <- if (identical(method, "kernel")) ce$delta_boot_ci_lower else NULL
    ci_hi_stored <- if (identical(method, "kernel")) ce$delta_boot_ci_upper else NULL
    ci_method_stored <- if (identical(method, "kernel")) "bootstrap_percentile" else NA_character_
    if (!identical(method, "kernel")) {
      ci_lo_stored <- ce$delta_ci_lower
      ci_hi_stored <- ce$delta_ci_upper
      ci_method_stored <- "second_stage_mapping_pointwise"
      gaussian_crit <- .fdid_scalar_number(ce$gaussian_crit_delta)
      simultaneous_crit <- .fdid_scalar_number(ce$simultaneous_crit_delta)
      curve_vcov_available <- !is.null(ce$delta_vcov) && is.matrix(ce$delta_vcov)
    }
    band_lo_stored <- if (identical(method, "kernel")) ce$delta_band_lower %||% x$band_lower else ce$delta_band_lower
    band_hi_stored <- if (identical(method, "kernel")) ce$delta_band_upper %||% x$band_upper else ce$delta_band_upper
    band_method <- if (identical(method, "kernel")) {
      ce$delta_uniform_method %||% ce$simultaneous_band_method
    } else {
      ce$simultaneous_band_method %||% "stored_simultaneous_band"
    }
    band_zeta <- ce$delta_uniform_zeta %||% ce$bootstrap_zeta %||% NA_real_
    band_coverage <- ce$delta_uniform_coverage %||% NA_real_
    n_boot_eff <- ce$delta_uniform_n_boot_eff %||% NA_integer_
    ylab_default <- "Derivative DID effect delta(g) = d theta(g)/dG"
    curve_label <- "Derivative curve"
  }

  if (is.null(band_method) || length(band_method) == 0L) {
    band_method <- NA_character_
  }
  band_method <- as.character(band_method[1L])
  pointwise_method <- as.character(pointwise_method[1L])

  if (length(estimate) != n) {
    stop("The selected curve is missing or does not match eval_g.")
  }

  has_boot_ci <- .fdid_curve_has_values(ci_lo_stored, n) &&
    .fdid_curve_has_values(ci_hi_stored, n)
  if (has_boot_ci) {
    pointwise_lower <- as.numeric(ci_lo_stored)
    pointwise_upper <- as.numeric(ci_hi_stored)
    pointwise_method <- ci_method_stored
    pointwise_label <- if (identical(pointwise_method, "bootstrap_percentile")) {
      "Bootstrap CI"
    } else if (identical(method, "dml_flex")) {
      "Mapping CI"
    } else {
      "Pointwise CI"
    }
  } else {
    pointwise_lower <- estimate - z * se
    pointwise_upper <- estimate + z * se
  }

  has_true_band <- .fdid_curve_has_values(band_lo_stored, n) &&
    .fdid_curve_has_values(band_hi_stored, n)
  if (identical(method, "kernel")) {
    has_true_band <- has_true_band &&
      (identical(band_method, "bootstrap_quantile_envelope") ||
         (identical(vartype, "bootstrap") && !identical(band_method, "normal_pointwise")))
    if (has_true_band && (is.na(band_method) || !nzchar(band_method))) {
      band_method <- "bootstrap_quantile_envelope"
    }
  }
  band_lower <- if (has_true_band) as.numeric(band_lo_stored) else rep(NA_real_, n)
  band_upper <- if (has_true_band) as.numeric(band_hi_stored) else rep(NA_real_, n)
  if (!has_true_band) {
    band_method <- NA_character_
    band_zeta <- NA_real_
    band_coverage <- NA_real_
    n_boot_eff <- NA_integer_
    gaussian_crit <- NA_real_
    simultaneous_crit <- NA_real_
  }
  band_label <- if (has_true_band) {
    .fdid_curve_band_label(method, band_method, practical_band_method)
  } else {
    NA_character_
  }
  band_caveat <- if (has_true_band) {
    .fdid_curve_band_caveat(method, band_method)
  } else {
    NA_character_
  }

  out <- data.frame(
    g = evg,
    estimate = estimate,
    se = se,
    pointwise_lower = pointwise_lower,
    pointwise_upper = pointwise_upper,
    band_lower = band_lower,
    band_upper = band_upper,
    curve = curve,
    curve_label = curve_label,
    method = method,
    vartype = vartype,
    target_estimand = target_estimand,
    inference_scope = inference_scope,
    inference_method = inference_method,
    signal_map = signal_map,
    density_method = density_method,
    pointwise_method = pointwise_method,
    pointwise_label = pointwise_label,
    band_method = band_method,
    band_label = band_label,
    band_caveat = band_caveat,
    practical_band_method = practical_band_method,
    gaussian_crit = gaussian_crit,
    simultaneous_crit = simultaneous_crit,
    curve_vcov_available = curve_vcov_available,
    band_zeta = band_zeta,
    band_coverage = band_coverage,
    n_boot_eff = n_boot_eff,
    ylab_default = ylab_default,
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$g), , drop = FALSE]
  if (!is.null(xlim)) {
    out <- out[out$g >= xlim[1] & out$g <= xlim[2], , drop = FALSE]
  }
  if (nrow(out) == 0L) {
    stop("No curve evaluation points fall within xlim.")
  }
  out <- out[order(out$g), , drop = FALSE]
  flags <- .fdid_curve_interval_flags(out, interval, ci, show.uniform.CI)
  if (!flags$pointwise) {
    out$pointwise_lower <- out$pointwise_upper <- NA_real_
  }
  if (!flags$uniform) {
    out$band_lower <- out$band_upper <- NA_real_
    out$band_label <- out$band_caveat <- NA_character_
  }
  out
}

.fdid_contiguous_segments <- function(data, cols, x_col = "g") {
  if (nrow(data) == 0L) return(list())
  cols <- c(x_col, cols)
  ok <- rep(TRUE, nrow(data))
  for (cc in cols) {
    if (!(cc %in% names(data))) return(list())
    ok <- ok & is.finite(data[[cc]])
  }
  valid_idx <- which(ok)
  if (length(valid_idx) == 0L) return(list())
  run_id <- cumsum(c(TRUE, diff(valid_idx) != 1L))
  lapply(split(valid_idx, run_id), function(ii) data[ii, , drop = FALSE])
}

.fdid_curve_interval_flags <- function(curve_df, interval, ci, show.uniform.CI) {
  has_pointwise <- any(is.finite(curve_df$pointwise_lower) &
                         is.finite(curve_df$pointwise_upper))
  has_uniform <- any(is.finite(curve_df$band_lower) &
                       is.finite(curve_df$band_upper))
  if (!isTRUE(ci) || identical(interval, "none")) {
    return(list(pointwise = FALSE, uniform = FALSE))
  }
  pointwise <- interval %in% c("auto", "pointwise", "both") && has_pointwise
  uniform <- interval %in% c("auto", "uniform", "both") &&
    has_uniform && isTRUE(show.uniform.CI)
  list(pointwise = pointwise, uniform = uniform)
}

.fdid_support_values <- function(x, plot_xlim) {
  g_support <- as.numeric(x$G)
  g_support <- g_support[is.finite(g_support)]
  g_support[g_support >= plot_xlim[1] & g_support <= plot_xlim[2]]
}

.fdid_support_available <- function(g_support, Xdistr) {
  if (identical(Xdistr, "none")) return(FALSE)
  if (identical(Xdistr, "rug")) return(length(g_support) >= 1L)
  length(g_support) >= 2L && length(unique(g_support)) >= 2L
}

.fdid_draw_support <- function(g_support, Xdistr, plot_xlim, bottom, top,
                               hist.color, density.color, support.alpha) {
  if (!.fdid_support_available(g_support, Xdistr)) return(invisible(NULL))
  if (Xdistr %in% c("histogram", "hist")) {
    h <- graphics::hist(g_support, breaks = 80, plot = FALSE)
    if (length(h$counts) > 0L && max(h$counts, na.rm = TRUE) > 0) {
      ymax <- bottom + h$counts / max(h$counts, na.rm = TRUE) * (top - bottom)
      graphics::rect(
        pmax(h$breaks[-length(h$breaks)], plot_xlim[1]),
        bottom,
        pmin(h$breaks[-1L], plot_xlim[2]),
        ymax,
        col = grDevices::adjustcolor(hist.color, support.alpha),
        border = grDevices::adjustcolor(hist.color, min(1, support.alpha + 0.3))
      )
    }
  }
  if (identical(Xdistr, "density")) {
    dens <- tryCatch(
      stats::density(g_support, from = plot_xlim[1], to = plot_xlim[2]),
      error = function(e) NULL
    )
    if (!is.null(dens) && length(dens$x) > 0L &&
        max(dens$y, na.rm = TRUE) > 0) {
      dens_y <- bottom + dens$y / max(dens$y, na.rm = TRUE) * (top - bottom)
      graphics::polygon(
        c(dens$x, rev(dens$x)),
        c(dens_y, rep(bottom, length(dens_y))),
        col = grDevices::adjustcolor(density.color, support.alpha),
        border = grDevices::adjustcolor(density.color, min(1, support.alpha + 0.3))
      )
    }
  }
  if (identical(Xdistr, "rug")) {
    graphics::rug(g_support, side = 1,
                  col = grDevices::adjustcolor(hist.color, support.alpha),
                  ticksize = 0.04)
  }
  invisible(NULL)
}

.fdid_plot_support_panel <- function(g_support, Xdistr, plot_xlim, xlab,
                                     hist.color, density.color, support.alpha,
                                     eval_g = NULL, eval.color = "steelblue") {
  graphics::plot(NA, xlim = plot_xlim, ylim = c(0, 1), type = "n",
                 axes = FALSE, xlab = xlab, ylab = "")
  graphics::axis(1)
  graphics::box()
  graphics::mtext("Observed G", side = 2, line = 2.2, cex = 0.8)
  .fdid_draw_support(g_support, Xdistr, plot_xlim, 0, 1,
                     hist.color, density.color, support.alpha)
  if (!is.null(eval_g) && length(eval_g) > 0L) {
    graphics::rug(eval_g, side = 3,
                  col = grDevices::adjustcolor(eval.color, 0.35),
                  ticksize = 0.08)
  }
  invisible(NULL)
}

.fdid_curve_marker_values <- function(diff.values = NULL, contrast.values = NULL) {
  vals <- c(diff.values, contrast.values)
  vals <- suppressWarnings(as.numeric(vals))
  vals[is.finite(vals)]
}

.fdid_unique_numeric <- function(x) {
  out <- numeric(0L)
  for (val in x) {
    if (!any(abs(out - val) <= sqrt(.Machine$double.eps))) {
      out <- c(out, val)
    }
  }
  out
}

.fdid_format_g <- function(x) {
  format(signif(as.numeric(x), 4L), trim = TRUE)
}

.fdid_legend_margin_lines <- function(labels) {
  labels <- as.character(labels)
  max_label <- if (length(labels) > 0L) max(nchar(labels), na.rm = TRUE) else 0L
  label_factor <- if (is.finite(max_label) && max_label > 24L) 1.1 else 0.95
  max(4.4, 2.2 + label_factor * max(1L, length(labels)))
}

.fdid_mar_with_bottom_legend <- function(mar, labels) {
  mar[1L] <- max(mar[1L], mar[1L] + .fdid_legend_margin_lines(labels))
  mar
}

.fdid_draw_bottom_legend <- function(legend, col = NULL, lty = NULL,
                                     lwd = NULL, pch = NULL,
                                     fill = NULL, pt.cex = NULL,
                                     cex = 0.82) {
  usr <- graphics::par("usr")
  x_mid <- mean(usr[1:2])
  y_gap <- 0.34 * diff(usr[3:4])
  old_xpd <- graphics::par("xpd")
  graphics::par(xpd = NA)
  on.exit(graphics::par(xpd = old_xpd), add = TRUE)
  legend_args <- list(
    x = x_mid,
    y = usr[3] - y_gap,
    legend = legend,
    xjust = 0.5,
    yjust = 1,
    bg = "white",
    box.col = "gray70",
    cex = cex
  )
  if (!is.null(col)) legend_args$col <- col
  if (!is.null(lty)) legend_args$lty <- lty
  if (!is.null(lwd)) legend_args$lwd <- lwd
  if (!is.null(pch)) legend_args$pch <- pch
  if (!is.null(fill)) legend_args$fill <- fill
  if (!is.null(pt.cex)) legend_args$pt.cex <- pt.cex
  do.call(graphics::legend, legend_args)
  invisible(NULL)
}

.fdid_draw_legend_panel <- function(legend, col = NULL, lty = NULL,
                                    lwd = NULL, pch = NULL,
                                    fill = NULL, pt.cex = NULL,
                                    cex = 0.82) {
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
  legend_args <- list(
    x = 0.5,
    y = 0.55,
    legend = legend,
    xjust = 0.5,
    yjust = 0.5,
    bg = "white",
    box.col = "gray70",
    cex = cex
  )
  if (!is.null(col)) legend_args$col <- col
  if (!is.null(lty)) legend_args$lty <- lty
  if (!is.null(lwd)) legend_args$lwd <- lwd
  if (!is.null(pch)) legend_args$pch <- pch
  if (!is.null(fill)) legend_args$fill <- fill
  if (!is.null(pt.cex)) legend_args$pt.cex <- pt.cex
  do.call(graphics::legend, legend_args)
  invisible(NULL)
}

.fdid_nearest_grid <- function(evg, value) {
  evg[which.min(abs(evg - value))]
}

.fdid_resolve_ref_g <- function(ref.g, evg, targets = NULL,
                                target.type = "level") {
  if (is.null(ref.g)) {
    if (identical(target.type, "quantile") && length(targets) > 0L) {
      return(as.numeric(targets[1L]))
    }
    return(.fdid_nearest_grid(evg, stats::median(evg)))
  }
  if (is.character(ref.g)) {
    if (length(ref.g) != 1L || is.na(ref.g) || !nzchar(ref.g)) {
      stop("ref.g must be a single reference rule or a single finite number.")
    }
    rule <- tolower(ref.g)
    if (rule %in% c("median", "med")) {
      return(.fdid_nearest_grid(evg, stats::median(evg)))
    }
    if (rule %in% c("mean", "average", "avg")) {
      return(.fdid_nearest_grid(evg, mean(evg)))
    }
    if (rule %in% c("min", "minimum", "lowest")) {
      return(min(evg))
    }
    if (rule %in% c("max", "maximum", "highest")) {
      return(max(evg))
    }
    if (rule %in% c("first", "q1")) {
      if (length(targets) == 0L) {
        stop("ref.g='first' requires resolved target values.")
      }
      return(as.numeric(targets[1L]))
    }
    stop("Unknown ref.g rule. Use a finite number, 'median', 'mean', ",
         "'min', 'max', or 'first'.")
  }
  if (!is.numeric(ref.g) || length(ref.g) != 1L || !is.finite(ref.g)) {
    stop("ref.g must be a single reference rule or a single finite number.")
  }
  as.numeric(ref.g)
}

.fdid_contrast_quantile_source <- function(x, evg,
                                           quantile.source = "support") {
  if (identical(quantile.source, "support")) {
    g_support <- suppressWarnings(as.numeric(x$G))
    g_support <- g_support[is.finite(g_support)]
    g_support <- g_support[g_support >= min(evg) & g_support <= max(evg)]
    if (length(unique(g_support)) >= 2L) {
      return(g_support)
    }
  }
  evg
}

.fdid_contrast_quantile_probs <- function(spec, n.contrast) {
  if (is.null(spec)) {
    spec <- n.contrast
  }
  if (is.character(spec)) {
    if (length(spec) != 1L || is.na(spec)) {
      stop("Quantile target shorthand must be a single value such as 'q5' or '5'.")
    }
    low <- tolower(spec)
    if (grepl("^q[0-9]+$", low)) {
      spec <- as.numeric(sub("^q", "", low))
    } else if (grepl("^[0-9]+$", low)) {
      spec <- as.numeric(low)
    } else {
      stop("Unknown quantile target shorthand. Use target.g='q5', target.g='5', or ",
           "target.type='quantile', target.g=5.")
    }
  }
  if (!is.numeric(spec) || length(spec) == 0L || any(!is.finite(spec))) {
    stop("Quantile targets must be a finite numeric vector or shorthand like 'q5' or '5'.")
  }
  if (length(spec) == 1L && spec > 1) {
    k <- as.integer(floor(spec))
    if (k < 2L || abs(k - spec) > sqrt(.Machine$double.eps)) {
      stop("A scalar quantile target count must be an integer greater than 1.")
    }
    return(seq(0, 1, length.out = k))
  }
  if (all(spec >= 0 & spec <= 1)) {
    probs <- spec
  } else if (all(spec >= 1 & spec <= 100)) {
    probs <- spec / 100
  } else {
    stop("Quantile probabilities must be in [0, 1], percent values in ",
         "[1, 100], or a scalar count greater than 1.")
  }
  sort(.fdid_unique_numeric(as.numeric(probs)))
}

.fdid_resolve_contrast_targets <- function(x, evg, supplied_targets,
                                           n.contrast, max.contrast,
                                           include.ref,
                                           target.type = "level",
                                           quantile.source = "support") {
  if (identical(target.type, "quantile")) {
    probs <- .fdid_contrast_quantile_probs(supplied_targets, n.contrast)
    q_source <- .fdid_contrast_quantile_source(x, evg, quantile.source)
    vals <- as.numeric(stats::quantile(q_source, probs = probs,
                                       names = FALSE, na.rm = TRUE,
                                       type = 7))
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0L) {
      stop("No finite quantile target values are available.")
    }
    vals_unique <- .fdid_unique_numeric(vals)
    if (length(vals_unique) < length(vals)) {
      warning("Duplicate quantile target values were removed.", call. = FALSE)
    }
    labels <- paste0("Q", seq_along(vals_unique), "\n",
                     .fdid_format_g(vals_unique))
    return(list(values = vals_unique, labels = labels,
                supplied = !is.null(supplied_targets), auto = FALSE,
                target_type = "quantile"))
  }

  if (is.null(supplied_targets)) {
    list(values = numeric(0L), labels = character(0L), supplied = FALSE,
         auto = TRUE, target_type = "level")
  } else {
    if (!is.numeric(supplied_targets) || length(supplied_targets) == 0L ||
        any(!is.finite(supplied_targets))) {
      stop("target.g and contrast.values must be finite numeric vectors.")
    }
    vals <- .fdid_unique_numeric(as.numeric(supplied_targets))
    list(values = vals, labels = .fdid_format_g(vals), supplied = TRUE,
         auto = FALSE, target_type = "level")
  }
}

.fdid_reference_contrast_row <- function(x, ref.g, target_estimand, label) {
  out <- .fdid_result_row(
    target = "level contrast",
    method = x$method,
    curve = "level",
    g0 = ref.g,
    g1 = ref.g,
    estimate = 0,
    std_error = 0,
    conf_low = 0,
    conf_high = 0,
    inference = "reference",
    note = "Reference row; contrast of fitted curve value with itself."
  )
  out$target_g <- as.numeric(ref.g)
  out$ref_g <- as.numeric(ref.g)
  out$reference_g <- as.numeric(ref.g)
  out$is_reference <- TRUE
  out$target_estimand <- target_estimand
  out$label <- label
  out
}

.fdid_contrast_plot_data <- function(x, ref.g = NULL, contrast.values = NULL,
                                     target.g = NULL, n.contrast = 7L,
                                     max.contrast = 10L,
                                     include.ref = TRUE,
                                     target.type = c("level", "quantile"),
                                     quantile.source = c("support", "eval_g"),
                                     inference = c("auto", "replicate", "vcov",
                                                   "band", "pointwise"),
                                     interpolate = c("linear", "none"),
                                     alpha = NULL) {
  inference <- match.arg(inference)
  interpolate <- match.arg(interpolate)
  target.type <- match.arg(target.type)
  quantile.source <- match.arg(quantile.source)
  if (!is.logical(include.ref) || length(include.ref) != 1L ||
      is.na(include.ref)) {
    stop("include.ref must be TRUE or FALSE.")
  }
  if (!(x$method %in% c("kernel", "dml_flex"))) {
    stop("contrast plots are available for method='kernel' and method='dml_flex' only.")
  }
  if (is.null(x$eval_g) || is.null(x$curve_event)) {
    stop("contrast plots require stored continuous-G curve data.")
  }
  if (!is.null(target.g) && !is.null(contrast.values)) {
    stop("Supply only one of target.g or contrast.values for contrast plots.")
  }
  if (identical(target.type, "level") && is.character(target.g) &&
      length(target.g) == 1L &&
      grepl("^(q[0-9]+|[0-9]+)$", tolower(target.g))) {
    target.type <- "quantile"
  }
  if (identical(target.type, "quantile") && !is.null(contrast.values)) {
    stop("contrast.values supplies level targets; use target.g with ",
         "target.type='quantile'.")
  }
  if (!is.numeric(n.contrast) || length(n.contrast) != 1L ||
      !is.finite(n.contrast) || n.contrast < 1) {
    stop("n.contrast must be a positive finite number.")
  }
  if (!is.numeric(max.contrast) || length(max.contrast) != 1L ||
      !is.finite(max.contrast) || max.contrast < 1) {
    stop("max.contrast must be a positive finite number.")
  }
  n.contrast <- as.integer(floor(n.contrast))
  max.contrast <- as.integer(floor(max.contrast))
  display_cap <- min(n.contrast, max.contrast)

  evg <- as.numeric(x$eval_g)
  evg <- evg[is.finite(evg)]
  if (length(evg) == 0L) {
    stop("No finite eval_g values are available for contrast plotting.")
  }
  evg <- sort(.fdid_unique_numeric(evg))

  supplied_targets <- target.g %||% contrast.values
  targets <- .fdid_resolve_contrast_targets(
    x = x,
    evg = evg,
    supplied_targets = supplied_targets,
    n.contrast = display_cap,
    max.contrast = max.contrast,
    include.ref = include.ref,
    target.type = target.type,
    quantile.source = quantile.source
  )
  ref.g <- .fdid_resolve_ref_g(ref.g, evg, targets$values, target.type)

  if (targets$auto) {
    candidates <- evg[abs(evg - ref.g) > sqrt(.Machine$double.eps)]
    n_nonref <- if (isTRUE(include.ref)) display_cap - 1L else display_cap
    n_nonref <- min(max(0L, n_nonref), length(candidates))
    if (length(candidates) == 0L && !isTRUE(include.ref)) {
      stop("No target eval_g values remain after excluding the reference.")
    }
    idx <- if (n_nonref > 0L) {
      unique(round(seq(1, length(candidates), length.out = n_nonref)))
    } else {
      integer(0L)
    }
    target.g <- candidates[idx]
    target_labels <- .fdid_format_g(target.g)
  } else {
    target.g <- targets$values
    target_labels <- targets$labels
  }

  requested_ref <- abs(target.g - ref.g) <= sqrt(.Machine$double.eps)
  add_ref <- isTRUE(include.ref) || any(requested_ref)
  target.g_nonref <- target.g[!requested_ref]
  target_labels_nonref <- target_labels[!requested_ref]
  total_rows <- length(target.g_nonref) + as.integer(add_ref)
  if (!targets$auto && total_rows > max.contrast) {
    warning("More than max.contrast contrast rows supplied; plotting all rows.",
            call. = FALSE)
  }

  rows <- lapply(target.g_nonref, function(g1) {
    fdid_contrast(x, g0 = ref.g, g1 = g1, type = "level",
                  inference = inference, interpolate = interpolate,
                  alpha = alpha)
  })
  target_estimand <- if (identical(x$method, "kernel")) {
    "mu(g) - mu(g_ref)"
  } else {
    "theta(g) - theta(g_ref)"
  }
  if (length(rows) > 0L) {
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out$target_g <- as.numeric(out$g1)
    out$ref_g <- as.numeric(out$g0)
    out$reference_g <- out$ref_g
    out$is_reference <- FALSE
    out$target_estimand <- target_estimand
    out$label <- target_labels_nonref
  } else {
    out <- NULL
  }
  if (add_ref) {
    ref_label <- if (identical(target.type, "quantile") && any(requested_ref)) {
      target_labels[which.min(abs(target.g - ref.g))]
    } else {
      paste0("ref\n", .fdid_format_g(ref.g))
    }
    ref_row <- .fdid_reference_contrast_row(x, ref.g, target_estimand, ref_label)
    out <- if (is.null(out)) ref_row else rbind(ref_row, out)
  }
  if (is.null(out) || nrow(out) == 0L) {
    stop("No finite contrasts are available to plot.")
  }
  ord <- order(out$target_g, !out$is_reference)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}

.fdid_plot_contrast_base <- function(contrast_df, xlab = NULL, ylab = NULL,
                                     main = NULL, ylim = NULL,
                                     line.color = "steelblue",
                                     ci.color = "steelblue",
                                     ci = TRUE,
                                     show.grid = FALSE, ...) {
  if (nrow(contrast_df) == 0L ||
      !any(is.finite(contrast_df$estimate))) {
    stop("No finite contrasts are available to plot.")
  }
  xpos <- seq_len(nrow(contrast_df))
  all_y <- contrast_df$estimate
  nonref_ci <- if ("is_reference" %in% names(contrast_df)) {
    !contrast_df$is_reference
  } else {
    rep(TRUE, nrow(contrast_df))
  }
  has_ci <- isTRUE(ci) & is.finite(contrast_df$conf.low) &
    is.finite(contrast_df$conf.high) & nonref_ci
  if (any(has_ci)) {
    all_y <- c(all_y, contrast_df$conf.low[has_ci], contrast_df$conf.high[has_ci])
  }
  all_y <- all_y[is.finite(all_y)]
  if (is.null(ylim)) {
    rng <- range(c(all_y, 0), na.rm = TRUE)
    pad <- if (is.finite(diff(rng)) && diff(rng) > 0) 0.08 * diff(rng) else 0.5
    ylim <- c(rng[1] - pad, rng[2] + pad)
  }
  label_lines <- max(lengths(strsplit(as.character(contrast_df$label), "\n",
                                      fixed = TRUE)), 1L)
  old_mar <- graphics::par("mar")
  needed_bottom <- max(4.8, 3.8 + 0.65 * label_lines)
  if (old_mar[1L] < needed_bottom) {
    new_mar <- old_mar
    new_mar[1L] <- needed_bottom
    graphics::par(mar = new_mar)
    on.exit(graphics::par(mar = old_mar), add = TRUE)
  }
  graphics::plot(
    xpos, contrast_df$estimate, type = "n", xaxt = "n",
    xlim = c(0.5, length(xpos) + 0.5), ylim = ylim,
    xlab = xlab %||% "Target G",
    ylab = ylab %||% contrast_df$target_estimand[1L],
    main = main %||% paste0("Level contrasts vs G = ",
                            format(signif(contrast_df$reference_g[1L], 4L),
                                   trim = TRUE)),
    ...
  )
  if (isTRUE(show.grid)) {
    graphics::grid(col = "gray90", lty = 1)
    graphics::box()
  }
  graphics::abline(h = 0, col = "gray50", lwd = 1.5, lty = 2)
  if (any(has_ci)) {
    graphics::arrows(xpos[has_ci], contrast_df$conf.low[has_ci],
                     xpos[has_ci], contrast_df$conf.high[has_ci],
                     angle = 90, code = 3, length = 0.05,
                     col = ci.color)
  }
  ref_row <- "is_reference" %in% names(contrast_df) &
    isTRUE(any(contrast_df$is_reference))
  if (ref_row) {
    nonref <- !contrast_df$is_reference
    if (any(nonref)) {
      graphics::points(xpos[nonref], contrast_df$estimate[nonref], pch = 16,
                       col = line.color, cex = 1.05)
    }
    graphics::points(xpos[contrast_df$is_reference],
                     contrast_df$estimate[contrast_df$is_reference],
                     pch = 21, bg = "white", col = "gray20", cex = 1.1)
  } else {
    graphics::points(xpos, contrast_df$estimate, pch = 16,
                     col = line.color, cex = 1.05)
  }
  graphics::axis(1, at = xpos, labels = FALSE, tcl = -0.25)
  label_cex <- if (nrow(contrast_df) > 7L || label_lines > 1L) 0.82 else 0.9
  graphics::axis(1, at = xpos, labels = contrast_df$label,
                 tick = FALSE, line = 0.35, cex.axis = label_cex)
  graphics::box()
  inf_rows <- if ("is_reference" %in% names(contrast_df)) {
    !contrast_df$is_reference
  } else {
    rep(TRUE, nrow(contrast_df))
  }
  inf <- unique(stats::na.omit(contrast_df$inference[inf_rows]))
  if (length(inf) == 1L) {
    graphics::mtext(paste("Inference:", inf), side = 3, line = 0.2,
                    adj = 1, cex = 0.75, col = "gray35")
  } else if (length(inf) > 1L) {
    graphics::mtext("Inference varies by target; see fdid_contrast()",
                    side = 3, line = 0.2, adj = 1, cex = 0.75,
                    col = "gray35")
  }
  invisible(NULL)
}

.fdid_auto_curve_legend_position <- function(curve_df, plot_xlim, plot_ylim,
                                             embedded_support = FALSE) {
  candidates <- c("topright", "topleft", "bottomright", "bottomleft")
  xmid <- mean(plot_xlim)
  ymid <- mean(plot_ylim)
  score <- setNames(rep(0, length(candidates)), candidates)
  y_cols <- c("estimate", "pointwise_lower", "pointwise_upper",
              "band_lower", "band_upper")
  y_cols <- y_cols[y_cols %in% names(curve_df)]
  pts <- data.frame(
    x = rep(curve_df$g, length(y_cols)),
    y = unlist(curve_df[y_cols], use.names = FALSE)
  )
  pts <- pts[is.finite(pts$x) & is.finite(pts$y), , drop = FALSE]
  if (nrow(pts) == 0L) return("topright")
  for (pos in candidates) {
    right <- grepl("right", pos)
    top <- grepl("top", pos)
    in_x <- if (right) pts$x >= xmid else pts$x <= xmid
    in_y <- if (top) pts$y >= ymid else pts$y <= ymid
    score[pos] <- sum(in_x & in_y)
    if (embedded_support && !top) score[pos] <- score[pos] + 100
  }
  candidates[which.min(score[candidates])]
}

.fdid_plot_curve_base <- function(x, curve_df, interval, ci, show.uniform.CI,
                                  Xdistr, support.panel, show.eval_g,
                                  diff.values, contrast.values,
                                  xlab, ylab, main, ylim, xlim,
                                  line.color, line.size, ci.color, ci.alpha,
                                  band.color, band.lty, hist.color,
                                  density.color, support.alpha,
                                  show.grid = FALSE,
                                  legend.position = "auto", ...) {
  flags <- .fdid_curve_interval_flags(curve_df, interval, ci, show.uniform.CI)
  finite_g <- curve_df$g[is.finite(curve_df$g)]
  finite_y <- curve_df$estimate[is.finite(curve_df$estimate)]
  if (length(finite_y) == 0L) {
    stop("The selected curve has no finite estimates to plot.")
  }
  plot_xlim <- xlim %||% range(finite_g, na.rm = TRUE)
  if (!all(is.finite(plot_xlim)) || plot_xlim[1] == plot_xlim[2]) {
    pad <- if (is.finite(plot_xlim[1])) max(abs(plot_xlim[1]) * 0.04, 0.5) else 0.5
    plot_xlim <- c(plot_xlim[1] - pad, plot_xlim[2] + pad)
  }

  all_y <- finite_y
  if (flags$pointwise) {
    all_y <- c(all_y, curve_df$pointwise_lower, curve_df$pointwise_upper)
  }
  if (flags$uniform) {
    all_y <- c(all_y, curve_df$band_lower, curve_df$band_upper)
  }
  all_y <- all_y[is.finite(all_y)]
  if (length(all_y) == 0L) all_y <- finite_y
  if (is.null(ylim)) {
    rng <- range(all_y, na.rm = TRUE)
    pad <- if (is.finite(diff(rng)) && diff(rng) > 0) 0.05 * diff(rng) else 0.5
    curve_ylim <- c(rng[1] - pad, rng[2] + pad)
  } else {
    curve_ylim <- ylim
  }

  g_support <- .fdid_support_values(x, plot_xlim)
  support_available <- .fdid_support_available(g_support, Xdistr)
  separate_support <- identical(support.panel, "separate") &&
    support_available && !identical(Xdistr, "none")
  embedded_support <- identical(support.panel, "embedded") &&
    support_available && !identical(Xdistr, "rug")

  legend_labels <- "Estimate"
  legend_col <- line.color
  legend_lty <- 1
  legend_lwd <- line.size
  legend_pch <- NA
  legend_pt_cex <- 1
  if (flags$pointwise) {
    point_idx <- which(!is.na(curve_df$pointwise_label) &
                         nzchar(curve_df$pointwise_label))
    point_label <- if (length(point_idx) > 0L) {
      curve_df$pointwise_label[point_idx[1L]]
    } else {
      "Pointwise CI"
    }
    legend_labels <- c(legend_labels, point_label)
    legend_col <- c(legend_col, grDevices::adjustcolor(ci.color, ci.alpha))
    legend_lty <- c(legend_lty, NA)
    legend_lwd <- c(legend_lwd, NA)
    legend_pch <- c(legend_pch, 15)
    legend_pt_cex <- c(legend_pt_cex, 1.8)
  }
  if (flags$uniform) {
    band_idx <- which(!is.na(curve_df$band_label) &
                        nzchar(curve_df$band_label))
    band_label <- if (length(band_idx) > 0L) {
      curve_df$band_label[band_idx[1L]]
    } else {
      "Uniform band"
    }
    legend_labels <- c(legend_labels, band_label)
    legend_col <- c(legend_col, band.color)
    legend_lty <- c(legend_lty, band.lty)
    legend_lwd <- c(legend_lwd, 1.2)
    legend_pch <- c(legend_pch, NA)
    legend_pt_cex <- c(legend_pt_cex, 1)
  }
  draw_legend <- !identical(legend.position, "none") &&
    (length(legend_labels) > 1L || !identical(legend.position, "auto"))
  legend_below <- draw_legend && identical(legend.position, "auto")

  if (separate_support) {
    oldpar <- graphics::par(no.readonly = TRUE)
    if (legend_below) {
      graphics::layout(matrix(c(1, 2, 3), ncol = 1), heights = c(4, 1, 0.85))
    } else {
      graphics::layout(matrix(c(1, 2), ncol = 1), heights = c(4, 1))
    }
    graphics::par(mar = c(2.2, 3.8, 2.2, 0.8))
    on.exit({
      try(graphics::layout(1), silent = TRUE)
      try(graphics::par(oldpar), silent = TRUE)
    }, add = TRUE)
  } else if (legend_below) {
    old_mar <- graphics::par("mar")
    graphics::par(mar = .fdid_mar_with_bottom_legend(old_mar, legend_labels))
    on.exit(graphics::par(mar = old_mar), add = TRUE)
  }

  plot_ylim <- curve_ylim
  if (embedded_support && is.null(ylim)) {
    y_span <- diff(curve_ylim)
    if (!is.finite(y_span) || y_span <= 0) y_span <- 1
    plot_ylim <- c(curve_ylim[1] - y_span / 5, curve_ylim[2])
  }

  curve_xlab <- if (separate_support) "" else xlab %||% "Baseline factor G"
  graphics::plot(curve_df$g, curve_df$estimate, type = "n",
                 xlim = plot_xlim, ylim = plot_ylim,
                 xlab = curve_xlab,
                 ylab = ylab %||% curve_df$ylab_default[1L],
                 main = main, ...)
  if (isTRUE(show.grid)) {
    graphics::grid(col = "gray90", lty = 1)
    graphics::box()
  }

  if (embedded_support) {
    y_span <- diff(curve_ylim)
    if (!is.finite(y_span) || y_span <= 0) y_span <- 1
    support_bottom <- if (is.null(ylim)) plot_ylim[1] else curve_ylim[1]
    support_top <- if (is.null(ylim)) curve_ylim[1] else curve_ylim[1] + y_span / 5
    .fdid_draw_support(g_support, Xdistr, plot_xlim, support_bottom, support_top,
                       hist.color, density.color, support.alpha)
  }
  if (identical(support.panel, "embedded") && identical(Xdistr, "rug") &&
      support_available) {
    .fdid_draw_support(g_support, Xdistr, plot_xlim, curve_ylim[1],
                       curve_ylim[1], hist.color, density.color, support.alpha)
  }

  graphics::abline(h = 0, col = "gray50", lwd = 1.5, lty = 2)
  markers <- .fdid_curve_marker_values(diff.values, contrast.values)
  markers <- markers[markers >= plot_xlim[1] & markers <= plot_xlim[2]]
  if (length(markers) > 0L) {
    graphics::abline(v = markers,
                     col = grDevices::adjustcolor("firebrick3", 0.75),
                     lwd = 1, lty = 2)
  }

  if (flags$pointwise) {
    ci_segments <- .fdid_contiguous_segments(curve_df,
                                             c("pointwise_lower", "pointwise_upper"))
    for (seg in ci_segments) {
      if (nrow(seg) >= 2L) {
        graphics::polygon(
          c(seg$g, rev(seg$g)),
          c(seg$pointwise_lower, rev(seg$pointwise_upper)),
          col = grDevices::adjustcolor(ci.color, ci.alpha),
          border = NA
        )
      }
    }
  }
  if (flags$uniform) {
    band_segments <- .fdid_contiguous_segments(curve_df,
                                               c("band_lower", "band_upper"))
    for (seg in band_segments) {
      if (nrow(seg) >= 2L) {
        graphics::lines(seg$g, seg$band_lower,
                        col = band.color, lwd = 1.2, lty = band.lty)
        graphics::lines(seg$g, seg$band_upper,
                        col = band.color, lwd = 1.2, lty = band.lty)
      }
    }
  }

  estimate_segments <- .fdid_contiguous_segments(curve_df, "estimate")
  for (seg in estimate_segments) {
    if (nrow(seg) >= 2L) {
      graphics::lines(seg$g, seg$estimate, lwd = line.size, col = line.color)
    } else {
      graphics::points(seg$g, seg$estimate, pch = 16, col = line.color)
    }
  }

  if (isTRUE(show.eval_g)) {
    graphics::rug(curve_df$g, side = 1,
                  col = grDevices::adjustcolor(line.color, 0.35),
                  ticksize = 0.025)
  }

  if (draw_legend) {
    legend_cex <- if (max(nchar(legend_labels), na.rm = TRUE) > 18L ||
                      length(legend_labels) > 2L) 0.78 else 0.85
    if (legend_below && !separate_support) {
      .fdid_draw_bottom_legend(
        legend = legend_labels,
        col = legend_col,
        lty = legend_lty,
        lwd = legend_lwd,
        pch = legend_pch,
        pt.cex = legend_pt_cex,
        cex = legend_cex
      )
    } else if (!legend_below) {
      graphics::legend(legend.position, legend = legend_labels, col = legend_col,
                       lty = legend_lty, lwd = legend_lwd, pch = legend_pch,
                       pt.cex = legend_pt_cex, bty = "n", cex = legend_cex)
    }
  }

  if (separate_support) {
    graphics::par(mar = c(2.8, 3.8, 0.2, 0.8))
    .fdid_plot_support_panel(
      g_support, Xdistr, plot_xlim, xlab %||% "Baseline factor G",
      hist.color, density.color, support.alpha,
      eval_g = if (isTRUE(show.eval_g)) curve_df$g else NULL,
      eval.color = line.color
    )
    if (legend_below && draw_legend) {
      legend_cex <- if (max(nchar(legend_labels), na.rm = TRUE) > 18L ||
                        length(legend_labels) > 2L) 0.78 else 0.85
      graphics::par(mar = c(0, 0, 0, 0))
      .fdid_draw_legend_panel(
        legend = legend_labels,
        col = legend_col,
        lty = legend_lty,
        lwd = legend_lwd,
        pch = legend_pch,
        pt.cex = legend_pt_cex,
        cex = legend_cex
      )
    }
  }

  invisible(NULL)
}

#' Plot Results from FDID Analysis
#'
#' Provides visualisations for FDID results, including raw means, dynamic
#' effects, propensity-score overlap, continuous-G curves, and continuous-G
#' level contrasts. The comparison plot of multiple methods has been removed;
#' use \code{plot.fdid_list()} for that.
#'
#' @param x   An \code{fdid} object.
#' @param type One of \code{"raw"}, \code{"dynamic"}, \code{"overlap"},
#'   \code{"curve"}, or \code{"contrast"}. The \code{"curve"} and
#'   \code{"contrast"} types are available for \code{method="kernel"} and
#'   \code{method="dml_flex"} objects with stored curve data.
#' @param curve Which continuous-G curve to plot when \code{type="curve"}.
#'   Options are \code{"auto"}, \code{"level"}, \code{"derivative"}, and
#'   \code{"both"}.
#'   \code{"auto"} preserves method defaults: derivative for \code{"kernel"}
#'   and level for \code{"dml_flex"}.
#' @param connected Logical; if \code{TRUE}, connects points with lines in
#'        the "raw" and "dynamic" plots.  Default is \code{FALSE}.
#' @param ci Logical; if \code{TRUE}, draw 95\% CIs when available. Default is
#'   \code{TRUE}. For continuous-G curves, use \code{interval} to choose
#'   pointwise intervals, simultaneous bands, both, or neither. For
#'   continuous-G contrast plots, controls coefficient CI bars.
#' @param interval For continuous-G curve plots, one of \code{"auto"},
#'   \code{"pointwise"}, \code{"uniform"}, \code{"both"}, or \code{"none"}.
#'   \code{"auto"} shows pointwise intervals when available and stored
#'   simultaneous bands when available. Kernel bootstrap bands are finite-grid
#'   bootstrap quantile-envelope bands. Kernel robust pointwise normal intervals
#'   are not shown as uniform bands. For \code{method="dml_flex"}, stored bands
#'   are labeled as practical signal-mapping bands or, for
#'   \code{signal_map="blp_spline"}, finite-grid Gaussian max-t BLP-spline
#'   covariance bands.
#' @param show.uniform.CI Logical; if \code{FALSE}, suppress stored
#'   simultaneous bands in continuous-G curve plots.
#' @param shade_periods Shaded intervals on the time axis. Default uses \code{x$tr_period}, i.e. event periods.
#'        Set to \code{NULL} to remove shaded area.
#' @param alpha_shade Transparency for shading the treatment period.
#' @param palette A palette name from \strong{RColorBrewer}.  Default \code{"Set2"}.
#' @param group_labels Labels for the two groups.
#' @param xlab,ylab,main Axis labels and main title.
#' @param ylim Y-axis limits.  Default \code{NULL} (computed automatically).
#' @param xlim X-axis limits. For \code{type="curve"}, curve rows and support
#'   overlays are restricted to this range before drawing.
#' @param Xdistr How to display the observed support of continuous \code{G}
#'   under a curve plot. One of \code{"histogram"}, \code{"hist"},
#'   \code{"density"}, \code{"rug"}, or \code{"none"}. Default is
#'   \code{"histogram"}.
#' @param support.panel For continuous-G curve plots, \code{"embedded"} draws
#'   the support display in a lower strip of the main plot, while
#'   \code{"separate"} draws it in a small aligned panel below the curve. When
#'   \code{curve="both"}, level and derivative panels are drawn side-by-side
#'   and support displays are embedded in each panel.
#' @param show.eval_g Logical; if \code{TRUE}, draw ticks for the stored
#'   evaluation grid in continuous-G curve plots.
#' @param line.color,line.size,ci.color,ci.alpha,band.color,band.lty,hist.color,density.color,support.alpha
#'   Conservative base-graphics style controls for continuous-G curve plots.
#'   \code{line.color} and \code{ci.color} also style continuous-G contrast
#'   coefficient plots.
#' @param diff.values,contrast.values Optional numeric \code{G} values to mark
#'   with vertical reference lines in continuous-G curve plots. For
#'   \code{type="contrast"}, \code{contrast.values} supplies target
#'   \code{G} values and cannot be combined with \code{target.g}.
#' @param ref.g Reference \code{G} value for \code{type="contrast"}. May be a
#'   finite number or one of \code{"median"}, \code{"mean"}, \code{"min"},
#'   \code{"max"}, or \code{"first"}. For level targets, the default
#'   \code{NULL} is equivalent to \code{"median"}. For quantile targets, the
#'   default uses the first selected quantile target.
#' @param target.g Target \code{G} values for \code{type="contrast"}. With
#'   \code{target.type="level"}, numeric values are interpreted as absolute
#'   \code{G} levels. With \code{target.type="quantile"}, a scalar integer
#'   requests that many quantile points, a vector in \code{[0,1]} supplies
#'   quantile probabilities, and shorthand such as \code{"q5"} or \code{"5"}
#'   requests five quantile points.
#' @param n.contrast,max.contrast Automatic contrast target count and cap.
#'   In automatic and quantile modes these counts include the reference row
#'   when \code{include.ref=TRUE}. User-supplied targets above
#'   \code{max.contrast} are allowed with a warning.
#' @param include.ref Logical; if \code{TRUE}, include a reference row with
#'   coefficient zero in \code{type="contrast"} plots. Default is \code{TRUE}.
#' @param target.type For \code{type="contrast"}, whether \code{target.g}
#'   denotes absolute \code{G} levels or quantile targets.
#' @param quantile.source For quantile contrast targets, use observed
#'   continuous \code{G} support within the fitted \code{eval_g} range
#'   (\code{"support"}) or the fitted \code{eval_g} grid itself
#'   (\code{"eval_g"}).
#' @param inference Inference rule passed to \code{\link{fdid_contrast}} for
#'   \code{type="contrast"}.
#' @param interpolate Interpolation rule passed to \code{\link{fdid_contrast}}
#'   for \code{type="contrast"}.
#' @param show.grid Logical; if \code{TRUE}, draw a light background grid for
#'   continuous-G curve and contrast plots.
#' @param legend.position Legend placement for plots that draw a legend. Use
#'   \code{"auto"} to draw a boxed legend below the plot, one of the standard
#'   base graphics corners for inside-plot placement, or \code{"none"} to
#'   suppress the legend.
#' @param theme.bw Reserved for compatibility with ggplot-style plotting
#'   controls; currently ignored by the base graphics method.
#' @param ...  Additional graphics parameters.
#'
#' @return Produces a plot; invisibly returns \code{NULL}.
#' @examples
#' \donttest{
#' data(fdid)
#' mortality$uniqueid <- paste(mortality$provid, mortality$countyid, sep = "-")
#' mortality$G <- ifelse(mortality$pczupu >= median(mortality$pczupu, na.rm = TRUE), 1, 0)
#' s <- fdid_prepare(
#'   data = mortality, Y_label = "mortality",
#'   X_labels = c("avggrain", "lnpop"),
#'   G_label = "G", unit_label = "uniqueid", time_label = "year"
#' )
#' result <- fdid(s, tr_period = 1958:1961, ref_period = 1957)
#' plot(result, type = "raw")
#' plot(result, type = "dynamic")
#' }
#' @author Rivka Lipkovitz, Enhan Liu
#' @importFrom rlang %||%
#' @export
plot.fdid <- function(x,
                      type = c("raw", "dynamic", "overlap", "curve", "contrast"),
                      connected = FALSE,
                      ci = TRUE,
                      shade_periods = x$tr_period,
                      alpha_shade = 0.2,
                      palette = "Set2",
                      group_labels = c("Group 0", "Group 1"),
                      xlab = NULL,
                      ylab = NULL,
                      main = NULL,
                      ylim = NULL,
                      xlim = NULL,
                      curve = c("auto", "level", "derivative", "both"),
                      Xdistr = c("histogram", "hist", "density", "rug", "none"),
                      interval = c("auto", "pointwise", "uniform", "both", "none"),
                      show.uniform.CI = TRUE,
                      support.panel = c("embedded", "separate"),
                      show.eval_g = FALSE,
                      line.color = "steelblue",
                      line.size = 2,
                      ci.color = "steelblue",
                      ci.alpha = 0.25,
                      band.color = "gray45",
                      band.lty = 3,
                      hist.color = "gray50",
                      density.color = "gray50",
                      support.alpha = 0.25,
                      diff.values = NULL,
                      contrast.values = NULL,
                      ref.g = NULL,
                      target.g = NULL,
                      n.contrast = 7L,
                      max.contrast = 10L,
                      include.ref = TRUE,
                      target.type = c("level", "quantile"),
                      quantile.source = c("support", "eval_g"),
                      inference = c("auto", "replicate", "vcov", "band",
                                    "pointwise"),
                      interpolate = c("linear", "none"),
                      show.grid = FALSE,
                      legend.position = c("auto", "topright", "topleft",
                                          "bottomright", "bottomleft", "none"),
                      theme.bw = FALSE,
                      ...) {

  type <- match.arg(type)
  curve <- match.arg(curve)
  Xdistr <- match.arg(Xdistr)
  interval <- match.arg(interval)
  support.panel <- match.arg(support.panel)
  inference <- match.arg(inference)
  interpolate <- match.arg(interpolate)
  legend.position <- match.arg(legend.position)

  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Please install the 'RColorBrewer' package.")
  }

  available_palettes <- rownames(RColorBrewer::brewer.pal.info)
  if (!(palette %in% available_palettes)) {
    stop("Palette not found. Choose from: ",
         paste(available_palettes, collapse = ", "))
  }

  group_colors <- RColorBrewer::brewer.pal(3, palette)[1:2]

  if (!is.null(xlim)) {
    if (!is.numeric(xlim) || length(xlim) != 2L || any(!is.finite(xlim)) ||
        xlim[1] >= xlim[2]) {
      stop("xlim must be a numeric vector of length 2 with xlim[1] < xlim[2].")
    }
  }
  if (!is.null(ylim)) {
    if (!is.numeric(ylim) || length(ylim) != 2L || any(!is.finite(ylim)) ||
        ylim[1] >= ylim[2]) {
      stop("ylim must be a numeric vector of length 2 with ylim[1] < ylim[2].")
    }
  }

 shade_treatment <- function(tr_period, yrange, col = "gray") {

    # remove shading: user passes NULL/empty
    if (is.null(tr_period) || length(tr_period) == 0 ||
        !is.finite(alpha_shade) || alpha_shade <= 0) {
      return(invisible(NULL))
    }

    draw_one <- function(tp) {
      tp <- as.numeric(tp)
      tp <- tp[!is.na(tp)]
      if (length(tp) == 0) return(invisible(NULL))

      rect(min(tp) - 0.5,
           yrange[1] - 10,
           max(tp) + 0.5,
           yrange[2] + 10,
           col = adjustcolor(col, alpha.f = alpha_shade),
           border = NA)
      invisible(NULL)
    }

    if (is.list(tr_period)) {
      for (tp in tr_period) draw_one(tp)
    } else {
      draw_one(tr_period)
    }

    invisible(NULL)
  }

  ## ------------------------------------------------------------------
  ## RAW means
  ## ------------------------------------------------------------------
  if (type == "raw") {
    rawdf <- x$raw_means
    if (is.null(rawdf)) stop("No 'raw_means' found in the 'fdid' object.")
    if (!all(c("time", "group", "meanY") %in% names(rawdf))) {
      stop("'raw_means' is missing required columns.")
    }

    times <- sort(unique(rawdf$time))

    # If CI columns exist (or seY exists), prepare 95% CI for plotting
    has_ci <- ci && all(c("CI_Lower", "CI_Upper") %in% names(rawdf))

    if (is.null(ylim)) {
      if (has_ci) {
        tmp <- range(rawdf$CI_Lower, rawdf$CI_Upper, na.rm = TRUE)
      } else {
        tmp <- range(rawdf$meanY, na.rm = TRUE)
      }
      ylim <- c(min(tmp[1], 0.99 * tmp[1]), max(tmp[2], 1.01 * tmp[2]))
    }

    raw_legend_below <- !identical(legend.position, "none") &&
      identical(legend.position, "auto")
    if (raw_legend_below) {
      old_mar <- graphics::par("mar")
      graphics::par(mar = .fdid_mar_with_bottom_legend(old_mar, group_labels))
      on.exit(graphics::par(mar = old_mar), add = TRUE)
    }

    plot(NULL, xlim = xlim %||% range(times), ylim = ylim,
         xlab = xlab %||% "Time",
         ylab = ylab %||% "Mean outcome",
         main = main, mgp = c(2, .8, 0), ...)
    shade_treatment(shade_periods, ylim)

    g1 <- rawdf[rawdf$group == "Group 1", ]
    g0 <- rawdf[rawdf$group == "Group 0", ]

    if (connected) {
      lines(g1$time, g1$meanY, col = group_colors[2], lwd = 2)
      lines(g0$time, g0$meanY, col = group_colors[1], lwd = 2)
    }

    # 95% CI error bars (default on when available)
    if (exists("has_ci") && isTRUE(has_ci)) {
      ok1 <- !is.na(g1$CI_Lower) & !is.na(g1$CI_Upper)
      ok0 <- !is.na(g0$CI_Lower) & !is.na(g0$CI_Upper)
      if (any(ok1)) {
        arrows(g1$time[ok1], g1$CI_Lower[ok1], g1$time[ok1], g1$CI_Upper[ok1],
               angle = 90, code = 3, length = 0.05, col = group_colors[2])
      }
      if (any(ok0)) {
        arrows(g0$time[ok0], g0$CI_Lower[ok0], g0$time[ok0], g0$CI_Upper[ok0],
               angle = 90, code = 3, length = 0.05, col = group_colors[1])
      }
    }

    points(g1$time, g1$meanY, pch = 16, col = group_colors[2])
    points(g0$time, g0$meanY, pch = 16, col = group_colors[1])

    if (!identical(legend.position, "none")) {
      if (raw_legend_below) {
        .fdid_draw_bottom_legend(
          legend = group_labels,
          col = group_colors,
          pch = 16
        )
      } else {
        legend(legend.position, legend = group_labels,
               col = group_colors, pch = 16, bty = "n")
      }
    }
    return(invisible(NULL))
  }

  ## ------------------------------------------------------------------
  ## DYNAMIC effects
  ## ------------------------------------------------------------------
  if (type == "dynamic") {
    dyn <- x$dynamic
    if (is.null(dyn) || nrow(dyn) == 0) stop("No dynamic data found.")
    yrs <- as.numeric(rownames(dyn))
    est <- dyn$Estimate
    lo  <- dyn$CI_Lower
    hi  <- dyn$CI_Upper

    if (is.null(ylim)) {
      if (isTRUE(ci)) {
        rng  <- range(lo, hi, na.rm = TRUE)
      } else {
        rng  <- range(est, na.rm = TRUE)
      }
      ylim <- c(min(rng[1], 0), max(rng[2], 0)) * 1.05
    }

    plot(yrs, est, type = "n", xlim = xlim %||% range(yrs), ylim = ylim,
         xlab = xlab %||% "Time",
         ylab = ylab %||% "Coefficients",
         main = main, mgp = c(2, .8, 0), ...)

    shade_treatment(shade_periods, ylim)
    abline(h = 0, col = "gray50", lwd = 2, lty = 2)

    if (connected) lines(yrs, est, lwd = 2)

    if (isTRUE(ci)) {
      ok <- lo != hi
      if (any(ok)) {
        arrows(yrs[ok], lo[ok], yrs[ok], hi[ok],
               angle = 90, code = 3, length = 0.05)
      }
    }
    points(yrs, est, pch = 16)
    return(invisible(NULL))
  }

  ## ------------------------------------------------------------------
  ## CONTRAST (level contrast against a reference G)
  ## ------------------------------------------------------------------
  if (type == "contrast") {
    contrast_df <- .fdid_contrast_plot_data(
      x,
      ref.g = ref.g,
      contrast.values = contrast.values,
      target.g = target.g,
      n.contrast = n.contrast,
      max.contrast = max.contrast,
      include.ref = include.ref,
      target.type = target.type,
      quantile.source = quantile.source,
      inference = inference,
      interpolate = interpolate
    )
    .fdid_plot_contrast_base(
      contrast_df = contrast_df,
      xlab = xlab,
      ylab = ylab,
      main = main,
      ylim = ylim,
      line.color = line.color,
      ci.color = ci.color,
      ci = ci,
      show.grid = show.grid,
      ...
    )
    return(invisible(NULL))
  }

  ## ------------------------------------------------------------------
  ## CURVE (level theta/mu or derivative delta)
  ## ------------------------------------------------------------------
  if (type == "curve") {
    if (is.null(x$eval_g) || is.null(x$curve_event)) {
      stop("No curve data in this fdid object. ",
           "Use type='dynamic' for scalar estimates, or use method='kernel' ",
           "or 'dml_flex'.")
    }

    if (!(x$method %in% c("kernel", "dml_flex"))) {
      stop("curve plots are available for method='kernel' and method='dml_flex' only.")
    }

    curves_to_plot <- if (identical(curve, "auto")) {
      if (identical(x$method, "dml_flex")) "level" else "derivative"
    } else if (identical(curve, "both")) {
      c("level", "derivative")
    } else {
      curve
    }

    if (length(curves_to_plot) > 1L) {
      oldpar <- graphics::par(no.readonly = TRUE)
      on.exit({
        graphics::layout(1)
        graphics::par(oldpar)
      }, add = TRUE)
      graphics::layout(matrix(seq_along(curves_to_plot), nrow = 1L))
      graphics::par(mar = c(4.0, 4.0, 2.8, 0.8), oma = c(0, 0, 0, 0))
    }

    for (cc in curves_to_plot) {
      curve_df <- .fdid_curve_plot_data(
        x, curve = cc, interval = interval, xlim = xlim,
        ci = ci, show.uniform.CI = show.uniform.CI
      )
      panel_main <- main
      if (length(curves_to_plot) > 1L && is.null(main)) {
        panel_main <- curve_df$curve_label[1L]
      }
      panel_support <- if (length(curves_to_plot) > 1L) "embedded" else support.panel
      .fdid_plot_curve_base(
        x = x,
        curve_df = curve_df,
        interval = interval,
        ci = ci,
        show.uniform.CI = show.uniform.CI,
        Xdistr = Xdistr,
        support.panel = panel_support,
        show.eval_g = show.eval_g,
        diff.values = diff.values,
        contrast.values = contrast.values,
        xlab = xlab,
        ylab = ylab,
        main = panel_main,
        ylim = ylim,
        xlim = xlim,
        line.color = line.color,
        line.size = line.size,
        ci.color = ci.color,
        ci.alpha = ci.alpha,
        band.color = band.color,
        band.lty = band.lty,
        hist.color = hist.color,
        density.color = density.color,
        support.alpha = support.alpha,
        show.grid = show.grid,
        legend.position = legend.position,
        ...
      )
    }
    return(invisible(NULL))
  }

  ## ------------------------------------------------------------------
  ## OVERLAP
  ## ------------------------------------------------------------------
  ps <- x$ps; G <- x$G
  if (is.null(ps) || is.null(G)) stop("No propensity scores found.")

  brks <- seq(0, 1, 0.02)
  h0 <- hist(ps[G == 0], breaks = brks, plot = FALSE)
  h1 <- hist(ps[G == 1], breaks = brks, plot = FALSE)
  h0$density <- h0$counts / sum(h0$counts)
  h1$density <- h1$counts / sum(h1$counts)

  maxd <- max(h0$density, h1$density)
  ylim2 <- ylim %||% c(-maxd, maxd) * 1.1

  overlap_legend_below <- !identical(legend.position, "none") &&
    identical(legend.position, "auto")
  if (overlap_legend_below) {
    old_mar <- graphics::par("mar")
    graphics::par(mar = .fdid_mar_with_bottom_legend(old_mar, group_labels))
    on.exit(graphics::par(mar = old_mar), add = TRUE)
  }

  plot(h0, freq = FALSE, col = adjustcolor(group_colors[1], .6), ylim = ylim2,
       xlab = xlab %||% "Propensity Score",
       main = main %||% "Overlap of Propensity Scores")
  h1$density <- -h1$density
  plot(h1, freq = FALSE, add = TRUE,
       col = adjustcolor(group_colors[2], .6))
  abline(h = 0, lty = 2)
  if (!identical(legend.position, "none")) {
    if (overlap_legend_below) {
      .fdid_draw_bottom_legend(
        legend = group_labels,
        fill = adjustcolor(group_colors, .6)
      )
    } else {
      legend(legend.position, legend = group_labels,
             fill = adjustcolor(group_colors, .6), bty = "n")
    }
  }

  invisible(NULL)
}
