#' Report a continuous-G curve contrast
#'
#' Extract a level-curve contrast, such as \eqn{\theta(g_1)-\theta(g_0)}
#' for \code{method = "dml_flex"} or \eqn{\mu(g_1)-\mu(g_0)} for
#' \code{method = "kernel"}. The helper uses stored curve values and, when
#' available, stored curve covariance matrices or replicate curves.
#'
#' @param x An object returned by \code{\link{fdid}} with \code{method} equal
#'   to \code{"kernel"} or \code{"dml_flex"}.
#' @param g0,g1 Numeric scalar endpoints for the contrast.
#' @param type Curve type. \code{"level"} is the usual Notes-style contrast;
#'   \code{"derivative"} reports a contrast between derivative values.
#' @param inference Inference rule. \code{"auto"} uses stored curve covariance
#'   matrices when available, then replicate curves, then conservative
#'   band-implied intervals, then approximate pointwise intervals.
#'   \code{"pointwise"} ignores endpoint covariance and is labeled approximate.
#'   \code{"band"} uses stored simultaneous or bootstrap bands when available.
#'   \code{"replicate"} uses stored bootstrap or multiplier replicate curves.
#'   \code{"vcov"} uses the covariance identity
#'   \eqn{V_{11} + V_{22} - 2V_{12}} when a curve covariance matrix is stored.
#' @param interpolate If \code{"linear"}, linearly interpolate between stored
#'   \code{eval_g} rows. If \code{"none"}, \code{g0} and \code{g1} must be on
#'   the fitted grid.
#' @param alpha Significance level. Defaults to the fitted object's
#'   \code{alpha}, or 0.05 for older objects.
#'
#' @return A one-row data frame with estimate, standard error when available,
#'   confidence limits, inference label, and diagnostic note.
#' @export
fdid_contrast <- function(x, g0, g1,
                          type = c("level", "derivative"),
                          inference = c("auto", "replicate", "vcov", "band", "pointwise"),
                          interpolate = c("linear", "none"),
                          alpha = NULL) {
  type <- match.arg(type)
  inference <- match.arg(inference)
  interpolate <- match.arg(interpolate)
  alpha <- .fdid_alpha(x, alpha)
  .fdid_check_g(g0, "g0")
  .fdid_check_g(g1, "g1")

  comp <- .fdid_curve_components(x, type)
  est0 <- .fdid_interp(comp$g, comp$estimate, g0, interpolate)
  est1 <- .fdid_interp(comp$g, comp$estimate, g1, interpolate)
  estimate <- est1 - est0

  candidates <- if (inference == "auto") {
    c("vcov", "replicate", "band", "pointwise")
  } else {
    inference
  }

  out <- NULL
  for (rule in candidates) {
    out <- switch(
      rule,
	      replicate = .fdid_contrast_replicate(comp, g0, g1, estimate, interpolate,
	                                           alpha),
	      vcov = .fdid_contrast_vcov(comp, g0, g1, estimate, interpolate, alpha),
	      band = .fdid_contrast_band(comp, g0, g1, estimate, interpolate),
      pointwise = .fdid_contrast_pointwise(comp, g0, g1, estimate, interpolate,
                                           alpha)
    )
    if (!is.null(out)) break
  }
  if (is.null(out)) {
    stop("No usable inference information is available for this contrast.")
  }

  .fdid_result_row(
    target = if (type == "level") "level contrast" else "derivative contrast",
    method = x$method,
    curve = type,
    g0 = g0,
    g1 = g1,
    estimate = estimate,
    std_error = out$std_error,
    conf_low = out$conf_low,
    conf_high = out$conf_high,
    inference = out$inference,
    note = out$note
  )
}

#' Report a continuous-G derivative value
#'
#' Extract \eqn{\theta'(g_0)} for \code{method = "dml_flex"} or
#' \eqn{\delta(g_0)} for \code{method = "kernel"} from a fitted continuous-G
#' curve object.
#'
#' @inheritParams fdid_contrast
#' @param g0 Numeric scalar evaluation point.
#'
#' @return A one-row data frame with estimate, standard error when available,
#'   confidence limits, inference label, and diagnostic note.
#' @export
fdid_derivative <- function(x, g0,
	                            inference = c("auto", "replicate", "vcov", "band", "pointwise"),
                            interpolate = c("linear", "none"),
                            alpha = NULL) {
  inference <- match.arg(inference)
  interpolate <- match.arg(interpolate)
  alpha <- .fdid_alpha(x, alpha)
  .fdid_check_g(g0, "g0")

  comp <- .fdid_curve_components(x, "derivative")
  estimate <- .fdid_interp(comp$g, comp$estimate, g0, interpolate)

  candidates <- if (inference == "auto") {
    c("vcov", "replicate", "band", "pointwise")
  } else {
    inference
  }

  out <- NULL
  for (rule in candidates) {
    out <- switch(
      rule,
	      replicate = .fdid_value_replicate(comp, g0, estimate, interpolate, alpha),
	      vcov = .fdid_value_vcov(comp, g0, estimate, interpolate, alpha),
	      band = .fdid_value_band(comp, g0, estimate, interpolate),
      pointwise = .fdid_value_pointwise(comp, g0, estimate, interpolate, alpha)
    )
    if (!is.null(out)) break
  }
  if (is.null(out)) {
    stop("No usable inference information is available for this derivative value.")
  }

  .fdid_result_row(
    target = "derivative value",
    method = x$method,
    curve = "derivative",
    g0 = g0,
    g1 = NA_real_,
    estimate = estimate,
    std_error = out$std_error,
    conf_low = out$conf_low,
    conf_high = out$conf_high,
    inference = out$inference,
    note = out$note
  )
}

.fdid_alpha <- function(x, alpha) {
  if (is.null(alpha)) alpha <- x$alpha
  if (is.null(alpha)) alpha <- 0.05
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single finite number strictly between 0 and 1.")
  }
  as.numeric(alpha)
}

.fdid_check_g <- function(g, nm) {
  if (!is.numeric(g) || length(g) != 1L || !is.finite(g)) {
    stop(nm, " must be a single finite number.")
  }
}

.fdid_curve_components <- function(x, type) {
  if (!inherits(x, "fdid")) {
    stop("x must be an object returned by fdid().")
  }
  if (!(x$method %in% c("kernel", "dml_flex"))) {
    stop("Continuous-G curve reporting is available only for method='kernel' ",
         "and method='dml_flex'.")
  }
  if (is.null(x$eval_g) || is.null(x$curve_event)) {
    stop("x does not contain stored continuous-G curve data.")
  }

  ce <- x$curve_event
  g <- as.numeric(x$eval_g)
  if (type == "level") {
	    if (identical(x$method, "kernel")) {
	      estimate <- ce$mu_hat
	      se <- ce$se_mu
		      lower <- ce$mu_band_lower
		      upper <- ce$mu_band_upper
		      reps <- ce$mu_bootstrap
		      vcov <- if (identical(x$vartype, "bootstrap") &&
		                  !is.null(ce$mu_boot_vcov)) ce$mu_boot_vcov else ce$mu_vcov
    } else {
      estimate <- ce$theta_hat
      se <- ce$se_theta
	      lower <- ce$theta_band_lower
	      upper <- ce$theta_band_upper
	      reps <- ce$theta_multiplier_reps
	      vcov <- ce$theta_vcov
    }
  } else {
    estimate <- ce$delta_hat
    se <- ce$se_delta
	    if (identical(x$method, "kernel")) {
		      lower <- ce$delta_band_lower
		      upper <- ce$delta_band_upper
		      reps <- ce$delta_bootstrap
		      vcov <- if (identical(x$vartype, "bootstrap") &&
		                  !is.null(ce$delta_boot_vcov)) ce$delta_boot_vcov else ce$delta_vcov
    } else {
	      lower <- ce$delta_band_lower
	      upper <- ce$delta_band_upper
	      reps <- ce$delta_multiplier_reps
	      vcov <- ce$delta_vcov
    }
  }

  if (is.null(estimate) || length(estimate) != length(g)) {
    stop("The requested curve is missing or does not match eval_g.")
  }
  list(
    g = g,
    estimate = as.numeric(estimate),
    se = if (is.null(se)) rep(NA_real_, length(g)) else as.numeric(se),
    lower = if (is.null(lower)) NULL else as.numeric(lower),
	    upper = if (is.null(upper)) NULL else as.numeric(upper),
	    reps = reps,
	    vcov = vcov
	  )
	}

.fdid_interp <- function(g, y, g0, interpolate) {
  ok <- is.finite(g) & is.finite(y)
  if (!any(ok)) stop("No finite curve values are available.")
  g <- g[ok]
  y <- y[ok]
  if (g0 < min(g) || g0 > max(g)) {
    stop("Requested G value is outside the fitted eval_g range.")
  }
  exact <- which(abs(g - g0) <= sqrt(.Machine$double.eps))
  if (length(exact) > 0L) return(y[exact[1L]])
  if (interpolate == "none") {
    stop("Requested G value is not on the fitted eval_g grid; use ",
         "interpolate='linear' or refit with this G value in eval_g.")
  }
  if (length(g) < 2L) {
    stop("Linear interpolation needs at least two finite curve grid points.")
  }
  stats::approx(g, y, xout = g0, ties = mean)$y
}

.fdid_interp_reps <- function(g, reps, g0, interpolate) {
  if (is.null(reps) || !is.matrix(reps) || nrow(reps) < 2L ||
      ncol(reps) != length(g)) {
    return(NULL)
  }
  vals <- apply(reps, 1L, function(row) {
    tryCatch(.fdid_interp(g, as.numeric(row), g0, interpolate),
             error = function(e) NA_real_)
  })
  vals[is.finite(vals)]
}

.fdid_contrast_replicate <- function(comp, g0, g1, estimate, interpolate, alpha) {
  r0 <- .fdid_interp_reps(comp$g, comp$reps, g0, interpolate)
  r1 <- .fdid_interp_reps(comp$g, comp$reps, g1, interpolate)
  if (is.null(r0) || is.null(r1)) return(NULL)
  n <- min(length(r0), length(r1))
  if (n < 2L) return(NULL)
  rr <- r1[seq_len(n)] - r0[seq_len(n)]
  rr <- rr[is.finite(rr)]
  if (length(rr) < 2L) return(NULL)
  ci <- stats::quantile(rr, probs = c(alpha / 2, 1 - alpha / 2),
                        names = FALSE, na.rm = TRUE)
  list(
    std_error = stats::sd(rr, na.rm = TRUE),
    conf_low = ci[1L],
    conf_high = ci[2L],
    inference = "replicate",
    note = paste0("Uses stored bootstrap or multiplier replicate curves; ",
                  "replicate count = ", length(rr), ".")
  )
}

.fdid_grid_weights <- function(g, g0, interpolate) {
  ok <- is.finite(g)
  if (!any(ok)) stop("No finite eval_g values are available.")
  idx <- which(ok)
  g_ok <- g[ok]
  if (g0 < min(g_ok) || g0 > max(g_ok)) {
    stop("Requested G value is outside the fitted eval_g range.")
  }
  exact <- which(abs(g_ok - g0) <= sqrt(.Machine$double.eps))
  w <- numeric(length(g))
  if (length(exact) > 0L) {
    w[idx[exact[1L]]] <- 1
    return(w)
  }
  if (interpolate == "none") {
    stop("Requested G value is not on the fitted eval_g grid; use ",
         "interpolate='linear' or refit with this G value in eval_g.")
  }
  if (length(g_ok) < 2L) {
    stop("Linear interpolation needs at least two finite curve grid points.")
  }
  ord <- order(g_ok)
  g_ord <- g_ok[ord]
  idx_ord <- idx[ord]
  j <- findInterval(g0, g_ord)
  if (j <= 0L) j <- 1L
  if (j >= length(g_ord)) j <- length(g_ord) - 1L
  g_left <- g_ord[j]
  g_right <- g_ord[j + 1L]
  if (!is.finite(g_right - g_left) || abs(g_right - g_left) < 1e-12) {
    stop("Linear interpolation is undefined for duplicate eval_g values.")
  }
  lambda <- (g0 - g_left) / (g_right - g_left)
  w[idx_ord[j]] <- 1 - lambda
  w[idx_ord[j + 1L]] <- lambda
  w
}

.fdid_contrast_vcov <- function(comp, g0, g1, estimate, interpolate, alpha) {
  if (is.null(comp$vcov) || !is.matrix(comp$vcov)) return(NULL)
  if (nrow(comp$vcov) != length(comp$g) || ncol(comp$vcov) != length(comp$g)) {
    return(NULL)
  }
  w0 <- .fdid_grid_weights(comp$g, g0, interpolate)
  w1 <- .fdid_grid_weights(comp$g, g1, interpolate)
  w <- w1 - w0
  v <- as.numeric(t(w) %*% comp$vcov %*% w)
  if (!is.finite(v)) return(NULL)
  se <- sqrt(pmax(0, v))
  z <- stats::qnorm(1 - alpha / 2)
  list(
    std_error = se,
    conf_low = estimate - z * se,
    conf_high = estimate + z * se,
    inference = "vcov",
    note = "Uses the stored curve covariance matrix; exact-grid contrasts reduce to V11 + V22 - 2*V12."
  )
}

.fdid_contrast_band <- function(comp, g0, g1, estimate, interpolate) {
  if (is.null(comp$lower) || is.null(comp$upper)) return(NULL)
  lo0 <- .fdid_interp(comp$g, comp$lower, g0, interpolate)
  hi0 <- .fdid_interp(comp$g, comp$upper, g0, interpolate)
  lo1 <- .fdid_interp(comp$g, comp$lower, g1, interpolate)
  hi1 <- .fdid_interp(comp$g, comp$upper, g1, interpolate)
  if (!all(is.finite(c(lo0, hi0, lo1, hi1)))) return(NULL)
  list(
    std_error = NA_real_,
    conf_low = lo1 - hi0,
    conf_high = hi1 - lo0,
    inference = "band",
    note = "Conservative band-implied interval: lower(g1)-upper(g0), upper(g1)-lower(g0)."
  )
}

.fdid_contrast_pointwise <- function(comp, g0, g1, estimate, interpolate, alpha) {
  se0 <- .fdid_interp(comp$g, comp$se, g0, interpolate)
  se1 <- .fdid_interp(comp$g, comp$se, g1, interpolate)
  if (!all(is.finite(c(se0, se1)))) return(NULL)
  se <- sqrt(se0^2 + se1^2)
  z <- stats::qnorm(1 - alpha / 2)
  list(
    std_error = se,
    conf_low = estimate - z * se,
    conf_high = estimate + z * se,
    inference = "pointwise_approx",
    note = "Approximate endpoint-independent interval; endpoint covariance is not used."
  )
}

.fdid_value_replicate <- function(comp, g0, estimate, interpolate, alpha) {
  rr <- .fdid_interp_reps(comp$g, comp$reps, g0, interpolate)
  if (is.null(rr) || length(rr) < 2L) return(NULL)
  ci <- stats::quantile(rr, probs = c(alpha / 2, 1 - alpha / 2),
                        names = FALSE, na.rm = TRUE)
  list(
    std_error = stats::sd(rr, na.rm = TRUE),
    conf_low = ci[1L],
    conf_high = ci[2L],
    inference = "replicate",
    note = paste0("Uses stored bootstrap or multiplier replicate curves; ",
                  "replicate count = ", length(rr), ".")
  )
}

.fdid_value_vcov <- function(comp, g0, estimate, interpolate, alpha) {
  if (is.null(comp$vcov) || !is.matrix(comp$vcov)) return(NULL)
  if (nrow(comp$vcov) != length(comp$g) || ncol(comp$vcov) != length(comp$g)) {
    return(NULL)
  }
  w <- .fdid_grid_weights(comp$g, g0, interpolate)
  v <- as.numeric(t(w) %*% comp$vcov %*% w)
  if (!is.finite(v)) return(NULL)
  se <- sqrt(pmax(0, v))
  z <- stats::qnorm(1 - alpha / 2)
  list(
    std_error = se,
    conf_low = estimate - z * se,
    conf_high = estimate + z * se,
    inference = "vcov",
    note = "Uses the stored curve covariance matrix at the requested G value."
  )
}

.fdid_value_band <- function(comp, g0, estimate, interpolate) {
  if (is.null(comp$lower) || is.null(comp$upper)) return(NULL)
  lo <- .fdid_interp(comp$g, comp$lower, g0, interpolate)
  hi <- .fdid_interp(comp$g, comp$upper, g0, interpolate)
  if (!all(is.finite(c(lo, hi)))) return(NULL)
  list(
    std_error = NA_real_,
    conf_low = lo,
    conf_high = hi,
    inference = "band",
    note = "Uses stored curve band at the requested G value."
  )
}

.fdid_value_pointwise <- function(comp, g0, estimate, interpolate, alpha) {
  se <- .fdid_interp(comp$g, comp$se, g0, interpolate)
  if (!is.finite(se)) return(NULL)
  z <- stats::qnorm(1 - alpha / 2)
  list(
    std_error = se,
    conf_low = estimate - z * se,
    conf_high = estimate + z * se,
    inference = "pointwise",
    note = "Pointwise normal interval at the requested G value."
  )
}

.fdid_result_row <- function(target, method, curve, g0, g1, estimate,
                             std_error, conf_low, conf_high, inference, note) {
  data.frame(
    target = target,
    method = method,
    curve = curve,
    g0 = g0,
    g1 = g1,
    estimate = estimate,
    std.error = std_error,
    conf.low = conf_low,
    conf.high = conf_high,
    inference = inference,
    note = note,
    stringsAsFactors = FALSE
  )
}
