utils::globalVariables("tr_period")

#' Summary Method for FDID Objects
#'
#' @param object An object of class \code{fdid}.
#' @param ... Additional arguments (not used).
#'
#' @return Prints a summary of the \code{fdid} object.
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
#' summary(result)
#' }
#' @author Rivka Lipkovitz, Enhan Liu
#' @export
summary.fdid <- function(object, ...) {
  if (!inherits(object, "fdid")) {
    stop("Object must be of class 'fdid'.")
  }


  # --- Helper: format one row ---
  fmt_row <- function(label, est, se, ci_lo, ci_hi, lw = 16, nw = 11) {
    sprintf("  %-*s %*s  %*s  [%s, %s]",
            lw, label,
            nw, formatC(est, format = "f", digits = 4),
            nw, formatC(se, format = "f", digits = 4),
            formatC(ci_lo, format = "f", digits = 4),
            formatC(ci_hi, format = "f", digits = 4))
  }

  fmt_header <- function(lw = 16, nw = 11) {
    sprintf("  %-*s %*s  %*s  %s", lw, "", nw, "Estimate", nw, "Std.Error", "   95% CI")
  }

  # Keep summary output ASCII-safe for non-UTF-8 render environments.
  sep        <- paste(rep("-", 72), collapse = "")
  sep_double <- paste(rep("=", 72), collapse = "")
  dots       <- paste(rep(".", 68), collapse = "")

  # --- Period labels ---
  pre_times  <- object$entire_period[object$entire_period < min(object$tr_period)]
  post_times <- object$entire_period[object$entire_period > max(object$tr_period)]

  # === Header ===
  cat("\n")
  cat("  Factorial Difference-in-Differences (FDID) Summary\n")
  cat(" ", sep_double, "\n")
  cat("  Method:            ", object$method, "\n")
  cat("  Variance Type:     ", object$vartype, "\n")
  cat("  Reference Period:  ", paste(object$ref_period, collapse = ", "), "\n")
  if (length(pre_times) > 0)
    cat("  Pre-Event Period:  ", paste(pre_times, collapse = ", "), "\n")
  if (length(object$tr_period) > 0)
    cat("  Event Period:      ", paste(object$tr_period, collapse = ", "), "\n")
  if (length(post_times) > 0)
    cat("  Post-Event Period: ", paste(post_times, collapse = ", "), "\n")
  if (!is.null(object$cluster)) {
    n_clusters <- length(unique(object$cluster[!is.na(object$cluster)]))
    cat("  Clusters:          ", n_clusters, "\n")
  }
  if (!is.null(object$dml_metadata)) {
    if (!is.null(object$dml_metadata$target_estimand))
      cat("  DML Target:        ", object$dml_metadata$target_estimand, "\n")
    if (!is.null(object$dml_metadata$inference_method))
      cat("  DML Inference:     ", object$dml_metadata$inference_method, "\n")
  }
  cat(" ", sep_double, "\n\n")

  # === Aggregate Estimates ===
  cat("  Aggregate Estimates\n")
  cat(" ", sep, "\n")
  cat(fmt_header(), "\n")
  cat(" ", sep, "\n")

  # Helper: check if an estimate row is valid (not all NA)
  is_valid <- function(df) {
    !is.null(df) && !all(is.na(df$Estimate))
  }

  pre  <- object$est$pre
  evt  <- object$est$event
  post <- object$est$post

  if (is_valid(pre))
    cat(fmt_row("Pre-Event",  pre$Estimate,  pre$Std.Error,  pre$CI_Lower,  pre$CI_Upper),  "\n")
  if (is_valid(evt) && any(object$tr_period %in% object$entire_period))
    cat(fmt_row("Event",      evt$Estimate,  evt$Std.Error,  evt$CI_Lower,  evt$CI_Upper),  "\n")
  if (is_valid(post))
    cat(fmt_row("Post-Event", post$Estimate, post$Std.Error, post$CI_Lower, post$CI_Upper), "\n")
  cat(" ", sep, "\n\n")

  # === Dynamic Estimates ===
  cat("  Dynamic Estimates\n")
  cat(" ", sep, "\n")
  cat(fmt_header(), "\n")
  cat(" ", sep, "\n")

  dyn <- object$dynamic
  earliest_event <- min(object$tr_period)
  latest_event   <- max(object$tr_period)

  for (i in seq_len(nrow(dyn))) {
    yr   <- rownames(dyn)[i]
    yr_n <- as.numeric(yr)

    if (yr_n == earliest_event && i > 1) cat("  ", dots, "\n")
    if (yr_n == latest_event + 1)        cat("  ", dots, "\n")

    cat(fmt_row(yr, dyn$Estimate[i], dyn$Std.Error[i],
                dyn$CI_Lower[i], dyn$CI_Upper[i]), "\n")
  }

  cat(" ", sep, "\n")

  # Continuous-G curve summary (kernel / dml_flex)
  if (!is.null(object$eval_g) && length(object$eval_g) > 0) {
    cat("\n  Moderation Curve Summary (G-domain)\n")
    cat(" ", sep, "\n")
    cat(sprintf("  Method: %s  |  Eval points: %d  |  G range: [%.3f, %.3f]\n",
                object$method, length(object$eval_g),
                min(object$eval_g), max(object$eval_g)))
    if (identical(object$method, "dml_flex"))
      cat("  Target: level curve theta(g) and derivative curve theta'(g)\n")
    if (identical(object$method, "kernel"))
      cat("  Target: level curve mu(g) and derivative curve delta(g)\n")
	    if (!is.null(object$h0))
	      cat(sprintf("  Bandwidth h0: %.4f\n", object$h0))
	    if (identical(object$method, "kernel") &&
	        !is.null(object$bootstrap_resampling) &&
	        identical(object$vartype, "bootstrap")) {
	      cat(sprintf("  Kernel bootstrap: %s resampling", object$bootstrap_resampling))
	      if (!is.null(object$n_clusters) && is.finite(object$n_clusters))
	        cat(sprintf("  |  clusters: %d", object$n_clusters))
	      cat("\n")
	    }
    if (!is.null(object$learner))
      cat(sprintf("  Learner: %s  |  K folds: %d  |  S splits: %d\n",
                  object$learner, object$K, object$S))
    if (!is.null(object$signal_map))
      cat(sprintf("  Signal map: %s", object$signal_map))
    if (!is.null(object$density_method))
      cat(sprintf("  |  Density: %s", object$density_method))
    if (!is.null(object$signal_map) || !is.null(object$density_method))
      cat("\n")
	    if (!is.null(object$curve_event$simultaneous_band_method)) {
	      cat(sprintf("  Simultaneous bands: %s\n",
	                  object$curve_event$simultaneous_band_method))
	    }
	    if (identical(object$method, "kernel") &&
	        identical(object$vartype, "bootstrap") &&
	        !is.null(object$curve_event$delta_uniform_zeta)) {
	      cat(sprintf("  Kernel delta band zeta: %.4f",
	                  object$curve_event$delta_uniform_zeta))
	      if (!is.null(object$curve_event$delta_uniform_coverage) &&
	          is.finite(object$curve_event$delta_uniform_coverage)) {
	        cat(sprintf("  |  empirical coverage: %.3f",
	                    object$curve_event$delta_uniform_coverage))
	      }
	      cat("\n")
	    }
    if (identical(object$method, "dml_flex")) {
      cat("  (Full theta(g) and theta'(g) curves in $curve_event; scalar $est$event is grid-average theta)\n")
    } else {
      cat("  (Full mu(g) and delta(g) curves in $curve_event; scalar $est$event is interval average)\n")
    }
    cat(" ", sep, "\n")
  }

  if (identical(object$method, "dml_incremental")) {
    cat("\n  DML Incremental Signal\n")
    cat(" ", sep, "\n")
    cat("  Target: observed-population average derivative E[partial_g mu(G, X)]\n")
    if (!is.null(object$learner))
      cat(sprintf("  Learner: %s  |  K folds: %d  |  S splits: %d\n",
                  object$learner, object$K, object$S))
    if (!is.null(object$density_method))
      cat(sprintf("  Density: %s\n", object$density_method))
    cat("  (Scalar estimate uses the cross-fitted orthogonal average-derivative signal)\n")
    cat(" ", sep, "\n")
  }

  if (!is.null(object$dml_diagnostics)) {
    diag <- object$dml_diagnostics
    fmt_q <- function(x) {
      if (is.null(x) || all(is.na(x))) return("not available")
      keep <- intersect(c("p05", "p50", "p95"), names(x))
      paste(sprintf("%s=%.4f", keep, x[keep]), collapse = ", ")
    }
    cat("\n  DML Diagnostics\n")
    cat(" ", sep, "\n")
    if (!is.null(diag$learner))
      cat(sprintf("  Learner: %s\n", diag$learner))
    if (!is.null(diag$density_method) && !is.na(diag$density_method))
      cat(sprintf("  Density method: %s\n", diag$density_method))
    if (!is.null(diag$n) && !is.na(diag$n))
      cat(sprintf("  Complete observations: %d\n", diag$n))
    if (!is.null(diag$n_clusters) && !is.na(diag$n_clusters))
      cat(sprintf("  Clusters used for scalar DML SEs: %d\n", diag$n_clusters))
    if (!is.null(diag$min_density) && is.finite(diag$min_density)) {
      cat(sprintf("  Min s_hat(G|X): %.5f  |  share < 0.01: %.3f\n",
                  diag$min_density, diag$share_density_below_1e_2))
    }
    cat("  Signal quantiles: ", fmt_q(diag$signal_quantiles), "\n", sep = "")
    if (!is.null(diag$d_mu_quantiles) && !all(is.na(diag$d_mu_quantiles)))
      cat("  partial_g mu quantiles: ", fmt_q(diag$d_mu_quantiles), "\n", sep = "")
    if (!is.null(diag$d_log_density_quantiles) && !all(is.na(diag$d_log_density_quantiles)))
      cat("  partial_g log s quantiles: ", fmt_q(diag$d_log_density_quantiles), "\n", sep = "")
    if (!is.null(diag$correction_quantiles) && !all(is.na(diag$correction_quantiles)))
      cat("  Residual-correction quantiles: ", fmt_q(diag$correction_quantiles), "\n", sep = "")
    cat(" ", sep, "\n")
  }

  invisible(object)
}
