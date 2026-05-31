#' Expand Prepared FDID Covariates
#'
#' Replaces the prepared \code{x1}, \code{x2}, ... covariates with a richer
#' design matrix formed from raw covariates, polynomial or B-spline basis
#' functions, and optional pairwise interactions.
#'
#' @param s A data frame produced by \code{fdid_prepare}.
#' @param basis_type Basis expansion for covariates. Options are
#'   \code{"none"} (raw covariates), \code{"polynomial"}, and \code{"bspline"}.
#' @param include_interactions Logical; if \code{TRUE}, add all pairwise
#'   interactions among the expanded main-effect basis columns.
#' @param poly_degree Integer polynomial degree used when
#'   \code{basis_type = "polynomial"}.
#' @param spline_df Integer degrees of freedom used when
#'   \code{basis_type = "bspline"}.
#' @param spline_degree Integer spline degree used when
#'   \code{basis_type = "bspline"}.
#' @param standardize_basis Logical; if \code{TRUE}, center and scale expanded
#'   covariate columns. This is useful for regularized and ML learners.
#'
#' @return A prepared FDID data frame with expanded covariates renamed to
#'   \code{x1}, \code{x2}, ...
#'
#' @export
fdid_expand_covariates <- function(s,
                                   basis_type = c("none", "polynomial", "bspline"),
                                   include_interactions = FALSE,
                                   poly_degree = 2L,
                                   spline_df = 4L,
                                   spline_degree = 2L,
                                   standardize_basis = TRUE) {
  stopifnot(is.data.frame(s))
  basis_type <- match.arg(basis_type)
  include_interactions <- isTRUE(include_interactions)

  x_cols <- grep("^x[0-9]+$", names(s), value = TRUE)
  if (length(x_cols) == 0L) return(s)
  x_cols <- x_cols[order(as.integer(sub("^x", "", x_cols)))]

  if (basis_type == "none" && !include_interactions) return(s)

  X <- as.matrix(s[, x_cols, drop = FALSE])
  storage.mode(X) <- "double"

  make_poly <- function(z, nm) {
    deg <- max(1L, as.integer(poly_degree))
    out <- sapply(seq_len(deg), function(d) z^d)
    colnames(out) <- paste0(nm, "_p", seq_len(deg))
    out
  }

  make_bs <- function(z, nm) {
    finite_unique <- unique(z[is.finite(z)])
    if (length(finite_unique) <= spline_degree + 1L ||
        stats::sd(z, na.rm = TRUE) < 1e-8) {
      out <- matrix(z, ncol = 1L)
      colnames(out) <- nm
      return(out)
    }
    out <- tryCatch(
      {
        bs <- splines::bs(
          z,
          df = max(1L, as.integer(spline_df)),
          degree = max(1L, as.integer(spline_degree)),
          intercept = FALSE
        )
        bs <- as.matrix(bs)
        colnames(bs) <- paste0(nm, "_bs", seq_len(ncol(bs)))
        bs
      },
      error = function(e) {
        raw <- matrix(z, ncol = 1L)
        colnames(raw) <- nm
        raw
      }
    )
    out
  }

  main_list <- lapply(seq_along(x_cols), function(j) {
    z <- X[, j]
    if (basis_type == "polynomial") {
      make_poly(z, x_cols[j])
    } else if (basis_type == "bspline") {
      make_bs(z, x_cols[j])
    } else {
      out <- matrix(z, ncol = 1L)
      colnames(out) <- x_cols[j]
      out
    }
  })
  X_exp <- do.call(cbind, main_list)

  if (include_interactions && ncol(X_exp) > 1L) {
    pair_idx <- utils::combn(seq_len(ncol(X_exp)), 2L)
    X_int <- matrix(NA_real_, nrow(X_exp), ncol(pair_idx))
    int_names <- character(ncol(pair_idx))
    for (j in seq_len(ncol(pair_idx))) {
      a <- pair_idx[1L, j]
      b <- pair_idx[2L, j]
      X_int[, j] <- X_exp[, a] * X_exp[, b]
      int_names[j] <- paste0(colnames(X_exp)[a], "_x_", colnames(X_exp)[b])
    }
    colnames(X_int) <- int_names
    X_exp <- cbind(X_exp, X_int)
  }

  keep <- vapply(seq_len(ncol(X_exp)), function(j) {
    z <- X_exp[, j]
    any(is.finite(z)) && is.finite(stats::sd(z, na.rm = TRUE)) &&
      stats::sd(z, na.rm = TRUE) > 1e-10
  }, logical(1L))
  X_exp <- X_exp[, keep, drop = FALSE]

  if (ncol(X_exp) > 0L && isTRUE(standardize_basis)) {
    centers <- colMeans(X_exp, na.rm = TRUE)
    scales <- apply(X_exp, 2L, stats::sd, na.rm = TRUE)
    scales[!is.finite(scales) | scales < 1e-8] <- 1
    X_exp <- sweep(sweep(X_exp, 2L, centers, "-"), 2L, scales, "/")
  }

  non_x <- s[, setdiff(names(s), x_cols), drop = FALSE]
  if (ncol(X_exp) == 0L) {
    out <- non_x
  } else {
    X_df <- as.data.frame(X_exp)
    names(X_df) <- paste0("x", seq_len(ncol(X_df)))
    out <- cbind(non_x, X_df)
  }
  attr(out, "fdid_basis") <- list(
    basis_type = basis_type,
    include_interactions = include_interactions,
    poly_degree = as.integer(poly_degree),
    spline_df = as.integer(spline_df),
    spline_degree = as.integer(spline_degree),
    standardize_basis = isTRUE(standardize_basis),
    raw_covariates = x_cols,
    expanded_terms = colnames(X_exp)
  )
  out
}

#' Prepare Data for Factorial Difference-in-Differences Analysis
#'
#' Prepares a dataset for factorial difference-in-differences (FDID) analysis by reshaping the data into a wide format,
#' averaging time-varying covariates, and renaming columns for consistency in subsequent analysis.
#'
#' @param data A data frame containing the dataset to be processed.
#' @param Y_label A string specifying the column name of the outcome variable.
#' @param X_labels A character vector specifying the column names of the time-varying covariates.
#' @param G_label A string specifying the column name of the group variable (e.g., treatment vs. control).
#' @param unit_label A string specifying the column name of the unit identifier (e.g., individual or entity).
#' @param time_label A string specifying the column name of the time variable.
#' @param cluster_label An optional string specifying the column name of the clustering variable. Default is `NULL`.
#' @param basis_type Optional basis expansion for covariates before estimation.
#'   Options are \code{"none"} (default), \code{"polynomial"}, and
#'   \code{"bspline"}.
#' @param include_interactions Logical; if \code{TRUE}, add pairwise
#'   interactions among the expanded covariate basis columns.
#' @param poly_degree Integer polynomial degree used when
#'   \code{basis_type = "polynomial"}.
#' @param spline_df Integer degrees of freedom used when
#'   \code{basis_type = "bspline"}.
#' @param spline_degree Integer spline degree used when
#'   \code{basis_type = "bspline"}.
#' @param standardize_basis Logical; if \code{TRUE}, center and scale expanded
#'   covariate columns. Ignored when no expansion or interactions are requested.
#'
#' @return A data frame in wide format with the following:
#' - Outcome variable pivoted to wide format with time columns.
#' - Time-varying covariates averaged across time.
#' - Columns renamed:
#'   - Unit identifier -> `unit`
#'   - Raw or expanded covariates -> `x1`, `x2`, ...
#'   - Group variable -> `G`
#'   - Clustering variable (if provided) -> `c`
#'
#' @examples
#' data <- data.frame(
#'   id = rep(1:3, each = 4),
#'   time = rep(1:4, times = 3),
#'   outcome = rnorm(12),
#'   covar1 = runif(12),
#'   covar2 = runif(12),
#'   group = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1)
#' )
#' fdid_data <- fdid_prepare(
#'   data = data,
#'   Y_label = "outcome",
#'   X_labels = c("covar1", "covar2"),
#'   G_label = "group",
#'   unit_label = "id",
#'   time_label = "time"
#' )
#' head(fdid_data)
#'
#' @import dplyr
#' @import tidyr
#' @author Rivka Lipkovitz
#' @export
fdid_prepare <- function(data,
                             Y_label,
                             X_labels = NULL,
                             G_label,
                             unit_label,
                             time_label,
                             cluster_label = NULL,
                             basis_type = c("none", "polynomial", "bspline"),
                             include_interactions = FALSE,
                             poly_degree = 2L,
                             spline_df = 4L,
                             spline_degree = 2L,
                             standardize_basis = TRUE) {

  basis_type <- match.arg(basis_type)

  # 1. Pivot the outcome to wide
  wide_data <- data %>%
    tidyr::pivot_wider(
      id_cols     = tidyselect::all_of(unit_label),
      names_from  = tidyselect::all_of(time_label),
      values_from = tidyselect::all_of(Y_label),
      names_prefix = "Y_"  # keep the time values as column names
    )

  # 2. Average the time-varying covariates using older dplyr syntax
  #    Summarize only the X_labels
  if (!is.null(X_labels) && length(X_labels) > 0L) {
    covar_data <- data %>%
      group_by(!!rlang::sym(unit_label)) %>%
      dplyr::summarise_at(
        .vars = X_labels,
        .funs = ~ mean(.x, na.rm = TRUE)
      ) %>%
      ungroup()
  } else {
    covar_data <- data %>%
      distinct(!!rlang::sym(unit_label))
  }

  # 3. Extract G_label and cluster_label (if present).
  #    We'll take the first occurrence for each unit
  #    (assuming they do not vary over time).
  #    Then we'll join these back to covar_data.
  if (!is.null(cluster_label)) {
    group_cluster_data <- data %>%
      distinct(
        !!rlang::sym(unit_label),
        !!rlang::sym(G_label),
        !!rlang::sym(cluster_label)
      )
    covar_data <- covar_data %>%
      left_join(group_cluster_data, by = unit_label)
  } else {
    group_data <- data %>%
      distinct(
        !!rlang::sym(unit_label),
        !!rlang::sym(G_label)
      )
    covar_data <- covar_data %>%
      left_join(group_data, by = unit_label)
  }

  # 4. Join these averaged covariates + group info back to wide outcome data
  wide_data <- wide_data %>%
    left_join(covar_data, by = unit_label)

  # 5. Rename columns:
  #    - First column => "unit"
  #    - Covariates => x1, x2, ...
  #    - G_label => "G"
  #    - cluster_label => "c"
  #    (using older rename_at / rename_with might need caution in older dplyr)
  colnames(wide_data)[1] <- "unit"

  # figure out which columns are the newly-averaged covariates
  # (the same names as X_labels), so we can rename them to x1, x2, ...
  x_cols <- X_labels
  # we also rename them in the same order they appear in X_labels
  if (!is.null(x_cols) && length(x_cols) > 0L) {
    for (i in seq_along(x_cols)) {
      old_col <- x_cols[i]
      new_col <- paste0("x", i)
      wide_data <- wide_data %>%
        rename(!!new_col := !!rlang::sym(old_col))
    }
  }

  # rename G_label => "G"
  wide_data <- wide_data %>%
    rename(G = !!rlang::sym(G_label))

  # rename cluster_label => "c"
  if (!is.null(cluster_label)) {
    wide_data <- wide_data %>%
      rename(c = !!rlang::sym(cluster_label))
  }

  if ((!identical(basis_type, "none") || isTRUE(include_interactions)) &&
      !is.null(X_labels) && length(X_labels) > 0L) {
    wide_data <- fdid_expand_covariates(
      wide_data,
      basis_type = basis_type,
      include_interactions = include_interactions,
      poly_degree = poly_degree,
      spline_df = spline_df,
      spline_degree = spline_degree,
      standardize_basis = standardize_basis
    )
  }

  return(wide_data)
}
