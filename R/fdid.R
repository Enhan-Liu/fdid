#' Factorial Difference-in-Differences Estimation
#'
#' Performs factorial difference-in-differences (FDID) estimation using various methods and variance estimation techniques.
#'
#' @param s A data frame prepared using \code{fdid_prepare}.
#' @param tr_period A numeric vector specifying the treatment periods.
#' @param ref_period A numeric scalar specifying the reference period.
#' @param entire_period A numeric vector specifying the total range of time periods.
#'   If \code{NULL}, estimation is performed on all available time periods.
#'   Example: \code{c(1958, 1959, 1960, 1961)}.
#' @param method A string specifying the estimation method.
#'   Options: \code{"ols1"}, \code{"ols2"}, \code{"did"}, \code{"ebal"}, \code{"ipw"}, \code{"aipw"},
#'   \code{"kernel"}, \code{"dml_binary"}, \code{"dml_plr"}, \code{"dml_flex"},
#'   \code{"dml_incremental"}.
#'   Default is \code{"ols1"}.
#' @param vartype A string specifying the variance estimation type.
#'   Options: \code{"robust"}, \code{"bootstrap"}, \code{"jackknife"}.
#'   For \code{method="kernel"}, currently only \code{"robust"} and
#'   \code{"bootstrap"} are supported. DML methods currently support
#'   \code{"robust"} only; scalar DML robust standard errors are
#'   cluster-aware when \code{fdid_prepare(cluster_label = ...)} supplies a
#'   cluster column.
#'   Default is \code{"robust"}.
#' @param missing_data How to handle missing data. Two options:
#' \itemize{
#'   \item \code{"listwise"}: Drop any row missing \strong{any} relevant column (including
#'     outcomes in the periods used).
#'   \item \code{"available"}: Drop rows only if they are missing in group/covariates/cluster
#'     columns, but allow partial usage of outcomes.
#' }
#' Default is \code{"listwise"}.
#' @param nsims Number of simulations for bootstrap variance estimation.
#'   Default is \code{1000}.
#' @param parallel Logical; whether to perform parallel computations.
#'   Default is \code{FALSE}.
#' @param cores Number of cores for parallel computations.
#'   Default is \code{2}.
#' @param target.pop Character; the target population for averaging: \code{"all"}, \code{"1"}, or \code{"0"}.
#'   \code{"all"} corresponds to the full sample. \code{"1"} targets the \code{G=1} population.
#'   \code{"0"} targets the \code{G=0} population. Default is \code{"all"}.
#' @param eval_g Numeric vector of evaluation points for continuous-G curve methods
#'   (\code{"kernel"}, \code{"dml_flex"}).
#'   Default \code{NULL} auto-selects a trimmed grid.
#' @param h0 Numeric baseline bandwidth for \code{method="kernel"}.
#'   Default \code{NULL} selects by 10-fold LSCV.
#' @param K_folds Integer; number of CV folds for kernel bandwidth selection. Default 10.
#' @param boot Integer; bootstrap or multiplier replications for curve bands.
#'   For \code{method="kernel"} with \code{vartype="bootstrap"}, this controls
#'   the kernel bootstrap level and derivative bands and stored replicate
#'   curves. For \code{method="dml_flex"}, this controls the practical
#'   multiplier band and stored multiplier replicate curves over
#'   signal-mapping residuals, and Gaussian max-\eqn{t} draws for BLP-spline
#'   covariance bands when available.
#'   Default 500.
#' @param alpha Numeric significance level for confidence intervals. Default 0.05.
#' @param trim Numeric tail-trim fraction for \code{eval_g} construction. Default 0.05.
#' @param learner Learner for DML nuisance functions.
#'   Options: \code{"linear"} (default), \code{"glmnet"} / \code{"lasso"},
#'   \code{"ridge"}, \code{"elasticnet"}, \code{"ranger"} / \code{"rf"},
#'   \code{"grf"}, \code{"gam"}, \code{"nnet"}, and
#'   \code{"xgboost"} / \code{"boosting"}. Optional learners are checked
#'   lazily when requested.
#' @param K Integer; number of DML cross-fitting folds. Default 5.
#' @param S Integer; number of random splits for median aggregation. Default 3.
#' @param signal_map Signal-to-\code{G} mapping for \code{method="dml_flex"}.
#'   Options are \code{"local_poly"} (default), \code{"spline"},
#'   \code{"blp_spline"}, \code{"gam"}, and \code{"kernel"}. The
#'   \code{"blp_spline"} option stores a spline-basis covariance matrix used
#'   for covariance-aware contrasts and Gaussian max-\eqn{t} curve bands.
#' @param map_degree Integer polynomial degree for \code{signal_map="local_poly"}
#'   and spline degree for \code{signal_map="spline"}. Default 2.
#' @param map_df Optional integer degrees of freedom or basis size for
#'   \code{signal_map="spline"} and \code{signal_map="gam"}. Default \code{NULL}
#'   chooses a small data-dependent value.
#' @param density_method Conditional-density estimator for \code{method="dml_flex"}
#'   and \code{method="dml_incremental"}. Options are \code{"residual_kde"}
#'   (default), \code{"location_scale"}, and \code{"local_kde"}.
#' @param dml_inference Scalar DML interval rule. \code{"score"} (default)
#'   reports analytical score-based pointwise intervals. \code{"score_multiplier"}
#'   uses the stored score multiplier critical value for scalar DML intervals
#'   when available.
#' @param dml_boot Integer; multiplier draws for scalar DML score multiplier
#'   inference. Default 500.
#' @param dml_multiplier Multiplier law for scalar DML score inference:
#'   \code{"normal"} (default), \code{"wild"}, or \code{"exponential"}.
#'
#' @return A list with the following components:
#'   \item{est}{A list with three elements:
#'              \code{$pre}, \code{$event}, and \code{$post} containing
#'              aggregated pre-treatment, overall event, and post-treatment
#'              FDID estimates, respectively.}
#'   \item{dynamic}{Dynamic FDID estimates for each time in \code{entire_period}.}
#'   \item{raw_means}{Raw mean outcomes by group for each time in \code{entire_period}.}
#'   \item{tr_period}{Treatment periods used.}
#'   \item{ref_period}{Reference period used.}
#'   \item{entire_period}{All time periods for dynamic estimation.}
#'   \item{method}{Method used.}
#'   \item{vartype}{Variance type used.}
#'   \item{alpha}{Significance level used for confidence intervals.}
#'   \item{times}{All numeric time columns found.}
#'   \item{G}{Baseline factor; binary for binary methods and continuous for continuous-G methods.}
#'   \item{ps}{Propensity scores (if \code{ipw}, \code{aipw}, or \code{dml_binary} method used).}
#'   \item{curve_event}{Continuous-G curve fields for \code{method="kernel"} and
#'              \code{method="dml_flex"}. For \code{dml_flex}, this contains
#'              both \code{theta_hat} and \code{delta_hat}, pointwise
#'              confidence intervals, curve covariance matrices when available,
#'              and simultaneous bands when available. For reporting level contrasts and
#'              fixed-\code{g} derivatives, use \code{\link{fdid_contrast}} and
#'              \code{\link{fdid_derivative}}.}
#'   \item{incremental_event}{Scalar observed-average-derivative signal fields
#'              for \code{method="dml_incremental"}.}
#'   \item{scalar_event}{Scalar DML score, influence-function, and multiplier
#'              inference fields for scalar DML methods.}
#'   \item{dml_metadata}{DML target and inference metadata when available.}
#'   \item{dml_diagnostics}{DML signal, density, derivative, correction, and
#'              cluster diagnostics when available.}
#'   \item{call}{The matched call.}
#'   \item{target.pop}{Character indicating the target population used.}
#'
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
#' @import dplyr
#' @import estimatr
#' @importFrom grf causal_forest
#' @importFrom rlang sym
#' @importFrom tidyselect all_of
#' @importFrom doFuture registerDoFuture `%dofuture%`
#' @importFrom future plan multisession sequential
#' @importFrom stats cov density approxfun dnorm qnorm sd quantile median mad bw.SJ lm glm predict coef fitted complete.cases setNames
#' @importFrom grDevices adjustcolor
#' @importFrom graphics abline polygon lines
#' @author Rivka Lipkovitz, Enhan Liu
#' @export
fdid <- function(s,
                 tr_period,
                 ref_period,
                 entire_period = NULL,
                 method   = "ols1",
                 vartype  = "robust",
                 missing_data = c("listwise", "available"),
                 nsims    = 1000,
                 parallel = FALSE,
                 cores    = 2,
                 target.pop = c("all","1","0"),
                 eval_g   = NULL,
                 h0       = NULL,
                 K_folds  = 10L,
                 boot     = 500L,
                 alpha    = 0.05,
                 trim     = 0.05,
                 learner  = "linear",
                 K        = 5L,
                 S        = 3L,
                 signal_map = "local_poly",
                 map_degree = 2L,
                 map_df = NULL,
                 density_method = "residual_kde",
                 dml_inference = "score",
                 dml_boot = 500L,
                 dml_multiplier = "normal") {

  # Helper: convert a numeric time value to its Y_-prefixed column name
  ycol <- function(t) paste0("Y_", t)

  # Helper to suppress ebal output:
  silent_ebalance <- function(Treatment, X) {
    tf_con <- file(tempfile(), open = "wt")
    on.exit({
      sink(type = "message"); sink(type = "output")
      close(tf_con)
    }, add = TRUE)
    sink(tf_con, type = "output"); sink(tf_con, type = "message")

    out <- suppressWarnings(ebal::ebalance(Treatment = Treatment, X = X))

    return(out)
  }

  # Match arg:
  missing_data <- match.arg(missing_data)
  target.pop <- match.arg(target.pop)

  # --- Basic checks ---
  if (!requireNamespace("estimatr", quietly=TRUE)) {
    stop("Package 'estimatr' not installed. Please install it.")
  }
  if (!requireNamespace("car", quietly=TRUE)) {
    stop("Package 'car' not installed. Please install it.")
  }

  # If there's a cluster column, store its name
  cluster <- if ("c" %in% names(s)) "c" else NULL

  stopifnot(is.data.frame(s))

  valid_methods <- c("ols1", "ols2", "did", "ebal", "aipw", "ipw",
                     "kernel", "dml_binary", "dml_plr", "dml_flex",
                     "dml_incremental")
  if (!(method %in% valid_methods)) {
    stop("method must be one of: ", paste(valid_methods, collapse=", "))
  }

  valid_vtypes <- c("robust","bootstrap","jackknife")
  if (!(vartype %in% valid_vtypes)) {
    stop("vartype must be one of: ", paste(valid_vtypes, collapse=", "))
  }
  if (method == "kernel" && !(vartype %in% c("robust", "bootstrap"))) {
    stop("method='kernel' currently supports vartype='robust' or vartype='bootstrap' only.")
  }
  if (method %in% c("dml_binary", "dml_plr", "dml_flex", "dml_incremental") &&
      vartype != "robust") {
    stop("DML methods currently support vartype='robust' only.")
  }
  if (!"G" %in% names(s)) {
    stop("No column named 'G' found in s (the baseline factor).")
  }

  # Continuous-G and DML tuning parameters should fail early when impossible.
  # Method-specific dispatchers repeat sample-size checks after missing-data
  # filtering, but these checks catch malformed inputs before any estimation.
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single finite number strictly between 0 and 1.")
  }
  if (!is.numeric(trim) || length(trim) != 1L || !is.finite(trim) ||
      trim < 0 || trim >= 0.5) {
    stop("trim must be a single finite number in [0, 0.5).")
  }
  if (!is.null(eval_g)) {
    if (!is.numeric(eval_g) || length(eval_g) == 0L ||
        any(!is.finite(eval_g))) {
      stop("eval_g must be a nonempty finite numeric vector when supplied.")
    }
  }
  if (!is.null(h0) &&
      (!is.numeric(h0) || length(h0) != 1L || !is.finite(h0) || h0 <= 0)) {
    stop("h0 must be a single positive finite number when supplied.")
  }
  if (!is.numeric(boot) || length(boot) != 1L || !is.finite(boot) ||
      boot < 1) {
    stop("boot must be a positive integer-like value.")
  }
  if (!is.numeric(dml_boot) || length(dml_boot) != 1L ||
      !is.finite(dml_boot) || dml_boot < 1) {
    stop("dml_boot must be a positive integer-like value.")
  }
  if (!is.numeric(K_folds) || length(K_folds) != 1L || !is.finite(K_folds) ||
      K_folds < 2) {
    stop("K_folds must be an integer-like value of at least 2.")
  }
  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) || K < 2) {
    stop("K must be an integer-like value of at least 2.")
  }
  if (!is.numeric(S) || length(S) != 1L || !is.finite(S) || S < 1) {
    stop("S must be an integer-like value of at least 1.")
  }
  K_folds <- as.integer(K_folds)
  boot <- as.integer(boot)
  dml_boot <- as.integer(dml_boot)
  K <- as.integer(K)
  S <- as.integer(S)
  map_degree <- as.integer(map_degree)

  ref_period_check <- if (length(ref_period) > 1) max(ref_period) else ref_period
  Y_tr_cols_check <- ycol(tr_period)
  Y_ref_col_check <- ycol(ref_period_check)
  if (any(!(Y_tr_cols_check %in% names(s)))) {
    stop("Some treatment-time columns not found in data: ",
         paste(Y_tr_cols_check[!(Y_tr_cols_check %in% names(s))], collapse = ", "))
  }
  if (!(Y_ref_col_check %in% names(s))) {
    stop("Reference-time column '", Y_ref_col_check, "' not found in data.")
  }

  # ── Continuous-G methods: dispatch early and return complete fdid object ──
	  if (method == "kernel") {
	    return(run_kernel_method(
	      s = s, tr_period = tr_period, ref_period = ref_period,
	      entire_period = entire_period, eval_g = eval_g, h0 = h0,
	      K_folds = K_folds, boot = boot, alpha = alpha, trim = trim,
	      vartype = vartype, missing_data = match.arg(missing_data),
	      cluster = cluster,
	      the_call = match.call()
	    ))
	  }
  if (method %in% c("dml_binary", "dml_plr", "dml_flex", "dml_incremental")) {
    return(run_dml_method(
      s = s, tr_period = tr_period, ref_period = ref_period,
      entire_period = entire_period, method = method, learner = learner,
      K = K, S = S, eval_g = eval_g, trim = trim, boot = boot, alpha = alpha,
      signal_map = signal_map, map_degree = map_degree, map_df = map_df,
      density_method = density_method,
      dml_inference = dml_inference, dml_boot = dml_boot,
      dml_multiplier = dml_multiplier,
      missing_data = match.arg(missing_data), the_call = match.call()
    ))
  }
  if (method %in% c("did","ols1") && target.pop != "all"){
    warning("For method = '", method,
            "', target.pop != 'all' does not change the estimate; returning the same estimate for all target.pop values.")
  }

  # Covariates are those named x1, x2, ... (per your original code)
  covar <- grep("^x\\d+$", names(s), value=TRUE)

  # If user gave multiple ref_period, pick last
  if (length(ref_period) > 1) {
    ref_period <- max(ref_period)
  }

  # Identify outcome columns with Y_ prefix (e.g. Y_2021, Y_1958)
  all_y_cols <- grep("^Y_[-]?[0-9]+(\\.[0-9]+)?$", names(s), value = TRUE)
  numeric_times <- suppressWarnings(as.numeric(sub("^Y_", "", all_y_cols)))
  numeric_times <- numeric_times[!is.na(numeric_times)]
  numeric_times <- sort(unique(numeric_times))

  # Check user-supplied time columns
  Y_tr_cols <- ycol(tr_period)
  Y_ref_col <- ycol(ref_period)
  missing_tr <- setdiff(Y_tr_cols, names(s))
  if (length(missing_tr) > 0) {
    stop("Some treatment-time columns not found in data: ",
         paste(missing_tr, collapse=", "))
  }
  if (!(Y_ref_col %in% names(s))) {
    stop("Reference-time column '", Y_ref_col, "' not found in data.")
  }

  # Missing-data approach
  needed_cols_minimal <- c("G", covar)
  if (!is.null(cluster)) needed_cols_minimal <- c(needed_cols_minimal, cluster)

  # If entire_period is NULL, use all numeric_times
  if (!is.null(entire_period)) {
    all_times <- sort(intersect(numeric_times, entire_period))
  } else {
    all_times <- numeric_times
  }

  if (missing_data == "listwise") {
    needed_cols_all <- c(needed_cols_minimal, ycol(all_times))
    s <- s[complete.cases(s[, needed_cols_all, drop=FALSE]), ]
  } else {
    # "available": only remove missing in G, covars, cluster
    s <- s[complete.cases(s[, needed_cols_minimal, drop=FALSE]), ]
  }

  # --- AUTOMATICALLY ENCODE ANY FACTOR COVARIATES FOR IPW & AIPW ---
  if (method %in% c("aipw", "ipw") && length(covar) > 0) {
    # Check if any of these covars are not numeric
    is_non_num <- sapply(s[, covar, drop=FALSE], function(z) !is.numeric(z))
    if (any(is_non_num)) {
      # Build a formula for model.matrix
      form_covar <- as.formula(paste("~", paste(covar, collapse = " + ")))
      mm <- model.matrix(form_covar, data = s)
      mm <- mm[, -1, drop = FALSE]  # drop the intercept column

      # Drop old covar columns
      keep_cols <- setdiff(names(s), covar)
      # Combine
      s <- cbind(s[, keep_cols, drop=FALSE], mm)
      # Update covar to match new dummy columns
      covar <- colnames(mm)
    }
  }

  # Create single FDID outcome: average of treat periods minus reference
  Y_tr_cols_eff <- setdiff(unique(Y_tr_cols), Y_ref_col)
  tmp <- rowMeans(s[, Y_tr_cols_eff, drop = FALSE], na.rm = TRUE)
  s$tempY <- tmp - s[[Y_ref_col]]

  # ----- Parallel backend selector -----
  # mclapply (fork) is fast on Unix/Mac; multisession is too slow for small
  # tasks (~4ms/iter) due to IPC overhead. On Windows, fall back to sequential.
  use_mclapply <- parallel && (.Platform$OS.type == "unix")
  use_future   <- parallel && !use_mclapply &&
                    vartype %in% c("bootstrap", "jackknife")
  if (use_future) {
    if (!requireNamespace("foreach",  quietly = TRUE) ||
        !requireNamespace("doFuture", quietly = TRUE) ||
        !requireNamespace("future",   quietly = TRUE)) {
      stop("Parallel requires 'foreach', 'doFuture', and 'future' packages.")
    }
    future::plan(future::multisession, workers = cores)
    doFuture::registerDoFuture()
    on.exit(future::plan(future::sequential), add = TRUE)
  }

  # ----- Estimation Subroutines -----

  est_did <- function(d, covar) {
    fml <- tempY ~ G
    if (is.null(cluster)) {
      mod <- estimatr::lm_robust(fml, data=d, se_type="stata")
    } else {
      mod <- estimatr::lm_robust(fml, data=d, se_type="stata",
                                 clusters=d[[cluster]])
    }
    cfit <- summary(mod)$coefficients["G", c(1,2,5,6)]
    names(cfit) <- c("Estimate","Std.Error","CI_Lower","CI_Upper")
    cfit
  }

  est_ols1 <- function(d, covar) {
    # Y ~ G + X (no interactions)
    if (length(covar) > 0) {
      fml <- as.formula(paste("tempY ~ G +", paste(covar, collapse=" + ")))
    } else {
      fml <- tempY ~ G
    }
    if (is.null(cluster)) {
      mod <- estimatr::lm_robust(fml, data=d, se_type="stata")
    } else {
      mod <- estimatr::lm_robust(fml, data=d, se_type="stata",
                                 clusters=d[[cluster]])
    }
    cfit <- summary(mod)$coefficients["G", c(1,2,5,6)]
    names(cfit) <- c("Estimate","Std.Error","CI_Lower","CI_Upper")
    cfit
  }

  est_ols2 <- function(d, covar) {
    # Y ~ G * X (with interactions).
    if (length(covar) == 0) {
      return(est_did(d, covar))
    }

    # Prep data
    X <- scale(d[, covar, drop = FALSE])
    Y <- d$tempY
    Z <- d$G
    n <- nrow(X)

    # Fit a plain lm to get coefficients
    linmod <- lm(Y ~ Z * X)
    coefs  <- coef(linmod)

    # Names of the interaction terms ("Z:X1", "Z:X2", …)
    int_names <- paste0("Z:X", colnames(X))

    # Mean of each (scaled) covariate in each group
    X_bar1 <- colMeans(X[Z == 1, , drop = FALSE])
    X_bar0 <- colMeans(X[Z == 0, , drop = FALSE])

    xbar <- switch(target.pop,
                   "1"   = X_bar1,
                   "0"   = X_bar0,
                   "all" = rep(0, ncol(X))
    )

    # — Point estimate at chosen xbar:
    beta_Z  <- coefs["Z"]
    beta_ZX <- coefs[int_names]
    est_super <- beta_Z + sum(beta_ZX * xbar)

    # — Variance matrix
    if (is.null(cluster)) {
      V <- car::hccm(linmod)                      # HC0
    } else {
      # require sandwich for cluster‐robust
      if (!requireNamespace("sandwich", quietly = TRUE)) {
        stop("Please install the 'sandwich' package for clustered SEs")
      }
      V <- sandwich::vcovCL(linmod, cluster = d[[cluster]])
    }

      all_nms <- names(coefs)
      g <- numeric(length(coefs)); names(g) <- all_nms
      g["Z"] <- 1
      g[int_names] <- xbar

      var_super <- as.numeric(t(g) %*% V %*% g)

      # Super-population correction (Lin 2013): add γ'Σ_X γ / n
      # This term is non-zero even when target.pop = "all" (xbar = 0)
      # and ensures SEs match the paper's reported results.
      if (target.pop == "all") {
        gamma      <- coefs[int_names]
        sp_correction <- as.numeric(t(gamma) %*% cov(X) %*% gamma) / n
        var_super  <- var_super + sp_correction
      }

      se_super  <- sqrt(var_super)

      out <- c(
        Estimate  = est_super,
        Std.Error = se_super,
        CI_Lower  = est_super - 1.96 * se_super,
        CI_Upper  = est_super + 1.96 * se_super
      )

    names(out) <- c("Estimate","Std.Error","CI_Lower","CI_Upper")
    out
  }

  est_ebal <- function(d, covar) {
    if (!requireNamespace("ebal", quietly=TRUE)) {
      stop("Method 'ebal' requires the 'ebal' package.")
    }
    # ebal code that reweights controls to match the distribution of treated => ATT
    if (target.pop != "1") {
      stop("Currently, 'ebal' only implemented for target.pop='1'. Please set target.pop='1' or choose another method.")
    }
    if (length(covar) == 0) {
      return(est_did(d, covar))
    }
    group1 <- d$G == 1
    X <- d[, covar, drop=FALSE]
    eout <- silent_ebalance(Treatment = group1, X = X)

    wts <- rep(0, nrow(d))
    wts[group1] <- 1
    wts[!group1] <- eout$w

    fml <- tempY ~ G
    if (is.null(cluster)) {
      mod <- estimatr::lm_robust(fml, data=d, weights=wts, se_type="stata")
    } else {
      mod <- estimatr::lm_robust(fml, data=d, weights=wts, se_type="stata",
                                 clusters=d[[cluster]])
    }
    cfit <- summary(mod)$coefficients["G", c(1,2,5,6)]
    names(cfit) <- c("Estimate","Std.Error","CI_Lower","CI_Upper")
    cfit
  }

  # IPW estimation using grf::probability_forest
  est_ipw <- function(d, covar) {
    if (!requireNamespace("grf", quietly = TRUE)) {
      stop("Method 'ipw' requires the 'grf' package.")
    }

    X <- d[, covar, drop = FALSE]
    Y <- d$tempY
    W <- d$G

    # 1) propensity score
    if (length(covar) == 0) {
      ps <- rep(mean(W, na.rm = TRUE), nrow(d))
    } else {
      pf <- grf::probability_forest(X = X, Y = as.factor(W), seed = 1234)
      ps <- pf$predictions[, 2]
    }

    # avoid 0/1 probabilities (infinite weights)
    eps <- 1e-6
    ps <- pmin(pmax(as.numeric(ps), eps), 1 - eps)

    # 2) weights: ATT vs ATE
    wts <- switch(target.pop,
                  "1"   = ifelse(W == 1, 1, ps / (1 - ps)),          # ATT
                  "0"   = ifelse(W == 1, (1 - ps) / ps, 1),          # ATC
                  "all" = ifelse(W == 1, 1 / ps, 1 / (1 - ps))       # ATE
    )

    # 3) weighted regression (robust SE, optional cluster like other est_*)
    fml <- tempY ~ G
    if (is.null(cluster)) {
      mod <- estimatr::lm_robust(fml, data = d, weights = wts, se_type = "stata")
    } else {
      mod <- estimatr::lm_robust(
        fml, data = d, weights = wts, se_type = "stata",
        clusters = d[[cluster]]
      )
    }

    # align with aipw style: construct a named vector 'ipw'
    cfit <- summary(mod)$coefficients["G", c(1, 2)]
    ipw  <- c(estimate = as.numeric(cfit[1]), std.err = as.numeric(cfit[2]))

    est <- as.numeric(ipw["estimate"])
    se_ <- as.numeric(ipw["std.err"])
    c(
      Estimate  = est,
      Std.Error = se_,
      CI_Lower  = est - 1.96 * se_,
      CI_Upper  = est + 1.96 * se_
    )
  }

  # AIPW estimation using grf::causal_forest
  est_aipw <- function(d, covar) {
    if (!requireNamespace("grf", quietly = TRUE)) {
      stop("Method 'aipw' requires the 'grf' package.")
    }
    X <- d[, covar, drop = FALSE]
    Y <- d$tempY
    W <- d$G

    c.forest <- grf::causal_forest(X, Y, W)
    att      <- grf::average_treatment_effect(
      c.forest,
      target.sample = switch(target.pop, "1" = "treated", "0" = "control", "all" = "all"),
      method        = "AIPW"
    )
    est <- as.numeric(att["estimate"])
    se_ <- as.numeric(att["std.err"])
    c(
      Estimate  = est,
      Std.Error = se_,
      CI_Lower  = est - 1.96 * se_,
      CI_Upper  = est + 1.96 * se_
    )
  }


  run_method <- function(d, method, covar) {
    if (method=="did")  return(est_did(d, covar))
    if (method=="ols1") return(est_ols1(d, covar))
    if (method=="ols2") return(est_ols2(d, covar))
    if (method=="ebal") return(est_ebal(d, covar))
    if (method=="ipw")  return(est_ipw(d, covar))
    if (method=="aipw") return(est_aipw(d, covar))
    stop("Unknown method: ", method)
  }

  # ----- Variance wrappers -----

  do_robust <- function(d, method, covar) {
    run_method(d, method, covar)
  }

  do_bootstrap <- function(d, method, covar, B = nsims) {
    main_est <- run_method(d, method, covar)["Estimate"]
    one_boot <- function(i) {
      idx <- sample.int(nrow(d), replace = TRUE)
      db  <- d[idx, ]
      run_method(db, method, covar)["Estimate"]
    }
    if (!parallel) {
      vals <- replicate(B, one_boot(1), simplify = TRUE)
    } else if (use_mclapply) {
      vals <- unlist(parallel::mclapply(1:B, one_boot, mc.cores = cores,
                                        mc.set.seed = TRUE))
    } else {
      vals <- foreach::foreach(
        i = 1:B,
        .combine = c,
        .options.future = list(seed = TRUE)
      ) %dofuture% {
        one_boot(i)
      }
    }
    se_        <- sd(vals, na.rm = TRUE)
    boot_mean  <- mean(vals, na.rm = TRUE)
    q          <- quantile(vals, probs = c(0.025, 0.975), na.rm = TRUE)
    ci_lower   <- main_est - (q[2] - boot_mean)   # bias-corrected
    ci_upper   <- main_est + (boot_mean - q[1])
    c(Estimate=main_est, Std.Error=se_, CI_Lower=ci_lower, CI_Upper=ci_upper)
  }

  do_jackknife <- function(d, method, covar) {
    n <- nrow(d)
    main_est <- run_method(d, method, covar)["Estimate"]
    one_jack <- function(i) {
      d_j <- d[-i, ]
      run_method(d_j, method, covar)["Estimate"]
    }
    if (!parallel) {
      jvals <- numeric(n)
      for (i in seq_len(n)) {
        jvals[i] <- one_jack(i)
      }
    } else if (use_mclapply) {
      jvals <- unlist(parallel::mclapply(1:n, one_jack, mc.cores = cores,
                                         mc.set.seed = TRUE))
    } else {
      jvals <- foreach::foreach(
        i = 1:n,
        .combine = c,
        .options.future = list(seed = TRUE)
      ) %dofuture% {
        one_jack(i)
      }
    }
    est_mean <- mean(jvals, na.rm=TRUE)
    se_ <- sqrt((n - 1)/n * sum((jvals - est_mean)^2))
    #ci_ <- quantile(jvals, probs=c(0.025, 0.975), na.rm=TRUE)
    alpha <- 0.05
    z_crit <- qnorm(1 - alpha/2)
    ci_lower <- main_est - z_crit * se_
    ci_upper <- main_est + z_crit * se_
    c(Estimate=main_est, Std.Error=se_, CI_Lower=ci_lower, CI_Upper=ci_upper)
  }

  get_estimate <- function(d, method, vartype, covar) {
    if (vartype == "robust") {
      do_robust(d, method, covar)
    } else if (vartype == "bootstrap") {
      do_bootstrap(d, method, covar, nsims)
    } else if (vartype == "jackknife") {
      do_jackknife(d, method, covar)
    }
  }

  # ----- Single FDID estimate (static) -----
  static_result_vec <- get_estimate(s, method, vartype, covar)
  static_result <- as.data.frame(t(static_result_vec))

  # ----- Dynamic FDID for each time column -----
  dynamic_times <- all_times
  dynamic_df <- data.frame(
    Estimate  = numeric(length(dynamic_times)),
    Std.Error = numeric(length(dynamic_times)),
    CI_Lower  = numeric(length(dynamic_times)),
    CI_Upper  = numeric(length(dynamic_times)),
    row.names = as.character(dynamic_times)
  )
  last_ref_num <- as.numeric(ref_period)

  for (i in seq_along(dynamic_times)) {
    t_yr <- dynamic_times[i]
    if (t_yr == last_ref_num) {
      # By definition, the reference period effect is 0
      dynamic_df[i, ] <- 0
    } else {
      t_col <- ycol(t_yr)
      ref_col <- ycol(last_ref_num)
      d_dyn <- s
      d_dyn$tempY <- d_dyn[[t_col]] - d_dyn[[ref_col]]
      d_dyn <- d_dyn[!is.na(d_dyn$tempY), , drop=FALSE]
      est_vec <- get_estimate(d_dyn, method, vartype, covar)
      dynamic_df[i, ] <- est_vec
    }
  }

  # ----- Raw means by group & time (with 95% CI) -----
  zval <- stats::qnorm(0.975)

  rawdf <- data.frame(
    time     = numeric(0),
    group    = character(0),
    meanY    = numeric(0),
    n        = integer(0),
    seY      = numeric(0),
    CI_Lower = numeric(0),
    CI_Upper = numeric(0)
  )



  for (yr in dynamic_times) {
    col_yr <- ycol(yr)
    if (!col_yr %in% names(s)) {
      warning("Column '", col_yr, "' not found in data. Skipping raw_means.")
      next
    }

    y1 <- s[[col_yr]][s$G == 1]
    y0 <- s[[col_yr]][s$G == 0]

    n1 <- sum(!is.na(y1)); m1 <- mean(y1, na.rm = TRUE)
    n0 <- sum(!is.na(y0)); m0 <- mean(y0, na.rm = TRUE)

    se1 <- if (n1 <= 1) NA_real_ else stats::sd(y1, na.rm = TRUE) / sqrt(n1)
    se0 <- if (n0 <= 1) NA_real_ else stats::sd(y0, na.rm = TRUE) / sqrt(n0)

    rawdf <- rbind(
      rawdf,
      data.frame(
        time     = yr,
        group    = "Group 1",
        meanY    = m1,
        n        = n1,
        seY      = se1,
        CI_Lower = m1 - zval * se1,
        CI_Upper = m1 + zval * se1
      ),
      data.frame(
        time     = yr,
        group    = "Group 0",
        meanY    = m0,
        n        = n0,
        seY      = se0,
        CI_Lower = m0 - zval * se0,
        CI_Upper = m0 + zval * se0
      )
    )
  }

  # ----- Propensity Scores (only for IPW & AIPW) -----
  ps_vec <- NULL
  if (method %in% c("aipw", "ipw")) {
    if (length(covar) == 0) {
      ps_vec <- rep(mean(s$G, na.rm=TRUE), nrow(s))
    } else if (method == "ipw") {
      # match IPW estimation ps (GRF)
      if (!requireNamespace("grf", quietly = TRUE)) {
        stop("Method 'ipw' requires the 'grf' package.")
      }
      X <- s[, covar, drop = FALSE]
      pf <- grf::probability_forest(X = X, Y = as.factor(s$G), seed = 1234)
      ps_vec <- as.numeric(pf$predictions[, 2])
    } else {
      # AIPW output ps (GLM)
      ps_formula <- as.formula(paste("G ~", paste(covar, collapse=" + ")))
      ps_fit <- glm(ps_formula, data=s, family=binomial)
      ps_vec <- predict(ps_fit, newdata=s, type="response")
    }
  }

  # ----- Pre-event & Post-event aggregates -----
  earliest_event_time <- min(tr_period)
  latest_event_time   <- max(tr_period)

  get_aggregate_estimate <- function(sdata, time_vec, ref_col) {
    if (length(time_vec) == 0) {
      return(c(Estimate=NA, Std.Error=NA, CI_Lower=NA, CI_Upper=NA))
    }

    time_cols <- setdiff(unique(ycol(time_vec)), ref_col)
    if (length(time_cols) == 0L) {
      return(c(Estimate=NA, Std.Error=NA, CI_Lower=NA, CI_Upper=NA))
    }

    d_agg <- sdata

    tmp <- rowMeans(d_agg[, time_cols, drop = FALSE], na.rm = TRUE)
    tmp[is.nan(tmp)] <- NA_real_
    d_agg$tempY <- tmp - d_agg[[ref_col]]

    d_agg <- d_agg[!is.na(d_agg$tempY), , drop=FALSE]
    get_estimate(d_agg, method, vartype, covar)
  }

  pre_times <- dynamic_times[dynamic_times < earliest_event_time]
  post_times <- dynamic_times[dynamic_times > latest_event_time]
  pre_event_vec  <- get_aggregate_estimate(s, pre_times, Y_ref_col)
  post_event_vec <- get_aggregate_estimate(s, post_times, Y_ref_col)

  pre_event_result  <- as.data.frame(t(pre_event_vec))
  event_result      <- static_result
  post_event_result <- as.data.frame(t(post_event_vec))

  forced_cols <- c("Estimate", "Std.Error", "CI_Lower", "CI_Upper")
  colnames(pre_event_result)   <- forced_cols
  colnames(event_result)       <- forced_cols
  colnames(post_event_result)  <- forced_cols

  est_list <- list(
    pre   = pre_event_result,
    event = event_result,
    post  = post_event_result
  )

  # Compile output
  out <- list(
    est           = est_list,
    dynamic       = dynamic_df,
    raw_means     = rawdf,
    tr_period     = tr_period,
    ref_period    = ref_period,
    entire_period = dynamic_times,
    method        = method,
    vartype       = vartype,
    alpha         = alpha,
    times         = numeric_times,
    G             = s$G,
    ps            = ps_vec,
    call          = match.call(),
    target.pop    = target.pop
  )
  class(out) <- "fdid"
  return(out)
}
