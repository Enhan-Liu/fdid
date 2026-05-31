# DML estimators for continuous-G FDID
# Algorithms: Chernozhukov et al. (2018) §3 (PLR), §5 (IRM/AIPW),
#             Semenova & Chernozhukov (2021) Example 2.1 (flex)
# Called internally by fdid() when method ∈ {"dml_binary","dml_plr",
#                                            "dml_flex","dml_incremental"}

# ── Learner helpers ───────────────────────────────────────────────────────────

.valid_dml_learners <- paste(
  "'linear', 'glmnet'/'lasso', 'ridge', 'elasticnet', 'ranger'/'rf',",
  "'grf', 'gam', 'nnet', or 'xgboost'/'boosting'"
)

.normalize_learner <- function(learner) {
  if (!is.character(learner) || length(learner) != 1L || is.na(learner)) {
    stop("learner must be a single string. Choose: ", .valid_dml_learners, ".",
         call. = FALSE)
  }
  ll <- gsub("[[:space:]_-]+", "", tolower(learner))
  if (ll %in% c("default", "linear", "logistic", "lm", "ols", "l", "d")) {
    return("linear")
  }
  if (ll %in% c("glmnet", "lasso", "regularization", "regularisation", "r")) {
    return("lasso")
  }
  if (ll %in% c("ridge")) {
    return("ridge")
  }
  if (ll %in% c("elasticnet", "enet")) {
    return("elasticnet")
  }
  if (ll %in% c("ranger", "rf", "randomforest", "forest")) {
    return("ranger")
  }
  if (ll %in% c("grf", "regressionforest", "probabilityforest")) {
    return("grf")
  }
  if (ll %in% c("gam", "mgcv")) {
    return("gam")
  }
  if (ll %in% c("nnet", "nn", "network", "neuralnetwork")) {
    return("nnet")
  }
  if (ll %in% c("xgboost", "xgb", "boosting", "boost", "gradientboosting",
                "gradientboost", "histgradientboosting", "hgb")) {
    return("xgboost")
  }
  stop("Unknown learner: '", learner, "'. Choose: ", .valid_dml_learners, ".",
       call. = FALSE)
}

.glmnet_alpha <- function(type) {
  switch(type,
         lasso = 1,
         ridge = 0,
         elasticnet = 0.5,
         1)
}

.cv_nfolds <- function(n) {
  max(3L, min(5L, floor(n / 2L)))
}

.as_feature_df <- function(X, nm = NULL) {
  if (is.null(nm)) nm <- paste0("v", seq_len(ncol(X)))
  setNames(as.data.frame(X), nm)
}

.glmnet_design <- function(X) {
  X <- as.matrix(X)
  if (ncol(X) == 1L) cbind(X, .fdid_zero = 0) else X
}

.scale_train <- function(X) {
  center <- colMeans(X, na.rm = TRUE)
  scale <- apply(X, 2, stats::sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale < 1e-8] <- 1
  Xs <- sweep(sweep(X, 2, center, "-"), 2, scale, "/")
  list(X = Xs, center = center, scale = scale)
}

.scale_new <- function(X, center, scale) {
  sweep(sweep(X, 2, center, "-"), 2, scale, "/")
}

.gam_formula <- function(response, df, family = c("gaussian", "binomial")) {
  family <- match.arg(family)
  vars <- setdiff(names(df), response)
  if (length(vars) == 0L) return(stats::as.formula(paste(response, "~ 1")))
  terms <- vapply(vars, function(v) {
    n_unique <- length(unique(df[[v]][is.finite(df[[v]])]))
    if (n_unique >= 6L) {
      paste0("s(", v, ", k = ", min(5L, n_unique - 1L), ")")
    } else {
      v
    }
  }, character(1L))
  stats::as.formula(paste(response, "~", paste(terms, collapse = " + ")))
}

.fit_reg_model <- function(Y_tr, X_tr, learner) {
  default <- mean(Y_tr, na.rm = TRUE)
  if (!is.finite(default)) default <- 0
  if (ncol(X_tr) == 0L || nrow(X_tr) < 5L) {
    return(list(type = "mean", value = default, p = ncol(X_tr)))
  }

  X_tr <- as.matrix(X_tr)
  nm <- paste0("v", seq_len(ncol(X_tr)))
  df_tr <- data.frame(Y = Y_tr, .as_feature_df(X_tr, nm))
  ll <- .normalize_learner(learner)

  if (ll == "linear") {
    return(list(type = "lm", model = stats::lm(Y ~ ., data = df_tr),
                nm = nm, p = ncol(X_tr), value = default))
  }

  if (ll %in% c("lasso", "ridge", "elasticnet")) {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("learner='", learner, "' requires the 'glmnet' package.",
           call. = FALSE)
    }
    X_fit <- .glmnet_design(X_tr)
    mod <- glmnet::cv.glmnet(
      X_fit, Y_tr, alpha = .glmnet_alpha(ll), nfolds = .cv_nfolds(nrow(X_tr))
    )
    return(list(type = "glmnet", model = mod, p = ncol(X_tr), value = default,
                glmnet_pad = ncol(X_tr) == 1L))
  }

  if (ll == "ranger") {
    if (!requireNamespace("ranger", quietly = TRUE)) {
      stop("learner='", learner, "' requires the 'ranger' package.",
           call. = FALSE)
    }
    mod <- ranger::ranger(Y ~ ., data = df_tr, num.trees = 200,
                          min.node.size = 5)
    return(list(type = "ranger", model = mod, nm = nm, p = ncol(X_tr),
                value = default))
  }

  if (ll == "grf") {
    if (!requireNamespace("grf", quietly = TRUE)) {
      stop("learner='grf' requires the 'grf' package.", call. = FALSE)
    }
    mod <- grf::regression_forest(as.matrix(X_tr), Y_tr, num.trees = 500)
    return(list(type = "grf_regression", model = mod, p = ncol(X_tr),
                value = default))
  }

  if (ll == "gam") {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop("learner='gam' requires the 'mgcv' package.", call. = FALSE)
    }
    fml <- .gam_formula("Y", df_tr, family = "gaussian")
    mod <- mgcv::gam(fml, data = df_tr, method = "REML")
    return(list(type = "gam", model = mod, nm = nm, p = ncol(X_tr),
                value = default))
  }

  if (ll == "nnet") {
    if (!requireNamespace("nnet", quietly = TRUE)) {
      stop("learner='nnet' requires the 'nnet' package.", call. = FALSE)
    }
    sc <- .scale_train(X_tr)
    size <- min(5L, max(1L, ncol(X_tr)))
    mod <- nnet::nnet(x = sc$X, y = Y_tr, size = size, linout = TRUE,
                      decay = 0.01, maxit = 500, trace = FALSE,
                      MaxNWts = 10000)
    return(list(type = "nnet_regression", model = mod, center = sc$center,
                scale = sc$scale, p = ncol(X_tr), value = default))
  }

  if (ll == "xgboost") {
    if (!requireNamespace("xgboost", quietly = TRUE)) {
      stop("learner='", learner, "' requires the 'xgboost' package.",
           call. = FALSE)
    }
    mod <- .fit_xgboost_model(X_tr, Y_tr, objective = "reg:squarederror")
    return(list(type = "xgboost", model = mod, p = ncol(X_tr),
                value = default))
  }

  stop("Unknown learner: '", learner, "'. Choose: ", .valid_dml_learners, ".",
       call. = FALSE)
}

.predict_reg_model <- function(model, X_new) {
  X_new <- as.matrix(X_new)
  n_new <- nrow(X_new)
  out <- switch(model$type,
    mean = rep(model$value, n_new),
    lm = as.numeric(stats::predict(model$model,
                                   newdata = .as_feature_df(X_new, model$nm))),
    glmnet = as.numeric(stats::predict(model$model, newx = .glmnet_design(X_new),
                                       s = "lambda.min")),
    ranger = as.numeric(stats::predict(model$model,
                                       data = .as_feature_df(X_new, model$nm))$predictions),
    grf_regression = as.numeric(stats::predict(model$model,
                                               as.matrix(X_new))$predictions),
    gam = as.numeric(stats::predict(model$model,
                                    newdata = .as_feature_df(X_new, model$nm),
                                    type = "response")),
    nnet_regression = as.numeric(stats::predict(
      model$model, .scale_new(X_new, model$center, model$scale), type = "raw"
    )),
    xgboost = as.numeric(stats::predict(model$model, as.matrix(X_new))),
    stop("Internal error: unknown regression model type '", model$type, "'.",
         call. = FALSE)
  )
  out[!is.finite(out)] <- model$value
  out
}

.fit_reg <- function(Y_tr, X_tr, X_new, learner) {
  .predict_reg_model(.fit_reg_model(Y_tr, X_tr, learner), X_new)
}

.fit_xgboost_model <- function(X_tr, y_tr, objective) {
  dtrain <- xgboost::xgb.DMatrix(as.matrix(X_tr), label = y_tr)
  xgboost::xgb.train(
    params = list(
      objective = objective,
      eta = 0.05,
      max_depth = 3L,
      min_child_weight = 5,
      subsample = 0.9,
      colsample_bytree = 0.9,
      nthread = 1L,
      verbosity = 0L
    ),
    data = dtrain,
    nrounds = 100L
  )
}

.fit_classif_model <- function(G_tr, X_tr, learner) {
  ps0 <- mean(as.integer(G_tr == 1), na.rm = TRUE)
  if (!is.finite(ps0)) ps0 <- 0.5
  if (ncol(X_tr) == 0L || nrow(X_tr) < 5L || length(unique(G_tr)) < 2L) {
    return(list(type = "mean_classif", value = ps0, p = ncol(X_tr)))
  }

  X_tr <- as.matrix(X_tr)
  y01 <- as.integer(G_tr == 1)
  nm <- paste0("v", seq_len(ncol(X_tr)))
  df_tr <- data.frame(G = y01, .as_feature_df(X_tr, nm))
  ll <- .normalize_learner(learner)

  if (ll == "linear") {
    return(list(type = "glm_binomial",
                model = stats::glm(G ~ ., data = df_tr,
                                   family = stats::binomial),
                nm = nm, p = ncol(X_tr), value = ps0))
  }

  if (ll %in% c("lasso", "ridge", "elasticnet")) {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("learner='", learner, "' requires the 'glmnet' package.",
           call. = FALSE)
    }
    X_fit <- .glmnet_design(X_tr)
    mod <- glmnet::cv.glmnet(
      X_fit, y01, family = "binomial", alpha = .glmnet_alpha(ll),
      nfolds = .cv_nfolds(nrow(X_tr))
    )
    return(list(type = "glmnet_binomial", model = mod, p = ncol(X_tr),
                value = ps0, glmnet_pad = ncol(X_tr) == 1L))
  }

  if (ll == "ranger") {
    if (!requireNamespace("ranger", quietly = TRUE)) {
      stop("learner='", learner, "' requires the 'ranger' package.",
           call. = FALSE)
    }
    df_rf <- df_tr
    df_rf$G <- as.factor(df_rf$G)
    mod <- ranger::ranger(G ~ ., data = df_rf, probability = TRUE,
                          num.trees = 200)
    return(list(type = "ranger_probability", model = mod, nm = nm,
                p = ncol(X_tr), value = ps0))
  }

  if (ll == "grf") {
    if (!requireNamespace("grf", quietly = TRUE)) {
      stop("learner='grf' requires the 'grf' package.", call. = FALSE)
    }
    mod <- grf::regression_forest(as.matrix(X_tr), y01, num.trees = 500)
    return(list(type = "grf_probability", model = mod, p = ncol(X_tr),
                value = ps0))
  }

  if (ll == "gam") {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop("learner='gam' requires the 'mgcv' package.", call. = FALSE)
    }
    fml <- .gam_formula("G", df_tr, family = "binomial")
    mod <- mgcv::gam(fml, data = df_tr, family = stats::binomial,
                     method = "REML")
    return(list(type = "gam_binomial", model = mod, nm = nm,
                p = ncol(X_tr), value = ps0))
  }

  if (ll == "nnet") {
    if (!requireNamespace("nnet", quietly = TRUE)) {
      stop("learner='nnet' requires the 'nnet' package.", call. = FALSE)
    }
    sc <- .scale_train(X_tr)
    size <- min(5L, max(1L, ncol(X_tr)))
    mod <- nnet::nnet(x = sc$X, y = y01, size = size, entropy = TRUE,
                      decay = 0.01, maxit = 500, trace = FALSE,
                      MaxNWts = 10000)
    return(list(type = "nnet_probability", model = mod, center = sc$center,
                scale = sc$scale, p = ncol(X_tr), value = ps0))
  }

  if (ll == "xgboost") {
    if (!requireNamespace("xgboost", quietly = TRUE)) {
      stop("learner='", learner, "' requires the 'xgboost' package.",
           call. = FALSE)
    }
    mod <- .fit_xgboost_model(X_tr, y01, objective = "binary:logistic")
    return(list(type = "xgboost_probability", model = mod, p = ncol(X_tr),
                value = ps0))
  }

  stop("Unknown learner: '", learner, "'. Choose: ", .valid_dml_learners, ".",
       call. = FALSE)
}

.predict_classif_model <- function(model, X_new) {
  X_new <- as.matrix(X_new)
  n_new <- nrow(X_new)
  ps <- switch(model$type,
    mean_classif = rep(model$value, n_new),
    glm_binomial = as.numeric(stats::predict(model$model,
                                             newdata = .as_feature_df(X_new, model$nm),
                                             type = "response")),
    glmnet_binomial = as.numeric(stats::predict(model$model,
                                                newx = .glmnet_design(X_new),
                                                s = "lambda.min",
                                                type = "response")),
    ranger_probability = {
      pmat <- stats::predict(model$model,
                             data = .as_feature_df(X_new, model$nm))$predictions
      lvl1 <- which(colnames(pmat) == "1")
      if (length(lvl1) == 0L) lvl1 <- ncol(pmat)
      as.numeric(pmat[, lvl1])
    },
    grf_probability = as.numeric(stats::predict(model$model,
                                                as.matrix(X_new))$predictions),
    gam_binomial = as.numeric(stats::predict(model$model,
                                             newdata = .as_feature_df(X_new, model$nm),
                                             type = "response")),
    nnet_probability = as.numeric(stats::predict(
      model$model, .scale_new(X_new, model$center, model$scale), type = "raw"
    )),
    xgboost_probability = as.numeric(stats::predict(model$model,
                                                    as.matrix(X_new))),
    stop("Internal error: unknown classification model type '", model$type, "'.",
         call. = FALSE)
  )
  ps[!is.finite(ps)] <- model$value
  pmin(pmax(ps, 0.01), 0.99)
}

.fit_classif <- function(G_tr, X_tr, X_new, learner) {
  .predict_classif_model(.fit_classif_model(G_tr, X_tr, learner), X_new)
}

.dml_score_variance <- function(score, center = 0, jacobian = 1, cluster = NULL) {
  score <- as.numeric(score)
  ok <- is.finite(score)
  if (!any(ok)) return(NA_real_)
  score <- score[ok] - center
  n <- length(score)
  jacobian <- as.numeric(jacobian)
  if (!is.finite(jacobian) || abs(jacobian) < 1e-10) return(Inf)
  infl <- score / jacobian
  if (is.null(cluster)) {
    return(sum(infl^2, na.rm = TRUE) / (n^2))
  }
  cl <- cluster[ok]
  if (length(cl) != n || all(is.na(cl))) {
    return(sum(infl^2, na.rm = TRUE) / (n^2))
  }
  cl <- as.factor(cl)
  sums <- rowsum(infl, cl, reorder = FALSE)[, 1]
  g <- length(sums)
  if (g <= 1L) return(sum(infl^2, na.rm = TRUE) / (n^2))
  (g / (g - 1)) * sum(sums^2, na.rm = TRUE) / (n^2)
}

.diagnostic_quantiles <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(c(min = NA_real_, p01 = NA_real_, p05 = NA_real_,
             p50 = NA_real_, p95 = NA_real_, p99 = NA_real_, max = NA_real_))
  }
  c(
    min = min(x),
    p01 = unname(stats::quantile(x, 0.01, names = FALSE)),
    p05 = unname(stats::quantile(x, 0.05, names = FALSE)),
    p50 = unname(stats::quantile(x, 0.50, names = FALSE)),
    p95 = unname(stats::quantile(x, 0.95, names = FALSE)),
    p99 = unname(stats::quantile(x, 0.99, names = FALSE)),
    max = max(x)
  )
}

.dml_diagnostics <- function(signal = NULL, density = NULL, d_log_density = NULL,
                             d_mu = NULL, correction = NULL, cluster = NULL,
                             density_method = NULL, learner = NULL) {
  list(
    learner = learner,
    density_method = density_method,
    n = if (!is.null(signal)) length(signal) else if (!is.null(density)) length(density) else NA_integer_,
    n_clusters = if (!is.null(cluster)) length(unique(cluster[!is.na(cluster)])) else NA_integer_,
    signal_quantiles = .diagnostic_quantiles(signal),
    density_quantiles = .diagnostic_quantiles(density),
    d_log_density_quantiles = .diagnostic_quantiles(d_log_density),
    d_mu_quantiles = .diagnostic_quantiles(d_mu),
    correction_quantiles = .diagnostic_quantiles(correction),
    min_density = if (!is.null(density)) suppressWarnings(min(density, na.rm = TRUE)) else NA_real_,
    share_density_below_1e_3 = if (!is.null(density)) mean(density < 1e-3, na.rm = TRUE) else NA_real_,
    share_density_below_1e_2 = if (!is.null(density)) mean(density < 1e-2, na.rm = TRUE) else NA_real_
  )
}

# ── Algorithm A: dml_binary — IRM / AIPW ─────────────────────────────────────
# Chernozhukov et al. (2018) §5, DML2 with median aggregation over S splits

.est_dml_binary <- function(DeltaY, G, X, K, S, alpha, learner, cluster = NULL,
                            dml_inference = "score", dml_boot = 500L,
                            dml_multiplier = "normal") {
  n      <- length(DeltaY)
  g_vals <- sort(unique(G))
  if (!identical(g_vals, c(0, 1)))
    stop("method='dml_binary' requires binary G in {0, 1}.")

  tau_s  <- numeric(S)
  var_s  <- numeric(S)
  ps_mat <- matrix(NA_real_, n, S)
  psi_mat <- matrix(NA_real_, n, S)
  influence_mat <- matrix(NA_real_, n, S)
  mu0_global <- mean(DeltaY[G == 0], na.rm = TRUE)
  mu1_global <- mean(DeltaY[G == 1], na.rm = TRUE)
  pi_global  <- mean(G == 1, na.rm = TRUE)

  for (ss in seq_len(S)) {
    fold_id <- sample(rep(seq_len(K), length.out = n))
    mu0_hat <- numeric(n)
    mu1_hat <- numeric(n)
    pi_hat  <- numeric(n)

    for (k in seq_len(K)) {
      tr  <- which(fold_id != k)
      ho  <- which(fold_id == k)
      X_tr <- if (ncol(X) > 0) X[tr, , drop = FALSE] else matrix(0, length(tr), 0)
      X_ho <- if (ncol(X) > 0) X[ho, , drop = FALSE] else matrix(0, length(ho), 0)

      idx0 <- tr[G[tr] == 0]; idx1 <- tr[G[tr] == 1]
      X0   <- if (ncol(X) > 0) X[idx0, , drop = FALSE] else matrix(0, length(idx0), 0)
      X1   <- if (ncol(X) > 0) X[idx1, , drop = FALSE] else matrix(0, length(idx1), 0)

      mu0_hat[ho] <- if (length(idx0) >= 5) .fit_reg(DeltaY[idx0], X0, X_ho, learner) else
                       rep(mu0_global, length(ho))
      mu1_hat[ho] <- if (length(idx1) >= 5) .fit_reg(DeltaY[idx1], X1, X_ho, learner) else
                       rep(mu1_global, length(ho))
      pi_hat[ho]  <- if (length(unique(G[tr])) == 2 && ncol(X_tr) > 0)
                       .fit_classif(G[tr], X_tr, X_ho, learner) else
                       rep(pi_global, length(ho))
    }

    pi_hat <- pmin(pmax(pi_hat, 0.01), 0.99)
    psi_i  <- mu1_hat - mu0_hat +
              as.numeric(G == 1) / pi_hat * (DeltaY - mu1_hat) -
              as.numeric(G == 0) / (1 - pi_hat) * (DeltaY - mu0_hat)

    tau_s[ss] <- mean(psi_i)
    var_s[ss] <- .dml_score_variance(psi_i, center = tau_s[ss], cluster = cluster)
    ps_mat[, ss] <- pi_hat
    psi_mat[, ss] <- psi_i
    influence_mat[, ss] <- psi_i - tau_s[ss]
  }

  tau_f <- median(tau_s)
  se_f  <- sqrt(max(0, median(var_s + (tau_s - tau_f)^2)))
  result <- c(Estimate  = tau_f,
              Std.Error = se_f,
              CI_Lower  = tau_f - stats::qnorm(1 - alpha/2) * se_f,
              CI_Upper  = tau_f + stats::qnorm(1 - alpha/2) * se_f)
  inf <- rowMeans(influence_mat, na.rm = TRUE)
  out <- .dml_attach_scalar_inference(
    result = result,
    estimate = tau_f,
    se = se_f,
    alpha = alpha,
    influence = inf,
    split_estimates = tau_s,
    dml_inference = dml_inference,
    dml_boot = dml_boot,
    dml_multiplier = dml_multiplier,
    cluster = cluster
  )
  out$ps <- rowMeans(ps_mat, na.rm = TRUE)
  out$signal <- rowMeans(psi_mat, na.rm = TRUE)
  out$metadata <- .dml_metadata("binary", dml_inference = out$scalar_inference)
  out
}

# ── Algorithm B: dml_plr — Partially Linear Regression ───────────────────────
# Chernozhukov et al. (2018) §3, Robinson (1988) partialling-out

.est_dml_plr <- function(DeltaY, G, X, K, S, alpha, learner, cluster = NULL,
                         dml_inference = "score", dml_boot = 500L,
                         dml_multiplier = "normal") {
  n <- length(DeltaY)
  beta_s <- numeric(S)
  var_s  <- numeric(S)
  influence_mat <- matrix(NA_real_, n, S)

  for (ss in seq_len(S)) {
    fold_id <- sample(rep(seq_len(K), length.out = n))
    V_tilde <- numeric(n)
    U_tilde <- numeric(n)

    for (k in seq_len(K)) {
      tr  <- which(fold_id != k)
      ho  <- which(fold_id == k)
      X_tr <- if (ncol(X) > 0) X[tr, , drop = FALSE] else matrix(0, length(tr), 0)
      X_ho <- if (ncol(X) > 0) X[ho, , drop = FALSE] else matrix(0, length(ho), 0)
      l_hat    <- .fit_reg(DeltaY[tr], X_tr, X_ho, learner)
      m_hat    <- .fit_reg(G[tr],      X_tr, X_ho, learner)
      V_tilde[ho] <- G[ho]       - m_hat
      U_tilde[ho] <- DeltaY[ho]  - l_hat
    }

    denom <- sum(V_tilde^2)
    if (denom < 1e-10) { beta_s[ss] <- 0; var_s[ss] <- Inf; next }
    beta_s[ss] <- sum(V_tilde * U_tilde) / denom

    psi_hat    <- V_tilde * (U_tilde - beta_s[ss] * V_tilde)
    J_hat      <- mean(V_tilde^2)
    var_s[ss]  <- .dml_score_variance(psi_hat, center = 0,
                                      jacobian = J_hat, cluster = cluster)
    influence_mat[, ss] <- psi_hat / J_hat
  }

  b_f  <- median(beta_s)
  se_f <- sqrt(max(0, median(pmax(0, var_s) + (beta_s - b_f)^2)))
  result <- c(Estimate  = b_f,
              Std.Error = se_f,
              CI_Lower  = b_f - stats::qnorm(1 - alpha/2) * se_f,
              CI_Upper  = b_f + stats::qnorm(1 - alpha/2) * se_f)
  out <- .dml_attach_scalar_inference(
    result = result,
    estimate = b_f,
    se = se_f,
    alpha = alpha,
    influence = rowMeans(influence_mat, na.rm = TRUE),
    split_estimates = beta_s,
    dml_inference = dml_inference,
    dml_boot = dml_boot,
    dml_multiplier = dml_multiplier,
    cluster = cluster
  )
  out$metadata <- .dml_metadata("plr", dml_inference = out$scalar_inference)
  out
}

# ── Algorithm C: dml_flex — Semenova & Chernozhukov (2021) dose-response ──────

.local_poly_fit2 <- function(g0, y, x, bw, degree = 2L) {
  .local_poly_eval2(g0, y, x, bw, degree)[["theta"]]
}

.local_poly_eval2 <- function(g0, y, x, bw, degree = 2L) {
  u  <- (x - g0) / bw
  W  <- dnorm(u)
  ok <- W > 1e-12
  if (sum(ok) < degree + 2L) {
    return(c(theta = NA_real_, delta = NA_real_,
             se_theta = NA_real_, se_delta = NA_real_))
  }
  dx <- x[ok] - g0
  P  <- outer(dx, 0L:degree, "^")
  WP <- P * W[ok]
  A  <- crossprod(WP, P)
  bv <- crossprod(WP, y[ok])
  A_inv <- tryCatch(solve(A), error = function(e) matrix(NA_real_, degree + 1L, degree + 1L))
  if (anyNA(A_inv)) {
    return(c(theta = NA_real_, delta = NA_real_,
             se_theta = NA_real_, se_delta = NA_real_))
  }
  th <- as.numeric(A_inv %*% bv)
  e_hat <- as.numeric(y[ok] - P %*% th)
  meat <- crossprod(P * (W[ok]^2 * e_hat^2), P)
  V <- A_inv %*% meat %*% A_inv
  delta <- if (degree >= 1L) th[2L] else NA_real_
  se_delta <- if (degree >= 1L) sqrt(pmax(0, V[2L, 2L])) else NA_real_
  c(theta = th[1L],
    delta = delta,
    se_theta = sqrt(pmax(0, V[1L, 1L])),
    se_delta = se_delta)
}

.select_bw_loocv2 <- function(y, x, n_grid = 15L, degree = 2L) {
  n_x     <- length(x)
  sg      <- stats::sd(x); if (sg < 1e-10) sg <- 1
  h_min   <- max(0.05 * sg, (degree + 2L) / n_x)
  h_max   <- 2.0 * sg
  if (h_min >= h_max) h_max <- 4.0 * h_min
  bw_grid <- exp(seq(log(h_min), log(h_max), length.out = n_grid))
  loo_err <- vapply(bw_grid, function(bw) {
    errs <- vapply(seq_along(y), function(i) {
      yhat <- .local_poly_fit2(x[i], y[-i], x[-i], bw, degree = degree)
      if (is.na(yhat)) (y[i] - mean(y[-i]))^2 else (y[i] - yhat)^2
    }, numeric(1L))
    mean(errs, na.rm = TRUE)
  }, numeric(1L))
  bw_grid[which.min(loo_err)]
}

.normalize_signal_map <- function(signal_map) {
  if (is.null(signal_map)) return("local_poly")
  if (!is.character(signal_map) || length(signal_map) != 1L || is.na(signal_map)) {
    stop("signal_map must be one of 'local_poly', 'spline', 'blp_spline', 'gam', or 'kernel'.",
         call. = FALSE)
  }
  sm <- gsub("[[:space:]_-]+", "", tolower(signal_map))
  if (sm %in% c("localpoly", "localpolynomial", "localquad", "localquadratic")) {
    return("local_poly")
  }
  if (sm %in% c("spline", "bspline", "bs")) return("spline")
  if (sm %in% c("blpspline", "blp", "sievecovariance", "sieve")) return("blp_spline")
  if (sm %in% c("gam", "mgcv")) return("gam")
  if (sm %in% c("kernel", "locallinear", "nadarayawatson", "nw")) return("kernel")
  stop("Unknown signal_map: '", signal_map,
       "'. Choose 'local_poly', 'spline', 'blp_spline', 'gam', or 'kernel'.",
       call. = FALSE)
}

.normalize_density_method <- function(density_method) {
  if (is.null(density_method)) return("residual_kde")
  if (!is.character(density_method) || length(density_method) != 1L ||
      is.na(density_method)) {
    stop("density_method must be one of 'residual_kde', 'location_scale', or 'local_kde'.",
         call. = FALSE)
  }
  dm <- gsub("[[:space:]_-]+", "", tolower(density_method))
  if (dm %in% c("residualkde", "residual", "kde")) return("residual_kde")
  if (dm %in% c("locationscale", "locscale", "heteroskedastic")) return("location_scale")
  if (dm %in% c("localkde", "localdensity", "conditionalkde")) return("local_kde")
  stop("Unknown density_method: '", density_method,
       "'. Choose 'residual_kde', 'location_scale', or 'local_kde'.",
       call. = FALSE)
}

.normalize_dml_inference <- function(dml_inference) {
  if (is.null(dml_inference)) return("score")
  if (!is.character(dml_inference) || length(dml_inference) != 1L ||
      is.na(dml_inference)) {
    stop("dml_inference must be 'score' or 'score_multiplier'.",
         call. = FALSE)
  }
  di <- gsub("[[:space:]_-]+", "", tolower(dml_inference))
  if (di %in% c("score", "analytical", "analytic", "robust")) return("score")
  if (di %in% c("scoremultiplier", "multiplier", "multiplierscore")) {
    return("score_multiplier")
  }
  stop("Unknown dml_inference: '", dml_inference,
       "'. Choose 'score' or 'score_multiplier'.", call. = FALSE)
}

.normalize_dml_multiplier <- function(dml_multiplier) {
  if (is.null(dml_multiplier)) return("normal")
  if (!is.character(dml_multiplier) || length(dml_multiplier) != 1L ||
      is.na(dml_multiplier)) {
    stop("dml_multiplier must be 'normal', 'wild', or 'exponential'.",
         call. = FALSE)
  }
  dm <- gsub("[[:space:]_-]+", "", tolower(dml_multiplier))
  if (dm %in% c("normal", "gaussian")) return("normal")
  if (dm %in% c("wild", "rademacher")) return("wild")
  if (dm %in% c("exponential", "bayes", "bayesian")) return("exponential")
  stop("Unknown dml_multiplier: '", dml_multiplier,
       "'. Choose 'normal', 'wild', or 'exponential'.", call. = FALSE)
}

.dml_metadata <- function(method, dml_inference = "score",
                          uniform_band_method = NA_character_) {
  method <- sub("^dml_", "", method)
  scalar_inference_type <- if (identical(dml_inference, "score_multiplier")) {
    "score_multiplier"
  } else if (identical(dml_inference, "second_stage_mapping")) {
    "second_stage_mapping"
  } else {
    "analytical_score"
  }
  base <- list(
    target_family = "continuous_G_DML",
    target_population = "all_complete_observations",
    inference_method = dml_inference,
    uniform_band_method = uniform_band_method
  )
  if (method == "binary") {
    return(c(base, list(
      target_estimand = "binary_AIPW_IRM_contrast",
      target_scale = "scalar",
      inference_type = scalar_inference_type,
      inference_scope = "pointwise_scalar",
      orthogonal_score = "binary_AIPW_IRM"
    )))
  }
  if (method == "plr") {
    return(c(base, list(
      target_estimand = "partially_linear_slope",
      target_scale = "scalar",
      inference_type = scalar_inference_type,
      inference_scope = "pointwise_scalar",
      orthogonal_score = "partialling_out"
    )))
  }
  if (method == "flex") {
    return(c(base, list(
      target_estimand = "fixed_G_level_curve",
      target_scale = "curve",
      inference_type = scalar_inference_type,
      inference_scope = "pointwise_curve_and_uniform_band",
      orthogonal_score = "continuous_treatment_level_curve_signal"
    )))
  }
  if (method == "incremental") {
    return(c(base, list(
      target_estimand = "observed_population_average_derivative",
      target_scale = "scalar",
      inference_type = scalar_inference_type,
      inference_scope = "pointwise_scalar",
      orthogonal_score = "average_derivative_density_score"
    )))
  }
  base
}

.dml_multiplier_weights <- function(n, B, method) {
  method <- .normalize_dml_multiplier(method)
  if (method == "normal") {
    return(matrix(stats::rnorm(n * B), nrow = n, ncol = B))
  }
  if (method == "wild") {
    return(matrix(sample(c(-1, 1), n * B, replace = TRUE), nrow = n, ncol = B))
  }
  matrix(stats::rexp(n * B, rate = 1) - 1, nrow = n, ncol = B)
}

.dml_scalar_multiplier <- function(influence, se, alpha, cluster = NULL,
                                   B = 500L, method = "normal") {
  influence <- as.numeric(influence)
  ok <- is.finite(influence)
  if (!any(ok)) return(NULL)
  influence <- influence[ok]
  n <- length(influence)
  B <- as.integer(B)
  if (!is.finite(B) || B < 1L) return(NULL)
  method <- .normalize_dml_multiplier(method)
  if (!is.finite(se) || se <= 0) return(NULL)

  if (!is.null(cluster)) {
    cl <- cluster[ok]
    if (length(cl) == n && !all(is.na(cl))) {
      cl <- as.factor(cl)
      infl_use <- rowsum(influence, cl, reorder = FALSE)[, 1]
      denom_n <- n
    } else {
      infl_use <- influence
      denom_n <- n
    }
  } else {
    infl_use <- influence
    denom_n <- n
  }
  scale <- denom_n * se
  if (!is.finite(scale) || scale <= 0) return(NULL)
  W <- .dml_multiplier_weights(length(infl_use), B, method)
  t_vals <- as.numeric(crossprod(infl_use, W) / scale)
  t_vals <- t_vals[is.finite(t_vals)]
  if (length(t_vals) < 1L) return(NULL)
  crit <- as.numeric(stats::quantile(abs(t_vals), probs = 1 - alpha,
                                     names = FALSE, na.rm = TRUE))
  if (!is.finite(crit)) return(NULL)
  list(
    t = t_vals,
    crit = crit,
    method = paste0(method, " score multiplier"),
    B = length(t_vals)
  )
}

.dml_attach_scalar_inference <- function(result, estimate, se, alpha,
                                         influence, split_estimates,
                                         dml_inference, dml_boot,
                                         dml_multiplier, cluster = NULL) {
  dml_inference <- .normalize_dml_inference(dml_inference)
  out <- list(
    result = result,
    influence = influence,
    split_estimates = split_estimates,
    scalar_inference = "analytical_score",
    scalar_multiplier_reps = NULL,
    scalar_multiplier_crit = NA_real_,
    scalar_multiplier_method = NA_character_
  )
  mult <- .dml_scalar_multiplier(
    influence = influence,
    se = se,
    alpha = alpha,
    cluster = cluster,
    B = dml_boot,
    method = dml_multiplier
  )
  if (!is.null(mult)) {
    out$scalar_multiplier_reps <- mult$t
    out$scalar_multiplier_crit <- mult$crit
    out$scalar_multiplier_method <- mult$method
  }
  if (dml_inference == "score_multiplier" && !is.null(mult)) {
    result["CI_Lower"] <- estimate - mult$crit * se
    result["CI_Upper"] <- estimate + mult$crit * se
    out$result <- result
    out$scalar_inference <- "score_multiplier"
  }
  out
}

.dml_gaussian_max_t_crit <- function(V, alpha = 0.05, B = 1000L) {
  if (is.null(V) || !is.matrix(V) || nrow(V) == 0L || nrow(V) != ncol(V)) {
    return(NA_real_)
  }
  V <- (V + t(V)) / 2
  se <- sqrt(pmax(diag(V), 0))
  ok <- is.finite(se) & se > 1e-10
  if (!any(ok)) return(NA_real_)
  V_ok <- V[ok, ok, drop = FALSE]
  eg <- eigen(V_ok, symmetric = TRUE)
  vals <- pmax(eg$values, 0)
  root <- eg$vectors %*% diag(sqrt(vals), nrow = length(vals), ncol = length(vals))
  B <- as.integer(B)
  if (!is.finite(B) || B < 1L) B <- 1000L
  Z <- matrix(stats::rnorm(B * nrow(V_ok)), nrow = B) %*% t(root)
  max_t_stat <- apply(abs(sweep(Z, 2, se[ok], "/")), 1L, max,
                      na.rm = TRUE)
  crit <- as.numeric(stats::quantile(max_t_stat[is.finite(max_t_stat)],
                                     probs = 1 - alpha, names = FALSE,
                                     na.rm = TRUE))
  if (!is.finite(crit)) NA_real_ else crit
}

.dml_average_vcov <- function(vcov_list, dim_expected) {
  valid <- Filter(function(V) {
    is.matrix(V) &&
      nrow(V) == dim_expected &&
      ncol(V) == dim_expected &&
      all(is.finite(V))
  }, vcov_list)
  if (length(valid) == 0L) return(NULL)
  V <- Reduce("+", valid) / length(valid)
  (V + t(V)) / 2
}

.dml_rescale_vcov_diag <- function(V, se_target) {
  if (is.null(V) || !is.matrix(V)) return(NULL)
  V <- (V + t(V)) / 2
  se_old <- sqrt(pmax(diag(V), 0))
  se_target <- as.numeric(se_target)
  ok <- is.finite(se_old) & se_old > 1e-10 &
    is.finite(se_target) & se_target >= 0
  scale <- rep(1, length(se_old))
  scale[ok] <- se_target[ok] / se_old[ok]
  V <- sweep(sweep(V, 1L, scale, "*"), 2L, scale, "*")
  diag_ok <- is.finite(se_target) & se_target >= 0
  diag(V)[diag_ok] <- se_target[diag_ok]^2
  (V + t(V)) / 2
}

.safe_kde_bw <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  sx <- stats::sd(x)
  if (!is.finite(sx) || sx < 1e-8) sx <- 1
  bw <- tryCatch(stats::bw.SJ(x),
                 error = function(e) 1.06 * sx * n^(-1/5))
  if (!is.finite(bw) || bw <= 0) bw <- 1.06 * sx * n^(-1/5)
  pmax(bw, 1e-6)
}

.fit_cond_density <- function(G_tr, X_tr, learner, density_method) {
  density_method <- .normalize_density_method(density_method)
  n_tr <- length(G_tr)
  density_floor <- 0.005

  if (density_method == "residual_kde") {
    m_model  <- .fit_reg_model(G_tr, X_tr, learner)
    m_hat_tr <- .predict_reg_model(m_model, X_tr)
    V_tr     <- G_tr - m_hat_tr
    kde_V    <- stats::density(V_tr, bw = .safe_kde_bw(V_tr), n = 1024)
    return(list(method = "residual_kde",
                m_model = m_model,
                kde_fn = stats::approxfun(kde_V, rule = 2),
                floor = density_floor))
  }

  if (density_method == "location_scale") {
    m_model  <- .fit_reg_model(G_tr, X_tr, learner)
    m_hat_tr <- .predict_reg_model(m_model, X_tr)
    V_tr     <- G_tr - m_hat_tr
    sigma_floor <- max(0.05 * stats::sd(G_tr), 1e-4)
    logv_tr <- log(V_tr^2 + sigma_floor^2)
    sigma_model <- .fit_reg_model(logv_tr, X_tr, learner)
    sigma_hat_tr <- sqrt(pmax(exp(.predict_reg_model(sigma_model, X_tr)),
                              sigma_floor^2))
    Z_tr <- V_tr / sigma_hat_tr
    kde_Z <- stats::density(Z_tr, bw = .safe_kde_bw(Z_tr), n = 1024)
    return(list(method = "location_scale",
                m_model = m_model,
                sigma_model = sigma_model,
                kde_fn = stats::approxfun(kde_Z, rule = 2),
                sigma_floor = sigma_floor,
                floor = density_floor))
  }

  # local_kde: product Gaussian kernel over scaled X and Gaussian KDE over G.
  bw_g <- .safe_kde_bw(G_tr)
  if (ncol(X_tr) == 0L) {
    kde_G <- stats::density(G_tr, bw = bw_g, n = 1024)
    return(list(method = "local_kde",
                marginal = TRUE,
                kde_fn = stats::approxfun(kde_G, rule = 2),
                floor = density_floor))
  }

  sc <- .scale_train(as.matrix(X_tr))
  p_x <- ncol(X_tr)
  bw_x <- max(n_tr^(-1 / (4 + p_x)), 0.25)
  list(method = "local_kde",
       marginal = FALSE,
       G_tr = G_tr,
       X_scaled = sc$X,
       center = sc$center,
       scale = sc$scale,
       bw_g = bw_g,
       bw_x = bw_x,
       floor = density_floor)
}

.predict_cond_density <- function(model, g_new, X_new) {
  g_new <- as.numeric(g_new)
  X_new <- as.matrix(X_new)
  n_new <- length(g_new)

  if (model$method == "residual_kde") {
    m_hat <- .predict_reg_model(model$m_model, X_new)
    out <- as.numeric(model$kde_fn(g_new - m_hat))
    return(pmax(out, model$floor))
  }

  if (model$method == "location_scale") {
    m_hat <- .predict_reg_model(model$m_model, X_new)
    logv_hat <- .predict_reg_model(model$sigma_model, X_new)
    sigma_hat <- sqrt(pmax(exp(logv_hat), model$sigma_floor^2))
    out <- as.numeric(model$kde_fn((g_new - m_hat) / sigma_hat)) / sigma_hat
    return(pmax(out, model$floor))
  }

  if (model$method == "local_kde" && isTRUE(model$marginal)) {
    return(pmax(as.numeric(model$kde_fn(g_new)), model$floor))
  }

  if (model$method == "local_kde") {
    X_scaled_new <- .scale_new(X_new, model$center, model$scale)
    out <- vapply(seq_len(n_new), function(i) {
      dx <- sweep(model$X_scaled, 2, X_scaled_new[i, ], "-") / model$bw_x
      logw <- -0.5 * rowSums(dx^2)
      logw <- logw - max(logw)
      w <- exp(logw)
      sw <- sum(w)
      if (!is.finite(sw) || sw <= .Machine$double.eps) {
        w <- rep(1, length(model$G_tr))
        sw <- length(model$G_tr)
      }
      kg <- stats::dnorm((g_new[i] - model$G_tr) / model$bw_g) / model$bw_g
      sum(w * kg) / sw
    }, numeric(1L))
    return(pmax(out, model$floor))
  }

  stop("Internal error: unknown conditional density model.", call. = FALSE)
}

.finite_diff_points <- function(eval_g, x) {
  rg <- range(x, eval_g, finite = TRUE)
  span <- diff(rg)
  if (!is.finite(span) || span <= 0) span <- 1
  eps <- max(1e-5 * span, .Machine$double.eps^0.25)
  x_plus <- pmin(eval_g + eps, rg[2L])
  x_minus <- pmax(eval_g - eps, rg[1L])
  bad <- x_plus <= x_minus
  if (any(bad)) {
    x_plus[bad] <- eval_g[bad] + eps
    x_minus[bad] <- eval_g[bad] - eps
  }
  list(plus = x_plus, minus = x_minus, denom = x_plus - x_minus)
}

.map_signal_lm <- function(fit, eval_g, x) {
  nd <- data.frame(g = eval_g)
  pred <- stats::predict(fit, newdata = nd, se.fit = TRUE)
  fd <- .finite_diff_points(eval_g, x)
  Xp <- stats::model.matrix(stats::delete.response(stats::terms(fit)),
                            data.frame(g = fd$plus))
  Xm <- stats::model.matrix(stats::delete.response(stats::terms(fit)),
                            data.frame(g = fd$minus))
  D <- sweep(Xp - Xm, 1, fd$denom, "/")
  co <- stats::coef(fit)
  V <- stats::vcov(fit)
  delta <- as.numeric(D %*% co)
  se_delta <- sqrt(pmax(0, rowSums((D %*% V) * D)))
  list(theta = as.numeric(pred$fit),
       delta = delta,
       se_theta = as.numeric(pred$se.fit),
       se_delta = se_delta)
}

.map_signal_gam <- function(fit, eval_g, x) {
  nd <- data.frame(g = eval_g)
  pred <- stats::predict(fit, newdata = nd, se.fit = TRUE, type = "response")
  fd <- .finite_diff_points(eval_g, x)
  Xp <- stats::predict(fit, newdata = data.frame(g = fd$plus), type = "lpmatrix")
  Xm <- stats::predict(fit, newdata = data.frame(g = fd$minus), type = "lpmatrix")
  D <- sweep(Xp - Xm, 1, fd$denom, "/")
  co <- stats::coef(fit)
  V <- fit$Vp
  delta <- as.numeric(D %*% co)
  se_delta <- sqrt(pmax(0, rowSums((D %*% V) * D)))
  list(theta = as.numeric(pred$fit),
       delta = delta,
       se_theta = as.numeric(pred$se.fit),
       se_delta = se_delta)
}

.safe_solve <- function(A) {
  out <- tryCatch(solve(A), error = function(e) NULL)
  if (!is.null(out) && all(is.finite(out))) return(out)
  eg <- eigen((A + t(A)) / 2, symmetric = TRUE)
  tol <- max(dim(A)) * max(abs(eg$values), na.rm = TRUE) * .Machine$double.eps
  vals_inv <- ifelse(abs(eg$values) > tol, 1 / eg$values, 0)
  eg$vectors %*% diag(vals_inv, nrow = length(vals_inv), ncol = length(vals_inv)) %*%
    t(eg$vectors)
}

.map_signal_blp_spline <- function(xi_hat, G, eval_g, map_degree, map_df) {
  n_unique <- length(unique(G))
  degree_use <- min(max(as.integer(map_degree), 1L), 3L)
  df_use <- if (is.null(map_df)) {
    min(max(degree_use + 2L, 5L), max(degree_use + 1L, n_unique - 1L))
  } else {
    min(max(as.integer(map_df), degree_use + 1L), max(degree_use + 1L, n_unique - 1L))
  }
  dat <- data.frame(y = xi_hat, g = G)
  boundary <- range(G, finite = TRUE)
  fml <- y ~ splines::bs(g, df = df_use, degree = degree_use,
                         Boundary.knots = boundary)
  fit <- stats::lm(fml, data = dat)
  terms_obj <- stats::delete.response(stats::terms(fit))
  B_train <- stats::model.matrix(fit)
  B_grid <- stats::model.matrix(terms_obj, data.frame(g = eval_g))
  beta <- stats::coef(fit)
  beta[!is.finite(beta)] <- 0
  resid <- as.numeric(stats::residuals(fit))
  XtX_inv <- .safe_solve(crossprod(B_train))
  meat <- crossprod(B_train * resid)
  V_beta <- XtX_inv %*% meat %*% XtX_inv
  V_beta <- (V_beta + t(V_beta)) / 2
  theta <- as.numeric(B_grid %*% beta)
  theta_vcov <- B_grid %*% V_beta %*% t(B_grid)
  theta_vcov <- (theta_vcov + t(theta_vcov)) / 2

  fd <- .finite_diff_points(eval_g, G)
  Bp <- stats::model.matrix(terms_obj, data.frame(g = fd$plus))
  Bm <- stats::model.matrix(terms_obj, data.frame(g = fd$minus))
  D_grid <- sweep(Bp - Bm, 1L, fd$denom, "/")
  delta <- as.numeric(D_grid %*% beta)
  delta_vcov <- D_grid %*% V_beta %*% t(D_grid)
  delta_vcov <- (delta_vcov + t(delta_vcov)) / 2
  level_derivative_vcov <- B_grid %*% V_beta %*% t(D_grid)

  list(
    theta = theta,
    delta = delta,
    se_theta = sqrt(pmax(0, diag(theta_vcov))),
    se_delta = sqrt(pmax(0, diag(delta_vcov))),
    theta_vcov = theta_vcov,
    delta_vcov = delta_vcov,
    level_derivative_vcov = level_derivative_vcov,
    map_info = list(df = df_use, degree = degree_use, type = "blp_spline")
  )
}

.map_signal_to_g <- function(xi_hat, G, eval_g, signal_map, map_degree, map_df) {
  signal_map <- .normalize_signal_map(signal_map)
  map_degree <- as.integer(map_degree)
  if (length(map_degree) != 1L || is.na(map_degree) ||
      map_degree < 1L || map_degree > 5L) {
    stop("map_degree must be a single integer between 1 and 5.", call. = FALSE)
  }
  if (!is.null(map_df)) {
    map_df <- as.integer(map_df)
    if (length(map_df) != 1L || is.na(map_df) || map_df < 3L) {
      stop("map_df must be NULL or a single integer >= 3.", call. = FALSE)
    }
  }

  n_unique <- length(unique(G))
  if (signal_map %in% c("local_poly", "kernel")) {
    degree_use <- if (signal_map == "kernel") 1L else map_degree
    bw_2nd <- tryCatch(
      .select_bw_loocv2(xi_hat, G, degree = degree_use),
      error = function(e) {
        sg <- stats::sd(G)
        if (!is.finite(sg) || sg < 1e-10) sg <- 1
        0.5 * sg
      }
    )
    lp <- vapply(eval_g, function(g0)
      .local_poly_eval2(g0, xi_hat, G, bw_2nd, degree = degree_use),
      numeric(4L))
    return(list(theta = lp["theta", ],
                delta = lp["delta", ],
                se_theta = lp["se_theta", ],
                se_delta = lp["se_delta", ],
                map_info = list(bandwidth = bw_2nd, degree = degree_use)))
  }

  if (signal_map == "spline") {
    degree_use <- min(max(map_degree, 1L), 3L)
    df_use <- if (is.null(map_df)) {
      min(max(degree_use + 2L, 5L), max(degree_use + 1L, n_unique - 1L))
    } else {
      min(max(map_df, degree_use + 1L), max(degree_use + 1L, n_unique - 1L))
    }
    dat <- data.frame(y = xi_hat, g = G)
    fit <- stats::lm(
      y ~ splines::bs(g, df = df_use, degree = degree_use,
                      Boundary.knots = range(G, finite = TRUE)),
      data = dat
    )
    mp <- .map_signal_lm(fit, eval_g, G)
    mp$map_info <- list(df = df_use, degree = degree_use)
    return(mp)
  }

  if (signal_map == "blp_spline") {
    return(.map_signal_blp_spline(xi_hat, G, eval_g, map_degree, map_df))
  }

  if (signal_map == "gam") {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop("signal_map='gam' requires the 'mgcv' package.", call. = FALSE)
    }
    k_use <- if (is.null(map_df)) min(10L, n_unique - 1L) else map_df
    k_use <- min(max(3L, as.integer(k_use)), max(3L, n_unique - 1L))
    dat <- data.frame(y = xi_hat, g = G)
    fit <- mgcv::gam(stats::as.formula(paste0("y ~ s(g, k = ", k_use, ")")),
                     data = dat, method = "REML")
    mp <- .map_signal_gam(fit, eval_g, G)
    mp$map_info <- list(k = k_use)
    return(mp)
  }

  stop("Internal error: unknown signal_map.", call. = FALSE)
}

.predict_mu_g_derivative <- function(mu_model, G_new, X_new, G_support) {
  X_new <- as.matrix(X_new)
  fd <- .finite_diff_points(G_new, G_support)
  feat_plus <- cbind(fd$plus, X_new)
  feat_minus <- cbind(fd$minus, X_new)
  pred_plus <- .predict_reg_model(mu_model, feat_plus)
  pred_minus <- .predict_reg_model(mu_model, feat_minus)
  (pred_plus - pred_minus) / fd$denom
}

.predict_log_density_g_derivative <- function(s_model, G_new, X_new, G_support) {
  X_new <- as.matrix(X_new)
  fd <- .finite_diff_points(G_new, G_support)
  s_plus <- .predict_cond_density(s_model, fd$plus, X_new)
  s_minus <- .predict_cond_density(s_model, fd$minus, X_new)
  (log(pmax(s_plus, 1e-8)) - log(pmax(s_minus, 1e-8))) / fd$denom
}

.est_dml_incremental_scalar <- function(DeltaY, G, X, K, S, alpha, learner,
                                        density_method = "residual_kde",
                                        cluster = NULL,
                                        dml_inference = "score",
                                        dml_boot = 500L,
                                        dml_multiplier = "normal") {
  density_method <- .normalize_density_method(density_method)
  n <- length(DeltaY)
  tau_s <- numeric(S)
  var_s <- numeric(S)
  phi_mat <- matrix(NA_real_, n, S)
  density_mat <- matrix(NA_real_, n, S)
  d_log_s_mat <- matrix(NA_real_, n, S)
  d_mu_mat <- matrix(NA_real_, n, S)
  correction_mat <- matrix(NA_real_, n, S)
  influence_mat <- matrix(NA_real_, n, S)

  for (ss in seq_len(S)) {
    fold_id <- sample(rep(seq_len(K), length.out = n))
    phi_hat <- numeric(n)
    s_hat_all <- numeric(n)
    d_log_s_all <- numeric(n)
    d_mu_all <- numeric(n)
    correction_all <- numeric(n)

    for (k in seq_len(K)) {
      tr <- which(fold_id != k)
      ho <- which(fold_id == k)
      X_tr <- if (ncol(X) > 0) X[tr, , drop = FALSE] else matrix(0, length(tr), 0)
      X_ho <- if (ncol(X) > 0) X[ho, , drop = FALSE] else matrix(0, length(ho), 0)

      feat_tr <- cbind(G[tr], X_tr)
      mu_model <- .fit_reg_model(DeltaY[tr], feat_tr, learner)
      s_model <- .fit_cond_density(G[tr], X_tr, learner, density_method)

      feat_ho <- cbind(G[ho], X_ho)
      mu_hat <- .predict_reg_model(mu_model, feat_ho)
      d_mu_hat <- .predict_mu_g_derivative(mu_model, G[ho], X_ho, G[tr])
      d_log_s_hat <- .predict_log_density_g_derivative(s_model, G[ho], X_ho, G[tr])

      s_hat_all[ho] <- .predict_cond_density(s_model, G[ho], X_ho)
      d_mu_all[ho] <- d_mu_hat
      d_log_s_all[ho] <- d_log_s_hat
      correction_all[ho] <- -d_log_s_hat * (DeltaY[ho] - mu_hat)
      phi_hat[ho] <- d_mu_hat + correction_all[ho]
    }

    tau_s[ss] <- mean(phi_hat, na.rm = TRUE)
    var_s[ss] <- .dml_score_variance(phi_hat, center = tau_s[ss], cluster = cluster)
    phi_mat[, ss] <- phi_hat
    influence_mat[, ss] <- phi_hat - tau_s[ss]
    density_mat[, ss] <- s_hat_all
    d_log_s_mat[, ss] <- d_log_s_all
    d_mu_mat[, ss] <- d_mu_all
    correction_mat[, ss] <- correction_all
  }

  tau_f <- stats::median(tau_s, na.rm = TRUE)
  se_f <- sqrt(max(0, stats::median(var_s + (tau_s - tau_f)^2, na.rm = TRUE)))
  phi_avg <- rowMeans(phi_mat, na.rm = TRUE)
  density_avg <- rowMeans(density_mat, na.rm = TRUE)
  d_log_s_avg <- rowMeans(d_log_s_mat, na.rm = TRUE)
  d_mu_avg <- rowMeans(d_mu_mat, na.rm = TRUE)
  correction_avg <- rowMeans(correction_mat, na.rm = TRUE)
  result <- c(
    Estimate  = tau_f,
    Std.Error = se_f,
    CI_Lower  = tau_f - stats::qnorm(1 - alpha/2) * se_f,
    CI_Upper  = tau_f + stats::qnorm(1 - alpha/2) * se_f
  )
  out <- .dml_attach_scalar_inference(
    result = result,
    estimate = tau_f,
    se = se_f,
    alpha = alpha,
    influence = rowMeans(influence_mat, na.rm = TRUE),
    split_estimates = tau_s,
    dml_inference = dml_inference,
    dml_boot = dml_boot,
    dml_multiplier = dml_multiplier,
    cluster = cluster
  )
  out$signal <- phi_avg
  out$density_method <- density_method
  out$diagnostics <- .dml_diagnostics(
    signal = phi_avg,
    density = density_avg,
    d_log_density = d_log_s_avg,
    d_mu = d_mu_avg,
    correction = correction_avg,
    cluster = cluster,
    density_method = density_method,
    learner = learner
  )
  out$metadata <- .dml_metadata("incremental", dml_inference = out$scalar_inference)
  out
}

.est_dml_flex <- function(DeltaY, G, X, eval_g, K, S, alpha, learner,
                          signal_map = "local_poly", map_degree = 2L,
                          map_df = NULL, density_method = "residual_kde",
                          boot = 100L) {
  signal_map <- .normalize_signal_map(signal_map)
  density_method <- .normalize_density_method(density_method)
  n  <- length(DeltaY)
  m  <- length(eval_g)
  all_theta <- matrix(NA_real_, S, m)
  all_delta <- matrix(NA_real_, S, m)
  all_se_theta <- matrix(NA_real_, S, m)
  all_se_delta <- matrix(NA_real_, S, m)
  xi_mat <- matrix(NA_real_, n, S)
  density_mat <- matrix(NA_real_, n, S)
  map_info <- vector("list", S)
  theta_vcov_list <- vector("list", S)
  delta_vcov_list <- vector("list", S)
  level_derivative_vcov_list <- vector("list", S)

  for (ss in seq_len(S)) {
    fold_id <- sample(rep(seq_len(K), length.out = n))
    xi_hat  <- numeric(n)

    # Build per-fold nuisance models
    mu_models   <- vector("list", K)
    s_models    <- vector("list", K)
    omega_fns   <- vector("list", K)

    for (k in seq_len(K)) {
      tr   <- which(fold_id != k)
      X_tr <- if (ncol(X) > 0) X[tr, , drop = FALSE] else matrix(0, length(tr), 0)

      # 2a: outcome regression μ(g,x) = E[ΔY | G=g, X=x], treat G as feature
      feat_tr <- cbind(G[tr], X_tr)
      mu_models[[k]] <- .fit_reg_model(DeltaY[tr], feat_tr, learner)
      mu_models[[k]]$p_x <- ncol(X_tr)

      # 2b: conditional density s(g|x)
      s_models[[k]] <- .fit_cond_density(G[tr], X_tr, learner, density_method)

      # 2c: marginal density ω̂(g)
      bw_g   <- .safe_kde_bw(G[tr])
      kde_G  <- stats::density(G[tr], bw = bw_g, n = 1024)
      omega_fns[[k]] <- stats::approxfun(kde_G, rule = 2)
    }

    # Step 3: pseudo-outcomes
    for (k in seq_len(K)) {
      ho       <- which(fold_id == k)
      n_ho     <- length(ho)
      mu_info  <- mu_models[[k]]
      s_info   <- s_models[[k]]
      omega_fn <- omega_fns[[k]]

      X_ho <- if (ncol(X) > 0) X[ho, , drop = FALSE] else matrix(0, n_ho, 0)
      s_pred_ho <- .predict_cond_density(s_info, G[ho], X_ho)
      density_mat[ho, ss] <- s_pred_ho

      # Predict μ̂(G[i], X[i]) for held-out observations (using stored model/coefs)
      feat_ho <- cbind(G[ho], X_ho)
      mu_pred_ho <- .predict_reg_model(mu_info, feat_ho)

      # LOO within-fold mean: E_{X}[μ̂(G_i, X)] for i in held-out
      if (n_ho <= 1L) {
        m_loo_vec <- mu_pred_ho
      } else if (mu_info$type == "lm") {
        # Exact LOO using linear structure: Σ_j μ̂(G_i, X_j) = intercept*n_ho + β_G*G_i*n_ho + β_X' ΣX_j
        co <- stats::coef(mu_info$model)
        b0 <- if (!is.na(co["(Intercept)"])) co["(Intercept)"] else 0
        bG <- if (!is.na(co["v1"])) co["v1"] else 0
        p_x <- if (is.null(mu_info$p_x)) 0L else mu_info$p_x
        bX_names <- if (p_x > 0L) paste0("v", seq.int(2L, p_x + 1L)) else character(0)
        bX <- if (p_x > 0L) co[bX_names] else numeric(0)
        bX[is.na(bX)] <- 0
        X_sum_ho <- if (ncol(X_ho) > 0) colSums(X_ho) else numeric(0)
        mu_sum_over_ho <- n_ho * b0 + n_ho * bG * G[ho] +
                          if (length(bX) > 0) as.numeric(bX %*% X_sum_ho) else 0
        # Subtract own prediction and divide by (n_ho - 1)
        mu_own <- b0 + bG * G[ho] +
                  if (ncol(X_ho) > 0) as.numeric(X_ho %*% bX) else 0
        m_loo_vec <- (mu_sum_over_ho - mu_own) / (n_ho - 1L)
      } else {
        # Exact within-fold LOO average for non-linear learners.
        m_loo_vec <- vapply(seq_len(n_ho), function(ii) {
          g_i <- G[ho[ii]]
          mate_idx <- setdiff(seq_len(n_ho), ii)
          if (length(mate_idx) == 0L) return(mu_pred_ho[ii])
          X_mates <- if (ncol(X_ho) > 0L) {
            X_ho[mate_idx, , drop = FALSE]
          } else {
            matrix(0, length(mate_idx), 0)
          }
          feat_i <- cbind(rep(g_i, length(mate_idx)), X_mates)
          mean(.predict_reg_model(mu_info, feat_i), na.rm = TRUE)
        }, numeric(1L))
      }

      # Assemble pseudo-outcome
      for (ii in seq_len(n_ho)) {
        i        <- ho[ii]
        g_i      <- G[i]
        s_i      <- s_pred_ho[ii]
        omega_i  <- pmax(omega_fn(g_i), 1e-8)
        ipw_term <- (DeltaY[i] - mu_pred_ho[ii]) / s_i * omega_i
        xi_hat[i] <- ipw_term + m_loo_vec[ii]
      }
    }

    # Step 4: second-stage mapping from orthogonal signals back to G
    mp <- .map_signal_to_g(xi_hat, G, eval_g, signal_map,
                           map_degree = map_degree, map_df = map_df)
    xi_mat[, ss] <- xi_hat
    all_theta[ss, ]    <- mp$theta
    all_delta[ss, ]    <- mp$delta
    all_se_theta[ss, ] <- mp$se_theta
    all_se_delta[ss, ] <- mp$se_delta
    map_info[[ss]] <- mp$map_info
    theta_vcov_list[ss] <- list(mp$theta_vcov)
    delta_vcov_list[ss] <- list(mp$delta_vcov)
    level_derivative_vcov_list[ss] <- list(mp$level_derivative_vcov)
  }

  # Median aggregation across splits
  theta_final <- apply(all_theta, 2, stats::median, na.rm = TRUE)
  delta_final <- apply(all_delta, 2, stats::median, na.rm = TRUE)
  theta_split_dev <- sweep(all_theta, 2, theta_final, "-")
  delta_split_dev <- sweep(all_delta, 2, delta_final, "-")
  se_theta <- sqrt(apply(all_se_theta^2 + theta_split_dev^2, 2, stats::median, na.rm = TRUE))
  se_delta <- sqrt(apply(all_se_delta^2 + delta_split_dev^2, 2, stats::median, na.rm = TRUE))
  theta_vcov <- .dml_rescale_vcov_diag(.dml_average_vcov(theta_vcov_list, m),
                                       se_theta)
  delta_vcov <- .dml_rescale_vcov_diag(.dml_average_vcov(delta_vcov_list, m),
                                       se_delta)
  level_derivative_vcov <- .dml_average_vcov(level_derivative_vcov_list, m)

  g_range <- diff(range(eval_g))
  interval_avg <- if (is.finite(g_range) && g_range > 0) {
    (theta_final[m] - theta_final[1L]) / g_range
  } else {
    delta_final[1L]
  }
  interval_se <- if (is.finite(g_range) && g_range > 0) {
    sqrt(se_theta[1L]^2 + se_theta[m]^2) / g_range
  } else {
    se_delta[1L]
  }
  mean_delta <- mean(delta_final, na.rm = TRUE)
  n_delta_se <- sum(is.finite(se_delta))
  mean_delta_se <- if (n_delta_se > 0L) {
    sqrt(mean(se_delta^2, na.rm = TRUE) / n_delta_se)
  } else {
    NA_real_
  }

  # Multiplier bands for the second-stage signal mapping. This is a practical
  # simultaneous-band diagnostic for the curve object; it is still reported as
  # pointwise-mapping based inference in the docs until the theory is finalized.
  boot <- as.integer(boot)
  boot <- if (is.finite(boot)) max(1L, boot) else 100L
  max_t_theta <- rep(NA_real_, boot)
  max_t_delta <- rep(NA_real_, boot)
  theta_multiplier_reps <- matrix(NA_real_, boot, m)
  delta_multiplier_reps <- matrix(NA_real_, boot, m)
  se_theta_safe <- pmax(se_theta, 1e-8)
  se_delta_safe <- pmax(se_delta, 1e-8)
  theta_obs_list <- vector("list", S)
  resid_list <- vector("list", S)
  for (ss in seq_len(S)) {
    obs_map <- tryCatch(
      .map_signal_to_g(xi_mat[, ss], G, G, signal_map,
                       map_degree = map_degree, map_df = map_df),
      error = function(e) NULL
    )
    if (!is.null(obs_map)) {
      theta_obs_list[[ss]] <- obs_map$theta
      resid_list[[ss]] <- xi_mat[, ss] - obs_map$theta
    }
  }
  for (bb in seq_len(boot)) {
    theta_b <- matrix(NA_real_, S, m)
    delta_b <- matrix(NA_real_, S, m)
    for (ss in seq_len(S)) {
      if (is.null(theta_obs_list[[ss]]) || is.null(resid_list[[ss]])) next
      xi_star <- theta_obs_list[[ss]] + resid_list[[ss]] * stats::rnorm(n)
      mp_b <- tryCatch(
        .map_signal_to_g(xi_star, G, eval_g, signal_map,
                         map_degree = map_degree, map_df = map_df),
        error = function(e) NULL
      )
      if (!is.null(mp_b)) {
        theta_b[ss, ] <- mp_b$theta
        delta_b[ss, ] <- mp_b$delta
      }
    }
    theta_b_final <- apply(theta_b, 2, stats::median, na.rm = TRUE)
    delta_b_final <- apply(delta_b, 2, stats::median, na.rm = TRUE)
    theta_multiplier_reps[bb, ] <- theta_b_final
    delta_multiplier_reps[bb, ] <- delta_b_final
    max_t_theta[bb] <- max(abs(theta_b_final - theta_final) / se_theta_safe,
                           na.rm = TRUE)
    max_t_delta[bb] <- max(abs(delta_b_final - delta_final) / se_delta_safe,
                           na.rm = TRUE)
  }
  colnames(theta_multiplier_reps) <- colnames(delta_multiplier_reps) <-
    paste0("g", seq_along(eval_g))
  crit_theta <- stats::quantile(max_t_theta[is.finite(max_t_theta)],
                                probs = 1 - alpha, names = FALSE, na.rm = TRUE)
  crit_delta <- stats::quantile(max_t_delta[is.finite(max_t_delta)],
                                probs = 1 - alpha, names = FALSE, na.rm = TRUE)
  practical_crit_theta <- crit_theta
  practical_crit_delta <- crit_delta
  if (!is.finite(practical_crit_theta)) practical_crit_theta <- stats::qnorm(1 - alpha / 2)
  if (!is.finite(practical_crit_delta)) practical_crit_delta <- stats::qnorm(1 - alpha / 2)

  gaussian_crit_theta <- .dml_gaussian_max_t_crit(theta_vcov, alpha = alpha, B = boot)
  gaussian_crit_delta <- .dml_gaussian_max_t_crit(delta_vcov, alpha = alpha, B = boot)
  crit_theta <- if (is.finite(gaussian_crit_theta)) gaussian_crit_theta else practical_crit_theta
  crit_delta <- if (is.finite(gaussian_crit_delta)) gaussian_crit_delta else practical_crit_delta
  uniform_band_method <- if (is.finite(gaussian_crit_theta) || is.finite(gaussian_crit_delta)) {
    "gaussian_max_t_blp_spline_covariance"
  } else {
    "practical_signal_residual_multiplier"
  }
  density_avg <- rowMeans(density_mat, na.rm = TRUE)

  list(
    theta_hat = theta_final,
    delta_hat = delta_final,
    se_theta  = se_theta,
    se_delta  = se_delta,
    theta_ci_lower = theta_final - stats::qnorm(1 - alpha/2) * se_theta,
    theta_ci_upper = theta_final + stats::qnorm(1 - alpha/2) * se_theta,
    delta_ci_lower = delta_final - stats::qnorm(1 - alpha/2) * se_delta,
    delta_ci_upper = delta_final + stats::qnorm(1 - alpha/2) * se_delta,
    theta_band_lower = theta_final - crit_theta * se_theta,
    theta_band_upper = theta_final + crit_theta * se_theta,
    delta_band_lower = delta_final - crit_delta * se_delta,
    delta_band_upper = delta_final + crit_delta * se_delta,
    theta_multiplier_reps = theta_multiplier_reps,
    delta_multiplier_reps = delta_multiplier_reps,
    simultaneous_crit_theta = crit_theta,
    simultaneous_crit_delta = crit_delta,
    gaussian_crit_theta = gaussian_crit_theta,
    gaussian_crit_delta = gaussian_crit_delta,
    practical_simultaneous_crit_theta = practical_crit_theta,
    practical_simultaneous_crit_delta = practical_crit_delta,
    simultaneous_band_method = uniform_band_method,
    practical_band_method = "practical residual multiplier over cross-fitted signal-mapping residuals",
    theta_vcov = theta_vcov,
    delta_vcov = delta_vcov,
    level_derivative_vcov = level_derivative_vcov,
    ci_lower  = theta_final - stats::qnorm(1 - alpha/2) * se_theta,
    ci_upper  = theta_final + stats::qnorm(1 - alpha/2) * se_theta,
    interval_average = c(
      Estimate  = interval_avg,
      Std.Error = interval_se,
      CI_Lower  = interval_avg - stats::qnorm(1 - alpha/2) * interval_se,
      CI_Upper  = interval_avg + stats::qnorm(1 - alpha/2) * interval_se
    ),
    mean_delta = c(
      Estimate  = mean_delta,
      Std.Error = mean_delta_se,
      CI_Lower  = mean_delta - stats::qnorm(1 - alpha/2) * mean_delta_se,
      CI_Upper  = mean_delta + stats::qnorm(1 - alpha/2) * mean_delta_se
    ),
    signal_map = signal_map,
    map_degree = map_degree,
    map_df = map_df,
    density_method = density_method,
    map_info = map_info,
    metadata = .dml_metadata("flex", dml_inference = "second_stage_mapping",
                             uniform_band_method = uniform_band_method),
    diagnostics = .dml_diagnostics(
      signal = rowMeans(xi_mat, na.rm = TRUE),
      density = density_avg,
      cluster = NULL,
      density_method = density_method,
      learner = learner
    )
  )
}

# ── Main DML dispatcher ───────────────────────────────────────────────────────

run_dml_method <- function(
  s,
  tr_period,
  ref_period,
  entire_period = NULL,
  method        = "dml_plr",
  learner       = "linear",
  K             = 5L,
  S             = 3L,
  eval_g        = NULL,
  trim          = 0.05,
  boot          = 500L,
  alpha         = 0.05,
  signal_map    = "local_poly",
  map_degree    = 2L,
  map_df        = NULL,
  density_method = "residual_kde",
  dml_inference = "score",
  dml_boot      = 500L,
  dml_multiplier = "normal",
  missing_data  = "listwise",
  the_call      = NULL
) {

  ycol <- function(t) paste0("Y_", t)
  dml_type <- sub("^dml_", "", method)
  signal_map <- if (dml_type == "flex") .normalize_signal_map(signal_map) else NULL
  density_method <- if (dml_type %in% c("flex", "incremental")) {
    .normalize_density_method(density_method)
  } else {
    NULL
  }
  dml_inference <- .normalize_dml_inference(dml_inference)
  dml_multiplier <- .normalize_dml_multiplier(dml_multiplier)

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

  needed_cols <- c("G", covar)
  keep <- if (missing_data == "listwise")
    complete.cases(s[, c(needed_cols, ycol(all_times)), drop = FALSE])
  else
    complete.cases(s[, needed_cols, drop = FALSE])
  s <- s[keep, , drop = FALSE]

  n     <- nrow(s)
  p     <- length(covar)
  G_vec <- as.numeric(s$G)
  X_mat <- if (p > 0) as.matrix(s[, covar, drop = FALSE]) else matrix(0, n, 0)
  cluster_vec <- if ("c" %in% names(s)) s$c else NULL

  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) ||
      K < 2L || K > n) {
    stop("K must be between 2 and the number of complete observations (", n, ").")
  }
  if (!is.numeric(S) || length(S) != 1L || !is.finite(S) || S < 1L) {
    stop("S must be a positive integer-like value.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single finite number strictly between 0 and 1.")
  }
  if (!is.numeric(trim) || length(trim) != 1L || !is.finite(trim) ||
      trim < 0 || trim >= 0.5) {
    stop("trim must be a single finite number in [0, 0.5).")
  }
  if (!is.numeric(boot) || length(boot) != 1L || !is.finite(boot) || boot < 1L) {
    stop("boot must be a positive integer-like value.")
  }
  if (!is.numeric(dml_boot) || length(dml_boot) != 1L ||
      !is.finite(dml_boot) || dml_boot < 1L) {
    stop("dml_boot must be a positive integer-like value.")
  }
  K <- as.integer(K)
  S <- as.integer(S)
  boot <- as.integer(boot)
  dml_boot <- as.integer(dml_boot)

  # Validate G type
  if (dml_type == "binary") {
    g_unique <- sort(unique(G_vec))
    if (!identical(g_unique, c(0, 1)))
      stop("method='dml_binary' requires binary G in {0, 1}.")
  } else {
    if (length(unique(G_vec)) <= 2)
      stop("method='", method, "' requires continuous G (more than 2 unique values).")
  }

  # eval_g for flexible level and derivative curves
  if (dml_type == "flex") {
    if (is.null(eval_g)) {
      g_lo   <- as.numeric(quantile(G_vec, trim))
      g_hi   <- as.numeric(quantile(G_vec, 1 - trim))
      eval_g <- seq(g_lo, g_hi, length.out = 50L)
    } else {
      eval_g <- sort(as.numeric(eval_g))
    }
  }
  m_eval <- if (!is.null(eval_g)) length(eval_g) else 0L

  # ── Helper: run DML on one tempY vector ──
  run_one <- function(DeltaY, G, X, cluster = NULL) {
    if (dml_type == "binary") {
      res <- .est_dml_binary(
        DeltaY, G, X, K, S, alpha, learner,
        cluster = cluster,
        dml_inference = dml_inference,
        dml_boot = dml_boot,
        dml_multiplier = dml_multiplier
      )
      list(result = res$result, ps = res$ps, curve = NULL, scalar = res)
    } else if (dml_type == "plr") {
      res <- .est_dml_plr(
        DeltaY, G, X, K, S, alpha, learner,
        cluster = cluster,
        dml_inference = dml_inference,
        dml_boot = dml_boot,
        dml_multiplier = dml_multiplier
      )
      list(result = res$result, ps = NULL, curve = NULL, scalar = res)
    } else if (dml_type == "flex") {
      cv <- .est_dml_flex(
        DeltaY, G, X, eval_g, K, S, alpha, learner,
        signal_map = signal_map,
        map_degree = map_degree,
        map_df = map_df,
        density_method = density_method,
        boot = boot
      )
      scalar_est <- mean(cv$theta_hat, na.rm = TRUE)
      scalar_se  <- mean(cv$se_theta, na.rm = TRUE)
      r <- c(Estimate  = scalar_est,
             Std.Error = scalar_se,
             CI_Lower  = scalar_est - stats::qnorm(1 - alpha/2) * scalar_se,
             CI_Upper  = scalar_est + stats::qnorm(1 - alpha/2) * scalar_se)
      list(result = r, ps = NULL, curve = cv)
    } else if (dml_type == "incremental") {
      inc <- .est_dml_incremental_scalar(
        DeltaY, G, X, K, S, alpha, learner,
        density_method = density_method,
        cluster = cluster,
        dml_inference = dml_inference,
        dml_boot = dml_boot,
        dml_multiplier = dml_multiplier
      )
      list(result = inc$result, ps = NULL, curve = NULL,
           scalar = inc, incremental = inc)
    }
  }

  # ── Event-period tempY ──
  Y_ref_col  <- ycol(ref_period)
  Y_tr_cols  <- setdiff(ycol(tr_period), Y_ref_col)
  tempY_event <- rowMeans(s[, Y_tr_cols, drop = FALSE], na.rm = TRUE) - s[[Y_ref_col]]

  event_res <- run_one(tempY_event, G_vec, X_mat, cluster = cluster_vec)
  result_event <- event_res$result
  stored_ps    <- event_res$ps
  curve_event  <- event_res$curve
  scalar_event <- if (!is.null(event_res$scalar)) event_res$scalar else NULL
  incremental_event <- if (!is.null(event_res$incremental)) event_res$incremental else NULL

  est_event_df <- as.data.frame(t(result_event))

  # ── Dynamic loop ──
  dynamic_df <- data.frame(
    Estimate = numeric(length(all_times)), Std.Error = numeric(length(all_times)),
    CI_Lower = numeric(length(all_times)), CI_Upper  = numeric(length(all_times)),
    row.names = as.character(all_times)
  )

  for (i in seq_along(all_times)) {
    t_yr <- all_times[i]
    if (t_yr == ref_period) { dynamic_df[i, ] <- 0; next }
    t_col <- ycol(t_yr)
    if (!t_col %in% names(s)) next
    tempY_t <- s[[t_col]] - s[[Y_ref_col]]
    ok      <- !is.na(tempY_t)
    if (sum(ok) < max(10L, K)) next
    dres <- run_one(tempY_t[ok], G_vec[ok],
                    if (p > 0) X_mat[ok, , drop = FALSE] else X_mat,
                    cluster = if (!is.null(cluster_vec)) cluster_vec[ok] else NULL)
    dynamic_df[i, ] <- dres$result[c("Estimate", "Std.Error", "CI_Lower", "CI_Upper")]
  }

  # ── Pre/post aggregates ──
  pre_times  <- all_times[all_times < min(tr_period)]
  post_times <- all_times[all_times > max(tr_period)]

  dml_aggregate <- function(time_vec) {
    na_df <- data.frame(Estimate=NA, Std.Error=NA, CI_Lower=NA, CI_Upper=NA)
    if (length(time_vec) == 0) return(na_df)
    t_cols <- setdiff(unique(ycol(time_vec)), Y_ref_col)
    if (length(t_cols) == 0) return(na_df)
    tmp <- rowMeans(s[, t_cols, drop = FALSE], na.rm = TRUE)
    tmp[is.nan(tmp)] <- NA_real_
    DY_a <- tmp - s[[Y_ref_col]]
    ok   <- !is.na(DY_a)
    if (sum(ok) < max(10L, K)) return(na_df)
    r <- run_one(DY_a[ok], G_vec[ok],
                 if (p > 0) X_mat[ok, , drop = FALSE] else X_mat,
                 cluster = if (!is.null(cluster_vec)) cluster_vec[ok] else NULL)
    as.data.frame(t(r$result[c("Estimate","Std.Error","CI_Lower","CI_Upper")]))
  }

  pre_event_result  <- dml_aggregate(pre_times)
  post_event_result <- dml_aggregate(post_times)

  # ── raw_means: median-split approximation ──
  zval  <- stats::qnorm(0.975)
  G_med <- stats::median(G_vec, na.rm = TRUE)
  rawdf <- do.call(rbind, lapply(all_times, function(yr) {
    col_yr <- ycol(yr)
    if (!col_yr %in% names(s)) return(NULL)
    y1 <- s[[col_yr]][G_vec >= G_med]; y0 <- s[[col_yr]][G_vec < G_med]
    n1 <- sum(!is.na(y1)); m1 <- mean(y1, na.rm = TRUE)
    n0 <- sum(!is.na(y0)); m0 <- mean(y0, na.rm = TRUE)
    se1 <- if (n1 <= 1) NA_real_ else stats::sd(y1, na.rm = TRUE) / sqrt(n1)
    se0 <- if (n0 <= 1) NA_real_ else stats::sd(y0, na.rm = TRUE) / sqrt(n0)
    rbind(
      data.frame(time=yr, group="Group 1", meanY=m1, n=n1, seY=se1,
                 CI_Lower=m1-zval*se1, CI_Upper=m1+zval*se1, stringsAsFactors=FALSE),
      data.frame(time=yr, group="Group 0", meanY=m0, n=n0, seY=se0,
                 CI_Lower=m0-zval*se0, CI_Upper=m0+zval*se0, stringsAsFactors=FALSE)
    )
  }))
  if (is.null(rawdf)) rawdf <- data.frame(
    time=numeric(0), group=character(0), meanY=numeric(0), n=integer(0),
    seY=numeric(0), CI_Lower=numeric(0), CI_Upper=numeric(0), stringsAsFactors=FALSE)

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
    method        = method,
    vartype       = "robust",
    alpha         = alpha,
    times         = numeric_times,
    G             = G_vec,
    ps            = stored_ps,
    cluster       = cluster_vec,
    call          = the_call,
    target.pop    = "all",
    dml_method    = dml_type,
    K             = K,
    S             = S,
    learner       = learner,
    signal_map    = signal_map,
    map_degree    = map_degree,
    map_df        = map_df,
    density_method = density_method,
    dml_inference = dml_inference,
    dml_boot      = dml_boot,
    dml_multiplier = dml_multiplier,
    eval_g        = eval_g,
    curve_event   = curve_event,
    scalar_event  = scalar_event,
    incremental_event = incremental_event,
    dml_metadata = if (!is.null(curve_event) && !is.null(curve_event$metadata)) {
      curve_event$metadata
    } else if (!is.null(scalar_event) && !is.null(scalar_event$metadata)) {
      scalar_event$metadata
    } else {
      NULL
    },
    dml_diagnostics = if (!is.null(curve_event) && !is.null(curve_event$diagnostics)) {
      curve_event$diagnostics
    } else if (!is.null(incremental_event) && !is.null(incremental_event$diagnostics)) {
      incremental_event$diagnostics
    } else {
      NULL
    }
  )
  class(out) <- "fdid"
  out
}
