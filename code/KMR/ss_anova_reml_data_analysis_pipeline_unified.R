# REML score test script for carrying out data analysis

data_analysis_pipeline <- function(Y, X, Z_cov, init_params_full=c(1,1), alpha=.05) {
  
  require(PairwiseSims)
  Z_cov <- as.matrix(Z_cov)
  p <- ncol(Z_cov)
  
  # 1. Estimate global model to get the anchor for plots and interactions[cite: 4, 5]
  fit_full <- fit_unified_full_model_reml(Y, X, Z_cov, init_params = init_params_full)
  
  if (is.null(fit_full) || fit_full$convergence != 0) {
    warning("Global unified REML model failed to converge. Proceeding with default init.")
    optimal_log_rho <- init_params_full[1]
  } else {
    optimal_log_rho <- fit_full$kernel_params[1]
  }
  
  # 2. Test for variable importance for each coordinate (Sequential)[cite: 4, 5]
  var.importance.pvals <- rep(NA, p)
  names(var.importance.pvals) <- if(!is.null(colnames(Z_cov))) colnames(Z_cov) else paste0("Var", 1:p)
  
  for (j in 1:p) {
    score_results <- tryCatch({
      run_unified_var_imp_pipeline_reml(Y, X, Z_cov, test_idx_vec = j, global_log_rho = optimal_log_rho)
    }, error = function(e) {
      message("Error on Var ", j, ": ", e$message)
      return(NULL)
    })
    
    if (!is.null(score_results)) {
      var.importance.pvals[j] <- score_results$p_value
    }
  }
  
  # 3. Test for interaction between any variables with p-value < alpha (Sequential)[cite: 4, 5]
  sig.pair <- which(var.importance.pvals < alpha)
  
  if (length(sig.pair) > 1) {
    interaction.tests <- t(combn(sig.pair, 2))
    interaction.tests <- cbind(interaction.tests, rep(NA, nrow(interaction.tests)))
    row.names(interaction.tests) <- 1:nrow(interaction.tests)
    colnames(interaction.tests) <- c("Var1", "Var2", "p-value")
    
    var_names <- names(var.importance.pvals)
    for (i in 1:nrow(interaction.tests)) {
      row.names(interaction.tests)[i] <- paste(var_names[interaction.tests[i, 1]],
                                               var_names[interaction.tests[i, 2]], 
                                               sep=" x ")
    }
    
    for (i in 1:nrow(interaction.tests)) {
      score_res <- tryCatch({
        run_unified_interaction_pipeline_reml(Y, X, Z_cov, 
                                              idx1 = interaction.tests[i, 1], 
                                              idx2 = interaction.tests[i, 2], 
                                              global_log_rho = optimal_log_rho)
      }, error = function(e) {
        message("Error on Interaction ", i, ": ", e$message)
        return(NULL)
      })
      
      if (!is.null(score_res)) {
        interaction.tests[i, 3] <- score_res$p_value
      }
    }
    
  } else {
    interaction.tests <- NULL
  }
  
  return(list(
    var.importance.pvals = var.importance.pvals,
    interaction.pvals = interaction.tests,
    global_log_rho = optimal_log_rho 
  ))
}

# Runs one score test and returns its p-value, or "ERROR: <message>" if it fails
# (including a NULL result), so worker failures are visible to the caller.
run_test_or_error <- function(Y, X, Z_cov, init_params_full, test_type, idx_args, global_rho) {
  res <- tryCatch(
    run_unified_pipeline_reml(Y, X, Z_cov, init_params_full = init_params_full,
                              test_type = test_type, idx_args = idx_args, global_rho = global_rho),
    error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(res)) return(res)
  if (is.null(res) || is.null(res$p_value) || is.na(res$p_value)) {
    return(paste("ERROR: test returned no p-value for", test_type, paste(idx_args, collapse = ",")))
  }
  res$p_value
}

# Stops with a list of every failed test, instead of passing NA p-values on.
stop_if_test_errors <- function(results, labels) {
  failed <- grepl("^ERROR", results)
  if (any(failed)) {
    stop(sprintf("%d of %d tests failed:\n%s", sum(failed), length(results),
                 paste0("  ", labels[failed], " -> ", results[failed], collapse = "\n")),
         call. = FALSE)
  }
  invisible(TRUE)
}

data_analysis_pipeline_parallel <- function(Y, X, Z_cov, init_params_full=c(1,1), alpha=.05) {
  require(foreach) 
  require(PairwiseSims)
  
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  Z_cov <- as.matrix(Z_cov)
  p <- ncol(Z_cov)
  
  fit_full <- fit_unified_full_model_reml(Y, X, Z_cov, init_params = init_params_full)

  # Fail loudly: silently falling back to the starting value of rho would give
  # tests based on an arbitrary lengthscale.
  if (is.null(fit_full) || fit_full$convergence != 0) {
    stop("Global unified REML model failed to converge (convergence code: ",
         if (is.null(fit_full)) "NULL fit" else fit_full$convergence,
         "). Not proceeding with the default rho.")
  }
  optimal_log_rho <- fit_full$kernel_params[1]

  # FIXED: Explicitly export the modern master functions to the worker nodes
  funcs_to_export <- c("fit_unified_full_model_reml",
                       "fit_unified_H0_fixed_rho_reml",
                       "run_unified_score_test_reml",
                       "run_unified_pipeline_reml",
                       "run_test_or_error",
                       "eval_unified_test_kernels",
                       "sum_prod_blocked",
                       "tr_prod_blocked")

  var_names <- if(!is.null(colnames(Z_cov))) colnames(Z_cov) else paste0("Var", 1:p)

  # Each worker returns either a p-value or an "ERROR: ..." string, so failures
  # (e.g. std::bad_alloc) are reported instead of silently becoming NA.
  var_imp_res <- foreach(j = 1:p, .combine = c,
                         .packages = c("MASS", "PairwiseSims"),
                         .export = funcs_to_export) %dopar% {
                           run_test_or_error(Y, X, Z_cov, init_params_full, "var_imp", j, optimal_log_rho)
                         }
  stop_if_test_errors(var_imp_res, paste("variable importance:", var_names))
  var.importance.pvals <- as.numeric(var_imp_res)

  names(var.importance.pvals) <- var_names
  
  sig.pair <- which(var.importance.pvals < alpha)
  
  if (length(sig.pair) > 1) {
    interaction.tests <- t(combn(sig.pair, 2))
    var_names <- names(var.importance.pvals)
    test_names <- paste(var_names[interaction.tests[, 1]], var_names[interaction.tests[, 2]], sep=" x ")
    
    int_res <- foreach(i = 1:nrow(interaction.tests), .combine = c,
                       .packages = c("MASS", "PairwiseSims"), .export = funcs_to_export) %dopar% {
                         run_test_or_error(Y, X, Z_cov, init_params_full, "interaction",
                                           c(interaction.tests[i, 1], interaction.tests[i, 2]), optimal_log_rho)
                       }
    stop_if_test_errors(int_res, paste("interaction:", test_names))
    int_pvals <- as.numeric(int_res)

    interaction_df <- data.frame(
      Var1 = interaction.tests[, 1],
      Var2 = interaction.tests[, 2],
      `p-value` = int_pvals,
      row.names = test_names,
      check.names = FALSE
    )
  } else {
    interaction_df <- NULL
  }
  
  return(list(
    var.importance.pvals = var.importance.pvals,
    interaction.pvals = interaction_df,
    global_log_rho = optimal_log_rho 
  ))
}

# ---------------------------------------------------------
# The Core Prediction Engine and Plot Wrappers (Updated for Raw Matrices)
# ---------------------------------------------------------
predict_unified_surface_with_se <- function(fit_full, Y_train, X_train, Z_train, Z_new, X_new) {
  
  Z_train <- as.matrix(Z_train)
  Z_new <- as.matrix(Z_new)
  X_train <- as.matrix(X_train)
  X_new <- as.matrix(X_new)
  
  n_train <- nrow(Z_train)
  n_new <- nrow(Z_new)
  p_fixed <- ncol(X_train)
  P <- ncol(Z_train)
  L_total <- P + (P * (P - 1)) / 2
  
  log_rho <- fit_full$kernel_params[1]
  log_lam <- fit_full$kernel_params[2]
  rho_a <- exp(log_rho)
  lam_a <- exp(log_lam)
  
  beta_hat <- fit_full$beta_hat
  sigma2_hat <- fit_full$sigma2_hat
  
  cpp_train <- eval_unified_full_model(Z_train, log_rho, log_lam, n_train)
  
  U <- try(chol(cpp_train$W), silent = TRUE)
  if (inherits(U, "try-error")) stop("W matrix computationally singular.")
  W_inv <- chol2inv(U)
  
  alpha <- W_inv %*% (Y_train - (X_train %*% beta_hat))
  
  # --- Pure R implementation of Cross-Covariance ---
  S_1D_cross <- matrix(0, n_new, n_train)
  S_sq_cross <- matrix(0, n_new, n_train)
  S_1D_diag  <- rep(0, n_new)
  S_sq_diag  <- rep(0, n_new)
  
  for (j in 1:P) {
    # Train distances for centering[cite: 5]
    D_train <- outer(Z_train[, j], Z_train[, j], "-")^2
    R_train <- exp(-D_train / rho_a)
    col_means_train <- colMeans(R_train)
    grand_mean <- mean(R_train)
    
    # Cross distances[cite: 5]
    D_cross <- outer(Z_new[, j], Z_train[, j], "-")^2
    R_cross <- exp(-D_cross / rho_a)
    row_means_cross <- rowMeans(R_cross)
    
    # Orthogonal projection[cite: 5]
    K_cross_j <- sweep(R_cross, 1, row_means_cross, "-")
    K_cross_j <- sweep(K_cross_j, 2, col_means_train, "-")
    K_cross_j <- K_cross_j + grand_mean
    
    S_1D_cross <- S_1D_cross + K_cross_j
    S_sq_cross <- S_sq_cross + K_cross_j^2
    
    # Diagonal components[cite: 5]
    K_diag_j <- 1.0 - 2.0 * row_means_cross + grand_mean
    S_1D_diag <- S_1D_diag + K_diag_j
    S_sq_diag <- S_sq_diag + K_diag_j^2
  }
  
  K_anova_cross <- (S_1D_cross + 0.5 * (S_1D_cross^2 - S_sq_cross)) / L_total
  K_anova_diag  <- (S_1D_diag + 0.5 * (S_1D_diag^2 - S_sq_diag)) / L_total
  
  scaled_K_cross <- lam_a * K_anova_cross
  scaled_K_diag  <- lam_a * K_anova_diag
  
  Y_pred <- (X_new %*% beta_hat) + (scaled_K_cross %*% alpha)
  
  var_reduction <- rowSums((scaled_K_cross %*% W_inv) * scaled_K_cross)
  
  X_tilde <- forwardsolve(t(U), X_train)
  cov_beta <- chol2inv(qr.R(qr(X_tilde))) 
  
  R_mat <- X_new - (scaled_K_cross %*% W_inv %*% X_train)
  var_fixed <- rowSums((R_mat %*% cov_beta) * R_mat)
  
  se2 <- sigma2_hat * (scaled_K_diag - var_reduction + var_fixed)
  se <- sqrt(pmax(se2, 1e-10)) 
  t_crit <- qt(0.975, df = n_train - p_fixed)
  
  return(data.frame(
    Predicted_Y = as.numeric(Y_pred),
    SE = as.numeric(se),
    CI_Lower = as.numeric(Y_pred - t_crit * se),
    CI_Upper = as.numeric(Y_pred + t_crit * se)
  ))
}

library(ggplot2)

plot_main_effect_slice <- function(target_idx, fit_full, Y, X, Z_cov, 
                                   orig_min = 0, orig_max = 1, unlog = FALSE, 
                                   grid_points = 100) {
  
  Z_cov <- as.matrix(Z_cov)
  X <- as.matrix(X)
  p <- ncol(Z_cov)
  
  var_name <- if(!is.null(colnames(Z_cov))) colnames(Z_cov)[target_idx] else paste("Variable", target_idx)
  
  Z_new <- matrix(apply(Z_cov, 2, median), nrow = grid_points, ncol = p, byrow = TRUE)
  x_seq <- seq(min(Z_cov[, target_idx]), max(Z_cov[, target_idx]), length.out = grid_points)
  Z_new[, target_idx] <- x_seq
  
  X_new <- matrix(apply(X, 2, median), nrow = grid_points, ncol = ncol(X), byrow = TRUE)
  if (all(X[, 1] == 1)) X_new[, 1] <- 1 
  
  pred_df <- predict_unified_surface_with_se(fit_full, Y, X, Z_cov, Z_new, X_new)
  
  x_real <- (x_seq * (orig_max - orig_min)) + orig_min
  if (unlog) x_real <- exp(x_real)
  
  plot_data <- cbind(data.frame(X_val = x_real), pred_df)
  
  ggplot(plot_data, aes(x = X_val, y = Predicted_Y)) +
    geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), fill = "blue", alpha = 0.2) +
    geom_line(color = "blue", linewidth = 1.2) +
    theme_minimal() +
    labs(title = paste("Conditional Main Effect Slice for", var_name),
         subtitle = "Shaded region represents 95% pointwise credible interval",
         x = var_name,
         y = "Predicted Y")
}

plot_interaction_slice <- function(idx_main, idx_quant, fit_full, Y, X, Z_cov, 
                                   orig_min_main = 0, orig_max_main = 1, unlog_main = FALSE,
                                   orig_min_quant = 0, orig_max_quant = 1, unlog_quant = FALSE,
                                   grid_points = 100) {
  
  Z_cov <- as.matrix(Z_cov)
  X <- as.matrix(X)
  p <- ncol(Z_cov)
  
  name_main <- if(!is.null(colnames(Z_cov))) colnames(Z_cov)[idx_main] else paste("Variable", idx_main)
  name_quant <- if(!is.null(colnames(Z_cov))) colnames(Z_cov)[idx_quant] else paste("Variable", idx_quant)
  
  quantiles <- quantile(Z_cov[, idx_quant], probs = c(0.10, 0.50, 0.90))[cite: 4, 5]
  
  quant_real <- (quantiles * (orig_max_quant - orig_min_quant)) + orig_min_quant
  if (unlog_quant) quant_real <- exp(quant_real)
  
  quant_labels <- paste0(c("10th (", "50th (", "90th ("), round(quant_real, 2), ")")
  
  plot_data <- data.frame()
  
  for (i in seq_along(quantiles)) {
    q_val <- quantiles[i]
    q_label <- quant_labels[i]
    
    Z_new <- matrix(apply(Z_cov, 2, median), nrow = grid_points, ncol = p, byrow = TRUE)
    x_seq <- seq(min(Z_cov[, idx_main]), max(Z_cov[, idx_main]), length.out = grid_points)
    
    Z_new[, idx_main] <- x_seq
    Z_new[, idx_quant] <- q_val
    
    X_new <- matrix(apply(X, 2, median), nrow = grid_points, ncol = ncol(X), byrow = TRUE)
    if (all(X[, 1] == 1)) X_new[, 1] <- 1
    
    pred_df <- predict_unified_surface_with_se(fit_full, Y, X, Z_cov, Z_new, X_new)
    
    x_real <- (x_seq * (orig_max_main - orig_min_main)) + orig_min_main
    if (unlog_main) x_real <- exp(x_real)
    
    plot_data <- rbind(plot_data, cbind(
      data.frame(X_main = x_real, Quantile_Level = factor(q_label, levels = quant_labels)), 
      pred_df
    ))
  }
  
  ggplot(plot_data, aes(x = X_main, y = Predicted_Y, color = Quantile_Level, fill = Quantile_Level, group = Quantile_Level)) +
    geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.2) +
    theme_minimal() +
    scale_color_viridis_d(name = paste(name_quant, "\n(Quantiles)")) +
    scale_fill_viridis_d(name = paste(name_quant, "\n(Quantiles)")) +
    labs(title = paste("Interaction:", name_main, "conditional on", name_quant),
         subtitle = "Shaded regions represent 95% pointwise credible intervals",
         x = name_main,
         y = "Predicted Y")
}