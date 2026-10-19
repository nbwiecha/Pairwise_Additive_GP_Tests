# =========================================================
# REML score test script for carrying out data analysis
# =========================================================

data_analysis_pipeline_sp <- function(Y, X, Z_cov, coords_sp, init_params_full=c(0,0,0,0), alpha=.05) {
  
  require(PairwiseSims)
  Z_cov <- as.matrix(Z_cov)
  coords_sp <- as.matrix(coords_sp)
  p <- ncol(Z_cov)
  
  # 1. Estimate global model to get the anchor for plots and interactions[cite: 5]
  fit_full <- fit_unified_full_model_reml_sp(Y, X, Z_cov, coords_sp, init_params = init_params_full)
  
  if (is.null(fit_full) || fit_full$convergence != 0) {
    warning("Global unified REML model failed to converge. Proceeding with default init.")
    optimal_rhos <- c(init_params_full[1], init_params_full[3])
  } else {
    optimal_rhos <- c(fit_full$kernel_params[1], fit_full$kernel_params[3])
  }
  
  # 2. Test for variable importance for each coordinate (Sequential)[cite: 5]
  var.importance.pvals <- rep(NA, p)
  names(var.importance.pvals) <- if(!is.null(colnames(Z_cov))) colnames(Z_cov) else paste0("Var", 1:p)
  
  for (j in 1:p) {
    score_results <- tryCatch({
      run_unified_pipeline_reml_sp(Y, X, Z_cov, coords_sp, 
                                   init_params_full = init_params_full,
                                   test_type = "var_imp", idx_args = j, 
                                   global_rhos = optimal_rhos)
    }, error = function(e) {
      message("Error on Var ", j, ": ", e$message)
      return(NULL)
    })
    
    if (!is.null(score_results)) {
      var.importance.pvals[j] <- score_results$p_value
    }
  }
  
  # 3. Test for interaction between any variables with p-value < alpha (Sequential)[cite: 5]
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
        run_unified_pipeline_reml_sp(Y, X, Z_cov, coords_sp, 
                                     init_params_full = init_params_full,
                                     test_type = "interaction", 
                                     idx_args = c(interaction.tests[i, 1], interaction.tests[i, 2]), 
                                     global_rhos = optimal_rhos)
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
    global_rhos = optimal_rhos 
  ))
}

# Runs one spatial score test and returns its p-value, or "ERROR: <message>" if it
# fails (including a NULL result), so worker failures are visible to the caller.
run_test_or_error_sp <- function(Y, X, Z_cov, coords_sp, init_params_full, test_type, idx_args, global_rhos) {
  res <- tryCatch(
    run_unified_pipeline_reml_sp(Y, X, Z_cov, coords_sp, init_params_full = init_params_full,
                                 test_type = test_type, idx_args = idx_args, global_rhos = global_rhos),
    error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(res)) return(res)
  if (is.null(res) || is.null(res$p_value) || is.na(res$p_value)) {
    return(paste("ERROR: test returned no p-value for", test_type, paste(idx_args, collapse = ",")))
  }
  res$p_value
}

# Stops with a list of every failed test, instead of passing NA p-values on.
stop_if_test_errors_sp <- function(results, labels) {
  failed <- grepl("^ERROR", results)
  if (any(failed)) {
    stop(sprintf("%d of %d tests failed:\n%s", sum(failed), length(results),
                 paste0("  ", labels[failed], " -> ", results[failed], collapse = "\n")),
         call. = FALSE)
  }
  invisible(TRUE)
}

data_analysis_pipeline_parallel_sp <- function(Y, X, Z_cov, coords_sp, init_params_full=c(0,0,0,0), alpha=.05) {
  
  require(PairwiseSims)
  require(foreach) 
  
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  Z_cov <- as.matrix(Z_cov)
  coords_sp <- as.matrix(coords_sp)
  p <- ncol(Z_cov)
  
  # 1. Estimate global model to get the anchor for plots and interactions[cite: 5]
  fit_full <- fit_unified_full_model_reml_sp(Y, X, Z_cov, coords_sp, init_params = init_params_full)
  
  # Fail loudly rather than testing at the starting values of the lengthscales.
  if (is.null(fit_full) || fit_full$convergence != 0) {
    stop("Global unified REML model failed to converge (convergence code: ",
         if (is.null(fit_full)) "NULL fit" else fit_full$convergence,
         "). Not proceeding with the default rhos.")
  }
  optimal_rhos <- c(fit_full$kernel_params[1], fit_full$kernel_params[3])

  funcs_to_export <- c("fit_unified_full_model_reml_sp",
                       "fit_unified_H0_fixed_rho_sp_reml",
                       "run_unified_score_test_reml_sp",
                       "run_unified_pipeline_reml_sp",
                       "run_test_or_error_sp")

  var_names <- if(!is.null(colnames(Z_cov))) colnames(Z_cov) else paste0("Var", 1:p)

  # 2. Variable Importance Tests (Parallelized). Each worker returns a p-value or
  # an "ERROR: ..." string so that failures are reported, not turned into NA.
  var_imp_res <- foreach(j = 1:p, .combine = c,
                         .packages = c("PairwiseSims", "MASS"),
                         .export = funcs_to_export) %dopar% {
                           run_test_or_error_sp(Y, X, Z_cov, coords_sp, init_params_full,
                                                "var_imp", j, optimal_rhos)
                         }
  stop_if_test_errors_sp(var_imp_res, paste("variable importance:", var_names))
  var.importance.pvals <- as.numeric(var_imp_res)

  names(var.importance.pvals) <- var_names
  
  # 3. Interaction Tests (Parallelized)[cite: 5]
  sig.pair <- which(var.importance.pvals < alpha)
  
  if (length(sig.pair) > 1) {
    interaction.tests <- t(combn(sig.pair, 2))
    
    var_names <- names(var.importance.pvals)
    test_names <- paste(var_names[interaction.tests[, 1]], 
                        var_names[interaction.tests[, 2]], sep=" x ")
    
    int_res <- foreach(i = 1:nrow(interaction.tests),
                       .combine = c, .packages = c("PairwiseSims", "MASS"),
                       .export = funcs_to_export) %dopar% {
                         run_test_or_error_sp(Y, X, Z_cov, coords_sp, init_params_full, "interaction",
                                              c(interaction.tests[i, 1], interaction.tests[i, 2]), optimal_rhos)
                       }
    stop_if_test_errors_sp(int_res, paste("interaction:", test_names))
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
    global_rhos = optimal_rhos 
  ))
}

# ---------------------------------------------------------
# The Core Prediction Engine and Plot Wrappers (Updated for SP)
# ---------------------------------------------------------
predict_unified_surface_with_se_sp <- function(fit_full, Y_train, X_train, Z_train, coords_train, 
                                               Z_new, coords_new, X_new) {
  
  Z_train <- as.matrix(Z_train)
  coords_train <- as.matrix(coords_train)
  Z_new <- as.matrix(Z_new)
  coords_new <- as.matrix(coords_new)
  X_train <- as.matrix(X_train)
  X_new <- as.matrix(X_new)
  
  n_train <- nrow(Z_train)
  n_new <- nrow(Z_new)
  p_fixed <- ncol(X_train)
  P <- ncol(Z_train)
  L_total <- P + (P * (P - 1)) / 2
  
  rho_a  <- exp(fit_full$kernel_params[1])
  lam_a  <- exp(fit_full$kernel_params[2])
  rho_sp <- exp(fit_full$kernel_params[3])
  lam_sp <- exp(fit_full$kernel_params[4])
  
  beta_hat <- fit_full$beta_hat
  sigma2_hat <- fit_full$sigma2_hat
  
  # Fetch the original W matrix to get alpha and W_inv
  cpp_train <- eval_unified_full_model_sp(Z_train, coords_train, 
                                          fit_full$kernel_params[1], fit_full$kernel_params[2], 
                                          fit_full$kernel_params[3], fit_full$kernel_params[4], n_train)
  
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
    # Train distances for centering
    D_train <- outer(Z_train[, j], Z_train[, j], "-")^2
    R_train <- exp(-D_train / rho_a)
    col_means_train <- colMeans(R_train)
    grand_mean <- mean(R_train)
    
    # Cross distances
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
  
  # Spatial Cross-Covariance
  D2_sp_cross <- outer(coords_new[, 1], coords_train[, 1], "-")^2 + 
    outer(coords_new[, 2], coords_train[, 2], "-")^2
  K_sp_cross <- exp(-D2_sp_cross / rho_sp)
  K_sp_diag  <- rep(1, n_new)
  
  # Combine Scaled Components
  scaled_K_cross <- lam_a * K_anova_cross + lam_sp * K_sp_cross
  scaled_K_diag  <- lam_a * K_anova_diag  + lam_sp * K_sp_diag
  
  Y_pred <- (X_new %*% beta_hat) + (scaled_K_cross %*% alpha)
  
  # Standard Errors[cite: 5]
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

plot_main_effect_slice_sp <- function(target_idx, fit_full, Y, X, Z_cov, coords_sp, 
                                      orig_min = 0, orig_max = 1, unlog = FALSE, 
                                      grid_points = 100) {
  
  Z_cov <- as.matrix(Z_cov)
  coords_sp <- as.matrix(coords_sp)
  X <- as.matrix(X)
  p <- ncol(Z_cov)
  
  var_name <- if(!is.null(colnames(Z_cov))) colnames(Z_cov)[target_idx] else paste("Variable", target_idx)
  
  # Set non-target ANOVA variables to their median[cite: 5]
  Z_new <- matrix(apply(Z_cov, 2, median), nrow = grid_points, ncol = p, byrow = TRUE)
  x_seq <- seq(min(Z_cov[, target_idx]), max(Z_cov[, target_idx]), length.out = grid_points)
  Z_new[, target_idx] <- x_seq
  
  # Set spatial coordinates to their median for the slice
  coords_new <- matrix(apply(coords_sp, 2, median), nrow = grid_points, ncol = ncol(coords_sp), byrow = TRUE)
  
  X_new <- matrix(apply(X, 2, median), nrow = grid_points, ncol = ncol(X), byrow = TRUE)
  if (all(X[, 1] == 1)) X_new[, 1] <- 1 
  
  pred_df <- predict_unified_surface_with_se_sp(fit_full, Y, X, Z_cov, coords_sp, Z_new, coords_new, X_new)
  
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

plot_interaction_slice_sp <- function(idx_main, idx_quant, fit_full, Y, X, Z_cov, coords_sp, 
                                      orig_min_main = 0, orig_max_main = 1, unlog_main = FALSE,
                                      orig_min_quant = 0, orig_max_quant = 1, unlog_quant = FALSE,
                                      grid_points = 100) {
  
  Z_cov <- as.matrix(Z_cov)
  coords_sp <- as.matrix(coords_sp)
  X <- as.matrix(X)
  p <- ncol(Z_cov)
  
  name_main <- if(!is.null(colnames(Z_cov))) colnames(Z_cov)[idx_main] else paste("Variable", idx_main)
  name_quant <- if(!is.null(colnames(Z_cov))) colnames(Z_cov)[idx_quant] else paste("Variable", idx_quant)
  
  quantiles <- quantile(Z_cov[, idx_quant], probs = c(0.10, 0.50, 0.90))[cite: 5]
  
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
    
    coords_new <- matrix(apply(coords_sp, 2, median), nrow = grid_points, ncol = ncol(coords_sp), byrow = TRUE)
    
    X_new <- matrix(apply(X, 2, median), nrow = grid_points, ncol = ncol(X), byrow = TRUE)
    if (all(X[, 1] == 1)) X_new[, 1] <- 1
    
    pred_df <- predict_unified_surface_with_se_sp(fit_full, Y, X, Z_cov, coords_sp, Z_new, coords_new, X_new)
    
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