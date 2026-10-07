# =========================================================
# FILE: ss_anova_reml_unified.R
# SS-ANOVA REML Approach with Full-State Caching & Analysis Pipelines
# =========================================================

# 1. Fit Global Unified Model (2 Parameters)
fit_unified_full_model_reml <- function(Y, X, Z_cov, init_params) {
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  Z_cov <- as.matrix(Z_cov)
  n <- length(Y)
  p <- ncol(X)
  
  # --- FULL-STATE CACHE ---
  last_params <- NULL
  last_state <- NULL
  
  get_state <- function(params) {
    if (!is.null(last_params) && all(params == last_params)) return(last_state)
    
    cpp_res <- eval_unified_full_model(Z_cov, params[1], params[2], n)
    U <- try(chol(cpp_res$W), silent = TRUE)
    
    if (inherits(U, "try-error")) {
      state <- list(error = TRUE)
    } else {
      X_tilde <- forwardsolve(t(U), X)
      Y_tilde <- forwardsolve(t(U), Y)
      
      beta_hat <- tryCatch({ qr.solve(X_tilde, Y_tilde) }, error = function(e) MASS::ginv(X_tilde) %*% Y_tilde)
      sigma2_hat <- sum((Y_tilde - (X_tilde %*% beta_hat))^2) / (n - p)
      qr_X <- qr(X_tilde)
      
      state <- list(error = FALSE, cpp_res = cpp_res, U = U, 
                    X_tilde = X_tilde, Y_tilde = Y_tilde, 
                    beta_hat = beta_hat, sigma2_hat = sigma2_hat, qr_X = qr_X)
    }
    
    last_params <<- params
    last_state <<- state
    
    if (runif(1) < 0.1) gc(verbose = FALSE) 
    return(state)
  }
  
  obj_REML <- function(params) {
    st <- get_state(params)
    if (st$error || st$sigma2_hat <= 1e-12) return(1e10)
    return((n - p) * log(st$sigma2_hat) + 2 * sum(log(diag(st$U))) + 2 * sum(log(abs(diag(qr.R(st$qr_X))))))
  }
  
  grad_REML <- function(params) {
    st <- get_state(params)
    if (st$error || st$sigma2_hat <= 1e-12) return(c(0, 0))
    
    W_inv <- chol2inv(st$U)
    W_inv_X <- W_inv %*% X
    
    inv_R <- tryCatch(chol2inv(qr.R(st$qr_X)), error = function(e) MASS::ginv(t(st$X_tilde) %*% st$X_tilde))
    P_mat <- W_inv - (W_inv_X %*% inv_R %*% t(W_inv_X))
    
    alpha <- W_inv %*% (Y - (X %*% st$beta_hat))
    
    neg_v_ell <- numeric(2)
    dW_list <- list(st$cpp_res$dW_dlogrho, st$cpp_res$dW_dloglam)
    
    for (k in 1:2) {
      term1 <- sum(P_mat * dW_list[[k]])
      term2 <- as.numeric(t(alpha) %*% dW_list[[k]] %*% alpha) / st$sigma2_hat
      neg_v_ell[k] <- term1 - term2
    }
    return(neg_v_ell)
  }
  
  opt_res <- optim(par = init_params, fn = obj_REML, gr = grad_REML, 
                   method = "L-BFGS-B", lower = c(-20, -20), upper = c(20, 20),
                   control = list(pgtol = 1e-5, factr = 1e7))
  
  final_st <- get_state(opt_res$par)
  
  return(list(
    kernel_params = opt_res$par,
    beta_hat = final_st$beta_hat,
    sigma2_hat = final_st$sigma2_hat,
    neg_reml_val = opt_res$value,
    convergence = opt_res$convergence
  ))
}

# 2. Fit Null Model Template (Lengthscale Fixed)
fit_unified_H0_fixed_rho_reml <- function(Y, X, Z_cov, init_log_lam, fixed_log_rho, test_type, idx_args = NULL) {
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  Z_cov <- as.matrix(Z_cov)
  n <- length(Y)
  p <- ncol(X)
  I_n <- diag(n)
  
  # --- ONE-TIME PRECOMPUTATION ---
  if (test_type == "var_imp") {
    cpp_res <- eval_unified_var_imp_fixed_rho(Z_cov, as.integer(c(idx_args) - 1), 
                                              fixed_log_rho, fixed_log_rho, 0, n)
  } else if (test_type == "interaction") {
    cpp_res <- eval_unified_interaction_fixed_rho(Z_cov, as.integer(idx_args[1] - 1), as.integer(idx_args[2] - 1), 
                                                  fixed_log_rho, fixed_log_rho, 0, n)
  } else if (test_type == "all_interactions") {
    cpp_res <- eval_unified_all_interactions_fixed_rho(Z_cov, fixed_log_rho, 0, n)
  }
  
  K_null <- cpp_res$dW0_dloglam 
  
  # --- FULL-STATE CACHE FOR NULL MODEL ---
  last_lam <- NULL
  last_state <- NULL
  
  get_state <- function(log_lam) {
    if (!is.null(last_lam) && log_lam == last_lam) return(last_state)
    
    lam_val <- exp(log_lam)
    W0 <- I_n + lam_val * K_null
    U <- try(chol(W0), silent = TRUE)
    
    if (inherits(U, "try-error")) {
      state <- list(error = TRUE)
    } else {
      X_tilde <- forwardsolve(t(U), X)
      Y_tilde <- forwardsolve(t(U), Y)
      
      beta_hat <- tryCatch({ qr.solve(X_tilde, Y_tilde) }, error = function(e) MASS::ginv(X_tilde) %*% Y_tilde)
      sigma2_hat <- sum((Y_tilde - (X_tilde %*% beta_hat))^2) / (n - p)
      qr_X <- qr(X_tilde)
      
      state <- list(error = FALSE, U = U, X_tilde = X_tilde, Y_tilde = Y_tilde, 
                    beta_hat = beta_hat, sigma2_hat = sigma2_hat, qr_X = qr_X,
                    lam_val = lam_val)
    }
    
    last_lam <<- log_lam
    last_state <<- state
    if (runif(1) < 0.1) gc(verbose = FALSE) 
    return(state)
  }
  
  obj_REML <- function(log_lam) {
    st <- get_state(log_lam)
    if (st$error || st$sigma2_hat <= 1e-12) return(1e10)
    return((n - p) * log(st$sigma2_hat) + 2 * sum(log(diag(st$U))) + 2 * sum(log(abs(diag(qr.R(st$qr_X))))))
  }
  
  grad_REML <- function(log_lam) {
    st <- get_state(log_lam)
    if (st$error || st$sigma2_hat <= 1e-12) return(0)
    
    dW0_dloglam <- st$lam_val * K_null
    
    W_inv <- chol2inv(st$U)
    W_inv_X <- W_inv %*% X
    
    inv_R <- tryCatch(chol2inv(qr.R(st$qr_X)), error = function(e) MASS::ginv(t(st$X_tilde) %*% st$X_tilde))
    P_mat <- W_inv - (W_inv_X %*% inv_R %*% t(W_inv_X))
    
    alpha <- W_inv %*% (Y - (X %*% st$beta_hat))
    
    term1 <- sum(P_mat * dW0_dloglam) 
    term2 <- as.numeric(t(alpha) %*% dW0_dloglam %*% alpha) / st$sigma2_hat
    return(term1 - term2)
  }
  
  opt_res <- optim(par = init_log_lam, fn = obj_REML, gr = grad_REML, 
                   method = "L-BFGS-B", lower = -20, upper = 20, control = list(pgtol = 1e-5, factr = 1e7))
  
  final_st <- get_state(opt_res$par)
  
  return(list(
    kernel_params = c(fixed_log_rho, opt_res$par), 
    beta_hat = final_st$beta_hat,
    sigma2_hat = final_st$sigma2_hat,
    neg_reml_val = opt_res$value,
    test_type = test_type,
    idx_args = idx_args
  ))
}

# 3. Exact Score Test Template (2x2 Information Matrix)
run_unified_score_test_reml <- function(fit_H0, Y, X, Z_cov) {
  require(MASS)
  n <- nrow(X)
  p <- ncol(X)
  I_n <- diag(n)
  
  fixed_log_rho  <- fit_H0$kernel_params[1]
  log_lam        <- fit_H0$kernel_params[2]
  sigma2         <- fit_H0$sigma2_hat
  
  lam_val <- exp(log_lam)
  
  if (fit_H0$test_type == "var_imp") {
    cpp_res <- eval_unified_var_imp_fixed_rho(Z_cov, as.integer(c(fit_H0$idx_args) - 1), 
                                              fixed_log_rho, fixed_log_rho, log_lam, n)
  } else if (fit_H0$test_type == "interaction") {
    cpp_res <- eval_unified_interaction_fixed_rho(Z_cov, as.integer(fit_H0$idx_args[1] - 1), as.integer(fit_H0$idx_args[2] - 1), 
                                                  fixed_log_rho, fixed_log_rho, log_lam, n)
  } else if (fit_H0$test_type == "all_interactions") {
    cpp_res <- eval_unified_all_interactions_fixed_rho(Z_cov, fixed_log_rho, log_lam, n)
  }
  
  U <- chol(cpp_res$W0)
  W_inv <- chol2inv(U)
  
  X_tilde <- forwardsolve(t(U), X)
  Y_tilde <- forwardsolve(t(U), Y)
  qr_X <- qr(X_tilde)
  
  W_inv_X <- W_inv %*% X
  inv_R <- tryCatch(chol2inv(qr.R(qr_X)), error = function(e) MASS::ginv(t(X_tilde) %*% X_tilde))
  P_mat <- W_inv - (W_inv_X %*% inv_R %*% t(W_inv_X))
  
  beta_hat <- tryCatch({ qr.solve(X_tilde, Y_tilde) }, error = function(e) MASS::ginv(X_tilde) %*% Y_tilde)
  alpha <- W_inv %*% (Y - (X %*% beta_hat))
  
  Q <- as.numeric((t(alpha) %*% cpp_res$K_test %*% alpha) / sigma2)
  
  tr_prod <- function(A, B) sum(A * t(B))
  
  PK <- P_mat %*% cpp_res$K_test
  PV_lam <- P_mat %*% cpp_res$dW0_dloglam
  
  tr_PK <- sum(diag(PK))
  tr_PV_lam <- sum(diag(PV_lam))
  
  I_nuis <- matrix(0, 2, 2)
  I_nuis[1, 1] <- (n - p) / 2
  I_nuis[1, 2] <- 0.5 * tr_PV_lam
  I_nuis[2, 1] <- I_nuis[1, 2]
  I_nuis[2, 2] <- 0.5 * tr_prod(PV_lam, PV_lam)
  
  I_cross <- matrix(0, 2, 1)
  I_cross[1, 1] <- 0.5 * tr_PK
  I_cross[2, 1] <- 0.5 * tr_prod(PV_lam, PK) 
  
  I_test <- 0.5 * tr_prod(PK, PK)            
  
  I_nuis_inv <- tryCatch(solve(I_nuis), error = function(e) MASS::ginv(I_nuis))
  V_schur <- as.numeric(I_test - (t(I_cross) %*% I_nuis_inv %*% I_cross))
  
  if (V_schur <= 1e-12 || tr_PK <= 1e-12) {
    return(list(Q_stat = Q, p_value = 1.0, Satterthwaite_Scale = NA, Satterthwaite_DF = NA))
  }
  
  mu_Q <- tr_PK 
  sigma2_Q <- 4 * V_schur 
  
  scale_a <- sigma2_Q / (2 * mu_Q)
  df_d    <- (2 * mu_Q^2) / sigma2_Q
  
  p_val <- pchisq(Q / scale_a, df = df_d, lower.tail = FALSE)
  
  return(list(Q_stat = Q, p_value = p_val, Satterthwaite_Scale = scale_a, Satterthwaite_DF = df_d))
}

# 4. Master Pipeline Wrapper
run_unified_pipeline_reml <- function(Y, X, Z_cov, init_params_full, 
                                      test_type = "var_imp", idx_args = NULL, 
                                      global_rho = NULL) {
  
  if (is.null(global_rho)) {
    fit_full <- fit_unified_full_model_reml(Y, X, Z_cov, init_params = init_params_full)
    if (is.null(fit_full) || fit_full$convergence != 0) return(NULL) 
    
    fixed_log_rho  <- fit_full$kernel_params[1]
    init_log_lam   <- fit_full$kernel_params[2]
  } else {
    fixed_log_rho  <- global_rho
    init_log_lam   <- init_params_full[2]
  }
  
  fit_H0 <- fit_unified_H0_fixed_rho_reml(Y, X, Z_cov, init_log_lam, fixed_log_rho, test_type, idx_args)
  score_res <- run_unified_score_test_reml(fit_H0, Y, X, Z_cov)
  
  idx_args_str <- if(is.null(idx_args)) "" else paste(idx_args, collapse=",")
  
  return(list(Test_Type = test_type, Index_Args = idx_args_str,
              Fixed_Log_Rho = fixed_log_rho, H0_Log_Lam = fit_H0$kernel_params[2],
              Q_stat = score_res$Q_stat, p_value = score_res$p_value,
              Satterthwaite_Scale = score_res$Satterthwaite_Scale, Satterthwaite_DF = score_res$Satterthwaite_DF))
}

# 5. Backward Compatibility Wrappers
run_unified_var_imp_pipeline_reml <- function(Y, X, Z_cov, init_params_full, test_idx_vec, global_log_rho=NULL) {
  run_unified_pipeline_reml(Y, X, Z_cov, init_params_full, test_type="var_imp", idx_args=test_idx_vec, global_rho=global_log_rho)
}
run_unified_interaction_pipeline_reml <- function(Y, X, Z_cov, init_params_full, idx1, idx2, global_log_rho=NULL) {
  run_unified_pipeline_reml(Y, X, Z_cov, init_params_full, test_type="interaction", idx_args=c(idx1, idx2), global_rho=global_log_rho)
}
run_unified_all_interactions_pipeline_reml <- function(Y, X, Z_cov, init_params_full, global_log_rho=NULL) {
  run_unified_pipeline_reml(Y, X, Z_cov, init_params_full, test_type="all_interactions", idx_args=NULL, global_rho=global_log_rho)
}


# =========================================================
# 6. High-Level Data Analysis Pipelines (Sequential & Parallel)
# =========================================================

data_analysis_pipeline <- function(Y, X, Z_cov, init_params_full=c(0,0), alpha=.05) {
  require(PairwiseSims)
  Z_cov <- as.matrix(Z_cov)
  p <- ncol(Z_cov)
  
  fit_full <- fit_unified_full_model_reml(Y, X, Z_cov, init_params = init_params_full)
  
  if (is.null(fit_full) || fit_full$convergence != 0) {
    warning("Global unified REML model failed to converge. Proceeding with default init.")
    optimal_log_rho <- init_params_full[1]
  } else {
    optimal_log_rho <- fit_full$kernel_params[1]
  }
  
  var.importance.pvals <- rep(NA, p)
  names(var.importance.pvals) <- if(!is.null(colnames(Z_cov))) colnames(Z_cov) else paste0("Var", 1:p)
  
  for (j in 1:p) {
    score_results <- tryCatch({
      run_unified_pipeline_reml(Y, X, Z_cov, init_params_full = init_params_full,
                                test_type = "var_imp", idx_args = j, global_rho = optimal_log_rho)
    }, error = function(e) {
      message("Error on Var ", j, ": ", e$message)
      return(NULL)
    })
    
    if (!is.null(score_results)) {
      var.importance.pvals[j] <- score_results$p_value
    }
  }
  
  sig.pair <- which(var.importance.pvals < alpha)
  
  if (length(sig.pair) > 1) {
    interaction.tests <- t(combn(sig.pair, 2))
    interaction.tests <- cbind(interaction.tests, rep(NA, nrow(interaction.tests)))
    row.names(interaction.tests) <- 1:nrow(interaction.tests)
    colnames(interaction.tests) <- c("Var1", "Var2", "p-value")
    
    var_names <- names(var.importance.pvals)
    for (i in 1:nrow(interaction.tests)) {
      row.names(interaction.tests)[i] <- paste(var_names[interaction.tests[i, 1]],
                                               var_names[interaction.tests[i, 2]], sep=" x ")
    }
    
    for (i in 1:nrow(interaction.tests)) {
      score_res <- tryCatch({
        run_unified_pipeline_reml(Y, X, Z_cov, init_params_full = init_params_full,
                                  test_type = "interaction", 
                                  idx_args = c(interaction.tests[i, 1], interaction.tests[i, 2]), 
                                  global_rho = optimal_log_rho)
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

data_analysis_pipeline_parallel <- function(Y, X, Z_cov, init_params_full=c(0,0), alpha=.05) {
  require(foreach) 
  require(PairwiseSims)
  
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  Z_cov <- as.matrix(Z_cov)
  p <- ncol(Z_cov)
  
  fit_full <- fit_unified_full_model_reml(Y, X, Z_cov, init_params = init_params_full)
  
  if (is.null(fit_full) || fit_full$convergence != 0) {
    warning("Global unified REML model failed to converge. Proceeding with default init.")
    optimal_log_rho <- init_params_full[1]
  } else {
    optimal_log_rho <- fit_full$kernel_params[1]
  }
  
  funcs_to_export <- c("fit_unified_full_model_reml",
                       "fit_unified_H0_fixed_rho_reml",
                       "run_unified_score_test_reml",
                       "run_unified_pipeline_reml")
  
  var.importance.pvals <- foreach(j = 1:p, .combine = c, 
                                  .packages = c("MASS", "PairwiseSims"),
                                  .export = funcs_to_export) %dopar% {
                                    score_results <- tryCatch({
                                      run_unified_pipeline_reml(Y, X, Z_cov, init_params_full = init_params_full,
                                                                test_type = "var_imp", idx_args = j, global_rho = optimal_log_rho)
                                    }, error = function(e) {
                                      message("Worker Error on Var ", j, ": ", e$message)
                                      return(NULL)
                                    })
                                    if (!is.null(score_results)) score_results$p_value else NA
                                  }
  
  names(var.importance.pvals) <- if(!is.null(colnames(Z_cov))) colnames(Z_cov) else paste0("Var", 1:p)
  
  sig.pair <- which(var.importance.pvals < alpha)
  
  if (length(sig.pair) > 1) {
    interaction.tests <- t(combn(sig.pair, 2))
    var_names <- names(var.importance.pvals)
    test_names <- paste(var_names[interaction.tests[, 1]], var_names[interaction.tests[, 2]], sep=" x ")
    
    int_pvals <- foreach(i = 1:nrow(interaction.tests), .combine = c, 
                         .packages = c("MASS", "PairwiseSims"), .export = funcs_to_export) %dopar% {
                           score_res <- tryCatch({
                             run_unified_pipeline_reml(Y, X, Z_cov, init_params_full = init_params_full,
                                                       test_type = "interaction", 
                                                       idx_args = c(interaction.tests[i, 1], interaction.tests[i, 2]), 
                                                       global_rho = optimal_log_rho)
                           }, error = function(e) {
                             message("Worker Error on Interaction ", i, ": ", e$message)
                             return(NULL)
                           })
                           if (!is.null(score_res)) score_res$p_value else NA
                         }
    
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

# =========================================================
# 7. Core Prediction Engine and Plot Wrappers
# =========================================================

predict_unified_surface_with_se <- function(fit_full, Y_train, X_train, Z_train, Z_new, X_new) {
  Z_train <- as.matrix(Z_train)
  Z_new <- as.matrix(Z_new)
  X_train <- as.matrix(X_train)
  X_new <- as.matrix(X_new)
  
  n <- nrow(Z_train)
  p_fixed <- ncol(X_train)
  
  log_rho <- fit_full$kernel_params[1]
  log_lam <- fit_full$kernel_params[2]
  lambda <- exp(log_lam)
  
  beta_hat <- fit_full$beta_hat
  sigma2_hat <- fit_full$sigma2_hat
  
  cpp_train <- eval_unified_full_model(Z_train, log_rho, log_lam, n)
  U <- try(chol(cpp_train$W), silent = TRUE)
  if (inherits(U, "try-error")) stop("W matrix computationally singular.")
  W_inv <- chol2inv(U)
  
  alpha <- W_inv %*% (Y_train - (X_train %*% beta_hat))
  
  cpp_cross <- eval_unified_K_cross_and_diag(Z_new, Z_train, log_rho)
  K_cross <- cpp_cross$K_cross
  K_diag  <- as.numeric(cpp_cross$K_diag)
  
  Y_pred <- (X_new %*% beta_hat) + (lambda * K_cross %*% alpha)
  
  scaled_K <- lambda * K_cross
  var_reduction <- rowSums((scaled_K %*% W_inv) * scaled_K)
  
  X_tilde <- forwardsolve(t(U), X_train)
  cov_beta <- chol2inv(qr.R(qr(X_tilde))) 
  
  R_mat <- X_new - (scaled_K %*% W_inv %*% X_train)
  var_fixed <- rowSums((R_mat %*% cov_beta) * R_mat)
  
  se2 <- sigma2_hat * (lambda * K_diag - var_reduction + var_fixed)
  se <- sqrt(pmax(se2, 1e-10)) 
  t_crit <- qt(0.975, df = n - p_fixed)
  
  return(data.frame(
    Predicted_Y = as.numeric(Y_pred),
    SE = as.numeric(se),
    CI_Lower = as.numeric(Y_pred - t_crit * se),
    CI_Upper = as.numeric(Y_pred + t_crit * se)
  ))
}

library(ggplot2)

plot_main_effect_slice <- function(target_idx, fit_full, Y, X, Z_train, 
                                   orig_min = 0, orig_max = 1, unlog = FALSE, 
                                   grid_points = 100) {
  Z_train <- as.matrix(Z_train)
  X <- as.matrix(X)
  p <- ncol(Z_train)
  
  var_name <- if(!is.null(colnames(Z_train))) colnames(Z_train)[target_idx] else paste("Variable", target_idx)
  
  Z_new <- matrix(apply(Z_train, 2, median), nrow = grid_points, ncol = p, byrow = TRUE)
  x_seq <- seq(min(Z_train[, target_idx]), max(Z_train[, target_idx]), length.out = grid_points)
  Z_new[, target_idx] <- x_seq
  
  X_new <- matrix(apply(X, 2, median), nrow = grid_points, ncol = ncol(X), byrow = TRUE)
  if (all(X[, 1] == 1)) X_new[, 1] <- 1 
  
  pred_df <- predict_unified_surface_with_se(fit_full, Y, X, Z_train, Z_new, X_new)
  
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

plot_interaction_slice <- function(idx_main, idx_quant, fit_full, Y, X, Z_train, 
                                   orig_min_main = 0, orig_max_main = 1, unlog_main = FALSE,
                                   orig_min_quant = 0, orig_max_quant = 1, unlog_quant = FALSE,
                                   grid_points = 100) {
  Z_train <- as.matrix(Z_train)
  X <- as.matrix(X)
  p <- ncol(Z_train)
  
  name_main <- if(!is.null(colnames(Z_train))) colnames(Z_train)[idx_main] else paste("Variable", idx_main)
  name_quant <- if(!is.null(colnames(Z_train))) colnames(Z_train)[idx_quant] else paste("Variable", idx_quant)
  
  quantiles <- quantile(Z_train[, idx_quant], probs = c(0.10, 0.50, 0.90))
  
  quant_real <- (quantiles * (orig_max_quant - orig_min_quant)) + orig_min_quant
  if (unlog_quant) quant_real <- exp(quant_real)
  
  quant_labels <- paste0(c("10th (", "50th (", "90th ("), round(quant_real, 2), ")")
  
  plot_data <- data.frame()
  
  for (i in seq_along(quantiles)) {
    q_val <- quantiles[i]
    q_label <- quant_labels[i]
    
    Z_new <- matrix(apply(Z_train, 2, median), nrow = grid_points, ncol = p, byrow = TRUE)
    x_seq <- seq(min(Z_train[, idx_main]), max(Z_train[, idx_main]), length.out = grid_points)
    
    Z_new[, idx_main] <- x_seq
    Z_new[, idx_quant] <- q_val
    
    X_new <- matrix(apply(X, 2, median), nrow = grid_points, ncol = ncol(X), byrow = TRUE)
    if (all(X[, 1] == 1)) X_new[, 1] <- 1
    
    pred_df <- predict_unified_surface_with_se(fit_full, Y, X, Z_train, Z_new, X_new)
    
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


# =========================================================
# Backward Compatibility Wrappers (Append to Section 5)
# =========================================================

# --- 1. Null Model Fitting Wrappers ---
fit_unified_var_imp_H0_fixed_rho_reml <- function(Y, X, Z_cov, init_log_lam, test_idx_vec, fixed_log_rho) {
  fit_unified_H0_fixed_rho_reml(Y, X, Z_cov, init_log_lam, fixed_log_rho, 
                                test_type = "var_imp", idx_args = test_idx_vec)
}

fit_unified_interaction_H0_fixed_rho_reml <- function(Y, X, Z_cov, init_log_lam, idx1, idx2, fixed_log_rho) {
  fit_unified_H0_fixed_rho_reml(Y, X, Z_cov, init_log_lam, fixed_log_rho, 
                                test_type = "interaction", idx_args = c(idx1, idx2))
}

fit_unified_all_interactions_H0_fixed_rho_reml <- function(Y, X, Z_cov, init_log_lam, fixed_log_rho) {
  fit_unified_H0_fixed_rho_reml(Y, X, Z_cov, init_log_lam, fixed_log_rho, 
                                test_type = "all_interactions", idx_args = NULL)
}

# --- 2. Score Test Evaluation Wrappers ---
run_unified_var_imp_score_test_reml <- function(fit_H0, Y, X, Z_cov) {
  run_unified_score_test_reml(fit_H0, Y, X, Z_cov)
}

run_unified_interaction_score_test_reml <- function(fit_H0, Y, X, Z_cov) {
  run_unified_score_test_reml(fit_H0, Y, X, Z_cov)
}

run_unified_all_interactions_score_test_reml <- function(fit_H0, Y, X, Z_cov) {
  run_unified_score_test_reml(fit_H0, Y, X, Z_cov)
}