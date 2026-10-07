# =========================================================
# FILE: ss_anova_reml_unified_sp.R
# SS-ANOVA REML Approach with Full-State Caching
# =========================================================

# 1. Fit Global Unified Model (4 Parameters)
fit_unified_full_model_reml_sp <- function(Y, X, Z_cov, coords_sp, init_params) {
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  Z_cov <- as.matrix(Z_cov)
  coords_sp <- as.matrix(coords_sp)
  n <- length(Y)
  p <- ncol(X)
  
  # --- FULL-STATE CACHE ---
  # Stores C++ matrices and the expensive O(N^3) Cholesky decompositions
  last_params <- NULL
  last_state <- NULL
  
  get_state <- function(params) {
    if (!is.null(last_params) && all(params == last_params)) return(last_state)
    
    cpp_res <- eval_unified_full_model_sp(Z_cov, coords_sp, params[1], params[2], params[3], params[4], n)
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
    
    # Keep memory clean during heavy optimization
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
    if (st$error || st$sigma2_hat <= 1e-12) return(c(0, 0, 0, 0))
    
    W_inv <- chol2inv(st$U)
    W_inv_X <- W_inv %*% X
    
    inv_R <- tryCatch(chol2inv(qr.R(st$qr_X)), error = function(e) MASS::ginv(t(st$X_tilde) %*% st$X_tilde))
    P_mat <- W_inv - (W_inv_X %*% inv_R %*% t(W_inv_X))
    
    alpha <- W_inv %*% (Y - (X %*% st$beta_hat))
    
    neg_v_ell <- numeric(4)
    dW_list <- list(st$cpp_res$dW_dlogrho_a, st$cpp_res$dW_dloglam_a, 
                    st$cpp_res$dW_dlogrho_sp, st$cpp_res$dW_dloglam_sp)
    
    for (k in 1:4) {
      term1 <- sum(P_mat * dW_list[[k]])
      term2 <- as.numeric(t(alpha) %*% dW_list[[k]] %*% alpha) / st$sigma2_hat
      neg_v_ell[k] <- term1 - term2
    }
    return(neg_v_ell)
  }
  
  opt_res <- optim(par = init_params, fn = obj_REML, gr = grad_REML, 
                   method = "L-BFGS-B", lower = rep(-20, 4), upper = rep(20, 4),
                   control = list(pgtol = 1e-5, factr = 1e7))
  
  # Retrieve final state directly from cache (0 compute cost)
  final_st <- get_state(opt_res$par)
  
  return(list(
    kernel_params = opt_res$par,
    beta_hat = final_st$beta_hat,
    sigma2_hat = final_st$sigma2_hat,
    neg_reml_val = opt_res$value,
    convergence = opt_res$convergence
  ))
}

# 2. Fit Null Model Template (Lengthscales Fixed)
fit_unified_H0_fixed_rho_sp_reml <- function(Y, X, Z_cov, coords_sp, init_lams, 
                                             fixed_log_rho_a, fixed_log_rho_sp, test_type, idx_args = NULL) {
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  Z_cov <- as.matrix(Z_cov)
  coords_sp <- as.matrix(coords_sp)
  n <- length(Y)
  p <- ncol(X)
  I_n <- diag(n)
  
  # --- ONE-TIME PRECOMPUTATION ---
  if (test_type == "var_imp") {
    cpp_res <- eval_unified_var_imp_fixed_rho_sp(Z_cov, coords_sp, 
                                                 as.integer(idx_args - 1), fixed_log_rho_a, fixed_log_rho_sp, n)
  } else if (test_type == "interaction") {
    cpp_res <- eval_unified_interaction_fixed_rho_sp(Z_cov, coords_sp, 
                                                     as.integer(idx_args[1] - 1), as.integer(idx_args[2] - 1), 
                                                     fixed_log_rho_a, fixed_log_rho_sp, n)
  } else if (test_type == "all_interactions") {
    cpp_res <- eval_unified_all_interactions_fixed_rho_sp(Z_cov, coords_sp, 
                                                          fixed_log_rho_a, fixed_log_rho_sp, n)
  }
  
  K_null_a <- cpp_res$K_null_a
  K_sp <- cpp_res$K_sp
  
  # --- FULL-STATE CACHE FOR NULL MODEL ---
  last_lams <- NULL
  last_state <- NULL
  
  get_state <- function(log_lams) {
    if (!is.null(last_lams) && all(log_lams == last_lams)) return(last_state)
    
    lam_a <- exp(log_lams[1])
    lam_sp <- exp(log_lams[2])
    W0 <- I_n + lam_a * K_null_a + lam_sp * K_sp
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
                    lam_a = lam_a, lam_sp = lam_sp)
    }
    
    last_lams <<- log_lams
    last_state <<- state
    if (runif(1) < 0.1) gc(verbose = FALSE) 
    return(state)
  }
  
  obj_REML <- function(log_lams) {
    st <- get_state(log_lams)
    if (st$error || st$sigma2_hat <= 1e-12) return(1e10)
    return((n - p) * log(st$sigma2_hat) + 2 * sum(log(diag(st$U))) + 2 * sum(log(abs(diag(qr.R(st$qr_X))))))
  }
  
  grad_REML <- function(log_lams) {
    st <- get_state(log_lams)
    if (st$error || st$sigma2_hat <= 1e-12) return(c(0, 0))
    
    dW0_dloglam_a <- st$lam_a * K_null_a
    dW0_dloglam_sp <- st$lam_sp * K_sp
    
    W_inv <- chol2inv(st$U)
    W_inv_X <- W_inv %*% X
    
    inv_R <- tryCatch(chol2inv(qr.R(st$qr_X)), error = function(e) MASS::ginv(t(st$X_tilde) %*% st$X_tilde))
    P_mat <- W_inv - (W_inv_X %*% inv_R %*% t(W_inv_X))
    
    alpha <- W_inv %*% (Y - (X %*% st$beta_hat))
    
    neg_v_ell <- numeric(2)
    dW_list <- list(dW0_dloglam_a, dW0_dloglam_sp)
    
    for (k in 1:2) {
      term1 <- sum(P_mat * dW_list[[k]])
      term2 <- as.numeric(t(alpha) %*% dW_list[[k]] %*% alpha) / st$sigma2_hat
      neg_v_ell[k] <- term1 - term2
    }
    return(neg_v_ell)
  }
  
  opt_res <- optim(par = init_lams, fn = obj_REML, gr = grad_REML, 
                   method = "L-BFGS-B", lower = c(-20, -20), upper = c(20, 20), 
                   control = list(pgtol = 1e-5, factr = 1e7))
  
  final_st <- get_state(opt_res$par)
  
  return(list(
    kernel_params = c(fixed_log_rho_a, opt_res$par[1], fixed_log_rho_sp, opt_res$par[2]), 
    beta_hat = final_st$beta_hat,
    sigma2_hat = final_st$sigma2_hat,
    neg_reml_val = opt_res$value,
    test_type = test_type,
    idx_args = idx_args
  ))
}

# 3. Exact Score Test Template (3x3 Information Matrix)
run_unified_score_test_reml_sp <- function(fit_H0, Y, X, Z_cov, coords_sp) {
  require(MASS)
  n <- nrow(X)
  p <- ncol(X)
  I_n <- diag(n)
  
  fixed_log_rho_a  <- fit_H0$kernel_params[1]
  log_lam_a        <- fit_H0$kernel_params[2]
  fixed_log_rho_sp <- fit_H0$kernel_params[3]
  log_lam_sp       <- fit_H0$kernel_params[4]
  sigma2           <- fit_H0$sigma2_hat
  
  lam_a <- exp(log_lam_a)
  lam_sp <- exp(log_lam_sp)
  
  # Fetch test space geometry
  if (fit_H0$test_type == "var_imp") {
    cpp_res <- eval_unified_var_imp_fixed_rho_sp(Z_cov, coords_sp, 
                                                 as.integer(fit_H0$idx_args - 1), fixed_log_rho_a, fixed_log_rho_sp, n)
  } else if (fit_H0$test_type == "interaction") {
    cpp_res <- eval_unified_interaction_fixed_rho_sp(Z_cov, coords_sp, 
                                                     as.integer(fit_H0$idx_args[1] - 1), as.integer(fit_H0$idx_args[2] - 1), 
                                                     fixed_log_rho_a, fixed_log_rho_sp, n)
  } else if (fit_H0$test_type == "all_interactions") {
    cpp_res <- eval_unified_all_interactions_fixed_rho_sp(Z_cov, coords_sp, 
                                                          fixed_log_rho_a, fixed_log_rho_sp, n)
  }
  
  W0 <- I_n + lam_a * cpp_res$K_null_a + lam_sp * cpp_res$K_sp
  U <- chol(W0)
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
  
  dW0_dloglam_a <- lam_a * cpp_res$K_null_a
  dW0_dloglam_sp <- lam_sp * cpp_res$K_sp
  
  PK <- P_mat %*% cpp_res$K_test
  PV_lam_a <- P_mat %*% dW0_dloglam_a
  PV_lam_sp <- P_mat %*% dW0_dloglam_sp
  
  tr_PK <- sum(diag(PK))
  tr_PV_lam_a <- sum(diag(PV_lam_a))
  tr_PV_lam_sp <- sum(diag(PV_lam_sp))
  
  I_nuis <- matrix(0, 3, 3)
  I_nuis[1, 1] <- (n - p) / 2
  I_nuis[1, 2] <- 0.5 * tr_PV_lam_a
  I_nuis[1, 3] <- 0.5 * tr_PV_lam_sp
  
  I_nuis[2, 1] <- I_nuis[1, 2]
  I_nuis[2, 2] <- 0.5 * tr_prod(PV_lam_a, PV_lam_a)
  I_nuis[2, 3] <- 0.5 * tr_prod(PV_lam_a, PV_lam_sp)
  
  I_nuis[3, 1] <- I_nuis[1, 3]
  I_nuis[3, 2] <- I_nuis[2, 3]
  I_nuis[3, 3] <- 0.5 * tr_prod(PV_lam_sp, PV_lam_sp)
  
  I_cross <- matrix(0, 3, 1)
  I_cross[1, 1] <- 0.5 * tr_PK
  I_cross[2, 1] <- 0.5 * tr_prod(PV_lam_a, PK)
  I_cross[3, 1] <- 0.5 * tr_prod(PV_lam_sp, PK) 
  
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
run_unified_pipeline_reml_sp <- function(Y, X, Z_cov, coords_sp, init_params_full, 
                                         test_type = "var_imp", idx_args = NULL, 
                                         global_rhos = NULL) {
  
  if (is.null(global_rhos)) {
    fit_full <- fit_unified_full_model_reml_sp(Y, X, Z_cov, coords_sp, init_params = init_params_full)
    if (is.null(fit_full) || fit_full$convergence != 0) return(NULL) 
    
    fixed_log_rho_a  <- fit_full$kernel_params[1]
    init_log_lam_a   <- fit_full$kernel_params[2]
    fixed_log_rho_sp <- fit_full$kernel_params[3]
    init_log_lam_sp  <- fit_full$kernel_params[4]
  } else {
    fixed_log_rho_a  <- global_rhos[1]
    fixed_log_rho_sp <- global_rhos[2]
    init_log_lam_a   <- init_params_full[2]
    init_log_lam_sp  <- init_params_full[4]
  }
  
  fit_H0 <- fit_unified_H0_fixed_rho_sp_reml(Y, X, Z_cov, coords_sp, 
                                             init_lams = c(init_log_lam_a, init_log_lam_sp), 
                                             fixed_log_rho_a, fixed_log_rho_sp, test_type, idx_args)
  
  score_res <- run_unified_score_test_reml_sp(fit_H0, Y, X, Z_cov, coords_sp)
  
  # Ensure idx_args isn't dropped into nothingness when returning NULL
  idx_args_str <- if(is.null(idx_args)) "" else paste(idx_args, collapse=",")
  
  return(list(Test_Type = test_type, Index_Args = idx_args_str,
              Fixed_Log_Rho_A = fixed_log_rho_a, Fixed_Log_Rho_Sp = fixed_log_rho_sp, 
              H0_Log_Lam_A = fit_H0$kernel_params[2], H0_Log_Lam_Sp = fit_H0$kernel_params[4],
              Q_stat = score_res$Q_stat, p_value = score_res$p_value,
              Satterthwaite_Scale = score_res$Satterthwaite_Scale, Satterthwaite_DF = score_res$Satterthwaite_DF))
}