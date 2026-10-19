# ---------------------------------------------------------
# FILE: code/simulation_global_interaction.R
# Standalone Simulation: Test for ANY Interaction
# ---------------------------------------------------------
rm(list=ls())

# 1. Define Environment (1 for Cluster, 0 for Local)
hpc <- 1 

# 2. Source Master Framework (Automatically handles cluster, libraries, and sub-scripts)
if(hpc) {
  source("code/simulation_framework.R")
} else {
  source("~/GitHub/Pairwise_Additive_GP_Tests/code/final_simulations/simulation_framework.R")
}

# =========================================================
# CORE EVALUATION WRAPPER (Runs one replicate)
# =========================================================
evaluate_global_interaction <- function(n, p, rho_corr, mean_fn, snr) {
  
  # 1. Generate Correlated Predictors via Gaussian Copula
  Sigma <- matrix(rho_corr, nrow = p, ncol = p)
  diag(Sigma) <- 1
  Z <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  S <- pnorm(Z) 
  colnames(S) <- paste0("S", 1:p)
  
  D_list <- lapply(1:p, function(j) as.matrix(dist(S[,j])))
  coords_list <- list(D_1D_list = D_list)
  
  X_int <- matrix(1, nrow = n) 
  
  # 2. Generate Signal and Apply Exact SNR (Training)
  EY <- mean_fn(S)
  var_signal <- var(EY)
  var_noise <- var_signal / snr
  
  EY_train_raw <- EY + rnorm(n, sd = sqrt(var_noise))
  mean_Y <- mean(EY_train_raw)
  sd_Y <- sd(EY_train_raw)
  
  Y <- as.numeric((EY_train_raw - mean_Y) / sd_Y)
  
  # 3. Generate Held-Out Test Data for MSE (Strictly No Noise)
  Z_test <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  S_test <- pnorm(Z_test)
  colnames(S_test) <- paste0("S", 1:p)
  X_int_test <- matrix(1, nrow = n)
  
  # Scale test true signal exactly like Y for valid MSE measurement
  EY_test_raw <- mean_fn(S_test)
  EY_test <- as.numeric((EY_test_raw - mean_Y) / sd_Y)
  
  # Initialize output using standardized names for the report generator
  res <- c(
    Unified_REML_pval = NA,
    Unified_REML_MSE  = NA
  )
  
  # 4. Execute the Global Interaction Test & Prediction Pipeline
  # Mirrors the var_imp/interaction tests in simulation_framework.R: fit the full
  # (alternative) model once to estimate rho, use it for prediction, then run the
  # all-interactions score test with rho fixed at that estimate.
  tryCatch({
    fit_full <- fit_unified_full_model_reml(Y, X_int, Z_cov = S, init_params = c(0, 0))
    if (!is.null(fit_full) && fit_full$convergence == 0) {
      # Prediction accuracy (MSE) of the full model
      pred_uni <- predict_unified_surface_with_se(fit_full, Y_train = Y, X_train = X_int,
                                                  Z_train = S, Z_new = S_test, X_new = X_int_test)
      res["Unified_REML_MSE"] <- as.numeric(mean((pred_uni$Predicted_Y - EY_test)^2))

      # Score test for any pairwise interaction
      test_out <- run_unified_pipeline_reml(
        Y, X_int, Z_cov = S, init_params_full = fit_full$kernel_params,
        test_type = "all_interactions", idx_args = NULL,
        global_rho = fit_full$kernel_params[1]
      )
      if (!is.null(test_out)) {
        res["Unified_REML_pval"] <- as.numeric(test_out$p_value)
      }
    }

  }, error = function(e) {
    warning(paste("Worker Global Interaction Test Error:", e$message))
  })
  
  return(res)
}

# =========================================================
# MASTER SIMULATION LOOP
# =========================================================
run_global_interaction_simulation <- function(B = 1000, n = 300, p = 10, 
                                              rho_corr = 0.0, snr = 1.0, 
                                              mean_fn, is_null_sim = FALSE) {
  
  cat(sprintf("\nStarting Global Interaction Simulation: B=%d, n=%d, p=%d, rho=%.2f, SNR=%.2f\n", 
              B, n, p, rho_corr, snr))
  cat(sprintf("Condition: %s\n", ifelse(is_null_sim, "NULL (Type I Error)", "ALTERNATIVE (Power)")))
  
  # Export necessary functions for global test (plus framework predictors)
  explicit_exports <- c(
    "evaluate_global_interaction",
    "fit_unified_full_model_reml",
    "fit_unified_H0_fixed_rho_reml",
    "run_unified_score_test_reml",
    "run_unified_pipeline_reml",
    "fit_unified_all_interactions_H0_fixed_rho_reml",
    "run_unified_all_interactions_score_test_reml",
    "run_unified_all_interactions_pipeline_reml",
    "eval_unified_full_model",
    "eval_unified_all_interactions_fixed_rho",
    "predict_unified_surface_with_se"
  )
  funcs_to_export <- unique(c(explicit_exports, as.character(lsf.str())))
  
  results_matrix <- foreach(b = 1:B, 
                            .combine = rbind, 
                            .packages = c("MASS", "PairwiseSims", "mgcv", "bkmr", "fields"), 
                            .export = funcs_to_export) %dopar% {
                              
                              evaluate_global_interaction(n = n, p = p, 
                                                          rho_corr = rho_corr, 
                                                          mean_fn = mean_fn, 
                                                          snr = snr)
                            }
  
  rownames(results_matrix) <- NULL
  return(as.data.frame(results_matrix))
}


# =========================================================
# EXECUTION SCRIPT
# =========================================================

# 1. DEFINE FUNCTIONAL FORMS (Matched to Pure Interaction)
# Null: S1 has no effect (Interaction term vanishes)
mean_fn_null <- function(S) {
  sin(2 * pi * S[,3]) + sin(2 * pi * S[,4])
}

# Alternative: S1 is present via pure interaction with S2
mean_fn_alt <- function(S) {
  sin(2 * pi * S[,3]) + sin(2 * pi * S[,4]) + 15 * (S[,1] - 0.5) * (S[,2] - 0.5)
}

B_reps <- 2000
n_val <- 100
p_vars <- 10
snr_val <- 1.0
rho_val <- 0.0

set.seed(123)

cat("\n--- Running Type I Error Evaluation ---\n")
results_null <- run_global_interaction_simulation(
  B = B_reps, n = n_val, p = p_vars, rho_corr = rho_val, snr = snr_val, 
  mean_fn = mean_fn_null, is_null_sim = TRUE
)

cat("\n--- Running Power Evaluation ---\n")
results_alt <- run_global_interaction_simulation(
  B = B_reps, n = n_val, p = p_vars, rho_corr = rho_val, snr = snr_val, 
  mean_fn = mean_fn_alt, is_null_sim = FALSE
)

# Clean up cluster (initialized in framework)
stopCluster(cl)

# =========================================================
# SUMMARIZE AND SAVE OUTPUTS
# =========================================================
alpha_thresh <- 0.10

payload <- list(
  Type_I_Error = results_null,
  Power        = results_alt
)

# Fire the agnostic summarizer from the framework
study_results <- generate_study_report(
  results_list = payload, 
  alpha = alpha_thresh, 
  save_path = "outputs/SimResults_AnyInteraction.Rdata" # file name read by output_table.R
)