# ---------------------------------------------------------
# FILE: code/simulation_framework.R
# Master Framework & Parallelization Setup
# ---------------------------------------------------------
library(MASS)
library(mgcv)
library(bkmr)
library(foreach)
library(doSNOW)
library(PairwiseSims)

# =========================================================
# GLOBAL SETUP & CLUSTER INITIALIZATION
# =========================================================
if (!exists("hpc")) hpc <- 1

if (hpc) {
  source("code/ss_anova_reml_unified.R")
  source("code/KMR.R")
  source("code/KMR_pairwise.R") 
  
  n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
  if (is.na(n_cores)) {
    n_cores <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE", unset = "1"))
  }
} else {
  n_cores <- parallel::detectCores() - 1 
  setwd("~/GitHub/Pairwise_Additive_GP_Tests")
  source("code/KMR/ss_anova_reml_unified.R")
  source("code/KMR/KMR.R")
  source("code/KMR/KMR_pairwise.R")
}

cl <- makeCluster(n_cores)
registerDoSNOW(cl)

cat(paste0("\n[System Setup] Initialized cluster with ", n_cores, " cores. HPC Mode: ", as.logical(hpc), "\n"))

# ---------------------------------------------------------
# CORE EVALUATION WRAPPER (Runs one replicate)
# ---------------------------------------------------------
evaluate_all_methods <- function(n, p, rho_corr, mean_fn, snr, run_bkmr = FALSE, target_idx = 1, test_type = "main") {
  
  # 1. Generate Correlated Predictors via Gaussian Copula
  Sigma <- matrix(rho_corr, nrow = p, ncol = p)
  diag(Sigma) <- 1
  Z <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  S <- pnorm(Z) 
  colnames(S) <- paste0("S", 1:p)
  
  # Kept for backward compatibility with older methods
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
  S_scaled <- scale(S)
  
  # 3. Generate Held-Out Test Data for MSE (Strictly No Noise)
  Z_test <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  S_test <- pnorm(Z_test)
  colnames(S_test) <- paste0("S", 1:p)
  X_int_test <- matrix(1, nrow = n)
  
  EY_test_raw <- mean_fn(S_test)
  EY_test <- as.numeric((EY_test_raw - mean_Y) / sd_Y)
  S_test_scaled <- scale(S_test, center = attr(S_scaled, "scaled:center"), scale = attr(S_scaled, "scaled:scale"))
  
  # --- CALCULATE GLOBAL RHO ANCHOR ---
  fit_global <- tryCatch({
    fit_unified_full_model_reml(Y, X_int, Z_cov = S, init_params = c(0, 0))
  }, error = function(e) NULL)
  
  global_log_rho <- if (!is.null(fit_global) && fit_global$convergence == 0) fit_global$kernel_params[1] else 0
  
  res <- c(
    Unified_REML_pval       = NA,
    Unified_REML_MSE        = NA,
    Full_GP_pval            = NA,
    Full_GP_MSE             = NA,
    Old_Pairwise_GP_pval    = NA,
    Old_Pairwise_GP_MSE     = NA,
    GAM_Smooth_pval         = NA,
    GAM_Smooth_MSE          = NA,
    OLS_pval                = NA,
    OLS_MSE                 = NA,
    BKMR_PIP                = NA,
    BKMR_MSE                = NA
  )
  
  # --- Method 1: Unified REML ---
  tryCatch({
    fit_full <- fit_unified_full_model_reml(Y, X_int, Z_cov = S, init_params = c(0,0))
    if (!is.null(fit_full) && fit_full$convergence == 0) {
      pred_uni <- predict_unified_surface_with_se(fit_full, Y_train = Y, X_train = X_int, Z_train = S, Z_new = S_test, X_new = X_int_test)
      res["Unified_REML_MSE"] <- mean((pred_uni$Predicted_Y - EY_test)^2)
      
      if (test_type == "interaction" && length(target_idx) == 2) {
        tmp <- run_unified_pipeline_reml(
          Y, X_int, Z_cov = S, init_params_full = fit_full$kernel_params, 
          test_type = "interaction", idx_args = target_idx, global_rho = fit_full$kernel_params[1]
        )
        if(!is.null(tmp)) res["Unified_REML_pval"] <- as.numeric(tmp$p_value)
      } else {
        tmp <- run_unified_pipeline_reml(
          Y, X_int, Z_cov = S, init_params_full = fit_full$kernel_params, 
          test_type = "var_imp", idx_args = target_idx, global_rho = fit_full$kernel_params[1]
        )
        if(!is.null(tmp)) res["Unified_REML_pval"] <- as.numeric(tmp$p_value)
      }
    }
  }, error = function(e) { warning(paste("Unified REML Error:", e$message)) })
  
  # --- Method 5: Full GP (KMR) ---
  tryCatch({
    fit.full.reg <- Kernel.reg.Gauss(yy = Y, zz = t(S), xx = t(X_int), zz.p = t(S_test), lam.min = .0001, rho.min = .0001, sig.min = .0001, maxiter = 500, rho.max = 500, lam.max = 500)
    pred_full <- as.numeric(X_int_test %*% fit.full.reg$bet + fit.full.reg$h.hat)
    res["Full_GP_MSE"] <- mean((pred_full - EY_test)^2)
    
    if (!(length(target_idx) > 1 || test_type == "interaction")) {
      S_reord <- S[, c(target_idx, setdiff(1:p, target_idx))]
      test.full <- Kernel.test.Gauss(yy = Y, zz = t(S_reord), xx = t(X_int), lam.min = .0001, rho.min = .0001, sig.min = .0001, maxiter = 500, rho.max = 500, lam.max = 500, jitter = FALSE)
      if (length(test.full) == 6) res["Full_GP_pval"] <- as.numeric(test.full$pvalue[1])
    }
  }, error = function(e) { warning(paste("Full GP Error:", e$message)) })
  
  # --- Method 6: Old Pairwise GP ---
  tryCatch({
    fit.pw.reg <- Kernel.reg.Gauss.pairwise(yy = Y, zz = t(S), xx = t(X_int), zz.p = t(S_test), lam.min = .0001, rho.min = .0001, sig.min = .0001, maxiter = 500, rho.max = 500, lam.max = 500)
    pred_pw <- as.numeric(X_int_test %*% fit.pw.reg$bet + fit.pw.reg$h.hat)
    res["Old_Pairwise_GP_MSE"] <- mean((pred_pw - EY_test)^2)
    
    if (!(length(target_idx) > 1 || test_type == "interaction")) {
      S_reord <- S[, c(target_idx, setdiff(1:p, target_idx))]
      test.pw <- Kernel.test.Gauss.pairwise(yy = Y, zz = t(S_reord), xx = t(X_int), lam.min = .0001, rho.min = .0001, sig.min = .0001, maxiter = 500, rho.max = 500, lam.max = 500, jitter = FALSE)
      if ("pvalue" %in% names(test.pw)) res["Old_Pairwise_GP_pval"] <- as.numeric(test.pw$pvalue[1])
    }
  }, error = function(e) { warning(paste("Old Pairwise GP Error:", e$message)) })
  
  # --- Method 4: GAM Smooth ---
  tryCatch({
    df <- data.frame(Y = Y, S)
    df_test <- data.frame(S_test)
    
    if (test_type == "interaction" && length(target_idx) == 2) {
      bg_vars <- setdiff(1:p, target_idx)
      bg_terms <- if (length(bg_vars) > 0) paste(paste0("s(S", bg_vars, ")"), collapse = " + ") else ""
      main_ti <- paste0("ti(S", target_idx[1], ") + ti(S", target_idx[2], ")")
      int_ti <- paste0("ti(S", target_idx[1], ", S", target_idx[2], ")")
      
      form <- as.formula(if(bg_terms != "") paste("Y ~", bg_terms, "+", main_ti, "+", int_ti) else paste("Y ~", main_ti, "+", int_ti))
      fit_gam <- gam(form, data = df, method = "REML")
      
      s_table <- summary(fit_gam)$s.table
      
      # Strict substring match for both variables and 'ti'
      row_idx <- which(grepl("ti", rownames(s_table), fixed = TRUE) & 
                         grepl(paste0("S", target_idx[1]), rownames(s_table), fixed = TRUE) & 
                         grepl(paste0("S", target_idx[2]), rownames(s_table), fixed = TRUE))
      
      if(length(row_idx) > 0) res["GAM_Smooth_pval"] <- as.numeric(s_table[row_idx[1], "p-value"])
      res["GAM_Smooth_MSE"] <- as.numeric(mean((predict(fit_gam, newdata = df_test) - EY_test)^2))
      
    } else if (length(target_idx) == 1) {
      form <- as.formula(paste("Y ~", paste(paste0("s(S", 1:p, ")"), collapse = " + ")))
      fit_gam <- gam(form, data = df, method = "REML")
      
      s_table <- summary(fit_gam)$s.table
      
      # Exact match bypassing regex processing
      main_name <- paste0("s(S", target_idx, ")")
      row_idx <- grep(main_name, rownames(s_table), fixed = TRUE)
      
      if(length(row_idx) > 0) res["GAM_Smooth_pval"] <- as.numeric(s_table[row_idx[1], "p-value"])
      res["GAM_Smooth_MSE"] <- as.numeric(mean((predict(fit_gam, newdata = df_test) - EY_test)^2))
      
    } else {
      null_vars <- setdiff(1:p, target_idx)
      form_alt <- as.formula(paste("Y ~", paste(paste0("s(S", 1:p, ")"), collapse = " + ")))
      form_null <- as.formula(if(length(null_vars)==0) "Y ~ 1" else paste("Y ~", paste(paste0("s(S", null_vars, ")"), collapse = " + ")))
      
      m_null <- gam(form_null, data = df, method = "ML")
      m_alt  <- gam(form_alt, data = df, method = "ML")
      
      a_table <- anova(m_null, m_alt, test = "F")
      res["GAM_Smooth_pval"] <- as.numeric(a_table[["Pr(>F)"]][2])
      res["GAM_Smooth_MSE"] <- as.numeric(mean((predict(m_alt, newdata = df_test) - EY_test)^2))
    }
  }, error = function(e) { warning(paste("GAM Error:", e$message)) })
  
  # --- Method 5: Ordinary Least Squares ---
  tryCatch({
    df <- data.frame(Y = Y, S)
    df_test <- data.frame(S_test)
    
    if (test_type == "interaction" && length(target_idx) == 2) {
      form_null <- as.formula("Y ~ .")
      form_alt <- as.formula(paste0("Y ~ . + S", target_idx[1], ":S", target_idx[2]))
      
      fit_lm_null <- lm(form_null, data = df)
      fit_lm_alt <- lm(form_alt, data = df)
      
      a_table <- anova(fit_lm_null, fit_lm_alt)
      res["OLS_pval"] <- as.numeric(a_table[["Pr(>F)"]][2])
      res["OLS_MSE"] <- as.numeric(mean((predict(fit_lm_alt, newdata = df_test) - EY_test)^2))
      
    } else if (length(target_idx) == 1) {
      fit_lm <- lm(Y ~ ., data = df)
      target_name <- paste0("S", target_idx)
      coef_table <- summary(fit_lm)$coefficients
      
      if (target_name %in% rownames(coef_table)) {
        res["OLS_pval"] <- as.numeric(coef_table[target_name, "Pr(>|t|)"])
      }
      res["OLS_MSE"] <- as.numeric(mean((predict(fit_lm, newdata = df_test) - EY_test)^2))
      
    } else {
      null_vars <- setdiff(1:p, target_idx)
      form_alt <- as.formula("Y ~ .")
      form_null <- as.formula(if(length(null_vars)==0) "Y ~ 1" else paste("Y ~", paste(colnames(S)[null_vars], collapse = " + ")))
      
      fit_lm_null <- lm(form_null, data = df)
      fit_lm_alt <- lm(form_alt, data = df)
      
      a_table <- anova(fit_lm_null, fit_lm_alt)
      res["OLS_pval"] <- as.numeric(a_table[["Pr(>F)"]][2])
      res["OLS_MSE"] <- as.numeric(mean((predict(fit_lm_alt, newdata = df_test) - EY_test)^2))
    }
  }, error = function(e) { warning(paste("OLS Error:", e$message)) })
  
  # --- Method 6: BKMR ---
  if (run_bkmr) {
    tryCatch({
      if (test_type == "interaction" || length(target_idx) > 1) {
        bkmr_groups <- rep(0, p)
        bkmr_groups[target_idx] <- 1
        bkmr_groups[-target_idx] <- 2:(p - length(target_idx) + 1)
        fit_bkmr <- kmbayes(y = Y, Z = S, X = X_int, iter = 10000, varsel = TRUE, groups = bkmr_groups, verbose = FALSE)
        if (test_type != "interaction") {
          pips <- ExtractPIPs(fit_bkmr)
          res["BKMR_PIP"] <- as.numeric(pips$groupPIP[pips$group == 1][1])
        }
      } else {
        fit_bkmr <- kmbayes(y = Y, Z = S, X = X_int, iter = 10000, varsel = TRUE, verbose = FALSE)
        pips <- ExtractPIPs(fit_bkmr)
        target_name <- paste0("S", target_idx)
        res["BKMR_PIP"] <- as.numeric(pips$PIP[pips$variable == target_name])
      }
      
      sel_burnin <- seq(floor(fit_bkmr$iter / 2) + 1, fit_bkmr$iter)
      # Predict at the held-out test points (S_test), matching the other methods' MSE
      h_new <- bkmr::ComputePostmeanHnew(fit_bkmr, Znew = S_test, sel = sel_burnin, method = "approx")
      beta_mean <- colMeans(fit_bkmr$beta[sel_burnin, , drop = FALSE])
      preds <- as.numeric(X_int_test %*% beta_mean + h_new$postmean)
      res["BKMR_MSE"] <- as.numeric(mean((preds - EY_test)^2))
      
    }, error = function(e) { warning(paste("BKMR Error:", e$message)) })
  }
  
  return(res)
}

# ---------------------------------------------------------
# MASTER SIMULATION LOOP
# ---------------------------------------------------------
run_formal_simulation <- function(B = 1000, n = 300, p = 10, rho_corr = 0.0, snr = 1.0, 
                                  mean_fn, is_null_sim = FALSE, run_bkmr = FALSE, 
                                  target_idx = 1, test_type = "main") {
  
  cat(paste0("\nStarting simulation: B=", B, ", n=", n, ", p=", p, ", rho=", rho_corr, ", SNR=", snr, "\n"))
  cat(paste0("Condition: ", ifelse(is_null_sim, "NULL", "ALTERNATIVE"), "\n"))
  cat(paste0("Testing: ", toupper(test_type), "\n"))
  
  explicit_exports <- c(
    "evaluate_all_methods", 
    "fit_unified_full_model_reml", "fit_unified_H0_fixed_rho_reml",
    "run_unified_score_test_reml", "run_unified_pipeline_reml",
    "eval_unified_full_model", "eval_unified_var_imp_fixed_rho", 
    "eval_unified_interaction_fixed_rho", "eval_unified_all_interactions_fixed_rho",
    "eval_unified_K_cross_and_diag", "predict_unified_surface_with_se",
    "Kernel.test.Gauss", "Kernel.reg.Gauss", "kernel.Gauss",
    "obj.loglik.se", "grad.loglik.se", "Trace",
    "Kernel.test.Gauss.pairwise", "Kernel.reg.Gauss.pairwise", 
    "create_d_mat_list_1d", "V_fn", "obj.loglik.fast", 
    "grad.loglik.fast", "kernel.Gauss.pairwise", "d.rho.kernel.Gauss.pairwise"
  )
  
  funcs_to_export <- unique(c(explicit_exports, as.character(lsf.str())))
  
  results_matrix <- foreach(b = 1:B, .combine = rbind, 
                            .packages = c("MASS", "PairwiseSims", "mgcv", "bkmr", "fields"), 
                            .export = funcs_to_export) %dopar% {
                              evaluate_all_methods(n = n, p = p, rho_corr = rho_corr, mean_fn = mean_fn, snr = snr, run_bkmr = run_bkmr, target_idx = target_idx, test_type = test_type)
                            }
  rownames(results_matrix) <- NULL
  return(as.data.frame(results_matrix))
}

# ---------------------------------------------------------
# AGNOSTIC STUDY REPORT GENERATOR 
# ---------------------------------------------------------
generate_study_report <- function(results_list, alpha = 0.10, save_path = NULL) {
  summary_tables <- list()
  cat(paste0("\n=== GLOBAL SIMULATION SUMMARY (alpha = ", alpha, ") ===\n"))
  
  for (scenario_name in names(results_list)) {
    res_df <- results_list[[scenario_name]]
    
    pval_cols <- grep("pval", colnames(res_df), value = TRUE)
    report_df <- data.frame()
    if (length(pval_cols) > 0) {
      rejection_rates <- colMeans(res_df[, pval_cols, drop = FALSE] < alpha, na.rm = TRUE)
      na_counts_pval  <- colSums(is.na(res_df[, pval_cols, drop = FALSE]))
      report_df <- data.frame(Method = pval_cols, Rejection_Rate = rejection_rates, NAs = na_counts_pval, row.names = NULL)
      cat(paste0("\n--- ", toupper(scenario_name), ": HYPOTHESIS TESTING ---\n"))
      print(report_df, digits = 3)
    }
    
    mse_cols <- grep("MSE", colnames(res_df), value = TRUE)
    mse_df <- data.frame()
    if (length(mse_cols) > 0) {
      mse_means <- colMeans(res_df[, mse_cols, drop = FALSE], na.rm = TRUE)
      na_counts_mse <- colSums(is.na(res_df[, mse_cols, drop = FALSE]))
      mse_df <- data.frame(Method = mse_cols, Mean_MSE = mse_means, NAs = na_counts_mse, row.names = NULL)
      cat(paste0("\n--- ", toupper(scenario_name), ": PREDICTIVE ACCURACY ---\n"))
      print(mse_df, digits = 4)
    }
    
    pip_cols <- grep("PIP", colnames(res_df), value = TRUE)
    bkmr_report <- data.frame()
    if (length(pip_cols) > 0) {
      bkmr_report <- data.frame(Method=character(0), Rate_50=numeric(0), Rate_alpha=numeric(0), NAs=numeric(0))
      for (pip_col in pip_cols) {
        if (!all(is.na(res_df[[pip_col]]))) {
          rate_50 <- mean(res_df[[pip_col]] > 0.50, na.rm = TRUE)
          rate_a  <- mean(res_df[[pip_col]] > (1 - alpha), na.rm = TRUE)
          pip_nas <- sum(is.na(res_df[[pip_col]]))
        } else {
          rate_50 <- NA
          rate_a  <- NA
          pip_nas <- nrow(res_df) 
        }
        bkmr_report <- rbind(bkmr_report, data.frame(Method=pip_col, Rate_50=rate_50, Rate_alpha=rate_a, NAs=pip_nas))
      }
      cat(paste0("\n--- ", toupper(scenario_name), ": BKMR INCLUSION PROBABILITIES ---\n"))
      print(bkmr_report, digits = 3)
    }
    summary_tables[[scenario_name]] <- list(Testing = report_df, Accuracy = mse_df, PIP = bkmr_report)
  }
  
  if (!is.null(save_path)) save(results_list, summary_tables, file = save_path)
  return(summary_tables)
}