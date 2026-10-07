rm(list=ls())

# 1. Define Environment (1 for Cluster, 0 for Local)
hpc <- 1 

# 2. Source Master Framework (This automatically handles the cluster, libraries, and sub-scripts)
if(hpc) {
  source("code/simulation_framework.R")
} else {
  source("~/GitHub/Pairwise_Additive_GP_Tests/code/final_simulations/simulation_framework.R")
}

# =========================================================
# 1. DEFINE FUNCTIONAL FORMS (Scenario 2)
# Comparing Cosine vs. Sine functional forms
# =========================================================

# --- FUNCTION 1: COSINE ---
# Null: S1 is absent
mean_fn_cos_null <- function(S) {
  cos(2 * pi * S[,2])+ 15 * (S[,3] - 0.5) * (S[,4] - 0.5)
}

# Alternative: S1 is present
mean_fn_cos_alt <- function(S) {
  cos(2 * pi * S[,1]) + cos(2 * pi * S[,2])+ 15 * (S[,3] - 0.5) * (S[,4] - 0.5)
}

# --- FUNCTION 2: SINE ---
# Null: S1 is absent
mean_fn_sin_null <- function(S) {
  sin(2 * pi * S[,2])+ 15 * (S[,3] - 0.5) * (S[,4] - 0.5)
}

# Alternative: S1 is present 
mean_fn_sin_alt <- function(S) {
  sin(2 * pi * S[,1]) + sin(2 * pi * S[,2])+ 15 * (S[,3] - 0.5) * (S[,4] - 0.5)
}

# =========================================================
# 2. SIMULATION PARAMETERS
# =========================================================
B_reps <- 2000
n_val <- 300
p_vars <- 10
snr_val <- 1.0
rho_val <- 0.0

# I left this as TRUE per your upload, but you can toggle to FALSE if you want to skip MCMC
run_bkmr_sim <- TRUE 

# =========================================================
# 3. RUN SIMULATIONS
# =========================================================

cat("\n=========================================================\n")
cat("EVALUATING FUNCTION 1: COSINE\n")
cat("=========================================================\n")

res_cos_null <- run_formal_simulation(B = B_reps, n = n_val, p = p_vars, 
                                      rho_corr = rho_val, snr = snr_val, 
                                      mean_fn = mean_fn_cos_null, is_null_sim = TRUE,
                                      run_bkmr = run_bkmr_sim)

res_cos_alt  <- run_formal_simulation(B = B_reps, n = n_val, p = p_vars, 
                                      rho_corr = rho_val, snr = snr_val, 
                                      mean_fn = mean_fn_cos_alt, is_null_sim = FALSE,
                                      run_bkmr = run_bkmr_sim)


cat("\n=========================================================\n")
cat("EVALUATING FUNCTION 2: SINE\n")
cat("=========================================================\n")

res_sin_null <- run_formal_simulation(B = B_reps, n = n_val, p = p_vars, 
                                      rho_corr = rho_val, snr = snr_val, 
                                      mean_fn = mean_fn_sin_null, is_null_sim = TRUE,
                                      run_bkmr = run_bkmr_sim)

res_sin_alt  <- run_formal_simulation(B = B_reps, n = n_val, p = p_vars, 
                                      rho_corr = rho_val, snr = snr_val, 
                                      mean_fn = mean_fn_sin_alt, is_null_sim = FALSE,
                                      run_bkmr = run_bkmr_sim)

# Clean up cluster
stopCluster(cl)

# =========================================================
# 4. SUMMARIZE AND SAVE OUTPUTS
# =========================================================
alpha_thresh <- 0.10

# 1. Pack the payload with all 4 evaluations
payload <- list(
  Cosine_Type_I_Error = res_cos_null,
  Cosine_Power        = res_cos_alt,
  Sine_Type_I_Error   = res_sin_null,
  Sine_Power          = res_sin_alt
)

# 2. Fire the agnostic summarizer
study_results <- generate_study_report(
  results_list = payload, 
  alpha = alpha_thresh, 
  save_path = "outputs/SimResults_Scenario2_Trig_p10_n300.Rdata"
)
