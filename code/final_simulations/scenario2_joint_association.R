# ---------------------------------------------------------
# SIMULATION SCENARIO 2: Joint Association / Group Testing
# Testing the joint significance of a block of variables (S1 - S5)
# ---------------------------------------------------------

rm(list=ls())

# 1. Define Environment (1 for Cluster, 0 for Local)
hpc <- 1 

if(hpc) {
  source("code/simulation_framework.R")
} else {
  source("~/GitHub/Pairwise_Additive_GP_Tests/code/final_simulations/simulation_framework.R")
}

# =========================================================
# 1. DEFINE FUNCTIONAL FORMS (Scenario 2)
# =========================================================
# We are testing the JOINT hypothesis for the block: S1, S2, S3, S4, S5.
# H0: f(S1, S2, S3, S4, S5) = 0
# H1: f(S1, S2, S3, S4, S5) != 0

# Null Function: The entire target block (S1-S5) is absent. 
# Only background variables (S6, S7) drive the response.
mean_fn_null <- function(S) {
  sin(2 * pi * S[,6]) + sin(2 * pi * S[,7])
}

# Alternative Function: The block has an effect (e.g., S1 and S2 are active)
# Background variables (S6, S7) are still present.
mean_fn_alt <- function(S) {
  sin(2 * pi * S[,1]) + sin(2 * pi * S[,2]) + sin(2 * pi * S[,6]) + sin(2 * pi * S[,7])
}

# =========================================================
# 2. SIMULATION PARAMETERS
# =========================================================
B_reps <- 2000
n_val <- 100
p_total <- 10            
target_block <- 1:5      # The block of variables we are testing jointly
snr_val <- 1.0
run_bkmr_sim <- TRUE    # Toggle to TRUE if you want to run Group PIPs

# Evaluate across different correlation structures
rho_sequence <- c(0.0, 0.5)

results_null_list <- list()
results_alt_list <- list()

# =========================================================
# 3. RUN SIMULATION LOOP OVER RHO (Correlation)
# =========================================================
for (r in rho_sequence) {
  cat(sprintf("\n=========================================================\n"))
  cat(sprintf("EVALUATING JOINT CORRELATION: rho = %.2f\n", r))
  cat(sprintf("=========================================================\n"))
  
  # ---------------------------------------------------------
  # NULL STUDY (Type I Error for the Block)
  # ---------------------------------------------------------
  res_null <- run_formal_simulation(B = B_reps, 
                                    n = n_val, 
                                    p = p_total, 
                                    rho_corr = r, 
                                    snr = snr_val, 
                                    mean_fn = mean_fn_null, 
                                    is_null_sim = TRUE,
                                    run_bkmr = run_bkmr_sim,
                                    target_idx = target_block) 
  
  res_null$rho_val <- r
  results_null_list[[paste0("rho_", r)]] <- res_null
  
  # ---------------------------------------------------------
  # ALTERNATIVE STUDY (Power for the Block)
  # ---------------------------------------------------------
  res_alt <- run_formal_simulation(B = B_reps, 
                                   n = n_val, 
                                   p = p_total, 
                                   rho_corr = r, 
                                   snr = snr_val, 
                                   mean_fn = mean_fn_alt, 
                                   is_null_sim = FALSE,
                                   run_bkmr = run_bkmr_sim,
                                   target_idx = target_block) 
  
  res_alt$rho_val <- r
  results_alt_list[[paste0("rho_", r)]] <- res_alt
}

# Clean up
stopCluster(cl)

# =========================================================
# 4. AGGREGATE AND SAVE OUTPUTS
# =========================================================
alpha_thresh <- 0.10

payload <- list()
for (r in rho_sequence) {
  payload[[paste0("Joint_Type_I_Error_rho", r)]] <- results_null_list[[paste0("rho_", r)]]
  payload[[paste0("Joint_Power_rho", r)]]        <- results_alt_list[[paste0("rho_", r)]]
}

# Fire the agnostic summarizer
study_results <- generate_study_report(
  results_list = payload, 
  alpha = alpha_thresh, 
  save_path = "outputs/SimResults_Scenario2_JointAssociation.Rdata"
)
