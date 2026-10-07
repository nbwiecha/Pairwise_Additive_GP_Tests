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
# 1. DEFINE FUNCTIONAL FORMS
# S1 and S2 enter ONLY through the pure interaction
# =========================================================

# Null: S1 has no effect (Interaction term vanishes)
mean_fn_null <- function(S) {
  sin(2 * pi * S[,3]) + sin(2 * pi * S[,4])
}

# Alternative: S1 is present via pure interaction with S2
mean_fn_alt <- function(S) {
  sin(2 * pi * S[,3]) + sin(2 * pi * S[,4]) + 15 * (S[,1] - 0.5) * (S[,2] - 0.5)
}

# =========================================================
# 2. SIMULATION PARAMETERS
# =========================================================
B_reps <- 2000
n_val <- 100
p_vars <- 10
snr_val <- 1.0
rho_val <- 0.0
run_bkmr_sim <- FALSE 

# =========================================================
# 3. RUN SIMULATIONS
# =========================================================

cat("\n=========================================================\n")
cat(sprintf("EVALUATING PURE INTERACTION TEST (n = %d, p = %d, rho = %.2f)\n", n_val, p_vars, rho_val))
cat("=========================================================\n")

# Run Null Simulation (Type I Error)
res_null <- run_formal_simulation(B = B_reps, 
                                  n = n_val, 
                                  p = p_vars, 
                                  rho_corr = rho_val, 
                                  snr = snr_val, 
                                  mean_fn = mean_fn_null, 
                                  is_null_sim = TRUE,
                                  run_bkmr = run_bkmr_sim,
                                  target_idx = c(1, 2),
                                  test_type = "interaction")

# Run Alternative Simulation (Power)
res_alt  <- run_formal_simulation(B = B_reps, 
                                  n = n_val, 
                                  p = p_vars, 
                                  rho_corr = rho_val, 
                                  snr = snr_val, 
                                  mean_fn = mean_fn_alt, 
                                  is_null_sim = FALSE,
                                  run_bkmr = run_bkmr_sim,
                                  target_idx = c(1, 2),
                                  test_type = "interaction")

# Clean up cluster
stopCluster(cl)

# =========================================================
# 4. SUMMARIZE AND SAVE OUTPUTS
# =========================================================
alpha_thresh <- 0.10

# 1. Pack the payload using highly descriptive names for the report
payload <- list(
  Type_I_Error = res_null,
  Power        = res_alt
)

# 2. Fire the agnostic summarizer
# It will print to the console and save everything to the .Rdata file in one step
study_results <- generate_study_report(
  results_list = payload, 
  alpha = alpha_thresh, 
  save_path = "outputs/SimResults_Scenario5_TestInteraction_n100_p10.Rdata"
)
