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
# 1. DEFINE FUNCTIONAL FORMS (Scenario 1)
# Linear in S1 and S2, with pure interaction of S3 and S4
# =========================================================

mean_fn <- function(S) {
  sin(5 * pi * S[,1]) + sin(5 * pi * S[,2]) + 15 * sin(2*S[,3]) * sin(3*S[,4])
}

# =========================================================
# 2. SIMULATION PARAMETERS
# =========================================================
B_reps <- 2000
n_val <- 200
snr_val <- 1.0
rho_val <- 0
run_bkmr_sim <- FALSE

# Define the sequence of p values to evaluate
p_sequence <- c(5, 7, 10, 13, 15, 17, 20)
# p_sequence <- c(5, 10, 13, 15, 17, 20, 25, 30) 

# Lists to store the raw results of each p-dimension
results_null_list <- list()
results_alt_list <- list()

# =========================================================
# 3. RUN SIMULATION LOOP OVER p
# =========================================================
for (p_vars in p_sequence) {
  cat(paste0("\n=========================================================\n"))
  cat(paste0("EVALUATING DIMENSION: p = ", p_vars, "\n"))
  cat(paste0("=========================================================\n"))
  
  # ---------------------------------------------------------
  # NULL STUDY (Type I Error)
  # ---------------------------------------------------------
  # res_null <- run_formal_simulation(B = B_reps, 
  #                                   n = n_val, 
  #                                   p = p_vars, 
  #                                   rho_corr = rho_val, 
  #                                   snr = snr_val, 
  #                                   mean_fn = mean_fn_null, 
  #                                   is_null_sim = TRUE,
  #                                   run_bkmr = run_bkmr_sim)
  # 
  # res_null$p_dim <- p_vars
  # results_null_list[[paste0("p_", p_vars)]] <- res_null
  # 
  # ---------------------------------------------------------
  # ALTERNATIVE STUDY (Power)
  # ---------------------------------------------------------
  res_alt <- run_formal_simulation(B = B_reps, 
                                   n = n_val, 
                                   p = p_vars, 
                                   rho_corr = rho_val, 
                                   snr = snr_val, 
                                   mean_fn = mean_fn, 
                                   is_null_sim = FALSE,
                                   run_bkmr = run_bkmr_sim)
  
  res_alt$p_dim <- p_vars
  results_alt_list[[paste0("p_", p_vars)]] <- res_alt
}

# Clean up
stopCluster(cl)

# =========================================================
# 4. AGGREGATE AND SAVE OUTPUTS
# =========================================================
alpha_thresh <- 0.10

# 1. Dynamically pack the payload using descriptive names for the report
payload <- list()
for (p_vars in p_sequence) {
  payload[[paste0("Type_I_Error_p", p_vars)]] <- results_null_list[[paste0("p_", p_vars)]]
  payload[[paste0("Power_p", p_vars)]]        <- results_alt_list[[paste0("p_", p_vars)]]
}

# 2. Fire the agnostic summarizer
# It will print 10 clean tables to the console and save everything in one step
study_results <- generate_study_report(
  results_list = payload, 
  alpha = alpha_thresh, 
  save_path = "outputs/SimResults_Scenario1_Dimensionality_n200.Rdata"
)
