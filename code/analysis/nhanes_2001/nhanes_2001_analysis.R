# NHANES 2001 analysis using PAGP

rm(list=ls())
library(tidyverse)
library(foreach)
library(doSNOW) 
library(MASS)
library(fields)
library(Rcpp)
library(RcppArmadillo)
library(PairwiseSims) # Holds compiled C++ code

# =========================================================
# SYSTEM SETUP & PATHS
# =========================================================
hpc <- 1

if(hpc){
  # FIXED: Point to the new separated script
  source("code/ss_anova_reml_unified.R")
  source("code/ss_anova_reml_data_analysis_pipeline_unified.R")
  
  n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
  if (is.na(n_cores)) {
    n_cores <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE", unset = "1"))
  }
  
  load('data/nhanes_2001.Rdata')
  
} else {
  n_cores <- parallel::detectCores() - 1 
  setwd("~/GitHub/Pairwise_Additive_GP_Tests")
  
  # FIXED: Point to the new separated script
  source("code/KMR/ss_anova_reml_unified.R")
  source("code/KMR/ss_anova_reml_data_analysis_pipeline_unified.R")
  load('data/nhanes_2001/nhanes_2001.Rdata')
  
}

# Clean up environment if necessary
suppressWarnings(rm(xx, zz, yy))

# Setup the parallel cluster
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

# =========================================================
# DATA PREPARATION
# =========================================================
debug <- FALSE
if(debug){
  dat <- dat[1:300,]
}
# 1. Format Response
Y <- as.numeric(scale(dat$ln_ltl))

# 2. Format Nonparametric Predictors (S)
S <- dat %>%
  dplyr::select(-SEQN, -ln_ltl, -age, -age2, -female, -high, -some_coll, -college,
                -nh_black, -mexam, -othrace, -white_ct, -lymph_pct, -mon_pct,
                -neut_pct, -eos_pct, -bas_pct, -under, -over, -ln_cot) %>%
  as.matrix()

# 3. Format Fixed Covariates (X)
X_raw <- dat %>%
  dplyr::select(age, age2, female, high, some_coll, college,
                nh_black, mexam, othrace, white_ct, lymph_pct, mon_pct,
                neut_pct, eos_pct, bas_pct, under, over, ln_cot) %>%
  as.matrix()

# FIXED: Explicitly add the intercept column for the matrix algebra
X <- cbind(Intercept = 1, X_raw)

# 4. Transformations (Min-Max Scaling on S)
for(j in 1:ncol(S)){
  S[,j] <- S[,j] - min(S[,j])
  if(max(S[,j]) > 0) {
    S[,j] <- S[,j] / max(S[,j])
  }
}

n <- nrow(S)
p <- ncol(S)

# =========================================================
# EXECUTION
# =========================================================

cat("\n=========================================\n")
cat(sprintf("RUNNING NHANES PARALLEL PIPELINE (n=%d, p=%d)\n", n, p))
cat("=========================================\n")

nhanes_analysis <- data_analysis_pipeline_parallel(Y, X, Z_cov = S)

# Group tests
global_log_rho <- nhanes_analysis$global_log_rho
colnames(S)
group1 <- 1:8
group2 <- 9:10
group3 <- 11:18

group1_test <- run_unified_pipeline_reml(
  Y = Y, 
  X = X, 
  Z_cov = S,
  init_params_full = c(1,1),
  test_type = "var_imp",
  idx_args = group1, 
  global_rho = global_log_rho
)

group2_test <- run_unified_pipeline_reml(
  Y = Y, 
  X = X, 
  Z_cov = S,
  init_params_full = c(1,1),
  test_type = "var_imp",
  idx_args = group2, 
  global_rho = global_log_rho
)

group3_test <- run_unified_pipeline_reml(
  Y = Y, 
  X = X, 
  Z_cov = S,
  init_params_full = c(1,1),
  test_type = "var_imp",
  idx_args = group3, 
  global_rho = global_log_rho
)

# Print high-level results to console
print(nhanes_analysis)
print(group1_test)
print(group2_test)
print(group3_test)
# Save output and close cluster
stopCluster(cl)
save(nhanes_analysis, group1_test, group2_test, group3_test, file = "outputs/nhanes_analysis_results.Rdata")
cat("\nPipeline complete. Results saved to outputs/nhanes_analysis_results.Rdata\n")