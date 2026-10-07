# NCWell analysis using competitors (with spatial residuals) 

rm(list=ls())
library(tidyverse)

hpc <- 0
library(foreach)
library(doSNOW) # Using doSNOW instead of doParallel to support the progress bar
library(MASS)
library(fields)
library(Rcpp)
library(RcppArmadillo)
library(PairwiseSims) # Holds compiled C++ code
library(GpGp)

if(hpc){
  source("code/ss_anova_reml_unified.R")
  source("code/ss_anova_reml_data_analysis_pipeline_unified.R")
  
  # Swapped to Spatial versions
  source("code/KMR_sp.R") 
  source("code/nc_tract_analysis.R")
  source("code/competitor_analyses_sp.R")
  
  data <- read.csv("data//NC_tract_all_outcomes.csv")
  
  n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
  # fallback 
  if (is.na(n_cores)) {
    n_cores <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE", unset = "1"))
  }
}else{
  n_cores <- parallel::detectCores() - 1 
  setwd("~/GitHub/Pairwise_Additive_GP_Tests")
  
  source("code/KMR/ss_anova_reml_unified.R")
  source("code/KMR/ss_anova_reml_data_analysis_pipeline_unified.R")
  
  # Swapped to Spatial versions
  source("code/KMR/KMR_sp.R")
  source("code/analysis/nc_tract/nc_tract_analysis.R")
  source("code/analysis/competitor_analyses_sp.R")
  
  data <- read.csv("data//NC_tract_all_outcomes.csv")
}

library(Rcpp)
library(RcppArmadillo)
library(PairwiseSims)
library(GpGp)


sub_data <- data

outcome_name <- "obesity_pct"
Y <- sub_data[,outcome_name]
Y <- scale(Y)

# --- Extract Coordinates Separately ---
coords_sp <- sub_data %>%
  dplyr::select(lat, lon) %>%
  as.matrix()

# --- Extract S without lat/lon ---
S <- sub_data %>%
  dplyr::select(Age, Median.Income, Female.Pct, Black.Pct, 
                Hisp.Pct, Pop, pct_well, ends_with(".Mean_avg")) %>%
  as.matrix()

S.as <- sub_data %>%
  dplyr::select(Age, Median.Income, Female.Pct, Black.Pct, 
                Hisp.Pct, Pop, pct_well, "Arsenic.Mean_avg") %>%
  as.matrix()

# Log transformations
cols_to_log <- grep("\\.Mean_avg$|NatAm\\.Pct", colnames(S))
for (idx in cols_to_log) {
  if (min(S[, idx], na.rm = TRUE) > 0) {
    S[, idx] <- log(S[, idx])
  }
}

cols_to_log_as <- grep("\\.Mean_avg$|NatAm\\.Pct", colnames(S.as))
for (idx in cols_to_log_as) {
  if (min(S.as[, idx], na.rm = TRUE) > 0) {
    S.as[, idx] <- log(S.as[, idx])
  }
}

X <- matrix(1, nrow=length(Y))

# Min-Max Scaling for S matrices
for(j in 1:ncol(S)){
  S[,j] <- S[,j] - min(S[,j])
  S[,j] <- S[,j]/max(S[,j])
}

for(j in 1:ncol(S.as)){
  S.as[,j] <- S.as[,j] - min(S.as[,j])
  S.as[,j] <- S.as[,j]/max(S.as[,j])
}

# Min-Max Scaling for Coordinates (maintains lengthscale stability)
for(j in 1:ncol(coords_sp)){
  coords_sp[,j] <- coords_sp[,j] - min(coords_sp[,j])
  coords_sp[,j] <- coords_sp[,j]/max(coords_sp[,j])
}


# =========================================================
# Spatial GLS Joint Test for Metals (Replaces OLS)
# =========================================================
# S[, 1:7] contains covariates (Age, Median.Income, Female.Pct, Black.Pct, Hisp.Pct, Pop, pct_well)
# S[, 8:ncol(S)] contains the metals
X_gls_null <- cbind(X, S[, 1:7])
X_gls_alt <- cbind(X, S)

set.seed(123) # GpGp's ordering non-deterministic
fit_gls_null <- fit_model(y = as.numeric(Y), locs = coords_sp, X = X_gls_null, covfun_name = "matern_isotropic", silent = TRUE)
fit_gls_alt  <- fit_model(y = as.numeric(Y), locs = coords_sp, X = X_gls_alt, covfun_name = "matern_isotropic", silent = TRUE)

betas <- fit_gls_alt$betahat
se_betas <- fit_gls_alt$sebeta
p_vals_gls <- 2 * pnorm(-abs(betas / se_betas))

# Likelihood Ratio Test
lrt_stat <- 2 * (fit_gls_alt$loglik - fit_gls_null$loglik)
df_diff <- ncol(X_gls_alt) - ncol(X_gls_null)
p_val_joint <- pchisq(lrt_stat, df = df_diff, lower.tail = FALSE)

cat("\n==================================================\n")
cat(sprintf("Spatial GLS Joint Test for Metals (LRT p-value): %e\n", p_val_joint))
cat("==================================================\n\n")

print(round(p_vals_gls,3)[9:19])
# =========================================================
# Run Competitor Analyses
# =========================================================
iters <- 1

compet_analysis_full <- run_competitor_analyses(Y = Y, 
                                                X = X, 
                                                S = S, 
                                                coords_sp = coords_sp,
                                                bkmr_iter = iters, 
                                                bkmr_knots = 1000,
                                                cl = NULL)

save.image("outputs//nc_obesity_analysis_competitors_sp.Rdata")
