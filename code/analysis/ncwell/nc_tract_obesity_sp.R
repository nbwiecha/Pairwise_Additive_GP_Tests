# NCWell analysis using PAGP + spatial residual


rm(list=ls())
library(tidyverse)

hpc <- 0
library(foreach)
library(doSNOW) # Using doSNOW instead of doParallel to support the progress bar
library(MASS)
library(fields)
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(PairwiseSims) # Holds compiled C++ code

if(hpc){
  source("code/ss_anova_reml_unified.R")
  source("code/ss_anova_reml_data_analysis_pipeline_unified.R")
  source("code/KMR.R")
  source("code/nc_tract_analysis.R")
  data <- read.csv("data//NC_tract_all_outcomes.csv")
  
  n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
  # fallback 
  if (is.na(n_cores)) {
    n_cores <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE", unset = "1"))
  }
  # n_cores <- 20
}else{
  # Cap background BLAS threads (e.g., to 4 or half your physical cores)
  Sys.setenv(OPENBLAS_NUM_THREADS = 1)
  Sys.setenv(OMP_NUM_THREADS = 1)
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
  }
  
  n_cores <- parallel::detectCores() - 1 
  setwd("~/GitHub/Pairwise_Additive_GP_Tests")
  # setwd("H:/Pairwise_Additive_GP_Tests")
  # sourceCpp("code/KMR/ss_anova_kernel_unified_sp.cpp")
  source("code/KMR/ss_anova_reml_unified_sp.R")
  source("code/KMR/ss_anova_reml_data_analysis_pipeline_unified_sp.R")
  # source("code/KMR/KMR.R")
  source("code/analysis/nc_tract/nc_tract_analysis.R")
  data <- read.csv("data//NC_tract_all_outcomes.csv")
  
}

# Setup the parallel cluster
# cl <- makeCluster(n_cores)
# registerDoSNOW(cl)

Y <- scale(data$obesity_pct)

S <- data %>%
  dplyr::select(lat, lon, Age, Median.Income, Female.Pct, Black.Pct, 
                Hisp.Pct, Pop, pct_well, ends_with(".Mean_avg")) %>%
  as.matrix()

S.as <- data %>%
  dplyr::select(lat, lon, Age, Median.Income, Female.Pct, Black.Pct, 
                Hisp.Pct, Pop, pct_well, "Arsenic.Mean_avg") %>%
  as.matrix()


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

for(j in 1:ncol(S)){
  
  S[,j] <- S[,j] - min(S[,j])
  S[,j] <- S[,j]/max(S[,j])
  
}

for(j in 1:ncol(S.as)){
  
  S.as[,j] <- S.as[,j] - min(S.as[,j])
  S.as[,j] <- S.as[,j]/max(S.as[,j])
  
}
head(S)

n <- nrow(S)
p <- ncol(S)
p.as <- ncol(S.as)



init_params <- c(0, 0, 0, 0)

obesity_spat <- data_analysis_pipeline_parallel_sp(Y, X, S, S[,1:2])

round(obesity_spat$var.importance.pvals, 3)
round(obesity_spat$interaction.pvals, 3)

save(obesity_spat,file="outputs//nc_obesity_analysis_sp.Rdata")
