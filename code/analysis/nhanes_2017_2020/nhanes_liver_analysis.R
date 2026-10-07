# Analysis of NHANES liver data using PAGP

library(tidyverse)
library(PairwiseSims) # Holds compiled C++ code
library(foreach)
library(doSNOW) # Using doSNOW instead of doParallel to support the progress bar

hpc <- 1

if(hpc){
  source("code/ss_anova_reml_unified.R")
  source("code/ss_anova_reml_data_analysis_pipeline_unified.R")
  source("code/KMR.R")
  load("data/NHANES_2017_2020_Liver.Rdata")
  
  n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
  # fallback 
  if (is.na(n_cores)) {
    n_cores <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE", unset = "1"))
  }
  # n_cores <- 20
}else{
  n_cores <- parallel::detectCores() - 1 
  setwd("~/GitHub/Pairwise_Additive_GP_Tests")
  # setwd("H:/Pairwise_Additive_GP_Tests")
  source("code/KMR/ss_anova_reml_unified.R")
  source("code/KMR/ss_anova_reml_data_analysis_pipeline_unified.R")
  load("data/NHANES_2017_2020_Liver.Rdata")
}

cl <- makeCluster(n_cores)
registerDoSNOW(cl)

debug <- FALSE
if(debug){
  data <- data[1:100,]
}

Y_ALT <- as.numeric(scale(log(data$ALT)))
Y_AST <- as.numeric(scale(log(data$AST)))

S <- data %>%
  mutate(Log_Mercury = log(Mercury),
         Log_Manganese = log(Manganese),
         Log_Lead = log(Lead),
         Log_Cadmium = log(Cadmium),
         Log_Selenium = log(Selenium)) %>%
  dplyr::select(Log_Mercury, Log_Manganese, Log_Lead, Log_Cadmium, Log_Selenium, 
                Male, Age, Poverty_Ratio, 
                Alcohol_User, Smoker, HDL_C, Total_Cholesterol,
                Energy_Intake_kcal, Hypertension, Diabetes, Race_Ethnicity) %>%
  mutate(Black = Race_Ethnicity=="Non-Hispanic Black") %>%
  dplyr::select(-Race_Ethnicity)

for(j in 1:ncol(S)){
  
  S[,j] <- S[,j] - min(S[,j])
  S[,j] <- S[,j]/max(S[,j])
  
}

X <- model.matrix(ALT~Education_Level + Sedentary_Cat, data=data)
head(X)

p <- ncol(S)

for(j in 1:p){
  if(is.numeric(S[,j])){
    hist(S[,j], main=colnames(S)[j], breaks=100)    
  }
  
}

# Convert S to a numeric matrix for the updated pipeline functions
S_mat <- as.matrix(S)

ALT_analysis_adj <- tryCatch(
  {
    # Attempt to run the pipeline with the full X matrix
    data_analysis_pipeline_parallel(Y_ALT, X, S_mat)
  },
  error = function(e) {
    # If the full X fails, print the error and fall back to intercept-only
    message("Full X matrix failed. Falling back to intercept only. Error: ", e$message)
    
    # Using drop = FALSE ensures X remains a 1-column matrix rather than becoming a vector
    data_analysis_pipeline_parallel(Y_ALT, X[, 1, drop = FALSE], S_mat)
  }
)


ALT_analysis_adj
ALT_int_test <- run_unified_all_interactions_pipeline_reml(Y_ALT, X[,1,drop=FALSE], S_mat, 
                                                           c(1,1), ALT_analysis_adj$global_log_rho)
ALT_int_test

ALT_metals <- run_unified_var_imp_pipeline_reml(Y_ALT, X[,1,drop=FALSE], S_mat,
                                                init_params_full=c(1,1),
                                                test_idx_vec = 1:5, global_log_rho = ALT_analysis_adj$global_log_rho)
ALT_metals

AST_analysis_adj <- tryCatch(
  {
    # Attempt to run the pipeline with the full X matrix
    data_analysis_pipeline_parallel(Y_AST, X, S_mat)
  },
  error = function(e) {
    # If the full X fails, print the error and fall back to intercept-only
    message("Full X matrix failed. Falling back to intercept only. Error: ", e$message)
    
    # Using drop = FALSE ensures X remains a 1-column matrix rather than becoming a vector
    data_analysis_pipeline_parallel(Y_AST, X[, 1, drop = FALSE], S_mat)
  }
)

AST_analysis_adj
AST_int_test <- run_unified_all_interactions_pipeline_reml(Y_AST, X[,1,drop=FALSE], S_mat, 
                                                           c(1,1), AST_analysis_adj$global_log_rho)
AST_int_test

AST_metals <- run_unified_var_imp_pipeline_reml(Y_AST, X[,1,drop=FALSE], S_mat,
                                                init_params_full=c(1,1),
                                                test_idx_vec = 1:5, global_log_rho = AST_analysis_adj$global_log_rho)

AST_metals
save.image(file="outputs/NHANES_Liver_Analysis.Rdata")