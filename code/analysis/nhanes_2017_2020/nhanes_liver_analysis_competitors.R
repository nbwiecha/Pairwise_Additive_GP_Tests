# Analysis of NHANES liver data using competitors

library(tidyverse)
library(PairwiseSims) # Holds compiled C++ code
library(foreach)
library(doSNOW) # Using doSNOW instead of doParallel to support the progress bar
library(mgcv)
library(bkmr)

hpc <- 1

if(hpc){
  source("code/ss_anova_reml_unified.R")
  source("code/ss_anova_reml_data_analysis_pipeline_unified.R")
  source("code/KMR.R")
  source("code/competitor_analyses.R")
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
  source("code/KMR/KMR.R")
  source("code/KMR/ss_anova_reml_unified.R")
  source("code/KMR/ss_anova_reml_data_analysis_pipeline_unified.R")
  source("code/analysis/competitor_analyses.R")
  load("data/NHANES_2017_2020_Liver.Rdata")
}

# Guard against running on the old data file (n = 2,507 before never-smokers and
# never-drinkers were kept; n = 5,980 after the recode).
cat("Sample size:", nrow(data), "\n")
if (nrow(data) == 2507) stop("data/NHANES_2017_2020_Liver.Rdata is the old (pre-recode) file; rerun the cleaning script.")

# One BLAS thread per worker so the workers do not oversubscribe the allocated cores
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
cl <- makeCluster(n_cores)
registerDoSNOW(cl)
cat("Workers:", n_cores, "\n")

debug <- FALSE
if(debug){
  data <- data[1:100,]
  iters <- 500
}else{
  iters <- 10^5
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
                Race_Ethnicity, Male, Age, Poverty_Ratio, 
                Alcohol_User, Smoker, HDL_C, Total_Cholesterol,
                Energy_Intake_kcal, Hypertension, Diabetes) %>%
  mutate(Black = Race_Ethnicity=="Non-Hispanic Black") %>%
  dplyr::select(-Race_Ethnicity)

for(j in 1:ncol(S)){
  
  S[,j] <- S[,j] - min(S[,j])
  S[,j] <- S[,j]/max(S[,j])
  
}

X <- model.matrix(ALT~Education_Level + Sedentary_Cat, data=data)
head(X)

p <- ncol(S)

iters <- 1 # throws error for BKMR

compet_analysis_alt <- run_competitor_analyses(Y_ALT, X[,1,drop=FALSE], S, 
                                               bkmr_iter = iters, bkmr_knots = 1000)
compet_analysis_alt

compet_analysis_ast <- run_competitor_analyses(Y_AST, X[,1,drop=FALSE], S, 
                                               bkmr_iter = iters, bkmr_knots = 1000)
compet_analysis_ast
# 
# for(j in 1:p){
#   if(is.numeric(S[,j])){
#     hist(S[,j], main=colnames(S)[j], breaks=100)    
#   }
# 
# }

# GAMs
gam_alt <- gam(log(ALT) ~ s(log(Mercury)) + s(log(Manganese)) + s(log(Lead)) + s(log(Cadmium)) +
                 s(log(Selenium)) + Male + s(Age) + s(Poverty_Ratio) + Alcohol_User + 
                 Smoker + s(HDL_C) + s(Total_Cholesterol) + s(Energy_Intake_kcal) +
                 Hypertension + Diabetes + I(Race_Ethnicity=="Non-Hispanic Black"),
               method="REML", data=data)
summary(gam_alt)
gam.check(gam_alt)

gam_ast <- gam(log(AST) ~ s(log(Mercury)) + s(Manganese) + s(Lead) + s(Cadmium) +
                 s(Selenium) + Male + s(Age) + s(Poverty_Ratio) + Alcohol_User + 
                 Smoker + s(HDL_C) + s(Total_Cholesterol) + s(Energy_Intake_kcal) +
                 Hypertension + Diabetes + I(Race_Ethnicity=="Non-Hispanic Black"),
               method="REML", data=data)
summary(gam_ast)
gam.check(gam_ast)

save(compet_analysis_alt, compet_analysis_ast,
     gam_alt, gam_ast, file="outputs/NHANES_Liver_Analysis_Competitors.Rdata")
stopCluster(cl)
