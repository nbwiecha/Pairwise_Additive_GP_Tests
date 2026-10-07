# NHANES 2001 analysis using Competitors

rm(list=ls())
library(tidyverse)
library(foreach)
library(doSNOW) 
library(MASS)
library(fields)
library(Rcpp)
library(RcppArmadillo)
library(PairwiseSims) # Holds compiled C++ code
library(mgcv)
library(bkmr)
library(car)
# =========================================================
# SYSTEM SETUP & PATHS
# =========================================================
hpc <- 1

if(hpc){
  # FIXED: Point to the new separated script
  source("code/ss_anova_reml_unified.R")
  source("code/ss_anova_reml_data_analysis_pipeline_unified.R")
  source("code/KMR.R")
  source("code/competitor_analyses.R")

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
  source("code/KMR/KMR.R")
  source("code/analysis/competitor_analyses.R")
  
  load('data/nhanes_2001/nhanes_2001.Rdata')
  
}

# Clean up environment if necessary
suppressWarnings(rm(xx, zz, yy))

# Setup the parallel cluster
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

debug <- FALSE
if(debug){
  dat <- dat[1:100,]
  iters <- 500
}else{
  iters <- 10^5
}
# =========================================================
# DATA PREPARATION
# =========================================================

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

nhanes_analysis_compet <- run_competitor_analyses(Y, X, S, 
                                                  bkmr_iter = iters, 
                                                  bkmr_groups = c(rep(1,8),
                                                                  rep(2,2),
                                                                  rep(3,8)))
nhanes_analysis_compet
stopCluster(cl)
save(nhanes_analysis_compet, file = "outputs/nhanes_analysis_results_competitors.Rdata")

round(cor(S), 2)
round(cor(S[,1:8]), 2)

lm.fit.X <- lm(Y ~ X - 1)

lm.fit <- lm(Y~ X + S -1)
summary(lm.fit)

lm.fit.sig <- lm(Y ~ X + S[,c("ln_PCB180", "ln_pcndf23478")]-1)
summary(lm.fit.sig)
cor(S[,c("ln_PCB180", "ln_pcndf23478")])
confidenceEllipse(lm.fit.sig, which.coef=c(20,21))

lm.fit.153.180 <- lm(Y ~ X + S[,c("ln_PCB153", "ln_PCB180")]-1)
summary(lm.fit.153.180)
confint(lm.fit.153.180)
cor(S[,c("ln_PCB153", "ln_PCB180")])
confidenceEllipse(lm.fit.153.180, which.coef=c(20,21))
anova(lm.fit.X, lm.fit.153.180)

lm.fit.153 <- lm(Y ~ X + S[,c("ln_PCB153")]-1)
summary(lm.fit.153)

lm.fit.180 <- lm(Y ~ X + S[,c("ln_PCB180")]-1)
summary(lm.fit.180)


lm.fit.group1.out <- lm(Y ~ X + S[,9:ncol(S)]-1)
summary(lm.fit.group1.out)
anova(lm.fit, lm.fit.group1.out)

lm.fit.group2.out <- lm(Y ~ X + S[,-(9:10)]-1)
summary(lm.fit.group2.out)
anova(lm.fit, lm.fit.group2.out)

lm.fit.group3.out <- lm(Y ~ X + S[,-(11:ncol(S))]-1)
summary(lm.fit.group3.out)
anova(lm.fit, lm.fit.group3.out)
