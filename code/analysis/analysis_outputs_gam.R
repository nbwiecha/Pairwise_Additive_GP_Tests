# Analysis results: NC Well only with spatial PAGP model

rm(list=ls())
library(gratia)
library(tidyverse)
library(mgcv)
library(corrplot)
library(PairwiseSims)
setwd("~/GitHub/Pairwise_Additive_GP_Tests")
# setwd("H:/Pairwise_Additive_GP_Tests")


################################################################################
#                           NC Well                                            #
################################################################################
rm(list=ls())

load("outputs//nc_obesity_analysis_sp.Rdata")
source("code/KMR/ss_anova_reml_unified_sp.R")
source("code/KMR/ss_anova_reml_data_analysis_pipeline_unified_sp.R")

round(obesity_spat$var.importance.pvals, 4)

round(obesity_spat$var.importance.pvals, 4)[10:20]
# 
# round(compet_analysis_full$OLS[11:21,4],3)
# round(compet_analysis_full$GAM[10:20,4],3)
# round(compet_analysis_full$Previous_KMR[10:20],3)
# 
ncwell_table <- data.frame(
  Pairwise_GP = round(obesity_spat$var.importance.pvals, 4)[10:20]
,
  OLS = round(compet_analysis_full$OLS[10:20,4],4),
  KMR = round(compet_analysis_full$Previous_KMR[10:20],4),
GAM = round(compet_analysis_full$GAM[10:20,4],4)


)
write.csv(ncwell_table, file="outputs/ncwell_results.csv")

int.df <- round(obesity_spat$interaction.pvals, 4) %>%
  filter(`p-value` < 0.05)
write.csv(int.df, file="outputs/ncwell_interactions.csv")
round(analysis_full$ncwell_analysis$interaction.pvals, 4)
round(analysis_full$ncwell_analysis$interaction.pvals, 4) %>%
  filter(`p-value` < .05)

data <- read.csv("data//NC_tract_all_outcomes.csv")
sub_data <- data
outcome_name <- "obesity_pct"
Y <- sub_data[,outcome_name]
Y <- scale(Y)

# colnames(sub_data)
S <- sub_data %>%
  dplyr::select(lat, lon, Age, Median.Income, Female.Pct, Black.Pct, 
                Hisp.Pct, Pop, pct_well, ends_with(".Mean_avg")) %>%
  as.matrix()

cols_to_log <- grep("\\.Mean_avg$|NatAm\\.Pct", colnames(S))
for (idx in cols_to_log) {
  if (min(S[, idx], na.rm = TRUE) > 0) {
    S[, idx] <- log(S[, idx])
  }
}

colnames(S)
p <- ncol(S)
D_1D_list <- vector(mode = "list", length = p)
for (j in 1:p) {
  D_1D_list[[j]] <- as.matrix(dist(S[, j]))
}
names(D_1D_list) <- colnames(S)
coords_list <- list(D_1D_list = D_1D_list)

X <- matrix(1, nrow=length(Y))

for(j in 1:ncol(S)){
  
  S[,j] <- S[,j] - min(S[,j])
  S[,j] <- S[,j]/max(S[,j])
  
}

init_params <- c(0,0,0,0)

ncwell_int_test_sp <- run_unified_pipeline_reml_sp(
  Y = Y,
  X = X,
  Z_cov = S,
  coords_sp = S[,1:2],
  init_params_full = init_params,
  test_type = "all_interactions",
  idx_args = 11:20 # R-based index (1st covariate)
)
ncwell_int_test_sp


ncwell_metals <- run_unified_pipeline_reml_sp(
    Y = Y,
    X = X,
    Z_cov = S,
    coords_sp = S[,1:2],
    init_params_full = init_params,
    test_type = "var_imp",
    idx_args = 11:20 # R-based index (1st covariate)
  )
print(ncwell_metals)

save(ncwell_int_test, ncwell_metals, file="outputs/ncwell_obesity_int_joint.Rdata")

corrplot(cor(S))
round(cor(S),2)

data.gam <- data %>%
  mutate(Log_As = log(Arsenic.Mean_avg),
         Log_Pb = log(Lead.Mean_avg),
         Log_Cr = log(Chromium.Mean_avg),
         Log_Zn = log(Zinc.Mean_avg),
         Log_Ba = log(Barium.Mean_avg))

ncwell_form <-  "obesity_pct ~ ti(lat, k=100) +ti(lon,k=100) +
ti(lat, lon, k=20) + ti(Age) + ti(Median.Income) +
ti(Black.Pct) + ti(Hisp.Pct)  + ti(pct_well) + ti(Log_As) +
ti(Log_Pb) + ti(Log_As, Log_Pb) + 
ti(Age, Median.Income) +
ti(Age, Hisp.Pct) + ti(Age, pct_well) + 
ti(Median.Income, Black.Pct) + ti(Median.Income, Hisp.Pct) + ti(Median.Income, Log_As) +
ti(Black.Pct, Hisp.Pct)"

ncwell.gam <- gam(formula(ncwell_form), method="REML", data=data.gam)
plot(ncwell.gam, select=25)
summary(ncwell.gam)
gam.check(ncwell.gam)
draw(ncwell.gam, select=c(10,11,13, 14, 22, 25))
draw(ncwell.gam, select=c(25))

png("outputs/ncwell_main.png", width=900, height=600)
draw(ncwell.gam, select=c(9), contour=FALSE)
dev.off()


png("outputs/ncwell_int.png", width=900, height=600)
draw(ncwell.gam, select=c(11, 17))
dev.off()

write.csv(summary(ncwell.gam)$s.table, file="outputs/ncwell_gam_table.csv")

save(ncwell.gam, file="outputs/ncwell_obesity_analysis_gam.Rdata")

#### Competitors
load("outputs/nc_obesity_analysis_competitors_sp.Rdata")
round(compet_analysis_full$GLS,3)
round(compet_analysis_full$GAM,3)

load("outputs/ncwell_obesity_analysis_gam.Rdata")
