# Survey-weighted OLS sensitivity analysis for the NHANES 2017-2020 liver data
#
# ALT and the blood metals come from the exam/lab (MEC) component, so the MEC
# exam weight WTMECPRP is used. Following NHANES guidance, the survey design is
# built on the full P_DEMO sample and then subset to the analytic sample, so
# variance estimates reflect the complete design.
rm(list=ls())
library(tidyverse)
library(nhanesA)
library(survey)
setwd("~/GitHub/Pairwise_Additive_GP_Tests")
load("data/NHANES_2017_2020_Liver.Rdata")
cat("Analytic sample size:", nrow(data), "\n")

# Design variables and weights for everyone in NHANES 2017-March 2020
demo <- nhanes("P_DEMO", translated = FALSE) %>%
  select(SEQN, SDMVPSU, SDMVSTRA, WTMECPRP)

full <- demo %>%
  left_join(data %>% mutate(in_analysis = TRUE), by = "SEQN") %>%
  mutate(in_analysis = !is.na(in_analysis),
         Black = Race_Ethnicity == "Non-Hispanic Black") %>%
  filter(!is.na(WTMECPRP) & WTMECPRP > 0)   # interview-only participants have no MEC weight
stopifnot(sum(full$in_analysis) == nrow(data))

nhanesDesign <- svydesign(id = ~SDMVPSU,       # Primary Sampling Units (PSU)
                          strata  = ~SDMVSTRA, # Stratification used in the survey
                          weights = ~WTMECPRP, # MEC exam weight
                          nest    = TRUE,      # PSUs are nested within strata
                          data    = full)
nhanesDesign <- subset(nhanesDesign, in_analysis)

# Same covariates as the PAGP model (plus a quadratic age term)
form <- log(ALT) ~ log(Mercury) + log(Manganese) + log(Cadmium) + log(Selenium) + log(Lead) +
  Male + Hypertension + Age + I(Age^2) + Poverty_Ratio + Alcohol_User + HDL_C +
  Total_Cholesterol + Energy_Intake_kcal + Diabetes + Black + Smoker

# Unweighted OLS for comparison
fit.lm <- lm(form, data = full %>% filter(in_analysis))
summary(fit.lm)

fit.lm.svy <- svyglm(form, design = nhanesDesign)
summary(fit.lm.svy)

# Joint Wald test of the five metals
regTermTest(fit.lm.svy, ~ log(Mercury) + log(Manganese) + log(Lead) + log(Cadmium) + log(Selenium),
            method = "Wald")

metals <- c("log(Mercury)", "log(Manganese)", "log(Lead)", "log(Cadmium)", "log(Selenium)")
svy_table <- data.frame(
  Metal = metals,
  OLS_est = coef(fit.lm)[metals],
  OLS_p = summary(fit.lm)$coefficients[metals, 4],
  SVY_est = coef(fit.lm.svy)[metals],
  SVY_p = summary(fit.lm.svy)$coefficients[metals, 4],
  row.names = NULL
)
print(svy_table, digits = 3)
write.csv(svy_table, file = "outputs/nhanes_liver_svy_ols.csv", row.names = FALSE)
