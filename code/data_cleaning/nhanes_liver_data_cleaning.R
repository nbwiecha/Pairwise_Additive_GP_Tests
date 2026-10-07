# Install packages if you haven't already
# install.packages(c("nhanesA", "dplyr", "purrr"))

library(nhanesA)
library(tidyverse)
library(purrr)

setwd("~/GitHub/Pairwise_Additive_GP_Tests")
# 1. Download the required tables for 2017-March 2020 (Prefix "P_")
message("Downloading NHANES tables...")
tables_to_download <- c("P_DEMO", "P_BMX", "P_PAQ", "P_SMQ", "P_ALQ", 
                        "P_HDL", "P_TCHOL", "P_DR1TOT", "P_ALB_CR", 
                        "P_BIOPRO", "P_BPQ", "P_DIQ", "P_PBCD")

# Fetch all tables into a list
nhanes_list <- lapply(tables_to_download, nhanes)

# 2. Merge all tables using SEQN (Respondent Sequence Number)
merged_data <- nhanes_list %>% reduce(full_join, by = "SEQN")

# 3. Extract, clean, and derive variables
final_data <- merged_data %>%
  mutate(
    # Demographics
    Age = RIDAGEYR,
    Male = ifelse(RIAGENDR == "Male", TRUE, 0),
    
    # Education: 1-2 = <HS, 3 = HS, 4-5 = >HS
    # Education_Level = case_when(
    #   DMDEDUC2 %in% c(1, 2) ~ "Below high school",
    #   DMDEDUC2 == 3 ~ "High school",
    #   DMDEDUC2 %in% c(4, 5) ~ "Above high school",
    #   TRUE ~ NA_character_
    # ),
    Education_Level = DMDEDUC2,
    
    # Race/Ethnicity
    Race_Ethnicity = RIDRETH3,
    
    Poverty_Ratio = INDFMPIR,
    
    # Questionnaires
    Smoker = SMQ040 != "Not at all", # 1/2 = Current smoker everyday/some days
    Alcohol_User = case_when(
      ALQ121 %in% c("Every day", "Nearly every day", "3 to 4 times a week", 
                    "2 times a week", "Once a week", "2 to 3 times a month", 
                    "Once a month") ~ TRUE,
      
      ALQ121 %in% c("7 to 11 times in the last year", "3 to 6 times in the last year", 
                    "1 to 2 times in the last year", "Never in the last year") ~ FALSE,
      
      # Explicitly handling non-responses and any baseline NAs
      ALQ121 %in% c("Refused", "Don't know") ~ NA,
      TRUE ~ NA 
    ),
    Hypertension = ifelse(BPQ020 == "Yes", 1, 0),
    Diabetes = ifelse(DIQ010 == "Yes", 1, 0),
    
    # Physical/Diet
    BMI = BMXBMI,
    Energy_Intake_kcal = DR1TKCAL,
     
    # Sedentary Time (PAD680 is in minutes)
    Sedentary_Hrs = PAD680 / 60,
    Sedentary_Cat = case_when(
      Sedentary_Hrs < 3 ~ "<3 h",
      Sedentary_Hrs >= 3 & Sedentary_Hrs <= 6 ~ "3-6 h",
      Sedentary_Hrs > 6 ~ ">6 h",
      TRUE ~ NA_character_
    ),
    Sedentary_Cat = as.factor(Sedentary_Cat),
    
    # Routine Labs
    HDL_C = LBDHDD,
    Total_Cholesterol = LBXTC,
    Urinary_Albumin = URXUMA,
    Serum_Total_Protein = LBXSTP,
    Alkaline_Phosphatase = LBXSAPSI,
    AST = LBXSASSI,
    ALT = LBXSATSI,
    GGT = LBXSGTSI,
    Total_Bilirubin = LBXSTB,
    
    # Derived Liver Variables
    ALT_AST_Ratio = ALT / AST,
    HSI = (8 * ALT_AST_Ratio) + BMI + ifelse(Diabetes == 1, 2, 0) + ifelse(Male == 0, 2, 0),
    
    # Heavy Metals
    Mercury = LBXTHG,
    Manganese = LBXBMN,
    Lead = LBXBPB,
    Cadmium = LBXBCD,
    Selenium = LBXBSE
  ) %>%
  # Select only the needed variables to clean up the dataframe
  select(
    SEQN, Age, Male, Education_Level, Race_Ethnicity, Poverty_Ratio, Smoker, 
    Alcohol_User, BMI, Sedentary_Cat, HDL_C, Total_Cholesterol, Energy_Intake_kcal, 
    Urinary_Albumin, Serum_Total_Protein, Hypertension, Diabetes, 
    Alkaline_Phosphatase, AST, ALT, ALT_AST_Ratio, GGT, Total_Bilirubin, HSI, 
    Mercury, Manganese, Lead, Cadmium, Selenium
  )

# View the first few rows
head(final_data)
data <- final_data %>%
  drop_na()
summary(data)

lm.fit.alt <- lm(log(ALT) ~ log(Mercury) + log(Manganese) + log(Lead) + 
                   log(Cadmium) + log(Selenium) +
                   I(Race_Ethnicity=="Non-Hispanic Black") + Male + 
                   poly(Age, 2) + Poverty_Ratio + Alcohol_User + Smoker +
                   HDL_C + Total_Cholesterol + Energy_Intake_kcal +
                   Hypertension + Diabetes, data=data)
summary(lm.fit.alt)
hist(residuals(lm.fit.alt), breaks=50)
car::residualPlots(lm.fit.alt)

lm.fit.alt.2 <- lm(log(ALT) ~ Mercury + Manganese + Lead + Cadmium + Selenium +
                   Race_Ethnicity + Education_Level + Sedentary_Cat+ Male + 
                   poly(Age, 2) + Poverty_Ratio + Alcohol_User + Smoker +
                   HDL_C + Total_Cholesterol + Energy_Intake_kcal +
                   Hypertension + Diabetes, data=data)
summary(lm.fit.alt.2)

lm.fit.ast <- lm(log(AST) ~ Mercury + Manganese + Lead + Cadmium + Selenium +
                   I(Race_Ethnicity=="Non-Hispanic Black") + Male + 
                   poly(Age, 2) + Poverty_Ratio + Alcohol_User + Smoker +
                   HDL_C + Total_Cholesterol + Energy_Intake_kcal +
                   Hypertension + Diabetes, data=data)
summary(lm.fit.ast)
hist(residuals(lm.fit.ast), breaks=50)

lm.fit.ast.2 <- lm(log(AST) ~ Mercury + Manganese + Lead + Cadmium + Selenium +
                     Race_Ethnicity + Education_Level + Sedentary_Cat+ Male + 
                     poly(Age, 2) + Poverty_Ratio + Alcohol_User + Smoker +
                     HDL_C + Total_Cholesterol + Energy_Intake_kcal +
                     Hypertension + Diabetes, data=data)
summary(lm.fit.ast.2)

library(mgcv)
gam.fit.alt <- gam(log(ALT) ~ s(log(Mercury)) + s(Manganese) + s(Lead) + s(Cadmium) + s(Selenium) +
                    I(Race_Ethnicity=="Non-Hispanic Black") + Male + 
                    s(Age) + Poverty_Ratio + Alcohol_User + Smoker +
                    HDL_C + Total_Cholesterol + Energy_Intake_kcal +
                    Hypertension + Diabetes, data=data)
summary(gam.fit.alt)
plot(gam.fit.alt)
hist(residuals(gam.fit.alt), breaks=50)

gam.fit.ast <- gam(log(ALT) ~ s(log(Mercury)) + s(Manganese) + s(Lead) + s(Cadmium) + s(Selenium) +
                     I(Race_Ethnicity=="Non-Hispanic Black") + Male + 
                     s(Age) + Poverty_Ratio + Alcohol_User + Smoker +
                     HDL_C + Total_Cholesterol + Energy_Intake_kcal +
                     Hypertension + Diabetes, data=data)
summary(gam.fit.ast)
plot(gam.fit.ast)
hist(residuals(gam.fit.ast), breaks=50)

save(data, file="data/NHANES_2017_2020_Liver.Rdata")
