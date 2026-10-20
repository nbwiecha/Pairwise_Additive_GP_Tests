# NHANES 2017-2020 liver analysis: results tables, follow-up GAM, and figures
# (n = 5,980 after the smoking/drinking recode). Run after nhanes_liver_analysis.R
# and nhanes_liver_analysis_competitors.R.

rm(list = ls())
library(tidyverse)
library(mgcv)
library(gratia)
library(patchwork)
library(knitr)
setwd("~/GitHub/Pairwise_Additive_GP_Tests")
dir.create("outputs/figs", showWarnings = FALSE, recursive = TRUE)

# Plotting and summary-table helpers: load only the function definitions
load_functions <- function(file) {
  exprs <- parse(file)
  for (e in exprs) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) {
      eval(e, envir = globalenv())
    }
  }
}
load_functions("code/analysis/addl_gam_plots.R")
load_functions("code/analysis/data_summary_tables.R")

load("outputs/NHANES_Liver_Analysis.Rdata")              # PAGP results, data, S, Y_ALT
pagp <- ALT_analysis_adj
load("outputs/NHANES_Liver_Analysis_Competitors.Rdata")  # compet_analysis_alt, gam_alt
cat("n =", nrow(data), "\n")

metals <- c("Log_Mercury", "Log_Manganese", "Log_Lead", "Log_Cadmium", "Log_Selenium")

# ---------------------------------------------------------
# 1. Per-metal p-values (main-text table)
# ---------------------------------------------------------
gam_add <- summary(gam_alt)$s.table
liver_table <- data.frame(
  Metal = c("Mercury", "Manganese", "Lead", "Cadmium", "Selenium"),
  PAGP = pagp$var.importance.pvals[metals],
  OLS  = compet_analysis_alt$OLS[metals, 4],
  GKM  = compet_analysis_alt$Previous_KMR[metals],
  GAM  = gam_add[paste0("s(log(", c("Mercury", "Manganese", "Lead", "Cadmium", "Selenium"), "))"), 4],
  row.names = NULL
)
print(liver_table, digits = 3)
write.csv(liver_table, "outputs/nhanes_liver_metal_pvals.csv", row.names = FALSE)

# Joint tests of the five metals
cat("PAGP joint test of metals: p =", signif(ALT_metals$p_value, 3), "\n")
S_df <- as.data.frame(S)
ols_full <- lm(Y_ALT ~ ., data = S_df)
ols_red  <- lm(Y_ALT ~ ., data = S_df[, setdiff(colnames(S_df), metals)])
ols_joint_p <- anova(ols_red, ols_full)$`Pr(>F)`[2]
cat("OLS joint test of metals: p =", signif(ols_joint_p, 3), "\n")
cat("PAGP test for any interaction: p =", signif(ALT_int_test$p_value, 3), "\n")

# ---------------------------------------------------------
# 2. Significant PAGP interactions (supplement table)
# ---------------------------------------------------------
int_tab <- pagp$interaction.pvals %>%
  rownames_to_column("Interaction") %>%
  rename(p = `p-value`) %>%
  filter(p < 0.05) %>%
  arrange(p) %>%
  mutate(involves_metal = grepl("Log_", Interaction))
print(int_tab, digits = 3)
cat("Interaction tests run:", nrow(pagp$interaction.pvals), "\n")
write.csv(int_tab, "outputs/nhanes_liver_sig_interactions.csv", row.names = FALSE)

# ---------------------------------------------------------
# 3. Follow-up GAM: metals and PAGP-significant interactions, plus all covariates
# ---------------------------------------------------------
data.gam <- data %>%
  mutate(Log_Mercury = log(Mercury), Log_Manganese = log(Manganese), Log_Lead = log(Lead),
         Log_Cadmium = log(Cadmium), Log_Selenium = log(Selenium),
         Male = as.numeric(Male), Black = as.numeric(Race_Ethnicity == "Non-Hispanic Black"),
         Alcohol_User = as.numeric(Alcohol_User), Smoker = as.numeric(Smoker),
         Hypertension = as.numeric(Hypertension), Diabetes = as.numeric(Diabetes))

form_alt <- "log(ALT) ~ ti(Log_Mercury) + ti(Log_Manganese) + ti(Log_Lead) + ti(Log_Cadmium) +
ti(Log_Selenium) +
ti(Log_Mercury, Age) + ti(Log_Mercury, by = Black) +
ti(Log_Manganese, Log_Selenium) + ti(Log_Manganese, by = Male) + ti(Log_Manganese, Age) +
ti(Log_Manganese, HDL_C) +
ti(Log_Lead, Age) + ti(Log_Lead, by = Alcohol_User) + ti(Log_Lead, HDL_C) +
ti(Log_Cadmium, HDL_C) +
ti(Log_Selenium, HDL_C) + ti(Log_Selenium, Energy_Intake_kcal) + ti(Log_Selenium, by = Hypertension) +
ti(Age, by = Male) + ti(Total_Cholesterol, by = Male) + ti(Energy_Intake_kcal, by = Male) +
ti(Age, HDL_C) + ti(Age, Total_Cholesterol) +
ti(Age) + ti(HDL_C) + ti(Total_Cholesterol) + ti(Energy_Intake_kcal) + ti(Poverty_Ratio) +
Male + Black + Smoker + Alcohol_User + Hypertension + Diabetes"

gam_follow <- gam(formula(form_alt), method = "REML", data = data.gam)
gam_tab <- summary(gam_follow)$s.table
print(round(gam_tab, 4))
write.csv(gam_tab, "outputs/nhanes_liver_gam_followup_pvals.csv")
save(gam_follow, data.gam, form_alt, liver_table, int_tab, ols_joint_p,
     file = "outputs/nhanes_liver_followup.Rdata")

# ---------------------------------------------------------
# 4. Figures
# ---------------------------------------------------------
p_hg <- plot_custom_main_effect(gam_follow, "ti(Log_Mercury)", "Mercury", x_label = "ln(Mercury)", show_rug = TRUE)
p_mn <- plot_custom_main_effect(gam_follow, "ti(Log_Manganese)", "Manganese", x_label = "ln(Manganese)", show_rug = TRUE)
p_cd <- plot_custom_main_effect(gam_follow, "ti(Log_Cadmium)", "Cadmium", x_label = "ln(Cadmium)", show_rug = TRUE)
p_se <- plot_custom_main_effect(gam_follow, "ti(Log_Selenium)", "Selenium", x_label = "ln(Selenium)", show_rug = TRUE)
p_main <- wrap_plots(p_hg, p_mn, p_cd, p_se, ncol = 2)
ggsave("outputs/figs/nhanes_main.png", plot = p_main, width = 7, height = 5)

# Bonferroni-significant PAGP interactions for selenium and manganese, which are
# also significant in the follow-up GAM.
# Each line is the metal's main effect plus its interaction with the modifier.
p_se_htn <- plot_partial_effect_1d(gam_follow, data.gam, x_var = "Log_Selenium", mod_var = "Hypertension",
                                   title = "Selenium by Hypertension", x_label = "ln(Selenium)",
                                   legend_title = "Hypertension", quantiles = c(0, 1),
                                   quantile_labels = c("No", "Yes"))
p_mn_male <- plot_partial_effect_1d(gam_follow, data.gam, x_var = "Log_Manganese", mod_var = "Male",
                                    title = "Manganese by Sex", x_label = "ln(Manganese)",
                                    legend_title = "Sex", quantiles = c(0, 1),
                                    quantile_labels = c("Female", "Male"))
p_int <- (p_se_htn | p_mn_male) & labs(y = "Partial Effect") &
  theme(legend.position = "bottom", legend.title.position = "top", legend.key.width = unit(0.8, "cm"),
        plot.title = element_text(face = "bold"))
ggsave("outputs/figs/nhanes_int.png", plot = p_int, width = 8, height = 3.6)

# ---------------------------------------------------------
# 5. Summary statistics (supplement tables)
# ---------------------------------------------------------
data2 <- data %>%
  select(-Education_Level, -Sedentary_Cat, -SEQN, -AST, -ALT_AST_Ratio,
         -GGT, -Total_Bilirubin, -HSI)
tabs <- generate_summary_tables(data2)
write.csv(tabs[[1]], file = "outputs/liver_cont.csv")
write.csv(tabs[[2]], file = "outputs/liver_fact.csv")
