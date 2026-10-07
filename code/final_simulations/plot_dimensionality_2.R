# ---------------------------------------------------------
# VISUALIZE DIMENSIONALITY SIMULATION RESULTS
# Generates summary tables and ggplot line charts for Testing and MSE
# ---------------------------------------------------------

rm(list=ls())
library(tidyverse)
library(ggplot2)
library(knitr)
setwd("~/GitHub/Pairwise_Additive_GP_Tests")

# =========================================================
# 1. LOAD AND AGGREGATE DATA
# =========================================================
load("outputs/SimResults_Scenario1_Dimensionality_n100.Rdata")

alpha_thresh <- 0.10
summary_df_test <- data.frame()
summary_df_mse <- data.frame()

# Helper function to map old method names to new requested labels
rename_methods <- function(name) {
  if (name == "Full GP") return("GKM")
  if (name == "GAM Smooth") return("GAM")
  if (name == "Old Pairwise GP") return("Pairwise GKM")
  if (name == "Unified REML") return("PAGP")
  return(name) # Return original name (e.g., OLS, BKMR) if not in the list
}

# Loop through the list of results
for (scenario_name in names(results_list)) {
  
  res_df <- results_list[[scenario_name]]
  p_val <- res_df$p_dim[1]
  
  # Identify if this scenario evaluated Power or Type I Error
  sim_type <- ifelse(grepl("Type_I_Error", scenario_name), "Type I Error", "Power")
  
  # Separate P-value and MSE columns
  pval_cols <- grep("pval", colnames(res_df), value = TRUE)
  mse_cols <- grep("MSE", colnames(res_df), value = TRUE)
  
  # Calculate Rejection Rates for Hypothesis Testing
  for (m in pval_cols) {
    if (!all(is.na(res_df[[m]]))) {
      rate <- mean(res_df[[m]] < alpha_thresh, na.rm = TRUE)
      
      # Clean and rename method for the plot
      clean_name <- gsub("_pval", "", m)
      clean_name <- gsub("_", " ", clean_name)
      clean_name <- rename_methods(clean_name)
      
      summary_df_test <- rbind(summary_df_test, data.frame(
        p = p_val,
        Method = clean_name,
        Type = sim_type,
        Rate = rate
      ))
    }
  }
  
  # Calculate Average Mean Squared Error
  for (m in mse_cols) {
    if (!all(is.na(res_df[[m]]))) {
      mean_mse <- mean(res_df[[m]], na.rm = TRUE)
      
      # Clean and rename method for the plot
      clean_name <- gsub("_MSE", "", m)
      clean_name <- gsub("_", " ", clean_name)
      clean_name <- rename_methods(clean_name)
      
      summary_df_mse <- rbind(summary_df_mse, data.frame(
        p = p_val,
        Method = clean_name,
        Type = sim_type,
        MSE = mean_mse
      ))
    }
  }
}

# Factor the simulation type
summary_df_test$Type <- factor(summary_df_test$Type, levels = c("Type I Error", "Power"))
summary_df_mse$Type <- factor(summary_df_mse$Type, levels = c("Type I Error", "Power"))

# =========================================================
# 2. GENERATE GGPLOT2 LINE CHARTS
# =========================================================

# --- PLOT 1: Hypothesis Testing (Power / Type I Error) ---
p_test <- ggplot(filter(summary_df_test, p<25), aes(x = p, y = Rate, color = Method, linetype = Type, shape = Method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = alpha_thresh, linetype = "dashed", color = "black", linewidth = 0.8, alpha = 0.6) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_x_continuous(breaks = unique(summary_df_test$p)) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "right",
    legend.key.width = unit(3, "lines"), 
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Effect of Input Dimensionality on Test Performance",
    subtitle = sprintf("Sample Size n = 100, SNR = 1.0, Nominal Alpha = %.2f (Dashed Line)", alpha_thresh),
    x = "Total Input Variables (p)",
    y = "Empirical Rejection Rate",
    color = "Method",
    linetype = "Scenario",
    shape = "Method"
  )

# --- PLOT 2: Predictive Accuracy (MSE) ---
# Filter out the null simulations
mse_plot_data <- summary_df_mse %>% filter(Type == "Power")

p_mse <- ggplot(filter(mse_plot_data, p<25), aes(x = p, y = MSE, color = Method, shape = Method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = unique(mse_plot_data$p)) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "right",
    legend.key.width = unit(3, "lines"), 
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Effect of Input Dimensionality on Predictive Accuracy",
    subtitle = "Sample Size n = 100, SNR = 1.0 (Alternative Scenario Only)",
    x = "Total Input Variables (p)",
    y = "Mean Squared Error (MSE)",
    color = "Method",
    shape = "Method"
  )

# =========================================================
# 3. SAVE OUTPUTS
# =========================================================
ggsave("outputs/Plot_Dimensionality_Testing_n100.png", plot = p_test, width = 11, height = 7, dpi = 300)
ggsave("outputs/Plot_Dimensionality_MSE_n100.png", plot = p_mse, width = 11, height = 7, dpi = 300)

cat("Plots saved successfully to outputs/ directory.\n")
