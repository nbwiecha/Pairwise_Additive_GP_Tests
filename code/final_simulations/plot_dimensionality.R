# ---------------------------------------------------------
# VISUALIZE DIMENSIONALITY SIMULATION RESULTS
# Generates a summary table and a ggplot line chart
# ---------------------------------------------------------

rm(list=ls())
library(tidyverse)
library(ggplot2)
library(knitr)
setwd("~/GitHub/Pairwise_Additive_GP_Tests")
# =========================================================
# 1. LOAD AND AGGREGATE DATA
# =========================================================
# Load the saved payload from the simulation script
# This loads 'results_list' and 'final_output' into the environment
load("outputs/SimResults_Scenario1_Dimensionality_n100.Rdata")

alpha_thresh <- 0.10
summary_df <- data.frame()

# Loop through the list of results to extract the empirical rejection rates
for (scenario_name in names(results_list)) {
  
  res_df <- results_list[[scenario_name]]
  p_val <- res_df$p_dim[1]
  
  # Identify if this scenario evaluated Power or Type I Error
  sim_type <- ifelse(grepl("Type_I_Error", scenario_name), "Type I Error", "Power")
  
  # Identify frequentist methods (Ignore PIP and p_dim columns)
  methods <- grep("PIP|p_dim", colnames(res_df), invert = TRUE, value = TRUE)
  
  # Calculate rejection rates (p_value < alpha)
  for (m in methods) {
    if (m %in% colnames(res_df)) {
      # Only process if the method wasn't completely NA (e.g., if skipped)
      if (!all(is.na(res_df[[m]]))) {
        rate <- mean(res_df[[m]] < alpha_thresh, na.rm = TRUE)
        summary_df <- rbind(summary_df, data.frame(
          p = p_val,
          Method = m,
          Type = sim_type,
          Rate = rate
        ))
      }
    }
  }
}

# Factor the simulation type to ensure Type I Error plots/prints before Power
summary_df$Type <- factor(summary_df$Type, levels = c("Type I Error", "Power"))

# Clean up method names for the plot/table
summary_df$Method <- gsub("_", " ", summary_df$Method)

# =========================================================
# 2. GENERATE CLEAN SUMMARY TABLE
# =========================================================
# Pivot the dataframe to make a wide, readable table
table_df <- summary_df %>%
  pivot_wider(names_from = Method, values_from = Rate) %>%
  arrange(Type, p)

cat("\n=========================================================\n")
cat(sprintf("   SIMULATION RESULTS TABLE (Alpha = %.2f)\n", alpha_thresh))
cat("=========================================================\n")
print(kable(table_df, digits = 3, format = "markdown"))
cat("\n")

# =========================================================
# 3. GENERATE GGPLOT2 LINE CHART
# =========================================================
p_plot <- ggplot(summary_df, aes(x = p, y = Rate, color = Method, linetype = Type, shape = Method)) +
  # Add lines and points
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  
  # Add a nominal alpha reference line to easily check Type I Error control
  geom_hline(yintercept = alpha_thresh, linetype = "dashed", color = "black", linewidth = 0.8, alpha = 0.6) +
  
  # Formatting axes
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_x_continuous(breaks = unique(summary_df$p)) +
  
  # Clean publication theme
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "right",
    legend.key.width = unit(3, "lines"), # Make lines in legend wide enough to see the linetype
    panel.grid.minor = element_blank()
  ) +
  
  # Labels
  labs(
    title = "Effect of Input Dimensionality (p) on Test Performance",
    subtitle = sprintf("Sample Size n = 100, SNR = 1.0, Nominal Alpha = %.2f (Dashed Line)", alpha_thresh),
    x = "Total Input Variables (p)",
    y = "Empirical Rejection Rate",
    color = "Testing Method",
    linetype = "Simulation Type",
    shape = "Testing Method"
  )

# Save the high-resolution plot
ggsave("outputs/Plot_Dimensionality_n100.png", plot = p_plot, width = 11, height = 7, dpi = 300)
cat("Plot saved successfully to outputs/Plot_Dimensionality_n100.png\n")
