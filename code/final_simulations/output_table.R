# ---------------------------------------------------------
# MASTER SIMULATION SUMMARY SCRIPT
# Aggregates Type I Error, Power, and MSE across all scenarios
# ---------------------------------------------------------

rm(list=ls())
library(tidyverse)
library(knitr)

# Set working directory to your repository
setwd("~/GitHub/Pairwise_Additive_GP_Tests")

# 1. Define the simulation files based on the outputs directory
sim_files <- c(
  "SimResults_AnyInteraction.Rdata",
  "SimResults_AnyInteraction_n300.Rdata",
  "SimResults_Scenario2_JointAssociation.Rdata",
  "SimResults_Scenario2_JointAssociation_n300.Rdata",
  "SimResults_Scenario2_Trig_p10.Rdata",
  "SimResults_Scenario2_Trig_p10_n300.Rdata",
  "SimResults_Scenario5_PureInteraction_n100_p10.Rdata",
  "SimResults_Scenario5_PureInteraction_n300_p10.Rdata",
  "SimResults_Scenario5_TestInteraction_n100_p10.Rdata",
  "SimResults_Scenario5_TestInteraction_n300_p10.Rdata"
)

# Set nominal alpha for hypothesis testing (0.10)
alpha_thresh <- 0.10 

# Helper function to map old method names to requested labels
rename_methods <- function(name) {
  # Strip suffixes
  name <- gsub("_pval|_MSE|_PIP", "", name)
  name <- gsub("_", " ", name)
  
  # Map to new names
  if (name == "Full GP") return("GKM (Maity)")
  if (name == "GAM Smooth") return("GAM")
  if (name == "Old Pairwise GP") return("Pairwise GKM (Maity)")
  if (name == "Unified REML") return("Pairwise GP")
  return(name)
}

# Initialize empty dataframe to collect results
all_results <- data.frame()

# 2. Loop through each file and extract metrics
for (file in sim_files) {
  file_path <- file.path("outputs", file)
  
  if (!file.exists(file_path)) {
    warning(sprintf("File not found: %s", file_path))
    next
  }
  
  # Clear results_list to prevent carry-over from previous iterations
  if (exists("results_list")) rm(results_list)
  
  # Load the Rdata file
  load(file_path)
  
  # --- FORMAT ADAPTER FOR 'ANY INTERACTION' SIMULATIONS ---
  if (grepl("AnyInteraction", file)) {
    # If the workspace contains the standalone null/alt objects, format them
    if (exists("results_null") && exists("results_alt")) {
      
      # Rename 'p_value' to include '_pval' so grep catches it downstream
      if ("p_value" %in% colnames(results_null)) {
        colnames(results_null)[colnames(results_null) == "p_value"] <- "Global_Interaction_Test_pval"
      }
      if ("p_value" %in% colnames(results_alt)) {
        colnames(results_alt)[colnames(results_alt) == "p_value"] <- "Global_Interaction_Test_pval"
      }
      
      # Package into the expected list structure
      results_list <- list(
        "null" = results_null,
        "alt"  = results_alt
      )
      
      # Clean up the global variables to keep memory tidy
      rm(results_null, results_alt)
    }
  }
  # --------------------------------------------------------
  
  # Ensure the list loaded successfully before proceeding
  if (!exists("results_list")) {
    warning(sprintf("Could not find 'results_list' or valid interaction data in: %s", file))
    next
  }
  
  # Create a clean base name for the table
  base_scenario <- gsub("SimResults_|.Rdata", "", file)
  
  # Loop through sub-scenarios
  for (scen_name in names(results_list)) {
    res_df <- results_list[[scen_name]]
    
    # Determine if this sub-scenario represents Type I Error or Power
    is_null <- grepl("null|type_i", tolower(scen_name))
    stat_type <- ifelse(is_null, "Type I Error", "Power")
    mse_type <- ifelse(is_null, "MSE (Null)", "MSE (Alt)")
    
    # Strip null/alt identifiers to find the underlying condition
    condition <- gsub("(?i)type_i_error|power|null|alternative|alt", "", scen_name)
    condition <- gsub("^_+|_+$", "", condition) # Clean up leftover underscores
    if (condition == "") condition <- "Overall"
    
    # Identify testing, MSE, and PIP columns
    pval_cols <- grep("pval", colnames(res_df), value = TRUE)
    mse_cols <- grep("MSE", colnames(res_df), value = TRUE)
    pip_cols <- grep("PIP", colnames(res_df), value = TRUE)
    
    # Process Hypothesis Testing (Rejection Rates)
    for (m in pval_cols) {
      if (!all(is.na(res_df[[m]]))) {
        rate <- mean(res_df[[m]] < alpha_thresh, na.rm = TRUE)
        clean_name <- rename_methods(m)
        
        all_results <- rbind(all_results, data.frame(
          Scenario = base_scenario,
          Condition = condition,
          Method = clean_name,
          Metric = stat_type,
          Value = rate
        ))
      }
    }
    
    # Process BKMR Hypothesis Testing (PIP Cutoffs)
    for (m in pip_cols) {
      if (!all(is.na(res_df[[m]]))) {
        # Calculate rates for both cutoffs
        rate_50 <- mean(res_df[[m]] > 0.50, na.rm = TRUE)
        rate_90 <- mean(res_df[[m]] > 0.90, na.rm = TRUE)
        
        base_name <- rename_methods(m)
        
        # Add row for PIP > 0.5
        all_results <- rbind(all_results, data.frame(
          Scenario = base_scenario,
          Condition = condition,
          Method = paste0(base_name, " (PIP > 0.5)"),
          Metric = stat_type,
          Value = rate_50
        ))
        
        # Add row for PIP > 0.9
        all_results <- rbind(all_results, data.frame(
          Scenario = base_scenario,
          Condition = condition,
          Method = paste0(base_name, " (PIP > 0.9)"),
          Metric = stat_type,
          Value = rate_90
        ))
      }
    }
    
    # Process Predictive Accuracy (Mean MSE)
    for (m in mse_cols) {
      if (!all(is.na(res_df[[m]]))) {
        avg_mse <- mean(res_df[[m]], na.rm = TRUE)
        clean_name <- rename_methods(m)
        
        all_results <- rbind(all_results, data.frame(
          Scenario = base_scenario,
          Condition = condition,
          Method = clean_name,
          Metric = mse_type,
          Value = avg_mse
        ))
      }
    }
  }
}

# 3. Pivot and format the final table
final_table <- all_results %>%
  pivot_wider(names_from = Metric, values_from = Value) %>%
  arrange(Scenario, Condition, Method)

# Ensure columns are ordered logically
col_order <- c("Scenario", "Condition", "Method", "Type I Error", "Power", "MSE (Null)", "MSE (Alt)")
final_table <- final_table[, intersect(col_order, colnames(final_table))]

# Print the formatted markdown table
cat("\n=========================================================\n")
cat(sprintf("   MASTER SIMULATION SUMMARY TABLE (Alpha = %.2f)\n", alpha_thresh))
cat("=========================================================\n")
print(kable(final_table, digits = 3, format = "markdown", na = "-"))
cat("\n")

# Save it out as a CSV for easy sharing/viewing in Excel
write.csv(final_table, "outputs/Master_Simulation_Summary.csv", row.names = FALSE)