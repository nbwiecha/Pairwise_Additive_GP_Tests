# Analysis Table 1's

rm(list=ls())

library(tidyverse)
library(knitr)

setwd("~/GitHub/Pairwise_Additive_GP_Tests")

library(tidyverse)
library(knitr)

generate_summary_tables <- function(df) {
  
  # Helper functions to separate variable types
  is_cat_or_bin <- function(x) {
    is.character(x) || is.factor(x) || is.logical(x) || (is.numeric(x) && length(unique(na.omit(x))) <= 2)
  }
  
  is_strictly_cont <- function(x) {
    is.numeric(x) && length(unique(na.omit(x))) > 2
  }
  
  # ---------------------------------------------------------
  # 1. Continuous Variables Summary
  # ---------------------------------------------------------
  cont_vars <- df %>% select(where(is_strictly_cont))
  
  if (ncol(cont_vars) > 0) {
    cont_table <- cont_vars %>%
      pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
      group_by(Variable) %>%
      summarise(
        Mean = mean(Value, na.rm = TRUE),
        Median = median(Value, na.rm = TRUE),
        SD = sd(Value, na.rm = TRUE)
      ) %>%
      ungroup()
    
    cat("\n### Continuous Variables\n")
    print(kable(cont_table, digits = 3, format = "markdown"))
  } else {
    cont_table <- NULL
  }
  
  # ---------------------------------------------------------
  # 2. Categorical/Binary Variables Summary
  # ---------------------------------------------------------
  cat_vars <- df %>% select(where(is_cat_or_bin))
  
  if (ncol(cat_vars) > 0) {
    cat_table <- cat_vars %>%
      mutate(across(everything(), as.character)) %>%
      pivot_longer(cols = everything(), names_to = "Variable", values_to = "Category") %>%
      filter(!is.na(Category)) %>% 
      group_by(Variable, Category) %>%
      summarise(Count = n(), .groups = "drop") %>%
      # Group by variable again to calculate the percentage of the valid total
      group_by(Variable) %>%
      mutate(
        Percent = (Count / sum(Count)) * 100,
        `Percent (%)` = sprintf("%.1f%%", Percent) # Format with a percent sign and 1 decimal
      ) %>%
      select(-Percent) %>% # Drop the raw numeric column in favor of the formatted string
      ungroup() %>%
      arrange(Variable, Category)
    
    cat("\n### Categorical and Binary Variables\n")
    print(kable(cat_table, format = "markdown", align = c("l", "l", "r", "r")))
  } else {
    cat_table <- NULL
  }
  
  # Return the underlying dataframes invisibly for further use if needed
  invisible(list(continuous = cont_table, categorical = cat_table))
}

# ==========================================
# Example Usage:
# ==========================================
# dummy_data <- data.frame(
#   Age = rnorm(100, 50, 10),
#   BMI = rnorm(100, 25, 4),
#   Sex = sample(c("Male", "Female"), 100, replace = TRUE),
#   Treatment_Group = sample(c(0, 1), 100, replace = TRUE) # Numeric binary
# )
# 
# out <- generate_summary_tables(dummy_data)

# NHANES 2017-2020 liver data: summary tables are produced by
# code/analysis/nhanes_2017_2020/nhanes_liver_followup.R

###################### NC Well ##########################
rm(list=ls())
data <- read.csv("data//NC_tract_all_outcomes.csv")
head(data)
data2 <- data %>%
  select(-tract_fips, -diabetes_pct, -lat, -lon,-stroke_pct, -cancer_pct,
         -coronary_pct, -depression_pct, -NatAm.Pct, -numtests)

tabs_ncwell <- generate_summary_tables(data2)

write.csv(tabs_ncwell[[1]], file="outputs//ncwell_cont.csv")
write.csv(tabs_ncwell[[2]], file="outputs//ncwell_fact.csv")
