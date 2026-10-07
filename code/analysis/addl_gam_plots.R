# Additional GAM plots for Pairwise GP MS
# run after analysis_outputs_gam.R

# Utility functions
library(mgcv)
library(gratia)
library(ggplot2)
library(dplyr)
library(patchwork)

plot_custom_main_effect <- function(model, smooth_term, clean_title, x_label = NULL, show_rug = FALSE) {
  
  sm <- smooth_estimates(model, smooth = smooth_term) %>%
    add_confint()
  
  var_name <- gsub("^ti\\(|^s\\(|\\)$", "", smooth_term)
  if (is.null(x_label)) x_label <- clean_title
  
  p <- ggplot(sm, aes(x = .data[[var_name]], y = .estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "steelblue") +
    geom_line(linewidth = 1, color = "midnightblue") +
    labs(
      title = clean_title,
      x = x_label,
      y = "Partial Effect"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  # Optionally add the rug plot using the original model data
  if (show_rug) {
    p <- p + geom_rug(
      data = model$model, 
      aes(x = .data[[var_name]]), 
      inherit.aes = FALSE, # Prevents ggplot from looking for 'y = .estimate' in raw data
      alpha = 0.5, 
      sides = "b"
    )
  }
  
  return(p)
}

plot_custom_interaction <- function(model, smooth_term, clean_title, x_label = NULL, y_label = NULL, 
                                    mask_too_far = TRUE, show_points = FALSE, point_size = 1.5) {
  
  sm <- smooth_estimates(model, smooth = smooth_term)
  
  if (mask_too_far) {
    sm <- sm %>%
      add_too_far(model$model) %>%
      filter(!.too_far)
  }
  
  clean_vars <- gsub("^ti\\(|^te\\(|^s\\(|\\)$", "", smooth_term) 
  clean_vars <- gsub(" ", "", clean_vars)                         
  vars <- strsplit(clean_vars, ",")[[1]]                          
  
  var_x <- vars[1]
  var_y <- vars[2]
  
  if (is.null(x_label)) x_label <- var_x
  if (is.null(y_label)) y_label <- var_y
  
  p <- ggplot(sm, aes(x = .data[[var_x]], y = .data[[var_y]])) +
    geom_tile(aes(fill = .estimate)) +
    geom_contour(aes(z = .estimate), color = "black", alpha = 0.3, bins = 10) +
    
    scale_fill_gradient2(
      low = "midnightblue", 
      mid = "white", 
      high = "darkred", 
      midpoint = 0,
      name = "Effect"
    ) +
    
    labs(
      title = clean_title,
      x = x_label,
      y = y_label
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid = element_blank() 
    )
  
  if (show_points) {
    p <- p + geom_point(
      data = model$model, 
      aes(x = .data[[var_x]], y = .data[[var_y]]), 
      inherit.aes = FALSE, 
      alpha = 0.3, 
      size = point_size, # Now controlled by the function argument
      color = "black"
    )
  }
  
  return(p)
}
plot_custom_interaction_sig <- function(model, smooth_term, clean_title, x_label = NULL, y_label = NULL, 
                                        mask_too_far = TRUE, show_points = FALSE, point_size = 1.5) {
  
  # Extract smooth estimates (using gratia)
  sm <- smooth_estimates(model, smooth = smooth_term)
  
  # Apply data support masking
  if (mask_too_far) {
    sm <- sm %>%
      add_too_far(model$model) %>%
      filter(!.too_far)
  }
  
  # Calculate the Z-score (Estimate / Standard Error)
  sm <- sm %>%
    mutate(.z_score = .estimate / .se)
  
  # Extract variable names
  clean_vars <- gsub("^ti\\(|^te\\(|^s\\(|\\)$", "", smooth_term) 
  clean_vars <- gsub(" ", "", clean_vars)                         
  vars <- strsplit(clean_vars, ",")[[1]]                          
  
  var_x <- vars[1]
  var_y <- vars[2]
  
  if (is.null(x_label)) x_label <- var_x
  if (is.null(y_label)) y_label <- var_y
  
  # Build the significance plot
  p <- ggplot(sm, aes(x = .data[[var_x]], y = .data[[var_y]])) +
    geom_tile(aes(fill = .z_score)) +
    
    # Highlight the boundaries of significance at alpha = 0.05 (|Z| = 1.96)
    geom_contour(aes(z = .z_score), breaks = c(-1.96, 1.96), color = "black", linewidth = 1) +
    
    # Add a subtle dashed line where the estimated effect is exactly zero
    geom_contour(aes(z = .z_score), breaks = 0, color = "grey50", linetype = "dashed", linewidth = 0.5) +
    
    # Diverging color scale centered at 0
    scale_fill_gradient2(
      low = "#2166AC", 
      mid = "white", 
      high = "#B2182B", 
      midpoint = 0,
      limits = c(-4, 4),           # Cap limits so extreme Z-scores don't wash out the gradient
      oob = scales::squish,        # Requires the 'scales' package
      name = "Z-Score\n(Est/SE)"
    ) +
    
    labs(
      title = clean_title,
      subtitle = "Solid black contours bound regions significant at the 95% level (|Z| > 1.96)",
      x = x_label,
      y = y_label
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid = element_blank() 
    )
  
  # Add raw data points if requested
  if (show_points) {
    p <- p + geom_point(
      data = model$model, 
      aes(x = .data[[var_x]], y = .data[[var_y]]), 
      inherit.aes = FALSE, 
      alpha = 0.3, 
      size = point_size, 
      color = "black"
    )
  }
  
  return(p)
}

library(mgcv)
library(gratia)
library(ggplot2)
library(dplyr)

plot_custom_varying_effect <- function(model, smooth_term, clean_title, x_label = NULL, show_rug = FALSE) {
  
  # 1. Extract estimates and CI
  sm <- smooth_estimates(model, smooth = smooth_term) %>%
    add_confint()
  
  # 2. Extract the main variable name and the 'by' variable name
  # For "ti(Log_Manganese):Male":
  # This regex strips "^ti(", "^s(", and anything from "):" to the end
  var_name <- gsub("^ti\\(|^te\\(|^s\\(|\\):.*$", "", smooth_term) 
  
  # This regex extracts just the text after the colon
  by_var <- sub(".*:(.*)$", "\\1", smooth_term)
  
  if (is.null(x_label)) x_label <- clean_title
  
  # 3. Build the plot
  p <- ggplot(sm, aes(x = .data[[var_name]], y = .estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "steelblue") +
    geom_line(linewidth = 1, color = "midnightblue") +
    labs(
      title = clean_title,
      x = x_label,
      y = paste("Partial Effect", paste0("(", by_var, ")")) # e.g., "Partial Effect (Male)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  # 4. Smart Rug Plot
  if (show_rug) {
    rug_data <- model$model
    
    # If the modifying variable is a binary dummy column (like a 'Male' 0/1 column), 
    # we filter the rug to only show observations for that specific group.
    if (by_var %in% colnames(rug_data)) {
      # Assumes dummy variable format where the active group is 1
      rug_data <- rug_data %>% filter(.data[[by_var]] == 1)
    }
    
    p <- p + geom_rug(
      data = rug_data, 
      aes(x = .data[[var_name]]), 
      inherit.aes = FALSE, 
      alpha = 0.5, 
      sides = "b"
    )
  }
  
  return(p)
}

plot_partial_effect_1d <- function(model, data, x_var, mod_var, 
                                   title = NULL, x_label = NULL, legend_title = NULL,
                                   quantiles = c(0.1, 0.5, 0.9), 
                                   quantile_labels = NULL) {
  
  # 1. Auto-generate grammatically correct ordinal labels if none provided
  if (is.null(quantile_labels)) {
    percentages <- quantiles * 100
    suffixes <- rep("th", length(percentages))
    suffixes[percentages %% 10 == 1 & percentages %% 100 != 11] <- "st"
    suffixes[percentages %% 10 == 2 & percentages %% 100 != 12] <- "nd"
    suffixes[percentages %% 10 == 3 & percentages %% 100 != 13] <- "rd"
    quantile_labels <- paste0(percentages, suffixes)
  }
  
  # 2. Create a baseline dataframe 
  base_data <- data %>%
    summarise(across(everything(), \(x) if(is.numeric(x)) mean(x, na.rm = TRUE) else first(x)))
  
  # 3. Create the specific prediction grid dynamically based on provided quantiles
  x_seq <- seq(min(data[[x_var]], na.rm = TRUE), max(data[[x_var]], na.rm = TRUE), length.out = 100)
  mod_quantiles <- quantile(data[[mod_var]], probs = quantiles, na.rm = TRUE)
  
  grid_vals <- expand.grid(x_seq, mod_quantiles)
  colnames(grid_vals) <- c(x_var, mod_var)
  
  newdata <- base_data %>%
    slice(rep(1, nrow(grid_vals)))
  
  newdata[[x_var]] <- grid_vals[[x_var]]
  newdata[[mod_var]] <- grid_vals[[mod_var]]
  
  # 4. Extract the linear predictor matrix
  Xp <- predict(model, newdata = newdata, type = "lpmatrix")
  
  # 5. Isolate the target terms dynamically
  # Escape any periods in variable names for regex (e.g., Median.Income -> Median\\.Income)
  x_regex <- gsub("\\.", "\\\\.", x_var)
  mod_regex <- gsub("\\.", "\\\\.", mod_var)
  
  # Match the main effect
  idx_main <- grepl(paste0("^ti\\(", x_regex, "\\)"), colnames(Xp))
  
  # Match the interaction regardless of order: ti(x,mod) or ti(mod,x)
  idx_interaction <- grepl(paste0("^ti\\((", x_regex, ",", mod_regex, "|", mod_regex, ",", x_regex, ")\\)"), colnames(Xp))
  
  terms_to_keep <- idx_main | idx_interaction
  
  # Zero out all other nuisance covariates
  Xp[, !terms_to_keep] <- 0
  
  # 6. Calculate the isolated fit and SE
  grid_vals$fit <- as.vector(Xp %*% coef(model))
  grid_vals$se  <- sqrt(rowSums((Xp %*% vcov(model)) * Xp))
  
  # 7. Apply default labels if not provided
  if (is.null(title)) title <- paste("Estimated Partial Effect of", x_var, "Modified by", mod_var)
  if (is.null(x_label)) x_label <- x_var
  if (is.null(legend_title)) legend_title <- paste(mod_var, "\n(Quantiles)")
  
  # 8. Setup dynamic linetypes to accommodate any number of quantiles
  base_linetypes <- c("solid", "dashed", "longdash", "dotdash", "dotted", "twodash")
  use_linetypes <- rep(base_linetypes, length.out = length(quantiles))
  
  # 9. Plot the isolated partial effect
  p <- ggplot(grid_vals, aes(
    x = .data[[x_var]], 
    y = fit, 
    color = as.factor(.data[[mod_var]]), 
    fill = as.factor(.data[[mod_var]]),
    linetype = as.factor(.data[[mod_var]]) 
  )) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "darkgrey", linewidth = 0.8) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se), alpha = 0.2, color = NA) +
    scale_color_viridis_d(
      name = legend_title, 
      labels = quantile_labels,
      option = "cividis"
    ) +
    scale_fill_viridis_d(
      name = legend_title, 
      labels = quantile_labels,
      option = "cividis"
    ) +
    scale_linetype_manual(
      name = legend_title, 
      labels = quantile_labels,
      values = use_linetypes
    ) +
    labs(
      title = title,
      x = x_label,
      y = "Partial Effect (Centered)"
    ) +
    theme_minimal() +
    theme(
      legend.key.width = unit(1.5, "cm") 
    )
  
  return(p)
}
################### NC Well #################
load("outputs/ncwell_obesity_analysis_gam.Rdata")
data <- read.csv("data//NC_tract_all_outcomes.csv")
data.gam <- data %>%
  mutate(Log_As = log(Arsenic.Mean_avg),
         Log_Pb = log(Lead.Mean_avg),
         Log_Cr = log(Chromium.Mean_avg),
         Log_Zn = log(Zinc.Mean_avg),
         Log_Ba = log(Barium.Mean_avg))
round(summary(ncwell.gam)$s.table[,4], 3)
library(dplyr)
library(ggplot2)
library(mgcv)


################################################################################
#                    demographic line plots                                    #
################################################################################

as_income <- plot_partial_effect_1d(ncwell.gam, data.gam, x_var="Log_As",
                                    x_label="ln(As)", mod_var="Median.Income",
                       title="Estimated Partial Effect of Arsenic Modified by Median Income",
                       legend_title = "Median Income (Quantiles)") +
  ylab("Partial Effect")

ggsave("outputs/figs/ncwell_as_income_line.png", as_income)

as_pb <- plot_partial_effect_1d(ncwell.gam, data.gam, x_var="Log_As", mod_var="Log_Pb",
                       title="Estimated Partial Effect of Arsenic Modified by Lead",
                       legend_title = "Median Income (Quantiles)") #OK, Nothing. Good.

plot(data.gam$Log_As, data.gam$Log_Pb)

ggsave("outputs/figs/ncwell_as_pb_line.png", as_pb)


plot_partial_effect_1d(ncwell.gam, data.gam, x_var="Black.Pct", mod_var="Median.Income",
                       title="Estimated Partial Effect of Black Pct Modified by Median Income",
                       legend_title = "Median Income (Quantiles)")

plot(data.gam$Black.Pct, data.gam$Median.Income)

income.quant.90 <- quantile(data.gam$Median.Income, probs=.9)
black.income.90 <- data.gam %>%
  filter(Median.Income > income.quant.90) %>%
  select(Black.Pct)
hist(as.matrix(black.income.90))

hisp.age.int <- plot_partial_effect_1d(ncwell.gam, data.gam, x_var="Hisp.Pct", mod_var="Age",
                       title="Estimated Partial Effect of Hispanic Pct on Obesity Pct, Modified by Age",
                       legend_title = "Age (Quantiles)", quantiles=c(.1, .25, .5))
ggsave(filename="outputs/figs/ncwell_hisp_age.png", hisp.age.int)

plot(data.gam$Hisp.Pct, data.gam$Age)

age.quant.90 <- quantile(data.gam$Age, probs=.9)
hisp.age.90 <- data.gam %>%
  filter(Age > age.quant.90) %>%
  select(Hisp.Pct)
hist(as.matrix(hisp.age.90))

age.quant.50 <- quantile(data.gam$Age, probs=.5)
hisp.age.50 <- data.gam %>%
  filter(Age > age.quant.50) %>%
  select(Hisp.Pct)
hist(as.matrix(hisp.age.50))

age.quant.05 <- quantile(data.gam$Age, probs=.05)
age.quant.15 <- quantile(data.gam$Age, probs=.15)
hisp.age.10 <- data.gam %>%
  filter(Age > age.quant.05,
         Age < age.quant.15) %>%
  select(Hisp.Pct)
hist(as.matrix(hisp.age.10))

age.quant.2 <- quantile(data.gam$Age, probs=.2)
age.quant.3 <- quantile(data.gam$Age, probs=.3)
hisp.age.25 <- data.gam %>%
  filter(Age > age.quant.2,
         Age < age.quant.3) %>%
  select(Hisp.Pct)
hist(as.matrix(hisp.age.25))

age.quant.45 <- quantile(data.gam$Age, probs=.45)
age.quant.55 <- quantile(data.gam$Age, probs=.55)
hisp.age.5 <- data.gam %>%
  filter(Age > age.quant.45,
         Age < age.quant.55) %>%
  select(Hisp.Pct)
hist(as.matrix(hisp.age.5))


# Main eff's: nonlinear?
income.main <- plot_custom_main_effect(ncwell.gam, "ti(Median.Income)", "Estimated Partial Effect of Median Income on Obesity Pct", "Median Income") # Great
age.main <- plot_custom_main_effect(ncwell.gam, "ti(Age)","Estimated Partial Effect of Median Age on Obesity Pct", "Age") # Great

ggsave(filename="outputs/figs/ncwell_income.png", income.main)
ggsave(filename="outputs/figs/ncwell_age.png", age.main)
################################################################################
#                        Main and interaction effect plots                     #
################################################################################

# Generate the individual plots with your clean labels
p_as <- plot_custom_main_effect(ncwell.gam, "ti(Log_As)", "Arsenic", 
                                x_label="ln(Arsenic)", show_rug = TRUE) +
  ylim(-2, 1.5)
p_pb <- plot_custom_main_effect(ncwell.gam, "ti(Log_Pb)", "Lead",
                                x_label="ln(Lead)", show_rug = TRUE) +
  ylim(-2, 1.5)
# p_cr <- plot_custom_main_effect(ncwell.gam, "ti(Log_Cr)", "Chromium", x_label="ln(Chromium)", show_rug = TRUE)
# p_zn <- plot_custom_main_effect(ncwell.gam, "ti(Log_Zn)", "Zinc", x_label="ln(Zinc)", show_rug = TRUE)
# p_ba <- plot_custom_main_effect(ncwell.gam, "ti(Log_Ba)", "Barium", x_label="ln(Barium)", show_rug = TRUE)
# Combine them into a single layout using patchwork
p_main <- wrap_plots(p_as, p_pb)
ggsave("outputs/figs/ncwell_main.png", plot=p_main, width=7, height=5)
# Generate the individual 2D plots


p_as_income <- plot_custom_interaction(
  model = ncwell.gam, 
  smooth_term = "ti(Median.Income,Log_As)", 
  clean_title = "Income & Arsenic Interaction",
  x_label = "Median Income",
  y_label = "Ln(Arsenic)",
  mask_too_far = FALSE,
  show_points = TRUE
)

# p_as_pb <- plot_custom_interaction(
#   model = ncwell.gam, 
#   smooth_term = "ti(Log_As,Log_Pb)", 
#   clean_title = "Lead & Arsenic Interaction",
#   x_label = "ln(Arsenic)",
#   y_label = "ln(Lead)",
#   mask_too_far = FALSE,
#   show_points = TRUE
# )

# Display side-by-side
# p_int <- p_age_cr | p_as_income
ggsave("outputs/figs/ncwell_int.png", plot=p_as_income, width=9, height=5)



# 
# p_age_cr_sig <- plot_custom_interaction_sig(
#   model = ncwell.gam, 
#   smooth_term = "ti(Age,Log_Cr)", 
#   clean_title = "Age & Chromium Interaction",
#   x_label = "Age",
#   y_label = "Ln(Chromium)",
#   mask_too_far = FALSE,
#   show_points = TRUE)

p_as_income_sig <- plot_custom_interaction_sig(
  model = ncwell.gam, 
  smooth_term = "ti(Median.Income,Log_As)", 
  clean_title = "Income & Arsenic Interaction",
  x_label = "Median Income",
  y_label = "Ln(Arsenic)",
  mask_too_far = FALSE,
  show_points = TRUE)

# p_age_cr_sig | p_as_income_sig
# ggsave("outputs/figs/int_age_cr_sig.png", plot=p_age_cr_sig, width=9, height=5)
ggsave("outputs/figs/int_as_income_sig.png", plot=p_as_income_sig, width=9, height=5)

#### NHANES Liver ####
gam.alt.nohg
p_mn_alt <- plot_custom_main_effect(gam.alt.nohg, "ti(Log_Manganese)", "Manganese", 
                                x_label="ln(Manganese)", show_rug = TRUE)
p_mn_male_alt <- plot_custom_varying_effect(gam.alt.nohg, "ti(Log_Manganese):Male", "Manganese, Male Only", 
                                x_label="ln(Manganese)", show_rug = TRUE)
p_pb_alt <- plot_custom_main_effect(gam.alt.nohg, "ti(Log_Lead)", "Lead", 
                                x_label="ln(Lead)", show_rug = TRUE)
p_cd_alt <- plot_custom_main_effect(gam.alt.nohg, "ti(Log_Cadmium)", "Cadmium", 
                                x_label="ln(Cadmium)", show_rug = TRUE)
p_se_alt <- plot_custom_main_effect(gam.alt.nohg, "ti(Log_Selenium)", "Selenium", 
                                x_label="ln(Selenium)", show_rug = TRUE)

p_univariate_nhanes <- wrap_plots(p_mn_alt, p_mn_male_alt, p_pb_alt, p_cd_alt, p_se_alt, ncol = 3)
ggsave("outputs/figs/nhanes_main.png", plot=p_univariate_nhanes, width=7, height=5)

draw(gam.alt.nohg, select=c(6,9))

p_mn_hdl <- plot_custom_interaction(
  model = gam.alt.nohg, 
  smooth_term = "ti(Log_Manganese,HDL_C)", 
  clean_title = "Manganese & HDL-C Interaction",
  x_label = "Manganese",
  y_label = "HDL-C",
  mask_too_far = FALSE,
  show_points = TRUE
)

p_cd_hdl <- plot_custom_interaction(
  model = gam.alt.nohg, 
  smooth_term = "ti(Log_Cadmium,HDL_C)", 
  clean_title = "Cadmium & HDL-C Interaction",
  x_label = "Cadmium",
  y_label = "HDL-C",
  mask_too_far = FALSE,
  show_points = TRUE
)

p_int_nhanes <- p_mn_hdl | p_cd_hdl
ggsave("outputs/figs/nhanes_int.png", plot=p_int_nhanes, width=9, height=3)

p_mn_hdl_sig <- plot_custom_interaction_sig(
  model = gam.alt.nohg, 
  smooth_term = "ti(Log_Manganese,HDL_C)", 
  clean_title = "Manganese & HDL-C Interaction",
  x_label = "Manganese",
  y_label = "HDL-C",
  mask_too_far = FALSE,
  show_points = TRUE
)

p_cd_hdl_sig <- plot_custom_interaction_sig(
  model = gam.alt.nohg, 
  smooth_term = "ti(Log_Cadmium,HDL_C)", 
  clean_title = "Cadmium & HDL-C Interaction",
  x_label = "Cadmium",
  y_label = "HDL-C",
  mask_too_far = FALSE,
  show_points = TRUE
)

p_int_nhanes_sig <- p_mn_hdl_sig | p_cd_hdl_sig
ggsave("outputs/figs/nhanes_int.png", plot=p_int_nhanes, width=9, height=3)

ggsave("outputs/figs/int_mn_hdl_sig.png", plot=p_mn_hdl_sig, width=9, height=5)
ggsave("outputs/figs/int_cd_hdl_sig.png", plot=p_cd_hdl_sig, width=9, height=5)

############# NHANES Telomere ###############
p_pcb_126_telo <- plot_custom_main_effect(gam.telo, "s(log_PCB126)", "PCB 126", 
                                    x_label="ln(PCB 126)", show_rug = TRUE)
p_furan_1_telo <- plot_custom_main_effect(gam.telo, "s(log_Furan1)", "Furan 1", 
                                         x_label="ln(Furan 1)", show_rug = TRUE)
p_main_telo <- p_pcb_126_telo | p_furan_1_telo
ggsave("outputs/figs/telo_main.png", plot=p_main_telo, width=6, height=3)
