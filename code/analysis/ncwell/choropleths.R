# Load required libraries
library(tidyverse)
library(sf)
library(tigris)      # For downloading census tract shapefiles
library(viridis)     # For color scales
setwd("~/GitHub/Pairwise_Additive_GP_Tests")

data <- read.csv("data//NC_tract_all_outcomes.csv")
# 1. Download North Carolina Census Tract Geometries
# Set cb = TRUE for a generalized (lighter) boundary file, great for plotting
options(tigris_use_cache = TRUE)
nc_tracts <- tracts(state = "NC", cb = TRUE, year = 2010)

# 2. Prepare Your Raw Data (Example structure)
# Replace this block with your actual dataset. 
# Ensure your data has a 'GEOID' column that matches the 11-digit tract GEOIDs.
# set.seed(123)

my_raw_data <- data %>%
  select(tract_fips, obesity_pct, Arsenic.Mean_avg) %>%
  mutate(tract_fips=as.character(tract_fips))

# 3. Join the Data to the Geometries
# Strip the ACS prefix from the tigris GEOID to match your 11-digit FIPS
nc_tracts_clean <- nc_tracts %>%
  mutate(tract_fips = str_remove(GEO_ID, "1400000US"))

# Now join using your data's column name (assuming your column is named 'tract_fips')
map_data <- nc_tracts_clean %>%
  left_join(my_raw_data, by = c("tract_fips" = "tract_fips"))

# 4. Plot the Choropleth Map using ggplot2
nc_choropleth_obesity <- ggplot(data = map_data) +
  geom_sf(aes(fill = obesity_pct, color=after_scale(fill)), 
          # color = "grey", 
          linewidth=0.0001) + # Thin white borders
  scale_fill_viridis_c(
    option = "magma",    # "viridis", "magma", "plasma", or "cividis"
    na.value = "grey80", # Color for tracts with missing data
    name = "Obesity Pct"
  ) +
  labs(
    title = "Estimated Obesity Prevalence by Census Tract",
    subtitle = "North Carolina",
    caption = "Source: CDC PLACES"
  ) +
  theme_void() + # Removes axes, background grid
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "right"
  )

# 5. Display the map
print(nc_choropleth_obesity)

# Optional: Save for your Beamer slides / paper
ggsave("outputs/figs/Obesity_NC_Tract_Choropleth.png", plot = nc_choropleth_obesity, width = 8, height = 5, dpi = 300)



# 4. Plot the Choropleth Map using ggplot2
nc_choropleth_as <- ggplot(data = map_data) +
  geom_sf(aes(fill = log(Arsenic.Mean_avg), color=after_scale(fill)), 
          # color = "grey", 
          linewidth=0.0001) + # Thin white borders
  scale_fill_viridis_c(
    option = "magma",    # "viridis", "magma", "plasma", or "cividis"
    na.value = "grey80", # Color for tracts with missing data
    name = "iAs, ppb\n(log)"
  ) +
  labs(
    title = "Average Inorganic Arsenic (iAs) in Private Wells by Census Tract",
    subtitle = "North Carolina",
    caption = "Source: NCWell Database"
  ) +
  theme_void() + # Removes axes, background grid
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "right"
  )

# 5. Display the map
print(nc_choropleth_as)

# Optional: Save for your Beamer slides / paper
ggsave("outputs/figs/iAs_NC_Tract_Choropleth.png", plot = nc_choropleth_as, width = 8, height = 5, dpi = 300)
