# ==================================================
# SWOR HABITAT SUITABILITY CLIMATOLOGY PIPELINE
# ==================================================
#
# DESCRIPTION
# This script processes daily swordfish habitat suitability
# predictions and produces publication-quality summaries:
#
# 1. Loads daily prediction rasters (1993–2018)
# 2. Builds a multi-temporal raster stack
# 3. Computes:
#      - Long-term grand mean habitat suitability
#      - Monthly climatology (seasonal cycle)
# 4. Crops outputs to the California Current region
# 5. Converts rasters to data frames for ggplot mapping
# 6. Generates publication-quality figures:
#      - Mean habitat suitability map
#      - Monthly climatology (facet maps)
#
# OUTPUTS
# - Mean_All_Days_0.05.png
# - Monthly_Climatology_0.05.png
# - Processed raster summaries (in memory / optional export)
#

# ==================================================

rm(list = ls())

# --------------------------------------------------
# LIBRARIES
# --------------------------------------------------

library(terra)
library(stringr)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(scales)
library(grid)

# --------------------------------------------------
# CUSTOM COLOR PALETTE
# --------------------------------------------------

custom_colors <- colorRampPalette(c(
  "#9b59b6",
           "#3498db",
           "#1abc9c",
           "#f1c40f",
           "#e74c3c"
))(100)

# --------------------------------------------------
# PATHS
# --------------------------------------------------

setwd("/Volumes/Triple_Bottom_Line/Nerea_working/Forecast_swor/FINAL/")

data_dir <- "Gridded/5_Predictions_ROMs_0.05/swor/"
out_dir  <- "Gridded/6_Average/"

# --------------------------------------------------
# LIST FILES
# --------------------------------------------------

files <- list.files(
  data_dir,
  pattern = "\\.tif$",
  full.names = TRUE
)

files <- sort(files)

# --------------------------------------------------
# EXTRACT DATES
# --------------------------------------------------

dates <- str_extract(basename(files),
                     "\\d{4}-\\d{2}-\\d{2}")

dates <- ymd(dates)

# --------------------------------------------------
# LOAD RASTER STACK
# --------------------------------------------------

r_stack <- rast(files)

# --------------------------------------------------
# LONGITUDE CORRECTION
# --------------------------------------------------

r_stack <- rotate(r_stack)

# --------------------------------------------------
# SUBSET TIME RANGE
# --------------------------------------------------

keep <- which(
  dates >= as.Date("1993-01-01") &
    dates <= as.Date("2018-12-31")
)

r_stack <- r_stack[[keep]]
dates <- dates[keep]

time(r_stack) <- dates

# --------------------------------------------------
# GRAND MEAN HABITAT SUITABILITY
# --------------------------------------------------

mean_all_days <- mean(r_stack, na.rm = TRUE)

# --------------------------------------------------
# MONTHLY CLIMATOLOGY
# --------------------------------------------------

month_index <- month(dates)

monthly_mean <- tapp(
  r_stack,
  index = month_index,
  fun = mean,
  na.rm = TRUE
)

names(monthly_mean) <- month.abb

# --------------------------------------------------
# CROP REGION (CALIFORNIA CURRENT)
# --------------------------------------------------

ext_crop <- ext(-130, -117, 30, 47)

mean_all_days <- shift(mean_all_days, dx = -360)
monthly_mean  <- shift(monthly_mean, dx = -360)

mean_crop <- crop(mean_all_days, ext_crop)
monthly_crop <- crop(monthly_mean, ext_crop)

# --------------------------------------------------
# CONVERT TO DATA FRAMES
# --------------------------------------------------

mean_df <- as.data.frame(mean_crop, xy = TRUE, na.rm = FALSE)
colnames(mean_df)[3] <- "prediction"

monthly_df <- as.data.frame(monthly_crop, xy = TRUE, na.rm = FALSE)

monthly_long <- monthly_df %>%
  pivot_longer(
    cols = -c(x, y),
    names_to = "month",
    values_to = "prediction"
  )

monthly_long$month <- factor(monthly_long$month,
                             levels = month.abb)

# --------------------------------------------------
# COASTLINE
# --------------------------------------------------

world <- ne_countries(scale = "medium", returnclass = "sf")

world_crop <- st_crop(
  world,
  xmin = -130, xmax = -117,
  ymin = 30, ymax = 47
)

# --------------------------------------------------
# COLOR SCALE
# --------------------------------------------------

fill_limits <- c(0, 1)

# --------------------------------------------------
# PLOT THEME
# --------------------------------------------------

theme_map <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "gray85", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#dceaf4", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black", size = 10),
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.width = unit(8, "cm")
  )

# --------------------------------------------------
# GRAND MEAN MAP
# --------------------------------------------------

p_mean <- ggplot() +
  geom_tile(
    data = mean_df,
    aes(x = x, y = y, fill = prediction)
  ) +
  geom_sf(
    data = world_crop,
    fill = "gray70",
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_gradientn(
    colours = custom_colors,
    limits = c(0, 1),
    oob = squish,
    na.value = "transparent",
    breaks = seq(0, 1, 0.2),
    name = "Habitat\nSuitability"
  ) +
  coord_sf(
    xlim = c(-130, -117),
    ylim = c(30, 47),
    expand = FALSE
  ) +
  labs(
    title = "Swordfish Habitat Suitability",
    subtitle = "Mean Across All Days (1993–2018)",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_map

# --------------------------------------------------
# SAVE GRAND MEAN MAP
# --------------------------------------------------

ggsave(
  filename = file.path(out_dir, "Mean_All_Days_0.05.png"),
  plot = p_mean,
  width = 9,
  height = 7,
  dpi = 500,
  bg = "white"
)

# --------------------------------------------------
# MONTHLY CLIMATOLOGY MAP
# --------------------------------------------------

p_monthly <- ggplot() +
  geom_tile(
    data = monthly_long,
    aes(x = x, y = y, fill = prediction)
  ) +
  geom_sf(
    data = world_crop,
    fill = "gray70",
    color = "black",
    linewidth = 0.15
  ) +
  scale_fill_gradientn(
    colours = custom_colors,
    limits = c(0, 1),
    oob = squish,
    na.value = "transparent",
    breaks = seq(0, 1, 0.2),
    name = "Suitability"
  ) +
  coord_sf(
    xlim = c(-130, -117),
    ylim = c(30, 47),
    expand = FALSE
  ) +
  facet_wrap(~month, ncol = 4) +
  labs(
    title = "Monthly Swordfish Habitat Climatology",
    subtitle = "1993–2018 Mean",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_map

# --------------------------------------------------
# SAVE MONTHLY CLIMATOLOGY
# --------------------------------------------------

ggsave(
  filename = file.path(out_dir, "Monthly_Climatology_0.05.png"),
  plot = p_monthly,
  width = 14,
  height = 10,
  dpi = 500,
  bg = "white"
)

# --------------------------------------------------
# OPTIONAL DISPLAY
# --------------------------------------------------

print(p_mean)
print(p_monthly)

# ==================================================
# END
# ==================================================