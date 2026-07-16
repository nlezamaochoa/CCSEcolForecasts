############################################################
# Swordfish Habitat Suitability
# Daily Predictions → Grand Mean + Monthly Climatology
# Publication-Quality Maps
############################################################

rm(list = ls())

############################################################
# 1. Libraries
############################################################

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


############################################################
# Custom Color Palette
############################################################

custom_colors <- colorRampPalette(c(
  "#9b59b6",  # purple
  "#3498db",  # blue
  "#1abc9c",  # turquoise/green
  "#f1c40f",  # yellow
  "#e74c3c"   # red
))(100)


############################################################
# 2. Directories
############################################################

setwd("/Volumes/Triple_Bottom_Line/Data/EcoROMS/species_predictions/")

data_dir <- "swor/"
out_dir  <- "/Volumes/Triple_Bottom_Line/Nerea_working/Forecast_swor/FINAL/Gridded/6_Average/"

############################################################
# 3. List Files
############################################################dd

files <- list.files(
  data_dir,
  pattern = "mean\\.grd$",
  full.names = TRUE
)

files <- sort(files)

############################################################
# 4. Extract Dates
############################################################

# Expected format:
# prediction_day_2000-02-05.grd

dates <- str_extract(
  basename(files),
  "\\d{4}-\\d{2}-\\d{2}"
)

dates <- ymd(dates)

############################################################
# 5. Load Raster Stack
############################################################

r_stack <- rast(files)

############################################################
# 6. Convert Longitude
############################################################
# IMPORTANT:
# Converts 0–360 longitude to -180–180
############################################################

r_stack <- rotate(r_stack)

############################################################
# 7. Subset Dates
############################################################

keep <- which(
  dates >= as.Date("1993-01-01") &
    dates <= as.Date("2018-12-31")
)

r_stack <- r_stack[[keep]]
dates <- dates[keep]

time(r_stack) <- dates

############################################################
# 8. Compute Grand Mean
############################################################

mean_all_days <- mean(
  r_stack,
  na.rm = TRUE
)

############################################################
# 9. Monthly Climatology
############################################################

month_index <- month(dates)

monthly_mean <- tapp(
  r_stack,
  index = month_index,
  fun = mean,
  na.rm = TRUE
)

names(monthly_mean) <- month.abb

############################################################
# 10. Crop Region
############################################################
# Southern California Current
############################################################

ext_crop <- ext(  -130,
  -117,
  30,
  47
)

library(terra)

mean_all_days <- shift(mean_all_days, dx = -360)
monthly_mean  <- shift(monthly_mean,  dx = -360)

ext(mean_all_days)
ext(monthly_mean)

mean_crop <- crop(mean_all_days, ext_crop)

monthly_crop <- crop(monthly_mean, ext_crop)

############################################################
# 11. Convert to Data Frames
############################################################

# Grand mean
mean_df <- as.data.frame(
  mean_crop,
  xy = TRUE,
  na.rm = FALSE
)

colnames(mean_df)[3] <- "prediction"

# Monthly climatology
monthly_df <- as.data.frame(
  monthly_crop,
  xy = TRUE,
  na.rm = FALSE
)

monthly_long <- monthly_df %>%
  pivot_longer(
    cols = -c(x, y),
    names_to = "month",
    values_to = "prediction"
  )

monthly_long$month <- factor(
  monthly_long$month,
  levels = month.abb
)

############################################################
# 12. Coastline
############################################################

world <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)

world_crop <- st_crop(
  world,
  xmin = -130,
  xmax = -117,
  ymin = 30,
  ymax = 47
)

############################################################
# 13. Common Color Limits
############################################################
fill_limits <- c(0, 1)

############################################################
# 14. Publication Theme
############################################################

theme_map <- theme_minimal(base_size = 13) +
  
  theme(
    
    panel.grid.major = element_line(
      color = "gray85",
      linewidth = 0.2
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.background = element_rect(
      fill = "#dceaf4",
      color = NA
    ),
    
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 12
    ),
    
    axis.text = element_text(
      color = "black",
      size = 10
    ),
    
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    
    legend.position = "bottom",
    
    legend.title = element_text(
      face = "bold"
    ),
    
    legend.width = unit(8, "cm")
  )

############################################################
# 15. Grand Mean Map
############################################################

p_mean <- ggplot() +
  
  geom_tile(
    data = mean_df,
    aes(
      x = x,
      y = y,
      fill = prediction
    )
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
  )+
  
  coord_sf(
    xlim = c(-130, -117),
    ylim = c(30, 47),
    expand = FALSE
  ) +
  
  labs(
    title = "Swordfish Habitat Suitability",
    subtitle = "Mean old model",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_map

############################################################
# 16. Save Grand Mean Map
############################################################

ggsave(
  filename = file.path(
    out_dir,
    "Mean_All_Days_old.png"
  ),
  plot = p_mean,
  width = 9,
  height = 7,
  dpi = 500,
  bg = "white"
)

############################################################
# 17. Monthly Climatology Map
############################################################

p_monthly <- ggplot() +
  
  geom_tile(
    data = monthly_long,
    aes(
      x = x,
      y = y,
      fill = prediction
    )
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
  
  facet_wrap(
    ~month,
    ncol = 4
  ) +
  
  labs(
    title = "Monthly Swordfish Habitat Climatology",
    subtitle = "Old model",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_map

############################################################
# 18. Save Monthly Climatology
############################################################

ggsave(
  filename = file.path(
    out_dir,
    "Monthly_Climatology_old.png"
  ),
  plot = p_monthly,
  width = 14,
  height = 10,
  dpi = 500,
  bg = "white"
)

############################################################
# 19. Optional Quick Plot
############################################################

print(p_mean)

print(p_monthly)

############################################################
# END


############################################################
