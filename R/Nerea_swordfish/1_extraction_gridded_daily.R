#==================================================
# SWOR ENVIRONMENTAL DATA EXTRACTION PIPELINE (MOM6)
# ==================================================
#
# DESCRIPTION
# This script extracts daily environmental variables from
# MOM6 ocean model NetCDF outputs at SWOR observation locations.
#
# It builds a matched observation–environment dataset by
# combining spatial and temporal extraction from large gridded
# ocean model outputs.
#
# MAIN STEPS
# 1. Load SWOR observation data and format dates/coordinates
# 2. Convert longitude to 0–360 to match ocean model grid
# 3. Extract daily environmental variables from NetCDF files:
#      - bbv
#      - ild
#      - sos
#      - ssh
#      - ssu_rotate
#      - ssv_rotate
#      - tos
# 4. Extract static bathymetry (deptho)
# 5. Compute spatial focal standard deviation (5x5 window):
#      - tos_sd
#      - ssh_sd
#      - deptho_sd
# 6. Compute lunar phase for each observation date
# 7. Compute derived metric:
#      - EKE (Eddy Kinetic Energy)
# 8. Clean dataset and save final CSV output
#
# OUTPUT
# swor_env_MOM6_gridded_daily.csv
# A merged dataset of SWOR observations + environmental covariates
#
# NOTES
# - Designed for large NetCDF time series (1993–2017)
# - Requires exact date matching between observations and model output
# - Uses terra for spatial extraction and raster handling
#
# ==================================================

library(terra)
library(dplyr)
library(lubridate)
library(lunar)

# --------------------------------------------------
# USER SETTINGS (EDIT FOR YOUR SYSTEM)
# --------------------------------------------------

project_dir <- "/YOUR/PROJECT/ROOT"
data_dir <- file.path(project_dir, "data")
output_dir <- file.path(project_dir, "output")

# --------------------------------------------------
# FILE PATHS
# --------------------------------------------------

nc_files <- list(
  bbv  = file.path(data_dir, "ocean_daily.19930101-20191231.bbv.nc"),
  ild  = file.path(data_dir, "ocean_daily.19930101-20191231.ild.nc"),
  sos  = file.path(data_dir, "ocean_daily.19930101-20191231.sos.nc"),
  ssh  = file.path(data_dir, "ocean_daily.19930101-20191231.ssh.nc"),
  ssu_rotate = file.path(data_dir, "ocean_daily.19930101-20191231.ssu_rotate.nc"),
  ssv_rotate = file.path(data_dir, "ocean_daily.19930101-20191231.ssv_rotate.nc"),
  tos  = file.path(data_dir, "ocean_daily.19930101-20191231.tos.nc"),
  deptho = file.path(data_dir, "deptho_field_RG.nc")
)

# --------------------------------------------------
# LOAD SWOR DATA
# --------------------------------------------------

swor <- read.csv(file.path(data_dir, "swor_dataframe.csv"))

swor$date <- as.Date(swor$dt, format = "%m/%d/%Y")

swor <- swor %>%
  filter(date >= as.Date("1993-01-01"),
         date <= as.Date("2017-01-31"))

# convert longitude to 0–360 system
swor$lon <- ifelse(swor$lon < 0, swor$lon + 360, swor$lon)

pts <- vect(swor, geom = c("lon", "lat"), crs = "EPSG:4326")

final_df <- swor

# --------------------------------------------------
# DAILY EXTRACTION FUNCTION
# --------------------------------------------------

extract_by_date <- function(nc_path, var, swor, pts) {
  
  cat("Loading:", var, "\n")
  
  r <- rast(nc_path, subds = var)
  r_dates <- as.Date(time(r))
  
  pts2 <- project(pts, crs(r))
  
  out <- rep(NA_real_, nrow(swor))
  udates <- sort(unique(swor$date))
  
  cat("Unique dates:", length(udates), "\n")
  
  for (i in seq_along(udates)) {
    
    d <- udates[i]
    
    if (i %% 50 == 0) {
      cat("  Date", i, "of", length(udates), "\n")
    }
    
    lyr <- match(d, r_dates)
    if (is.na(lyr)) next
    
    rows <- which(swor$date == d)
    
    vals <- terra::extract(r[[lyr]], pts2[rows])
    
    out[rows] <- vals[, 2]
  }
  
  out
}

# --------------------------------------------------
# RUN DAILY VARIABLES
# --------------------------------------------------

daily_vars <- c("bbv","ild","sos","ssh","ssu_rotate","ssv_rotate","tos")

for (v in daily_vars) {
  
  cat("\n====================\n")
  cat("Processing:", v, "\n")
  cat("====================\n")
  
  final_df[[v]] <- extract_by_date(
    nc_files[[v]],
    v,
    swor,
    pts
  )
}

# --------------------------------------------------
# STATIC DEPTH
# --------------------------------------------------

deptho_r <- rast(nc_files$deptho)
pts_depth <- project(pts, crs(deptho_r))

final_df$deptho <- terra::extract(deptho_r, pts_depth)[, 2]

# --------------------------------------------------
# FOCAL SD FUNCTION (SPATIAL VARIABILITY)
# --------------------------------------------------

extract_focal_sd_by_date <- function(nc_path, var, swor, pts, window = 5) {
  
  cat("Loading focal SD:", var, "\n")
  
  r <- rast(nc_path, subds = var)
  r_dates <- as.Date(time(r))
  
  pts2 <- project(pts, crs(r))
  
  out <- rep(NA_real_, nrow(swor))
  w <- matrix(1, window, window)
  
  udates <- sort(unique(swor$date))
  
  cat("Unique dates:", length(udates), "\n")
  
  for (i in seq_along(udates)) {
    
    d <- udates[i]
    
    if (i %% 50 == 0) {
      cat("  Date", i, "of", length(udates), "\n")
    }
    
    lyr <- match(d, r_dates)
    if (is.na(lyr)) next
    
    rows <- which(swor$date == d)
    
    r_sd <- focal(
      r[[lyr]],
      w = w,
      fun = sd,
      na.rm = TRUE,
      pad = TRUE
    )
    
    vals <- terra::extract(r_sd, pts2[rows])
    out[rows] <- vals[, 2]
    
    rm(r_sd)
    gc()
  }
  
  out
}

# --------------------------------------------------
# SD VARIABLES
# --------------------------------------------------

final_df$tos_sd <- extract_focal_sd_by_date(
  nc_files$tos, "tos", swor, pts, window = 5
)

final_df$ssh_sd <- extract_focal_sd_by_date(
  nc_files$ssh, "ssh", swor, pts, window = 5
)

# --------------------------------------------------
# DEPTH SD
# --------------------------------------------------

deptho_r <- rast(nc_files$deptho)

deptho_sd <- focal(
  deptho_r,
  w = matrix(1, 5, 5),
  fun = sd,
  na.rm = TRUE,
  pad = TRUE
)

final_df$deptho_sd <- terra::extract(deptho_sd, pts_depth)[, 2]

# --------------------------------------------------
# MOON PHASE
# --------------------------------------------------

final_df$moon_phase <- lunar.phase(final_df$date)

# --------------------------------------------------
# CLEAN DATASET
# --------------------------------------------------

final_df <- final_df %>%
  filter(!is.na(tos) | !is.na(ssh))

# --------------------------------------------------
# SAVE OUTPUT
# --------------------------------------------------

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_file <- file.path(output_dir, "swor_env_MOM6_gridded_daily.csv")

write.csv(final_df, output_file, row.names = FALSE)

cat("\nDONE ->", output_file, "\n")

# ==================================================
# DERIVED METRIC: EDDY KINETIC ENERGY (EKE)
# ==================================================

swor_out <- read.csv(output_file)

swor_out$eke <- 0.5 * (
  swor_out$ssu_rotate^2 +
    swor_out$ssv_rotate^2
)

swor_out <- swor_out %>%
  filter(!is.na(eke))

write.csv(swor_out, output_file, row.names = FALSE)

cat("DONE: EKE added\n")