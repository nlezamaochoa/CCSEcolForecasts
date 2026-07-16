# ==================================================
# ROMS / MOM6 RASTER CROPPING PIPELINE (CALIFORNIA CURRENT)
# ==================================================
#
# DESCRIPTION
# This script crops global/regional ocean model rasters
# to the California Current System (CCS) domain for use in
# SWOR habitat modelling and downstream prediction workflows.
#
# The workflow:
#
# 1. Reads all .grd raster files recursively from prediction directory
# 2. Converts longitude system if needed (0–360 → -180–180)
# 3. Crops rasters to California Current spatial extent
# 4. Preserves folder structure in output directory
# 5. Writes cropped rasters to new directory
#
# OUTPUT
# - Cropped .grd rasters in identical folder hierarchy
# - Restricted to CCS domain:
#   lon: -134 to -115.5
#   lat: 30 to 48
#
# USE CASE
# - Reducing data volume for modelling
# - Standardizing spatial domain across variables
# - Preparing inputs for habitat suitability models
#
# ==================================================

library(terra)

# --------------------------------------------------
# INPUT / OUTPUT DIRECTORIES
# --------------------------------------------------

pred_dir <- "/Volumes/Triple_Bottom_Line/Nerea_working/Forecast_swor/FINAL/Gridded/4_Pred"
out_dir  <- "/Volumes/Triple_Bottom_Line/Nerea_working/Forecast_swor/FINAL/Gridded/4_Pred_ROMs"

# --------------------------------------------------
# STUDY AREA EXTENT (CALIFORNIA CURRENT SYSTEM)
# --------------------------------------------------

ccs_ext <- ext(-134, -115.5, 30, 48)

# --------------------------------------------------
# FIND ALL RASTER FILES
# --------------------------------------------------

raster_files <- list.files(
  pred_dir,
  pattern = "\\.grd$",
  recursive = TRUE,
  full.names = TRUE
)

# --------------------------------------------------
# PROCESS EACH RASTER
# --------------------------------------------------

for (f in raster_files) {
  
  cat("Processing:", basename(f), "\n")
  
  # Read raster
  r <- rast(f)
  
  # ------------------------------------------------
  # LONGITUDE CORRECTION
  # ------------------------------------------------
  # Convert 0–360 grids to -180–180 if needed
  # ------------------------------------------------
  
  if (xmin(r) > 0) {
    r <- rotate(r)
  }
  
  # ------------------------------------------------
  # CROP TO STUDY AREA
  # ------------------------------------------------
  
  r_crop <- crop(r, ccs_ext)
  
  # ------------------------------------------------
  # PRESERVE FOLDER STRUCTURE
  # ------------------------------------------------
  
  rel_path <- sub(paste0("^", pred_dir, "/"), "", f)
  out_file <- file.path(out_dir, rel_path)
  
  dir.create(
    dirname(out_file),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # ------------------------------------------------
  # SAVE OUTPUT
  # ------------------------------------------------
  
  writeRaster(
    r_crop,
    out_file,
    overwrite = TRUE
  )
}

# --------------------------------------------------
# FINISH
# --------------------------------------------------

cat("Finished cropping all rasters to CCS domain\n")