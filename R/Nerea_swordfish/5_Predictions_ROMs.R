# ==================================================
# SWOR SPATIOTEMPORAL PREDICTION PIPELINE (ROMS/MOM6)
# ==================================================
#
# DESCRIPTION
# This script generates daily habitat suitability predictions
# for SWOR using a trained Boosted Regression Tree (BRT) model
# applied to ROMS/MOM6 environmental raster stacks.
#
# The workflow includes:
#
# 1. Clean workspace and remove auxiliary GIS artifacts
#    - Deletes unwanted .aux.xml files from prediction directories
#
# 2. Load required modelling and spatial libraries
#
# 3. Load trained BRT model ensemble
#
# 4. Build daily environmental raster stacks
#    - Matches predictor variables to model training order
#
# 5. Run spatial prediction for each day
#    - Applies BRT model to gridded environmental data
#
# 6. Spatial post-processing
#    - 5x5 focal smoothing for coastal gap filling
#
# 7. Export outputs
#    - GeoTIFF rasters
#    - PNG visualization products (optional function call)
#
# OUTPUT STRUCTURE
# For each day:
#   swor_YYYY-MM-DD_mean.tif
#   + optional PNG visualization
#
# ==================================================

# --------------------------------------------------
# CLEAN AUX FILES
# --------------------------------------------------

files <- list.files(
  path = "/Volumes/Triple_Bottom_Line/Nerea_working/Forecast_swor/FINAL/Gridded/4_Pred_ROMsb",
  pattern = "aux\\.xml$",
  recursive = TRUE,
  full.names = TRUE
)

file.remove(files)

# --------------------------------------------------
# LOAD LIBRARIES
# --------------------------------------------------

library(raster)
library(dismo)
library(rgdal)
library(maptools)
library(maps)
library(mapdata)
library(ncdf4)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(terra)

# --------------------------------------------------
# OPTIONAL NORMALIZATION FUNCTION
# --------------------------------------------------

range01 <- function(r){
  r.min <- cellStats(r, "min")
  r.max <- cellStats(r, "max")
  (r - r.min) / (r.max - r.min)
}

# --------------------------------------------------
# PATHS
# --------------------------------------------------

setwd("/Volumes/Triple_Bottom_Line/Nerea_working/Forecast_swor/FINAL/Gridded/4_Pred_ROMsb/")

template <- raster("2005-01-01/ild.grd")

outDir <- "/Volumes/Triple_Bottom_Line/Nerea_working/Forecast_swor/FINAL/Gridded/5_Predictions_ROMs_0.05/"

# --------------------------------------------------
# LOAD MODEL
# --------------------------------------------------

Species <- "SWOR"

modrep <- readRDS(
  "/Volumes/Triple_Bottom_Line/Nerea_working/Forecast_swor/FINAL/Gridded/3_Output/0.05/SWOR.res1.tc3.lr03.single_gridded_0.05.rds"
)

spname <- "swor"

# --------------------------------------------------
# PREDICTION FUNCTION
# --------------------------------------------------

predCIs_ROMS <- function(get_date,
                         spname,
                         modrep,
                         stackfile,
                         template,
                         outDir,
                         studyarea,
                         droppath){
  
  spDir <- file.path(outDir, spname)
  
  if(!dir.exists(spDir)){
    dir.create(spDir, recursive = TRUE)
  }
  
  # ------------------------------------------------
  # ALIGN VARIABLE NAMES (MUST MATCH MODEL TRAINING)
  # ------------------------------------------------
  
  names(stackfile)[1]  <- "bbv"
  names(stackfile)[2]  <- "deptho_sd"
  names(stackfile)[3]  <- "deptho"
  names(stackfile)[4]  <- "eke"
  names(stackfile)[5]  <- "ild"
  names(stackfile)[6]  <- "moon_phase"
  names(stackfile)[7]  <- "sos"
  names(stackfile)[8]  <- "ssh_sd"
  names(stackfile)[9]  <- "ssh"
  names(stackfile)[10] <- "ssu_rotate"
  names(stackfile)[11] <- "ssv_rotate"
  names(stackfile)[12] <- "tos_sd"
  names(stackfile)[13] <- "tos"
  
  # ------------------------------------------------
  # MODEL PREDICTION
  # ------------------------------------------------
  
  pred.raster <- predict(stackfile, modrep, type = "response")
  
  plot(pred.raster)
  maps::map("world2", add = TRUE, col = grey(0.7), fill = TRUE)
  
  # ------------------------------------------------
  # SPATIAL SMOOTHING (5x5 FOCAL FILTER)
  # ------------------------------------------------
  
  w <- matrix(1, 5, 5)
  
  pred.raster <- focal(
    pred.raster,
    w = w,
    fun = function(x){
      if(all(is.na(x))){
        NA
      } else {
        mean(x, na.rm = TRUE)
      }
    },
    NAonly = TRUE,
    pad = TRUE
  )
  
  # ------------------------------------------------
  # OUTPUT NAME
  # ------------------------------------------------
  
  out_file <- file.path(
    spDir,
    paste0(spname, "_", get_date, "_mean")
  )
  
  # ------------------------------------------------
  # SAVE RASTER
  # ------------------------------------------------
  
  writeRaster(
    pred.raster,
    filename = out_file,
    overwrite = TRUE,
    format = "GTiff"
  )
  
  # ------------------------------------------------
  # OPTIONAL PNG EXPORT
  # ------------------------------------------------
  
  make_png_operationalization(
    r = pred.raster,
    spname = spname,
    get_date = get_date,
    outDir = outDir,
    type = "mean"
  )
}

# --------------------------------------------------
# RUN LOOP OVER DATES
# --------------------------------------------------

get_date <- list.dirs(getwd(),
                      recursive = FALSE,
                      full.names = FALSE)

for(i in seq_along(get_date)){
  
  flist <- list.files(
    get_date[i],
    recursive = TRUE,
    full.names = TRUE,
    pattern = paste(
      "tos.grd|tos_sd.grd|ssh.grd|ssh_sd.grd|sos.grd|",
      "deptho.grd|deptho_sd.grd|ssu_rotate.grd|ssv_rotate.grd|",
      "ild.grd|bbv.grd|EKE.grd|lunar_illumination.grd$",
      sep = ""
    )
  )
  
  stackfile <- stack(flist)
  
  predCIs_ROMS(
    get_date = get_date[i],
    modrep = modrep,
    spname = spname,
    stackfile = stackfile,
    template = template,
    outDir = outDir,
    studyarea = studyarea,
    droppath = droppath
  )
  
  print(get_date[i])
}