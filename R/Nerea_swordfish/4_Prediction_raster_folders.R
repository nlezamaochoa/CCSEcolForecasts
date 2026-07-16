# ==================================================
# ROMS / MOM6 DAILY PREDICTION & DERIVED VARIABLES PIPELINE
# ==================================================
#
# DESCRIPTION
# This script processes ocean model NetCDF outputs and generates
# daily geospatial rasters for environmental predictors used in
# SWOR habitat and distribution modelling.
#
# It performs the following operations:
#
# 1. Convert NetCDF variables into daily GeoTIFF/GRD rasters
#    - Interpolates time dimension into daily folders
#    - Handles coordinate system adjustments (0–360 → -180–180)
#
# 2. Creates static environmental layers
#    - Depth (deptho) aligned to model grid
#    - Replicated across all daily folders
#
# 3. Computes spatial variability layers
#    - Focal standard deviation (5x5 window)
#    - Applied to deptho (deptho_sd)
#
# 4. Computes derived dynamic variables
#    - Eddy Kinetic Energy (EKE)
#    - log(EKE)
#
# 5. Computes astronomical covariate
#    - Lunar illumination per day
#
# OUTPUT STRUCTURE
# Each daily folder contains:
#   - ild.grd, ssh.grd, tos.grd, etc.
#   - deptho.grd
#   - deptho_sd.grd
#   - EKE.grd
#   - logEKE.grd
#   - lunar_illumination.grd
#
# USE CASE
# - Spatiotemporal ecological modelling
# - Fisheries habitat prediction (SWOR)
# - High-resolution environmental reconstruction
#
# ==================================================

# --------------------------------------------------
# LOAD LIBRARIES
# --------------------------------------------------

library(ncdf4)
library(raster)
library(lunar)

# --------------------------------------------------
# INPUT / OUTPUT PATHS
# --------------------------------------------------

pred_folder <- "/Volumes/Triple_Bottom_Line/Nerea_working/Forecast_swor/FINAL/Gridded/4_Pred/"
data_folder <- "/Volumes/Triple_Bottom_Line/Nerea_working/Forecast_swor/FINAL/Gridded/7_Regridded/"

# --------------------------------------------------
# CREATE DAILY FOLDERS LIST
# --------------------------------------------------

daily_folders <- list.dirs(
  pred_folder,
  recursive = FALSE,
  full.names = TRUE
)

daily_folders <- daily_folders[
  grepl("^\\d{4}-\\d{2}-\\d{2}$", basename(daily_folders))
]

cat("Found", length(daily_folders), "daily folders\n")

# ==================================================
# PART 1 — CONVERT NETCDF VARIABLES TO DAILY RASTERS
# ==================================================

nc_files <- list.files(data_folder,
                       pattern = "\\.nc$",
                       full.names = TRUE)

for (nc_file in nc_files) {
  
  cat("\nProcessing:", basename(nc_file), "\n")
  
  nc <- nc_open(nc_file)
  var_name <- names(nc$var)[1]
  
  # --------------------------------------------------
  # TIME
  # --------------------------------------------------
  
  time_vals <- ncvar_get(nc, "time")
  time_units <- ncatt_get(nc, "time", "units")$value
  
  origin <- substr(sub(".*since\\s+", "", time_units), 1, 10)
  dates <- as.Date(time_vals, origin = origin)
  
  cat("Time steps:", length(dates), "\n")
  
  # --------------------------------------------------
  # COORDINATES
  # --------------------------------------------------
  
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  # convert 0–360 → -180–180 if needed
  lon <- ifelse(lon > 180, lon - 360, lon)
  
  lon_order <- order(lon)
  lon <- lon[lon_order]
  
  # --------------------------------------------------
  # DATA
  # --------------------------------------------------
  
  var_data_all <- ncvar_get(nc, var_name)
  nc_close(nc)
  
  nlon <- length(lon)
  nlat <- length(lat)
  
  # --------------------------------------------------
  # TEMPLATE RASTER
  # --------------------------------------------------
  
  r_template <- raster(
    nrows = nlat,
    ncols = nlon,
    xmn = min(lon),
    xmx = max(lon),
    ymn = min(lat),
    ymx = max(lat),
    crs = CRS("+proj=longlat +datum=WGS84")
  )
  
  # --------------------------------------------------
  # LOOP OVER TIME
  # --------------------------------------------------
  
  for (tt in seq_along(dates)) {
    
    date_str <- format(dates[tt], "%Y-%m-%d")
    daily_folder <- file.path(pred_folder, date_str)
    
    if (!dir.exists(daily_folder)) next
    
    slice <- var_data_all[, , tt]
    slice <- slice[lon_order, ]
    
    if (lat[1] < lat[length(lat)]) {
      slice <- slice[, ncol(slice):1]
    }
    
    r <- setValues(r_template, as.vector(slice))
    
    writeRaster(
      r,
      filename = file.path(daily_folder, paste0(var_name, ".grd")),
      format = "raster",
      overwrite = TRUE
    )
    
    if (tt %% 100 == 0) {
      cat("Completed", tt, "of", length(dates), "\n")
    }
  }
  
  cat("Finished:", var_name, "\n")
}

cat("\nNetCDF → daily rasters completed\n")

# ==================================================
# PART 2 — DEPTHO STATIC LAYER
# ==================================================

ref <- raster(file.path(pred_folder, "1993-01-01", "ild.grd"))

nc_static <- list.files(
  "/Volumes/Triple_Bottom_Line/Nerea_working/Forecast_swor/FINAL/Gridded/7_static/",
  pattern = "static\\.nc$",
  full.names = TRUE
)[1]

nc <- nc_open(nc_static)
deptho <- ncvar_get(nc, "deptho")
nc_close(nc)

deptho <- drop(deptho)

if (ncell(ref) == length(deptho)) {
  values(ref) <- as.vector(deptho)
} else if (ncell(ref) == length(t(deptho))) {
  values(ref) <- as.vector(t(deptho))
} else {
  stop("deptho grid mismatch")
}

r_deptho <- ref

for (f in daily_folders) {
  writeRaster(
    r_deptho,
    filename = file.path(f, "deptho.grd"),
    format = "raster",
    overwrite = TRUE
  )
}

cat("Deptho written\n")

# ==================================================
# PART 3 — FOCAL SD (DEPTHO)
# ==================================================

w <- matrix(1, 5, 5)

for (folder in daily_folders) {
  
  r_file <- file.path(folder, "deptho.grd")
  
  if (!file.exists(r_file)) next
  
  r <- raster(r_file)
  
  r_sd <- focal(
    r,
    w = w,
    fun = sd,
    na.rm = TRUE,
    pad = TRUE
  )
  
  writeRaster(
    r_sd,
    filename = file.path(folder, "deptho_sd.grd"),
    format = "raster",
    overwrite = TRUE
  )
}

cat("Deptho SD completed\n")

# ==================================================
# PART 4 — EKE + LUNAR VARIABLES
# ==================================================

for (folder in daily_folders) {
  
  date_str <- basename(folder)
  current_date <- as.Date(date_str)
  
  u_file <- file.path(folder, "ssu_rotate.grd")
  v_file <- file.path(folder, "ssv_rotate.grd")
  
  if (!file.exists(u_file) | !file.exists(v_file)) next
  
  u <- raster(u_file)
  v <- raster(v_file)
  
  # -------------------------
  # EKE
  # -------------------------
  
  eke <- (u^2 + v^2) / 2
  logeke <- calc(eke, function(x) { x[x <= 0] <- NA; log(x) })
  
  # -------------------------
  # LUNAR
  # -------------------------
  
  illum <- lunar.illumination(current_date)
  
  illum_r <- raster(u)
  values(illum_r) <- illum
  
  # -------------------------
  # SAVE
  # -------------------------
  
  writeRaster(eke,
              file.path(folder, "EKE.grd"),
              overwrite = TRUE)
  
  writeRaster(logeke,
              file.path(folder, "logEKE.grd"),
              overwrite = TRUE)
  
  writeRaster(illum_r,
              file.path(folder, "lunar_illumination.grd"),
              overwrite = TRUE)
  
  cat("Done:", date_str, "\n")
}

cat("\nALL PROCESSING COMPLETE\n")