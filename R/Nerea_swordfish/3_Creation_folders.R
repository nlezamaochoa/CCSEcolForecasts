 #==================================================
# DAILY FOLDER STRUCTURE CREATION FOR SWOR PREDICTIONS
# ==================================================
#
# DESCRIPTION
# This script creates a complete directory structure of daily folders
# for storing model predictions or outputs derived from environmental
# and fisheries modelling workflows.
#
# Each folder corresponds to a single day between 1993 and 2019,
# enabling structured storage of time-explicit prediction outputs
# (e.g., habitat suitability maps, BRT predictions, diagnostics).
#
# MAIN STEPS
# 1. Load required libraries
# 2. Identify input variable directories (dynamic regridded data)
# 3. Define output prediction folder structure
# 4. Generate full sequence of daily dates (1993–2019)
# 5. Create one folder per date if it does not already exist
#
# OUTPUT
# A hierarchical folder structure:
#
# Pred/
#   ├── 1993-01-01/
#   ├── 1993-01-02/
#   ├── ...
#   └── 2019-12-31/
#
# USE CASE
# - Storing daily model predictions
# - Organizing large spatiotemporal outputs
# - Supporting batch processing pipelines
#
# ==================================================

library(ncdf4)

# --------------------------------------------------
# INPUT VARIABLE DIRECTORIES
# --------------------------------------------------

variable_folders <- list.files(
  "/Volumes/Triple_Bottom_Line/Nerea_working/Forecast_swor/FINAL/Gridded/6_dynamic_regridded/",
  full.names = TRUE
)

# --------------------------------------------------
# OUTPUT DIRECTORY
# --------------------------------------------------

pred_folder <- "/Volumes/Triple_Bottom_Line/Nerea_working/Forecast_swor/FINAL/Gridded/Pred/"

if (!dir.exists(pred_folder)) {
  dir.create(pred_folder, recursive = TRUE)
}

# --------------------------------------------------
# CREATE DAILY DATE SEQUENCE
# --------------------------------------------------

dates <- format(
  seq(
    from = as.Date("1993-01-01"),
    to   = as.Date("2019-12-31"),
    by   = "day"
  ),
  "%Y-%m-%d"
)

# --------------------------------------------------
# CREATE DAILY FOLDERS
# --------------------------------------------------

for (date_str in dates) {
  
  daily_folder <- file.path(pred_folder, date_str)
  
  if (!dir.exists(daily_folder)) {
    dir.create(
      daily_folder,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
}

# --------------------------------------------------
# SUMMARY MESSAGE
# --------------------------------------------------

cat(
  length(dates),
  "daily folders created successfully in:\n",
  pred_folder,
  "\n"
)