# ============================================================================
# 01_ingest_vnrp.R
# UCFR Filamentous Algae Project
# Stage 1: Ingest VNRP biological and water chemistry data from WQP
#
# Source:   Water Quality Portal (WQP)
#           Organization: TSWQC_WQX (Tri-State Water Quality Council)
# Output:   0_data/vnrp_raw.csv
#
# Notes:
#   - Pulls full period of record on every run
#   - Returns long format with 4 columns of interest
#   - Unzips programmatically, no manual steps required
#   - Unit standardization and variable harmonization handled in later stages
# ============================================================================

library(httr)
library(readr)

# ----------------------------------------------------------------------------
# 1. Define output path
# ----------------------------------------------------------------------------

out_dir  <- "0_data"
out_file <- file.path(out_dir, "vnrp_raw.csv")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ----------------------------------------------------------------------------
# 2. Define WQP query URL
# ----------------------------------------------------------------------------

# WQP supports GET with query parameters — simpler and more reliable than POST
# providers must be repeated as separate parameters, not comma-separated

wqp_url <- paste0(
  "https://www.waterqualitydata.us/data/Result/search",
  "?mimeType=csv",
  "&zip=yes",
  "&organization=TSWQC_WQX",
  "&dataProfile=biological",
  "&providers=NWIS",
  "&providers=STORET"
)

cat("Request URL:\n", wqp_url, "\n\n")

# ----------------------------------------------------------------------------
# 3. Pull data from WQP
# ----------------------------------------------------------------------------

cat("Querying Water Quality Portal...\n")

tmp_zip <- tempfile(fileext = ".zip")

response <- GET(
  url = wqp_url,
  write_disk(tmp_zip, overwrite = TRUE),
  progress()
)

if (http_error(response)) {
  stop("WQP query failed with status: ", status_code(response),
       "\nResponse: ", content(response, as = "text", encoding = "UTF-8"))
}

cat("Download complete.\n")

# ----------------------------------------------------------------------------
# 4. Unzip and read
# ----------------------------------------------------------------------------

cat("Unzipping response...\n")

tmp_dir   <- tempdir()
zip_files <- unzip(tmp_zip, exdir = tmp_dir)

# WQP zips typically contain one CSV — find it
csv_file <- zip_files[grepl("\\.csv$", zip_files, ignore.case = TRUE)]

if (length(csv_file) == 0) {
  stop("No CSV found in WQP zip file. Files found: ",
       paste(basename(zip_files), collapse = ", "))
}

if (length(csv_file) > 1) {
  warning("Multiple CSVs found in zip — using first: ", basename(csv_file[1]))
  csv_file <- csv_file[1]
}

cat("Reading CSV:", basename(csv_file), "\n")

vnrp_raw <- read_csv(
  csv_file,
  col_types = cols(.default = col_character()),  # read all as character first
  show_col_types = FALSE
)

cat("Raw dimensions:", nrow(vnrp_raw), "rows x", ncol(vnrp_raw), "columns\n")

# ----------------------------------------------------------------------------
# 5. Retain only columns of interest
# ----------------------------------------------------------------------------

cols_keep <- c(
  "ActivityStartDate",
  "MonitoringLocationIdentifier",
  "CharacteristicName",
  "ResultMeasureValue",
  "ResultMeasure/MeasureUnitCode"
)

missing_cols <- setdiff(cols_keep, names(vnrp_raw))
if (length(missing_cols) > 0) {
  stop("Expected columns not found in WQP data: ",
       paste(missing_cols, collapse = ", "))
}

vnrp_slim <- vnrp_raw[ , cols_keep]

# Rename unit column to something friendlier
names(vnrp_slim)[names(vnrp_slim) == "ResultMeasure/MeasureUnitCode"] <- "UnitCode"

cat("Columns retained:", paste(names(vnrp_slim), collapse = ", "), "\n")

# ----------------------------------------------------------------------------
# 6. Basic ingestion diagnostics
# ----------------------------------------------------------------------------

cat("\n--- Ingestion Summary ---\n")
cat("Total records:       ", nrow(vnrp_slim), "\n")
cat("Unique sites:        ", length(unique(vnrp_slim$MonitoringLocationIdentifier)), "\n")
cat("Unique variables:    ", length(unique(vnrp_slim$CharacteristicName)), "\n")
cat("Date range:          ",
    min(vnrp_slim$ActivityStartDate, na.rm = TRUE), "to",
    max(vnrp_slim$ActivityStartDate, na.rm = TRUE), "\n")
cat("Missing values in ResultMeasureValue: ",
    sum(is.na(vnrp_slim$ResultMeasureValue) |
          vnrp_slim$ResultMeasureValue == ""), "\n")

# ----------------------------------------------------------------------------
# 7. Write to 0_data
# ----------------------------------------------------------------------------

write_csv(vnrp_slim, out_file)
cat("\nSaved to:", out_file, "\n")

# ----------------------------------------------------------------------------
# 8. Cleanup temp files
# ----------------------------------------------------------------------------

unlink(tmp_zip)
unlink(zip_files)

cat("Done.\n")