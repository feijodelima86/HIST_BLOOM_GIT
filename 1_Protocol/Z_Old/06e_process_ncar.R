# ============================================================================
# 06e_process_ncar.R
# UCFR Filamentous Algae Project
# Stage 6: Download and process NCAR GDEX d010014 climate scenario streamflow
#
# Inputs:
#   0_data/ncar_folders.csv       (ESM-scenario inventory, from 06a/06d)
#
# Downloads (one-time, cached):
#   0_data/ncar_site_nc/{ESM_scen}_mizuRoute_daily_site.nc  (~73 MB × 27)
#
# Outputs:
#   0_data/site_reach_lookup.csv        (site-to-reach mapping, for provenance)
#   2_incremental/ncar_processed.csv    (annual metrics, keyed by study site)
#
# Metrics (per site × ESM × scenario × year):
#   summer_mean_q_cfs  — mean daily Q over Jun–Aug
#   Q_peak_cfs         — annual max daily Q
#   DOY_peak           — day of year (1–366) of annual max (freshet timing)
#   Q_baseflow_cfs     — mean daily Q over baseline period (single value per
#                        site × ESM, repeated across years)
#   anomaly            — (Q_peak / Q_baseflow) ^ (1/3)
#
# Metric definitions match 05_process_usgs.R exactly (except Q_obs_cfs and
# Days_Since_Freshet, which are observation-specific and inapplicable to
# projected scenarios). summer_mean_q_cfs is new — it serves as the
# projection-mode analog of Q_obs_cfs.
#
# NetCDF details (confirmed by 06d inspection of CCSM4_rcp85):
#   Variable:  streamflow [seg, time]   units: m3/s
#   Reach ID:  reachID [seg]            8-byte int (MERIT Hydro pfaf IDs)
#   Time:      days since 1950-01-01    calendar: proleptic_gregorian
#   Span:      1950-01-01 to ~2099-09-30  (~54695 daily time steps)
#   File size: ~73 MB per ESM-scenario combo (414 sites × full time)
#
# Run from project root.
# ============================================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(readr)
  library(dplyr)
  library(tidyr)
})

# ============================================================================
# 0. Configuration
# ============================================================================

BASELINE_START <- 1998   # first year of baseline period
BASELINE_END   <- 2025   # last year of baseline period (matches USGS script)

M3S_TO_CFS <- 35.31467  # conversion factor: m³/s -> cfs

TDS_BASE  <- "https://tds.gdex.ucar.edu/thredds"
FILE_BASE <- file.path(TDS_BASE, "fileServer/files/d010014/pnw_hydrology")

CACHE_DIR <- "0_data/ncar_site_nc"
OUT_DIR   <- "2_incremental"
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR,   showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 1. Site-reach mapping (locked)
# ============================================================================
# 4 unique MERIT Hydro reaches on the Clark Fork mainstem cover all 7 UCFR
# study sites. Mapping established in 06b; confirmed present in _site.nc by 06d.

SITE_REACH <- data.frame(
  site       = c("DL",      "GR",      "BN",      "MS",      "BM",
                 "HU",      "FH"),
  reach_code = c("CLALO",   "CLALO",   "CLADR",   "CLABE",   "CLABE",
                 "CLABE",   "CLAPL"),
  reach_id   = c(78017501,  78017501,  78017483,  78017471,  78017471,
                 78017471,  78015509),
  stringsAsFactors = FALSE
)

UNIQUE_REACHES <- SITE_REACH[!duplicated(SITE_REACH$reach_id),
                             c("reach_code", "reach_id")]

write_csv(SITE_REACH, "0_data/site_reach_lookup.csv")
cat("Site-reach lookup saved -> 0_data/site_reach_lookup.csv\n")
cat("  Study sites: ", paste(SITE_REACH$site, collapse = ", "), "\n")
cat("  Unique reaches: ", paste(UNIQUE_REACHES$reach_code, collapse = ", "), "\n\n")

# ============================================================================
# 2. Read ESM-scenario folder inventory
# ============================================================================

folders_df <- read_csv("0_data/ncar_folders.csv", show_col_types = FALSE)
cat("ESM-scenario combos: ", nrow(folders_df), "\n\n")

# ============================================================================
# 3. Helper functions
# ============================================================================

download_retry <- function(url, dest, n_try = 3, wait = 10) {
  # Skip if already cached (>1 MB = real file, not an error page)
  if (file.exists(dest) && file.size(dest) > 1e6) return(TRUE)
  for (i in seq_len(n_try)) {
    ok <- tryCatch({
      suppressWarnings(download.file(url, dest, mode = "wb", quiet = TRUE))
      TRUE
    }, error = function(e) FALSE)
    if (ok && file.exists(dest) && file.size(dest) > 1e6) return(TRUE)
    if (file.exists(dest)) file.remove(dest)
    cat(sprintf("    retry %d/%d ... ", i, n_try))
    if (i < n_try) Sys.sleep(wait * i)
  }
  FALSE
}

compute_reach_metrics <- function(q_cfs, dates, years, months,
                                  baseline_start, baseline_end) {
  # Computes annual metrics for one reach from one ESM-scenario file.
  # Returns a data.frame with one row per year.
  
  # --- Baseflow: single value over baseline period ---
  bl_mask  <- years >= baseline_start & years <= baseline_end
  baseflow <- mean(q_cfs[bl_mask], na.rm = TRUE)
  
  # --- Annual peak Q ---
  annual_peak <- tapply(q_cfs, years, function(x) {
    if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  })
  
  # --- DOY of peak (day of year, 1–366) ---
  annual_doy <- tapply(seq_along(q_cfs), years, function(idx) {
    qsub <- q_cfs[idx]; dsub <- dates[idx]
    if (all(is.na(qsub))) return(NA_integer_)
    as.integer(format(dsub[which.max(qsub)], "%j"))
  })
  
  # --- Summer mean Q (JJA, require >=60 valid days) ---
  jja       <- months %in% 6:8
  jja_q     <- q_cfs[jja]
  jja_years <- years[jja]
  annual_summer <- tapply(jja_q, jja_years, function(x) {
    if (sum(!is.na(x)) < 60) NA_real_ else mean(x, na.rm = TRUE)
  })
  
  # --- Anomaly ---
  annual_anomaly <- (as.numeric(annual_peak) / baseflow) ^ (1/3)
  
  # --- Assemble ---
  yrs <- as.integer(names(annual_peak))
  data.frame(
    year              = yrs,
    summer_mean_q_cfs = as.numeric(annual_summer[as.character(yrs)]),
    Q_peak_cfs        = as.numeric(annual_peak),
    DOY_peak          = as.integer(annual_doy[as.character(yrs)]),
    Q_baseflow_cfs    = baseflow,
    anomaly           = as.numeric(annual_anomaly),
    stringsAsFactors  = FALSE
  )
}

# ============================================================================
# Phase 1: Download all _site.nc files
# ============================================================================

cat(strrep("=", 64), "\n")
cat(" Phase 1: Downloading _mizuRoute_daily_site.nc files\n")
cat(strrep("=", 64), "\n\n")

download_status <- character(nrow(folders_df))
for (i in seq_len(nrow(folders_df))) {
  folder <- folders_df$folder[i]
  fname  <- paste0(folder, "_mizuRoute_daily_site.nc")
  url    <- file.path(FILE_BASE, folder, fname)
  dest   <- file.path(CACHE_DIR, fname)
  
  cat(sprintf("[%2d/%2d] %-45s ", i, nrow(folders_df), folder))
  if (download_retry(url, dest)) {
    cat(sprintf("OK  (%4.0f MB)\n", file.size(dest) / 1e6))
    download_status[i] <- "ok"
  } else {
    cat("FAILED\n")
    download_status[i] <- "failed"
  }
}

folders_df$download_status <- download_status
n_ok     <- sum(download_status == "ok")
n_failed <- sum(download_status == "failed")
cat(sprintf("\nDownload summary: %d OK, %d failed out of %d total\n",
            n_ok, n_failed, nrow(folders_df)))

if (n_failed > 0) {
  cat("Failed folders:\n")
  cat(paste("  ", folders_df$folder[download_status == "failed"],
            collapse = "\n"), "\n")
}
if (n_ok == 0) stop("No files downloaded — check network / GDEX status.")
cat("\n")

# ============================================================================
# Phase 2: Extract streamflow and compute annual metrics
# ============================================================================

cat(strrep("=", 64), "\n")
cat(" Phase 2: Computing annual metrics per reach\n")
cat(strrep("=", 64), "\n\n")

all_results <- vector("list", n_ok * nrow(UNIQUE_REACHES))
result_idx  <- 0
files_processed <- 0

for (i in seq_len(nrow(folders_df))) {
  if (folders_df$download_status[i] != "ok") next
  
  folder   <- folders_df$folder[i]
  esm      <- folders_df$esm[i]
  scenario <- folders_df$scenario[i]
  cmip     <- folders_df$cmip[i]
  fname    <- paste0(folder, "_mizuRoute_daily_site.nc")
  path     <- file.path(CACHE_DIR, fname)
  
  cat(sprintf("[%2d/%2d] %s\n", i, nrow(folders_df), folder))
  
  nc <- nc_open(path)
  
  # Time axis: proleptic_gregorian, days since 1950-01-01
  time_vals <- nc$dim$time$vals
  origin    <- as.Date("1950-01-01")
  dates     <- origin + time_vals
  years     <- as.integer(format(dates, "%Y"))
  months    <- as.integer(format(dates, "%m"))
  
  # Reach IDs (read once per file)
  rids <- ncvar_get(nc, "reachID")
  
  for (r in seq_len(nrow(UNIQUE_REACHES))) {
    rc  <- UNIQUE_REACHES$reach_code[r]
    rid <- UNIQUE_REACHES$reach_id[r]
    seg_idx <- which(rids == rid)
    
    if (length(seg_idx) == 0) {
      warning(sprintf("  Reach %s (%d) not found in %s — skipping", rc, rid, fname))
      next
    }
    
    # Extract: streamflow[seg, time] — seg is dim 1, time is dim 2
    q_m3s <- ncvar_get(nc, "streamflow",
                       start = c(seg_idx, 1),
                       count = c(1, -1))
    q_cfs <- as.numeric(q_m3s) * M3S_TO_CFS
    
    # Compute annual metrics
    metrics_df <- compute_reach_metrics(q_cfs, dates, years, months,
                                        BASELINE_START, BASELINE_END)
    
    # Tag with ESM metadata
    metrics_df$reach_code <- rc
    metrics_df$reach_id   <- rid
    metrics_df$esm        <- esm
    metrics_df$scenario   <- scenario
    metrics_df$cmip       <- cmip
    
    result_idx <- result_idx + 1
    all_results[[result_idx]] <- metrics_df
  }
  
  nc_close(nc)
  files_processed <- files_processed + 1
}

cat(sprintf("\nProcessed %d files, %d reach-file combinations\n",
            files_processed, result_idx))

# ============================================================================
# Phase 3: Assemble, map to study sites, write CSV
# ============================================================================

cat("\n", strrep("=", 64), "\n", sep = "")
cat(" Phase 3: Mapping to study sites and writing output\n")
cat(strrep("=", 64), "\n\n")

# Bind reach-level results
reach_df <- do.call(rbind, all_results[seq_len(result_idx)])
cat("Reach-level rows: ", nrow(reach_df), "\n")

# Join to study sites (expands 4 reaches to 7 sites)
site_df <- merge(SITE_REACH[, c("site", "reach_code", "reach_id")],
                 reach_df,
                 by = c("reach_code", "reach_id"),
                 all.x = TRUE)

# Reorder columns for readability
col_order <- c("site", "reach_code", "reach_id", "esm", "scenario", "cmip",
               "year", "summer_mean_q_cfs", "Q_peak_cfs", "DOY_peak",
               "Q_baseflow_cfs", "anomaly")
site_df <- site_df[, col_order]

# Sort
site_df <- site_df[order(site_df$site, site_df$esm, site_df$scenario,
                         site_df$year), ]
rownames(site_df) <- NULL

cat("Site-level rows:  ", nrow(site_df), "\n")
cat("Unique sites:     ", paste(sort(unique(site_df$site)), collapse = ", "), "\n")
cat("Unique ESMs:      ", length(unique(site_df$esm)), "\n")
cat("Year range:       ", min(site_df$year), "–", max(site_df$year), "\n\n")

# ============================================================================
# Summary statistics — sanity check
# ============================================================================

cat("--- Baseflow by reach (mean across ESMs, cfs) ---\n")
bf_summary <- site_df %>%
  filter(!duplicated(paste(reach_code, esm, scenario))) %>%
  group_by(reach_code) %>%
  summarise(
    mean_baseflow = round(mean(Q_baseflow_cfs, na.rm = TRUE), 1),
    min_baseflow  = round(min(Q_baseflow_cfs, na.rm = TRUE), 1),
    max_baseflow  = round(max(Q_baseflow_cfs, na.rm = TRUE), 1),
    .groups = "drop"
  )
print(as.data.frame(bf_summary))

cat("\n--- Anomaly range by reach (across all ESMs & years) ---\n")
anom_summary <- site_df %>%
  group_by(reach_code) %>%
  summarise(
    min_anom  = round(min(anomaly, na.rm = TRUE), 3),
    mean_anom = round(mean(anomaly, na.rm = TRUE), 3),
    max_anom  = round(max(anomaly, na.rm = TRUE), 3),
    .groups = "drop"
  )
print(as.data.frame(anom_summary))

cat("\n--- DOY_peak summary by reach (across all ESMs & years) ---\n")
doy_summary <- site_df %>%
  group_by(reach_code) %>%
  summarise(
    median_DOY = round(median(DOY_peak, na.rm = TRUE), 0),
    q25_DOY    = round(quantile(DOY_peak, 0.25, na.rm = TRUE), 0),
    q75_DOY    = round(quantile(DOY_peak, 0.75, na.rm = TRUE), 0),
    .groups = "drop"
  )
print(as.data.frame(doy_summary))

cat("\n--- Rows with NA by metric ---\n")
metric_cols <- c("summer_mean_q_cfs", "Q_peak_cfs", "DOY_peak",
                 "Q_baseflow_cfs", "anomaly")
for (m in metric_cols) {
  n_na <- sum(is.na(site_df[[m]]))
  cat(sprintf("  %-22s  %d NA  (%d total)\n", m, n_na, nrow(site_df)))
}

# ============================================================================
# Write output
# ============================================================================

out_file <- file.path(OUT_DIR, "ncar_processed.csv")
write_csv(site_df, out_file)
cat("\nSaved -> ", out_file, "\n", sep = "")
cat(sprintf("  %d rows × %d columns\n", nrow(site_df), ncol(site_df)))
cat("\nDone.\n")