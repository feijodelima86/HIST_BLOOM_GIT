# =============================================================================
# 09_climate_deltas.R
# UCFR Cladophora Bloom Prediction Pipeline
#
# PURPOSE:
#   Derive climate deltas for all projection predictors from NCAR daily Q
#   and NCCV monthly temperature. Outputs a delta table used by Script 10
#   to build the projection grid.
#
# INPUTS:
#   0_data/ncar_daily_q.csv                      — NCAR SUMMA/mizuRoute daily Q
#   0_data/UCF_HUC17010201_MMM_english.csv       — NCCV upper Clark Fork (monthly)
#   0_data/MCF_HUC17010204_MMM_english.csv       — NCCV middle Clark Fork (monthly)
#   2_incremental/ucfr_model_ready.csv           — observed data (baseline means)
#
# OUTPUTS:
#   2_incremental/climate_deltas.csv             — delta table (one row per
#                                                  site x scenario x horizon x month)
#
# DELTA CONVENTIONS:
#   Temp_oC          — additive (°C), from NCCV mean temperature F→C
#   Q_obs_cfs        — multiplicative (% change as proportion), from NCAR summer mean Q
#   anomaly          — additive (dimensionless), recomputed from NCAR peak/baseflow Q
#   Days_Since_Freshet — additive (days), from NCAR DOY of peak Q
#
# REACH-TO-SITE LOOKUP:
#   CLALO → DL
#   CLADR → GR, BN
#   CLABE → MS, BM
#   CLAPL → HU, FH
#
# HUC-TO-SITE LOOKUP (NCCV temperature):
#   UCF (HUC17010201) → DL, GR, BN, MS
#   MCF (HUC17010204) → BM, HU, FH
#
# FILTERS:
#   NCAR: cmip == "CMIP6"; scenarios ssp245 and ssp585 only
#   NCCV: scenarios ssp245 and ssp585 only (ssp370 dropped)
#
# AUTHOR: [Rafa]
# =============================================================================


# -----------------------------------------------------------------------------
# 0. USER PARAMETERS
# -----------------------------------------------------------------------------

baseline_start <- 1998
baseline_end   <- 2022

horizons <- list(
  "2050" = c(2040, 2060),
  "2080" = c(2070, 2090)
)

scenarios <- c("ssp245", "ssp585")

# Growing season months retained for delta computation
# (June=6 included for chemistry window but not a date point itself)
growing_months <- 6:9

# Date points (DOY): early July, August, mid-September
# Used downstream in Script 10; carried here for reference
date_points <- c(188, 213, 258)

# Site order (longitudinal, upstream to downstream)
site_order <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")


# -----------------------------------------------------------------------------
# 1. REACH-TO-SITE AND HUC-TO-SITE LOOKUPS
# -----------------------------------------------------------------------------

# NCAR reach columns → UCFR sites
reach_site <- list(
  CLALO = "DL",
  CLADR = c("GR", "BN"),
  CLABE = c("MS", "BM"),
  CLAPL = c("HU", "FH")
)

# NCCV HUC file → UCFR sites
huc_site <- list(
  UCF = c("DL", "GR", "BN", "MS"),
  MCF = c("BM", "HU", "FH")
)


# -----------------------------------------------------------------------------
# 2. LOAD AND FILTER NCAR DAILY Q
# -----------------------------------------------------------------------------

cat("Loading NCAR daily Q...\n")

ncar_raw <- read.csv("2_incremental/ncar_daily_q.csv", stringsAsFactors = FALSE)

# Confirm required columns present
required_cols <- c("esm", "scenario", "cmip", "date", "CLALO", "CLADR", "CLABE", "CLAPL")
missing_cols  <- setdiff(required_cols, names(ncar_raw))
if (length(missing_cols) > 0) {
  stop("Missing columns in ncar_daily_q.csv: ", paste(missing_cols, collapse = ", "))
}

# Filter to CMIP6 and target scenarios
ncar <- ncar_raw[ncar_raw$cmip == "CMIP6" & ncar_raw$scenario %in% scenarios, ]
cat("  CMIP6 rows retained:", nrow(ncar), "\n")
cat("  ESMs present:", paste(unique(ncar$esm), collapse = ", "), "\n")
cat("  Scenarios present:", paste(unique(ncar$scenario), collapse = ", "), "\n")

# Parse date
ncar$date_parsed <- as.Date(ncar$date, format = "%Y-%m-%d")
ncar$year        <- as.integer(format(ncar$date_parsed, "%Y"))
ncar$month       <- as.integer(format(ncar$date_parsed, "%m"))

# Drop year 1950 (spin-up; consistent with pipeline convention)
ncar <- ncar[ncar$year > 1950, ]

head(ncar_raw$date, 20)
class(ncar_raw$date)


# -----------------------------------------------------------------------------
# 3. COMPUTE NCAR ANNUAL METRICS PER ESM × SCENARIO × REACH
# -----------------------------------------------------------------------------
# For each water year (calendar year used here; freshet typically May-June):
#   - Q_peak_cfs:     annual maximum daily Q
#   - Q_baseflow_cfs: mean of August-September daily Q (low-flow proxy)
#   - DOY_peak:       day of year of annual peak Q
#   - summer_mean_q:  mean of June-September daily Q
#
# Units: NCAR values are in m³/s based on magnitude; convert to cfs (× 35.3147)

CMS_TO_CFS <- 35.3147

reach_cols <- c("CLALO", "CLADR", "CLABE", "CLAPL")

# Convert reach columns to cfs
ncar[, reach_cols] <- ncar[, reach_cols] * CMS_TO_CFS

cat("Computing NCAR annual metrics...\n")

# Split by ESM × scenario for annual metric computation
ncar_list <- split(ncar, list(ncar$esm, ncar$scenario), drop = TRUE)

ncar_annual_list <- lapply(ncar_list, function(df) {
  
  esm_val <- unique(df$esm)
  scen_val <- unique(df$scenario)
  
  # Annual metrics per year per reach
  years <- sort(unique(df$year))
  
  results <- lapply(years, function(yr) {
    yr_df <- df[df$year == yr, ]
    
    row_out <- data.frame(
      esm      = esm_val,
      scenario = scen_val,
      year     = yr,
      stringsAsFactors = FALSE
    )
    
    for (rc in reach_cols) {
      q <- yr_df[[rc]]
      doy <- as.integer(format(yr_df$date_parsed, "%j"))
      
      # Annual peak
      peak_idx <- which.max(q)
      row_out[[paste0(rc, "_Q_peak")]]     <- q[peak_idx]
      row_out[[paste0(rc, "_DOY_peak")]]   <- doy[peak_idx]
      
      # Baseflow: mean of Aug-Sep (months 8-9)
      bf_idx <- yr_df$month %in% c(8, 9)
      row_out[[paste0(rc, "_Q_baseflow")]] <- if (sum(bf_idx) > 0) mean(q[bf_idx], na.rm = TRUE) else NA
      
      # Summer mean: Jun-Sep (months 6-9)
      sm_idx <- yr_df$month %in% 6:9
      row_out[[paste0(rc, "_summer_mean")]] <- if (sum(sm_idx) > 0) mean(q[sm_idx], na.rm = TRUE) else NA
    }
    
    row_out
  })
  
  do.call(rbind, results)
})

ncar_annual <- do.call(rbind, ncar_annual_list)
rownames(ncar_annual) <- NULL

cat("  Annual metrics computed:", nrow(ncar_annual), "rows\n")


# -----------------------------------------------------------------------------
# 4. COMPUTE ANOMALY FROM NCAR ANNUAL METRICS
# -----------------------------------------------------------------------------
# anomaly = (Q_peak / Q_baseflow)^(1/3)
# Computed per ESM × scenario × year × reach

for (rc in reach_cols) {
  peak_col <- paste0(rc, "_Q_peak")
  bf_col   <- paste0(rc, "_Q_baseflow")
  anom_col <- paste0(rc, "_anomaly")
  ncar_annual[[anom_col]] <- (ncar_annual[[peak_col]] / ncar_annual[[bf_col]])^(1/3)
}


# -----------------------------------------------------------------------------
# 5. COMPUTE NCAR DELTAS PER REACH × SCENARIO × HORIZON
# -----------------------------------------------------------------------------
# For each metric, delta = horizon_mean - baseline_mean (additive)
# For Q_peak and summer_mean, also compute % delta = (horizon_mean / baseline_mean) - 1

metrics_additive <- c("DOY_peak", "anomaly")
metrics_pct      <- c("Q_peak", "Q_baseflow", "summer_mean")

cat("Computing NCAR deltas...\n")

ncar_delta_list <- list()

for (scen in scenarios) {
  scen_df <- ncar_annual[ncar_annual$scenario == scen, ]
  
  for (hz_name in names(horizons)) {
    hz_range <- horizons[[hz_name]]
    
    base_df <- scen_df[scen_df$year >= baseline_start & scen_df$year <= baseline_end, ]
    hz_df   <- scen_df[scen_df$year >= hz_range[1]    & scen_df$year <= hz_range[2], ]
    
    for (rc in reach_cols) {
      sites_for_reach <- reach_site[[rc]]
      
      for (site in sites_for_reach) {
        
        row_delta <- data.frame(
          site     = site,
          scenario = scen,
          horizon  = hz_name,
          stringsAsFactors = FALSE
        )
        
        # Additive deltas (median, q25, q75 across ESMs)
        for (met in metrics_additive) {
          col <- paste0(rc, "_", met)
          base_vals <- tapply(base_df[[col]], base_df$esm, mean, na.rm = TRUE)
          hz_vals   <- tapply(hz_df[[col]],   hz_df$esm,   mean, na.rm = TRUE)
          
          # Match ESMs present in both
          common_esm <- intersect(names(base_vals), names(hz_vals))
          deltas <- hz_vals[common_esm] - base_vals[common_esm]
          
          row_delta[[paste0("delta_", met, "_med")]] <- median(deltas, na.rm = TRUE)
          row_delta[[paste0("delta_", met, "_q25")]] <- quantile(deltas, 0.25, na.rm = TRUE)
          row_delta[[paste0("delta_", met, "_q75")]] <- quantile(deltas, 0.75, na.rm = TRUE)
        }
        
        # Multiplicative deltas (% change)
        for (met in metrics_pct) {
          col <- paste0(rc, "_", met)
          base_vals <- tapply(base_df[[col]], base_df$esm, mean, na.rm = TRUE)
          hz_vals   <- tapply(hz_df[[col]],   hz_df$esm,   mean, na.rm = TRUE)
          
          common_esm <- intersect(names(base_vals), names(hz_vals))
          pct_deltas <- (hz_vals[common_esm] / base_vals[common_esm]) - 1
          
          row_delta[[paste0("delta_", met, "_pct_med")]] <- median(pct_deltas, na.rm = TRUE)
          row_delta[[paste0("delta_", met, "_pct_q25")]] <- quantile(pct_deltas, 0.25, na.rm = TRUE)
          row_delta[[paste0("delta_", met, "_pct_q75")]] <- quantile(pct_deltas, 0.75, na.rm = TRUE)
        }
        
        ncar_delta_list[[length(ncar_delta_list) + 1]] <- row_delta
      }
    }
  }
}

ncar_deltas <- do.call(rbind, ncar_delta_list)
rownames(ncar_deltas) <- NULL
cat("  NCAR delta rows:", nrow(ncar_deltas), "\n")


# -----------------------------------------------------------------------------
# 6. LOAD AND PROCESS NCCV TEMPERATURE
# -----------------------------------------------------------------------------

cat("Loading NCCV temperature files...\n")

F_to_C <- function(f) (f - 32) * 5 / 9

load_nccv <- function(path, huc_label) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Parse date
  df$date_parsed <- as.Date(trimws(df$Date), format = "%m/%d/%Y")
  df$year        <- as.integer(format(df$date_parsed, "%Y"))
  df$month       <- as.integer(format(df$date_parsed, "%m"))
  df$huc         <- huc_label
  
  # Extract mean temperature columns for ssp245 and ssp585 only
  # Convert F to C
  for (scen in c("ssp245", "ssp585")) {
    col_name <- paste0(scen, " Mean temperature (deg_F)")
    if (!col_name %in% names(df)) {
      stop("Expected column not found in ", path, ": ", col_name)
    }
    df[[paste0("temp_C_", scen)]] <- F_to_C(as.numeric(df[[col_name]]))
  }
  
  df[, c("date_parsed", "year", "month", "huc", "temp_C_ssp245", "temp_C_ssp585")]
}

nccv_ucf <- load_nccv("0_data/UCF_HUC17010201_MMM_english.csv", "UCF")
nccv_mcf <- load_nccv("0_data/MCF_HUC17010204_MMM_english.csv", "MCF")

nccv <- rbind(nccv_ucf, nccv_mcf)

# Filter to growing season months only
nccv <- nccv[nccv$month %in% growing_months, ]


# -----------------------------------------------------------------------------
# 7. COMPUTE NCCV TEMPERATURE DELTAS PER HUC × SCENARIO × HORIZON × MONTH
# -----------------------------------------------------------------------------
# Additive delta: horizon mean temp - baseline mean temp (°C)
# Month-specific (Jun, Jul, Aug, Sep separately)

cat("Computing NCCV temperature deltas...\n")

nccv_delta_list <- list()

for (huc_label in c("UCF", "MCF")) {
  huc_df    <- nccv[nccv$huc == huc_label, ]
  sites_for_huc <- huc_site[[huc_label]]
  
  for (scen in scenarios) {
    temp_col  <- paste0("temp_C_", scen)
    base_df   <- huc_df[huc_df$year >= baseline_start & huc_df$year <= baseline_end, ]
    
    for (hz_name in names(horizons)) {
      hz_range <- horizons[[hz_name]]
      hz_df    <- huc_df[huc_df$year >= hz_range[1] & huc_df$year <= hz_range[2], ]
      
      for (mo in growing_months) {
        base_mean <- mean(base_df[[temp_col]][base_df$month == mo], na.rm = TRUE)
        hz_mean   <- mean(hz_df[[temp_col]][hz_df$month == mo],     na.rm = TRUE)
        delta_T   <- hz_mean - base_mean
        
        for (site in sites_for_huc) {
          nccv_delta_list[[length(nccv_delta_list) + 1]] <- data.frame(
            site          = site,
            scenario      = scen,
            horizon       = hz_name,
            month         = mo,
            delta_Temp_C  = delta_T,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

nccv_deltas <- do.call(rbind, nccv_delta_list)
rownames(nccv_deltas) <- NULL
cat("  NCCV delta rows:", nrow(nccv_deltas), "\n")


# -----------------------------------------------------------------------------
# 8. MERGE NCAR AND NCCV DELTAS INTO FINAL DELTA TABLE
# -----------------------------------------------------------------------------
# NCAR deltas are site × scenario × horizon (no month dimension — applied uniformly)
# NCCV deltas are site × scenario × horizon × month
#
# Final table: one row per site × scenario × horizon × month
# Columns:
#   site, scenario, horizon, month
#   delta_Temp_C                      (from NCCV, month-specific)
#   delta_summer_mean_pct_med/q25/q75 (from NCAR, applied to Q_obs_cfs)
#   delta_anomaly_med/q25/q75         (from NCAR, additive)
#   delta_DOY_peak_med/q25/q75        (from NCAR, additive → Days_Since_Freshet)

cat("Merging NCAR and NCCV deltas...\n")

# Expand NCAR deltas to monthly (same value repeated for each growing month)
ncar_deltas_monthly <- ncar_deltas[rep(seq_len(nrow(ncar_deltas)), each = length(growing_months)), ]
ncar_deltas_monthly$month <- rep(growing_months, times = nrow(ncar_deltas))
rownames(ncar_deltas_monthly) <- NULL

# Merge on site × scenario × horizon × month
climate_deltas <- merge(
  nccv_deltas,
  ncar_deltas_monthly,
  by = c("site", "scenario", "horizon", "month"),
  all = TRUE
)

# Enforce site order and sort
climate_deltas$site <- factor(climate_deltas$site, levels = site_order)
climate_deltas <- climate_deltas[order(
  climate_deltas$site,
  climate_deltas$scenario,
  climate_deltas$horizon,
  climate_deltas$month
), ]
climate_deltas$site <- as.character(climate_deltas$site)

cat("  Final delta table rows:", nrow(climate_deltas), "\n")
cat("  Expected rows: 7 sites × 2 scenarios × 2 horizons × 4 months =",
    7 * 2 * 2 * 4, "\n")


# -----------------------------------------------------------------------------
# 9. SANITY CHECKS
# -----------------------------------------------------------------------------

cat("\n--- Sanity checks ---\n")

# Check for NAs
na_counts <- colSums(is.na(climate_deltas))
if (any(na_counts > 0)) {
  cat("  WARNING: NAs detected in columns:\n")
  print(na_counts[na_counts > 0])
} else {
  cat("  No NAs in delta table.\n")
}

# Print temperature deltas for August by site and scenario (spot check)
cat("\n  Temperature deltas (Aug, °C) by site × scenario × horizon:\n")
aug_check <- climate_deltas[climate_deltas$month == 8, 
                            c("site", "scenario", "horizon", "delta_Temp_C")]
print(aug_check)

# Print anomaly deltas for DL (spot check)
cat("\n  Anomaly deltas for DL:\n")
dl_check <- climate_deltas[climate_deltas$site == "DL",
                           c("scenario", "horizon", "month",
                             "delta_anomaly_med", "delta_anomaly_q25", "delta_anomaly_q75")]
print(dl_check[!duplicated(dl_check[, c("scenario", "horizon")]), ])

# Print summer Q % deltas for MS (spot check)
cat("\n  Summer Q % deltas for MS:\n")
ms_check <- climate_deltas[climate_deltas$site == "MS",
                           c("scenario", "horizon", "month",
                             "delta_summer_mean_pct_med")]
print(ms_check[!duplicated(ms_check[, c("scenario", "horizon")]), ])


# -----------------------------------------------------------------------------
# 10. WRITE OUTPUT
# -----------------------------------------------------------------------------

out_path <- "2_incremental/climate_deltas.csv"
write.csv(climate_deltas, out_path, row.names = FALSE)
cat("\nOutput written to:", out_path, "\n")
cat("Done.\n")