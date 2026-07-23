# =============================================================================
# 12_time_slices.R
# UCFR Cladophora Bloom Prediction Pipeline
#
# PURPOSE:
#   Aggregate the full projection trajectory from Script 11 into time slice
#   summaries for the baseline period, 2050 horizon, and 2080 horizon.
#   These summaries are the primary input for the paper figures in Script 13.
#
# INPUTS:
#   2_incremental/projections_monthly.csv  — from Script 11
#   2_incremental/ucfr_model_ready.csv     — observed data (for baseline
#                                            observed bloom values)
#
# OUTPUTS:
#   2_incremental/time_slice_summaries.csv — aggregated time slice table
#                                            one row per site × scenario ×
#                                            date point × horizon
#   2_incremental/baseline_observed.csv    — observed bloom means per site ×
#                                            month for validation overlay
#
# TIME SLICES:
#   baseline: 1998–2022 (from projection grid, delta = 0 period)
#   2050:     2040–2060 mean
#   2080:     2070–2090 mean
#
# COLUMNS IN TIME SLICE SUMMARIES:
#   Identifiers: site, scenario, date_point, month, horizon
#   Predicted bloom: pred_CHLa_med, pred_CHLa_q25, pred_CHLa_q75 (µg/L)
#   Predicted log bloom: pred_logCHLa_med, pred_logCHLa_q25, pred_logCHLa_q75
#   Mean projected predictors: Temp_oC, Q_obs_cfs_med, anomaly_med, DSF_med
#   Extrapolation: pct_flagged (% of years in slice with extrapolation flag)
#
# AUTHOR: [Rafa]
# =============================================================================


# -----------------------------------------------------------------------------
# 0. USER PARAMETERS
# -----------------------------------------------------------------------------

baseline_start <- 1998
baseline_end   <- 2022

horizons <- list(
  "baseline" = c(1998, 2022),
  "2050"     = c(2040, 2060),
  "2080"     = c(2070, 2090)
)

scenarios  <- c("ssp245", "ssp585")

site_order <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")

date_point_order <- c("early_july", "august", "mid_september")


# -----------------------------------------------------------------------------
# 1. LOAD INPUTS
# -----------------------------------------------------------------------------

cat("Loading inputs...\n")

proj     <- read.csv("2_incremental/projections_monthly.csv", stringsAsFactors = FALSE)
obs_data <- read.csv("2_incremental/ucfr_model_ready.csv",   stringsAsFactors = FALSE)

cat("  Projection rows:", nrow(proj), "\n")
cat("  Observed rows:", nrow(obs_data), "\n")


# -----------------------------------------------------------------------------
# 2. AGGREGATE PROJECTIONS TO TIME SLICES
# -----------------------------------------------------------------------------

cat("Aggregating to time slices...\n")

# Columns to average across years within each slice
mean_cols <- c(
  "Temp_oC", "Q_obs_cfs_med", "Q_obs_cfs_q25", "Q_obs_cfs_q75",
  "anomaly_med", "anomaly_q25", "anomaly_q75",
  "DSF_med", "DSF_q25", "DSF_q75",
  "pred_SPC_med", "pred_SPC_q25", "pred_SPC_q75",
  "pred_logTP_med", "pred_logTP_q25", "pred_logTP_q75",
  "pred_logTN_med", "pred_logTN_q25", "pred_logTN_q75",
  "pred_logCHLa_med", "pred_logCHLa_q25", "pred_logCHLa_q75",
  "pred_logCHLa_lo", "pred_logCHLa_hi",
  "pred_CHLa_med", "pred_CHLa_q25", "pred_CHLa_q75",
  "pred_CHLa_lo", "pred_CHLa_hi"
)

slice_list <- list()

for (scen in scenarios) {
  scen_proj <- proj[proj$scenario == scen, ]
  
  for (hz_name in names(horizons)) {
    hz_range <- horizons[[hz_name]]
    
    hz_proj <- scen_proj[
      scen_proj$year >= hz_range[1] & scen_proj$year <= hz_range[2],
    ]
    
    # Mean across years per site × date_point
    hz_means <- aggregate(
      hz_proj[, mean_cols],
      by  = list(
        site       = hz_proj$site,
        date_point = hz_proj$date_point,
        month      = hz_proj$month,
        doy        = hz_proj$doy
      ),
      FUN = mean,
      na.rm = TRUE
    )
    
    # Extrapolation: % of years flagged per site × date_point
    extrap_pct <- aggregate(
      extrapolation_flag ~ site + date_point,
      data = hz_proj,
      FUN  = function(x) round(100 * mean(x), 1)
    )
    names(extrap_pct)[3] <- "pct_flagged"
    
    hz_means <- merge(hz_means, extrap_pct, by = c("site", "date_point"), all.x = TRUE)
    
    hz_means$scenario <- scen
    hz_means$horizon  <- hz_name
    
    slice_list[[length(slice_list) + 1]] <- hz_means
  }
}

time_slices <- do.call(rbind, slice_list)
rownames(time_slices) <- NULL

# Enforce ordering
time_slices$site       <- factor(time_slices$site,       levels = site_order)
time_slices$date_point <- factor(time_slices$date_point, levels = date_point_order)
time_slices$horizon    <- factor(time_slices$horizon,
                                 levels = c("baseline", "2050", "2080"))

time_slices <- time_slices[order(
  time_slices$scenario,
  time_slices$horizon,
  time_slices$date_point,
  time_slices$site
), ]

time_slices$site       <- as.character(time_slices$site)
time_slices$date_point <- as.character(time_slices$date_point)
time_slices$horizon    <- as.character(time_slices$horizon)

cat("  Time slice rows:", nrow(time_slices), "\n")
cat("  Expected rows: 2 scenarios × 3 horizons × 3 date points × 7 sites =",
    2 * 3 * 3 * 7, "\n")


# -----------------------------------------------------------------------------
# 3. BASELINE OBSERVED BLOOM — FOR VALIDATION OVERLAY
# -----------------------------------------------------------------------------
# Compute mean observed log10(CHLa) and CHLa per site × month from the
# actual sampling record. Months 7, 8, 9 retained (matching date points).
# This is used in Script 13 to overlay observed baseline onto projections.

cat("Computing observed baseline bloom means...\n")

obs_baseline <- obs_data[
  obs_data$Year >= baseline_start & obs_data$Year <= baseline_end &
    obs_data$Month %in% c(7, 8, 9),
]

# Month → date_point label mapping
month_to_dp <- c("7" = "early_july", "8" = "august", "9" = "mid_september")

obs_means <- aggregate(
  cbind(logCHLa, CHLa) ~ Site + Month,
  data = obs_baseline,
  FUN  = mean,
  na.rm = TRUE
)
names(obs_means)[1] <- "site"
names(obs_means)[2] <- "month"

obs_means$date_point   <- month_to_dp[as.character(obs_means$month)]
obs_means$obs_CHLa_mean    <- obs_means$CHLa
obs_means$obs_logCHLa_mean <- obs_means$logCHLa

# Also compute SD and N for error bars
obs_sd <- aggregate(
  cbind(logCHLa, CHLa) ~ Site + Month,
  data = obs_baseline,
  FUN  = sd,
  na.rm = TRUE
)
obs_n <- aggregate(
  logCHLa ~ Site + Month,
  data = obs_baseline,
  FUN  = length
)
names(obs_sd)[1:2] <- c("site", "month")
names(obs_n)[1:2]  <- c("site", "month")
names(obs_sd)[3:4] <- c("obs_logCHLa_sd", "obs_CHLa_sd")
names(obs_n)[3]    <- "obs_n"

obs_summary <- merge(obs_means[, c("site", "month", "date_point",
                                   "obs_CHLa_mean", "obs_logCHLa_mean")],
                     obs_sd[, c("site", "month", "obs_logCHLa_sd", "obs_CHLa_sd")],
                     by = c("site", "month"))
obs_summary <- merge(obs_summary, obs_n, by = c("site", "month"))

# SE for plotting
obs_summary$obs_logCHLa_se <- obs_summary$obs_logCHLa_sd / sqrt(obs_summary$obs_n)

cat("  Observed baseline rows:", nrow(obs_summary), "\n")


# -----------------------------------------------------------------------------
# 4. SANITY CHECKS
# -----------------------------------------------------------------------------

cat("\n--- Sanity checks ---\n")

# NA check
na_ts <- colSums(is.na(time_slices))
if (any(na_ts > 0)) {
  cat("  WARNING: NAs in time_slices:\n")
  print(na_ts[na_ts > 0])
} else {
  cat("  No NAs in time_slices.\n")
}

# Baseline vs 2080 contrast: August, median CHLa by site
cat("\n  August bloom (µg/L): baseline vs 2050 vs 2080 — ssp585 median\n")
aug_contrast <- time_slices[
  time_slices$date_point == "august" & time_slices$scenario == "ssp585",
  c("site", "horizon", "pred_CHLa_med", "pred_CHLa_q25", "pred_CHLa_q75",
    "Temp_oC", "pct_flagged")
]
print(aug_contrast[order(aug_contrast$site, aug_contrast$horizon), ])

# Scenario contrast at 2080: ssp245 vs ssp585
cat("\n  2080 August: ssp245 vs ssp585 (median CHLa µg/L)\n")
scen_contrast <- time_slices[
  time_slices$date_point == "august" & time_slices$horizon == "2080",
  c("site", "scenario", "pred_CHLa_med", "pct_flagged")
]
print(scen_contrast[order(scen_contrast$site, scen_contrast$scenario), ])

# Date point contrast: bloom across growing season at 2080, ssp585
cat("\n  2080 ssp585: bloom by date point across sites (median CHLa µg/L)\n")
dp_contrast <- time_slices[
  time_slices$horizon == "2080" & time_slices$scenario == "ssp585",
  c("site", "date_point", "pred_CHLa_med")
]
# Reshape wide for readability
dp_wide <- reshape(dp_contrast,
                   idvar     = "site",
                   timevar   = "date_point",
                   direction = "wide")
names(dp_wide) <- gsub("pred_CHLa_med\\.", "", names(dp_wide))
print(dp_wide[order(match(dp_wide$site, site_order)), ])

# Observed baseline spot check
cat("\n  Observed baseline bloom means (µg/L) by site × month:\n")
print(obs_summary[order(obs_summary$month,
                        match(obs_summary$site, site_order)),
                  c("site", "month", "date_point", "obs_CHLa_mean",
                    "obs_logCHLa_mean", "obs_n")])


# -----------------------------------------------------------------------------
# 5. WRITE OUTPUTS
# -----------------------------------------------------------------------------

write.csv(time_slices,  "2_incremental/time_slice_summaries.csv", row.names = FALSE)
write.csv(obs_summary,  "2_incremental/baseline_observed.csv",    row.names = FALSE)

cat("\nOutputs written:\n")
cat("  2_incremental/time_slice_summaries.csv\n")
cat("  2_incremental/baseline_observed.csv\n")
cat("Done.\n")