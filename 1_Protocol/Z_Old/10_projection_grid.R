# =============================================================================
# 10_projection_grid.R
# UCFR Cladophora Bloom Prediction Pipeline
#
# PURPOSE:
#   Build the full projection grid for climate scenario analysis.
#   For each site × year × scenario × date point, generates a complete row
#   of environmental predictors by linearly interpolating climate deltas
#   between the baseline anchor year and the 2050/2080 horizon midpoints.
#
# INPUTS:
#   2_incremental/climate_deltas.csv       — delta table from Script 09
#   2_incremental/ucfr_model_ready.csv     — observed data (baseline means)
#
# OUTPUTS:
#   2_incremental/projection_grid.csv      — full projection grid
#                                            (site × year × scenario × date point)
#
# PREDICTOR CONSTRUCTION:
#   Temp_oC            = baseline_mean + interpolated additive delta (NCCV, month-specific)
#   Q_obs_cfs          = baseline_mean × (1 + interpolated % delta) (NCAR summer mean)
#   anomaly            = baseline_mean + interpolated additive delta (NCAR)
#   Days_Since_Freshet = baseline_mean + interpolated additive delta (NCAR DOY_peak)
#
# DATE POINT → MONTH MAPPING:
#   DOY 188 (early July)    → month 7
#   DOY 213 (August)        → month 8
#   DOY 258 (mid-September) → month 9
#
# INTERPOLATION:
#   Monotone Hermite spline through three anchors:
#     baseline_anchor_year → delta = 0
#     2050 horizon midpoint (2050) → delta_2050
#     2080 horizon midpoint (2080) → delta_2080
#   Produces smooth continuous trajectory with no kinks at anchor points.
#   Years before baseline_anchor_year carry delta = 0 (baseline conditions)
#   Years after 2080 carry delta = delta_2080 (held flat)
#
# AUTHOR: [Rafa]
# =============================================================================


# -----------------------------------------------------------------------------
# 0. USER PARAMETERS
# -----------------------------------------------------------------------------

baseline_start       <- 1998
baseline_end         <- 2022
baseline_anchor_year <- 2010   # <-- adjust to check trajectory shape iteratively

projection_start     <- 2000
projection_end       <- 2090

scenarios  <- c("ssp245", "ssp585")

horizons <- list(
  "2050" = c(2040, 2060),
  "2080" = c(2070, 2090)
)

# Horizon midpoints used as interpolation anchors
horizon_midpoints <- c("2050" = 2050, "2080" = 2080)

# Date points: DOY → month mapping
date_points <- data.frame(
  doy   = c(188, 213, 258),
  month = c(7,   8,   9),
  label = c("early_july", "august", "mid_september"),
  stringsAsFactors = FALSE
)

# Site order (longitudinal, upstream to downstream)
site_order <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")


# -----------------------------------------------------------------------------
# 1. LOAD INPUTS
# -----------------------------------------------------------------------------

cat("Loading inputs...\n")

deltas   <- read.csv("2_incremental/climate_deltas.csv",   stringsAsFactors = FALSE)
obs_data <- read.csv("2_incremental/ucfr_model_ready.csv", stringsAsFactors = FALSE)

cat("  Delta rows:", nrow(deltas), "\n")
cat("  Observed rows:", nrow(obs_data), "\n")


# -----------------------------------------------------------------------------
# 2. COMPUTE BASELINE MEANS PER SITE × MONTH
# -----------------------------------------------------------------------------
# For each predictor, compute the mean over the baseline period
# grouped by site × month. Only months 7, 8, 9 are needed (date point months).

cat("Computing baseline means...\n")

baseline_obs <- obs_data[
  obs_data$Year >= baseline_start & obs_data$Year <= baseline_end &
    obs_data$Month %in% date_points$month,
]

predictors <- c("Q_obs_cfs", "Temp_oC", "anomaly", "Days_Since_Freshet")

# Compute means per site × month
baseline_means <- aggregate(
  baseline_obs[, predictors],
  by  = list(site = baseline_obs$Site, month = baseline_obs$Month),
  FUN = mean,
  na.rm = TRUE
)

# Rename for clarity
names(baseline_means)[names(baseline_means) %in% predictors] <-
  paste0("base_", predictors)

cat("  Baseline mean rows:", nrow(baseline_means), "\n")
cat("  Sites:", paste(sort(unique(baseline_means$site)), collapse = ", "), "\n")
cat("  Months:", paste(sort(unique(baseline_means$month)), collapse = ", "), "\n")

# Sanity check: confirm all site × month combinations present
expected_combos <- length(site_order) * nrow(date_points)
if (nrow(baseline_means) != expected_combos) {
  warning("Expected ", expected_combos, " baseline mean rows but got ",
          nrow(baseline_means), ". Check for missing site × month combinations.")
}


# -----------------------------------------------------------------------------
# 3. BUILD INTERPOLATION FUNCTION — MONOTONE HERMITE SPLINE
# -----------------------------------------------------------------------------
# Fits a monotone-preserving spline (Fritsch-Carlson) through three anchors:
#   (baseline_anchor_year, 0), (2050, delta_2050), (2080, delta_2080)
# Evaluated at each projection year.
# Years <= baseline_anchor_year: delta = 0 (held flat)
# Years > 2080: delta = delta_2080 (held flat at last anchor)
#
# monoH.FC prevents overshoot between anchors — critical when both deltas
# have the same sign and a linear kink at 2050 would be ecologically misleading.

interpolate_delta <- function(year, delta_2050, delta_2080,
                              anchor = baseline_anchor_year) {
  mid_2050 <- horizon_midpoints["2050"]
  mid_2080 <- horizon_midpoints["2080"]
  
  # Three anchor points
  x_anchors <- c(anchor,  mid_2050,    mid_2080)
  y_anchors <- c(0,       delta_2050,  delta_2080)
  
  # Build monotone spline through anchors
  spline_fn <- splinefun(x_anchors, y_anchors, method = "monoH.FC")
  
  # Evaluate, clamping outside anchor range
  delta_out          <- numeric(length(year))
  pre                <- year <= anchor
  interior           <- year > anchor & year <= mid_2080
  post               <- year > mid_2080
  
  delta_out[pre]      <- 0
  delta_out[interior] <- spline_fn(year[interior])
  delta_out[post]     <- delta_2080
  
  delta_out
}


# -----------------------------------------------------------------------------
# 4. BUILD PROJECTION GRID
# -----------------------------------------------------------------------------
# For each site × scenario × date point, generate one row per projection year.
# Apply interpolated deltas to baseline means.

cat("Building projection grid...\n")

years <- projection_start:projection_end
grid_list <- list()

for (site in site_order) {
  for (scen in scenarios) {
    for (dp_idx in seq_len(nrow(date_points))) {
      
      dp_doy   <- date_points$doy[dp_idx]
      dp_month <- date_points$month[dp_idx]
      dp_label <- date_points$label[dp_idx]
      
      # Pull baseline means for this site × month
      bm_row <- baseline_means[
        baseline_means$site == site & baseline_means$month == dp_month,
      ]
      
      if (nrow(bm_row) == 0) {
        warning("No baseline mean for site=", site, " month=", dp_month, " — skipping.")
        next
      }
      
      # Pull deltas for this site × scenario × month
      # Temperature: month-specific
      # Q, anomaly, DSF: month-invariant (same value repeated across months in delta table)
      delta_row_2050 <- deltas[
        deltas$site == site & deltas$scenario == scen &
          deltas$horizon == "2050" & deltas$month == dp_month,
      ]
      delta_row_2080 <- deltas[
        deltas$site == site & deltas$scenario == scen &
          deltas$horizon == "2080" & deltas$month == dp_month,
      ]
      
      if (nrow(delta_row_2050) == 0 || nrow(delta_row_2080) == 0) {
        warning("Missing delta rows for site=", site, " scenario=", scen,
                " month=", dp_month, " — skipping.")
        next
      }
      
      # Build year-by-year rows
      n_years <- length(years)
      
      # Interpolated deltas for each predictor across all years
      delta_Temp_med <- interpolate_delta(years,
                                          delta_row_2050$delta_Temp_C,
                                          delta_row_2080$delta_Temp_C)
      
      delta_Q_med    <- interpolate_delta(years,
                                          delta_row_2050$delta_summer_mean_pct_med,
                                          delta_row_2080$delta_summer_mean_pct_med)
      delta_Q_q25    <- interpolate_delta(years,
                                          delta_row_2050$delta_summer_mean_pct_q25,
                                          delta_row_2080$delta_summer_mean_pct_q25)
      delta_Q_q75    <- interpolate_delta(years,
                                          delta_row_2050$delta_summer_mean_pct_q75,
                                          delta_row_2080$delta_summer_mean_pct_q75)
      
      delta_anom_med <- interpolate_delta(years,
                                          delta_row_2050$delta_anomaly_med,
                                          delta_row_2080$delta_anomaly_med)
      delta_anom_q25 <- interpolate_delta(years,
                                          delta_row_2050$delta_anomaly_q25,
                                          delta_row_2080$delta_anomaly_q25)
      delta_anom_q75 <- interpolate_delta(years,
                                          delta_row_2050$delta_anomaly_q75,
                                          delta_row_2080$delta_anomaly_q75)
      
      delta_DSF_med  <- interpolate_delta(years,
                                          delta_row_2050$delta_DOY_peak_med,
                                          delta_row_2080$delta_DOY_peak_med)
      delta_DSF_q25  <- interpolate_delta(years,
                                          delta_row_2050$delta_DOY_peak_q25,
                                          delta_row_2080$delta_DOY_peak_q25)
      delta_DSF_q75  <- interpolate_delta(years,
                                          delta_row_2050$delta_DOY_peak_q75,
                                          delta_row_2080$delta_DOY_peak_q75)
      
      # Apply deltas to baseline means
      # Temperature: additive (single delta — NCCV has no ESM spread)
      Temp_med <- bm_row$base_Temp_oC + delta_Temp_med
      
      # Q_obs_cfs: multiplicative
      Q_med    <- bm_row$base_Q_obs_cfs * (1 + delta_Q_med)
      Q_q25    <- bm_row$base_Q_obs_cfs * (1 + delta_Q_q25)
      Q_q75    <- bm_row$base_Q_obs_cfs * (1 + delta_Q_q75)
      
      # anomaly: additive
      anom_med <- bm_row$base_anomaly + delta_anom_med
      anom_q25 <- bm_row$base_anomaly + delta_anom_q25
      anom_q75 <- bm_row$base_anomaly + delta_anom_q75
      
      # Days_Since_Freshet: additive (earlier freshet → larger DSF)
      # DSF = sampling_DOY - DOY_peak; earlier peak (negative DOY delta) → more days since
      DSF_med  <- bm_row$base_Days_Since_Freshet - delta_DSF_med
      DSF_q25  <- bm_row$base_Days_Since_Freshet - delta_DSF_q25
      DSF_q75  <- bm_row$base_Days_Since_Freshet - delta_DSF_q75
      
      grid_block <- data.frame(
        site          = site,
        year          = years,
        scenario      = scen,
        date_point    = dp_label,
        doy           = dp_doy,
        month         = dp_month,
        
        # Baseline means (constant across years — useful for validation)
        base_Temp_oC           = bm_row$base_Temp_oC,
        base_Q_obs_cfs         = bm_row$base_Q_obs_cfs,
        base_anomaly           = bm_row$base_anomaly,
        base_Days_Since_Freshet = bm_row$base_Days_Since_Freshet,
        
        # Projected predictors (median)
        Temp_oC            = Temp_med,
        Q_obs_cfs_med      = Q_med,
        Q_obs_cfs_q25      = Q_q25,
        Q_obs_cfs_q75      = Q_q75,
        anomaly_med        = anom_med,
        anomaly_q25        = anom_q25,
        anomaly_q75        = anom_q75,
        Days_Since_Freshet_med = DSF_med,
        Days_Since_Freshet_q25 = DSF_q25,
        Days_Since_Freshet_q75 = DSF_q75,
        
        stringsAsFactors = FALSE
      )
      
      grid_list[[length(grid_list) + 1]] <- grid_block
    }
  }
}

projection_grid <- do.call(rbind, grid_list)
rownames(projection_grid) <- NULL

cat("  Projection grid rows:", nrow(projection_grid), "\n")
cat("  Expected rows:", length(site_order), "sites ×",
    length(years), "years ×",
    length(scenarios), "scenarios ×",
    nrow(date_points), "date points =",
    length(site_order) * length(years) * length(scenarios) * nrow(date_points), "\n")


# -----------------------------------------------------------------------------
# 5. SANITY CHECKS
# -----------------------------------------------------------------------------

cat("\n--- Sanity checks ---\n")

# Check for NAs
na_counts <- colSums(is.na(projection_grid))
if (any(na_counts > 0)) {
  cat("  WARNING: NAs in columns:\n")
  print(na_counts[na_counts > 0])
} else {
  cat("  No NAs in projection grid.\n")
}

# Check interpolation: DL, August, ssp585 — Temp should rise smoothly
cat("\n  Temp_oC trajectory: DL / august / ssp585 (every 10 years)\n")
check <- projection_grid[
  projection_grid$site == "DL" &
    projection_grid$date_point == "august" &
    projection_grid$scenario == "ssp585" &
    projection_grid$year %% 10 == 0,
  c("year", "Temp_oC", "Q_obs_cfs_med", "anomaly_med", "Days_Since_Freshet_med")
]
print(check)

# Check that baseline period rows carry delta = 0 (projected == baseline)
cat("\n  Baseline period check: projected Temp should equal base_Temp for years <=",
    baseline_anchor_year, "\n")
pre_anchor <- projection_grid[
  projection_grid$year <= baseline_anchor_year &
    projection_grid$site == "DL" &
    projection_grid$date_point == "august" &
    projection_grid$scenario == "ssp585",
  c("year", "base_Temp_oC", "Temp_oC")
]
print(pre_anchor)

# Range checks: no negative Q or anomaly
if (any(projection_grid$Q_obs_cfs_med < 0, na.rm = TRUE)) {
  cat("  WARNING: negative Q_obs_cfs_med values detected.\n")
} else {
  cat("  Q_obs_cfs_med: all non-negative. OK\n")
}

if (any(projection_grid$anomaly_med < 0, na.rm = TRUE)) {
  cat("  WARNING: negative anomaly_med values detected.\n")
} else {
  cat("  anomaly_med: all non-negative. OK\n")
}


# -----------------------------------------------------------------------------
# 6. WRITE OUTPUT
# -----------------------------------------------------------------------------

out_path <- "2_incremental/projection_grid.csv"
write.csv(projection_grid, out_path, row.names = FALSE)
cat("\nOutput written to:", out_path, "\n")
cat("Done.\n")