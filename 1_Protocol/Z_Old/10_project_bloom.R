# ============================================================================
# 10_project_bloom.R
# UCFR Filamentous Algae Project
# Stage 10: Project bloom biomass — annual time series via delta method
#
# Inputs:
#   2_incremental/ucfr_model_ready.csv       (observed 1998-2023)
#   2_incremental/ncar_processed.csv         (NCAR CMIP6 annual metrics)
#   0_data/UCF_HUC17010201_MMM_english.csv   (NCCV upper Clark Fork)
#   0_data/MCF_HUC17010204_MMM_english.csv   (NCCV middle Clark Fork)
#   3_models/brt_SPC.rds
#   3_models/brt_TP.rds
#   3_models/brt_TN.rds
#   3_models/brt_bloom_fitted.rds
#
# Outputs:
#   2_incremental/projections_annual.csv
#     One row per site x ESM x scenario x year (1998-2099)
#     source: "observed" / "observed_interp" / "projected"
#
#   2_incremental/projections_timeslice.csv
#     Median + IQR across ESMs x years
#     One row per site x scenario x horizon (2050, 2080)
#
# Design (Option B):
#   1998-2023 (observed period):
#     - Predictors from ucfr_model_ready.csv, aggregated annual per site
#     - Same values for all scenarios (no divergence)
#     - Gaps interpolated within site, flagged
#
#   2024-2099 (projected period):
#     - anomaly, Q_obs_cfs (=summer_mean_q_cfs), DOY_peak from NCAR
#     - Temp_oC from NCCV JJA mean (F -> C)
#     - Days_Since_Freshet = 213 - DOY_peak (August 1)
#     - Diverges by scenario and ESM
#
# Notes:
#   - Year 1950 discarded from NCAR (corrupt)
#   - CMIP6 only: 5 ESMs
#   - SSP370 dropped
#   - Extrapolation flagged when predictors exceed training envelope
# ============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(dismo)
  library(gbm)
})

# ============================================================================
# 0. Configuration
# ============================================================================

OBS_START <- 1998
OBS_END   <- 2023
PROJ_START <- 2024
PROJ_END   <- 2099

HORIZON_2050 <- c(2040, 2060)
HORIZON_2080 <- c(2070, 2090)

CMIP6_ESMS <- c("CanESM5", "CMCC-CM2-SR5", "MIROC-ES2L",
                "MPI-M.MPI-ESM1-2-LR", "NorESM2-MM")

SCENARIO_MAP <- c(ssp245 = "RCP4.5", ssp585 = "RCP8.5")

UCF_SITES  <- c("DL", "GR", "BN", "MS")
MCF_SITES  <- c("BM", "HU", "FH")
SITE_ORDER <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")

SUMMER_MONTHS <- 6:8
SAMPLE_DOY    <- 213

OUT_DIR <- "2_incremental"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Helper functions
# ============================================================================

f_to_c <- function(f) (f - 32) * 5 / 9

run_chain <- function(newdat, brt_SPC, brt_TP, brt_TN, brt_bloom) {
  pred_SPC   <- predict(brt_SPC,   newdata = newdat,
                        n.trees = brt_SPC$gbm.call$best.trees)
  pred_logTP <- predict(brt_TP,    newdata = newdat,
                        n.trees = brt_TP$gbm.call$best.trees)
  pred_logTN <- predict(brt_TN,    newdata = newdat,
                        n.trees = brt_TN$gbm.call$best.trees)
  
  bloom_dat <- data.frame(
    pred_SPC           = pred_SPC,
    pred_logTP         = pred_logTP,
    pred_logTN         = pred_logTN,
    anomaly            = newdat$anomaly,
    Q_obs_cfs          = newdat$Q_obs_cfs,
    Temp_oC            = newdat$Temp_oC,
    Days_Since_Freshet = newdat$Days_Since_Freshet
  )
  
  pred_logCHLa <- predict(brt_bloom, newdata = bloom_dat,
                          n.trees = brt_bloom$gbm.call$best.trees)
  
  data.frame(
    pred_SPC     = pred_SPC,
    pred_logTP   = pred_logTP,
    pred_logTN   = pred_logTN,
    pred_logCHLa = pred_logCHLa,
    pred_CHLa    = 10^pred_logCHLa
  )
}

# Linear interpolation within each site for missing predictor values
interp_within_site <- function(df, value_col, year_col = "year",
                               site_col = "site") {
  out <- df
  out$interp_flag <- FALSE
  for (s in unique(out[[site_col]])) {
    idx <- which(out[[site_col]] == s)
    yrs <- out[[year_col]][idx]
    vals <- out[[value_col]][idx]
    if (sum(!is.na(vals)) >= 2) {
      filled <- approx(yrs[!is.na(vals)], vals[!is.na(vals)],
                       xout = yrs, rule = 2)$y
      na_mask <- is.na(vals) & !is.na(filled)
      out[[value_col]][idx[na_mask]] <- filled[na_mask]
      out$interp_flag[idx[na_mask]] <- TRUE
    }
  }
  out
}

# ============================================================================
# 1. Load models and training envelope
# ============================================================================

cat("Loading BRT model objects...\n")
brt_SPC   <- readRDS("3_models/brt_SPC.rds")
brt_TP    <- readRDS("3_models/brt_TP.rds")
brt_TN    <- readRDS("3_models/brt_TN.rds")
brt_bloom <- readRDS("3_models/brt_bloom_fitted.rds")

cat("Loading observed data...\n")
obs <- as.data.frame(read_csv("2_incremental/ucfr_model_ready.csv",
                              show_col_types = FALSE))

pred_vars <- c("anomaly", "Q_obs_cfs", "Temp_oC", "Days_Since_Freshet")

envelope <- list()
for (v in pred_vars) {
  envelope[[paste0(v, "_min")]] <- min(obs[[v]], na.rm = TRUE)
  envelope[[paste0(v, "_max")]] <- max(obs[[v]], na.rm = TRUE)
}

cat("Training envelope:\n")
for (v in pred_vars) {
  cat(sprintf("  %-25s  [%8.2f, %8.2f]\n",
              v, envelope[[paste0(v, "_min")]], envelope[[paste0(v, "_max")]]))
}
cat("\n")

check_extrapolation <- function(newdat) {
  flag <- rep(FALSE, nrow(newdat))
  for (v in pred_vars) {
    flag <- flag | (newdat[[v]] < envelope[[paste0(v, "_min")]]) |
      (newdat[[v]] > envelope[[paste0(v, "_max")]])
  }
  flag
}

# ============================================================================
# 2. Build observed annual table (1998-2023)
# ============================================================================

cat(strrep("=", 60), "\n")
cat(" Phase 1: Observed annual aggregation (1998-2023)\n")
cat(strrep("=", 60), "\n\n")

obs_filt <- obs[obs$Year >= OBS_START & obs$Year <= OBS_END, ]

# Aggregate to annual mean per site x year across available months
obs_annual <- obs_filt %>%
  group_by(Site, Year) %>%
  summarise(
    anomaly            = mean(anomaly,            na.rm = TRUE),
    Q_obs_cfs          = mean(Q_obs_cfs,          na.rm = TRUE),
    Temp_oC            = mean(Temp_oC,            na.rm = TRUE),
    Days_Since_Freshet = mean(Days_Since_Freshet, na.rm = TRUE),
    n_months           = n(),
    .groups = "drop"
  ) %>%
  as.data.frame()

# Replace NaN with NA
for (v in pred_vars) {
  obs_annual[[v]][is.nan(obs_annual[[v]])] <- NA_real_
}

names(obs_annual)[names(obs_annual) == "Site"] <- "site"
names(obs_annual)[names(obs_annual) == "Year"] <- "year"

# Build a complete site x year grid then merge
full_grid <- expand.grid(site = SITE_ORDER, year = OBS_START:OBS_END,
                         stringsAsFactors = FALSE)
obs_annual <- merge(full_grid, obs_annual, by = c("site", "year"), all.x = TRUE)
obs_annual$n_months[is.na(obs_annual$n_months)] <- 0

cat(sprintf("Observed annual grid: %d rows\n", nrow(obs_annual)))
cat("NA counts before interpolation:\n")
for (v in pred_vars) {
  cat(sprintf("  %-25s  %d\n", v, sum(is.na(obs_annual[[v]]))))
}

# Interpolate gaps within site
obs_annual$any_interp <- FALSE
for (v in pred_vars) {
  obs_annual <- interp_within_site(obs_annual, v)
  obs_annual$any_interp <- obs_annual$any_interp | obs_annual$interp_flag
  obs_annual$interp_flag <- NULL
}

cat("\nNA counts after interpolation:\n")
for (v in pred_vars) {
  cat(sprintf("  %-25s  %d\n", v, sum(is.na(obs_annual[[v]]))))
}
cat(sprintf("Rows with any interpolation: %d\n\n", sum(obs_annual$any_interp)))

obs_annual$source <- ifelse(obs_annual$any_interp,
                            "observed_interp", "observed")

# ============================================================================
# 3. Build NCCV annual JJA temperature
# ============================================================================

cat(strrep("=", 60), "\n")
cat(" Phase 2: NCCV annual JJA temperature\n")
cat(strrep("=", 60), "\n\n")

read_nccv_temp <- function(path, region_label) {
  cat(sprintf("Reading %s (%s)...\n", basename(path), region_label))
  
  raw <- read_csv(path, show_col_types = FALSE)
  names(raw) <- trimws(names(raw))
  
  raw$Date  <- as.Date(trimws(raw[[1]]), format = "%m/%d/%Y")
  raw$Year  <- as.integer(format(raw$Date, "%Y"))
  raw$Month <- as.integer(format(raw$Date, "%m"))
  
  jja <- raw[raw$Month %in% SUMMER_MONTHS, ]
  
  results <- data.frame()
  for (sc in c("ssp245", "ssp585")) {
    col_f <- paste0(sc, " Mean temperature (deg_F)")
    if (!col_f %in% names(jja)) next
    
    jja_c <- f_to_c(jja[[col_f]])
    annual <- tapply(jja_c, jja$Year, mean, na.rm = TRUE)
    
    df <- data.frame(
      year     = as.integer(names(annual)),
      Temp_oC  = as.numeric(annual),
      scenario = SCENARIO_MAP[sc],
      region   = region_label,
      stringsAsFactors = FALSE
    )
    results <- rbind(results, df)
  }
  results
}

ucf_temp <- read_nccv_temp("0_data/UCF_HUC17010201_MMM_english.csv", "UCF")
mcf_temp <- read_nccv_temp("0_data/MCF_HUC17010204_MMM_english.csv", "MCF")

nccv_temp <- rbind(ucf_temp, mcf_temp)

# Expand to sites
nccv_temp_sites <- data.frame()
for (s in SITE_ORDER) {
  region <- if (s %in% UCF_SITES) "UCF" else "MCF"
  d <- nccv_temp[nccv_temp$region == region, ]
  d$site <- s
  nccv_temp_sites <- rbind(nccv_temp_sites, d)
}

cat(sprintf("\nNCCV annual JJA temp: %d rows\n", nrow(nccv_temp_sites)))
cat(sprintf("Year range: %d to %d\n\n",
            min(nccv_temp_sites$year), max(nccv_temp_sites$year)))

# ============================================================================
# 4. Build NCAR annual table (2024-2099)
# ============================================================================

cat(strrep("=", 60), "\n")
cat(" Phase 3: NCAR annual hydrology (2024-2099)\n")
cat(strrep("=", 60), "\n\n")

ncar <- as.data.frame(read_csv("2_incremental/ncar_processed.csv",
                               show_col_types = FALSE))

ncar_filt <- ncar[
  ncar$esm      %in% CMIP6_ESMS &
    ncar$scenario %in% names(SCENARIO_MAP) &
    ncar$year     >= PROJ_START &
    ncar$year     <= PROJ_END,
]

ncar_filt$scenario <- SCENARIO_MAP[ncar_filt$scenario]

cat(sprintf("NCAR projected rows: %d\n", nrow(ncar_filt)))
cat(sprintf("Year range: %d to %d\n",
            min(ncar_filt$year), max(ncar_filt$year)))
cat(sprintf("Sites: %s\n",
            paste(sort(unique(ncar_filt$site)), collapse = ", ")))
cat(sprintf("ESMs: %s\n\n",
            paste(sort(unique(ncar_filt$esm)), collapse = ", ")))

# Rename for chain compatibility
names(ncar_filt)[names(ncar_filt) == "summer_mean_q_cfs"] <- "Q_obs_cfs"

# Days_Since_Freshet
ncar_filt$Days_Since_Freshet <- SAMPLE_DOY - ncar_filt$DOY_peak

# Merge NCCV temperature onto NCAR
ncar_proj <- merge(
  ncar_filt[ , c("site", "esm", "scenario", "year",
                 "anomaly", "Q_obs_cfs", "Days_Since_Freshet",
                 "DOY_peak", "Q_peak_cfs", "Q_baseflow_cfs")],
  nccv_temp_sites[ , c("site", "scenario", "year", "Temp_oC")],
  by = c("site", "scenario", "year"),
  all.x = TRUE
)

ncar_proj$source <- "projected"

cat(sprintf("Projected table: %d rows\n", nrow(ncar_proj)))
cat(sprintf("Missing Temp_oC: %d rows\n", sum(is.na(ncar_proj$Temp_oC))))
n_neg <- sum(ncar_proj$Days_Since_Freshet < 0, na.rm = TRUE)
if (n_neg > 0) {
  cat(sprintf("NOTE: %d projected rows have negative DSF (freshet after Aug 1)\n",
              n_neg))
}
cat("\n")

# ============================================================================
# 5. Expand observed across scenarios and ESMs
# ============================================================================

cat(strrep("=", 60), "\n")
cat(" Phase 4: Expanding observed across scenarios x ESMs\n")
cat(strrep("=", 60), "\n\n")

obs_expanded <- data.frame()
for (sc in c("RCP4.5", "RCP8.5")) {
  for (esm in CMIP6_ESMS) {
    d <- obs_annual
    d$scenario <- sc
    d$esm      <- esm
    obs_expanded <- rbind(obs_expanded, d)
  }
}

obs_expanded$DOY_peak       <- SAMPLE_DOY - obs_expanded$Days_Since_Freshet
obs_expanded$Q_peak_cfs     <- NA_real_
obs_expanded$Q_baseflow_cfs <- NA_real_

cat(sprintf("Observed expanded: %d rows\n", nrow(obs_expanded)))

common_cols <- c("site", "esm", "scenario", "year",
                 "anomaly", "Q_obs_cfs", "Temp_oC",
                 "Days_Since_Freshet", "DOY_peak",
                 "Q_peak_cfs", "Q_baseflow_cfs", "source")

obs_final  <- obs_expanded[ , common_cols]
proj_final <- ncar_proj[ , common_cols]

pred_table <- rbind(obs_final, proj_final)
pred_table <- pred_table[order(pred_table$scenario, pred_table$esm,
                               pred_table$site, pred_table$year), ]
rownames(pred_table) <- NULL

cat(sprintf("Combined predictor table: %d rows\n", nrow(pred_table)))
cat(sprintf("  Observed: %d  Interp: %d  Projected: %d\n\n",
            sum(pred_table$source == "observed"),
            sum(pred_table$source == "observed_interp"),
            sum(pred_table$source == "projected")))

pred_table$extrap_flag <- check_extrapolation(pred_table)
cat(sprintf("Extrapolation flagged: %d rows (%.1f%%)\n\n",
            sum(pred_table$extrap_flag, na.rm = TRUE),
            100 * mean(pred_table$extrap_flag, na.rm = TRUE)))

# ============================================================================
# 6. Run prediction chain
# ============================================================================

cat(strrep("=", 60), "\n")
cat(" Phase 5: Running prediction chain\n")
cat(strrep("=", 60), "\n\n")

complete_mask <- complete.cases(pred_table[ , pred_vars])
pred_complete <- pred_table[complete_mask, ]

cat(sprintf("Complete cases: %d of %d (%.1f%%)\n\n",
            nrow(pred_complete), nrow(pred_table),
            100 * nrow(pred_complete) / nrow(pred_table)))

chain_out <- run_chain(
  newdat    = pred_complete[ , pred_vars],
  brt_SPC   = brt_SPC,
  brt_TP    = brt_TP,
  brt_TN    = brt_TN,
  brt_bloom = brt_bloom
)

annual_out <- cbind(pred_complete, chain_out)
annual_out <- annual_out[order(annual_out$scenario, annual_out$esm,
                               annual_out$site, annual_out$year), ]
rownames(annual_out) <- NULL

cat(sprintf("Annual projection rows: %d\n", nrow(annual_out)))

# ============================================================================
# 7. Hindcast check
# ============================================================================

cat("\n--- Hindcast check (1998-2023 median, RCP4.5 + ESM 1) ---\n")
cat("(Observed period is identical across scenarios x ESMs)\n\n")
cat(sprintf("  %-6s  %10s  %10s\n",
            "Site", "med_logCHLa", "med_CHLa"))
cat(paste(rep("-", 32), collapse = ""), "\n")

hindcast <- annual_out[annual_out$year >= OBS_START &
                         annual_out$year <= OBS_END &
                         annual_out$scenario == "RCP4.5" &
                         annual_out$esm == CMIP6_ESMS[1], ] %>%
  group_by(site) %>%
  summarise(
    med_logCHLa = median(pred_logCHLa, na.rm = TRUE),
    med_CHLa    = median(pred_CHLa,    na.rm = TRUE),
    .groups = "drop"
  )

for (i in seq_len(nrow(hindcast))) {
  d <- hindcast[i, ]
  cat(sprintf("  %-6s  %10.3f  %10.2f\n",
              d$site, d$med_logCHLa, d$med_CHLa))
}

# ============================================================================
# 8. Time slice summaries
# ============================================================================

cat("\n", strrep("=", 60), "\n", sep = "")
cat(" Phase 6: Time slice summaries\n")
cat(strrep("=", 60), "\n\n")

label_horizon <- function(year) {
  ifelse(year >= HORIZON_2050[1] & year <= HORIZON_2050[2], "2050",
         ifelse(year >= HORIZON_2080[1] & year <= HORIZON_2080[2], "2080", NA))
}

annual_out$horizon <- label_horizon(annual_out$year)
timeslice_in <- annual_out[!is.na(annual_out$horizon), ]

timeslice_out <- timeslice_in %>%
  group_by(site, scenario, horizon) %>%
  summarise(
    logCHLa_med  = median(pred_logCHLa,  na.rm = TRUE),
    logCHLa_q25  = quantile(pred_logCHLa, 0.25, na.rm = TRUE),
    logCHLa_q75  = quantile(pred_logCHLa, 0.75, na.rm = TRUE),
    CHLa_med     = median(pred_CHLa,      na.rm = TRUE),
    CHLa_q25     = quantile(pred_CHLa,    0.25, na.rm = TRUE),
    CHLa_q75     = quantile(pred_CHLa,    0.75, na.rm = TRUE),
    SPC_med      = median(pred_SPC,       na.rm = TRUE),
    logTP_med    = median(pred_logTP,     na.rm = TRUE),
    logTN_med    = median(pred_logTN,     na.rm = TRUE),
    anomaly_med  = median(anomaly,        na.rm = TRUE),
    Q_obs_med    = median(Q_obs_cfs,      na.rm = TRUE),
    Temp_med     = median(Temp_oC,        na.rm = TRUE),
    DSF_med      = median(Days_Since_Freshet, na.rm = TRUE),
    pct_extrap   = 100 * mean(extrap_flag, na.rm = TRUE),
    n_obs        = n(),
    .groups = "drop"
  ) %>%
  as.data.frame()

for (col in names(timeslice_out)) {
  if (is.numeric(timeslice_out[[col]])) {
    digits <- if (grepl("logCHLa|logTP|logTN|anomaly", col)) 3 else 1
    timeslice_out[[col]] <- round(timeslice_out[[col]], digits)
  }
}

timeslice_out <- timeslice_out[order(timeslice_out$scenario,
                                     timeslice_out$horizon,
                                     match(timeslice_out$site, SITE_ORDER)), ]
rownames(timeslice_out) <- NULL

cat("--- Time Slice Bloom Predictions ---\n\n")
cat(sprintf("  %-6s  %-8s  %-6s  %10s  %10s  %10s  %8s\n",
            "Site", "Scenario", "Horiz",
            "CHLa_med", "CHLa_q25", "CHLa_q75", "pct_extrap"))
cat(paste(rep("-", 68), collapse = ""), "\n")
for (i in seq_len(nrow(timeslice_out))) {
  d <- timeslice_out[i, ]
  cat(sprintf("  %-6s  %-8s  %-6s  %10.2f  %10.2f  %10.2f  %8.1f%%\n",
              d$site, d$scenario, d$horizon,
              d$CHLa_med, d$CHLa_q25, d$CHLa_q75, d$pct_extrap))
}

# ============================================================================
# 9. Baseline vs scenario change
# ============================================================================

cat("\n--- Baseline vs Scenario Change ---\n\n")
cat(sprintf("  %-6s  %-8s  %-6s  %10s  %10s  %8s\n",
            "Site", "Scenario", "Horiz",
            "Base_CHLa", "Proj_CHLa", "Delta%"))
cat(paste(rep("-", 58), collapse = ""), "\n")

baseline_vals <- hindcast
names(baseline_vals)[names(baseline_vals) == "med_CHLa"] <- "baseline_CHLa"

ts_compare <- merge(
  timeslice_out[ , c("site", "scenario", "horizon", "CHLa_med")],
  baseline_vals[ , c("site", "baseline_CHLa")],
  by = "site"
)
ts_compare$delta_pct <- 100 * (ts_compare$CHLa_med - ts_compare$baseline_CHLa) /
  ts_compare$baseline_CHLa
ts_compare <- ts_compare[order(ts_compare$scenario, ts_compare$horizon,
                               match(ts_compare$site, SITE_ORDER)), ]

for (i in seq_len(nrow(ts_compare))) {
  d <- ts_compare[i, ]
  cat(sprintf("  %-6s  %-8s  %-6s  %10.2f  %10.2f  %+8.1f%%\n",
              d$site, d$scenario, d$horizon,
              d$baseline_CHLa, d$CHLa_med, d$delta_pct))
}

# ============================================================================
# 10. Write outputs
# ============================================================================

cat("\n", strrep("=", 60), "\n", sep = "")
cat(" Writing outputs\n")
cat(strrep("=", 60), "\n\n")

write_csv(annual_out,    file.path(OUT_DIR, "projections_annual.csv"))
write_csv(timeslice_out, file.path(OUT_DIR, "projections_timeslice.csv"))

cat(sprintf("Saved: 2_incremental/projections_annual.csv\n"))
cat(sprintf("  %d rows x %d cols\n", nrow(annual_out), ncol(annual_out)))
cat(sprintf("Saved: 2_incremental/projections_timeslice.csv\n"))
cat(sprintf("  %d rows x %d cols\n", nrow(timeslice_out), ncol(timeslice_out)))

cat("\nNOTE: 1998-2023 uses observed predictor values, identical across\n")
cat("scenarios and ESMs. 2024-2099 uses NCAR CMIP6 hydrology + NCCV temp.\n")
cat("Scenarios diverge starting in 2024.\n")
cat("\nDone.\n")