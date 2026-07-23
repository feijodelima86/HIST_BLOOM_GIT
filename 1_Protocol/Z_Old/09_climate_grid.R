# ============================================================================
# 09_climate_grid.R
# UCFR Filamentous Algae Project
# Stage 9: Synthesize climate scenario grid from NCCV and NCAR mizuRoute
#
# Inputs:
#   0_data/UCF_HUC17010201_MMM_english.csv   (NCCV upper Clark Fork)
#   0_data/MCF_HUC17010204_MMM_english.csv   (NCCV middle Clark Fork)
#   2_incremental/ncar_processed.csv         (NCAR CMIP6 annual metrics)
#
# Output:
#   2_incremental/climate_scenario_grid.csv
#     One row per scenario x horizon x site (2 x 2 x 7 = 28 rows)
#     Columns: site, scenario, horizon,
#              temp_C,                  (absolute, from NCCV)
#              delta_summer_q_pct,      (% change from NCCV runoff, JJA)
#              Q_peak_cfs_med,          (median across ESMs, from NCAR)
#              Q_peak_cfs_q25,
#              Q_peak_cfs_q75,
#              Q_baseflow_cfs_med,
#              Q_baseflow_cfs_q25,
#              Q_baseflow_cfs_q75,
#              DOY_peak_med,
#              DOY_peak_q25,
#              DOY_peak_q75,
#              anomaly_med,
#              anomaly_q25,
#              anomaly_q75
#
# Scenario definitions:
#   ssp245 -> "RCP4.5"  (low emissions)
#   ssp585 -> "RCP8.5"  (high emissions)
#   ssp370 dropped — no clean CMIP5 analog, not used
#
# CMIP generation:
#   CMIP6 only: CanESM5, CMCC-CM2-SR5, MIROC-ES2L,
#               MPI-M.MPI-ESM1-2-LR, NorESM2-MM
#   CMIP5 dropped for internal consistency with NCCV (also CMIP6)
#
# Baseline period: 1998-2022 (matches biological observation record)
# Horizon windows: 2040-2060 -> "2050", 2070-2090 -> "2080"
#
# Site-to-NCCV region assignment:
#   Upper Clark Fork (UCF): DL, GR, BN, MS
#   Middle Clark Fork (MCF): BM, HU, FH
#
# Notes:
#   - NCAR year 1950 discarded (corrupt values)
#   - Temperature converted from deg_F to deg_C
#   - NCCV runoff used as % change proxy for summer mean Q (JJA)
#   - NCAR provides peak Q, baseflow Q, DOY_peak, anomaly
#   - All delta values expressed relative to 1998-2022 baseline
# ============================================================================

library(readr)
library(dplyr)

# ============================================================================
# 0. Configuration
# ============================================================================

BASELINE_START <- 1998
BASELINE_END   <- 2022

HORIZON_2050   <- c(2040, 2060)
HORIZON_2080   <- c(2070, 2090)

# CMIP6 ESMs only
CMIP6_ESMS <- c("CanESM5", "CMCC-CM2-SR5", "MIROC-ES2L",
                "MPI-M.MPI-ESM1-2-LR", "NorESM2-MM")

# Scenario mapping: NCAR scenario label -> output label
SCENARIO_MAP <- c(ssp245 = "RCP4.5", ssp585 = "RCP8.5")

# Site-to-region assignment
UCF_SITES <- c("DL", "GR", "BN", "MS")
MCF_SITES <- c("BM", "HU", "FH")
SITE_ORDER <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")

# Summer months for runoff (JJA)
SUMMER_MONTHS <- 6:8

# Output paths
OUT_DIR  <- "2_incremental"
OUT_FILE <- file.path(OUT_DIR, "climate_scenario_grid.csv")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Helper: Fahrenheit to Celsius
# ============================================================================

f_to_c <- function(f) (f - 32) * 5 / 9

# ============================================================================
# Helper: quantile summary for a numeric vector
# ============================================================================

qsum <- function(x, na.rm = TRUE) {
  list(
    med = median(x, na.rm = na.rm),
    q25 = quantile(x, 0.25, na.rm = na.rm),
    q75 = quantile(x, 0.75, na.rm = na.rm)
  )
}

# ============================================================================
# Phase 1: Process NCCV data
# ============================================================================

cat(strrep("=", 60), "\n")
cat(" Phase 1: Processing NCCV data\n")
cat(strrep("=", 60), "\n\n")

read_nccv <- function(path, region_label) {
  cat(sprintf("Reading %s (%s)...\n", basename(path), region_label))
  
  raw <- read_csv(path, show_col_types = FALSE)
  
  # Parse date — format is " 1/15/1950" (leading space, month/day/year)
  raw$Date <- as.Date(trimws(raw$Date), format = "%m/%d/%Y")
  raw$Year  <- as.integer(format(raw$Date, "%Y"))
  raw$Month <- as.integer(format(raw$Date, "%m"))
  
  cat(sprintf("  Rows: %d  |  Date range: %s to %s\n",
              nrow(raw),
              min(raw$Date, na.rm = TRUE),
              max(raw$Date, na.rm = TRUE)))
  
  # Standardize column names: remove leading/trailing spaces
  names(raw) <- trimws(names(raw))
  
  # Rename columns to short names for clarity
  # Pattern: "ssp245 Mean temperature (deg_F)" -> "ssp245_temp_mean_F"
  col_map <- c(
    "ssp245 Mean temperature (deg_F)"  = "ssp245_temp_mean_F",
    "ssp245 Runoff (in/mo)"            = "ssp245_runoff",
    "ssp370 Mean temperature (deg_F)"  = "ssp370_temp_mean_F",
    "ssp370 Runoff (in/mo)"            = "ssp370_runoff",
    "ssp585 Mean temperature (deg_F)"  = "ssp585_temp_mean_F",
    "ssp585 Runoff (in/mo)"            = "ssp585_runoff"
  )
  
  for (old in names(col_map)) {
    if (old %in% names(raw)) {
      names(raw)[names(raw) == old] <- col_map[old]
    }
  }
  
  raw$region <- region_label
  raw
}

ucf <- read_nccv("0_data/UCF_HUC17010201_MMM_english.csv", "UCF")
mcf <- read_nccv("0_data/MCF_HUC17010204_MMM_english.csv", "MCF")

# ============================================================================
# Phase 1a: Compute NCCV baseline and horizon summaries
# ============================================================================

summarize_nccv <- function(df, region_label) {
  
  # Baseline: annual mean temperature and annual sum of summer runoff
  bl <- df[df$Year >= BASELINE_START & df$Year <= BASELINE_END, ]
  
  # Summer runoff baseline (JJA mean across years)
  bl_summer <- bl[bl$Month %in% SUMMER_MONTHS, ]
  
  # Annual summer runoff per year (sum of JJA months)
  bl_summer_annual_245 <- tapply(bl_summer$ssp245_runoff, bl_summer$Year, sum, na.rm = TRUE)
  bl_summer_annual_585 <- tapply(bl_summer$ssp585_runoff, bl_summer$Year, sum, na.rm = TRUE)
  
  baseline_summer_runoff_245 <- mean(bl_summer_annual_245, na.rm = TRUE)
  baseline_summer_runoff_585 <- mean(bl_summer_annual_585, na.rm = TRUE)
  
  # Annual mean temperature baseline (all months)
  baseline_temp_245 <- mean(bl$ssp245_temp_mean_F, na.rm = TRUE)
  baseline_temp_585 <- mean(bl$ssp585_temp_mean_F, na.rm = TRUE)
  
  cat(sprintf("\n  %s baseline (%d-%d):\n", region_label,
              BASELINE_START, BASELINE_END))
  cat(sprintf("    SSP245 mean temp: %.2f F (%.2f C)\n",
              baseline_temp_245, f_to_c(baseline_temp_245)))
  cat(sprintf("    SSP585 mean temp: %.2f F (%.2f C)\n",
              baseline_temp_585, f_to_c(baseline_temp_585)))
  cat(sprintf("    SSP245 summer runoff: %.3f in/mo\n",
              baseline_summer_runoff_245))
  cat(sprintf("    SSP585 summer runoff: %.3f in/mo\n",
              baseline_summer_runoff_585))
  
  # Horizon summaries
  results <- data.frame()
  
  for (h in list(list(label = "2050", yrs = HORIZON_2050),
                 list(label = "2080", yrs = HORIZON_2080))) {
    
    hor <- df[df$Year >= h$yrs[1] & df$Year <= h$yrs[2], ]
    hor_summer <- hor[hor$Month %in% SUMMER_MONTHS, ]
    
    for (sc in c("ssp245", "ssp585")) {
      
      # Skip ssp370
      temp_col   <- paste0(sc, "_temp_mean_F")
      runoff_col <- paste0(sc, "_runoff")
      
      if (!temp_col %in% names(hor)) next
      
      # Horizon mean temperature
      hor_temp_F <- mean(hor[[temp_col]], na.rm = TRUE)
      hor_temp_C <- f_to_c(hor_temp_F)
      
      # Horizon summer runoff (annual sum per year, then mean)
      hor_summer_annual <- tapply(hor_summer[[runoff_col]], hor_summer$Year,
                                  sum, na.rm = TRUE)
      hor_summer_runoff <- mean(hor_summer_annual, na.rm = TRUE)
      
      # % change in summer runoff relative to baseline
      bl_runoff <- if (sc == "ssp245") baseline_summer_runoff_245 else baseline_summer_runoff_585
      delta_summer_q_pct <- 100 * (hor_summer_runoff - bl_runoff) / bl_runoff
      
      cat(sprintf("    %s %s: temp=%.2f C  summer_runoff=%.3f in/mo  delta=%.1f%%\n",
                  sc, h$label, hor_temp_C, hor_summer_runoff, delta_summer_q_pct))
      
      results <- rbind(results, data.frame(
        region             = region_label,
        scenario_ncar      = sc,
        scenario           = SCENARIO_MAP[sc],
        horizon            = h$label,
        temp_C             = round(hor_temp_C, 3),
        delta_summer_q_pct = round(delta_summer_q_pct, 2),
        stringsAsFactors   = FALSE
      ))
    }
  }
  results
}

cat("\n")
ucf_summary <- summarize_nccv(ucf, "UCF")
mcf_summary <- summarize_nccv(mcf, "MCF")

nccv_summary <- rbind(ucf_summary, mcf_summary)

cat("\n--- NCCV Summary ---\n")
print(as.data.frame(nccv_summary))

# ============================================================================
# Phase 2: Process NCAR mizuRoute data
# ============================================================================

cat("\n", strrep("=", 60), "\n", sep = "")
cat(" Phase 2: Processing NCAR mizuRoute data\n")
cat(strrep("=", 60), "\n\n")

cat("Reading ncar_processed.csv...\n")
ncar <- as.data.frame(read_csv("2_incremental/ncar_processed.csv",
                               show_col_types = FALSE))

cat(sprintf("  Rows: %d\n", nrow(ncar)))
cat(sprintf("  ESMs: %s\n", paste(sort(unique(ncar$esm)), collapse = ", ")))
cat(sprintf("  Scenarios: %s\n", paste(sort(unique(ncar$scenario)), collapse = ", ")))
cat(sprintf("  Year range: %d to %d\n", min(ncar$year), max(ncar$year)))
cat(sprintf("  Sites: %s\n\n", paste(sort(unique(ncar$site)), collapse = ", ")))

# Filter: CMIP6 only, drop ssp370, discard year 1950
ncar_filt <- ncar[
  ncar$esm      %in% CMIP6_ESMS &
    ncar$scenario %in% names(SCENARIO_MAP) &
    ncar$year     != 1950,
]

cat(sprintf("After filtering (CMIP6, ssp245/585, no 1950): %d rows\n\n",
            nrow(ncar_filt)))

# Verify ESM coverage
cat("ESM x scenario coverage after filtering:\n")
xt <- table(ncar_filt$esm, ncar_filt$scenario)
print(xt)
cat("\n")

# ============================================================================
# Phase 2a: Baseline from NCAR (1998-2022)
# ============================================================================

cat("Computing NCAR baselines (site x ESM x scenario)...\n")

ncar_baseline <- ncar_filt[
  ncar_filt$year >= BASELINE_START & ncar_filt$year <= BASELINE_END,
] %>%
  group_by(site, esm, scenario) %>%
  summarise(
    baseline_Q_peak_cfs    = median(Q_peak_cfs,        na.rm = TRUE),
    baseline_Q_baseflow    = median(Q_baseflow_cfs,     na.rm = TRUE),
    baseline_DOY_peak      = median(DOY_peak,           na.rm = TRUE),
    baseline_anomaly       = median(anomaly,            na.rm = TRUE),
    baseline_n             = sum(!is.na(Q_peak_cfs)),
    .groups = "drop"
  )

cat(sprintf("  Baseline records: %d\n\n", nrow(ncar_baseline)))

# ============================================================================
# Phase 2b: Horizon summaries from NCAR
# ============================================================================

cat("Computing NCAR horizon summaries...\n")

ncar_results <- data.frame()

for (h in list(list(label = "2050", yrs = HORIZON_2050),
               list(label = "2080", yrs = HORIZON_2080))) {
  
  hor <- ncar_filt[ncar_filt$year >= h$yrs[1] & ncar_filt$year <= h$yrs[2], ]
  
  # Summarize per site x ESM x scenario for this horizon
  hor_sum <- hor %>%
    group_by(site, esm, scenario) %>%
    summarise(
      Q_peak_cfs   = median(Q_peak_cfs,    na.rm = TRUE),
      Q_baseflow   = median(Q_baseflow_cfs, na.rm = TRUE),
      DOY_peak     = median(DOY_peak,       na.rm = TRUE),
      anomaly      = median(anomaly,        na.rm = TRUE),
      .groups = "drop"
    )
  
  hor_sum$horizon <- h$label
  
  ncar_results <- rbind(ncar_results, as.data.frame(hor_sum))
}

# Now aggregate across ESMs: median + IQR per site x scenario x horizon
ncar_agg <- ncar_results %>%
  group_by(site, scenario, horizon) %>%
  summarise(
    Q_peak_cfs_med    = median(Q_peak_cfs,   na.rm = TRUE),
    Q_peak_cfs_q25    = quantile(Q_peak_cfs, 0.25, na.rm = TRUE),
    Q_peak_cfs_q75    = quantile(Q_peak_cfs, 0.75, na.rm = TRUE),
    Q_baseflow_med    = median(Q_baseflow,    na.rm = TRUE),
    Q_baseflow_q25    = quantile(Q_baseflow,  0.25, na.rm = TRUE),
    Q_baseflow_q75    = quantile(Q_baseflow,  0.75, na.rm = TRUE),
    DOY_peak_med      = median(DOY_peak,      na.rm = TRUE),
    DOY_peak_q25      = quantile(DOY_peak,    0.25, na.rm = TRUE),
    DOY_peak_q75      = quantile(DOY_peak,    0.75, na.rm = TRUE),
    anomaly_med       = median(anomaly,       na.rm = TRUE),
    anomaly_q25       = quantile(anomaly,     0.25, na.rm = TRUE),
    anomaly_q75       = quantile(anomaly,     0.75, na.rm = TRUE),
    n_esm             = n(),
    .groups = "drop"
  )

# Remap scenario labels
ncar_agg$scenario <- SCENARIO_MAP[ncar_agg$scenario]

cat(sprintf("  NCAR aggregated rows: %d\n\n", nrow(ncar_agg)))

cat("--- NCAR Horizon Summary (median across ESMs) ---\n")
cat(sprintf("  %-6s  %-8s  %-6s  %10s  %10s  %8s  %8s\n",
            "Site", "Scenario", "Horiz", "Q_peak_med", "Q_base_med",
            "DOY_med", "anom_med"))
cat(paste(rep("-", 66), collapse = ""), "\n")

for (i in seq_len(nrow(ncar_agg))) {
  d <- ncar_agg[i, ]
  cat(sprintf("  %-6s  %-8s  %-6s  %10.1f  %10.1f  %8.1f  %8.3f\n",
              d$site, d$scenario, d$horizon,
              d$Q_peak_cfs_med, d$Q_baseflow_med,
              d$DOY_peak_med, d$anomaly_med))
}

# ============================================================================
# Phase 2c: Compute NCAR baseline summary for comparison
# ============================================================================

ncar_bl_agg <- ncar_baseline %>%
  group_by(site, scenario) %>%
  summarise(
    baseline_Q_peak_med = median(baseline_Q_peak_cfs, na.rm = TRUE),
    baseline_Q_base_med = median(baseline_Q_baseflow,  na.rm = TRUE),
    baseline_DOY_med    = median(baseline_DOY_peak,    na.rm = TRUE),
    baseline_anom_med   = median(baseline_anomaly,     na.rm = TRUE),
    .groups = "drop"
  )

ncar_bl_agg$scenario <- SCENARIO_MAP[ncar_bl_agg$scenario]

cat("\n--- NCAR Baseline Summary (1998-2022 median) ---\n")
cat(sprintf("  %-6s  %-8s  %10s  %10s  %8s  %8s\n",
            "Site", "Scenario", "Q_peak_med", "Q_base_med",
            "DOY_med", "anom_med"))
cat(paste(rep("-", 58), collapse = ""), "\n")

for (i in seq_len(nrow(ncar_bl_agg))) {
  d <- ncar_bl_agg[i, ]
  cat(sprintf("  %-6s  %-8s  %10.1f  %10.1f  %8.1f  %8.3f\n",
              d$site, d$scenario,
              d$baseline_Q_peak_med, d$baseline_Q_base_med,
              d$baseline_DOY_med, d$baseline_anom_med))
}

# ============================================================================
# Phase 3: Assemble final scenario grid
# ============================================================================

cat("\n", strrep("=", 60), "\n", sep = "")
cat(" Phase 3: Assembling climate scenario grid\n")
cat(strrep("=", 60), "\n\n")

# Expand NCCV summaries to site level
nccv_sites <- data.frame()

for (s in SITE_ORDER) {
  region <- if (s %in% UCF_SITES) "UCF" else "MCF"
  d <- nccv_summary[nccv_summary$region == region, ]
  d$site <- s
  nccv_sites <- rbind(nccv_sites, d)
}

# Merge NCCV and NCAR
grid <- merge(
  nccv_sites[ , c("site", "scenario", "horizon",
                  "temp_C", "delta_summer_q_pct")],
  ncar_agg[ , c("site", "scenario", "horizon",
                "Q_peak_cfs_med", "Q_peak_cfs_q25", "Q_peak_cfs_q75",
                "Q_baseflow_med", "Q_baseflow_q25", "Q_baseflow_q75",
                "DOY_peak_med", "DOY_peak_q25", "DOY_peak_q75",
                "anomaly_med", "anomaly_q25", "anomaly_q75",
                "n_esm")],
  by = c("site", "scenario", "horizon"),
  all.x = TRUE
)

# Order columns and rows
grid <- grid[order(grid$scenario, grid$horizon,
                   match(grid$site, SITE_ORDER)), ]
rownames(grid) <- NULL

# Round numeric columns
num_cols <- c("temp_C", "delta_summer_q_pct",
              "Q_peak_cfs_med", "Q_peak_cfs_q25", "Q_peak_cfs_q75",
              "Q_baseflow_med", "Q_baseflow_q25", "Q_baseflow_q75",
              "DOY_peak_med", "DOY_peak_q25", "DOY_peak_q75",
              "anomaly_med", "anomaly_q25", "anomaly_q75")

for (col in num_cols) {
  if (col %in% names(grid)) {
    digits <- if (grepl("anomaly|temp|delta", col)) 3 else 1
    grid[[col]] <- round(grid[[col]], digits)
  }
}

# ============================================================================
# Summary table
# ============================================================================

cat("--- Final Climate Scenario Grid ---\n\n")
cat(sprintf("  %-6s  %-8s  %-6s  %6s  %8s  %10s  %10s  %8s  %8s\n",
            "Site", "Scenario", "Horiz", "Temp_C", "dSumQ%",
            "Q_peak_med", "Q_base_med", "DOY_med", "anom_med"))
cat(paste(rep("-", 82), collapse = ""), "\n")

for (i in seq_len(nrow(grid))) {
  d <- grid[i, ]
  cat(sprintf("  %-6s  %-8s  %-6s  %6.2f  %8.1f  %10.1f  %10.1f  %8.1f  %8.3f\n",
              d$site, d$scenario, d$horizon,
              d$temp_C, d$delta_summer_q_pct,
              d$Q_peak_cfs_med, d$Q_baseflow_med,
              d$DOY_peak_med, d$anomaly_med))
}

cat(sprintf("\nGrid dimensions: %d rows x %d cols\n", nrow(grid), ncol(grid)))
cat(sprintf("Scenarios: %s\n", paste(sort(unique(grid$scenario)), collapse = ", ")))
cat(sprintf("Horizons:  %s\n", paste(sort(unique(grid$horizon)), collapse = ", ")))
cat(sprintf("Sites:     %s\n", paste(SITE_ORDER, collapse = ", ")))
cat(sprintf("ESMs per scenario x horizon: %d (CMIP6 only)\n",
            unique(grid$n_esm[!is.na(grid$n_esm)])[1]))

# ============================================================================
# Sanity checks
# ============================================================================

cat("\n--- Sanity Checks ---\n")

# Check expected row count
expected_rows <- length(SITE_ORDER) * 2 * 2  # 7 sites x 2 scenarios x 2 horizons
if (nrow(grid) == expected_rows) {
  cat(sprintf("  [OK] Row count: %d (expected %d)\n", nrow(grid), expected_rows))
} else {
  cat(sprintf("  [WARN] Row count: %d (expected %d)\n", nrow(grid), expected_rows))
}

# Check no NA in key columns
key_cols <- c("temp_C", "delta_summer_q_pct", "Q_peak_cfs_med",
              "Q_baseflow_med", "DOY_peak_med", "anomaly_med")
for (col in key_cols) {
  n_na <- sum(is.na(grid[[col]]))
  status <- if (n_na == 0) "[OK]" else "[WARN]"
  cat(sprintf("  %s %s: %d NA\n", status, col, n_na))
}

# Check direction of change: RCP8.5 should be warmer than RCP4.5
for (s in SITE_ORDER) {
  for (h in c("2050", "2080")) {
    t45 <- grid$temp_C[grid$site == s & grid$scenario == "RCP4.5" & grid$horizon == h]
    t85 <- grid$temp_C[grid$site == s & grid$scenario == "RCP8.5" & grid$horizon == h]
    if (length(t45) > 0 && length(t85) > 0 && !is.na(t45) && !is.na(t85)) {
      if (t85 <= t45) {
        cat(sprintf("  [WARN] %s %s: RCP8.5 temp (%.2f) not warmer than RCP4.5 (%.2f)\n",
                    s, h, t85, t45))
      }
    }
  }
}
cat("  [OK] Temperature direction check complete\n")

# ============================================================================
# Write output
# ============================================================================

write_csv(grid, OUT_FILE)
cat(sprintf("\nSaved -> %s\n", OUT_FILE))
cat(sprintf("  %d rows x %d columns\n", nrow(grid), ncol(grid)))
cat("\nDone.\n")