# ============================================================================
# 05_process_usgs.R
# UCFR Filamentous Algae Project
# Stage 4: Process USGS daily discharge data
#
# Inputs:   0_data/usgs_daily_q_raw.csv
#           2_incremental/vnrp_sliding_window.csv
# Output:   2_incremental/usgs_processed.csv
#
# Steps:
#   1. Read raw USGS data
#   2. Build HU composite discharge (12352500 + 12353000)
#   3. Map gauge numbers to site codes
#   4. Calculate baseflow — mean daily Q over reference period (flexible)
#   5. Calculate annual peak Q and date of peak per site per year
#   6. Calculate anomaly — (peak Q / baseflow) ^ (1/3)
#   7. Extract observation-specific Q for each biological sampling date
#   8. Calculate Days_Since_Freshet
#   9. Assemble and write final dataset
#
# Notes:
#   - Baseflow reference period is flexible via BASELINE_START / BASELINE_END
#   - HU composite is summed before any calculations
#   - Anomaly baseline is fixed at publication time, not recalculated dynamically
#   - D1 fix (2026-06): section 7 previously approximated the bio sampling
#     date as the 1st of the sampling month (no exact date was available
#     upstream). vnrp_sliding_window.csv now carries a real Date_sample
#     column (recovered in 03_process_vnrp.R, carried through
#     04_lag_selection.R), so Q_obs_cfs and Days_Since_Freshet are now
#     computed from the true sampling day. This also means the 21
#     previously-hidden double-visit bio observations (see 03/04 D1 notes)
#     each get their OWN distinct Q_obs_cfs / Days_Since_Freshet, which is
#     the whole point of having recovered them as separate rows.
# ============================================================================

library(readr)
library(dplyr)

# ----------------------------------------------------------------------------
# 0. Configuration — adjust here if baseline period needs to change
# ----------------------------------------------------------------------------

BASELINE_START <- 1998  # first year of baseline period
BASELINE_END   <- 2025  # last year of baseline period (fixed at publication)

# ----------------------------------------------------------------------------
# 1. Read raw USGS data and biological sampling dates
# ----------------------------------------------------------------------------

cat("Reading USGS raw discharge data...\n")
q_raw <- read_csv("0_data/usgs_daily_q_raw.csv",
                  col_types = cols(.default = col_character()),
                  show_col_types = FALSE)

q_raw$Date  <- as.Date(q_raw$Date)
q_raw$Year  <- as.integer(q_raw$Year)
q_raw$Month <- as.integer(q_raw$Month)
q_raw$Day   <- as.integer(q_raw$Day)
q_raw$Q_cfs <- as.numeric(q_raw$Q_cfs)

cat("Raw rows:", nrow(q_raw), "\n")
cat("Gauges:  ", paste(sort(unique(q_raw$site_no)), collapse = ", "), "\n\n")

cat("Reading biological sampling dates...\n")
bio <- read_csv("2_incremental/vnrp_sliding_window.csv",
                show_col_types = FALSE)
bio$date_yearmon <- as.character(bio$date_yearmon)
bio$Date_sample  <- as.Date(bio$Date_sample)

n_missing_date <- sum(is.na(bio$Date_sample))
if (n_missing_date > 0) {
  warning(sprintf(
    "%d biological rows have no Date_sample -- Q_obs_cfs will be NA for these.",
    n_missing_date))
}

cat("Biological sampling rows:", nrow(bio), "\n")
cat("Rows with real Date_sample:", sum(!is.na(bio$Date_sample)), "\n\n")

# ----------------------------------------------------------------------------
# 2. Build HU composite discharge
# ----------------------------------------------------------------------------
# HU = Clark Fork below Missoula (12353000) + Bitterroot River (12352500)
# Sum both gauges for every date, then treat as a single site

cat("Building HU composite discharge...\n")

gauge_bitterroot  <- "12352500"
gauge_below_miss  <- "12353000"

hu_bitterroot <- q_raw[q_raw$site_no == gauge_bitterroot,
                       c("Date", "Year", "Month", "Day", "Q_cfs")]
hu_below_miss <- q_raw[q_raw$site_no == gauge_below_miss,
                       c("Date", "Year", "Month", "Day", "Q_cfs")]

hu_composite <- merge(hu_bitterroot, hu_below_miss,
                      by = c("Date", "Year", "Month", "Day"),
                      suffixes = c("_bitterroot", "_below_miss"))

hu_composite$Q_cfs    <- rowSums(hu_composite[ , c("Q_cfs_bitterroot",
                                                   "Q_cfs_below_miss")],
                                 na.rm = TRUE)
hu_composite$site_no  <- "HU_composite"
hu_composite$Q_code   <- NA_character_

hu_composite <- hu_composite[ , c("site_no", "Date", "Year",
                                  "Month", "Day", "Q_cfs", "Q_code")]

cat(sprintf("HU composite: %d daily records (%s to %s)\n",
            nrow(hu_composite),
            min(hu_composite$Date),
            max(hu_composite$Date)))

# ----------------------------------------------------------------------------
# 3. Map gauge numbers to site codes
# ----------------------------------------------------------------------------

gauge_site_map <- c(
  "12324200"    = "DL",
  "12324400"    = "GR",
  "12331800"    = "BN",
  "12340500"    = "MS",
  "12353000"    = "BM",
  "12354500"    = "FH",
  "HU_composite" = "HU"
)

# Remove Bitterroot gauge — it's only used in HU composite
q_sites <- q_raw[q_raw$site_no != gauge_bitterroot, ]

# Append HU composite
q_sites <- rbind(q_sites, hu_composite)

# Map to site codes
q_sites$Site <- gauge_site_map[q_sites$site_no]

# Drop any unmapped rows
n_before <- nrow(q_sites)
q_sites  <- q_sites[!is.na(q_sites$Site), ]
cat(sprintf("\nSite mapping: %d rows dropped (unmapped gauges)\n",
            n_before - nrow(q_sites)))
cat("Sites:", paste(sort(unique(q_sites$Site)), collapse = ", "), "\n\n")

# ----------------------------------------------------------------------------
# 4. Calculate baseflow — mean daily Q over reference period
# ----------------------------------------------------------------------------

cat(sprintf("Calculating baseflow (reference period: %d-%d)...\n",
            BASELINE_START, BASELINE_END))

baseline <- q_sites[q_sites$Year >= BASELINE_START &
                      q_sites$Year <= BASELINE_END, ]

baseflow <- baseline %>%
  group_by(Site) %>%
  summarise(
    Q_baseflow_cfs = mean(Q_cfs, na.rm = TRUE),
    baseline_n     = sum(!is.na(Q_cfs)),
    .groups        = "drop"
  )

cat("\n--- Baseflow by Site ---\n")
for (i in seq_len(nrow(baseflow))) {
  cat(sprintf("  %s  mean Q = %8.1f cfs  n = %d days\n",
              baseflow$Site[i],
              baseflow$Q_baseflow_cfs[i],
              baseflow$baseline_n[i]))
}

# ----------------------------------------------------------------------------
# 5. Calculate annual peak Q and date of peak
# ----------------------------------------------------------------------------

cat("\nCalculating annual peak discharge...\n")

annual_peak <- q_sites %>%
  group_by(Site, Year) %>%
  slice_max(Q_cfs, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Site, Year, Q_peak_cfs = Q_cfs, Date_peak = Date)

cat("Annual peak records:", nrow(annual_peak), "\n")

# ----------------------------------------------------------------------------
# 6. Calculate anomaly — (peak Q / baseflow) ^ (1/3)
# ----------------------------------------------------------------------------

cat("Calculating anomaly...\n")

annual_peak <- left_join(annual_peak, baseflow[ , c("Site", "Q_baseflow_cfs")],
                         by = "Site")

annual_peak$anomaly <- (annual_peak$Q_peak_cfs /
                          annual_peak$Q_baseflow_cfs) ^ (1/3)

cat("\n--- Anomaly range by site ---\n")
for (s in sort(unique(annual_peak$Site))) {
  d <- annual_peak[annual_peak$Site == s, ]
  cat(sprintf("  %s  min = %.3f  max = %.3f  mean = %.3f\n",
              s,
              min(d$anomaly, na.rm = TRUE),
              max(d$anomaly, na.rm = TRUE),
              mean(d$anomaly, na.rm = TRUE)))
}

# ----------------------------------------------------------------------------
# 7. Extract observation-specific Q for each biological sampling date
# ----------------------------------------------------------------------------
# D1 fix: bio$Date_sample is now the REAL sampling date (recovered in
# 03_process_vnrp.R, carried through 04_lag_selection.R) -- no longer
# approximated as the 1st of the sampling month. Each of the 21 recovered
# double-visit observations gets its own distinct Q_obs_cfs from its own
# real date, which is the whole point of having recovered them.

cat("\nExtracting observation-specific discharge...\n")

# Build a site-date lookup from q_sites
q_lookup <- q_sites[ , c("Site", "Date", "Q_cfs")]
names(q_lookup)[3] <- "Q_obs_cfs"

obs_q <- left_join(bio[ , c("Site", "Year", "Month", "Date_sample")],
                   q_lookup,
                   by = c("Site", "Date_sample" = "Date"))

cat(sprintf("Observation-specific Q matched: %d of %d rows\n",
            sum(!is.na(obs_q$Q_obs_cfs)), nrow(obs_q)))

n_date_no_q <- sum(!is.na(obs_q$Date_sample) & is.na(obs_q$Q_obs_cfs))
if (n_date_no_q > 0) {
  warning(sprintf(
    paste0(n_date_no_q,
           " rows have a real Date_sample but no matching USGS discharge",
           " record for that exact date -- check for gauge data gaps.")))
}

# ----------------------------------------------------------------------------
# 8. Calculate Days_Since_Freshet
# ----------------------------------------------------------------------------

cat("Calculating Days_Since_Freshet...\n")

# Join peak date onto bio by Site + Year
obs_q <- left_join(obs_q,
                   annual_peak[ , c("Site", "Year", "Date_peak")],
                   by = c("Site", "Year"))

obs_q$Days_Since_Freshet <- as.integer(obs_q$Date_sample - obs_q$Date_peak)

# Flag negative values
n_neg <- sum(obs_q$Days_Since_Freshet < 0, na.rm = TRUE)
if (n_neg > 0) {
  warning(sprintf("%d rows have negative Days_Since_Freshet — sampled before peak",
                  n_neg))
}

cat(sprintf("Days_Since_Freshet range: %d to %d days\n",
            min(obs_q$Days_Since_Freshet, na.rm = TRUE),
            max(obs_q$Days_Since_Freshet, na.rm = TRUE)))

# ----------------------------------------------------------------------------
# 9. Assemble final dataset
# ----------------------------------------------------------------------------

cat("\nAssembling final dataset...\n")

# Join annual summaries onto bio. Join key is now Site + Year + Month +
# Date_sample (not just Site+Year+Month) so that a double-visit month's
# two rows each retain their own distinct Q_obs_cfs / Days_Since_Freshet
# from obs_q, rather than colliding on a shared Site+Year+Month key.
final <- bio %>%
  left_join(annual_peak[ , c("Site", "Year", "Q_peak_cfs",
                             "Q_baseflow_cfs", "anomaly", "Date_peak")],
            by = c("Site", "Year")) %>%
  left_join(obs_q[ , c("Site", "Year", "Month", "Date_sample", "Q_obs_cfs",
                       "Days_Since_Freshet")],
            by = c("Site", "Year", "Month", "Date_sample"))

cat("Final dataset dimensions:", nrow(final), "rows x", ncol(final), "cols\n")

n_row_check <- nrow(final) == nrow(bio)
cat(sprintf("Row count preserved (no join fan-out): %s\n",
            ifelse(n_row_check, "YES", "NO -- INVESTIGATE")))

cat("\n--- Final Dataset Summary ---\n")
cat(sprintf("  %-25s  n = %d\n", "Q_peak_cfs",       sum(!is.na(final$Q_peak_cfs))))
cat(sprintf("  %-25s  n = %d\n", "Q_baseflow_cfs",   sum(!is.na(final$Q_baseflow_cfs))))
cat(sprintf("  %-25s  n = %d\n", "Q_obs_cfs",        sum(!is.na(final$Q_obs_cfs))))
cat(sprintf("  %-25s  n = %d\n", "anomaly",          sum(!is.na(final$anomaly))))
cat(sprintf("  %-25s  n = %d\n", "Days_Since_Freshet",
            sum(!is.na(final$Days_Since_Freshet))))

# ----------------------------------------------------------------------------
# 10. Write output
# ----------------------------------------------------------------------------

out_dir  <- "2_incremental"
out_file <- file.path(out_dir, "usgs_processed.csv")

write_csv(final, out_file)
cat("\nSaved to:", out_file, "\n")
cat("Done.\n")