# ============================================================================
# 07c_ncar_discharge_envelope.R
# UCFR Filamentous Algae Project
# MA20 discharge envelope — production version
#
# Inputs:
#   2_incremental/ncar_daily_q.csv          (from 06f; spin-up trimmed)
#   2_incremental/ncar_water_year_peaks.csv (from 06g; peak Q + timing)
#
# Output:
#   2_incremental/ncar_discharge_envelope.csv
#     One row per site x esm x scenario x water_year.
#     MA20-smoothed metrics ready for delta computation against USGS baseline.
#
# Columns in output:
#   site, esm, scenario, cmip, water_year,
#   peak_q_cfs_raw, baseflow_q_cfs_raw, mean_q_cfs_raw,
#   days_since_wy_start_raw,
#   peak_q_cfs_ma20, baseflow_q_cfs_ma20, mean_q_cfs_ma20,
#   days_since_wy_start_ma20,
#   anomaly_ma20          <- smooth-then-ratio: (peak_ma20 / base_ma20)^(1/3)
#   anomaly_raw           <- ratio-then-smooth: kept for diagnostic contrast only
#
# MA window: 20 years, centered, computed independently per metric per
#   site x member. Coverage guard: require >= 10 non-NA years within the
#   window before returning a smoothed value (otherwise NA).
#
# Smooth-then-ratio (D2 resolution):
#   peak_q_cfs and baseflow_q_cfs are each MA20-smoothed independently.
#   anomaly_ma20 = (peak_q_cfs_ma20 / baseflow_q_cfs_ma20)^(1/3).
#   This differs from the observed pipeline (05_process_usgs.R), where
#   baseflow is a single fixed long-term mean over 1998-2025. The two
#   definitions are analogous but not numerically equivalent — the NCAR
#   pipeline uses a rolling smoothed reference so the ratio tracks
#   decadal shifts in the disturbance regime rather than departures
#   from a fixed historical anchor.
#
# NCAR dataset constraint (Coxon et al.):
#   Individual daily/annual ESM series are NOT valid for event analysis.
#   Dataset mandates decadal+ climatological means computed per ESM
#   individually before any ensemble summary. MA20 satisfies this.
#   Never pool across ESM members before computing per-member metrics.
#
# Baseflow definition:
#   Annual minimum daily Q within the water year (Oct-Sep boundary
#   inherited from 06g's variable-boundary water years).
#   Computed here directly from ncar_daily_q.csv keyed to each
#   water year's wy_start_date / wy_end_date from ncar_water_year_peaks.csv.
#
# Delta computation (downstream, not here):
#   This script produces smoothed absolute values per member.
#   Deltas (future / historical ratio or difference) are computed in the
#   projection script by comparing future MA20 values against the
#   MA20 mean over the 1998-2025 baseline window, per member.
#   Never apply ensemble mean before per-member delta; never use raw
#   annual values as absolute discharge inputs to M1.
# ============================================================================

# ============================================================================
# CONFIGURATION
# ============================================================================
MA_WINDOW    <- 20L     # years; centered; do not change without re-validating
MA_MIN_OBS   <- 10L     # minimum non-NA years within window to return value

IN_DAILY     <- "2_incremental/ncar_daily_q.csv"
IN_PEAKS     <- "2_incremental/ncar_water_year_peaks.csv"
OUT_FILE     <- "2_incremental/ncar_discharge_envelope.csv"

SITES        <- c("CLALO", "CLADR", "CLABE", "CLAPL")
# ============================================================================


# ============================================================================
# 1. LOAD INPUTS
# ============================================================================
daily_df <- read.csv(IN_DAILY,  stringsAsFactors = FALSE)
peaks_df <- read.csv(IN_PEAKS,  stringsAsFactors = FALSE)

daily_df$date          <- as.Date(daily_df$date)
peaks_df$wy_start_date <- as.Date(peaks_df$wy_start_date)
peaks_df$wy_end_date   <- as.Date(peaks_df$wy_end_date)

# Drop flagged water years (degenerate boundaries from 06g)
peaks_df <- peaks_df[!peaks_df$flag_length | is.na(peaks_df$flag_length), ]

members <- unique(peaks_df[, c("esm", "scenario", "cmip")])
members <- members[order(members$scenario, members$esm), ]

cat("Members loaded:", nrow(members), "\n")
cat("Sites:", paste(SITES, collapse = ", "), "\n")
cat("Water year range:",
    min(peaks_df$water_year, na.rm = TRUE), "-",
    max(peaks_df$water_year, na.rm = TRUE), "\n\n")


# ============================================================================
# 2. CENTERED MOVING AVERAGE HELPER
# ============================================================================
# Centered MA over a vector of annual values (one per water year).
# half  = floor(MA_WINDOW / 2)
# For each position i, average values in [i - half, i + half].
# Returns NA where fewer than MA_MIN_OBS non-NA values exist in window.

centered_ma <- function(x, window = MA_WINDOW, min_obs = MA_MIN_OBS) {
  n    <- length(x)
  out  <- rep(NA_real_, n)
  half <- floor(window / 2L)
  for (i in seq_len(n)) {
    lo  <- max(1L, i - half)
    hi  <- min(n,  i + half)
    win <- x[lo:hi]
    ok  <- sum(!is.na(win))
    if (ok >= min_obs) out[i] <- mean(win, na.rm = TRUE)
  }
  out
}


# ============================================================================
# 3. MAIN LOOP: compute raw annual metrics + MA20 per site x member
# ============================================================================
# For each site x member:
#   (a) From peaks_df: inherit peak_q_cfs and days_since_wy_start (already
#       computed by 06g from the variable water-year boundaries).
#   (b) From daily_df: compute per-water-year baseflow (annual min) and
#       mean_q_cfs (annual mean) by subsetting daily data to each WY window.
#   (c) MA20-smooth all four raw metrics independently.
#   (d) Compute anomaly_ma20 = (peak_ma20 / base_ma20)^(1/3).
#   (e) Compute anomaly_raw  = MA20(raw_peak / raw_base)^(1/3) — diagnostic
#       contrast retained per D2 resolution; not used in projections.

all_out <- vector("list", nrow(members) * length(SITES))
k       <- 0L

for (m in seq_len(nrow(members))) {
  esm_i      <- members$esm[m]
  scenario_i <- members$scenario[m]
  cmip_i     <- members$cmip[m]
  
  # Daily series for this member, sorted by date
  daily_m <- daily_df[daily_df$esm      == esm_i &
                        daily_df$scenario == scenario_i, ]
  daily_m <- daily_m[order(daily_m$date), ]
  
  # Peaks/timing for this member
  peaks_m <- peaks_df[peaks_df$esm      == esm_i &
                        peaks_df$scenario == scenario_i, ]
  
  for (site in SITES) {
    k <- k + 1L
    
    site_peaks <- peaks_m[peaks_m$site == site, ]
    site_peaks <- site_peaks[order(site_peaks$water_year), ]
    if (nrow(site_peaks) == 0) next
    
    n_wy <- nrow(site_peaks)
    
    # Vectors for raw annual metrics
    baseflow_raw <- rep(NA_real_, n_wy)
    mean_q_raw   <- rep(NA_real_, n_wy)
    
    daily_q <- daily_m[[site]]
    daily_d <- daily_m$date
    
    for (i in seq_len(n_wy)) {
      wy_start <- site_peaks$wy_start_date[i]
      wy_end   <- site_peaks$wy_end_date[i]
      
      idx <- which(daily_d > wy_start & daily_d <= wy_end)
      if (length(idx) == 0) next
      
      q_win <- daily_q[idx]
      if (all(is.na(q_win))) next
      
      baseflow_raw[i] <- min(q_win,  na.rm = TRUE)
      mean_q_raw[i]   <- mean(q_win, na.rm = TRUE)
    }
    
    # Raw peak and timing from 06g
    peak_raw   <- site_peaks$peak_q_cfs
    timing_raw <- as.numeric(site_peaks$days_since_wy_start)
    
    # Raw anomaly (ratio-then-smooth) — diagnostic only
    anom_ratio_raw <- (peak_raw / baseflow_raw)^(1/3)
    anom_raw_ma20  <- centered_ma(anom_ratio_raw)
    
    # MA20 smooth each metric independently
    peak_ma20     <- centered_ma(peak_raw)
    base_ma20     <- centered_ma(baseflow_raw)
    mean_ma20     <- centered_ma(mean_q_raw)
    timing_ma20   <- centered_ma(timing_raw)
    
    # Smooth-then-ratio anomaly (production)
    anom_ma20 <- (peak_ma20 / base_ma20)^(1/3)
    
    all_out[[k]] <- data.frame(
      site                  = site,
      esm                   = esm_i,
      scenario              = scenario_i,
      cmip                  = cmip_i,
      water_year            = site_peaks$water_year,
      # Raw annual values (for diagnostics and baseline computation)
      peak_q_cfs_raw        = round(peak_raw,     2),
      baseflow_q_cfs_raw    = round(baseflow_raw,  2),
      mean_q_cfs_raw        = round(mean_q_raw,    2),
      days_since_wy_start_raw = timing_raw,
      # MA20-smoothed values (production inputs for delta computation)
      peak_q_cfs_ma20       = round(peak_ma20,    2),
      baseflow_q_cfs_ma20   = round(base_ma20,    2),
      mean_q_cfs_ma20       = round(mean_ma20,    2),
      days_since_wy_start_ma20 = round(timing_ma20, 2),
      # Anomaly: smooth-then-ratio (production) and ratio-then-smooth (diagnostic)
      anomaly_ma20          = round(anom_ma20,    4),
      anomaly_raw           = round(anom_raw_ma20, 4),
      stringsAsFactors = FALSE
    )
  }
}

envelope <- do.call(rbind, all_out)
envelope <- envelope[order(envelope$site, envelope$scenario,
                           envelope$esm, envelope$water_year), ]
rownames(envelope) <- NULL


# ============================================================================
# 4. SAVE OUTPUT
# ============================================================================
write.csv(envelope, OUT_FILE, row.names = FALSE)


# ============================================================================
# 5. SCORECARD
# ============================================================================
cat("============================================================\n")
cat("DISCHARGE ENVELOPE SCORECARD — MA20 production\n")
cat("============================================================\n\n")

# --- 5a. Coverage summary ---
coverage_sc <- data.frame(
  Metric  = c("Members (esm x scenario)",
              "Sites",
              "Total rows",
              "Water year range",
              "MA20 NA rate (peak)",
              "MA20 NA rate (baseflow)",
              "MA20 NA rate (anomaly)"),
  Value   = c(
    nrow(members),
    length(SITES),
    nrow(envelope),
    paste(min(envelope$water_year, na.rm = TRUE), "-",
          max(envelope$water_year, na.rm = TRUE)),
    paste0(round(mean(is.na(envelope$peak_q_cfs_ma20))    * 100, 1), "%"),
    paste0(round(mean(is.na(envelope$baseflow_q_cfs_ma20)) * 100, 1), "%"),
    paste0(round(mean(is.na(envelope$anomaly_ma20))        * 100, 1), "%")
  ),
  stringsAsFactors = FALSE
)
cat("--- Coverage ---\n")
print(coverage_sc, row.names = FALSE)
cat("\n")

# --- 5b. Per-site MA20 means across full record (sanity check vs USGS) ---
site_means <- do.call(rbind, lapply(SITES, function(s) {
  d <- envelope[envelope$site == s, ]
  data.frame(
    Site             = s,
    peak_ma20_mean   = round(mean(d$peak_q_cfs_ma20,    na.rm = TRUE), 0),
    base_ma20_mean   = round(mean(d$baseflow_q_cfs_ma20, na.rm = TRUE), 0),
    mean_ma20_mean   = round(mean(d$mean_q_cfs_ma20,    na.rm = TRUE), 0),
    anomaly_ma20_mean = round(mean(d$anomaly_ma20,       na.rm = TRUE), 3),
    stringsAsFactors = FALSE
  )
}))
cat("--- Per-site ensemble-mean MA20 values (all years, all members) ---\n")
cat("(Sanity check: compare peak/mean against known USGS basin-scale Q)\n\n")
print(site_means, row.names = FALSE)
cat("\n")

# --- 5c. Anomaly drift check: early (1960-1979) vs late (2060-2079) ---
anom_drift <- do.call(rbind, lapply(SITES, function(s) {
  d     <- envelope[envelope$site == s & !is.na(envelope$anomaly_ma20), ]
  early <- mean(d$anomaly_ma20[d$water_year >= 1960 & d$water_year <= 1979],
                na.rm = TRUE)
  late  <- mean(d$anomaly_ma20[d$water_year >= 2060 & d$water_year <= 2079],
                na.rm = TRUE)
  data.frame(
    Site         = s,
    anomaly_1960_79 = round(early, 3),
    anomaly_2060_79 = round(late,  3),
    delta        = round(late - early, 3),
    direction    = ifelse(late > early, "increasing", "declining"),
    stringsAsFactors = FALSE
  )
}))
cat("--- Anomaly drift (early vs late century, ensemble mean) ---\n")
cat("(Should show mild decline per D2 ecological finding)\n\n")
print(anom_drift, row.names = FALSE)
cat("\n")

# --- 5d. Output file ---
sc_files <- data.frame(
  File     = OUT_FILE,
  Contents = paste0("MA20 discharge envelope. ",
                    nrow(envelope), " rows x ", ncol(envelope), " cols. ",
                    "Feed to delta computation in projection script."),
  stringsAsFactors = FALSE
)
cat("--- Output ---\n")
print(sc_files, row.names = FALSE)
cat("\n")
cat("============================================================\n")
cat("Done. Canonical discharge envelope for projection pipeline.\n")
cat("Next: 07d_ncar_temperature_envelope.R (NorthWEST pairing).\n")
cat("============================================================\n")