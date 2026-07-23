# ============================================================================
# explore_ma_window_sensitivity.R
# UCFR Filamentous Algae Project
# Exploratory: derive water-year hydroclimate metrics (peak magnitude/timing,
# baseflow, anomaly, mean Q) per NCAR ESM member, then compare several
# centered moving-average window lengths against fixed 30-year bin means,
# across all four NCAR reaches (CLALO, CLADR, CLABE, CLAPL).
#
# Rationale: GDEX d010014 documentation (Coxon et al., Nature Sci Data)
# explicitly warns against using individual daily/annual values, and against
# averaging across ESMs before computing metrics. This script respects both:
# all smoothing is done WITHIN one ESM member's own trajectory, never across
# members, and the whole point is to find a window length that turns
# individually-unreliable annual values into a trustworthy long-term signal.
#
# Reads:
#   2_incremental/ncar_daily_q.csv          (from 06f)
#   2_incremental/ncar_water_year_peaks.csv (from 06g; peak date/magnitude)
#
# Output:
#   2_incremental/wy_metrics_all_sites.csv         (per-wy derived metrics)
#   4_products/diag_ma_window_sensitivity_upper.png (CLALO, CLADR)
#   4_products/diag_ma_window_sensitivity_lower.png (CLABE, CLAPL)
#   scorecard (printed, no file)
# ============================================================================

suppressPackageStartupMessages({
  library(readr)
})

# ============================================================================
# Configuration
# ============================================================================

SITES   <- c("CLALO", "CLADR", "CLABE", "CLAPL")
MA_WINDOWS <- c(10, 15, 20, 30)       # centered moving-average windows (years)
BIN_WIDTH  <- 30                      # fixed bin width for comparison

daily_df <- read_csv("2_incremental/ncar_daily_q.csv", show_col_types = FALSE)
daily_df$date <- as.Date(daily_df$date)

wy_peaks <- read_csv("2_incremental/ncar_water_year_peaks.csv", show_col_types = FALSE)
wy_peaks$wy_start_date <- as.Date(wy_peaks$wy_start_date)
wy_peaks$wy_end_date   <- as.Date(wy_peaks$wy_end_date)
wy_peaks$peak_date     <- as.Date(wy_peaks$peak_date)

members <- unique(wy_peaks[, c("esm", "scenario", "cmip")])

# ============================================================================
# Step 1: derive baseflow_q_cfs, mean_q_cfs, anomaly per water year, per
# site, per member -- joining onto the existing peak/timing table from 06g
# ============================================================================

wy_metrics_list <- vector("list", nrow(members) * length(SITES))
k <- 0L

for (m in seq_len(nrow(members))) {
  esm_i <- members$esm[m]; scn_i <- members$scenario[m]; cmip_i <- members$cmip[m]
  
  sub_daily <- daily_df[daily_df$esm == esm_i & daily_df$scenario == scn_i, ]
  sub_daily <- sub_daily[order(sub_daily$date), ]
  
  sub_peaks <- wy_peaks[wy_peaks$esm == esm_i & wy_peaks$scenario == scn_i, ]
  
  for (site in SITES) {
    k <- k + 1L
    site_peaks <- sub_peaks[sub_peaks$site == site & !sub_peaks$flag_length, ]
    if (nrow(site_peaks) == 0) { wy_metrics_list[[k]] <- NULL; next }
    
    q_full <- sub_daily[[site]]
    d_full <- sub_daily$date
    
    baseflow <- numeric(nrow(site_peaks))
    meanq    <- numeric(nrow(site_peaks))
    
    for (i in seq_len(nrow(site_peaks))) {
      idx <- which(d_full > site_peaks$wy_start_date[i] &
                     d_full <= site_peaks$wy_end_date[i])
      baseflow[i] <- if (length(idx) == 0) NA_real_ else min(q_full[idx], na.rm = TRUE)
      meanq[i]    <- if (length(idx) == 0) NA_real_ else mean(q_full[idx], na.rm = TRUE)
    }
    
    chunk <- data.frame(
      site = site, esm = esm_i, scenario = scn_i, cmip = cmip_i,
      water_year = site_peaks$water_year,
      peak_q_cfs = site_peaks$peak_q_cfs,
      days_since_wy_start = site_peaks$days_since_wy_start,
      baseflow_q_cfs = baseflow,
      mean_q_cfs = meanq,
      stringsAsFactors = FALSE
    )
    chunk$anomaly <- (chunk$peak_q_cfs / chunk$baseflow_q_cfs)^(1/3)
    
    wy_metrics_list[[k]] <- chunk
  }
}

wy_metrics <- do.call(rbind, wy_metrics_list)
wy_metrics <- wy_metrics[order(wy_metrics$site, wy_metrics$scenario,
                               wy_metrics$esm, wy_metrics$water_year), ]
rownames(wy_metrics) <- NULL

write_csv(wy_metrics, "2_incremental/wy_metrics_all_sites.csv")

centered_ma <- function(x, yrs, window) {
  half <- floor(window / 2)
  out <- rep(NA_real_, length(x))
  ord <- order(yrs)
  x_o <- x[ord]; y_o <- yrs[ord]
  for (i in seq_along(x_o)) {
    lo <- y_o[i] - half
    hi <- y_o[i] + half
    w_idx <- which(y_o >= lo & y_o <= hi)
    if (length(w_idx) >= max(3, floor(window / 2))) {
      out[ord[i]] <- mean(x_o[w_idx], na.rm = TRUE)
    }
  }
  out
}

bin_mean <- function(x, yrs, bin_width, bin_origin) {
  bin_id <- floor((yrs - bin_origin) / bin_width)
  bin_center <- bin_origin + bin_id * bin_width + bin_width / 2
  ave(x, bin_id, FUN = function(v) mean(v, na.rm = TRUE))
}

# ============================================================================
# Step 3: apply smoothing per site x member x metric
# ============================================================================

METRICS <- c("peak_q_cfs", "days_since_wy_start", "baseflow_q_cfs",
             "mean_q_cfs", "anomaly")
BIN_ORIGIN <- min(wy_metrics$water_year, na.rm = TRUE)

for (w in MA_WINDOWS) {
  for (metric in METRICS) {
    new_col <- sprintf("%s_ma%d", metric, w)
    wy_metrics[[new_col]] <- NA_real_
  }
}
for (metric in METRICS) {
  wy_metrics[[sprintf("%s_bin%d", metric, BIN_WIDTH)]] <- NA_real_
}

for (site in SITES) {
  for (m in seq_len(nrow(members))) {
    esm_i <- members$esm[m]; scn_i <- members$scenario[m]
    idx <- which(wy_metrics$site == site & wy_metrics$esm == esm_i &
                   wy_metrics$scenario == scn_i)
    if (length(idx) < 5) next
    
    yrs <- wy_metrics$water_year[idx]
    
    for (metric in METRICS) {
      x <- wy_metrics[[metric]][idx]
      
      for (w in MA_WINDOWS) {
        wy_metrics[[sprintf("%s_ma%d", metric, w)]][idx] <- centered_ma(x, yrs, w)
      }
      wy_metrics[[sprintf("%s_bin%d", metric, BIN_WIDTH)]][idx] <-
        bin_mean(x, yrs, BIN_WIDTH, BIN_ORIGIN)
    }
  }
}

write_csv(wy_metrics, "2_incremental/wy_metrics_all_sites.csv")

# ============================================================================
# Step 4: comparison figure -- one representative member per site, all
# window lengths + bin mean overlaid against the raw annual values
# ============================================================================

REP_ESM <- "NorESM2-MM"
REP_SCN <- "ssp585"

ma_colors <- c("10" = "orange", "15" = "forestgreen",
               "20" = "blue", "30" = "purple")

SITE_GROUPS <- list(
  upper = c("CLALO", "CLADR"),
  lower = c("CLABE", "CLAPL")
)

dir.create("4_products", showWarnings = FALSE, recursive = TRUE)

for (grp_name in names(SITE_GROUPS)) {
  grp_sites <- SITE_GROUPS[[grp_name]]
  
  fig_file <- sprintf("4_products/diag_ma_window_sensitivity_%s.png", grp_name)
  png(fig_file, width = 2400, height = 3200, res = 200)
  par(mfrow = c(length(METRICS), length(grp_sites)), mar = c(4, 4, 3, 1))
  
  for (metric in METRICS) {
    for (site in grp_sites) {
      sub <- wy_metrics[wy_metrics$site == site & wy_metrics$esm == REP_ESM &
                          wy_metrics$scenario == REP_SCN, ]
      sub <- sub[order(sub$water_year), ]
      
      plot(sub$water_year, sub[[metric]], type = "p", pch = 16, cex = 0.5,
           col = "gray70",
           xlab = "Water year", ylab = metric,
           main = sprintf("%s: %s (raw + smoothed)", site, metric))
      
      for (w in MA_WINDOWS) {
        lines(sub$water_year, sub[[sprintf("%s_ma%d", metric, w)]],
              col = ma_colors[as.character(w)], lwd = 1.8)
      }
      lines(sub$water_year, sub[[sprintf("%s_bin%d", metric, BIN_WIDTH)]],
            col = "black", lwd = 2, lty = 2)
      
      if (metric == METRICS[1] && site == grp_sites[1]) {
        legend("topright",
               legend = c(paste0("MA", MA_WINDOWS), "30yr bin", "raw"),
               col = c(ma_colors[as.character(MA_WINDOWS)], "black", "gray70"),
               lty = c(rep(1, length(MA_WINDOWS)), 2, NA),
               pch = c(rep(NA, length(MA_WINDOWS)), NA, 16),
               lwd = 1.8, bty = "n", cex = 0.65)
      }
    }
  }
  
  dev.off()
}

# ============================================================================
# Scorecard
# ============================================================================

n_na_by_metric <- sapply(METRICS, function(mt) sum(is.na(wy_metrics[[mt]])))

scorecard_summary <- data.frame(
  metric = c(
    "n_members", "n_sites", "n_wy_rows",
    "water_year_min", "water_year_max",
    paste0("n_missing_", METRICS)
  ),
  value = c(
    nrow(members), length(SITES), nrow(wy_metrics),
    min(wy_metrics$water_year, na.rm = TRUE),
    max(wy_metrics$water_year, na.rm = TRUE),
    n_na_by_metric
  ),
  stringsAsFactors = FALSE
)

scorecard_summary