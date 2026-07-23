# ============================================================================
# 06g_ncar_water_year_peaks.R
# UCFR Filamentous Algae Project
# Stage 6: Water-year segmentation and peak discharge timing/magnitude
#          across the NCAR ESM-scenario ensemble
#
# Purpose:
#   Define water years by the low-flow date in a Sep15-Nov15 search window
#   (rather than a fixed Oct 1 boundary), then extract the date and
#   magnitude of peak discharge within each bounded water year. This lets
#   us look directly at how freshet timing and peak magnitude drift over
#   1952-2099 across all 26 ESM-scenario members and all 4 NCAR reaches,
#   before committing to a baseline-window choice for delta computation.
#
# Input:
#   2_incremental/ncar_daily_q.csv   (from 06f; spin-up already trimmed)
#
# Output:
#   2_incremental/ncar_water_year_peaks.csv
#     One row per site x esm x scenario x water_year:
#       wy_start_date, wy_end_date, peak_date, peak_q_cfs,
#       days_since_wy_start, wy_length_days
#   4_products/fig_freshet_timing_<SITE>.png   (one per site, 4 total)
#
# Run from project root. Requires 06f to have run first.
# ============================================================================

suppressPackageStartupMessages({
  library(readr)
})

# ============================================================================
# Configuration
# ============================================================================

IN_FILE   <- "2_incremental/ncar_daily_q.csv"
OUT_DIR   <- "2_incremental"
FIG_DIR   <- "4_products"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

SITES <- c("CLALO", "CLADR", "CLABE", "CLAPL")

# Low-flow search window (month-day, applied within each calendar year)
LOWFLOW_WINDOW_START <- c(month = 9,  day = 15)
LOWFLOW_WINDOW_END   <- c(month = 11, day = 15)

# Sanity bound on water-year length (days) - flags degenerate boundaries
WY_LENGTH_MIN <- 300
WY_LENGTH_MAX <- 430

# ============================================================================
# Read data
# ============================================================================

daily_df <- read_csv(IN_FILE, show_col_types = FALSE)
daily_df$date <- as.Date(daily_df$date)

members <- unique(daily_df[, c("esm", "scenario", "cmip")])

# ============================================================================
# Helper: find low-flow boundary dates for one site's series within one
# member, by searching every Sep15-Nov15 window across the years spanned
# ============================================================================

find_wy_boundaries <- function(dates, q) {
  yrs <- unique(format(dates, "%Y"))
  yrs <- as.integer(yrs)
  
  boundaries <- as.Date(character(0))
  
  for (yr in yrs) {
    win_start <- as.Date(sprintf("%d-%02d-%02d", yr,
                                 LOWFLOW_WINDOW_START["month"],
                                 LOWFLOW_WINDOW_START["day"]))
    win_end   <- as.Date(sprintf("%d-%02d-%02d", yr,
                                 LOWFLOW_WINDOW_END["month"],
                                 LOWFLOW_WINDOW_END["day"]))
    
    idx <- which(dates >= win_start & dates <= win_end)
    if (length(idx) == 0) next
    if (all(is.na(q[idx]))) next
    
    min_idx <- idx[which.min(q[idx])]
    boundaries <- c(boundaries, dates[min_idx])
  }
  
  sort(unique(boundaries))
}

# ============================================================================
# Helper: given boundaries, extract peak date/magnitude per bounded
# water year, with a length sanity check
# ============================================================================

extract_wy_peaks <- function(dates, q, boundaries) {
  n_wy <- length(boundaries) - 1
  if (n_wy < 1) {
    return(data.frame(
      wy_start_date = as.Date(character(0)),
      wy_end_date   = as.Date(character(0)),
      wy_length_days = integer(0),
      peak_date     = as.Date(character(0)),
      peak_q_cfs    = numeric(0),
      days_since_wy_start = integer(0),
      flag_length   = logical(0)
    ))
  }
  
  out <- vector("list", n_wy)
  
  for (i in seq_len(n_wy)) {
    wy_start <- boundaries[i]
    wy_end   <- boundaries[i + 1]
    wy_len   <- as.integer(wy_end - wy_start)
    
    idx <- which(dates > wy_start & dates <= wy_end)
    if (length(idx) == 0 || all(is.na(q[idx]))) {
      out[[i]] <- data.frame(
        wy_start_date = wy_start, wy_end_date = wy_end,
        wy_length_days = wy_len, peak_date = as.Date(NA),
        peak_q_cfs = NA_real_, days_since_wy_start = NA_integer_,
        flag_length = (wy_len < WY_LENGTH_MIN | wy_len > WY_LENGTH_MAX)
      )
      next
    }
    
    peak_idx  <- idx[which.max(q[idx])]
    peak_date <- dates[peak_idx]
    
    out[[i]] <- data.frame(
      wy_start_date = wy_start,
      wy_end_date   = wy_end,
      wy_length_days = wy_len,
      peak_date     = peak_date,
      peak_q_cfs    = q[peak_idx],
      days_since_wy_start = as.integer(peak_date - wy_start),
      flag_length   = (wy_len < WY_LENGTH_MIN | wy_len > WY_LENGTH_MAX)
    )
  }
  
  do.call(rbind, out)
}

# ============================================================================
# Main loop: site x member
# ============================================================================

all_peaks <- vector("list", nrow(members) * length(SITES))
k <- 0L

for (m in seq_len(nrow(members))) {
  esm_i      <- members$esm[m]
  scenario_i <- members$scenario[m]
  cmip_i     <- members$cmip[m]
  
  sub <- daily_df[daily_df$esm == esm_i & daily_df$scenario == scenario_i, ]
  sub <- sub[order(sub$date), ]
  
  for (site in SITES) {
    k <- k + 1L
    q <- sub[[site]]
    dates <- sub$date
    
    boundaries <- find_wy_boundaries(dates, q)
    peaks <- extract_wy_peaks(dates, q, boundaries)
    
    if (nrow(peaks) > 0) {
      peaks$site     <- site
      peaks$esm      <- esm_i
      peaks$scenario <- scenario_i
      peaks$cmip     <- cmip_i
      peaks$water_year <- as.integer(format(peaks$wy_start_date, "%Y"))
    }
    
    all_peaks[[k]] <- peaks
  }
}

wy_peaks <- do.call(rbind, all_peaks)
wy_peaks <- wy_peaks[, c("site", "esm", "scenario", "cmip", "water_year",
                         "wy_start_date", "wy_end_date", "wy_length_days",
                         "peak_date", "peak_q_cfs", "days_since_wy_start",
                         "flag_length")]
wy_peaks <- wy_peaks[order(wy_peaks$site, wy_peaks$scenario, wy_peaks$esm,
                           wy_peaks$water_year), ]
rownames(wy_peaks) <- NULL

out_file <- file.path(OUT_DIR, "ncar_water_year_peaks.csv")
write_csv(wy_peaks, out_file)

# ============================================================================
# Figures: one per site, all members overlaid, peak timing vs water year
# ============================================================================

member_colors <- setNames(
  rainbow(nrow(members), s = 0.6, v = 0.8),
  paste(members$esm, members$scenario)
)

for (site in SITES) {
  site_dat <- wy_peaks[wy_peaks$site == site & !wy_peaks$flag_length, ]
  
  fig_file <- file.path(FIG_DIR, sprintf("fig_freshet_timing_%s.png", site))
  png(fig_file, width = 2000, height = 1400, res = 200)
  
  plot(NA, NA,
       xlim = range(site_dat$water_year, na.rm = TRUE),
       ylim = range(site_dat$days_since_wy_start, na.rm = TRUE),
       xlab = "Water year", ylab = "Days since water-year start",
       main = sprintf("%s: peak discharge timing, 1952-2099 (all ESM members)", site))
  
  for (mem in unique(paste(site_dat$esm, site_dat$scenario))) {
    md <- site_dat[paste(site_dat$esm, site_dat$scenario) == mem, ]
    md <- md[order(md$water_year), ]
    lines(md$water_year, md$days_since_wy_start,
          col = member_colors[mem], lwd = 1)
  }
  
  dev.off()
}

# ============================================================================
# Scorecard
# ============================================================================

n_flagged <- sum(wy_peaks$flag_length, na.rm = TRUE)
n_na_peak <- sum(is.na(wy_peaks$peak_q_cfs))

scorecard_summary <- data.frame(
  metric = c(
    "n_members",
    "n_sites",
    "n_water_year_rows",
    "n_flagged_wy_length",
    "n_missing_peak",
    "water_year_range_min",
    "water_year_range_max",
    paste0("median_days_since_start_", SITES)
  ),
  value = c(
    nrow(members),
    length(SITES),
    nrow(wy_peaks),
    n_flagged,
    n_na_peak,
    min(wy_peaks$water_year, na.rm = TRUE),
    max(wy_peaks$water_year, na.rm = TRUE),
    sapply(SITES, function(s) {
      median(wy_peaks$days_since_wy_start[wy_peaks$site == s], na.rm = TRUE)
    })
  ),
  stringsAsFactors = FALSE
)

scorecard_summary