# ============================================================================
# 07a_ncar_daily_wide.R
# UCFR Filamentous Algae Project
# Stage 6: Produce wide-format daily streamflow from cached NCAR NetCDFs
#
# Inputs:
#   0_data/ncar_folders.csv                  (ESM-scenario inventory)
#   0_data/ncar_site_nc/*_site.nc            (cached downloads from 06e)
#
# Output:
#   2_incremental/ncar_daily_q.csv
#     Columns: esm, scenario, cmip, date, CLALO, CLADR, CLABE, CLAPL
#     One row per ESM-scenario-date. Reach columns are daily Q in cfs.
#     Sorted by scenario, esm, date.
#     Spin-up window (pre-1952) dropped at extraction — NCAR routed Q
#     during 1950-1951 is an unreliable model-initialization artifact.
#
# Reach-to-site mapping (for reference — not applied here):
#   CLALO -> DL, GR    CLADR -> BN    CLABE -> MS, BM, HU    CLAPL -> FH
#
# Run from project root. Requires 06e to have run first (cached .nc files).
# ============================================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(readr)
})

# ============================================================================
# Configuration
# ============================================================================

M3S_TO_CFS  <- 35.31467
SPINUP_CUTOFF <- as.Date("1952-01-01")  # drop pre-1952 NCAR spin-up artifact

CACHE_DIR  <- "0_data/ncar_site_nc"
OUT_DIR    <- "2_incremental"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 4 unique MERIT reaches in upstream-to-downstream order
REACHES <- data.frame(
  code     = c("CLALO",   "CLADR",   "CLABE",   "CLAPL"),
  reach_id = c(78017501,  78017483,  78017471,  78015509),
  stringsAsFactors = FALSE
)

# ============================================================================
# Read folder inventory, filter to successfully cached files
# ============================================================================

folders_df <- read_csv("0_data/ncar_folders.csv", show_col_types = FALSE)

folders_df$nc_path <- file.path(CACHE_DIR,
                                paste0(folders_df$folder,
                                       "_mizuRoute_daily_site.nc"))
folders_df$cached <- file.exists(folders_df$nc_path) &
  file.size(folders_df$nc_path) > 1e6

folders_ok <- folders_df[folders_df$cached, ]

# ============================================================================
# Extract daily Q for 4 reaches from each file, dropping pre-1952 spin-up
# ============================================================================

all_chunks   <- vector("list", nrow(folders_ok))
scorecard_rows <- vector("list", nrow(folders_ok))

for (i in seq_len(nrow(folders_ok))) {
  folder   <- folders_ok$folder[i]
  esm      <- folders_ok$esm[i]
  scenario <- folders_ok$scenario[i]
  cmip     <- folders_ok$cmip[i]
  path     <- folders_ok$nc_path[i]
  
  nc <- nc_open(path)
  
  # Time axis -> dates
  time_vals <- nc$dim$time$vals
  dates_full <- as.Date("1950-01-01") + time_vals
  
  # Drop NCAR spin-up window: pre-1952 values are unreliable model artifacts
  keep   <- dates_full >= SPINUP_CUTOFF
  dates  <- dates_full[keep]
  n_days <- length(dates)
  
  # Reach IDs (read once)
  rids <- ncvar_get(nc, "reachID")
  
  # Pre-allocate matrix: rows = days (post-cutoff), cols = 4 reaches
  q_mat <- matrix(NA_real_, nrow = n_days, ncol = nrow(REACHES))
  colnames(q_mat) <- REACHES$code
  
  n_missing_reach <- 0L
  
  for (r in seq_len(nrow(REACHES))) {
    seg_idx <- which(rids == REACHES$reach_id[r])
    if (length(seg_idx) == 0) {
      n_missing_reach <- n_missing_reach + 1L
      next
    }
    # streamflow[seg, time] -> start = c(seg_idx, 1), count = c(1, all)
    q_m3s <- ncvar_get(nc, "streamflow",
                       start = c(seg_idx, 1),
                       count = c(1, -1))
    q_mat[, r] <- as.numeric(q_m3s)[keep] * M3S_TO_CFS
  }
  
  nc_close(nc)
  
  # Build data.frame for this ESM-scenario
  chunk <- data.frame(
    esm      = esm,
    scenario = scenario,
    cmip     = cmip,
    date     = dates,
    q_mat,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  all_chunks[[i]] <- chunk
  
  scorecard_rows[[i]] <- data.frame(
    folder          = folder,
    esm             = esm,
    scenario        = scenario,
    cmip            = cmip,
    n_days_raw      = length(dates_full),
    n_days_kept     = n_days,
    date_min        = min(dates),
    date_max        = max(dates),
    n_missing_reach = n_missing_reach,
    stringsAsFactors = FALSE
  )
}

# ============================================================================
# Bind, sort, write
# ============================================================================

daily_df <- do.call(rbind, all_chunks)
daily_df <- daily_df[order(daily_df$scenario, daily_df$esm, daily_df$date), ]
rownames(daily_df) <- NULL

out_file <- file.path(OUT_DIR, "ncar_daily_q.csv")
write_csv(daily_df, out_file)

# ============================================================================
# Scorecard
# ============================================================================

scorecard <- do.call(rbind, scorecard_rows)
rownames(scorecard) <- NULL

mean_q_by_reach <- sapply(REACHES$code, function(rc) mean(daily_df[[rc]], na.rm = TRUE))

scorecard_summary <- data.frame(
  metric = c(
    "n_files_listed",
    "n_files_cached",
    "n_files_processed",
    "n_files_with_missing_reach",
    "total_rows_written",
    "date_range_min",
    "date_range_max",
    "n_unique_esm_scenarios",
    paste0("mean_Q_cfs_", REACHES$code)
  ),
  value = c(
    nrow(folders_df),
    sum(folders_df$cached),
    nrow(folders_ok),
    sum(scorecard$n_missing_reach > 0),
    nrow(daily_df),
    as.character(min(daily_df$date)),
    as.character(max(daily_df$date)),
    length(unique(paste(daily_df$esm, daily_df$scenario))),
    sprintf("%.1f", mean_q_by_reach)
  ),
  stringsAsFactors = FALSE
)

scorecard
scorecard_summary