# ============================================================================
# 02_ingest_usgs.R
# UCFR Filamentous Algae Project
# Stage 1: Ingest USGS daily discharge data from NWIS
#
# Sources:  USGS Water Data API via dataRetrieval::read_waterdata_daily()
#           Seven gauges covering six study sites plus Huron composite
#
# Gauge list:
#   12324200  Clark Fork at Deer Lodge MT                        -> DL
#   12324400  Clark Fork above Little Blackfoot River MT         -> GR
#   12331800  Clark Fork at Bonita MT                            -> BN
#   12340500  Clark Fork above Missoula MT                       -> MS
#   12352500  Bitterroot River near Missoula MT                  -> HU (component)
#   12353000  Clark Fork below Missoula MT                       -> BM + HU (component)
#   12354500  Clark Fork at St. Regis MT                         -> FH
#
# HU composite discharge = 12352500 + 12353000
#
# Output:   0_data/usgs_daily_q_raw.csv
#
# Notes:
#   - Pulls full period of record on every run
#   - Gauges pulled one at a time in a loop to avoid server timeouts
#   - Uses read_waterdata_daily() — replaces deprecated readNWISdv()
#   - HU composite is NOT calculated here — handled in later pipeline stage
#   - Observation-specific Q extraction handled at join stage
# ============================================================================

library(dataRetrieval)
library(readr)

# ----------------------------------------------------------------------------
# 1. Define output path
# ----------------------------------------------------------------------------

out_dir  <- "0_data"
out_file <- file.path(out_dir, "usgs_daily_q_raw.csv")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ----------------------------------------------------------------------------
# 2. Define gauge list
# ----------------------------------------------------------------------------

# read_waterdata_daily requires "USGS-" prefix on site IDs
gauges <- c(
  "USGS-12324200",  # Clark Fork at Deer Lodge MT                -> DL
  "USGS-12324400",  # Clark Fork above Little Blackfoot River MT -> GR
  "USGS-12331800",  # Clark Fork at Bonita MT                    -> BN
  "USGS-12340500",  # Clark Fork above Missoula MT               -> MS
  "USGS-12352500",  # Bitterroot River near Missoula MT          -> HU component
  "USGS-12353000",  # Clark Fork below Missoula MT               -> BM + HU component
  "USGS-12354500"   # Clark Fork at St. Regis MT                 -> FH
)

param_cd <- "00060"  # discharge, cubic feet per second
stat_id  <- "00003"  # mean daily value

# ----------------------------------------------------------------------------
# 3. Pull full period of record gauge by gauge
# ----------------------------------------------------------------------------
# Pulling all gauges in one call times out — loop with a pause between requests

cat("Querying USGS Water Data API...\n")
cat("Parameter: 00060 (discharge, cfs) | Statistic: 00003 (mean daily)\n\n")

q_list <- vector("list", length(gauges))

for (i in seq_along(gauges)) {
  g <- gauges[i]
  cat(sprintf("  Pulling %s (%d of %d)...", g, i, length(gauges)))
  
  result <- tryCatch({
    read_waterdata_daily(
      monitoring_location_id = g,
      parameter_code         = param_cd,
      statistic_id           = stat_id,
      time                   = c(NA, NA)  # full period of record
    )
  }, error = function(e) {
    cat(" FAILED:", conditionMessage(e), "\n")
    NULL
  })
  
  if (!is.null(result) && nrow(result) > 0) {
    cat(sprintf(" %d records\n", nrow(result)))
    q_list[[i]] <- result
  } else {
    cat(" no data returned\n")
  }
  
  Sys.sleep(2)  # pause between requests to avoid rate limiting
}

# ----------------------------------------------------------------------------
# 4. Combine and standardize columns
# ----------------------------------------------------------------------------

q_raw <- do.call(rbind, q_list)

# Drop sf geometry — read_waterdata_daily returns a spatial dataframe
# which causes tibble column renaming to fail downstream
q_raw <- as.data.frame(q_raw)
q_raw$geometry <- NULL

cat("\nRaw combined dimensions:", nrow(q_raw), "rows x", ncol(q_raw), "columns\n")
cat("Column names:", paste(names(q_raw), collapse = ", "), "\n\n")

# Note: USGS-12354700 (Clark Fork above Flathead / original FH gauge) is down.
# Replaced with USGS-12354500 (Clark Fork at St. Regis MT) as FH proxy.

cols_expected <- c("monitoring_location_id", "time", "value", "approval_status")
missing_cols  <- setdiff(cols_expected, names(q_raw))

if (length(missing_cols) > 0) {
  stop("Expected columns not found: ",
       paste(missing_cols, collapse = ", "),
       "\nActual columns: ", paste(names(q_raw), collapse = ", "))
}

q_slim <- q_raw[ , cols_expected]

names(q_slim) <- c("site_no", "Date", "Q_cfs", "Q_code")

# Strip "USGS-" prefix from site_no to match gauge numbers used elsewhere
q_slim$site_no <- sub("USGS-", "", q_slim$site_no)

# ----------------------------------------------------------------------------
# 5. Parse date and derive year / month / day
# ----------------------------------------------------------------------------

q_slim$Date  <- as.Date(q_slim$Date)
q_slim$Year  <- as.integer(format(q_slim$Date, "%Y"))
q_slim$Month <- as.integer(format(q_slim$Date, "%m"))
q_slim$Day   <- as.integer(format(q_slim$Date, "%d"))

# Ensure Q_cfs is numeric
q_slim$Q_cfs <- as.numeric(q_slim$Q_cfs)

# ----------------------------------------------------------------------------
# 6. Ingestion diagnostics
# ----------------------------------------------------------------------------

cat("\n--- Ingestion Summary by Gauge ---\n")

for (g in sub("USGS-", "", gauges)) {
  d <- q_slim[q_slim$site_no == g, ]
  cat(sprintf(
    "  %s  n=%-6d  %s to %s  missing Q: %d\n",
    g,
    nrow(d),
    as.character(min(d$Date, na.rm = TRUE)),
    as.character(max(d$Date, na.rm = TRUE)),
    sum(is.na(d$Q_cfs))
  ))
}

cat("\nTotal records:          ", nrow(q_slim), "\n")
cat("Total missing Q values: ", sum(is.na(q_slim$Q_cfs)), "\n")

# ----------------------------------------------------------------------------
# 7. Write to 0_data
# ----------------------------------------------------------------------------

write_csv(q_slim, out_file)
cat("\nSaved to:", out_file, "\n")
cat("Done.\n")