# ============================================================================
# 06_final_join.R
# UCFR Filamentous Algae Project
# Stage 5: Finalization, TP transform, QC, and clean output
#
# Input:    2_incremental/usgs_processed.csv
# Output:   2_incremental/ucfr_model_ready.csv
#           2_incremental/qc_flags.csv (only if any flags raised)
#
# Steps:
#   1. Read assembled dataset
#   2. Select and order final modeling columns
#   3. Derive corrected logTP_mg_L (replaces any pre-existing column of the
#      same name) under an explicit TP_SCALE toggle
#   4. Completeness report by site and year (captured into scorecard)
#   5. QC checks (flag, never silently drop)
#   6. Write final model-ready dataset
#   7. Print end-of-script scorecard (no mid-script cat())
#
# Notes:
#   - The join between VNRP and USGS was completed in 05_process_usgs.R
#   - This script is for finalization, column ordering, TP/Q transforms, and QC
#   - TP fix (D3, 2026-06): logTP_mg_L previously stored log10(1+TP_mg_L),
#     which is ~linear over the observed TP range and barely variance-
#     stabilizing. Corrected default is log10(TP_mg_L * 1000), i.e. log10 of
#     concentration in ug/L (Suplee 24 ug/L threshold -> log10(24) = 1.38).
#     TP_SCALE toggle below preserves the old transform for comparison only.
#     This replaces logTP_mg_L in place -- any script reading
#     ucfr_model_ready.csv picks up the corrected column automatically.
#   - logQ_obs_cfs fix (2026-06): this column was referenced by
#     08_temporal_validation.R (and presumably the M1 bloom model formula)
#     but was never derived anywhere in the coded pipeline -- it only ever
#     existed as a manually-added column outside the pipeline. Now derived
#     here as plain log10(Q_obs_cfs); no log1p/offset needed since observed
#     discharge here is always well above zero (hundreds-to-thousands cfs).
# ============================================================================

library(readr)
library(dplyr)

# ----------------------------------------------------------------------------
# 0. Configuration
# ----------------------------------------------------------------------------

TP_SCALE <- "log10_ugL"   # "log10_ugL" (recommended/default) or "log1p_mgL"

in_file  <- "2_incremental/usgs_processed.csv"
out_dir  <- "2_incremental"
out_file <- file.path(out_dir, "ucfr_model_ready.csv")

# ----------------------------------------------------------------------------
# 1. Read assembled dataset
# ----------------------------------------------------------------------------

dat <- read_csv(in_file, show_col_types = FALSE)

# ----------------------------------------------------------------------------
# 2. Select and order final modeling columns
# ----------------------------------------------------------------------------

final_cols <- c(
  # Identifiers
  "Site", "Year", "Month", "date_yearmon",
  
  # Response variables
  "CHLa", "logCHLa", "AFDM",
  
  # Hydrological predictors
  "anomaly", "Days_Since_Freshet",
  "Q_peak_cfs", "Q_baseflow_cfs", "Q_obs_cfs",
  
  # Nutrient predictors
  "TP_mg_L", "TN_mg_L", "SRP_mg_L",
  "NH4_mg_L", "NO3_mg_L", "DIN_mg_L",
  
  # Physical / water quality predictors
  "pH", "Temp_oC", "SPC", "TDS", "TURBIDITY"
)

missing_cols <- setdiff(final_cols, names(dat))
if (length(missing_cols) > 0) {
  warning("Expected columns not found -- will be omitted:\n  ",
          paste(missing_cols, collapse = ", "))
  final_cols <- intersect(final_cols, names(dat))
}

dat_final <- dat[ , final_cols]

# ----------------------------------------------------------------------------
# 3. Derive corrected logTP_mg_L (replaces in place; see header note)
# ----------------------------------------------------------------------------

if (!"TP_mg_L" %in% names(dat_final)) {
  stop("TP_mg_L not present in assembled dataset -- cannot derive logTP_mg_L.")
}

if (TP_SCALE == "log10_ugL") {
  dat_final$logTP_mg_L <- log10(dat_final$TP_mg_L * 1000)
} else if (TP_SCALE == "log1p_mgL") {
  dat_final$logTP_mg_L <- log10(1 + dat_final$TP_mg_L)
} else {
  stop("TP_SCALE must be 'log10_ugL' or 'log1p_mgL'")
}

# Ensure logTP_mg_L is positioned with the other nutrient predictors
if ("TN_mg_L" %in% names(dat_final)) {
  dat_final <- dat_final %>%
    relocate(logTP_mg_L, .after = TP_mg_L)
}

# ----------------------------------------------------------------------------
# 3b. Derive logQ_obs_cfs (was never coded anywhere upstream -- see header)
# ----------------------------------------------------------------------------

if (!"Q_obs_cfs" %in% names(dat_final)) {
  stop("Q_obs_cfs not present in assembled dataset -- cannot derive logQ_obs_cfs.")
}

n_nonpos_q <- sum(!is.na(dat_final$Q_obs_cfs) & dat_final$Q_obs_cfs <= 0)
if (n_nonpos_q > 0) {
  warning(sprintf(
    "%d rows have Q_obs_cfs <= 0 -- log10 will produce NaN/Inf for these.",
    n_nonpos_q))
}

dat_final$logQ_obs_cfs <- log10(dat_final$Q_obs_cfs)

dat_final <- dat_final %>%
  relocate(logQ_obs_cfs, .after = Q_obs_cfs)

# ----------------------------------------------------------------------------
# 4. Completeness report by site and year (captured, not printed mid-script)
# ----------------------------------------------------------------------------

key_vars <- c("CHLa", "TP_mg_L", "logTP_mg_L", "TN_mg_L", "anomaly",
              "Days_Since_Freshet", "Q_obs_cfs", "logQ_obs_cfs")
key_vars <- intersect(key_vars, names(dat_final))

completeness_by_site <- dat_final %>%
  group_by(Site) %>%
  summarise(
    n         = n(),
    year_min  = min(Year, na.rm = TRUE),
    year_max  = max(Year, na.rm = TRUE),
    across(all_of(key_vars), ~ sum(!is.na(.x)), .names = "n_{.col}"),
    .groups   = "drop"
  )

all_vars <- c("CHLa", "logCHLa", "AFDM", "TP_mg_L", "logTP_mg_L", "TN_mg_L",
              "SRP_mg_L", "NH4_mg_L", "NO3_mg_L", "DIN_mg_L", "pH", "Temp_oC",
              "SPC", "TDS", "TURBIDITY", "anomaly", "Days_Since_Freshet",
              "Q_peak_cfs", "Q_baseflow_cfs", "Q_obs_cfs", "logQ_obs_cfs")
all_vars <- intersect(all_vars, names(dat_final))

completeness_overall <- data.frame(
  Variable     = all_vars,
  n_obs        = sapply(all_vars, function(v) sum(!is.na(dat_final[[v]]))),
  n_missing    = sapply(all_vars, function(v) sum(is.na(dat_final[[v]]))),
  pct_complete = sapply(all_vars, function(v)
    round(100 * sum(!is.na(dat_final[[v]])) / nrow(dat_final))),
  row.names    = NULL,
  stringsAsFactors = FALSE
)

# ----------------------------------------------------------------------------
# 5. QC CHECKS -- flag, never silently drop
# ----------------------------------------------------------------------------

qc_flags <- data.frame(
  Site = character(), Year = integer(), Month = integer(),
  Flag = character(), stringsAsFactors = FALSE
)

add_flags <- function(qc_flags, flagged, label) {
  if (nrow(flagged) == 0) return(qc_flags)
  rbind(qc_flags, data.frame(
    Site = flagged$Site, Year = flagged$Year, Month = flagged$Month,
    Flag = label, stringsAsFactors = FALSE
  ))
}

# QC1: negative Days_Since_Freshet
if ("Days_Since_Freshet" %in% names(dat_final)) {
  neg_dsf <- dat_final[!is.na(dat_final$Days_Since_Freshet) &
                         dat_final$Days_Since_Freshet < 0, ]
  qc_flags <- add_flags(qc_flags, neg_dsf, "Negative_Days_Since_Freshet")
}

# QC2: logCHLa outliers (> 3 SD from mean)
if ("logCHLa" %in% names(dat_final)) {
  chla_mean <- mean(dat_final$logCHLa, na.rm = TRUE)
  chla_sd   <- sd(dat_final$logCHLa, na.rm = TRUE)
  chla_out  <- dat_final[!is.na(dat_final$logCHLa) &
                           abs(dat_final$logCHLa - chla_mean) > 3 * chla_sd, ]
  qc_flags <- add_flags(qc_flags, chla_out, "CHLa_outlier_3SD")
}

# QC3: implausible TP values (raw mg/L > 1 -- likely unit error)
if ("TP_mg_L" %in% names(dat_final)) {
  tp_high <- dat_final[!is.na(dat_final$TP_mg_L) & dat_final$TP_mg_L > 1, ]
  qc_flags <- add_flags(qc_flags, tp_high, "TP_gt_1mgL")
}

# QC4: anomaly out of plausible range (<= 0)
if ("anomaly" %in% names(dat_final)) {
  anom_bad <- dat_final[!is.na(dat_final$anomaly) & dat_final$anomaly <= 0, ]
  qc_flags <- add_flags(qc_flags, anom_bad, "Anomaly_lte_zero")
}

# -- ADD FUTURE QC CHECKS BELOW THIS LINE ------------------------------------
# Template:
# flagged <- dat_final[condition, ]
# qc_flags <- add_flags(qc_flags, flagged, "FLAG_LABEL")
# -----------------------------------------------------------------------------

qc_log_path <- NA_character_
if (nrow(qc_flags) > 0) {
  qc_log_path <- file.path(out_dir, "qc_flags.csv")
  write_csv(qc_flags, qc_log_path)
}

# ----------------------------------------------------------------------------
# 6. Write final model-ready dataset
# ----------------------------------------------------------------------------

write_csv(dat_final, out_file)

# ----------------------------------------------------------------------------
# 7. SCORECARD
# ----------------------------------------------------------------------------

scorecard_summary <- data.frame(
  metric = c(
    "input_file",
    "output_file",
    "n_rows",
    "n_cols",
    "n_sites",
    "year_min",
    "year_max",
    "TP_SCALE",
    "logTP_mg_L_formula",
    "logQ_obs_cfs_formula",
    "n_Q_obs_cfs_nonpositive",
    "n_qc_flags_total",
    "qc_flag_log_path"
  ),
  value = c(
    in_file,
    out_file,
    nrow(dat_final),
    ncol(dat_final),
    length(unique(dat_final$Site)),
    min(dat_final$Year, na.rm = TRUE),
    max(dat_final$Year, na.rm = TRUE),
    TP_SCALE,
    ifelse(TP_SCALE == "log10_ugL", "log10(TP_mg_L * 1000)", "log10(1 + TP_mg_L)"),
    "log10(Q_obs_cfs)",
    n_nonpos_q,
    nrow(qc_flags),
    ifelse(is.na(qc_log_path), "none", qc_log_path)
  ),
  stringsAsFactors = FALSE
)

scorecard_summary
completeness_by_site
completeness_overall
if (nrow(qc_flags) > 0) qc_flags