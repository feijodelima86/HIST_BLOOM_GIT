# ============================================================================
# 04_lag_selection.R
# UCFR Filamentous Algae Project
# Stage 3: Temporal lag selection for chemistry predictors
#
# Input:    2_incremental/vnrp_processed.csv
# Outputs:  2_incremental/lag_selection.csv
#           4_products/diagnostics/lag_scatterplots.pdf
#
# Steps:
#   1. Read processed VNRP data
#   2. Split biological and chemistry streams
#   3. For each variable x month lag, merge onto biological stream by Site+Year
#   4. Compute correlations of CHLa vs each variable x lag combination
#   5. Produce scatterplots (one page per variable, four lags side by side)
#   6. Output correlation summary and lag selection table
#
# Notes:
#   - CHLa is log10-transformed for correlations and plots (response variable)
#   - Months 6-9 (June-September) tested as lags
#   - Lag selection is documented explicitly for pipeline reproducibility
#   - Backcheck against existing n=221 dataset to be done manually by user
#   - D1 fix (2026-06): vnrp_processed.csv now contains visit-level bio rows
#     (some Site-Year-Months have 2 real bio visits, each its own row, with
#     that month's chemistry broadcast onto both -- see 03_process_vnrp.R).
#     Section 3's chem_month lookup must therefore be de-duplicated to one
#     row per Site-Year before the Site+Year join below, otherwise a
#     double-visit month produces two IDENTICAL chemistry rows and the join
#     fans out many-to-many. De-duplication is a safe no-op: duplicate rows
#     for the same Site-Year-Month carry identical broadcast chemistry
#     values, so keeping the first occurrence loses no information.
#     Section 7's get_chem_month()/mean() lookup is unaffected -- mean() of
#     duplicate identical values is already a no-op, confirmed separately.
# ============================================================================

library(readr)
library(dplyr)
library(lubridate)

# ----------------------------------------------------------------------------
# 1. Read processed VNRP data
# ----------------------------------------------------------------------------

in_file      <- "2_incremental/vnrp_processed.csv"
out_dir_inc  <- "2_incremental"
out_dir_diag <- "4_products/diagnostics"

for (d in c(out_dir_inc, out_dir_diag)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

cat("Reading processed VNRP data...\n")
vnrp <- read_csv(in_file, show_col_types = FALSE)
cat("Rows:", nrow(vnrp), "\n\n")

# ----------------------------------------------------------------------------
# 2. Split biological and chemistry streams
# ----------------------------------------------------------------------------

# Biological stream: rows where CHLa is not NA
# Log10-transform CHLa — this is the response variable in the model
bio <- vnrp[!is.na(vnrp$CHLa), ]
bio$logCHLa <- log10(bio$CHLa)

cat("Biological stream (CHLa non-NA):", nrow(bio), "rows\n")
cat("Year range:", min(bio$Year), "to", max(bio$Year), "\n")
cat("Sites:", paste(sort(unique(bio$Site)), collapse = ", "), "\n\n")

# Chemistry stream: all rows, split by month
chem_vars <- c("TP_mg_L", "TN_mg_L", "SRP_mg_L", "NH4_mg_L",
               "NO3_mg_L", "DIN_mg_L", "pH", "Temp_oC",
               "SPC", "TDS", "TURBIDITY")

lag_months      <- c(6, 7, 8, 9)
lag_month_names <- c("June", "July", "August", "September")

# ----------------------------------------------------------------------------
# 3. Build merged dataset for each variable x lag combination
# ----------------------------------------------------------------------------

# For each lag month, extract chemistry and merge onto bio by Site + Year
# Suffix each chemistry column with the month name to avoid collisions

cat("Building lag-merged datasets...\n")

# Start with the biological stream
merged <- bio[ , c("Site", "Year", "Month", "logCHLa", "CHLa")]

for (i in seq_along(lag_months)) {
  m      <- lag_months[i]
  m_name <- lag_month_names[i]
  
  chem_month <- vnrp[vnrp$Month == m, c("Site", "Year", chem_vars)]
  
  # D1 fix: de-duplicate to one row per Site-Year before the join. A
  # double-visit month now produces 2 identical chemistry rows here (same
  # month-level value broadcast to both bio visits in 03) -- without this,
  # the Site+Year join below fans out many-to-many and silently duplicates
  # merged's rows. Safe: duplicate rows carry identical chemistry values,
  # so keeping the first occurrence drops no information.
  chem_month <- chem_month[!duplicated(chem_month[c("Site", "Year")]), ]
  
  # Rename chemistry columns to include lag month
  names(chem_month)[names(chem_month) %in% chem_vars] <-
    paste0(chem_vars, "_", m_name)
  
  merged <- left_join(merged, chem_month, by = c("Site", "Year"))
}

cat("Merged dimensions:", nrow(merged), "rows x", ncol(merged), "cols\n\n")

# ----------------------------------------------------------------------------
# 4. Compute correlations: logCHLa vs each variable x lag
# ----------------------------------------------------------------------------

cat("Computing correlations...\n")

cor_results <- data.frame(
  Variable  = character(),
  Lag_Month = character(),
  r         = numeric(),
  n         = integer(),
  stringsAsFactors = FALSE
)

for (v in chem_vars) {
  for (m_name in lag_month_names) {
    col_name <- paste0(v, "_", m_name)
    if (col_name %in% names(merged)) {
      valid <- merged[!is.na(merged[[col_name]]) & !is.na(merged$logCHLa), ]
      if (nrow(valid) > 5) {
        r <- cor(valid$logCHLa, valid[[col_name]], use = "complete.obs")
      } else {
        r <- NA
      }
      cor_results <- rbind(cor_results, data.frame(
        Variable  = v,
        Lag_Month = m_name,
        r         = round(r, 3),
        n         = nrow(valid),
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Print correlation table
cat("\n--- Correlation Table: log10(CHLa) vs Chemistry by Lag Month ---\n")
cat(sprintf("  %-15s  %-10s  %6s  %5s\n", "Variable", "Lag", "r", "n"))
cat(paste(rep("-", 45), collapse = ""), "\n")
for (i in seq_len(nrow(cor_results))) {
  cat(sprintf("  %-15s  %-10s  %6.3f  %5d\n",
              cor_results$Variable[i],
              cor_results$Lag_Month[i],
              cor_results$r[i],
              cor_results$n[i]))
}

# ----------------------------------------------------------------------------
# 5. Identify best lag per variable (highest absolute correlation)
# ----------------------------------------------------------------------------

best_lag <- cor_results %>%
  group_by(Variable) %>%
  filter(!is.na(r)) %>%
  slice_max(abs(r), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(Best_Lag = Lag_Month, Best_r = r, Best_n = n)

cat("\n--- Best Lag per Variable ---\n")
cat(sprintf("  %-15s  %-10s  %6s  %5s\n", "Variable", "Best Lag", "r", "n"))
cat(paste(rep("-", 45), collapse = ""), "\n")
for (i in seq_len(nrow(best_lag))) {
  cat(sprintf("  %-15s  %-10s  %6.3f  %5d\n",
              best_lag$Variable[i],
              best_lag$Best_Lag[i],
              best_lag$Best_r[i],
              best_lag$Best_n[i]))
}

# ----------------------------------------------------------------------------
# 6. Scatterplots: one page per variable, four lags side by side
# ----------------------------------------------------------------------------

cat("\nGenerating scatterplots...\n")

pdf(file.path(out_dir_diag, "lag_scatterplots.pdf"),
    width = 11, height = 8)

for (v in chem_vars) {
  
  # Determine x-axis range across all lags for consistent scaling
  all_vals <- unlist(lapply(lag_month_names, function(m_name) {
    col <- paste0(v, "_", m_name)
    if (col %in% names(merged)) merged[[col]] else NULL
  }))
  x_range <- range(all_vals, na.rm = TRUE)
  y_range <- range(merged$logCHLa, na.rm = TRUE)
  
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
  
  for (i in seq_along(lag_month_names)) {
    m_name   <- lag_month_names[i]
    col_name <- paste0(v, "_", m_name)
    
    if (!col_name %in% names(merged)) next
    
    valid <- merged[!is.na(merged[[col_name]]) & !is.na(merged$logCHLa), ]
    r_val <- cor_results$r[cor_results$Variable == v &
                             cor_results$Lag_Month == m_name]
    n_val <- cor_results$n[cor_results$Variable == v &
                             cor_results$Lag_Month == m_name]
    
    # Color points by site
    site_cols <- c(DL = "#E41A1C", GR = "#377EB8", BN = "#4DAF4A",
                   MS = "#984EA3", BM = "#FF7F00", HU = "#A65628",
                   FH = "#F781BF")
    pt_col <- site_cols[valid$Site]
    
    plot(valid[[col_name]], valid$logCHLa,
         xlim  = x_range,
         ylim  = y_range,
         xlab  = v,
         ylab  = "log10(CHLa)",
         main  = sprintf("%s lag  |  r = %.3f  |  n = %d", m_name, r_val, n_val),
         pch   = 16,
         col   = pt_col,
         cex   = 0.8)
    
    # Add regression line if enough points
    if (nrow(valid) > 5) {
      abline(lm(logCHLa ~ valid[[col_name]], data = valid),
             col = "grey40", lty = 2)
    }
    
    # Add Suplee TP threshold line for TP variable
    if (v == "TP_mg_L") {
      abline(v = 0.024, col = "red", lty = 3, lwd = 1.5)
    }
  }
  
  # Add site legend and overall title
  mtext(paste("log10(CHLa) vs", v, "— Lag Month Comparison"),
        outer = TRUE, cex = 1.1, font = 2)
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
      mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottomright",
         legend = names(site_cols),
         col    = site_cols,
         pch    = 16,
         horiz  = TRUE,
         cex    = 0.7,
         bty    = "n")
}

dev.off()
cat("Scatterplots saved to:", file.path(out_dir_diag, "lag_scatterplots.pdf"), "\n")

# ----------------------------------------------------------------------------
# 7. Sliding window chemistry assignment
# ----------------------------------------------------------------------------
# For each biological sampling event, average chemistry from the two months
# bracketing or leading up to the sampling month:
#   Sampling in June or July -> average(June, July)
#   Sampling in August       -> average(July, August)
#   Sampling in September    -> average(August, September)
#
# This matches the original manual Excel approach and is more ecologically
# defensible than a single fixed lag month.
#
# Note: get_chem_month() below averages over whatever rows match a given
# Site-Year-Month -- for a double-visit bio month this includes duplicate
# IDENTICAL chemistry rows (broadcast in 03), so mean() here is already a
# safe no-op and needs no de-duplication fix (unlike section 3 above).

cat("\nBuilding sliding window chemistry dataset...\n")

# We need chemistry in a format we can look up by Site + Year + Month
# Use the long-format vnrp data before pivoting

# For each chem variable, build a site-year-month lookup
chem_long <- vnrp[ , c("Site", "Year", "Month", chem_vars)]

# Function to get chemistry value for a given site, year, and month
get_chem_month <- function(site, year, month, var, chem_df) {
  row <- chem_df[chem_df$Site == site &
                   chem_df$Year == year &
                   chem_df$Month == month, var, drop = TRUE]
  if (length(row) == 0 || all(is.na(row))) return(NA_real_)
  mean(row, na.rm = TRUE)
}

# Define the two-month window for each sampling month
window_months <- function(sampling_month) {
  if (sampling_month %in% c(6, 7)) return(c(6, 7))
  if (sampling_month == 8)          return(c(7, 8))
  if (sampling_month == 9)          return(c(8, 9))
  return(c(NA, NA))
}

# Build output row by row
cat("Processing", nrow(bio), "biological sampling rows...\n")

sw_list <- vector("list", nrow(bio))

for (i in seq_len(nrow(bio))) {
  row         <- bio[i, ]
  s           <- row$Site
  y           <- row$Year
  m           <- row$Month
  win         <- window_months(m)
  
  out_row <- data.frame(
    Site          = s,
    Year          = y,
    Month         = m,
    date_yearmon  = row$date_yearmon,
    Date_sample   = row$Date_sample,
    logCHLa       = row$logCHLa,
    CHLa          = row$CHLa,
    AFDM          = row$AFDM,
    window_m1     = win[1],
    window_m2     = win[2],
    stringsAsFactors = FALSE
  )
  
  for (v in chem_vars) {
    v1 <- if (!is.na(win[1])) get_chem_month(s, y, win[1], v, chem_long) else NA_real_
    v2 <- if (!is.na(win[2])) get_chem_month(s, y, win[2], v, chem_long) else NA_real_
    out_row[[v]] <- mean(c(v1, v2), na.rm = TRUE)
    # If both months are NA, result should be NA not NaN
    if (is.nan(out_row[[v]])) out_row[[v]] <- NA_real_
  }
  
  sw_list[[i]] <- out_row
}

sw_data <- do.call(rbind, sw_list)

cat("Sliding window dataset dimensions:",
    nrow(sw_data), "rows x", ncol(sw_data), "cols\n")

# Recalculate DIN from windowed NH4 + NO3
if (all(c("NH4_mg_L", "NO3_mg_L") %in% names(sw_data))) {
  sw_data$DIN_mg_L <- sw_data$NH4_mg_L + sw_data$NO3_mg_L
  sw_data$DIN_mg_L[is.nan(sw_data$DIN_mg_L)] <- NA_real_
}

cat("\n--- Sliding Window Dataset: Non-NA counts per variable ---\n")
for (v in c("CHLa", "AFDM", chem_vars, "DIN_mg_L", "Date_sample")) {
  if (v %in% names(sw_data)) {
    cat(sprintf("  %-15s  n = %d\n", v, sum(!is.na(sw_data[[v]]))))
  }
}

# ----------------------------------------------------------------------------
# 8. Write outputs
# ----------------------------------------------------------------------------

# Full correlation table
write_csv(cor_results,
          file.path(out_dir_inc, "lag_correlations_full.csv"))

# Best lag per variable
write_csv(best_lag,
          file.path(out_dir_inc, "lag_selection.csv"))

# Sliding window chemistry dataset — input to join step
write_csv(sw_data,
          file.path(out_dir_inc, "vnrp_sliding_window.csv"))

cat("\nSaved:\n")
cat("  2_incremental/lag_correlations_full.csv\n")
cat("  2_incremental/lag_selection.csv\n")
cat("  2_incremental/vnrp_sliding_window.csv\n")
cat("\nDone.\n")
cat("\nNOTE: Review lag_scatterplots.pdf and compare lag_selection.csv\n")
cat("against existing n=221 dataset before proceeding to join step.\n")