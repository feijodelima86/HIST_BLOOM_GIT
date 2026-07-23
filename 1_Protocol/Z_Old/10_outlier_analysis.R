# ============================================================================
# 11_outlier_analysis.R
# UCFR Filamentous Algae Project
# Stage 10: Outlier and influence analysis
#
# Inputs:   2_incremental/ucfr_model_ready.csv
#           2_incremental/brt_bloom_fitted.csv
#           2_incremental/sem_jackknife_predictions.csv
#           3_models/brt_bloom_fitted.rds
# Outputs:  2_incremental/outlier_summary.csv
#           4_products/diagnostics/outlier_analysis.pdf
#
# Steps:
#   1. Full model residuals — which obs does trained model struggle with
#   2. Jackknife residuals — which obs are hardest to predict when withheld
#   3. Year-level influence — years with consistently high residuals
#   4. Site x year influence — specific cells driving model behavior
#   5. Predictor space outliers — Mahalanobis distance
#   6. Cross-reference full model vs jackknife outliers
# ============================================================================

library(readr)
library(dismo)
library(gbm)

# ----------------------------------------------------------------------------
# 1. Read data
# ----------------------------------------------------------------------------

cat("Reading data...\n")
dat    <- as.data.frame(read_csv("2_incremental/ucfr_model_ready.csv",
                                 show_col_types = FALSE))
fitted <- as.data.frame(read_csv("2_incremental/brt_bloom_fitted.csv",
                                 show_col_types = FALSE))
jk     <- as.data.frame(read_csv("2_incremental/sem_jackknife_predictions.csv",
                                 show_col_types = FALSE))

for (d in c("2_incremental", "4_products/diagnostics")) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

dat$logCHLa <- log10(dat$CHLa)

# Merge full model fitted values onto main dataset
full_mod <- merge(dat[ , c("Site", "Year", "Month", "logCHLa",
                           "anomaly", "Q_obs_cfs", "Temp_oC",
                           "Days_Since_Freshet", "SPC",
                           "TP_mg_L", "TN_mg_L")],
                  fitted[ , c("Site", "Year", "Month",
                              "pred_logCHLa", "resid_logCHLa")],
                  by = c("Site", "Year", "Month"), all.x = FALSE)

# Jackknife residuals
jk$resid_jk <- jk$Observed - jk$Predicted

cat(sprintf("Full model observations: %d\n", nrow(full_mod)))
cat(sprintf("Jackknife observations:  %d\n\n", nrow(jk)))

site_order <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")
site_cols  <- c(DL = "#E41A1C", GR = "#377EB8", BN = "#4DAF4A",
                MS = "#984EA3", BM = "#FF7F00", HU = "#A65628",
                FH = "#F781BF")

# ----------------------------------------------------------------------------
# 2. Define outlier thresholds
# ----------------------------------------------------------------------------

# Threshold: absolute residual > 2 SD from mean residual
full_resid_sd  <- sd(full_mod$resid_logCHLa, na.rm = TRUE)
full_resid_mn  <- mean(full_mod$resid_logCHLa, na.rm = TRUE)
full_thresh    <- 2 * full_resid_sd

jk_resid_sd    <- sd(jk$resid_jk, na.rm = TRUE)
jk_resid_mn    <- mean(jk$resid_jk, na.rm = TRUE)
jk_thresh      <- 2 * jk_resid_sd

full_mod$outlier_full <- abs(full_mod$resid_logCHLa - full_resid_mn) > full_thresh
jk$outlier_jk         <- abs(jk$resid_jk - jk_resid_mn) > jk_thresh

cat(sprintf("Full model outlier threshold: |resid| > %.3f log units\n", full_thresh))
cat(sprintf("Jackknife outlier threshold:  |resid| > %.3f log units\n\n", jk_thresh))

cat(sprintf("Full model outliers: %d of %d (%.1f%%)\n",
            sum(full_mod$outlier_full, na.rm = TRUE), nrow(full_mod),
            100 * mean(full_mod$outlier_full, na.rm = TRUE)))
cat(sprintf("Jackknife outliers:  %d of %d (%.1f%%)\n\n",
            sum(jk$outlier_jk, na.rm = TRUE), nrow(jk),
            100 * mean(jk$outlier_jk, na.rm = TRUE)))

# ----------------------------------------------------------------------------
# 3. Year-level influence
# ----------------------------------------------------------------------------

cat("--- Year-level Mean Absolute Residual ---\n")
cat(sprintf("  %-6s  %10s  %10s  %8s\n",
            "Year", "MAR_full", "MAR_jk", "n"))
cat(paste(rep("-", 42), collapse = ""), "\n")

years_full <- sort(unique(full_mod$Year))
years_jk   <- sort(unique(jk$Year))
all_yrs    <- sort(union(years_full, years_jk))

yr_influence <- data.frame(
  Year     = all_yrs,
  MAR_full = NA_real_,
  MAR_jk   = NA_real_,
  n_full   = NA_integer_,
  n_jk     = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_along(all_yrs)) {
  yr <- all_yrs[i]
  df <- full_mod[full_mod$Year == yr, ]
  dj <- jk[jk$Year == yr, ]
  
  mar_f <- if (nrow(df) > 0) mean(abs(df$resid_logCHLa), na.rm = TRUE) else NA
  mar_j <- if (nrow(dj) > 0) mean(abs(dj$resid_jk),      na.rm = TRUE) else NA
  
  yr_influence$MAR_full[i] <- round(mar_f, 3)
  yr_influence$MAR_jk[i]   <- round(mar_j, 3)
  yr_influence$n_full[i]   <- nrow(df)
  yr_influence$n_jk[i]     <- nrow(dj)
  
  cat(sprintf("  %-6d  %10s  %10s  %8d\n",
              yr,
              if (!is.na(mar_f)) sprintf("%.3f", mar_f) else "—",
              if (!is.na(mar_j)) sprintf("%.3f", mar_j) else "—",
              max(nrow(df), nrow(dj))))
}

# Flag high-influence years
high_inf_full <- yr_influence$Year[!is.na(yr_influence$MAR_full) &
                                     yr_influence$MAR_full >
                                     mean(yr_influence$MAR_full, na.rm=TRUE) +
                                     sd(yr_influence$MAR_full, na.rm=TRUE)]
high_inf_jk   <- yr_influence$Year[!is.na(yr_influence$MAR_jk) &
                                     yr_influence$MAR_jk >
                                     mean(yr_influence$MAR_jk, na.rm=TRUE) +
                                     sd(yr_influence$MAR_jk, na.rm=TRUE)]

cat(sprintf("\nHigh-influence years (full model): %s\n",
            paste(high_inf_full, collapse = ", ")))
cat(sprintf("High-influence years (jackknife):  %s\n\n",
            paste(high_inf_jk, collapse = ", ")))

# ----------------------------------------------------------------------------
# 4. Site x Year outlier cells
# ----------------------------------------------------------------------------

cat("--- Site x Year Outlier Cells ---\n\n")

cat("Full model outliers:\n")
full_out <- full_mod[full_mod$outlier_full & !is.na(full_mod$outlier_full), ]
full_out <- full_out[order(abs(full_out$resid_logCHLa), decreasing = TRUE), ]
cat(sprintf("  %-6s  %-6s  %-6s  %10s  %10s  %10s\n",
            "Site", "Year", "Month", "Observed", "Predicted", "Residual"))
cat(paste(rep("-", 58), collapse = ""), "\n")
for (i in seq_len(nrow(full_out))) {
  cat(sprintf("  %-6s  %-6d  %-6d  %10.3f  %10.3f  %10.3f\n",
              full_out$Site[i], full_out$Year[i], full_out$Month[i],
              full_out$logCHLa[i],
              full_out$logCHLa[i] - full_out$resid_logCHLa[i],
              full_out$resid_logCHLa[i]))
}

cat("\nJackknife outliers:\n")
jk_out <- jk[jk$outlier_jk & !is.na(jk$outlier_jk), ]
jk_out <- jk_out[order(abs(jk_out$resid_jk), decreasing = TRUE), ]
cat(sprintf("  %-6s  %-6s  %-6s  %10s  %10s  %10s\n",
            "Site", "Year", "Month", "Observed", "Predicted", "Residual"))
cat(paste(rep("-", 58), collapse = ""), "\n")
for (i in seq_len(nrow(jk_out))) {
  cat(sprintf("  %-6s  %-6d  %-6d  %10.3f  %10.3f  %10.3f\n",
              jk_out$Site[i], jk_out$Year[i], jk_out$Month[i],
              jk_out$Observed[i], jk_out$Predicted[i],
              jk_out$resid_jk[i]))
}

# Cross-reference — outliers in both
cat("\nOutliers in BOTH full model and jackknife:\n")
full_keys <- paste(full_out$Site, full_out$Year, full_out$Month)
jk_keys   <- paste(jk_out$Site,   jk_out$Year,   jk_out$Month)
both_keys <- intersect(full_keys, jk_keys)
if (length(both_keys) > 0) {
  cat(paste(" ", both_keys, collapse = "\n"), "\n")
} else {
  cat("  None\n")
}

# ----------------------------------------------------------------------------
# 5. Mahalanobis distance — predictor space outliers
# ----------------------------------------------------------------------------

cat("\n--- Predictor Space Outliers (Mahalanobis Distance) ---\n")

pred_vars <- c("anomaly", "Q_obs_cfs", "Temp_oC", "Days_Since_Freshet")
pred_dat  <- full_mod[ , pred_vars]
pred_dat  <- pred_dat[complete.cases(pred_dat), ]

mah_dist  <- mahalanobis(pred_dat,
                         center = colMeans(pred_dat),
                         cov    = cov(pred_dat))

# Chi-squared threshold with 4 df, p=0.001
mah_thresh <- qchisq(0.999, df = length(pred_vars))
mah_flag   <- mah_dist > mah_thresh

full_mod_cc           <- full_mod[complete.cases(full_mod[ , pred_vars]), ]
full_mod_cc$mah_dist  <- mah_dist
full_mod_cc$mah_flag  <- mah_flag

cat(sprintf("Mahalanobis threshold (chi² p=0.001, df=%d): %.2f\n",
            length(pred_vars), mah_thresh))
cat(sprintf("Predictor space outliers: %d of %d (%.1f%%)\n\n",
            sum(mah_flag), length(mah_flag),
            100 * mean(mah_flag)))

mah_out <- full_mod_cc[full_mod_cc$mah_flag, ]
mah_out <- mah_out[order(mah_out$mah_dist, decreasing = TRUE), ]
if (nrow(mah_out) > 0) {
  cat(sprintf("  %-6s  %-6s  %8s  %8s  %8s  %8s  %10s\n",
              "Site", "Year", "anomaly", "Q_obs", "Temp", "DSF", "Mah_dist"))
  cat(paste(rep("-", 65), collapse = ""), "\n")
  for (i in seq_len(nrow(mah_out))) {
    cat(sprintf("  %-6s  %-6d  %8.3f  %8.1f  %8.1f  %8.1f  %10.2f\n",
                mah_out$Site[i], mah_out$Year[i],
                mah_out$anomaly[i], mah_out$Q_obs_cfs[i],
                mah_out$Temp_oC[i], mah_out$Days_Since_Freshet[i],
                mah_out$mah_dist[i]))
  }
}

# ----------------------------------------------------------------------------
# 6. Plots
# ----------------------------------------------------------------------------

cat("\nGenerating outlier plots...\n")

pdf("4_products/diagnostics/outlier_analysis.pdf",
    width = 12, height = 14)

par(mfrow = c(4, 2), mar = c(4, 4, 3, 1))

# 1. Full model residuals by year
boxplot(resid_logCHLa ~ Year, data = full_mod,
        xlab = "Year", ylab = "Residual (log10 units)",
        main = "Full Model Residuals by Year",
        col  = "lightblue", border = "steelblue",
        las = 2, cex.axis = 0.7)
abline(h = 0,            col = "red",    lty = 2)
abline(h =  full_thresh, col = "orange", lty = 3)
abline(h = -full_thresh, col = "orange", lty = 3)

# 2. Jackknife residuals by year
boxplot(resid_jk ~ Year, data = jk,
        xlab = "Year", ylab = "Residual (log10 units)",
        main = "Jackknife Residuals by Year",
        col  = "lightyellow", border = "darkorange",
        las = 2, cex.axis = 0.7)
abline(h = 0,          col = "red",    lty = 2)
abline(h =  jk_thresh, col = "orange", lty = 3)
abline(h = -jk_thresh, col = "orange", lty = 3)

# 3. Full model residuals by site
boxplot(resid_logCHLa ~ factor(Site, levels = site_order),
        data = full_mod,
        xlab = "Site", ylab = "Residual (log10 units)",
        main = "Full Model Residuals by Site",
        col  = site_cols[site_order], border = "grey30",
        las = 1)
abline(h = 0,            col = "red",    lty = 2)
abline(h =  full_thresh, col = "orange", lty = 3)
abline(h = -full_thresh, col = "orange", lty = 3)

# 4. Jackknife residuals by site
boxplot(resid_jk ~ factor(Site, levels = site_order),
        data = jk,
        xlab = "Site", ylab = "Residual (log10 units)",
        main = "Jackknife Residuals by Site",
        col  = site_cols[site_order], border = "grey30",
        las = 1)
abline(h = 0,          col = "red",    lty = 2)
abline(h =  jk_thresh, col = "orange", lty = 3)
abline(h = -jk_thresh, col = "orange", lty = 3)

# 5. Year MAR comparison
yr_both <- yr_influence[!is.na(yr_influence$MAR_full) &
                          !is.na(yr_influence$MAR_jk), ]
y_rng <- range(c(yr_both$MAR_full, yr_both$MAR_jk), na.rm = TRUE)
plot(yr_both$Year, yr_both$MAR_full,
     type = "b", pch = 16, col = "steelblue",
     xlab = "Year", ylab = "Mean Absolute Residual",
     main = "Year-level Influence: Full Model vs Jackknife",
     ylim = c(0, y_rng[2] * 1.1))
lines(yr_both$Year, yr_both$MAR_jk,
      type = "b", pch = 17, col = "darkorange")
legend("topleft", legend = c("Full model", "Jackknife"),
       col = c("steelblue", "darkorange"), pch = c(16, 17),
       lty = 1, bty = "n", cex = 0.85)
# Highlight high-influence years
for (yr in union(high_inf_full, high_inf_jk)) {
  abline(v = yr, col = "red", lty = 3, lwd = 0.8)
}

# 6. Mahalanobis distance
plot(seq_len(nrow(full_mod_cc)), sort(full_mod_cc$mah_dist, decreasing = TRUE),
     type = "h", col = ifelse(sort(full_mod_cc$mah_dist, decreasing=TRUE) >
                                mah_thresh, "red", "steelblue"),
     xlab = "Observation rank",
     ylab = "Mahalanobis distance",
     main = "Predictor Space Outliers (Mahalanobis)")
abline(h = mah_thresh, col = "red", lty = 2)
text(sum(mah_flag) + 2, mah_thresh + 0.5,
     sprintf("p=0.001 threshold\n(%.1f)", mah_thresh),
     col = "red", cex = 0.7, adj = 0)

# 7. Full model: observed vs predicted, outliers highlighted
rng_f <- range(c(full_mod$logCHLa, full_mod$logCHLa - full_mod$resid_logCHLa),
               na.rm = TRUE)
plot(full_mod$logCHLa - full_mod$resid_logCHLa, full_mod$logCHLa,
     xlim = rng_f, ylim = rng_f,
     xlab = "Predicted log10(CHLa)",
     ylab = "Observed log10(CHLa)",
     main = "Full Model: Outliers Highlighted",
     pch  = ifelse(full_mod$outlier_full, 17, 16),
     col  = ifelse(full_mod$outlier_full, "red",
                   adjustcolor(site_cols[full_mod$Site], 0.6)),
     cex  = ifelse(full_mod$outlier_full, 1.2, 0.7))
abline(0, 1, col = "grey40", lty = 2)
legend("topleft", legend = c("Normal", "Outlier"),
       pch = c(16, 17), col = c("grey50", "red"),
       bty = "n", cex = 0.8)

# 8. Jackknife: observed vs predicted, outliers highlighted
rng_j <- range(c(jk$Observed, jk$Predicted), na.rm = TRUE)
plot(jk$Predicted, jk$Observed,
     xlim = rng_j, ylim = rng_j,
     xlab = "Predicted log10(CHLa)",
     ylab = "Observed log10(CHLa)",
     main = "Jackknife: Outliers Highlighted",
     pch  = ifelse(jk$outlier_jk, 17, 16),
     col  = ifelse(jk$outlier_jk, "red",
                   adjustcolor(site_cols[jk$Site], 0.6)),
     cex  = ifelse(jk$outlier_jk, 1.2, 0.7))
abline(0, 1, col = "grey40", lty = 2)
legend("topleft", legend = c("Normal", "Outlier"),
       pch = c(16, 17), col = c("grey50", "red"),
       bty = "n", cex = 0.8)

dev.off()
cat("Plots saved to 4_products/diagnostics/outlier_analysis.pdf\n")

# ----------------------------------------------------------------------------
# 7. Save summary
# ----------------------------------------------------------------------------

outlier_summary <- merge(
  full_mod[ , c("Site", "Year", "Month", "logCHLa",
                "resid_logCHLa", "outlier_full")],
  jk[ , c("Site", "Year", "Month", "resid_jk", "outlier_jk")],
  by = c("Site", "Year", "Month"), all = TRUE
)

write_csv(outlier_summary, "2_incremental/outlier_summary.csv")
cat("Summary saved to 2_incremental/outlier_summary.csv\n")
cat("Done.\n")
