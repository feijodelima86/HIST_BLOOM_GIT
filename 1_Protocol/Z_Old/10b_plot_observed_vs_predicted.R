# ============================================================================
# plot_observed_vs_predicted.R
# UCFR Filamentous Algae Project
# Diagnostic plot: observed log10(CHLa) vs predicted annual trajectory by site
#
# Inputs:
#   2_incremental/ucfr_model_ready.csv       (observed monthly CHLa)
#   2_incremental/projections_annual.csv     (full chain predictions)
#
# Output:
#   4_products/diagnostics/observed_vs_predicted_by_site.pdf
#
# Layout:
#   7 panels stacked vertically (DL, GR, BN, MS, BM, HU, FH)
#   X-axis: year
#   Y-axis: log10(CHLa) in log10(mg/m²) — linear scale, log10-transformed values
#   Black points: observed log10(CHLa)
#   Lines: median predicted log10(CHLa) across ESMs, one per scenario
#   Ribbons: IQR across ESMs
#   Vertical line at 2023.5: observed -> projected transition
# ============================================================================

library(readr)
library(dplyr)

# ----------------------------------------------------------------------------
# 1. Configuration
# ----------------------------------------------------------------------------

SITE_ORDER <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")

COL_OBS    <- "black"
COL_RCP45  <- "#1A5276"
COL_RCP85  <- "#922B21"
COL_RIBBON_45 <- adjustcolor(COL_RCP45, alpha.f = 0.2)
COL_RIBBON_85 <- adjustcolor(COL_RCP85, alpha.f = 0.2)

OUT_DIR  <- "4_products/diagnostics"
OUT_FILE <- file.path(OUT_DIR, "observed_vs_predicted_by_site.pdf")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# 2. Read data
# ----------------------------------------------------------------------------

cat("Reading observed data...\n")
obs <- as.data.frame(read_csv("2_incremental/ucfr_model_ready.csv",
                              show_col_types = FALSE))

cat("Reading annual projections...\n")
proj <- as.data.frame(read_csv("2_incremental/projections_annual.csv",
                               show_col_types = FALSE))

cat(sprintf("Observed rows: %d\n", nrow(obs)))
cat(sprintf("Projection rows: %d\n\n", nrow(proj)))

# ----------------------------------------------------------------------------
# 3. Summarize predictions on log10 scale: median + IQR across ESMs
# ----------------------------------------------------------------------------

# Work directly with pred_logCHLa (already log10) from projections file
pred_summary <- proj %>%
  group_by(site, scenario, year) %>%
  summarise(
    logCHLa_med = median(pred_logCHLa, na.rm = TRUE),
    logCHLa_q25 = quantile(pred_logCHLa, 0.25, na.rm = TRUE),
    logCHLa_q75 = quantile(pred_logCHLa, 0.75, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  as.data.frame()

cat(sprintf("Summarized predictions: %d rows\n\n", nrow(pred_summary)))

# ----------------------------------------------------------------------------
# 4. X range
# ----------------------------------------------------------------------------

x_min <- min(c(obs$Year, pred_summary$year), na.rm = TRUE)
x_max <- max(c(obs$Year, pred_summary$year), na.rm = TRUE)

# ----------------------------------------------------------------------------
# 5. Generate plot
# ----------------------------------------------------------------------------

cat("Generating plot...\n")

pdf(OUT_FILE, width = 11, height = 14)

par(mfrow = c(7, 1), mar = c(2.5, 4.5, 2, 1), oma = c(3, 0, 3, 0))

for (s in SITE_ORDER) {
  
  # Observed: filter to positive CHLa, then log10 transform
  obs_s   <- obs[obs$Site == s & !is.na(obs$CHLa) & obs$CHLa > 0, ]
  obs_log <- log10(obs_s$CHLa)
  
  # Predictions for this site, on log10 scale already
  pred_s  <- pred_summary[pred_summary$site == s, ]
  pred_45 <- pred_s[pred_s$scenario == "RCP4.5", ]
  pred_85 <- pred_s[pred_s$scenario == "RCP8.5", ]
  pred_45 <- pred_45[order(pred_45$year), ]
  pred_85 <- pred_85[order(pred_85$year), ]
  
  # Y range: combine all log10 values
  y_vals <- c(obs_log, pred_s$logCHLa_q25, pred_s$logCHLa_q75)
  y_vals <- y_vals[is.finite(y_vals)]
  y_max  <- max(y_vals) + 0.15
  y_min  <- min(y_vals) - 0.15
  
  # Empty plot (linear axis, log10-transformed values)
  plot(NA, NA,
       xlim = c(x_min, x_max),
       ylim = c(y_min, y_max),
       xlab = "",
       ylab = expression(log[10] ~ CHLa ~ (mg/m^2)),
       main = sprintf("Site %s", s),
       las  = 1,
       cex.main = 1.1,
       font.main = 2)
  
  # Observed/projected boundary
  abline(v = 2023.5, col = "grey60", lty = 2, lwd = 1)
  
  # RCP4.5 ribbon and median line
  if (nrow(pred_45) > 0) {
    polygon(c(pred_45$year, rev(pred_45$year)),
            c(pred_45$logCHLa_q25, rev(pred_45$logCHLa_q75)),
            col = COL_RIBBON_45, border = NA)
    lines(pred_45$year, pred_45$logCHLa_med,
          col = COL_RCP45, lwd = 1.8)
  }
  
  # RCP8.5 ribbon and median line
  if (nrow(pred_85) > 0) {
    polygon(c(pred_85$year, rev(pred_85$year)),
            c(pred_85$logCHLa_q25, rev(pred_85$logCHLa_q75)),
            col = COL_RIBBON_85, border = NA)
    lines(pred_85$year, pred_85$logCHLa_med,
          col = COL_RCP85, lwd = 1.8)
  }
  
  # Observed points
  points(obs_s$Year, obs_log,
         pch = 16, col = COL_OBS, cex = 0.85)
  
  # Sample size
  text(x_min + 2, y_max - 0.05,
       sprintf("n obs = %d", nrow(obs_s)),
       cex = 0.75, col = "grey40", adj = 0)
}

# X-axis label
mtext("Year", side = 1, line = 2.5, cex = 0.95)

# Overall title
mtext("Observed vs predicted log10(CHLa) by site",
      outer = TRUE, side = 3, line = 0.8, cex = 1.15, font = 2)

# Legend at the top
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
    mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n",
     xlab = "", ylab = "")
legend("top",
       legend = c("Observed log10(CHLa)",
                  "RCP4.5 median",
                  "RCP4.5 IQR",
                  "RCP8.5 median",
                  "RCP8.5 IQR",
                  "Obs / Proj boundary"),
       col    = c(COL_OBS, COL_RCP45, COL_RIBBON_45,
                  COL_RCP85, COL_RIBBON_85, "grey60"),
       pch    = c(16, NA, 15, NA, 15, NA),
       lty    = c(NA, 1, NA, 1, NA, 2),
       lwd    = c(NA, 1.8, NA, 1.8, NA, 1),
       horiz  = TRUE,
       cex    = 0.75,
       bty    = "n",
       inset  = c(0, 0.005))

dev.off()
cat(sprintf("Saved: %s\n", OUT_FILE))
cat("Done.\n")