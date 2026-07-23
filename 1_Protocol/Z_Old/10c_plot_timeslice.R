# ============================================================================
# plot_timeslice.R
# UCFR Filamentous Algae Project
# Figure: time slice bloom predictions on log10 scale
#
# Input:    2_incremental/projections_timeslice.csv
# Output:   4_products/diagnostics/timeslice_bloom_by_site.pdf
#
# Layout:
#   X-axis: site (DL, GR, BN, MS, BM, HU, FH — upstream to downstream)
#   Y-axis: log10(CHLa mg/m²) — linear axis, values pre-transformed
#   Four groups per site: RCP4.5 2050, RCP4.5 2080, RCP8.5 2050, RCP8.5 2080
#   Point = median, vertical bar = IQR
#   Optional reference line at baseline (hindcast median) per site
# ============================================================================

library(readr)
library(dplyr)

# ----------------------------------------------------------------------------
# 1. Configuration
# ----------------------------------------------------------------------------

SITE_ORDER <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")

# Colors
COL_RCP45_2050 <- "#5DADE2"   # light blue
COL_RCP45_2080 <- "#1A5276"   # dark blue
COL_RCP85_2050 <- "#E59866"   # light red
COL_RCP85_2080 <- "#922B21"   # dark red

OUT_DIR  <- "4_products/diagnostics"
OUT_FILE <- file.path(OUT_DIR, "timeslice_bloom_by_site.pdf")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# 2. Read data
# ----------------------------------------------------------------------------

cat("Reading time slice projections...\n")
ts <- as.data.frame(read_csv("2_incremental/projections_timeslice.csv",
                             show_col_types = FALSE))

cat("Reading annual projections (for baseline)...\n")
ann <- as.data.frame(read_csv("2_incremental/projections_annual.csv",
                              show_col_types = FALSE))

cat("Reading observed data...\n")
obs <- as.data.frame(read_csv("2_incremental/ucfr_model_ready.csv",
                              show_col_types = FALSE))

cat(sprintf("Timeslice rows: %d\n", nrow(ts)))
cat(sprintf("Annual rows: %d\n", nrow(ann)))
cat(sprintf("Observed rows: %d\n\n", nrow(obs)))

# Compute observed log10(CHLa) summary per site (1998-2023)
obs$logCHLa <- log10(obs$CHLa)
obs_summary <- obs %>%
  filter(Year >= 1998 & Year <= 2023 &
           !is.na(CHLa) & CHLa > 0) %>%
  group_by(Site) %>%
  summarise(
    obs_med = median(logCHLa, na.rm = TRUE),
    obs_q25 = quantile(logCHLa, 0.25, na.rm = TRUE),
    obs_q75 = quantile(logCHLa, 0.75, na.rm = TRUE),
    n_obs   = n(),
    .groups = "drop"
  ) %>%
  as.data.frame()

cat("Observed log10(CHLa) per site (1998-2023):\n")
for (i in seq_len(nrow(obs_summary))) {
  cat(sprintf("  %-6s  med=%.3f  IQR=[%.3f, %.3f]  n=%d\n",
              obs_summary$Site[i], obs_summary$obs_med[i],
              obs_summary$obs_q25[i], obs_summary$obs_q75[i],
              obs_summary$n_obs[i]))
}
cat("\n")

# ----------------------------------------------------------------------------
# 3. Compute hindcast baseline (1998-2023 median log10CHLa)
# ----------------------------------------------------------------------------

baseline <- ann %>%
  filter(year >= 1998 & year <= 2023 &
           scenario == "RCP4.5" & esm == "CanESM5") %>%
  group_by(site) %>%
  summarise(
    baseline_logCHLa = median(pred_logCHLa, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  as.data.frame()

cat("Baseline log10(CHLa) per site:\n")
for (i in seq_len(nrow(baseline))) {
  cat(sprintf("  %-6s  %.3f\n",
              baseline$site[i], baseline$baseline_logCHLa[i]))
}
cat("\n")

# ----------------------------------------------------------------------------
# 4. Build group labels and x-axis positions
# ----------------------------------------------------------------------------

# 5 groups per site, slightly offset: observed + 4 scenario/horizon combos
site_x_center <- setNames(seq_along(SITE_ORDER), SITE_ORDER)
offsets <- c("Observed"    = -0.28,
             "RCP4.5_2050" = -0.14,
             "RCP4.5_2080" = -0.00,
             "RCP8.5_2050" =  0.14,
             "RCP8.5_2080" =  0.28)

group_cols <- c("Observed"    = "#2C3E50",
                "RCP4.5_2050" = COL_RCP45_2050,
                "RCP4.5_2080" = COL_RCP45_2080,
                "RCP8.5_2050" = COL_RCP85_2050,
                "RCP8.5_2080" = COL_RCP85_2080)

ts$group <- paste(ts$scenario, ts$horizon, sep = "_")
ts$x     <- site_x_center[ts$site] + offsets[ts$group]

# ----------------------------------------------------------------------------
# 5. Y range
# ----------------------------------------------------------------------------

y_vals <- c(ts$logCHLa_med, ts$logCHLa_q25, ts$logCHLa_q75,
            baseline$baseline_logCHLa,
            obs_summary$obs_med, obs_summary$obs_q25, obs_summary$obs_q75)
y_vals <- y_vals[is.finite(y_vals)]
y_min  <- min(y_vals) - 0.15
y_max  <- max(y_vals) + 0.15

# ----------------------------------------------------------------------------
# 6. Generate plot
# ----------------------------------------------------------------------------

cat("Generating plot...\n")

pdf(OUT_FILE, width = 10, height = 6)

par(mar = c(4, 4.5, 3, 1), oma = c(0, 0, 0, 0))

plot(NA, NA,
     xlim = c(0.5, length(SITE_ORDER) + 0.5),
     ylim = c(y_min, y_max),
     xlab = "",
     ylab = expression(log[10] ~ CHLa ~ (mg/m^2)),
     xaxt = "n",
     main = "Time slice bloom predictions by site",
     las  = 1,
     font.main = 2)

# Site labels on x-axis
axis(1, at = seq_along(SITE_ORDER), labels = SITE_ORDER, las = 1)
mtext("Site (upstream -> downstream)", side = 1, line = 2.5, cex = 0.95)

# Light vertical dividers between sites
for (i in seq_along(SITE_ORDER)[-length(SITE_ORDER)]) {
  abline(v = i + 0.5, col = "grey85", lty = 1, lwd = 0.5)
}

# Baseline reference lines per site (horizontal segments)
for (i in seq_len(nrow(baseline))) {
  s <- baseline$site[i]
  x_center <- site_x_center[s]
  segments(x0 = x_center - 0.35, x1 = x_center + 0.35,
           y0 = baseline$baseline_logCHLa[i],
           col = "grey40", lty = 3, lwd = 1.2)
}

# Plot observed point + IQR per site first (so projection points draw on top if overlap)
for (i in seq_len(nrow(obs_summary))) {
  s        <- obs_summary$Site[i]
  x_center <- site_x_center[s]
  if (is.na(x_center)) next
  x_pos    <- x_center + offsets["Observed"]
  col      <- group_cols["Observed"]
  
  # IQR bar
  segments(x0 = x_pos, x1 = x_pos,
           y0 = obs_summary$obs_q25[i], y1 = obs_summary$obs_q75[i],
           col = col, lwd = 2.5)
  # Whisker caps
  segments(x0 = x_pos - 0.025, x1 = x_pos + 0.025,
           y0 = obs_summary$obs_q25[i], y1 = obs_summary$obs_q25[i],
           col = col, lwd = 2.5)
  segments(x0 = x_pos - 0.025, x1 = x_pos + 0.025,
           y0 = obs_summary$obs_q75[i], y1 = obs_summary$obs_q75[i],
           col = col, lwd = 2.5)
  # Median point
  points(x_pos, obs_summary$obs_med[i],
         pch = 21, bg = col, col = "black",
         cex = 1.4, lwd = 0.8)
}

# Plot IQR bars and median points for scenario projections
for (i in seq_len(nrow(ts))) {
  d <- ts[i, ]
  col <- group_cols[d$group]
  
  # IQR bar
  segments(x0 = d$x, x1 = d$x,
           y0 = d$logCHLa_q25, y1 = d$logCHLa_q75,
           col = col, lwd = 2.5)
  
  # Whisker caps
  segments(x0 = d$x - 0.025, x1 = d$x + 0.025,
           y0 = d$logCHLa_q25, y1 = d$logCHLa_q25,
           col = col, lwd = 2.5)
  segments(x0 = d$x - 0.025, x1 = d$x + 0.025,
           y0 = d$logCHLa_q75, y1 = d$logCHLa_q75,
           col = col, lwd = 2.5)
  
  # Median point
  points(d$x, d$logCHLa_med,
         pch = 21, bg = col, col = "black",
         cex = 1.4, lwd = 0.8)
}

# Legend
legend("topright",
       legend = c("RCP4.5 2050", "RCP4.5 2080",
                  "RCP8.5 2050", "RCP8.5 2080",
                  "1998-2023 baseline"),
       pch    = c(21, 21, 21, 21, NA),
       pt.bg  = c(COL_RCP45_2050, COL_RCP45_2080,
                  COL_RCP85_2050, COL_RCP85_2080, NA),
       col    = c("black", "black", "black", "black", "grey40"),
       lty    = c(NA, NA, NA, NA, 3),
       lwd    = c(NA, NA, NA, NA, 1.2),
       pt.cex = c(1.4, 1.4, 1.4, 1.4, NA),
       cex    = 0.75,
       bty    = "n",
       inset  = c(0.01, 0.01),
       ncol   = 1)

dev.off()
cat(sprintf("Saved: %s\n", OUT_FILE))
cat("Done.\n")