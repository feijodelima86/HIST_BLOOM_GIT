# ============================================================================
# 15_results_figures.R
# Results-section figures (main text).
#
# Section A: Figure 1 -- GAM smooth-term panel, all six M1 partial effects.
#   Panel order mirrors the M1 formula (not an influence ranking): lag_y,
#   anomaly, logQ, Days_Since_Freshet, logTP, Temp_oC. select = 1:6 excludes
#   s(Site, bs="re"), the 7th smooth in the model object.
#   No main= per panel -- mgcv's default y-axis label already reads
#   "s(term, edf)", which identifies the panel; a title would repeat it.
#   seWithMean = TRUE: CI bands include uncertainty about the overall mean,
#   not just the smooth's shape around it (Wood, mgcv) -- flip to FALSE for
#   the narrower shape-only band if you want it to match the summary table's
#   edf/F-test framing instead.
# Section B: Site/reach map -- not yet built.
# Section C: Figure 4 -- hydrology driver drift by reach, 4 NCAR reaches,
#   ensemble mean + p10-p90 across all members, no bracket split. Three
#   panels: anomaly, mean discharge, peak timing -- the raw envelope
#   columns validated in 04a/07c (D2 resolution), not the site-anchored,
#   delta-scaled predictors M1 sees during projection. Temp_oC deliberately
#   excluded -- site-level not reach-level, deterministic per bracket (no
#   ensemble spread), monotonic warming everywhere. Belongs in 3.4 instead,
#   paired with the calibration-ceiling crossing story.
# Section D: Supplementary -- temperature projections by site (3.4). Single
#   panel, 7 sites (not grouped by reach -- MS/BM/HU share CLABE but diverge
#   sharply on ceiling-crossing, so reach grouping would be misleading here).
#   Both brackets per site, calibration ceiling marked, 2026-2098 window.
# Section E: Supplementary -- ensemble variance-by-year diagnostic (3.4).
#   SD of pred_logCHLa across members within each Site x bracket x year
#   cell, then averaged across the 7 sites per bracket x year. Growing SD
#   (not collapsing) supports the recursion-relaxation caveat.
#
# Inputs : 3_models/bloom_model_M1.rds
#          2_incremental/ncar_discharge_envelope.csv (Section C)
#          2_incremental/ncar_temperature_envelope.csv (Section D)
#          2_incremental/bloom_projections_members.csv (Section E)
# Outputs: 4_products/gam_smooth_panel.pdf
#          4_products/hydrology_drift_by_reach.pdf
#          4_products/temp_projections_by_site.pdf
#          4_products/ensemble_variance_by_year.pdf
# ============================================================================

library(mgcv)

m_M1 <- readRDS("3_models/bloom_model_M1.rds")

# ----------------------------------------------------------------------------
# Section A: Figure 1 -- GAM smooth-term panel
# ----------------------------------------------------------------------------

xlab_labels <- c(
  "Lag-year max log10(CHLa)",
  "(Peak Q / baseflow Q)^(1/3)",
  "log10(Q, cfs)",
  "Days since freshet",
  "log10(TP, mg/L)",
  "Temperature (deg C)"
)

if (!dir.exists("4_products")) dir.create("4_products", recursive = TRUE)
pdf("4_products/gam_smooth_panel.pdf", width = 9, height = 6, family = "Helvetica")
par(mfrow = c(2, 3), mar = c(4, 4.5, 2, 1))

for (i in 1:6) {
  plot(m_M1,
       select     = i,
       residuals  = TRUE,
       shade      = TRUE,
       shade.col  = "gray85",
       seWithMean = TRUE,
       pch        = 1,
       cex        = 0.6,
       col        = "steelblue4",
       lwd        = 2,
       xlab       = xlab_labels[i])
  abline(h = 0, lty = 2, col = "gray60")
}

dev.off()
cat("Wrote 4_products/gam_smooth_panel.pdf\n")

# ----------------------------------------------------------------------------
# Section B: Site/reach map -- not yet built
# ----------------------------------------------------------------------------
# Needs site lat/lon (not currently an input to this script) and sf for the
# base layer, per project package conventions. Add here once scoped.

# ----------------------------------------------------------------------------
# Section C: Figure 4 -- hydrology driver drift by reach
# ----------------------------------------------------------------------------

env <- read.csv("2_incremental/ncar_discharge_envelope.csv", stringsAsFactors = FALSE)

reach_order  <- c("CLALO", "CLADR", "CLABE", "CLAPL")
reach_colors <- c(CLALO = "#08519c", CLADR = "#3182bd",
                  CLABE = "#6baed6", CLAPL = "#bdd7e7")

# Ensemble mean + p10-p90 across all members, per reach x water_year.
drift_summary <- function(col) {
  x    <- env[[col]]
  keys <- list(site = env$site, water_year = env$water_year)
  m    <- aggregate(x, by = keys, FUN = mean, na.rm = TRUE)
  p10  <- aggregate(x, by = keys, FUN = function(v) quantile(v, 0.10, na.rm = TRUE))
  p90  <- aggregate(x, by = keys, FUN = function(v) quantile(v, 0.90, na.rm = TRUE))
  names(m)[3] <- "mean"; names(p10)[3] <- "p10"; names(p90)[3] <- "p90"
  merge(merge(m, p10, by = c("site", "water_year")), p90, by = c("site", "water_year"))
}

plot_drift_panel <- function(d, ylab, log = "") {
  xr <- range(d$water_year, na.rm = TRUE)
  yr <- range(c(d$p10, d$p90), na.rm = TRUE)
  plot(NA, xlim = xr, ylim = yr, xlab = "Water year", ylab = ylab, log = log)
  for (r in reach_order) {
    dd <- d[d$site == r, ]; dd <- dd[order(dd$water_year), ]
    polygon(c(dd$water_year, rev(dd$water_year)), c(dd$p90, rev(dd$p10)),
            col = adjustcolor(reach_colors[r], alpha.f = 0.15), border = NA)
  }
  for (r in reach_order) {
    dd <- d[d$site == r, ]; dd <- dd[order(dd$water_year), ]
    lines(dd$water_year, dd$mean, col = reach_colors[r], lwd = 2)
  }
}

pdf("4_products/hydrology_drift_by_reach.pdf", width = 12, height = 4.5, family = "Helvetica")
par(mfrow = c(1, 3), mar = c(4, 4.5, 2, 1))

plot_drift_panel(drift_summary("anomaly_ma20"),
                 "Anomaly, (peak Q / baseflow Q)^(1/3), MA20")
plot_drift_panel(drift_summary("mean_q_cfs_ma20"),
                 "Mean discharge, MA20 (cfs)", log = "y")
plot_drift_panel(drift_summary("days_since_wy_start_ma20"),
                 "Peak timing, MA20 (days since WY start)")

legend("topright", legend = reach_order, col = reach_colors[reach_order],
       lwd = 2, bty = "n", cex = 0.8)

dev.off()
cat("Wrote 4_products/hydrology_drift_by_reach.pdf\n")

# ----------------------------------------------------------------------------
# Section D: Supplementary -- temperature projections by site
# ----------------------------------------------------------------------------

temp_env <- read.csv("2_incremental/ncar_temperature_envelope.csv", stringsAsFactors = FALSE)
temp_env <- unique(temp_env[, c("site", "water_year", "Temp_oC_low", "Temp_oC_high")])
temp_env <- temp_env[temp_env$water_year >= 2026 & temp_env$water_year <= 2098, ]

site_order  <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")
site_colors <- setNames(
  c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7"),
  site_order
)
CAL_CEILING <- 22.2   # observed calibration ceiling, deg C -- see 3.4 text

pdf("4_products/temp_projections_by_site.pdf", width = 7, height = 5, family = "Helvetica")
par(mar = c(4, 4.5, 2, 1))

yr <- range(c(temp_env$Temp_oC_low, temp_env$Temp_oC_high), na.rm = TRUE)
plot(NA, xlim = c(2026, 2098), ylim = yr,
     xlab = "Water year", ylab = "Stream temperature, deg C")
abline(h = CAL_CEILING, lty = 2, col = "gray50")

for (s in site_order) {
  d <- temp_env[temp_env$site == s, ]; d <- d[order(d$water_year), ]
  lines(d$water_year, d$Temp_oC_high, col = site_colors[s], lwd = 2, lty = 1)
  lines(d$water_year, d$Temp_oC_low,  col = site_colors[s], lwd = 2, lty = 3)
}

legend("topleft", legend = site_order, col = site_colors[site_order],
       lwd = 2, bty = "n", cex = 0.8, ncol = 2)
legend("bottomright", legend = c("High bracket (RCP8.5/SSP585)", "Low bracket (RCP4.5/SSP245)"),
       lty = c(1, 3), lwd = 2, col = "gray30", bty = "n", cex = 0.75)

dev.off()
cat("Wrote 4_products/temp_projections_by_site.pdf\n")

# ----------------------------------------------------------------------------
# Section E: Supplementary -- ensemble variance-by-year diagnostic
# ----------------------------------------------------------------------------

members <- read.csv("2_incremental/bloom_projections_members.csv", stringsAsFactors = FALSE)

site_sd <- aggregate(pred_logCHLa ~ Site + bracket + year, data = members,
                     FUN = function(x) sd(x, na.rm = TRUE))
names(site_sd)[4] <- "sd_logCHLa"

bracket_sd <- aggregate(sd_logCHLa ~ bracket + year, data = site_sd,
                        FUN = mean, na.rm = TRUE)

bracket_colors <- c(low = "#56B4E9", high = "#D55E00")

pdf("4_products/ensemble_variance_by_year.pdf", width = 7, height = 5, family = "Helvetica")
par(mar = c(4, 4.5, 2, 1))

xr <- range(bracket_sd$year, na.rm = TRUE)
yr <- range(bracket_sd$sd_logCHLa, na.rm = TRUE)
plot(NA, xlim = xr, ylim = yr,
     xlab = "Water year", ylab = "SD of predicted log10(CHLa) across members")

for (b in c("low", "high")) {
  d <- bracket_sd[bracket_sd$bracket == b, ]; d <- d[order(d$year), ]
  lines(d$year, d$sd_logCHLa, col = bracket_colors[b], lwd = 2)
}

legend("topleft", legend = c("Low bracket", "High bracket"),
       col = bracket_colors[c("low", "high")], lwd = 2, bty = "n", cex = 0.85)

dev.off()
cat("Wrote 4_products/ensemble_variance_by_year.pdf\n")