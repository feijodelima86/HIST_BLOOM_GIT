# ============================================================================
# 15b_plot_projections_reworked.R
# UCFR Cladophora Pipeline -- projection figure
#
# Per-site bloom trajectories, low vs high scenario bracket:
#   - median logCHLa line
#   - 10th-90th percentile ribbon (across ensemble members)
#   - observed historical points (grey) + lowess smooth through them
#
# Panels ordered upstream -> downstream; shared x/y-axis for cross-site
# comparability, now spanning the observed record through the projection
# horizon.
#
# Observed layer: ucfr_model_ready.csv only carries Year/Month (no exact
# sample date), so double-visit months plot two points at nearly the same
# x position -- expected, not a bug.
#
# Input :  2_incremental/bloom_projections.csv
#          2_incremental/ucfr_model_ready.csv       (observed historical)
# Output:  4_products/bloom_trajectories.pdf
# ============================================================================

PATH_SUMMARY <- "2_incremental/bloom_projections.csv"
PATH_OBS     <- "2_incremental/ucfr_model_ready.csv"
OUT_PDF      <- "4_products/bloom_trajectories.pdf"

site_order <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")   # upstream -> downstream

col_low   <- "#377EB8"                       # blue  = low (RCP4.5 / SSP245)
col_high  <- "#E41A1C"                        # red   = high (RCP8.5 / SSP585)
fill_low  <- adjustcolor(col_low,  alpha.f = 0.22)
fill_high <- adjustcolor(col_high, alpha.f = 0.22)

col_obs_pt <- adjustcolor("grey45", alpha.f = 0.6)   # observed points
col_obs_sm <- "grey20"                                # observed lowess smooth

# ----------------------------------------------------------------------------
s <- read.csv(PATH_SUMMARY, stringsAsFactors = FALSE)

obs <- read.csv(PATH_OBS, stringsAsFactors = FALSE)
obs <- obs[!is.na(obs$logCHLa), ]
obs$decimal_year <- obs$Year + (obs$Month - 0.5) / 12

xlim <- range(c(obs$decimal_year, s$year), na.rm = TRUE)
ylim <- range(c(obs$logCHLa, s$p10_logCHLa, s$p90_logCHLa), na.rm = TRUE)
ylim <- ylim + c(-0.03, 0.03) * diff(ylim)

ribbon <- function(d, fill) {
  d <- d[order(d$year), ]
  polygon(c(d$year, rev(d$year)),
          c(d$p10_logCHLa, rev(d$p90_logCHLa)),
          col = fill, border = NA)
}
med_line <- function(d, col) {
  d <- d[order(d$year), ]
  lines(d$year, d$median_logCHLa, col = col, lwd = 2)
}

if (!dir.exists("4_products")) dir.create("4_products", recursive = TRUE)

pdf(OUT_PDF, width = 9, height = 11)
par(mfrow = c(4, 2), mar = c(3.2, 3.6, 2.2, 0.8),
    mgp = c(2.1, 0.6, 0), oma = c(2.2, 1.5, 2.4, 0.5))

for (st in site_order) {
  ds <- s[s$Site == st, ]
  lo <- ds[ds$bracket == "low", ]
  hi <- ds[ds$bracket == "high", ]
  do <- obs[obs$Site == st, ]
  
  plot(NA, xlim = xlim, ylim = ylim,
       xlab = "Year", ylab = expression(log[10] * " Chl " * italic(a)), main = st)
  
  ribbon(lo, fill_low);  ribbon(hi, fill_high)
  
  points(do$decimal_year, do$logCHLa, pch = 16, cex = 0.5, col = col_obs_pt)
  if (nrow(do) >= 4) {
    sm <- lowess(do$decimal_year, do$logCHLa, f = 2/3)
    lines(sm$x, sm$y, col = col_obs_sm, lwd = 1.5)
  }
  
  med_line(lo, col_low); med_line(hi, col_high)
  box()
}

# legend panel (8th cell)
plot.new()
leg <- c("Observed", "Observed smooth (lowess)",
         "Low median (RCP4.5/SSP245)", "High median (RCP8.5/SSP585)",
         "10th-90th percentile")
lcol <- c(col_obs_pt, col_obs_sm, col_low, col_high, "grey50")
llwd <- c(NA, 1.5, 2, 2, NA); lpch <- c(16, NA, NA, NA, 15)
legend("center", legend = leg, col = lcol, lwd = llwd, pch = lpch,
       pt.cex = c(1.2, NA, NA, NA, 2.2), bty = "n", cex = 1.0)

mtext("Projected Cladophora bloom trajectories (NCAR ensemble)",
      outer = TRUE, cex = 1.05, font = 2, line = 0.6)
mtext("panels upstream -> downstream  |  shaded median band = 10th-90th pct",
      outer = TRUE, side = 1, cex = 0.75, line = 0.6)

dev.off()
cat("Wrote", OUT_PDF, "\n")