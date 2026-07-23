# ============================================================================
# 13b_plot_projections.R
# UCFR Cladophora Pipeline -- projection figure
#
# Per-site bloom trajectories, low vs high scenario bracket:
#   - median logCHLa line
#   - 10th-90th percentile ribbon (across ensemble members)
#
# Panels ordered upstream -> downstream; shared y-axis for cross-site
# comparability.
#
# Input :  2_incremental/bloom_projections.csv
# Output:  4_products/bloom_trajectories.pdf
# ============================================================================

PATH_SUMMARY <- "2_incremental/bloom_projections.csv"
OUT_PDF      <- "4_products/bloom_trajectories.pdf"

site_order <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")   # upstream -> downstream

col_low   <- "#377EB8"                       # blue  = low (RCP4.5 / SSP245)
col_high  <- "#E41A1C"                        # red   = high (RCP8.5 / SSP585)
fill_low  <- adjustcolor(col_low,  alpha.f = 0.22)
fill_high <- adjustcolor(col_high, alpha.f = 0.22)

# ----------------------------------------------------------------------------
s <- read.csv(PATH_SUMMARY, stringsAsFactors = FALSE)

xlim <- range(s$year, na.rm = TRUE)
ylim <- range(c(s$p10_logCHLa, s$p90_logCHLa), na.rm = TRUE)
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
  
  plot(NA, xlim = xlim, ylim = ylim,
       xlab = "Year", ylab = expression(log[10] * " Chl " * italic(a)), main = st)
  
  ribbon(lo, fill_low);  ribbon(hi, fill_high)
  med_line(lo, col_low); med_line(hi, col_high)
  box()
}

# legend panel (8th cell)
plot.new()
leg <- c("Low median (RCP4.5/SSP245)", "High median (RCP8.5/SSP585)",
         "10th-90th percentile")
lcol <- c(col_low, col_high, "grey50")
llwd <- c(2, 2, NA); lpch <- c(NA, NA, 15)
legend("center", legend = leg, col = lcol, lwd = llwd, pch = lpch,
       pt.cex = 2.2, bty = "n", cex = 1.05)

mtext("Projected Cladophora bloom trajectories (NCAR ensemble)",
      outer = TRUE, cex = 1.05, font = 2, line = 0.6)
mtext("panels upstream -> downstream  |  shaded median band = 10th-90th pct",
      outer = TRUE, side = 1, cex = 0.75, line = 0.6)

dev.off()
cat("Wrote", OUT_PDF, "\n")