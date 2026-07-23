# ============================================================================
# diag_projection_uncertainty_check.R
# UCFR Filamentous Algae Project
# Diagnostic: does the projection ribbon understate real predictive
# uncertainty? Not a manuscript deliverable -- this answers whether the
# parked TP Monte Carlo / uncertainty-propagation work is warranted; it
# does not run it.
#
# Section A: observed annual-max overlay on the projection trajectories --
#   apples-to-apples fix for 13b_plot_projections.R, which overlaid ALL
#   visits (within-season noise included) against a projected series that
#   is effectively an annual-max equivalent (13_project_bloom.R's
#   recursion seeds next year's lag_y from each yearly prediction, the
#   same way observed annual max seeds it during fitting).
# Section B: ribbon half-width (p90-p10)/2 from bloom_projections.csv, per
#   site, vs M1's residual SD (in-sample and LOYO), per site. The ribbon
#   only reflects spread across NCAR ensemble members run through M1 as
#   point predictions -- no residual/model uncertainty is added. Turns
#   "looks too narrow" into a fold-multiplier. LOYO is the honest number
#   (out-of-sample); in-sample is shown for contrast and is optimistically
#   biased low, since those are the residuals M1 was fit to minimize.
#
# Inputs:  2_incremental/ucfr_model_ready.csv         (Section A)
#          2_incremental/bloom_projections.csv        (Section A, B)
#          2_incremental/m1_predictions.csv           (Section B, in-sample)
#          2_incremental/temporal_val_predictions.csv (Section B, LOYO)
#
# Outputs: 4_products/diagnostics/observed_annual_max_overlay.pdf
#          4_products/diagnostics/ribbon_vs_residual_sd.pdf
#          4_products/diagnostics/ribbon_vs_residual_sd_table.csv
#          console: comparison table
# ============================================================================

site_order <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")

if (!dir.exists("4_products/diagnostics")) {
  dir.create("4_products/diagnostics", recursive = TRUE)
}

# ============================================================================
# Section A: observed annual-max overlay
# ============================================================================

s <- read.csv("2_incremental/bloom_projections.csv", stringsAsFactors = FALSE)

obs <- read.csv("2_incremental/ucfr_model_ready.csv", stringsAsFactors = FALSE)
obs <- obs[!is.na(obs$logCHLa), ]
obs_annmax <- aggregate(logCHLa ~ Site + Year, data = obs, FUN = max, na.rm = TRUE)

col_low    <- "#377EB8"
col_high   <- "#E41A1C"
fill_low   <- adjustcolor(col_low,  alpha.f = 0.22)
fill_high  <- adjustcolor(col_high, alpha.f = 0.22)
col_obs_pt <- adjustcolor("grey45", alpha.f = 0.7)
col_obs_sm <- "grey20"

xlim <- range(c(obs_annmax$Year, s$year), na.rm = TRUE)
ylim <- range(c(obs_annmax$logCHLa, s$p10_logCHLa, s$p90_logCHLa), na.rm = TRUE)
ylim <- ylim + c(-0.03, 0.03) * diff(ylim)

ribbon <- function(d, fill) {
  d <- d[order(d$year), ]
  polygon(c(d$year, rev(d$year)), c(d$p10_logCHLa, rev(d$p90_logCHLa)),
          col = fill, border = NA)
}
med_line <- function(d, col) {
  d <- d[order(d$year), ]
  lines(d$year, d$median_logCHLa, col = col, lwd = 2)
}

pdf("4_products/diagnostics/observed_annual_max_overlay.pdf",
    width = 9, height = 11, family = "Helvetica")
par(mfrow = c(4, 2), mar = c(3.2, 3.6, 2.2, 0.8),
    mgp = c(2.1, 0.6, 0), oma = c(2.2, 1.5, 2.4, 0.5))

for (st in site_order) {
  ds <- s[s$Site == st, ]
  lo <- ds[ds$bracket == "low", ]
  hi <- ds[ds$bracket == "high", ]
  do <- obs_annmax[obs_annmax$Site == st, ]
  do <- do[order(do$Year), ]
  
  plot(NA, xlim = xlim, ylim = ylim,
       xlab = "Year", ylab = expression(log[10] * " Chl " * italic(a)), main = st)
  
  ribbon(lo, fill_low); ribbon(hi, fill_high)
  
  points(do$Year, do$logCHLa, pch = 16, cex = 0.8, col = col_obs_pt)
  if (nrow(do) >= 4) {
    sm <- lowess(do$Year, do$logCHLa, f = 2/3)
    lines(sm$x, sm$y, col = col_obs_sm, lwd = 1.5)
  }
  
  med_line(lo, col_low); med_line(hi, col_high)
  box()
}

plot.new()
leg <- c("Observed annual max", "Observed smooth (lowess)",
         "Low median (RCP4.5/SSP245)", "High median (RCP8.5/SSP585)",
         "10th-90th percentile")
lcol <- c(col_obs_pt, col_obs_sm, col_low, col_high, "grey50")
llwd <- c(NA, 1.5, 2, 2, NA); lpch <- c(16, NA, NA, NA, 15)
legend("center", legend = leg, col = lcol, lwd = llwd, pch = lpch,
       pt.cex = c(1.4, NA, NA, NA, 2.2), bty = "n", cex = 1.0)

mtext("Diagnostic: observed ANNUAL MAX vs projections (apples-to-apples)",
      outer = TRUE, cex = 1.05, font = 2, line = 0.6)
mtext("compare against bloom_trajectories.pdf, which overlays all visits",
      outer = TRUE, side = 1, cex = 0.75, line = 0.6)

dev.off()
cat("Wrote 4_products/diagnostics/observed_annual_max_overlay.pdf\n")

# ============================================================================
# Section B: ribbon half-width vs residual SD, by site
# ============================================================================

s$hw <- (s$p90_logCHLa - s$p10_logCHLa) / 2

ribbon_mean <- aggregate(hw ~ Site, data = s, FUN = mean, na.rm = TRUE)
ribbon_min  <- aggregate(hw ~ Site, data = s, FUN = min,  na.rm = TRUE)
ribbon_max  <- aggregate(hw ~ Site, data = s, FUN = max,  na.rm = TRUE)
names(ribbon_mean)[2] <- "ribbon_hw_mean"
names(ribbon_min)[2]  <- "ribbon_hw_min"
names(ribbon_max)[2]  <- "ribbon_hw_max"
ribbon_by_site <- merge(merge(ribbon_mean, ribbon_min, by = "Site"),
                        ribbon_max, by = "Site")

# In-sample residual SD, per site
insamp <- read.csv("2_incremental/m1_predictions.csv", stringsAsFactors = FALSE)
insamp <- insamp[insamp$scheme == "in_sample", ]
sd_insamp <- aggregate(Resid ~ Site, data = insamp,
                       FUN = function(x) sd(x, na.rm = TRUE))
names(sd_insamp) <- c("Site", "resid_sd_insample")

# LOYO residual SD, per site
loyo <- read.csv("2_incremental/temporal_val_predictions.csv", stringsAsFactors = FALSE)
loyo <- loyo[loyo$scheme == "LOYO", ]
loyo$Resid <- loyo$Observed - loyo$Predicted
sd_loyo <- aggregate(Resid ~ Site, data = loyo,
                     FUN = function(x) sd(x, na.rm = TRUE))
names(sd_loyo) <- c("Site", "resid_sd_loyo")

tab <- merge(ribbon_by_site, sd_insamp, by = "Site")
tab <- merge(tab, sd_loyo, by = "Site")
tab$fold_insample <- round(tab$resid_sd_insample / tab$ribbon_hw_mean, 1)
tab$fold_loyo      <- round(tab$resid_sd_loyo / tab$ribbon_hw_mean, 1)

num_cols <- c("ribbon_hw_mean", "ribbon_hw_min", "ribbon_hw_max",
              "resid_sd_insample", "resid_sd_loyo")
tab[, num_cols] <- round(tab[, num_cols], 3)
tab <- tab[match(site_order, tab$Site), ]

overall <- data.frame(
  Site              = "ALL",
  ribbon_hw_mean    = round(mean(s$hw, na.rm = TRUE), 3),
  ribbon_hw_min     = round(min(s$hw, na.rm = TRUE), 3),
  ribbon_hw_max     = round(max(s$hw, na.rm = TRUE), 3),
  resid_sd_insample = round(sd(insamp$Resid, na.rm = TRUE), 3),
  resid_sd_loyo     = round(sd(loyo$Resid, na.rm = TRUE), 3),
  fold_insample     = NA,
  fold_loyo         = NA
)
overall$fold_insample <- round(overall$resid_sd_insample / overall$ribbon_hw_mean, 1)
overall$fold_loyo      <- round(overall$resid_sd_loyo / overall$ribbon_hw_mean, 1)
tab <- rbind(tab, overall)

cat("\n--- Ribbon half-width vs residual SD, by site ---\n")
cat("(fold_loyo = resid_sd_loyo / ribbon_hw_mean -- the honest comparison;\n")
cat(" fold_insample uses fitting residuals, optimistically biased low)\n\n")
print(tab, row.names = FALSE)

write.csv(tab, "4_products/diagnostics/ribbon_vs_residual_sd_table.csv",
          row.names = FALSE)

# Companion plot: fold_loyo by site, headline number at a glance
tab_sites <- tab[tab$Site != "ALL", ]
tab_sites <- tab_sites[match(site_order, tab_sites$Site), ]

pdf("4_products/diagnostics/ribbon_vs_residual_sd.pdf",
    width = 7, height = 5, family = "Helvetica")
par(mar = c(4, 4.5, 2, 1))
bp <- barplot(tab_sites$fold_loyo, names.arg = tab_sites$Site,
              ylab = "LOYO residual SD / ribbon half-width (fold)",
              xlab = "Site", col = "grey70", border = NA)
abline(h = 1, lty = 2, col = "gray40")
text(bp, tab_sites$fold_loyo, labels = tab_sites$fold_loyo, pos = 3, cex = 0.8)
dev.off()
cat("Wrote 4_products/diagnostics/ribbon_vs_residual_sd.pdf\n")