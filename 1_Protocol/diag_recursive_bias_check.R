# ============================================================================
# diag_recursive_bias_check.R
# UCFR Cladophora Bloom Prediction Pipeline -- QAQC diagnostic
# NOT manuscript-facing. Not part of the numbered 1_protocol/ sequence.
#
# Question: does feeding a GAM conditional-mean prediction forward as next
# year's lag_y (13_project_bloom.R's recursion; lag_y is defined in
# 10_bloom_model_M1.R as the prior year's observed ANNUAL MAX) introduce a
# systematic, growing downward bias -- distinct from the already-quantified
# ribbon-width (variance) problem?
#
# Two independent checks, both read-only against existing outputs. Neither
# refits M1 or the TP submodel, so no mgcv dependency here.
#
#   Part 1 -- Signed error by step (bias check)
#     11_temporal_validation.R's Scheme D (recursive-mode validation) already
#     writes recursive_val_predictions.csv with a signed
#     Error = Observed - Predicted_rec column. Its own scorecard only ever
#     aggregates RMSE/MAE (magnitudes), which hides sign. This re-aggregates
#     the same file for mean signed error by step. A positive mean that grows
#     with step count is exactly the direction the mean-vs-max mismatch
#     predicts (observed annual max > fed-forward conditional mean).
#     This does NOT require re-running 11_temporal_validation.R -- it reads
#     the CSV that script already writes.
#
#   Part 2 -- Annual-mean overlay on the projection ribbons
#     The mid-session audit rebuilt the observed overlay on
#     bloom_trajectories.pdf using annual MAX (mitigated but did not close
#     the discontinuity) and proposed, but did not build, an annual MEAN
#     version to isolate the effect. This builds that overlay: both observed
#     annual max and annual mean plotted per site against the low/high
#     projection ribbons.
#     CAVEAT (added after first real run): this panel runs 75 years under
#     drifting climate drivers (anomaly/temp/TP), so any lag_y-scale effect
#     is confounded with genuine driver-trend effects -- it is NOT a clean
#     test of whether the recursion is mean-scale or max-scale. Part 1 is
#     the clean test (validation years stay in the observed driver range).
#     Treat Part 2 as a visual sanity check only (e.g. "is there still a
#     visible discontinuity at the seed year") -- do not read a mean-vs-max
#     scale conclusion off it.
#
# Inputs
#   2_incremental/recursive_val_predictions.csv   (11_temporal_validation.R, Scheme D)
#   2_incremental/ucfr_model_ready.csv            (observed visit-level data)
#   2_incremental/bloom_projections.csv           (13_project_bloom.R ensemble summary)
#
# Outputs
#   4_products/diagnostics/recursive_signed_error_by_step.csv
#   4_products/diagnostics/recursive_bias_check.pdf
#   console: scorecard data frame
# ============================================================================

# ============================================================================
# CONFIGURATION -- edit here only
# ============================================================================
PATH_REC   <- "2_incremental/recursive_val_predictions.csv"
PATH_OBS   <- "2_incremental/ucfr_model_ready.csv"
PATH_PROJ  <- "2_incremental/bloom_projections.csv"

OUT_DIR    <- "4_products/diagnostics"
OUT_CSV    <- file.path(OUT_DIR, "recursive_signed_error_by_step.csv")
OUT_PDF    <- file.path(OUT_DIR, "recursive_bias_check.pdf")

site_order <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")   # upstream -> downstream

col_low   <- "#377EB8"                        # blue = low (RCP4.5/SSP245)
col_high  <- "#E41A1C"                        # red  = high (RCP8.5/SSP585)
fill_low  <- adjustcolor(col_low,  alpha.f = 0.22)
fill_high <- adjustcolor(col_high, alpha.f = 0.22)
# ============================================================================

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

require_cols <- function(df, cols, what) {
  miss <- setdiff(cols, names(df))
  if (length(miss))
    stop(sprintf("%s is missing column(s): %s\n  Available: %s",
                 what, paste(miss, collapse = ", "),
                 paste(names(df), collapse = ", ")), call. = FALSE)
}

# ============================================================================
# PART 1 -- SIGNED ERROR BY STEP
# ============================================================================
rec <- read.csv(PATH_REC, stringsAsFactors = FALSE)
require_cols(rec, c("Site", "Year", "step", "Observed", "Predicted_rec", "Error"),
             "recursive_val_predictions.csv")

max_step <- max(rec$step, na.rm = TRUE)
err_by_step <- do.call(rbind, lapply(seq_len(max_step), function(k) {
  d <- rec[rec$step == k & is.finite(rec$Error), ]
  n <- nrow(d)
  data.frame(
    step    = k,
    n       = n,
    MeanErr = round(mean(d$Error), 4),
    SD_Err  = round(if (n > 1) sd(d$Error) else NA_real_, 4),
    SE_Err  = round(if (n > 1) sd(d$Error) / sqrt(n) else NA_real_, 4),
    RMSE    = round(sqrt(mean(d$Error^2)), 4)
  )
}))

# Descriptive linear trend of MeanErr vs step -- not a formal inferential
# model (n = number of steps, usually small), just a slope for the scorecard.
trend <- if (nrow(err_by_step) >= 3 && var(err_by_step$step) > 0) {
  tryCatch(lm(MeanErr ~ step, data = err_by_step), error = function(e) NULL)
} else NULL
trend_slope <- if (!is.null(trend)) round(coef(trend)[["step"]], 5) else NA_real_

write.csv(err_by_step, OUT_CSV, row.names = FALSE)


# ============================================================================
# PART 2 -- ANNUAL-MEAN vs ANNUAL-MAX OVERLAY ON PROJECTION RIBBONS
# ============================================================================
obs  <- read.csv(PATH_OBS,  stringsAsFactors = FALSE)
proj <- read.csv(PATH_PROJ, stringsAsFactors = FALSE)

require_cols(obs,  c("Site", "Year", "logCHLa"), "ucfr_model_ready.csv")
require_cols(proj, c("Site", "bracket", "year", "median_logCHLa",
                     "p10_logCHLa", "p90_logCHLa"), "bloom_projections.csv")

obs_clean <- obs[is.finite(obs$logCHLa), ]

ann_max  <- aggregate(logCHLa ~ Site + Year, data = obs_clean, FUN = max)
names(ann_max)[3] <- "ann_max"
ann_mean <- aggregate(logCHLa ~ Site + Year, data = obs_clean, FUN = mean)
names(ann_mean)[3] <- "ann_mean"
ann <- merge(ann_max, ann_mean, by = c("Site", "Year"))

xlim <- range(c(proj$year, ann$Year), na.rm = TRUE)
ylim <- range(c(proj$p10_logCHLa, proj$p90_logCHLa,
                ann$ann_max, ann$ann_mean), na.rm = TRUE)
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

pdf(OUT_PDF, width = 10, height = 11)

# --- Page 1: signed error by step (bias check) ---
par(mfrow = c(1, 1), mar = c(4.2, 4.2, 2, 1))
has_se <- any(is.finite(err_by_step$SE_Err))
ylim_e <- range(c(err_by_step$MeanErr,
                  if (has_se) err_by_step$MeanErr - err_by_step$SE_Err else NULL,
                  if (has_se) err_by_step$MeanErr + err_by_step$SE_Err else NULL,
                  0), na.rm = TRUE)
ylim_e <- ylim_e + c(-0.05, 0.05) * diff(ylim_e)

plot(err_by_step$step, err_by_step$MeanErr, type = "b", pch = 16,
     col = "darkorange", lwd = 2, ylim = ylim_e,
     xlab = "Steps ahead (years)",
     ylab = "Mean signed error, Observed minus Predicted (recursive)")
if (has_se) {
  segments(err_by_step$step, err_by_step$MeanErr - err_by_step$SE_Err,
           err_by_step$step, err_by_step$MeanErr + err_by_step$SE_Err,
           col = "darkorange")
}
abline(h = 0, lty = 2, col = "grey40")
if (!is.null(trend)) {
  abline(trend, lty = 3, col = "steelblue")
  legend("topleft", legend = paste0("linear trend: ", trend_slope, " / step"),
         bty = "n", text.col = "steelblue", cex = 0.85)
}
text(err_by_step$step, err_by_step$MeanErr, paste0("n=", err_by_step$n),
     pos = 3, cex = 0.65, col = "grey40")
mtext("positive = model underpredicts (annual-max-vs-conditional-mean direction)",
      side = 3, line = 0.3, cex = 0.8, col = "grey30")

# --- Pages 2+: per-site overlay ---
par(mfrow = c(4, 2), mar = c(3.2, 3.6, 2.2, 0.8),
    mgp = c(2.1, 0.6, 0), oma = c(2.2, 1.5, 2.4, 0.5))

for (st in site_order) {
  ds <- proj[proj$Site == st, ]
  lo <- ds[ds$bracket == "low", ]
  hi <- ds[ds$bracket == "high", ]
  da <- ann[ann$Site == st, ]
  da <- da[order(da$Year), ]
  
  plot(NA, xlim = xlim, ylim = ylim,
       xlab = "Year", ylab = expression(log[10] * " Chl " * italic(a)), main = st)
  
  if (nrow(lo)) { ribbon(lo, fill_low);  med_line(lo, col_low) }
  if (nrow(hi)) { ribbon(hi, fill_high); med_line(hi, col_high) }
  
  if (nrow(da)) {
    lines(da$Year, da$ann_max,  col = "black",  lty = 1, lwd = 1)
    points(da$Year, da$ann_max, col = "black",  pch = 17, cex = 0.8)
    lines(da$Year, da$ann_mean,  col = "grey40", lty = 2, lwd = 1)
    points(da$Year, da$ann_mean, col = "grey40", pch = 16, cex = 0.8)
  }
  box()
}

plot.new()
leg  <- c("Low median (RCP4.5/SSP245)", "High median (RCP8.5/SSP585)",
          "10th-90th percentile", "Observed annual max", "Observed annual mean")
lcol <- c(col_low, col_high, "grey50", "black", "grey40")
llwd <- c(2, 2, NA, 1, 1)
lpch <- c(NA, NA, 15, 17, 16)
llty <- c(1, 1, NA, 1, 2)
legend("center", legend = leg, col = lcol, lwd = llwd, pch = lpch, lty = llty,
       pt.cex = c(2.2, 2.2, 2.2, 0.9, 0.9), bty = "n", cex = 0.9)

mtext("QAQC: annual max/mean vs. projection ribbons (visual reference only)",
      outer = TRUE, cex = 1.0, font = 2, line = 0.6)
mtext("not a test of recursion scale -- conflates lag_y effects with 75yr driver drift; see Page 1 / signed-error CSV for the actual test",
      outer = TRUE, side = 1, cex = 0.65, line = 0.6)

dev.off()


# ============================================================================
# SCORECARD
# ============================================================================
cat("\n")
cat("============================================================\n")
cat("QAQC -- RECURSIVE-MODE BIAS CHECK\n")
cat("Not manuscript-facing.\n")
cat("============================================================\n\n")

cat("--- Signed error by step (Observed - Predicted_rec) ---\n")
print(err_by_step, row.names = FALSE)
cat("\n")
cat("Linear trend of MeanErr vs step:",
    ifelse(is.na(trend_slope), "not estimable (< 3 steps or no step variance)",
           paste0(trend_slope, " per step")), "\n")
cat("Interpretation: a trend near zero is consistent with no growing bias.\n")
cat("A positive, growing trend is consistent with the annual-max-vs-\n")
cat("conditional-mean mismatch described in the audit -- distinct from,\n")
cat("and additive to, the already-quantified ribbon-width (variance) gap.\n\n")

sc_files <- data.frame(
  File = c(OUT_CSV, OUT_PDF),
  Contents = c(
    "Signed mean error, SD/SE, RMSE by recursive step.",
    "Page 1: mean signed error by step. Pages 2+: per-site annual max/mean overlay on projection ribbons."
  ),
  stringsAsFactors = FALSE
)
cat("--- Output files ---\n")
print(sc_files, row.names = FALSE)
cat("\n")
cat("============================================================\n")
cat("Done.\n")
cat("============================================================\n")