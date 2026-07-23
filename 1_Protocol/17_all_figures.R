# ============================================================================
# 17_all_figures.R
# UCFR Cladophora Bloom Prediction Pipeline
# ----------------------------------------------------------------------------
# ONE script that generates every manuscript figure, in canonical order.
#
# Design intent (per Rafa, 2026-07): a single place where all figures live so
# graphical parameters can be adjusted in one file, and any edit is visible
# across the figure set. Redundancy with the source scripts (15, 15b, 13c,
# diag_*) is INTENTIONAL and accepted for now -- housekeeping/consolidation
# (retiring the source scripts, moving the exceedance computation into a
# numbered 13d) is deferred. Each section is self-contained: it reads its own
# inputs, sets its own colours/params, opens its own device, and closes it, so
# any single section can be run on its own by selecting and executing it.
#
# Canonical figure order (final numbers are authoritative in the .Rmd
# captions, NOT here -- section labels below are organisational only):
#
#   Sec 0  Methods  Site/reach map ............... PLACEHOLDER (not built)
#   Sec 1  3.1 Fig  GAM smooth-term panel ........ from 15 Section A
#   Sec 2  3.2 Fig  Bloom trajectories ........... from 15b (canonical: superset of 13b)
#   Sec 3  3.3 Fig  Exceedance probability ....... HEADLINE; computed inline (ex diag_exceedance)
#   Sec 4  3.4 Fig  Driver decomposition ......... from 13c
#   Sec 5  3.4 Fig  Hydrology driver drift ....... from 15 Section C  [main-vs-SI placement OPEN]
#   Sec 6  SI  Fig  Recursive-mode validation .... from diag_recursive_bias_check
#   Sec 7  SI  Fig  Temperature projections ...... from 15 Section D
#   Sec 8  SI  Fig  Ensemble variance by year .... from 15 Section E
#
# Inputs (all pre-existing pipeline intermediates -- nothing new required):
#   3_models/bloom_model_M1.rds                     (Sec 1, 4)
#   2_incremental/bloom_projections.csv             (Sec 2, 6)
#   2_incremental/ucfr_model_ready.csv              (Sec 2, 6)
#   2_incremental/bloom_projections_members.csv     (Sec 3, 4, 8)
#   2_incremental/recursive_val_predictions.csv     (Sec 3, 6)
#   2_incremental/ncar_discharge_envelope.csv       (Sec 5)
#   2_incremental/ncar_temperature_envelope.csv     (Sec 7)
#
# Outputs -> 4_products/Final_figures/  (one figure per file, for the .Rmd):
#   site_reach_map_PLACEHOLDER.pdf
#   gam_smooth_panel.pdf
#   bloom_trajectories.pdf
#   exceedance_trajectories.pdf            (headline: P(exceed) by site through 2098)
#   exceedance_naive_vs_corrected.pdf      (methods visual: BN ~0 -> ~0.50 story)
#   exceedance_residual_pool.pdf           (transparency: the bootstrap pool)
#   driver_decomposition.pdf
#   hydrology_drift_by_reach.pdf
#   SI_recursive_signed_error.pdf          (the actual bias test)
#   SI_recursive_maxmean_overlay.pdf       (visual reference only -- NOT a test)
#   temp_projections_by_site.pdf
#   ensemble_variance_by_year.pdf
# Also writes (data product feeding the .Rmd threshold table):
#   2_incremental/exceedance_probability.csv
# ============================================================================

library(mgcv)

OUT_DIR <- "4_products/Final_figures"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Shared palette (kept per-section too, so sections stay independent) --------
COL_LOW  <- "#377EB8"   # low  bracket (RCP4.5 / SSP245)
COL_HIGH <- "#E41A1C"   # high bracket (RCP8.5 / SSP585)
SITE_ORDER <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")   # upstream -> downstream

# small helper used by a couple of sections
.require_cols <- function(df, cols, what) {
  miss <- setdiff(cols, names(df))
  if (length(miss))
    stop(sprintf("%s is missing column(s): %s\n  Available: %s",
                 what, paste(miss, collapse = ", "),
                 paste(names(df), collapse = ", ")), call. = FALSE)
}


# ============================================================================
# SECTION 0 -- Methods: Site/reach map  (PLACEHOLDER, not yet built)
# ----------------------------------------------------------------------------
# Reserves the figure slot for the co-author report. Real build needs site
# lat/lon (not currently an input here) + sf for the base layer, per project
# package conventions (see 15 Section B). Replace this block when scoped.
# ============================================================================
pdf(file.path(OUT_DIR, "1_site_reach_map_PLACEHOLDER.pdf"),
    width = 7, height = 5, family = "Helvetica")
par(mar = c(1, 1, 1, 1))
plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
box(col = "grey60")
text(0.5, 0.62, "Site / reach map", cex = 1.7, font = 2)
text(0.5, 0.46, "PLACEHOLDER -- not yet built", cex = 1.1, col = "grey40")
text(0.5, 0.30, "needs site lat/lon + sf base layer (see 15 Section B)",
     cex = 0.8, col = "grey55")
dev.off()
cat("Wrote", file.path(OUT_DIR, "1_site_reach_map_PLACEHOLDER.pdf"), "\n")


# ============================================================================
# SECTION 1 -- 3.1  Figure: GAM smooth-term panel  (from 15 Section A)
# ----------------------------------------------------------------------------
# All six M1 partial effects. Panel order mirrors the M1 formula (not an
# influence ranking). select = 1:6 excludes s(Site, bs="re") (7th smooth).
# No per-panel main= -- mgcv's default y-axis label "s(term, edf)" identifies
# each panel; a title would repeat it. seWithMean = TRUE: CI bands include
# uncertainty about the overall mean (flip to FALSE for shape-only bands).
# ============================================================================
s1_m1 <- readRDS("3_models/bloom_model_M1.rds")

s1_xlab <- c(
  "Lag-year max log10(CHLa)",
  "(Peak Q / baseflow Q)^(1/3)",
  "log10(Q, cfs)",
  "Days since freshet",
  "log10(TP, mg/L)",
  "Temperature (deg C)"
)

pdf(file.path(OUT_DIR, "2_gam_smooth_panel.pdf"),
    width = 9, height = 6, family = "Helvetica")
par(mfrow = c(2, 3), mar = c(4, 4.5, 2, 1))
for (i in 1:6) {
  plot(s1_m1,
       select     = i,
       residuals  = TRUE,
       shade      = TRUE,
       shade.col  = "gray85",
       seWithMean = TRUE,
       pch        = 1,
       cex        = 0.6,
       col        = "steelblue4",
       lwd        = 2,
       xlab       = s1_xlab[i])
  abline(h = 0, lty = 2, col = "gray60")
}
dev.off()
cat("Wrote", file.path(OUT_DIR, "2_gam_smooth_panel.pdf"), "\n")


# ============================================================================
# SECTION 2 -- 3.2  Figure: Bloom trajectories  (from 15b, canonical)
# ----------------------------------------------------------------------------
# Per-site low vs high, median line + 10th-90th ribbon (across members).
# Observed historical points + lowess smooth REMOVED (2026-07, per Rafa) --
# projections-only view. Panels upstream -> downstream, shared axes.
# Per-panel main=st title and legend REMOVED (2026-07, per Rafa) --
# publication figure, non-panel elements to be handled via caption.
#
# Panel H (2026-07, per Rafa): general diagnostic slot, filled with pooled
# recursive-mode validation -- mean signed error by step-ahead, all sites
# pooled, from recursive_val_predictions.csv. Same computation as Section 6 /
# SI_recursive_signed_error.pdf, but text annotations (n= labels, trend-slope
# legend) dropped to match this figure's panel-only convention; the fitted
# trend line itself is kept as a plain (unlabeled) graphical element.
#
# Panel letters A-H (2026-07, per Rafa): added top-right of every panel via
# s2_panel_letter(), to make external legend/caption writing easier. A-G =
# SITE_ORDER (upstream -> downstream), H = the diagnostic panel above.
#
# NOTE (not acted on): Section 2 now also reads recursive_val_predictions.csv,
# which isn't yet listed in the master input list at the top of this script
# (line ~34) -- worth adding there when the whole script is next touched.
# ============================================================================
s2_s <- read.csv("2_incremental/bloom_projections.csv", stringsAsFactors = FALSE)

s2_fill_low  <- adjustcolor(COL_LOW,  alpha.f = 0.22)
s2_fill_high <- adjustcolor(COL_HIGH, alpha.f = 0.22)

s2_xlim <- range(s2_s$year, na.rm = TRUE)
s2_ylim <- range(c(s2_s$p10_logCHLa, s2_s$p90_logCHLa), na.rm = TRUE)
s2_ylim <- s2_ylim + c(-0.03, 0.03) * diff(s2_ylim)

s2_ribbon <- function(d, fill) {
  d <- d[order(d$year), ]
  polygon(c(d$year, rev(d$year)), c(d$p10_logCHLa, rev(d$p90_logCHLa)),
          col = fill, border = NA)
}
s2_med <- function(d, col) {
  d <- d[order(d$year), ]
  lines(d$year, d$median_logCHLa, col = col, lwd = 2)
}
s2_panel_letter <- function(letter) {
  usr <- par("usr")
  text(usr[2], usr[4], labels = letter, xpd = NA,
       adj = c(1.5, 1.5), font = 2, cex = 1.1)
}
s2_letters <- LETTERS[1:8]

pdf(file.path(OUT_DIR, "3_bloom_trajectories.pdf"), width = 9, height = 11)
par(mfrow = c(4, 2), mar = c(3.2, 3.6, 2.2, 0.8),
    mgp = c(2.1, 0.6, 0), oma = c(2.2, 1.5, 2.4, 0.5))

for (i in seq_along(SITE_ORDER)) {
  st <- SITE_ORDER[i]
  ds <- s2_s[s2_s$Site == st, ]
  lo <- ds[ds$bracket == "low",  ]
  hi <- ds[ds$bracket == "high", ]
  plot(NA, xlim = s2_xlim, ylim = s2_ylim,
       xlab = "Year",
       ylab = expression(log[10] * " Chl " * italic(a) * " (mg/m"^2 * ")"),
       main = st)
  s2_ribbon(lo, s2_fill_low);  s2_ribbon(hi, s2_fill_high)
  s2_med(lo, COL_LOW); s2_med(hi, COL_HIGH)
  box(lwd=2)
  s2_panel_letter(s2_letters[i])
}

# --- panel H: pooled recursive-mode validation (general diagnostic) --------
s2_rec <- read.csv("2_incremental/recursive_val_predictions.csv", stringsAsFactors = FALSE)
.require_cols(s2_rec, c("Site", "step", "Error"), "recursive_val_predictions.csv")

s2_maxstep <- max(s2_rec$step, na.rm = TRUE)
s2_ebs <- do.call(rbind, lapply(seq_len(s2_maxstep), function(k) {
  d <- s2_rec[s2_rec$step == k & is.finite(s2_rec$Error), ]
  n <- nrow(d)
  data.frame(step    = k,
             MeanErr = mean(d$Error),
             SE_Err  = if (n > 1) sd(d$Error) / sqrt(n) else NA_real_)
}))

s2_h_has_se <- any(is.finite(s2_ebs$SE_Err))
s2_h_ylim <- range(c(s2_ebs$MeanErr,
                     if (s2_h_has_se) s2_ebs$MeanErr - s2_ebs$SE_Err else NULL,
                     if (s2_h_has_se) s2_ebs$MeanErr + s2_ebs$SE_Err else NULL, 0),
                   na.rm = TRUE)
s2_h_ylim <- s2_h_ylim + c(-0.05, 0.05) * diff(s2_h_ylim)

s2_h_trend <- if (nrow(s2_ebs) >= 3 && var(s2_ebs$step) > 0) {
  tryCatch(lm(MeanErr ~ step, data = s2_ebs), error = function(e) NULL)
} else NULL

plot(s2_ebs$step, s2_ebs$MeanErr, type = "b", pch = 16, col = "darkorange", lwd = 2,
     xlim = c(1, s2_maxstep), ylim = s2_h_ylim,
     xlab = "Steps ahead (years)",
     ylab = expression("Mean signed error, " * log[10] * " Chl " * italic(a) * " (mg/m"^2 * ")"))
if (s2_h_has_se)
  segments(s2_ebs$step, s2_ebs$MeanErr - s2_ebs$SE_Err,
           s2_ebs$step, s2_ebs$MeanErr + s2_ebs$SE_Err, col = "darkorange")
abline(h = 0, lty = 2, col = "grey40")
if (!is.null(s2_h_trend)) abline(s2_h_trend, lty = 3, col = "steelblue")
box(lwd=2)
s2_panel_letter(s2_letters[8])

mtext("Projected Cladophora bloom trajectories (NCAR ensemble)",
      outer = TRUE, cex = 1.05, font = 2, line = 0.6)
mtext("panels upstream -> downstream  |  shaded band = 10th-90th pct",
      outer = TRUE, side = 1, cex = 0.75, line = 0.6)
dev.off()
cat("Wrote", file.path(OUT_DIR, "3_bloom_trajectories.pdf"), "\n")

# ============================================================================
# SECTION 3 -- 3.3  HEADLINE: Exceedance probability  (computation ex diag_exceedance)
# ----------------------------------------------------------------------------
# P(CHLa > 150 mg/m^2) (Suplee et al. 2009, JAWRA) under honest error
# propagation ("C+"): NCAR ensemble members (outer) x bootstrap resample of
# Scheme D recursive signed errors (inner). The recursive Error pool is
# already annual-max-calibrated and is real out-of-sample predictive error, so
# it satisfies both audit conditions (full predictive uncertainty; annual-max
# target) by construction. NAIVE (ensemble spread only) is reported alongside
# CORRECTED to show what carrying real uncertainty changes.
#
# NOTE: computation duplicated from diag_exceedance_probability_check.R for
# self-containment. Same SEED / draws -> reproduces its numbers exactly. When
# housekeeping promotes this to a numbered 13d, move the computation + CSV
# there and have this section read the CSV instead.
#
# LIMITATION (state if this goes to Methods): the residual pool is drawn from
# validation years inside the observed driver envelope; temp_extrapolated
# projection years extend beyond it, so corrected P there is a floor, not a
# full accounting -- flagged, not fixed (consistent with the HU divergence).
# ============================================================================
S3_THRESHOLD_MG_M2 <- 150
S3_THRESHOLD_LOG10 <- log10(S3_THRESHOLD_MG_M2)
S3_B_DRAWS         <- 2000          # residual bootstrap draws per member
S3_POOL_ACROSS     <- TRUE          # TRUE: pooled ~90-row pool (primary). FALSE: per-site (sensitivity only)
S3_SEED            <- 20260701
S3_BENCH_YEAR      <- NULL          # NULL = max(year) present (end-century headline)

set.seed(S3_SEED)

s3_mem <- read.csv("2_incremental/bloom_projections_members.csv", stringsAsFactors = FALSE)
s3_rec <- read.csv("2_incremental/recursive_val_predictions.csv", stringsAsFactors = FALSE)
.require_cols(s3_mem, c("Site", "bracket", "year", "pred_logCHLa"),
              "bloom_projections_members.csv")
.require_cols(s3_rec, c("Site", "Error"), "recursive_val_predictions.csv")
if (!"temp_extrapolated" %in% names(s3_mem)) {
  warning("temp_extrapolated not found in members file -- extrapolation flag will be NA.")
  s3_mem$temp_extrapolated <- NA
}

# residual bootstrap pools
s3_resid_all <- s3_rec$Error[is.finite(s3_rec$Error)]
if (length(s3_resid_all) < 20)
  warning("Pooled residual pool has only ", length(s3_resid_all),
          " values -- bootstrap tails will be unstable.")
s3_resid_by_site <- split(s3_rec$Error[is.finite(s3_rec$Error)],
                          s3_rec$Site[is.finite(s3_rec$Error)])
s3_get_pool <- function(site) {
  if (S3_POOL_ACROSS) return(s3_resid_all)
  p <- s3_resid_by_site[[site]]
  if (is.null(p) || length(p) < 5) {
    warning("Site ", site, " has < 5 residuals -- falling back to pooled.")
    return(s3_resid_all)
  }
  p
}

# P(exceed) per Site x bracket x year
s3_groups <- unique(s3_mem[, c("Site", "bracket", "year")])
s3_groups <- s3_groups[order(s3_groups$Site, s3_groups$bracket, s3_groups$year), ]
s3_out <- data.frame(
  Site = s3_groups$Site, bracket = s3_groups$bracket, year = s3_groups$year,
  n_members = NA_integer_, n_draws = NA_integer_,
  P_exceed_naive = NA_real_, P_exceed_corrected = NA_real_,
  frac_temp_extrap = NA_real_
)
for (i in seq_len(nrow(s3_groups))) {
  s <- s3_groups$Site[i]; b <- s3_groups$bracket[i]; y <- s3_groups$year[i]
  d <- s3_mem[s3_mem$Site == s & s3_mem$bracket == b & s3_mem$year == y, ]
  if (!nrow(d)) next
  preds <- d$pred_logCHLa[is.finite(d$pred_logCHLa)]
  n_mem <- length(preds)
  if (!n_mem) next
  s3_out$P_exceed_naive[i] <- mean(preds > S3_THRESHOLD_LOG10)   # ensemble spread only
  pool <- s3_get_pool(s)
  resid_draws <- matrix(sample(pool, n_mem * S3_B_DRAWS, replace = TRUE),
                        nrow = n_mem, ncol = S3_B_DRAWS)
  composite <- preds + resid_draws       # composite[i,j] = preds[i] + draw[i,j]
  s3_out$n_members[i] <- n_mem
  s3_out$n_draws[i]   <- length(composite)
  s3_out$P_exceed_corrected[i] <- mean(composite > S3_THRESHOLD_LOG10)
  if (!all(is.na(d$temp_extrapolated)))
    s3_out$frac_temp_extrap[i] <- mean(d$temp_extrapolated, na.rm = TRUE)
}
write.csv(s3_out, "2_incremental/exceedance_probability.csv", row.names = FALSE)
cat("Wrote 2_incremental/exceedance_probability.csv\n")

s3_bench_year <- if (is.null(S3_BENCH_YEAR)) max(s3_out$year, na.rm = TRUE) else S3_BENCH_YEAR
s3_bench <- s3_out[s3_out$year == s3_bench_year, ]

# --- 3a. HEADLINE figure: per-site P(exceed) through 2098 --------------------
pdf(file.path(OUT_DIR, "4_exceedance_trajectories.pdf"), width = 9, height = 11)
par(mfrow = c(4, 2), mar = c(3.2, 3.6, 2.2, 0.8),
    mgp = c(2.1, 0.6, 0), oma = c(2.2, 1.5, 2.4, 0.5))
s3_yrs <- sort(unique(s3_out$year))
for (st in SITE_ORDER) {
  d  <- s3_out[s3_out$Site == st, ]
  lo <- d[d$bracket == "low",  ]; lo <- lo[order(lo$year), ]
  hi <- d[d$bracket == "high", ]; hi <- hi[order(hi$year), ]
  plot(NA, xlim = range(s3_yrs), ylim = c(0, 1),
       xlab = "Year", ylab = "P(exceed)", main = st)
  extrap_yrs <- d$year[!is.na(d$frac_temp_extrap) & d$frac_temp_extrap > 0.5]
  if (length(extrap_yrs))
    rect(min(extrap_yrs), -0.05, max(s3_yrs) + 1, 1.05,
         col = adjustcolor("grey50", alpha.f = 0.12), border = NA)
  if (nrow(lo)) {
    lines(lo$year, lo$P_exceed_naive,     col = COL_LOW,  lty = 2, lwd = 1)
    lines(lo$year, lo$P_exceed_corrected, col = COL_LOW,  lty = 1, lwd = 2)
  }
  if (nrow(hi)) {
    lines(hi$year, hi$P_exceed_naive,     col = COL_HIGH, lty = 2, lwd = 1)
    lines(hi$year, hi$P_exceed_corrected, col = COL_HIGH, lty = 1, lwd = 2)
  }
  box()
}
plot.new()
legend("center",
       legend = c("Low, corrected", "Low, naive", "High, corrected", "High, naive",
                  "temp-extrapolated (majority)"),
       col = c(COL_LOW, COL_LOW, COL_HIGH, COL_HIGH, "grey50"),
       lty = c(1, 2, 1, 2, NA), lwd = c(2, 1, 2, 1, NA), pch = c(NA, NA, NA, NA, 15),
       pt.cex = 2, bty = "n", cex = 0.9)
mtext("Projected exceedance probability, P(CHLa > 150 mg/m^2)",
      outer = TRUE, cex = 1.0, font = 2, line = 0.6)
mtext("thick = corrected (+ residual bootstrap) | thin dashed = naive (ensemble only) | Suplee et al. 2009",
      outer = TRUE, side = 1, cex = 0.7, line = 0.6)
dev.off()
cat("Wrote", file.path(OUT_DIR, "4_exceedance_trajectories.pdf"), "\n")

# --- 3a. HEADLINE figure: per-site P(exceed) through 2098 --------------------
# Temp-extrapolated shaded region REMOVED (2026-07, per Rafa) -- frac_temp_extrap
# is still computed and written to exceedance_probability.csv upstream (feeds
# the SI extrapolation table); only this plot's shading was cut.
#
# Naive lines REMOVED (2026-07, per Rafa) -- headline now shows corrected
# P(exceed) only. The naive-vs-corrected comparison (why the residual
# bootstrap can't be skipped) stays in exceedance_naive_vs_corrected.pdf
# (Section 3b) as a methods-justification / supplementary figure.
#
# Per-panel legend REMOVED, panel letters A-H added -- consistent with
# bloom_trajectories.pdf (Section 2). main = st kept for site identification,
# same as Section 2's final state.
#
# Panel H: general diagnostic slot, filled with the pooled residual bootstrap
# pool -- the same error distribution behind every corrected P(exceed) value
# in panels A-G. Analogous role to Section 2's recursive-validation panel H:
# it's the empirical basis for trusting the numbers next to it, not a new
# result. Adapted from Section 3c / exceedance_residual_pool.pdf; the n=/
# mean/SD legend text is dropped to match the panel-only convention, but the
# zero and mean reference lines stay as plain graphical elements.
#
# Bottom caption edited (not just left stale) to drop "thin dashed = naive",
# since that line style no longer appears in the plot.
# ============================================================================
s3_panel_letter <- function(letter) {
  usr <- par("usr")
  text(usr[2], usr[4], labels = letter, xpd = NA,
       adj = c(1.5, 1.5), font = 2, cex = 1.1)
}
s3_letters <- LETTERS[1:8]

pdf(file.path(OUT_DIR, "4_exceedance_trajectories.pdf"), width = 9, height = 11)
par(mfrow = c(4, 2), mar = c(3.2, 3.6, 2.2, 0.8),
    mgp = c(2.1, 0.6, 0), oma = c(2.2, 1.5, 2.4, 0.5))
s3_yrs <- sort(unique(s3_out$year))

for (i in seq_along(SITE_ORDER)) {
  st <- SITE_ORDER[i]
  d  <- s3_out[s3_out$Site == st, ]
  lo <- d[d$bracket == "low",  ]; lo <- lo[order(lo$year), ]
  hi <- d[d$bracket == "high", ]; hi <- hi[order(hi$year), ]
  plot(NA, xlim = range(s3_yrs), ylim = c(0, 1),
       xlab = "Year", ylab = "P(exceed)", main = st)
  if (nrow(lo)) lines(lo$year, lo$P_exceed_corrected, col = COL_LOW,  lwd = 2)
  if (nrow(hi)) lines(hi$year, hi$P_exceed_corrected, col = COL_HIGH, lwd = 2)
  box()
  s3_panel_letter(s3_letters[i])
}

# --- panel H: pooled residual bootstrap pool (general diagnostic) ----------
hist(s3_resid_all, breaks = 20, col = "grey70", border = "white", main = "",
     xlab = "Historical recursive error, Observed annual max - Predicted (log10 CHLa)",
     ylab = "Frequency")
abline(v = 0, lty = 2, col = "grey30")
abline(v = mean(s3_resid_all), lty = 2, col = "darkorange", lwd = 2)
box()
s3_panel_letter(s3_letters[8])

mtext("Projected exceedance probability, P(CHLa > 150 mg/m^2)",
      outer = TRUE, cex = 1.0, font = 2, line = 0.6)
mtext("corrected (+ residual bootstrap) | Suplee et al. 2009",
      outer = TRUE, side = 1, cex = 0.7, line = 0.6)
dev.off()
cat("Wrote", file.path(OUT_DIR, "4_exceedance_trajectories.pdf"), "\n")

# ============================================================================
# SECTION 4 -- 3.4  Figure: Driver decomposition  (from 13c)
# ----------------------------------------------------------------------------
# M1 is additive, so predict(type="terms") partitions logCHLa exactly into
# per-smooth contributions. Each term's change from early (2026-2035) to late
# (2089-2098) window, as share of |total change|. lag_y = persistence/
# inherited, not a climate forcing.
#
# BOTH SCENARIOS (2026-07, per Rafa) -- previously high bracket only, now two
# stacked panels (high on top, low below), same % share method for each.
# Cross-checked against a freeze-one-driver counterfactual re-projection
# (diag_freeze_one_driver.R): top-ranked driver agreed at all 7 sites under
# both scenarios, magnitudes matched within ~0.5 pct points where compared --
# this static method isn't misattributing shared (concurvity) variance.
# ============================================================================
s4_m1   <- readRDS("3_models/bloom_model_M1.rds")
s4_grid <- read.csv("2_incremental/bloom_projections_members.csv", stringsAsFactors = FALSE)
s4_grid$Site <- factor(s4_grid$Site, levels = levels(s4_m1$model$Site))

s4_tt  <- predict(s4_m1, newdata = s4_grid, type = "terms")
s4_drv <- c("s(lag_y)", "s(anomaly)", "s(logQ_obs_cfs)",
            "s(Days_Since_Freshet)", "s(logTP_mg_L)", "s(Temp_oC)")
s4_tt <- s4_tt[, s4_drv]
colnames(s4_tt) <- c("lag", "anomaly", "logQ", "DSF", "logTP", "Temp")

s4_g <- data.frame(Site = s4_grid$Site, bracket = s4_grid$bracket,
                   year = s4_grid$year, s4_tt)

s4_pct_by_bracket <- function(bracket) {
  g     <- s4_g[s4_g$bracket == bracket, ]
  early <- g[g$year <= 2035, ]
  late  <- g[g$year >= 2089, ]
  delta <- sapply(SITE_ORDER, function(st)
    colMeans(late[late$Site == st, 4:9]) - colMeans(early[early$Site == st, 4:9]))
  apply(abs(delta), 2, function(x) 100 * x / sum(x))   # drivers x sites, cols sum to 100
}
s4_pct_high <- s4_pct_by_bracket("high")
s4_pct_low  <- s4_pct_by_bracket("low")

s4_cols <- c(lag = "grey75", anomaly = "#1b9e77", logQ = "#377EB8",
             DSF = "#7570b3", logTP = "#d95f02", Temp = "#e7298a")

pdf(file.path(OUT_DIR, "5_driver_decomposition.pdf"), width = 8, height = 9.5)
par(mfrow = c(2, 1), mar = c(4, 4.2, 3, 6.5), xpd = TRUE)

barplot(s4_pct_high, col = s4_cols[rownames(s4_pct_high)], border = "white", las = 1,
        ylab = "% of projected end-century change (|change| share)",
        main = "High scenario (RCP8.5/SSP585)")
legend(x = ncol(s4_pct_high) * 1.25, y = 100, legend = rownames(s4_pct_high),
       fill = s4_cols[rownames(s4_pct_high)], bty = "n", cex = 0.9)

barplot(s4_pct_low, col = s4_cols[rownames(s4_pct_low)], border = "white", las = 1,
        ylab = "% of projected end-century change (|change| share)",
        main = "Low scenario (RCP4.5/SSP245)")
legend(x = ncol(s4_pct_low) * 1.25, y = 100, legend = rownames(s4_pct_low),
       fill = s4_cols[rownames(s4_pct_low)], bty = "n", cex = 0.9)

dev.off()
cat("Wrote", file.path(OUT_DIR, "5_driver_decomposition.pdf"), "\n")


# ============================================================================
# SECTION 5 -- 3.4  Figure: Hydrology driver drift by reach  (from 15 Section C)
# ----------------------------------------------------------------------------
# PLACEMENT OPEN: main text (memory: "Fig 4, locked") vs 3.4/SI (this figure's
# own note). Included here regardless; resolve number/placement in the .Rmd.
# 4 NCAR reaches, ensemble mean + p10-p90 across all members, no bracket split.
# Temp_oC excluded (site- not reach-level, deterministic per bracket) --
# belongs in 3.4 with the ceiling story.
#
# Panel C REDESIGNED (2026-07, per Rafa): the single overlapping-4-reach
# panel let some reach lines go invisible where two reaches' semi-transparent
# ribbons stacked and the blended fill landed close to a line color. Fixed
# structurally, not by tweaking alpha: C is now 4 stacked single-reach
# sub-panels (one ribbon each, nothing left to overlap with), occupying the
# same total footprint as one of panels A/B via layout(), sharing one x-axis
# at the bottom. Ribbons kept in each sub-panel (per Rafa, 2026-07).
#
# Panel letters A/B/C added. C labels the 4-stack once (top of the CLALO
# sub-panel), not once per reach -- the 4 strips are one conceptual panel.
#
# Separate reach-color legend DROPPED: each reach sub-panel in C now carries
# its own reach-code label (mtext, right margin, colored to match) -- and
# colors are identical across A/B/C, so C's labels already serve as the
# figure's color key. Flag if you'd rather have an explicit legend back --
# e.g. for a reader who scans A/B before reaching C.
# ============================================================================
s5_env <- read.csv("2_incremental/ncar_discharge_envelope.csv", stringsAsFactors = FALSE)
s5_reach_order  <- c("CLALO", "CLADR", "CLABE", "CLAPL")
s5_reach_colors <- c(CLALO = "#08519c", CLADR = "#3182bd",
                     CLABE = "#6baed6", CLAPL = "#bdd7e7")

s5_drift <- function(col) {
  x    <- s5_env[[col]]
  keys <- list(site = s5_env$site, water_year = s5_env$water_year)
  m    <- aggregate(x, by = keys, FUN = mean, na.rm = TRUE)
  p10  <- aggregate(x, by = keys, FUN = function(v) quantile(v, 0.10, na.rm = TRUE))
  p90  <- aggregate(x, by = keys, FUN = function(v) quantile(v, 0.90, na.rm = TRUE))
  names(m)[3] <- "mean"; names(p10)[3] <- "p10"; names(p90)[3] <- "p90"
  merge(merge(m, p10, by = c("site", "water_year")), p90, by = c("site", "water_year"))
}

s5_panel_letter <- function(letter) {
  usr <- par("usr")
  text(usr[2], usr[4], labels = letter, xpd = NA,
       adj = c(1.5, 1.5), font = 2, cex = 1.1)
}

# --- panels A, B: unchanged plotting logic, letter added ------------------
s5_panel <- function(d, ylab, log = "", letter) {
  xr <- range(d$water_year, na.rm = TRUE)
  yr <- range(c(d$p10, d$p90), na.rm = TRUE)
  plot(NA, xlim = xr, ylim = yr, xlab = "Water year", ylab = ylab, log = log)
  for (r in s5_reach_order) {
    dd <- d[d$site == r, ]; dd <- dd[order(dd$water_year), ]
    polygon(c(dd$water_year, rev(dd$water_year)), c(dd$p90, rev(dd$p10)),
            col = adjustcolor(s5_reach_colors[r], alpha.f = 0.15), border = NA)
  }
  for (r in s5_reach_order) {
    dd <- d[d$site == r, ]; dd <- dd[order(dd$water_year), ]
    lines(dd$water_year, dd$mean, col = s5_reach_colors[r], lwd = 2)
  }
  s5_panel_letter(letter)
}

# --- panel C: 4 stacked single-reach sub-panels, one x-axis shared at bottom
s5_reach_stack <- function(d, ylab_title) {
  xr <- range(d$water_year, na.rm = TRUE)
  n  <- length(s5_reach_order)
  for (i in seq_along(s5_reach_order)) {
    r  <- s5_reach_order[i]
    dd <- d[d$site == r, ]; dd <- dd[order(dd$water_year), ]
    yr <- range(c(dd$p10, dd$p90), na.rm = TRUE)
    is_top    <- i == 1
    is_bottom <- i == n
    par(mar = c(if (is_bottom) 4   else 0.6,
                4.5,
                if (is_top)    2   else 0.6,
                1.8))
    plot(NA, xlim = xr, ylim = yr, xlab = "", ylab = "",
         xaxt = if (is_bottom) "s" else "n")
    polygon(c(dd$water_year, rev(dd$water_year)), c(dd$p90, rev(dd$p10)),
            col = adjustcolor(s5_reach_colors[r], alpha.f = 0.15), border = NA)
    lines(dd$water_year, dd$mean, col = s5_reach_colors[r], lwd = 2)
    mtext(r, side = 4, line = 0.3, cex = 0.65, col = s5_reach_colors[r], font = 2)
    if (i == 2) mtext(ylab_title, side = 2, line = 3, cex = 0.7)
    if (is_bottom) mtext("Water year", side = 1, line = 2.4, cex = 0.85)
    if (is_top) s5_panel_letter("C")
  }
}

pdf(file.path(OUT_DIR, "6_hydrology_drift_by_reach.pdf"),
    width = 12, height = 6, family = "Helvetica")
layout(matrix(c(1, 2, 3,
                1, 2, 4,
                1, 2, 5,
                1, 2, 6), nrow = 4, byrow = TRUE),
       widths = c(1, 1, 1))

par(mar = c(4, 4.5, 2, 1))
s5_panel(s5_drift("anomaly_ma20"),
         "Anomaly, (peak Q / baseflow Q)^(1/3), MA20", letter = "A")
s5_panel(s5_drift("mean_q_cfs_ma20"),
         "Mean discharge, MA20 (cfs)", log = "y", letter = "B")
s5_reach_stack(s5_drift("days_since_wy_start_ma20"),
               "Peak timing MA20 (days)")

dev.off()
cat("Wrote", file.path(OUT_DIR, "6_hydrology_drift_by_reach.pdf"), "\n")
# ============================================================================
# SECTION 6 -- SI  Figure: Recursive-mode validation  (from diag_recursive_bias_check)
# ----------------------------------------------------------------------------
# Split into two single-page files:
#   6a SI_recursive_signed_error.pdf   -- THE actual bias test (validation
#      years stay in observed driver range). Flat/near-zero trend + steady
#      ~0.08 offset = no growing bias, distinct from the ribbon-width gap.
#   6b SI_recursive_maxmean_overlay.pdf-- visual reference ONLY. Runs 75yr
#      under drifting drivers, so it CONFLATES lag_y-scale effects with driver
#      trends -- NOT a test of recursion scale. Do not read a scale conclusion
#      off it. Kept for the "is there a visible seed-year discontinuity" check.
# ============================================================================
s6_rec  <- read.csv("2_incremental/recursive_val_predictions.csv", stringsAsFactors = FALSE)
.require_cols(s6_rec, c("Site", "Year", "step", "Observed", "Predicted_rec", "Error"),
              "recursive_val_predictions.csv")

s6_maxstep <- max(s6_rec$step, na.rm = TRUE)
s6_ebs <- do.call(rbind, lapply(seq_len(s6_maxstep), function(k) {
  d <- s6_rec[s6_rec$step == k & is.finite(s6_rec$Error), ]
  n <- nrow(d)
  data.frame(step = k, n = n,
             MeanErr = round(mean(d$Error), 4),
             SD_Err  = round(if (n > 1) sd(d$Error) else NA_real_, 4),
             SE_Err  = round(if (n > 1) sd(d$Error) / sqrt(n) else NA_real_, 4),
             RMSE    = round(sqrt(mean(d$Error^2)), 4))
}))
s6_trend <- if (nrow(s6_ebs) >= 3 && var(s6_ebs$step) > 0) {
  tryCatch(lm(MeanErr ~ step, data = s6_ebs), error = function(e) NULL)
} else NULL
s6_slope <- if (!is.null(s6_trend)) round(coef(s6_trend)[["step"]], 5) else NA_real_

# --- 6a. signed error by step (the test) ------------------------------------
pdf(file.path(OUT_DIR, "SI_recursive_signed_error.pdf"),
    width = 8, height = 6, family = "Helvetica")
par(mfrow = c(1, 1), mar = c(4.2, 4.2, 2, 1))
s6_has_se <- any(is.finite(s6_ebs$SE_Err))
s6_ye <- range(c(s6_ebs$MeanErr,
                 if (s6_has_se) s6_ebs$MeanErr - s6_ebs$SE_Err else NULL,
                 if (s6_has_se) s6_ebs$MeanErr + s6_ebs$SE_Err else NULL, 0), na.rm = TRUE)
s6_ye <- s6_ye + c(-0.05, 0.05) * diff(s6_ye)
plot(s6_ebs$step, s6_ebs$MeanErr, type = "b", pch = 16, col = "darkorange", lwd = 2,
     ylim = s6_ye, xlab = "Steps ahead (years)",
     ylab = "Mean signed error, Observed minus Predicted (recursive)")
if (s6_has_se)
  segments(s6_ebs$step, s6_ebs$MeanErr - s6_ebs$SE_Err,
           s6_ebs$step, s6_ebs$MeanErr + s6_ebs$SE_Err, col = "darkorange")
abline(h = 0, lty = 2, col = "grey40")
if (!is.null(s6_trend)) {
  abline(s6_trend, lty = 3, col = "steelblue")
  legend("topleft", legend = paste0("linear trend: ", s6_slope, " / step"),
         bty = "n", text.col = "steelblue", cex = 0.85)
}
text(s6_ebs$step, s6_ebs$MeanErr, paste0("n=", s6_ebs$n), pos = 3, cex = 0.65, col = "grey40")
mtext("positive = model underpredicts (annual-max-vs-conditional-mean direction)",
      side = 3, line = 0.3, cex = 0.8, col = "grey30")
dev.off()
cat("Wrote", file.path(OUT_DIR, "SI_recursive_signed_error.pdf"), "\n")

# --- 6b. annual max/mean overlay (visual reference only) --------------------
s6_obs  <- read.csv("2_incremental/ucfr_model_ready.csv", stringsAsFactors = FALSE)
s6_proj <- read.csv("2_incremental/bloom_projections.csv", stringsAsFactors = FALSE)
.require_cols(s6_obs,  c("Site", "Year", "logCHLa"), "ucfr_model_ready.csv")
.require_cols(s6_proj, c("Site", "bracket", "year", "median_logCHLa",
                         "p10_logCHLa", "p90_logCHLa"), "bloom_projections.csv")

s6_oc <- s6_obs[is.finite(s6_obs$logCHLa), ]
s6_amax  <- aggregate(logCHLa ~ Site + Year, data = s6_oc, FUN = max);  names(s6_amax)[3]  <- "ann_max"
s6_amean <- aggregate(logCHLa ~ Site + Year, data = s6_oc, FUN = mean); names(s6_amean)[3] <- "ann_mean"
s6_ann <- merge(s6_amax, s6_amean, by = c("Site", "Year"))

s6_xlim <- range(c(s6_proj$year, s6_ann$Year), na.rm = TRUE)
s6_ylim <- range(c(s6_proj$p10_logCHLa, s6_proj$p90_logCHLa,
                   s6_ann$ann_max, s6_ann$ann_mean), na.rm = TRUE)
s6_ylim <- s6_ylim + c(-0.03, 0.03) * diff(s6_ylim)
s6_fill_low  <- adjustcolor(COL_LOW,  alpha.f = 0.22)
s6_fill_high <- adjustcolor(COL_HIGH, alpha.f = 0.22)
s6_ribbon <- function(d, fill) {
  d <- d[order(d$year), ]
  polygon(c(d$year, rev(d$year)), c(d$p10_logCHLa, rev(d$p90_logCHLa)), col = fill, border = NA)
}
s6_med <- function(d, col) { d <- d[order(d$year), ]; lines(d$year, d$median_logCHLa, col = col, lwd = 2) }

pdf(file.path(OUT_DIR, "SI_recursive_maxmean_overlay.pdf"), width = 9, height = 11)
par(mfrow = c(4, 2), mar = c(3.2, 3.6, 2.2, 0.8),
    mgp = c(2.1, 0.6, 0), oma = c(2.2, 1.5, 2.4, 0.5))
for (st in SITE_ORDER) {
  ds <- s6_proj[s6_proj$Site == st, ]
  lo <- ds[ds$bracket == "low",  ]
  hi <- ds[ds$bracket == "high", ]
  da <- s6_ann[s6_ann$Site == st, ]; da <- da[order(da$Year), ]
  plot(NA, xlim = s6_xlim, ylim = s6_ylim,
       xlab = "Year", ylab = expression(log[10] * " Chl " * italic(a)), main = st)
  if (nrow(lo)) { s6_ribbon(lo, s6_fill_low);  s6_med(lo, COL_LOW) }
  if (nrow(hi)) { s6_ribbon(hi, s6_fill_high); s6_med(hi, COL_HIGH) }
  if (nrow(da)) {
    lines(da$Year, da$ann_max,  col = "black",  lty = 1, lwd = 1)
    points(da$Year, da$ann_max, col = "black",  pch = 17, cex = 0.8)
    lines(da$Year, da$ann_mean,  col = "grey40", lty = 2, lwd = 1)
    points(da$Year, da$ann_mean, col = "grey40", pch = 16, cex = 0.8)
  }
  box()
}
plot.new()
legend("center",
       legend = c("Low median (RCP4.5/SSP245)", "High median (RCP8.5/SSP585)",
                  "10th-90th percentile", "Observed annual max", "Observed annual mean"),
       col = c(COL_LOW, COL_HIGH, "grey50", "black", "grey40"),
       lwd = c(2, 2, NA, 1, 1), pch = c(NA, NA, 15, 17, 16), lty = c(1, 1, NA, 1, 2),
       pt.cex = c(2.2, 2.2, 2.2, 0.9, 0.9), bty = "n", cex = 0.9)
mtext("SI: annual max/mean vs projection ribbons (visual reference only)",
      outer = TRUE, cex = 1.0, font = 2, line = 0.6)
mtext("NOT a test of recursion scale -- conflates lag_y with 75yr driver drift; see SI_recursive_signed_error",
      outer = TRUE, side = 1, cex = 0.65, line = 0.6)
dev.off()
cat("Wrote", file.path(OUT_DIR, "SI_recursive_maxmean_overlay.pdf"), "\n")


# ============================================================================
# SECTION 7 -- SI  Figure: Temperature projections by site  (from 15 Section D)
# ----------------------------------------------------------------------------
# Single panel, 7 sites (not grouped by reach -- MS/BM/HU share CLABE but
# diverge on ceiling-crossing). Both brackets per site, calibration ceiling
# marked, 2026-2098.
# ============================================================================
s7_te <- read.csv("2_incremental/ncar_temperature_envelope.csv", stringsAsFactors = FALSE)
s7_te <- unique(s7_te[, c("site", "water_year", "Temp_oC_low", "Temp_oC_high")])
s7_te <- s7_te[s7_te$water_year >= 2026 & s7_te$water_year <= 2098, ]

s7_colors <- setNames(
  c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
  SITE_ORDER)
S7_CEILING <- 22.2   # observed calibration ceiling, deg C (see 3.4 text)

pdf(file.path(OUT_DIR, "temp_projections_by_site.pdf"),
    width = 7, height = 5, family = "Helvetica")
par(mar = c(4, 4.5, 2, 1))
s7_yr <- range(c(s7_te$Temp_oC_low, s7_te$Temp_oC_high), na.rm = TRUE)
plot(NA, xlim = c(2026, 2098), ylim = s7_yr,
     xlab = "Water year", ylab = "Stream temperature, deg C")
abline(h = S7_CEILING, lty = 2, col = "gray50")
for (s in SITE_ORDER) {
  d <- s7_te[s7_te$site == s, ]; d <- d[order(d$water_year), ]
  lines(d$water_year, d$Temp_oC_high, col = s7_colors[s], lwd = 2, lty = 1)
  lines(d$water_year, d$Temp_oC_low,  col = s7_colors[s], lwd = 2, lty = 3)
}
legend("topleft", legend = SITE_ORDER, col = s7_colors[SITE_ORDER],
       lwd = 2, bty = "n", cex = 0.8, ncol = 2)
legend("bottomright", legend = c("High bracket (RCP8.5/SSP585)", "Low bracket (RCP4.5/SSP245)"),
       lty = c(1, 3), lwd = 2, col = "gray30", bty = "n", cex = 0.75)
dev.off()
cat("Wrote", file.path(OUT_DIR, "temp_projections_by_site.pdf"), "\n")


# ============================================================================
# SECTION 8 -- SI  Figure: Ensemble variance by year  (from 15 Section E)
# ----------------------------------------------------------------------------
# SD of pred_logCHLa across members within each Site x bracket x year cell,
# averaged across the 7 sites per bracket x year. Growing SD (not collapsing)
# supports the recursion-relaxation caveat.
# ============================================================================
s8_mem <- read.csv("2_incremental/bloom_projections_members.csv", stringsAsFactors = FALSE)
s8_site_sd <- aggregate(pred_logCHLa ~ Site + bracket + year, data = s8_mem,
                        FUN = function(x) sd(x, na.rm = TRUE))
names(s8_site_sd)[4] <- "sd_logCHLa"
s8_bsd <- aggregate(sd_logCHLa ~ bracket + year, data = s8_site_sd, FUN = mean, na.rm = TRUE)
s8_bcol <- c(low = "#56B4E9", high = "#D55E00")

pdf(file.path(OUT_DIR, "ensemble_variance_by_year.pdf"),
    width = 7, height = 5, family = "Helvetica")
par(mar = c(4, 4.5, 2, 1))
s8_xr <- range(s8_bsd$year, na.rm = TRUE)
s8_yr <- range(s8_bsd$sd_logCHLa, na.rm = TRUE)
plot(NA, xlim = s8_xr, ylim = s8_yr,
     xlab = "Water year", ylab = "SD of predicted log10(CHLa) across members")
for (b in c("low", "high")) {
  d <- s8_bsd[s8_bsd$bracket == b, ]; d <- d[order(d$year), ]
  lines(d$year, d$sd_logCHLa, col = s8_bcol[b], lwd = 2)
}
legend("topleft", legend = c("Low bracket", "High bracket"),
       col = s8_bcol[c("low", "high")], lwd = 2, bty = "n", cex = 0.85)
dev.off()
cat("Wrote", file.path(OUT_DIR, "ensemble_variance_by_year.pdf"), "\n")


# ============================================================================
# DONE
# ============================================================================
cat("\nAll figures written to ", OUT_DIR, "/\n", sep = "")

