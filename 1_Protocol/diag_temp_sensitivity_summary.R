# ============================================================================
# diag_temp_sensitivity_summary.R
# UCFR Cladophora Bloom Prediction Pipeline -- Temperature-smooth sensitivity
# (defensive/reviewer-response analysis, NOT manuscript-facing)
#
# Cross-variant comparison, reading the outputs of
# diag_temp_sensitivity_projection.R (already run for all 5 variants).
# Does not touch M1, tp_submodel, or any predict() call -- pure read/
# summarize/plot over the CSVs already on disk.
#
# Inputs (all in 4_products/diagnostics/temp_sensitivity/):
#   bloom_projections_members_<variant>.csv   x5 (V0, V1, V2, V3mild, V3strong)
#   temp_sensitivity_eoc.csv                  end-of-century summary, all variants
#   temp_sensitivity_scorecard.csv            per-variant run summary
#
# Outputs (same folder):
#   temp_sensitivity_comparison.pdf   per-site end-of-century P(exceed) across
#                                      variants (both brackets) + driver-share
#                                      panel by variant
#   temp_sensitivity_summary.csv      numeric table for a reviewer-response
#                                      letter: end-century median logCHLa,
#                                      P(exceed 150 mg/m^2), Temp-driver share,
#                                      by variant x site x bracket
# ============================================================================

library(mgcv)

VARIANTS <- c("V0", "V1", "V2", "V3mild", "V3strong")
DIAG_DIR <- "4_products/diagnostics/temp_sensitivity/"

EXCEEDANCE_THRESHOLD_MG_M2 <- 150     # Suplee et al. 2009 JAWRA, log10(mg/m^2) scale
LOG_EXCEEDANCE_THRESHOLD   <- log10(EXCEEDANCE_THRESHOLD_MG_M2)

PATH_M1 <- "3_models/bloom_model_M1.rds"
PATH_RESID_POOL <- "2_incremental/recursive_val_predictions.csv"  # Scheme D residuals (Error column)
SEED <- 20260701
B_DRAWS_PER_MEMBER <- 2000

OUT_PDF     <- paste0(DIAG_DIR, "temp_sensitivity_comparison.pdf")
OUT_SUMMARY <- paste0(DIAG_DIR, "temp_sensitivity_summary.csv")


# ----------------------------------------------------------------------------
# 1. LOAD -- member-level projections for all 5 variants
# ----------------------------------------------------------------------------
members_list <- setNames(vector("list", length(VARIANTS)), VARIANTS)
for (v in VARIANTS) {
  f <- paste0(DIAG_DIR, "bloom_projections_members_", v, ".csv")
  if (!file.exists(f)) stop("Missing members file for variant ", v, ": ", f)
  df <- read.csv(f, stringsAsFactors = FALSE)
  df$temp_variant <- v
  members_list[[v]] <- df
}
members_all <- do.call(rbind, members_list)
rownames(members_all) <- NULL

require(mgcv)
m1 <- readRDS(PATH_M1)


# ----------------------------------------------------------------------------
# 2. RESIDUAL BOOTSTRAP -- reuse Scheme D's already-calibrated recursive
#    residuals (unchanged pool: only member predictions differ across
#    variants, not the error distribution being resampled onto them).
# ----------------------------------------------------------------------------
if (!file.exists(PATH_RESID_POOL))
  stop("Residual pool not found at ", PATH_RESID_POOL,
       " -- this script reuses Scheme D's recursive-validation errors ",
       "and does not regenerate them.")

resid_pool <- read.csv(PATH_RESID_POOL, stringsAsFactors = FALSE)
if (!"Error" %in% names(resid_pool))
  stop("Expected an 'Error' column (Observed - Predicted_rec) in ", PATH_RESID_POOL,
       " but found: ", paste(names(resid_pool), collapse = ", "))
residuals_vec <- resid_pool$Error
residuals_vec <- residuals_vec[is.finite(residuals_vec)]

set.seed(SEED)

exceed_by_group <- function(pred_vec, n_draws_per_member) {
  draws <- sample(residuals_vec, size = length(pred_vec) * n_draws_per_member, replace = TRUE)
  draws <- matrix(draws, nrow = length(pred_vec), ncol = n_draws_per_member)
  boot_pred <- pred_vec + draws
  mean(boot_pred > LOG_EXCEEDANCE_THRESHOLD)
}

eoc_year <- max(members_all$year)

exceed_tab_rows <- list()
i <- 1
for (v in VARIANTS) {
  d_v <- members_all[members_all$temp_variant == v & members_all$year == eoc_year, ]
  for (st in sort(unique(d_v$Site))) {
    for (br in sort(unique(d_v$bracket))) {
      dd <- d_v[d_v$Site == st & d_v$bracket == br, ]
      if (!nrow(dd)) next
      p_exceed <- exceed_by_group(dd$pred_logCHLa, B_DRAWS_PER_MEMBER)
      exceed_tab_rows[[i]] <- data.frame(
        temp_variant = v, Site = st, bracket = br,
        eoc_year = eoc_year,
        median_logCHLa = median(dd$pred_logCHLa),
        n_members = nrow(dd),
        p_exceed_150 = p_exceed,
        stringsAsFactors = FALSE
      )
      i <- i + 1
    }
  }
}
exceed_tab <- do.call(rbind, exceed_tab_rows)


# ----------------------------------------------------------------------------
# 3. DRIVER-SHARE BY VARIANT -- CHANGED to match diag_freeze_one_driver.R's
#    actual method exactly (confirmed by reading that script): freeze
#    Temp_oC at each (Site, member)'s 2026 value across the WHOLE recursive
#    trajectory (not a single-year cross-member snapshot -- Temp_oC is
#    member-invariant within a Site x bracket x year, so freezing across
#    members at one year is always a no-op, which is what produced the
#    all-zero bug), re-run the SAME recursion engine used for the real
#    projection, and compare early-window (<=2035) vs late-window (>=2089)
#    mean logCHLa between the real and frozen trajectories. Applied per
#    variant so the Temp share is re-estimated under each Temp-treatment
#    assumption, on member-level data already on disk from
#    diag_temp_sensitivity_projection.R -- no re-fitting, no new predict()
#    logic beyond what that script already validated.
# ----------------------------------------------------------------------------
is_early <- function(y) y <= 2035
is_late  <- function(y) y >= 2089

run_recursion_temp_frozen <- function(grid) {
  years <- sort(unique(grid$year))
  proj_start <- min(years)
  grid$pred_logCHLa <- NA_real_
  
  seed_row <- grid[grid$year == proj_start, c("uid", "lag_y")]
  seed_row <- seed_row[!duplicated(seed_row$uid), ]
  seed_lag <- setNames(seed_row$lag_y, seed_row$uid)
  lag_state <- seed_lag[unique(grid$uid)]
  
  anchor <- grid[grid$year == proj_start, c("uid", "Temp_oC")]
  anchor <- anchor[!duplicated(anchor$uid), ]
  fixed_temp <- setNames(anchor$Temp_oC, anchor$uid)
  grid$Temp_oC <- fixed_temp[grid$uid]
  
  for (t in years) {
    idx <- which(grid$year == t)
    nd  <- grid[idx, ]
    nd$lag_y <- lag_state[nd$uid]
    p <- as.numeric(predict(m1, newdata = nd[, c("lag_y", "anomaly",
                                                 "logQ_obs_cfs", "Days_Since_Freshet",
                                                 "logTP_mg_L", "Temp_oC", "Site")]))
    grid$pred_logCHLa[idx] <- p
    lag_state[nd$uid] <- p
  }
  grid$pred_logCHLa
}

window_mean <- function(grid, pred, keep_fun) {
  d <- data.frame(Site = grid$Site, bracket = grid$bracket, year = grid$year, pred = pred)
  d <- d[keep_fun(d$year), ]
  ag <- aggregate(pred ~ Site + bracket, data = d, FUN = mean, na.rm = TRUE)
  names(ag)[3] <- "mean_logCHLa"
  ag
}

temp_share_rows <- list()
for (v in VARIANTS) {
  
  d_v <- members_all[members_all$temp_variant == v, ]
  d_v$Site <- factor(d_v$Site, levels = levels(m1$model$Site))
  d_v$uid  <- paste(d_v$Site, d_v$member, sep = "|")
  
  grid_v <- d_v[, c("Site", "uid", "bracket", "year", "anomaly",
                    "logQ_obs_cfs", "Days_Since_Freshet", "logTP_mg_L",
                    "Temp_oC", "lag_y")]
  
  frozen_pred <- run_recursion_temp_frozen(grid_v)
  
  actual_early <- window_mean(d_v, d_v$pred_logCHLa, is_early)
  actual_late  <- window_mean(d_v, d_v$pred_logCHLa, is_late)
  actual <- merge(actual_early, actual_late, by = c("Site", "bracket"),
                  suffixes = c("_early", "_late"))
  actual$actual_delta <- actual$mean_logCHLa_late - actual$mean_logCHLa_early
  
  frozen_early <- window_mean(grid_v, frozen_pred, is_early)
  frozen_late  <- window_mean(grid_v, frozen_pred, is_late)
  frozen <- merge(frozen_early, frozen_late, by = c("Site", "bracket"),
                  suffixes = c("_early", "_late"))
  frozen$frozen_delta <- frozen$mean_logCHLa_late - frozen$mean_logCHLa_early
  
  m <- merge(actual[, c("Site", "bracket", "actual_delta")],
             frozen[, c("Site", "bracket", "frozen_delta")],
             by = c("Site", "bracket"))
  m$temp_contribution <- m$actual_delta - m$frozen_delta
  m$temp_driver_share <- ifelse(m$actual_delta != 0,
                                m$temp_contribution / m$actual_delta, NA_real_)
  m$temp_variant <- v
  
  temp_share_rows[[v]] <- m[, c("temp_variant", "Site", "bracket",
                                "actual_delta", "frozen_delta",
                                "temp_contribution", "temp_driver_share")]
}
temp_share_tab <- do.call(rbind, temp_share_rows)


# ----------------------------------------------------------------------------
# 4. UNIFIED SUMMARY TABLE
# ----------------------------------------------------------------------------
summary_tab <- merge(exceed_tab, temp_share_tab,
                     by = c("temp_variant", "Site", "bracket"))
summary_tab <- summary_tab[order(summary_tab$Site, summary_tab$bracket,
                                 factor(summary_tab$temp_variant, levels = VARIANTS)), ]

write.csv(summary_tab, OUT_SUMMARY, row.names = FALSE)


# ----------------------------------------------------------------------------
# 5. FIGURE -- base R, two panels per page: (a) P(exceed) by site across
#    variants, split by bracket; (b) Temp driver share by site across variants.
#    No main= titles (axis labels already identify the panel).
# ----------------------------------------------------------------------------
variant_pch <- c(V0 = 16, V1 = 17, V2 = 15, V3mild = 3, V3strong = 4)
variant_col <- c(V0 = "black", V1 = "#1b9e77", V2 = "#7570b3",
                 V3mild = "#d95f02", V3strong = "#e7298a")

sites <- sort(unique(summary_tab$Site))

pdf(OUT_PDF, width = 10, height = 7)

for (br in c("low", "high")) {
  d_br <- summary_tab[summary_tab$bracket == br, ]
  plot(NA, xlim = c(1, length(sites)), ylim = c(0, 1),
       xaxt = "n", xlab = "Site", ylab = sprintf("P(exceed 150 mg/m^2), %s bracket, EOC %d", br, eoc_year))
  axis(1, at = seq_along(sites), labels = sites)
  abline(h = c(0.25, 0.5, 0.75), col = "grey85", lty = 3)
  for (v in VARIANTS) {
    dv <- d_br[d_br$temp_variant == v, ]
    dv <- dv[match(sites, dv$Site), ]
    points(seq_along(sites), dv$p_exceed_150, pch = variant_pch[v], col = variant_col[v])
    lines(seq_along(sites), dv$p_exceed_150, col = variant_col[v], lty = 2)
  }
  legend("topleft", legend = VARIANTS, pch = variant_pch[VARIANTS],
         col = variant_col[VARIANTS], bty = "n", cex = 0.8)
}

for (br in c("low", "high")) {
  d_br <- summary_tab[summary_tab$bracket == br, ]
  plot(NA, xlim = c(1, length(sites)), ylim = range(d_br$temp_driver_share, na.rm = TRUE),
       xaxt = "n", xlab = "Site", ylab = sprintf("Temp driver share, %s bracket, EOC %d", br, eoc_year))
  axis(1, at = seq_along(sites), labels = sites)
  for (v in VARIANTS) {
    dv <- d_br[d_br$temp_variant == v, ]
    dv <- dv[match(sites, dv$Site), ]
    points(seq_along(sites), dv$temp_driver_share, pch = variant_pch[v], col = variant_col[v])
    lines(seq_along(sites), dv$temp_driver_share, col = variant_col[v], lty = 2)
  }
  legend("topleft", legend = VARIANTS, pch = variant_pch[VARIANTS],
         col = variant_col[VARIANTS], bty = "n", cex = 0.8)
}

dev.off()


# ----------------------------------------------------------------------------
# 6. SCORECARD
# ----------------------------------------------------------------------------
cat("\n============ diag_temp_sensitivity_summary.R SCORECARD ============\n")
cat("\n-- eoc_year used:", eoc_year, "--\n")
cat("\n-- rows per variant in members_all --\n")
print(table(members_all$temp_variant))
cat("\n-- summary table (head) --\n")
print(head(summary_tab, 10), row.names = FALSE)
cat("\n-- outputs written --\n")
print(data.frame(file = c(OUT_PDF, OUT_SUMMARY)), row.names = FALSE)
cat("\n=====================================================================\n")