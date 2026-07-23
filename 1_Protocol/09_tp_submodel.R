# ============================================================================
# 09_tp_submodel.R
# UCFR Filamentous Algae Project
# TP submodel — predicts log TP from climate-projectable hydroclimate drivers
#
# Role: the single chemistry bridge in the GAM projection chain. The bloom
#   model M1 consumes a log-TP predictor; this submodel supplies it under
#   climate scenarios. (TN excluded, SPC excluded from M1, so TP is the only
#   chemistry prediction the pipeline needs.)
#
# RESPONSE SCALE (see TP_SCALE below):
#   The model-ready column `logTP_mg_L` is log10(1 + TP_mg_L) — a log1p applied
#   as a blanket across nutrient columns (sensible for SRP/NH4/NO3, which have
#   zeros). For TP this is ~linear over the observed range (~0.004-0.06 mg/L),
#   so it does little to stabilize variance. We therefore derive the response
#   here directly from raw TP_mg_L rather than trusting that column:
#     "log10_ugL"  : log10(TP_mg_L * 1000)  -- plain log10 of concentration.
#                    Positive-valued; Suplee 24 ug/L threshold = log10(24)=1.38.
#                    (Identical to log10(mg/L) up to a +3 intercept shift.)
#     "log1p_mgL"  : log10(1 + TP_mg_L)     -- reproduces the original column,
#                    for side-by-side comparison.
#
# Structure (decided):
#   <logTP> ~ s(anomaly) + s(logQ_obs_cfs) + s(Days_Since_Freshet)
#             + s(Temp_oC) + s(Site, bs="re")
#   - GAM (mgcv), not BRT — extrapolation consistency with the pivot.
#   - Site random intercept carries site-level TP structure (geological P at
#     GR/BN, WWTP point source at BM) hydroclimate cannot see. Projectable
#     because site identity is fixed; assumes site-baseline stationarity
#     (watershed-wide caveat; WWTP-upgrade knob deferred).
#   - Per-variable k via k_spec.
#
# Validation:
#   - In-sample R2 / deviance explained.
#   - LOSO (leave-one-site-out, RE EXCLUDED -> population-level prediction):
#     strict spatial test; isolates the hydroclimate-only signal. Expected
#     weak, esp. BM/FH/GR. Reported as pooled and within-site R2.
#   - LOYO (leave-one-year-out, RE RETAINED): the operationally relevant test
#     — in projection the site is always one of the known 7, only the year is
#     out of sample.
#   - Per-site residuals to quantify the BM gap explicitly.
#
# Propagation hook: fitted object saved as .rds so Vp and residual SD are
#   retrievable later for Monte Carlo posterior propagation. Point-estimate
#   build only.
#
# Input:    2_incremental/ucfr_model_ready.csv
# Outputs:  3_models/tp_submodel.rds
#           2_incremental/tp_submodel_predictions.csv
#           4_products/diagnostics/tp_submodel.pdf
# ============================================================================

library(mgcv)

# --- CONFIGURATION ---------------------------------------------------------
TP_SCALE       <- "log10_ugL"    # "log10_ugL" (recommended) or "log1p_mgL"
PREDICTORS     <- c("anomaly", "logQ_obs_cfs", "Days_Since_Freshet", "Temp_oC")
k_spec         <- c(anomaly = 5, logQ_obs_cfs = 5,
                    Days_Since_Freshet = 5, Temp_oC = 5)
OUTLIER_SD     <- 2.0
OUTLIER_ROUNDS <- 2
# ---------------------------------------------------------------------------

# --- HELPER: fit GAM from predictor list + k_spec (self-contained) ---------
fit_gam <- function(predictors, k_spec, response, data, site_re = TRUE) {
  terms <- vapply(predictors,
                  function(p) paste0("s(", p, ", k=", k_spec[[p]], ")"),
                  character(1))
  if (site_re) terms <- c(terms, "s(Site, bs='re')")
  f <- as.formula(paste(response, "~", paste(terms, collapse = " + ")))
  gam(f, data = data, method = "REML")
}

# --- METRIC HELPERS --------------------------------------------------------
r2_pooled <- function(obs, pred)
  1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
rmse_fn   <- function(obs, pred) sqrt(mean((obs - pred)^2))
within_r2 <- function(df) {           # vs site means (harder, removes offset)
  ss_res <- sum((df$Observed - df$Predicted)^2)
  ss_tot <- sum((df$Observed - ave(df$Observed, df$Site))^2)
  1 - ss_res / ss_tot
}

# --- 1. READ & DERIVE RESPONSE ---------------------------------------------
dat <- read.csv("2_incremental/ucfr_model_ready.csv",
                stringsAsFactors = FALSE)

# derive TP response directly from raw TP_mg_L (do not trust logTP_mg_L column)
if (TP_SCALE == "log10_ugL") {
  dat$TP_response <- log10(dat$TP_mg_L * 1000)   # = log10(mg/L) + 3
} else if (TP_SCALE == "log1p_mgL") {
  dat$TP_response <- log10(1 + dat$TP_mg_L)      # reproduces logTP_mg_L
} else {
  stop("TP_SCALE must be 'log10_ugL' or 'log1p_mgL'")
}
TP_RESPONSE <- "TP_response"

keep_cols <- c(TP_RESPONSE, "TP_mg_L", PREDICTORS, "Site", "Year", "Month")
mdat <- dat[complete.cases(dat[, c(TP_RESPONSE, PREDICTORS)]), keep_cols]
mdat$Site <- factor(mdat$Site)

# --- 2. OUTLIER REMOVAL (+/- SD residuals, same logic as core scripts) -----
outlier_log <- data.frame()
for (r in seq_len(OUTLIER_ROUNDS)) {
  m_tmp     <- fit_gam(PREDICTORS, k_spec, TP_RESPONSE, mdat)
  resid_tmp <- residuals(m_tmp)
  thr       <- OUTLIER_SD * sd(resid_tmp)
  out       <- which(abs(resid_tmp) > thr)
  if (length(out) == 0) break
  log_r       <- mdat[out, c("Site", "Year", "Month", "TP_mg_L", TP_RESPONSE)]
  log_r$resid <- round(resid_tmp[out], 3)
  log_r$round <- r
  outlier_log <- rbind(outlier_log, log_r)
  mdat <- mdat[-out, ]
  mdat$Site <- factor(mdat$Site)      # drop any emptied level
}

# --- 3. FULL MODEL ---------------------------------------------------------
m_full <- fit_gam(PREDICTORS, k_spec, TP_RESPONSE, mdat)
mdat$resid_full <- residuals(m_full)
saveRDS(m_full, "3_models/tp_submodel.rds")   # propagation hook

# concurvity: overall (each term vs all others) + pairwise (anomaly vs logQ)
cc_overall  <- round(concurvity(m_full, full = TRUE), 3)
cc_pairwise <- round(concurvity(m_full, full = FALSE)$estimate, 3)

# --- 4. LOSO (leave-one-site-out, RE EXCLUDED) -----------------------------
sites     <- levels(mdat$Site)
loso_list <- list()
for (s in sites) {
  train <- mdat[mdat$Site != s, ]; train$Site <- factor(train$Site)
  test  <- mdat[mdat$Site == s, ]
  m_s <- tryCatch(fit_gam(PREDICTORS, k_spec, TP_RESPONSE, train),
                  error = function(e) NULL)
  if (is.null(m_s)) next
  test2 <- test                       # dummy training level so design builds;
  test2$Site <- factor(levels(train$Site)[1],   # exclude= then zeroes it out
                       levels = levels(train$Site))
  pred <- predict(m_s, newdata = test2, exclude = "s(Site)")
  loso_list[[s]] <- data.frame(scheme = "LOSO", Site = test$Site,
                               Year = test$Year, Month = test$Month,
                               Observed = test[[TP_RESPONSE]],
                               Predicted = as.numeric(pred),
                               stringsAsFactors = FALSE)
}
loso_df <- do.call(rbind, loso_list); row.names(loso_df) <- NULL

# --- 5. LOYO (leave-one-year-out, RE RETAINED) -----------------------------
years     <- sort(unique(mdat$Year))
loyo_list <- list()
for (y in years) {
  train <- mdat[mdat$Year != y, ]; train$Site <- factor(train$Site)
  test  <- mdat[mdat$Year == y, ]
  if (nrow(test) == 0) next
  m_y <- tryCatch(fit_gam(PREDICTORS, k_spec, TP_RESPONSE, train),
                  error = function(e) NULL)
  if (is.null(m_y)) next
  pred <- predict(m_y, newdata = test)          # site known -> RE used
  loyo_list[[as.character(y)]] <-
    data.frame(scheme = "LOYO", Site = test$Site,
               Year = test$Year, Month = test$Month,
               Observed = test[[TP_RESPONSE]],
               Predicted = as.numeric(pred),
               stringsAsFactors = FALSE)
}
loyo_df <- do.call(rbind, loyo_list); row.names(loyo_df) <- NULL

# --- 6. SAVE PREDICTIONS ---------------------------------------------------
all_preds <- rbind(loso_df, loyo_df)
all_preds$TP_scale <- TP_SCALE
write.csv(all_preds, "2_incremental/tp_submodel_predictions.csv",
          row.names = FALSE)

# --- 7. DIAGNOSTIC PDF -----------------------------------------------------
pdf("4_products/diagnostics/tp_submodel.pdf", width = 10, height = 8)

rng <- range(c(all_preds$Observed, all_preds$Predicted))   # shared axes
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

plot(loso_df$Predicted, loso_df$Observed, xlim = rng, ylim = rng,
     pch = 16, cex = 0.7, col = adjustcolor("steelblue", 0.6), bty = "o",
     xlab = "Predicted logTP", ylab = "Observed logTP",
     main = paste0("LOSO (RE excl)  R2 = ",
                   round(r2_pooled(loso_df$Observed, loso_df$Predicted), 3),
                   "  [", TP_SCALE, "]"))
abline(0, 1, lty = 2)

plot(loyo_df$Predicted, loyo_df$Observed, xlim = rng, ylim = rng,
     pch = 16, cex = 0.7, col = adjustcolor("darkgreen", 0.6), bty = "o",
     xlab = "Predicted logTP", ylab = "Observed logTP",
     main = paste0("LOYO (RE ret)  R2 = ",
                   round(r2_pooled(loyo_df$Observed, loyo_df$Predicted), 3),
                   "  [", TP_SCALE, "]"))
abline(0, 1, lty = 2)

boxplot(resid_full ~ Site, data = mdat, bty = "o",
        ylab = "Full-model residual", xlab = "Site",
        main = "Residuals by site", col = "grey85")
abline(h = 0, lty = 2)

loso_site_r2 <- sapply(split(loso_df, loso_df$Site),
                       function(d) r2_pooled(d$Observed, d$Predicted))
barplot(loso_site_r2, las = 2, ylab = "LOSO R2 (RE excl)",
        main = "LOSO skill by held-out site", col = "grey85")
abline(h = 0, lty = 2)

par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
plot(m_full, pages = 0, shade = TRUE,
     shade.col = adjustcolor("steelblue", 0.3), seWithMean = TRUE)

dev.off()

# --- 8. SCORECARD + SUPPORTING TABLES (printed at end) ---------------------
scorecard <- data.frame(
  Scheme    = c("Full model (in-sample)", "LOSO (RE excluded)",
                "LOYO (RE retained)"),
  n         = c(nrow(mdat), nrow(loso_df), nrow(loyo_df)),
  n_sites   = c(nlevels(mdat$Site), length(unique(loso_df$Site)),
                length(unique(loyo_df$Site))),
  n_years   = c(length(unique(mdat$Year)), length(unique(loso_df$Year)),
                length(unique(loyo_df$Year))),
  R2_pooled = c(round(summary(m_full)$r.sq, 4),
                round(r2_pooled(loso_df$Observed, loso_df$Predicted), 4),
                round(r2_pooled(loyo_df$Observed, loyo_df$Predicted), 4)),
  R2_within = c(NA, round(within_r2(loso_df), 4), NA),
  RMSE      = c(round(rmse_fn(mdat[[TP_RESPONSE]], fitted(m_full)), 4),
                round(rmse_fn(loso_df$Observed, loso_df$Predicted), 4),
                round(rmse_fn(loyo_df$Observed, loyo_df$Predicted), 4)),
  stringsAsFactors = FALSE
)

sites_ord   <- names(loso_site_r2)
site_detail <- data.frame(
  Site       = sites_ord,
  n          = as.integer(table(loso_df$Site)[sites_ord]),
  LOSO_R2    = round(loso_site_r2[sites_ord], 3),
  resid_mean = round(tapply(mdat$resid_full, mdat$Site, mean)[sites_ord], 3),
  resid_sd   = round(tapply(mdat$resid_full, mdat$Site, sd)[sites_ord], 3),
  stringsAsFactors = FALSE
)
row.names(site_detail) <- NULL

dev_expl <- round(summary(m_full)$dev.expl * 100, 1)

cat("=== TP SUBMODEL SCORECARD  [scale:", TP_SCALE, "] ===\n")
cat("Deviance explained (full model):", dev_expl, "%\n")
cat("N outliers removed:", nrow(outlier_log), "of",
    nrow(mdat) + nrow(outlier_log), "\n\n")
print(scorecard, row.names = FALSE)
cat("\n=== PER-SITE DETAIL (LOSO R2 + full-model residuals) ===\n")
print(site_detail, row.names = FALSE)
cat("\n=== CONCURVITY overall (rows: worst/observed/estimate) ===\n")
print(cc_overall)
cat("\n=== CONCURVITY pairwise estimate (watch anomaly vs logQ; >0.8 concern) ===\n")
print(cc_pairwise)
if (nrow(outlier_log) > 0) {
  cat("\n=== OUTLIERS REMOVED (TP_mg_L shown alongside response) ===\n")
  print(outlier_log, row.names = FALSE)
}
cat("\nSaved: 3_models/tp_submodel.rds |",
    "2_incremental/tp_submodel_predictions.csv |",
    "4_products/diagnostics/tp_submodel.pdf\n")