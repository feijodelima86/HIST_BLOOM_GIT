## ==========================================================================
## 07d_gam_chla_NP_ratio.R
## --------------------------------------------------------------------------
## Purpose : Add TN:TP ratio as a predictor and rerun full GAMs with
##           shrinkage selection for two site configurations:
##             Config A: 6 sites (drop FH only)
##             Config B: 5 sites (drop FH + HU)
##           Both include the AIC single-term-drop permutation.
##
## Notes   : TN_mg_L and TP_mg_L are already as N and P, so simple
##           division gives molar-equivalent ratio.  Raw ratio used
##           (not log) to avoid perfect collinearity with logTN - logTP.
##
## Inputs  : 2_incremental/ucfr_model_ready.csv
## Outputs : Console diagnostics, base-R plots
## ==========================================================================

library(mgcv)

# --------------------------------------------------------------------------
# 1.  Load & prepare
# --------------------------------------------------------------------------

dat <- read.csv("2_incremental/ucfr_model_ready.csv", stringsAsFactors = FALSE)

dat$logTP   <- log10(dat$TP_mg_L)
dat$logTN   <- log10(dat$TN_mg_L)
dat$NP_ratio <- dat$TN_mg_L / dat$TP_mg_L

cat("NP_ratio summary:\n")
print(summary(dat$NP_ratio))
cat("\n")

# Predictors now include NP_ratio
predictors <- c("anomaly", "Q_obs_cfs", "Temp_oC",
                "Days_Since_Freshet", "SPC", "logTP", "logTN", "NP_ratio")

vars <- c("logCHLa", predictors)

# ==========================================================================
# Helper function: fit GAM, LOSO, temporal, drop-AIC for a given subset
# ==========================================================================

run_gam_battery <- function(dat_sub, config_label, site_levels) {
  
  dat_sub$Site <- factor(dat_sub$Site, levels = site_levels)
  dat_sub <- dat_sub[complete.cases(dat_sub[, vars]), ]
  
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("  CONFIG: ", config_label, "\n")
  cat(strrep("=", 70), "\n")
  cat("Complete cases:", nrow(dat_sub), "\n")
  cat("Sites:", paste(levels(droplevels(dat_sub$Site)), collapse = ", "), "\n")
  cat("Years:", paste(range(dat_sub$Year), collapse = "–"), "\n")
  cat("Obs per site:\n")
  print(table(dat_sub$Site))
  cat("\n")
  
  # --- Correlations ---
  cat("=== Pairwise correlations ===\n")
  cor_mat <- cor(dat_sub[, predictors], use = "complete.obs")
  print(round(cor_mat, 2))
  
  high <- which(abs(cor_mat) > 0.7 & upper.tri(cor_mat), arr.ind = TRUE)
  if (nrow(high) > 0) {
    cat("\nPairs with |r| > 0.7:\n")
    for (i in seq_len(nrow(high))) {
      cat("  ", predictors[high[i, 1]], " – ", predictors[high[i, 2]],
          " : r =", round(cor_mat[high[i, 1], high[i, 2]], 3), "\n")
    }
  } else {
    cat("\nNo predictor pairs exceed |r| = 0.7\n")
  }
  cat("\n")
  
  # --- Fit full GAM ---
  m1 <- gam(logCHLa ~ s(anomaly,            bs = "ts", k = 5) +
              s(Q_obs_cfs,          bs = "ts", k = 5) +
              s(Temp_oC,            bs = "ts", k = 5) +
              s(Days_Since_Freshet, bs = "ts", k = 5) +
              s(SPC,                bs = "ts", k = 5) +
              s(logTP,              bs = "ts", k = 5) +
              s(logTN,              bs = "ts", k = 5) +
              s(NP_ratio,           bs = "ts", k = 5),
            data   = dat_sub,
            method = "REML",
            select = TRUE)
  
  cat("=== GAM summary ===\n")
  print(summary(m1))
  
  # --- Concurvity ---
  cat("\n=== Concurvity (worst case) ===\n")
  print(round(concurvity(m1, full = TRUE), 3))
  
  cat("\n=== Concurvity (pairwise, estimate) ===\n")
  cc <- concurvity(m1, full = FALSE)
  print(round(cc$estimate, 3))
  
  # --- gam.check ---
  par(mfrow = c(2, 2))
  gam.check(m1)
  title(sub = config_label, line = -1, outer = TRUE)
  
  # --- Smooth plots ---
  par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))
  plot(m1, shade = TRUE, shade.col = "lightblue",
       residuals = TRUE, pch = 16, cex = 0.4, col = "grey40",
       pages = 0)
  mtext(config_label, side = 3, line = -1.5, outer = TRUE, cex = 0.9)
  
  # --- LOSO ---
  sites <- levels(droplevels(dat_sub$Site))
  loso_preds <- data.frame()
  
  for (s in sites) {
    train <- dat_sub[dat_sub$Site != s, ]
    test  <- dat_sub[dat_sub$Site == s, ]
    
    m_cv <- gam(logCHLa ~ s(anomaly,            bs = "ts", k = 5) +
                  s(Q_obs_cfs,          bs = "ts", k = 5) +
                  s(Temp_oC,            bs = "ts", k = 5) +
                  s(Days_Since_Freshet, bs = "ts", k = 5) +
                  s(SPC,                bs = "ts", k = 5) +
                  s(logTP,              bs = "ts", k = 5) +
                  s(logTN,              bs = "ts", k = 5) +
                  s(NP_ratio,           bs = "ts", k = 5),
                data   = train,
                method = "REML",
                select = TRUE)
    
    test$pred <- predict(m_cv, newdata = test, type = "response")
    loso_preds <- rbind(loso_preds,
                        data.frame(Site = test$Site, Year = test$Year,
                                   obs = test$logCHLa, pred = test$pred))
  }
  
  ss_res <- sum((loso_preds$obs - loso_preds$pred)^2)
  ss_tot <- sum((loso_preds$obs - mean(loso_preds$obs))^2)
  loso_r2 <- 1 - ss_res / ss_tot
  
  cat("\n=== LOSO cross-validation ===\n")
  cat("Overall LOSO R²:", round(loso_r2, 3), "\n")
  cat("Overall LOSO RMSE:",
      round(sqrt(mean((loso_preds$obs - loso_preds$pred)^2)), 3), "\n\n")
  
  cat("Per-site LOSO:\n")
  for (s in sites) {
    sub <- loso_preds[loso_preds$Site == s, ]
    if (nrow(sub) < 3) { cat("  ", s, ": n <3, skipped\n"); next }
    r2_s <- 1 - sum((sub$obs - sub$pred)^2) / sum((sub$obs - mean(sub$obs))^2)
    cat(sprintf("  %s : n = %3d, R² = %6.3f, RMSE = %.3f\n",
                s, nrow(sub), r2_s,
                sqrt(mean((sub$obs - sub$pred)^2))))
  }
  
  # --- Temporal jackknife ---
  years <- sort(unique(dat_sub$Year))
  temp_preds <- data.frame()
  
  for (y in years) {
    train <- dat_sub[dat_sub$Year != y, ]
    test  <- dat_sub[dat_sub$Year == y, ]
    
    m_ty <- gam(logCHLa ~ s(anomaly,            bs = "ts", k = 5) +
                  s(Q_obs_cfs,          bs = "ts", k = 5) +
                  s(Temp_oC,            bs = "ts", k = 5) +
                  s(Days_Since_Freshet, bs = "ts", k = 5) +
                  s(SPC,                bs = "ts", k = 5) +
                  s(logTP,              bs = "ts", k = 5) +
                  s(logTN,              bs = "ts", k = 5) +
                  s(NP_ratio,           bs = "ts", k = 5),
                data   = train,
                method = "REML",
                select = TRUE)
    
    test$pred <- predict(m_ty, newdata = test, type = "response")
    temp_preds <- rbind(temp_preds,
                        data.frame(Site = test$Site, Year = test$Year,
                                   obs = test$logCHLa, pred = test$pred))
  }
  
  ss_res_t <- sum((temp_preds$obs - temp_preds$pred)^2)
  ss_tot_t <- sum((temp_preds$obs - mean(temp_preds$obs))^2)
  temp_r2  <- 1 - ss_res_t / ss_tot_t
  
  cat("\n=== Temporal jackknife ===\n")
  cat("Overall temporal R²:", round(temp_r2, 3), "\n")
  cat("Overall temporal RMSE:",
      round(sqrt(mean((temp_preds$obs - temp_preds$pred)^2)), 3), "\n")
  
  # --- Single-term drop AIC ---
  cat("\n=== Single-term drop AIC ===\n")
  cat(sprintf("Full model AIC: %.1f\n", AIC(m1)))
  
  for (d in predictors) {
    f_drop <- as.formula(
      paste("logCHLa ~",
            paste(sprintf('s(%s, bs = "ts", k = 5)',
                          setdiff(predictors, d)),
                  collapse = " + ")))
    m_drop <- gam(f_drop, data = dat_sub, method = "REML", select = TRUE)
    cat(sprintf("  Drop %-20s  AIC = %7.1f  (delta = %+.1f)\n",
                d, AIC(m_drop), AIC(m_drop) - AIC(m1)))
  }
  
  # --- Return key metrics for comparison ---
  invisible(list(loso_r2 = loso_r2, temp_r2 = temp_r2,
                 aic = AIC(m1), model = m1,
                 loso_preds = loso_preds))
}

# ==========================================================================
# 2.  Config A: 6 sites (drop FH only)
# ==========================================================================

dat_A <- dat
res_A <- run_gam_battery(dat_A, "FULL",
                         c("DL","GR","BN","MS","BM","HU","FH"))

# ==========================================================================
# 3.  Config B: 5 sites (drop FH + HU)
# ==========================================================================

dat_B <- dat[!(dat$Site %in% c("FH", "HU")), ]
res_B <- run_gam_battery(dat_B, "5 sites (no FH, no HU)",
                         c("DL","GR","BN","MS","BM"))

# ==========================================================================
# 4.  Summary comparison
# ==========================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("  CROSS-CONFIGURATION COMPARISON\n")
cat(strrep("=", 70), "\n\n")
cat("                     7-site(07a)  7-site red(07b)  6-site(A)  5-site(B)\n")
cat(sprintf("  LOSO R²            0.451        0.426            %.3f      %.3f\n",
            res_A$loso_r2, res_B$loso_r2))
cat(sprintf("  Temporal R²        0.566        0.567            %.3f      %.3f\n",
            res_A$temp_r2, res_B$temp_r2))
cat("  BRT LOSO (ref)     0.279\n")
cat("  BRT Temporal (ref) 0.404\n")

cat("\n--- Done. Compare EDF tables and drop-AIC across configs. ---\n")
cat("--- Key question: does NP_ratio survive shrinkage where logTN did not? ---\n")