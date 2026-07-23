## ==========================================================================
## 07f_gam_submodels_no_baseflow.R
## --------------------------------------------------------------------------
## Purpose : Rerun GAM shrinkage selection for SPC and logTP submodels
##           after removing Q_baseflow_cfs from candidates.  Baseflow is
##           already embedded in anomaly = (peak/baseflow)^(1/3), so
##           including it separately is circular.
##
##           Candidate predictors:
##             anomaly, Q_obs_cfs, Temp_oC, Days_Since_Freshet, Q_peak_cfs
##
##           All 7 sites retained.
##
## Inputs  : 2_incremental/ucfr_model_ready.csv
## Outputs : Console diagnostics, base-R plots
## ==========================================================================

library(mgcv)

# --------------------------------------------------------------------------
# 1.  Load & prepare
# --------------------------------------------------------------------------

dat <- read.csv("2_incremental/ucfr_model_ready.csv", stringsAsFactors = FALSE)

dat$logTP <- log10(dat$TP_mg_L)
dat$Site  <- factor(dat$Site,
                    levels = c("DL","GR","BN","MS","BM","HU","FH"))

# Candidate predictors — no Q_baseflow_cfs
hydro_preds <- c("anomaly", "Q_obs_cfs", "Temp_oC",
                 "Days_Since_Freshet", "Q_peak_cfs")

# ==========================================================================
# Helper (same as 07e but reproduced for standalone use)
# ==========================================================================

run_submodel_battery <- function(dat, response, response_label, predictors) {
  
  vars <- c(response, predictors)
  dat_mod <- dat[complete.cases(dat[, vars]), ]
  
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("  SUBMODEL: ", response_label, "\n")
  cat(strrep("=", 70), "\n")
  cat("Complete cases:", nrow(dat_mod), "\n")
  cat("Sites:", paste(levels(droplevels(dat_mod$Site)), collapse = ", "), "\n")
  cat("Years:", paste(range(dat_mod$Year), collapse = "–"), "\n\n")
  
  cat("Response summary:\n")
  print(summary(dat_mod[[response]]))
  cat("\n")
  
  # --- Correlations ---
  cat("=== Pairwise correlations among predictors ===\n")
  cor_mat <- cor(dat_mod[, predictors], use = "complete.obs")
  print(round(cor_mat, 2))
  
  high <- which(abs(cor_mat) > 0.7 & upper.tri(cor_mat), arr.ind = TRUE)
  if (nrow(high) > 0) {
    cat("\nPairs with |r| > 0.7:\n")
    for (i in seq_len(nrow(high))) {
      cat("  ", predictors[high[i,1]], " – ", predictors[high[i,2]],
          " : r =", round(cor_mat[high[i,1], high[i,2]], 3), "\n")
    }
  } else {
    cat("\nNo predictor pairs exceed |r| = 0.7\n")
  }
  cat("\n")
  
  cat("=== Correlation with response ===\n")
  r_resp <- sapply(predictors, function(p) {
    cor(dat_mod[[p]], dat_mod[[response]], use = "complete.obs")
  })
  print(round(sort(abs(r_resp), decreasing = TRUE), 3))
  cat("\n")
  
  # --- Fit full GAM ---
  smooth_terms <- paste(sprintf('s(%s, bs = "ts", k = 5)', predictors),
                        collapse = " + ")
  f <- as.formula(paste(response, "~", smooth_terms))
  
  m1 <- gam(f, data = dat_mod, method = "REML", select = TRUE)
  
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
  
  # --- Smooth plots ---
  n_preds <- length(predictors)
  nr <- ceiling(n_preds / 3)
  par(mfrow = c(nr, 3), mar = c(4, 4, 2, 1))
  plot(m1, shade = TRUE, shade.col = "lightblue",
       residuals = TRUE, pch = 16, cex = 0.4, col = "grey40",
       pages = 0)
  mtext(response_label, side = 3, line = -1.5, outer = TRUE, cex = 0.9)
  
  # --- LOSO ---
  sites <- levels(droplevels(dat_mod$Site))
  loso_preds <- data.frame()
  
  for (s in sites) {
    train <- dat_mod[dat_mod$Site != s, ]
    test  <- dat_mod[dat_mod$Site == s, ]
    
    m_cv <- gam(f, data = train, method = "REML", select = TRUE)
    test$pred <- predict(m_cv, newdata = test, type = "response")
    
    loso_preds <- rbind(loso_preds,
                        data.frame(Site = test$Site, Year = test$Year,
                                   obs = test[[response]], pred = test$pred))
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
  years <- sort(unique(dat_mod$Year))
  temp_preds <- data.frame()
  
  for (y in years) {
    train <- dat_mod[dat_mod$Year != y, ]
    test  <- dat_mod[dat_mod$Year == y, ]
    
    m_ty <- gam(f, data = train, method = "REML", select = TRUE)
    test$pred <- predict(m_ty, newdata = test, type = "response")
    
    temp_preds <- rbind(temp_preds,
                        data.frame(Site = test$Site, Year = test$Year,
                                   obs = test[[response]], pred = test$pred))
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
      paste(response, "~",
            paste(sprintf('s(%s, bs = "ts", k = 5)',
                          setdiff(predictors, d)),
                  collapse = " + ")))
    m_drop <- gam(f_drop, data = dat_mod, method = "REML", select = TRUE)
    cat(sprintf("  Drop %-20s  AIC = %7.1f  (delta = %+.1f)\n",
                d, AIC(m_drop), AIC(m_drop) - AIC(m1)))
  }
  
  # --- LOSO plot ---
  par(mfrow = c(1, 1), mar = c(4.5, 4.5, 2, 1))
  site_cols <- setNames(
    c("firebrick","darkorange","forestgreen","steelblue",
      "gold3","mediumpurple","grey40"),
    c("DL","GR","BN","MS","BM","HU","FH"))
  
  plot(loso_preds$obs, loso_preds$pred,
       col  = site_cols[as.character(loso_preds$Site)],
       pch  = 16, cex = 0.9,
       xlab = paste("Observed", response_label),
       ylab = paste("LOSO Predicted", response_label),
       main = sprintf("%s LOSO  (R² = %.3f)", response_label, loso_r2))
  abline(0, 1, lty = 2)
  legend("topleft", legend = names(site_cols), col = site_cols,
         pch = 16, cex = 0.7, bty = "n", ncol = 2)
  
  # --- Residuals by site ---
  loso_preds$resid <- loso_preds$obs - loso_preds$pred
  boxplot(resid ~ Site, data = loso_preds,
          col = site_cols[levels(droplevels(loso_preds$Site))],
          xlab = "Site", ylab = "LOSO Residual",
          main = paste(response_label, ": LOSO residuals by site"))
  abline(h = 0, lty = 2)
  
  invisible(list(loso_r2 = loso_r2, temp_r2 = temp_r2,
                 aic = AIC(m1), model = m1,
                 loso_preds = loso_preds, data = dat_mod))
}

# ==========================================================================
# 2.  SPC submodel
# ==========================================================================

res_spc <- run_submodel_battery(dat, "SPC", "SPC (Specific Conductance)",
                                hydro_preds)

# ==========================================================================
# 3.  logTP submodel
# ==========================================================================

res_tp <- run_submodel_battery(dat, "logTP", "log10 TP",
                               hydro_preds)

# ==========================================================================
# 4.  Summary
# ==========================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("  SUBMODEL SUMMARY (no Q_baseflow_cfs)\n")
cat(strrep("=", 70), "\n\n")
cat("                     SPC          logTP\n")
cat(sprintf("  LOSO R²            %.3f        %.3f\n",
            res_spc$loso_r2, res_tp$loso_r2))
cat(sprintf("  Temporal R²        %.3f        %.3f\n",
            res_spc$temp_r2, res_tp$temp_r2))
cat(sprintf("  Full model AIC     %.1f       %.1f\n",
            res_spc$aic, res_tp$aic))

cat("\n--- Compare to 07e (with Q_baseflow_cfs): ---\n")
cat("  07e SPC:   LOSO = 0.299, Temporal = 0.827\n")
cat("  07e logTP: LOSO = -0.246, Temporal = 0.425\n")

cat("\n--- Review drop-AIC and decide reduced submodels. ---\n")