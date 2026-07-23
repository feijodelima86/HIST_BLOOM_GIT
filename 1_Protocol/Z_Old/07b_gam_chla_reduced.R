## ==========================================================================
## 07b_gam_chla_reduced.R
## --------------------------------------------------------------------------
## Purpose : Refit GAM with only the 4 predictors that survived shrinkage
##           selection in 07a: SPC, anomaly, logTP, Days_Since_Freshet.
##           Compare LOSO and temporal jackknife to full 7-predictor model.
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
                    levels = c("DL","GR","BM","BN","MS","HU","FH"))

vars <- c("logCHLa", "anomaly", "Days_Since_Freshet", "SPC", "logTP")
dat_mod <- dat[complete.cases(dat[, vars]), ]

cat("Complete cases:", nrow(dat_mod), "of", nrow(dat), "\n")
cat("Sites:", paste(levels(droplevels(dat_mod$Site)), collapse = ", "), "\n")
cat("Years:", paste(range(dat_mod$Year), collapse = "–"), "\n\n")

# --------------------------------------------------------------------------
# 2.  Fit reduced GAM
# --------------------------------------------------------------------------

m_red <- gam(logCHLa ~ s(SPC,                bs = "ts", k = 5) +
               s(anomaly,            bs = "ts", k = 5) +
               s(logTP,              bs = "ts", k = 5) +
               s(Days_Since_Freshet, bs = "ts", k = 5),
             data   = dat_mod,
             method = "REML",
             select = TRUE)

cat("=== Reduced GAM summary ===\n")
print(summary(m_red))

# --------------------------------------------------------------------------
# 3.  Diagnostics
# --------------------------------------------------------------------------

cat("\n=== Concurvity (worst case) ===\n")
print(round(concurvity(m_red, full = TRUE), 3))

cat("\n=== Concurvity (pairwise, estimate) ===\n")
cc <- concurvity(m_red, full = FALSE)
print(round(cc$estimate, 3))

par(mfrow = c(2, 2))
gam.check(m_red)

# --------------------------------------------------------------------------
# 4.  Smooth partial effects
# --------------------------------------------------------------------------

par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(m_red, shade = TRUE, shade.col = "lightblue",
     residuals = TRUE, pch = 16, cex = 0.4, col = "grey40",
     pages = 0)

# --------------------------------------------------------------------------
# 5.  LOSO cross-validation
# --------------------------------------------------------------------------

sites <- levels(dat_mod$Site)
loso_preds <- data.frame()

for (s in sites) {
  train <- dat_mod[dat_mod$Site != s, ]
  test  <- dat_mod[dat_mod$Site == s, ]
  
  m_cv <- gam(logCHLa ~ s(SPC,                bs = "ts", k = 5) +
                s(anomaly,            bs = "ts", k = 5) +
                s(logTP,              bs = "ts", k = 5) +
                s(Days_Since_Freshet, bs = "ts", k = 5),
              data   = train,
              method = "REML",
              select = TRUE)
  
  test$pred <- predict(m_cv, newdata = test, type = "response")
  
  loso_preds <- rbind(loso_preds,
                      data.frame(Site = test$Site,
                                 Year = test$Year,
                                 obs  = test$logCHLa,
                                 pred = test$pred))
}

ss_res <- sum((loso_preds$obs - loso_preds$pred)^2)
ss_tot <- sum((loso_preds$obs - mean(loso_preds$obs))^2)
loso_r2 <- 1 - ss_res / ss_tot

cat("\n=== LOSO cross-validation (reduced model) ===\n")
cat("Overall LOSO R²:", round(loso_r2, 3), "\n")
cat("Overall LOSO RMSE:", round(sqrt(mean((loso_preds$obs - loso_preds$pred)^2)), 3), "\n\n")

cat("Per-site LOSO:\n")
for (s in sites) {
  sub <- loso_preds[loso_preds$Site == s, ]
  if (nrow(sub) < 3) { cat("  ", s, ": n <3, skipped\n"); next }
  r2_site <- 1 - sum((sub$obs - sub$pred)^2) / sum((sub$obs - mean(sub$obs))^2)
  cat(sprintf("  %s : n = %3d, R² = %6.3f, RMSE = %.3f\n",
              s, nrow(sub), r2_site,
              sqrt(mean((sub$obs - sub$pred)^2))))
}

# --------------------------------------------------------------------------
# 6.  Temporal leave-one-year-out
# --------------------------------------------------------------------------

years <- sort(unique(dat_mod$Year))
temp_preds <- data.frame()

for (y in years) {
  train <- dat_mod[dat_mod$Year != y, ]
  test  <- dat_mod[dat_mod$Year == y, ]
  
  m_ty <- gam(logCHLa ~ s(SPC,                bs = "ts", k = 5) +
                s(anomaly,            bs = "ts", k = 5) +
                s(logTP,              bs = "ts", k = 5) +
                s(Days_Since_Freshet, bs = "ts", k = 5),
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

cat("\n=== Temporal jackknife (reduced model) ===\n")
cat("Overall temporal R²:", round(temp_r2, 3), "\n")
cat("Overall temporal RMSE:", round(sqrt(mean((temp_preds$obs - temp_preds$pred)^2)), 3), "\n")

# --------------------------------------------------------------------------
# 7.  Comparison table
# --------------------------------------------------------------------------

cat("\n=== Model comparison ===\n")
cat("                     Full (7 pred)    Reduced (4 pred)\n")
cat(sprintf("  LOSO R²            0.451            %.3f\n", loso_r2))
cat(sprintf("  Temporal R²        0.566            %.3f\n", temp_r2))
cat(sprintf("  BRT LOSO (ref)     0.279\n"))
cat(sprintf("  BRT Temporal (ref) 0.404\n"))

# --------------------------------------------------------------------------
# 8.  LOSO obs vs pred plot
# --------------------------------------------------------------------------

par(mfrow = c(1, 1), mar = c(4.5, 4.5, 2, 1))
site_cols <- setNames(
  c("firebrick","darkorange","gold3","forestgreen",
    "steelblue","mediumpurple","grey40"),
  c("DL","GR","BM","BN","MS","HU","FH"))

plot(loso_preds$obs, loso_preds$pred,
     col  = site_cols[as.character(loso_preds$Site)],
     pch  = 16, cex = 0.9,
     xlab = "Observed log10 CHLa",
     ylab = "LOSO Predicted log10 CHLa",
     main = sprintf("Reduced GAM LOSO  (R² = %.3f)", loso_r2))
abline(0, 1, lty = 2)
abline(h = 2, v = 2, col = "red", lty = 3)
legend("topleft", legend = names(site_cols), col = site_cols,
       pch = 16, cex = 0.7, bty = "n", ncol = 2)

# --------------------------------------------------------------------------
# 9.  Residuals by site (boxplot)
# --------------------------------------------------------------------------

loso_preds$resid <- loso_preds$obs - loso_preds$pred

par(mfrow = c(1, 1), mar = c(4.5, 4.5, 2, 1))
boxplot(resid ~ Site, data = loso_preds,
        col = site_cols[levels(loso_preds$Site)],
        xlab = "Site", ylab = "LOSO Residual (obs – pred)",
        main = "Reduced GAM: LOSO residuals by site")
abline(h = 0, lty = 2)

cat("\n--- Done. Compare to 07a results. ---\n")