## ==========================================================================
## 07a_gam_chla_feasibility.R
## --------------------------------------------------------------------------
## Purpose : Can a GAM predict log10 CHLa from observed environmental
##           predictors?  If not, there is no point building the full
##           cascade.  Uses shrinkage smooths (select = TRUE, bs = "ts")
##           so the fitting procedure itself performs variable selection.
##
## Inputs  : 2_incremental/ucfr_model_ready.csv
## Outputs : Console diagnostics, base-R plots
##
## Notes   : Binomial nuisance-bloom model (threshold 100 mg/m²,
##           logCHLa >= 2) kept in back pocket — flagged but not run.
## ==========================================================================

library(mgcv)

# --------------------------------------------------------------------------
# 1.  Load & prepare
# --------------------------------------------------------------------------

dat <- read.csv("2_incremental/ucfr_model_ready.csv", stringsAsFactors = FALSE)

# Log-transform TP and TN (matching pipeline convention)
dat$logTP <- log10(dat$TP_mg_L)
dat$logTN <- log10(dat$TN_mg_L)

# Site as factor (ordered longitudinally for plotting convenience)
dat$Site <- factor(dat$Site,
                   levels = c("DL","GR","BM","BN","MS","HU","FH"))

# Drop rows with any NA in the predictor / response set
vars <- c("logCHLa", "anomaly", "Q_obs_cfs", "Temp_oC",
          "Days_Since_Freshet", "SPC", "logTP", "logTN")
dat_mod <- dat[complete.cases(dat[, vars]), ]

cat("Complete cases:", nrow(dat_mod), "of", nrow(dat), "\n")
cat("Sites represented:", paste(levels(droplevels(dat_mod$Site)),
                                collapse = ", "), "\n")
cat("Year range:", range(dat_mod$Year), "\n\n")

# --------------------------------------------------------------------------
# 2.  Pairwise correlations & concurvity preview
# --------------------------------------------------------------------------

predictors <- c("anomaly","Q_obs_cfs","Temp_oC",
                "Days_Since_Freshet","SPC","logTP","logTN")

cat("=== Pairwise Pearson correlations among predictors ===\n")
cor_mat <- cor(dat_mod[, predictors], use = "complete.obs")
print(round(cor_mat, 2))
cat("\n")

# Flag pairs with |r| > 0.7
high <- which(abs(cor_mat) > 0.7 & upper.tri(cor_mat), arr.ind = TRUE)
if (nrow(high) > 0) {
  cat("Pairs with |r| > 0.7:\n")
  for (i in seq_len(nrow(high))) {
    cat("  ", predictors[high[i,1]], " – ", predictors[high[i,2]],
        " : r =", round(cor_mat[high[i,1], high[i,2]], 3), "\n")
  }
} else {
  cat("No predictor pairs exceed |r| = 0.7\n")
}
cat("\n")

# --------------------------------------------------------------------------
# 3.  Fit GAM — shrinkage selection
# --------------------------------------------------------------------------

#  bs = "ts"  : thin-plate with shrinkage (penalises null space)
#  select = TRUE : extra penalty allows full removal of terms
#  method = "REML" : preferred for penalised selection
#  k = 5 : conservative basis dimension given ~N site-years

m1 <- gam(logCHLa ~ s(anomaly,            bs = "ts", k = 5) +
            s(Q_obs_cfs,          bs = "ts", k = 5) +
            s(Temp_oC,            bs = "ts", k = 5) +
            s(Days_Since_Freshet, bs = "ts", k = 5) +
            s(SPC,                bs = "ts", k = 5) +
            s(logTP,              bs = "ts", k = 5) +
            s(logTN,              bs = "ts", k = 5),
          data   = dat_mod,
          method = "REML",
          select = TRUE)

cat("=== GAM summary ===\n")
print(summary(m1))

# --------------------------------------------------------------------------
# 4.  Concurvity (GAM analog of multicollinearity)
# --------------------------------------------------------------------------

cat("\n=== Concurvity (worst case) ===\n")
print(round(concurvity(m1, full = TRUE), 3))

cat("\n=== Concurvity (pairwise, estimate) ===\n")
cc <- concurvity(m1, full = FALSE)
print(round(cc$estimate, 3))

# --------------------------------------------------------------------------
# 5.  Diagnostic plots
# --------------------------------------------------------------------------

par(mfrow = c(2, 2))
gam.check(m1)

# --------------------------------------------------------------------------
# 6.  Smooth partial effects
# --------------------------------------------------------------------------

par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))
plot(m1, shade = TRUE, shade.col = "lightblue",
     residuals = TRUE, pch = 16, cex = 0.4, col = "grey40",
     pages = 0)

# --------------------------------------------------------------------------
# 7.  LOSO cross-validation  (the key test)
# --------------------------------------------------------------------------

sites <- levels(dat_mod$Site)
loso_preds <- data.frame()

for (s in sites) {
  train <- dat_mod[dat_mod$Site != s, ]
  test  <- dat_mod[dat_mod$Site == s, ]
  
  m_cv <- gam(logCHLa ~ s(anomaly,            bs = "ts", k = 5) +
                s(Q_obs_cfs,          bs = "ts", k = 5) +
                s(Temp_oC,            bs = "ts", k = 5) +
                s(Days_Since_Freshet, bs = "ts", k = 5) +
                s(SPC,                bs = "ts", k = 5) +
                s(logTP,              bs = "ts", k = 5) +
                s(logTN,              bs = "ts", k = 5),
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

# Overall LOSO R²
ss_res <- sum((loso_preds$obs - loso_preds$pred)^2)
ss_tot <- sum((loso_preds$obs - mean(loso_preds$obs))^2)
loso_r2 <- 1 - ss_res / ss_tot

cat("\n=== LOSO cross-validation ===\n")
cat("Overall LOSO R²:", round(loso_r2, 3), "\n")
cat("Overall LOSO RMSE:", round(sqrt(mean((loso_preds$obs - loso_preds$pred)^2)), 3), "\n\n")

# Per-site LOSO
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
# 8.  LOSO observed vs predicted plot
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
     main = sprintf("GAM LOSO  (R² = %.3f)", loso_r2))
abline(0, 1, lty = 2)
# Nuisance bloom threshold for reference
abline(h = 2, v = 2, col = "red", lty = 3)
legend("topleft", legend = names(site_cols), col = site_cols,
       pch = 16, cex = 0.7, bty = "n", ncol = 2)

# --------------------------------------------------------------------------
# 9.  Temporal leave-one-year-out (secondary check)
# --------------------------------------------------------------------------

years <- sort(unique(dat_mod$Year))
temp_preds <- data.frame()

for (y in years) {
  train <- dat_mod[dat_mod$Year != y, ]
  test  <- dat_mod[dat_mod$Year == y, ]
  
  m_ty <- gam(logCHLa ~ s(anomaly,            bs = "ts", k = 5) +
                s(Q_obs_cfs,          bs = "ts", k = 5) +
                s(Temp_oC,            bs = "ts", k = 5) +
                s(Days_Since_Freshet, bs = "ts", k = 5) +
                s(SPC,                bs = "ts", k = 5) +
                s(logTP,              bs = "ts", k = 5) +
                s(logTN,              bs = "ts", k = 5),
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
cat("Overall temporal RMSE:", round(sqrt(mean((temp_preds$obs - temp_preds$pred)^2)), 3), "\n")

# --------------------------------------------------------------------------
# 10. Quick benchmark: how much does each predictor contribute?
# --------------------------------------------------------------------------

#  Refit dropping one smooth at a time; compare AIC
cat("\n=== Single-term drop AIC ===\n")
cat(sprintf("Full model AIC: %.1f\n", AIC(m1)))

drop_terms <- predictors
for (d in drop_terms) {
  f_drop <- as.formula(
    paste("logCHLa ~",
          paste(sprintf('s(%s, bs = "ts", k = 5)',
                        setdiff(predictors, d)),
                collapse = " + ")))
  m_drop <- gam(f_drop, data = dat_mod, method = "REML", select = TRUE)
  cat(sprintf("  Drop %-20s  AIC = %7.1f  (delta = %+.1f)\n",
              d, AIC(m_drop), AIC(m_drop) - AIC(m1)))
}

# ==========================================================================
# BACK POCKET: Binomial nuisance-bloom model
# --------------------------------------------------------------------------
# bloom_yn <- ifelse(dat_mod$logCHLa >= 2, 1, 0)
# m_binom <- gam(bloom_yn ~ s(anomaly, bs="ts", k=5) + ... ,
#                family = binomial(link = "logit"),
#                data = dat_mod, method = "REML", select = TRUE)
# ==========================================================================

cat("\n--- Done. Review EDF table, LOSO R², and smooth plots. ---\n")
cat("--- If LOSO R² substantially beats 0.279, GAMs are worth pursuing. ---\n")

