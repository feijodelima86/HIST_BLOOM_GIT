# ============================================================================
# 09_2_sem_temporal_jackknife.R
# UCFR Filamentous Algae Project
# Stage 8.2: Temporal jackknife validation of full BRT-SEM chain
#
# Input:    2_incremental/ucfr_model_ready.csv
# Outputs:  2_incremental/sem_jackknife_predictions.csv
#           2_incremental/sem_jackknife_performance.csv
#           4_products/diagnostics/sem_jackknife_fit.pdf
#
# Procedure:
#   For each year in the record:
#     1. Withhold all observations from that year (all sites)
#     2. Refit all four submodels on remaining years
#     3. Predict withheld year through full chain
#     4. Compare to observed logCHLa
#
# Full chain (no observed chemistry):
#   hydrology -> SPC
#   hydrology -> log10(TP)
#   hydrology -> log10(TN)
#   pred_SPC + pred_logTP + pred_logTN + hydrology -> log10(CHLa)
#
# Notes:
#   - All sites retained in both train and test — only year withheld
#   - No site identity leakage issue — cleaner than LOSO
#   - Tests temporal generalization rather than spatial
#   - Directly relevant for climate scenario projections into future years
#   - BRT settings: tc=3 for SPC/TP/TN, tc=4 for bloom
# ============================================================================

library(readr)
library(dismo)
library(gbm)

# ----------------------------------------------------------------------------
# 1. Read and prepare data
# ----------------------------------------------------------------------------

cat("Reading model-ready dataset...\n")
dat <- as.data.frame(
  read_csv("2_incremental/ucfr_model_ready.csv", show_col_types = FALSE)
)

dat$logCHLa   <- log10(dat$CHLa)
dat$logTP_obs <- log10(dat$TP_mg_L)
dat$logTN_obs <- log10(dat$TN_mg_L)

for (d in c("2_incremental", "4_products/diagnostics")) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

hydro_preds <- c("anomaly", "Q_obs_cfs", "Temp_oC", "Days_Since_Freshet")
bloom_preds <- c("pred_SPC", "pred_logTP", "pred_logTN", hydro_preds)

chain_vars  <- c("Site", "Year", "Month",
                 hydro_preds, "SPC", "logTP_obs", "logTN_obs", "logCHLa")

dat_model <- as.data.frame(dat[complete.cases(dat[ , chain_vars]), ])

years <- sort(unique(dat_model$Year))
cat(sprintf("Complete cases: %d rows\n", nrow(dat_model)))
cat(sprintf("Years in record: %d (%d to %d)\n\n",
            length(years), min(years), max(years)))

# ----------------------------------------------------------------------------
# 2. BRT fitting function
# ----------------------------------------------------------------------------

fit_brt <- function(data, response, predictors, tc, lr = 0.01, bag = 0.75) {
  data_sub <- as.data.frame(data[ , c(predictors, response), drop = FALSE])
  tryCatch(
    gbm.step(
      data            = data_sub,
      gbm.x           = which(names(data_sub) %in% predictors),
      gbm.y           = which(names(data_sub) == response),
      family          = "gaussian",
      tree.complexity = tc,
      learning.rate   = lr,
      bag.fraction    = bag,
      verbose         = FALSE,
      plot.main       = FALSE
    ),
    error = function(e) {
      cat(sprintf("    FIT FAILED: %s\n", conditionMessage(e)))
      NULL
    }
  )
}

# ----------------------------------------------------------------------------
# 3. Temporal jackknife loop
# ----------------------------------------------------------------------------

cat("Running temporal jackknife...\n")
cat(sprintf("(%d iterations — one per year)\n\n", length(years)))

jk_out <- vector("list", length(years))

for (i in seq_along(years)) {
  yr    <- years[i]
  train <- as.data.frame(dat_model[dat_model$Year != yr, ])
  test  <- as.data.frame(dat_model[dat_model$Year == yr, ])
  
  cat(sprintf("--- Leaving out %d (train n=%d, test n=%d) ---\n",
              yr, nrow(train), nrow(test)))
  
  # -- Step 1: fit intermediate submodels ----------------------------------
  brt_SPC_jk <- fit_brt(train, "SPC",       hydro_preds, tc = 3)
  brt_TP_jk  <- fit_brt(train, "logTP_obs", hydro_preds, tc = 3)
  brt_TN_jk  <- fit_brt(train, "logTN_obs", hydro_preds, tc = 3)
  
  if (any(sapply(list(brt_SPC_jk, brt_TP_jk, brt_TN_jk), is.null))) {
    cat(sprintf("  SKIPPING %d — intermediate model fit failed\n\n", yr))
    next
  }
  
  # -- Step 2: predict intermediates for withheld year ---------------------
  test$pred_SPC   <- predict(brt_SPC_jk, newdata = test[ , hydro_preds],
                             n.trees = brt_SPC_jk$gbm.call$best.trees)
  test$pred_logTP <- predict(brt_TP_jk,  newdata = test[ , hydro_preds],
                             n.trees = brt_TP_jk$gbm.call$best.trees)
  test$pred_logTN <- predict(brt_TN_jk,  newdata = test[ , hydro_preds],
                             n.trees = brt_TN_jk$gbm.call$best.trees)
  
  # -- Step 3: predict intermediates for training years -------------------
  train$pred_SPC   <- predict(brt_SPC_jk, newdata = train[ , hydro_preds],
                              n.trees = brt_SPC_jk$gbm.call$best.trees)
  train$pred_logTP <- predict(brt_TP_jk,  newdata = train[ , hydro_preds],
                              n.trees = brt_TP_jk$gbm.call$best.trees)
  train$pred_logTN <- predict(brt_TN_jk,  newdata = train[ , hydro_preds],
                              n.trees = brt_TN_jk$gbm.call$best.trees)
  
  # -- Step 4: fit bloom submodel on training years -----------------------
  brt_bloom_jk <- fit_brt(train, "logCHLa", bloom_preds, tc = 4)
  
  if (is.null(brt_bloom_jk)) {
    cat(sprintf("  SKIPPING %d — bloom model fit failed\n\n", yr))
    next
  }
  
  # -- Step 5: predict bloom for withheld year ----------------------------
  test$pred_logCHLa <- predict(brt_bloom_jk,
                               newdata = test[ , bloom_preds],
                               n.trees = brt_bloom_jk$gbm.call$best.trees)
  
  # -- Step 6: evaluate ---------------------------------------------------
  r_yr   <- cor(test$logCHLa, test$pred_logCHLa, use = "complete.obs")
  rmse   <- sqrt(mean((test$logCHLa - test$pred_logCHLa)^2, na.rm = TRUE))
  
  cat(sprintf("  r = %.3f  RMSE = %.3f  n = %d\n\n",
              r_yr, rmse, nrow(test)))
  
  jk_out[[i]] <- data.frame(
    Site          = as.character(test$Site),
    Year          = test$Year,
    Month         = test$Month,
    Observed      = test$logCHLa,
    Predicted     = test$pred_logCHLa,
    pred_SPC      = test$pred_SPC,
    pred_logTP    = test$pred_logTP,
    pred_logTN    = test$pred_logTN,
    stringsAsFactors = FALSE
  )
}

# ----------------------------------------------------------------------------
# 4. Compile and summarize
# ----------------------------------------------------------------------------

jk_all <- do.call(rbind, jk_out)

r_overall    <- cor(jk_all$Observed, jk_all$Predicted, use = "complete.obs")
rmse_overall <- sqrt(mean((jk_all$Observed - jk_all$Predicted)^2,
                          na.rm = TRUE))

cat("══════════════════════════════════════════\n")
cat("Temporal Jackknife Performance Summary\n")
cat("══════════════════════════════════════════\n")
cat(sprintf("  Overall r    = %.3f\n", r_overall))
cat(sprintf("  Overall R²   = %.3f\n", r_overall^2))
cat(sprintf("  Overall RMSE = %.3f log units\n\n", rmse_overall))

# Per-year performance
cat("Per-year Performance:\n")
cat(sprintf("  %-6s  %6s  %6s  %8s  %5s\n", "Year", "r", "R²", "RMSE", "n"))
cat(paste(rep("-", 42), collapse = ""), "\n")

perf_list <- vector("list", length(years))

for (i in seq_along(years)) {
  yr <- years[i]
  d  <- jk_all[jk_all$Year == yr, ]
  if (nrow(d) < 3) next
  r    <- cor(d$Observed, d$Predicted, use = "complete.obs")
  rmse <- sqrt(mean((d$Observed - d$Predicted)^2, na.rm = TRUE))
  cat(sprintf("  %-6d  %6.3f  %6.3f  %8.3f  %5d\n",
              yr, r, r^2, rmse, nrow(d)))
  perf_list[[i]] <- data.frame(Year = yr, r = r, R2 = r^2,
                               RMSE = rmse, n = nrow(d))
}

perf_df <- do.call(rbind, perf_list)

# Per-site performance across all jackknife folds
cat("\nPer-site Performance (across all withheld years):\n")
cat(sprintf("  %-6s  %6s  %6s  %8s  %5s\n", "Site", "r", "R²", "RMSE", "n"))
cat(paste(rep("-", 42), collapse = ""), "\n")

for (s in c("DL","GR","BN","MS","BM","HU","FH")) {
  d <- jk_all[jk_all$Site == s, ]
  if (nrow(d) < 3) next
  r    <- cor(d$Observed, d$Predicted, use = "complete.obs")
  rmse <- sqrt(mean((d$Observed - d$Predicted)^2, na.rm = TRUE))
  cat(sprintf("  %-6s  %6.3f  %6.3f  %8.3f  %5d\n",
              s, r, r^2, rmse, nrow(d)))
}

# Context
cat("\n--- Context ---\n")
cat("  Training R² (V1 all observed):  0.830\n")
cat("  Training R² (V2 full chain):    0.709\n")
cat(sprintf("  LOSO R²   (spatial):            0.279\n"))
cat(sprintf("  Jackknife R² (temporal):        %.3f\n", r_overall^2))

# ----------------------------------------------------------------------------
# 5. Plots
# ----------------------------------------------------------------------------

cat("\nGenerating jackknife plots...\n")

site_cols <- c(DL = "#E41A1C", GR = "#377EB8", BN = "#4DAF4A",
               MS = "#984EA3", BM = "#FF7F00", HU = "#A65628",
               FH = "#F781BF")

pdf("4_products/diagnostics/sem_jackknife_fit.pdf",
    width = 12, height = 10)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# 1. Overall observed vs predicted
rng <- range(c(jk_all$Observed, jk_all$Predicted), na.rm = TRUE)
plot(jk_all$Observed, jk_all$Predicted,
     xlim = rng, ylim = rng,
     xlab = "Observed log10(CHLa)",
     ylab = "Jackknife Predicted log10(CHLa)",
     main = sprintf("All years  |  r=%.3f  R²=%.3f",
                    r_overall, r_overall^2),
     pch  = 16,
     col  = site_cols[jk_all$Site],
     cex  = 0.8)
abline(0, 1, col = "grey40", lty = 2)
legend("topleft", legend = names(site_cols), col = site_cols,
       pch = 16, cex = 0.65, bty = "n")

# 2. Per-year r values over time
plot(perf_df$Year, perf_df$r,
     type = "b", pch = 16,
     xlab = "Year withheld",
     ylab = "Pearson r",
     main = "Jackknife r by Year",
     ylim = c(-0.2, 1),
     col  = "steelblue")
abline(h = 0,          col = "grey40", lty = 2)
abline(h = r_overall,  col = "red",    lty = 2, lwd = 1.5)
text(max(perf_df$Year), r_overall + 0.04,
     sprintf("Overall r=%.3f", r_overall),
     col = "red", cex = 0.75, adj = 1)

# 3. Per-year RMSE over time
plot(perf_df$Year, perf_df$RMSE,
     type = "b", pch = 16,
     xlab = "Year withheld",
     ylab = "RMSE (log10 units)",
     main = "Jackknife RMSE by Year",
     col  = "darkorange")
abline(h = rmse_overall, col = "red", lty = 2, lwd = 1.5)
text(max(perf_df$Year), rmse_overall + 0.01,
     sprintf("Overall RMSE=%.3f", rmse_overall),
     col = "red", cex = 0.75, adj = 1)

# 4. Residuals over time
jk_all$resid <- jk_all$Observed - jk_all$Predicted
boxplot(resid ~ Year, data = jk_all,
        xlab = "Year",
        ylab = "Residual (obs - pred)",
        main = "Residuals by Year",
        col  = "lightblue",
        border = "steelblue",
        las = 2,
        cex.axis = 0.7)
abline(h = 0, col = "red", lty = 2)

dev.off()
cat("Plots saved to 4_products/diagnostics/sem_jackknife_fit.pdf\n")

# ----------------------------------------------------------------------------
# 6. Save outputs
# ----------------------------------------------------------------------------

write_csv(jk_all,   "2_incremental/sem_jackknife_predictions.csv")
write_csv(perf_df,  "2_incremental/sem_jackknife_performance.csv")

cat("\nSaved:\n")
cat("  2_incremental/sem_jackknife_predictions.csv\n")
cat("  2_incremental/sem_jackknife_performance.csv\n")
cat("Done.\n")