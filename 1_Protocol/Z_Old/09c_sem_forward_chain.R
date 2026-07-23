# ============================================================================
# 09_2_sem_forward_chain.R
# UCFR Filamentous Algae Project
# Stage 8.2: Forward chaining temporal validation of full BRT-SEM chain
#
# Input:    2_incremental/ucfr_model_ready.csv
# Outputs:  2_incremental/sem_fwdchain_predictions.csv
#           2_incremental/sem_fwdchain_performance.csv
#           4_products/diagnostics/sem_fwdchain_fit.pdf
#
# Procedure:
#   For each year Y from (min_year + min_train_window) to max_year:
#     1. Train all four submodels on years < Y only
#     2. Predict year Y through full chain
#     3. Compare to observed logCHLa
#
# Notes:
#   - Respects temporal ordering — never uses future data to predict past
#   - Minimum training window = 10 years before first test year
#   - First test year = min_year + min_train_window
#   - Most realistic analog to climate scenario projection
#   - BRT settings: tc=3 for SPC/TP/TN, tc=4 for bloom
# ============================================================================

library(readr)
library(dismo)
library(gbm)

# ----------------------------------------------------------------------------
# 0. Configuration
# ----------------------------------------------------------------------------

MIN_TRAIN_YEARS <- 3  # minimum years of data before first prediction

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

chain_vars <- c("Site", "Year", "Month",
                hydro_preds, "SPC", "logTP_obs", "logTN_obs", "logCHLa")

dat_model <- as.data.frame(dat[complete.cases(dat[ , chain_vars]), ])

all_years  <- sort(unique(dat_model$Year))
min_year   <- min(all_years)
first_test <- min_year + MIN_TRAIN_YEARS
test_years <- all_years[all_years >= first_test]

cat(sprintf("Complete cases:      %d rows\n", nrow(dat_model)))
cat(sprintf("Full year range:     %d to %d\n", min_year, max(all_years)))
cat(sprintf("Min training window: %d years\n", MIN_TRAIN_YEARS))
cat(sprintf("First test year:     %d\n", first_test))
cat(sprintf("Test years:          %d to %d (%d iterations)\n\n",
            min(test_years), max(test_years), length(test_years)))

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
# 3. Forward chaining loop
# ----------------------------------------------------------------------------

cat("Running forward chaining validation...\n\n")

fc_out <- vector("list", length(test_years))

for (i in seq_along(test_years)) {
  yr    <- test_years[i]
  train <- as.data.frame(dat_model[dat_model$Year < yr, ])
  test  <- as.data.frame(dat_model[dat_model$Year == yr, ])
  
  cat(sprintf("--- Predicting %d (trained on %d-%d, n_train=%d, n_test=%d) ---\n",
              yr, min(train$Year), max(train$Year), nrow(train), nrow(test)))
  
  # -- Step 1: fit intermediate submodels on all prior years ---------------
  brt_SPC_fc <- fit_brt(train, "SPC",       hydro_preds, tc = 3)
  brt_TP_fc  <- fit_brt(train, "logTP_obs", hydro_preds, tc = 3)
  brt_TN_fc  <- fit_brt(train, "logTN_obs", hydro_preds, tc = 3)
  
  if (any(sapply(list(brt_SPC_fc, brt_TP_fc, brt_TN_fc), is.null))) {
    cat(sprintf("  SKIPPING %d — intermediate model fit failed\n\n", yr))
    next
  }
  
  # -- Step 2: predict intermediates for test year -------------------------
  test$pred_SPC   <- predict(brt_SPC_fc, newdata = test[ , hydro_preds],
                             n.trees = brt_SPC_fc$gbm.call$best.trees)
  test$pred_logTP <- predict(brt_TP_fc,  newdata = test[ , hydro_preds],
                             n.trees = brt_TP_fc$gbm.call$best.trees)
  test$pred_logTN <- predict(brt_TN_fc,  newdata = test[ , hydro_preds],
                             n.trees = brt_TN_fc$gbm.call$best.trees)
  
  # -- Step 3: predict intermediates for training years -------------------
  train$pred_SPC   <- predict(brt_SPC_fc, newdata = train[ , hydro_preds],
                              n.trees = brt_SPC_fc$gbm.call$best.trees)
  train$pred_logTP <- predict(brt_TP_fc,  newdata = train[ , hydro_preds],
                              n.trees = brt_TP_fc$gbm.call$best.trees)
  train$pred_logTN <- predict(brt_TN_fc,  newdata = train[ , hydro_preds],
                              n.trees = brt_TN_fc$gbm.call$best.trees)
  
  # -- Step 4: fit bloom submodel on all prior years ----------------------
  brt_bloom_fc <- fit_brt(train, "logCHLa", bloom_preds, tc = 4)
  
  if (is.null(brt_bloom_fc)) {
    cat(sprintf("  SKIPPING %d — bloom model fit failed\n\n", yr))
    next
  }
  
  # -- Step 5: predict bloom for test year --------------------------------
  test$pred_logCHLa <- predict(brt_bloom_fc,
                               newdata = test[ , bloom_preds],
                               n.trees = brt_bloom_fc$gbm.call$best.trees)
  
  # -- Step 6: evaluate ---------------------------------------------------
  r_yr <- cor(test$logCHLa, test$pred_logCHLa, use = "complete.obs")
  rmse <- sqrt(mean((test$logCHLa - test$pred_logCHLa)^2, na.rm = TRUE))
  
  cat(sprintf("  r = %.3f  RMSE = %.3f  n_train_years = %d\n\n",
              r_yr, rmse, length(unique(train$Year))))
  
  fc_out[[i]] <- data.frame(
    Site          = as.character(test$Site),
    Year          = test$Year,
    Month         = test$Month,
    Observed      = test$logCHLa,
    Predicted     = test$pred_logCHLa,
    pred_SPC      = test$pred_SPC,
    pred_logTP    = test$pred_logTP,
    pred_logTN    = test$pred_logTN,
    n_train_years = length(unique(train$Year)),
    stringsAsFactors = FALSE
  )
}

# ----------------------------------------------------------------------------
# 4. Compile and summarize
# ----------------------------------------------------------------------------

fc_all <- do.call(rbind, fc_out)

r_overall    <- cor(fc_all$Observed, fc_all$Predicted, use = "complete.obs")
rmse_overall <- sqrt(mean((fc_all$Observed - fc_all$Predicted)^2,
                          na.rm = TRUE))

cat("══════════════════════════════════════════\n")
cat("Forward Chaining Performance Summary\n")
cat("══════════════════════════════════════════\n")
cat(sprintf("  Overall r    = %.3f\n", r_overall))
cat(sprintf("  Overall R²   = %.3f\n", r_overall^2))
cat(sprintf("  Overall RMSE = %.3f log units\n\n", rmse_overall))

# Per-year performance
cat("Per-year Performance:\n")
cat(sprintf("  %-6s  %6s  %6s  %8s  %5s  %10s\n",
            "Year", "r", "R²", "RMSE", "n", "train_yrs"))
cat(paste(rep("-", 52), collapse = ""), "\n")

perf_list <- vector("list", length(test_years))

for (i in seq_along(test_years)) {
  yr <- test_years[i]
  d  <- fc_all[fc_all$Year == yr, ]
  if (nrow(d) < 3) next
  r      <- cor(d$Observed, d$Predicted, use = "complete.obs")
  rmse   <- sqrt(mean((d$Observed - d$Predicted)^2, na.rm = TRUE))
  n_tr   <- unique(d$n_train_years)
  cat(sprintf("  %-6d  %6.3f  %6.3f  %8.3f  %5d  %10d\n",
              yr, r, r^2, rmse, nrow(d), n_tr))
  perf_list[[i]] <- data.frame(Year = yr, r = r, R2 = r^2,
                               RMSE = rmse, n = nrow(d),
                               n_train_years = n_tr)
}

perf_df <- do.call(rbind, perf_list)

# Per-site performance
cat("\nPer-site Performance (across all test years):\n")
cat(sprintf("  %-6s  %6s  %6s  %8s  %5s\n", "Site", "r", "R²", "RMSE", "n"))
cat(paste(rep("-", 42), collapse = ""), "\n")

for (s in c("DL","GR","BN","MS","BM","HU","FH")) {
  d <- fc_all[fc_all$Site == s, ]
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
cat(sprintf("  Forward chain R² (temporal):    %.3f\n", r_overall^2))

# ----------------------------------------------------------------------------
# 5. Plots
# ----------------------------------------------------------------------------

cat("\nGenerating forward chaining plots...\n")

site_cols <- c(DL = "#E41A1C", GR = "#377EB8", BN = "#4DAF4A",
               MS = "#984EA3", BM = "#FF7F00", HU = "#A65628",
               FH = "#F781BF")

pdf("4_products/diagnostics/sem_fwdchain_fit.pdf",
    width = 12, height = 10)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# 1. Overall observed vs predicted
rng <- range(c(fc_all$Observed, fc_all$Predicted), na.rm = TRUE)
plot(fc_all$Observed, fc_all$Predicted,
     xlim = rng, ylim = rng,
     xlab = "Observed log10(CHLa)",
     ylab = "Forward Chain Predicted log10(CHLa)",
     main = sprintf("All test years  |  r=%.3f  R²=%.3f",
                    r_overall, r_overall^2),
     pch  = 16,
     col  = site_cols[fc_all$Site],
     cex  = 0.8)
abline(0, 1, col = "grey40", lty = 2)
legend("topleft", legend = names(site_cols), col = site_cols,
       pch = 16, cex = 0.65, bty = "n")

# 2. Per-year r over time with training window size
par(mar = c(4, 4, 3, 4))
plot(perf_df$Year, perf_df$r,
     type = "b", pch = 16,
     xlab = "Test Year",
     ylab = "Pearson r",
     main = "Forward Chain r by Year",
     ylim = c(-0.5, 1),
     col  = "steelblue")
abline(h = 0,         col = "grey40", lty = 2)
abline(h = r_overall, col = "red",    lty = 2, lwd = 1.5)
# Add training window size on secondary axis
par(new = TRUE)
plot(perf_df$Year, perf_df$n_train_years,
     type = "l", lty = 3, col = "grey60",
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(4, col = "grey60", col.axis = "grey60")
mtext("Training years", side = 4, line = 2.5, col = "grey60", cex = 0.8)
par(mar = c(4, 4, 3, 1))

# 3. Per-year RMSE over time
plot(perf_df$Year, perf_df$RMSE,
     type = "b", pch = 16,
     xlab = "Test Year",
     ylab = "RMSE (log10 units)",
     main = "Forward Chain RMSE by Year",
     col  = "darkorange")
abline(h = rmse_overall, col = "red", lty = 2, lwd = 1.5)

# 4. Residuals by year
fc_all$resid <- fc_all$Observed - fc_all$Predicted
boxplot(resid ~ Year, data = fc_all,
        xlab = "Year",
        ylab = "Residual (obs - pred)",
        main = "Residuals by Test Year",
        col  = "lightblue",
        border = "steelblue",
        las = 2,
        cex.axis = 0.65)
abline(h = 0, col = "red", lty = 2)

dev.off()
cat("Plots saved to 4_products/diagnostics/sem_fwdchain_fit.pdf\n")

# ----------------------------------------------------------------------------
# 6. Save outputs
# ----------------------------------------------------------------------------

write_csv(fc_all,  "2_incremental/sem_fwdchain_predictions.csv")
write_csv(perf_df, "2_incremental/sem_fwdchain_performance.csv")

cat("\nSaved:\n")
cat("  2_incremental/sem_fwdchain_predictions.csv\n")
cat("  2_incremental/sem_fwdchain_performance.csv\n")
cat("Done.\n")