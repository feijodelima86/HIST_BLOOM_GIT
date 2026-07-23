# ============================================================================
# 12_sensitivity_no_BM.R
# UCFR Filamentous Algae Project
# Sensitivity analysis: V1, V2, LOSO, Jackknife excluding BM site
# ============================================================================

library(readr)
library(dismo)
library(gbm)

cat("Reading data...\n")
dat <- as.data.frame(read_csv("2_incremental/ucfr_model_ready.csv",
                              show_col_types = FALSE))
nut <- as.data.frame(read_csv("2_incremental/brt_nutrients_fitted.csv",
                              show_col_types = FALSE))

dat <- merge(dat, nut[ , c("Site", "Year", "Month",
                           "pred_SPC", "pred_logTP", "pred_logTN")],
             by = c("Site", "Year", "Month"), all.x = TRUE)

dat$logCHLa   <- log10(dat$CHLa)
dat$logTP_obs <- log10(dat$TP_mg_L)
dat$logTN_obs <- log10(dat$TN_mg_L)

# Drop BM
dat <- dat[dat$Site != "BM", ]
cat(sprintf("Rows after dropping BM: %d\n\n", nrow(dat)))

hydro_preds <- c("anomaly", "Q_obs_cfs", "Temp_oC", "Days_Since_Freshet")
preds_v1    <- c("SPC", "logTP_obs", "logTN_obs", hydro_preds)
preds_v2    <- c("pred_SPC", "pred_logTP", "pred_logTN", hydro_preds)
bloom_preds <- c("pred_SPC", "pred_logTP", "pred_logTN", hydro_preds)

vars_v1 <- c("Site", "Year", "Month", preds_v1, "logCHLa")
vars_v2 <- c("Site", "Year", "Month", preds_v2, "logCHLa")
chain_vars <- c("Site", "Year", "Month", hydro_preds,
                "SPC", "logTP_obs", "logTN_obs", "logCHLa")

dat_v1 <- as.data.frame(dat[complete.cases(dat[ , vars_v1]), ])
dat_v2 <- as.data.frame(dat[complete.cases(dat[ , vars_v2]), ])
dat_chain <- as.data.frame(dat[complete.cases(dat[ , chain_vars]), ])

cat(sprintf("Complete cases V1: %d  V2: %d  Chain: %d\n\n",
            nrow(dat_v1), nrow(dat_v2), nrow(dat_chain)))

# ----------------------------------------------------------------------------
# BRT fitting function
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
    error = function(e) { cat(sprintf("FIT FAILED: %s\n", conditionMessage(e))); NULL }
  )
}

# ----------------------------------------------------------------------------
# V1 — all observed
# ----------------------------------------------------------------------------

cat("Fitting V1 (all observed, no BM)...\n")
brt_v1 <- fit_brt(dat_v1, "logCHLa", preds_v1, tc = 4)
pred_v1 <- predict(brt_v1, newdata = dat_v1,
                   n.trees = brt_v1$gbm.call$best.trees)
r_v1 <- cor(dat_v1$logCHLa, pred_v1)
cat(sprintf("  V1: r=%.3f  R²=%.3f\n\n", r_v1, r_v1^2))

# ----------------------------------------------------------------------------
# V2 — full chain
# ----------------------------------------------------------------------------

cat("Fitting V2 (full chain, no BM)...\n")
brt_v2 <- fit_brt(dat_v2, "logCHLa", preds_v2, tc = 4)
pred_v2 <- predict(brt_v2, newdata = dat_v2,
                   n.trees = brt_v2$gbm.call$best.trees)
r_v2 <- cor(dat_v2$logCHLa, pred_v2)
cat(sprintf("  V2: r=%.3f  R²=%.3f\n\n", r_v2, r_v2^2))

# ----------------------------------------------------------------------------
# LOSO
# ----------------------------------------------------------------------------

cat("Running LOSO (no BM)...\n")
sites    <- sort(unique(dat_chain$Site))
loso_out <- vector("list", length(sites))

for (i in seq_along(sites)) {
  s     <- sites[i]
  train <- as.data.frame(dat_chain[dat_chain$Site != s, ])
  test  <- as.data.frame(dat_chain[dat_chain$Site == s, ])
  cat(sprintf("  Leaving out %s...\n", s))
  
  brt_SPC_l  <- fit_brt(train, "SPC",       hydro_preds, tc = 3)
  brt_TP_l   <- fit_brt(train, "logTP_obs", hydro_preds, tc = 3)
  brt_TN_l   <- fit_brt(train, "logTN_obs", hydro_preds, tc = 3)
  if (any(sapply(list(brt_SPC_l, brt_TP_l, brt_TN_l), is.null))) next
  
  test$pred_SPC   <- predict(brt_SPC_l, newdata = test[ , hydro_preds],
                             n.trees = brt_SPC_l$gbm.call$best.trees)
  test$pred_logTP <- predict(brt_TP_l,  newdata = test[ , hydro_preds],
                             n.trees = brt_TP_l$gbm.call$best.trees)
  test$pred_logTN <- predict(brt_TN_l,  newdata = test[ , hydro_preds],
                             n.trees = brt_TN_l$gbm.call$best.trees)
  
  train$pred_SPC   <- predict(brt_SPC_l, newdata = train[ , hydro_preds],
                              n.trees = brt_SPC_l$gbm.call$best.trees)
  train$pred_logTP <- predict(brt_TP_l,  newdata = train[ , hydro_preds],
                              n.trees = brt_TP_l$gbm.call$best.trees)
  train$pred_logTN <- predict(brt_TN_l,  newdata = train[ , hydro_preds],
                              n.trees = brt_TN_l$gbm.call$best.trees)
  
  brt_bloom_l <- fit_brt(train, "logCHLa", bloom_preds, tc = 4)
  if (is.null(brt_bloom_l)) next
  
  test$pred_logCHLa <- predict(brt_bloom_l,
                               newdata = test[ , bloom_preds],
                               n.trees = brt_bloom_l$gbm.call$best.trees)
  
  r_s <- cor(test$logCHLa, test$pred_logCHLa, use = "complete.obs")
  cat(sprintf("    r = %.3f\n", r_s))
  
  loso_out[[i]] <- data.frame(
    Site = test$Site, Year = test$Year,
    Observed = test$logCHLa, Predicted = test$pred_logCHLa
  )
}

loso_all <- do.call(rbind, loso_out)
r_loso   <- cor(loso_all$Observed, loso_all$Predicted, use = "complete.obs")
cat(sprintf("\n  LOSO: r=%.3f  R²=%.3f\n\n", r_loso, r_loso^2))

# ----------------------------------------------------------------------------
# Jackknife
# ----------------------------------------------------------------------------

cat("Running Jackknife (no BM)...\n")
years   <- sort(unique(dat_chain$Year))
jk_out  <- vector("list", length(years))

for (i in seq_along(years)) {
  yr    <- years[i]
  train <- as.data.frame(dat_chain[dat_chain$Year != yr, ])
  test  <- as.data.frame(dat_chain[dat_chain$Year == yr, ])
  
  brt_SPC_j  <- fit_brt(train, "SPC",       hydro_preds, tc = 3)
  brt_TP_j   <- fit_brt(train, "logTP_obs", hydro_preds, tc = 3)
  brt_TN_j   <- fit_brt(train, "logTN_obs", hydro_preds, tc = 3)
  if (any(sapply(list(brt_SPC_j, brt_TP_j, brt_TN_j), is.null))) next
  
  test$pred_SPC   <- predict(brt_SPC_j, newdata = test[ , hydro_preds],
                             n.trees = brt_SPC_j$gbm.call$best.trees)
  test$pred_logTP <- predict(brt_TP_j,  newdata = test[ , hydro_preds],
                             n.trees = brt_TP_j$gbm.call$best.trees)
  test$pred_logTN <- predict(brt_TN_j,  newdata = test[ , hydro_preds],
                             n.trees = brt_TN_j$gbm.call$best.trees)
  
  train$pred_SPC   <- predict(brt_SPC_j, newdata = train[ , hydro_preds],
                              n.trees = brt_SPC_j$gbm.call$best.trees)
  train$pred_logTP <- predict(brt_TP_j,  newdata = train[ , hydro_preds],
                              n.trees = brt_TP_j$gbm.call$best.trees)
  train$pred_logTN <- predict(brt_TN_j,  newdata = train[ , hydro_preds],
                              n.trees = brt_TN_j$gbm.call$best.trees)
  
  brt_bloom_j <- fit_brt(train, "logCHLa", bloom_preds, tc = 4)
  if (is.null(brt_bloom_j)) next
  
  test$pred_logCHLa <- predict(brt_bloom_j,
                               newdata = test[ , bloom_preds],
                               n.trees = brt_bloom_j$gbm.call$best.trees)
  
  jk_out[[i]] <- data.frame(
    Site = test$Site, Year = test$Year,
    Observed = test$logCHLa, Predicted = test$pred_logCHLa
  )
}

jk_all  <- do.call(rbind, jk_out)
r_jk    <- cor(jk_all$Observed, jk_all$Predicted, use = "complete.obs")
cat(sprintf("\n  Jackknife: r=%.3f  R²=%.3f\n\n", r_jk, r_jk^2))

# ----------------------------------------------------------------------------
# Summary comparison
# ----------------------------------------------------------------------------

cat("══════════════════════════════════════════════════════\n")
cat("Performance Summary — 6 sites (BM excluded)\n")
cat("══════════════════════════════════════════════════════\n")
cat(sprintf("  %-35s  %6s  %6s\n", "Model", "r", "R²"))
cat(paste(rep("-", 52), collapse = ""), "\n")
cat(sprintf("  %-35s  %6.3f  %6.3f\n", "V1 training (all observed)",   r_v1,  r_v1^2))
cat(sprintf("  %-35s  %6.3f  %6.3f\n", "V2 training (full chain)",     r_v2,  r_v2^2))
cat(sprintf("  %-35s  %6.3f  %6.3f\n", "LOSO (spatial)",               r_loso, r_loso^2))
cat(sprintf("  %-35s  %6.3f  %6.3f\n", "Jackknife (temporal)",         r_jk,  r_jk^2))
cat("\nAll sites comparison:\n")
cat("  V1:        R²=0.830\n")
cat("  V2:        R²=0.709\n")
cat("  LOSO:      R²=0.279\n")
cat("  Jackknife: R²=0.404\n")
cat("Done.\n")

# ============================================================================
# 12_sensitivity_no_FH.R
# UCFR Filamentous Algae Project
# Sensitivity analysis: V1, V2, LOSO, Jackknife excluding FH site
# ============================================================================

library(readr)
library(dismo)
library(gbm)

cat("Reading data...\n")
dat <- as.data.frame(read_csv("2_incremental/ucfr_model_ready.csv",
                              show_col_types = FALSE))
nut <- as.data.frame(read_csv("2_incremental/brt_nutrients_fitted.csv",
                              show_col_types = FALSE))

dat <- merge(dat, nut[ , c("Site", "Year", "Month",
                           "pred_SPC", "pred_logTP", "pred_logTN")],
             by = c("Site", "Year", "Month"), all.x = TRUE)

dat$logCHLa   <- log10(dat$CHLa)
dat$logTP_obs <- log10(dat$TP_mg_L)
dat$logTN_obs <- log10(dat$TN_mg_L)

# Drop FH
dat <- dat[dat$Site != "FH", ]
cat(sprintf("Rows after dropping FH: %d\n\n", nrow(dat)))

hydro_preds <- c("anomaly", "Q_obs_cfs", "Temp_oC", "Days_Since_Freshet")
preds_v1    <- c("SPC", "logTP_obs", "logTN_obs", hydro_preds)
preds_v2    <- c("pred_SPC", "pred_logTP", "pred_logTN", hydro_preds)
bloom_preds <- c("pred_SPC", "pred_logTP", "pred_logTN", hydro_preds)

vars_v1 <- c("Site", "Year", "Month", preds_v1, "logCHLa")
vars_v2 <- c("Site", "Year", "Month", preds_v2, "logCHLa")
chain_vars <- c("Site", "Year", "Month", hydro_preds,
                "SPC", "logTP_obs", "logTN_obs", "logCHLa")

dat_v1 <- as.data.frame(dat[complete.cases(dat[ , vars_v1]), ])
dat_v2 <- as.data.frame(dat[complete.cases(dat[ , vars_v2]), ])
dat_chain <- as.data.frame(dat[complete.cases(dat[ , chain_vars]), ])

cat(sprintf("Complete cases V1: %d  V2: %d  Chain: %d\n\n",
            nrow(dat_v1), nrow(dat_v2), nrow(dat_chain)))

# ----------------------------------------------------------------------------
# BRT fitting function
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
    error = function(e) { cat(sprintf("FIT FAILED: %s\n", conditionMessage(e))); NULL }
  )
}

# ----------------------------------------------------------------------------
# V1 — all observed
# ----------------------------------------------------------------------------

cat("Fitting V1 (all observed, no FH)...\n")
brt_v1 <- fit_brt(dat_v1, "logCHLa", preds_v1, tc = 4)
pred_v1 <- predict(brt_v1, newdata = dat_v1,
                   n.trees = brt_v1$gbm.call$best.trees)
r_v1 <- cor(dat_v1$logCHLa, pred_v1)
cat(sprintf("  V1: r=%.3f  R²=%.3f\n\n", r_v1, r_v1^2))

# ----------------------------------------------------------------------------
# V2 — full chain
# ----------------------------------------------------------------------------

cat("Fitting V2 (full chain, no FH)...\n")
brt_v2 <- fit_brt(dat_v2, "logCHLa", preds_v2, tc = 4)
pred_v2 <- predict(brt_v2, newdata = dat_v2,
                   n.trees = brt_v2$gbm.call$best.trees)
r_v2 <- cor(dat_v2$logCHLa, pred_v2)
cat(sprintf("  V2: r=%.3f  R²=%.3f\n\n", r_v2, r_v2^2))

# ----------------------------------------------------------------------------
# LOSO
# ----------------------------------------------------------------------------

cat("Running LOSO (no FH)...\n")
sites    <- sort(unique(dat_chain$Site))
loso_out <- vector("list", length(sites))

for (i in seq_along(sites)) {
  s     <- sites[i]
  train <- as.data.frame(dat_chain[dat_chain$Site != s, ])
  test  <- as.data.frame(dat_chain[dat_chain$Site == s, ])
  cat(sprintf("  Leaving out %s...\n", s))
  
  brt_SPC_l  <- fit_brt(train, "SPC",       hydro_preds, tc = 3)
  brt_TP_l   <- fit_brt(train, "logTP_obs", hydro_preds, tc = 3)
  brt_TN_l   <- fit_brt(train, "logTN_obs", hydro_preds, tc = 3)
  if (any(sapply(list(brt_SPC_l, brt_TP_l, brt_TN_l), is.null))) next
  
  test$pred_SPC   <- predict(brt_SPC_l, newdata = test[ , hydro_preds],
                             n.trees = brt_SPC_l$gbm.call$best.trees)
  test$pred_logTP <- predict(brt_TP_l,  newdata = test[ , hydro_preds],
                             n.trees = brt_TP_l$gbm.call$best.trees)
  test$pred_logTN <- predict(brt_TN_l,  newdata = test[ , hydro_preds],
                             n.trees = brt_TN_l$gbm.call$best.trees)
  
  train$pred_SPC   <- predict(brt_SPC_l, newdata = train[ , hydro_preds],
                              n.trees = brt_SPC_l$gbm.call$best.trees)
  train$pred_logTP <- predict(brt_TP_l,  newdata = train[ , hydro_preds],
                              n.trees = brt_TP_l$gbm.call$best.trees)
  train$pred_logTN <- predict(brt_TN_l,  newdata = train[ , hydro_preds],
                              n.trees = brt_TN_l$gbm.call$best.trees)
  
  brt_bloom_l <- fit_brt(train, "logCHLa", bloom_preds, tc = 4)
  if (is.null(brt_bloom_l)) next
  
  test$pred_logCHLa <- predict(brt_bloom_l,
                               newdata = test[ , bloom_preds],
                               n.trees = brt_bloom_l$gbm.call$best.trees)
  
  r_s <- cor(test$logCHLa, test$pred_logCHLa, use = "complete.obs")
  cat(sprintf("    r = %.3f\n", r_s))
  
  loso_out[[i]] <- data.frame(
    Site = test$Site, Year = test$Year,
    Observed = test$logCHLa, Predicted = test$pred_logCHLa
  )
}

loso_all <- do.call(rbind, loso_out)
r_loso   <- cor(loso_all$Observed, loso_all$Predicted, use = "complete.obs")
cat(sprintf("\n  LOSO: r=%.3f  R²=%.3f\n\n", r_loso, r_loso^2))

# ----------------------------------------------------------------------------
# Jackknife
# ----------------------------------------------------------------------------

cat("Running Jackknife (no FH)...\n")
years   <- sort(unique(dat_chain$Year))
jk_out  <- vector("list", length(years))

for (i in seq_along(years)) {
  yr    <- years[i]
  train <- as.data.frame(dat_chain[dat_chain$Year != yr, ])
  test  <- as.data.frame(dat_chain[dat_chain$Year == yr, ])
  
  brt_SPC_j  <- fit_brt(train, "SPC",       hydro_preds, tc = 3)
  brt_TP_j   <- fit_brt(train, "logTP_obs", hydro_preds, tc = 3)
  brt_TN_j   <- fit_brt(train, "logTN_obs", hydro_preds, tc = 3)
  if (any(sapply(list(brt_SPC_j, brt_TP_j, brt_TN_j), is.null))) next
  
  test$pred_SPC   <- predict(brt_SPC_j, newdata = test[ , hydro_preds],
                             n.trees = brt_SPC_j$gbm.call$best.trees)
  test$pred_logTP <- predict(brt_TP_j,  newdata = test[ , hydro_preds],
                             n.trees = brt_TP_j$gbm.call$best.trees)
  test$pred_logTN <- predict(brt_TN_j,  newdata = test[ , hydro_preds],
                             n.trees = brt_TN_j$gbm.call$best.trees)
  
  train$pred_SPC   <- predict(brt_SPC_j, newdata = train[ , hydro_preds],
                              n.trees = brt_SPC_j$gbm.call$best.trees)
  train$pred_logTP <- predict(brt_TP_j,  newdata = train[ , hydro_preds],
                              n.trees = brt_TP_j$gbm.call$best.trees)
  train$pred_logTN <- predict(brt_TN_j,  newdata = train[ , hydro_preds],
                              n.trees = brt_TN_j$gbm.call$best.trees)
  
  brt_bloom_j <- fit_brt(train, "logCHLa", bloom_preds, tc = 4)
  if (is.null(brt_bloom_j)) next
  
  test$pred_logCHLa <- predict(brt_bloom_j,
                               newdata = test[ , bloom_preds],
                               n.trees = brt_bloom_j$gbm.call$best.trees)
  
  jk_out[[i]] <- data.frame(
    Site = test$Site, Year = test$Year,
    Observed = test$logCHLa, Predicted = test$pred_logCHLa
  )
}

jk_all  <- do.call(rbind, jk_out)
r_jk    <- cor(jk_all$Observed, jk_all$Predicted, use = "complete.obs")
cat(sprintf("\n  Jackknife: r=%.3f  R²=%.3f\n\n", r_jk, r_jk^2))

# ----------------------------------------------------------------------------
# Summary comparison
# ----------------------------------------------------------------------------

cat("══════════════════════════════════════════════════════\n")
cat("Performance Summary — 6 sites (FH excluded)\n")
cat("══════════════════════════════════════════════════════\n")
cat(sprintf("  %-35s  %6s  %6s\n", "Model", "r", "R²"))
cat(paste(rep("-", 52), collapse = ""), "\n")
cat(sprintf("  %-35s  %6.3f  %6.3f\n", "V1 training (all observed)",   r_v1,  r_v1^2))
cat(sprintf("  %-35s  %6.3f  %6.3f\n", "V2 training (full chain)",     r_v2,  r_v2^2))
cat(sprintf("  %-35s  %6.3f  %6.3f\n", "LOSO (spatial)",               r_loso, r_loso^2))
cat(sprintf("  %-35s  %6.3f  %6.3f\n", "Jackknife (temporal)",         r_jk,  r_jk^2))
cat("\nAll sites comparison:\n")
cat("  V1:        R²=0.830\n")
cat("  V2:        R²=0.709\n")
cat("  LOSO:      R²=0.279\n")
cat("  Jackknife: R²=0.404\n")
cat("Done.\n")


# ============================================================================
# 14_sensitivity_no_outliers.R
# UCFR Filamentous Algae Project
# Sensitivity analysis: V1, V2, LOSO, Jackknife excluding the 11
# observations flagged as outliers in both full model and jackknife
#
# Flagged observations (Site Year Month):
#   FH 2008 8, BM 2009 9, FH 2016 9, FH 2022 8, DL 2009 9,
#   BM 2016 7, BM 2013 8, BM 2015 7, FH 2021 8, MS 2023 7, BM 2019 9
# ============================================================================

library(readr)
library(dismo)
library(gbm)

# ----------------------------------------------------------------------------
# 1. Define flagged observations and load data
# ----------------------------------------------------------------------------

flagged <- data.frame(
  Site  = c("FH","BM","FH","FH","DL","BM","BM","BM","FH","MS","BM"),
  Year  = c(2008,2009,2016,2022,2009,2016,2013,2015,2021,2023,2019),
  Month = c(   8,   9,   9,   8,   9,   7,   8,   7,   8,   7,   9),
  stringsAsFactors = FALSE
)

cat("Flagged observations to exclude:\n")
for (i in seq_len(nrow(flagged))) {
  cat(sprintf("  %s %d %s\n", flagged$Site[i], flagged$Year[i],
              month.abb[flagged$Month[i]]))
}
cat("\n")

cat("Reading data...\n")
dat <- as.data.frame(read_csv("2_incremental/ucfr_model_ready.csv",
                              show_col_types = FALSE))
nut <- as.data.frame(read_csv("2_incremental/brt_nutrients_fitted.csv",
                              show_col_types = FALSE))

dat <- merge(dat, nut[ , c("Site", "Year", "Month",
                           "pred_SPC", "pred_logTP", "pred_logTN")],
             by = c("Site", "Year", "Month"), all.x = TRUE)

dat$logCHLa   <- log10(dat$CHLa)
dat$logTP_obs <- log10(dat$TP_mg_L)
dat$logTN_obs <- log10(dat$TN_mg_L)

# Flag and remove outlier observations
dat$flag_key     <- paste(dat$Site, dat$Year, dat$Month)
flagged$flag_key <- paste(flagged$Site, flagged$Year, flagged$Month)
dat$outlier_flag <- dat$flag_key %in% flagged$flag_key

cat(sprintf("Total rows before exclusion: %d\n", nrow(dat)))
cat(sprintf("Flagged rows excluded:       %d\n", sum(dat$outlier_flag)))
dat <- dat[!dat$outlier_flag, ]
cat(sprintf("Rows after exclusion:        %d\n\n", nrow(dat)))

# ----------------------------------------------------------------------------
# 2. Setup
# ----------------------------------------------------------------------------

hydro_preds <- c("anomaly", "Q_obs_cfs", "Temp_oC", "Days_Since_Freshet")
preds_v1    <- c("SPC", "logTP_obs", "logTN_obs", hydro_preds)
preds_v2    <- c("pred_SPC", "pred_logTP", "pred_logTN", hydro_preds)
bloom_preds <- c("pred_SPC", "pred_logTP", "pred_logTN", hydro_preds)

vars_v1    <- c("Site", "Year", "Month", preds_v1, "logCHLa")
vars_v2    <- c("Site", "Year", "Month", preds_v2, "logCHLa")
chain_vars <- c("Site", "Year", "Month", hydro_preds,
                "SPC", "logTP_obs", "logTN_obs", "logCHLa")

dat_v1    <- as.data.frame(dat[complete.cases(dat[ , vars_v1]), ])
dat_v2    <- as.data.frame(dat[complete.cases(dat[ , vars_v2]), ])
dat_chain <- as.data.frame(dat[complete.cases(dat[ , chain_vars]), ])

cat(sprintf("Complete cases V1: %d  V2: %d  Chain: %d\n\n",
            nrow(dat_v1), nrow(dat_v2), nrow(dat_chain)))

# ----------------------------------------------------------------------------
# 3. BRT fitting function
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
    error = function(e) { cat(sprintf("FIT FAILED: %s\n", conditionMessage(e))); NULL }
  )
}

# ----------------------------------------------------------------------------
# 4. V1 — all observed
# ----------------------------------------------------------------------------

cat("Fitting V1 (all observed, outliers removed)...\n")
brt_v1  <- fit_brt(dat_v1, "logCHLa", preds_v1, tc = 4)
pred_v1 <- predict(brt_v1, newdata = dat_v1,
                   n.trees = brt_v1$gbm.call$best.trees)
r_v1    <- cor(dat_v1$logCHLa, pred_v1)
cat(sprintf("  V1: r=%.3f  R²=%.3f\n\n", r_v1, r_v1^2))

# ----------------------------------------------------------------------------
# 5. V2 — full chain
# ----------------------------------------------------------------------------

cat("Fitting V2 (full chain, outliers removed)...\n")
brt_v2  <- fit_brt(dat_v2, "logCHLa", preds_v2, tc = 4)
pred_v2 <- predict(brt_v2, newdata = dat_v2,
                   n.trees = brt_v2$gbm.call$best.trees)
r_v2    <- cor(dat_v2$logCHLa, pred_v2)
cat(sprintf("  V2: r=%.3f  R²=%.3f\n\n", r_v2, r_v2^2))

# ----------------------------------------------------------------------------
# 6. LOSO
# ----------------------------------------------------------------------------

cat("Running LOSO (outliers removed)...\n")
sites    <- sort(unique(dat_chain$Site))
loso_out <- vector("list", length(sites))

for (i in seq_along(sites)) {
  s     <- sites[i]
  train <- as.data.frame(dat_chain[dat_chain$Site != s, ])
  test  <- as.data.frame(dat_chain[dat_chain$Site == s, ])
  cat(sprintf("  Leaving out %s...\n", s))
  
  brt_SPC_l <- fit_brt(train, "SPC",       hydro_preds, tc = 3)
  brt_TP_l  <- fit_brt(train, "logTP_obs", hydro_preds, tc = 3)
  brt_TN_l  <- fit_brt(train, "logTN_obs", hydro_preds, tc = 3)
  if (any(sapply(list(brt_SPC_l, brt_TP_l, brt_TN_l), is.null))) next
  
  test$pred_SPC   <- predict(brt_SPC_l, newdata = test[ , hydro_preds],
                             n.trees = brt_SPC_l$gbm.call$best.trees)
  test$pred_logTP <- predict(brt_TP_l,  newdata = test[ , hydro_preds],
                             n.trees = brt_TP_l$gbm.call$best.trees)
  test$pred_logTN <- predict(brt_TN_l,  newdata = test[ , hydro_preds],
                             n.trees = brt_TN_l$gbm.call$best.trees)
  
  train$pred_SPC   <- predict(brt_SPC_l, newdata = train[ , hydro_preds],
                              n.trees = brt_SPC_l$gbm.call$best.trees)
  train$pred_logTP <- predict(brt_TP_l,  newdata = train[ , hydro_preds],
                              n.trees = brt_TP_l$gbm.call$best.trees)
  train$pred_logTN <- predict(brt_TN_l,  newdata = train[ , hydro_preds],
                              n.trees = brt_TN_l$gbm.call$best.trees)
  
  brt_bloom_l <- fit_brt(train, "logCHLa", bloom_preds, tc = 4)
  if (is.null(brt_bloom_l)) next
  
  test$pred_logCHLa <- predict(brt_bloom_l,
                               newdata = test[ , bloom_preds],
                               n.trees = brt_bloom_l$gbm.call$best.trees)
  
  r_s <- cor(test$logCHLa, test$pred_logCHLa, use = "complete.obs")
  cat(sprintf("    r = %.3f\n", r_s))
  
  loso_out[[i]] <- data.frame(
    Site = test$Site, Year = test$Year,
    Observed = test$logCHLa, Predicted = test$pred_logCHLa
  )
}

loso_all <- do.call(rbind, loso_out)
r_loso   <- cor(loso_all$Observed, loso_all$Predicted, use = "complete.obs")
cat(sprintf("\n  LOSO: r=%.3f  R²=%.3f\n\n", r_loso, r_loso^2))

# ----------------------------------------------------------------------------
# 7. Jackknife
# ----------------------------------------------------------------------------

cat("Running Jackknife (outliers removed)...\n")
years  <- sort(unique(dat_chain$Year))
jk_out <- vector("list", length(years))

for (i in seq_along(years)) {
  yr    <- years[i]
  train <- as.data.frame(dat_chain[dat_chain$Year != yr, ])
  test  <- as.data.frame(dat_chain[dat_chain$Year == yr, ])
  
  brt_SPC_j <- fit_brt(train, "SPC",       hydro_preds, tc = 3)
  brt_TP_j  <- fit_brt(train, "logTP_obs", hydro_preds, tc = 3)
  brt_TN_j  <- fit_brt(train, "logTN_obs", hydro_preds, tc = 3)
  if (any(sapply(list(brt_SPC_j, brt_TP_j, brt_TN_j), is.null))) next
  
  test$pred_SPC   <- predict(brt_SPC_j, newdata = test[ , hydro_preds],
                             n.trees = brt_SPC_j$gbm.call$best.trees)
  test$pred_logTP <- predict(brt_TP_j,  newdata = test[ , hydro_preds],
                             n.trees = brt_TP_j$gbm.call$best.trees)
  test$pred_logTN <- predict(brt_TN_j,  newdata = test[ , hydro_preds],
                             n.trees = brt_TN_j$gbm.call$best.trees)
  
  train$pred_SPC   <- predict(brt_SPC_j, newdata = train[ , hydro_preds],
                              n.trees = brt_SPC_j$gbm.call$best.trees)
  train$pred_logTP <- predict(brt_TP_j,  newdata = train[ , hydro_preds],
                              n.trees = brt_TP_j$gbm.call$best.trees)
  train$pred_logTN <- predict(brt_TN_j,  newdata = train[ , hydro_preds],
                              n.trees = brt_TN_j$gbm.call$best.trees)
  
  brt_bloom_j <- fit_brt(train, "logCHLa", bloom_preds, tc = 4)
  if (is.null(brt_bloom_j)) next
  
  test$pred_logCHLa <- predict(brt_bloom_j,
                               newdata = test[ , bloom_preds],
                               n.trees = brt_bloom_j$gbm.call$best.trees)
  
  jk_out[[i]] <- data.frame(
    Site = test$Site, Year = test$Year,
    Observed = test$logCHLa, Predicted = test$pred_logCHLa
  )
}

jk_all <- do.call(rbind, jk_out)
r_jk   <- cor(jk_all$Observed, jk_all$Predicted, use = "complete.obs")
cat(sprintf("\n  Jackknife: r=%.3f  R²=%.3f\n\n", r_jk, r_jk^2))

# ----------------------------------------------------------------------------
# 8. Summary
# ----------------------------------------------------------------------------

cat("══════════════════════════════════════════════════════\n")
cat("Performance Summary — 11 outlier observations removed\n")
cat("══════════════════════════════════════════════════════\n")
cat(sprintf("  %-35s  %6s  %6s\n", "Model", "r", "R²"))
cat(paste(rep("-", 52), collapse = ""), "\n")
cat(sprintf("  %-35s  %6.3f  %6.3f\n", "V1 training (all observed)",  r_v1,  r_v1^2))
cat(sprintf("  %-35s  %6.3f  %6.3f\n", "V2 training (full chain)",    r_v2,  r_v2^2))
cat(sprintf("  %-35s  %6.3f  %6.3f\n", "LOSO (spatial)",              r_loso, r_loso^2))
cat(sprintf("  %-35s  %6.3f  %6.3f\n", "Jackknife (temporal)",        r_jk,  r_jk^2))
cat("\nAll observations comparison:\n")
cat("  V1:        R²=0.830\n")
cat("  V2:        R²=0.709\n")
cat("  LOSO:      R²=0.279\n")
cat("  Jackknife: R²=0.404\n")
cat("Done.\n")

