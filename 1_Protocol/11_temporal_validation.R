# ============================================================================
# 11_temporal_validation.R
# UCFR Filamentous Algae Project
# Temporal validation of M1 bloom model — four schemes in one pass
#
# Inputs:  3_models/bloom_model_M1.rds            (fitted M1 object)
#          2_incremental/m1_predictions.csv        (clean obs_ids from M1)
#          2_incremental/ucfr_model_ready.csv      (full dataset)
#
# Outputs: 2_incremental/temporal_val_predictions.csv
#          4_products/diagnostics/temporal_validation.pdf
#          console: scorecard data frame
#
# Four validation schemes:
#
#   A. LOYO jackknife — leave one year out, predict with RE retained.
#      Operationally relevant: known sites, unknown future year.
#      Primary validation claim for the paper.
#
#   B. Forward-chaining — train on years 1:t, predict t+1, expanding window.
#      Respects temporal ordering. Key secondary test because the lag makes
#      this a genuine conditional forecast.
#
#   C. Lag-only baseline — fit s(lag_y) + s(Site, bs="re") only.
#      Persistence benchmark. Partitions M1 skill into:
#        - autocorrelation component (lag-only R²)
#        - genuine mechanism (M1 R² minus lag-only R²)
#      Without this a referee will assume most skill is bloom persistence.
#
#   D. Recursive-mode validation — seed one true lag_y per site, feed
#      model predictions forward as lag_y for subsequent years.
#      Reports RMSE as a function of steps forward (error-growth curve).
#      This curve empirically defines the trustworthy projection horizon
#      and converts the recursive-lag problem from a hidden flaw into a
#      stated, quantified finding.
#
# LOSO framing note (from 08_bloom_model_M1.R):
#   LOSO is NOT reported here. LOSO tests generalization to a novel site
#   (RE excluded), which is never the actual projection task — all seven
#   UCFR sites are known. LOSO lives in 08_bloom_model_M1.R as a
#   supplementary transparency item only. LOYO is the primary claim.
#
# Clean dataset dependency:
#   This script uses exactly the 221 observations M1 was trained on.
#   Row identity comes from obs_id in m1_predictions.csv (in_sample rows).
#   This guarantees consistency with M1 — no independent outlier removal.
# ============================================================================

library(mgcv)

# ============================================================================
# CONFIGURATION
# ============================================================================
MIN_TRAIN_YRS <- 10    # minimum training years for forward-chaining
# ============================================================================


# ============================================================================
# 1. LOAD MODEL AND CLEAN DATASET
# ============================================================================
m_M1   <- readRDS("3_models/bloom_model_M1.rds")
f_M1   <- formula(m_M1)
RESPONSE <- as.character(f_M1[[2]])

# Identify which observations M1 was trained on
preds_csv <- read.csv("2_incremental/m1_predictions.csv",
                      stringsAsFactors = FALSE)
clean_ids <- preds_csv$obs_id[preds_csv$scheme == "in_sample"]

# Full dataset — needed to reconstruct predictor columns
dat_full <- read.csv("2_incremental/ucfr_model_ready.csv",
                     stringsAsFactors = FALSE)
dat_full$obs_id <- seq_len(nrow(dat_full))

# Derive logAFDM if needed
if (RESPONSE == "logAFDM" && !"logAFDM" %in% names(dat_full)) {
  dat_full$logAFDM[dat_full$AFDM <= 0] <- NA
  dat_full$logAFDM <- log10(dat_full$AFDM)
}

# Derive lag_y (same logic as 08_bloom_model_M1.R)
ann_max <- aggregate(dat_full[[RESPONSE]] ~ Site + Year,
                     data = dat_full, FUN = max, na.rm = TRUE)
names(ann_max) <- c("Site", "Year", "annual_max")
lag_df <- ann_max
lag_df$Year <- lag_df$Year + 1
names(lag_df)[names(lag_df) == "annual_max"] <- "lag_y"
dat_full <- merge(dat_full,
                  lag_df[, c("Site", "Year", "lag_y")],
                  by = c("Site", "Year"), all.x = TRUE)
dat_full <- dat_full[order(dat_full$obs_id), ]

# Extract the clean M1 training set by obs_id
PREDICTORS <- c("lag_y", "anomaly", "logQ_obs_cfs",
                "Days_Since_Freshet", "logTP_mg_L", "Temp_oC")
keep_cols  <- c("obs_id", "Site", "Year", "Month", RESPONSE, PREDICTORS)
mdat       <- dat_full[dat_full$obs_id %in% clean_ids, keep_cols]
mdat$Site  <- factor(mdat$Site)

# Sanity check
if (nrow(mdat) != length(clean_ids)) {
  stop("obs_id mismatch: expected ", length(clean_ids),
       " rows, got ", nrow(mdat),
       ". Re-run 08_bloom_model_M1.R to regenerate m1_predictions.csv.")
}

years <- sort(unique(mdat$Year))
sites <- levels(mdat$Site)


# ============================================================================
# A. LEAVE-ONE-YEAR-OUT (LOYO) JACKKNIFE
# ============================================================================
# RE retained: we are predicting at known sites in unknown years.
# This is the operationally relevant test and the primary paper claim.

loyo_list <- vector("list", length(years))
names(loyo_list) <- as.character(years)

for (y in years) {
  train <- mdat[mdat$Year != y, ]
  test  <- mdat[mdat$Year == y, ]
  if (nrow(test) == 0) next
  
  m_loyo <- tryCatch(
    gam(f_M1, data = train, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(m_loyo)) next
  
  # For sites present in training: predict with RE
  # For sites absent from training: predict at population level (exclude RE)
  # In practice all 7 sites appear in every year except early years —
  # handle generically so the loop is robust.
  train_sites <- levels(droplevels(train$Site))
  test2       <- test
  test2$Site  <- factor(as.character(test2$Site), levels = train_sites)
  
  in_train  <- !is.na(test2$Site)
  pred_loyo <- rep(NA_real_, nrow(test))
  
  if (any(in_train)) {
    pred_loyo[in_train] <- predict(m_loyo,
                                   newdata = test2[in_train, ])
  }
  if (any(!in_train)) {
    test_oos        <- test[!in_train, ]
    test_oos$Site   <- factor(train_sites[1], levels = train_sites)
    pred_loyo[!in_train] <- predict(m_loyo, newdata = test_oos,
                                    exclude = 's(Site, bs="re")',
                                    newdata.guaranteed = TRUE)
  }
  
  obs_loyo <- test[[RESPONSE]]
  loyo_list[[as.character(y)]] <- data.frame(
    scheme    = "LOYO",
    obs_id    = test$obs_id,
    Year      = test$Year,
    Site      = as.character(test$Site),
    Month     = test$Month,
    Observed  = obs_loyo,
    Predicted = pred_loyo,
    stringsAsFactors = FALSE
  )
}

loyo_df       <- do.call(rbind, loyo_list)
row.names(loyo_df) <- NULL

loyo_ss_res   <- sum((loyo_df$Observed - loyo_df$Predicted)^2, na.rm = TRUE)
loyo_ss_tot   <- sum((loyo_df$Observed - mean(loyo_df$Observed, na.rm = TRUE))^2,
                     na.rm = TRUE)
r2_loyo       <- 1 - loyo_ss_res / loyo_ss_tot
rmse_loyo     <- sqrt(mean((loyo_df$Observed - loyo_df$Predicted)^2, na.rm = TRUE))

# Per-year LOYO summary
loyo_yr <- do.call(rbind, lapply(split(loyo_df, loyo_df$Year), function(d) {
  ss_r <- sum((d$Observed - d$Predicted)^2, na.rm = TRUE)
  ss_t <- sum((d$Observed - mean(d$Observed, na.rm = TRUE))^2, na.rm = TRUE)
  data.frame(
    Year = d$Year[1],
    n    = nrow(d),
    R2   = round(ifelse(ss_t > 0, 1 - ss_r / ss_t, NA), 4),
    RMSE = round(sqrt(mean((d$Observed - d$Predicted)^2, na.rm = TRUE)), 4)
  )
}))
row.names(loyo_yr) <- NULL


# ============================================================================
# B. FORWARD-CHAINING (expanding window)
# ============================================================================
# Train on years 1:t, predict t+1. RE retained.
# Key secondary test — respects temporal ordering.

first_test <- years[1] + MIN_TRAIN_YRS
test_years <- years[years >= first_test]

fc_list <- vector("list", length(test_years))
names(fc_list) <- as.character(test_years)

for (y in test_years) {
  train <- mdat[mdat$Year < y, ]
  test  <- mdat[mdat$Year == y, ]
  if (nrow(test) == 0) next
  
  # Sites in test but not in training get population-level prediction
  train_sites <- levels(droplevels(train$Site))
  test2       <- test
  test2$Site  <- factor(as.character(test2$Site), levels = train_sites)
  
  m_fc <- tryCatch(
    gam(f_M1, data = train, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(m_fc)) next
  
  in_train <- !is.na(test2$Site)
  pred_fc  <- rep(NA_real_, nrow(test))
  
  if (any(in_train)) {
    pred_fc[in_train] <- predict(m_fc, newdata = test2[in_train, ])
  }
  if (any(!in_train)) {
    test_oos       <- test[!in_train, ]
    test_oos$Site  <- factor(train_sites[1], levels = train_sites)
    pred_fc[!in_train] <- predict(m_fc, newdata = test_oos,
                                  exclude = 's(Site, bs="re")',
                                  newdata.guaranteed = TRUE)
  }
  
  obs_fc <- test[[RESPONSE]]
  fc_list[[as.character(y)]] <- data.frame(
    scheme    = "FwdChain",
    obs_id    = test$obs_id,
    Year      = test$Year,
    Site      = as.character(test$Site),
    Month     = test$Month,
    Observed  = obs_fc,
    Predicted = pred_fc,
    stringsAsFactors = FALSE
  )
}

fc_df         <- do.call(rbind, fc_list)
row.names(fc_df) <- NULL

fc_ss_res     <- sum((fc_df$Observed - fc_df$Predicted)^2, na.rm = TRUE)
fc_ss_tot     <- sum((fc_df$Observed - mean(fc_df$Observed, na.rm = TRUE))^2,
                     na.rm = TRUE)
r2_fc         <- 1 - fc_ss_res / fc_ss_tot
rmse_fc       <- sqrt(mean((fc_df$Observed - fc_df$Predicted)^2, na.rm = TRUE))

# Per-year forward-chain summary
fc_yr <- do.call(rbind, lapply(split(fc_df, fc_df$Year), function(d) {
  y    <- d$Year[1]
  ss_r <- sum((d$Observed - d$Predicted)^2, na.rm = TRUE)
  ss_t <- sum((d$Observed - mean(d$Observed, na.rm = TRUE))^2, na.rm = TRUE)
  data.frame(
    Year      = y,
    n_train   = sum(mdat$Year < y),
    n_test    = nrow(d),
    R2        = round(ifelse(ss_t > 0, 1 - ss_r / ss_t, NA), 4),
    RMSE      = round(sqrt(mean((d$Observed - d$Predicted)^2, na.rm = TRUE)), 4)
  )
}))
row.names(fc_yr) <- NULL


# ============================================================================
# C. LAG-ONLY BASELINE MODEL
# ============================================================================
# Fits s(lag_y, k=5) + s(Site, bs="re") on the same clean dataset.
# Partitions M1 skill:
#   - lag-only R²       = persistence component (bloom autocorrelation)
#   - M1 R² - lag R²   = mechanistic component (genuine environmental signal)
#
# Also runs LOYO on the lag-only model so the persistence vs mechanism
# decomposition holds at the validation level, not just in-sample.

f_lag <- as.formula(paste0(RESPONSE, " ~ s(lag_y, k=5) + s(Site, bs=\"re\")"))

m_lag       <- gam(f_lag, data = mdat, method = "REML")
r2_lag_insamp <- summary(m_lag)$r.sq

# M1 in-sample R² for comparison (recompute from fitted values on mdat)
m_M1_refit  <- gam(f_M1, data = mdat, method = "REML")
r2_M1_insamp <- summary(m_M1_refit)$r.sq
r2_mechanism_insamp <- r2_M1_insamp - r2_lag_insamp

# LOYO on lag-only model
lag_loyo_list <- vector("list", length(years))
for (y in years) {
  train <- mdat[mdat$Year != y, ]
  test  <- mdat[mdat$Year == y, ]
  if (nrow(test) == 0) next
  
  train_sites <- levels(droplevels(train$Site))
  test2       <- test
  test2$Site  <- factor(as.character(test2$Site), levels = train_sites)
  
  m_lag_loyo <- tryCatch(
    gam(f_lag, data = train, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(m_lag_loyo)) next
  
  in_train <- !is.na(test2$Site)
  pred_lag <- rep(NA_real_, nrow(test))
  if (any(in_train))
    pred_lag[in_train] <- predict(m_lag_loyo, newdata = test2[in_train, ])
  if (any(!in_train)) {
    test_oos      <- test[!in_train, ]
    test_oos$Site <- factor(train_sites[1], levels = train_sites)
    pred_lag[!in_train] <- predict(m_lag_loyo, newdata = test_oos,
                                   exclude = 's(Site, bs="re")',
                                   newdata.guaranteed = TRUE)
  }
  
  lag_loyo_list[[as.character(y)]] <- data.frame(
    Observed  = test[[RESPONSE]],
    Predicted = pred_lag,
    stringsAsFactors = FALSE
  )
}

lag_loyo_df   <- do.call(rbind, lag_loyo_list)
lag_ss_res    <- sum((lag_loyo_df$Observed - lag_loyo_df$Predicted)^2, na.rm = TRUE)
lag_ss_tot    <- sum((lag_loyo_df$Observed - mean(lag_loyo_df$Observed, na.rm = TRUE))^2,
                     na.rm = TRUE)
r2_lag_loyo   <- 1 - lag_ss_res / lag_ss_tot
r2_mechanism_loyo <- r2_loyo - r2_lag_loyo


# ============================================================================
# D. RECURSIVE-MODE VALIDATION
# ============================================================================
# Tests error propagation under the self-feeding lag architecture used
# in climate projections: seed year t with observed lag_y, predict t+1,
# feed predicted value as lag_y for t+2, and so on.
#
# Run per site so error-growth curves are interpretable per monitoring location.
# Aggregate across sites for a pooled RMSE-by-step summary.
#
# Uses M1 fit on full clean dataset (m_M1_refit above).
# Prediction uses RE for all sites (known sites — same as LOYO).
#
# Output: for each site, a sequence of (step, observed, predicted_recursive).
# Step = number of years ahead from the seed year.

# For recursive mode we need one observation per site per year (annual max,
# since lag_y is defined as the prior year's annual max).
ann_obs <- aggregate(mdat[[RESPONSE]] ~ Site + Year, data = mdat,
                     FUN = max, na.rm = TRUE)
names(ann_obs) <- c("Site", "Year", "y_obs")

# Representative predictor values per site-year: use the observation with the
# higher response value within each Site-Year (i.e., the one that defines
# the annual max — its predictor context is what matters for the max bloom).
mdat_ann <- do.call(rbind, lapply(split(mdat, list(mdat$Site, mdat$Year),
                                        drop = TRUE), function(d) {
                                          d[which.max(d[[RESPONSE]]), ]
                                        }))
mdat_ann <- mdat_ann[order(mdat_ann$Site, mdat_ann$Year), ]

rec_list <- vector("list", length(sites))
names(rec_list) <- sites

for (s in sites) {
  site_dat  <- mdat_ann[mdat_ann$Site == s, ]
  site_yrs  <- sort(site_dat$Year)
  if (length(site_yrs) < 3) next   # need at least seed + 2 steps
  
  # Seed: first year with observed lag_y (i.e., second year in sequence)
  # Use the observed lag_y for the seed year, then propagate recursively.
  seed_idx <- which(!is.na(site_dat$lag_y))[1]
  if (is.na(seed_idx)) next
  
  n_steps   <- nrow(site_dat) - seed_idx
  if (n_steps < 1) next
  
  rec_rows  <- vector("list", n_steps)
  lag_carry <- site_dat$lag_y[seed_idx]   # true observed lag for seed year
  
  for (step in seq_len(n_steps)) {
    row_idx   <- seed_idx + step - 1
    pred_row  <- site_dat[row_idx, ]
    pred_row$lag_y <- lag_carry
    
    pred_val  <- tryCatch(
      as.numeric(predict(m_M1_refit, newdata = pred_row)),
      error = function(e) NA_real_
    )
    
    rec_rows[[step]] <- data.frame(
      Site            = s,
      Year            = pred_row$Year,
      step            = step,
      lag_y_used      = round(lag_carry, 4),
      Observed        = round(pred_row[[RESPONSE]], 4),
      Predicted_rec   = round(pred_val, 4),
      Error           = round(pred_row[[RESPONSE]] - pred_val, 4),
      stringsAsFactors = FALSE
    )
    
    # Feed prediction forward as next year's lag_y
    lag_carry <- pred_val
  }
  
  rec_list[[s]] <- do.call(rbind, rec_rows)
}

rec_df <- do.call(rbind, rec_list)
row.names(rec_df) <- NULL

# RMSE by step (pooled across sites)
max_steps    <- max(rec_df$step, na.rm = TRUE)
rmse_by_step <- do.call(rbind, lapply(seq_len(max_steps), function(k) {
  d <- rec_df[rec_df$step == k & !is.na(rec_df$Error), ]
  data.frame(
    step     = k,
    n_sites  = nrow(d),
    RMSE     = round(sqrt(mean(d$Error^2)), 4),
    MAE      = round(mean(abs(d$Error)), 4)
  )
}))


# ============================================================================
# 5. SAVE OUTPUTS
# ============================================================================
if (!dir.exists("2_incremental")) dir.create("2_incremental", recursive = TRUE)

# Combined predictions (LOYO + FwdChain; recursive saved separately)
all_preds <- rbind(
  loyo_df[, c("scheme","obs_id","Year","Site","Month","Observed","Predicted")],
  fc_df[,   c("scheme","obs_id","Year","Site","Month","Observed","Predicted")]
)
write.csv(all_preds, "2_incremental/temporal_val_predictions.csv",
          row.names = FALSE)

write.csv(rec_df,        "2_incremental/recursive_val_predictions.csv",
          row.names = FALSE)
write.csv(rmse_by_step,  "2_incremental/recursive_rmse_by_step.csv",
          row.names = FALSE)


# ============================================================================
# 6. DIAGNOSTIC PLOTS
# ============================================================================
if (!dir.exists("4_products/diagnostics"))
  dir.create("4_products/diagnostics", recursive = TRUE)

pdf("4_products/diagnostics/temporal_validation.pdf", width = 10, height = 8)

# --- Page 1: Obs vs pred, R² by year ---
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

plot(loyo_df$Predicted, loyo_df$Observed,
     pch = 16, cex = 0.7, col = adjustcolor("steelblue", 0.6),
     xlab = "Predicted", ylab = "Observed",
     main = paste0("LOYO (R² = ", round(r2_loyo, 3), ")"))
abline(0, 1, lty = 2, col = "red")

plot(fc_df$Predicted, fc_df$Observed,
     pch = 16, cex = 0.7, col = adjustcolor("darkgreen", 0.6),
     xlab = "Predicted", ylab = "Observed",
     main = paste0("Forward-Chain (R² = ", round(r2_fc, 3), ")"))
abline(0, 1, lty = 2, col = "red")

ylim_loyo <- c(min(c(loyo_yr$R2, 0), na.rm = TRUE),
               max(c(loyo_yr$R2, 1), na.rm = TRUE))
plot(loyo_yr$Year, loyo_yr$R2, type = "b", pch = 16,
     xlab = "Year", ylab = "R²", ylim = ylim_loyo,
     main = "LOYO R² by Year")
abline(h = 0,         lty = 3, col = "grey50")
abline(h = r2_loyo,   lty = 2, col = "steelblue")
text(max(loyo_yr$Year), r2_loyo, "pooled", pos = 2,
     col = "steelblue", cex = 0.8)

ylim_fc <- c(min(c(fc_yr$R2, 0), na.rm = TRUE),
             max(c(fc_yr$R2, 1), na.rm = TRUE))
plot(fc_yr$Year, fc_yr$R2, type = "b", pch = 16,
     xlab = "Year", ylab = "R²", ylim = ylim_fc,
     main = "Forward-Chain R² by Year")
abline(h = 0,       lty = 3, col = "grey50")
abline(h = r2_fc,   lty = 2, col = "darkgreen")
text(max(fc_yr$Year), r2_fc, "pooled", pos = 2,
     col = "darkgreen", cex = 0.8)

# --- Page 2: Residual time series ---
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

loyo_df$Resid <- loyo_df$Observed - loyo_df$Predicted
plot(jitter(loyo_df$Year, amount = 0.2), loyo_df$Resid,
     pch = 16, cex = 0.6, col = adjustcolor("steelblue", 0.5),
     xlab = "Year", ylab = "Residual", main = "LOYO Residuals by Year")
abline(h = 0, lty = 2)
loyo_yr_mean <- aggregate(Resid ~ Year, data = loyo_df, FUN = mean)
lines(loyo_yr_mean$Year, loyo_yr_mean$Resid, col = "red", lwd = 2)

fc_df$Resid <- fc_df$Observed - fc_df$Predicted
plot(jitter(fc_df$Year, amount = 0.2), fc_df$Resid,
     pch = 16, cex = 0.6, col = adjustcolor("darkgreen", 0.5),
     xlab = "Year", ylab = "Residual", main = "Forward-Chain Residuals by Year")
abline(h = 0, lty = 2)
fc_yr_mean <- aggregate(Resid ~ Year, data = fc_df, FUN = mean)
lines(fc_yr_mean$Year, fc_yr_mean$Resid, col = "red", lwd = 2)

# --- Page 3: Recursive RMSE growth curve ---
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))

# Also show LOYO RMSE as a horizontal reference (non-recursive baseline)
plot(rmse_by_step$step, rmse_by_step$RMSE,
     type = "b", pch = 16, col = "darkorange", lwd = 2,
     xlab = "Steps ahead (years)", ylab = "RMSE",
     main = "Recursive-Mode Error Growth",
     ylim = c(0, max(rmse_by_step$RMSE, rmse_loyo) * 1.15))
abline(h = rmse_loyo, lty = 2, col = "steelblue")
text(max(rmse_by_step$step), rmse_loyo,
     paste0("LOYO RMSE = ", round(rmse_loyo, 3)),
     pos = 2, col = "steelblue", cex = 0.85)
# Add n_sites labels
text(rmse_by_step$step, rmse_by_step$RMSE,
     paste0("n=", rmse_by_step$n_sites),
     pos = 3, cex = 0.7, col = "grey40")

# --- Page 4: Recursive predictions per site ---
n_sites <- length(sites)
n_cols  <- 2
n_rows  <- ceiling(n_sites / n_cols)
par(mfrow = c(n_rows, n_cols), mar = c(3, 3, 2, 1))

for (s in sites) {
  d <- rec_df[rec_df$Site == s, ]
  if (nrow(d) == 0) next
  ylim_s <- range(c(d$Observed, d$Predicted_rec), na.rm = TRUE)
  plot(d$Year, d$Observed,
       type = "b", pch = 16, col = "black", lwd = 1.5,
       xlab = "", ylab = RESPONSE,
       main = s, ylim = ylim_s)
  lines(d$Year, d$Predicted_rec,
        type = "b", pch = 1, col = "darkorange", lty = 2, lwd = 1.5)
  legend("topleft", legend = c("Observed", "Recursive pred"),
         col = c("black", "darkorange"), lty = c(1, 2), pch = c(16, 1),
         cex = 0.7, bty = "n")
}

dev.off()


# ============================================================================
# 7. SCORECARD
# ============================================================================
cat("\n")
cat("============================================================\n")
cat("TEMPORAL VALIDATION SCORECARD — M1\n")
cat("Response:", RESPONSE, "\n")
cat("Clean n (from 08_bloom_model_M1.R):", nrow(mdat), "\n")
cat("============================================================\n\n")

# --- A. Summary table ---
sc_main <- data.frame(
  Scheme = c(
    "M1 in-sample (refit on clean n)",
    "A. LOYO jackknife          [PRIMARY]",
    "B. Forward-chaining",
    "C. Lag-only baseline (LOYO)",
    "   Mechanism component (A minus C)"
  ),
  n = c(
    nrow(mdat),
    nrow(loyo_df),
    nrow(fc_df),
    nrow(lag_loyo_df),
    NA
  ),
  R2_or_delta = c(
    round(r2_M1_insamp,         4),
    round(r2_loyo,              4),
    round(r2_fc,                4),
    round(r2_lag_loyo,          4),
    round(r2_mechanism_loyo,    4)
  ),
  RMSE = c(
    round(sqrt(mean(residuals(m_M1_refit)^2)), 4),
    round(rmse_loyo, 4),
    round(rmse_fc,   4),
    round(sqrt(lag_ss_res / nrow(lag_loyo_df)), 4),
    NA
  ),
  stringsAsFactors = FALSE
)
cat("--- Overall performance ---\n")
print(sc_main, row.names = FALSE)
cat("\n")
cat("Note: 'Mechanism component' = LOYO R²(M1) - LOYO R²(lag-only).\n")
cat("      Fraction of held-out skill attributable to environmental\n")
cat("      predictors beyond pure bloom autocorrelation.\n\n")

# --- B. LOYO by year ---
cat("--- A. LOYO by year ---\n")
print(loyo_yr, row.names = FALSE)
cat("\n")

# --- C. Forward-chain by year ---
cat("--- B. Forward-chain by year ---\n")
print(fc_yr, row.names = FALSE)
cat("\n")

# --- D. Recursive RMSE by step ---
cat("--- D. Recursive-mode RMSE by steps ahead ---\n")
cat("(LOYO RMSE =", round(rmse_loyo, 4),
    "— shown as non-recursive reference)\n")
print(rmse_by_step, row.names = FALSE)
cat("\n")
cat("Interpretation: the step at which recursive RMSE exceeds LOYO RMSE\n")
cat("is the empirical limit of the trustworthy projection horizon.\n\n")

# --- E. LOSO note ---
cat("--- LOSO note ---\n")
cat("LOSO is not reported here. LOSO tests generalization to a novel\n")
cat("site (RE excluded) — never the actual projection task. All seven\n")
cat("UCFR sites are known. LOSO lives in 08_bloom_model_M1.R as a\n")
cat("supplementary transparency item. LOYO is the primary claim.\n\n")

# --- F. Output files ---
sc_files <- data.frame(
  File = c(
    "2_incremental/temporal_val_predictions.csv",
    "2_incremental/recursive_val_predictions.csv",
    "2_incremental/recursive_rmse_by_step.csv",
    "4_products/diagnostics/temporal_validation.pdf"
  ),
  Contents = c(
    "LOYO and forward-chain obs/predicted rows.",
    "Per-site recursive predictions with step index.",
    "Pooled RMSE and MAE by steps-ahead.",
    "Diagnostic plots: obs/pred, R² by year, recursive growth curve."
  ),
  stringsAsFactors = FALSE
)
cat("--- Output files ---\n")
print(sc_files, row.names = FALSE)
cat("\n")
cat("============================================================\n")
cat("Done.\n")
cat("============================================================\n")