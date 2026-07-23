# ============================================================================
# 08_temporal_validation.R
# UCFR Filamentous Algae Project
# Temporal validation of lagged GAM bloom model
#
# Input:    2_incremental/ucfr_model_ready.csv
# Outputs:  2_incremental/temporal_val_predictions.csv
#           4_products/diagnostics/temporal_validation.pdf
#
# Two validation schemes:
#   (1) Leave-one-year-out (LOYO) jackknife — all years except Y used to
#       train, predict held-out Y. Symmetric — every year gets a turn.
#   (2) Forward-chaining — train on years 1:t, predict t+1. Expanding
#       window with minimum training period. Respects temporal ordering.
#       Key test because the lag makes this a genuine conditional forecast.
#
# Lag note: lag_y = previous year's observed annual max at same site.
#   This is honestly available at prediction time (you know what happened
#   last year), so both schemes use observed lag. No predicted lag chaining.
#
# Same response toggle and outlier removal logic as core model scripts.
# ============================================================================

library(mgcv)

# --- CONFIGURATION ---------------------------------------------------------
RESPONSE       <- "logCHLa"    # toggle: "logCHLa" or "logAFDM"
K              <- 5            # basis dimension for smooths
OUTLIER_SD     <- 2.0          # flag residuals beyond ±SD
OUTLIER_ROUNDS <- 2            # number of outlier removal rounds
MIN_TRAIN_YRS  <- 10           # minimum training years for forward-chaining
# ---------------------------------------------------------------------------

# --- 1. READ & PREP -------------------------------------------------------
dat <- read.csv("2_incremental/ucfr_model_ready.csv",
                stringsAsFactors = FALSE)
names(dat) <- gsub("/", "_over_", names(dat))

# Derive lag: previous year's annual max at same site
ann_max <- aggregate(dat[[RESPONSE]] ~ Site + Year, data = dat,
                     FUN = max, na.rm = TRUE)
names(ann_max) <- c("Site", "Year", "annual_max")
ann_max$lag_Year <- ann_max$Year + 1
lag_df <- ann_max[, c("Site", "lag_Year", "annual_max")]
names(lag_df) <- c("Site", "Year", "lag_y")
dat <- merge(dat, lag_df, by = c("Site", "Year"), all.x = TRUE)

# Model formula — matches M1 from iterative selection
PREDICTORS <- c("lag_y", "anomaly", "logQ_obs_cfs",
                "Days_Since_Freshet", "logTP_mg_L", "Temp_oC")
f1 <- as.formula(paste0(RESPONSE, " ~ ",
                        paste0("s(", PREDICTORS, ", k=", K, ")",
                               collapse = " + ")))

# Complete cases only
keep_cols <- c(RESPONSE, PREDICTORS, "Site", "Year", "Month")
mdat <- dat[complete.cases(dat[, c(RESPONSE, PREDICTORS)]), keep_cols]
mdat$Site <- factor(mdat$Site)

cat("Response:", RESPONSE, "\n")
cat("Formula:", deparse(f1, width.cutoff = 200), "\n")
cat("Complete cases before outlier removal:", nrow(mdat), "\n")
cat("Year range:", min(mdat$Year), "-", max(mdat$Year), "\n")
cat("Sites:", paste(levels(mdat$Site), collapse = ", "), "\n\n")

# --- 2. OUTLIER REMOVAL (same logic as core scripts) ----------------------
for (r in seq_len(OUTLIER_ROUNDS)) {
  m_tmp <- gam(f1, data = mdat, method = "REML")
  resid_tmp <- residuals(m_tmp)
  sd_tmp <- sd(resid_tmp)
  outliers <- which(abs(resid_tmp) > OUTLIER_SD * sd_tmp)
  
  if (length(outliers) == 0) {
    cat("Outlier round", r, ": none found. Dataset stable.\n")
    break
  }
  
  cat("Outlier round", r, ": removing", length(outliers), "observations",
      "(SD threshold:", round(OUTLIER_SD * sd_tmp, 3), ")\n")
  
  # Log which obs are removed
  removed <- mdat[outliers, c("Site", "Year", "Month", RESPONSE)]
  removed$resid <- round(resid_tmp[outliers], 3)
  print(removed, row.names = FALSE)
  cat("\n")
  
  mdat <- mdat[-outliers, ]
}

cat("Clean dataset:", nrow(mdat), "observations\n")
cat("Years retained:", length(unique(mdat$Year)),
    "(", min(mdat$Year), "-", max(mdat$Year), ")\n\n")

# Refit on clean data for reference R²
m_full <- gam(f1, data = mdat, method = "REML")
cat("Full-data R² (adj):", round(summary(m_full)$r.sq, 4), "\n")
cat("Full-data dev expl:", round(summary(m_full)$dev.expl * 100, 1), "%\n\n")


# ============================================================================
# 3. LEAVE-ONE-YEAR-OUT (LOYO) JACKKNIFE
# ============================================================================
cat("====================================================\n")
cat("LEAVE-ONE-YEAR-OUT JACKKNIFE\n")
cat("====================================================\n\n")

years <- sort(unique(mdat$Year))
loyo_preds <- list()

for (y in years) {
  train <- mdat[mdat$Year != y, ]
  test  <- mdat[mdat$Year == y, ]
  
  if (nrow(test) == 0) next
  
  m_loyo <- tryCatch(
    gam(f1, data = train, method = "REML"),
    error = function(e) NULL
  )
  
  if (is.null(m_loyo)) {
    cat("  Year", y, ": GAM failed to converge. Skipping.\n")
    next
  }
  
  pred <- predict(m_loyo, newdata = test)
  obs  <- test[[RESPONSE]]
  
  loyo_preds[[as.character(y)]] <- data.frame(
    scheme    = "LOYO",
    Year      = test$Year,
    Site      = test$Site,
    Month     = test$Month,
    Observed  = obs,
    Predicted = pred,
    stringsAsFactors = FALSE
  )
  
  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  r2 <- ifelse(ss_tot > 0, 1 - ss_res / ss_tot, NA)
  cat("  Year", y, ": n =", nrow(test),
      " R² =", round(r2, 3),
      " RMSE =", round(sqrt(mean((obs - pred)^2)), 3), "\n")
}

loyo_df <- do.call(rbind, loyo_preds)
row.names(loyo_df) <- NULL

# Pooled LOYO R² (vs grand mean of observed)
loyo_ss_res <- sum((loyo_df$Observed - loyo_df$Predicted)^2)
loyo_ss_tot <- sum((loyo_df$Observed - mean(loyo_df$Observed))^2)
loyo_r2_pooled <- 1 - loyo_ss_res / loyo_ss_tot
loyo_rmse <- sqrt(mean((loyo_df$Observed - loyo_df$Predicted)^2))

cat("\nLOYO pooled:  R² =", round(loyo_r2_pooled, 4),
    "  RMSE =", round(loyo_rmse, 4),
    "  n =", nrow(loyo_df), "\n\n")


# ============================================================================
# 4. FORWARD-CHAINING (expanding window)
# ============================================================================
cat("====================================================\n")
cat("FORWARD-CHAINING (min training window:", MIN_TRAIN_YRS, "years)\n")
cat("====================================================\n\n")

# First test year: min_year + MIN_TRAIN_YRS
first_test <- years[1] + MIN_TRAIN_YRS
test_years <- years[years >= first_test]

if (length(test_years) == 0) {
  cat("ERROR: not enough years for forward-chaining with",
      MIN_TRAIN_YRS, "year minimum.\n")
} else {
  cat("Training starts:", years[1], "\n")
  cat("First test year:", first_test, "\n")
  cat("Test years:", paste(test_years, collapse = ", "), "\n\n")
}

fc_preds <- list()

for (y in test_years) {
  train <- mdat[mdat$Year < y, ]
  test  <- mdat[mdat$Year == y, ]
  
  if (nrow(test) == 0) next
  
  n_train_yrs <- length(unique(train$Year))
  
  m_fc <- tryCatch(
    gam(f1, data = train, method = "REML"),
    error = function(e) NULL
  )
  
  if (is.null(m_fc)) {
    cat("  Year", y, ": GAM failed (", nrow(train), "train obs,",
        n_train_yrs, "years). Skipping.\n")
    next
  }
  
  pred <- predict(m_fc, newdata = test)
  obs  <- test[[RESPONSE]]
  
  fc_preds[[as.character(y)]] <- data.frame(
    scheme    = "FwdChain",
    Year      = test$Year,
    Site      = test$Site,
    Month     = test$Month,
    Observed  = obs,
    Predicted = pred,
    stringsAsFactors = FALSE
  )
  
  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  r2 <- ifelse(ss_tot > 0, 1 - ss_res / ss_tot, NA)
  cat("  Year", y, ": train =", nrow(train), "obs (",
      n_train_yrs, "yrs)  test =", nrow(test),
      "  R² =", round(r2, 3),
      "  RMSE =", round(sqrt(mean((obs - pred)^2)), 3), "\n")
}

fc_df <- do.call(rbind, fc_preds)
row.names(fc_df) <- NULL

# Pooled forward-chaining R²
fc_ss_res <- sum((fc_df$Observed - fc_df$Predicted)^2)
fc_ss_tot <- sum((fc_df$Observed - mean(fc_df$Observed))^2)
fc_r2_pooled <- 1 - fc_ss_res / fc_ss_tot
fc_rmse <- sqrt(mean((fc_df$Observed - fc_df$Predicted)^2))

cat("\nFwd-chain pooled:  R² =", round(fc_r2_pooled, 4),
    "  RMSE =", round(fc_rmse, 4),
    "  n =", nrow(fc_df), "\n\n")


# ============================================================================
# 5. SAVE PREDICTIONS
# ============================================================================
all_preds <- rbind(loyo_df, fc_df)
write.csv(all_preds, "2_incremental/temporal_val_predictions.csv",
          row.names = FALSE)
cat("Predictions saved: 2_incremental/temporal_val_predictions.csv\n\n")


# ============================================================================
# 6. DIAGNOSTIC PLOTS
# ============================================================================
pdf("4_products/diagnostics/temporal_validation.pdf",
    width = 10, height = 8)

# --- Panel 1: LOYO obs vs pred ---
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

plot(loyo_df$Predicted, loyo_df$Observed,
     pch = 16, cex = 0.7, col = adjustcolor("steelblue", 0.6),
     xlab = "Predicted", ylab = "Observed",
     main = paste0("LOYO Jackknife (R² = ",
                   round(loyo_r2_pooled, 3), ")"))
abline(0, 1, lty = 2, col = "red")

# --- Panel 2: Forward-chain obs vs pred ---
plot(fc_df$Predicted, fc_df$Observed,
     pch = 16, cex = 0.7, col = adjustcolor("darkgreen", 0.6),
     xlab = "Predicted", ylab = "Observed",
     main = paste0("Forward-Chaining (R² = ",
                   round(fc_r2_pooled, 3), ")"))
abline(0, 1, lty = 2, col = "red")

# --- Panel 3: LOYO R² by year ---
loyo_yr <- do.call(rbind, lapply(split(loyo_df, loyo_df$Year), function(d) {
  ss_r <- sum((d$Observed - d$Predicted)^2)
  ss_t <- sum((d$Observed - mean(d$Observed))^2)
  data.frame(Year = d$Year[1], R2 = ifelse(ss_t > 0, 1 - ss_r / ss_t, NA),
             n = nrow(d))
}))

plot(loyo_yr$Year, loyo_yr$R2, type = "b", pch = 16,
     xlab = "Year", ylab = "R²",
     main = "LOYO R² by Year", ylim = c(min(loyo_yr$R2, 0, na.rm = TRUE),
                                        max(loyo_yr$R2, 1, na.rm = TRUE)))
abline(h = 0, lty = 3, col = "grey50")
abline(h = loyo_r2_pooled, lty = 2, col = "steelblue")
text(max(loyo_yr$Year), loyo_r2_pooled, "pooled", pos = 2,
     col = "steelblue", cex = 0.8)

# --- Panel 4: Forward-chain R² by year ---
fc_yr <- do.call(rbind, lapply(split(fc_df, fc_df$Year), function(d) {
  ss_r <- sum((d$Observed - d$Predicted)^2)
  ss_t <- sum((d$Observed - mean(d$Observed))^2)
  data.frame(Year = d$Year[1], R2 = ifelse(ss_t > 0, 1 - ss_r / ss_t, NA),
             n = nrow(d))
}))

plot(fc_yr$Year, fc_yr$R2, type = "b", pch = 16,
     xlab = "Year", ylab = "R²",
     main = "Forward-Chain R² by Year",
     ylim = c(min(fc_yr$R2, 0, na.rm = TRUE),
              max(fc_yr$R2, 1, na.rm = TRUE)))
abline(h = 0, lty = 3, col = "grey50")
abline(h = fc_r2_pooled, lty = 2, col = "darkgreen")
text(max(fc_yr$Year), fc_r2_pooled, "pooled", pos = 2,
     col = "darkgreen", cex = 0.8)

# --- Page 2: Residual time series ---
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

loyo_df$Resid <- loyo_df$Observed - loyo_df$Predicted
plot(loyo_df$Year + runif(nrow(loyo_df), -0.2, 0.2),
     loyo_df$Resid,
     pch = 16, cex = 0.6, col = adjustcolor("steelblue", 0.5),
     xlab = "Year", ylab = "Residual",
     main = "LOYO Residuals by Year")
abline(h = 0, lty = 2)
loyo_yr_mean <- aggregate(Resid ~ Year, data = loyo_df, FUN = mean)
lines(loyo_yr_mean$Year, loyo_yr_mean$Resid, col = "red", lwd = 2)

fc_df$Resid <- fc_df$Observed - fc_df$Predicted
plot(fc_df$Year + runif(nrow(fc_df), -0.2, 0.2),
     fc_df$Resid,
     pch = 16, cex = 0.6, col = adjustcolor("darkgreen", 0.5),
     xlab = "Year", ylab = "Residual",
     main = "Forward-Chain Residuals by Year")
abline(h = 0, lty = 2)
fc_yr_mean <- aggregate(Resid ~ Year, data = fc_df, FUN = mean)
lines(fc_yr_mean$Year, fc_yr_mean$Resid, col = "red", lwd = 2)

dev.off()
cat("Plots saved: 4_products/diagnostics/temporal_validation.pdf\n\n")


# ============================================================================
# 7. SCORECARD
# ============================================================================
cat("====================================================\n")
cat("TEMPORAL VALIDATION SCORECARD\n")
cat("====================================================\n\n")

scorecard <- data.frame(
  Scheme          = c("Full model (in-sample)",
                      "LOYO Jackknife",
                      "Forward-Chaining"),
  n               = c(nrow(mdat),
                      nrow(loyo_df),
                      nrow(fc_df)),
  n_years         = c(length(unique(mdat$Year)),
                      length(unique(loyo_df$Year)),
                      length(unique(fc_df$Year))),
  R2_adj          = c(round(summary(m_full)$r.sq, 4),
                      NA, NA),
  R2_pooled       = c(NA,
                      round(loyo_r2_pooled, 4),
                      round(fc_r2_pooled, 4)),
  RMSE            = c(round(sqrt(mean(residuals(m_full)^2)), 4),
                      round(loyo_rmse, 4),
                      round(fc_rmse, 4)),
  stringsAsFactors = FALSE
)
print(scorecard, row.names = FALSE)

cat("\n--- Per-Year Detail ---\n\n")

cat("LOYO by year:\n")
loyo_yr$RMSE <- sapply(split(loyo_df, loyo_df$Year), function(d)
  round(sqrt(mean((d$Observed - d$Predicted)^2)), 3))
print(loyo_yr, row.names = FALSE)

cat("\nForward-chain by year:\n")
fc_yr$RMSE <- sapply(split(fc_df, fc_df$Year), function(d)
  round(sqrt(mean((d$Observed - d$Predicted)^2)), 3))
fc_yr$train_yrs <- sapply(split(fc_df, fc_df$Year), function(d) {
  y <- d$Year[1]
  length(unique(mdat$Year[mdat$Year < y]))
})
print(fc_yr, row.names = FALSE)

# Context: spatial validation for comparison
cat("\n--- Context ---\n")
cat("LOSO pooled R² (from core model script): 0.624\n")
cat("Full-data R² (adj):                     ",
    round(summary(m_full)$r.sq, 4), "\n")
cat("Expectation: LOYO > LOSO; FwdChain <= LOYO\n")
cat("If FwdChain ~ LOYO, temporal generalization is solid.\n")
cat("If FwdChain << LOYO, model may be overfitting to future structure.\n")

cat("\n--- Done. ---\n")