# ============================================================================
# 08_brt_bloom.R
# UCFR Filamentous Algae Project
# Stage 7: BRT bloom submodel with observed SPC, predicted TP and TN
#
# Inputs:   2_incremental/ucfr_model_ready.csv
#           2_incremental/brt_nutrients_fitted.csv
# Outputs:  3_models/brt_bloom_obs.rds      (observed nutrients + SPC)
#           3_models/brt_bloom_fitted.rds   (predicted nutrients + observed SPC)
#           3_models/brt_bloom_simple.rds   (simplified model)
#           2_incremental/brt_bloom_fitted.csv
#           4_products/diagnostics/brt_bloom_pdp.pdf
#           4_products/diagnostics/brt_bloom_fit.pdf
#           4_products/diagnostics/brt_bloom_simplify.pdf
#
# Model structure:
#   Version 1 — all observed (baseline):
#   log10(CHLa) ~ SPC + logTP_obs + logTN_obs + anomaly + Q_obs_cfs +
#                 Temp_oC + Days_Since_Freshet
#
#   Version 2 — predicted nutrients, observed SPC (chain model):
#   log10(CHLa) ~ SPC + pred_logTP + pred_logTN + anomaly + Q_obs_cfs +
#                 Temp_oC + Days_Since_Freshet
#
#   Version 3 — simplified via gbm.simplify
#
# BRT settings: tc=4, lr=0.01, bag=0.75
#
# Notes:
#   - SPC used as observed covariate — hydrology -> SPC chain weak at lower sites
#   - gbm.simplify used for data-driven variable selection
#   - No Site predictor anywhere — fully transferable model
# ============================================================================

library(readr)
library(dismo)
library(gbm)

# ----------------------------------------------------------------------------
# 1. Read and prepare data
# ----------------------------------------------------------------------------

cat("Reading data...\n")
dat <- as.data.frame(read_csv("2_incremental/ucfr_model_ready.csv",
                              show_col_types = FALSE))
nut <- as.data.frame(read_csv("2_incremental/brt_nutrients_fitted.csv",
                              show_col_types = FALSE))

dat <- merge(dat, nut[ , c("Site", "Year", "Month",
                           "pred_SPC", "pred_logTP", "pred_logTN")],
             by = c("Site", "Year", "Month"), all.x = TRUE)

for (d in c("3_models", "2_incremental", "4_products/diagnostics")) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

dat$logCHLa   <- log10(dat$CHLa)
dat$logTP_obs <- log10(dat$TP_mg_L)
dat$logTN_obs <- log10(dat$TN_mg_L)

# ----------------------------------------------------------------------------
# 2. Define predictor sets
# ----------------------------------------------------------------------------

hydro      <- c("anomaly", "Q_obs_cfs", "Temp_oC", "Days_Since_Freshet")

# Version 1 — all observed (ceiling / benchmark)
preds_v1   <- c("SPC", "logTP_obs", "logTN_obs", hydro)

# Version 2 — full chain: predicted SPC, TP, TN from hydrology
# This is the climate-ready model — no observed chemistry required
preds_v2   <- c("pred_SPC", "pred_logTP", "pred_logTN", hydro)

vars_v1    <- c("Site", "Year", "Month", preds_v1, "logCHLa")
vars_v2    <- c("Site", "Year", "Month", preds_v2, "logCHLa")

dat_v1     <- as.data.frame(dat[complete.cases(dat[ , vars_v1]), ])
dat_v2     <- as.data.frame(dat[complete.cases(dat[ , vars_v2]), ])

cat(sprintf("Complete cases V1 (all observed):    %d rows\n", nrow(dat_v1)))
cat(sprintf("Complete cases V2 (full chain):      %d rows\n\n", nrow(dat_v2)))

cat(sprintf("log10(CHLa)  mean=%.3f  sd=%.3f  range=[%.3f, %.3f]\n\n",
            mean(dat_v1$logCHLa), sd(dat_v1$logCHLa),
            min(dat_v1$logCHLa), max(dat_v1$logCHLa)))

# ----------------------------------------------------------------------------
# 3. BRT fitting function
# ----------------------------------------------------------------------------

fit_brt <- function(data, response, predictors, tc = 4, lr = 0.01, bag = 0.75) {
  data_sub <- as.data.frame(data[ , c(predictors, response), drop = FALSE])
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
  )
}

print_diagnostics <- function(model, data, response, predictors, label) {
  pred <- predict(model, newdata = data, n.trees = model$gbm.call$best.trees)
  r    <- cor(data[[response]], pred, use = "complete.obs")
  
  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("  r = %.3f  R² = %.3f  trees = %d\n\n",
              r, r^2, model$gbm.call$best.trees))
  
  cat("  Variable Importance:\n")
  imp <- model$contributions
  for (i in seq_len(nrow(imp))) {
    cat(sprintf("    %-25s  %.1f%%\n", imp$var[i], imp$rel.inf[i]))
  }
  
  cat("\n  Per-site residuals:\n")
  cat(sprintf("    %-6s  %6s  %8s  %5s\n", "Site", "r", "RMSE", "n"))
  cat("  ", paste(rep("-", 32), collapse = ""), "\n")
  data$pred_ <- pred
  for (s in c("DL","GR","BN","MS","BM","HU","FH")) {
    d <- data[data$Site == s, ]
    if (nrow(d) == 0) next
    rs   <- cor(d[[response]], d$pred_, use = "complete.obs")
    rmse <- sqrt(mean((d[[response]] - d$pred_)^2, na.rm = TRUE))
    cat(sprintf("    %-6s  %6.3f  %8.3f  %5d\n", s, rs, rmse, nrow(d)))
  }
  invisible(pred)
}

# ----------------------------------------------------------------------------
# 4. Fit Version 1 — all observed
# ----------------------------------------------------------------------------

cat("Fitting Version 1 (all observed)...\n")
brt_bloom_obs <- fit_brt(dat_v1, "logCHLa", preds_v1)
pred_v1 <- print_diagnostics(brt_bloom_obs, dat_v1,
                             "logCHLa", preds_v1,
                             "Version 1: All Observed")

# ----------------------------------------------------------------------------
# 5. Fit Version 2 — predicted nutrients, observed SPC
# ----------------------------------------------------------------------------

cat("\nFitting Version 2 (full chain: predicted SPC + nutrients)...\n")
brt_bloom_fit <- fit_brt(dat_v2, "logCHLa", preds_v2)
pred_v2 <- print_diagnostics(brt_bloom_fit, dat_v2,
                             "logCHLa", preds_v2,
                             "Version 2: Full Chain (predicted SPC + nutrients)")

# ----------------------------------------------------------------------------
# 6. Chain performance comparison
# ----------------------------------------------------------------------------

r_v1 <- cor(dat_v1$logCHLa, pred_v1)
r_v2 <- cor(dat_v2$logCHLa, pred_v2)

cat("\n--- Chain Performance Comparison ---\n")
cat(sprintf("  Version 1 (all observed):        r = %.3f  R² = %.3f\n",
            r_v1, r_v1^2))
cat(sprintf("  Version 2 (full chain):          r = %.3f  R² = %.3f\n",
            r_v2, r_v2^2))
cat(sprintf("  R² cost of full chain:           %.3f\n", r_v1^2 - r_v2^2))

# ----------------------------------------------------------------------------
# 7. gbm.simplify on Version 2
# ----------------------------------------------------------------------------

cat("\nRunning gbm.simplify on Version 2...\n")

data_v2_sub <- as.data.frame(dat_v2[ , c(preds_v2, "logCHLa"), drop = FALSE])

simp <- gbm.simplify(
  brt_bloom_fit,
  n.drops   = length(preds_v2) - 1,
  plot      = FALSE
)

cat("\nSimplification results:\n")
cat(sprintf("  %-30s  %s\n", "Variables dropped", "CV deviance"))
cat(paste(rep("-", 50), collapse = ""), "\n")
for (i in seq_along(simp$deviance.summary$mean)) {
  cat(sprintf("  %-30s  %.4f (±%.4f)\n",
              paste(simp$pred.list[[i]], collapse = ", "),
              simp$deviance.summary$mean[i],
              simp$deviance.summary$se[i]))
}

# Identify optimal number of predictors
best_n_drop <- which.min(simp$deviance.summary$mean)
vars_to_drop <- simp$pred.list[[best_n_drop]]
preds_simple <- preds_v2[!preds_v2 %in% vars_to_drop]

cat(sprintf("\nOptimal: drop %d variable(s): %s\n",
            length(vars_to_drop),
            paste(vars_to_drop, collapse = ", ")))
cat(sprintf("Simplified predictor set: %s\n\n",
            paste(preds_simple, collapse = ", ")))

# Refit simplified model
cat("Fitting simplified model...\n")
brt_bloom_simple <- fit_brt(dat_v2, "logCHLa", preds_simple)
pred_simple <- print_diagnostics(brt_bloom_simple, dat_v2,
                                 "logCHLa", preds_simple,
                                 "Version 3: Simplified")

# Save simplify plot
pdf("4_products/diagnostics/brt_bloom_simplify.pdf", width = 8, height = 5)
plot(simp$deviance.summary$mean,
     type = "b", pch = 16,
     xlab = "Number of variables dropped",
     ylab = "CV deviance",
     main = "gbm.simplify — bloom model")
abline(v = best_n_drop, col = "red", lty = 2)
dev.off()

# ----------------------------------------------------------------------------
# 8. Save fitted values
# ----------------------------------------------------------------------------

dat_v2$pred_logCHLa  <- pred_v2
dat_v2$resid_logCHLa <- dat_v2$logCHLa - pred_v2

fitted_out <- dat_v2[ , c("Site", "Year", "Month", "logCHLa",
                          "pred_SPC", "pred_logTP", "pred_logTN",
                          "pred_logCHLa", "resid_logCHLa")]

write_csv(fitted_out, "2_incremental/brt_bloom_fitted.csv")
cat("\nFitted values saved to 2_incremental/brt_bloom_fitted.csv\n")

# ----------------------------------------------------------------------------
# 9. Partial dependence plots — Version 2
# ----------------------------------------------------------------------------

cat("\nGenerating partial dependence plots...\n")

pred_labels <- c(
  pred_SPC           = "predicted SPC (umho/cm)",
  pred_logTP         = "predicted log10(TP)",
  pred_logTN         = "predicted log10(TN)",
  anomaly            = "Anomaly",
  Q_obs_cfs          = "Q obs (cfs)",
  Temp_oC            = "Temperature (°C)",
  Days_Since_Freshet = "Days Since Freshet"
)

pdf("4_products/diagnostics/brt_bloom_pdp.pdf", width = 12, height = 8)
par(mfrow = c(2, 4), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))

for (p in preds_v2) {
  if (p %in% brt_bloom_fit$var.names) {
    gbm.plot(brt_bloom_fit,
             variable.no = which(brt_bloom_fit$var.names == p),
             smooth      = TRUE,
             rug         = TRUE,
             plot.layout = c(1, 1),
             write.title = FALSE,
             y.label     = "log10(CHLa) marginal effect",
             x.label     = pred_labels[p])
    if (p == "pred_logTP") {
      abline(v = log10(0.024), col = "red", lty = 3, lwd = 1.5)
    }
  }
}
mtext("Partial Dependence: log10(CHLa) — Version 2",
      outer = TRUE, cex = 1.1, font = 2)
dev.off()
cat("PDP saved to 4_products/diagnostics/brt_bloom_pdp.pdf\n")

# ----------------------------------------------------------------------------
# 10. Observed vs fitted plots
# ----------------------------------------------------------------------------

cat("Generating observed vs fitted plots...\n")

site_cols <- c(DL = "#E41A1C", GR = "#377EB8", BN = "#4DAF4A",
               MS = "#984EA3", BM = "#FF7F00", HU = "#A65628",
               FH = "#F781BF")

pdf("4_products/diagnostics/brt_bloom_fit.pdf", width = 12, height = 5)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

for (info in list(
  list(obs = dat_v1$logCHLa, pred = pred_v1,   site = dat_v1$Site, label = "V1: All Observed"),
  list(obs = dat_v2$logCHLa, pred = pred_v2,   site = dat_v2$Site, label = "V2: Full Chain"),
  list(obs = dat_v2$logCHLa, pred = pred_simple, site = dat_v2$Site, label = "V3: Simplified")
)) {
  r   <- cor(info$obs, info$pred, use = "complete.obs")
  rng <- range(c(info$obs, info$pred), na.rm = TRUE)
  plot(info$obs, info$pred,
       xlim = rng, ylim = rng,
       xlab = "Observed log10(CHLa)",
       ylab = "Fitted log10(CHLa)",
       main = sprintf("%s\nr=%.3f  R²=%.3f", info$label, r, r^2),
       pch = 16, col = site_cols[as.character(info$site)], cex = 0.9)
  abline(0, 1, col = "grey40", lty = 2)
  legend("topleft", legend = names(site_cols), col = site_cols,
         pch = 16, cex = 0.6, bty = "n")
}
dev.off()
cat("Fit plots saved to 4_products/diagnostics/brt_bloom_fit.pdf\n")

# ----------------------------------------------------------------------------
# 11. Save model objects
# ----------------------------------------------------------------------------

saveRDS(brt_bloom_obs,    "3_models/brt_bloom_obs.rds")
saveRDS(brt_bloom_fit,    "3_models/brt_bloom_fitted.rds")
saveRDS(brt_bloom_simple, "3_models/brt_bloom_simple.rds")
cat("\nModel objects saved to 3_models/\n")
cat("Done.\n")