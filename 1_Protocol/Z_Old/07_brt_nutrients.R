# ============================================================================
# 07_brt_nutrients.R
# UCFR Filamentous Algae Project
# Stage 6: BRT submodels for SPC, TP and TN — no Site predictor
#
# Input:    2_incremental/ucfr_model_ready.csv
# Outputs:  3_models/brt_SPC.rds
#           3_models/brt_TP.rds
#           3_models/brt_TN.rds
#           2_incremental/brt_nutrients_fitted.csv
#           4_products/diagnostics/brt_nutrients_pdp.pdf
#           4_products/diagnostics/brt_nutrients_fit.pdf
#
# Model structure:
#   SPC      ~ anomaly + Q_obs_cfs + Temp_oC + Days_Since_Freshet
#   log10(TP) ~ anomaly + Q_obs_cfs + Temp_oC + Days_Since_Freshet
#   log10(TN) ~ anomaly + Q_obs_cfs + Temp_oC + Days_Since_Freshet
#
# BRT settings: tc=3, lr=0.01, bag=0.75
#
# Notes:
#   - No Site predictor — fully transferable to new monitoring locations
#   - SPC is parallel to nutrients in causal chain
#   - Q_obs_cfs added as dilution signal for all three submodels
#   - Partial dependence plots for all four predictors
# ============================================================================

library(readr)
library(dismo)
library(gbm)

# ----------------------------------------------------------------------------
# 1. Read and prepare data
# ----------------------------------------------------------------------------

cat("Reading model-ready dataset...\n")
dat <- as.data.frame(read_csv("2_incremental/ucfr_model_ready.csv",
                              show_col_types = FALSE))

for (d in c("3_models", "2_incremental", "4_products/diagnostics")) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

dat$logTP <- log10(dat$TP_mg_L)
dat$logTN <- log10(dat$TN_mg_L)

predictors <- c("anomaly", "Q_obs_cfs", "Temp_oC", "Days_Since_Freshet")

vars_needed <- c("Site", "Year", "Month", predictors,
                 "SPC", "logTP", "logTN")
dat_model   <- as.data.frame(dat[complete.cases(dat[ , vars_needed]), ])

cat(sprintf("Complete cases: %d of %d rows\n\n", nrow(dat_model), nrow(dat)))

cat("Response variable summaries:\n")
cat(sprintf("  SPC        mean=%.1f   sd=%.1f   range=[%.1f, %.1f]\n",
            mean(dat_model$SPC), sd(dat_model$SPC),
            min(dat_model$SPC), max(dat_model$SPC)))
cat(sprintf("  log10(TP)  mean=%.3f  sd=%.3f  range=[%.3f, %.3f]\n",
            mean(dat_model$logTP), sd(dat_model$logTP),
            min(dat_model$logTP), max(dat_model$logTP)))
cat(sprintf("  log10(TN)  mean=%.3f  sd=%.3f  range=[%.3f, %.3f]\n\n",
            mean(dat_model$logTN), sd(dat_model$logTN),
            min(dat_model$logTN), max(dat_model$logTN)))

cat("Observed predictor ranges by site:\n")
cat(sprintf("  %-6s  %-18s  %-18s  %-18s  %-18s\n",
            "Site", "anomaly", "Q_obs_cfs", "Temp_oC", "Days_Since_Freshet"))
for (s in c("DL","GR","BN","MS","BM","HU","FH")) {
  d <- dat_model[dat_model$Site == s, ]
  if (nrow(d) == 0) next
  cat(sprintf("  %-6s  [%5.2f,%5.2f]    [%6.0f,%6.0f]  [%4.1f,%4.1f]      [%3.0f,%3.0f]\n",
              s,
              min(d$anomaly), max(d$anomaly),
              min(d$Q_obs_cfs), max(d$Q_obs_cfs),
              min(d$Temp_oC), max(d$Temp_oC),
              min(d$Days_Since_Freshet), max(d$Days_Since_Freshet)))
}
cat("\n")

# ----------------------------------------------------------------------------
# 2. BRT fitting function
# ----------------------------------------------------------------------------

fit_brt <- function(data, response, predictors, tc = 3, lr = 0.01, bag = 0.75) {
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

# ----------------------------------------------------------------------------
# 3. Fit submodels
# ----------------------------------------------------------------------------

cat("Fitting BRT submodel for SPC...\n")
brt_SPC <- fit_brt(dat_model, "SPC", predictors)

cat("Fitting BRT submodel for log10(TP)...\n")
brt_TP <- fit_brt(dat_model, "logTP", predictors)

cat("Fitting BRT submodel for log10(TN)...\n")
brt_TN <- fit_brt(dat_model, "logTN", predictors)

# ----------------------------------------------------------------------------
# 4. Training performance and variable importance
# ----------------------------------------------------------------------------

pred_SPC <- predict(brt_SPC, newdata = dat_model,
                    n.trees = brt_SPC$gbm.call$best.trees)
pred_TP  <- predict(brt_TP,  newdata = dat_model,
                    n.trees = brt_TP$gbm.call$best.trees)
pred_TN  <- predict(brt_TN,  newdata = dat_model,
                    n.trees = brt_TN$gbm.call$best.trees)

r_SPC <- cor(dat_model$SPC,   pred_SPC)
r_TP  <- cor(dat_model$logTP, pred_TP)
r_TN  <- cor(dat_model$logTN, pred_TN)

cat("\n--- Training Performance ---\n")
cat(sprintf("  SPC        r = %.3f  R² = %.3f  trees = %d\n",
            r_SPC, r_SPC^2, brt_SPC$gbm.call$best.trees))
cat(sprintf("  log10(TP)  r = %.3f  R² = %.3f  trees = %d\n",
            r_TP, r_TP^2, brt_TP$gbm.call$best.trees))
cat(sprintf("  log10(TN)  r = %.3f  R² = %.3f  trees = %d\n",
            r_TN, r_TN^2, brt_TN$gbm.call$best.trees))

cat("\n--- Variable Importance (%) ---\n")
for (nm in c("SPC", "log10(TP)", "log10(TN)")) {
  brt_obj <- list(SPC = brt_SPC, "log10(TP)" = brt_TP,
                  "log10(TN)" = brt_TN)[[nm]]
  cat(sprintf("  %s:\n", nm))
  imp <- brt_obj$contributions
  for (i in seq_len(nrow(imp))) {
    cat(sprintf("    %-25s  %.1f%%\n", imp$var[i], imp$rel.inf[i]))
  }
}

# ----------------------------------------------------------------------------
# 5. Per-site residual diagnostics
# ----------------------------------------------------------------------------

dat_model$pred_SPC  <- pred_SPC
dat_model$pred_logTP <- pred_TP
dat_model$pred_logTN <- pred_TN

cat("\n--- Per-site Residual Diagnostics ---\n")
cat(sprintf("  %-6s  %6s  %6s  %6s  %8s  %8s  %8s  %5s\n",
            "Site", "r_SPC", "r_TP", "r_TN",
            "RMSE_SPC", "RMSE_TP", "RMSE_TN", "n"))
cat(paste(rep("-", 72), collapse = ""), "\n")

for (s in c("DL","GR","BN","MS","BM","HU","FH")) {
  d <- dat_model[dat_model$Site == s, ]
  if (nrow(d) == 0) next
  r_spc  <- cor(d$SPC,   d$pred_SPC,   use = "complete.obs")
  r_tp   <- cor(d$logTP, d$pred_logTP, use = "complete.obs")
  r_tn   <- cor(d$logTN, d$pred_logTN, use = "complete.obs")
  rmse_spc <- sqrt(mean((d$SPC   - d$pred_SPC)^2,   na.rm = TRUE))
  rmse_tp  <- sqrt(mean((d$logTP - d$pred_logTP)^2, na.rm = TRUE))
  rmse_tn  <- sqrt(mean((d$logTN - d$pred_logTN)^2, na.rm = TRUE))
  cat(sprintf("  %-6s  %6.3f  %6.3f  %6.3f  %8.1f  %8.3f  %8.3f  %5d\n",
              s, r_spc, r_tp, r_tn,
              rmse_spc, rmse_tp, rmse_tn, nrow(d)))
}

# ----------------------------------------------------------------------------
# 6. Save fitted values
# ----------------------------------------------------------------------------

fitted_out <- dat_model[ , c("Site", "Year", "Month",
                             "SPC", "logTP", "logTN",
                             "pred_SPC", "pred_logTP", "pred_logTN")]

write_csv(fitted_out, "2_incremental/brt_nutrients_fitted.csv")
cat("\nFitted values saved to 2_incremental/brt_nutrients_fitted.csv\n")

# ----------------------------------------------------------------------------
# 7. Partial dependence plots
# ----------------------------------------------------------------------------

cat("\nGenerating partial dependence plots...\n")

pred_labels <- c(
  anomaly            = "Anomaly",
  Q_obs_cfs          = "Q obs (cfs)",
  Temp_oC            = "Temperature (°C)",
  Days_Since_Freshet = "Days Since Freshet"
)

pdf("4_products/diagnostics/brt_nutrients_pdp.pdf",
    width = 12, height = 10)

responses <- list(
  list(model = brt_SPC, label = "SPC"),
  list(model = brt_TP,  label = "log10(TP)"),
  list(model = brt_TN,  label = "log10(TN)")
)

for (resp in responses) {
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
  for (p in predictors) {
    if (p %in% resp$model$var.names) {
      gbm.plot(resp$model,
               variable.no = which(resp$model$var.names == p),
               smooth      = TRUE,
               rug         = TRUE,
               plot.layout = c(1, 1),
               write.title = FALSE,
               y.label     = paste(resp$label, "marginal effect"),
               x.label     = pred_labels[p])
      if (resp$label == "log10(TP)" && p == "Q_obs_cfs") {
        mtext("dilution effect", side = 3, line = 0.3, cex = 0.7, col = "red")
      }
    }
  }
  mtext(paste("Partial Dependence:", resp$label),
        outer = TRUE, cex = 1.1, font = 2)
}

dev.off()
cat("PDP saved to 4_products/diagnostics/brt_nutrients_pdp.pdf\n")

# ----------------------------------------------------------------------------
# 8. Observed vs fitted plots
# ----------------------------------------------------------------------------

cat("Generating observed vs fitted plots...\n")

site_cols <- c(DL = "#E41A1C", GR = "#377EB8", BN = "#4DAF4A",
               MS = "#984EA3", BM = "#FF7F00", HU = "#A65628",
               FH = "#F781BF")

pdf("4_products/diagnostics/brt_nutrients_fit.pdf",
    width = 12, height = 5)

par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

plot_fit <- function(obs, pred, label, site) {
  r   <- cor(obs, pred, use = "complete.obs")
  rng <- range(c(obs, pred), na.rm = TRUE)
  plot(obs, pred,
       xlim = rng, ylim = rng,
       xlab = paste("Observed", label),
       ylab = paste("Fitted", label),
       main = sprintf("%s  r=%.3f  R²=%.3f", label, r, r^2),
       pch  = 16, col = site_cols[as.character(site)], cex = 0.9)
  abline(0, 1, col = "grey40", lty = 2)
  legend("topleft", legend = names(site_cols), col = site_cols,
         pch = 16, cex = 0.6, bty = "n")
}

plot_fit(dat_model$SPC,   dat_model$pred_SPC,   "SPC",       dat_model$Site)
plot_fit(dat_model$logTP, dat_model$pred_logTP, "log10(TP)", dat_model$Site)
plot_fit(dat_model$logTN, dat_model$pred_logTN, "log10(TN)", dat_model$Site)

dev.off()
cat("Fit plots saved to 4_products/diagnostics/brt_nutrients_fit.pdf\n")

# ----------------------------------------------------------------------------
# 9. Save model objects
# ----------------------------------------------------------------------------

saveRDS(brt_SPC, "3_models/brt_SPC.rds")
saveRDS(brt_TP,  "3_models/brt_TP.rds")
saveRDS(brt_TN,  "3_models/brt_TN.rds")
cat("\nModel objects saved to 3_models/\n")
cat("Done.\n")