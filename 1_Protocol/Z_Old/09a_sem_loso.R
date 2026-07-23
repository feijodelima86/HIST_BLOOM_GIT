# ============================================================================
# 09_1_sem_loso.R
# UCFR Filamentous Algae Project
# Stage 8.1: LOSO validation of full BRT-SEM chain
#
# Input:    2_incremental/ucfr_model_ready.csv
# Outputs:  2_incremental/sem_loso_predictions.csv
#           2_incremental/sem_loso_performance.csv
#           4_products/diagnostics/sem_loso_fit.pdf
#
# Full chain (no observed chemistry anywhere):
#   hydrology -> SPC
#   hydrology -> log10(TP)
#   hydrology -> log10(TN)
#   pred_SPC + pred_logTP + pred_logTN + hydrology -> log10(CHLa)
#
# LOSO procedure:
#   For each left-out site:
#     1. Refit all four submodels on 6 training sites
#     2. Predict left-out site through full chain
#     3. No site identity used anywhere
#
# BRT settings: tc=3 for SPC/TP/TN submodels, tc=4 for bloom
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

# Hydrological predictors — same for all three intermediate submodels
hydro_preds <- c("anomaly", "Q_obs_cfs", "Temp_oC", "Days_Since_Freshet")

# Bloom predictors — predicted intermediates + hydrology
bloom_preds <- c("pred_SPC", "pred_logTP", "pred_logTN", hydro_preds)

# Complete cases needed
chain_vars <- c("Site", "Year", "Month",
                hydro_preds, "SPC", "logTP_obs", "logTN_obs", "logCHLa")

dat_model <- as.data.frame(dat[complete.cases(dat[ , chain_vars]), ])
cat(sprintf("Complete cases: %d rows\n\n", nrow(dat_model)))

# ----------------------------------------------------------------------------
# 2. BRT fitting function — subsets data to avoid column leakage
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
# 3. LOSO loop
# ----------------------------------------------------------------------------

sites    <- sort(unique(dat_model$Site))
loso_out <- vector("list", length(sites))

cat("Running LOSO chain validation...\n")
cat("(Refitting 4 submodels per iteration — this will take several minutes)\n\n")

for (i in seq_along(sites)) {
  s     <- sites[i]
  train <- as.data.frame(dat_model[dat_model$Site != s, ])
  test  <- as.data.frame(dat_model[dat_model$Site == s, ])
  
  cat(sprintf("--- Leaving out %s (train n=%d, test n=%d) ---\n",
              s, nrow(train), nrow(test)))
  
  # -- Step 1: fit intermediate submodels on training sites ----------------
  cat("  Fitting brt_SPC...\n")
  brt_SPC_cv <- fit_brt(train, "SPC", hydro_preds, tc = 3)
  
  cat("  Fitting brt_TP...\n")
  brt_TP_cv  <- fit_brt(train, "logTP_obs", hydro_preds, tc = 3)
  
  cat("  Fitting brt_TN...\n")
  brt_TN_cv  <- fit_brt(train, "logTN_obs", hydro_preds, tc = 3)
  
  if (any(sapply(list(brt_SPC_cv, brt_TP_cv, brt_TN_cv), is.null))) {
    cat(sprintf("  SKIPPING %s — intermediate model fit failed\n\n", s))
    next
  }
  
  # -- Step 2: predict intermediates for left-out site ---------------------
  test$pred_SPC    <- predict(brt_SPC_cv, newdata = test[ , hydro_preds],
                              n.trees = brt_SPC_cv$gbm.call$best.trees)
  test$pred_logTP  <- predict(brt_TP_cv,  newdata = test[ , hydro_preds],
                              n.trees = brt_TP_cv$gbm.call$best.trees)
  test$pred_logTN  <- predict(brt_TN_cv,  newdata = test[ , hydro_preds],
                              n.trees = brt_TN_cv$gbm.call$best.trees)
  
  cat(sprintf("  Intermediates: SPC [%.0f,%.0f]  TP [%.3f,%.3f]  TN [%.3f,%.3f]\n",
              min(test$pred_SPC), max(test$pred_SPC),
              min(test$pred_logTP), max(test$pred_logTP),
              min(test$pred_logTN), max(test$pred_logTN)))
  
  # -- Step 3: predict intermediates for training sites -------------------
  train$pred_SPC   <- predict(brt_SPC_cv, newdata = train[ , hydro_preds],
                              n.trees = brt_SPC_cv$gbm.call$best.trees)
  train$pred_logTP <- predict(brt_TP_cv,  newdata = train[ , hydro_preds],
                              n.trees = brt_TP_cv$gbm.call$best.trees)
  train$pred_logTN <- predict(brt_TN_cv,  newdata = train[ , hydro_preds],
                              n.trees = brt_TN_cv$gbm.call$best.trees)
  
  # -- Step 4: fit bloom submodel on training sites -----------------------
  cat("  Fitting brt_bloom...\n")
  brt_bloom_cv <- fit_brt(train, "logCHLa", bloom_preds, tc = 4)
  
  if (is.null(brt_bloom_cv)) {
    cat(sprintf("  SKIPPING %s — bloom model fit failed\n\n", s))
    next
  }
  
  # -- Step 5: predict bloom for left-out site ----------------------------
  test$pred_logCHLa <- predict(brt_bloom_cv,
                               newdata = test[ , bloom_preds],
                               n.trees = brt_bloom_cv$gbm.call$best.trees)
  
  # -- Step 6: evaluate ---------------------------------------------------
  r_site <- cor(test$logCHLa, test$pred_logCHLa, use = "complete.obs")
  rmse   <- sqrt(mean((test$logCHLa - test$pred_logCHLa)^2, na.rm = TRUE))
  
  cat(sprintf("  Result: r = %.3f  RMSE = %.3f  n = %d\n\n",
              r_site, rmse, nrow(test)))
  
  loso_out[[i]] <- data.frame(
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

loso_all <- do.call(rbind, loso_out)

r_overall    <- cor(loso_all$Observed, loso_all$Predicted, use = "complete.obs")
rmse_overall <- sqrt(mean((loso_all$Observed - loso_all$Predicted)^2,
                          na.rm = TRUE))

cat("══════════════════════════════════════════\n")
cat("LOSO Full Chain Performance Summary\n")
cat("══════════════════════════════════════════\n")
cat(sprintf("  Overall r    = %.3f\n", r_overall))
cat(sprintf("  Overall R²   = %.3f\n", r_overall^2))
cat(sprintf("  Overall RMSE = %.3f log units\n\n", rmse_overall))

cat("Per-site LOSO Performance:\n")
cat(sprintf("  %-6s  %6s  %6s  %8s  %5s\n", "Site", "r", "R²", "RMSE", "n"))
cat(paste(rep("-", 42), collapse = ""), "\n")

perf_list <- vector("list", length(sites))

for (i in seq_along(sites)) {
  s  <- sites[i]
  d  <- loso_all[loso_all$Site == s, ]
  if (nrow(d) < 3) next
  r    <- cor(d$Observed, d$Predicted, use = "complete.obs")
  rmse <- sqrt(mean((d$Observed - d$Predicted)^2, na.rm = TRUE))
  cat(sprintf("  %-6s  %6.3f  %6.3f  %8.3f  %5d\n",
              s, r, r^2, rmse, nrow(d)))
  perf_list[[i]] <- data.frame(Site = s, r = r, R2 = r^2,
                               RMSE = rmse, n = nrow(d))
}

perf_df <- do.call(rbind, perf_list)

# Context vs training performance
cat("\n--- Context ---\n")
cat("  Training R² (V1 all observed):  0.830\n")
cat("  Training R² (V2 full chain):    0.709\n")
cat(sprintf("  LOSO R²   (full chain):         %.3f\n", r_overall^2))

# ----------------------------------------------------------------------------
# 5. LOSO plots
# ----------------------------------------------------------------------------

site_cols <- c(DL = "#E41A1C", GR = "#377EB8", BN = "#4DAF4A",
               MS = "#984EA3", BM = "#FF7F00", HU = "#A65628",
               FH = "#F781BF")

# ----------------------------------------------------------------------------
# 5. Path diagram with variable importance from full trained models
# ----------------------------------------------------------------------------

cat("\nGenerating path diagram...\n")

brt_SPC   <- readRDS("3_models/brt_SPC.rds")
brt_TP    <- readRDS("3_models/brt_TP.rds")
brt_TN    <- readRDS("3_models/brt_TN.rds")
brt_bloom <- readRDS("3_models/brt_bloom_fitted.rds")

get_imp <- function(contributions, varname) {
  val <- contributions$rel.inf[contributions$var == varname]
  if (length(val) == 0) return(0)
  round(val, 1)
}

imp_SPC   <- brt_SPC$contributions
imp_TP    <- brt_TP$contributions
imp_TN    <- brt_TN$contributions
imp_bloom <- brt_bloom$contributions

pdf("4_products/diagnostics/sem_path_diagram.pdf", width = 13, height = 7)
par(mar = c(1, 1, 2, 1), bg = "white")
plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(0, 8),
     xaxt = "n", yaxt = "n", bty = "n",
     main = "BRT-SEM Path Diagram — UCFR Filamentous Algae",
     cex.main = 1.1)

nodes <- list(
  anomaly = c(1.2, 7.0),
  qobs    = c(1.2, 5.5),
  temp    = c(1.2, 4.0),
  dsf     = c(1.2, 2.5),
  spc     = c(4.5, 6.5),
  tp      = c(4.5, 4.5),
  tn      = c(4.5, 2.5),
  bloom   = c(8.2, 4.5)
)

node_labels <- list(
  anomaly = "Anomaly",
  qobs    = "Q obs (cfs)",
  temp    = "Temperature",
  dsf     = "Days Since\nFreshet",
  spc     = "pred SPC",
  tp      = "pred log10(TP)",
  tn      = "pred log10(TN)",
  bloom   = "log10(CHLa)\nBloom Biomass"
)

node_cols <- list(
  anomaly = "#AED6F1", qobs = "#AED6F1",
  temp    = "#AED6F1", dsf  = "#AED6F1",
  spc     = "#A9DFBF", tp   = "#A9DFBF",
  tn      = "#A9DFBF", bloom = "#F9E79F"
)

box_w <- 1.4
box_h = 0.55

draw_node <- function(name) {
  x <- nodes[[name]][1]; y <- nodes[[name]][2]
  rect(x - box_w/2, y - box_h/2, x + box_w/2, y + box_h/2,
       col = node_cols[[name]], border = "grey40", lwd = 1.5)
  text(x, y, node_labels[[name]], cex = 0.72, font = 2)
}
for (nm in names(nodes)) draw_node(nm)

draw_arrow <- function(from, to, label, col = "grey30", lwd = 1.5, offset = 0) {
  x0 <- nodes[[from]][1] + box_w/2
  y0 <- nodes[[from]][2] + offset
  x1 <- nodes[[to]][1]   - box_w/2
  y1 <- nodes[[to]][2]   + offset
  arrows(x0, y0, x1, y1, length = 0.1, angle = 20, col = col, lwd = lwd)
  mx <- (x0 + x1) / 2; my <- (y0 + y1) / 2 + 0.18
  text(mx, my, paste0(label, "%"), cex = 0.60, col = col, font = 2)
}

# Hydrology -> SPC
draw_arrow("anomaly", "spc", get_imp(imp_SPC, "anomaly"),   col="#1A5276", offset= 0.15)
draw_arrow("qobs",    "spc", get_imp(imp_SPC, "Q_obs_cfs"), col="#1A5276", offset= 0.05)
draw_arrow("temp",    "spc", get_imp(imp_SPC, "Temp_oC"),   col="#1A5276", offset=-0.05)
draw_arrow("dsf",     "spc", get_imp(imp_SPC, "Days_Since_Freshet"), col="#1A5276", offset=-0.15)

# Hydrology -> TP
draw_arrow("anomaly", "tp", get_imp(imp_TP, "anomaly"),   col="#1A5276", offset= 0.15)
draw_arrow("qobs",    "tp", get_imp(imp_TP, "Q_obs_cfs"), col="#1A5276", offset= 0.05)
draw_arrow("temp",    "tp", get_imp(imp_TP, "Temp_oC"),   col="#1A5276", offset=-0.05)
draw_arrow("dsf",     "tp", get_imp(imp_TP, "Days_Since_Freshet"), col="#1A5276", offset=-0.15)

# Hydrology -> TN
draw_arrow("anomaly", "tn", get_imp(imp_TN, "anomaly"),   col="#1A5276", offset= 0.15)
draw_arrow("qobs",    "tn", get_imp(imp_TN, "Q_obs_cfs"), col="#1A5276", offset= 0.05)
draw_arrow("temp",    "tn", get_imp(imp_TN, "Temp_oC"),   col="#1A5276", offset=-0.05)
draw_arrow("dsf",     "tn", get_imp(imp_TN, "Days_Since_Freshet"), col="#1A5276", offset=-0.15)

# Intermediates -> bloom
draw_arrow("spc", "bloom", get_imp(imp_bloom, "pred_SPC"),    col="#1E8449", lwd=2)
draw_arrow("tp",  "bloom", get_imp(imp_bloom, "pred_logTP"),  col="#1E8449", lwd=2)
draw_arrow("tn",  "bloom", get_imp(imp_bloom, "pred_logTN"),  col="#1E8449", lwd=2)

# Direct hydrology -> bloom
draw_arrow("anomaly", "bloom", get_imp(imp_bloom, "anomaly"),           col="#922B21", lwd=1.2, offset= 0.2)
draw_arrow("qobs",    "bloom", get_imp(imp_bloom, "Q_obs_cfs"),         col="#922B21", lwd=1.2, offset= 0.1)
draw_arrow("temp",    "bloom", get_imp(imp_bloom, "Temp_oC"),           col="#922B21", lwd=1.2, offset=-0.1)
draw_arrow("dsf",     "bloom", get_imp(imp_bloom, "Days_Since_Freshet"),col="#922B21", lwd=1.2, offset=-0.2)

# Performance annotations
text(8.2, 1.8, sprintf("Training R² (V1): 0.830"), cex = 0.75, col = "grey20")
text(8.2, 1.4, sprintf("Training R² (V2): 0.709"), cex = 0.75, col = "grey20")
text(8.2, 1.0, sprintf("LOSO R²:  %.3f", r_overall^2), cex = 0.75, col = "grey20", font = 2)

legend(0.1, 1.2,
       legend = c("Hydrology → Intermediates",
                  "Intermediates → Bloom",
                  "Hydrology → Bloom (direct)"),
       col    = c("#1A5276", "#1E8449", "#922B21"),
       lwd    = 2, bty = "n", cex = 0.72)

dev.off()
cat("Path diagram saved to 4_products/diagnostics/sem_path_diagram.pdf\n\n")

# ----------------------------------------------------------------------------
# 6. LOSO fit plots
# ----------------------------------------------------------------------------

pdf("4_products/diagnostics/sem_loso_fit.pdf", width = 12, height = 9)
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))

rng <- range(c(loso_all$Observed, loso_all$Predicted), na.rm = TRUE)
plot(loso_all$Observed, loso_all$Predicted,
     xlim = rng, ylim = rng,
     xlab = "Observed log10(CHLa)",
     ylab = "LOSO Predicted log10(CHLa)",
     main = sprintf("All sites  |  r=%.3f  R²=%.3f",
                    r_overall, r_overall^2),
     pch = 16, col = site_cols[loso_all$Site], cex = 0.85)
abline(0, 1, col = "grey40", lty = 2)
legend("topleft", legend = names(site_cols), col = site_cols,
       pch = 16, cex = 0.6, bty = "n")

for (s in sites) {
  d <- loso_all[loso_all$Site == s, ]
  if (nrow(d) < 3) next
  r     <- cor(d$Observed, d$Predicted, use = "complete.obs")
  rng_s <- range(c(d$Observed, d$Predicted), na.rm = TRUE)
  plot(d$Observed, d$Predicted,
       xlim = rng_s, ylim = rng_s,
       xlab = "Observed log10(CHLa)",
       ylab = "LOSO Predicted log10(CHLa)",
       main = sprintf("%s  |  r=%.3f  R²=%.3f  n=%d",
                      s, r, r^2, nrow(d)),
       pch = 16, col = site_cols[s], cex = 0.9)
  abline(0, 1, col = "grey40", lty = 2)
  text(d$Observed, d$Predicted,
       labels = substr(d$Year, 3, 4),
       cex = 0.5, pos = 3, col = "grey50")
}

dev.off()
cat("\nLOSO fit plots saved to 4_products/diagnostics/sem_loso_fit.pdf\n")

# ----------------------------------------------------------------------------
# 6. Save outputs
# ----------------------------------------------------------------------------

write_csv(loso_all, "2_incremental/sem_loso_predictions.csv")
write_csv(perf_df,  "2_incremental/sem_loso_performance.csv")

cat("Saved:\n")
cat("  2_incremental/sem_loso_predictions.csv\n")
cat("  2_incremental/sem_loso_performance.csv\n")
cat("Done.\n")