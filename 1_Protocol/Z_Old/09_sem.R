# ============================================================================
# 09_sem.R
# UCFR Filamentous Algae Project
# Stage 8: Full SEM chain assembly, evaluation, and path diagram
#
# Inputs:   3_models/brt_TP.rds
#           3_models/brt_TN.rds
#           3_models/brt_bloom_fitted.rds
#           2_incremental/ucfr_model_ready.csv
# Outputs:  2_incremental/sem_predictions.csv
#           4_products/diagnostics/sem_chain_fit.pdf
#           4_products/diagnostics/sem_path_diagram.pdf
#
# Steps:
#   1. Load model objects and data
#   2. Run full prediction chain: hydrology -> nutrients -> bloom
#   3. Overall and per-site chain performance
#   4. Path diagram with variable importance as path weights
#   5. Observed vs chain-predicted plots by site
#   6. Save predictions
#
# Notes:
#   - Path weights are variable importance (%) from each BRT submodel
#   - Standardized path coefficients not used — BRTs are non-parametric
#   - Chain constrained to observed predictor ranges per site
# ============================================================================

library(readr)
library(dismo)
library(gbm)

# ----------------------------------------------------------------------------
# 1. Load model objects and data
# ----------------------------------------------------------------------------

cat("Loading model objects...\n")
brt_TP    <- readRDS("3_models/brt_TP.rds")
brt_TN    <- readRDS("3_models/brt_TN.rds")
brt_bloom <- readRDS("3_models/brt_bloom_fitted.rds")

cat("Loading model-ready dataset...\n")
dat <- as.data.frame(
  read_csv("2_incremental/ucfr_model_ready.csv", show_col_types = FALSE)
)

dat$logCHLa <- log10(dat$CHLa)
dat$Site    <- as.factor(dat$Site)

# Complete cases for full chain
chain_vars <- c("Site", "Year", "anomaly", "Temp_oC",
                "Days_Since_Freshet", "logCHLa")
dat_chain  <- dat[complete.cases(dat[ , chain_vars]), ]

cat(sprintf("Rows available for chain prediction: %d\n\n", nrow(dat_chain)))

# ----------------------------------------------------------------------------
# 2. Run full prediction chain
# ----------------------------------------------------------------------------

cat("Running prediction chain...\n")

# Step 1a: predict log10(TP) from hydrology + Site
dat_chain$pred_logTP <- predict(
  brt_TP,
  newdata = dat_chain,
  n.trees = brt_TP$gbm.call$best.trees
)

# Step 1b: predict log10(TN) from hydrology + Site
dat_chain$pred_logTN <- predict(
  brt_TN,
  newdata = dat_chain,
  n.trees = brt_TN$gbm.call$best.trees
)

# Step 2: predict log10(CHLa) from predicted nutrients + hydrology + Site
dat_chain$pred_logCHLa <- predict(
  brt_bloom,
  newdata = dat_chain,
  n.trees = brt_bloom$gbm.call$best.trees
)

cat("Chain predictions complete.\n\n")

# ----------------------------------------------------------------------------
# 3. Chain performance
# ----------------------------------------------------------------------------

# Residuals
dat_chain$resid_logCHLa <- dat_chain$logCHLa - dat_chain$pred_logCHLa

# Overall
r_chain   <- cor(dat_chain$logCHLa, dat_chain$pred_logCHLa,
                 use = "complete.obs")
rmse_chain <- sqrt(mean(dat_chain$resid_logCHLa^2, na.rm = TRUE))

cat("--- Full Chain Performance ---\n")
cat(sprintf("  Overall r    = %.3f\n", r_chain))
cat(sprintf("  Overall R²   = %.3f\n", r_chain^2))
cat(sprintf("  Overall RMSE = %.3f log units\n\n", rmse_chain))

# Per-site
cat("--- Per-site Chain Performance ---\n")
cat(sprintf("  %-6s  %6s  %6s  %8s  %5s\n", "Site", "r", "R²", "RMSE", "n"))
cat(paste(rep("-", 42), collapse = ""), "\n")

site_stats <- data.frame(
  Site = character(), r = numeric(), R2 = numeric(),
  RMSE = numeric(), n = integer(), stringsAsFactors = FALSE
)

for (s in sort(levels(dat_chain$Site))) {
  d    <- dat_chain[dat_chain$Site == s, ]
  r    <- cor(d$logCHLa, d$pred_logCHLa, use = "complete.obs")
  rmse <- sqrt(mean(d$resid_logCHLa^2, na.rm = TRUE))
  cat(sprintf("  %-6s  %6.3f  %6.3f  %8.3f  %5d\n",
              s, r, r^2, rmse, nrow(d)))
  site_stats <- rbind(site_stats, data.frame(
    Site = s, r = r, R2 = r^2, RMSE = rmse, n = nrow(d)
  ))
}

# ----------------------------------------------------------------------------
# 4. Path diagram
# ----------------------------------------------------------------------------

cat("\nGenerating path diagram...\n")

# Extract variable importance from each submodel
imp_TP    <- brt_TP$contributions
imp_TN    <- brt_TN$contributions
imp_bloom <- brt_bloom$contributions

# Helper to get importance for a variable
get_imp <- function(contributions, varname) {
  val <- contributions$rel.inf[contributions$var == varname]
  if (length(val) == 0) return(0)
  round(val, 1)
}

# Hydrology -> TP importances
ht_anomaly <- get_imp(imp_TP, "anomaly")
ht_temp    <- get_imp(imp_TP, "Temp_oC")
ht_dsf     <- get_imp(imp_TP, "Days_Since_Freshet")
ht_site    <- get_imp(imp_TP, "Site")

# Hydrology -> TN importances
hn_anomaly <- get_imp(imp_TN, "anomaly")
hn_temp    <- get_imp(imp_TN, "Temp_oC")
hn_dsf     <- get_imp(imp_TN, "Days_Since_Freshet")
hn_site    <- get_imp(imp_TN, "Site")

# Nutrients + hydrology -> bloom importances
bc_tp      <- get_imp(imp_bloom, "pred_logTP")
bc_tn      <- get_imp(imp_bloom, "pred_logTN")
bc_anomaly <- get_imp(imp_bloom, "anomaly")
bc_temp    <- get_imp(imp_bloom, "Temp_oC")
bc_dsf     <- get_imp(imp_bloom, "Days_Since_Freshet")
bc_site    <- get_imp(imp_bloom, "Site")

pdf("4_products/diagnostics/sem_path_diagram.pdf",
    width = 12, height = 7)

par(mar = c(1, 1, 2, 1), bg = "white")
plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(0, 7),
     xaxt = "n", yaxt = "n", bty = "n",
     main = "BRT-SEM Path Diagram — UCFR Filamentous Algae",
     cex.main = 1.1)

# ── Node positions ──────────────────────────────────────────────────────────
# Hydrology nodes (left column)
nodes <- list(
  anomaly = c(1.2, 6.0),
  temp    = c(1.2, 4.5),
  dsf     = c(1.2, 3.0),
  site    = c(1.2, 1.5),
  tp      = c(4.5, 5.5),
  tn      = c(4.5, 3.0),
  bloom   = c(8.0, 4.2)
)

node_labels <- list(
  anomaly = "Anomaly",
  temp    = "Temperature",
  dsf     = "Days Since\nFreshet",
  site    = "Site",
  tp      = "log10(TP)",
  tn      = "log10(TN)",
  bloom   = "log10(CHLa)\nBloom Biomass"
)

node_cols <- list(
  anomaly = "#AED6F1", temp = "#AED6F1", dsf = "#AED6F1", site = "#D5DBDB",
  tp = "#A9DFBF", tn = "#A9DFBF", bloom = "#F9E79F"
)

# Draw nodes
box_w <- 1.3
box_h <- 0.55

draw_node <- function(name) {
  x   <- nodes[[name]][1]
  y   <- nodes[[name]][2]
  col <- node_cols[[name]]
  rect(x - box_w/2, y - box_h/2, x + box_w/2, y + box_h/2,
       col = col, border = "grey40", lwd = 1.5)
  text(x, y, node_labels[[name]], cex = 0.75, font = 2)
}

for (nm in names(nodes)) draw_node(nm)

# ── Arrow drawing helper ─────────────────────────────────────────────────────
draw_arrow <- function(from, to, label, col = "grey30", lwd = 1.5) {
  x0 <- nodes[[from]][1] + box_w/2
  y0 <- nodes[[from]][2]
  x1 <- nodes[[to]][1]   - box_w/2
  y1 <- nodes[[to]][2]
  # Offset slightly for overlapping arrows
  arrows(x0, y0, x1, y1, length = 0.1, angle = 20,
         col = col, lwd = lwd)
  mx <- (x0 + x1) / 2
  my <- (y0 + y1) / 2 + 0.18
  text(mx, my, paste0(label, "%"), cex = 0.62, col = col, font = 2)
}

# ── Hydrology -> TP arrows ───────────────────────────────────────────────────
draw_arrow("anomaly", "tp",  ht_anomaly, col = "#1A5276")
draw_arrow("temp",    "tp",  ht_temp,    col = "#1A5276")
draw_arrow("dsf",     "tp",  ht_dsf,     col = "#1A5276")
draw_arrow("site",    "tp",  ht_site,    col = "#7F8C8D")

# ── Hydrology -> TN arrows ───────────────────────────────────────────────────
draw_arrow("anomaly", "tn",  hn_anomaly, col = "#1A5276")
draw_arrow("temp",    "tn",  hn_temp,    col = "#1A5276")
draw_arrow("dsf",     "tn",  hn_dsf,     col = "#1A5276")
draw_arrow("site",    "tn",  hn_site,    col = "#7F8C8D")

# ── Nutrients -> bloom arrows ─────────────────────────────────────────────────
draw_arrow("tp",      "bloom", bc_tp,      col = "#1E8449")
draw_arrow("tn",      "bloom", bc_tn,      col = "#1E8449")

# ── Direct hydrology -> bloom arrows ─────────────────────────────────────────
draw_arrow("anomaly", "bloom", bc_anomaly, col = "#922B21", lwd = 1.2)
draw_arrow("dsf",     "bloom", bc_dsf,     col = "#922B21", lwd = 1.2)
draw_arrow("temp",    "bloom", bc_temp,    col = "#922B21", lwd = 1.2)
draw_arrow("site",    "bloom", bc_site,    col = "#7F8C8D", lwd = 1.2)

# ── Legend ───────────────────────────────────────────────────────────────────
legend(0.1, 0.9,
       legend = c("Hydrology → Nutrients", "Nutrients → Bloom",
                  "Hydrology → Bloom (direct)", "Site effect"),
       col    = c("#1A5276", "#1E8449", "#922B21", "#7F8C8D"),
       lwd    = 2, bty = "n", cex = 0.75)

# ── Performance annotation ───────────────────────────────────────────────────
text(8.0, 1.2,
     sprintf("Chain R² = %.3f", r_chain^2),
     cex = 0.85, font = 2, col = "grey20")

dev.off()
cat("Path diagram saved to 4_products/diagnostics/sem_path_diagram.pdf\n")

# ----------------------------------------------------------------------------
# 5. Observed vs chain-predicted plots by site
# ----------------------------------------------------------------------------

cat("Generating chain fit plots...\n")

site_cols <- c(DL = "#E41A1C", GR = "#377EB8", BN = "#4DAF4A",
               MS = "#984EA3", BM = "#FF7F00", HU = "#A65628",
               FH = "#F781BF")

pdf("4_products/diagnostics/sem_chain_fit.pdf",
    width = 12, height = 9)

par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))

# Overall plot
rng <- range(c(dat_chain$logCHLa, dat_chain$pred_logCHLa), na.rm = TRUE)
plot(dat_chain$logCHLa, dat_chain$pred_logCHLa,
     xlim = rng, ylim = rng,
     xlab = "Observed log10(CHLa)",
     ylab = "Chain-predicted log10(CHLa)",
     main = sprintf("All sites  |  r=%.3f  R²=%.3f", r_chain, r_chain^2),
     pch  = 16,
     col  = site_cols[as.character(dat_chain$Site)],
     cex  = 0.85)
abline(0, 1, col = "grey40", lty = 2)
legend("topleft", legend = names(site_cols), col = site_cols,
       pch = 16, cex = 0.6, bty = "n")

# Per-site plots
for (s in sort(levels(dat_chain$Site))) {
  d    <- dat_chain[dat_chain$Site == s, ]
  r    <- cor(d$logCHLa, d$pred_logCHLa, use = "complete.obs")
  rng_s <- range(c(d$logCHLa, d$pred_logCHLa), na.rm = TRUE)
  
  plot(d$logCHLa, d$pred_logCHLa,
       xlim = rng_s, ylim = rng_s,
       xlab = "Observed log10(CHLa)",
       ylab = "Predicted log10(CHLa)",
       main = sprintf("%s  |  r=%.3f  R²=%.3f  n=%d",
                      s, r, r^2, nrow(d)),
       pch  = 16,
       col  = site_cols[s],
       cex  = 0.9)
  abline(0, 1, col = "grey40", lty = 2)
  
  # Add year labels for context
  text(d$logCHLa, d$pred_logCHLa,
       labels = substr(d$Year, 3, 4),
       cex = 0.5, pos = 3, col = "grey50")
}

dev.off()
cat("Chain fit plots saved to 4_products/diagnostics/sem_chain_fit.pdf\n")

# ----------------------------------------------------------------------------
# 6. Save predictions
# ----------------------------------------------------------------------------

out_cols <- c("Site", "Year", "Month",
              "logCHLa", "pred_logTP", "pred_logTN",
              "pred_logCHLa", "resid_logCHLa",
              "anomaly", "Temp_oC", "Days_Since_Freshet")

write_csv(dat_chain[ , out_cols],
          "2_incremental/sem_predictions.csv")

cat("\nPredictions saved to 2_incremental/sem_predictions.csv\n")

# Save per-site performance table
write_csv(site_stats, "2_incremental/sem_site_performance.csv")
cat("Site performance saved to 2_incremental/sem_site_performance.csv\n")

cat("Done.\n")