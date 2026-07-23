# ============================================================================
# 10_climate_scenarios.R
# UCFR Filamentous Algae Project
# Stage 9: Climate scenario predictions
#
# Inputs:   3_models/brt_SPC.rds
#           3_models/brt_TP.rds
#           3_models/brt_TN.rds
#           3_models/brt_bloom_fitted.rds
#           2_incremental/ucfr_model_ready.csv
# Outputs:  2_incremental/scenario_predictions.csv
#           4_products/diagnostics/scenario_results.pdf
#
# Scenario matrix:
#   4 pathways x 3 time slices = 12 scenarios per site
#   HS: High precip, Snow  — smaller freshet, less flow
#   LS: Low precip,  Snow  — smaller freshet, much less flow
#   HR: High precip, Rain  — earlier/flashier freshet, more flow
#   LR: Low precip,  Rain  — earlier/flashier freshet, more flow
#
# Baseline: 2020 site means for hydrology, interpolated temperature
#
# Notes:
#   - Temperature linearly interpolated from NorthWEST projections
#   - DSF shifts from Larson snow scenarios; rain scenarios use same
#     shifts as corresponding snow scenario + uncertainty flag
#   - Q changes applied equally to peak Q (anomaly) and Q_obs_cfs
#   - Predictions constrained to observed training range per site
#   - All predictions on log10 scale, back-transformed for reporting
# ============================================================================

library(readr)
library(dismo)
library(gbm)

# ----------------------------------------------------------------------------
# 1. Load models and data
# ----------------------------------------------------------------------------

cat("Loading models...\n")
brt_SPC   <- readRDS("3_models/brt_SPC.rds")
brt_TP    <- readRDS("3_models/brt_TP.rds")
brt_TN    <- readRDS("3_models/brt_TN.rds")
brt_bloom <- readRDS("3_models/brt_bloom_fitted.rds")

cat("Loading model-ready dataset...\n")
dat <- as.data.frame(
  read_csv("2_incremental/ucfr_model_ready.csv", show_col_types = FALSE)
)

for (d in c("2_incremental", "4_products/diagnostics")) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Site order upstream to downstream
site_order <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")

# ----------------------------------------------------------------------------
# 2. NorthWEST temperature projections (site-specific)
# ----------------------------------------------------------------------------

northwest_temp <- data.frame(
  Site = c("BM","BM","BM","BN","BN","BN","DL","DL","DL",
           "FH","FH","FH","GR","GR","GR","HU","HU","HU","MS","MS","MS"),
  Year = rep(c(2011, 2040, 2080), 7),
  Temp = c(18.74, 20.31, 21.45,
           14.59, 15.99, 17.01,
           15.85, 17.30, 18.36,
           18.78, 20.35, 21.49,
           16.67, 18.15, 19.23,
           20.24, 21.87, 23.05,
           18.03, 19.57, 20.69)
)

# Linear interpolation function for a given site and target year
interp_temp <- function(site, target_year) {
  d <- northwest_temp[northwest_temp$Site == site, ]
  approx(d$Year, d$Temp, xout = target_year)$y
}

# ----------------------------------------------------------------------------
# 3. Scenario definitions
# ----------------------------------------------------------------------------

# Time slices
time_slices <- c(2020, 2050, 2080)

# DSF shifts in days (negative = earlier freshet)
# Snow scenarios from Larson
dsf_shifts <- data.frame(
  Scenario  = c("HS", "LS", "HR", "LR"),
  Year_2020 = c(-14,  -19,  -14,  -19),   # rain uses same as snow + flag
  Year_2050 = c(-32,  -40,  -32,  -40),
  Year_2080 = c(-46,  -65,  -46,  -65),
  Rain_flag = c(FALSE, FALSE, TRUE, TRUE)
)

# Q changes as % from baseline (positive = increase)
q_changes <- data.frame(
  Scenario  = c("HS",   "LS",    "HR",   "LR"),
  Year_2020 = c(-7.5,  -11.1,   12.5,   75.0),
  Year_2050 = c(-11.7, -17.6,   48.8,   45.1),
  Year_2080 = c(-16.8, -26.1,   50.7,   50.9)
)

scenario_labels <- c(
  HS = "High Precip / Snow",
  LS = "Low Precip / Snow",
  HR = "High Precip / Rain",
  LR = "Low Precip / Rain"
)

# ----------------------------------------------------------------------------
# 4. Compute baseline conditions per site (2020)
# ----------------------------------------------------------------------------

cat("Computing site baseline conditions (2020)...\n\n")

baseline <- data.frame(
  Site               = site_order,
  anomaly_base       = NA_real_,
  Q_obs_base         = NA_real_,
  DSF_base           = NA_real_,
  Temp_base_2020     = NA_real_,
  anomaly_min        = NA_real_,
  anomaly_max        = NA_real_,
  Q_obs_min          = NA_real_,
  Q_obs_max          = NA_real_,
  DSF_min            = NA_real_,
  DSF_max            = NA_real_,
  Temp_min           = NA_real_,
  Temp_max           = NA_real_,
  stringsAsFactors   = FALSE
)

for (i in seq_along(site_order)) {
  s <- site_order[i]
  d <- dat[dat$Site == s & !is.na(dat$anomaly) &
             !is.na(dat$Q_obs_cfs) & !is.na(dat$Days_Since_Freshet), ]
  
  baseline$anomaly_base[i]   <- mean(d$anomaly,           na.rm = TRUE)
  baseline$Q_obs_base[i]     <- mean(d$Q_obs_cfs,         na.rm = TRUE)
  baseline$DSF_base[i]       <- mean(d$Days_Since_Freshet, na.rm = TRUE)
  baseline$Temp_base_2020[i] <- interp_temp(s, 2020)
  
  # Observed ranges for extrapolation flagging
  baseline$anomaly_min[i] <- min(d$anomaly,            na.rm = TRUE)
  baseline$anomaly_max[i] <- max(d$anomaly,            na.rm = TRUE)
  baseline$Q_obs_min[i]   <- min(d$Q_obs_cfs,          na.rm = TRUE)
  baseline$Q_obs_max[i]   <- max(d$Q_obs_cfs,          na.rm = TRUE)
  baseline$DSF_min[i]     <- min(d$Days_Since_Freshet,  na.rm = TRUE)
  baseline$DSF_max[i]     <- max(d$Days_Since_Freshet,  na.rm = TRUE)
  baseline$Temp_min[i]    <- min(d$Temp_oC,             na.rm = TRUE)
  baseline$Temp_max[i]    <- max(d$Temp_oC,             na.rm = TRUE)
}

cat("Baseline conditions (2020):\n")
cat(sprintf("  %-6s  %8s  %10s  %8s  %8s\n",
            "Site", "anomaly", "Q_obs_cfs", "DSF", "Temp_2020"))
cat(paste(rep("-", 50), collapse = ""), "\n")
for (i in seq_along(site_order)) {
  cat(sprintf("  %-6s  %8.3f  %10.1f  %8.1f  %8.2f\n",
              baseline$Site[i],
              baseline$anomaly_base[i],
              baseline$Q_obs_base[i],
              baseline$DSF_base[i],
              baseline$Temp_base_2020[i]))
}

# ----------------------------------------------------------------------------
# 5. Run prediction chain helper
# ----------------------------------------------------------------------------

run_chain <- function(anomaly, Q_obs_cfs, Temp_oC, Days_Since_Freshet) {
  newdat <- data.frame(
    anomaly            = anomaly,
    Q_obs_cfs          = Q_obs_cfs,
    Temp_oC            = Temp_oC,
    Days_Since_Freshet = Days_Since_Freshet
  )
  
  pred_SPC   <- predict(brt_SPC,   newdata = newdat,
                        n.trees = brt_SPC$gbm.call$best.trees)
  pred_logTP <- predict(brt_TP,    newdata = newdat,
                        n.trees = brt_TP$gbm.call$best.trees)
  pred_logTN <- predict(brt_TN,    newdata = newdat,
                        n.trees = brt_TN$gbm.call$best.trees)
  
  bloom_dat <- data.frame(
    pred_SPC           = pred_SPC,
    pred_logTP         = pred_logTP,
    pred_logTN         = pred_logTN,
    anomaly            = anomaly,
    Q_obs_cfs          = Q_obs_cfs,
    Temp_oC            = Temp_oC,
    Days_Since_Freshet = Days_Since_Freshet
  )
  
  pred_logCHLa <- predict(brt_bloom, newdata = bloom_dat,
                          n.trees = brt_bloom$gbm.call$best.trees)
  
  list(
    pred_SPC     = pred_SPC,
    pred_logTP   = pred_logTP,
    pred_logTN   = pred_logTN,
    pred_logCHLa = pred_logCHLa,
    pred_CHLa    = 10^pred_logCHLa
  )
}

# ----------------------------------------------------------------------------
# 6. Run scenarios
# ----------------------------------------------------------------------------

cat("\nRunning scenarios...\n\n")

results <- data.frame(
  Site              = character(),
  Scenario          = character(),
  Scenario_label    = character(),
  Year              = integer(),
  Rain_uncertainty  = logical(),
  anomaly           = numeric(),
  Q_obs_cfs         = numeric(),
  Temp_oC           = numeric(),
  Days_Since_Freshet = numeric(),
  pred_SPC          = numeric(),
  pred_logTP        = numeric(),
  pred_logTN        = numeric(),
  pred_logCHLa      = numeric(),
  pred_CHLa         = numeric(),
  extrapolation_flag = logical(),
  stringsAsFactors  = FALSE
)

for (s in site_order) {
  b <- baseline[baseline$Site == s, ]
  
  # Baseline prediction (2020)
  base_pred <- run_chain(b$anomaly_base, b$Q_obs_base,
                         b$Temp_base_2020, b$DSF_base)
  
  results <- rbind(results, data.frame(
    Site               = s,
    Scenario           = "Baseline",
    Scenario_label     = "Baseline (2020)",
    Year               = 2020,
    Rain_uncertainty   = FALSE,
    anomaly            = b$anomaly_base,
    Q_obs_cfs          = b$Q_obs_base,
    Temp_oC            = b$Temp_base_2020,
    Days_Since_Freshet = b$DSF_base,
    pred_SPC           = base_pred$pred_SPC,
    pred_logTP         = base_pred$pred_logTP,
    pred_logTN         = base_pred$pred_logTN,
    pred_logCHLa       = base_pred$pred_logCHLa,
    pred_CHLa          = base_pred$pred_CHLa,
    extrapolation_flag = FALSE,
    stringsAsFactors   = FALSE
  ))
  
  # Climate scenarios
  for (sc in c("HS", "LS", "HR", "LR")) {
    sc_dsf <- dsf_shifts[dsf_shifts$Scenario == sc, ]
    sc_q   <- q_changes[q_changes$Scenario == sc, ]
    rain_flag <- dsf_shifts$Rain_flag[dsf_shifts$Scenario == sc]
    
    for (yr in time_slices) {
      yr_col <- paste0("Year_", yr)
      
      # Apply changes
      q_pct    <- sc_q[[yr_col]] / 100
      dsf_shift <- sc_dsf[[yr_col]]
      temp_proj <- interp_temp(s, yr)
      
      new_anomaly <- b$anomaly_base * (1 + q_pct)
      new_Q_obs   <- b$Q_obs_base   * (1 + q_pct)
      new_DSF     <- b$DSF_base     + dsf_shift
      new_Temp    <- temp_proj
      
      # Check extrapolation
      extrap <- new_anomaly < b$anomaly_min | new_anomaly > b$anomaly_max |
        new_Q_obs   < b$Q_obs_min   | new_Q_obs   > b$Q_obs_max   |
        new_DSF     < b$DSF_min     | new_DSF     > b$DSF_max     |
        new_Temp    < b$Temp_min    | new_Temp    > b$Temp_max
      
      pred <- run_chain(new_anomaly, new_Q_obs, new_Temp, new_DSF)
      
      results <- rbind(results, data.frame(
        Site               = s,
        Scenario           = sc,
        Scenario_label     = scenario_labels[sc],
        Year               = yr,
        Rain_uncertainty   = rain_flag,
        anomaly            = new_anomaly,
        Q_obs_cfs          = new_Q_obs,
        Temp_oC            = new_Temp,
        Days_Since_Freshet = new_DSF,
        pred_SPC           = pred$pred_SPC,
        pred_logTP         = pred$pred_logTP,
        pred_logTN         = pred$pred_logTN,
        pred_logCHLa       = pred$pred_logCHLa,
        pred_CHLa          = pred$pred_CHLa,
        extrapolation_flag = extrap,
        stringsAsFactors   = FALSE
      ))
    }
  }
}

# ----------------------------------------------------------------------------
# 7. Summary table
# ----------------------------------------------------------------------------

cat("--- Scenario Predictions: pred_CHLa (mg/m²) ---\n\n")

for (s in site_order) {
  cat(sprintf("Site: %s\n", s))
  base_chla <- results$pred_CHLa[results$Site == s &
                                   results$Scenario == "Baseline"]
  cat(sprintf("  Baseline 2020: %.2f mg/m²  [log10=%.3f]\n\n",
              base_chla, log10(base_chla)))
  cat(sprintf("  %-6s  %-5s  %10s  %10s  %6s  %s\n",
              "Scen", "Year", "CHLa mg/m²", "log10CHLa", "Δ%", "Extrap"))
  cat("  ", paste(rep("-", 52), collapse = ""), "\n")
  
  sc_dat <- results[results$Site == s & results$Scenario != "Baseline", ]
  for (i in seq_len(nrow(sc_dat))) {
    delta_pct <- 100 * (sc_dat$pred_CHLa[i] - base_chla) / base_chla
    extrap_flag <- if (sc_dat$extrapolation_flag[i]) "!" else ""
    rain_flag   <- if (sc_dat$Rain_uncertainty[i]) "*" else ""
    cat(sprintf("  %-6s  %-5d  %10.2f  %10.3f  %+6.1f%%  %s%s\n",
                sc_dat$Scenario[i],
                sc_dat$Year[i],
                sc_dat$pred_CHLa[i],
                sc_dat$pred_logCHLa[i],
                delta_pct,
                extrap_flag,
                rain_flag))
  }
  cat("\n")
}

cat("! = extrapolation beyond observed training range\n")
cat("* = rain scenario — DSF uncertainty higher\n\n")

# ----------------------------------------------------------------------------
# 8. Plots
# ----------------------------------------------------------------------------

cat("Generating scenario plots...\n")

site_cols <- c(DL = "#E41A1C", GR = "#377EB8", BN = "#4DAF4A",
               MS = "#984EA3", BM = "#FF7F00", HU = "#A65628",
               FH = "#F781BF")

scen_cols <- c(HS = "#2471A3", LS = "#1A5276",
               HR = "#CB4335", LR = "#922B21")

scen_lty  <- c(HS = 1, LS = 2, HR = 1, LR = 2)

pdf("4_products/diagnostics/scenario_results.pdf",
    width = 14, height = 10)

par(mfrow = c(3, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))

for (s in site_order) {
  d_site <- results[results$Site == s, ]
  base   <- d_site$pred_logCHLa[d_site$Scenario == "Baseline"]
  
  # Y range across all scenarios
  y_rng <- range(d_site$pred_logCHLa, na.rm = TRUE)
  y_rng <- c(min(y_rng[1], base) - 0.1, max(y_rng[2], base) + 0.1)
  
  plot(0, 0, type = "n",
       xlim = c(2018, 2082),
       ylim = y_rng,
       xlab = "Year",
       ylab = "log10(CHLa mg/m²)",
       main = sprintf("%s — Baseline: %.2f mg/m²",
                      s, 10^base))
  
  abline(h = base, col = "grey60", lty = 3, lwd = 1.5)
  
  for (sc in c("HS", "LS", "HR", "LR")) {
    d_sc <- d_site[d_site$Scenario == sc, ]
    d_sc <- d_sc[order(d_sc$Year), ]
    
    # Add baseline point
    yvals <- c(base, d_sc$pred_logCHLa)
    xvals <- c(2020, d_sc$Year)
    
    lines(xvals, yvals,
          col = scen_cols[sc], lty = scen_lty[sc], lwd = 1.8)
    points(xvals, yvals,
           col  = scen_cols[sc],
           pch  = ifelse(d_sc$Rain_uncertainty, 2, 16),
           cex  = 0.9)
    
    # Mark extrapolation points
    extrap_idx <- which(d_sc$extrapolation_flag)
    if (length(extrap_idx) > 0) {
      points(d_sc$Year[extrap_idx], d_sc$pred_logCHLa[extrap_idx],
             pch = 4, col = "red", cex = 1.2, lwd = 2)
    }
  }
}

# Legend panel
plot(0, 0, type = "n", xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "")
legend("center",
       legend = c("HS: High Precip/Snow",
                  "LS: Low Precip/Snow",
                  "HR: High Precip/Rain",
                  "LR: Low Precip/Rain",
                  "Baseline (2020)",
                  "Rain scenario (DSF uncertain)",
                  "Extrapolation warning"),
       col    = c(scen_cols, "grey60", scen_cols["HR"], "red"),
       lty    = c(scen_lty, 3, NA, NA),
       pch    = c(rep(16, 4), NA, 2, 4),
       lwd    = c(rep(1.8, 4), 1.5, NA, NA),
       bty    = "n", cex = 0.85)

mtext("Climate Scenario Predictions — UCFR Filamentous Algae",
      outer = TRUE, cex = 1.1, font = 2)

dev.off()
cat("Scenario plots saved to 4_products/diagnostics/scenario_results.pdf\n")

# ----------------------------------------------------------------------------
# 9. Save results
# ----------------------------------------------------------------------------

write_csv(results, "2_incremental/scenario_predictions.csv")
cat("Predictions saved to 2_incremental/scenario_predictions.csv\n")
cat("Done.\n")