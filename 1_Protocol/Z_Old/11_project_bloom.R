# =============================================================================
# 11_project_bloom.R
# UCFR Cladophora Bloom Prediction Pipeline
#
# PURPOSE:
#   Apply the full BRT prediction chain (SPC → TP → TN → bloom) to the
#   projection grid from Script 10. Produces predicted log10(CHLa) for each
#   site × year × scenario × date point, with median and uncertainty bounds
#   propagated from ESM q25/q75 hydrological deltas.
#
# INPUTS:
#   2_incremental/projection_grid.csv      — from Script 10
#   3_models/brt_SPC.rds                   — BRT submodel for SPC
#   3_models/brt_TP.rds                    — BRT submodel for log10(TP)
#   3_models/brt_TN.rds                    — BRT submodel for log10(TN)
#   3_models/brt_bloom_fitted.rds          — BRT bloom model
#
# OUTPUTS:
#   2_incremental/projections_monthly.csv  — full projection output
#                                            one row per site × year × scenario
#                                            × date point, with predicted chain
#                                            values and uncertainty bounds
#
# PREDICTION CHAIN:
#   Pass 1 (median):  Q_obs_cfs_med, anomaly_med, DSF_med, Temp_oC
#     → pred_SPC_med, pred_logTP_med, pred_logTN_med → pred_logCHLa_med
#   Pass 2 (q25):     Q_obs_cfs_q25, anomaly_q25, DSF_q25, Temp_oC
#     → pred_SPC_q25, pred_logTP_q25, pred_logTN_q25 → pred_logCHLa_q25
#   Pass 3 (q75):     Q_obs_cfs_q75, anomaly_q75, DSF_q75, Temp_oC
#     → pred_SPC_q75, pred_logTP_q75, pred_logTN_q75 → pred_logCHLa_q75
#
# NOTE: Temperature has no ESM spread (NCCV single ensemble mean), so Temp_oC
#       is the same across all three passes.
#
# EXTRAPOLATION FLAGGING:
#   Each predictor is checked against the training envelope (min/max from
#   ucfr_model_ready.csv). Rows with any predictor outside the envelope
#   are flagged with extrapolation_flag = TRUE.
#
# AUTHOR: [Rafa]
# =============================================================================


# -----------------------------------------------------------------------------
# 0. USER PARAMETERS
# -----------------------------------------------------------------------------

baseline_start <- 1998
baseline_end   <- 2022

site_order <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")


# -----------------------------------------------------------------------------
# 1. LOAD INPUTS
# -----------------------------------------------------------------------------

cat("Loading inputs...\n")

library(gbm)
library(dismo)

grid     <- read.csv("2_incremental/projection_grid.csv",   stringsAsFactors = FALSE)
obs_data <- read.csv("2_incremental/ucfr_model_ready.csv",  stringsAsFactors = FALSE)

brt_SPC   <- readRDS("3_models/brt_SPC.rds")
brt_TP    <- readRDS("3_models/brt_TP.rds")
brt_TN    <- readRDS("3_models/brt_TN.rds")
brt_bloom <- readRDS("3_models/brt_bloom_fitted.rds")

cat("  Projection grid rows:", nrow(grid), "\n")


# -----------------------------------------------------------------------------
# 2. DEFINE TRAINING ENVELOPE FOR EXTRAPOLATION FLAGGING
# -----------------------------------------------------------------------------
# Use the full observed dataset (not just baseline) as the training envelope,
# since the models were trained on all available observations.

cat("Defining training envelope...\n")

env_predictors <- c("Q_obs_cfs", "anomaly", "Days_Since_Freshet", "Temp_oC")

training_env <- lapply(env_predictors, function(col) {
  if (col %in% names(obs_data)) {
    c(min(obs_data[[col]], na.rm = TRUE),
      max(obs_data[[col]], na.rm = TRUE))
  } else {
    NULL
  }
})
names(training_env) <- env_predictors

cat("  Training envelope defined for:", paste(env_predictors, collapse = ", "), "\n")


# -----------------------------------------------------------------------------
# 3. PREDICTION FUNCTION
# -----------------------------------------------------------------------------
# Takes a data frame with columns:
#   Q_obs_cfs, anomaly, Days_Since_Freshet, Temp_oC
# Returns a data frame with full chain predictions:
#   pred_SPC, pred_logTP, pred_logTN, pred_logCHLa

run_chain <- function(input_df) {
  
  # --- Intermediate submodels ---
  # Each submodel was trained on: anomaly, Q_obs_cfs, Temp_oC, Days_Since_Freshet
  
  sub_predictors <- c("anomaly", "Q_obs_cfs", "Temp_oC", "Days_Since_Freshet")
  
  pred_SPC   <- predict(brt_SPC,   newdata = input_df[, sub_predictors],
                        n.trees = brt_SPC$gbm.call$best.trees,   type = "response")
  pred_logTP <- predict(brt_TP,    newdata = input_df[, sub_predictors],
                        n.trees = brt_TP$gbm.call$best.trees,    type = "response")
  pred_logTN <- predict(brt_TN,    newdata = input_df[, sub_predictors],
                        n.trees = brt_TN$gbm.call$best.trees,    type = "response")
  
  # --- Bloom model ---
  # Trained on: pred_SPC, pred_logTP, pred_logTN,
  #             anomaly, Q_obs_cfs, Temp_oC, Days_Since_Freshet
  
  bloom_input <- data.frame(
    pred_SPC           = pred_SPC,
    pred_logTP         = pred_logTP,
    pred_logTN         = pred_logTN,
    anomaly            = input_df$anomaly,
    Q_obs_cfs          = input_df$Q_obs_cfs,
    Temp_oC            = input_df$Temp_oC,
    Days_Since_Freshet = input_df$Days_Since_Freshet
  )
  
  pred_logCHLa <- predict(brt_bloom, newdata = bloom_input,
                          n.trees = brt_bloom$gbm.call$best.trees, type = "response")
  
  data.frame(
    pred_SPC     = pred_SPC,
    pred_logTP   = pred_logTP,
    pred_logTN   = pred_logTN,
    pred_logCHLa = pred_logCHLa
  )
}


# -----------------------------------------------------------------------------
# 4. BUILD INPUT FRAMES FOR EACH PASS
# -----------------------------------------------------------------------------

cat("Building prediction input frames...\n")

# Pass 1: median hydrological predictors
input_med <- data.frame(
  Q_obs_cfs          = grid$Q_obs_cfs_med,
  anomaly            = grid$anomaly_med,
  Days_Since_Freshet = grid$Days_Since_Freshet_med,
  Temp_oC            = grid$Temp_oC
)

# Pass 2: q25 hydrological predictors (lower bound)
input_q25 <- data.frame(
  Q_obs_cfs          = grid$Q_obs_cfs_q25,
  anomaly            = grid$anomaly_q25,
  Days_Since_Freshet = grid$Days_Since_Freshet_q25,
  Temp_oC            = grid$Temp_oC
)

# Pass 3: q75 hydrological predictors (upper bound)
input_q75 <- data.frame(
  Q_obs_cfs          = grid$Q_obs_cfs_q75,
  anomaly            = grid$anomaly_q75,
  Days_Since_Freshet = grid$Days_Since_Freshet_q75,
  Temp_oC            = grid$Temp_oC
)


# -----------------------------------------------------------------------------
# 5. RUN PREDICTION CHAIN FOR ALL THREE PASSES
# -----------------------------------------------------------------------------

cat("Running prediction chain...\n")
cat("  Pass 1: median...\n")
chain_med <- run_chain(input_med)

cat("  Pass 2: q25...\n")
chain_q25 <- run_chain(input_q25)

cat("  Pass 3: q75...\n")
chain_q75 <- run_chain(input_q75)

cat("  Prediction chain complete.\n")


# -----------------------------------------------------------------------------
# 6. FLAG EXTRAPOLATION
# -----------------------------------------------------------------------------
# Flag any row where a projected predictor (median pass) falls outside the
# training envelope. Conservative: flags on median only; q25/q75 may also
# exceed envelope but the flag is row-level.

cat("Flagging extrapolation...\n")

flag_matrix <- data.frame(
  Q_out    = input_med$Q_obs_cfs < training_env$Q_obs_cfs[1] |
    input_med$Q_obs_cfs > training_env$Q_obs_cfs[2],
  anom_out = input_med$anomaly < training_env$anomaly[1] |
    input_med$anomaly > training_env$anomaly[2],
  DSF_out  = input_med$Days_Since_Freshet < training_env$Days_Since_Freshet[1] |
    input_med$Days_Since_Freshet > training_env$Days_Since_Freshet[2],
  Temp_out = input_med$Temp_oC < training_env$Temp_oC[1] |
    input_med$Temp_oC > training_env$Temp_oC[2]
)

extrapolation_flag <- rowSums(flag_matrix) > 0

cat("  Rows flagged for extrapolation:", sum(extrapolation_flag),
    "of", nrow(grid), "\n")
cat("  Extrapolation by predictor:\n")
print(colSums(flag_matrix))


# -----------------------------------------------------------------------------
# 7. ASSEMBLE OUTPUT
# -----------------------------------------------------------------------------

cat("Assembling output...\n")

projections <- data.frame(
  # Identifiers
  site       = grid$site,
  year       = grid$year,
  scenario   = grid$scenario,
  date_point = grid$date_point,
  doy        = grid$doy,
  month      = grid$month,
  
  # Baseline means (for reference)
  base_Temp_oC            = grid$base_Temp_oC,
  base_Q_obs_cfs          = grid$base_Q_obs_cfs,
  base_anomaly            = grid$base_anomaly,
  base_Days_Since_Freshet = grid$base_Days_Since_Freshet,
  
  # Projected predictors (median)
  Temp_oC                = grid$Temp_oC,
  Q_obs_cfs_med          = grid$Q_obs_cfs_med,
  Q_obs_cfs_q25          = grid$Q_obs_cfs_q25,
  Q_obs_cfs_q75          = grid$Q_obs_cfs_q75,
  anomaly_med            = grid$anomaly_med,
  anomaly_q25            = grid$anomaly_q25,
  anomaly_q75            = grid$anomaly_q75,
  DSF_med                = grid$Days_Since_Freshet_med,
  DSF_q25                = grid$Days_Since_Freshet_q25,
  DSF_q75                = grid$Days_Since_Freshet_q75,
  
  # Chain predictions — median pass
  pred_SPC_med     = chain_med$pred_SPC,
  pred_logTP_med   = chain_med$pred_logTP,
  pred_logTN_med   = chain_med$pred_logTN,
  pred_logCHLa_med = chain_med$pred_logCHLa,
  
  # Chain predictions — q25 pass (raw, from q25 hydro inputs)
  pred_SPC_q25     = chain_q25$pred_SPC,
  pred_logTP_q25   = chain_q25$pred_logTP,
  pred_logTN_q25   = chain_q25$pred_logTN,
  pred_logCHLa_q25 = chain_q25$pred_logCHLa,
  
  # Chain predictions — q75 pass (raw, from q75 hydro inputs)
  pred_SPC_q75     = chain_q75$pred_SPC,
  pred_logTP_q75   = chain_q75$pred_logTP,
  pred_logTN_q75   = chain_q75$pred_logTN,
  pred_logCHLa_q75 = chain_q75$pred_logCHLa,
  
  # Envelope-corrected uncertainty bands for plotting
  # BRT nonlinearity means q25/q75 inputs don't map monotonically to
  # q25/q75 outputs (e.g. high flow suppresses bloom via dilution).
  # Band is defined as min/max across all three passes, guaranteeing
  # the median line always falls within the ribbon.
  pred_logCHLa_lo = pmin(chain_med$pred_logCHLa,
                         chain_q25$pred_logCHLa,
                         chain_q75$pred_logCHLa),
  pred_logCHLa_hi = pmax(chain_med$pred_logCHLa,
                         chain_q25$pred_logCHLa,
                         chain_q75$pred_logCHLa),
  
  # Back-transformed CHLa (µg/L) for reporting
  pred_CHLa_med = 10^chain_med$pred_logCHLa,
  pred_CHLa_q25 = 10^chain_q25$pred_logCHLa,
  pred_CHLa_q75 = 10^chain_q75$pred_logCHLa,
  pred_CHLa_lo  = 10^pmin(chain_med$pred_logCHLa,
                          chain_q25$pred_logCHLa,
                          chain_q75$pred_logCHLa),
  pred_CHLa_hi  = 10^pmax(chain_med$pred_logCHLa,
                          chain_q25$pred_logCHLa,
                          chain_q75$pred_logCHLa),
  
  # Extrapolation flag
  extrapolation_flag = extrapolation_flag,
  
  stringsAsFactors = FALSE
)

# Enforce site order
projections$site <- factor(projections$site, levels = site_order)
projections <- projections[order(
  projections$site,
  projections$scenario,
  projections$date_point,
  projections$year
), ]
projections$site <- as.character(projections$site)


# -----------------------------------------------------------------------------
# 8. SANITY CHECKS
# -----------------------------------------------------------------------------

cat("\n--- Sanity checks ---\n")

# NA check
na_counts <- colSums(is.na(projections))
if (any(na_counts > 0)) {
  cat("  WARNING: NAs in columns:\n")
  print(na_counts[na_counts > 0])
} else {
  cat("  No NAs in projections.\n")
}

# Predicted CHLa range
cat("\n  pred_CHLa_med range (µg/L):",
    round(min(projections$pred_CHLa_med), 2), "to",
    round(max(projections$pred_CHLa_med), 2), "\n")

# Spot check: DL / august / ssp585 bloom trajectory every 10 years
cat("\n  Bloom trajectory: DL / august / ssp585 (every 10 years)\n")
traj_check <- projections[
  projections$site == "DL" &
    projections$date_point == "august" &
    projections$scenario == "ssp585" &
    projections$year %% 10 == 0,
  c("year", "Temp_oC", "Q_obs_cfs_med", "pred_SPC_med",
    "pred_logTP_med", "pred_logCHLa_med", "pred_CHLa_med")
]
print(traj_check)

# Spot check: scenario contrast at 2080 across sites for August
cat("\n  2080 August bloom by site and scenario (median, µg/L):\n")
slice_2080 <- projections[
  projections$year == 2080 &
    projections$date_point == "august",
  c("site", "scenario", "pred_CHLa_med", "pred_CHLa_q25", "pred_CHLa_q75",
    "extrapolation_flag")
]
print(slice_2080[order(slice_2080$scenario, slice_2080$site), ])

# Extrapolation summary by site and scenario
cat("\n  Extrapolation flags by site × scenario:\n")
extrap_summary <- aggregate(
  extrapolation_flag ~ site + scenario,
  data = projections,
  FUN  = sum
)
print(extrap_summary[order(extrap_summary$scenario, extrap_summary$site), ])


# -----------------------------------------------------------------------------
# 9. WRITE OUTPUT
# -----------------------------------------------------------------------------

out_path <- "2_incremental/projections_monthly.csv"
write.csv(projections, out_path, row.names = FALSE)
cat("\nOutput written to:", out_path, "\n")
cat("Done.\n")