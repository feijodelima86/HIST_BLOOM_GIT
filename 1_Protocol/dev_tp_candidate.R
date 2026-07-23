# ============================================================================
# dev_tp_candidate.R
# Model-selection sandbox for the TP submodel.
# Same idea as dev_chla_candidate.R -- edit tp_k_spec, re-run.
# ============================================================================

source("1_protocol/dev_model_selection_utils.R")

dat <- read.csv("2_incremental/ucfr_model_ready.csv", stringsAsFactors = FALSE)

# --- Response -------------------------------------------------------------
# Derived directly from raw TP_mg_L, matching 09_tp_submodel.R's own
# reasoning: the model-ready logTP_mg_L column (log10(1+TP_mg_L)) is close
# to linear over the observed TP range and barely variance-stabilizing.
TP_SCALE <- "log10_ugL"   # "log10_ugL" (recommended) or "log1p_mgL"
if (TP_SCALE == "log10_ugL") {
  dat$TP_response <- log10(dat$TP_mg_L * 1000)
} else if (TP_SCALE == "log1p_mgL") {
  dat$TP_response <- log10(1 + dat$TP_mg_L)
} else {
  stop("TP_SCALE must be 'log10_ugL' or 'log1p_mgL'")
}
RESPONSE <- "TP_response"

# ============================================================================
# CANDIDATE PREDICTOR SET -- this is the block you edit
# ----------------------------------------------------------------------------
# Active = current 09_tp_submodel.R set. Everything else that exists in
# ucfr_model_ready.csv is listed commented-out below.
# ============================================================================
tp_k_spec <- c(
  anomaly            = 5,
 # logQ_obs_cfs       = 5,
  Days_Since_Freshet = 5,
  Temp_oC            = 5     # NOTE: current TP submodel uses k=5 here, but
  # M1 uses k=10 for Temp_oC. Worth deciding
  # whether that difference is intentional.
  
  # --- available, not currently in TP submodel -- uncomment to test ---
  # , pH                 = 5
  # , SPC                = 5
  # , TDS                = 5
  # , TURBIDITY          = 5
  # , TN_mg_L            = 5
  # , SRP_mg_L           = 5
  # , NH4_mg_L           = 5
  # , NO3_mg_L           = 5
  # , DIN_mg_L           = 5
  # , Q_peak_cfs         = 5
  # , Q_baseflow_cfs     = 5
)

TP_SITE_RE        <- F
TP_OUTLIER_SD     <- 2.0
TP_OUTLIER_ROUNDS <- 2
TP_SELECT_PENALTY <- T
TP_DO_LOSO        <- TRUE
TP_DO_LOYO        <- TRUE
TP_DO_DROP_AIC    <- TRUE

result_tp <- run_trial(
  label          = "TP candidate",
  response       = RESPONSE,
  k_spec         = tp_k_spec,
  data_full      = dat,
  site_re        = TP_SITE_RE,
  outlier_sd     = TP_OUTLIER_SD,
  outlier_rounds = TP_OUTLIER_ROUNDS,
  select_penalty = TP_SELECT_PENALTY,
  do_loso        = TP_DO_LOSO,
  do_loyo        = TP_DO_LOYO,
  do_drop_aic    = TP_DO_DROP_AIC
)

