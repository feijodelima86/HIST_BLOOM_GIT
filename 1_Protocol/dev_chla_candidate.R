# ============================================================================
# dev_chla_candidate.R
# Model-selection sandbox for the CHLa (M1) submodel.
# NOT manuscript-facing, NOT the pipeline -- see dev_model_selection_utils.R
# for what this actually runs.
#
# HOW TO USE
#   - To drop a variable: delete or comment out its line in chla_k_spec.
#   - To add a variable: add a line "varname = k". k=5 is a reasonable
#     default guess for a monotone-looking relationship; k=10 if you expect
#     a hump/inflection, matching the logTP/Temp convention already in M1.
#   - Re-run the whole script. Everything downstream reads chla_k_spec.
#
# Assumes this file sits in the same folder as dev_model_selection_utils.R
# and is run with the project root as the working directory (same
# convention as every other pipeline script). Adjust the source() path
# below if you put it somewhere else.
# ============================================================================

source("1_protocol/dev_model_selection_utils.R")

dat <- read.csv("2_incremental/ucfr_model_ready.csv", stringsAsFactors = FALSE)

# --- Response -----------------------------------------------------------
# Swap to "logAFDM" the same way 10_bloom_model_M1.R does if you ever want to
# try that side; not built out here since AFDM isn't the flagship response.
RESPONSE <- "logCHLa"

# --- lag_y: previous year's site-level annual max of RESPONSE -----------
# To test the model WITHOUT the persistence term: comment out this whole
# block AND remove lag_y from chla_k_spec below.
ann_max <- aggregate(dat[[RESPONSE]] ~ Site + Year, data = dat,
                     FUN = max, na.rm = TRUE)
names(ann_max) <- c("Site", "Year", "annual_max")
ann_max$Year <- ann_max$Year + 1
names(ann_max)[names(ann_max) == "annual_max"] <- "lag_y"
dat <- merge(dat, ann_max, by = c("Site", "Year"), all.x = TRUE)

# --- TN:TP ratio (nutrient stoichiometry) --------------------------------
# Classic Redfield-type nutrient-limitation indicator. Raw ratio is what's
# offered as a candidate below -- keeps the Redfield-16 reference point
# directly readable. (log10 version also computed here in case you want it
# later; not offered as a candidate.) Guarded against TP_mg_L <= 0.
n_bad_tp <- sum(!is.finite(dat$TN_mg_L / dat$TP_mg_L))
if (n_bad_tp > 0) {
  warning(n_bad_tp, " rows with non-finite TN:TP (TP_mg_L <= 0 or NA) -- set to NA.")
}
dat$TN_TP_ratio <- dat$TN_mg_L / dat$TP_mg_L
dat$TN_TP_ratio[!is.finite(dat$TN_TP_ratio)] <- NA_real_
dat$logTN_TP_ratio <- log10(dat$TN_TP_ratio)

# ============================================================================
# CANDIDATE PREDICTOR SET -- this is the block you edit
# ----------------------------------------------------------------------------
# Active = current M1. Everything else that exists in ucfr_model_ready.csv
# is listed commented-out below, ready to uncomment.
# ============================================================================
chla_k_spec <- c(
  lag_y              = 5,
  anomaly            = 5,
  #logQ_obs_cfs       = 5,
  Days_Since_Freshet = 5,
  logTP_mg_L         = 10,
  Temp_oC            = 10
  
  # --- available, not currently in M1 -- uncomment to test ---
  # , pH                 = 5
   , SPC                = 5
  # , TDS                = 5
  # , TURBIDITY          = 5
  # , TN_mg_L            = 5
  # , TN_TP_ratio        = 5   # N:P stoichiometry, raw ratio
  # , SRP_mg_L           = 5
  # , NH4_mg_L           = 5
  # , NO3_mg_L           = 5
  # , DIN_mg_L           = 5
  # , Q_peak_cfs         = 5   # raw component of `anomaly` -- expect concurvity if both included
  # , Q_baseflow_cfs     = 5   # raw component of `anomaly` -- expect concurvity if both included
)

CHLA_SITE_RE        <- F   # include s(Site, bs="re")?
CHLA_OUTLIER_SD     <- 2.0
CHLA_OUTLIER_ROUNDS <- 2      # set to 0 to skip outlier removal (fixed-data comparisons)
CHLA_SELECT_PENALTY <- T  # TRUE = mgcv's select=TRUE extra shrinkage -- see utils header note
CHLA_DO_LOSO        <- TRUE
CHLA_DO_LOYO        <- TRUE
CHLA_DO_DROP_AIC    <- TRUE

result_chla <- run_trial(
  label          = "CHLa candidate",
  response       = RESPONSE,
  k_spec         = chla_k_spec,
  data_full      = dat,
  site_re        = CHLA_SITE_RE,
  outlier_sd     = CHLA_OUTLIER_SD,
  outlier_rounds = CHLA_OUTLIER_ROUNDS,
  select_penalty = CHLA_SELECT_PENALTY,
  do_loso        = CHLA_DO_LOSO,
  do_loyo        = CHLA_DO_LOYO,
  do_drop_aic    = CHLA_DO_DROP_AIC
)