# ============================================================================
# 10_bloom_model_M1.R
# UCFR Filamentous Algae Project
# Standalone fitting script for M1 — six-predictor additive GAM bloom model
#
# Inputs:  2_incremental/ucfr_model_ready.csv   (n=284 visit-level obs)
# Outputs: 3_models/bloom_model_M1.rds          (fitted mgcv::gam object)
#          2_incremental/m1_predictions.csv      (fitted values + residuals)
#          console: scorecard data frame
#
# This script is the single source of truth for what M1 is.
# All downstream scripts (11_temporal_validation.R, projection pipeline)
# load 3_models/bloom_model_M1.rds rather than re-specifying the formula.
#
# Model: RESPONSE ~ s(lag_y) + s(anomaly) + s(logQ_obs_cfs) +
#                   s(Days_Since_Freshet) + s(logTP_mg_L) + s(Temp_oC)
#
# REVISION 2026-07-23: Site random effect s(Site, bs="re") REMOVED from the
# default configuration (see SITE_RE toggle below). A concurvity audit
# (12_concurvity.R) found s(Site) <-> s(logQ_obs_cfs) worst-case concurvity
# of 0.96 -- Site was absorbing nearly all of logQ's between-site signal,
# which also explains why lag_y and logQ looked non-significant in-sample
# despite carrying real held-out predictive skill. Dropping Site clears the
# concurvity flag entirely and every remaining term becomes significant, at
# a real but modest LOYO cost (0.710 -> ~0.67-0.69 depending on shrinkage;
# re-verify under this script's own defaults, not just the sandbox trial).
# See dev_chla_candidate.R trial history for the full comparison.
#
# This is a DIFFERENT call than the TP submodel, which KEEPS Site
# (09_tp_submodel.R unchanged) -- there, Site carries two specific,
# non-substitutable mechanisms (BM's WWTP point source, GR/BN geologic P),
# confirmed by LOSO collapsing to R2 < -2 without it. Here, no comparably
# specific mechanism was identified for what Site was capturing beyond a
# generic downstream/discharge-magnitude gradient that logQ_obs_cfs already
# tracks.
#
# IMPORTANT DOWNSTREAM CONSEQUENCE, not yet resolved: 13_project_bloom.R
# currently gets per-site projection differentiation from Site's fitted
# intercept. Without it, two sites with similar projected
# anomaly/logQ/DSF/logTP/Temp trajectories will get near-identical
# projected bloom trajectories, even if they have real, unexplained
# baseline differences Site RE used to carry. Understand this before
# regenerating any projection-stage figure or table.
#
# lag_y = previous year's site-level annual max of RESPONSE.
#         Represents Cladophora propagule bank legacy (overwinters as
#         basal filaments; prior bloom magnitude seeds next-year establishment).
#         Derived in-script from the response column; not in ucfr_model_ready.
#
# logAFDM is derived in-script (log10(AFDM)) if RESPONSE = "logAFDM".
#         Not present as a column in ucfr_model_ready.csv.
#
# LOSO validation: when SITE_RE = TRUE, uses RE-excluded population-level
# predictions for held-out sites (same mgcv::exclude pattern as
# 09_tp_submodel.R). When SITE_RE = FALSE (current default), no exclude is
# needed -- the held-out site's rows predict exactly like any other row,
# since Site was never a model term. This also means LOSO now tests the
# actual production model's spatial generalization directly, rather than a
# deliberately RE-crippled version of a different model.
#
# Row identity: obs_id (1:n in input row order) is the stable key for the
# predictions CSV. Site+Year+Month alone is not unique — n=284 includes
# double-visit months where two sampling events in the same month have
# distinct Days_Since_Freshet and Q_obs_cfs values.
#
# After refitting: re-run 12_concurvity.R against the new
# bloom_model_M1.rds to confirm the concurvity picture looks as expected.
# ============================================================================

library(mgcv)

# ============================================================================
# CONFIGURATION — edit here only
# ============================================================================
RESPONSE       <- "logCHLa"   # "logCHLa" or "logAFDM"
OUTLIER_SD     <- 2.0         # residual SD threshold for outlier removal
OUTLIER_ROUNDS <- 2           # number of outlier removal passes

# Site random effect -- see REVISION note above. Toggle back to TRUE for a
# sensitivity re-check; LOSO logic below (Section 8) handles both cases.
SITE_RE        <- FALSE

# Per-variable basis dimensions.
# Higher k for logTP and Temp: both showed non-monotonic partial responses
# during variable selection (TP: calcite co-precipitation inflection;
# Temp: unimodal growth optimum). Hydrology predictors expected monotone — k=5.
k_spec <- c(
  lag_y              = 5,
  anomaly            = 5,
  logQ_obs_cfs       = 5,
  Days_Since_Freshet = 5,
  logTP_mg_L         = 10,
  Temp_oC            = 10
)

PREDICTORS <- names(k_spec)   # order defines formula term order
# ============================================================================


# ============================================================================
# 1. READ & VALIDATE INPUT
# ============================================================================
dat <- read.csv("2_incremental/ucfr_model_ready.csv",
                stringsAsFactors = FALSE)

# Hard stop: no gsub on column names. All downstream references use exact names.
required_cols <- c("Site", "Year", "Month",
                   "CHLa", "logCHLa", "AFDM",
                   PREDICTORS[PREDICTORS != "lag_y"])  # lag_y derived below
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Stable row identifier — Site+Year+Month not unique in double-visit months
dat$obs_id <- seq_len(nrow(dat))

# ============================================================================
# 2. DERIVE RESPONSE AND lag_y
# ============================================================================

# --- logAFDM (if toggled) ---
# logCHLa is already in ucfr_model_ready.csv as log10(CHLa).
# logAFDM is not present; derive as log10(AFDM) with zero guard.
if (RESPONSE == "logAFDM") {
  n_zero_afdm <- sum(dat$AFDM <= 0, na.rm = TRUE)
  if (n_zero_afdm > 0) {
    warning(n_zero_afdm, " rows with AFDM <= 0 set to NA before log10.")
    dat$AFDM[dat$AFDM <= 0] <- NA
  }
  dat$logAFDM <- log10(dat$AFDM)
} else if (RESPONSE != "logCHLa") {
  stop("RESPONSE must be 'logCHLa' or 'logAFDM'. Got: '", RESPONSE, "'")
}

# --- lag_y: previous year's site-level annual max of RESPONSE ---
# Aggregated at annual level, merged back to observation level.
# Both within-year visits share the same lag_y (previous year's max).
# Years with no prior-year data get NA lag_y (correctly excluded downstream).
ann_max <- aggregate(dat[[RESPONSE]] ~ Site + Year,
                     data = dat, FUN = max, na.rm = TRUE)
names(ann_max) <- c("Site", "Year", "annual_max")
lag_df <- ann_max
lag_df$Year <- lag_df$Year + 1          # lag: this row's year gets last year's max
names(lag_df)[names(lag_df) == "annual_max"] <- "lag_y"

dat <- merge(dat, lag_df[, c("Site", "Year", "lag_y")],
             by = c("Site", "Year"), all.x = TRUE)

# Restore original row order (merge can shuffle)
dat <- dat[order(dat$obs_id), ]

# ============================================================================
# 3. BUILD FORMULA
# ============================================================================
smooth_terms <- mapply(
  function(v, k) paste0("s(", v, ", k=", k, ")"),
  PREDICTORS, k_spec[PREDICTORS]
)
all_terms <- smooth_terms
if (SITE_RE) all_terms <- c(all_terms, 's(Site, bs="re")')

f_M1 <- as.formula(
  paste(RESPONSE, "~",
        paste(all_terms, collapse = " + "))
)

# ============================================================================
# 4. COMPLETE CASES
# ============================================================================
keep_cols <- c("obs_id", "Site", "Year", "Month",
               RESPONSE, PREDICTORS)
mdat <- dat[complete.cases(dat[, c(RESPONSE, PREDICTORS)]), keep_cols]
mdat$Site <- factor(mdat$Site)

n_raw        <- nrow(dat)
n_complete   <- nrow(mdat)
n_dropped_na <- n_raw - n_complete

# ============================================================================
# 5. OUTLIER REMOVAL — two rounds at ± OUTLIER_SD * SD(residuals)
# ============================================================================
# Fit on complete cases, flag residuals beyond threshold, refit.
# Outlier log accumulated for scorecard transparency.
outlier_log <- list()

for (r in seq_len(OUTLIER_ROUNDS)) {
  m_tmp    <- gam(f_M1, data = mdat, method = "REML")
  resid_r  <- residuals(m_tmp)
  sd_r     <- sd(resid_r)
  flag     <- which(abs(resid_r) > OUTLIER_SD * sd_r)
  
  if (length(flag) == 0) break
  
  outlier_log[[r]] <- data.frame(
    round   = r,
    obs_id  = mdat$obs_id[flag],
    Site    = mdat$Site[flag],
    Year    = mdat$Year[flag],
    Month   = mdat$Month[flag],
    y_obs   = round(mdat[[RESPONSE]][flag], 4),
    resid   = round(resid_r[flag], 4),
    stringsAsFactors = FALSE
  )
  mdat <- mdat[-flag, ]
}

n_clean        <- nrow(mdat)
n_outliers_rm  <- n_complete - n_clean
outlier_detail <- if (length(outlier_log) > 0) do.call(rbind, outlier_log) else
  data.frame()

# ============================================================================
# 6. FIT M1 ON CLEAN DATA
# ============================================================================
m_M1 <- gam(f_M1, data = mdat, method = "REML")
s_M1 <- summary(m_M1)

fitted_vals  <- fitted(m_M1)
resid_vals   <- residuals(m_M1)
r2_adj       <- s_M1$r.sq
dev_expl     <- s_M1$dev.expl
rmse_insamp  <- sqrt(mean(resid_vals^2))
n_model      <- nrow(mdat)

# ============================================================================
# 7. SAVE FITTED OBJECT
# ============================================================================
if (!dir.exists("3_models")) dir.create("3_models", recursive = TRUE)
saveRDS(m_M1, "3_models/bloom_model_M1.rds")

# ============================================================================
# 8. LOSO VALIDATION
# ============================================================================
# SITE_RE = TRUE : RE excluded for held-out sites, predicting at population
#                  level (no site-specific intercept for a novel site).
# SITE_RE = FALSE: no exclude needed -- Site was never a model term, so the
#                  held-out site's rows predict exactly like any other row.
#                  LOSO now tests the actual production model directly.
#
# Reported as:
#   (a) pooled R² vs grand mean of all held-out observations
#   (b) per-site R² and RMSE for site-level diagnostic transparency

sites     <- levels(mdat$Site)
loso_list <- vector("list", length(sites))
names(loso_list) <- sites

for (s in sites) {
  train <- mdat[mdat$Site != s, ]
  test  <- mdat[mdat$Site == s, ]
  
  if (nrow(test) == 0) next
  
  m_loso <- tryCatch(
    gam(f_M1, data = train, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(m_loso)) next
  
  if (SITE_RE) {
    # Assign held-out site a dummy level from training data, then exclude RE
    test2        <- test
    test2$Site   <- factor(levels(train$Site)[1], levels = levels(train$Site))
    pred_loso    <- predict(m_loso, newdata = test2,
                            exclude = 's(Site, bs="re")',
                            newdata.guaranteed = TRUE)
  } else {
    pred_loso <- predict(m_loso, newdata = test)
  }
  
  obs_loso     <- test[[RESPONSE]]
  ss_res       <- sum((obs_loso - pred_loso)^2)
  ss_tot       <- sum((obs_loso - mean(obs_loso))^2)
  r2_site      <- ifelse(ss_tot > 0, 1 - ss_res / ss_tot, NA)
  rmse_site    <- sqrt(mean((obs_loso - pred_loso)^2))
  
  loso_list[[s]] <- data.frame(
    Site      = s,
    n         = nrow(test),
    R2        = round(r2_site, 4),
    RMSE      = round(rmse_site, 4),
    obs_id    = test$obs_id,
    Observed  = obs_loso,
    Predicted = pred_loso,
    stringsAsFactors = FALSE
  )
}

loso_df <- do.call(rbind, loso_list)
row.names(loso_df) <- NULL

# Pooled LOSO R² — all held-out predictions vs grand mean
loso_ss_res    <- sum((loso_df$Observed - loso_df$Predicted)^2)
loso_ss_tot    <- sum((loso_df$Observed - mean(loso_df$Observed))^2)
r2_loso_pooled <- 1 - loso_ss_res / loso_ss_tot
rmse_loso      <- sqrt(mean((loso_df$Observed - loso_df$Predicted)^2))

# Per-site summary (separate from row-level loso_df)
loso_site_summary <- do.call(rbind, lapply(loso_list, function(d) {
  d[1, c("Site", "n", "R2", "RMSE")]
}))
row.names(loso_site_summary) <- NULL

# ============================================================================
# 9. SAVE PREDICTIONS CSV
# ============================================================================
# In-sample fitted values on clean dataset
insamp_df <- data.frame(
  obs_id    = mdat$obs_id,
  Site      = mdat$Site,
  Year      = mdat$Year,
  Month     = mdat$Month,
  Observed  = mdat[[RESPONSE]],
  Fitted    = round(fitted_vals, 6),
  Resid     = round(resid_vals, 6),
  scheme    = "in_sample",
  stringsAsFactors = FALSE
)

# LOSO held-out predictions
loso_out_df <- data.frame(
  obs_id    = loso_df$obs_id,
  Site      = loso_df$Site,
  Year      = mdat$Year[match(loso_df$obs_id, mdat$obs_id)],
  Month     = mdat$Month[match(loso_df$obs_id, mdat$obs_id)],
  Observed  = loso_df$Observed,
  Fitted    = round(loso_df$Predicted, 6),
  Resid     = round(loso_df$Observed - loso_df$Predicted, 6),
  scheme    = "LOSO",
  stringsAsFactors = FALSE
)

pred_out <- rbind(insamp_df, loso_out_df)
if (!dir.exists("2_incremental")) dir.create("2_incremental", recursive = TRUE)
write.csv(pred_out, "2_incremental/m1_predictions.csv", row.names = FALSE)

# ============================================================================
# 10. SCORECARD
# ============================================================================
cat("\n")
cat("============================================================\n")
cat("M1 BLOOM MODEL — SCORECARD\n")
cat("Response:", RESPONSE, "\n")
cat("Site random effect included:", SITE_RE, "\n")
cat("Formula:", deparse(f_M1, width.cutoff = 120), "\n")
cat("============================================================\n\n")

# --- 10a. Sample flow ---
flow_sc <- data.frame(
  Stage              = c("Raw input rows",
                         "Complete cases (all predictors present)",
                         "Dropped (NA in any predictor)",
                         "Outliers removed (2 rounds ±2 SD)",
                         "Clean n for model fitting"),
  n                  = c(n_raw, n_complete, n_dropped_na,
                         n_outliers_rm, n_clean),
  stringsAsFactors   = FALSE
)
cat("--- Sample flow ---\n")
print(flow_sc, row.names = FALSE)
cat("\n")

# --- 10b. In-sample performance ---
insamp_sc <- data.frame(
  Metric  = c("R² (adj)", "Deviance explained (%)", "RMSE"),
  Value   = c(round(r2_adj, 4),
              round(dev_expl * 100, 2),
              round(rmse_insamp, 4)),
  stringsAsFactors = FALSE
)
cat("--- In-sample fit (n =", n_model, ") ---\n")
print(insamp_sc, row.names = FALSE)
cat("\n")

# --- 10c. Smooth term summary ---
sm_names <- rownames(s_M1$s.table)
smooth_sc <- data.frame(
  Term    = sm_names,
  edf     = round(s_M1$s.table[, "edf"], 3),
  F_stat  = round(s_M1$s.table[, "F"], 3),
  p_value = formatC(s_M1$s.table[, "p-value"],
                    format = "e", digits = 2),
  stringsAsFactors = FALSE
)
cat("--- Smooth term summary ---\n")
print(smooth_sc, row.names = FALSE)
cat("\n")

# --- 10d. LOSO spatial validation ---
loso_sc <- data.frame(
  Scheme   = paste0("LOSO", if (SITE_RE) " (RE excluded)" else " (no Site term)"),
  n        = nrow(loso_df),
  R2_pooled = round(r2_loso_pooled, 4),
  RMSE     = round(rmse_loso, 4),
  stringsAsFactors = FALSE
)
cat("--- LOSO spatial validation ---\n")
print(loso_sc, row.names = FALSE)
cat("\n")

cat("--- LOSO per-site breakdown ---\n")
print(loso_site_summary, row.names = FALSE)
cat("\n")

# --- 10e. Outlier detail ---
if (nrow(outlier_detail) > 0) {
  cat("--- Outliers removed ---\n")
  print(outlier_detail, row.names = FALSE)
  cat("\n")
} else {
  cat("--- Outliers removed: none ---\n\n")
}

# --- 10f. Output file locations ---
files_sc <- data.frame(
  File    = c("3_models/bloom_model_M1.rds",
              "2_incremental/m1_predictions.csv"),
  Contents = c("Fitted mgcv::gam object (M1). Load with readRDS().",
               "In-sample fitted values + LOSO held-out predictions."),
  stringsAsFactors = FALSE
)
cat("--- Output files ---\n")
print(files_sc, row.names = FALSE)
cat("\n")

cat("============================================================\n")
cat("Done. M1 is the canonical bloom model.\n")
cat("Downstream scripts load 3_models/bloom_model_M1.rds.\n")
cat("Do not re-specify the formula elsewhere.\n")
cat("Re-run 12_concurvity.R next to confirm the concurvity picture.\n")
cat("============================================================\n")