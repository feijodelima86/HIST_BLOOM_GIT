# ============================================================================
# diag_temp_sensitivity_projection.R
# UCFR Cladophora Bloom Prediction Pipeline -- Temperature-smooth sensitivity
# (defensive/reviewer-response analysis, NOT manuscript-facing)
#
# This is 13_project_bloom.R with a minimal diff, now wrapped to loop over
# all 5 TEMP_VARIANT values in one run instead of one manual edit per pass.
# Nothing upstream of Step 7 moves: M1, tp_submodel, and all validation
# scripts are untouched and unread by this script beyond loading the saved
# .rds objects. All data loading (Sections 1-5) that does NOT depend on the
# treatment happens once, outside the loop; only Sections 6-10 (TP predict,
# M1 recursion, outputs, scorecard) re-run per variant.
#
# Changes relative to the single-variant version already validated:
#   1. VARIANTS_TO_RUN -- run all 5 in sequence instead of one at a time.
#   2. Sections 1-5 (load models, observed climatology, discharge envelope,
#      site expansion, temperature envelope) hoisted OUTSIDE the loop --
#      identical computation every variant, no reason to repeat it 5x.
#   3. Sections 6-10 wrapped in a for-loop over VARIANTS_TO_RUN, each
#      producing its own bloom_projections_members_<variant>.csv /
#      bloom_projections_<variant>.csv, exactly as before.
#   4. Per-variant scorecards are captured into a single data frame
#      (unified_scorecard) instead of printed separately per run, plus one
#      combined extrapolation table (unified_extrap) and one combined
#      end-of-century table (unified_eoc) across all variants. Written to
#      temp_sensitivity_scorecard.csv / temp_sensitivity_extrap.csv /
#      temp_sensitivity_eoc.csv in the same diagnostics folder.
#
# V0 is still a pure pass-through and was already confirmed byte-identical
# (post column-reorder) to bloom_projections_members.csv -- that check is
# not repeated here since it's already passed; see the commented block at
# the bottom if you want to re-run it.
# ============================================================================

library(mgcv)


# ============================================================================
# CONFIGURATION -- edit here only
# ============================================================================

# --- CHANGED: run all variants in one pass instead of picking one ----------
VARIANTS_TO_RUN <- c("V0", "V1", "V2", "V3mild", "V3strong")

TREATMENT <- list(
  V0       = list(type = "none"),
  V1       = list(type = "cap",     cap = 22.2),
  V2       = list(type = "cap",     cap = 28),
  V3mild   = list(type = "decline", threshold = 23, slope = 0.05),
  V3strong = list(type = "decline", threshold = 23, slope = 0.15)
)

unknown_variants <- setdiff(VARIANTS_TO_RUN, names(TREATMENT))
if (length(unknown_variants)) {
  stop("Unknown TEMP_VARIANT(s): ", paste(unknown_variants, collapse = ", "),
       " -- must be among ", paste(names(TREATMENT), collapse = ", "))
}

# --- file paths (unchanged from 13_project_bloom.R) -------------------------
PATH_M1        <- "3_models/bloom_model_M1.rds"
PATH_TP        <- "3_models/tp_submodel.rds"
PATH_OBS       <- "2_incremental/ucfr_model_ready.csv"
PATH_DISCH     <- "2_incremental/ncar_discharge_envelope.csv"
PATH_TEMP      <- "2_incremental/ncar_temperature_envelope.csv"

# --- outputs -> diagnostics/temp_sensitivity/, variant-suffixed ------------
OUT_DIR <- "4_products/diagnostics/temp_sensitivity/"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# --- CHANGED: unified scorecard outputs (one row per variant) --------------
OUT_SCORECARD <- paste0(OUT_DIR, "temp_sensitivity_scorecard.csv")
OUT_EXTRAP    <- paste0(OUT_DIR, "temp_sensitivity_extrap.csv")
OUT_EOC       <- paste0(OUT_DIR, "temp_sensitivity_eoc.csv")

# --- windows (unchanged) -----------------------------------------------------
BASELINE_START <- 1981
BASELINE_END   <- 2010
PROJ_START     <- 2026
PROJ_END       <- 2099
SEED_YEAR      <- 2025

# --- observed-data column names (unchanged) ---------------------------------
OBS_SITE   <- "Site"
OBS_YEAR   <- "Year"
OBS_CHLA   <- "logCHLa"
OBS_LOGQ   <- "logQ_obs_cfs"
OBS_ANOM   <- "anomaly"
OBS_DSF    <- "Days_Since_Freshet"

# --- NCAR discharge envelope column names (unchanged) -----------------------
DQ_REACH   <- "site"
DQ_ESM     <- "esm"
DQ_SCEN    <- "scenario"
DQ_WY      <- "water_year"
DQ_MEANQ   <- "mean_q_cfs_ma20"
DQ_ANOM    <- "anomaly_ma20"
DQ_DAYS    <- "days_since_wy_start_ma20"

# --- NCAR temperature envelope column names (unchanged) ---------------------
DT_SITE    <- "site"
DT_WY      <- "water_year"
DT_TLOW    <- "Temp_oC_low"
DT_THIGH   <- "Temp_oC_high"

# --- reach -> site mapping (unchanged) --------------------------------------
reach_site_map <- data.frame(
  reach = c("CLALO", "CLADR", "CLADR", "CLABE", "CLABE", "CLABE", "CLAPL"),
  Site  = c("DL",    "GR",    "BN",    "MS",    "BM",    "HU",    "FH"),
  stringsAsFactors = FALSE
)

VERBOSE <- FALSE
# ============================================================================


# ----------------------------------------------------------------------------
# small helpers (unchanged)
# ----------------------------------------------------------------------------
require_cols <- function(df, cols, what) {
  miss <- setdiff(cols, names(df))
  if (length(miss))
    stop(sprintf("%s is missing column(s): %s\n  Available: %s",
                 what, paste(miss, collapse = ", "),
                 paste(names(df), collapse = ", ")), call. = FALSE)
  invisible(TRUE)
}

assign_bracket <- function(s) {
  key <- toupper(gsub("[^A-Za-z0-9]", "", s))
  out <- rep(NA_character_, length(key))
  out[grepl("370", key)]                              <- "DROP"
  out[is.na(out) & (grepl("245", key) | grepl("RCP45", key))] <- "low"
  out[is.na(out) & (grepl("585", key) | grepl("RCP85", key))] <- "high"
  out
}

# ----------------------------------------------------------------------------
# apply_temp_treatment() -- unchanged from the validated version.
#
# Same function serves both tp and m1 -- pass whichever model and its own
# newdata; the treatment is applied to that model's own s(Temp_oC) smooth,
# so both models see internally-consistent Temp_oC behavior.
#
#   "none"    (V0)            : predict(..., type="response"), unchanged
#                                from the original script's direct calls.
#   "cap"     (V1, V2)        : Way A -- freeze the input at `cap` before
#                                predict(). Two-line change, no re-summing.
#   "decline" (V3mild/strong) : Way B -- predict(type="terms"), freeze the
#                                Temp_oC term at `threshold`, subtract
#                                slope * (Temp_oC - threshold) above it,
#                                re-sum, re-apply the link.
#
# Uses model$family$linkinv() rather than assuming identity link, so the
# decline branch is correct regardless of which family either GAM was fit
# with (a harmless no-op in the identity-link case).
# ----------------------------------------------------------------------------
apply_temp_treatment <- function(model, newdata, treatment) {
  
  if (treatment$type == "none") {
    return(as.numeric(predict(model, newdata = newdata, type = "response")))
  }
  
  if (treatment$type == "cap") {
    nd <- newdata
    nd$Temp_oC <- pmin(newdata$Temp_oC, treatment$cap)
    return(as.numeric(predict(model, newdata = nd, type = "response")))
  }
  
  if (treatment$type == "decline") {
    
    term_pred <- predict(model, newdata = newdata, type = "terms")
    intercept <- attr(term_pred, "constant")
    
    temp_col <- grep("Temp_oC", colnames(term_pred), value = TRUE)
    if (length(temp_col) != 1) {
      stop("apply_temp_treatment: expected exactly one Temp_oC smooth column, found: ",
           paste(temp_col, collapse = ", "))
    }
    
    temp_term <- term_pred[, temp_col]
    temp_val  <- newdata$Temp_oC
    
    nd_frozen <- newdata
    nd_frozen$Temp_oC <- treatment$threshold
    frozen_val <- predict(model, newdata = nd_frozen, type = "terms")[, temp_col]
    
    new_temp_term <- ifelse(
      temp_val > treatment$threshold,
      frozen_val - treatment$slope * (temp_val - treatment$threshold),
      temp_term
    )
    
    term_pred[, temp_col] <- new_temp_term
    linear_pred <- rowSums(term_pred) + intercept
    return(as.numeric(model$family$linkinv(linear_pred)))
  }
  
  stop("apply_temp_treatment: unknown treatment$type: ", treatment$type)
}


# ============================================================================
# 1. LOAD MODELS (unchanged, hoisted outside the loop -- identical every run)
# ============================================================================
m1 <- readRDS(PATH_M1)
tp <- readRDS(PATH_TP)

site_levels <- levels(m1$model$Site)
if (is.null(site_levels))
  stop("Could not recover Site factor levels from M1 model frame.", call. = FALSE)


# ============================================================================
# 2. OBSERVED CLIMATOLOGY + 2025 SEED (unchanged, hoisted outside the loop)
# ============================================================================
obs <- read.csv(PATH_OBS, stringsAsFactors = FALSE)
require_cols(obs, c(OBS_SITE, OBS_YEAR, OBS_LOGQ, OBS_ANOM, OBS_DSF), "ucfr_model_ready.csv")

if (OBS_CHLA %in% names(obs)) {
  obs$.logCHLa <- obs[[OBS_CHLA]]
} else if ("CHLa" %in% names(obs)) {
  obs$.logCHLa <- log10(obs$CHLa)
} else {
  stop("Need a logCHLa column (or CHLa to derive it) in ucfr_model_ready.csv.",
       call. = FALSE)
}

clim_src <- obs[stats::complete.cases(obs[, c(OBS_LOGQ, OBS_ANOM, OBS_DSF)]), ]
clim <- aggregate(
  cbind(obs_mean_logQ = clim_src[[OBS_LOGQ]],
        obs_mean_anom = clim_src[[OBS_ANOM]],
        obs_mean_dsf  = clim_src[[OBS_DSF]]) ~ clim_src[[OBS_SITE]],
  FUN = mean)
names(clim)[1] <- "Site"
clim_n <- as.data.frame(table(Site = clim_src[[OBS_SITE]]))

seed_src <- obs[obs[[OBS_YEAR]] == SEED_YEAR & is.finite(obs$.logCHLa), ]
seed_lag <- tapply(seed_src$.logCHLa, seed_src[[OBS_SITE]], max)
seed_lag <- seed_lag[is.finite(seed_lag)]

obs_range <- sapply(c(OBS_LOGQ, OBS_ANOM, OBS_DSF), function(cc)
  range(obs[[cc]], na.rm = TRUE))

obs_range_temp  <- range(obs$Temp_oC, na.rm = TRUE)
obs_range_logtp <- range(tp$model[[1]], na.rm = TRUE)


# ============================================================================
# 3. DISCHARGE ENVELOPE -> per (reach, member) baselines + future deltas
#    (unchanged, hoisted outside the loop)
# ============================================================================
disch <- read.csv(PATH_DISCH, stringsAsFactors = FALSE)
require_cols(disch, c(DQ_REACH, DQ_ESM, DQ_SCEN, DQ_WY, DQ_MEANQ, DQ_ANOM, DQ_DAYS),
             "ncar_discharge_envelope.csv")

disch$reach    <- disch[[DQ_REACH]]
disch$esm      <- disch[[DQ_ESM]]
disch$scenario <- disch[[DQ_SCEN]]
disch$year     <- disch[[DQ_WY]]
disch$mq       <- disch[[DQ_MEANQ]]
disch$anm      <- disch[[DQ_ANOM]]
disch$dys      <- disch[[DQ_DAYS]]

disch$bracket <- assign_bracket(disch$scenario)
scen_audit <- unique(disch[, c("scenario", "bracket")])
disch <- disch[!is.na(disch$bracket) & disch$bracket != "DROP", ]

disch$member <- paste(disch$esm, disch$scenario, sep = "|")
member_tab <- unique(disch[, c("member", "esm", "scenario", "bracket")])

base_src <- disch[disch$year >= BASELINE_START & disch$year <= BASELINE_END, ]
baselines <- aggregate(
  cbind(base_mq = base_src$mq, base_anm = base_src$anm, base_dys = base_src$dys) ~
    reach + member, data = base_src, FUN = mean, na.rm = TRUE)

fut <- disch[disch$year >= PROJ_START & disch$year <= PROJ_END,
             c("reach", "member", "esm", "scenario", "bracket", "year",
               "mq", "anm", "dys")]
fut <- merge(fut, baselines, by = c("reach", "member"), all.x = TRUE)


# ============================================================================
# 4. EXPAND reach -> site, apply deltas to site-specific observed climatology
#    (unchanged, hoisted outside the loop -- base_grid has everything that
#    does NOT depend on temp treatment: anomaly, logQ_obs_cfs, DSF, Temp_oC,
#    uid. logTP_mg_L and pred_logCHLa are treatment-dependent and get
#    (re)computed per variant inside the loop.)
# ============================================================================
base_grid <- merge(fut, reach_site_map, by = "reach")
base_grid <- merge(base_grid, clim, by = "Site", all.x = TRUE)

base_grid$anomaly            <- (base_grid$anm / base_grid$base_anm) * base_grid$obs_mean_anom
base_grid$logQ_obs_cfs       <- base_grid$obs_mean_logQ + log10(base_grid$mq / base_grid$base_mq)
base_grid$Days_Since_Freshet <- base_grid$obs_mean_dsf - (base_grid$dys - base_grid$base_dys)


# ============================================================================
# 5. TEMPERATURE -> Temp_oC by bracket (low/high trajectory)
#    (unchanged, hoisted outside the loop)
# ============================================================================
temp <- read.csv(PATH_TEMP, stringsAsFactors = FALSE)
require_cols(temp, c(DT_SITE, DT_WY, DT_TLOW, DT_THIGH),
             "ncar_temperature_envelope.csv")

temp_lk <- aggregate(
  cbind(Temp_oC_low = temp[[DT_TLOW]], Temp_oC_high = temp[[DT_THIGH]]) ~
    Site + year,
  data = data.frame(Site = temp[[DT_SITE]], year = temp[[DT_WY]],
                    Temp_oC_low = temp[[DT_TLOW]], Temp_oC_high = temp[[DT_THIGH]]),
  FUN = mean, na.rm = TRUE)

base_grid <- merge(base_grid, temp_lk, by = c("Site", "year"), all.x = TRUE)
base_grid$Temp_oC <- ifelse(base_grid$bracket == "low",
                            base_grid$Temp_oC_low, base_grid$Temp_oC_high)
base_grid$Site <- factor(base_grid$Site, levels = site_levels)
base_grid$uid  <- paste(base_grid$Site, base_grid$member, sep = "|")

units <- unique(base_grid[, c("Site", "member", "uid")])
units$seed <- seed_lag[as.character(units$Site)]
seed_by_uid <- setNames(units$seed, units$uid)


# ============================================================================
# CHANGED: loop over all variants. Sections 6-10 from the single-variant
# script run once per TEMP_VARIANT, starting from the shared base_grid.
# ============================================================================
scorecard_rows <- vector("list", length(VARIANTS_TO_RUN))
extrap_rows    <- vector("list", length(VARIANTS_TO_RUN))
eoc_rows       <- vector("list", length(VARIANTS_TO_RUN))

for (v in VARIANTS_TO_RUN) {
  
  treatment <- TREATMENT[[v]]
  message("=== running TEMP_VARIANT = ", v, " ===")
  
  grid <- base_grid   # fresh copy each variant; base_grid itself untouched
  
  OUT_SUMMARY <- paste0(OUT_DIR, "bloom_projections_", v, ".csv")
  OUT_MEMBERS <- paste0(OUT_DIR, "bloom_projections_members_", v, ".csv")
  
  # -- 6. PROJECTED logTP via TP submodel -----------------------------------
  grid$logTP_mg_L <- apply_temp_treatment(
    tp,
    grid[, c("anomaly", "logQ_obs_cfs", "Days_Since_Freshet", "Temp_oC", "Site")],
    treatment
  )
  
  # -- 7. RECURSION ENGINE -- year-step M1 with self-fed lag_y --------------
  grid$pred_logCHLa <- NA_real_
  lag_state <- seed_by_uid
  
  years <- PROJ_START:PROJ_END
  for (t in years) {
    idx <- which(grid$year == t)
    if (!length(idx)) next
    nd <- grid[idx, ]
    nd$lag_y <- lag_state[nd$uid]
    
    p <- apply_temp_treatment(
      m1,
      nd[, c("lag_y", "anomaly", "logQ_obs_cfs",
             "Days_Since_Freshet", "logTP_mg_L", "Temp_oC", "Site")],
      treatment
    )
    
    grid$pred_logCHLa[idx] <- p
    lag_state[nd$uid] <- p
    if (VERBOSE) message(sprintf("  year %d done (%d units)", t, length(idx)))
  }
  
  # -- 8. ENSEMBLE SUMMARY ---------------------------------------------------
  summ_fun <- function(x) {
    x <- x[is.finite(x)]
    c(median = if (length(x)) median(x) else NA_real_,
      p10    = if (length(x)) as.numeric(quantile(x, 0.10)) else NA_real_,
      p90    = if (length(x)) as.numeric(quantile(x, 0.90)) else NA_real_,
      n      = length(x))
  }
  agg <- aggregate(pred_logCHLa ~ Site + bracket + year, data = grid,
                   FUN = summ_fun, na.action = na.pass)
  summary_tab <- data.frame(agg[, c("Site", "bracket", "year")], agg$pred_logCHLa)
  names(summary_tab)[4:7] <- c("median_logCHLa", "p10_logCHLa", "p90_logCHLa", "n_members")
  summary_tab <- summary_tab[order(summary_tab$Site, summary_tab$bracket,
                                   summary_tab$year), ]
  
  # -- 9. WRITE PER-VARIANT OUTPUTS ------------------------------------------
  members_tab <- grid[order(grid$Site, grid$member, grid$year),
                      c("Site", "reach", "member", "esm", "scenario", "bracket",
                        "year", "anomaly", "logQ_obs_cfs",
                        "Days_Since_Freshet", "logTP_mg_L", "Temp_oC",
                        "pred_logCHLa")]
  
  members_tab$temp_extrapolated <- members_tab$Temp_oC > obs_range_temp[2]
  
  members_tab$lag_y <- NA_real_
  mu <- paste(members_tab$Site, members_tab$member, sep = "|")
  for (u in unique(mu)) {
    ri <- which(mu == u)
    ri <- ri[order(members_tab$year[ri])]
    lags <- c(seed_by_uid[[u]], members_tab$pred_logCHLa[ri][-length(ri)])
    members_tab$lag_y[ri] <- lags
  }
  
  ex_agg <- aggregate(temp_extrapolated ~ Site + bracket + year,
                      data = members_tab, FUN = mean)
  names(ex_agg)[4] <- "frac_temp_extrap"
  summary_tab <- merge(summary_tab, ex_agg,
                       by = c("Site", "bracket", "year"), all.x = TRUE)
  
  write.csv(summary_tab, OUT_SUMMARY, row.names = FALSE)
  write.csv(members_tab, OUT_MEMBERS, row.names = FALSE)
  
  # -- 10. PER-VARIANT ROWS FOR THE UNIFIED SCORECARD ------------------------
  scorecard_rows[[v]] <- data.frame(
    temp_variant        = v,
    baseline_window     = sprintf("%d-%d", BASELINE_START, BASELINE_END),
    projection_window   = sprintf("%d-%d", PROJ_START, PROJ_END),
    seed_year           = SEED_YEAR,
    n_members           = length(unique(member_tab$member)),
    n_low_members       = sum(member_tab$bracket == "low"),
    n_high_members      = sum(member_tab$bracket == "high"),
    n_sites             = length(unique(grid$Site)),
    n_member_site_units = nrow(units),
    grid_rows           = nrow(grid),
    n_na_logTP          = sum(is.na(grid$logTP_mg_L)),
    n_na_pred_logCHLa   = sum(is.na(grid$pred_logCHLa)),
    min_pred_logCHLa    = min(grid$pred_logCHLa, na.rm = TRUE),
    median_pred_logCHLa = median(grid$pred_logCHLa, na.rm = TRUE),
    max_pred_logCHLa    = max(grid$pred_logCHLa, na.rm = TRUE),
    n_temp_extrapolated = sum(members_tab$temp_extrapolated),
    frac_temp_extrapolated = mean(members_tab$temp_extrapolated),
    members_file        = OUT_MEMBERS,
    summary_file        = OUT_SUMMARY,
    stringsAsFactors = FALSE
  )
  
  extrap_rows[[v]] <- data.frame(
    temp_variant = v,
    predictor    = c("logQ_obs_cfs", "anomaly", "Days_Since_Freshet", "Temp_oC", "logTP_mg_L"),
    obs_min      = c(obs_range[1, OBS_LOGQ], obs_range[1, OBS_ANOM], obs_range[1, OBS_DSF],
                     obs_range_temp[1], obs_range_logtp[1]),
    obs_max      = c(obs_range[2, OBS_LOGQ], obs_range[2, OBS_ANOM], obs_range[2, OBS_DSF],
                     obs_range_temp[2], obs_range_logtp[2]),
    proj_min     = c(min(grid$logQ_obs_cfs, na.rm = TRUE),
                     min(grid$anomaly, na.rm = TRUE),
                     min(grid$Days_Since_Freshet, na.rm = TRUE),
                     min(grid$Temp_oC, na.rm = TRUE),
                     min(grid$logTP_mg_L, na.rm = TRUE)),
    proj_max     = c(max(grid$logQ_obs_cfs, na.rm = TRUE),
                     max(grid$anomaly, na.rm = TRUE),
                     max(grid$Days_Since_Freshet, na.rm = TRUE),
                     max(grid$Temp_oC, na.rm = TRUE),
                     max(grid$logTP_mg_L, na.rm = TRUE)),
    n_below      = c(sum(grid$logQ_obs_cfs < obs_range[1, OBS_LOGQ], na.rm = TRUE),
                     sum(grid$anomaly      < obs_range[1, OBS_ANOM], na.rm = TRUE),
                     sum(grid$Days_Since_Freshet < obs_range[1, OBS_DSF], na.rm = TRUE),
                     sum(grid$Temp_oC      < obs_range_temp[1], na.rm = TRUE),
                     sum(grid$logTP_mg_L   < obs_range_logtp[1], na.rm = TRUE)),
    n_above      = c(sum(grid$logQ_obs_cfs > obs_range[2, OBS_LOGQ], na.rm = TRUE),
                     sum(grid$anomaly      > obs_range[2, OBS_ANOM], na.rm = TRUE),
                     sum(grid$Days_Since_Freshet > obs_range[2, OBS_DSF], na.rm = TRUE),
                     sum(grid$Temp_oC      > obs_range_temp[2], na.rm = TRUE),
                     sum(grid$logTP_mg_L   > obs_range_logtp[2], na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  
  eoc_this <- summary_tab[summary_tab$year == max(summary_tab$year),
                          c("Site", "bracket", "median_logCHLa", "p10_logCHLa", "p90_logCHLa")]
  eoc_this$temp_variant <- v
  eoc_rows[[v]] <- eoc_this[, c("temp_variant", "Site", "bracket",
                                "median_logCHLa", "p10_logCHLa", "p90_logCHLa")]
  
  message("  wrote ", nrow(members_tab), " rows to ", OUT_MEMBERS)
}


# ============================================================================
# CHANGED: unified scorecard across all variants -- one row per variant,
# plus combined extrapolation and end-of-century tables. Printed once at
# the end instead of once per run, and written to CSV for the comparison
# script.
# ============================================================================
unified_scorecard <- do.call(rbind, scorecard_rows)
unified_extrap    <- do.call(rbind, extrap_rows)
unified_eoc       <- do.call(rbind, eoc_rows)

write.csv(unified_scorecard, OUT_SCORECARD, row.names = FALSE)
write.csv(unified_extrap, OUT_EXTRAP, row.names = FALSE)
write.csv(unified_eoc, OUT_EOC, row.names = FALSE)

cat("\n================ UNIFIED TEMP-SENSITIVITY SCORECARD ================\n")
cat("\n-- scenario -> bracket audit (shared across all variants) --\n")
print(scen_audit, row.names = FALSE)
cat("\n-- per-variant run summary --\n")
print(unified_scorecard[, c("temp_variant", "grid_rows", "n_na_pred_logCHLa",
                            "min_pred_logCHLa", "median_pred_logCHLa", "max_pred_logCHLa",
                            "n_temp_extrapolated", "frac_temp_extrapolated")],
      row.names = FALSE)
cat("\n-- Temp_oC extrapolation vs observed range, by variant --\n")
print(unified_extrap[unified_extrap$predictor == "Temp_oC", ], row.names = FALSE)
cat("\n-- end-of-century (", PROJ_END, ") median logCHLa, all variants --\n", sep = "")
print(unified_eoc[order(unified_eoc$Site, unified_eoc$bracket, unified_eoc$temp_variant), ],
      row.names = FALSE)
cat("\n-- outputs written --\n")
print(data.frame(file = c(OUT_SCORECARD, OUT_EXTRAP, OUT_EOC,
                          unified_scorecard$members_file, unified_scorecard$summary_file)),
      row.names = FALSE)
cat("\n======================================================================\n")


# ============================================================================
# FIDELITY CHECK -- V0 only, already confirmed in the single-variant run
# (byte-identical to bloom_projections_members.csv after column reorder).
# Re-run manually if you want to re-verify after this restructuring.
# ============================================================================
# original <- read.csv("2_incremental/bloom_projections_members.csv")
# new_v0   <- read.csv(paste0(OUT_DIR, "bloom_projections_members_V0.csv"))
# o <- original[do.call(order, original[c("Site","member","year")]), ]
# n <- new_v0[do.call(order, new_v0[c("Site","member","year")]), ]
# rownames(o) <- NULL; rownames(n) <- NULL
# n <- n[names(o)]
# identical(o, n)