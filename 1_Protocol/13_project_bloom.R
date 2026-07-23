# ============================================================================
# 13_project_bloom.R
# UCFR Cladophora Bloom Prediction Pipeline -- Step 7: Future Projections
#
# Runs the bloom model M1 forward recursively, 2026-2099, for 21 NCAR
# ESM-scenario members (SSP370 excluded) x 7 sites x 2 scenario brackets
# (low = RCP4.5/SSP245, high = RCP8.5/SSP585). One annual prediction per
# (site, member, year); the prediction becomes next year's lag_y.
#
# Inputs
#   3_models/bloom_model_M1.rds         M1 GAM (6 predictors + s(Site, bs="re"))
#   3_models/tp_submodel.rds            TP submodel (logTP on log10(ug/L) scale)
#   2_incremental/ucfr_model_ready.csv  observed data (climatology + 2025 seed)
#   2_incremental/ncar_discharge_envelope.csv   MA20 hydrology, per reach
#   2_incremental/ncar_temperature_envelope.csv stream temp, per site (low/high)
#
# Outputs
#   2_incremental/bloom_projections.csv          ensemble summary (median/p10/p90
#                                                 per site x bracket x year)
#   2_incremental/bloom_projections_members.csv  full member-level trajectories
#
# --------------------------------------------------------------------------
# PREDICTOR CONSTRUCTION (delta-on-observed; never substitute raw NCAR into M1)
#   lag_y        recursive. Seed = observed 2025 site-level annual max logCHLa;
#                thereafter prior-year predicted logCHLa (single event/yr => the
#                annual max is that prediction).
#   anomaly      multiplicative delta (anomaly is linear; within-NCAR ratio
#                cancels the annual-min baseflow offset):
#                  anomaly_fut = (future_anomaly_ma20 / base_anomaly_ma20)
#                                 * obs_mean_anomaly[site]
#   logQ_obs_cfs LOG-ADDITIVE delta (a 10% flow change adds ~0.041 log units,
#                it does NOT scale the ~3.4 logQ value by 1.1):
#                  logQ_fut = obs_mean_logQ[site]
#                             + log10(future_mean_q_ma20 / base_mean_q_ma20)
#   Days_Since_Freshet  additive delta, MINUS sign (NCAR days_since_wy_start is
#                timing from Oct 1 [earlier freshet = smaller]; observed DSF is
#                sample-date minus freshet-peak [earlier freshet = LARGER] --
#                inversely related):
#                  DSF_fut = obs_mean_DSF[site]
#                            - (future_days_ma20 - base_days_ma20)
#   logTP_mg_L   TP submodel, direct feed. Submodel response is log10(ug/L),
#                identical scale to M1's logTP_mg_L predictor -- no conversion.
#                Inputs: projected anomaly/logQ/DSF/Temp + Site.
#   Temp_oC      direct from temperature envelope (already stream-temp scale).
#                Low-bracket members use the low trajectory, high use high.
#
#   Baseline window: 1981-2010 (WMO 30-yr climatological normal), per member.
# ============================================================================

library(mgcv)


# ============================================================================
# CONFIGURATION -- edit here only
# ============================================================================

# --- file paths -------------------------------------------------------------
PATH_M1        <- "3_models/bloom_model_M1.rds"
PATH_TP        <- "3_models/tp_submodel.rds"
PATH_OBS       <- "2_incremental/ucfr_model_ready.csv"
PATH_DISCH     <- "2_incremental/ncar_discharge_envelope.csv"
PATH_TEMP      <- "2_incremental/ncar_temperature_envelope.csv"
OUT_SUMMARY    <- "2_incremental/bloom_projections.csv"
OUT_MEMBERS    <- "2_incremental/bloom_projections_members.csv"

# --- windows ----------------------------------------------------------------
BASELINE_START <- 1981
BASELINE_END   <- 2010
PROJ_START     <- 2026
PROJ_END       <- 2099
SEED_YEAR      <- 2025          # observed annual-max year that seeds lag_y

# --- observed-data column names (ucfr_model_ready.csv) ----------------------
OBS_SITE   <- "Site"
OBS_YEAR   <- "Year"
OBS_CHLA   <- "logCHLa"         # falls back to log10(CHLa) if column absent
OBS_LOGQ   <- "logQ_obs_cfs"
OBS_ANOM   <- "anomaly"
OBS_DSF    <- "Days_Since_Freshet"

# --- NCAR discharge envelope column names -----------------------------------
DQ_REACH   <- "site"            # reach code column (CLALO/CLADR/CLABE/CLAPL)
DQ_ESM     <- "esm"
DQ_SCEN    <- "scenario"
DQ_WY      <- "water_year"
DQ_MEANQ   <- "mean_q_cfs_ma20"
DQ_ANOM    <- "anomaly_ma20"
DQ_DAYS    <- "days_since_wy_start_ma20"

# --- NCAR temperature envelope column names ---------------------------------
DT_SITE    <- "site"            # actual site (DL..FH)
DT_WY      <- "water_year"
DT_TLOW    <- "Temp_oC_low"
DT_THIGH   <- "Temp_oC_high"

# --- reach -> site mapping (delta shared across sites in a reach) ------------
reach_site_map <- data.frame(
  reach = c("CLALO", "CLADR", "CLADR", "CLABE", "CLABE", "CLABE", "CLAPL"),
  Site  = c("DL",    "GR",    "BN",    "MS",    "BM",    "HU",    "FH"),
  stringsAsFactors = FALSE
)

VERBOSE <- FALSE                # progress messages during recursion
# ============================================================================


# ----------------------------------------------------------------------------
# small helpers
# ----------------------------------------------------------------------------
require_cols <- function(df, cols, what) {
  miss <- setdiff(cols, names(df))
  if (length(miss))
    stop(sprintf("%s is missing column(s): %s\n  Available: %s",
                 what, paste(miss, collapse = ", "),
                 paste(names(df), collapse = ", ")), call. = FALSE)
  invisible(TRUE)
}

# scenario -> bracket. Auditable rule; the scenario->bracket table is emitted
# in the scorecard so the mapping can be verified against the actual strings.
assign_bracket <- function(s) {
  key <- toupper(gsub("[^A-Za-z0-9]", "", s))
  out <- rep(NA_character_, length(key))
  out[grepl("370", key)]                              <- "DROP"
  out[is.na(out) & (grepl("245", key) | grepl("RCP45", key))] <- "low"
  out[is.na(out) & (grepl("585", key) | grepl("RCP85", key))] <- "high"
  out
}


# ============================================================================
# 1. LOAD MODELS
# ============================================================================
m1 <- readRDS(PATH_M1)
tp <- readRDS(PATH_TP)

# Site factor levels must match the fitted models so s(Site, bs="re") returns
# the per-site fitted intercept (RE retained: projection is at known sites).
site_levels <- levels(m1$model$Site)
if (is.null(site_levels))
  stop("Could not recover Site factor levels from M1 model frame.", call. = FALSE)


# ============================================================================
# 2. OBSERVED CLIMATOLOGY + 2025 SEED
# ============================================================================
obs <- read.csv(PATH_OBS, stringsAsFactors = FALSE)
require_cols(obs, c(OBS_SITE, OBS_YEAR, OBS_LOGQ, OBS_ANOM, OBS_DSF), "ucfr_model_ready.csv")

# response on logCHLa scale: use the column if present, else derive from CHLa
if (OBS_CHLA %in% names(obs)) {
  obs$.logCHLa <- obs[[OBS_CHLA]]
} else if ("CHLa" %in% names(obs)) {
  obs$.logCHLa <- log10(obs$CHLa)
} else {
  stop("Need a logCHLa column (or CHLa to derive it) in ucfr_model_ready.csv.",
       call. = FALSE)
}

# per-site observed climatology (means of the three hydro anchors).
# Means over complete cases; robust to M1's residual-based outlier removal,
# which removed ~25 rows on residuals, not predictors.
clim_src <- obs[stats::complete.cases(obs[, c(OBS_LOGQ, OBS_ANOM, OBS_DSF)]), ]
clim <- aggregate(
  cbind(obs_mean_logQ = clim_src[[OBS_LOGQ]],
        obs_mean_anom = clim_src[[OBS_ANOM]],
        obs_mean_dsf  = clim_src[[OBS_DSF]]) ~ clim_src[[OBS_SITE]],
  FUN = mean)
names(clim)[1] <- "Site"
clim_n <- as.data.frame(table(Site = clim_src[[OBS_SITE]]))

# seed lag_y = observed SEED_YEAR site-level annual max logCHLa
seed_src <- obs[obs[[OBS_YEAR]] == SEED_YEAR & is.finite(obs$.logCHLa), ]
seed_lag <- tapply(seed_src$.logCHLa, seed_src[[OBS_SITE]], max)
seed_lag <- seed_lag[is.finite(seed_lag)]      # drop sites with no SEED_YEAR data

# observed predictor ranges (for extrapolation flagging in the scorecard)
obs_range <- sapply(c(OBS_LOGQ, OBS_ANOM, OBS_DSF), function(cc)
  range(obs[[cc]], na.rm = TRUE))

# Temp_oC and logTP ranges -- computed here (not down in the SCORECARD section
# where they conceptually belong) because temp_extrapolated in WRITE OUTPUTS
# (Sec. 9) needs obs_range_temp before the scorecard section ever runs.
obs_range_temp  <- range(obs$Temp_oC, na.rm = TRUE)
obs_range_logtp <- range(tp$model[[1]], na.rm = TRUE)   # tp$model[[1]] = fitted response


# ============================================================================
# 3. DISCHARGE ENVELOPE -> per (reach, member) baselines + future deltas
# ============================================================================
disch <- read.csv(PATH_DISCH, stringsAsFactors = FALSE)
require_cols(disch, c(DQ_REACH, DQ_ESM, DQ_SCEN, DQ_WY, DQ_MEANQ, DQ_ANOM, DQ_DAYS),
             "ncar_discharge_envelope.csv")

# normalise key names
disch$reach    <- disch[[DQ_REACH]]
disch$esm      <- disch[[DQ_ESM]]
disch$scenario <- disch[[DQ_SCEN]]
disch$year     <- disch[[DQ_WY]]
disch$mq       <- disch[[DQ_MEANQ]]
disch$anm      <- disch[[DQ_ANOM]]
disch$dys      <- disch[[DQ_DAYS]]

# bracket each member; drop SSP370
disch$bracket <- assign_bracket(disch$scenario)
scen_audit <- unique(disch[, c("scenario", "bracket")])
disch <- disch[!is.na(disch$bracket) & disch$bracket != "DROP", ]

# member = unique (esm, scenario)
disch$member <- paste(disch$esm, disch$scenario, sep = "|")
member_tab <- unique(disch[, c("member", "esm", "scenario", "bracket")])

# per (reach, member) baseline = mean of MA20 metrics over 1981-2010
base_src <- disch[disch$year >= BASELINE_START & disch$year <= BASELINE_END, ]
baselines <- aggregate(
  cbind(base_mq = base_src$mq, base_anm = base_src$anm, base_dys = base_src$dys) ~
    reach + member, data = base_src, FUN = mean, na.rm = TRUE)

# future rows over the projection window, joined to their member baseline
fut <- disch[disch$year >= PROJ_START & disch$year <= PROJ_END,
             c("reach", "member", "esm", "scenario", "bracket", "year",
               "mq", "anm", "dys")]
fut <- merge(fut, baselines, by = c("reach", "member"), all.x = TRUE)


# ============================================================================
# 4. EXPAND reach -> site, apply deltas to site-specific observed climatology
# ============================================================================
grid <- merge(fut, reach_site_map, by = "reach")          # many-to-one expansion
grid <- merge(grid, clim, by = "Site", all.x = TRUE)

grid$anomaly            <- (grid$anm / grid$base_anm) * grid$obs_mean_anom
grid$logQ_obs_cfs       <- grid$obs_mean_logQ + log10(grid$mq / grid$base_mq)
grid$Days_Since_Freshet <- grid$obs_mean_dsf - (grid$dys - grid$base_dys)


# ============================================================================
# 5. TEMPERATURE -> Temp_oC by bracket (low/high trajectory)
# ============================================================================
temp <- read.csv(PATH_TEMP, stringsAsFactors = FALSE)
require_cols(temp, c(DT_SITE, DT_WY, DT_TLOW, DT_THIGH),
             "ncar_temperature_envelope.csv")

# one low/high pair per (site, year): the two scenario trajectories are
# member-invariant by construction, so de-dup (mean is a no-op if unique).
temp_lk <- aggregate(
  cbind(Temp_oC_low = temp[[DT_TLOW]], Temp_oC_high = temp[[DT_THIGH]]) ~
    Site + year,
  data = data.frame(Site = temp[[DT_SITE]], year = temp[[DT_WY]],
                    t_low = temp[[DT_TLOW]], t_high = temp[[DT_THIGH]],
                    Temp_oC_low = temp[[DT_TLOW]], Temp_oC_high = temp[[DT_THIGH]]),
  FUN = mean, na.rm = TRUE)

grid <- merge(grid, temp_lk, by = c("Site", "year"), all.x = TRUE)
grid$Temp_oC <- ifelse(grid$bracket == "low", grid$Temp_oC_low, grid$Temp_oC_high)


# ============================================================================
# 6. PROJECTED logTP via TP submodel (direct feed, no conversion)
# ============================================================================
grid$Site <- factor(grid$Site, levels = site_levels)
grid$logTP_mg_L <- as.numeric(
  predict(tp, newdata = grid[, c("anomaly", "logQ_obs_cfs",
                                 "Days_Since_Freshet", "Temp_oC", "Site")]))


# ============================================================================
# 7. RECURSION ENGINE -- year-step M1 with self-fed lag_y
# ----------------------------------------------------------------------------
# Years are sequential (lag depends on the prior year), but sites x members are
# independent, so each year-step is a single vectorised predict across all
# (site, member) units.
# ============================================================================
grid$uid          <- paste(grid$Site, grid$member, sep = "|")
grid$pred_logCHLa <- NA_real_

units <- unique(grid[, c("Site", "member", "uid")])
units$seed <- seed_lag[as.character(units$Site)]
lag_state <- setNames(units$seed, units$uid)

years <- PROJ_START:PROJ_END
for (t in years) {
  idx <- which(grid$year == t)
  if (!length(idx)) next
  nd <- grid[idx, ]
  nd$lag_y <- lag_state[nd$uid]
  
  p <- as.numeric(predict(
    m1, newdata = nd[, c("lag_y", "anomaly", "logQ_obs_cfs",
                         "Days_Since_Freshet", "logTP_mg_L", "Temp_oC", "Site")]))
  
  grid$pred_logCHLa[idx] <- p
  lag_state[nd$uid] <- p          # single event/yr => this prediction is the
  # annual max and seeds next year's lag_y
  if (VERBOSE) message(sprintf("year %d done (%d units)", t, length(idx)))
}


# ============================================================================
# 8. ENSEMBLE SUMMARY (median / p10 / p90 across members, per bracket)
# ============================================================================
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


# ============================================================================
# 9. WRITE OUTPUTS
# ============================================================================
members_tab <- grid[order(grid$Site, grid$member, grid$year),
                    c("Site", "reach", "member", "esm", "scenario", "bracket",
                      "year", "anomaly", "logQ_obs_cfs",
                      "Days_Since_Freshet", "logTP_mg_L", "Temp_oC",
                      "pred_logCHLa")]
members_tab$lag_y <- NA_real_

# Temperature beyond the M1 calibration range -> additive projection is
# extrapolating the Temp smooth's tail; flag for the thermal-stress caveat.
members_tab$temp_extrapolated <- members_tab$Temp_oC > obs_range_temp[2]

# lag_y above is the last assigned value during recursion; recompute the
# per-row lag for an honest member-level record (lag at year t = pred at t-1).
members_tab$lag_y <- NA_real_
seed_by_uid <- setNames(units$seed, units$uid)
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


# ============================================================================
# 10. SCORECARD
# ============================================================================

# obs_range_temp / obs_range_logtp: computed earlier in Sec. 2 (needed by
# Sec. 9's temp_extrapolated flag before this scorecard section ever runs)

extrap <- data.frame(
  predictor = c("logQ_obs_cfs", "anomaly", "Days_Since_Freshet", "Temp_oC", "logTP_mg_L"),
  obs_min   = c(obs_range[1, OBS_LOGQ], obs_range[1, OBS_ANOM], obs_range[1, OBS_DSF],
                obs_range_temp[1], obs_range_logtp[1]),
  obs_max   = c(obs_range[2, OBS_LOGQ], obs_range[2, OBS_ANOM], obs_range[2, OBS_DSF],
                obs_range_temp[2], obs_range_logtp[2]),
  proj_min  = c(min(grid$logQ_obs_cfs, na.rm = TRUE),
                min(grid$anomaly, na.rm = TRUE),
                min(grid$Days_Since_Freshet, na.rm = TRUE),
                min(grid$Temp_oC, na.rm = TRUE),
                min(grid$logTP_mg_L, na.rm = TRUE)),
  proj_max  = c(max(grid$logQ_obs_cfs, na.rm = TRUE),
                max(grid$anomaly, na.rm = TRUE),
                max(grid$Days_Since_Freshet, na.rm = TRUE),
                max(grid$Temp_oC, na.rm = TRUE),
                max(grid$logTP_mg_L, na.rm = TRUE)),
  n_below   = c(sum(grid$logQ_obs_cfs < obs_range[1, OBS_LOGQ], na.rm = TRUE),
                sum(grid$anomaly      < obs_range[1, OBS_ANOM], na.rm = TRUE),
                sum(grid$Days_Since_Freshet < obs_range[1, OBS_DSF], na.rm = TRUE),
                sum(grid$Temp_oC      < obs_range_temp[1], na.rm = TRUE),
                sum(grid$logTP_mg_L   < obs_range_logtp[1], na.rm = TRUE)),
  n_above   = c(sum(grid$logQ_obs_cfs > obs_range[2, OBS_LOGQ], na.rm = TRUE),
                sum(grid$anomaly      > obs_range[2, OBS_ANOM], na.rm = TRUE),
                sum(grid$Days_Since_Freshet > obs_range[2, OBS_DSF], na.rm = TRUE),
                sum(grid$Temp_oC      > obs_range_temp[2], na.rm = TRUE),
                sum(grid$logTP_mg_L   > obs_range_logtp[2], na.rm = TRUE)),
  stringsAsFactors = FALSE)


run_config <- data.frame(
  item = c("baseline_window", "projection_window", "seed_year",
           "n_members", "n_low_members", "n_high_members", "n_sites",
           "n_member_site_units", "grid_rows"),
  value = c(sprintf("%d-%d", BASELINE_START, BASELINE_END),
            sprintf("%d-%d", PROJ_START, PROJ_END),
            as.character(SEED_YEAR),
            length(unique(member_tab$member)),
            sum(member_tab$bracket == "low"),
            sum(member_tab$bracket == "high"),
            length(unique(grid$Site)),
            nrow(units),
            nrow(grid)),
  stringsAsFactors = FALSE)

seed_table <- data.frame(Site = names(seed_lag),
                         seed_lag_y = as.numeric(seed_lag),
                         row.names = NULL)
missing_seed <- setdiff(site_levels, names(seed_lag))

clim_table <- merge(clim, clim_n, by = "Site", all.x = TRUE)
names(clim_table)[names(clim_table) == "Freq"] <- "n_obs"

bracket_summary <- aggregate(pred_logCHLa ~ bracket, data = grid,
                             FUN = function(x) {
                               x <- x[is.finite(x)]
                               c(min = min(x), median = median(x), max = max(x), n = length(x))
                             }, na.action = na.pass)
bracket_summary <- data.frame(bracket = bracket_summary$bracket,
                              bracket_summary$pred_logCHLa)

na_summary <- data.frame(
  field = c("anomaly", "logQ_obs_cfs", "Days_Since_Freshet", "Temp_oC",
            "logTP_mg_L", "pred_logCHLa"),
  n_na  = c(sum(is.na(grid$anomaly)), sum(is.na(grid$logQ_obs_cfs)),
            sum(is.na(grid$Days_Since_Freshet)), sum(is.na(grid$Temp_oC)),
            sum(is.na(grid$logTP_mg_L)), sum(is.na(grid$pred_logCHLa))),
  stringsAsFactors = FALSE)

output_files <- data.frame(
  file = c(OUT_SUMMARY, OUT_MEMBERS),
  rows = c(nrow(summary_tab), nrow(members_tab)),
  stringsAsFactors = FALSE)

# spot check: end-of-century (last year) median by site and bracket
eoc <- summary_tab[summary_tab$year == max(summary_tab$year),
                   c("Site", "bracket", "median_logCHLa", "p10_logCHLa", "p90_logCHLa")]

cat("\n================ 13_project_bloom.R SCORECARD ================\n")
cat("\n-- run config --\n");                 print(run_config, row.names = FALSE)
cat("\n-- scenario -> bracket audit --\n");  print(scen_audit, row.names = FALSE)
cat("\n-- members per bracket --\n")
print(as.data.frame(table(bracket = member_tab$bracket)), row.names = FALSE)
cat("\n-- 2025 seed lag_y by site --\n");     print(seed_table, row.names = FALSE)
if (length(missing_seed))
  cat("   !! sites with no", SEED_YEAR, "seed:", paste(missing_seed, collapse = ", "), "\n")
cat("\n-- observed climatology (delta anchors) --\n"); print(clim_table, row.names = FALSE)
cat("\n-- extrapolation vs observed predictor range --\n"); print(extrap, row.names = FALSE)
cat("\n-- NA counts in projection grid --\n"); print(na_summary, row.names = FALSE)
cat("\n-- predicted logCHLa range by bracket --\n"); print(bracket_summary, row.names = FALSE)
cat("\n-- end-of-century (", PROJ_END, ") median logCHLa --\n", sep = "")
print(eoc, row.names = FALSE)
cat("\n-- outputs written --\n");             print(output_files, row.names = FALSE)
cat("\n=============================================================\n")



hot <- grid[grid$Temp_oC > obs_range_temp[2], ]
aggregate(year ~ Site, data = hot, FUN = function(y) c(min = min(y), max = max(y), n = length(y)))