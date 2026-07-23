# ============================================================================
# 07d_ncar_temperature_envelope.R
# UCFR Filamentous Algae Project
# Temperature envelope — production version
#
# Input:   2_incremental/ncar_discharge_envelope.csv  (for member/scenario list)
# Output:  2_incremental/ncar_temperature_envelope.csv
#
# Produces annual site-level stream temperature trajectories for each
# NCAR ESM-scenario member across the full 1952-2099 projection window,
# under two scenario brackets.
#
# ARCHITECTURE:
#   Two scenario brackets only (symmetric with discharge pipeline):
#     Low  (RCP4.5 / SSP245): NorthWEST trajectory
#     High (RCP8.5 / SSP585): NorthWEST trajectory + linearly growing offset
#   SSP370 excluded from main pipeline. Data available in
#   ncar_discharge_envelope.csv for supplementary addition if requested.
#
# TEMPERATURE SOURCE — NorthWEST portal (hardcoded):
#   Site-specific mean summer stream temperature at 3 anchor years:
#     2011 = observed baseline (NorthWEST historical)
#     2040 = mid-century projection
#     2080 = end-of-century projection
#   Rows for 2020, 2050, 2099 were empty in the portal output.
#   FH has a row-order anomaly in the portal (2040/2050 swapped) —
#   values are read by Year, not row position.
#
# INTERPOLATION:
#   Linear between 2011→2040→2080 for each site.
#   Held flat at 2080 value from 2080→2099.
#   Pre-2011 years (historical spin-up 1952–2010): held flat at 2011 value.
#   This is appropriate — NorthWEST does not provide a pre-baseline
#   trajectory and we do not extrapolate backward.
#
# RCP8.5 OFFSET (applied on top of NorthWEST trajectory):
#   0.0°C at 2011 (no divergence at baseline)
#   +0.5°C at 2040 (mid-century scenario spread)
#   +1.5°C at 2080 (end-century scenario spread)
#   Held flat at +1.5°C from 2080→2099.
#   Anchored to MCA 2017 statewide summer warming spread between
#   RCP4.5 and RCP8.5 (~0.8°C mid-century, ~2.3°C end-century air temp)
#   × 0.7 air-to-stream transfer ratio for snowmelt-dominated systems.
#   Pre-2011: offset = 0.0°C (no divergence in historical period).
#
# SCENARIO MAPPING for NCAR members:
#   CMIP5 RCP4.5  → low  trajectory
#   CMIP5 RCP8.5  → high trajectory
#   CMIP6 SSP245  → low  trajectory
#   CMIP6 SSP585  → high trajectory
#   CMIP6 SSP370  → EXCLUDED (flagged in output; supp candidate)
#
# All members within a scenario bracket share the same temperature
# trajectory. No GCM-by-GCM temperature pairing attempted.
#
# Methods sentence:
#   "Stream temperature projections were derived from NorthWEST
#    site-specific estimates (low-emissions trajectory) with a
#    high-emissions offset anchored to MCA 2017 statewide summer
#    warming spread and a 0.7 air-to-stream transfer ratio."
# ============================================================================

# ============================================================================
# 1. HARDCODED NorthWEST ANCHOR VALUES
# ============================================================================
# Source: NorthWEST portal, manually extracted.
# Mean summer stream temperature (°C) at 3 anchor years per site.
# FH note: portal output had 2040/2050 rows swapped; values confirmed
# by Year field, not row position.

northwest_raw <- data.frame(
  Site = c(
    "DL","DL","DL",
    "GR","GR","GR",
    "BN","BN","BN",
    "MS","MS","MS",
    "BM","BM","BM",
    "HU","HU","HU",
    "FH","FH","FH"
  ),
  Year = c(
    2011, 2040, 2080,
    2011, 2040, 2080,
    2011, 2040, 2080,
    2011, 2040, 2080,
    2011, 2040, 2080,
    2011, 2040, 2080,
    2011, 2040, 2080
  ),
  Temp_oC = c(
    15.85, 17.30, 18.36,   # DL
    16.67, 18.15, 19.23,   # GR
    14.59, 15.99, 17.01,   # BN
    18.03, 19.57, 20.69,   # MS
    18.74, 20.31, 21.45,   # BM
    20.24, 21.87, 23.05,   # HU
    18.78, 20.35, 21.49    # FH
  ),
  stringsAsFactors = FALSE
)

SITES <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")

# ============================================================================
# 2. RCP8.5 OFFSET ANCHORS
# ============================================================================
# Linearly interpolated between these three points; held flat outside range.
# Pre-2011: 0.0°C. Post-2080: +1.5°C.

OFFSET_YEARS  <- c(2011, 2040, 2080)
OFFSET_VALUES <- c(0.0,  0.5,  1.5)   # °C above RCP4.5 trajectory

# ============================================================================
# 3. SCENARIO MAPPING
# ============================================================================
# Maps each NCAR scenario string to a temperature bracket.
# SSP370 → NA (excluded from main pipeline).

scenario_bracket <- c(
  "rcp45"  = "low",
  "rcp85"  = "high",
  "ssp245" = "low",
  "ssp585" = "high",
  "ssp370" = NA_character_   # excluded; supp candidate
)

# ============================================================================
# 4. LOAD MEMBER LIST FROM DISCHARGE ENVELOPE
# ============================================================================
env_df  <- read.csv("2_incremental/ncar_discharge_envelope.csv",
                    stringsAsFactors = FALSE)
members <- unique(env_df[, c("esm", "scenario", "cmip")])
members <- members[order(members$scenario, members$esm), ]

# Tag each member with its temperature bracket
members$temp_bracket <- scenario_bracket[members$scenario]

n_excluded <- sum(is.na(members$temp_bracket))
members_active <- members[!is.na(members$temp_bracket), ]

cat("Members total:", nrow(members), "\n")
cat("Members excluded (SSP370):", n_excluded, "\n")
cat("Members in main pipeline:", nrow(members_active), "\n\n")

# ============================================================================
# 5. INTERPOLATION HELPERS
# ============================================================================

# Linear interpolation between anchor points; flat outside range
interp_anchors <- function(year, anchor_years, anchor_vals) {
  if (year <= anchor_years[1]) return(anchor_vals[1])
  if (year >= anchor_years[length(anchor_years)]) {
    return(anchor_vals[length(anchor_vals)])
  }
  # Find surrounding segment
  i <- max(which(anchor_years <= year))
  j <- i + 1
  frac <- (year - anchor_years[i]) / (anchor_years[j] - anchor_years[i])
  anchor_vals[i] + frac * (anchor_vals[j] - anchor_vals[i])
}

# Apply to a vector of years
interp_vec <- function(years, anchor_years, anchor_vals) {
  sapply(years, interp_anchors,
         anchor_years = anchor_years, anchor_vals = anchor_vals)
}

# ============================================================================
# 6. BUILD NorthWEST TRAJECTORY PER SITE
# ============================================================================
# For each site: extract 3 anchor values, interpolate to full year sequence.
# Pre-2011: flat at 2011 value (no backward extrapolation).
# Post-2080: flat at 2080 value.

YEARS_FULL <- 1952:2099

nw_trajectories <- vector("list", length(SITES))
names(nw_trajectories) <- SITES

for (s in SITES) {
  anchors <- northwest_raw[northwest_raw$Site == s, ]
  anchors <- anchors[order(anchors$Year), ]
  
  if (nrow(anchors) != 3 || any(is.na(anchors$Temp_oC))) {
    stop("NorthWEST data incomplete for site: ", s,
         " — check hardcoded values.")
  }
  
  nw_trajectories[[s]] <- interp_vec(
    YEARS_FULL,
    anchor_years = anchors$Year,
    anchor_vals  = anchors$Temp_oC
  )
}

# ============================================================================
# 7. BUILD RCP8.5 OFFSET TRAJECTORY (shared across all sites/members)
# ============================================================================
offset_vec <- interp_vec(YEARS_FULL, OFFSET_YEARS, OFFSET_VALUES)

# ============================================================================
# 8. ASSEMBLE OUTPUT: one row per site x member x year
# ============================================================================
out_list <- vector("list", nrow(members_active) * length(SITES))
k <- 0L

for (m in seq_len(nrow(members_active))) {
  esm_i      <- members_active$esm[m]
  scenario_i <- members_active$scenario[m]
  cmip_i     <- members_active$cmip[m]
  bracket_i  <- members_active$temp_bracket[m]
  
  for (s in SITES) {
    k <- k + 1L
    
    base_traj <- nw_trajectories[[s]]
    
    temp_low  <- base_traj
    temp_high <- base_traj + offset_vec
    
    # Assign trajectory based on bracket
    temp_proj <- if (bracket_i == "low") temp_low else temp_high
    
    out_list[[k]] <- data.frame(
      site          = s,
      esm           = esm_i,
      scenario      = scenario_i,
      cmip          = cmip_i,
      temp_bracket  = bracket_i,
      water_year    = YEARS_FULL,
      Temp_oC_low   = round(temp_low,  3),   # RCP4.5/SSP245 trajectory
      Temp_oC_high  = round(temp_high, 3),   # RCP8.5/SSP585 trajectory
      Temp_oC       = round(temp_proj, 3),   # this member's assigned trajectory
      stringsAsFactors = FALSE
    )
  }
}

temp_env <- do.call(rbind, out_list)
row.names(temp_env) <- NULL
temp_env <- temp_env[order(temp_env$site, temp_env$scenario,
                           temp_env$esm, temp_env$water_year), ]

# ============================================================================
# 9. SAVE OUTPUT
# ============================================================================
if (!dir.exists("2_incremental")) dir.create("2_incremental", recursive = TRUE)
write.csv(temp_env, "2_incremental/ncar_temperature_envelope.csv",
          row.names = FALSE)

# ============================================================================
# 10. SCORECARD
# ============================================================================
cat("============================================================\n")
cat("TEMPERATURE ENVELOPE SCORECARD\n")
cat("============================================================\n\n")

# --- 10a. Coverage ---
cov_sc <- data.frame(
  Metric = c(
    "Sites",
    "Members in main pipeline",
    "Members excluded (SSP370)",
    "Years per member x site",
    "Total rows",
    "Scenario brackets"
  ),
  Value = c(
    length(SITES),
    nrow(members_active),
    n_excluded,
    length(YEARS_FULL),
    nrow(temp_env),
    "low (RCP4.5/SSP245), high (RCP8.5/SSP585)"
  ),
  stringsAsFactors = FALSE
)
cat("--- Coverage ---\n")
print(cov_sc, row.names = FALSE)
cat("\n")

# --- 10b. NorthWEST anchor values (sanity check) ---
cat("--- NorthWEST anchor values (hardcoded) ---\n")
print(northwest_raw, row.names = FALSE)
cat("\n")

# --- 10c. Trajectory at key years per site ---
key_years <- c(2011, 2040, 2080, 2099)
traj_sc <- do.call(rbind, lapply(SITES, function(s) {
  do.call(rbind, lapply(key_years, function(yr) {
    low_val  <- temp_env$Temp_oC_low[temp_env$site == s &
                                       temp_env$water_year == yr][1]
    high_val <- temp_env$Temp_oC_high[temp_env$site == s &
                                        temp_env$water_year == yr][1]
    data.frame(
      Site      = s,
      Year      = yr,
      Low_RCP45 = round(low_val,  2),
      High_RCP85 = round(high_val, 2),
      Offset    = round(high_val - low_val, 2),
      stringsAsFactors = FALSE
    )
  }))
}))
cat("--- Projected temperatures at key years ---\n")
cat("(Low = RCP4.5/SSP245; High = RCP8.5/SSP585; Offset = scenario spread)\n\n")
print(traj_sc, row.names = FALSE)
cat("\n")

# --- 10d. Offset trajectory (same for all sites) ---
offset_sc <- data.frame(
  Year   = OFFSET_YEARS,
  Offset = OFFSET_VALUES,
  Basis  = c(
    "Baseline — no scenario divergence",
    "MCA 2017 ~0.8°C air × 0.7 transfer ratio",
    "MCA 2017 ~2.3°C air × 0.7 transfer ratio"
  ),
  stringsAsFactors = FALSE
)
cat("--- RCP8.5 offset anchors ---\n")
print(offset_sc, row.names = FALSE)
cat("\n")

# --- 10e. Excluded members ---
excluded_members <- members[is.na(members$temp_bracket),
                            c("esm", "scenario", "cmip")]
cat("--- Excluded members (SSP370) ---\n")
cat("(Available in ncar_discharge_envelope.csv for supplementary use)\n\n")
print(excluded_members, row.names = FALSE)
cat("\n")

# --- 10f. Output ---
cat("--- Output ---\n")
sc_files <- data.frame(
  File = "2_incremental/ncar_temperature_envelope.csv",
  Contents = paste0(
    nrow(temp_env), " rows. ",
    "Columns: site, esm, scenario, cmip, temp_bracket, water_year, ",
    "Temp_oC_low, Temp_oC_high, Temp_oC (member-assigned trajectory)."
  ),
  stringsAsFactors = FALSE
)
print(sc_files, row.names = FALSE)
cat("\n")
cat("============================================================\n")
cat("Done. Temp envelope ready for projection pipeline (Step 7).\n")
cat("Join to discharge envelope on site x esm x scenario x water_year.\n")
cat("============================================================\n")