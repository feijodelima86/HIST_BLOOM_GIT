# ============================================================================
# 03_process_vnrp.R
# UCFR Filamentous Algae Project
# Stage 2: Process raw VNRP data from WQP
#
# Input:    0_data/vnrp_raw.csv
# Output:   2_incremental/vnrp_processed.csv
#
# Steps:
#   1. Read raw data
#   2. Extract site codes from MonitoringLocationIdentifier, drop non-study sites
#   3. Lowercase UnitCode, convert ug/l to mg/l where applicable
#   4. Apply CharacteristicName alias lookup to canonical variable names
#   5. Parse dates, derive Year / Month / Day
#   6. Aggregate chemistry to year-month means; keep bio at visit-date level
#   7. Pivot wide separately; join chemistry onto bio (fan-out by month)
#   8. Calculate DIN = NH4 + NO3
#   9. Write output
#
# Notes:
#   - Aliases do not overlap temporally so no priority logic needed
#   - NAs retained and not dropped — handled at join stage
#   - Chl-a and AFDM kept in same wide table as chemistry
#   - D1 fix, revised (2026-06): chemistry is sampled biweekly through the
#     season and is intentionally month-aggregated -- this is unchanged.
#     CHLa/AFDM are each a single real sampling VISIT, and some Site-Year-
#     Months contain two real visits (confirmed: 21 cases, ~3-4 weeks
#     apart, each with genuinely different Days_Since_Freshet and
#     chemistry context) that earlier pipeline versions silently averaged
#     into one row. Bio is now aggregated at the visit-date level (no
#     cross-visit averaging), and chemistry's month-level mean is broadcast
#     onto however many bio visits fall in that month. This recovers the
#     ~21 previously-hidden observations as real, distinct rows.
# ============================================================================

library(readr)
library(dplyr)
library(tidyr)
library(zoo)
library(lubridate)

# ----------------------------------------------------------------------------
# 1. Read raw data
# ----------------------------------------------------------------------------

in_file  <- "0_data/vnrp_raw.csv"
out_dir  <- "2_incremental"
out_file <- file.path(out_dir, "vnrp_processed.csv")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("Reading raw VNRP data...\n")
vnrp <- read_csv(in_file,
                 col_types = cols(.default = col_character()),
                 show_col_types = FALSE)

cat("Rows:", nrow(vnrp), "| Cols:", ncol(vnrp), "\n\n")

# ----------------------------------------------------------------------------
# 2. Extract site codes, drop non-study sites
# ----------------------------------------------------------------------------

# MonitoringLocationIdentifier format: TSWQC_WQX-CFRPO-{number}
# Extract trailing number and map to canonical site code

site_lookup <- c(
  "9"    = "DL",
  "10"   = "GR",
  "12"   = "BN",
  "15.5" = "MS",
  "18"   = "BM",
  "22"   = "HU",
  "25"   = "FH"
)

vnrp$site_num <- sub(".*CFRPO-", "", vnrp$MonitoringLocationIdentifier)
vnrp$Site     <- site_lookup[vnrp$site_num]

# Drop rows that don't match a study site
n_before <- nrow(vnrp)
vnrp <- vnrp[!is.na(vnrp$Site), ]
n_after <- nrow(vnrp)

cat(sprintf("Site filtering: %d rows dropped (%d non-study site records)\n",
            n_before - n_after, n_before - n_after))
cat("Sites retained:", paste(sort(unique(vnrp$Site)), collapse = ", "), "\n\n")

# ----------------------------------------------------------------------------
# 3. Unit standardization
# ----------------------------------------------------------------------------

# Lowercase all unit codes first
vnrp$UnitCode <- tolower(vnrp$UnitCode)

# Variables requiring ug/l -> mg/l conversion
convert_vars <- c(
  "Ammonia",
  "Ammonia-nitrogen",
  "Inorganic nitrogen (nitrate and nitrite) ***retired***use Nitrate + Nitrite",
  "Nitrate + Nitrite",
  "Nutrient-nitrogen***retired***use TOTAL NITROGEN, MIXED FORMS with speciation AS N",
  "Orthophosphate",
  "Phosphate-phosphorus",
  "Total Phosphorus, mixed forms",
  "Total Nitrogen, mixed forms"
)

vnrp$ResultMeasureValue <- as.numeric(vnrp$ResultMeasureValue)

mask <- vnrp$CharacteristicName %in% convert_vars & !is.na(vnrp$UnitCode) & vnrp$UnitCode == "ug/l"
n_converted <- sum(mask, na.rm = TRUE)

vnrp$ResultMeasureValue[mask] <- vnrp$ResultMeasureValue[mask] / 1000
vnrp$UnitCode[mask]            <- "mg/l"

cat(sprintf("Unit conversion: %d values converted from ug/l to mg/l\n\n",
            n_converted))

# ----------------------------------------------------------------------------
# 4. CharacteristicName alias lookup -> canonical variable names
# ----------------------------------------------------------------------------

alias_lookup <- c(
  # Ammonia / NH4
  "Ammonia"                                                                                      = "NH4_mg_L",
  "Ammonia-nitrogen"                                                                             = "NH4_mg_L",
  
  # Nitrate + Nitrite / NO3
  "Inorganic nitrogen (nitrate and nitrite) ***retired***use Nitrate + Nitrite"                  = "NO3_mg_L",
  "Nitrate + Nitrite"                                                                            = "NO3_mg_L",
  
  # SRP / Orthophosphate
  "Orthophosphate"                                                                               = "SRP_mg_L",
  
  # TP
  "Phosphate-phosphorus"                                                                         = "TP_mg_L",
  "Total Phosphorus, mixed forms"                                                                = "TP_mg_L",
  
  # TN
  "Kjeldahl nitrogen"                                                                            = "TN_mg_L",
  "Nutrient-nitrogen***retired***use TOTAL NITROGEN, MIXED FORMS with speciation AS N"          = "TN_mg_L",
  "Total Nitrogen, mixed forms"                                                                  = "TN_mg_L",
  
  # Physical / other
  "pH"                                                                                           = "pH",
  "Specific conductance"                                                                         = "SPC",
  "Temperature, water"                                                                           = "Temp_oC",
  "Total dissolved solids"                                                                       = "TDS",
  "Turbidity"                                                                                    = "TURBIDITY",
  
  # Biological
  "Weight"                                                                                       = "AFDM",
  "Chlorophyll a, corrected for pheophytin"                                                      = "CHLa"
)

vnrp$Variable <- alias_lookup[vnrp$CharacteristicName]

# Report variables not in lookup (will be dropped at pivot)
unmapped <- unique(vnrp$CharacteristicName[is.na(vnrp$Variable)])
if (length(unmapped) > 0) {
  cat("Variables not mapped (will be dropped):\n")
  cat(paste(" ", unmapped, collapse = "\n"), "\n\n")
}

# Keep only mapped variables
vnrp <- vnrp[!is.na(vnrp$Variable), ]

cat("Canonical variables retained:",
    paste(sort(unique(vnrp$Variable)), collapse = ", "), "\n\n")

# ----------------------------------------------------------------------------
# 5. Parse dates, derive Year / Month / Day
# ----------------------------------------------------------------------------

vnrp$ActivityStartDate <- as.Date(vnrp$ActivityStartDate)
vnrp$Year              <- year(vnrp$ActivityStartDate)
vnrp$Month             <- month(vnrp$ActivityStartDate)
vnrp$Day               <- day(vnrp$ActivityStartDate)
vnrp$date_yearmon      <- as.yearmon(vnrp$ActivityStartDate)

# ----------------------------------------------------------------------------
# 6. Aggregate chemistry to year-month means; keep bio at visit-date level
# ----------------------------------------------------------------------------
# Chemistry (TP, TN, SPC, etc.) is sampled roughly biweekly through the
# season -- month-level mean aggregation is an intentional, ecologically
# reasonable smoothing choice and is UNCHANGED here.
#
# CHLa and AFDM are different: each record is a single real bloom-sampling
# VISIT. Some site-year-months contain two real visits (confirmed: 21
# distinct Site-Year-Month groups, ~3-4 weeks apart, each with its own
# chemistry/discharge context) that were previously being silently averaged
# into one row. Bio is therefore aggregated at the visit-date level instead
# of being collapsed to year-month, so each real visit becomes its own row.
# Chemistry's month-level mean is then broadcast onto however many bio
# visits fall in that month (1 or 2) when the two streams are joined below.

cat("Aggregating chemistry to year-month means...\n")

chem_vars_all <- c("NH4_mg_L", "NO3_mg_L", "SPC", "SRP_mg_L", "TDS",
                   "TN_mg_L", "TP_mg_L", "Temp_oC", "pH", "TURBIDITY")
bio_vars      <- c("CHLa", "AFDM")

chem_agg <- vnrp[vnrp$Variable %in% chem_vars_all, ] %>%
  group_by(Site, Year, Month, date_yearmon, Variable) %>%
  summarise(
    ResultMeasureValue = mean(ResultMeasureValue, na.rm = TRUE),
    n_obs              = sum(!is.na(ResultMeasureValue)),
    .groups            = "drop"
  )

cat("Chemistry rows after month-level aggregation:", nrow(chem_agg), "\n\n")

cat("Building visit-level biological dataset (no date collapsing)...\n")

bio_agg <- vnrp[vnrp$Variable %in% bio_vars, ] %>%
  group_by(Site, Year, Month, date_yearmon, ActivityStartDate, Variable) %>%
  summarise(
    ResultMeasureValue = mean(ResultMeasureValue, na.rm = TRUE),
    # mean() here is a no-op safeguard for same-day replicate readings,
    # NOT a cross-visit average -- ActivityStartDate is in the grouping key
    n_obs              = sum(!is.na(ResultMeasureValue)),
    .groups            = "drop"
  )

n_bio_visits   <- nrow(bio_agg[bio_agg$Variable == "CHLa", ])
n_bio_visits_y <- length(unique(paste(vnrp$Site[vnrp$Variable == "CHLa"],
                                      vnrp$Year[vnrp$Variable == "CHLa"],
                                      vnrp$Month[vnrp$Variable == "CHLa"])))

cat(sprintf(
  "  CHLa visit-level rows: %d  (vs. %d distinct Site-Year-Month groups --\n",
  n_bio_visits, n_bio_visits_y))
cat(sprintf(
  "   difference of %d reflects the recovered double-visit months)\n\n",
  n_bio_visits - n_bio_visits_y))

# ----------------------------------------------------------------------------
# 7. Pivot wide: full month-level panel, with bio visit-dates expanded
# ----------------------------------------------------------------------------
# IMPORTANT: vnrp_wide must remain a full Site-Year-Month panel (one row
# per month for ANY variable present, chemistry-only months included) --
# 04_lag_selection.R's sliding-window chemistry lookup depends on having
# every month's chemistry available, not just months with a bio visit.
# The only structural change from the original script is that months with
# TWO real bio visits get TWO rows instead of being collapsed to one.

cat("Pivoting chemistry to wide format...\n")

chem_wide <- chem_agg %>%
  select(Site, Year, Month, date_yearmon, Variable, ResultMeasureValue) %>%
  pivot_wider(names_from = Variable, values_from = ResultMeasureValue)

cat("Chemistry wide dimensions:", nrow(chem_wide), "rows x",
    ncol(chem_wide), "cols\n\n")

cat("Pivoting biological data to wide format (visit-level)...\n")

bio_wide <- bio_agg %>%
  select(Site, Year, Month, date_yearmon, ActivityStartDate,
         Variable, ResultMeasureValue) %>%
  pivot_wider(names_from = Variable, values_from = ResultMeasureValue)

names(bio_wide)[names(bio_wide) == "ActivityStartDate"] <- "Date_sample"

cat("Bio wide dimensions:", nrow(bio_wide), "rows x", ncol(bio_wide), "cols\n\n")

cat("Building full month-level panel (chemistry-only months retained,",
    "bio double-visit months expanded)...\n")

# Months with a bio visit: join chemistry onto each visit (fan-out when
# there are 2 visits in the same month -- both correctly get that month's
# chemistry mean, not duplicated chemistry rows from chem_wide itself)
bio_with_chem <- left_join(
  bio_wide, chem_wide,
  by = c("Site", "Year", "Month", "date_yearmon")
)

# Months with chemistry but NO bio visit at all: these must still appear
# in the final panel (04's sliding-window lookup needs them), with bio
# columns as NA, exactly as the original (pre-D1) script produced.
chem_only <- anti_join(
  chem_wide, bio_wide,
  by = c("Site", "Year", "Month", "date_yearmon")
)
chem_only$Date_sample <- as.Date(NA)
if ("CHLa" %in% names(bio_with_chem) && !"CHLa" %in% names(chem_only)) {
  chem_only$CHLa <- NA_real_
}
if ("AFDM" %in% names(bio_with_chem) && !"AFDM" %in% names(chem_only)) {
  chem_only$AFDM <- NA_real_
}

# Align column order/sets before binding
all_cols  <- union(names(bio_with_chem), names(chem_only))
bio_with_chem[setdiff(all_cols, names(bio_with_chem))] <- NA
chem_only[setdiff(all_cols, names(chem_only))]         <- NA

vnrp_wide <- rbind(bio_with_chem[, all_cols], chem_only[, all_cols])
vnrp_wide <- vnrp_wide[order(vnrp_wide$Site, vnrp_wide$Year, vnrp_wide$Month,
                             vnrp_wide$Date_sample), ]
rownames(vnrp_wide) <- NULL

# Sanity checks:
#   (a) every row with non-NA CHLa should have a non-NA Date_sample
#   (b) total rows should equal bio visit-rows + chemistry-only-month rows,
#       with NO loss of chemistry-only months relative to a plain month
#       panel (this is the regression the original draft of this fix had)
n_chla_no_date <- sum(!is.na(vnrp_wide$CHLa) & is.na(vnrp_wide$Date_sample))
if (n_chla_no_date > 0) {
  warning(sprintf(
    "%d rows have CHLa but no Date_sample -- check step 6/7 grouping keys.",
    n_chla_no_date))
}

n_chem_only_months <- nrow(chem_only)

cat("Final panel dimensions:", nrow(vnrp_wide), "rows x",
    ncol(vnrp_wide), "cols\n")
cat(sprintf("  (%d bio-visit rows + %d chemistry-only-month rows)\n\n",
            nrow(bio_with_chem), n_chem_only_months))

# ----------------------------------------------------------------------------
# 8. Calculate DIN = NH4 + NO3
# ----------------------------------------------------------------------------

if (all(c("NH4_mg_L", "NO3_mg_L") %in% names(vnrp_wide))) {
  vnrp_wide$DIN_mg_L <- vnrp_wide$NH4_mg_L + vnrp_wide$NO3_mg_L
  cat("DIN calculated as NH4 + NO3\n")
} else {
  warning("NH4_mg_L or NO3_mg_L missing — DIN not calculated")
}

# ----------------------------------------------------------------------------
# 9. Diagnostics
# ----------------------------------------------------------------------------

cat("\n--- Processing Summary ---\n")
cat("Sites:       ", paste(sort(unique(vnrp_wide$Site)), collapse = ", "), "\n")
cat("Year range:  ", min(vnrp_wide$Year, na.rm = TRUE), "to",
    max(vnrp_wide$Year, na.rm = TRUE), "\n")
cat("Columns:     ", paste(names(vnrp_wide), collapse = ", "), "\n\n")

cat("Non-NA counts per variable:\n")
target_vars <- c("CHLa", "AFDM", "TP_mg_L", "TN_mg_L", "SRP_mg_L",
                 "NH4_mg_L", "NO3_mg_L", "DIN_mg_L", "pH",
                 "Temp_oC", "SPC", "TDS", "TURBIDITY", "Date_sample")

for (v in target_vars) {
  if (v %in% names(vnrp_wide)) {
    cat(sprintf("  %-15s  n = %d\n", v, sum(!is.na(vnrp_wide[[v]]))))
  }
}

cat("\n--- D1 fix: visit-level bio recovery ---\n")
cat(sprintf("  CHLa visit-level rows:                   %d\n", n_bio_visits))
cat(sprintf("  Distinct Site-Year-Month groups (CHLa):  %d\n", n_bio_visits_y))
cat(sprintf("  Recovered double-visit rows:             %d\n",
            n_bio_visits - n_bio_visits_y))
cat(sprintf("  Chemistry-only months (no bio visit):    %d\n", n_chem_only_months))
cat(sprintf("  Rows with CHLa but missing Date_sample:  %d\n", n_chla_no_date))
cat(sprintf("  Date_sample range:                       %s to %s\n",
            min(vnrp_wide$Date_sample, na.rm = TRUE),
            max(vnrp_wide$Date_sample, na.rm = TRUE)))

# ----------------------------------------------------------------------------
# 10. Write output
# ----------------------------------------------------------------------------

write_csv(vnrp_wide, out_file)
cat("\nSaved to:", out_file, "\n")
cat("Done.\n")