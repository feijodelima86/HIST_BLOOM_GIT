# ============================================================================
# diag_freeze_one_driver.R
# UCFR Cladophora Bloom Prediction Pipeline -- QAQC diagnostic
# NOT manuscript-facing yet. Reserve analysis (named in project instructions
# as "freeze-one-driver re-runs"), now run per request. To be merged into the
# numbered pipeline (likely a 13d-style script) during the housekeeping pass --
# kept in diag_ for now, per project convention.
#
# ----------------------------------------------------------------------------
# QUESTION
# Section 4 / 13c_driver_decomposition.R attributes end-century CHLa change to
# each of M1's 6 predictors using predict(type="terms") -- an EXACT partition
# of the additive linear predictor, but one that assumes each smooth's
# contribution is independently interpretable. Two things that glosses over:
#
#   (1) Concurvity is 0.46 max pairwise (12_concurvity.R), not zero -- some
#       attributed "share" could be shared explanatory power reassigned to
#       whichever term predict(type="terms") happens to credit. That script's
#       own comments name three ecologically-plausible high-overlap
#       candidates (not confirmed as THE 0.46 pair -- cross-check against
#       concurvity_worst.csv if you want the actual pair):
#         anomaly x logQ, Days_Since_Freshet x logQ, lag_y x logTP.
#   (2) The terms-based number is a static early-window-vs-late-window mean
#       difference, not a counterfactual -- it never asks "what would have
#       happened if this driver hadn't moved."
#
# This script answers (2) directly via freeze-one-driver counterfactual
# re-projection, and uses agreement/disagreement with the terms-based ranking
# as an indirect check on (1): concurvity-driven misattribution should show
# up as rank disagreement between the two methods, especially for whichever
# pair above is the real 0.46.
#
# ----------------------------------------------------------------------------
# METHOD
# For each of the 6 M1 predictors, re-run the SAME recursive projection engine
# as 13_project_bloom.R (its Section 7), holding that one predictor fixed at
# its 2026 (first projection year) value for the rest of the horizon, while
# the other 5 evolve exactly as in the real projection. Comparing the
# resulting end-century trajectory to the real one isolates that driver's
# counterfactual contribution.
#
# Three defaults, stated explicitly -- override if you want something else:
#
#   D1. FREEZE VALUE = each (Site, member)'s own 2026 value, held constant
#       through 2099. Cheapest, no new baseline computation, reuses a value
#       already in the real trajectory.
#
#   D2. lag_y's freeze is NOT an external-driver freeze -- it disables the
#       recursive feedback channel itself. Every year is predicted using the
#       constant 2025 seed instead of the prior year's prediction. This
#       isolates how much projected change is carried by persistence /
#       compounding, independent of any climate driver. Kept in the same
#       6-driver comparison as 13c for direct comparability, but flagged here
#       and in the output table so it isn't read as "just another driver."
#
#   D3. logTP_mg_L is frozen DIRECTLY as an M1 input, exactly like the other
#       5 -- no cascading re-run of the TP submodel under frozen
#       anomaly/logQ/DSF/Temp. Matches how 13c/Section 4 already treats
#       logTP, as one of 6 parallel M1 inputs, so the two methods stay
#       comparable. A structural TP-submodel cascade would be a further
#       extension, not built here.
#
# Both scenarios (low, high) are run this time -- 13c/Section 4 only ever
# ran high. Terms-based comparison is recomputed for low too, so both methods
# are checked on the same footing.
#
# ----------------------------------------------------------------------------
# INPUTS (no new upstream computation -- reuses 13_project_bloom.R's own output)
#   2_incremental/bloom_projections_members.csv  (real per-member trajectories:
#                                                  anomaly/logQ/DSF/logTP/Temp/
#                                                  lag_y/pred_logCHLa/bracket)
#   3_models/bloom_model_M1.rds
#
# OUTPUTS
#   4_products/diagnostics/freeze_one_driver.csv        (full detail table)
#   4_products/diagnostics/freeze_one_driver_check.pdf  (terms % vs freeze %,
#                                                         one page per bracket)
#   console: scorecard
# ============================================================================

library(mgcv)

# ============================================================================
# CONFIGURATION -- edit here only
# ============================================================================
PATH_MEMBERS <- "2_incremental/bloom_projections_members.csv"
PATH_M1      <- "3_models/bloom_model_M1.rds"

OUT_DIR <- "4_products/diagnostics"
OUT_CSV <- file.path(OUT_DIR, "freeze_one_driver.csv")
OUT_PDF <- file.path(OUT_DIR, "freeze_one_driver_check.pdf")

SITE_ORDER <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")
DRIVERS    <- c("lag_y", "anomaly", "logQ_obs_cfs",
                "Days_Since_Freshet", "logTP_mg_L", "Temp_oC")
DRIVER_LAB <- c(lag_y = "lag", anomaly = "anomaly", logQ_obs_cfs = "logQ",
                Days_Since_Freshet = "DSF", logTP_mg_L = "logTP", Temp_oC = "Temp")

is_early <- function(y) y <= 2035     # mirrors 13c/Section 4's own filters
is_late  <- function(y) y >= 2089     # (not a fixed year list, in case the
#  projection horizon's end year changes)
# ============================================================================

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

require_cols <- function(df, cols, what) {
  miss <- setdiff(cols, names(df))
  if (length(miss))
    stop(sprintf("%s is missing column(s): %s\n  Available: %s",
                 what, paste(miss, collapse = ", "),
                 paste(names(df), collapse = ", ")), call. = FALSE)
}

# ============================================================================
# 1. LOAD
# ============================================================================
mem <- read.csv(PATH_MEMBERS, stringsAsFactors = FALSE)
require_cols(mem, c("Site", "member", "bracket", "year", "anomaly",
                    "logQ_obs_cfs", "Days_Since_Freshet", "logTP_mg_L",
                    "Temp_oC", "pred_logCHLa", "lag_y"),
             "bloom_projections_members.csv")

m1 <- readRDS(PATH_M1)
site_levels <- levels(m1$model$Site)
if (is.null(site_levels))
  stop("Could not recover Site factor levels from M1 model frame.", call. = FALSE)

mem$Site <- factor(mem$Site, levels = site_levels)
mem$uid  <- paste(mem$Site, mem$member, sep = "|")

proj_start <- min(mem$year)

# each uid's lag_y at the first projection year == the observed seed that
# fed it (13_project_bloom.R Section 7: lag at year 2026 == seed_lag[Site]).
# Re-derived here rather than re-read from ucfr_model_ready.csv, so this
# script has exactly one data-file dependency beyond the model object.
seed_row <- mem[mem$year == proj_start, c("uid", "lag_y")]
seed_row <- seed_row[!duplicated(seed_row$uid), ]
seed_lag <- setNames(seed_row$lag_y, seed_row$uid)

# ============================================================================
# 2. RECURSION ENGINE -- mirrors 13_project_bloom.R Section 7 exactly,
#    parameterised by which driver (if any) to freeze at its proj_start value.
# ============================================================================
run_recursion <- function(grid, freeze = NULL) {
  years <- sort(unique(grid$year))
  grid$pred_logCHLa <- NA_real_
  lag_state <- seed_lag[unique(grid$uid)]
  
  if (!is.null(freeze) && freeze != "lag_y") {
    anchor <- grid[grid$year == proj_start, c("uid", freeze)]
    anchor <- anchor[!duplicated(anchor$uid), ]
    fixed  <- setNames(anchor[[freeze]], anchor$uid)
    grid[[freeze]] <- fixed[grid$uid]
  }
  
  for (t in years) {
    idx <- which(grid$year == t)
    nd  <- grid[idx, ]
    nd$lag_y <- if (!is.null(freeze) && freeze == "lag_y") {
      seed_lag[nd$uid]                # never updates -- recursion disabled
    } else {
      lag_state[nd$uid]
    }
    p <- as.numeric(predict(m1, newdata = nd[, c("lag_y", "anomaly",
                                                 "logQ_obs_cfs", "Days_Since_Freshet", "logTP_mg_L",
                                                 "Temp_oC", "Site")]))
    grid$pred_logCHLa[idx] <- p
    lag_state[nd$uid] <- p
  }
  grid$pred_logCHLa
}

# ============================================================================
# 3. RUN: sanity-check re-run (no freeze) + 6 frozen counterfactuals
# ============================================================================
base_grid <- mem[, c("Site", "uid", "bracket", "year", "anomaly",
                     "logQ_obs_cfs", "Days_Since_Freshet", "logTP_mg_L", "Temp_oC")]

check_pred  <- run_recursion(base_grid, freeze = NULL)
recheck_err <- max(abs(check_pred - mem$pred_logCHLa), na.rm = TRUE)

frozen <- list()
for (d in DRIVERS) {
  message("Freezing ", d, " ...")
  frozen[[d]] <- run_recursion(base_grid, freeze = d)
}

# ============================================================================
# 4. EARLY/LATE WINDOW MEANS, FREEZE-BASED CONTRIBUTION
# ============================================================================
window_mean <- function(pred, keep_fun) {
  d <- data.frame(Site = base_grid$Site, bracket = base_grid$bracket,
                  year = base_grid$year, pred = pred)
  d <- d[keep_fun(d$year), ]
  ag <- aggregate(pred ~ Site + bracket, data = d, FUN = mean, na.rm = TRUE)
  names(ag)[3] <- "mean_logCHLa"
  ag
}

actual_early <- window_mean(mem$pred_logCHLa, is_early)
actual_late  <- window_mean(mem$pred_logCHLa, is_late)
actual <- merge(actual_early, actual_late, by = c("Site", "bracket"),
                suffixes = c("_early", "_late"))
actual$actual_delta <- actual$mean_logCHLa_late - actual$mean_logCHLa_early

freeze_rows <- list()
for (d in DRIVERS) {
  fe <- window_mean(frozen[[d]], is_early)
  fl <- window_mean(frozen[[d]], is_late)
  f  <- merge(fe, fl, by = c("Site", "bracket"), suffixes = c("_early", "_late"))
  f$driver       <- d
  f$frozen_delta <- f$mean_logCHLa_late - f$mean_logCHLa_early
  freeze_rows[[d]] <- f[, c("Site", "bracket", "driver",
                            "mean_logCHLa_early", "mean_logCHLa_late", "frozen_delta")]
}
freeze_tab <- do.call(rbind, freeze_rows)
names(freeze_tab)[names(freeze_tab) == "mean_logCHLa_early"] <- "frozen_early"
names(freeze_tab)[names(freeze_tab) == "mean_logCHLa_late"]  <- "frozen_late"

out <- merge(freeze_tab, actual[, c("Site", "bracket", "actual_delta")],
             by = c("Site", "bracket"))
out$contribution <- out$actual_delta - out$frozen_delta   # change lost if D frozen
out <- out[order(out$Site, out$bracket, out$driver), ]
out$pct_share_freeze <- with(out, ave(abs(contribution), Site, bracket,
                                      FUN = function(x) 100 * x / sum(x)))

# ============================================================================
# 5. TERMS-BASED COMPARISON (13c / Section 4 method), BOTH scenarios this time
# ============================================================================
tt <- predict(m1, newdata = mem[, c("lag_y", "anomaly", "logQ_obs_cfs",
                                    "Days_Since_Freshet", "logTP_mg_L",
                                    "Temp_oC", "Site")], type = "terms")
term_cols <- c("s(lag_y)", "s(anomaly)", "s(logQ_obs_cfs)",
               "s(Days_Since_Freshet)", "s(logTP_mg_L)", "s(Temp_oC)")
tt <- tt[, term_cols]
colnames(tt) <- DRIVERS

tg <- data.frame(Site = mem$Site, bracket = mem$bracket, year = mem$year, tt)
term_rows <- list()
for (b in c("low", "high")) {
  gb <- tg[tg$bracket == b, ]
  te <- gb[is_early(gb$year), ]
  tl <- gb[is_late(gb$year), ]
  for (st in SITE_ORDER) {
    delta <- colMeans(tl[tl$Site == st, DRIVERS]) - colMeans(te[te$Site == st, DRIVERS])
    term_rows[[paste(st, b)]] <- data.frame(Site = st, bracket = b,
                                            driver = DRIVERS, term_delta = delta)
  }
}
term_tab <- do.call(rbind, term_rows)
term_tab$pct_share_terms <- with(term_tab, ave(abs(term_delta), Site, bracket,
                                               FUN = function(x) 100 * x / sum(x)))

out <- merge(out, term_tab[, c("Site", "bracket", "driver",
                               "term_delta", "pct_share_terms")],
             by = c("Site", "bracket", "driver"))
out$rank_freeze <- ave(-abs(out$contribution), out$Site, out$bracket, FUN = rank)
out$rank_terms  <- ave(-abs(out$term_delta),   out$Site, out$bracket, FUN = rank)
out$rank_agree  <- out$rank_freeze == out$rank_terms
out <- out[order(out$Site, out$bracket, out$driver), ]

write.csv(out, OUT_CSV, row.names = FALSE)
cat("Wrote", OUT_CSV, "\n")

# ============================================================================
# 6. FIGURE -- terms % vs freeze %, one page per bracket, 7-panel + legend
# ============================================================================
plot_bracket <- function(b, title) {
  d <- out[out$bracket == b, ]
  par(mfrow = c(4, 2), mar = c(5, 4, 2.2, 1), oma = c(1, 1, 2.4, 1))
  for (st in SITE_ORDER) {
    ds <- d[d$Site == st, ]
    ds <- ds[match(DRIVERS, ds$driver), ]
    barplot(rbind(ds$pct_share_terms, ds$pct_share_freeze),
            beside = TRUE, col = c("grey60", "darkorange"),
            names.arg = DRIVER_LAB[ds$driver], las = 2,
            ylim = c(0, max(c(ds$pct_share_terms, ds$pct_share_freeze), 100)),
            ylab = "% of |end-century change|", main = st)
  }
  plot.new()
  legend("center", legend = c("Terms-based (13c)", "Freeze-based (counterfactual)"),
         fill = c("grey60", "darkorange"), bty = "n", cex = 1.0)
  mtext(title, outer = TRUE, cex = 1.0, font = 2, line = 0.4)
}

pdf(OUT_PDF, width = 9, height = 11, family = "Helvetica")
plot_bracket("high", "Driver attribution: terms-based vs freeze-based (high scenario)")
plot_bracket("low",  "Driver attribution: terms-based vs freeze-based (low scenario)")
dev.off()
cat("Wrote", OUT_PDF, "\n")

# ============================================================================
# SCORECARD
# ============================================================================
cat("\n============================================================\n")
cat("QAQC -- FREEZE-ONE-DRIVER COUNTERFACTUAL CHECK\n")
cat("Not manuscript-facing. Reserve analysis, now run per request.\n")
cat("============================================================\n\n")

cat("--- Recursion fidelity check (no-freeze re-run vs real pipeline output) ---\n")
cat("  max abs diff =", round(recheck_err, 6),
    "\n  (should be ~0; a larger value means this script's recursion doesn't\n",
    "  match 13_project_bloom.R and nothing below should be trusted yet)\n\n")

cat("--- Rank agreement: does the top driver match between methods? ---\n")
disagree <- out[out$rank_freeze == 1 & out$rank_terms != 1,
                c("Site", "bracket", "driver", "pct_share_terms", "pct_share_freeze")]
if (nrow(disagree)) {
  cat("  Sites/scenarios where the top-ranked driver DIFFERS between methods:\n")
  print(disagree, row.names = FALSE)
} else {
  cat("  None -- top driver agrees in every Site x scenario.\n")
}

cat("\n--- Full comparison table (first 12 rows) ---\n")
print(head(out[, c("Site", "bracket", "driver", "pct_share_terms",
                   "pct_share_freeze", "rank_agree")], 12), row.names = FALSE)

cat("\n--- Candidate concurvity pairs to cross-check disagreements against ---\n")
cat("  (from 12_concurvity.R's own comments -- not confirmed as THE 0.46 pair)\n")
cat("  anomaly x logQ | Days_Since_Freshet x logQ | lag_y x logTP\n")

cat("\n--- Output files ---\n")
print(data.frame(file = c(OUT_CSV, OUT_PDF)), row.names = FALSE)
cat("\n============================================================\n")
cat("Done.\n")
cat("============================================================\n")