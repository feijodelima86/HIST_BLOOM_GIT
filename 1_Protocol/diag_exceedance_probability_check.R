# ============================================================================
# diag_exceedance_probability_check.R
# UCFR Cladophora Bloom Prediction Pipeline -- QAQC diagnostic
# NOT manuscript-facing yet. Purpose: find out, before submission, whether
# P(CHLa > 150 mg/m^2) survives honest error propagation or collapses the
# low-vs-high bracket contrast -- i.e. whether the reframing is safe to ship.
#
# METHOD ("C+" -- not the parked TP Monte Carlo / full simulation build):
#
#   Outer draws -- the NCAR ensemble members already computed in
#   13_project_bloom.R (bloom_projections_members.csv). Scenario/driver
#   uncertainty already in the pipeline; nothing new here.
#
#   Inner draws -- bootstrap resample (with replacement) from the row-level
#   signed Error = Observed - Predicted_rec values already computed by
#   11_temporal_validation.R's Scheme D (recursive_val_predictions.csv).
#
# Why this addresses both open conditions from the audit, essentially for
# free:
#   Condition 1 (carry real predictive uncertainty, not just ensemble
#     spread) -- satisfied because the bootstrap pool is real out-of-sample
#     recursive predictive error, not a model-internal se.fit.
#   Condition 2 (target the right thing -- annual max, not conditional
#     mean) -- satisfied automatically: Scheme D's Observed column IS the
#     site-year annual max (see 11_temporal_validation.R, mdat_ann /
#     which.max). Bootstrapping those errors onto a future point prediction
#     produces draws that are, by construction, on the annual-max scale.
#     No separate max-distribution model needed.
#   No normal/parametric assumption is imposed on the residual -- raw signed
#   errors are resampled directly, preserving whatever shape they have.
#
# This script also reports the NAIVE probability (ensemble spread only, no
# residual bootstrap -- i.e. what a plain reframing without Condition 1
# would report) alongside the CORRECTED one, so the effect of carrying real
# uncertainty is directly visible, not just asserted.
#
# Inputs
#   2_incremental/bloom_projections_members.csv   (13_project_bloom.R)
#   2_incremental/recursive_val_predictions.csv   (11_temporal_validation.R, Scheme D)
#
# Outputs
#   4_products/diagnostics/exceedance_probability.csv
#   4_products/diagnostics/exceedance_probability_check.pdf
#   console: scorecard
#
# Known limitation (state if this graduates to the manuscript): the residual
# pool is drawn from validation years within the observed driver envelope.
# Projection years flagged temp_extrapolated (see 13_project_bloom.R) extend
# beyond that envelope, so the corrected probability in those years is a
# floor, not a full accounting -- flagged, not fixed, per project convention
# for the HU thermal-stress divergence.
# ============================================================================

# ============================================================================
# CONFIGURATION -- edit here only
# ============================================================================
PATH_MEMBERS <- "2_incremental/bloom_projections_members.csv"
PATH_REC     <- "2_incremental/recursive_val_predictions.csv"

OUT_DIR <- "4_products/diagnostics"
OUT_CSV <- file.path(OUT_DIR, "exceedance_probability.csv")
OUT_PDF <- file.path(OUT_DIR, "exceedance_probability_check.pdf")

THRESHOLD_CHLA_MG_M2 <- 150             # Suplee et al. 2009 JAWRA, recreational
# nuisance threshold (mg Chl a / m^2)
THRESHOLD_LOG10 <- log10(THRESHOLD_CHLA_MG_M2)

B_DRAWS_PER_MEMBER <- 2000               # residual bootstrap draws per member
POOL_RESIDUALS_ACROSS_SITES <- TRUE     # TRUE: ~90-row pooled pool (stabler,
# primary). FALSE: per-site pool
# (~10-13 rows/site, noisier tails --
# use as a sensitivity check only).
SEED <- 20260701

BENCHMARK_YEAR <- NULL                  # NULL = use max(year) present, for
# the headline comparison page

site_order <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")
col_low  <- "#377EB8"
col_high <- "#E41A1C"
# ============================================================================

set.seed(SEED)
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
rec <- read.csv(PATH_REC,     stringsAsFactors = FALSE)

require_cols(mem, c("Site", "bracket", "year", "pred_logCHLa"),
             "bloom_projections_members.csv")
require_cols(rec, c("Site", "Error"), "recursive_val_predictions.csv")

if (!"temp_extrapolated" %in% names(mem)) {
  warning("temp_extrapolated column not found in members file -- ",
          "extrapolation flag will be NA in outputs.")
  mem$temp_extrapolated <- NA
}

# ============================================================================
# 2. RESIDUAL BOOTSTRAP POOLS
# ============================================================================
resid_all <- rec$Error[is.finite(rec$Error)]
if (length(resid_all) < 20)
  warning("Pooled residual pool has only ", length(resid_all),
          " values -- bootstrap tails will be unstable.")

resid_by_site <- split(rec$Error[is.finite(rec$Error)],
                       rec$Site[is.finite(rec$Error)])

get_pool <- function(site) {
  if (POOL_RESIDUALS_ACROSS_SITES) return(resid_all)
  p <- resid_by_site[[site]]
  if (is.null(p) || length(p) < 5) {
    warning("Site ", site, " has < 5 residuals -- falling back to pooled.")
    return(resid_all)
  }
  p
}

# ============================================================================
# 3. EXCEEDANCE PROBABILITY -- per Site x bracket x year
# ============================================================================
groups <- unique(mem[, c("Site", "bracket", "year")])
groups <- groups[order(groups$Site, groups$bracket, groups$year), ]

n_g <- nrow(groups)
out <- data.frame(
  Site = groups$Site, bracket = groups$bracket, year = groups$year,
  n_members = NA_integer_, n_draws = NA_integer_,
  P_exceed_naive     = NA_real_,   # ensemble spread only -- no residual bootstrap
  P_exceed_corrected = NA_real_,   # ensemble + residual bootstrap
  frac_temp_extrap   = NA_real_
)

for (i in seq_len(n_g)) {
  s <- groups$Site[i]; b <- groups$bracket[i]; y <- groups$year[i]
  d <- mem[mem$Site == s & mem$bracket == b & mem$year == y, ]
  if (!nrow(d)) next
  
  preds <- d$pred_logCHLa[is.finite(d$pred_logCHLa)]
  n_mem <- length(preds)
  if (!n_mem) next
  
  # naive: member spread only -- what a plain reframing without Condition 1 reports
  out$P_exceed_naive[i] <- mean(preds > THRESHOLD_LOG10)
  
  # corrected: each member x B bootstrap-resampled historical recursive errors
  pool <- get_pool(s)
  resid_draws <- matrix(sample(pool, n_mem * B_DRAWS_PER_MEMBER, replace = TRUE),
                        nrow = n_mem, ncol = B_DRAWS_PER_MEMBER)
  composite <- preds + resid_draws     # recycles preds down each column (by row)
  
  out$n_members[i] <- n_mem
  out$n_draws[i]    <- length(composite)
  out$P_exceed_corrected[i] <- mean(composite > THRESHOLD_LOG10)
  
  if (!all(is.na(d$temp_extrapolated)))
    out$frac_temp_extrap[i] <- mean(d$temp_extrapolated, na.rm = TRUE)
}

write.csv(out, OUT_CSV, row.names = FALSE)

bench_year <- if (is.null(BENCHMARK_YEAR)) max(out$year, na.rm = TRUE) else BENCHMARK_YEAR
bench <- out[out$year == bench_year, ]

# ============================================================================
# 4. PLOTS
# ============================================================================
pdf(OUT_PDF, width = 10, height = 11)

# --- Page 1: headline -- naive vs corrected, by site, at benchmark year ---
par(mfrow = c(1, 1), mar = c(5, 4.5, 3, 1))
bench_ord <- bench[order(match(bench$Site, site_order), bench$bracket), ]
x <- seq_along(site_order)
plot(NA, xlim = c(0.5, length(site_order) + 0.5), ylim = c(0, 1),
     xaxt = "n", xlab = "", ylab = paste0("P(CHLa > ", THRESHOLD_CHLA_MG_M2,
                                          " mg/m2)  --  year ", bench_year))
axis(1, at = x, labels = site_order)
abline(h = seq(0, 1, 0.25), col = "grey85", lty = 3)

off <- 0.15
for (k in seq_along(site_order)) {
  st <- site_order[k]
  rl <- bench_ord[bench_ord$Site == st & bench_ord$bracket == "low", ]
  rh <- bench_ord[bench_ord$Site == st & bench_ord$bracket == "high", ]
  if (nrow(rl)) {
    segments(x[k]-off, rl$P_exceed_naive, x[k]-off, rl$P_exceed_corrected,
             col = col_low, lwd = 1)
    points(x[k]-off, rl$P_exceed_naive,     pch = 1,  col = col_low, cex = 1.3)
    points(x[k]-off, rl$P_exceed_corrected, pch = 16, col = col_low, cex = 1.3)
  }
  if (nrow(rh)) {
    segments(x[k]+off, rh$P_exceed_naive, x[k]+off, rh$P_exceed_corrected,
             col = col_high, lwd = 1)
    points(x[k]+off, rh$P_exceed_naive,     pch = 1,  col = col_high, cex = 1.3)
    points(x[k]+off, rh$P_exceed_corrected, pch = 16, col = col_high, cex = 1.3)
  }
}
legend("topleft",
       legend = c("Low, naive (ensemble only)", "Low, corrected (+ residual bootstrap)",
                  "High, naive (ensemble only)", "High, corrected (+ residual bootstrap)"),
       col = c(col_low, col_low, col_high, col_high), pch = c(1, 16, 1, 16),
       bty = "n", cex = 0.85)
mtext("open = naive (ensemble spread only) | filled = corrected (+ historical predictive error)",
      side = 3, line = 0.3, cex = 0.8, col = "grey30")

# --- Page 2: residual bootstrap pool itself (transparency check) ---
par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3, 1))
hist(resid_all, breaks = 20, col = "grey70", border = "white",
     xlab = "Historical recursive error, Observed annual max - Predicted (log10 CHLa)",
     main = "Residual bootstrap pool (pooled across sites)")
abline(v = 0, lty = 2, col = "grey30")
abline(v = mean(resid_all), lty = 2, col = "darkorange", lwd = 2)
legend("topright",
       legend = c(paste0("n = ", length(resid_all)),
                  paste0("mean = ", round(mean(resid_all), 4)),
                  paste0("SD = ",   round(sd(resid_all), 4))),
       bty = "n", cex = 0.9)

# --- Pages 3+: per-site corrected (thick) vs naive (thin, dashed) over time ---
par(mfrow = c(4, 2), mar = c(3.2, 3.6, 2.2, 0.8),
    mgp = c(2.1, 0.6, 0), oma = c(2.2, 1.5, 2.4, 0.5))

yrs_all <- sort(unique(out$year))
for (st in site_order) {
  d <- out[out$Site == st, ]
  lo <- d[d$bracket == "low",  ]; lo <- lo[order(lo$year), ]
  hi <- d[d$bracket == "high", ]; hi <- hi[order(hi$year), ]
  
  plot(NA, xlim = range(yrs_all), ylim = c(0, 1),
       xlab = "Year", ylab = "P(exceed)", main = st)
  
  # shade years where a majority of members are temp-extrapolated
  extrap_yrs <- d$year[!is.na(d$frac_temp_extrap) & d$frac_temp_extrap > 0.5]
  if (length(extrap_yrs)) {
    rect(min(extrap_yrs), -0.05, max(yrs_all) + 1, 1.05,
         col = adjustcolor("grey50", alpha.f = 0.12), border = NA)
  }
  
  if (nrow(lo)) {
    lines(lo$year, lo$P_exceed_naive,     col = col_low, lty = 2, lwd = 1)
    lines(lo$year, lo$P_exceed_corrected, col = col_low, lty = 1, lwd = 2)
  }
  if (nrow(hi)) {
    lines(hi$year, hi$P_exceed_naive,     col = col_high, lty = 2, lwd = 1)
    lines(hi$year, hi$P_exceed_corrected, col = col_high, lty = 1, lwd = 2)
  }
  box()
}

plot.new()
leg  <- c("Low, corrected", "Low, naive", "High, corrected", "High, naive",
          "temp-extrapolated (majority)")
lcol <- c(col_low, col_low, col_high, col_high, "grey50")
llty <- c(1, 2, 1, 2, NA)
llwd <- c(2, 1, 2, 1, NA)
lpch <- c(NA, NA, NA, NA, 15)
legend("center", legend = leg, col = lcol, lty = llty, lwd = llwd, pch = lpch,
       pt.cex = 2, bty = "n", cex = 0.9)

mtext("Exceedance probability check -- ensemble-only (naive) vs +residual bootstrap (corrected)",
      outer = TRUE, cex = 1.0, font = 2, line = 0.6)
mtext(paste0("threshold = ", THRESHOLD_CHLA_MG_M2, " mg Chl a/m^2 (Suplee et al. 2009, JAWRA)"),
      outer = TRUE, side = 1, cex = 0.75, line = 0.6)

dev.off()


# ============================================================================
# SCORECARD
# ============================================================================
cat("\n")
cat("============================================================\n")
cat("QAQC -- EXCEEDANCE PROBABILITY CHECK\n")
cat("Not manuscript-facing.\n")
cat("============================================================\n\n")

cat("Threshold:", THRESHOLD_CHLA_MG_M2, "mg Chl a/m2  (log10 =",
    round(THRESHOLD_LOG10, 4), ")\n")
cat("Residual pool: ", ifelse(POOL_RESIDUALS_ACROSS_SITES, "pooled across sites",
                              "per-site"),
    ", n =", length(resid_all),
    ", mean =", round(mean(resid_all), 4),
    ", SD =", round(sd(resid_all), 4), "\n")
cat("Draws per group: n_members x", B_DRAWS_PER_MEMBER, "\n\n")

cat("--- Headline: P(exceed) at year", bench_year, "---\n")
print(bench_ord[, c("Site", "bracket", "P_exceed_naive", "P_exceed_corrected",
                    "frac_temp_extrap")], row.names = FALSE)
cat("\n")

# low-vs-high separation at benchmark year: does it survive correction?
sep_tab <- do.call(rbind, lapply(site_order, function(st) {
  rl <- bench[bench$Site == st & bench$bracket == "low",  ]
  rh <- bench[bench$Site == st & bench$bracket == "high", ]
  if (!nrow(rl) || !nrow(rh)) return(NULL)
  data.frame(
    Site = st,
    sep_naive     = round(rh$P_exceed_naive     - rl$P_exceed_naive,     4),
    sep_corrected = round(rh$P_exceed_corrected - rl$P_exceed_corrected, 4)
  )
}))
cat("--- High-minus-low separation at year", bench_year, "---\n")
cat("(does the bracket contrast survive carrying real uncertainty?)\n")
print(sep_tab, row.names = FALSE)
cat("\n")

sc_files <- data.frame(
  File = c(OUT_CSV, OUT_PDF),
  Contents = c(
    "P(exceed) naive and corrected, by Site x bracket x year, plus temp-extrapolation fraction.",
    "Page 1: headline naive-vs-corrected comparison. Page 2: residual pool. Pages 3+: per-site trajectories."
  ),
  stringsAsFactors = FALSE
)
cat("--- Output files ---\n")
print(sc_files, row.names = FALSE)
cat("\n")
cat("============================================================\n")
cat("Done.\n")
cat("============================================================\n")