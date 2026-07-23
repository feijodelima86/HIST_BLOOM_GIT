## ==========================================================================
## 07g_gam_submodels_kfix.R
## --------------------------------------------------------------------------
## Purpose : Fix k-index for anomaly and Q_peak_cfs (bumped to 10) in the
##           retained SPC and logTP submodels.  Run full diagnostics and
##           cross-validation, compare to 07f benchmarks (k = 5 throughout).
##
##           SPC retained:   Q_peak_cfs, Days_Since_Freshet, anomaly, Temp_oC
##           logTP retained: anomaly, Days_Since_Freshet, Q_peak_cfs
##
##           Q_obs_cfs was zeroed by shrinkage in both submodels (07f) and
##           is excluded here.
##
## Depends : utils_gam.R
## Inputs  : 2_incremental/ucfr_model_ready.csv
## Outputs : Console diagnostics, base-R plots
## ==========================================================================

source("1_protocol/utils_gam.R")   # loads mgcv

# --------------------------------------------------------------------------
# 1.  Load & prepare
# --------------------------------------------------------------------------

dat <- read.csv("2_incremental/ucfr_model_ready.csv", stringsAsFactors = FALSE)

dat$logTP <- log10(dat$TP_mg_L)
dat$Site  <- factor(dat$Site,
                    levels = c("DL", "GR", "BN", "MS", "BM", "HU", "FH"))

# ==========================================================================
# 2.  SPC submodel — k = 10 for anomaly, Q_peak_cfs
# ==========================================================================

spc_specs <- list(
  Q_peak_cfs         = list(k = 7),
  Days_Since_Freshet = list(k = 7),
  anomaly            = list(k = 7),
  Temp_oC            = list(k = 7)
)

dat_spc <- prep_model_data(dat, "SPC", spc_specs)
f_spc   <- build_gam_formula("SPC", spc_specs)
cat("  Formula:", deparse(f_spc), "\n\n")

m_spc <- fit_gam(f_spc, dat_spc)

# --- Diagnostics (check k-index in gam.check output) ---
gam_diagnostics(m_spc, "SPC  [k-fix: anomaly=10, Q_peak=10]")

# --- Cross-validation ---
loso_spc <- run_loso(f_spc, dat_spc, "SPC")
temp_spc <- run_temporal_jk(f_spc, dat_spc, "SPC")
print_validation(loso_spc, temp_spc, "SPC (k-fix)")

plot_loso(loso_spc, "SPC")
plot_loso_residuals(loso_spc, "SPC")

# --- Drop-AIC confirmation ---
drop_aic_table(m_spc, "SPC", spc_specs, dat_spc)

# ==========================================================================
# 3.  logTP submodel — k = 10 for anomaly, Q_peak_cfs
# ==========================================================================

tp_specs <- list(
  anomaly            = list(k = 7),
  Days_Since_Freshet = list(k = 7),
  Q_peak_cfs         = list(k = 7)
)

dat_tp <- prep_model_data(dat, "logTP", tp_specs)
f_tp   <- build_gam_formula("logTP", tp_specs)
cat("  Formula:", deparse(f_tp), "\n\n")

m_tp <- fit_gam(f_tp, dat_tp)

# --- Diagnostics ---
gam_diagnostics(m_tp, "logTP  [k-fix: anomaly=10, Q_peak=10]")

# --- Cross-validation ---
loso_tp <- run_loso(f_tp, dat_tp, "logTP")
temp_tp <- run_temporal_jk(f_tp, dat_tp, "logTP")
print_validation(loso_tp, temp_tp, "logTP (k-fix)")

plot_loso(loso_tp, "logTP")
plot_loso_residuals(loso_tp, "logTP")

# --- Drop-AIC confirmation ---
drop_aic_table(m_tp, "logTP", tp_specs, dat_tp)

# ==========================================================================
# 4.  Comparison to 07f benchmarks
# ==========================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("  k-FIX  vs  07f BENCHMARK  (k = 5 throughout)\n")
cat(strrep("=", 60), "\n\n")

spc_v  <- validation_stats(loso_spc)
spc_tv <- validation_stats(temp_spc)
tp_v   <- validation_stats(loso_tp)
tp_tv  <- validation_stats(temp_tp)

cat("                         07f (k=5)     k-fix\n")
cat(sprintf("  SPC   LOSO R²          0.716         %.3f\n",   spc_v["R2"]))
cat(sprintf("  SPC   Temporal R²      0.781         %.3f\n",   spc_tv["R2"]))
cat(sprintf("  logTP LOSO R²          0.195         %.3f\n",   tp_v["R2"]))
cat(sprintf("  logTP Temporal R²      0.431         %.3f\n\n", tp_tv["R2"]))

cat("Decision criteria:\n")
cat("  1. gam.check k-index for anomaly & Q_peak_cfs now adequate? (> 0.8)\n")
cat("  2. LOSO R² stable or improved vs 07f?\n")
cat("  3. edf values reasonable (not hitting k-1 ceiling)?\n")
cat("  If all yes → proceed to cascade assembly.\n")

