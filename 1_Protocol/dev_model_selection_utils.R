# ============================================================================
# dev_model_selection_utils.R
# UCFR Cladophora Bloom Prediction Pipeline -- model-selection sandbox
# ----------------------------------------------------------------------------
# NOT manuscript-facing. NOT the numbered pipeline. This is scratch tooling
# to freely try different predictor sets for EITHER the CHLa submodel (what
# 10_bloom_model_M1.R currently is) or the TP submodel (09_tp_submodel.R)
# before anything gets locked back into those numbered scripts.
#
# Source this file, then use dev_chla_candidate.R / dev_tp_candidate.R to
# actually try things -- those are the only files you should need to edit
# to add or remove a variable.
#
# Every trial prints, in order: sample flow, in-sample fit, per-term
# significance, concurvity (per-term AND full pairwise matrix -- this is
# deliberately always-on after today), drop-one-term AIC, LOSO, LOYO.
#
# Design notes / things worth knowing before you start swapping variables:
#
#   - Outlier removal (+/- OUTLIER_SD, OUTLIER_ROUNDS passes) is redone from
#     scratch for EVERY candidate set, matching how 09/10 currently work.
#     This means sample size and even which specific rows survive can shift
#     between trials -- that's real and appropriate (a better predictor set
#     may legitimately explain away what looked like an outlier under a
#     worse one), but it also means two trials' R2/RMSE aren't strictly
#     apples-to-apples unless you set outlier_rounds = 0 for a clean-data-
#     held-fixed comparison.
#
#   - select_penalty (off by default, matching current M1/TP practice) maps
#     to mgcv's own select=TRUE, which adds an extra shrinkage penalty that
#     can drive a weak term's edf toward ~0 on its own. There's an old,
#     orphaned utils file in this repo (0_utils_gam.R, used only by the
#     archived 07-series scripts in Z_Old) that had select=TRUE on by
#     default -- that's a real, different modeling choice, not just a
#     styling difference, so it isn't silently carried over here. Flip the
#     toggle if you want to see what automatic shrinkage does to lag_y/
#     logQ/logTP's edf.
#
#   - Drop-one-term AIC is refit under REML, same as everything else, to
#     stay consistent with project convention. Worth knowing: comparing AIC
#     across models fit by REML is most defensible when the smoothing
#     structure changes and the mean structure doesn't -- here we're
#     changing which predictors are in the model at all, which is the
#     harder case. Treat Delta_AIC as a second opinion alongside the
#     p-values and concurvity, not a standalone verdict.
# ============================================================================

library(mgcv)

# ----------------------------------------------------------------------------
# Build formula from a named k_spec vector (predictor = k) and fit.
# ----------------------------------------------------------------------------
fit_gam_from_kspec <- function(response, k_spec, data, site_re = TRUE,
                               select_penalty = FALSE) {
  predictors <- names(k_spec)
  terms <- vapply(predictors, function(p)
    paste0("s(", p, ", k=", k_spec[[p]], ")"), character(1))
  if (site_re) terms <- c(terms, 's(Site, bs="re")')
  f <- as.formula(paste(response, "~", paste(terms, collapse = " + ")))
  gam(f, data = data, method = "REML", select = select_penalty)
}

# ----------------------------------------------------------------------------
# Complete cases (dynamic on whatever's in k_spec right now) + outlier removal
# ----------------------------------------------------------------------------
prep_clean_data <- function(dat, response, k_spec, site_re = TRUE,
                            outlier_sd = 2.0, outlier_rounds = 2,
                            select_penalty = FALSE) {
  predictors <- names(k_spec)
  missing_needed <- setdiff(c(response, predictors), names(dat))
  if (length(missing_needed) > 0)
    stop("Not found in data: ", paste(missing_needed, collapse = ", "))
  
  dat$.obs_id <- seq_len(nrow(dat))
  keep_cols <- unique(intersect(
    c(".obs_id", "Site", "Year", "Month", response, predictors), names(dat)))
  
  cc   <- complete.cases(dat[, c(response, predictors)])
  mdat <- dat[cc, keep_cols]
  mdat$Site <- factor(mdat$Site)
  
  n_raw      <- nrow(dat)
  n_complete <- nrow(mdat)
  
  outlier_log <- list()
  if (outlier_rounds > 0) {
    for (r in seq_len(outlier_rounds)) {
      m_tmp   <- fit_gam_from_kspec(response, k_spec, mdat, site_re, select_penalty)
      resid_r <- residuals(m_tmp)
      sd_r    <- sd(resid_r)
      flag    <- which(abs(resid_r) > outlier_sd * sd_r)
      if (length(flag) == 0) break
      log_r <- mdat[flag, intersect(c(".obs_id", "Site", "Year", "Month"), names(mdat))]
      log_r$round <- r
      log_r$resid <- round(resid_r[flag], 4)
      outlier_log[[r]] <- log_r
      mdat <- mdat[-flag, ]
      mdat$Site <- factor(mdat$Site)
    }
  }
  
  list(data = mdat, n_raw = n_raw, n_complete = n_complete, n_clean = nrow(mdat),
       outliers = if (length(outlier_log) > 0) do.call(rbind, outlier_log) else data.frame())
}

# ----------------------------------------------------------------------------
# Concurvity: per-term worst/observed vs rest of model + full pairwise matrix
# with flagged (>0.8) pairs. This is the check that caught Site<->logQ today
# -- always run, not optional.
# ----------------------------------------------------------------------------
concurvity_report <- function(model) {
  cc_full   <- concurvity(model, full = TRUE)
  worst_mat <- concurvity(model, full = FALSE)$worst
  
  clean_names <- function(x) gsub('^s\\(|\\)$', "", x)
  colnames(cc_full)   <- clean_names(colnames(cc_full))
  colnames(worst_mat) <- clean_names(colnames(worst_mat))
  rownames(worst_mat) <- clean_names(rownames(worst_mat))
  
  summary_df <- data.frame(
    Term          = colnames(cc_full),
    Worst_vs_rest = round(cc_full["worst", ], 4),
    Obs_vs_rest   = round(cc_full["observed", ], 4),
    stringsAsFactors = FALSE
  )
  summary_df$Flag <- ifelse(summary_df$Worst_vs_rest > 0.9, "!!",
                            ifelse(summary_df$Worst_vs_rest > 0.8, "!", ""))
  row.names(summary_df) <- NULL
  
  worst_mat_r <- round(worst_mat, 4)
  high_idx <- which(worst_mat_r > 0.8, arr.ind = TRUE)
  if (nrow(high_idx) > 0) {
    pair_df <- data.frame(
      Term_A = rownames(worst_mat_r)[high_idx[, 1]],
      Term_B = colnames(worst_mat_r)[high_idx[, 2]],
      Worst  = worst_mat_r[high_idx],
      stringsAsFactors = FALSE
    )
    pair_df <- pair_df[pair_df$Term_A < pair_df$Term_B, ]
    row.names(pair_df) <- NULL
  } else {
    pair_df <- data.frame(Term_A = character(), Term_B = character(), Worst = numeric())
  }
  
  list(per_term = summary_df, pairwise = worst_mat_r, flagged_pairs = pair_df)
}

# ----------------------------------------------------------------------------
# Drop-one-term AIC: refit with each term removed, one at a time.
# Positive Delta_AIC = removing it makes the model worse (term earns its
# place). Negative = model is no worse (or better) without it.
# ----------------------------------------------------------------------------
drop_one_aic <- function(response, k_spec, data, site_re = TRUE, select_penalty = FALSE) {
  predictors <- names(k_spec)
  full_model <- fit_gam_from_kspec(response, k_spec, data, site_re, select_penalty)
  full_aic   <- AIC(full_model)
  
  rows <- lapply(predictors, function(p) {
    k_sub <- k_spec[setdiff(predictors, p)]
    if (length(k_sub) == 0) {
      f_sub <- if (site_re) as.formula(paste(response, '~ s(Site, bs="re")'))
      else as.formula(paste(response, "~ 1"))
      m_sub <- gam(f_sub, data = data, method = "REML", select = select_penalty)
    } else {
      m_sub <- fit_gam_from_kspec(response, k_sub, data, site_re, select_penalty)
    }
    data.frame(Term_dropped = p,
               AIC_without  = round(AIC(m_sub), 2),
               Delta_AIC    = round(AIC(m_sub) - full_aic, 2),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out <- out[order(-out$Delta_AIC), ]
  row.names(out) <- NULL
  out
}

# ----------------------------------------------------------------------------
# LOSO -- leave-one-site-out, RE excluded for the held-out site (population-
# level prediction, same exclude pattern as 09/10).
# ----------------------------------------------------------------------------
run_loso <- function(response, k_spec, data, site_re = TRUE, select_penalty = FALSE) {
  if (!"Site" %in% names(data)) stop("run_loso needs a Site column.")
  sites <- levels(factor(data$Site))
  rows  <- list()
  for (s in sites) {
    train <- data[data$Site != s, ]; train$Site <- factor(train$Site)
    test  <- data[data$Site == s, ]
    if (nrow(test) == 0) next
    m <- tryCatch(fit_gam_from_kspec(response, k_spec, train, site_re, select_penalty),
                  error = function(e) NULL)
    if (is.null(m)) next
    
    if (site_re) {
      test2 <- test
      test2$Site <- factor(levels(train$Site)[1], levels = levels(train$Site))
      pred <- as.numeric(predict(m, newdata = test2,
                                 exclude = 's(Site, bs="re")',
                                 newdata.guaranteed = TRUE))
    } else {
      pred <- as.numeric(predict(m, newdata = test))
    }
    
    rows[[s]] <- data.frame(Site = s, n = nrow(test),
                            Observed = test[[response]], Predicted = pred,
                            stringsAsFactors = FALSE)
  }
  detail <- do.call(rbind, rows); row.names(detail) <- NULL
  
  per_site <- do.call(rbind, lapply(split(detail, detail$Site), function(d) {
    ss_res <- sum((d$Observed - d$Predicted)^2)
    ss_tot <- sum((d$Observed - mean(d$Observed))^2)
    data.frame(Site = d$Site[1], n = nrow(d),
               R2   = ifelse(ss_tot > 0, round(1 - ss_res / ss_tot, 4), NA),
               RMSE = round(sqrt(mean((d$Observed - d$Predicted)^2)), 4))
  }))
  row.names(per_site) <- NULL
  
  pooled_ss_res <- sum((detail$Observed - detail$Predicted)^2)
  pooled_ss_tot <- sum((detail$Observed - mean(detail$Observed))^2)
  
  list(per_site    = per_site,
       pooled_R2   = round(1 - pooled_ss_res / pooled_ss_tot, 4),
       pooled_RMSE = round(sqrt(mean((detail$Observed - detail$Predicted)^2)), 4))
}

# ----------------------------------------------------------------------------
# LOYO -- leave-one-year-out, RE retained (site is always known at
# projection time; only the year is out of sample).
# ----------------------------------------------------------------------------
run_loyo <- function(response, k_spec, data, site_re = TRUE, select_penalty = FALSE) {
  if (!"Year" %in% names(data)) stop("run_loyo needs a Year column.")
  years <- sort(unique(data$Year))
  rows  <- list()
  for (y in years) {
    train <- data[data$Year != y, ]
    test  <- data[data$Year == y, ]
    if (nrow(test) == 0) next
    m <- tryCatch(fit_gam_from_kspec(response, k_spec, train, site_re, select_penalty),
                  error = function(e) NULL)
    if (is.null(m)) next
    pred <- as.numeric(predict(m, newdata = test))
    rows[[as.character(y)]] <- data.frame(Year = y, n = nrow(test),
                                          Observed = test[[response]], Predicted = pred,
                                          stringsAsFactors = FALSE)
  }
  detail <- do.call(rbind, rows); row.names(detail) <- NULL
  
  per_year <- do.call(rbind, lapply(split(detail, detail$Year), function(d) {
    ss_res <- sum((d$Observed - d$Predicted)^2)
    ss_tot <- sum((d$Observed - mean(d$Observed))^2)
    data.frame(Year = d$Year[1], n = nrow(d),
               R2   = ifelse(ss_tot > 0, round(1 - ss_res / ss_tot, 4), NA),
               RMSE = round(sqrt(mean((d$Observed - d$Predicted)^2)), 4))
  }))
  row.names(per_year) <- NULL
  per_year <- per_year[order(per_year$Year), ]
  
  pooled_ss_res <- sum((detail$Observed - detail$Predicted)^2)
  pooled_ss_tot <- sum((detail$Observed - mean(detail$Observed))^2)
  
  list(per_year    = per_year,
       pooled_R2   = round(1 - pooled_ss_res / pooled_ss_tot, 4),
       pooled_RMSE = round(sqrt(mean((detail$Observed - detail$Predicted)^2)), 4))
}

# ----------------------------------------------------------------------------
# Orchestrator: fit one candidate model and print the full scorecard.
# Returns everything invisibly so you can grab it (e.g. result$model) if you
# want to dig further without re-running.
# ----------------------------------------------------------------------------
run_trial <- function(label, response, k_spec, data_full, site_re = TRUE,
                      outlier_sd = 2.0, outlier_rounds = 2,
                      select_penalty = FALSE,
                      do_loso = TRUE, do_loyo = TRUE, do_drop_aic = TRUE) {
  
  predictors    <- names(k_spec)
  terms_display <- paste0("s(", predictors, ", k=", k_spec, ")")
  formula_str   <- paste(terms_display, collapse = " + ")
  if (site_re) formula_str <- paste(formula_str, '+ s(Site, bs="re")')
  
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("TRIAL:", label, "\n")
  cat("Formula:", response, "~", formula_str, "\n")
  cat("select penalty (extra shrinkage):", select_penalty, "\n")
  cat(strrep("=", 70), "\n\n")
  
  prep <- prep_clean_data(data_full, response, k_spec, site_re,
                          outlier_sd, outlier_rounds, select_penalty)
  mdat <- prep$data
  
  cat("--- Sample flow ---\n")
  print(data.frame(
    Stage = c("Raw rows supplied", "Complete cases",
              "Outliers removed", "Clean n for fitting"),
    n = c(prep$n_raw, prep$n_complete,
          prep$n_complete - prep$n_clean, prep$n_clean)
  ), row.names = FALSE)
  cat("\n")
  
  m <- fit_gam_from_kspec(response, k_spec, mdat, site_re, select_penalty)
  s <- summary(m)
  
  cat("--- In-sample fit (n =", nrow(mdat), ") ---\n")
  print(data.frame(
    Metric = c("R2(adj)", "Deviance explained (%)", "RMSE", "AIC"),
    Value  = c(round(s$r.sq, 4), round(s$dev.expl * 100, 2),
               round(sqrt(mean(residuals(m)^2)), 4), round(AIC(m), 2))
  ), row.names = FALSE)
  cat("\n")
  
  cat("--- Smooth term summary ---\n")
  print(data.frame(
    Term = rownames(s$s.table),
    edf  = round(s$s.table[, "edf"], 3),
    F    = round(s$s.table[, "F"], 3),
    p    = formatC(s$s.table[, "p-value"], format = "e", digits = 2),
    stringsAsFactors = FALSE
  ), row.names = FALSE)
  cat("\n")
  
  cc <- concurvity_report(m)
  cat("--- Concurvity: per-term worst/observed vs rest of model ---\n")
  cat("(> 0.8 flagged '!', > 0.9 flagged '!!')\n")
  print(cc$per_term, row.names = FALSE)
  cat("\n--- Concurvity: flagged pairwise (worst-case > 0.8) ---\n")
  if (nrow(cc$flagged_pairs) > 0) print(cc$flagged_pairs, row.names = FALSE) else cat("None.\n")
  cat("\n")
  
  if (do_drop_aic && length(predictors) > 1) {
    cat("--- Drop-one-term AIC ---\n")
    cat("(positive Delta_AIC = term earns its place; negative = model isn't\n")
    cat(" worse without it. Second opinion alongside p-values above, not a\n")
    cat(" standalone verdict -- see header note on REML AIC comparisons.)\n")
    print(drop_one_aic(response, k_spec, mdat, site_re, select_penalty), row.names = FALSE)
    cat("\n")
  }
  
  loso_res <- NULL
  if (do_loso) {
    loso_res <- run_loso(response, k_spec, mdat, site_re, select_penalty)
    cat("--- LOSO (leave-one-site-out, RE excluded) ---\n")
    cat("Pooled R2 =", loso_res$pooled_R2, " RMSE =", loso_res$pooled_RMSE, "\n")
    print(loso_res$per_site, row.names = FALSE)
    cat("\n")
  }
  
  loyo_res <- NULL
  if (do_loyo) {
    loyo_res <- run_loyo(response, k_spec, mdat, site_re, select_penalty)
    cat("--- LOYO (leave-one-year-out, RE retained) ---\n")
    cat("Pooled R2 =", loyo_res$pooled_R2, " RMSE =", loyo_res$pooled_RMSE, "\n")
    print(loyo_res$per_year, row.names = FALSE)
    cat("\n")
  }
  
  if (nrow(prep$outliers) > 0) {
    cat("--- Outliers removed ---\n")
    print(prep$outliers, row.names = FALSE)
    cat("\n")
  }
  
  cat(strrep("=", 70), "\n")
  cat("Done:", label, "\n")
  cat(strrep("=", 70), "\n")
  
  invisible(list(model = m, data = mdat, concurvity = cc, loso = loso_res, loyo = loyo_res))
}

