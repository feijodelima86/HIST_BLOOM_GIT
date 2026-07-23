## ==========================================================================
## utils_gam.R
## --------------------------------------------------------------------------
## Modular GAM utilities for UCFR Cladophora bloom prediction pipeline.
## Source this before running driver scripts (07g onward).
##
## Design notes:
##   - build_gam_formula() accepts per-predictor k via named list
##   - run_loso() returns .row_id for merging predictions into cascade
##   - All fitting uses REML + select = TRUE (shrinkage selection)
##   - Base R only
## ==========================================================================

library(mgcv)

# --------------------------------------------------------------------------
# Colour palette — consistent across all scripts
# --------------------------------------------------------------------------

site_colors <- function() {
  setNames(
    c("firebrick", "darkorange", "forestgreen", "steelblue",
      "gold3", "mediumpurple", "grey40"),
    c("DL", "GR", "BN", "MS", "BM", "HU", "FH"))
}

# --------------------------------------------------------------------------
# Formula builder
# --------------------------------------------------------------------------
# specs: named list, e.g.
#   list(anomaly = list(k = 10), Days_Since_Freshet = list(k = 5))
# Elements with no $k default to k = 5; $bs defaults to "ts".

build_gam_formula <- function(response, specs, bs_default = "ts") {
  terms <- vapply(names(specs), function(p) {
    s  <- specs[[p]]
    k  <- if (!is.null(s$k))  s$k  else 5
    bs <- if (!is.null(s$bs)) s$bs else bs_default
    sprintf('s(%s, bs = "%s", k = %d)', p, bs, k)
  }, character(1))
  as.formula(paste(response, "~", paste(terms, collapse = " + ")))
}

# --------------------------------------------------------------------------
# Fit GAM with project conventions
# --------------------------------------------------------------------------

fit_gam <- function(formula, data) {
  gam(formula, data = data, method = "REML", select = TRUE)
}

# --------------------------------------------------------------------------
# Prepare modelling data: complete cases for response + predictors
# --------------------------------------------------------------------------

prep_model_data <- function(dat, response, specs) {
  vars   <- c(response, names(specs))
  subset <- dat[complete.cases(dat[, vars]), ]
  cat(sprintf("  %s model data: n = %d  |  Sites: %s  |  Years: %s\n",
              response, nrow(subset),
              paste(levels(droplevels(subset$Site)), collapse = ", "),
              paste(range(subset$Year), collapse = "–")))
  subset
}

# --------------------------------------------------------------------------
# LOSO cross-validation
# --------------------------------------------------------------------------
# Returns data.frame with .row_id, Site, Year, Month, obs, pred
# .row_id is the row index within the supplied data (for merging).

run_loso <- function(formula, data, response, site_var = "Site") {
  data$.row_id <- seq_len(nrow(data))
  sites  <- levels(droplevels(data[[site_var]]))
  pieces <- vector("list", length(sites))
  
  for (i in seq_along(sites)) {
    s     <- sites[i]
    train <- data[data[[site_var]] != s, ]
    test  <- data[data[[site_var]] == s, ]
    m_cv  <- gam(formula, data = train, method = "REML", select = TRUE)
    
    pieces[[i]] <- data.frame(
      .row_id = test$.row_id,
      Site    = test[[site_var]],
      Year    = test$Year,
      Month   = test$Month,
      obs     = test[[response]],
      pred    = as.numeric(predict(m_cv, newdata = test, type = "response")),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, pieces)
  out[order(out$.row_id), ]
}

# --------------------------------------------------------------------------
# Temporal jackknife
# --------------------------------------------------------------------------

run_temporal_jk <- function(formula, data, response) {
  years  <- sort(unique(data$Year))
  pieces <- vector("list", length(years))
  
  for (i in seq_along(years)) {
    y     <- years[i]
    train <- data[data$Year != y, ]
    test  <- data[data$Year == y, ]
    m_ty  <- gam(formula, data = train, method = "REML", select = TRUE)
    
    pieces[[i]] <- data.frame(
      .row_id = seq_len(nrow(test)),
      Site    = test$Site,
      Year    = test$Year,
      Month   = test$Month,
      obs     = test[[response]],
      pred    = as.numeric(predict(m_ty, newdata = test, type = "response")),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, pieces)
}

# --------------------------------------------------------------------------
# Validation statistics
# --------------------------------------------------------------------------

validation_stats <- function(cv_preds) {
  ss_res <- sum((cv_preds$obs - cv_preds$pred)^2)
  ss_tot <- sum((cv_preds$obs - mean(cv_preds$obs))^2)
  r2     <- 1 - ss_res / ss_tot
  rmse   <- sqrt(mean((cv_preds$obs - cv_preds$pred)^2))
  c(R2 = r2, RMSE = rmse)
}

persite_stats <- function(cv_preds) {
  sites <- levels(droplevels(factor(cv_preds$Site,
                                    levels = c("DL","GR","BN","MS","BM","HU","FH"))))
  rows <- list()
  for (s in sites) {
    sub <- cv_preds[cv_preds$Site == s, ]
    if (nrow(sub) < 3) next
    ss_res <- sum((sub$obs - sub$pred)^2)
    ss_tot <- sum((sub$obs - mean(sub$obs))^2)
    rows[[length(rows) + 1]] <- data.frame(
      Site = s, n = nrow(sub),
      R2   = 1 - ss_res / ss_tot,
      RMSE = sqrt(mean((sub$obs - sub$pred)^2)),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# --------------------------------------------------------------------------
# Print validation summary
# --------------------------------------------------------------------------

print_validation <- function(loso_preds, temp_preds, label = "") {
  loso <- validation_stats(loso_preds)
  temp <- validation_stats(temp_preds)
  
  cat(sprintf("\n=== Validation: %s ===\n", label))
  cat(sprintf("  LOSO      R² = %.3f   RMSE = %.3f\n", loso["R2"], loso["RMSE"]))
  cat(sprintf("  Temporal  R² = %.3f   RMSE = %.3f\n", temp["R2"], temp["RMSE"]))
  
  cat("\nPer-site LOSO:\n")
  ps <- persite_stats(loso_preds)
  for (i in seq_len(nrow(ps))) {
    cat(sprintf("  %s : n = %3d,  R² = %6.3f,  RMSE = %.3f\n",
                ps$Site[i], ps$n[i], ps$R2[i], ps$RMSE[i]))
  }
}

# --------------------------------------------------------------------------
# Diagnostics: summary, concurvity, gam.check, smooth plots
# --------------------------------------------------------------------------

gam_diagnostics <- function(model, label = "") {
  cat("\n", strrep("=", 60), "\n")
  cat("  Diagnostics:", label, "\n")
  cat(strrep("=", 60), "\n\n")
  
  cat("=== Summary ===\n")
  print(summary(model))
  
  cat("\n=== Concurvity (worst) ===\n")
  print(round(concurvity(model, full = TRUE), 3))
  
  cat("\n=== Concurvity (pairwise, estimate) ===\n")
  cc <- concurvity(model, full = FALSE)
  print(round(cc$estimate, 3))
  
  cat("\n=== gam.check ===\n")
  par(mfrow = c(2, 2))
  gam.check(model)
  
  # Smooth partial-effect plots
  n_terms <- length(model$smooth)
  nr <- ceiling(n_terms / 3)
  par(mfrow = c(nr, min(n_terms, 3)), mar = c(4, 4, 2, 1))
  plot(model, shade = TRUE, shade.col = "lightblue",
       residuals = TRUE, pch = 16, cex = 0.4, col = "grey40",
       pages = 0)
  if (nzchar(label)) {
    mtext(label, side = 3, line = -1.5, outer = TRUE, cex = 0.9)
  }
}

# --------------------------------------------------------------------------
# Single-term drop-AIC table
# --------------------------------------------------------------------------

drop_aic_table <- function(model, response, specs, data) {
  predictors <- names(specs)
  full_aic   <- AIC(model)
  
  cat(sprintf("\n=== Drop-AIC (full = %.1f) ===\n", full_aic))
  
  for (d in predictors) {
    reduced_specs <- specs[setdiff(predictors, d)]
    f_drop <- build_gam_formula(response, reduced_specs)
    m_drop <- fit_gam(f_drop, data)
    delta  <- AIC(m_drop) - full_aic
    cat(sprintf("  Drop %-22s  AIC = %7.1f  (delta = %+.1f)\n",
                d, AIC(m_drop), delta))
  }
}

# --------------------------------------------------------------------------
# Plots
# --------------------------------------------------------------------------

plot_loso <- function(cv_preds, label = "") {
  cols  <- site_colors()
  stats <- validation_stats(cv_preds)
  
  par(mfrow = c(1, 1), mar = c(4.5, 4.5, 2, 1))
  plot(cv_preds$obs, cv_preds$pred,
       col  = cols[as.character(cv_preds$Site)],
       pch  = 16, cex = 0.9,
       xlab = paste("Observed", label),
       ylab = paste("LOSO Predicted", label),
       main = sprintf("%s LOSO  (R² = %.3f)", label, stats["R2"]))
  abline(0, 1, lty = 2)
  legend("topleft", legend = names(cols), col = cols,
         pch = 16, cex = 0.7, bty = "n", ncol = 2)
}

plot_loso_residuals <- function(cv_preds, label = "") {
  cols <- site_colors()
  cv_preds$resid <- cv_preds$obs - cv_preds$pred
  
  site_order <- intersect(c("DL","GR","BN","MS","BM","HU","FH"),
                          unique(as.character(cv_preds$Site)))
  cv_preds$Site <- factor(cv_preds$Site, levels = site_order)
  
  par(mfrow = c(1, 1), mar = c(4.5, 4.5, 2, 1))
  boxplot(resid ~ Site, data = cv_preds,
          col  = cols[site_order],
          xlab = "Site", ylab = "LOSO Residual",
          main = paste(label, ": LOSO residuals by site"))
  abline(h = 0, lty = 2)
}