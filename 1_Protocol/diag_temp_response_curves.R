# ============================================================================
# diag_temp_response_curves.R
# Actual M1 response to Temp_oC, per site, under each temp-sensitivity
# variant -- real predict() calls on your fitted model, not an illustration.
#
# Holds lag_y, anomaly, logQ_obs_cfs, Days_Since_Freshet, logTP_mg_L at each
# site's median value (computed from the real projection grid), sweeps
# Temp_oC across the actual observed+projected range, and applies the SAME
# apply_temp_treatment() used in diag_temp_sensitivity_projection.R -- so
# this is a direct readout of what that script is actually doing to the
# model, not a schematic.
#
# Inputs:
#   3_models/bloom_model_M1.rds
#   4_products/diagnostics/temp_sensitivity/bloom_projections_members_V0.csv
#     (used only to get realistic per-site median predictor values to hold
#     the other 5 drivers at -- V0 file chosen arbitrarily since medians of
#     anomaly/logQ/DSF/logTP/lag_y don't depend on the temp treatment)
#
# Output:
#   4_products/diagnostics/temp_sensitivity/temp_response_curves.pdf
#   one panel per site, 7 panels + legend
# ============================================================================

library(mgcv)

PATH_M1      <- "3_models/bloom_model_M1.rds"
PATH_MEMBERS <- "4_products/diagnostics/temp_sensitivity/bloom_projections_members_V0.csv"
OUT_PDF      <- "4_products/diagnostics/temp_sensitivity/temp_response_curves.pdf"

TREATMENT <- list(
  V0       = list(type = "none"),
  V1       = list(type = "cap",     cap = 22.2),
  V2       = list(type = "cap",     cap = 28),
  V3mild   = list(type = "decline", threshold = 23, slope = 0.05),
  V3strong = list(type = "decline", threshold = 23, slope = 0.15)
)

# --- same apply_temp_treatment() as diag_temp_sensitivity_projection.R -----
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
    temp_term <- term_pred[, temp_col]
    temp_val  <- newdata$Temp_oC
    nd_frozen <- newdata
    nd_frozen$Temp_oC <- treatment$threshold
    frozen_val <- predict(model, newdata = nd_frozen, type = "terms")[, temp_col]
    new_temp_term <- ifelse(temp_val > treatment$threshold,
                            frozen_val - treatment$slope * (temp_val - treatment$threshold),
                            temp_term)
    term_pred[, temp_col] <- new_temp_term
    linear_pred <- rowSums(term_pred) + intercept
    return(as.numeric(model$family$linkinv(linear_pred)))
  }
  stop("unknown treatment$type: ", treatment$type)
}

m1  <- readRDS(PATH_M1)
mem <- read.csv(PATH_MEMBERS, stringsAsFactors = FALSE)

sites <- levels(m1$model$Site)
temp_seq <- seq(min(mem$Temp_oC), max(mem$Temp_oC), length.out = 80)

variant_col <- c(V0 = "#2a78d6", V1 = "#1baf7a", V2 = "#eda100",
                 V3mild = "#e34948", V3strong = "#4a3aa7")
variant_lty <- c(V0 = 1, V1 = 2, V2 = 1, V3mild = 3, V3strong = 3)

pdf(OUT_PDF, width = 9, height = 11)
par(mfrow = c(4, 2), mar = c(4, 4, 2.2, 1))

for (st in sites) {
  d <- mem[mem$Site == st, ]
  if (!nrow(d)) { plot.new(); next }
  
  base_row <- data.frame(
    Site               = factor(st, levels = sites),
    lag_y              = median(d$lag_y, na.rm = TRUE),
    anomaly            = median(d$anomaly, na.rm = TRUE),
    logQ_obs_cfs       = median(d$logQ_obs_cfs, na.rm = TRUE),
    Days_Since_Freshet = median(d$Days_Since_Freshet, na.rm = TRUE),
    logTP_mg_L         = median(d$logTP_mg_L, na.rm = TRUE)
  )
  
  nd <- base_row[rep(1, length(temp_seq)), ]
  nd$Temp_oC <- temp_seq
  
  curves <- sapply(names(TREATMENT), function(v) apply_temp_treatment(m1, nd, TREATMENT[[v]]))
  
  plot(NA, xlim = range(temp_seq), ylim = range(curves),
       xlab = "Temp_oC", ylab = "predicted logCHLa")
  for (v in names(TREATMENT)) {
    lines(temp_seq, curves[, v], col = variant_col[v], lty = variant_lty[v], lwd = 2)
  }
  rug(d$Temp_oC, col = adjustcolor("grey40", 0.5))
  title(st, line = 0.5, cex.main = 1)
}

plot.new()
legend("center", legend = names(TREATMENT), col = variant_col, lty = variant_lty,
       lwd = 2, bty = "n", title = "variant")

dev.off()
cat("Wrote", OUT_PDF, "\n")