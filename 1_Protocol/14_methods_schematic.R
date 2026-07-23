# 14_methods_schematic.R
# Two-panel pipeline schematic for the manuscript Methods section.
# Panel A: observed data -> model fitting (M1 GAM, TP submodel) -> validation
# Panel B: NCAR climate envelopes -> delta-anchoring -> recursive projection -> output
# Base R graphics only. Purely illustrative -- no data dependencies.

OUT_PDF <- "4_products/methods_schematic.pdf"

# ---- palette (matched to the reviewed mockup: blue = observed, amber = ----
# ---- simulated/climate, teal = fitted model, gray = process/output)   ----
col_blue   <- list(fill = "#E6F1FB", border = "#185FA5", text = "#0C447C")
col_amber  <- list(fill = "#FAEEDA", border = "#854F0B", text = "#633806")
col_teal   <- list(fill = "#E1F5EE", border = "#0F6E56", text = "#085041")
col_gray   <- list(fill = "#F1EFE8", border = "#5F5E5A", text = "#2C2C2A")
col_line   <- "grey35"

# ---- helpers ----
# x, y = top-left corner (y increases DOWN the page, matching the design layout)
draw_box <- function(x, y, w, h, title, subtitle = NULL, pal,
                     cex.t = 0.80, cex.s = 0.66, lwd = 1) {
  rect(x, y, x + w, y + h, col = pal$fill, border = pal$border, lwd = lwd)
  if (is.null(subtitle)) {
    text(x + w / 2, y + h / 2, title, cex = cex.t, col = pal$text, font = 2)
  } else {
    text(x + w / 2, y + h * 0.36, title,    cex = cex.t, col = pal$text, font = 2)
    text(x + w / 2, y + h * 0.74, subtitle, cex = cex.s, col = pal$text, font = 1)
  }
}

# straight or L-shaped connector; xs/ys give the waypoints, arrowhead on the last segment
poly_arrow <- function(xs, ys, col = col_line, lwd = 1.1) {
  n <- length(xs)
  if (n > 2) {
    segments(xs[1:(n - 2)], ys[1:(n - 2)], xs[2:(n - 1)], ys[2:(n - 1)],
             col = col, lwd = lwd)
  }
  arrows(xs[n - 1], ys[n - 1], xs[n], ys[n],
         length = 0.07, angle = 25, col = col, lwd = lwd)
}

# swatch legend: items = list(list(pal=, label=)), laid out left to right at a given y
draw_legend <- function(y, items, x0 = 70, gap = 170, sw = 13) {
  for (i in seq_along(items)) {
    xi <- x0 + (i - 1) * gap
    rect(xi, y, xi + sw, y + sw, col = items[[i]]$pal$fill, border = items[[i]]$pal$border, lwd = 0.8)
    text(xi + sw + 6, y + sw / 2, items[[i]]$label, cex = 0.62, col = "grey25", adj = c(0, 0.5))
  }
}

new_panel <- function(height, label) {
  plot.new()
  plot.window(xlim = c(0, 680), ylim = c(height, 0), xaxs = "i", yaxs = "i")
  text(4, 14, label, cex = 1.05, font = 2, adj = c(0, 0.5))
}

# ================================================================
# Panel A -- data & model fitting
# ================================================================
draw_panel_a <- function() {
  new_panel(480, "(a)")
  
  draw_box(70, 40, 260, 56, "WQP + VNRP", "Chemistry + bio sampling", col_blue)
  draw_box(350, 40, 260, 56, "USGS daily Q", "Site-level discharge", col_blue)
  
  poly_arrow(c(200, 255), c(96, 140))
  poly_arrow(c(480, 425), c(96, 140))
  
  draw_box(180, 140, 320, 56, "Processing (03-06)", "Lag-match, hydrology join", col_gray)
  poly_arrow(c(340, 340), c(196, 220))
  
  draw_box(190, 220, 300, 56, "ucfr_model_ready.csv", "n=284 raw, 221 clean", col_blue)
  
  poly_arrow(c(340, 340, 200, 200), c(276, 290, 290, 300))
  poly_arrow(c(340, 340, 480, 480), c(276, 290, 290, 300))
  
  draw_box(70, 300, 260, 56, "M1 GAM (10)", "6 predictors + site RE", col_teal)
  draw_box(350, 300, 260, 56, "TP submodel (09)", "logTP, 4 drivers + RE", col_teal)
  
  poly_arrow(c(200, 200), c(356, 380))
  draw_box(70, 380, 260, 56, "Validation", "LOYO R2 = 0.710", col_teal)
  
  draw_legend(455, list(
    list(pal = col_blue, label = "Observed data"),
    list(pal = col_gray, label = "Process step"),
    list(pal = col_teal, label = "Fitted model")
  ))
}

# ================================================================
# Panel B -- projection engine, 2026-2098
# ================================================================
draw_panel_b <- function() {
  new_panel(880, "(b)")
  
  draw_box(70, 40, 260, 56, "Discharge envelope", "07c, MA20 smoothing", col_amber)
  draw_box(350, 40, 260, 56, "Temperature envelope", "07d, low/high bracket", col_amber)
  
  poly_arrow(c(200, 200), c(96, 140))
  poly_arrow(c(480, 480), c(96, 240))
  
  draw_box(70, 140, 260, 56, "Baseline vs future", "Delta over both windows", col_amber)
  poly_arrow(c(200, 200), c(196, 240))
  
  draw_box(70, 240, 260, 56, "Observed anchor", "2008-2025 site means", col_amber)
  draw_box(350, 240, 260, 56, "Temp_oC (bracket)", "Direct feed, no delta", col_amber)
  
  poly_arrow(c(200, 200), c(296, 340))
  
  draw_box(70, 340, 540, 150, NULL, NULL, col_amber, lwd = 1.4)
  text(94, 364, "Delta-anchoring", cex = 0.85, font = 2, col = col_amber$text, adj = c(0, 0.5))
  text(94, 390, "logQ = obs_logQ + log10(future_meanQ / baseline_meanQ)",
       cex = 0.66, col = col_amber$text, adj = c(0, 0.5), family = "mono")
  text(94, 410, "anomaly = (future_anom / baseline_anom) x obs_anomaly",
       cex = 0.66, col = col_amber$text, adj = c(0, 0.5), family = "mono")
  text(94, 430, "DSF = obs_DSF - (future_days - baseline_days)",
       cex = 0.66, col = col_amber$text, adj = c(0, 0.5), family = "mono")
  text(94, 450, "logTP = TP-submodel(anomaly, logQ, DSF, Temp), no delta",
       cex = 0.66, col = col_amber$text, adj = c(0, 0.5), family = "mono")
  
  poly_arrow(c(480, 480, 625, 625, 490), c(296, 320, 320, 578, 578))
  poly_arrow(c(480, 480, 625, 625, 490), c(296, 320, 320, 694, 694))
  
  poly_arrow(c(340, 340), c(490, 550))
  draw_box(190, 550, 300, 56, "TP submodel predict", "Produces future logTP", col_teal)
  
  poly_arrow(c(120, 120, 55, 55, 190), c(490, 510, 510, 694, 694))
  
  poly_arrow(c(340, 340), c(606, 666))
  draw_box(190, 666, 300, 56, "M1 predict", "6 predictors + site RE", col_teal)
  
  poly_arrow(c(340, 340), c(722, 782))
  text(360, 750, "recursive: logCHLa(t) -> lag_y(t+1)", cex = 0.64, col = col_line, adj = c(0, 0.5))
  draw_box(190, 782, 300, 56, "Projection outputs", "21 members x 7 sites", col_gray)
  
  draw_legend(850, list(
    list(pal = col_amber, label = "Climate / projection input"),
    list(pal = col_teal, label = "Model application"),
    list(pal = col_gray, label = "Output")
  ))
}

# ================================================================
# render
# ================================================================
dir.create(dirname(OUT_PDF), showWarnings = FALSE, recursive = TRUE)

pdf(OUT_PDF, width = 7, height = 14.0, pointsize = 10)
layout(matrix(1:2, nrow = 2), heights = c(480, 880))
par(mar = c(0.5, 0.5, 0.5, 0.5))

draw_panel_a()
draw_panel_b()

dev.off()

cat("Wrote", OUT_PDF, "\n")