# =============================================================================
# 13_paper_figures.R
# UCFR Cladophora Bloom Prediction Pipeline
#
# PURPOSE:
#   Publication-quality figures for climate projection results.
#   Figure 1: Time slice bloom predictions by site (3-panel, one per date point)
#             Sites on x-axis, 4 scenario×horizon colored points with error bars,
#             dotted baseline segment per site. Shared y-axis across panels.
#
# INPUTS:
#   2_incremental/time_slice_summaries.csv — from Script 12
#   2_incremental/baseline_observed.csv    — from Script 12
#
# OUTPUTS:
#   4_products/fig1_time_slice_by_site.pdf
#   4_products/fig1_time_slice_by_site.png
#
# DESIGN:
#   Colors: RCP4.5 2050 (light blue), RCP4.5 2080 (dark blue),
#           RCP8.5 2050 (light orange/red), RCP8.5 2080 (dark red)
#   Baseline: dotted horizontal segment per site (projected baseline mean)
#   Error bars: q25/q75 (not SE — these are ESM uncertainty bounds)
#   Y-axis: log10 CHLa (mg/m²), shared across panels
#   X-axis: sites upstream → downstream
#   BM flagged with asterisk — point-source caveat
#
# AUTHOR: [Rafa]
# =============================================================================


# -----------------------------------------------------------------------------
# 0. USER PARAMETERS
# -----------------------------------------------------------------------------

site_order      <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")
date_point_order <- c("early_july", "august", "mid_september")
date_point_labels <- c("Early July", "August", "Mid-September")

# Scenario × horizon combinations (order determines point x-offset)
combos <- data.frame(
  scenario   = c("ssp245", "ssp245", "ssp585", "ssp585"),
  horizon    = c("2050",   "2080",   "2050",   "2080"),
  label      = c("RCP4.5 2050", "RCP4.5 2080", "RCP8.5 2050", "RCP8.5 2080"),
  col        = c("#6BAED6",     "#08519C",     "#FD8D3C",     "#A50F15"),
  col_transp = c("#6BAED680",   "#08519C80",   "#FD8D3C80",   "#A50F1580"),
  pch        = c(21, 21, 21, 21),
  offset     = c(-0.3, -0.1, 0.1, 0.3),   # horizontal jitter within site
  stringsAsFactors = FALSE
)

# Point size and line width
cex_pt   <- 1.4
lwd_err  <- 1.8
lwd_base <- 1.2

# Y-axis limits (log10 scale) — shared across panels
y_lim <- c(1.0, 2.5)
y_at  <- seq(1.0, 2.5, by = 0.2)

# Figure dimensions
fig_width  <- 8    # inches
fig_height <- 10   # inches (3 stacked panels)

# Output paths
out_pdf <- "4_products/fig1_time_slice_by_site.pdf"
out_png <- "4_products/fig1_time_slice_by_site.png"


# -----------------------------------------------------------------------------
# 1. LOAD INPUTS
# -----------------------------------------------------------------------------

cat("Loading inputs...\n")

slices  <- read.csv("2_incremental/time_slice_summaries.csv", stringsAsFactors = FALSE)
obs_bl  <- read.csv("2_incremental/baseline_observed.csv",    stringsAsFactors = FALSE)

# Enforce site order as factor for x-axis positioning
slices$site  <- factor(slices$site,  levels = site_order)
obs_bl$site  <- factor(obs_bl$site,  levels = site_order)

n_sites <- length(site_order)
x_pos   <- seq_len(n_sites)   # integer x positions for sites


# -----------------------------------------------------------------------------
# 2. PLOT FUNCTION — SINGLE PANEL
# -----------------------------------------------------------------------------

plot_panel <- function(dp_label, dp_display, show_xaxis = FALSE,
                       show_legend = FALSE) {
  
  # Subset data for this date point
  sl_dp  <- slices[slices$date_point == dp_label, ]
  obs_dp <- obs_bl[obs_bl$date_point == dp_label, ]
  
  # Baseline slice (projected)
  base_dp <- sl_dp[sl_dp$horizon == "baseline", ]
  
  # Set up empty plot
  plot(NA,
       xlim = c(0.5, n_sites + 0.5),
       ylim = y_lim,
       xaxt = "n",
       yaxt = "n",
       xlab = "",
       ylab = "",
       bty  = "n")
  
  # Grid lines
  abline(h   = y_at,
         col = "grey90",
         lwd = 0.7)
  
  # Vertical site separators
  abline(v   = x_pos + 0.5,
         col = "grey85",
         lwd = 0.5,
         lty = 1)
  
  # Y-axis
  axis(2,
       at     = y_at,
       labels = formatC(y_at, digits = 1, format = "f"),
       las    = 2,
       cex.axis = 0.8,
       col.axis = "grey20",
       tck    = -0.02)
  
  # X-axis (only bottom panel)
  if (show_xaxis) {
    site_labels <- site_order
    site_labels[site_labels == "BM"] <- "BM*"   # flag point-source site
    axis(1,
         at       = x_pos,
         labels   = site_labels,
         cex.axis = 0.9,
         col.axis = "grey20",
         tck      = -0.02)
  }
  
  # Date point label (top-left of panel)
  mtext(dp_display,
        side = 3,
        adj  = 0.02,
        line = 0.3,
        cex  = 0.85,
        font = 2,
        col  = "grey25")
  
  # --- Baseline dotted segment per site ---
  for (i in seq_len(n_sites)) {
    site_i <- site_order[i]
    base_row <- base_dp[as.character(base_dp$site) == site_i, ]
    
    if (nrow(base_row) == 0) next
    
    # Use ssp245 baseline (identical for both scenarios by construction)
    bl_val <- base_row$pred_logCHLa_med[base_row$scenario == "ssp245"]
    if (length(bl_val) == 0) next
    
    segments(x0  = x_pos[i] - 0.45,
             x1  = x_pos[i] + 0.45,
             y0  = bl_val,
             y1  = bl_val,
             lty = 3,
             lwd = lwd_base,
             col = "grey40")
  }
  
  # --- Scenario × horizon points with error bars ---
  for (ci in seq_len(nrow(combos))) {
    scen  <- combos$scenario[ci]
    hz    <- combos$horizon[ci]
    col_i <- combos$col[ci]
    off_i <- combos$offset[ci]
    pch_i <- combos$pch[ci]
    
    sl_ci <- sl_dp[sl_dp$scenario == scen & sl_dp$horizon == hz, ]
    
    for (i in seq_len(n_sites)) {
      site_i <- site_order[i]
      row_i  <- sl_ci[as.character(sl_ci$site) == site_i, ]
      
      if (nrow(row_i) == 0) next
      
      xi   <- x_pos[i] + off_i
      ymed <- row_i$pred_logCHLa_med
      ylo  <- row_i$pred_logCHLa_lo
      yhi  <- row_i$pred_logCHLa_hi
      
      # Error bar
      segments(x0  = xi, x1  = xi,
               y0  = ylo, y1 = yhi,
               lwd = lwd_err,
               col = col_i)
      # Caps
      segments(x0  = xi - 0.04, x1 = xi + 0.04,
               y0  = ylo, y1 = ylo,
               lwd = lwd_err, col = col_i)
      segments(x0  = xi - 0.04, x1 = xi + 0.04,
               y0  = yhi, y1 = yhi,
               lwd = lwd_err, col = col_i)
      
      # Point
      points(xi, ymed,
             pch = pch_i,
             cex = cex_pt,
             col = col_i,
             bg  = col_i)
    }
  }
  
  # --- Legend (top panel only) ---
  if (show_legend) {
    legend("topright",
           legend = c(combos$label, "Projected baseline"),
           pch    = c(rep(21, nrow(combos)), NA),
           lty    = c(rep(NA, nrow(combos)), 3),
           lwd    = c(rep(NA, nrow(combos)), lwd_base),
           pt.bg  = c(combos$col, NA),
           col    = c(combos$col, "grey40"),
           pt.cex = cex_pt,
           bty    = "n",
           cex    = 0.78,
           x.intersp = 0.8)
  }
}


# -----------------------------------------------------------------------------
# 3. RENDER FIGURE — PDF AND PNG
# -----------------------------------------------------------------------------

render_figure <- function(file_path, device_fn, ...) {
  
  device_fn(file_path, ...)
  
  # Layout: 3 panels stacked, shared y-axis label
  par(mfrow  = c(3, 1),
      mar    = c(1.5, 4.5, 1.8, 1.2),   # bottom, left, top, right
      oma    = c(3.5, 0, 2.5, 0))        # outer margins for axis labels + title
  
  for (di in seq_along(date_point_order)) {
    dp     <- date_point_order[di]
    dp_lab <- date_point_labels[di]
    
    is_bottom <- di == length(date_point_order)
    is_top    <- di == 1
    
    plot_panel(dp_label    = dp,
               dp_display  = dp_lab,
               show_xaxis  = is_bottom,
               show_legend = is_top)
  }
  
  # Shared y-axis label
  mtext(expression(log[10]~CHLa~(mg~m^{-2})),
        side  = 2,
        outer = TRUE,
        line  = -1,
        cex   = 0.9,
        col   = "grey20")
  
  # Shared x-axis label (outer bottom)
  mtext("Site (upstream \u2192 downstream)",
        side  = 1,
        outer = TRUE,
        line  = 2,
        cex   = 0.9,
        col   = "grey20")
  
  # Main title
  mtext("Time slice bloom predictions by site",
        side  = 3,
        outer = TRUE,
        line  = 1,
        cex   = 1.05,
        font  = 2,
        col   = "grey15")
  
  # BM footnote
  mtext("* BM: below Missoula WWTP (point-source site; projections carry elevated uncertainty)",
        side  = 1,
        outer = TRUE,
        line  = 3.2,
        cex   = 0.65,
        adj   = 0,
        col   = "grey40")
  
  dev.off()
  cat("  Written:", file_path, "\n")
}

cat("Rendering figures...\n")

render_figure(out_pdf, pdf,
              width  = fig_width,
              height = fig_height,
              useDingbats = FALSE)

render_figure(out_png, png,
              width  = fig_width,
              height = fig_height,
              units  = "in",
              res    = 300)

cat("Done.\n")