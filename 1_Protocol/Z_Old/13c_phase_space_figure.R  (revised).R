# =============================================================================
# 13c_phase_space_figure.R  (revised)
# UCFR Cladophora Bloom Prediction Pipeline
#
# PURPOSE:
#   Figure B — Projected climate delta phase space.
#   X-axis: change in DOY of peak discharge relative to baseline (days;
#           negative = earlier freshet)
#   Y-axis: change in mean water temperature relative to baseline (degrees C)
#   Origin (0, 0) = baseline conditions.
#
#   Quadrant shading communicates ecological concern:
#     Upper left  (earlier + warmer) = salmon  -- cyanobacterial template
#     Upper right (later + warmer)   = light orange
#     Lower quadrants                = pale grey
#
#   Points: section mean per scenario x horizon x date point
#   Three panels: upper / mid / lower river
#
# INPUTS:
#   2_incremental/climate_deltas.csv
#   2_incremental/time_slice_summaries.csv
#
# OUTPUTS:
#   4_products/fig3_phase_space.pdf
#   4_products/fig3_phase_space.png
#
# RIVER SECTIONS:
#   Upper river: DL, GR
#   Mid-river:   BN, MS, BM
#   Lower river: HU, FH
#
# AUTHOR: [Rafa]
# =============================================================================


# -----------------------------------------------------------------------------
# 0. USER PARAMETERS
# -----------------------------------------------------------------------------

sections <- list(
  "Upper river\n(DL, GR)"    = c("DL", "GR"),
  "Mid-river\n(BN, MS, BM)"  = c("BN", "MS", "BM"),
  "Lower river\n(HU, FH)"    = c("HU", "FH")
)

scenarios  <- c("ssp245", "ssp585")
horizons   <- c("2050", "2080")

date_points       <- c("early_july", "august", "mid_september")
date_point_labels <- c("Early July", "August", "Mid-September")

col_45 <- "#2166AC"
col_85 <- "#B2182B"

dp_pch <- c("early_july" = 21, "august" = 22, "mid_september" = 23)

cex_2050 <- 1.6
cex_2080 <- 2.1

col_ul <- "#FDDBC7"   # upper left  -- earlier + warmer (salmon)
col_ur <- "#FEE8C8"   # upper right -- later + warmer   (light orange)
col_ll <- "#F7F7F7"   # lower left  -- earlier + cooler
col_lr <- "#F7F7F7"   # lower right -- later + cooler

lwd_ref <- 0.9
col_ref <- "grey40"

fig_width  <- 12
fig_height <- 5

out_pdf <- "4_products/fig3_phase_space.pdf"
out_png <- "4_products/fig3_phase_space.png"


# -----------------------------------------------------------------------------
# 1. LOAD INPUTS
# -----------------------------------------------------------------------------

cat("Loading inputs...\n")

deltas <- read.csv("2_incremental/climate_deltas.csv",       stringsAsFactors = FALSE)
slices <- read.csv("2_incremental/time_slice_summaries.csv", stringsAsFactors = FALSE)

cat("  Delta rows:", nrow(deltas), "\n")
cat("  Slice rows:", nrow(slices), "\n")


# -----------------------------------------------------------------------------
# 2. BUILD DELTA POINTS
# -----------------------------------------------------------------------------

delta_points_list <- list()

for (scen in scenarios) {
  for (hz in horizons) {
    for (dp in date_points) {
      
      dp_month <- switch(dp,
                         "early_july"    = 7,
                         "august"        = 8,
                         "mid_september" = 9
      )
      
      d_rows <- deltas[
        deltas$scenario == scen &
          deltas$horizon  == hz   &
          deltas$month    == dp_month,
      ]
      
      if (nrow(d_rows) == 0) next
      
      for (sec_name in names(sections)) {
        sec_sites <- sections[[sec_name]]
        sec_rows  <- d_rows[d_rows$site %in% sec_sites, ]
        if (nrow(sec_rows) == 0) next
        
        delta_points_list[[length(delta_points_list) + 1]] <- data.frame(
          section    = sec_name,
          scenario   = scen,
          horizon    = hz,
          date_point = dp,
          delta_x    = mean(sec_rows$delta_DOY_peak_med, na.rm = TRUE),
          delta_y    = mean(sec_rows$delta_Temp_C,       na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

delta_points <- do.call(rbind, delta_points_list)
rownames(delta_points) <- NULL

cat("  Delta point rows:", nrow(delta_points), "\n")
cat("  Expected:",
    length(scenarios) * length(horizons) * length(date_points) * length(sections), "\n")


# -----------------------------------------------------------------------------
# 3. GLOBAL AXIS LIMITS
# -----------------------------------------------------------------------------

x_range <- range(c(delta_points$delta_x, 0), na.rm = TRUE)
y_range <- range(c(delta_points$delta_y, 0), na.rm = TRUE)

x_pad <- max(diff(x_range) * 0.25, 2)
y_pad <- max(diff(y_range) * 0.25, 0.5)

x_lim <- c(x_range[1] - x_pad, x_range[2] + x_pad)
y_lim <- c(min(y_range[1] - y_pad, -0.5), y_range[2] + y_pad)

x_at <- pretty(x_lim, n = 6)
y_at <- pretty(y_lim, n = 5)

cat("  X lim:", round(x_lim, 1), "\n")
cat("  Y lim:", round(y_lim, 1), "\n")


# -----------------------------------------------------------------------------
# 4. HELPER: TRANSPARENT COLOR
# -----------------------------------------------------------------------------

make_transp <- function(hex_col, alpha_0_1) {
  rv <- col2rgb(hex_col) / 255
  rgb(rv[1], rv[2], rv[3], alpha = alpha_0_1)
}


# -----------------------------------------------------------------------------
# 5. SINGLE PANEL PLOT FUNCTION
# -----------------------------------------------------------------------------

plot_panel <- function(sec_name, show_yaxis, show_legend) {
  
  sec_df <- delta_points[delta_points$section == sec_name, ]
  
  plot(NA,
       xlim = x_lim,
       ylim = y_lim,
       xaxt = "n",
       yaxt = "n",
       xlab = "",
       ylab = "",
       bty  = "o")
  
  # Quadrant shading
  rect(x_lim[1], 0,        0,        y_lim[2], col = col_ul, border = NA)
  rect(0,        0,        x_lim[2], y_lim[2], col = col_ur, border = NA)
  rect(x_lim[1], y_lim[1], 0,        0,        col = col_ll, border = NA)
  rect(0,        y_lim[1], x_lim[2], 0,        col = col_lr, border = NA)
  
  box(col = "grey30")
  
  # Reference lines at origin
  abline(h = 0, lwd = lwd_ref, col = col_ref)
  abline(v = 0, lwd = lwd_ref, col = col_ref)
  
  # Grid
  abline(h = y_at[y_at != 0], col = "grey88", lwd = 0.4)
  abline(v = x_at[x_at != 0], col = "grey88", lwd = 0.4)
  
  # Points
  for (scen in scenarios) {
    scen_col <- if (scen == "ssp245") col_45 else col_85
    
    for (dp in date_points) {
      pch_i <- dp_pch[dp]
      
      for (hz in horizons) {
        cex_i <- if (hz == "2050") cex_2050 else cex_2080
        bg_i  <- if (hz == "2050") make_transp(scen_col, 0.5) else scen_col
        
        row_i <- sec_df[
          sec_df$scenario   == scen &
            sec_df$horizon    == hz   &
            sec_df$date_point == dp,
        ]
        if (nrow(row_i) == 0) next
        
        points(row_i$delta_x, row_i$delta_y,
               pch = pch_i,
               cex = cex_i,
               col = scen_col,
               bg  = bg_i,
               lwd = 1.2)
      }
    }
  }
  
  # Axes
  axis(1,
       at = x_at, labels = x_at,
       cex.axis = 0.75, col.axis = "grey20", tck = -0.025)
  
  if (show_yaxis) {
    axis(2,
         at = y_at, labels = y_at,
         las = 2, cex.axis = 0.75, col.axis = "grey20", tck = -0.025)
  } else {
    axis(2, at = y_at, labels = FALSE, tck = -0.02)
  }
  
  # Panel label
  mtext(sec_name, side = 3, line = 0.4, cex = 0.78, font = 2, col = "grey20")
  
  # Legend (first panel only)
  if (show_legend) {
    
    legend("topleft",
           legend  = c("RCP4.5 2050", "RCP4.5 2080",
                       "RCP8.5 2050", "RCP8.5 2080"),
           pch     = 21,
           pt.bg   = c(make_transp(col_45, 0.5), col_45,
                       make_transp(col_85, 0.5), col_85),
           col     = c(col_45, col_45, col_85, col_85),
           pt.cex  = c(cex_2050, cex_2080, cex_2050, cex_2080),
           bty     = "n",
           cex     = 0.68,
           title   = "Scenario x horizon",
           title.adj = 0)
    
    legend("bottomleft",
           legend    = date_point_labels,
           pch       = dp_pch,
           col       = "grey30",
           pt.bg     = "grey70",
           pt.cex    = 1.3,
           bty       = "n",
           cex       = 0.68,
           title     = "Date point",
           title.adj = 0)
  }
}


# -----------------------------------------------------------------------------
# 6. RENDER FIGURE
# -----------------------------------------------------------------------------

render_figure <- function(file_path, device_fn, ...) {
  
  device_fn(file_path, ...)
  
  n_sec <- length(sections)
  
  par(mfrow = c(1, n_sec),
      mar   = c(3.5, 3.8, 2.2, 0.8),
      oma   = c(2.0, 2.0, 2.5, 1.0))
  
  for (pi in seq_len(n_sec)) {
    plot_panel(
      sec_name    = names(sections)[pi],
      show_yaxis  = pi == 1,
      show_legend = pi == 1
    )
  }
  
  mtext("Change in peak discharge DOY (days; negative = earlier freshet)",
        side = 1, outer = TRUE, line = 0.5, cex = 0.82, col = "grey20")
  
  mtext("Change in water temperature (degrees C)",
        side = 2, outer = TRUE, line = 0.5, cex = 0.82, col = "grey20")
  
  mtext("Projected shifts in growing season environmental template by river section",
        side = 3, outer = TRUE, line = 1.2, cex = 0.92, font = 2, col = "grey15")
  
  dev.off()
  cat("  Written:", file_path, "\n")
}

cat("Rendering figures...\n")

render_figure(out_pdf, pdf,
              width = fig_width, height = fig_height, useDingbats = FALSE)

render_figure(out_png, png,
              width = fig_width, height = fig_height, units = "in", res = 300)

cat("Done.\n")