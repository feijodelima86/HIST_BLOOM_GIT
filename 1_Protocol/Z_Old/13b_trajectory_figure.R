# =============================================================================
# 13b_trajectory_figure.R
# UCFR Cladophora Bloom Prediction Pipeline
#
# PURPOSE:
#   Full trajectory figure: 7 rows (sites) × 3 columns (date points).
#   Each panel shows ssp245 and ssp585 scenario lines with shaded q25/q75
#   uncertainty ribbons, from 2000 to 2090.
#
# INPUTS:
#   2_incremental/projections_monthly.csv  — from Script 11
#
# OUTPUTS:
#   4_products/fig2_full_trajectory.pdf
#   4_products/fig2_full_trajectory.png
#
# DESIGN:
#   Lines:   ssp245 = blue (#2166AC), ssp585 = red (#B2182B)
#   Ribbons: same colors, alpha ~30% transparency
#   Baseline: dotted horizontal line (projected baseline mean, delta=0 period)
#   Vertical refs: grey lines at 2050 and 2080 (unlabeled in panels;
#                  labeled once in outer top margin)
#   Y-axis: log10 CHLa, shared; labels on left column (Early July) only
#   X-axis: year 2000–2090; labels on bottom row (FH) only; ticks all panels
#   Panel label: date point label top of each column; site label left of each row
#   Boxes: closed (bty = "o")
#
# AUTHOR: [Rafa]
# =============================================================================


# -----------------------------------------------------------------------------
# 0. USER PARAMETERS
# -----------------------------------------------------------------------------

site_order       <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")
date_point_order <- c("early_july", "august", "mid_september")
date_point_labels <- c("Early July", "August", "Mid-September")

year_range <- c(2000, 2090)
baseline_anchor_year <- 2010   # matches Script 10 — end of delta=0 period

# Colors
col_45     <- "#2166AC"   # ssp245 line
col_85     <- "#B2182B"   # ssp585 line
col_45_rgb <- col2rgb(col_45) / 255
col_85_rgb <- col2rgb(col_85) / 255
alpha_ribbon <- 0.25       # ribbon transparency

# Y-axis: shared log10 scale
y_lim <- c(1.0, 2.5)
y_at  <- seq(1.0, 2.5, by = 0.25)

# Line widths
lwd_line <- 1.6
lwd_base <- 1.0

# Figure dimensions — portrait (7 rows × 3 cols)
fig_width  <- 9    # inches
fig_height <- 14   # inches

# Output paths
out_pdf <- "4_products/fig2_full_trajectory.pdf"
out_png <- "4_products/fig2_full_trajectory.png"


# -----------------------------------------------------------------------------
# 1. LOAD INPUTS
# -----------------------------------------------------------------------------

cat("Loading inputs...\n")

proj <- read.csv("2_incremental/projections_monthly.csv", stringsAsFactors = FALSE)
cat("  Projection rows:", nrow(proj), "\n")


# -----------------------------------------------------------------------------
# 2. HELPER: TRANSPARENT COLOR FROM RGB
# -----------------------------------------------------------------------------

make_transp <- function(rgb_vec, alpha) {
  rgb(rgb_vec[1], rgb_vec[2], rgb_vec[3], alpha = alpha)
}

col_45_t <- make_transp(col_45_rgb, alpha_ribbon)
col_85_t <- make_transp(col_85_rgb, alpha_ribbon)


# -----------------------------------------------------------------------------
# 3. HELPER: POLYGON RIBBON
# -----------------------------------------------------------------------------

draw_ribbon <- function(x, y_lo, y_hi, col_transp) {
  # Remove NAs
  keep <- !is.na(x) & !is.na(y_lo) & !is.na(y_hi)
  x    <- x[keep]
  y_lo <- y_lo[keep]
  y_hi <- y_hi[keep]
  if (length(x) < 2) return(invisible(NULL))
  
  polygon(c(x, rev(x)),
          c(y_hi, rev(y_lo)),
          col    = col_transp,
          border = NA)
}


# -----------------------------------------------------------------------------
# 4. SINGLE PANEL PLOT FUNCTION
# -----------------------------------------------------------------------------

plot_panel <- function(site_i, dp_label,
                       show_yaxis, show_xaxis,
                       show_site_label, show_dp_label) {
  
  # Subset data
  pd <- proj[proj$site == site_i & proj$date_point == dp_label, ]
  pd <- pd[order(pd$year), ]
  
  pd_45 <- pd[pd$scenario == "ssp245", ]
  pd_85 <- pd[pd$scenario == "ssp585", ]
  
  # Baseline mean (delta = 0 period: years <= baseline_anchor_year)
  bl_45 <- mean(pd_45$pred_logCHLa_med[pd_45$year <= baseline_anchor_year],
                na.rm = TRUE)
  
  # Empty plot — closed box
  plot(NA,
       xlim = year_range,
       ylim = y_lim,
       xaxt = "n",
       yaxt = "n",
       xlab = "",
       ylab = "",
       bty  = "o")
  
  # Subtle grid
  abline(h   = y_at, col = "grey93", lwd = 0.5)
  
  # Vertical reference lines at 2050 and 2080
  abline(v   = c(2050, 2080),
         col = "grey75",
         lwd = 0.8,
         lty = 2)
  
  # Baseline dotted horizontal
  abline(h   = bl_45,
         lty = 3,
         lwd = lwd_base,
         col = "grey45")
  
  # --- ssp245 ribbon + line ---
  draw_ribbon(pd_45$year,
              pd_45$pred_logCHLa_lo,
              pd_45$pred_logCHLa_hi,
              col_45_t)
  lines(pd_45$year, pd_45$pred_logCHLa_med,
        col = col_45, lwd = lwd_line)
  
  # --- ssp585 ribbon + line ---
  draw_ribbon(pd_85$year,
              pd_85$pred_logCHLa_lo,
              pd_85$pred_logCHLa_hi,
              col_85_t)
  lines(pd_85$year, pd_85$pred_logCHLa_med,
        col = col_85, lwd = lwd_line)
  
  # Y-axis (left column only — Early July panels)
  if (show_yaxis) {
    axis(2,
         at       = y_at,
         labels   = formatC(y_at, digits = 2, format = "f"),
         las      = 2,
         cex.axis = 0.60,
         col.axis = "grey20",
         tck      = -0.025)
  } else {
    axis(2,
         at     = y_at,
         labels = FALSE,
         tck    = -0.02)
  }
  
  # X-axis (bottom row only — FH panels; ticks all rows)
  if (show_xaxis) {
    x_ticks <- seq(2000, 2090, by = 20)
    axis(1,
         at       = x_ticks,
         labels   = x_ticks,
         cex.axis = 0.60,
         col.axis = "grey20",
         tck      = -0.025)
  } else {
    axis(1,
         at     = seq(2000, 2090, by = 20),
         labels = FALSE,
         tck    = -0.02)
  }
  
  # Date point label — top of each column (top row only)
  if (show_dp_label) {
    mtext(date_point_labels[match(dp_label, date_point_order)],
          side = 3,
          line = 0.4,
          cex  = 0.75,
          font = 2,
          col  = "grey20")
  }
  
  # Site label — left of each row (left column only)
  if (show_site_label) {
    site_lab <- if (site_i == "BM") "BM*" else site_i
    mtext(site_lab,
          side = 2,
          line = 3.2,
          cex  = 0.75,
          font = 2,
          col  = "grey20",
          las  = 3)
  }
}


# -----------------------------------------------------------------------------
# 5. RENDER FIGURE
# -----------------------------------------------------------------------------

render_figure <- function(file_path, device_fn, ...) {
  
  device_fn(file_path, ...)
  
  n_dp   <- length(date_point_order)
  n_site <- length(site_order)
  
  # Layout: n_site rows × n_dp cols (sites as rows, date points as columns)
  par(mfrow  = c(n_site, n_dp),
      mar    = c(1.5, 2.2, 1.5, 0.5),
      oma    = c(4.0, 5.0, 3.5, 1.0))
  
  # Fill panels row by row (site by site, across date points)
  for (si in seq_len(n_site)) {
    site_i <- site_order[si]
    
    for (di in seq_len(n_dp)) {
      dp <- date_point_order[di]
      
      show_yaxis      <- di == 1            # left column: Early July
      show_xaxis      <- si == n_site       # bottom row: FH
      show_dp_label   <- si == 1            # date point label: top row only
      show_site_label <- di == 1            # site label: left column only
      
      plot_panel(site_i          = site_i,
                 dp_label        = dp,
                 show_yaxis      = show_yaxis,
                 show_xaxis      = show_xaxis,
                 show_site_label = show_site_label,
                 show_dp_label   = show_dp_label)
    }
  }
  
  # Shared y-axis label
  mtext(expression(log[10]~CHLa~(mg~m^{-2})),
        side  = 2,
        outer = TRUE,
        line  = 3.2,
        cex   = 0.85,
        col   = "grey20")
  
  # Shared x-axis label
  mtext("Year",
        side  = 1,
        outer = TRUE,
        line  = 2.5,
        cex   = 0.85,
        col   = "grey20")
  
  # Main title
  mtext("Projected bloom trajectory by site and growing season period",
        side  = 3,
        outer = TRUE,
        line  = 2.0,
        cex   = 0.95,
        font  = 2,
        col   = "grey15")
  
  # Horizon reference labels — placed above top row panels
  # 2050 ≈ 56% along x-axis (2000–2090); 2080 ≈ 89%
  # These are proportional positions within the plot area (not outer margin)
  mtext("2050",
        side  = 3,
        outer = TRUE,
        line  = 0.5,
        at    = 0.39,   # adjusted for 3-col layout
        cex   = 0.62,
        col   = "grey50")
  mtext("2080",
        side  = 3,
        outer = TRUE,
        line  = 0.5,
        at    = 0.73,
        cex   = 0.62,
        col   = "grey50")
  
  # Legend — horizontal, in outer bottom margin
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
      mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom",
         legend  = c("RCP4.5 (ssp245)", "RCP8.5 (ssp585)",
                     "ESM q25\u2013q75 range", "Projected baseline"),
         col     = c(col_45, col_85, "grey60", "grey45"),
         lty     = c(1, 1, NA, 3),
         lwd     = c(lwd_line, lwd_line, NA, lwd_base),
         fill    = c(NA, NA, "grey80", NA),
         border  = c(NA, NA, NA, NA),
         pt.cex  = 1.2,
         bty     = "n",
         horiz   = TRUE,
         cex     = 0.78,
         x.intersp = 0.6,
         xpd     = TRUE,
         inset   = c(0, 0.01))
  
  # BM footnote
  mtext("* BM: below Missoula WWTP (point-source site; hydrological predictors do not capture effluent dynamics)",
        side  = 1,
        outer = TRUE,
        line  = 3.2,
        cex   = 0.58,
        adj   = 0,
        col   = "grey45")
  
  dev.off()
  cat("  Written:", file_path, "\n")
}

cat("Rendering figures...\n")

render_figure(out_pdf, pdf,
              width       = fig_width,
              height      = fig_height,
              useDingbats = FALSE)

render_figure(out_png, png,
              width  = fig_width,
              height = fig_height,
              units  = "in",
              res    = 300)

cat("Done.\n")