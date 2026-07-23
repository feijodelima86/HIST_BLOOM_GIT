# ============================================================================
# 13c_driver_decomposition.R
# Per-site driver decomposition of the projected end-century bloom change.
# M1 is additive, so predict(type="terms") partitions logCHLa exactly into
# per-smooth contributions. We take each term's change from an early window
# (2026-2035) to a late window (2089-2098) and show its share of |total change|.
# High bracket. lag_y = persistence/inherited, not a climate forcing.
#
# Input : 2_incremental/bloom_projections_members.csv ; 3_models/bloom_model_M1.rds
# Output: 4_products/driver_decomposition.pdf
# ============================================================================

library(mgcv)

m1   <- readRDS("3_models/bloom_model_M1.rds")
grid <- read.csv("2_incremental/bloom_projections_members.csv", stringsAsFactors = FALSE)
grid$Site <- factor(grid$Site, levels = levels(m1$model$Site))

# exact additive term contributions
tt <- predict(m1, newdata = grid, type = "terms")
drv <- c("s(lag_y)", "s(anomaly)", "s(logQ_obs_cfs)",
         "s(Days_Since_Freshet)", "s(logTP_mg_L)", "s(Temp_oC)")
tt <- tt[, drv]
colnames(tt) <- c("lag", "anomaly", "logQ", "DSF", "logTP", "Temp")

g <- data.frame(Site = grid$Site, bracket = grid$bracket, year = grid$year, tt)
g <- g[g$bracket == "high", ]

# early vs late window means, then per-site delta per term
early <- g[g$year <= 2035, ]
late  <- g[g$year >= 2089, ]
site_order <- c("DL", "GR", "BN", "MS", "BM", "HU", "FH")
delta <- sapply(site_order, function(st)
  colMeans(late[late$Site == st, 4:9]) - colMeans(early[early$Site == st, 4:9]))

pct <- apply(abs(delta), 2, function(x) 100 * x / sum(x))   # drivers x sites, cols sum to 100

cols <- c(lag = "grey75", anomaly = "#1b9e77", logQ = "#377EB8",
          DSF = "#7570b3", logTP = "#d95f02", Temp = "#e7298a")

if (!dir.exists("4_products")) dir.create("4_products", recursive = TRUE)
pdf("4_products/driver_decomposition.pdf", width = 8, height = 5)
par(mar = c(4, 4.2, 3, 6.5), xpd = TRUE)
barplot(pct, col = cols[rownames(pct)], border = "white", las = 1,
        ylab = "% of projected end-century change (|change| share)",
        main = "Driver decomposition by site (high bracket)")
legend(x = ncol(pct) * 1.25, y = 100, legend = rownames(pct),
       fill = cols[rownames(pct)], bty = "n", cex = 0.9)
dev.off()
cat("Wrote 4_products/driver_decomposition.pdf\n")


m <- read.csv("2_incremental/bloom_projections_members.csv")

# spread across members, per site x bracket x year
v <- aggregate(pred_logCHLa ~ Site + bracket + year, data = m, FUN = sd)

cols <- c(DL="#1b9e77", GR="#d95f02", BN="#7570b3", MS="#e7298a",
          BM="#66a61e", HU="#e6ab02", FH="#a6761d")

plot(NA, xlim = range(v$year), ylim = range(v$pred_logCHLa, na.rm = TRUE),
     xlab = "Year", ylab = "SD of logCHLa across members")

for (st in names(cols))
  for (b in c("low","high")) {
    d <- v[v$Site == st & v$bracket == b, ]
    d <- d[order(d$year), ]
    lines(d$year, d$pred_logCHLa, col = cols[st], lty = if (b=="high") 1 else 2, lwd = 2)
  }

legend("topright", names(cols), col = cols, lwd = 2, bty = "n", cex = 0.8)
