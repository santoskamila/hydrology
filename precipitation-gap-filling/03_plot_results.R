# =============================================================================
# STEP 3: PLOTS OF THE VALIDATION RESULTS
# Requires: 02_validate_methods.R (results/validation_results.rds)
# Output:   PNG files in figures/
# Packages: ggplot2, gridExtra, igraph
# =============================================================================
rm(list = ls())
library(ggplot2); library(gridExtra); library(igraph)
source("R/utils.R")
source("R/taylor_diagram.R")

# ---- 1. Load -------------------------------------------------------------------------
if (!file.exists("results/validation_results.rds")) stop("File not found. Run 02_validate_methods.R first.")
V <- readRDS("results/validation_results.rds")
results <- V$results; gap_locations <- V$gap_locations
observed <- V$observed; predicted <- V$predicted
methods <- V$methods; station_names <- V$station_names; X_complete <- V$X_complete

out_dir <- "figures"
dir.create(out_dir, showWarnings = FALSE)
n_methods <- length(methods)
n_st <- length(station_names)
metrics <- c("RMSE", "MAE", "R2", "NSE", "Spearman", "Willmott_d")
method_colors <- matlab_lines(n_methods)
lab <- function(x) gsub("_", " ", x)

# ---- 2. Heatmap panel -----------------------------------------------------------------
cat("--- Heatmap panel ---\n")
plots <- lapply(metrics, function(met) {
  df <- do.call(rbind, lapply(station_names, function(st)
    data.frame(Gauge = st, Method = lab(methods), Value = results[[st]][[met]])))
  df$Gauge <- factor(df$Gauge, levels = station_names)
  df$Method <- factor(df$Method, levels = rev(lab(methods)))
  ggplot(df, aes(Gauge, Method, fill = Value)) +
    geom_tile(color = "white") +
    # cell values only when there are few gauges (otherwise unreadable)
    (if (n_st <= 8) geom_text(aes(label = formatC(Value, format = "g", digits = 3)), size = 2.5)) +
    scale_fill_gradientn(colours = matlab_jet(), name = NULL) +
    labs(title = met, x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"))
})
png(file.path(out_dir, "heatmap_panel.png"), width = 1600, height = 1000, res = 100)
grid.arrange(grobs = plots, ncol = 3,
             top = grid::textGrob("Performance Metrics Heatmaps",
                                  gp = grid::gpar(fontsize = 16, fontface = "bold")))
dev.off()

# ---- 3. Scatter panel (one subplot per method) ------------------------------------
cat("--- Scatter panel ---\n")
nr <- floor(sqrt(n_methods)); nc <- ceiling(n_methods / nr)
png(file.path(out_dir, "scatter_by_method.png"), width = 1200, height = 900, res = 100)
par(mfrow = c(nr, nc), oma = c(0, 0, 3, 0), mar = c(4, 4, 2, 1))
for (m in seq_len(n_methods)) {
  pred <- predicted[, m]
  v <- !is.na(pred) & pred >= 0 & !is.na(observed)
  lim <- max(c(observed[v], pred[v], 1))
  plot(observed[v], pred[v], pch = 16, col = rgb(0, 0, 1, 0.3), cex = 0.8,
       xlim = c(0, lim), ylim = c(0, lim), asp = 1,
       main = lab(methods[m]), xlab = "Observed (mm)", ylab = "Predicted (mm)")
  abline(0, 1); grid()
}
mtext("Observed vs. Predicted by Method (All Gauges)", outer = TRUE, cex = 1.3, font = 2)
dev.off()

# ---- 4. Taylor diagrams -------------------------------------------------------------
cat("--- Taylor diagrams ---\n")
taylor_stats <- function(obs, pred_m) {
  st <- matrix(0, n_methods + 1, 2)
  st[1, ] <- c(sd(obs, na.rm = TRUE), 1)
  for (m in seq_len(n_methods)) {
    pred <- pred_m[, m]
    v <- !is.na(pred) & pred >= 0 & !is.na(obs)
    if (sum(v) > 1) st[m + 1, ] <- c(sd(pred[v]), suppressWarnings(cor(obs[v], pred[v])))
  }
  st
}

n_tiles <- n_st + 1                    # +1 tile for the legend
nc <- ceiling(sqrt(n_tiles * 1.6)); nr <- ceiling(n_tiles / nc)
png(file.path(out_dir, "taylor_by_gauge.png"), width = 1500, height = 900, res = 100)
par(mfrow = c(nr, nc))
for (s in seq_len(n_st)) {
  idx <- which(gap_locations[, 1] == s)
  obs <- observed[idx]
  if (sum(!is.na(obs)) < 2) { plot.new(); title(station_names[s]); next }
  st <- taylor_stats(obs, predicted[idx, , drop = FALSE])
  taylor_diagram(st[, 1], st[, 2], showlabels = FALSE, title = station_names[s],
                 markersize = 1.3, colors = method_colors)
}
plot.new()
legend("center", legend = lab(methods), col = method_colors, pch = 16, title = "Methods", bty = "n")
dev.off()

# Aggregated Taylor diagram
st_all <- taylor_stats(observed, predicted)
png(file.path(out_dir, "taylor_aggregated.png"), width = 1000, height = 800, res = 100)
layout(matrix(1:2, 1), widths = c(4, 1.3))
taylor_diagram(st_all[, 1], st_all[, 2], showlabels = FALSE,
               title = "Aggregated Taylor Diagram (All Gauges)",
               markersize = 1.5, colors = method_colors)
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", legend = lab(methods), col = method_colors, pch = 16, bty = "n")
dev.off()

# ---- 5. Spearman correlation network --------------------------------------------
cat("--- Correlation network ---\n")
corr_threshold <- 0.6
cm <- suppressWarnings(cor(X_complete, method = "spearman", use = "pairwise.complete.obs"))
adj <- cm
adj[is.na(adj) | abs(adj) < corr_threshold] <- 0
diag(adj) <- 0
dimnames(adj) <- list(station_names, station_names)
G <- graph_from_adjacency_matrix(adj, mode = "upper", weighted = TRUE, diag = FALSE)
if (ecount(G) > 0) {
  w <- E(G)$weight
  pal <- matlab_jet(100)
  rng <- range(w); if (diff(rng) == 0) rng <- rng + c(-0.01, 0.01)
  edge_col <- pal[pmax(1, ceiling((w - rng[1]) / diff(rng) * 100))]
  set.seed(1)
  png(file.path(out_dir, "correlation_network.png"), width = 900, height = 900, res = 100)
  layout(matrix(1:2, 1), widths = c(6, 1))
  par(mar = c(1, 1, 3, 1))
  plot(G, layout = layout_with_fr(G), edge.width = 5 * abs(w), edge.color = edge_col,
       vertex.color = "red", vertex.size = 5, vertex.label.color = "black",
       vertex.label.cex = 0.8, vertex.label.dist = 1,
       main = sprintf("Spearman Correlation Network (rho > %.2f)", corr_threshold))
  par(mar = c(5, 1, 5, 4))
  image(1, seq(rng[1], rng[2], length.out = 100), t(seq_len(100)), col = pal,
        axes = FALSE, xlab = "", ylab = "")
  axis(4, las = 1); mtext("Correlation strength", side = 4, line = 2.5)
  dev.off()
}
cat("Plots saved to figures/.\n")
