taylor_diagram <- function(p_input, r_input, markers = 16, colors = NULL,
                           title = "Taylor Diagram", colobs = "black",
                           markersize = 1.4, showlabels = TRUE,
                           labels = c("Observed"), collabels = "black",
                           fontsize = 0.7, max_factor = 1.15) {
  # Draws a Taylor diagram with base R graphics.
  #   p_input - standard deviations (1st = observed, others = models)
  #   r_input - correlations        (1st = observed (1), others = models)
  # Invisibly returns the model colours, to build a legend.

  p0 <- p_input[1]
  r <- r_input[-1]; p <- p_input[-1]
  if (is.null(colors)) colors <- matlab_lines(length(p))
  colors <- rep(colors, length.out = length(p))
  label_obs <- labels[1]
  labels_mod <- labels[-1]
  if (length(labels_mod) != length(p)) labels_mod <- as.character(seq_along(p))

  max_p <- max_factor * max(c(p0, p), na.rm = TRUE)
  if (!is.finite(max_p) || max_p == 0) max_p <- 1

  op <- par(mar = c(3, 1, 3, 1)); on.exit(par(op))
  plot(NA, xlim = c(0, max_p), ylim = c(0, max_p), asp = 1,
       axes = FALSE, xlab = "", ylab = "", main = title)

  circle <- function(x, y, radius, ...) {
    th <- seq(0, 2 * pi, length.out = 101)
    lines(radius * cos(th) + x, radius * sin(th) + y, ...)
  }

  # Radial correlation lines
  tick_vals <- c(0, seq(0.2, 0.8, 0.2), 0.9, 0.95, 0.99, 1)
  ang <- acos(tick_vals)
  segments(0, 0, max_p * cos(ang), max_p * sin(ang), lty = 3)
  text(max_p * cos(ang), max_p * sin(ang), tick_vals, cex = 0.6, pos = 3, xpd = TRUE)
  a07 <- acos(0.7) - pi / 16
  text(max_p * cos(a07) * 1.06, max_p * sin(a07) * 1.06, "Correlation",
       srt = -45, cex = 0.75, xpd = TRUE)

  # Standard deviation arcs
  p_ticks <- seq(0, max_p, length.out = 6)
  for (pt in p_ticks[-1]) circle(0, 0, pt, lty = 3)
  text(p_ticks[-1], -0.03 * max_p, signif(p_ticks[-1], 3), cex = 0.6, xpd = TRUE)
  text(0.35 * max_p, -0.1 * max_p, "Standard Deviation", cex = 0.75, xpd = TRUE)

  # Centred RMSE arcs (around the observed point)
  for (pt in p_ticks[-1]) if (pt <= p0) circle(p0, 0, pt, col = "grey60", lty = 2)
  text(p0, -0.1 * max_p, "Centred RMSE", col = "grey50", cex = 0.75, xpd = TRUE)

  # Model points
  ok <- !is.na(p) & !is.na(r) & abs(r) <= 1
  x <- p * r; y <- p * sin(acos(pmin(pmax(r, -1), 1)))
  points(x[ok], y[ok], pch = markers, col = colors[ok], cex = markersize)

  # Observed point
  points(p0, 0, pch = 21, bg = "white", col = colobs, cex = markersize + 0.3, lwd = 1.5, xpd = TRUE)
  text(p0, 0.02 * max_p, paste0(" ", label_obs), col = colobs, font = 2, cex = 0.75, pos = 3, xpd = TRUE)

  if (isTRUE(showlabels)) {
    text(x[ok], y[ok], labels_mod[ok], col = collabels, cex = fontsize, pos = 2)
  }

  invisible(colors)
}
