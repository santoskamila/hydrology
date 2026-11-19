normalize_daily_units <- function(P) {
  # Normalises units when the column comes scaled (e.g. 18500000 -> 185)
  P <- as.numeric(P)
  x <- P[!is.na(P)]
  if (length(x) == 0) return(P)

  # 99th percentile (type = 5 reproduces MATLAB prctile)
  p99 <- unname(quantile(x, 0.99, type = 5))
  if (!is.finite(p99) || p99 == 0) return(P)

  # Try divisors until p99 <= 500
  for (d in c(10, 100, 1000, 10000, 100000, 1000000)) {
    if (p99 / d <= 500) { P <- P / d; break }
  }

  # Extra check: if the mean is still implausible, divide further
  if (mean(P, na.rm = TRUE) > 500) P <- P / 1000

  # No negatives
  P[!is.na(P) & P < 0] <- NA
  P
}
