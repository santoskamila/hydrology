# =============================================================================
# Shared helpers (R equivalents of MATLAB table/matrix operations)
# Used by all gap-filling functions.
# Monthly tables have the columns: Year, Month, GroupCount, MonthlyTotal
# =============================================================================

# MonthlyTotal of a table for a given year/month (numeric(0) if absent)
month_value <- function(df, year, month) {
  df$MonthlyTotal[which(df$Year == year & df$Month == month)]
}

# NA-safe checks (in MATLAB, NaN >= 0 is false)
is_valid_nonneg <- function(v) length(v) > 0 && !is.na(v[1]) && v[1] >= 0
is_valid_value  <- function(v) length(v) > 0 && !is.na(v[1])

# Equivalent of outerjoin(..., 'Keys', {'Year','Month'}, 'MergeKeys', true)
# between the target gauge and a set of neighbours. Returns Year, Month,
# MonthlyTotal and one MonthlyTotal_<id> column per neighbour, sorted by date.
join_neighbours <- function(target, ids, data) {
  T <- target[, c("Year", "Month", "MonthlyTotal")]
  for (id in ids) {
    d <- data[[id]][, c("Year", "Month", "MonthlyTotal")]
    names(d)[3] <- paste0("MonthlyTotal_", id)
    T <- merge(T, d, by = c("Year", "Month"), all = TRUE)
  }
  T <- T[order(T$Year, T$Month), , drop = FALSE]
  rownames(T) <- NULL
  T
}

# Equivalent of innerjoin of two monthly tables (columns _x and _y)
join_pair <- function(a, b) {
  m <- merge(a[, c("Year", "Month", "MonthlyTotal")],
             b[, c("Year", "Month", "MonthlyTotal")],
             by = c("Year", "Month"), suffixes = c("_x", "_y"))
  m[order(m$Year, m$Month), , drop = FALSE]
}

# Least squares (equivalent of MATLAB X \ y and regress).
# Pivoted QR with MATLAB's rank tolerance; redundant columns (rank
# deficiency) get a zero coefficient ("basic solution").
lsq <- function(X, y) {
  X <- as.matrix(X)
  q <- qr(X, LAPACK = TRUE)
  R <- qr.R(q)
  d <- abs(diag(R))
  b <- numeric(ncol(X))
  if (length(d) == 0 || d[1] == 0) return(b)
  tol <- max(dim(X)) * .Machine$double.eps * d[1]
  r <- sum(d > tol)
  qty <- qr.qty(q, y)[seq_len(r)]
  b[q$pivot[seq_len(r)]] <- backsolve(R[seq_len(r), seq_len(r), drop = FALSE], qty)
  b
}

# Equivalent of polyfit(x, y, 1): returns c(slope, intercept)
polyfit1 <- function(x, y) {
  b <- lsq(cbind(1, x), y)
  c(b[2], b[1])
}

# MATLAB 'lines' palette (7 colours, recycled)
matlab_lines <- function(n) {
  base <- c("#0072BD", "#D95319", "#EDB120", "#7E2F8E",
            "#77AC30", "#4DBEEE", "#A2142F")
  rep(base, length.out = n)
}

# MATLAB 'jet' palette
matlab_jet <- function(n = 64) {
  grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                "#7FFF7F", "yellow", "#FF7F00", "red",
                                "#7F0000"))(n)
}
