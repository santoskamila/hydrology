fill_pcr_nipals <- function(target_id, all_monthly, neighbours, gap_year, gap_month, num_components = 3) {
  # Gap filling by PCR with NIPALS
  # (regression without intercept; the gap month is left out of training)

  ids <- names(neighbours)

  # 1-2. Data matrix (Year, Month, MonthlyTotal, predictors)
  T <- join_neighbours(all_monthly[[target_id]], ids, all_monthly)

  # 3. Gap row
  is_gap <- T$Year == gap_year & T$Month == gap_month
  x_gap <- as.numeric(T[which(is_gap)[1], -(1:3)])

  # Enough data for the regression?
  if (sum(!is.na(x_gap)) < 2) return(-998)

  # 4. Historical base (without the gap row)
  X <- as.matrix(T[!is_gap, -(1:3), drop = FALSE])
  Y <- T$MonthlyTotal[!is_gap]

  # Rows with at least 2 predictors
  keep <- rowSums(!is.na(X)) >= 2
  X <- X[keep, , drop = FALSE]
  Y <- Y[keep]

  # Safety check before PCA
  if (nrow(X) < num_components || ncol(X) < 2) return(-997)

  Px <- tryCatch({
    np <- nipals_pca(X, num_components)
    # regress(Y, T_scores): no intercept, rows with NA are ignored
    ok <- !is.na(Y) & complete.cases(np$T)
    B <- lsq(np$T[ok, , drop = FALSE], Y[ok])
    as.numeric((x_gap %*% np$P) %*% B)   # NA if any neighbour is missing
  }, error = function(e) -996)           # internal error

  if (is.na(Px)) Px <- -995
  Px
}
