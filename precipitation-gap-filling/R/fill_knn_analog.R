fill_knn_analog <- function(target_id, all_monthly, neighbours, gap_year, gap_month, k = 5) {
  # k-NN ANALOGUES
  # Finds the k most similar historical months (normalised RMSE distance over
  # the neighbouring gauges) and takes an inverse-distance weighted mean.

  ids <- names(neighbours)

  # 1. Full data matrix (precipitation columns only)
  T <- join_neighbours(all_monthly[[target_id]], ids, all_monthly)
  X_hist <- as.matrix(T[, -(1:3), drop = FALSE])
  Y_hist <- T$MonthlyTotal

  # 2. Predictor vector for the gap month
  i_gap <- which(T$Year == gap_year & T$Month == gap_month)[1]
  x_gap <- X_hist[i_gap, ]

  valid_gap <- which(!is.na(x_gap))
  if (length(valid_gap) == 0) return(-998)

  min_common <- floor(length(valid_gap) * 0.5)
  if (min_common == 0) min_common <- 1

  # 3. NA-robust normalised distances
  n_hist <- nrow(X_hist)
  dist <- rep(Inf, n_hist)
  for (i in seq_len(n_hist)) {
    if (T$Year[i] == gap_year && T$Month[i] == gap_month) next
    xi <- X_hist[i, ]
    common <- which(!is.na(xi) & !is.na(x_gap))
    if (length(common) >= min_common) {
      dist[i] <- sqrt(sum((xi[common] - x_gap[common])^2) / length(common))
    }
  }

  # 4. k nearest neighbours
  ord <- order(dist)                 # stable sort, Inf at the end
  sorted_dist <- dist[ord]
  valid_Y <- which(!is.na(Y_hist[ord]))
  if (length(valid_Y) < k) return(-996)
  top <- valid_Y[1:k]

  top_values <- Y_hist[ord[top]]
  top_dist <- sorted_dist[top]

  # 5. Weighted mean
  epsilon <- 1e-6
  w <- 1 / (top_dist + epsilon)
  if (any(top_dist < epsilon)) {
    Px <- top_values[top_dist < epsilon][1]
  } else {
    Px <- sum(w * top_values) / sum(w)
  }

  if (is.na(Px)) Px <- -995
  if (Px < 0) Px <- 0
  Px
}
