fill_knn_simple <- function(target_id, all_monthly, k, dist_matrix, station_names, gap_year, gap_month) {
  # Fills one gap with the logic of R's impute.knn (simple mean of the k
  # closest gauges), using a pre-computed distance matrix.
  #
  # INPUTS:
  #   target_id     - gauge with the gap
  #   all_monthly   - list with all monthly tables
  #   k             - number of neighbouring gauges to use
  #   dist_matrix   - (n x n) distance matrix between gauges
  #   station_names - names of all gauges (for indexing)
  #   gap_year, gap_month - date of the gap

  i_target <- which(station_names == target_id)
  if (length(i_target) == 0) return(-999)  # gauge not found

  d_target <- dist_matrix[, i_target]

  # Gauges with valid data in the gap month
  idx <- c(); vals <- c()
  for (i in seq_along(station_names)) {
    if (i == i_target) next
    v <- month_value(all_monthly[[station_names[i]]], gap_year, gap_month)
    if (is_valid_nonneg(v)) { idx <- c(idx, i); vals <- c(vals, v[1]) }
  }

  if (length(idx) == 0) return(-998)  # no neighbour with valid data

  ord <- order(d_target[idx])
  n_use <- min(k, length(ord))
  if (n_use == 0) return(-997)

  Px <- mean(vals[ord[1:n_use]])  # simple mean
  if (Px < 0) Px <- 0
  Px
}
