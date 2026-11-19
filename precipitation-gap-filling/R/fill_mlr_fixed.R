fill_mlr_fixed <- function(target_id, all_monthly, neighbours, gap_year, gap_month, num_donors = 5) {
  # Multiple Linear Regression with the N most correlated donors that are
  # available in the gap month.

  target <- all_monthly[[target_id]]
  ids <- names(neighbours)

  # 1. N best candidate donors (by correlation)
  r <- rep(-1, length(ids))
  for (i in seq_along(ids)) {
    dc <- join_pair(target, neighbours[[ids[i]]])
    dc <- dc[!is.na(dc$MonthlyTotal_x) & !is.na(dc$MonthlyTotal_y), ]
    if (nrow(dc) > 12) r[i] <- suppressWarnings(cor(dc$MonthlyTotal_x, dc$MonthlyTotal_y))
  }
  # MATLAB sort 'descend' puts NaN first and is stable
  ord <- order(-ifelse(is.na(r), Inf, r))
  candidates <- ids[ord[seq_len(min(num_donors, length(ord)))]]

  # 2. Candidates AVAILABLE in the gap month
  donors <- c(); values <- c()
  for (id in candidates) {
    v <- month_value(all_monthly[[id]], gap_year, gap_month)
    if (is_valid_value(v)) { donors <- c(donors, id); values <- c(values, v[1]) }
  }
  n_donors <- length(donors)
  if (n_donors < 1) return(-998)

  # 3. Tailored model
  T <- join_neighbours(target, donors, all_monthly)
  T <- T[, -(1:2), drop = FALSE]                # MonthlyTotal + predictors
  train <- T[complete.cases(T), , drop = FALSE]
  if (nrow(train) < (n_donors + 2)) return(-997)

  X_reg <- cbind(1, as.matrix(train[, -1, drop = FALSE]))
  coef <- tryCatch(lsq(X_reg, train$MonthlyTotal), error = function(e) NULL)
  if (is.null(coef)) return(-996)

  # 4. Prediction
  Px <- sum(c(1, values) * coef)
  if (Px < 0) Px <- 0
  Px
}
