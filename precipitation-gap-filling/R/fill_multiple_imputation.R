fill_multiple_imputation <- function(target_id, all_monthly, neighbours, gap_year, gap_month) {
  # Fills a gap by Multiple Imputation with Stochastic Regression
  # (MLR in log space + Gaussian noise).

  m_imputations <- 5  # number of imputations
  log_const <- 1      # avoids log(0)

  target <- all_monthly[[target_id]]
  ids <- names(neighbours)
  n_avail <- length(ids)
  if (n_avail < 1) return(-998)

  # 1. Base model (MLR on the given neighbours)
  T <- join_neighbours(target, ids, all_monthly)
  T <- T[, -(1:2), drop = FALSE]
  train <- T[complete.cases(T), , drop = FALSE]
  if (nrow(train) < (n_avail + 2)) return(-997)

  Y_log <- log(train$MonthlyTotal + log_const)
  X_reg_log <- cbind(1, log(as.matrix(train[, -1, drop = FALSE]) + log_const))
  coef_log <- tryCatch(lsq(X_reg_log, Y_log), error = function(e) NULL)
  if (is.null(coef_log)) return(-996)

  # 2. Stochastic imputations
  vals <- sapply(ids, function(id) month_value(neighbours[[id]], gap_year, gap_month)[1])
  y_hat_log <- sum(c(1, log(vals + log_const)) * coef_log)

  sigma_e <- sd(Y_log - X_reg_log %*% coef_log)  # residual standard deviation
  imputed_log <- y_hat_log + rnorm(m_imputations) * sigma_e

  # 3. Back-transform and aggregate
  Px <- mean(exp(imputed_log) - log_const)
  if (Px < 0) Px <- 0
  Px
}
