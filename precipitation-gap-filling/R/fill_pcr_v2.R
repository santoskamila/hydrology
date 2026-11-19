fill_pcr_v2 <- function(target_id, all_monthly, neighbours, gap_year, gap_month) {
  # Principal Component Regression (PCR) with classical PCA (prcomp, the
  # equivalent of MATLAB pca()), trained on complete months only.

  ids <- names(neighbours)
  T <- join_neighbours(all_monthly[[target_id]], ids, all_monthly)
  T_clean <- T[, -(1:2), drop = FALSE]
  train <- T_clean[complete.cases(T_clean), , drop = FALSE]

  Y <- train$MonthlyTotal
  X <- as.matrix(train[, -1, drop = FALSE])

  if (nrow(train) < (length(ids) + 2)) return(-997)

  mean_Y <- mean(Y)
  sd_Y <- sd(Y)
  Y_std <- if (sd_Y > 0) (Y - mean_Y) / sd_Y else Y - mean_Y

  pc <- prcomp(X, center = TRUE, scale. = FALSE)
  coef <- lsq(cbind(1, pc$x), Y_std)

  i_gap <- which(T$Year == gap_year & T$Month == gap_month)[1]
  x_gap <- as.numeric(T_clean[i_gap, -1])
  if (any(is.na(x_gap))) return(-998)

  scores_gap <- (x_gap - pc$center) %*% pc$rotation
  Px <- sum(c(1, scores_gap) * coef) * sd_Y + mean_Y

  if (Px < 0) Px <- 0
  Px
}
