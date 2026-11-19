fill_linear_regression <- function(target_id, all_monthly, neighbours, year, month) {
  # Estimates precipitation (Px) by simple linear regression on the most
  # correlated neighbouring gauge (only non-negative months are used).

  target <- all_monthly[[target_id]]
  best_id <- ""
  best_r <- -1

  # 1. Most correlated neighbour
  for (id in names(neighbours)) {
    dc <- join_pair(target, neighbours[[id]])
    dc <- dc[which(dc$MonthlyTotal_x >= 0 & dc$MonthlyTotal_y >= 0), ]
    if (nrow(dc) > 12) {
      r <- suppressWarnings(cor(dc$MonthlyTotal_x, dc$MonthlyTotal_y))
      if (isTRUE(r > best_r)) { best_r <- r; best_id <- id }
    }
  }

  # 2. Regression on the best neighbour
  if (best_id != "") {
    best <- neighbours[[best_id]]
    dc <- join_pair(target, best)
    dc <- dc[which(dc$MonthlyTotal_x >= 0 & dc$MonthlyTotal_y >= 0), ]
    coef <- polyfit1(dc$MonthlyTotal_y, dc$MonthlyTotal_x)

    Pi <- month_value(best, year, month)
    Px <- if (is_valid_nonneg(Pi)) coef[1] * Pi[1] + coef[2] else -998  # best neighbour also missing
  } else {
    Px <- -999  # no correlated neighbour
  }

  # NOTE: as in the original MATLAB code, this line also turns the error
  # codes (-998/-999) into 0.
  if (Px < 0) Px <- 0
  Px
}
