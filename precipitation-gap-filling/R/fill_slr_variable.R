fill_slr_variable <- function(target_id, all_monthly, neighbours, gap_year, gap_month) {
  # "Simple Linear Regression (SLR) with Variable Donor Sites".
  # For each gap, picks the most correlated neighbour and applies a simple
  # linear regression. Same logic as fill_linear_regression, but filters on
  # "not NA" instead of ">= 0".

  target <- all_monthly[[target_id]]
  best_id <- ""
  best_r <- -1

  # 1. Most correlated neighbour (the "variable donor")
  for (id in names(neighbours)) {
    dc <- join_pair(target, neighbours[[id]])
    dc <- dc[!is.na(dc$MonthlyTotal_x) & !is.na(dc$MonthlyTotal_y), ]
    if (nrow(dc) > 12) {
      r <- suppressWarnings(cor(dc$MonthlyTotal_x, dc$MonthlyTotal_y))
      if (isTRUE(r > best_r)) { best_r <- r; best_id <- id }
    }
  }

  # 2. Regression on the best neighbour
  if (best_id != "") {
    best <- neighbours[[best_id]]
    dc <- join_pair(target, best)
    dc <- dc[!is.na(dc$MonthlyTotal_x) & !is.na(dc$MonthlyTotal_y), ]
    coef <- polyfit1(dc$MonthlyTotal_y, dc$MonthlyTotal_x)

    Pi <- month_value(best, gap_year, gap_month)
    Px <- if (is_valid_value(Pi)) coef[1] * Pi[1] + coef[2] else -998
  } else {
    Px <- -997
  }

  # NOTE: as in the original MATLAB code, error codes also become 0 here.
  if (Px < 0) Px <- 0
  Px
}
