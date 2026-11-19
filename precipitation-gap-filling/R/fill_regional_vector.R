fill_regional_vector <- function(target_id, all_monthly, neighbours, gap_year, gap_month, monthly_means) {
  # Estimates precipitation (Px) at the target gauge by the Regional Vector
  # method, based on deviations in anomalous months.

  ids <- names(neighbours)
  target <- all_monthly[[target_id]]

  # 1. Deviations
  dates <- unique(target[, c("Year", "Month")])
  dates <- dates[order(dates$Year, dates$Month), ]

  # Target deviation from its monthly mean
  target_dev <- target$MonthlyTotal - monthly_means[[target_id]][target$Month]

  # Mean deviation of the neighbours in a given month
  neighbour_mean_dev <- function(year, month) {
    s <- 0; n <- 0
    for (id in ids) {
      v <- month_value(neighbours[[id]], year, month)
      if (is_valid_nonneg(v)) {
        s <- s + (v[1] - monthly_means[[id]][month])
        n <- n + 1
      }
    }
    if (n > 0) s / n else NA  # NA: no neighbour in this month
  }
  neigh_dev <- mapply(neighbour_mean_dev, dates$Year, dates$Month)

  # 2. Anomalous months (target and neighbours both valid)
  valid <- !is.na(target$MonthlyTotal) & target$MonthlyTotal >= 0 & !is.na(neigh_dev)
  threshold <- sd(neigh_dev[valid])  # threshold = 1 standard deviation
  anomalous <- valid & !is.na(abs(neigh_dev) > threshold) & abs(neigh_dev) > threshold

  # 3. Linear regression on the anomalous months
  if (sum(anomalous) > 5) {
    x <- neigh_dev[anomalous]   # independent variable
    y <- target_dev[anomalous]  # dependent variable
    coef <- polyfit1(x, y)

    # 4. Estimate the gap
    dv_gap <- neighbour_mean_dev(gap_year, gap_month)
    if (!is.na(dv_gap)) {
      Px <- monthly_means[[target_id]][gap_month] + coef[1] * dv_gap + coef[2]
    } else {
      Px <- -998  # no neighbour with data in the gap month
    }
  } else {
    Px <- -997    # not enough anomalous months
  }

  # NOTE: as in the original MATLAB code, error codes also become 0 here.
  if (Px < 0) Px <- 0
  unname(Px)
}
