fill_arithmetic_mean <- function(target_id, all_monthly, neighbours, gap_year, gap_month) {
  # Fills a gap with the simple arithmetic mean of the given neighbouring gauges.

  ids <- names(neighbours)
  if (length(ids) == 0) return(-998)  # no neighbours available

  values <- c()
  for (id in ids) {
    v <- month_value(neighbours[[id]], gap_year, gap_month)
    if (is_valid_value(v)) values <- c(values, v[1])
  }

  if (length(values) == 0) return(-997)  # no neighbour with a valid value

  mean(values)
}
