fill_regional_weighting <- function(target_id, neighbours, year, month, annual_means) {
  # Estimates precipitation (Px) at the target gauge by Regional Weighting.
  # Px = (1/n) * SUM [ (Nx / Ni) * Pi ]

  Nx <- annual_means[[target_id]]  # long-term annual mean of the target gauge
  weighted_sum <- 0
  n <- 0                           # neighbours with valid data

  for (id in names(neighbours)) {
    Pi <- month_value(neighbours[[id]], year, month)
    if (is_valid_nonneg(Pi)) {             # neighbour has data and no gap
      Ni <- annual_means[[id]]             # long-term annual mean of neighbour
      if (Ni > 0) {                        # avoid division by zero
        weighted_sum <- weighted_sum + (Nx / Ni) * Pi[1]
        n <- n + 1
      }
    }
  }

  if (n > 0) weighted_sum / n else -999  # -999: no usable neighbour
}
