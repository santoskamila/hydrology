# =============================================================================
# Helpers for the validation script (02_validate_methods.R)
# =============================================================================

# Iterative 1-component PCA (L = monthly scores, C = gauge weights)
fill_iterative_pca_matrix <- function(P_gappy) {
  n_months <- nrow(P_gappy); n_st <- ncol(P_gappy)
  max_iterations <- 500; tolerance <- 1e-6
  L <- rowMeans(P_gappy, na.rm = TRUE); L[is.na(L)] <- 0
  C <- numeric(n_st)
  for (it in seq_len(max_iterations)) {
    L_old <- L; C <- numeric(n_st)
    for (j in seq_len(n_st)) {
      v <- !is.na(P_gappy[, j])
      if (any(v)) C[j] <- sum(L[v] * P_gappy[v, j]) / sum(L[v]^2)
    }
    for (i in seq_len(n_months)) {
      v <- !is.na(P_gappy[i, ])
      if (any(v)) L[i] <- sum(P_gappy[i, v] * C[v]) / sum(C[v]^2)
    }
    if (sqrt(sum((L - L_old)^2)) < tolerance) break
  }
  P_filled <- P_gappy
  idx <- which(is.na(P_gappy), arr.ind = TRUE)
  P_filled[idx] <- L[idx[, 1]] * C[idx[, 2]]
  P_filled
}

# Complete monthly calendar (Jan of first year to Dec of last year)
full_date_grid <- function(data, st_names) {
  y0 <- min(sapply(st_names, function(s) min(data[[s]]$Year)))
  y1 <- max(sapply(st_names, function(s) max(data[[s]]$Year)))
  expand.grid(Month = 1:12, Year = y0:y1)[, c("Year", "Month")]
}

# List of monthly tables -> matrix [months x gauges]
to_matrix <- function(data, st_names) {
  dates <- full_date_grid(data, st_names)
  key <- dates$Year * 100 + dates$Month
  X <- matrix(NA_real_, nrow(dates), length(st_names), dimnames = list(NULL, st_names))
  for (i in seq_along(st_names)) {
    d <- data[[st_names[i]]]
    loc <- match(d$Year * 100 + d$Month, key)
    X[loc[!is.na(loc)], i] <- d$MonthlyTotal[!is.na(loc)]
  }
  X
}

# Matrix [months x gauges] -> list of monthly tables (using a template)
to_list <- function(X, template, st_names) {
  dates <- full_date_grid(template, st_names)
  key <- dates$Year * 100 + dates$Month
  for (i in seq_along(st_names)) {
    d <- template[[st_names[i]]]
    loc <- match(d$Year * 100 + d$Month, key)
    d$MonthlyTotal[!is.na(loc)] <- X[loc[!is.na(loc)], i]
    template[[st_names[i]]] <- d
  }
  template
}

# Euclidean distance between gauges over their common months
distance_matrix <- function(data, st_names) {
  n <- length(st_names)
  D <- matrix(0, n, n, dimnames = list(st_names, st_names))
  for (i in seq_len(n - 1)) for (j in (i + 1):n) {
    dc <- join_pair(data[[st_names[i]]], data[[st_names[j]]])
    dc <- dc[!is.na(dc$MonthlyTotal_x) & !is.na(dc$MonthlyTotal_y), ]
    d <- if (nrow(dc) > 0) sqrt(sum((dc$MonthlyTotal_x - dc$MonthlyTotal_y)^2)) else Inf
    D[i, j] <- d; D[j, i] <- d
  }
  D
}

# Approximate equivalent of MATLAB fillmissing(X, 'knn') with k = 1:
# each missing value takes the value of the nearest row (month) that has it,
# using a NA-ignoring Euclidean distance scaled by the number of common columns.
fill_missing_knn <- function(X) {
  Xf <- X; p <- ncol(X)
  for (i in which(rowSums(is.na(X)) > 0)) {
    obs_i <- !is.na(X[i, ])
    if (!any(obs_i)) next
    d <- apply(X, 1, function(r) {
      com <- obs_i & !is.na(r)
      if (!any(com)) return(Inf)
      sqrt(sum((r[com] - X[i, com])^2) * p / sum(com))
    })
    d[i] <- Inf
    for (j in which(!obs_i)) {
      cand <- which(!is.na(X[, j]) & is.finite(d))
      if (length(cand) > 0) Xf[i, j] <- X[cand[which.min(d[cand])], j]
    }
  }
  Xf
}

# Performance metrics per gauge and method
compute_metrics <- function(methods, st_names, gap_locations, observed, predicted) {
  res <- list()
  for (s in seq_along(st_names)) {
    tab <- data.frame(RMSE = 0, MAE = 0, R2 = 0, NSE = 0, Spearman = 0, Willmott_d = 0)[rep(1, length(methods)), ]
    rownames(tab) <- methods
    idx <- which(gap_locations[, 1] == s)
    if (length(idx) > 0) {
      obs <- observed[idx]
      for (m in seq_along(methods)) {
        pred <- predicted[idx, m]
        v <- !is.na(pred) & pred >= 0
        o <- obs[v]; pv <- pred[v]
        if (length(o) > 1) {
          tab$RMSE[m] <- sqrt(mean((pv - o)^2))
          tab$MAE[m] <- mean(abs(pv - o))
          tab$R2[m] <- suppressWarnings(cor(pv, o))^2
          tab$NSE[m] <- 1 - sum((o - pv)^2) / sum((o - mean(o))^2)
          tab$Spearman[m] <- suppressWarnings(cor(pv, o, method = "spearman"))
          tab$Willmott_d[m] <- 1 - sum((pv - o)^2) / sum((abs(pv - mean(o)) + abs(o - mean(o)))^2)
        }
      }
    }
    res[[st_names[s]]] <- tab
  }
  res
}
