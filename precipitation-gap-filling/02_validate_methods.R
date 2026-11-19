# =============================================================================
# STEP 2: CROSS-VALIDATION OF GAP-FILLING METHODS
# 1. Builds a "complete" reference dataset (k-NN filling of the real gaps).
# 2. Removes a random share of the values (artificial gaps).
# 3. Refills them with each method and computes performance metrics.
#
# Requires: 01_prepare_data.R (results/prepared_data.rds)
# Output:   results/validation_results.rds, figures/cluster_dendrogram.png
# Packages: base R only (parallel is optional)
# =============================================================================
rm(list = ls())
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

# ---- 1. Settings ------------------------------------------------------------------
percent_artificial_gaps <- 5
max_neighbours <- 10
use_parallel <- FALSE   # TRUE = use several CPU cores (package parallel)
seed <- 1
dir.create("figures", showWarnings = FALSE)
cat(sprintf("Starting validation with %.0f%% artificial gaps.\n", percent_artificial_gaps))

# ---- 2. Load data -------------------------------------------------------------------
cat("Loading data...\n")
dp <- readRDS("results/prepared_data.rds")
monthly <- dp$monthly; annual_means <- dp$annual_means; monthly_means <- dp$monthly_means
station_names <- names(monthly)
n_st <- length(station_names)
for (st in station_names) {                     # gap months -> NA
  x <- monthly[[st]]$MonthlyTotal
  x[x < 0] <- NA
  monthly[[st]]$MonthlyTotal <- x
}

# ---- 3. Cluster analysis -------------------------------------------------------------
cat("Computing similarity and cluster analysis...\n")
dist_mat <- distance_matrix(monthly, station_names)
d_clust <- dist_mat
if (any(is.infinite(d_clust))) {
  warning("Some gauge pairs share no months; Inf distances replaced by 2x the largest finite distance in the dendrogram.")
  d_clust[is.infinite(d_clust)] <- 2 * max(d_clust[is.finite(d_clust)])
}
if (n_st > 2) {
  tree <- hclust(as.dist(d_clust), method = "ward.D2")   # = MATLAB linkage(..., 'ward')
  png("figures/cluster_dendrogram.png", width = 1000, height = 700)
  plot(tree, labels = station_names, xlab = "", sub = "",
       main = "Hierarchical Cluster Analysis of Rain Gauges",
       ylab = "Distance (similarity)")
  dev.off()
  cat("Dendrogram saved to figures/.\n\n")
}

# ---- 4. "Complete" reference dataset (ground truth) --------------------------------
cat("Building the complete reference dataset (k-NN filling)...\n")
X_complete <- fill_missing_knn(to_matrix(monthly, station_names))
complete_data <- to_list(X_complete, monthly, station_names)

# ---- 5. Artificial gaps -------------------------------------------------------------
cat("Introducing random artificial gaps...\n")
gappy <- complete_data
all_points <- do.call(rbind, lapply(seq_len(n_st), function(i)
  cbind(i, seq_len(nrow(gappy[[station_names[i]]])))))
set.seed(seed)
n_total <- nrow(all_points)
n_gaps <- round(percent_artificial_gaps / 100 * n_total)
gap_locations <- all_points[sample(n_total, n_gaps), , drop = FALSE]   # [gauge index, row]
observed <- numeric(n_gaps)
for (i in seq_len(n_gaps)) {
  st <- station_names[gap_locations[i, 1]]; r <- gap_locations[i, 2]
  observed[i] <- gappy[[st]]$MonthlyTotal[r]
  gappy[[st]]$MonthlyTotal[r] <- NA
}

# ---- 6. Fill the gaps ---------------------------------------------------------------
cat("Filling gaps...\n")
# 6a. Whole-matrix method: iterative PCA
P_pca <- fill_iterative_pca_matrix(to_matrix(gappy, station_names))
dates <- full_date_grid(gappy, station_names)
pred_pca <- sapply(seq_len(n_gaps), function(i) {
  d <- gappy[[station_names[gap_locations[i, 1]]]]
  r <- gap_locations[i, 2]
  P_pca[which(dates$Year == d$Year[r] & dates$Month == d$Month[r]), gap_locations[i, 1]]
})
cat("  > Iterative PCA done.\n")

# 6b. Gap-by-gap methods
gap_methods <- c("RW", "SLR_Var", "RV", "PCR_NIPALS", "KNN_Analog", "KNN_Simple", "MLR_Fixed", "Arithmetic_Mean")
#   RW = Regional Weighting, RV = Regional Vector

fill_one_gap <- function(i) {
  i_target <- gap_locations[i, 1]; row <- gap_locations[i, 2]
  st <- station_names[i_target]
  gy <- gappy[[st]]$Year[row]
  gm <- gappy[[st]]$Month[row]

  # Neighbours with data in the same month (matched by Year/Month),
  # ranked by distance, up to max_neighbours
  avail <- c()
  for (s in seq_along(station_names)) {
    if (s == i_target) next
    if (is_valid_value(month_value(gappy[[station_names[s]]], gy, gm))) avail <- c(avail, s)
  }
  nb <- list()
  if (length(avail) > 0) {
    ord <- order(dist_mat[avail, i_target])
    best <- avail[ord[seq_len(min(max_neighbours, length(ord)))]]
    for (s in best) nb[[station_names[s]]] <- gappy[[station_names[s]]]
  }

  safe <- function(expr) tryCatch(expr, error = function(e) NA_real_)
  c(safe(fill_regional_weighting(st, nb, gy, gm, annual_means)),
    safe(fill_slr_variable(st, gappy, nb, gy, gm)),
    safe(fill_regional_vector(st, gappy, nb, gy, gm, monthly_means)),
    safe(fill_pcr_nipals(st, gappy, nb, gy, gm)),
    safe(fill_knn_analog(st, gappy, nb, gy, gm, 5)),
    safe(fill_knn_simple(st, gappy, 5, dist_mat, station_names, gy, gm)),
    safe(fill_mlr_fixed(st, gappy, nb, gy, gm, 5)),
    safe(fill_arithmetic_mean(st, gappy, nb, gy, gm)))
  # safe(fill_multiple_imputation(st, gappy, nb, gy, gm))
}

if (use_parallel) {
  library(parallel)
  cl <- makeCluster(max(1, detectCores() - 1))
  clusterExport(cl, ls())
  res <- parLapply(cl, seq_len(n_gaps), fill_one_gap)
  stopCluster(cl)
} else {
  res <- lapply(seq_len(n_gaps), function(i) {
    if (i %% 20 == 0) cat(sprintf("    gap %d of %d\n", i, n_gaps))
    fill_one_gap(i)
  })
}
cat("Gap-by-gap methods done.\n")

# ---- 7. Metrics, display and save -----------------------------------------------
methods <- c(gap_methods, "Iterative_PCA")
predicted <- cbind(do.call(rbind, res), pred_pca)
colnames(predicted) <- methods
results <- compute_metrics(methods, station_names, gap_locations, observed, predicted)

cat("\n=======================================================\n")
cat("              VALIDATION RESULTS\n")
cat("=======================================================\n\n")
for (st in station_names) {
  cat(sprintf(">>> Method performance for gauge: %s\n", st))
  print(round(results[[st]], 4))
  cat("\n")
}
saveRDS(list(results = results, gap_locations = gap_locations, observed = observed,
             predicted = predicted, methods = methods, station_names = station_names,
             X_complete = X_complete),
        "results/validation_results.rds")
cat("Validation results saved to results/validation_results.rds\n")
