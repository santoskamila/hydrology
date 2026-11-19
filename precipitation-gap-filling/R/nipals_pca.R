nipals_pca <- function(X, num_components, tol = 1e-6, max_iter = 500) {
  # NIPALS algorithm: principal components of a matrix with missing data (NA).
  #
  # INPUTS:
  #   X              - data matrix (m observations x n variables), gaps = NA
  #   num_components - number of components to extract
  #   tol            - convergence tolerance (e.g. 1e-6)
  #   max_iter       - maximum iterations per component (e.g. 500)
  #
  # OUTPUT (list):
  #   T               - scores matrix
  #   P               - loadings matrix
  #   X_reconstructed - original X with the NAs filled by the model

  X <- as.matrix(X)
  m <- nrow(X); n <- ncol(X)
  X_current <- X

  T <- matrix(0, m, num_components)
  P <- matrix(0, n, num_components)

  for (h in seq_len(num_components)) {
    # Initialisation: first column with data is the initial score vector
    first_col <- which(colSums(!is.na(X_current)) > 0)[1]
    t_h <- X_current[, first_col]
    t_h[is.na(t_h)] <- 0

    for (iter in seq_len(max_iter)) {
      # Step 1: loadings (p) from scores (t), ignoring NAs
      p_h <- numeric(n)
      for (j in seq_len(n)) {
        ok <- !is.na(X_current[, j])
        p_h[j] <- sum(X_current[ok, j] * t_h[ok])
      }
      p_h <- p_h / sqrt(sum(p_h^2))   # unit length

      # Step 2: new scores (t) from loadings (p)
      t_h_old <- t_h
      t_h <- numeric(m)
      for (i in seq_len(m)) {
        ok <- !is.na(X_current[i, ])
        t_h[i] <- sum(X_current[i, ok] * p_h[ok])
      }

      # Step 3: convergence check
      if (sum((t_h - t_h_old)^2) < tol) break
    }

    T[, h] <- t_h
    P[, h] <- p_h

    # Step 4: deflate X
    X_current <- X_current - t_h %*% t(p_h)
  }

  # Reconstruct X and fill only the original NAs
  X_rec <- T %*% t(P)
  X_filled <- X
  X_filled[is.na(X)] <- X_rec[is.na(X)]

  list(T = T, P = P, X_reconstructed = X_filled)
}
