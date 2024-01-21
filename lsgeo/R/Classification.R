#' Return_Classification
#'
#' @param Y adjacent matrix
#' @param clique_size a number
#' @param K number of cliques
#' @param m clique_size - 1
#' @param mode current not use, default set to "Not Data"
#' @param computeRate current not use, default set to "False"
#' @param finalCliques list of doubles
#' @param dHat matrix
#' @param E_nu default set to 1
#'
#' @return list of classification_max, p_E, p_S, p_H, kappaHat_S, kappaHat_H, a, b
#' @export
#' @importFrom igraph graph.adjacency
Return_Classification <- function(Y, clique_size, K, m, mode, computeRate, finalCliques, dHat, E_nu = 1) {
  G <- graph.adjacency(Y, mode="undirected")
  n <- ncol(Y)
  J <- diag(K) - 1/K * outer(rep(1, K), rep(1, K))

  fHat <- -0.5 * J %*% (dHat^2) %*% J
  lambda_F <- sort(eigen(fHat, only.values=TRUE)$values)[1]

  print(dHat)
  upper_tri_values <- dHat[upper.tri(dHat)]
  # Calculate b
  b <- (pi / max(upper_tri_values))^2
  # Calculate a
  a <- (3 * min(upper_tri_values))^(-2)

  kappaHat_S <- mean(sapply(c(0, 1, 2), function(j) Estimate_Curvature_Sphere(dHat, K, a, b, j)))
  lambda_CS <- sort(eigen(cos(sqrt(kappaHat_S) * dHat), symmetric = TRUE)$values)[1]


  # Calculate kappaHat_H
  kappaHat_H <- mean(sapply(c(-2, -3, -4), function(j) Estimate_Curvature_Hyperbolic(dHat, K, -a, -b, j)))
  eigenvalues_CH <- eigen(cosh(sqrt(-kappaHat_H) * dHat), symmetric = TRUE)$values
  sorted_eigenvalues_CH <- sort(eigenvalues_CH, decreasing = TRUE)
  lambda_CH <- sorted_eigenvalues_CH[2]  # Second largest eigenvalue

  # Set rates
  rate_E <- 1/3
  rate_S <- 1/3
  rate_H <- 1/3

  # transform finalCliques to matrix
  finalCliques <- do.call(rbind, finalCliques)

  # Call Bootstrap_Romano function for Euclidean
  result_E <- Bootstrap_Romano(Y, finalCliques, clique_size, K, E_nu, J, lambda_F, 0.05, m, "E", 0, rate = rate_E, zeta = 1)
  rejectEuclidean <- result_E[[1]]
  p_E <- result_E[[2]]

  # Call Bootstrap_Romano function for Spherical
  result_S <- Bootstrap_Romano(Y, finalCliques, clique_size, K, E_nu, J, lambda_CS, 0.05, m, "S", kappaHat_S, rate = rate_S, zeta = 1)
  rejectSpherical <- result_S[[1]]
  p_S <- result_S[[2]]

  # Call Bootstrap_Romano function for Hyperbolic
  result_H <- Bootstrap_Romano(Y, finalCliques, clique_size, K, E_nu, J, lambda_CH, 0.05, m, "H", kappaHat_H, rate = rate_H, zeta = 1)
  rejectHyperbolic <- result_H[[1]]
  p_H <- result_H[[2]]


  # Find the maximum p value
  pMax <- max(p_E, p_S, p_H)

  # Determine classification_max based on the maximum p value
  classification_max <- ifelse(p_E == pMax, 0,
                               ifelse(p_S == pMax, 1, 2))


  # Set W_Hat, geometry, and kappaHat based on classification_max
  if (classification_max == 0) {
    W_Hat <- fHat
    geometry <- "Euclidean"
    kappaHat <- 0
  } else if (classification_max == 1) {
    W_Hat <- cos(sqrt(kappaHat_S) * dHat)
    geometry <- "Spherical"
    kappaHat <- kappaHat_S
  } else if (classification_max == 2) {
    W_Hat <- cosh(sqrt(-kappaHat_H) * dHat)
    geometry <- "Hyperbolic"
    kappaHat <- kappaHat_H
  }

  # Return the results
  return(list(classification_max, p_E, p_S, p_H, kappaHat_S, kappaHat_H, a, b))
}
