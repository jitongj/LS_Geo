## same outcome with python (distant matrix)

library(MASS) # for mvrnorm (multivariate normal distribution)

Make_Euclidean_D <- function(K, n, sigma, p) {
  # Initialize variables
  muVec <- matrix(0, nrow = K, ncol = p) # group centers
  V <- matrix(0, nrow = n, ncol = p)     # node locations
  c <- rep(0, n)                         # group memberships
  
  # Generate group centers
  for (i in 1:K) {
    muVec[i, ] <- mvrnorm(1, rep(0, p), sigma * diag(p))
  }
  
  # Compute distance matrix between groups
  DGroups <- as.matrix(dist(muVec))
  
  # Assign nodes to groups and generate locations
  for (i in 1:n) {
    c[i] <- floor(as.numeric(i - 1) * K / n) # same value with python
    V[i, ] <- mvrnorm(1, muVec[c[i] + 1, ], sigma / K * diag(p))
  }
  
  # Compute distance matrix between nodes
  D <- as.matrix(dist(V))
  
  # Return results
  return(list(DGroups = DGroups, D = D, c=c))
}