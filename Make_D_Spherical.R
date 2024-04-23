## outcome: same as python (matrix, contain 0)

Make_Separated_D_Spherical_3D <- function(K, n, kappa, ub_theta, ub_phi, delta) {
  DGroups <- matrix(0, nrow = K, ncol = K)
  theta <- runif(K, min = 0, max = ub_theta)
  phi <- runif(K, min = 0, max = ub_phi)
  Z <- matrix(0, nrow = K, ncol = 3)
  
  for (i in 1:K) {
    Z[i, ] <- c(1/sqrt(kappa) * sin(theta[i]) * cos(phi[i]), 
                1/sqrt(kappa) * sin(theta[i]) * sin(phi[i]), 
                1/sqrt(kappa) * cos(theta[i]))
  }
  
  for (i in 1:K) {
    for (j in 1:K) {
      if (i != j) {
        DGroups[i, j] <- 1/sqrt(kappa) * acos(kappa * sum(Z[i, ] * Z[j, ]))
      }
    }
  }
  
  theta_Pos <- numeric(n)
  phi_Pos <- numeric(n)
  c <- rep(0, n)                         # group memberships
  Z <- matrix(0, nrow = n, ncol = 3)
  for (i in 1:n) {
    # group_index is same as in python (may start from 0)
    group_index <- as.integer(floor((i - 1) * K / n))
    theta_Pos[i] <- runif(1, min = theta[group_index+1] - delta, max = theta[group_index+1] + delta)
    phi_Pos[i] <- runif(1, min = phi[group_index+1] - delta, max = phi[group_index+1] + delta)
    Z[i, ] <- c(1/sqrt(kappa) * sin(theta_Pos[i]) * cos(phi_Pos[i]), 
                1/sqrt(kappa) * sin(theta_Pos[i]) * sin(phi_Pos[i]), 
                1/sqrt(kappa) * cos(theta_Pos[i]))
    c[i] <- group_index
  }
  
  D <- matrix(0, nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        D[i, j] <- 1/sqrt(kappa) * acos(sum(Z[i, ] * Z[j, ]) * kappa)
      }
    }
  }
  diag(D) <- 0
  return(list(DGroups, D, c))
}
