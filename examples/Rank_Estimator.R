## pay attention to the dim of clique and indices1 and indices2
## output: the python index, which means it starts from 0.


# Loading necessary libraries
library(Matrix) # For matrix operations, equivalent to numpy
library(stats) # For statistical functions, similar to scipy.stats
library(graphics) # For plotting, similar to matplotlib

RankEstimate <- function(FHat, K, B, Y, clique, clique_size, E_nu, J, geometry, kappaHat) {
  
  # Eigenvalue Decomposition
  eig <- eigen(FHat)
  eigvalFHat <- eig$values # Eigenvalues
  eigvecFHat <- eig$vectors # Eigenvectors
  
  # Sorting eigenvalues and eigenvectors in decreasing order
  idx <- order(eigvalFHat, decreasing = TRUE)
  eigvalFHat <- eigvalFHat[idx]
  eigvecFHat <- eigvecFHat[, idx]
  
  # Bootstrap
  BStar <- array(0, dim = c(K, K, B))
  for (indexB in 1:B) {
    PB <- matrix(0, nrow = K, ncol = K)
    indices <- sample(0:(clique_size-1), clique_size, replace = TRUE)
    for (i1 in 1:K) {
      for (i2 in 1:K) {
        if (i1 < i2) {
          indices1 <- as.integer(clique[i1, indices + 1])
          indices2 <- as.integer(clique[i2, indices + 1])
          PB[i1, i2] <- max(1/clique_size^2, sum(Y[indices1 + 1, indices2 + 1]))
        }
      }
    }
    PB <- (PB + t(PB)) / clique_size^2
    diag(PB) <- E_nu
    DB <- -log(PB / E_nu)
    diag(DB) <- 0
    if (geometry == 'E') {
      FB <- -1/2 * J %*% (DB^2) %*% J
      eig <- eigen(FB)
    } else if (geometry == 'S') {
      CB <- cos(DB * sqrt(kappaHat))
      eig <- eigen(CB)
    } else if (geometry == 'H') {
      CB <- cosh(sqrt(-kappaHat) * DB)
      eig <- eigen(CB %*% CB)
    }
    #
    eigval <- eig$values
    eigvec <- eig$vectors
    idx <- order(eigval, decreasing = TRUE)
    eigval <- eigval[idx]
    eigvec <- eigvec[, idx]
    BStar[,,indexB] <- eigvec
  }
  
  phi <- function(k) {
    return(eigvalFHat[k] / sum(eigvalFHat))
  }
  
  
  f0 <- function(k) {
    if (k == 1) {
      return(0)
    }
    if (k > 1) {
      temp <- numeric(B) # create a numeric vector of length B
      for (i in 1:B) {
        temp[i] <- 1 - abs(det(t(eigvecFHat[,1:k]) %*% BStar[,1:k,i]))
      }
      return(mean(temp))
    }
  }
  
  fn <- function(j) {
    return(f0(j) / sum(sapply(1:(K-1), f0)))
  }
  
  g <- function(k) {
    return(fn(k) + phi(k))
  }
  
  valueC <- sapply(1:(K-1), g)
  # this is for the same python index, which may start from 0
  return(which.min(valueC) - 1)
  # for the R index: return(which.min(valueC))
}


