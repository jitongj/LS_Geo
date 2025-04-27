## dim and n must be 3 according to the codes
## outcome same as python (matrix, contain 0)

library(ggplot2)
library(pracma) # For matrix operations and acosh function

Generate_Hyperbolic <- function(n, dim, scale, kappa) {

  ## x is not used inside!!
  genHyperbolic <- function(x) {
    temp <- numeric(dim)
    temp <- runif(dim, min = -scale, max = scale)
    # according to python codes, dim must = 3
    # temp[1] <- runif(1, min = -scale, max = scale) # Generate first random number
    # temp[2] <- runif(1, min = -scale, max = scale) # Generate second random number
    temp[dim] <- sqrt(1/kappa + sum(temp[1:(dim-1)]^2))
    return(temp)
  }

  # Generate points
  pos <- matrix(0, nrow = n, ncol = dim)
  for (i in 1:n) {
    pos[i, ] <- genHyperbolic(n)
  }

  temp <- pos[1, ]

  B <- function(x, y) {
    return(x[length(x)] * y[length(y)] - sum(x[1:(dim-1)] * y[1:(dim-1)]))
  }

  D <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i < j) {
        D[i, j] <- acosh(B(pos[i, ], pos[j, ]) * kappa) / sqrt(kappa)
      }
    }
  }

  D <- D + t(D)
  diag(D) <- 0

  return(list(pos, D))
}


BForm <- function(x, y, p) {
  return(x[length(x)] * y[length(y)] - sum(x[1:(p-1)] * y[1:(p-1)]))
}




Make_D_FromMu_Hyperbolic <- function(n, p, K, muVec, kappa, delta) {
  V <- matrix(0, n, p)
  c <- rep(0, n)

  for (i in 1:n) {
    c[i] <- floor(as.numeric(i - 1) * K / n)
    for (j in 1:(p - 1)) {
      V[i, j] <- runif(1, muVec[c[i] + 1, j] - delta, muVec[c[i] + 1, j] + delta)
    }
    V[i, p] <- sqrt(1 / kappa + sum(V[i, 1:(p - 1)]^2))
  }

  D <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        D[i, j] <- acosh(max(1, BForm(V[i, ], V[j, ], p) * kappa)) / sqrt(kappa)
      }
    }
  }

  return(D, c)
}
