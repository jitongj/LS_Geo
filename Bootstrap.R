

library(MASS)  # For statistical methods
library(ggplot2)  # For plotting, similar to matplotlib in Python
library(matrixStats)
library(Matrix)


Bootstrap_Romano <- function(Y, cliques, clique_size, K, E_nu, J, lambdaHat, alpha, m, geometry, kappa, rate, zeta) {
  s <- 1000  # Number of bootstrap simulations
  temp <- numeric(s)  # Initialize a numeric vector of length s with zeroes
  lambdaStar <- numeric(s)  # Initialize a numeric vector of length s with zeroes
  
  cat("these are the cliques\n")
  print(cliques)
  
  cat("This is mstar\n")
  print(m)
  
  print(Y)
  print(clique_size)
  
  # Rest of the function logic goes here...
  for (i in 1:s) {
    pStar <- matrix(0, nrow = K, ncol = K)
    indices_1 <- sample(1:(clique_size - 1), size = m, replace = TRUE)
    indices_2 <- sample(1:(clique_size - 1), size = m, replace = TRUE)
    
    for (k in 1:K) {
      for (kPrime in 1:K) {
        if (k < kPrime) {

          # assumption：cliques not include 0) : Y_sub <- Y[cliques[k, ], cliques[kPrime, ]]
          # assumption：cliques include 0
          # if cliques is a matrix 
          Y_sub <- Y[cliques[k,] + 1, cliques[kPrime,]+1]
          # if cliques is a list
          #Y_sub <- Y[cliques[[k]] + 1, cliques[[kPrime]]+1]
          pStar[k, kPrime] <- max(1/clique_size^2, mean(Y_sub[indices_1 + 1, indices_2 + 1]))
        }
      }
    }
    pStar <- pStar + t(pStar)
    if (E_nu == 0) {
      E_nu <- 1
    }
    # Update Diagonal of pStar
    diag(pStar) <- E_nu
    
    # Calculate dStar
    dStar <- -log(pStar / E_nu)
    dStar <- dStar / zeta
    
    if (geometry == "E") {
      fStar <- -0.5 * J %*% (dStar^2) %*% J
      lambdaStar[i] <- min(eigen(fStar, symmetric = TRUE)$values) - lambdaHat
    }
    
    if (geometry == "S") {
      cStar <- cos(dStar * sqrt(kappa))
      lambdaStar[i] <- min(eigen(cStar, symmetric = TRUE)$values) - lambdaHat
    }
    
    if (geometry == "H") {
      cStar <- cosh(sqrt(-kappa) * dStar)
      lambdaStar[i]  <- sort(eigen(cStar, symmetric = TRUE)$values, decreasing = TRUE)[2] - lambdaHat
    }
    
  }
  

  if (geometry == "E") {
    ## type from 1-9 will generate same answer, except 7!!
    c_alpha <- quantile(m^(2*rate) * lambdaStar, probs = alpha, type = 8)
    pvalue <- sum(m^(2*rate) * lambdaStar < clique_size^(2*rate) * lambdaHat) / s
    return(list(c_alpha / (clique_size^(2*rate)) > lambdaHat, pvalue))
  }
  
  if (geometry == "S") {
    c_alpha <- quantile(m^(2*rate) * lambdaStar, probs = alpha, type = 8)
    pvalue <- sum(m^(2*rate) * lambdaStar < clique_size^(2*rate) * lambdaHat) / s
    return(list(c_alpha / (clique_size^(2*rate)) > lambdaHat, pvalue))
  }
  
  if (geometry == "H") {
    c_alpha <- quantile(m^(2*rate) * lambdaStar, probs = 1 - alpha,  type = 8)
    pvalue <- sum(m^(2*rate) * lambdaStar > clique_size^(2*rate) * lambdaHat) / s
    return(list(lambdaHat > c_alpha / (clique_size^(2*rate)), pvalue))
  }
  
}
