library(MASS)

Estimate_Rate <- function(Y, K, clique_size, E_nu, cliques, lambdaHat, geometry, kappa = 0) {
  J <- diag(K) - 1/K * matrix(1, nrow = K, ncol = K)
  b <- c(3^2, 5^2, 7^2)
  s <- seq(0.01, 0.49, length.out = 10)
  t <- seq(0.51, 0.99, length.out = 10)
  
  Return_Inv_CDF <- function(Y, b, s_i, t_i) {
    s <- 1000  # number of bootstrap samples
    lambdaStar <- numeric(s)
    for (i in 1:s) {
      pStar <- matrix(0, nrow = K, ncol = K)
      indices <- sample(1:(clique_size^2), size = b, replace = TRUE)
      for (k in 1:K) {
        for (kPrime in 1:K) {
          if (k < kPrime) {
            temp <- as.vector(Y[cliques[k, ], cliques[kPrime, ]])
            pStar[k, kPrime] <- max(1 / clique_size^2, mean(temp[indices]))
          }
        }
      }
      pStar <- pStar + t(pStar)
      diag(pStar) <- E_nu
      dStar <- -log(pStar / E_nu)
      if (geometry == 'Euclidean') {
        wStar <- -1/2 * J %*% (dStar^2) %*% J
        lambdaStar[i] <- eigen(wStar)$values[1] - lambdaHat
      }
      if (geometry == 'Spherical') {
        wStar <- cos(sqrt(kappa) * dStar)
        lambdaStar[i] <- eigen(wStar)$values[1] - lambdaHat
      }
      if (geometry == 'Hyperbolic') {
        wStar <- cosh(sqrt(-kappa) * dStar)
        lambdaStar[i] <- tail(eigen(wStar)$values, 1) - lambdaHat
      }
    }
    c_alpha_s <- quantile(lambdaStar, probs = s_i)
    c_alpha_t <- quantile(lambdaStar, probs = t_i)
    return(log(c_alpha_t - c_alpha_s))
  }
  
  Y_Matrix <- matrix(0, nrow = length(b), ncol = length(s))
  for (index1 in 1:length(b)) {
    for (index2 in 1:length(s)) {
      Y_Matrix[index1, index2] <- Return_Inv_CDF(Y, b[index1], s[index2], t[index2])
    }
  }
  Y_Dot <- rowMeans(Y_Matrix)
  Y_bar <- mean(Y_Matrix)
  log_bar <- mean(log(b))
  num <- sum((Y_Dot - Y_bar) * (log(b) - log_bar))
  denom <- sum((log(b) - log_bar)^2)
  rate <- -num / denom  # compute rate
  return(rate)
}
