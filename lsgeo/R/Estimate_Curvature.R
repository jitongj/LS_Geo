#' @ Estimate_Curvature_Sphere Function
#'
#' @param dHat a matrix
#' @param K number of cliques
#' @param a lower bound
#' @param b upper bound
#' @param index a index number
#'
#' @return minimum kappaVec
#' @export
#'
#' @importFrom stats quantile
Estimate_Curvature_Sphere <- function(dHat, K, a, b, index) {
  kappaVec <- seq(a, b, length.out = 10000)
  objFun <- numeric(length(kappaVec))

  for (i in seq_along(kappaVec)) {
    C <- cos(dHat * sqrt(kappaVec[i]))
    eigenvalues<- sort(eigen(C, symmetric = TRUE)$values, decreasing = FALSE)
    # if "index" is a R index: objFun[i] <- abs(eigenvalues[index])
    # if "index" is a python index:
    objFun[i] <- abs(eigenvalues[index+1])
  }

  index_min <- which.min(objFun)
  return(kappaVec[index_min])
}



#' @ Estimate_Curvature_Hyperbolic Function
#'
#' @param dHat a matrix
#' @param K number of cliques
#' @param a lower bound
#' @param b upper bound
#' @param index a index number
#'
#' @return minimum kappaVec
#' @export
#'
#' @importFrom stats quantile
Estimate_Curvature_Hyperbolic <- function(dHat, K, a, b, index) {
  kappaVec <- seq(a, b, length.out = 10000)
  objFun <- numeric(length(kappaVec))

  for (i in seq_along(kappaVec)) {
    C <- cosh(dHat * sqrt(-kappaVec[i]))
    eigenvalues <- eigen(C, symmetric = TRUE)$values
    n <- length(eigenvalues)
    index <- n + index + 1
    # if "index" is a R index: objFun[i] <- abs(eigenvalues[index])
    # if "index" is a python index:
    objFun[i] <- abs(eigenvalues[index+1])
  }

  index_min <- which.min(objFun)
  return(kappaVec[index_min])
}

