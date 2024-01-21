#' squeeze function
#'
#' @param arr array with one or more dimensions of size 1
#'
#' @return array with all single-dimensional (length 1) entries removed from its shape
#' @export
squeeze <- function(arr) {
  # Removing dimensions of length 1
  dims <- dim(arr) # Get dimensions of the array
  if (is.null(dims)) {
    # If 'arr' is a vector, return it as is
    return(arr)
  }

  new_dims <- dims[dims != 1] # Keep dimensions not equal to 1
  if (length(new_dims) == 0) {
    # If all dimensions are 1, return a single value
    return(arr[1])
  }

  dim(arr) <- new_dims # Set the new dimensions
  return(arr)
}



#' ComputePerformance
#'
#' @param P list of parameters
#'
#' @return list of reject, pvalue, kappaHat, D, DGroups,HatE_nu_vec
#' @export
#' @importFrom igraph cliques
#' @examples
#' \dontrun{
#' P <- list(
#'   n = 1200, # size of graph
#'   K_vec = c(5, 7, 10), # values of K (number of cliques)
#'   dim = 3, # dimension of latent space
#'   kappa = 0.75, # curvature of latent space
#'   numSim = 10, # number of sets of cliques to simulate
#'   beta = -0.01, # lower bound for distribution of fixed effects
#'   ell_vec = c(7), # values of ell (clique size)
#'   alpha = 0.05, # significance level of test
#'   zeta = 1, # scaling in front of distances (default = 1)
#'   scale = 2.5, # scale for generating points in hyperbolic space
#'   ub_theta = pi, # upper bound for theta in spherical LS
#'   ub_phi = 2 * pi, # upper bound for phi in spherical LS
#'   delta = 10^(-2), # spread of nodes around group centers in spherical LS
#'   rate = 1/3, # rate used when sub-sampling
#'   sigma = 0.8, # spread of points in Euclidean space
#'   geometry_true = 'E', # geometry of simulated graphs
#'   geometry_test = 'S', # geometry to test in latent space
#'   numSamples = 45, # nodes in one group to find cliques
#'   a = 1/3, # lower bound for searching kappa
#'   b = 1.5, # upper bound for searching kappa
#'   q = 2, # number of eigenvalues to minimize for kappa estimate
#'   m = c(6) # sub-sample rate, m = clique_size - 1
#' )
#' result <- ComputePerformance(P)
#' }
ComputePerformance <- function(P) {
  # 初始化变量
  wHat_array <- vector("list", length(P$K_vec))
  dHat_array <- vector("list", length(P$K_vec))
  for (i in seq_along(P$K_vec)) {
    wHat_array[[i]] <- vector("list", length(P$ell_vec))
    dHat_array[[i]] <- vector("list", length(P$ell_vec))
    for (j in seq_along(P$ell_vec)) {
      wHat_array[[i]][[j]] <- list()
      dHat_array[[i]][[j]] <- list()
    }
  }

  reject <- array(0, dim = c(length(P$K_vec), length(P$ell_vec), P$numSim))
  pvalue <- array(0, dim = c(length(P$K_vec), length(P$ell_vec), P$numSim))
  kappaHat <- array(0, dim = c(length(P$K_vec), length(P$ell_vec), P$numSim))
  HatE_nu_vec <- numeric(P$numSim)

  # 生成潜在空间位置
  # --------------------------------------------------------------------------------
  if (P$geometry_true == 'E') {
    kappa <- 0
    results <- Make_Euclidean_D(P$K_vec[length(P$K_vec)], P$n, P$sigma, P$dim)
    DGroups <- results$DGroups
    D <- results$D
  }
  if (P$geometry_true == "S") {
    results <- Make_Separated_D_Spherical_3D(P$K_vec[length(P$K_vec)], P$n, P$kappa, P$ub_theta, P$ub_phi, P$delta)
    DGroups <- results[[1]]   # 第一个元素
    D <- results[[2]]         # 第二个元素
  }
  if (P$geometry_true == "H") {
    results <- Generate_Hyperbolic(P$K_vec[length(P$K_vec)], P$dim, P$scale, P$kappa)
    muVec <- results[[1]]     # 第一个元素
    DGroups <- results[[2]]   # 第二个元素
    D <- Make_D_FromMu_Hyperbolic(P$n, P$dim, P$K_vec[length(P$K_vec)], muVec, P$kappa, P$delta)
  }
  # --------------------------------------------------------------------------------
  # 生成图
  nu <- runif(P$n, P$beta, 0) # 计算固定效应
  E_nu <- mean(exp(nu))^2
  edge_probabilities <- vector(mode = "numeric", length = choose(P$n, 2))
  iit <- 1
  for (i in 1:(P$n - 1)) {
    for (j in (i + 1):P$n) {
      edge_probabilities[iit] <- exp(nu[i] + nu[j] - P$zeta * D[i, j])
      iit <- iit + 1
    }
  }
  # --------------------------------------------------------------------------------
  dHat_array_changing_ell <- vector("list", length(P$ell_vec))
  finalCliques <- vector("list", length(P$ell_vec))
  # 对于每个元素，初始化为空列表
  for (i in 1:length(P$ell_vec)) {
    dHat_array_changing_ell[[i]] <- list()
    finalCliques[[i]] <- list()
  }

  # 假设 P 是一个列表，包含 numSim, ell_vec, K_vec, n, 和其他必要的参数

  for (indexSim in 1:P$numSim) {
    # the print output is R index, ie starts with 1
    cat("Simulation:", indexSim, "\n")
    t <- P$ell_vec[length(P$ell_vec)] - 1  # number of edges used in the "almost clique" to compute FE expectation.

    # 调用 GenerateY 函数，并将结果分配给相应的变量
    results <- GenerateY(P$K_vec[length(P$K_vec)], P$n, edge_probabilities, P$ell_vec[length(P$ell_vec)], E_nu, P$numSamples, t)
    Y <- results$Y
    pHat <- results$pHat
    dHat <- results$dHat
    cliques <- results$finalCliques
    HatE_nu_vec[indexSim] <- results$E_nu

    dHat <- dHat / P$zeta


    # 假设 P 是一个列表，包含 ell_vec, K_vec, n, 和其他必要的参数
    if ((length(P$ell_vec) - 1) > 0){
      for (indexTemp in 1:(length(P$ell_vec) - 1)) {
        # the print output is R index, ie starts with 1
        cat("this is index:", indexTemp, "\n")
        t <- P$ell_vec[indexTemp] - 1  # number of edges used in the "almost clique" to compute FE expectation.

        # 调用 findCliques_GivenY 函数，并将结果分配给相应的变量
        results <- findCliques_GivenY(Y, P$K_vec[length(P$K_vec)], P$n, P$ell_vec[indexTemp], E_nu, P$numSamples, t)
        dHat_array_changing_ell[[indexTemp]] <- results[[1]]
        finalCliques[[indexTemp]] <- results[[2]]
      }
    }


    dHat_array_changing_ell[[length(dHat_array_changing_ell)]] <- dHat
    finalCliques[[length(finalCliques)]] <- cliques

    for (indexEll in 1:length(P$ell_vec)) {
      temp_dHat <- squeeze(dHat_array_changing_ell[[indexEll]]) # 返回二维数组
      mStar <- P$m[[indexEll]]


      for (indexK in 1:length(P$K_vec)) {
        cat("indexEll:", indexEll, "\n")
        cat("indexK:", indexK, "\n")
        cat("this is d current", "\n")
        print(temp_dHat)

        K <- P$K_vec[indexK]
        J <- diag(K) - matrix(1/K, K, K)
        # indices_sampled is a R index, which = python index + 1
        indices_sampled <- sample(K, K, replace = FALSE)
        dHat_array[[indexK]][[indexEll]] <- temp_dHat[indices_sampled, indices_sampled]
        cat("this is J", "\n")
        print(J)
        print(dHat_array[[indexK]][[indexEll]])


        if (P$geometry_test == 'E') {
          wHat_array[[indexK]][[indexEll]] <- -1/2 * J %*% (dHat_array[[indexK]][[indexEll]]^2) %*% J
          kStar <- 0
          kappaHat[indexK, indexEll, indexSim] <- 0
        }

        if (P$geometry_test == "S") {
          kappaHat[indexK, indexEll, indexSim] <- mean(sapply(0:(P$q - 1), function(index) {
            Estimate_Curvature_Sphere(dHat_array[[indexK]][[indexEll]],
                                      P$K_vec[indexK], P$a, P$b, index)
          }))
          wHat_array[[indexK]][[indexEll]] <- cos(sqrt(kappaHat[indexK, indexEll, indexSim]) * dHat_array[[indexK]][[indexEll]])
          kStar <- 0
        }

        if (P$geometry_test == "H") {
          indices_eig <- sapply(0:(P$q - 1), function(temp) -2 - temp)
          kappaHat[indexK, indexEll, indexSim] <- mean(sapply(indices_eig, function(index) {
            Estimate_Curvature_Hyperbolic(dHat_array[[indexK]][[indexEll]],
                                          P$K_vec[indexK], -P$a, -P$b, index)
          }))
          wHat_array[[indexK]][[indexEll]] <- cosh(sqrt(-kappaHat[indexK, indexEll, indexSim]) * dHat_array[[indexK]][[indexEll]])
          kStar <- -2
        }


        eigWHat <- sort(eigen(wHat_array[[indexK]][[indexEll]])$values)

        #------------------------------
        # 需要将kStar转化成对应的R value
        if (kStar==0) {
          kStar = 1
        } else{
          kStar = length(eigWHat) - 1
        }
        #------------------------------

        temp <- Bootstrap_Romano(Y, finalCliques[[indexEll]], P$ell_vec[indexEll], P$K_vec[indexK],
                                 HatE_nu_vec[[indexSim]], J, eigWHat[kStar], P$alpha, mStar, P$geometry_test, kappaHat[indexK, indexEll, indexSim], P$rate, P$zeta)
        reject[indexK, indexEll, indexSim] <- temp[[1]]
        pvalue[indexK, indexEll, indexSim] <- temp[[2]]
      }
    }
  }


  # 使用自定义的 squeeze 函数移除单一维度
  LSreject <- squeeze(apply(reject, c(1, 2), mean)) # no use?

  # 构建返回的列表
  result <- list(
    reject = reject,
    pvalue = pvalue,
    kappaHat = kappaHat,
    D = D,
    DGroups = DGroups,
    HatE_nu_vec = HatE_nu_vec
  )

  # 返回结果
  return(result)
}

