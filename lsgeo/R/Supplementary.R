#' sort_cliques function
#'
#' @param cliques_adjusted a list of lists
#'
#' @return reordered cliques
#' @export
#'
sort_cliques <- function(cliques_adjusted) {
  lengths_of_sublists <- lapply(cliques_adjusted, length)
  flaten = unlist(lengths_of_sublists)
  uniqe_size = unique(flaten)
  uniqe_size = sort(uniqe_size)
  combined_list_of_lists <- list()
  for (u in uniqe_size) {
    indices<- which(sapply(cliques_adjusted, length) == u)
    sublist = cliques_adjusted[indices]
    matrix_of_lists <- do.call(rbind, sublist)
    # Ordering the matrix and then applying that order to the list of lists
    ordered_indices <- do.call(order, as.data.frame(matrix_of_lists))
    ordered_list_of_lists <- sublist[ordered_indices]

    combined_list_of_lists <- c(combined_list_of_lists, ordered_list_of_lists)
  }

  return(combined_list_of_lists)
}


#' GenerateGraph
#'
#' @param n siz of graph
#' @param probs probability to create an edge
#'
#' @return matrix Y
#' @export
#'
GenerateGraph <- function(n, probs) {
  Y <- matrix(0, n, n)

  # Get the indices of the upper triangle of the matrix, excluding the diagonal
  upper_tri_indices <- which(upper.tri(Y, diag = FALSE), arr.ind = TRUE)

  # 调整索引为行优先
  upper_tri_indices_row_major <- upper_tri_indices[order(upper_tri_indices[, 1], upper_tri_indices[, 2]), ]


  # Generate random numbers and compare with probability to create an edge
  Y[upper_tri_indices_row_major] <- runif(n = floor(n*(n-1)/2)) < probs

  return(Y + t(Y))
}


#' EstimatePhat
#'
#' @param k a number
#' @param kPrime a number
#' @param K number of cliques
#' @param n clique size times number of cliques
#' @param c clique
#' @param Y adjacency matrix for subgraph induced by cliques
#'
#' @return list of two number
#' @export
EstimatePhat <- function(k, kPrime, K, n, c, Y) {
  G <- vector("list", K)
  term <- 0

  # Initialize each list element as empty
  for (i in 1:K) {
    G[[i]] <- integer(0)
  }

  # Build the groups and calculate term
  for (i in 1:n) {
    G[[c[i] + 1]] <- c(G[[c[i] + 1]], i - 1)
    if (c[i] == k) {
      term <- term + sum(Y[i, ] * (c == kPrime))
    }
  }

  # Return values based on condition
  if (k != kPrime) {
    return(c(term, length(G[[k + 1]]) * length(G[[kPrime + 1]])))
  } else {
    return(c(term / 2, n * (n - 1) / 2))
  }
}



Compute_ReducedY_timed <- function(Y, finalCliques, t, clique_size, K_vec, K_current, n, max_ell) {
  time_elapsed <- 0
  while (time_elapsed < 100) {
    time1 <- Sys.time()
    time1 <- as.numeric(as.POSIXct(time1, origin="1970-01-01"))
    # python中的结果是微秒级别，下面代码可以达到同样级别
    # current_time <- Sys.time()
    # formatted_time <- format(current_time, "%Y-%m-%d %H:%M:%OS6")
    # posix_time <- as.POSIXct(formatted_time, format="%Y-%m-%d %H:%M:%OS", tz="UTC")
    # time1 <- as.numeric(format(posix_time, "%s.%OS6"))


    pHat <- matrix(0, nrow = K_current, ncol = K_current)
    current_cliques <- matrix(0, nrow = K_current, ncol = clique_size)

    # output of current_clique is same as python
    for (index in 1:K_current) {
      random_indices <- sample(max_ell, clique_size, replace = FALSE)
      current_cliques[index, ] <- finalCliques[index, random_indices]
    }

    tempIndices <- as.integer(as.vector(t(current_cliques)))

    # be careful of the index, outcome of A_c is same as in python
    A_c <- Y[tempIndices+1, tempIndices+1]
    c_clique <- rep(0, K_current * clique_size)
    for (i in 1:(K_current*clique_size)) {
      c_clique[i] <- floor(as.numeric(i-1) / clique_size)
    }

    pHatVec <- c()
    for (i in 0:(K_current-1)) {
      for (j in 0:(K_current-1)){
        if (i < j) {
          result <- EstimatePhat(i, j, K_current, clique_size * K_current, c_clique, A_c)
          pHatVec <- c(pHatVec, result[1] / result[2])
        }
      }
    }


    # Assign values in row-major order (the default assigning in R is col-major order)
    row_indices <- integer(0)
    col_indices <- integer(0)
    for (i in 1:(K_current-1)) {
      for (j in (i+1):K_current) {
        row_indices <- c(row_indices, i)
        col_indices <- c(col_indices, j)
      }
    }
    pHat[cbind(row_indices, col_indices)] <- pHatVec


    pHat <- pHat + t(pHat)

    # Computation of E_nu
    E_nu <- mean(sapply(t, function(j) approximateEnu(n, K_current, j, Y, current_cliques, clique_size)))

    # Updating Diagonal of pHat
    if (E_nu == 0) {
      E_nu <- 1
    }
    diag(pHat) <- E_nu

    # Print the Sum of Zeros in pHat
    print(sum(pHat == 0))

    if (all(pHat != 0)) {
      cat("found this great thing\n")
      dHat <- -log(pHat / E_nu)
      return(list(1, dHat))
    }

    time2 <- Sys.time()
    time2 <- as.numeric(as.POSIXct(time2, origin="1970-01-01"))
    time_elapsed <- time_elapsed + as.numeric(difftime(time2, time1, units = "secs"))
    cat(sprintf("this is time elapsed: %f", time_elapsed), "\n")

  }
  return(list(0, pHat))
}


Compute_ReducedY <- function(Y, finalCliques, t, clique_size, K_vec, K_current, n, max_ell) {
  while(TRUE) {
    pHat <- matrix(0, nrow = K_current, ncol = K_current)
    current_cliques <- matrix(0, nrow = K_current, ncol = clique_size)

    for (index in 1:K_current) {
      random_indices <- sample(max_ell, clique_size, replace = FALSE)
      current_cliques[index, ] <- finalCliques[index, random_indices]
    }

    tempIndices <- as.integer(c(as.vector(t(current_cliques))))
    A_c <- Y[tempIndices+1, tempIndices+1]
    c_clique <- rep(0, times = K_current * clique_size)

    for (i in 1:(K_current * clique_size)) {
      c_clique[i] <- floor((i - 1) / clique_size)
    }

    pHatVec <- c()
    for (i in 0:(K_current-1)) {
      for (j in 0:(K_current-1)){
        if (i < j) {
          result <- EstimatePhat(i, j, K_current, clique_size * K_current, c_clique, A_c)
          pHatVec <- c(pHatVec, result[1] / result[2])
        }
      }
    }


    # Assign values in row-major order (the default assigning in R is col-major order)
    row_indices <- integer(0)
    col_indices <- integer(0)
    for (i in 1:(K_current-1)) {
      for (j in (i+1):K_current) {
        row_indices <- c(row_indices, i)
        col_indices <- c(col_indices, j)
      }
    }
    pHat[cbind(row_indices, col_indices)] <- pHatVec
    pHat <- pHat + t(pHat)


    E_nu <- mean(sapply(t, function(j) approximateEnu(n, K_current, j, Y, current_cliques, clique_size)))

    if (E_nu == 0) {
      E_nu <- 1
    }

    diag(pHat) <- E_nu

    if (all(pHat != 0)) {
      dHat <- -log(pHat / E_nu)
      return(list(dHat))
    }
  }
}



#' findCliques_GivenY
#'
#' @param Y adjacency matrix
#' @param K number of cliques
#' @param n size of graph
#' @param clique_size clique size
#' @param E_nu a number
#' @param numTaken determines how many nodes to take in one group to find cliques
#' @param t number of edges used in the "almost clique" to compute FE expectation.
#'
#' @return list of matrix dHat and matrix finalCliques
#' @export
#' @importFrom igraph graph_from_adjacency_matrix
findCliques_GivenY <- function(Y, K, n, clique_size, E_nu, numTaken, t) {
  nodesPerClique <- as.integer(n / K)
  tempCliques <- vector("list", K)
  for (j in 1:K) {
    tempCliques[[j]] <- list()
  }
  finalCliques <- matrix(0, nrow = K, ncol = clique_size)
  cat("trying this\n")

  while (TRUE) {
    pHat <- matrix(0, nrow = K, ncol = K)
    found_empty <- FALSE


    for (indexClique in 1:K) {
      # tempIndices is same as in python
      tempIndices <- (nodesPerClique * (indexClique - 1)) + (0:(numTaken-1))
      Ytemp <- Y[tempIndices+1, tempIndices+1]

      # Convert matrix to graph and find cliques，（调整cliques顺序后和python一致，且内容一致，即含0）
      g <- graph_from_adjacency_matrix(as.matrix(Ytemp), mode = "undirected", diag = FALSE)
      cliques_with_ids <- cliques(g)
      cliques <- sort_cliques(cliques_with_ids)
      cliques <- lapply(cliques, function(clique) {
        sapply(clique, function(node) {
          node - 1
        })
      })

      # Filter cliques of the desired size
      tempCliques[[indexClique]] <- Filter(function(k) {length(k) == clique_size}, cliques)

      # Print the number of cliques found
      print(length(tempCliques[[indexClique]]))

      # Check for empty cliques and potentially break the loop
      if (length(tempCliques[[indexClique]]) == 0) {
        found_empty <- TRUE
        break
      }
    }

    if(found_empty) {
      next
    }

    # 循环遍历
    for (indexClique in 1:K) {
      #index_random is a R index, which is python index +1
      index_random <- sample(length(tempCliques[[indexClique]]), 1)
      selected_clique <- tempCliques[[indexClique]][[index_random]]
      finalCliques[indexClique, ] <- sapply(selected_clique, function(j) j + nodesPerClique*(indexClique - 1))
    }

    # same as in python
    tempIndices <- as.integer(c(t(finalCliques)))

    # 创建子图的邻接矩阵
    A_c <- Y[tempIndices + 1, tempIndices + 1]

    # 初始化 c_clique 向量
    c_clique <- rep(0, K * clique_size)

    # 填充 c_clique 向量
    for (i in 1:(K * clique_size)) {
      c_clique[i] <- floor((i - 1) / clique_size)
    }

    # 计算 pHatVec
    pHatVec <- c()  # 创建一个空的向量

    # 使用两个嵌套的 for 循环计算 pHatVec
    for (i in 0:(K-2)) {
      for (j in (i+1):(K-1)) {
        result <- EstimatePhat(i, j, K, clique_size * K, c_clique, A_c)[1] /
          EstimatePhat(i, j, K, clique_size * K, c_clique, A_c)[2]
        # Append the result to pHatVec
        pHatVec <- c(pHatVec, result)
      }
    }

    # 更新 pHat 矩阵
    # Assign values in row-major order (the default assigning in R is col-major order)
    row_indices <- integer(0)
    col_indices <- integer(0)
    for (i in 1:(K-1)) {
      for (j in (i+1):K) {
        row_indices <- c(row_indices, i)
        col_indices <- c(col_indices, j)
      }
    }
    pHat[cbind(row_indices, col_indices)] <- pHatVec
    pHat <- pHat + t(pHat)

    # 计算 E_nu 并更新 pHat 的对角线
    E_nu <- mean(sapply(t, function(j) approximateEnu(n, K, j, Y, finalCliques, clique_size)))
    diag(pHat) <- E_nu

    # check and calculate dHat
    print("this is how many pHat == 0")
    print(sum(pHat == 0))
    if (all(pHat != 0)) {
      dHat <- -log(pHat / E_nu)
      return(list(dHat, finalCliques))  # 返回 dHat 和 finalCliques
    }
  }
}


#' GenerateY
#'
#' @param K number of cliques
#' @param n size of graph
#' @param probs probability to create an edge
#' @param clique_size size of clique
#' @param E_nu a number
#' @param numTaken determines how many nodes to take in one group to find cliques
#' @param t number of edges used in the "almost clique" to compute FE expectation
#'
#' @return Y, pHat, dHat, finalCliques, E_nu
#' @export
#' @importFrom igraph graph_from_adjacency_matrix
GenerateY <- function(K, n, probs, clique_size, E_nu, numTaken, t) {
  nodesPerClique <- as.integer(n / K)
  tempCliques <- vector("list", K)
  for (j in 1:K) {
    tempCliques[[j]] <- list()
  }
  finalCliques <- matrix(0, nrow = K, ncol = clique_size)


  while(TRUE) {
    Y <- GenerateGraph(n, probs) # 假设 GenerateGraph 是另一个您有的函数
    pHat <- matrix(0, nrow = K, ncol = K)
    tempCliques <- vector("list", K)
    for (j in 1:K) {
      tempCliques[[j]] <- list()
    }
    found_empty <- FALSE

    for (indexClique in 1:K) {
      # tempIndices is same as in python
      tempIndices <- (nodesPerClique * (indexClique - 1)) + (0:(numTaken-1))
      Ytemp <- Y[tempIndices+1, tempIndices+1]

      # Convert matrix to graph and find cliques，（调整cliques顺序后和python一致，且内容一致，即含0）
      g <- graph_from_adjacency_matrix(as.matrix(Ytemp), mode = "undirected", diag = FALSE)
      cliques_with_ids <- cliques(g)
      cliques <- sort_cliques(cliques_with_ids)
      cliques <- lapply(cliques, function(clique) {
        sapply(clique, function(node) {
          node - 1
        })
      })

      # Filter cliques of the desired size

      tempCliques[[indexClique]] <- Filter(function(k) {length(k) == clique_size}, cliques)

      # Print the number of cliques found
      print(length(tempCliques[[indexClique]]))

      # Check for empty cliques and potentially break the loop
      if (length(tempCliques[[indexClique]]) == 0) {
        found_empty <- TRUE
        break
      }
    }

    if(found_empty) {
      next
    }


    # 循环遍历
    for (indexClique in 1:K) {
      # 更新 finalCliques 的每一行
      finalCliques[indexClique, ] <- sapply(tempCliques[[indexClique]][[1]], function(j) j + nodesPerClique * (indexClique - 1))
    }

    # same as in python
    tempIndices <- as.integer(c(t(finalCliques)))

    # 创建子图的邻接矩阵
    A_c <- Y[tempIndices + 1, tempIndices + 1]

    # 初始化 c_clique 向量
    c_clique <- rep(0, K * clique_size)

    # 填充 c_clique 向量
    for (i in 1:(K * clique_size)) {
      c_clique[i] <- floor((i - 1) / clique_size)
    }

    # 计算 pHatVec
    pHatVec <- c()  # 创建一个空的向量

    # 使用两个嵌套的 for 循环计算 pHatVec
    for (i in 0:(K-2)) {
      for (j in (i+1):(K-1)) {
        result <- EstimatePhat(i, j, K, clique_size * K, c_clique, A_c)[1] /
          EstimatePhat(i, j, K, clique_size * K, c_clique, A_c)[2]
        # Append the result to pHatVec
        pHatVec <- c(pHatVec, result)
      }
    }

    # 更新 pHat 矩阵
    # Assign values in row-major order (the default assigning in R is col-major order)
    row_indices <- integer(0)
    col_indices <- integer(0)
    for (i in 1:(K-1)) {
      for (j in (i+1):K) {
        row_indices <- c(row_indices, i)
        col_indices <- c(col_indices, j)
      }
    }
    pHat[cbind(row_indices, col_indices)] <- pHatVec
    pHat <- pHat + t(pHat)

    # 计算 E_nu 并更新 pHat 的对角线
    E_nu <- mean(sapply(t, function(j) approximateEnu(n, K, j, Y, finalCliques, clique_size)))
    diag(pHat) <- E_nu

    # 检查并计算 dHat
    print("this is how many pHat == 0")
    print(sum(pHat == 0))
    if (all(pHat != 0)) {
      dHat <- -log(pHat / E_nu)
      return(list(Y = Y, pHat = pHat, dHat = dHat, finalCliques = finalCliques, E_nu = E_nu))
    }
  }
}


#' approximateEnu
#'
#' @param n number of vertices in graph
#' @param K number of cliques
#' @param degCutoff a number
#' @param Y adjacent matrix
#' @param clique a matrix
#' @param clique_size a number
#'
#' @return a number
#' @export
#'
approximateEnu <- function(n, K, degCutoff, Y, clique, clique_size) {
  nodesPerClique <- as.integer(n / K)
  approx <- numeric(K)
  for (indexClique in 1:K) {
    # tempIndices is same as python
    tempIndices <- nodesPerClique * (indexClique - 1) + (0:(nodesPerClique-1))

    # if clique is a matrix, then use [indexClique,]
    cliqueIndices <- sapply(clique[indexClique,], function(x) as.integer(x))
    # if clique is a list， then use [[indexClique]], for PickCliqueFeasible function
    # cliqueIndices <- sapply(clique[[indexClique]], function(x) as.integer(x))

    # assuming clique is in python style, which means may contain 0
    YtempClique <- Y[cliqueIndices + 1, tempIndices + 1]
    # node is same as in python
    nodes <- numeric(0)
    for (j in 0:(nodesPerClique - 1)) {
      if (sum(YtempClique[, j + 1]) == degCutoff &&
          sum(YtempClique[, j + 1]) < clique_size &&
          !(j %in% cliqueIndices)) {
        nodes <- c(nodes, nodesPerClique * (indexClique-1) + j)
      }
    }

    if (length(nodes) < 3) {
      approx[indexClique] <- 0
    } else {
      approx[indexClique] <- sum(Y[nodes+1, nodes+1]) / (length(nodes) * (length(nodes) - 1))
    }
  }

  if (sum(approx) == 0) {
    return(0)
  } else {
    return(sum(approx) / length(which(approx > 0)))
  }
}


#=========================================================
## 目前假设cliques 是一个矩阵
## input: same as in python
## output: same as in python

ComputeDhat <- function(Y, K, clique_size, t, cliques, n) {
  pHat <- matrix(0, nrow = K, ncol = K)
  # 假设 cliques 是一个矩阵
  tempIndices <- as.integer(as.vector(t(cliques)))
  # 假设 cliques 是一个list of list
  # tempIndices <- as.integer(unlist(cliques))

  A_c <- Y[tempIndices + 1, tempIndices + 1]      # adjacency matrix for subgraph induced by cliques
  c_clique <- floor((0:(K*clique_size - 1)) / clique_size)


  pHatVec <- numeric()
  for (i in 0:(K-1)) {
    for (j in 0:(K-1)) {
      if(i < j ) {
        result <- EstimatePhat(i, j, K, clique_size * K, c_clique, A_c)
        pHatVec <- c(pHatVec, result[1] / result[2])
      }
    }
  }

  # Assign values in row-major order (the default assigning in R is col-major order)
  row_indices <- integer(0)
  col_indices <- integer(0)
  for (i in 1:(K-1)) {
    for (j in (i+1):K) {
      row_indices <- c(row_indices, i)
      col_indices <- c(col_indices, j)
    }
  }
  pHat[cbind(row_indices, col_indices)] <- pHatVec
  pHat <- pHat + t(pHat)

  results <- sapply(t, function(i) approximateEnu(n, K, i, Y, cliques, clique_size))
  positive_results <- results[results > 0]

  # 计算平均值，如果没有正数则返回NA
  E_nu <- if (length(positive_results) > 0) mean(positive_results) else NA
  if (is.na(E_nu)) {  # 如果E_nu是NA，则将E_nu设置为1
    E_nu <- 1
  }

  diag(pHat) <- E_nu

  dHat <- matrix(0, nrow = K, ncol = K)
  for (i in 1:K) {
    for (j in 1:K) {
      if (i != j) {
        dHat[i, j] <- -log(pHat[i, j] / E_nu)
      }
    }
  }

  cat(pHat, "\n")
  return(list(dHat, E_nu))
}


#' PickCliquesFeasible
#'
#' @param n number of vertices of the graph
#' @param K a number
#' @param Y adjacent matrix
#' @param clique list of lists
#' @param numSamples number of cliques
#' @param clique_size clique size
#' @param t number of edges used in the "almost clique" to compute FE expectation
#'
#' @return list of optimal indices, minimum objFun, True/False of whether min(objFun) < 100, dHat, cliques, Enu
#' @export
#'
PickCliquesFeasible <- function(n, K, Y, clique, numSamples, clique_size, t) {
  objFun <- rep(0, numSamples)
  indices <- matrix(0, nrow = numSamples, ncol = K)
  cliques_saved <- vector("list", numSamples)
  for (i in 1:numSamples) {
    cliques_saved [[i]] <- list()
  }
  dHat_saved <- vector("list", numSamples)
  for (i in 1:numSamples) {
    dHat_saved [[i]] <- list()
  }
  Enu_saved <- rep(0, numSamples)

  for (index in 1:numSamples) {
    print(index) # the print output is R index, for python output: print(index-1)
    pHat <- matrix(0, nrow = K, ncol = K)
    tempSum <- 0
    ## a is same as python
    a <- rep(1, K)
    while (length(a) > length(unique(a))) {
      a <- sample(1:(length(clique) - 1), K, replace = TRUE)
    }

    for (alpha in 1:K) {
      for (beta in 1:K) {
        if (alpha < beta) {
          tempSum <- tempSum + length(intersect(clique[[a[alpha]+1]], clique[[a[beta]+1]]))
        }
      }
    }

    cliques_sampled <- vector("list", length = K)
    for (j in 1:K) {
      cliques_sampled[[j]] <- clique[[a[j] + 1]]
    }


    indices[index, ] <- a
    cliques_saved[[index]] <- do.call(rbind, cliques_sampled)
    # tempIndices is same as in python
    tempIndices <- unlist(cliques_sampled)
    A_c <- Y[tempIndices+1, tempIndices+1]
    c_clique <- floor((0:(K*clique_size - 1)) / clique_size)

    pHatVec <- numeric()
    for (i in 0:(K-1)) {
      for (j in 0:(K-1)) {
        if (i < j) {
          estimate <- EstimatePhat(i, j, K, clique_size * K, c_clique, A_c)
          pHatVec <- c(pHatVec, estimate[1] / estimate[2])
        }
      }
    }

    # Assign values in row-major order (the default assigning in R is col-major order)
    row_indices <- integer(0)
    col_indices <- integer(0)
    for (i in 1:(K-1)) {
      for (j in (i+1):K) {
        row_indices <- c(row_indices, i)
        col_indices <- c(col_indices, j)
      }
    }
    pHat[cbind(row_indices, col_indices)] <- pHatVec
    pHat <- pHat + t(pHat)

    results <- sapply(t, function(i) approximateEnu(n, K, i, Y, do.call(rbind, cliques_sampled), clique_size))
    positive_results <- results[results > 0]
    # 计算平均值，如果没有正数则返回NA
    E_nu <- if (length(positive_results) > 0) mean(positive_results) else NA
    if (is.na(E_nu)) {  # 如果E_nu是NA，则将E_nu设置为1
      E_nu <- 1
    }

    Enu_saved[index] <- E_nu
    diag(pHat) <- E_nu
    dHat <- matrix(0, nrow = K, ncol = K)

    for (i in 1:K) {
      for (j in 1:K) {
        if (i != j) {
          dHat[i, j] <- -log(pHat[i, j] / E_nu)
        }
      }
    }

    dHat_saved[[index]] <- dHat  # 保存dHat
    objFun[index] <- tempSum + 100 * as.numeric(any(pHat == 0))

  }

  index_optim <- which.min(objFun)  # 找到objFun的最小值的索引
  return(list(indices[index_optim, ], min(objFun), min(objFun) < 100, dHat_saved[[index_optim]], cliques_saved[[index_optim]], Enu_saved[index_optim]))
}
