library(lsgeo)
# calculating Type1_Power
P <- list(
  n = 1200,                       # size of graph.
  K_vec = c(5, 7, 10),            # values of K (number of cliques).
  dim = 3,                        # dimension of latent space.
  kappa = 0.75,                   # curvature of latent space.
  numSim = 10, #100,                   # number of sets of cliques that we simulate.
  beta = -0.01,                   # lower bound for distribution of fixed effects.
  ell_vec = c(7),                 # values of ell (clique size).
  alpha = 0.05,                   # significance level of test.
  zeta = 1,                       # scaling in front of distances (default = 1).
  scale = 2.5,                    # scale used when generating points in hyperbolic space.
  ub_theta = pi,                  # upper bound for distribution of theta angle for points in spherical LS.
  ub_phi = 2 * pi,                # upper bound for distribution of phi angle for points in spherical LS.
  delta = 10^(-2),                # controls the spread of the nodes around their group centers in spherical LS.
  rate = 1/3,                     # rate used when sub-sampling. This is the exponent of tau_n = n^rate.
  sigma = 0.8,                    # controls the spread of the points in Euclidean space.
  geometry_true = 'E',            # we simulate graphs from geometry_true.
  geometry_test = 'S',            # we test if the latent space is geometry_test.
  numSamples = 45,                # determines how many nodes to take in one group to find cliques. (larger values slow down the code).
  a = 1/3,                        # lower bound used to search for kappa.
  b = 1.5,                        # upper bound used to search for kappa.
  q = 2,                          # number of eigenvalues we minimize to compute our estimate of kappa.
  m = c(6)                        # sub-sample rate. We take m = clique_size - 1.
)
result <- ComputePerformance(P)
results_dictionary <- list(parameters = P, results = result)

# calculating Type1_Power
library(igraph) 
K <- 12
Y <- as.matrix(read.csv("celegans131matrix.csv", header = FALSE, sep = ","))

# Create a graph from the adjacency matrix
G <- graph_from_adjacency_matrix(Y, mode = "undirected", diag = FALSE)
# Get the number of vertices
n <- vcount(G)
clique_number <- max(sapply( max_cliques(G), length))


clique_size <- 5
cat("this is clique size\n")
print(clique_size)
cliques_array <- cliques(G, min = clique_size, max = clique_size)
cliques_array <- lapply(cliques_array, as.numeric) 
cliques_array <- lapply(cliques_array, function(x) x - 1)
cliques <- sort_cliques(cliques_array)

numSamples <- 10**6
m <- 3 #clique_size - 1

t <- clique_number - 2
while(TRUE) {
  results <- PickCliquesFeasible(n, K, Y, cliques, numSamples, as.integer(clique_size), t)
  indices <- results[[1]]
  min_obj <- results[[2]]
  found_min <- results[[3]]
  dHat <- results[[4]]
  temp <- results[[5]] # is a list of matrix
  E_nu <- results[[6]]
  print(found_min)
  if (found_min) {
    break
  }
  
}

print("out of while loop")

# Function call
results <- Return_Classification(Y, clique_size, K, m, "Not Data", "False", cliques, dHat)

# Create a list (equivalent to a dictionary in Python)
data <- list(
  results = results,
  clique_size = clique_size,
  K = K,
  m = m,
  indices = indices,
  min = min_obj,
  numSamples = numSamples
)
