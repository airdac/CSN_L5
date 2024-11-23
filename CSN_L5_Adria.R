# UPC. FIB. 2024-25, fall semester.
# CSN. Lab work 5. Finding and assessing community structure
# Adrià Casanova, Dmitriy Chukhray

library(igraph)
library(igraphdata)
library(clustAnalytics)

#####################################################
# Load data
#####################################################

# karate
data(karate, package = "igraphdata")
karate <- upgrade_graph(karate)

# synthetic
b <- matrix(
  c(
    1, 0.1, 0.1, 0.1,
    0.1, 1, 0.1, 0.1,
    0.1, 0.1, 1, 0.1,
    0.1, 0.1, 0.1, 1
  ),
  ncol = 4
)
synthetic <- barabasi_albert_blocks(
  m = 4, p = rep(0.25, 4), B = b, t_max = 200,
  type = "Hajek", sample_with_replacement = FALSE
)
# check the generated network has 4 clusters
wc <- cluster_walktrap(synthetic)
unique(unname(membership(wc)))
plot(wc, synthetic)

# ENRON
# Transform to an undirected simple graph (remove weights)
# Use adjacency matrix
data("enron", package = "igraphdata")
enron <- upgrade_graph(enron)
M.enron <- as_adjacency_matrix(as_undirected(enron, mode = "each"))
M.enron[M.enron != 0] <- 1 # Replace all non-zero weights with 1
enron <- graph_from_adjacency_matrix(M.enron, mode = "undirected", weighted = NULL)

# A network of our choice
# Dolphins social network from
# https://github.com/balajisriraj/Network-Analysis-Key-Players-Community---Detection---Dolphins
dolphins <- read_graph("dolphins.gml", format = "gml")
plot(dolphins, vertex.size = 5, main = "Dolphin Social Network")


############################################################
# 1. Jaccard index between clusters of different clusterings
############################################################

# Function to compute Jaccard index between clusters of different clusterings
jaccard_sim <- function(clustering1, clustering2) {
  n_clusters1 <- length(unique(clustering1))
  n_clusters2 <- length(unique(clustering2))

  # Initialize an empty matrix to store the Jaccard indices
  jaccard_matrix <- matrix(0, nrow = n_clusters1, ncol = n_clusters2)
  rownames(jaccard_matrix) <- 1:n_clusters1
  colnames(jaccard_matrix) <- 1:n_clusters2

  # Loop over each pair of clusters
  for (i in 1:n_clusters1) {
    for (j in 1:n_clusters2) {
      # Get the nodes in each cluster
      cluster1_nodes <- which(clustering1 == i)
      cluster2_nodes <- which(clustering2 == j)

      # Compute the Jaccard index
      intersection <- length(intersect(cluster1_nodes, cluster2_nodes))
      union <- length(union(cluster1_nodes, cluster2_nodes))

      jaccard_matrix[i, j] <- intersection / union
    }
  }

  return(jaccard_matrix)
}

# Set up two clusterings on the karate dataset to test the function
wc <- cluster_walktrap(karate)
c1 <- unname(membership(wc))
c2 <- V(karate)$Faction

jaccard_table <- jaccard_sim(c1, c2)
print(jaccard_table)


#########################################################################
# 2. Most similar cluster in clustering 2 to each cluster of clustering 1
#########################################################################

# Function to find most similar cluster in clustering 2 to each cluster of clustering 1
match_clusters <- function(table, name1, name2) {
  clusters <- apply(table, 1, which.max)
  jaccard_indices <- table[cbind(1:nrow(table), clusters)]
  names(jaccard_indices) <- paste0(
    "(", name1, ".", 1:nrow(table), ",",
    name2, ".", clusters, ")"
  )
  return(jaccard_indices)
}

# Test
matched_clusters <- match_clusters(jaccard_table, "WC", "GT")
print(matched_clusters)


##############################
# 3. Global Jaccard similarity
##############################

# Weighted mean
Wmean <- function(v, w) {
  return(as.numeric(v %*% w))
}

# Average of jaccard indices weighted by fraction of number of nodes
# Compute fraction of number of nodes
w <- table(c1) / length(c1)
global_jaccard_sim <- Wmean(matched_clusters, w)
print(global_jaccard_sim)

# The weighted mean is a more reasonable similarity than the mean value because
# it accounts for the size of each cluster, giving more importance to larger clusters
# This reflects the idea that larger clusters have a greater impact on the overall
# structure of the network.

# Another way of combining the vector of Jaccard indices to quantify clusterings
# similarity is taking the harmonic mean, which mitigates the impact of large
# outliers and aggravates the impact of small ones. Hence, very similar clusters
# are not so important when using the harmonic mean, while very distinct clusters
# penalize a lot the global Jaccard index.
