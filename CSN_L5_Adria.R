# UPC. FIB. 2024-25, fall semester.
# CSN. Lab work 5. Finding and assessing community structure
# Adrià Casanova, Dmitriy Chukhray

library(igraph)
suppressPackageStartupMessages(library(igraph))
library(clustAnalytics)

set.seed(42)

#####################################################
# Load data
#####################################################

# karate
data(karate, package = "igraphdata")
karate <- upgrade_graph(karate)
plot(karate, main = "Karate Network")

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
plot(wc, synthetic, main = "Synthetic Network")

# ENRON
# Transform to an undirected simple weighted graph
data("enron", package = "igraphdata")
enron <- upgrade_graph(enron)
enron <- as.undirected(enron, mode = "each")

# Count the number of edges between each pair of nodes
edges <- ends(enron, E(enron)) # Get node pairs of all edges
edges_names <- apply(edges, 1, function(x) paste(sort(x), collapse = "-"))
edge_counts <- table(edges_names)

# Assign counts to edges' weights
E(enron)$weight <- as.integer(edge_counts[edges_names])

# Remove self-loops and multiedges
enron <- simplify(enron)

# Check code worked
head(E(enron)$weight)
plot(enron, vertex.size = 7, main = "ENRON Network")

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
  n_clusters1 <- unique(clustering1)
  n_clusters2 <- unique(clustering2)

  # Compute Jaccard index for a pair of clusters
  jaccard_index <- function(c1, c2) {
    cluster1_nodes <- which(clustering1 == c1)
    cluster2_nodes <- which(clustering2 == c2)

    intersection <- length(intersect(cluster1_nodes, cluster2_nodes))
    union <- length(union(cluster1_nodes, cluster2_nodes))

    return(intersection / union)
  }

  # Fill the output table
  jaccard_matrix <- outer(n_clusters1, n_clusters2, Vectorize(jaccard_index))
  rownames(jaccard_matrix) <- n_clusters1
  colnames(jaccard_matrix) <- n_clusters2

  return(jaccard_matrix)
}

lc <- cluster_louvain(karate)
lc1 <- unname(membership(lc))
lc2 <- V(karate)$Faction

jaccard_table_l <- jaccard_sim(lc1, lc2)
print("Jaccard index between clusters created by Louvain algorithm and ground truth clusters")
print(jaccard_table_l)

lpc <- cluster_label_prop(karate)
lpc1 <- unname(membership(lpc))
lpc2 <- V(karate)$Faction

jaccard_table_lp <- jaccard_sim(lpc1, lpc2)
print("Jaccard index between clusters created by Label Propagation algorithm and ground truth clusters")
print(jaccard_table_lp)


wc <- cluster_walktrap(karate)
c1 <- unname(membership(wc))
c2 <- V(karate)$Faction

jaccard_table_w <- jaccard_sim(c1, c2)
print("Jaccard index between clusters created by Walktrap algorithm and ground truth clusters")
print(jaccard_table_w)

ebc <- cluster_edge_betweenness(karate)
ebc1 <- unname(membership(ebc))
ebc2 <- V(karate)$Faction

jaccard_table_eb <- jaccard_sim(ebc1, ebc2)
print("Jaccard index between clusters created by Edge Betweenness algorithm and ground truth clusters")
print(jaccard_table_eb)


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

matched_clusters_l <- match_clusters(jaccard_table_l, "LC", "GT")
print("Most similar cluster in clustering 2 (GT) to each cluster of clustering 1 (LC)")
print(matched_clusters_l)

matched_clusters_lp <- match_clusters(jaccard_table_lp, "LPC", "GT")
print("Most similar cluster in clustering 2 (GT) to each cluster of clustering 1 (LPC)")
print(matched_clusters_lp)

matched_clusters_w <- match_clusters(jaccard_table_w, "WC", "GT")
print("Most similar cluster in clustering 2 (GT) to each cluster of clustering 1 (WC)")
print(matched_clusters_w)

matched_clusters_eb <- match_clusters(jaccard_table_eb, "EBC", "GT")
print("Most similar cluster in clustering 2 (GT) to each cluster of clustering 1 (EBC)")
print(matched_clusters_eb)


##############################
# 3. Global Jaccard similarity
##############################

# Weighted mean
Wmean <- function(v, w) {
  return(as.numeric(v %*% w))
}

# Average of jaccard indices weighted by fraction of number of nodes
# Compute fraction of number of nodes
w <- table(lc1) / length(lc1)
global_jaccard_sim_l <- Wmean(matched_clusters_l, w)
print("Global Jaccard similarity of clusters produced by Louvain algorithm")
print(global_jaccard_sim_l)

w <- table(lpc1) / length(lpc1)
global_jaccard_sim_lp <- Wmean(matched_clusters_lp, w)
print("Global Jaccard similarity of clusters produced by Label Propagation algorithm")
print(global_jaccard_sim_lp)

w <- table(c1) / length(c1)
global_jaccard_sim_w <- Wmean(matched_clusters_w, w)
print("Global Jaccard similarity of clusters produced by Walktrap algorithm")
print(global_jaccard_sim_w)

w <- table(ebc1) / length(ebc1)
global_jaccard_sim_eb <- Wmean(matched_clusters_eb, w)
print("Global Jaccard similarity of clusters produced by Edge Betweenness algorithm")
print(global_jaccard_sim_eb)

# The weighted mean is a more reasonable similarity than the mean value because
# it accounts for the size of each cluster, giving more importance to larger clusters
# This reflects the idea that larger clusters have a greater impact on the overall
# structure of the network.

# Another way of combining the vector of Jaccard indices to quantify clusterings
# similarity is taking the harmonic mean, which mitigates the impact of large
# outliers and aggravates the impact of small ones. Hence, very similar clusters
# are not so important when using the harmonic mean, while very distinct clusters
# penalize a lot the global Jaccard index. Basically, harmonic mean of is useful when
# we want to ensure that all clusters, regardless of size, are equally considered, 
# and when it's important to identify and emphasize any poor matches between clusters.


# Function to compute harmonic mean
Hmean <- function(v) {
  return(length(v) / sum(1 / v))
}

# Calculate Harmonic Mean Global Jaccard Similarity
harmonic_jaccard_sim_l <- Hmean(matched_clusters_l)
print("Harmonic Jaccard similarity of clusters produced by Louvain algorithm")
print(harmonic_jaccard_sim_l)

harmonic_jaccard_sim_lp <- Hmean(matched_clusters_lp)
print("Harmonic Jaccard similarity of clusters produced by Label Propagation algorithm")
print(harmonic_jaccard_sim_lp)

harmonic_jaccard_sim_w <- Hmean(matched_clusters_w)
print("Harmonic Jaccard similarity of clusters produced by Walktrap algorithm")
print(harmonic_jaccard_sim_w)

harmonic_jaccard_sim_eb <- Hmean(matched_clusters_eb)
print("Harmonic Jaccard similarity of clusters produced by Edge Betweenness algorithm")
print(harmonic_jaccard_sim_eb)



##############################
# 4.Evaluate the significance of clusterings produced by different community detection methods
##############################

# The first thing we have to do is decide which metrics from evaluate_significance are relevant.
# Then, based on these metrics we have to choose the best ranked clustering for data with no ground truth.
# We decided to choose modularity and coverage metrics.

# Define clustering algorithms
algorithms <- list(
  Louvain = cluster_louvain,
  LabelPropagation = cluster_label_prop,
  Walktrap = cluster_walktrap,
  EdgeBetweenness = cluster_edge_betweenness
)

# Networks without ground truth
networks_without_gt <- list(
  Enron = enron,
  Dolphins = dolphins
)

# Initialize a list to store modularity and coverage results
modularity_coverage_results <- list()

# Apply evaluate_significance once to each network without ground truth
for (network_name in names(networks_without_gt)) {
  graph <- networks_without_gt[[network_name]]
  
  cat("\nEvaluating", network_name, "Network:\n")
  
  # Apply evaluate_significance
  res <- evaluate_significance(graph, alg_list = algorithms)
  
  # Extract modularity and coverage for each algorithm
  modularity_scores <- res[14, ]  # Row 14 contains modularity
  coverage_scores <- res[19, ]   # Row 19 contains coverage
  
  # Store results in the main list
  modularity_coverage_results[[network_name]] <- list(
    modularity = modularity_scores,
    coverage = coverage_scores
  )
  
  # Print modularity and coverage for debugging
  cat("Modularity Scores:", modularity_scores, "\n")
  cat("Coverage Scores:", coverage_scores, "\n")
}

modularity_coverage_results$Enron$modularity
modularity_coverage_results$Enron$coverage

# For Enron data we have to go with Louvain produced clusters as the ground truth
# because its best modularity score (0.2606) is produced by the Louvain algorithm.
# Coverage metric is the best using Edge Betweenness algorithm (0.9945), but it is not much higher
# than under the Louvain algorithm (0.9777). If we compare the Edge Betweenness score in modularity 
# it is ~ 8.5 times lower than that of Louvain. Therefore, we can assume clusters produced by the Louvain algorithm
# are the closest to the ground truth

modularity_coverage_results$Dolphins$modularity
modularity_coverage_results$Dolphins$coverage

# For Dolphins data we have to go with Edge Betweenness produced clusters as the ground truth.
# In both modularity and coverage scores it ranks as the best (0.5193 and 0.7987 respectively).

# CAREFUL! LabelPropagation has a higher coverage than EdgeBetweenness, although it has a lower modularity as well...

# Define networks and assign specific ground truths
networks <- list(
  Karate = list(graph = karate, ground_truth = V(karate)$Faction),
  Synthetic = list(graph = synthetic, ground_truth = V(synthetic)$label),
  Enron = list(graph = enron, ground_truth = membership(cluster_louvain(enron))),  # Louvain as ground truth
  Dolphins = list(graph = dolphins, ground_truth = membership(cluster_edge_betweenness(dolphins)))  # Edge Betweenness as ground truth
)

# Function to compute Jaccard indices
compute_jaccard <- function(clustering, ground_truth) {
  jaccard_table <- jaccard_sim(clustering, ground_truth)  
  weights <- table(clustering) / length(clustering)       
  global_jaccard <- Wmean(apply(jaccard_table, 1, max), weights)  
  list(local = jaccard_table, global = global_jaccard)
}

# Initialize results
results <- list()

for (network_name in names(networks)) {
  network <- networks[[network_name]]
  graph <- network$graph
  ground_truth <- network$ground_truth
  
  cat("\nEvaluating", network_name, "Network:\n")
  
  # Evaluate clustering significance
  res <- evaluate_significance(graph, alg_list = algorithms, gt_clustering = ground_truth)
  
  # Compute Jaccard indices for each algorithm
  jaccard_indices <- lapply(algorithms, function(algo) {
    clustering <- membership(algo(graph))
    compute_jaccard(clustering, ground_truth)
  })
  
  # Store results
  results[[network_name]] <- list(significance = res, jaccard = jaccard_indices)
  print(res)
  
  # Print Jaccard indices
  cat("\nJaccard Indices:\n")
  for (algo_name in names(jaccard_indices)) {
    jaccard <- jaccard_indices[[algo_name]]
    cat("\nAlgorithm:", algo_name, "\n")
    cat("Global Jaccard Index:", jaccard$global, "\n")
    cat("Local Jaccard Table:\n")
    print(jaccard$local)
  }
}




























