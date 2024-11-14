# Introduction to igraph’s community detection algorithms

library(igraph)
library(igraphdata)
library(clustAnalytics)

# List of data sets in igraphdata
#data(package="igraphdata")

# Karate example
data(karate, package = "igraphdata")
karate <- upgrade_graph(karate)

plot(karate)

# Clustering
wc <- cluster_walktrap(karate)
modularity(wc)
membership(wc)
unname(membership(wc))
plot(wc, karate)

plot(karate, vertex.color=membership(wc))

# Hierarchical clustering
fc <- cluster_fast_greedy(karate)
plot_dendrogram(fc)

GN <- cluster_edge_betweenness(karate)
plot_dendrogram(GN)

modularity(GN); modularity(fc)

# Adjacency matrix
as_adjacency_matrix(as_undirected(karate,mode = "each"))

data("UKfaculty", package="igraphdata")
UKfaculty <- upgrade_graph(UKfaculty)
plot(UKfaculty)
as_adjacency_matrix(UKfaculty)
as_adjacency_matrix(as_undirected(UKfaculty, mode="each"))

# Evaluate significance of a few clusterings through various scoring functions
evaluate_significance(karate
                      , alg_list = list(Louvain = cluster_louvain
                                        , "label prop"= cluster_label_prop
                                        , walktrap=cluster_walktrap)
                      , gt_clustering=V(karate)$Faction)

# Generate a synthetic ground truth model
B <- matrix(c(1, 0.2, 0.2, 1), ncol=2)
G <- barabasi_albert_blocks(m=4, p=c(0.5, 0.5), B=B, t_max=100
                            , type="Hajek", sample_with_replacement = FALSE)
plot(G, vertex.color=(V(G)$label),vertex.label=NA,vertex.size=10)

