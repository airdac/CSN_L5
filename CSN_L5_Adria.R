library(igraph)
library(clustAnalytics)

data(karate, package = "igraphdata")

plot(karate)

wc <- walktrap.community(karate)
modularity(wc)
membership(wc)
unname(membership(wc))
plot(wc, karate)

plot(karate, vertex.color=membership(wc))

fc <- fastgreedy.community(karate)
dendPlot(fc)

GN <- edge.betweenness.community(karate)
dendPlot(GN)

modularity(GN); modularity(fc)

evaluate_significance(karate
                      , alg_list = list(Louvain = cluster_louvain
                                        , "label prop"= cluster_label_prop
                                        , walktrap=cluster_walktrap)
                      , gt_clustering=V(karate)$Faction)

B <- matrix(c(1, 0.2, 0.2, 1), ncol=2)
G <- barabasi_albert_blocks(m=4, p=c(0.5, 0.5), B=B, t_max=100
                            , type="Hajek", sample_with_replacement = FALSE)
plot(G, vertex.color=(V(G)$label),vertex.label=NA,vertex.size=10)

