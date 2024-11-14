# UPC. FIB. 2024-25, fall semester.
# CSN. Lab work 5. Finding and assessing community structure
# Adrià Casanova, Dmitriy Chukhray

library(igraph)
library(igraphdata)
library(clustAnalytics)

#####################################################
# Load data
#####################################################

data(karate, package = "igraphdata")
karate <- upgrade_graph(karate)

# Diagonal(B) = probability of connecting between clusters
# B has to be 4x4 instead
# G must have 200 nodes and 800 edges
# Code for introduction's example
B <- matrix(c(1, 0.2, 0.2, 1), ncol=2)
G <- barabasi_albert_blocks(m=4, p=c(0.5, 0.5), B=B, t_max=100
                            , type="Hajek", sample_with_replacement = FALSE)

# Transform to an undirected simple graph (remove weights)
# Use adjacency matrix
data("enron", package = "igraphdata")
enron <- upgrade_graph(enron)

# A network of our choice

#####################################################
# 1.
#####################################################

similarity(karate, method = "jaccard")
