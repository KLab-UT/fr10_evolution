library(ape)
library(adephylo)

treefile <- as.character(commandArgs(trailingOnly = TRUE)[1]))
distmatrix_file <- as.character(commandArgs(trailingOnly = TRUE)[2]))

# load tree
tree <- read.tree(treefile)

# get pairwise phylogenetic distance
pairwise_dist <- distTips(tree, tips='all', 'patristic')
pd_matrix <- as.matrix(pairwise_dist)

# write distances to csv file
pairwise_dist_df <- as.data.frame(pd_matrix)
write.csv(pairwise_dist_df, distmatrix_file, row.names = TRUE)