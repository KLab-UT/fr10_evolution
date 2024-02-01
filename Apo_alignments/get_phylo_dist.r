library(ape)
library(adephylo)

setwd('/Users/r_klabacka/OneDrive - Utah Tech University/KLab/Research/fr10_evolution/Apo_alignments')

# import tree
apoAII_tree <- read.tree('ApoA-II_aligned.contree')

# plot tree (just to check it)
plot(apoAII_tree)

# get pairwise phylogenetic distance
pairwise_dist <- distTips(apoAII_tree, tips='all', 'patristic')





