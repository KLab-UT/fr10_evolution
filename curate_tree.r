library(ape)
library(phytools)

fr10_tree <- read.tree('U44831_in_8342_fixed.fasta.contree')
fr10_tree_rooted <- root(fr10_tree, outgroup = "Pipa_carvalhoi")
nodelabels(fr10_tree_rooted$node.label,node=2:fr10_tree_rooted$Nnode+Ntip(fr10_tree_rooted),
    adj=c(1,-0.2),frame="none")

plotTree(fr10_tree_rooted, show.node.label = TRUE)
