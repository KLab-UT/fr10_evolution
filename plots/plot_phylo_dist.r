library(ape)
library(adephylo)

#treefile <- as.character(commandArgs(trailingOnly = TRUE)[1])
#distmatrix_file <- as.character(commandArgs(trailingOnly = TRUE)[2])
#query_species <- as.character(commandArgs(trailingOnly = TRUE)[3])
#sister_species <- as.character(commandArgs(trailingOnly = TRUE)[4])


distmatrix_file <- "apo-fr10_alignments/ApoA-I_aligned_distmatrix.csv"
query_species <- "Lithobates_sylvaticus"
sister_species <- "Xenopus_tropicalis"
query_gene <-"ApoA-I"


normalize_distances <- function(tree_file, query_species, sister_species){# load tree
  tree <- read.tree(tree_file)
  gene_name <- sub("_aligned.plottree", "", basename(tree_file))
  distmatrix_file <- sub(".treefile", "_distmatrix.csv", tree_file)
  print(tree_file)
  
  # get pairwise phylogenetic distance
  pairwise_dist <- distTips(tree, tips='all', 'patristic')
  pd_matrix <- as.matrix(pairwise_dist)
  
  # write distances to csv file
  pairwise_dist_df <- as.data.frame(pd_matrix)
  write.csv(pairwise_dist_df, distmatrix_file, row.names = TRUE)
  
  # read distmatrix
  distmatrix <- read.csv(distmatrix_file)
  
  # get distance between query species and sister species 
  qs_dist <- distmatrix[distmatrix$X == query_species, sister_species]
  
  # get the distances between the query and all other species not including query or sister
  q_distances <- distmatrix[distmatrix$X != query_species & distmatrix$X != sister_species, c("X", query_species), drop = FALSE]
  
  # Change the "X" column to row names
  rownames(q_distances) <- q_distances$X
  
  # Remove the "X" column
  q_distances <- q_distances[, -1, drop=FALSE]
  
  ##TODO: normalize q_distances
  normalized_distances <- q_distances / qs_dist 
  
  # Set row names for the normalized_distances dataframe
  rownames(normalized_distances) <- rownames(q_distances)
  
  # Transpose the dataframe
  normalized_distances <- t(normalized_distances)
  
  # Update row name
  rownames(normalized_distances)[rownames(normalized_distances) == query_species] <- query_gene
  return(normalized_distances)
}

tree_file <- "apo-fr10_alignments/ApoC-II_aligned.plottree"
query_species <- "Lithobates_sylvaticus"
sister_species <- "Xenopus_tropicalis"

result <- normalize_distances(tree_file, query_species, sister_species)




tree_files <- c("apo-fr10_alignments/ApoA-II_aligned.plottree", "apo-fr10_alignments/ApoA-V_aligned.plottree", "apo-fr10_alignments/ApoC-IV_aligned.plottree",
                "apo-fr10_alignments/ApoA-IV_aligned.plottree", "apo-fr10_alignments/ApoC-III_aligned.plottree", "apo-fr10_alignments/ApoC-I_aligned.plottree",
                "apo-fr10_alignments/ApoA-I_aligned.plottree", "apo-fr10_alignments/ApoC-II_aligned.plottree", "apo-fr10_alignments/АроЕ_aligned.plottree")
query_species <- "Lithobates_sylvaticus"
sister_species <- "Xenopus_tropicalis"

# Initialize an empty dataframe with column names
result_df <- data.frame()

# Loop through each pair and append to the result_df
for (i in seq_along(tree_files)) {
  pair_result <- normalize_distances(tree_files[i], query_species, sister_species)
  result_df <- rbind(result_df, pair_result)
}
