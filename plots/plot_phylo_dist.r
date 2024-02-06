library(ape)
library(adephylo)

normalize_distances <- function(tree_file, normalized_distances_fout, query_species, sister_species){# load tree
  # Load input and set up variables
  tree <- read.tree(tree_file)
  fout <- read.csv(normalized_distances_fout)
  gene_name <- sub("_aligned.plottree", "", basename(tree_file))
  distmatrix_file <- sub(".treefile", "_distmatrix.csv", tree_file)
  
  # get pairwise phylogenetic distance
  pairwise_dist <- distTips(tree, tips='all', 'patristic')
  pd_matrix <- as.matrix(pairwise_dist)
  
  # write distances to csv file
  pairwise_dist_df <- as.data.frame(pd_matrix)
  write.csv(pairwise_dist_df, distmatrix_file, row.names = TRUE)
  
  # read distmatrix
  distmatrix <- read.csv(distmatrix_file)
  
  # Check if sister_species is a column in distmatrix
  if (sister_species %in% colnames(distmatrix)) {
    # get distance between query species and sister species 
    qs_dist <- distmatrix[distmatrix$X == query_species, sister_species]
    
    # get the distances between the query and all other species not including query or sister
    q_distances <- distmatrix[distmatrix$X != query_species & distmatrix$X != sister_species, c("X", query_species), drop = FALSE]
    
    # Change the "X" column to row names
    rownames(q_distances) <- q_distances$X
    
    # Remove the "X" column
    q_distances <- q_distances[, -1, drop=FALSE]
    
    ## Normalize q_distances
    normalized_gene_distances <- q_distances / qs_dist 
    
    # Set row names for the normalized_distances dataframe
    rownames(normalized_gene_distances) <- rownames(q_distances)
    
    # Transpose the dataframe
    normalized_gene_distances <- t(normalized_gene_distances)
    
    # Update row name
    rownames(normalized_gene_distances)[rownames(normalized_gene_distances) == query_species] <- query_gene
    
    # Append normalized row to fout
    fout <- rbind(fout, normalized_gene_distances)
    
    # Identify columns in fout that are not in normalized_gene_distances
    missing_cols <- setdiff(colnames(fout), colnames(normalized_gene_distances))
    
    # Set values in the normalized row to NA for missing columns
    fout[nrow(fout), missing_cols] <- NA
    
    # Save the updated normalized_distances_fout
    write.csv(fout, normalized_distances_fout, row.names = TRUE)
    
  } 
  else {
    cat("sister_species is not a column in distmatrix. Skipping the normalization process.\n")
  }
}



treefile <- as.character(commandArgs(trailingOnly = TRUE)[1])
normalized_distances_fout <- as.character(commandArgs(trailingOnly = TRUE)[2])
query_species <- as.character(commandArgs(trailingOnly = TRUE)[3])
sister_species <- as.character(commandArgs(trailingOnly = TRUE)[4])

#tree_file <- "apo-fr10_alignments/ApoC-II_aligned.plottree"
#normalized_distances_fout <- "plots/normalized_distances_fr10.csv"
#query_species <- "Lithobates_sylvaticus"
#sister_species <- "Xenopus_tropicalis"
setwd("/home/reagan/bioinformatics/fr10_evolution/plots")

normalize_distances(treefile, normalized_distances_fout, query_species, sister_species)


