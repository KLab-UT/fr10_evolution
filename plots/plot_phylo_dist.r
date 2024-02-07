library(ape)
library(adephylo)

normalize_distances <- function(tree_file, normalized_distances_fout, query_species){# load tree
  # Load input and set up variables
  tree <- read.tree(tree_file)
  fout <- read.csv(normalized_distances_fout)
  gene_name <- sub("_aligned.plottree", "", basename(tree_file))
  distmatrix_file <- sub(".treefile", "_distmatrix.csv", tree_file)
  normalized_gene_fout <- sub("_aligned.plottree", "_normalized.csv", basename(tree_file))
  
  # get pairwise phylogenetic distance
  pairwise_dist <- distTips(tree, tips='all', 'patristic')
  pd_matrix <- as.matrix(pairwise_dist)
  
  # write distances to csv file
  pairwise_dist_df <- as.data.frame(pd_matrix)
  write.csv(pairwise_dist_df, distmatrix_file, row.names = TRUE)
  
  # read distmatrix
  distmatrix <- read.csv(distmatrix_file)
  
  # Check if Xenopus or Silurana
  if ("Xenopus_tropicalis" %in% colnames(distmatrix)) {
    sister_species <- "Xenopus_tropicalis"
  } else if ("Silurana_tropicalis" %in% colnames(distmatrix)){
    sister_species <- "Silurana_tropicalis"
  } else {
    print("Sister species not found in distmatirx")
  }
  
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
    
    # Sanitize distances dataframe
    sanitized_distances <- data.frame(matrix(NA, nrow = 1, ncol = ncol(fout)))
    
    cat("Species names in normalized distances output:", "\n")
    print(colnames(fout))
    colnames(sanitized_distances) <- colnames(fout)
    
    # Loop through species in normalized_gene_distances
    for (species in colnames(normalized_gene_distances)) {
      # Check if the species is present in sanitized_distances
      if (species %in% colnames(sanitized_distances)) {
        # Update the value in sanitized_distances with the value from normalized_gene_distances
        sanitized_distances[, species] <- normalized_gene_distances[, species]
      }
      else if (species == "Ornithorhynchus_anatini") {
        species <- "Ornithorhynchus_anatinus"
      }
      else if (species == "Chrysemys_picta") {
        species <- "picta_bellii"
      }
      else if (species == "Equus_caballus") {
        species <- "Equus_przewalskii"
      }
      
      
      else {
        # Print the species that is not present in sanitized_distances
        cat("Species not present in sanitized_distances:", species, "\n")
        
        cat("Species in gene normalized output:", "\n")
        print(colnames(sanitized_distances))
      }
    }
    

    rownames(sanitized_distances) <- gene_name
    
    write.csv(sanitized_distances, normalized_gene_fout, row.names = TRUE)
    
    
    # Update row name
    #rownames(normalized_gene_distances)[rownames(normalized_gene_distances) == query_species] <- gene_name
    
    # Append normalized row to fout
    #fout <- rbind(fout, sanitized_distances)
    
    # Save the updated normalized_distances_fout
    #write.csv(fout, normalized_distances_fout, row.names = TRUE)
    
  } 
  else {
    cat(paste(sister_species, " is not a column in distmatrix. Skipping the normalization process.\n"))
    print(colnames(distmatrix))
  }
}



treefile <- as.character(commandArgs(trailingOnly = TRUE)[1])
normalized_distances_fout <- as.character(commandArgs(trailingOnly = TRUE)[2])
query_species <- as.character(commandArgs(trailingOnly = TRUE)[3])

#treefile <- "../apo-fr10_alignments/ApoA-II_aligned.plottree"
#normalized_distances_fout <- "normalized_distances_fr10.csv"
#query_species <- "Lithobates_sylvaticus"

setwd("/uufs/chpc.utah.edu/common/home/u6052680/fr10_evolution/plots")

normalize_distances(treefile, normalized_distances_fout, query_species)




