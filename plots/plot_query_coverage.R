# Load the necessary libraries
library(ggplot2)
library(tidyr)
library(ggtree)

# Read the data from the CSV file
csv_file_path <- "hsp_logs/drp10_hsp_log_combined.csv"
data <- read.csv(csv_file_path, header = FALSE, stringsAsFactors = FALSE)
colnames(data) <- c("ID", "Species", "Start", "End")

# Extract unique species from the tree file
tree_file_path <- "searched_species.treefile"
tree <- read.tree(tree_file_path)
all_species <- tree$tip.label

# Create a data frame with all species and their ranges
all_data <- data.frame(Species = all_species, Start = NA, End = NA)
range_data <- bind_rows(list(all_data, data))

# Order the species based on their appearance in the tree
range_data$Species <- factor(range_data$Species, levels = tree$tip.label)

# Plot using ggplot
ggplot(range_data, aes(xmin = Start, xmax = End, ymin = Species, ymax = Species)) +
  geom_rect(fill = "skyblue", color = "skyblue", size = 1.5) +  # Adjust the size as needed
  scale_y_discrete(limits = rev(all_species), name = "Species Searched") +
  labs(x = "Range", title = expression(italic("drp10")~"query coverage")) +
  theme_minimal()
