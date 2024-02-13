# Read the data from CSV files
drp10_data <- read.csv("hsp_logs/drp10_hsp_log_combined.csv", header = FALSE, col.names = c("ID", "Species", "Start", "End"))
fr10_data <- read.csv("hsp_logs/fr10_hsp_log_combined.csv", header = FALSE, col.names = c("ID", "Species", "Start", "End"))

# Combine both datasets
combined_data <- bind_rows(
  mutate(fr10_data, Dataset = "fr10"),
  mutate(drp10_data, Dataset = "drp10")
)

# Extract unique species from the tree file
tree_file_path <- "searched_species.treefile"
tree <- read.tree(tree_file_path)
all_species <- tree$tip.label

# Create a data frame with all species and their ranges
all_data <- data.frame(
  Species = rep(all_species, 2),
  Start = -1, # Replace with a value that won't be plotted
  End = -1,   # Replace with a value that won't be plotted
  Dataset = rep(c("fr10", "drp10"), each = length(all_species))
)

range_data <- bind_rows(all_data, combined_data)

# Ensure Dataset is a factor with correct levels
range_data$Dataset <- factor(range_data$Dataset, levels = c("drp10", "fr10"))

# Plot using ggplot
ggplot(range_data, aes(xmin = Start, xmax = End, ymin = Species, ymax = Species, fill = Dataset, color = Dataset)) +
  geom_rect(linewidth = 1) +
  scale_y_discrete(limits = rev(all_species), name = "Species Searched") +
  labs(x = "Range", title = "Query Coverage Comparison") +
  theme_minimal() +
  scale_fill_manual(values = c("drp10" = "orange", "fr10" = "skyblue")) +
  scale_color_manual(values = c("drp10" = "orange", "fr10" = "skyblue"))
