library(ggplot2)
library(tidyr)
library(ggtree)
library(dplyr)

# Read the data from CSV files
drp10_data <- read.csv("hsp_logs/drp10_hsp_log_combined.csv", header = FALSE, col.names = c("ID", "Species", "Start", "End"))
fr10_data <- read.csv("hsp_logs/fr10_hsp_log_combined.csv", header = FALSE, col.names = c("ID", "Species", "Start", "End"))

# Combine both datasets
combined_data <- bind_rows(
  mutate(fr10_data, Query = "fr10"),
  mutate(drp10_data, Query = "drp10")
)


# Read species order
tree_order <- read.table("plots/tree_order.txt", header = FALSE, col.names = "Species")
all_species <- tree_order$Species

# Create a data frame with all species and their ranges
all_data <- data.frame(
  Species = rep(all_species, 2),
  Start = -1, # Replace with a value that won't be plotted
  End = -1,   # Replace with a value that won't be plotted
  Query = rep(c("fr10", "drp10"), each = length(all_species))
)

range_data <- bind_rows(all_data, combined_data)

# Ensure Dataset is a factor with correct levels
range_data$Query <- factor(range_data$Query, levels = c("drp10", "fr10"))

# Load and add overlap data
overlap_data <- read.csv("plots/coverage_overlap.csv", header = FALSE, col.names = c("Species", "Start", "End", "Query", "ID" ))
range_data <- bind_rows(range_data, overlap_data)

# Plot using ggplot
ggplot(range_data, aes(xmin = Start, xmax = End, ymin = Species, ymax = Species, fill = Query, color = Query)) +
  geom_rect(linewidth = 1) +
  scale_y_discrete(limits = all_species) +
  labs(x = "Query Range", title = "Query Coverage") +
  theme_minimal() +
  scale_fill_manual(values = c("drp10" = "brown1", "fr10" = "deepskyblue", "overlap" = "darkorchid")) +
  scale_color_manual(values = c("drp10" = "brown1", "fr10" = "deepskyblue", "overlap" = "darkorchid")) + 
  theme(text = element_text(size = 20))
