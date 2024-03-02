library(tidyverse)
library(RColorBrewer)

# Read csv
df <- read.csv("psi-blast_annotations.csv")

# Find the top 6 annotations by frequency
top_annotations <- head(names(sort(table(df$annotation), decreasing = TRUE)), 5)

# Filter the data frame to include only the top annotations
df_top <- df %>%
  filter(annotation %in% top_annotations) %>%
  group_by(annotation) %>%
  summarize(frequency = n())




# Create a data frame for other annotations
df_other <- df %>%
  filter(!(annotation %in% top_annotations)) %>%
  summarize(frequency = n()) %>%
  mutate(annotation = "Other")

# Combine top annotations and "Other" into a single data frame
df_combined <- bind_rows(df_top, df_other)

# Generate a color palette for each annotation
annotation_colors <- brewer.pal(length(unique(df_combined$annotation)), "Set3")

# Create a pie chart with annotations and pointers
ggplot(df_combined, aes(x = "", y = frequency, fill = annotation)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(label = frequency), position = position_stack(vjust = 0.5), color = "white") +
  scale_fill_manual(values = setNames(annotation_colors, levels(df_combined$annotation))) +
  labs(title = "Blastp hit Annotations", 
       fill = "Annotations") +
  theme_minimal() +
  theme(legend.position = "right", axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())

